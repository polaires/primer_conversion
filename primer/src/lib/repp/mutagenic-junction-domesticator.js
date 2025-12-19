/**
 * Mutagenic Junction Domestication for Golden Gate Assembly
 *
 * This module implements a hybrid domestication strategy that combines:
 * 1. Junction-based splitting (for sites in the middle of fragments)
 * 2. Silent mutations (to prevent site recreation after assembly)
 *
 * KEY INSIGHT: When we place a junction near an internal site, the Golden Gate
 * primers' homology regions will overlap with the recognition sequence. By
 * designing these primers with silent mutations, we:
 * - Introduce the mutation during the same PCR that generates fragments
 * - Choose overhangs freely for optimal fidelity (not constrained by original)
 * - Permanently remove the site (no re-cutting in one-pot reactions)
 *
 * This is the PREFERRED strategy for internal sites in the middle of fragments
 * where direct primer-based mutagenesis cannot reach.
 */

import { GOLDEN_GATE_ENZYMES, findInternalSites, calculateExperimentalFidelity } from './goldengate.js';
import { reverseComplement } from './enzymes.js';
import { CODON_TO_AA, CODON_TABLE, ECOLI_CODON_USAGE, YEAST_CODON_USAGE } from './silent-mutation-domesticator.js';
import { scanForFusionSites } from './fusion-site-scanner.js';

// ============================================================================
// CONFIGURATION
// ============================================================================

export const MUTAGENIC_JUNCTION_CONFIG = {
  // Primer design parameters
  homologyLength: 20,           // Length of primer homology region
  minHomologyLength: 15,        // Minimum acceptable homology
  maxHomologyLength: 30,        // Maximum homology length

  // Junction placement
  searchRadius: 30,             // Search this many bp around the internal site
  minDistanceFromSiteEdge: 0,   // Can place junction within the site

  // Overhang selection
  preferHighFidelityOverhangs: true,
  minOverhangScore: 50,         // Minimum acceptable overhang quality

  // Mutation preferences
  preferWobblePosition: true,   // Prefer 3rd codon position mutations
  avoidRareCodons: true,        // Penalize rare codons
  rareCodoThreshold: 10.0,      // Below this frequency is "rare"
};

// ============================================================================
// MAIN FUNCTIONS
// ============================================================================

/**
 * Design mutagenic junction to domesticate an internal site
 *
 * This finds the optimal junction position and primer mutations to:
 * 1. Split the fragment at the internal site
 * 2. Introduce silent mutations via primer homology regions
 * 3. Choose high-fidelity overhangs (not constrained by original sequence)
 *
 * @param {string} sequence - Full DNA sequence
 * @param {Object} site - Internal site object from findInternalSites
 * @param {string} enzyme - Enzyme name (default: 'BsaI')
 * @param {Object} options - Configuration options
 * @returns {Object} Mutagenic junction design
 */
export function designMutagenicJunction(sequence, site, enzyme = 'BsaI', options = {}) {
  const {
    frame = 0,
    organism = 'ecoli',
    homologyLength = MUTAGENIC_JUNCTION_CONFIG.homologyLength,
    searchRadius = MUTAGENIC_JUNCTION_CONFIG.searchRadius,
    existingOverhangs = [],
  } = options;

  const seq = sequence.toUpperCase();
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enz) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  const recognition = enz.recognition;
  const overhangLen = enz.overhangLength || 4;
  const codonUsage = organism === 'yeast' ? YEAST_CODON_USAGE : ECOLI_CODON_USAGE;

  const siteStart = site.position;
  const siteEnd = siteStart + recognition.length;

  // Find all candidate junction positions around the site
  const candidates = [];

  // Search range: from before the site to after the site
  const searchStart = Math.max(0, siteStart - searchRadius);
  const searchEnd = Math.min(seq.length - overhangLen, siteEnd + searchRadius);

  for (let junctionPos = searchStart; junctionPos <= searchEnd; junctionPos++) {
    // The overhang will be at sequence[junctionPos : junctionPos + overhangLen]
    const overhang = seq.slice(junctionPos, junctionPos + overhangLen);

    // Determine which parts of the recognition site fall in each fragment's primer region
    const fragment1End = junctionPos + overhangLen; // Fragment 1 ends here (includes overhang)
    const fragment2Start = junctionPos;             // Fragment 2 starts here (includes overhang)

    // Fragment 1's reverse primer homology region covers [fragment1End - homologyLength, fragment1End]
    const frag1HomologyStart = Math.max(0, fragment1End - homologyLength);
    const frag1HomologyEnd = fragment1End;

    // Fragment 2's forward primer homology region covers [fragment2Start, fragment2Start + homologyLength]
    const frag2HomologyStart = fragment2Start;
    const frag2HomologyEnd = Math.min(seq.length, fragment2Start + homologyLength);

    // Check which recognition site bases are covered by each primer
    const siteBasesInFrag1Primer = [];
    const siteBasesInFrag2Primer = [];

    for (let i = 0; i < recognition.length; i++) {
      const basePos = siteStart + i;
      if (basePos >= frag1HomologyStart && basePos < frag1HomologyEnd) {
        siteBasesInFrag1Primer.push({ position: basePos, indexInSite: i });
      }
      if (basePos >= frag2HomologyStart && basePos < frag2HomologyEnd) {
        siteBasesInFrag2Primer.push({ position: basePos, indexInSite: i });
      }
    }

    // We need at least ONE base of the site in primer reach to mutate
    const totalCoverage = new Set([
      ...siteBasesInFrag1Primer.map(b => b.position),
      ...siteBasesInFrag2Primer.map(b => b.position),
    ]).size;

    if (totalCoverage === 0) {
      // No part of the site is reachable by primers at this junction position
      continue;
    }

    // Find possible silent mutations for covered positions
    const mutations = findMutationsForCoveredBases(
      seq,
      site,
      siteBasesInFrag1Primer,
      siteBasesInFrag2Primer,
      frame,
      enzyme,
      codonUsage
    );

    if (mutations.length === 0) {
      // No valid silent mutations found for this junction position
      continue;
    }

    // Score the overhang (can be any sequence now, not constrained!)
    const overhangScore = scoreOverhangQuality(overhang, existingOverhangs, enzyme);

    // Score the mutation options
    const bestMutation = mutations[0]; // Already sorted by score

    candidates.push({
      junctionPosition: junctionPos,
      overhang,
      overhangScore,
      mutations,
      bestMutation,
      coverage: {
        fragment1Primer: siteBasesInFrag1Primer,
        fragment2Primer: siteBasesInFrag2Primer,
        totalBaseCoverage: totalCoverage,
      },
      // Combined score: overhang quality + mutation quality
      combinedScore: overhangScore.score * 0.6 + bestMutation.score * 0.4,
    });
  }

  // Sort by combined score
  candidates.sort((a, b) => b.combinedScore - a.combinedScore);

  if (candidates.length === 0) {
    return {
      success: false,
      site,
      error: 'NO_VALID_JUNCTION',
      message: `Could not find junction position with reachable silent mutation for site at ${siteStart}`,
    };
  }

  const best = candidates[0];

  // Design the actual primers with mutations
  const primerDesign = designMutagenicPrimers(
    seq,
    best.junctionPosition,
    best.overhang,
    best.bestMutation,
    enzyme,
    homologyLength
  );

  return {
    success: true,
    site,
    junctionPosition: best.junctionPosition,
    overhang: best.overhang,
    overhangScore: best.overhangScore,
    mutation: best.bestMutation,
    alternativeMutations: best.mutations.slice(1, 5),
    primers: primerDesign,
    coverage: best.coverage,
    combinedScore: best.combinedScore,
    candidates: candidates.slice(0, 10), // Top 10 alternatives
    message: `Junction at ${best.junctionPosition} with ${best.bestMutation.originalBase}→${best.bestMutation.newBase} ` +
             `mutation (${best.bestMutation.originalCodon}→${best.bestMutation.newCodon}, ${best.bestMutation.aminoAcid})`,
  };
}

/**
 * Find silent mutations for bases covered by primer homology regions
 */
function findMutationsForCoveredBases(
  sequence,
  site,
  frag1Coverage,
  frag2Coverage,
  frame,
  enzyme,
  codonUsage
) {
  const mutations = [];
  const recognition = site.sequence;
  const siteStart = site.position;
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  const enzRecognition = enz.recognition;
  const enzRecognitionRC = reverseComplement(enzRecognition);

  // Combine all covered positions
  const allCovered = [...frag1Coverage, ...frag2Coverage];

  // Remove duplicates (bases covered by both primers)
  const uniquePositions = new Map();
  for (const base of allCovered) {
    if (!uniquePositions.has(base.position)) {
      uniquePositions.set(base.position, {
        ...base,
        inFragment1Primer: frag1Coverage.some(b => b.position === base.position),
        inFragment2Primer: frag2Coverage.some(b => b.position === base.position),
      });
    } else {
      const existing = uniquePositions.get(base.position);
      existing.inFragment1Primer = existing.inFragment1Primer || frag1Coverage.some(b => b.position === base.position);
      existing.inFragment2Primer = existing.inFragment2Primer || frag2Coverage.some(b => b.position === base.position);
    }
  }

  // For each covered position, find silent mutations
  for (const [seqPos, baseInfo] of uniquePositions) {
    const originalBase = sequence[seqPos];
    const posInSite = baseInfo.indexInSite;

    // Determine codon context
    const adjustedPos = seqPos - frame;
    if (adjustedPos < 0) continue;

    const codonIndex = Math.floor(adjustedPos / 3);
    const codonStart = codonIndex * 3 + frame;
    const posInCodon = seqPos - codonStart;

    if (codonStart < 0 || codonStart + 3 > sequence.length) continue;
    if (posInCodon < 0 || posInCodon > 2) continue;

    const originalCodon = sequence.slice(codonStart, codonStart + 3);
    const originalAA = CODON_TO_AA[originalCodon];

    if (!originalAA) continue;

    // Try each alternative base
    for (const newBase of ['A', 'T', 'G', 'C']) {
      if (newBase === originalBase) continue;

      // Create new codon
      const newCodon =
        originalCodon.slice(0, posInCodon) +
        newBase +
        originalCodon.slice(posInCodon + 1);

      const newAA = CODON_TO_AA[newCodon];

      // Must be synonymous
      if (newAA !== originalAA) continue;

      // Check if mutation breaks the site
      const newSiteSeq =
        recognition.slice(0, posInSite) +
        newBase +
        recognition.slice(posInSite + 1);

      const breaksForward = newSiteSeq !== enzRecognition;
      const breaksReverse = newSiteSeq !== enzRecognitionRC;
      const siteIsBroken = site.orientation === 'forward' ? breaksForward : breaksReverse;

      if (!siteIsBroken) continue;

      // Score this mutation
      const codonFreq = codonUsage[newCodon] || 0;
      const originalFreq = codonUsage[originalCodon] || 0;

      let score = 80; // Base score

      // Prefer wobble position (3rd codon position)
      if (posInCodon === 2) {
        score += 10;
      }

      // Penalize rare codons
      if (codonFreq < MUTAGENIC_JUNCTION_CONFIG.rareCodoThreshold) {
        score -= 20;
      } else if (codonFreq >= originalFreq * 0.8) {
        score += 5; // Bonus for maintaining frequency
      }

      // Prefer mutations in middle of recognition site (more robust)
      const distFromMiddle = Math.abs(posInSite - (recognition.length - 1) / 2);
      score += (1 - distFromMiddle / (recognition.length / 2)) * 5;

      mutations.push({
        sequencePosition: seqPos,
        positionInSite: posInSite,
        positionInCodon: posInCodon,
        originalBase,
        newBase,
        originalCodon,
        newCodon,
        aminoAcid: originalAA,
        codonFrequency: codonFreq,
        originalCodonFrequency: originalFreq,
        inFragment1Primer: baseInfo.inFragment1Primer,
        inFragment2Primer: baseInfo.inFragment2Primer,
        score,
        isSynonymous: true,
        breaksSite: true,
      });
    }
  }

  // Sort by score
  mutations.sort((a, b) => b.score - a.score);

  return mutations;
}

/**
 * Score overhang quality for the junction
 */
function scoreOverhangQuality(overhang, existingOverhangs, enzyme) {
  let score = 70;
  const issues = [];
  const benefits = [];

  // Check for palindrome
  const rc = reverseComplement(overhang);
  if (overhang === rc) {
    score -= 30;
    issues.push('Palindromic - can self-ligate');
  }

  // Check for homopolymers
  if (/^(.)\1{3}$/.test(overhang)) {
    score -= 25;
    issues.push('Homopolymer - very low fidelity');
  } else if (/(.)\1{2}/.test(overhang)) {
    score -= 10;
    issues.push('Triplet repeat');
  }

  // Check GC content
  const gc = (overhang.match(/[GC]/g) || []).length / overhang.length;
  if (gc === 0.5) {
    score += 10;
    benefits.push('Optimal GC (50%)');
  } else if (gc < 0.25 || gc > 0.75) {
    score -= 10;
    issues.push('Extreme GC content');
  }

  // Check for conflicts with existing overhangs
  for (const existing of existingOverhangs) {
    if (overhang === existing || overhang === reverseComplement(existing)) {
      score -= 40;
      issues.push(`Conflicts with existing overhang ${existing}`);
      break;
    }
  }

  // TNNA pattern bonus
  if (/^T..A$/.test(overhang)) {
    score += 5;
    benefits.push('TNNA pattern - high efficiency');
  }

  return {
    score: Math.max(0, Math.min(100, score)),
    issues,
    benefits,
    gc: gc * 100,
    isPalindrome: overhang === rc,
  };
}

/**
 * Design the actual mutagenic primers
 */
function designMutagenicPrimers(sequence, junctionPos, overhang, mutation, enzyme, homologyLength) {
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  const recognition = enz.recognition;
  const overhangLen = overhang.length;

  // Fragment 1 reverse primer (binds to + strand, synthesizes toward 5' end)
  // Homology region: [junctionPos + overhangLen - homologyLength, junctionPos + overhangLen]
  const frag1HomologyStart = Math.max(0, junctionPos + overhangLen - homologyLength);
  const frag1HomologyEnd = junctionPos + overhangLen;
  let frag1Homology = sequence.slice(frag1HomologyStart, frag1HomologyEnd);

  // Fragment 2 forward primer (binds to - strand, synthesizes toward 3' end)
  // Homology region: [junctionPos, junctionPos + homologyLength]
  const frag2HomologyStart = junctionPos;
  const frag2HomologyEnd = Math.min(sequence.length, junctionPos + homologyLength);
  let frag2Homology = sequence.slice(frag2HomologyStart, frag2HomologyEnd);

  // Apply mutation to the appropriate primer(s)
  const mutPos = mutation.sequencePosition;

  if (mutation.inFragment1Primer && mutPos >= frag1HomologyStart && mutPos < frag1HomologyEnd) {
    const relPos = mutPos - frag1HomologyStart;
    frag1Homology = frag1Homology.slice(0, relPos) + mutation.newBase + frag1Homology.slice(relPos + 1);
  }

  if (mutation.inFragment2Primer && mutPos >= frag2HomologyStart && mutPos < frag2HomologyEnd) {
    const relPos = mutPos - frag2HomologyStart;
    frag2Homology = frag2Homology.slice(0, relPos) + mutation.newBase + frag2Homology.slice(relPos + 1);
  }

  // Build full primers
  // Forward primer (Fragment 2): 5'-[Flanking]-[Recognition]-[N]-[Overhang]-[Homology]-3'
  // Reverse primer (Fragment 1): 5'-[Flanking]-[Recognition]-[N]-[RC(Overhang)]-[RC(Homology)]-3'

  const flanking = 'GCGC'; // Standard flanking sequence
  const spacer = 'N';       // Single N spacer

  const fragment2ForwardPrimer = {
    sequence: flanking + recognition + spacer + overhang + frag2Homology,
    components: {
      flanking,
      recognition,
      spacer,
      overhang,
      homology: frag2Homology,
    },
    hasMutation: mutation.inFragment2Primer,
    mutationPosition: mutation.inFragment2Primer ? (mutPos - frag2HomologyStart) : null,
  };

  const fragment1ReversePrimer = {
    sequence: flanking + recognition + spacer + reverseComplement(overhang) + reverseComplement(frag1Homology),
    components: {
      flanking,
      recognition,
      spacer,
      overhang: reverseComplement(overhang),
      homology: reverseComplement(frag1Homology),
    },
    hasMutation: mutation.inFragment1Primer,
    mutationPosition: mutation.inFragment1Primer ? (frag1HomologyEnd - 1 - mutPos) : null,
  };

  return {
    fragment1: {
      reversePrimer: fragment1ReversePrimer,
      note: 'Reverse primer for Fragment 1 (upstream of junction)',
    },
    fragment2: {
      forwardPrimer: fragment2ForwardPrimer,
      note: 'Forward primer for Fragment 2 (downstream of junction)',
    },
    mutation: {
      applied: true,
      inFragment1Primer: mutation.inFragment1Primer,
      inFragment2Primer: mutation.inFragment2Primer,
      change: `${mutation.originalBase}→${mutation.newBase} at position ${mutPos}`,
      codonChange: `${mutation.originalCodon}→${mutation.newCodon} (${mutation.aminoAcid})`,
    },
  };
}

// ============================================================================
// BATCH PROCESSING
// ============================================================================

/**
 * Design mutagenic junctions for all internal sites in a sequence
 *
 * @param {string} sequence - DNA sequence
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Configuration options
 * @returns {Object} Complete domestication plan
 */
export function designAllMutagenicJunctions(sequence, enzyme = 'BsaI', options = {}) {
  const {
    frame = 0,
    organism = 'ecoli',
    existingOverhangs = [],
  } = options;

  const internalSites = findInternalSites(sequence, enzyme);

  if (!internalSites.hasSites) {
    return {
      success: true,
      needsDomestication: false,
      sites: [],
      junctions: [],
      message: `No internal ${enzyme} sites found`,
    };
  }

  const junctions = [];
  const failedSites = [];
  const allOverhangs = [...existingOverhangs];

  for (const site of internalSites.sites) {
    const result = designMutagenicJunction(sequence, site, enzyme, {
      frame,
      organism,
      existingOverhangs: allOverhangs,
    });

    if (result.success) {
      junctions.push(result);
      allOverhangs.push(result.overhang); // Add to existing for conflict checking
    } else {
      failedSites.push({
        site,
        error: result.error,
        message: result.message,
      });
    }
  }

  // Calculate overall fidelity with all new overhangs
  let fidelityResult = null;
  if (junctions.length > 0) {
    try {
      fidelityResult = calculateExperimentalFidelity(allOverhangs, enzyme);
    } catch (e) {
      // Fidelity calculation may not be available
    }
  }

  return {
    success: failedSites.length === 0,
    needsDomestication: true,
    totalSites: internalSites.count,
    sites: internalSites.sites,
    junctions,
    failedSites,
    allOverhangs,
    fidelity: fidelityResult?.assemblyFidelity,
    additionalFragments: junctions.length,
    message: `Designed ${junctions.length} mutagenic junction(s) for ${internalSites.count} internal site(s)` +
             (failedSites.length > 0 ? `. ${failedSites.length} site(s) failed.` : ''),
  };
}

// ============================================================================
// INTEGRATION WITH UNIFIED OPTIMIZER
// ============================================================================

/**
 * Compare all domestication strategies and recommend the best approach
 *
 * Strategy hierarchy:
 * 1. Direct primer mutation (site within ~25bp of existing junction)
 * 2. Mutagenic junction (site in middle of fragment)
 * 3. Alternative enzyme (if one has no internal sites)
 *
 * @param {string} sequence - DNA sequence
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Configuration options
 * @returns {Object} Strategy comparison and recommendation
 */
export function selectDomesticationStrategy(sequence, enzyme = 'BsaI', options = {}) {
  const {
    frame = 0,
    userJunctionPositions = [], // User's existing fragment boundaries
    requestedFragments = 2,
  } = options;

  const internalSites = findInternalSites(sequence, enzyme);

  if (!internalSites.hasSites) {
    return {
      needsDomestication: false,
      strategy: 'none',
      message: 'No internal sites found',
    };
  }

  // Calculate ideal junction positions for user's requested fragments
  const idealJunctions = calculateIdealJunctions(sequence.length, requestedFragments);
  const effectiveJunctions = userJunctionPositions.length > 0 ? userJunctionPositions : idealJunctions;

  // Analyze each site
  const siteAnalyses = internalSites.sites.map(site => {
    // Check distance to nearest junction
    const nearestJunction = findNearestJunction(site.position, effectiveJunctions);
    const distanceToJunction = nearestJunction ? nearestJunction.distance : Infinity;

    // Can we reach the site with a primer from an existing junction?
    const reachableByPrimer = distanceToJunction <= MUTAGENIC_JUNCTION_CONFIG.homologyLength;

    // Design mutagenic junction (for comparison)
    const mutagenicJunction = designMutagenicJunction(sequence, site, enzyme, options);

    return {
      site,
      distanceToNearestJunction: distanceToJunction,
      nearestJunction,
      strategies: {
        directPrimerMutation: {
          available: reachableByPrimer,
          fragmentIncrease: 0,
          description: 'Mutate via existing junction primer',
        },
        mutagenicJunction: {
          available: mutagenicJunction.success,
          fragmentIncrease: 1,
          junction: mutagenicJunction,
          description: 'Add junction with mutagenic primers',
        },
      },
      recommendedStrategy: reachableByPrimer ? 'directPrimerMutation' : 'mutagenicJunction',
    };
  });

  // Check for alternative enzymes
  const alternativeEnzymes = checkAlternativeEnzymes(sequence, enzyme);

  // Overall recommendation
  const directMutationCount = siteAnalyses.filter(a => a.strategies.directPrimerMutation.available).length;
  const mutagenicJunctionCount = siteAnalyses.filter(a =>
    !a.strategies.directPrimerMutation.available && a.strategies.mutagenicJunction.available
  ).length;
  const failedCount = siteAnalyses.filter(a =>
    !a.strategies.directPrimerMutation.available && !a.strategies.mutagenicJunction.available
  ).length;

  let overallStrategy = 'mutagenic_junction';
  let overallMessage = '';

  if (alternativeEnzymes.hasCompatible && failedCount > 0) {
    overallStrategy = 'alternative_enzyme';
    overallMessage = `Switch to ${alternativeEnzymes.best.enzyme} which has no internal sites`;
  } else if (directMutationCount === internalSites.count) {
    overallStrategy = 'direct_primer_mutation';
    overallMessage = 'All sites reachable by existing junction primers';
  } else if (failedCount === 0) {
    overallStrategy = 'mutagenic_junction';
    overallMessage = `Add ${mutagenicJunctionCount} mutagenic junction(s)`;
  } else {
    overallStrategy = 'hybrid';
    overallMessage = `${directMutationCount} direct, ${mutagenicJunctionCount} junctions, ${failedCount} need alternative enzyme`;
  }

  return {
    needsDomestication: true,
    totalSites: internalSites.count,
    siteAnalyses,
    alternativeEnzymes,
    summary: {
      directPrimerMutation: directMutationCount,
      mutagenicJunction: mutagenicJunctionCount,
      failed: failedCount,
      additionalFragments: mutagenicJunctionCount,
    },
    overallStrategy,
    message: overallMessage,
    onePotCompatible: true, // All our strategies work in one-pot!
  };
}

/**
 * Calculate ideal junction positions for N fragments
 */
function calculateIdealJunctions(seqLength, numFragments) {
  if (numFragments <= 1) return [];

  const junctions = [];
  const fragmentSize = seqLength / numFragments;

  for (let i = 1; i < numFragments; i++) {
    junctions.push(Math.round(fragmentSize * i));
  }

  return junctions;
}

/**
 * Find nearest junction to a position
 */
function findNearestJunction(position, junctions) {
  if (junctions.length === 0) return null;

  let nearest = null;
  let minDistance = Infinity;

  for (const junction of junctions) {
    const distance = Math.abs(position - junction);
    if (distance < minDistance) {
      minDistance = distance;
      nearest = { position: junction, distance };
    }
  }

  return nearest;
}

/**
 * Check for alternative enzymes without internal sites
 */
function checkAlternativeEnzymes(sequence, currentEnzyme) {
  const results = [];

  for (const [enzName, enzData] of Object.entries(GOLDEN_GATE_ENZYMES)) {
    if (enzName === currentEnzyme) continue;

    const sites = findInternalSites(sequence, enzName);
    results.push({
      enzyme: enzName,
      siteCount: sites.count,
      isCompatible: !sites.hasSites,
    });
  }

  results.sort((a, b) => a.siteCount - b.siteCount);

  return {
    alternatives: results,
    hasCompatible: results.some(r => r.isCompatible),
    best: results[0],
  };
}

// Functions are exported inline with their declarations above
