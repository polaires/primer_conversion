/**
 * ORF Detection and Reading Frame Validation for Golden Gate Domestication
 *
 * This module provides automatic detection of Open Reading Frames (ORFs) and
 * validates reading frames before silent mutations are applied. This is CRITICAL
 * for ensuring mutations remain synonymous.
 *
 * KEY SAFETY PRINCIPLE: Never auto-apply mutations without user confirmation
 * of the correct reading frame.
 */

import { CODON_TO_AA, CODON_TABLE, ECOLI_CODON_USAGE, YEAST_CODON_USAGE } from './silent-mutation-domesticator.js';

// ============================================================================
// CONFIGURATION
// ============================================================================

export const ORF_DETECTION_CONFIG = {
  // Minimum ORF length to consider (in amino acids)
  minOrfLength: 30,
  // Maximum number of ORFs to return
  maxOrfsToReturn: 10,
  // Consider alternative start codons (GTG, TTG in bacteria)
  useAlternativeStarts: false,
  // Standard start codons
  standardStartCodons: ['ATG'],
  // Alternative start codons (prokaryotes)
  alternativeStartCodons: ['ATG', 'GTG', 'TTG'],
  // Stop codons
  stopCodons: ['TAA', 'TAG', 'TGA'],
};

// ============================================================================
// ORF DETECTION
// ============================================================================

/**
 * Detect all Open Reading Frames (ORFs) in a sequence
 *
 * An ORF is defined as a sequence starting with a start codon (ATG) and
 * ending with a stop codon (TAA, TAG, TGA), with no internal stop codons.
 *
 * @param {string} sequence - DNA sequence to analyze
 * @param {Object} options - Detection options
 * @returns {Object} Detection result with ORFs and recommendations
 */
export function detectORFs(sequence, options = {}) {
  const {
    minLength = ORF_DETECTION_CONFIG.minOrfLength,
    useAlternativeStarts = ORF_DETECTION_CONFIG.useAlternativeStarts,
    organism = 'ecoli',
  } = options;

  const seq = sequence.toUpperCase();
  const orfs = [];

  const startCodons = useAlternativeStarts
    ? ORF_DETECTION_CONFIG.alternativeStartCodons
    : ORF_DETECTION_CONFIG.standardStartCodons;

  // Check all three reading frames on forward strand
  for (let frame = 0; frame < 3; frame++) {
    const frameOrfs = findOrfsInFrame(seq, frame, startCodons, minLength, organism);
    orfs.push(...frameOrfs);
  }

  // Check all three reading frames on reverse complement
  const rcSeq = reverseComplementSequence(seq);
  for (let frame = 0; frame < 3; frame++) {
    const frameOrfs = findOrfsInFrame(rcSeq, frame, startCodons, minLength, organism);
    // Mark these as reverse strand and adjust positions
    frameOrfs.forEach(orf => {
      orf.strand = 'reverse';
      // Convert positions to forward strand coordinates
      orf.forwardStart = seq.length - orf.end;
      orf.forwardEnd = seq.length - orf.start;
    });
    orfs.push(...frameOrfs);
  }

  // Sort by length (longest first) and score
  orfs.sort((a, b) => {
    if (b.score !== a.score) return b.score - a.score;
    return b.length - a.length;
  });

  // Limit results
  const limitedOrfs = orfs.slice(0, ORF_DETECTION_CONFIG.maxOrfsToReturn);

  // Generate recommendations
  const recommendation = generateOrfRecommendation(limitedOrfs, seq);

  return {
    orfs: limitedOrfs,
    totalFound: orfs.length,
    sequenceLength: seq.length,
    recommendation,
    hasOrfs: orfs.length > 0,
    bestOrf: limitedOrfs[0] || null,
  };
}

/**
 * Find ORFs in a specific reading frame
 */
function findOrfsInFrame(sequence, frame, startCodons, minLength, organism) {
  const orfs = [];
  const codonUsage = organism === 'yeast' ? YEAST_CODON_USAGE : ECOLI_CODON_USAGE;
  const stopCodons = ORF_DETECTION_CONFIG.stopCodons;

  // Scan for start codons
  for (let i = frame; i + 3 <= sequence.length; i += 3) {
    const codon = sequence.slice(i, i + 3);

    if (startCodons.includes(codon)) {
      // Found a start codon, look for stop codon
      const orf = extendOrf(sequence, i, stopCodons, codonUsage);

      if (orf && orf.proteinLength >= minLength) {
        orf.frame = frame;
        orf.strand = 'forward';
        orfs.push(orf);
      }
    }
  }

  return orfs;
}

/**
 * Extend from a start codon to find the complete ORF
 */
function extendOrf(sequence, startPos, stopCodons, codonUsage) {
  const codons = [];
  const protein = [];
  let codonUsageScore = 0;
  let hasRareCodon = false;

  for (let i = startPos; i + 3 <= sequence.length; i += 3) {
    const codon = sequence.slice(i, i + 3);
    const aa = CODON_TO_AA[codon];

    if (!aa) {
      // Invalid codon (contains N or other non-standard bases)
      return null;
    }

    codons.push(codon);
    protein.push(aa);

    // Track codon usage
    const usage = codonUsage[codon] || 0;
    codonUsageScore += usage;
    if (usage < 10) hasRareCodon = true;

    // Check for stop codon
    if (stopCodons.includes(codon)) {
      const endPos = i + 3;
      const proteinLength = protein.length - 1; // Exclude stop codon from count
      const avgCodonUsage = codonUsageScore / codons.length;

      // Calculate ORF score based on:
      // - Length (longer is better)
      // - Codon usage (higher average is better)
      // - Presence of rare codons (penalty)
      let score = 0;
      score += Math.min(proteinLength / 100, 1) * 40; // Length component (max 40)
      score += Math.min(avgCodonUsage / 30, 1) * 40;  // Codon usage component (max 40)
      if (hasRareCodon) score -= 10;
      score = Math.max(0, Math.min(100, score));

      return {
        start: startPos,
        end: endPos,
        length: endPos - startPos,
        proteinLength,
        proteinSequence: protein.slice(0, -1).join(''), // Exclude stop
        stopCodon: codon,
        codons,
        avgCodonUsage: avgCodonUsage.toFixed(1),
        hasRareCodon,
        score,
      };
    }
  }

  // No stop codon found - incomplete ORF
  // Still return it but mark as incomplete
  const proteinLength = protein.length;
  if (proteinLength >= 10) { // Only return if substantial
    const avgCodonUsage = codonUsageScore / codons.length;
    return {
      start: startPos,
      end: sequence.length,
      length: sequence.length - startPos,
      proteinLength,
      proteinSequence: protein.join(''),
      stopCodon: null,
      isIncomplete: true,
      codons,
      avgCodonUsage: avgCodonUsage.toFixed(1),
      hasRareCodon,
      score: 20, // Lower score for incomplete ORFs
    };
  }

  return null;
}

/**
 * Generate recommendation for ORF selection
 */
function generateOrfRecommendation(orfs, sequence) {
  if (orfs.length === 0) {
    return {
      type: 'NO_ORF',
      message: 'No significant ORFs detected. This may be a non-coding sequence.',
      action: 'Consider using junction-based domestication or alternative enzyme.',
      confidence: 'low',
    };
  }

  if (orfs.length === 1) {
    const orf = orfs[0];
    if (orf.start === 0 && orf.end === sequence.length) {
      return {
        type: 'FULL_SEQUENCE_ORF',
        message: `Single ORF covers entire sequence (frame ${orf.frame}, ${orf.proteinLength} aa).`,
        action: 'Use frame ' + orf.frame + ' for domestication.',
        recommendedFrame: orf.frame,
        confidence: 'high',
      };
    }
    return {
      type: 'SINGLE_ORF',
      message: `Single ORF detected at positions ${orf.start}-${orf.end} (frame ${orf.frame}, ${orf.proteinLength} aa).`,
      action: 'Confirm this is the coding region before domestication.',
      recommendedFrame: orf.frame,
      confidence: 'medium',
    };
  }

  // Multiple ORFs - more complex situation
  const bestOrf = orfs[0];
  const secondBest = orfs[1];

  // Check if best ORF is significantly better
  if (bestOrf.score > secondBest.score + 20) {
    return {
      type: 'CLEAR_BEST',
      message: `Best ORF: frame ${bestOrf.frame}, ${bestOrf.proteinLength} aa (score: ${bestOrf.score.toFixed(0)}). ${orfs.length - 1} other ORFs detected.`,
      action: 'Recommended frame ' + bestOrf.frame + ', but please verify.',
      recommendedFrame: bestOrf.frame,
      confidence: 'medium',
    };
  }

  return {
    type: 'AMBIGUOUS',
    message: `Multiple similar ORFs detected. Best: frame ${bestOrf.frame} (${bestOrf.proteinLength} aa), Second: frame ${secondBest.frame} (${secondBest.proteinLength} aa).`,
    action: 'USER MUST SELECT the correct reading frame.',
    recommendedFrame: null,
    confidence: 'low',
  };
}

// ============================================================================
// READING FRAME VALIDATION
// ============================================================================

/**
 * Validate that a specified reading frame produces a valid protein
 *
 * @param {string} sequence - DNA sequence
 * @param {number} frame - Reading frame (0, 1, or 2)
 * @param {Object} options - Validation options
 * @returns {Object} Validation result
 */
export function validateReadingFrame(sequence, frame, options = {}) {
  const {
    requireStartCodon = false,
    requireStopCodon = false,
    maxInternalStops = 0,
  } = options;

  const seq = sequence.toUpperCase();
  const validations = [];
  let isValid = true;

  // Check frame is valid
  if (frame < 0 || frame > 2) {
    return {
      isValid: false,
      validations: [{ check: 'VALID_FRAME', passed: false, message: 'Frame must be 0, 1, or 2' }],
    };
  }

  // Translate the sequence
  const translation = translateSequence(seq, frame);

  // Check for start codon
  if (requireStartCodon) {
    const firstCodon = seq.slice(frame, frame + 3);
    const hasStart = firstCodon === 'ATG';
    validations.push({
      check: 'START_CODON',
      passed: hasStart,
      message: hasStart ? 'Starts with ATG (Methionine)' : `First codon is ${firstCodon}, not ATG`,
    });
    if (!hasStart) isValid = false;
  }

  // Check for internal stop codons
  const internalStops = [];
  for (let i = 0; i < translation.protein.length - 1; i++) {
    if (translation.protein[i] === '*') {
      internalStops.push({
        position: i,
        codonPosition: frame + (i * 3),
        codon: translation.codons[i],
      });
    }
  }

  validations.push({
    check: 'INTERNAL_STOPS',
    passed: internalStops.length <= maxInternalStops,
    message: internalStops.length === 0
      ? 'No internal stop codons'
      : `${internalStops.length} internal stop codon(s) at positions: ${internalStops.map(s => s.position).join(', ')}`,
    details: internalStops,
  });

  if (internalStops.length > maxInternalStops) {
    isValid = false;
  }

  // Check for stop codon at end
  if (requireStopCodon) {
    const lastAA = translation.protein[translation.protein.length - 1];
    const hasStop = lastAA === '*';
    validations.push({
      check: 'STOP_CODON',
      passed: hasStop,
      message: hasStop ? 'Ends with stop codon' : 'No terminal stop codon',
    });
    if (!hasStop) isValid = false;
  }

  // Check for unknown codons
  const unknownCodons = translation.codons.filter(c => CODON_TO_AA[c] === undefined);
  if (unknownCodons.length > 0) {
    validations.push({
      check: 'VALID_CODONS',
      passed: false,
      message: `Contains ${unknownCodons.length} invalid codon(s): ${unknownCodons.join(', ')}`,
    });
    isValid = false;
  } else {
    validations.push({
      check: 'VALID_CODONS',
      passed: true,
      message: 'All codons are valid',
    });
  }

  return {
    isValid,
    frame,
    validations,
    translation,
    proteinLength: translation.protein.replace(/\*/g, '').length,
    internalStops,
  };
}

/**
 * Translate a DNA sequence to protein
 *
 * @param {string} sequence - DNA sequence
 * @param {number} frame - Reading frame (0, 1, or 2)
 * @returns {Object} Translation result with protein and codons
 */
export function translateSequence(sequence, frame = 0) {
  const seq = sequence.toUpperCase();
  const codons = [];
  const protein = [];
  const codonInfo = [];

  for (let i = frame; i + 3 <= seq.length; i += 3) {
    const codon = seq.slice(i, i + 3);
    const aa = CODON_TO_AA[codon] || '?';

    codons.push(codon);
    protein.push(aa);
    codonInfo.push({
      position: i,
      codon,
      aminoAcid: aa,
      isStop: aa === '*',
      isStart: codon === 'ATG',
    });
  }

  return {
    protein: protein.join(''),
    codons,
    codonInfo,
    frame,
    nucleotideLength: seq.length,
    codonCount: codons.length,
  };
}

// ============================================================================
// SITE CONTEXT ANALYSIS
// ============================================================================

/**
 * Analyze the codon context of an internal restriction site
 *
 * This determines which codons overlap with the site and what
 * synonymous mutation options are available.
 *
 * @param {string} sequence - DNA sequence
 * @param {Object} site - Site object with position and sequence
 * @param {number} frame - Reading frame
 * @returns {Object} Context analysis with mutation options
 */
export function analyzeSiteCodonContext(sequence, site, frame) {
  const seq = sequence.toUpperCase();
  const siteStart = site.position;
  const siteSeq = site.sequence;
  const siteEnd = siteStart + siteSeq.length;

  const overlappingCodons = [];
  const mutationOptions = [];

  // Find all codons that overlap with the site
  for (let codonStart = frame; codonStart + 3 <= seq.length; codonStart += 3) {
    const codonEnd = codonStart + 3;

    // Check if codon overlaps with site
    if (codonEnd > siteStart && codonStart < siteEnd) {
      const codon = seq.slice(codonStart, codonEnd);
      const aa = CODON_TO_AA[codon];

      // Determine which positions of the site are in this codon
      const sitePositionsInCodon = [];
      for (let i = 0; i < siteSeq.length; i++) {
        const sitePos = siteStart + i;
        if (sitePos >= codonStart && sitePos < codonEnd) {
          sitePositionsInCodon.push({
            siteIndex: i,
            codonPosition: sitePos - codonStart,
            base: siteSeq[i],
          });
        }
      }

      const codonInfo = {
        codonStart,
        codon,
        aminoAcid: aa,
        sitePositionsInCodon,
        hasSynonymousOptions: false,
        synonymousCodons: [],
      };

      // Find synonymous alternatives
      if (aa && CODON_TABLE[aa]) {
        const alternatives = CODON_TABLE[aa].filter(c => c !== codon);
        codonInfo.synonymousCodons = alternatives;
        codonInfo.hasSynonymousOptions = alternatives.length > 0;

        // Check which alternatives would break the site
        for (const altCodon of alternatives) {
          // Would this alternative break the site?
          let wouldBreakSite = false;
          for (const pos of sitePositionsInCodon) {
            if (altCodon[pos.codonPosition] !== codon[pos.codonPosition]) {
              wouldBreakSite = true;
              break;
            }
          }

          if (wouldBreakSite) {
            mutationOptions.push({
              codonStart,
              originalCodon: codon,
              newCodon: altCodon,
              aminoAcid: aa,
              changedPositions: sitePositionsInCodon.filter(
                pos => altCodon[pos.codonPosition] !== codon[pos.codonPosition]
              ),
            });
          }
        }
      }

      overlappingCodons.push(codonInfo);
    }
  }

  // Analyze mutation viability
  const canUseSilentMutation = mutationOptions.length > 0;
  const hasSingleCodonCoverage = overlappingCodons.length === 1;
  const allCodonsHaveAlternatives = overlappingCodons.every(c => c.hasSynonymousOptions);

  // Check for amino acids with no alternatives (Met, Trp)
  const immutableCodons = overlappingCodons.filter(c =>
    c.aminoAcid === 'M' || c.aminoAcid === 'W'
  );

  return {
    site,
    frame,
    overlappingCodons,
    mutationOptions,
    analysis: {
      canUseSilentMutation,
      hasSingleCodonCoverage,
      allCodonsHaveAlternatives,
      immutableCodons,
      totalOptions: mutationOptions.length,
    },
    recommendation: canUseSilentMutation
      ? `${mutationOptions.length} silent mutation option(s) available`
      : immutableCodons.length > 0
        ? `Site overlaps with ${immutableCodons.map(c => c.aminoAcid).join('/')} codon(s) - no synonymous alternatives`
        : 'No silent mutation options available',
  };
}

// ============================================================================
// PRE-FLIGHT VALIDATION
// ============================================================================

/**
 * Comprehensive pre-flight check before domestication
 *
 * This is the main entry point for validating that domestication
 * can proceed safely. It combines ORF detection, frame validation,
 * and site analysis.
 *
 * @param {string} sequence - DNA sequence
 * @param {Array} sites - Internal sites to domesticate
 * @param {Object} options - Validation options
 * @returns {Object} Pre-flight check result
 */
export function preFlightCheck(sequence, sites, options = {}) {
  const {
    frame = null, // If null, will try to detect
    organism = 'ecoli',
    requireUserConfirmation = true,
  } = options;

  const checks = [];
  const warnings = [];
  const errors = [];
  let recommendedFrame = frame;

  // Step 1: ORF Detection (if frame not specified)
  let orfDetection = null;
  if (frame === null) {
    orfDetection = detectORFs(sequence, { organism });

    if (orfDetection.recommendation.confidence === 'high') {
      recommendedFrame = orfDetection.recommendation.recommendedFrame;
      checks.push({
        name: 'ORF_DETECTION',
        status: 'passed',
        message: `Detected ORF with high confidence (frame ${recommendedFrame})`,
        requiresConfirmation: requireUserConfirmation,
      });
    } else if (orfDetection.recommendation.confidence === 'medium') {
      recommendedFrame = orfDetection.recommendation.recommendedFrame;
      warnings.push({
        name: 'ORF_DETECTION',
        message: 'ORF detected with medium confidence - user confirmation required',
        recommendedFrame,
      });
    } else {
      errors.push({
        name: 'ORF_DETECTION',
        message: 'Cannot determine reading frame automatically',
        details: orfDetection.recommendation,
      });
    }
  } else {
    recommendedFrame = frame;
    checks.push({
      name: 'FRAME_SPECIFIED',
      status: 'passed',
      message: `Using user-specified frame ${frame}`,
    });
  }

  // Step 2: Frame Validation (if we have a frame)
  let frameValidation = null;
  if (recommendedFrame !== null) {
    frameValidation = validateReadingFrame(sequence, recommendedFrame);

    if (frameValidation.isValid) {
      checks.push({
        name: 'FRAME_VALIDATION',
        status: 'passed',
        message: `Frame ${recommendedFrame} produces valid protein (${frameValidation.proteinLength} aa)`,
      });
    } else {
      const failedChecks = frameValidation.validations.filter(v => !v.passed);
      warnings.push({
        name: 'FRAME_VALIDATION',
        message: `Frame ${recommendedFrame} has issues: ${failedChecks.map(c => c.message).join('; ')}`,
        details: frameValidation,
      });
    }
  }

  // Step 3: Site Analysis (for each site)
  const siteAnalyses = [];
  let allSitesCanBeDomesticated = true;

  for (const site of sites) {
    if (recommendedFrame !== null) {
      const analysis = analyzeSiteCodonContext(sequence, site, recommendedFrame);
      siteAnalyses.push(analysis);

      if (!analysis.analysis.canUseSilentMutation) {
        allSitesCanBeDomesticated = false;
        warnings.push({
          name: 'SITE_ANALYSIS',
          message: `Site at ${site.position} cannot use silent mutation: ${analysis.recommendation}`,
          site,
          analysis,
        });
      }
    }
  }

  if (siteAnalyses.length > 0 && allSitesCanBeDomesticated) {
    checks.push({
      name: 'SITE_ANALYSIS',
      status: 'passed',
      message: `All ${sites.length} site(s) can be domesticated with silent mutations`,
    });
  }

  // Generate preview data for UI
  const preview = generateDomesticationPreview(sequence, sites, recommendedFrame, siteAnalyses);

  // Determine overall status
  const status = errors.length > 0
    ? 'error'
    : warnings.length > 0
      ? 'warning'
      : 'ready';

  const needsUserConfirmation = requireUserConfirmation ||
    (orfDetection?.recommendation.confidence !== 'high') ||
    warnings.length > 0;

  return {
    status,
    readyToProceed: status !== 'error' && (!needsUserConfirmation || frame !== null),
    needsUserConfirmation,
    recommendedFrame,
    checks,
    warnings,
    errors,
    orfDetection,
    frameValidation,
    siteAnalyses,
    preview,
    message: status === 'error'
      ? `Cannot proceed: ${errors.map(e => e.message).join('; ')}`
      : status === 'warning'
        ? `Review required: ${warnings.map(w => w.message).join('; ')}`
        : 'Ready for domestication',
  };
}

/**
 * Generate a preview of the domestication for user review
 */
function generateDomesticationPreview(sequence, sites, frame, siteAnalyses) {
  if (frame === null || sites.length === 0) {
    return null;
  }

  const seq = sequence.toUpperCase();
  const translation = translateSequence(seq, frame);

  // Build a visual representation
  const preview = {
    originalSequence: seq,
    frame,
    protein: translation.protein,
    sites: sites.map((site, i) => {
      const analysis = siteAnalyses[i];
      return {
        position: site.position,
        sequence: site.sequence,
        orientation: site.orientation,
        overlappingCodons: analysis?.overlappingCodons || [],
        bestMutation: analysis?.mutationOptions[0] || null,
        alternativeMutations: analysis?.mutationOptions.slice(1, 4) || [],
      };
    }),
    proteinLength: translation.protein.replace(/\*/g, '').length,
  };

  return preview;
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Reverse complement a DNA sequence
 */
function reverseComplementSequence(sequence) {
  const complement = {
    A: 'T', T: 'A', G: 'C', C: 'G',
    a: 't', t: 'a', g: 'c', c: 'g',
    N: 'N', n: 'n',
  };

  return sequence
    .split('')
    .reverse()
    .map(base => complement[base] || base)
    .join('');
}

/**
 * Format protein sequence for display (with position markers)
 */
export function formatProteinForDisplay(protein, lineLength = 60) {
  const lines = [];
  for (let i = 0; i < protein.length; i += lineLength) {
    const chunk = protein.slice(i, i + lineLength);
    lines.push({
      start: i + 1,
      end: Math.min(i + lineLength, protein.length),
      sequence: chunk,
    });
  }
  return lines;
}

/**
 * Compare two protein sequences and find differences
 */
export function compareProteins(original, domesticated) {
  const differences = [];
  const maxLen = Math.max(original.length, domesticated.length);

  for (let i = 0; i < maxLen; i++) {
    const origAA = original[i] || '-';
    const domAA = domesticated[i] || '-';

    if (origAA !== domAA) {
      differences.push({
        position: i + 1,
        original: origAA,
        domesticated: domAA,
        type: origAA === '-' ? 'insertion' : domAA === '-' ? 'deletion' : 'substitution',
      });
    }
  }

  return {
    identical: differences.length === 0,
    differences,
    originalLength: original.length,
    domesticatedLength: domesticated.length,
    identity: ((maxLen - differences.length) / maxLen * 100).toFixed(1),
  };
}
