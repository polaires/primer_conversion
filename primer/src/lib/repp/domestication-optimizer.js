/**
 * Unified Domestication Optimizer for Golden Gate Assembly
 *
 * This module provides a unified interface for domesticating sequences with
 * internal restriction sites. ALL strategies are one-pot compatible!
 *
 * Strategies ranked by preference:
 * 1. Direct primer mutation - Site within primer reach of existing junction (0 extra fragments)
 * 2. Mutagenic junction - Split at site with mutagenic primers (adds fragments, but one-pot safe)
 * 3. Alternative enzyme - If available enzyme has no internal sites (0 extra fragments)
 *
 * DEPRECATED: The old junction-based approach (without mutations) is NOT recommended
 * as it recreates the internal site in the assembled product, which gets re-cut in
 * one-pot reactions.
 */

import {
  domesticateWithSilentMutations,
  findAllSilentMutationCandidates,
  scoreMutationCandidates,
  verifyProteinSequence,
  validateDomestication,
  compareDomesticationStrategies,
  applyMutation,
  applyMutations,
  CODON_TABLE,
  CODON_TO_AA,
  ECOLI_CODON_USAGE,
  YEAST_CODON_USAGE,
} from './silent-mutation-domesticator.js';

import {
  designMutagenicJunction,
  designAllMutagenicJunctions,
  selectDomesticationStrategy,
  MUTAGENIC_JUNCTION_CONFIG,
} from './mutagenic-junction-domesticator.js';

import {
  analyzeForDomestication as analyzeForJunctionDomestication,
  findDomesticationJunctions,
  recommendAlternativeEnzymes,
  optimizeGlobalOverhangSet,
  validatePostDomestication,
  DOMESTICATION_DEFAULTS,
} from './auto-domestication-optimizer.js';

import { findInternalSites, GOLDEN_GATE_ENZYMES } from './goldengate.js';

// ============================================================================
// CONFIGURATION
// ============================================================================

export const DOMESTICATION_STRATEGY = {
  DIRECT_PRIMER_MUTATION: 'direct_primer_mutation', // Site near existing junction
  MUTAGENIC_JUNCTION: 'mutagenic_junction',         // Split with mutagenic primers
  ALTERNATIVE_ENZYME: 'alternative_enzyme',          // Switch enzyme
  LEGACY_JUNCTION: 'legacy_junction',               // Old method (NOT recommended)
  HYBRID: 'hybrid',
  NONE: 'none',
};

export const UNIFIED_DOMESTICATION_CONFIG = {
  // Strategy preferences (higher = more preferred)
  strategyPriority: {
    [DOMESTICATION_STRATEGY.DIRECT_PRIMER_MUTATION]: 100,  // Best: no extra fragments
    [DOMESTICATION_STRATEGY.MUTAGENIC_JUNCTION]: 90,       // Good: one-pot compatible
    [DOMESTICATION_STRATEGY.ALTERNATIVE_ENZYME]: 80,       // Good: no changes needed
    [DOMESTICATION_STRATEGY.LEGACY_JUNCTION]: 10,          // Bad: NOT one-pot compatible
  },

  // Warnings
  warnings: {
    legacyJunction: {
      severity: 'high',
      message: 'Legacy junction-based domestication recreates internal sites in assembled product. ' +
               'NOT compatible with one-pot Golden Gate reactions. ' +
               'Use mutagenic junction strategy instead.',
    },
    rareCodon: {
      severity: 'low',
      message: 'Silent mutation introduces a rare codon. Assembly will work but expression may be affected.',
    },
    newSiteCreated: {
      severity: 'medium',
      message: 'Silent mutation creates a new restriction site for a different enzyme.',
    },
  },
};

// ============================================================================
// MAIN UNIFIED DOMESTICATION FUNCTION
// ============================================================================

/**
 * Optimize domestication using the best available strategy
 *
 * Priority order:
 * 1. Silent mutations (preserves fragment count, one-pot compatible)
 * 2. Alternative enzyme (if one exists without internal sites)
 * 3. Junction-based (fallback, with warnings)
 *
 * @param {string} sequence - DNA sequence to domesticate
 * @param {string} enzyme - Enzyme name (default: 'BsaI')
 * @param {Object} options - Configuration options
 * @returns {Object} Comprehensive domestication result
 */
export function optimizeDomestication(sequence, enzyme = 'BsaI', options = {}) {
  const {
    frame = 0,
    isCodingSequence = true,
    organism = 'ecoli',
    allowJunctionFallback = true,
    preferAlternativeEnzyme = false,
    existingOverhangs = [],
    requestedFragments = null,
    checkAllEnzymes = true,
  } = options;

  const seq = sequence.toUpperCase();
  const result = {
    originalSequence: seq,
    domesticatedSequence: null,
    enzyme,
    strategy: null,
    strategies: {},
    recommendations: [],
    warnings: [],
    success: false,
  };

  // Step 1: Check if domestication is needed
  const internalSites = findInternalSites(seq, enzyme);

  if (!internalSites.hasSites) {
    return {
      ...result,
      domesticatedSequence: seq,
      strategy: DOMESTICATION_STRATEGY.NONE,
      success: true,
      message: `No internal ${enzyme} sites found - sequence is already compatible`,
    };
  }

  result.internalSites = internalSites;

  // Step 2: Try silent mutations (PREFERRED)
  if (isCodingSequence) {
    const silentResult = domesticateWithSilentMutations(seq, enzyme, {
      frame,
      isCodingSequence: true,
      organism,
      allowJunctionFallback: false,
      checkAllEnzymes,
    });

    result.strategies.silentMutation = silentResult;

    if (silentResult.success && silentResult.failedSites.length === 0) {
      // Perfect - all sites can be removed with silent mutations
      return {
        ...result,
        domesticatedSequence: silentResult.domesticatedSequence,
        strategy: DOMESTICATION_STRATEGY.SILENT_MUTATION,
        success: true,
        mutations: silentResult.mutations,
        verification: silentResult.verification,
        message: `Successfully domesticated ${internalSites.count} internal site(s) using silent mutations`,
        recommendations: [
          'Silent mutation domestication is compatible with one-pot Golden Gate assembly',
          'Verify protein expression if using rare codons',
        ],
      };
    }
  }

  // Step 3: Check alternative enzymes
  const alternativeEnzymes = recommendAlternativeEnzymes(seq, enzyme);
  const compatibleEnzymes = alternativeEnzymes.filter(e => e.isCompatible && !e.isCurrent);

  result.strategies.alternativeEnzyme = {
    available: compatibleEnzymes.length > 0,
    options: alternativeEnzymes,
    compatible: compatibleEnzymes,
  };

  if (compatibleEnzymes.length > 0 && preferAlternativeEnzyme) {
    const recommended = compatibleEnzymes[0];
    return {
      ...result,
      domesticatedSequence: seq, // No sequence change needed
      strategy: DOMESTICATION_STRATEGY.ALTERNATIVE_ENZYME,
      success: true,
      recommendedEnzyme: recommended,
      message: `Recommend switching to ${recommended.enzyme} which has no internal sites`,
      recommendations: [
        `Use ${recommended.enzyme} instead of ${enzyme}`,
        'No sequence modification required',
        'Compatible with one-pot Golden Gate assembly',
      ],
    };
  }

  // Step 4: Hybrid approach - silent mutations where possible, junction for the rest
  if (isCodingSequence && result.strategies.silentMutation) {
    const silentResult = result.strategies.silentMutation;

    if (silentResult.mutations.length > 0 && silentResult.failedSites.length > 0) {
      // Some sites handled by silent mutations, some need junction-based
      const partiallyDomesticated = silentResult.domesticatedSequence;

      if (allowJunctionFallback) {
        const junctionAnalysis = analyzeForJunctionDomestication(
          partiallyDomesticated,
          enzyme,
          { minFragmentSize: DOMESTICATION_DEFAULTS.minFragmentSize }
        );

        result.strategies.hybrid = {
          silentMutations: silentResult.mutations,
          junctionSites: silentResult.failedSites,
          junctionAnalysis,
        };

        result.warnings.push({
          ...UNIFIED_DOMESTICATION_CONFIG.warnings.junctionBased,
          affectedSites: silentResult.failedSites.length,
        });

        return {
          ...result,
          domesticatedSequence: partiallyDomesticated,
          strategy: DOMESTICATION_STRATEGY.HYBRID,
          success: true,
          mutations: silentResult.mutations,
          junctionRequired: silentResult.failedSites,
          message: `Hybrid domestication: ${silentResult.mutations.length} site(s) via silent mutations, ` +
                   `${silentResult.failedSites.length} site(s) require junction-based approach`,
          recommendations: [
            'Use sequential assembly protocol for junction-based sites',
            'Consider alternative enzyme if available',
          ],
        };
      }
    }
  }

  // Step 5: Junction-based fallback (LAST RESORT)
  if (allowJunctionFallback) {
    const junctionAnalysis = analyzeForJunctionDomestication(seq, enzyme, {
      minFragmentSize: DOMESTICATION_DEFAULTS.minFragmentSize,
    });

    result.strategies.junctionBased = junctionAnalysis;

    result.warnings.push(UNIFIED_DOMESTICATION_CONFIG.warnings.junctionBased);

    if (!junctionAnalysis.error) {
      return {
        ...result,
        domesticatedSequence: seq, // Sequence unchanged, fragments will be created
        strategy: DOMESTICATION_STRATEGY.JUNCTION_BASED,
        success: true,
        junctionAnalysis,
        additionalFragments: junctionAnalysis.additionalFragments,
        message: `Junction-based domestication: Will add ${junctionAnalysis.additionalFragments} fragment(s)`,
        recommendations: [
          'WARNING: NOT compatible with one-pot Golden Gate assembly',
          'Use sequential protocol: digest → heat-kill → ligate',
          'Consider using alternative enzyme: ' +
            (compatibleEnzymes.length > 0 ? compatibleEnzymes[0].enzyme : 'none available'),
        ],
      };
    }
  }

  // Step 6: No solution found
  return {
    ...result,
    domesticatedSequence: null,
    strategy: null,
    success: false,
    message: 'Unable to domesticate sequence with available methods',
    recommendations: [
      compatibleEnzymes.length > 0
        ? `Switch to ${compatibleEnzymes[0].enzyme} which has no internal sites`
        : 'Manual sequence modification may be required',
      'Consider synthesizing a domesticated version',
    ],
  };
}

// ============================================================================
// DETAILED ANALYSIS FUNCTIONS
// ============================================================================

/**
 * Get detailed analysis of all domestication options for a sequence
 *
 * @param {string} sequence - DNA sequence
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Analysis options
 * @returns {Object} Comprehensive analysis
 */
export function analyzeDomesticationOptions(sequence, enzyme = 'BsaI', options = {}) {
  const { frame = 0, organism = 'ecoli' } = options;

  const seq = sequence.toUpperCase();
  const internalSites = findInternalSites(seq, enzyme);

  if (!internalSites.hasSites) {
    return {
      needsDomestication: false,
      message: 'No internal sites found',
    };
  }

  // Analyze each site individually
  const siteAnalyses = internalSites.sites.map(site => {
    const silentCandidates = findAllSilentMutationCandidates(
      seq,
      site,
      frame,
      enzyme,
      organism === 'yeast' ? YEAST_CODON_USAGE : ECOLI_CODON_USAGE
    );

    const scoredCandidates = scoreMutationCandidates(silentCandidates, seq, enzyme, true);

    const junctionOptions = findDomesticationJunctions(seq, site, enzyme);

    return {
      site,
      silentMutation: {
        available: scoredCandidates.length > 0,
        candidateCount: scoredCandidates.length,
        bestCandidate: scoredCandidates[0] || null,
        topCandidates: scoredCandidates.slice(0, 5),
      },
      junctionBased: {
        available: junctionOptions.hasValidOption,
        candidateCount: junctionOptions.candidates.length,
        recommended: junctionOptions.recommended,
      },
    };
  });

  // Alternative enzymes
  const alternativeEnzymes = recommendAlternativeEnzymes(seq, enzyme);

  // Overall recommendation
  const allSilentAvailable = siteAnalyses.every(a => a.silentMutation.available);
  const somesilentAvailable = siteAnalyses.some(a => a.silentMutation.available);
  const hasCompatibleAlternative = alternativeEnzymes.some(e => e.isCompatible && !e.isCurrent);

  let recommendation;
  if (allSilentAvailable) {
    recommendation = {
      strategy: DOMESTICATION_STRATEGY.DIRECT_PRIMER_MUTATION,
      message: 'All sites can be removed with direct primer mutations (recommended)',
      onePotCompatible: true,
    };
  } else if (hasCompatibleAlternative) {
    const alt = alternativeEnzymes.find(e => e.isCompatible && !e.isCurrent);
    recommendation = {
      strategy: DOMESTICATION_STRATEGY.ALTERNATIVE_ENZYME,
      message: `Switch to ${alt.enzyme} which has no internal sites`,
      onePotCompatible: true,
      alternativeEnzyme: alt,
    };
  } else if (somesilentAvailable) {
    recommendation = {
      strategy: DOMESTICATION_STRATEGY.HYBRID,
      message: 'Some sites need mutagenic junction splitting - still one-pot compatible',
      onePotCompatible: true, // Now true with mutagenic junctions!
    };
  } else {
    recommendation = {
      strategy: DOMESTICATION_STRATEGY.MUTAGENIC_JUNCTION,
      message: 'Use mutagenic junction splitting (one-pot compatible)',
      onePotCompatible: true, // Mutagenic junctions ARE one-pot compatible!
    };
  }

  return {
    needsDomestication: true,
    siteCount: internalSites.count,
    siteAnalyses,
    alternativeEnzymes,
    recommendation,
    summary: {
      silentMutationSites: siteAnalyses.filter(a => a.silentMutation.available).length,
      junctionRequiredSites: siteAnalyses.filter(a => !a.silentMutation.available).length,
      compatibleAlternativeEnzymes: alternativeEnzymes.filter(e => e.isCompatible && !e.isCurrent).length,
    },
  };
}

/**
 * Generate a detailed report of domestication changes
 *
 * @param {Object} domesticationResult - Result from optimizeDomestication
 * @param {Object} options - Report options
 * @returns {Object} Detailed report
 */
export function generateDomesticationReport(domesticationResult, options = {}) {
  const { includeSequenceAlignment = true, frame = 0 } = options;

  const report = {
    summary: {
      success: domesticationResult.success,
      strategy: domesticationResult.strategy,
      enzyme: domesticationResult.enzyme,
    },
    changes: [],
    warnings: domesticationResult.warnings || [],
    recommendations: domesticationResult.recommendations || [],
  };

  // Document each mutation
  if (domesticationResult.mutations) {
    for (const mut of domesticationResult.mutations) {
      report.changes.push({
        type: 'silent_mutation',
        position: mut.mutation.sequencePosition,
        originalBase: mut.mutation.originalBase,
        newBase: mut.mutation.newBase,
        codon: {
          original: mut.mutation.originalCodon,
          new: mut.mutation.newCodon,
          aminoAcid: mut.mutation.aminoAcid,
        },
        siteRemoved: {
          position: mut.site.position,
          sequence: mut.site.sequence,
          orientation: mut.site.orientation,
        },
        score: mut.mutation.score,
      });
    }
  }

  // Sequence alignment
  if (includeSequenceAlignment && domesticationResult.originalSequence && domesticationResult.domesticatedSequence) {
    report.sequenceComparison = generateSequenceComparison(
      domesticationResult.originalSequence,
      domesticationResult.domesticatedSequence
    );
  }

  // Protein verification
  if (domesticationResult.domesticatedSequence) {
    report.proteinVerification = verifyProteinSequence(
      domesticationResult.originalSequence,
      domesticationResult.domesticatedSequence,
      frame
    );
  }

  return report;
}

/**
 * Generate a visual sequence comparison
 */
function generateSequenceComparison(original, domesticated) {
  const differences = [];
  const chunkSize = 60;

  for (let i = 0; i < original.length; i++) {
    if (original[i] !== domesticated[i]) {
      differences.push({
        position: i,
        original: original[i],
        domesticated: domesticated[i],
      });
    }
  }

  // Generate aligned chunks
  const chunks = [];
  for (let i = 0; i < original.length; i += chunkSize) {
    const origChunk = original.slice(i, i + chunkSize);
    const domChunk = domesticated.slice(i, i + chunkSize);

    let matchLine = '';
    for (let j = 0; j < origChunk.length; j++) {
      matchLine += origChunk[j] === domChunk[j] ? '|' : '*';
    }

    chunks.push({
      start: i,
      end: Math.min(i + chunkSize, original.length),
      original: origChunk,
      match: matchLine,
      domesticated: domChunk,
    });
  }

  return {
    totalDifferences: differences.length,
    differences,
    chunks,
  };
}

// ============================================================================
// EXPORTS
// ============================================================================

// Re-export from silent-mutation-domesticator
export {
  domesticateWithSilentMutations,
  findAllSilentMutationCandidates,
  scoreMutationCandidates,
  verifyProteinSequence,
  validateDomestication,
  compareDomesticationStrategies,
  applyMutation,
  applyMutations,
  CODON_TABLE,
  CODON_TO_AA,
  ECOLI_CODON_USAGE,
  YEAST_CODON_USAGE,
};

// Re-export mutagenic junction functions (PREFERRED for sites in middle of fragments)
export {
  designMutagenicJunction,
  designAllMutagenicJunctions,
  selectDomesticationStrategy,
  MUTAGENIC_JUNCTION_CONFIG,
};

// Re-export legacy junction-based functions (NOT recommended for one-pot assembly)
export {
  analyzeForJunctionDomestication,
  findDomesticationJunctions,
  recommendAlternativeEnzymes,
  optimizeGlobalOverhangSet,
  validatePostDomestication,
};
