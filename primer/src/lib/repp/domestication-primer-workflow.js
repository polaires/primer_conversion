/**
 * Domestication Primer Workflow
 *
 * Integrates Golden Gate domestication strategies with the Unified Primer Designer.
 * Provides a complete workflow from internal site detection to optimized primer design.
 *
 * Features:
 * 1. Automatic strategy selection (silent mutation vs mutagenic junction)
 * 2. Integration with NEB experimental fidelity data
 * 3. Comprehensive primer quality analysis
 * 4. Step-by-step workflow guidance
 * 5. Ready-to-order primer output with full QC
 *
 * This module bridges:
 * - Silent mutation domesticator
 * - Enhanced mutagenic junction optimizer
 * - Unified primer designer
 * - Thermodynamic analysis tools
 */

import {
  domesticateWithSilentMutations,
  findAllSilentMutationCandidates,
  scoreMutationCandidates,
  verifyProteinSequence,
  CODON_TO_AA,
  ECOLI_CODON_USAGE,
  YEAST_CODON_USAGE,
} from './silent-mutation-domesticator.js';

import {
  designEnhancedMutagenicJunction,
  designAllEnhancedJunctions,
  ENHANCED_JUNCTION_CONFIG,
  classifyQuality,
} from './enhanced-mutagenic-junction.js';

import {
  GOLDEN_GATE_ENZYMES,
  findInternalSites,
  calculateExperimentalFidelity,
  getEnzymeLigationData,
} from './goldengate.js';

import {
  OPTIMAL_FLANKING_SEQUENCES,
  OPTIMAL_SPACERS,
  GG_OPTIMIZER_DEFAULTS,
} from './goldengate-primer-optimizer.js';

import { recommendAlternativeEnzymes } from './auto-domestication-optimizer.js';
import { reverseComplement } from './enzymes.js';
import { calculateHairpinDG, calculateHomodimerDG, calculateHeterodimerDG } from '../equilibrium.js';
import { calculateTmQ5, calculate3primeTerminalDG } from '../tmQ5.js';
import {
  scoreTm,
  scoreGc,
  scoreHairpin,
  scoreHomodimer,
  scoreHeterodimer,
  scoreGcClamp,
  scoreHomopolymer,
  scoreTerminal3DG,
  scoreGQuadruplex,
  analyzeGQuadruplex,
} from '../scoring.js';

// ============================================================================
// CONFIGURATION
// ============================================================================

export const WORKFLOW_CONFIG = {
  // Strategy preferences
  preferSilentMutation: true,        // Prefer silent mutation over junction
  requireOnePotCompatible: true,     // Only one-pot compatible strategies

  // Quality thresholds
  minPrimerScore: 60,
  minAssemblyFidelity: 0.90,
  targetPrimerTm: 62,
  maxTmDifference: 5,

  // Primer design parameters
  homologyLength: { min: 15, optimal: 20, max: 30 },
  flankingLength: 6,

  // Output preferences
  includeOrderingInfo: true,
  includeThermoAnalysis: true,
  includeAlternatives: true,
};

// ============================================================================
// MAIN WORKFLOW FUNCTION
// ============================================================================

/**
 * Complete domestication workflow with integrated primer design
 *
 * This is the main entry point for state-of-the-art Golden Gate domestication.
 * It automatically:
 * 1. Detects internal restriction sites
 * 2. Selects optimal strategy (silent mutation or mutagenic junction)
 * 3. Designs high-quality primers with full thermodynamic analysis
 * 4. Validates the entire design
 * 5. Generates ready-to-order primer specifications
 *
 * @param {string} sequence - DNA sequence to domesticate
 * @param {string} enzyme - Enzyme name (default: 'BsaI')
 * @param {Object} options - Configuration options
 * @returns {Object} Complete workflow result with primers
 */
export function runDomesticationWorkflow(sequence, enzyme = 'BsaI', options = {}) {
  const {
    frame = 0,
    isCodingSequence = true,
    organism = 'ecoli',
    existingOverhangs = [],
    targetTm = WORKFLOW_CONFIG.targetPrimerTm,
    includeWorkflowGuide = true,
  } = options;

  const seq = sequence.toUpperCase();
  const startTime = Date.now();

  // Initialize result object
  const result = {
    success: false,
    sequence: seq,
    enzyme,
    organism,
    frame,

    // Analysis results
    analysis: null,
    strategy: null,
    domestication: null,

    // Primer designs
    primers: [],
    primerSummary: null,

    // Quality metrics
    quality: null,
    validation: null,

    // Workflow guide
    workflowGuide: null,

    // Timing
    processingTime: 0,
  };

  // ========== STEP 1: Analyze internal sites ==========
  const analysisResult = analyzeSequenceForDomestication(seq, enzyme, {
    frame,
    isCodingSequence,
    organism,
  });

  result.analysis = analysisResult;

  if (!analysisResult.needsDomestication) {
    result.success = true;
    result.strategy = 'none';
    result.message = 'No domestication needed - sequence is already compatible';

    if (includeWorkflowGuide) {
      result.workflowGuide = generateWorkflowGuide('none', result);
    }

    result.processingTime = Date.now() - startTime;
    return result;
  }

  // ========== STEP 2: Select and execute strategy ==========
  const strategyResult = selectAndExecuteStrategy(seq, enzyme, {
    frame,
    isCodingSequence,
    organism,
    existingOverhangs,
    analysisResult,
  });

  result.strategy = strategyResult.strategy;
  result.domestication = strategyResult;

  if (!strategyResult.success) {
    result.success = false;
    result.message = strategyResult.message;

    if (includeWorkflowGuide) {
      result.workflowGuide = generateWorkflowGuide('failed', result);
    }

    result.processingTime = Date.now() - startTime;
    return result;
  }

  // ========== STEP 3: Design integrated primers ==========
  const primerResult = designIntegratedPrimers(seq, strategyResult, enzyme, {
    frame,
    organism,
    targetTm,
    existingOverhangs,
  });

  result.primers = primerResult.primers;
  result.primerSummary = primerResult.summary;

  // ========== STEP 4: Validate complete design ==========
  const validationResult = validateCompleteDesign(result, enzyme, existingOverhangs);
  result.validation = validationResult;
  result.quality = validationResult.overallQuality;

  // ========== STEP 5: Generate workflow guide ==========
  if (includeWorkflowGuide) {
    result.workflowGuide = generateWorkflowGuide(result.strategy, result);
  }

  result.success = validationResult.isValid;
  result.message = generateResultMessage(result);
  result.processingTime = Date.now() - startTime;

  return result;
}

// ============================================================================
// ANALYSIS FUNCTIONS
// ============================================================================

/**
 * Comprehensive sequence analysis for domestication
 */
export function analyzeSequenceForDomestication(sequence, enzyme, options) {
  const { frame, isCodingSequence, organism } = options;

  const internalSites = findInternalSites(sequence, enzyme);

  if (!internalSites.hasSites) {
    return {
      needsDomestication: false,
      siteCount: 0,
      sites: [],
      message: `No internal ${enzyme} sites found`,
    };
  }

  const codonUsage = organism === 'yeast' ? YEAST_CODON_USAGE : ECOLI_CODON_USAGE;

  // Analyze each site
  const siteAnalyses = internalSites.sites.map(site => {
    // Check if silent mutation is possible
    const silentCandidates = isCodingSequence
      ? findAllSilentMutationCandidates(sequence, site, frame, enzyme, codonUsage)
      : [];

    const scoredCandidates = silentCandidates.length > 0
      ? scoreMutationCandidates(silentCandidates, sequence, enzyme, true)
      : [];

    const hasSilentOption = scoredCandidates.length > 0 && scoredCandidates[0].score >= 50;

    return {
      site,
      position: site.position,
      sequence: site.sequence,
      orientation: site.orientation,
      analysis: {
        silentMutationAvailable: hasSilentOption,
        silentMutationCount: scoredCandidates.length,
        bestSilentMutation: scoredCandidates[0] || null,
        requiresJunction: !hasSilentOption,
      },
    };
  });

  // Check alternative enzymes
  const alternativeEnzymes = recommendAlternativeEnzymes(sequence, enzyme);
  const compatibleAlternatives = alternativeEnzymes.filter(e => e.isCompatible && !e.isCurrent);

  // Determine recommended strategy
  const allSilentAvailable = siteAnalyses.every(s => s.analysis.silentMutationAvailable);
  const someSilentAvailable = siteAnalyses.some(s => s.analysis.silentMutationAvailable);

  let recommendedStrategy = 'mutagenic_junction';
  if (allSilentAvailable) {
    recommendedStrategy = 'silent_mutation';
  } else if (compatibleAlternatives.length > 0) {
    recommendedStrategy = 'alternative_enzyme';
  } else if (someSilentAvailable) {
    recommendedStrategy = 'hybrid';
  }

  return {
    needsDomestication: true,
    siteCount: internalSites.count,
    sites: internalSites.sites,
    siteAnalyses,
    alternativeEnzymes,
    compatibleAlternatives,
    recommendedStrategy,
    summary: {
      silentMutationSites: siteAnalyses.filter(s => s.analysis.silentMutationAvailable).length,
      junctionRequiredSites: siteAnalyses.filter(s => s.analysis.requiresJunction).length,
      totalSites: internalSites.count,
    },
  };
}

// ============================================================================
// STRATEGY SELECTION AND EXECUTION
// ============================================================================

/**
 * Select and execute the best domestication strategy
 */
export function selectAndExecuteStrategy(sequence, enzyme, options) {
  const {
    frame,
    isCodingSequence,
    organism,
    existingOverhangs,
    analysisResult,
  } = options;

  const { recommendedStrategy, siteAnalyses, sites } = analysisResult;

  // Strategy 1: All sites can use silent mutations
  if (recommendedStrategy === 'silent_mutation' && isCodingSequence) {
    const silentResult = domesticateWithSilentMutations(sequence, enzyme, {
      frame,
      isCodingSequence: true,
      organism,
      allowJunctionFallback: false,
    });

    if (silentResult.success && silentResult.failedSites.length === 0) {
      return {
        success: true,
        strategy: 'silent_mutation',
        domesticatedSequence: silentResult.domesticatedSequence,
        mutations: silentResult.mutations,
        additionalFragments: 0,
        onePotCompatible: true,
        message: `Applied ${silentResult.mutations.length} silent mutation(s)`,
      };
    }
  }

  // Strategy 2: Alternative enzyme available
  if (recommendedStrategy === 'alternative_enzyme') {
    const bestAlternative = analysisResult.compatibleAlternatives[0];
    return {
      success: true,
      strategy: 'alternative_enzyme',
      domesticatedSequence: sequence,
      recommendedEnzyme: bestAlternative.enzyme,
      currentEnzyme: enzyme,
      additionalFragments: 0,
      onePotCompatible: true,
      message: `Switch to ${bestAlternative.enzyme} - no internal sites`,
    };
  }

  // Strategy 3: Hybrid (some silent, some junction)
  if (recommendedStrategy === 'hybrid' && isCodingSequence) {
    // First apply silent mutations where possible
    const silentResult = domesticateWithSilentMutations(sequence, enzyme, {
      frame,
      isCodingSequence: true,
      organism,
      allowJunctionFallback: false,
    });

    const partiallyDomesticated = silentResult.domesticatedSequence;
    const remainingSites = findInternalSites(partiallyDomesticated, enzyme);

    if (remainingSites.hasSites) {
      // Design junctions for remaining sites
      const junctionResult = designAllEnhancedJunctions(partiallyDomesticated, enzyme, {
        frame,
        organism,
        existingOverhangs,
      });

      return {
        success: junctionResult.success || junctionResult.junctions.length > 0,
        strategy: 'hybrid',
        domesticatedSequence: partiallyDomesticated,
        mutations: silentResult.mutations,
        junctions: junctionResult.junctions,
        failedSites: junctionResult.failedSites,
        additionalFragments: junctionResult.junctions.length,
        fidelity: junctionResult.fidelity,
        onePotCompatible: true,
        message: `${silentResult.mutations.length} silent mutation(s), ${junctionResult.junctions.length} junction(s)`,
      };
    }

    return {
      success: true,
      strategy: 'silent_mutation',
      domesticatedSequence: partiallyDomesticated,
      mutations: silentResult.mutations,
      additionalFragments: 0,
      onePotCompatible: true,
      message: `Applied ${silentResult.mutations.length} silent mutation(s)`,
    };
  }

  // Strategy 4: Mutagenic junctions for all sites
  const junctionResult = designAllEnhancedJunctions(sequence, enzyme, {
    frame,
    organism,
    existingOverhangs,
    optimizeGlobalFidelity: true,
  });

  return {
    success: junctionResult.success || junctionResult.junctions.length > 0,
    strategy: 'mutagenic_junction',
    domesticatedSequence: sequence,
    junctions: junctionResult.junctions,
    failedSites: junctionResult.failedSites,
    additionalFragments: junctionResult.junctions.length,
    allOverhangs: junctionResult.allOverhangs,
    fidelity: junctionResult.fidelity,
    onePotCompatible: true,
    message: `Designed ${junctionResult.junctions.length} mutagenic junction(s)`,
  };
}

// ============================================================================
// INTEGRATED PRIMER DESIGN
// ============================================================================

/**
 * Design integrated primers with full thermodynamic analysis
 */
export function designIntegratedPrimers(sequence, strategyResult, enzyme, options) {
  const { frame, organism, targetTm, existingOverhangs } = options;

  const primers = [];
  const enz = GOLDEN_GATE_ENZYMES[enzyme];

  // Get optimal flanking and spacer
  const flankingConfig = OPTIMAL_FLANKING_SEQUENCES[enzyme] || OPTIMAL_FLANKING_SEQUENCES.BsaI;
  const spacerConfig = OPTIMAL_SPACERS[enzyme] || OPTIMAL_SPACERS.BsaI;
  const flanking = flankingConfig.default || 'GGTGCG';
  const spacer = spacerConfig.default || 'A';

  // Design primers based on strategy
  // Check both strategy and assemblyType since executeDomesticationPlan sets assemblyType
  const isSilentMutation = strategyResult.strategy === 'silent_mutation' ||
                           strategyResult.assemblyType === 'silent_mutation';

  if (isSilentMutation && strategyResult.mutations && strategyResult.mutations.length > 0) {
    // Silent Mutation Strategy - TWO STEP PROCESS:
    // Step 1: PCR mutagenesis primers to introduce mutations into template
    // Step 2: Golden Gate primers for final assembly

    const mutations = strategyResult.mutations || [];
    const originalSeq = strategyResult.originalSequence;
    const domesticatedSeq = strategyResult.domesticatedSequence;

    // STEP 1: Design mutagenesis PCR primers for each mutation site
    for (let i = 0; i < mutations.length; i++) {
      const mutation = mutations[i];
      const mutInfo = mutation.mutation || mutation;
      const sitePos = mutation.site?.position || mutInfo.position;

      // Design overlap extension PCR primers for site-directed mutagenesis
      const mutagenesisPrimers = designMutagenesisPrimers(
        originalSeq,
        mutInfo,
        sitePos,
        {
          targetTm: targetTm || 60,
          homologyLength: 20,
        }
      );

      primers.push({
        type: 'pcr_mutagenesis',
        step: 1,
        stepLabel: 'Step 1: Site-Directed Mutagenesis',
        name: `Mutagenesis Site ${i + 1}`,
        mutationIndex: i,
        site: mutation.site,
        mutation: mutInfo,
        codonChange: mutInfo.originalCodon
          ? `${mutInfo.originalCodon} → ${mutInfo.newCodon}`
          : `${mutInfo.originalBase} → ${mutInfo.newBase}`,
        ...mutagenesisPrimers,
        instructions: 'Use overlap extension PCR or QuikChange to introduce mutation',
      });
    }

    // STEP 2: Design Golden Gate primers for the domesticated template
    const ggPrimerPair = designFragmentPrimers(
      domesticatedSeq,
      0,
      domesticatedSeq.length,
      enzyme,
      {
        targetTm,
        flanking,
        spacer,
        overhangStart: existingOverhangs[0] || 'GGAG',
        overhangEnd: existingOverhangs[existingOverhangs.length - 1] || 'GCTT',
        mutations: [],
      }
    );

    primers.push({
      type: 'golden_gate',
      step: 2,
      stepLabel: 'Step 2: Golden Gate Assembly',
      name: 'Assembly Primers (for mutated template)',
      fragmentIndex: 0,
      ...ggPrimerPair,
      instructions: 'Use after obtaining mutated template from Step 1',
    });

  } else if (strategyResult.strategy === 'mutagenic_junction' ||
             strategyResult.strategy === 'hybrid' ||
             strategyResult.assemblyType === 'mutagenic_junction') {
    // For junction strategies, design primers for each fragment
    const junctions = strategyResult.junctions || [];

    // Sort junctions by position
    const sortedJunctions = [...junctions].sort((a, b) => a.junctionPosition - b.junctionPosition);

    // Design primers for each junction
    for (let i = 0; i < sortedJunctions.length; i++) {
      const junction = sortedJunctions[i];

      // Junction primers (from the enhanced design)
      const junctionPrimers = analyzeJunctionPrimers(junction, enzyme, targetTm);

      primers.push({
        type: 'junction',
        junctionIndex: i,
        junctionPosition: junction.junctionPosition,
        overhang: junction.overhang,
        mutation: junction.mutation,
        site: junction.site,
        ...junctionPrimers,
      });
    }

    // Add any silent mutation primers if hybrid
    if (strategyResult.mutations && strategyResult.mutations.length > 0) {
      for (const mutation of strategyResult.mutations) {
        primers.push({
          type: 'silent_mutation_info',
          mutation: mutation.mutation,
          site: mutation.site,
          note: 'Mutation applied directly to sequence - no additional primers needed',
        });
      }
    }
  }

  // Generate summary
  const summary = generatePrimerSummary(primers, enzyme);

  return {
    primers,
    summary,
  };
}

/**
 * Design primers for a fragment with specified overhangs
 */
function designFragmentPrimers(sequence, start, end, enzyme, options) {
  const { targetTm, flanking, spacer, overhangStart, overhangEnd, mutations } = options;
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  const recognition = enz.recognition;

  // Calculate homology length for target Tm
  const homologyLen = calculateHomologyLengthForTm(sequence, start, end, targetTm);

  // Forward primer homology (from start)
  const fwdHomology = sequence.slice(start, start + homologyLen.forward);

  // Reverse primer homology (from end)
  const revHomologyStart = Math.max(0, end - homologyLen.reverse);
  const revHomology = sequence.slice(revHomologyStart, end);

  // Build full primers
  const forwardPrimer = flanking + recognition + spacer + overhangStart + fwdHomology;
  const reversePrimer = flanking + recognition + spacer +
                        reverseComplement(overhangEnd) + reverseComplement(revHomology);

  // Analyze primers
  const fwdAnalysis = analyzePrimer(forwardPrimer, fwdHomology, true, targetTm);
  const revAnalysis = analyzePrimer(reversePrimer, revHomology, false, targetTm);

  return {
    forward: {
      sequence: forwardPrimer,
      ...fwdAnalysis,
      components: {
        flanking,
        recognition,
        spacer,
        overhang: overhangStart,
        homology: fwdHomology,
      },
    },
    reverse: {
      sequence: reversePrimer,
      ...revAnalysis,
      components: {
        flanking,
        recognition,
        spacer,
        overhang: reverseComplement(overhangEnd),
        homology: reverseComplement(revHomology),
      },
    },
    tmDifference: Math.abs(fwdAnalysis.homologyTm - revAnalysis.homologyTm),
    pairQuality: Math.round((fwdAnalysis.score + revAnalysis.score) / 2),
  };
}

/**
 * Design site-directed mutagenesis primers (for silent mutation strategy)
 * These primers are used in overlap extension PCR or QuikChange to introduce the mutation
 *
 * @param {string} sequence - Original DNA sequence
 * @param {Object} mutation - Mutation information (position, original/new base/codon)
 * @param {number} sitePosition - Position of the enzyme recognition site
 * @param {Object} options - Design options
 * @returns {Object} Mutagenesis primer pair
 */
function designMutagenesisPrimers(sequence, mutation, sitePosition, options = {}) {
  const {
    targetTm = 60,
    homologyLength = 20,  // Length of homology on each side of mutation
  } = options;

  // Get mutation position - could be in mutation object or calculated from site
  const mutPos = mutation.sequencePosition || mutation.position || sitePosition;

  // Design overlapping primers centered on the mutation
  // Forward mutagenic primer: ~15bp before mutation + mutation + ~15bp after
  const halfLength = Math.floor(homologyLength / 2);
  const fwdStart = Math.max(0, mutPos - halfLength);
  const fwdEnd = Math.min(sequence.length, mutPos + halfLength + 1);

  // Get the sequence and apply the mutation
  let fwdSequence = sequence.slice(fwdStart, fwdEnd);
  const mutPosInPrimer = mutPos - fwdStart;

  // Apply the mutation to the primer sequence
  if (mutation.newBase && mutPosInPrimer >= 0 && mutPosInPrimer < fwdSequence.length) {
    fwdSequence = fwdSequence.slice(0, mutPosInPrimer) +
                  mutation.newBase +
                  fwdSequence.slice(mutPosInPrimer + 1);
  }

  // Reverse primer is reverse complement of forward (for overlap extension)
  const revSequence = reverseComplement(fwdSequence);

  // Analyze primers
  const fwdAnalysis = analyzePrimer(fwdSequence, fwdSequence, true, targetTm);
  const revAnalysis = analyzePrimer(revSequence, revSequence, false, targetTm);

  // Calculate Tm for the overlap region
  const overlapTm = fwdAnalysis.homologyTm || calculateBasicTm(fwdSequence);

  return {
    forward: {
      sequence: fwdSequence,
      length: fwdSequence.length,
      ...fwdAnalysis,
      mutationPosition: mutPosInPrimer,
      hasMutation: true,
      components: {
        fivePrime: fwdSequence.slice(0, mutPosInPrimer),
        mutation: mutation.newBase || fwdSequence[mutPosInPrimer],
        threePrime: fwdSequence.slice(mutPosInPrimer + 1),
      },
    },
    reverse: {
      sequence: revSequence,
      length: revSequence.length,
      ...revAnalysis,
      hasMutation: true,
      components: {
        fivePrime: revSequence.slice(0, revSequence.length - mutPosInPrimer - 1),
        mutation: reverseComplement(mutation.newBase || 'N'),
        threePrime: revSequence.slice(revSequence.length - mutPosInPrimer),
      },
    },
    overlapTm,
    pcrProtocol: {
      annealingTemp: Math.round(overlapTm - 5),
      extensionTime: '30s per kb',
      method: 'Overlap Extension PCR or QuikChange',
    },
    pairQuality: Math.round((fwdAnalysis.score + revAnalysis.score) / 2),
  };
}

/**
 * Calculate basic Tm using simple 4+2 rule (fallback)
 */
function calculateBasicTm(sequence) {
  const at = (sequence.match(/[AT]/gi) || []).length;
  const gc = (sequence.match(/[GC]/gi) || []).length;
  if (sequence.length < 14) {
    return 2 * at + 4 * gc;
  }
  return 64.9 + 41 * (gc - 16.4) / sequence.length;
}

/**
 * Analyze junction primers with full thermodynamics
 */
function analyzeJunctionPrimers(junction, enzyme, targetTm) {
  const junctionPrimers = junction.primers;

  const fwdPrimer = junctionPrimers.fragment2.forwardPrimer;
  const revPrimer = junctionPrimers.fragment1.reversePrimer;

  // Analyze forward primer
  const fwdHomology = fwdPrimer.components.homology;
  const fwdAnalysis = analyzePrimer(fwdPrimer.sequence, fwdHomology, true, targetTm);

  // Analyze reverse primer
  const revHomology = revPrimer.components.homology;
  const revAnalysis = analyzePrimer(revPrimer.sequence, revHomology, false, targetTm);

  // Check heterodimer between the pair
  let heterodimerDG = 0;
  try {
    heterodimerDG = calculateHeterodimerDG(fwdPrimer.sequence, revPrimer.sequence);
  } catch (e) {}

  const heterodimerScore = scoreHeterodimer(heterodimerDG);

  return {
    forward: {
      sequence: fwdPrimer.sequence,
      ...fwdAnalysis,
      hasMutation: fwdPrimer.hasMutation,
      mutationPosition: fwdPrimer.mutationPosition,
      components: fwdPrimer.components,
    },
    reverse: {
      sequence: revPrimer.sequence,
      ...revAnalysis,
      hasMutation: revPrimer.hasMutation,
      mutationPosition: revPrimer.mutationPosition,
      components: revPrimer.components,
    },
    heterodimer: {
      dG: heterodimerDG,
      score: heterodimerScore,
      risk: heterodimerDG < -8 ? 'high' : heterodimerDG < -6 ? 'moderate' : 'low',
    },
    tmDifference: Math.abs(fwdAnalysis.homologyTm - revAnalysis.homologyTm),
    pairQuality: Math.round((fwdAnalysis.score + revAnalysis.score) / 2 * 0.8 + heterodimerScore * 20),
  };
}

/**
 * Analyze a single primer with full thermodynamics
 */
function analyzePrimer(fullSequence, homologyRegion, isForward, targetTm) {
  const seq = fullSequence.toUpperCase();
  const homology = homologyRegion.toUpperCase();

  // Calculate Tm for homology region
  let homologyTm = 62;
  try {
    homologyTm = calculateTmQ5(homology);
  } catch (e) {
    // Fallback
    const gc = (homology.match(/[GC]/g) || []).length;
    homologyTm = 64.9 + 41 * (gc - 16.4) / homology.length;
  }

  // GC content
  const gcCount = (homology.match(/[GC]/g) || []).length;
  const gc = gcCount / homology.length;

  // Thermodynamic properties
  let hairpinDG = 0, homodimerDG = 0, terminal3DG = -8;
  try {
    hairpinDG = calculateHairpinDG(seq);
    homodimerDG = calculateHomodimerDG(seq);
    const termResult = calculate3primeTerminalDG(homology);
    terminal3DG = termResult?.dG ?? -8;
  } catch (e) {}

  // G-quadruplex analysis
  const gQuadruplex = analyzeGQuadruplex(seq);

  // Score components
  const tmScore = scoreTm(homologyTm, { optimalTm: targetTm, tolerance: 3 });
  const gcScore = scoreGc(gc);
  const hairpinScore = scoreHairpin(hairpinDG);
  const homodimerScore = scoreHomodimer(homodimerDG);
  const terminal3DGScore = scoreTerminal3DG(terminal3DG);
  const gcClampScore = scoreGcClamp(homology);
  const homopolymerScore = scoreHomopolymer(seq);
  const gQuadScore = scoreGQuadruplex(seq);

  // Weighted composite score
  const score = Math.round(
    tmScore * 20 +
    gcScore * 10 +
    hairpinScore * 15 +
    homodimerScore * 10 +
    terminal3DGScore * 15 +
    gcClampScore * 10 +
    homopolymerScore * 5 +
    gQuadScore * 15
  );

  // Collect issues
  const issues = [];
  const warnings = [];

  if (tmScore < 0.7) issues.push(`Tm ${homologyTm.toFixed(1)}°C outside optimal range`);
  if (gcScore < 0.5) issues.push(`GC ${(gc * 100).toFixed(0)}% outside optimal range`);
  if (hairpinScore < 0.5) warnings.push(`Hairpin potential (ΔG=${hairpinDG.toFixed(1)})`);
  if (homodimerScore < 0.5) warnings.push(`Self-dimer potential (ΔG=${homodimerDG.toFixed(1)})`);
  if (terminal3DGScore < 0.5) warnings.push(`Weak 3' binding (ΔG=${terminal3DG.toFixed(1)})`);
  if (gcClampScore < 0.5) issues.push('No GC clamp at 3\' end');
  if (gQuadScore < 0.5) warnings.push('G-quadruplex risk');

  // Determine quality tier
  const qualityTier = score >= 85 ? 'excellent' : score >= 70 ? 'good' : score >= 55 ? 'acceptable' : 'poor';

  return {
    length: seq.length,
    homologyLength: homology.length,
    homologyTm,
    gc: gc * 100,
    gcPercent: `${(gc * 100).toFixed(1)}%`,

    thermodynamics: {
      hairpinDG,
      homodimerDG,
      terminal3DG,
    },

    scores: {
      tm: tmScore,
      gc: gcScore,
      hairpin: hairpinScore,
      homodimer: homodimerScore,
      terminal3DG: terminal3DGScore,
      gcClamp: gcClampScore,
      homopolymer: homopolymerScore,
      gQuadruplex: gQuadScore,
    },

    score,
    qualityTier,
    issues,
    warnings,

    hasGCClamp: gcClampScore >= 0.9,
    gQuadruplex,
  };
}

/**
 * Calculate optimal homology length for target Tm
 */
function calculateHomologyLengthForTm(sequence, start, end, targetTm) {
  const config = WORKFLOW_CONFIG.homologyLength;

  // Start with optimal length and adjust
  let fwdLen = config.optimal;
  let revLen = config.optimal;

  // Try to hit target Tm for forward primer
  for (let len = config.min; len <= config.max; len++) {
    const homology = sequence.slice(start, start + len);
    try {
      const tm = calculateTmQ5(homology);
      if (Math.abs(tm - targetTm) < Math.abs(calculateTmQ5(sequence.slice(start, start + fwdLen)) - targetTm)) {
        fwdLen = len;
      }
      if (Math.abs(tm - targetTm) < 2) break;
    } catch (e) {}
  }

  // Try to hit target Tm for reverse primer
  for (let len = config.min; len <= config.max; len++) {
    const homology = sequence.slice(end - len, end);
    try {
      const tm = calculateTmQ5(homology);
      if (Math.abs(tm - targetTm) < Math.abs(calculateTmQ5(sequence.slice(end - revLen, end)) - targetTm)) {
        revLen = len;
      }
      if (Math.abs(tm - targetTm) < 2) break;
    } catch (e) {}
  }

  return { forward: fwdLen, reverse: revLen };
}

/**
 * Generate primer summary for ordering
 * Creates cognitive names based on primer type and purpose
 */
export function generatePrimerSummary(primers, enzyme) {
  const summary = {
    totalPrimers: 0,
    junctionPrimers: 0,
    fragmentPrimers: 0,
    pcrPrimers: 0,
    goldenGatePrimers: 0,
    orderList: [],
    totalCost: 0,  // Estimated
  };

  // Track indices for each primer type for meaningful naming
  let sdmSiteIndex = 1;  // Site-Directed Mutagenesis
  let junctionIndex = 1;
  let ggIndex = 1;

  for (const primer of primers) {
    if (primer.type === 'silent_mutation_info') continue;

    // Determine cognitive name based on primer type
    let baseName;
    let notes;

    switch (primer.type) {
      case 'pcr_mutagenesis':
        baseName = `SDM_Site${primer.mutationIndex !== undefined ? primer.mutationIndex + 1 : sdmSiteIndex}`;
        notes = primer.codonChange ? `Mutagenesis: ${primer.codonChange}` : 'Site-directed mutagenesis';
        sdmSiteIndex++;
        summary.pcrPrimers += 2;
        break;
      case 'golden_gate':
        baseName = ggIndex === 1 ? 'GG_Assembly' : `GG_Frag${ggIndex}`;
        notes = 'Golden Gate assembly primer';
        ggIndex++;
        summary.goldenGatePrimers += 2;
        break;
      case 'junction':
        baseName = `Junction${primer.junctionIndex !== undefined ? primer.junctionIndex + 1 : junctionIndex}`;
        notes = primer.junctionPosition !== undefined
          ? `Junction at position ${primer.junctionPosition}`
          : 'Mutagenic junction primer';
        junctionIndex++;
        summary.junctionPrimers += 2;
        break;
      default:
        baseName = `Frag${summary.fragmentPrimers / 2 + 1}`;
        notes = 'Fragment primer';
        summary.fragmentPrimers += 2;
    }

    if (primer.forward) {
      summary.orderList.push({
        name: `${baseName}_Fwd`,
        sequence: primer.forward.sequence,
        length: primer.forward.length,
        tm: primer.forward.homologyTm,
        scale: '25nmol',
        purification: 'Standard',
        notes,
        type: primer.type,
        step: primer.step,
      });
      summary.totalPrimers++;
    }

    if (primer.reverse) {
      summary.orderList.push({
        name: `${baseName}_Rev`,
        sequence: primer.reverse.sequence,
        length: primer.reverse.length,
        tm: primer.reverse.homologyTm,
        scale: '25nmol',
        purification: 'Standard',
        notes,
        type: primer.type,
        step: primer.step,
      });
      summary.totalPrimers++;
    }
  }

  // Calculate total length
  summary.totalLength = summary.orderList.reduce((sum, p) => sum + (p.length || 0), 0);

  // Estimate cost (rough estimate: $0.30/base for standard primers)
  summary.totalCost = summary.orderList.reduce((sum, p) => sum + p.length * 0.30, 0);
  summary.estimatedCost = summary.totalCost;

  return summary;
}

// ============================================================================
// VALIDATION
// ============================================================================

/**
 * Validate the complete design
 */
export function validateCompleteDesign(result, enzyme, existingOverhangs) {
  const validations = [];
  let isValid = true;
  const warnings = [];
  const errors = [];

  // 1. Check strategy success
  if (!result.domestication || !result.domestication.success) {
    errors.push({ code: 'STRATEGY_FAILED', message: 'Domestication strategy failed' });
    isValid = false;
  }

  // 2. Check assembly fidelity
  if (result.domestication?.fidelity) {
    const fidelity = result.domestication.fidelity.assembly || result.domestication.fidelity;
    validations.push({
      check: 'ASSEMBLY_FIDELITY',
      passed: fidelity >= WORKFLOW_CONFIG.minAssemblyFidelity,
      value: fidelity,
      message: `Assembly fidelity: ${(fidelity * 100).toFixed(1)}%`,
    });

    if (fidelity < WORKFLOW_CONFIG.minAssemblyFidelity) {
      warnings.push({
        code: 'LOW_FIDELITY',
        message: `Assembly fidelity ${(fidelity * 100).toFixed(1)}% below target ${WORKFLOW_CONFIG.minAssemblyFidelity * 100}%`,
      });
    }
  }

  // 3. Check primer quality
  for (const primer of result.primers) {
    if (primer.forward && primer.forward.score < WORKFLOW_CONFIG.minPrimerScore) {
      warnings.push({
        code: 'LOW_PRIMER_QUALITY',
        message: `Forward primer score ${primer.forward.score} below threshold`,
        primer: primer.forward.sequence,
      });
    }
    if (primer.reverse && primer.reverse.score < WORKFLOW_CONFIG.minPrimerScore) {
      warnings.push({
        code: 'LOW_PRIMER_QUALITY',
        message: `Reverse primer score ${primer.reverse.score} below threshold`,
        primer: primer.reverse.sequence,
      });
    }

    // Check Tm difference
    if (primer.tmDifference > WORKFLOW_CONFIG.maxTmDifference) {
      warnings.push({
        code: 'HIGH_TM_DIFF',
        message: `Primer pair Tm difference ${primer.tmDifference.toFixed(1)}°C exceeds ${WORKFLOW_CONFIG.maxTmDifference}°C`,
      });
    }
  }

  // 4. Check for protein preservation (if applicable)
  if (result.domestication?.domesticatedSequence && result.frame !== undefined) {
    const proteinCheck = verifyProteinSequence(
      result.sequence,
      result.domestication.domesticatedSequence,
      result.frame
    );

    validations.push({
      check: 'PROTEIN_PRESERVED',
      passed: proteinCheck.identical,
      message: proteinCheck.message,
    });

    if (!proteinCheck.identical) {
      errors.push({
        code: 'PROTEIN_CHANGED',
        message: 'Domestication changed protein sequence!',
        details: proteinCheck.differences,
      });
      isValid = false;
    }
  }

  // Overall quality
  const avgPrimerScore = result.primers
    .filter(p => p.pairQuality)
    .reduce((sum, p) => sum + p.pairQuality, 0) / Math.max(1, result.primers.length);

  const overallQuality = {
    score: avgPrimerScore,
    tier: avgPrimerScore >= 85 ? 'excellent' : avgPrimerScore >= 70 ? 'good' : avgPrimerScore >= 55 ? 'acceptable' : 'poor',
  };

  return {
    isValid,
    validations,
    warnings,
    errors,
    overallQuality,
  };
}

// ============================================================================
// WORKFLOW GUIDE GENERATION
// ============================================================================

/**
 * Generate step-by-step workflow guide
 */
export function generateWorkflowGuide(strategy, result) {
  const steps = [];

  if (strategy === 'none') {
    steps.push({
      step: 1,
      title: 'No Domestication Required',
      description: 'Your sequence has no internal restriction sites for the selected enzyme.',
      action: 'Proceed directly to Golden Gate assembly',
      status: 'complete',
    });

    return { steps, strategy: 'none' };
  }

  if (strategy === 'failed') {
    steps.push({
      step: 1,
      title: 'Domestication Failed',
      description: 'Unable to find a suitable domestication strategy.',
      action: 'Consider using an alternative enzyme or manual sequence modification',
      status: 'error',
      recommendations: result.analysis?.compatibleAlternatives?.map(e => e.enzyme) || [],
    });

    return { steps, strategy: 'failed' };
  }

  // Successful domestication workflow
  steps.push({
    step: 1,
    title: 'Analysis Complete',
    description: `Found ${result.analysis?.siteCount || 0} internal site(s) requiring domestication.`,
    details: {
      sites: result.analysis?.sites?.map(s => ({ position: s.position, sequence: s.sequence })),
      strategy: result.strategy,
    },
    status: 'complete',
  });

  if (strategy === 'silent_mutation') {
    steps.push({
      step: 2,
      title: 'Silent Mutations Applied',
      description: 'Synonymous mutations introduced to remove internal sites while preserving protein sequence.',
      details: {
        mutations: result.domestication?.mutations?.map(m => ({
          position: m.mutation.sequencePosition,
          change: `${m.mutation.originalCodon}→${m.mutation.newCodon} (${m.mutation.aminoAcid})`,
        })),
      },
      status: 'complete',
    });
  } else if (strategy === 'mutagenic_junction' || strategy === 'hybrid') {
    steps.push({
      step: 2,
      title: 'Mutagenic Junctions Designed',
      description: 'Fragment boundaries with embedded mutations designed using NEB fidelity data.',
      details: {
        junctions: result.domestication?.junctions?.map(j => ({
          position: j.junctionPosition,
          overhang: j.overhang,
          fidelity: j.fidelity?.singleOverhang,
          quality: j.quality?.tier,
        })),
        fidelity: result.domestication?.fidelity,
      },
      status: 'complete',
    });

    if (strategy === 'hybrid' && result.domestication?.mutations?.length > 0) {
      steps.push({
        step: 3,
        title: 'Silent Mutations Applied',
        description: 'Additional silent mutations applied for sites near fragment boundaries.',
        details: {
          count: result.domestication.mutations.length,
        },
        status: 'complete',
      });
    }
  } else if (strategy === 'alternative_enzyme') {
    steps.push({
      step: 2,
      title: 'Alternative Enzyme Recommended',
      description: `Switch to ${result.domestication?.recommendedEnzyme} which has no internal sites.`,
      action: `Re-design assembly using ${result.domestication?.recommendedEnzyme}`,
      status: 'action_required',
    });

    return { steps, strategy };
  }

  // Primer ordering step
  steps.push({
    step: steps.length + 1,
    title: 'Order Primers',
    description: 'Optimized primers designed with full thermodynamic analysis.',
    details: {
      totalPrimers: result.primerSummary?.totalPrimers,
      estimatedCost: result.primerSummary?.totalCost?.toFixed(2),
      orderList: result.primerSummary?.orderList,
    },
    action: 'Copy primer sequences to order from your preferred supplier',
    status: 'ready',
  });

  // PCR step
  steps.push({
    step: steps.length + 1,
    title: 'PCR Amplification',
    description: 'Amplify fragments using the designed primers.',
    protocol: {
      enzyme: 'Q5 High-Fidelity DNA Polymerase (NEB M0491)',
      annealingTemp: `${result.primers[0]?.forward?.homologyTm?.toFixed(0) || 62}°C`,
      extensionTime: 'Calculate based on fragment size (30 sec/kb)',
      cycles: 30,
    },
    status: 'pending',
  });

  // Assembly step
  steps.push({
    step: steps.length + 1,
    title: 'Golden Gate Assembly',
    description: 'Assemble fragments using Golden Gate reaction.',
    protocol: {
      enzyme: result.enzyme,
      ligase: 'T4 DNA Ligase (NEB M0202)',
      buffer: 'T4 DNA Ligase Buffer',
      cycles: '30 cycles of 37°C (5 min) / 16°C (5 min), then 60°C (5 min)',
    },
    status: 'pending',
  });

  return {
    steps,
    strategy,
    summary: {
      totalSteps: steps.length,
      currentStep: steps.findIndex(s => s.status === 'ready' || s.status === 'pending') + 1,
      quality: result.quality,
    },
  };
}

// ============================================================================
// UTILITIES
// ============================================================================

/**
 * Generate result message
 */
function generateResultMessage(result) {
  const parts = [];

  if (result.success) {
    parts.push(`Successfully domesticated sequence using ${result.strategy.replace('_', ' ')} strategy.`);

    if (result.domestication?.additionalFragments > 0) {
      parts.push(`Created ${result.domestication.additionalFragments} additional fragment(s).`);
    }

    if (result.domestication?.fidelity) {
      const fidelity = result.domestication.fidelity.assembly || result.domestication.fidelity;
      if (typeof fidelity === 'number') {
        parts.push(`Assembly fidelity: ${(fidelity * 100).toFixed(1)}%.`);
      }
    }

    if (result.primerSummary) {
      parts.push(`${result.primerSummary.totalPrimers} primers designed.`);
    }
  } else {
    parts.push('Domestication failed.');
    if (result.validation?.errors?.length > 0) {
      parts.push(result.validation.errors[0].message);
    }
  }

  return parts.join(' ');
}

