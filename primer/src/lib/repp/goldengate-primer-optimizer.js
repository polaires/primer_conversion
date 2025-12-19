/**
 * Golden Gate Primer Optimizer
 *
 * State-of-the-art optimization for Golden Gate Assembly primer design.
 * Implements NEB best practices and Potapov et al. (2018) research.
 *
 * Key optimizations:
 * 1. Optimal 6bp flanking sequences (NEB recommendation)
 * 2. G:T mismatch detection and mitigation
 * 3. Automatic overhang optimization for maximum fidelity
 * 4. Comprehensive primer quality scoring
 *
 * References:
 * - NEB Golden Gate Assembly Guidelines (6bp flanking recommendation)
 * - Potapov et al. (2018) ACS Synth Biol - Overhang fidelity profiling
 * - Pryor et al. (2020) PLOS ONE - Data-optimized assembly design
 */

import { reverseComplement } from './enzymes.js';
import {
  GOLDEN_GATE_ENZYMES,
  OVERHANG_FIDELITY,
  HIGH_FIDELITY_SETS,
  designGoldenGatePrimers,
  getEnzymeLigationData,
  getOverhangFidelityExperimental,
  getRecommendedOverhangs,
  getOptimalOverhangSetExperimental,
  findOptimalOverhangSet,
  findInternalSites,
  suggestDomestication,
  findAlternativeEnzymes,
  calculateExperimentalFidelity,
} from './goldengate.js';
import { calculateHairpinDG, calculateHomodimerDG, calculateHeterodimerDG } from '../equilibrium.js';
import { calculate3primeTerminalDG } from '../tmQ5.js';
import {
  scoreTm,
  scoreGc,
  scoreHairpin,
  scoreHomodimer,
  scoreHeterodimer,
  scoreLength,
  scoreGcClamp,
  scoreHomopolymer,
  scoreTerminal3DG,
  score3PrimeComposition,
  scoreGQuadruplex,
  classifyQuality,
} from '../scoring.js';
import { DEFAULT_WEIGHTS } from '../weightCalibration.js';

// ============================================================================
// CONFIGURABLE DEFAULTS
// ============================================================================

/**
 * Default configuration for Golden Gate optimization
 * All values can be overridden via options parameters
 */
export const GG_OPTIMIZER_DEFAULTS = {
  // G:T mismatch handling
  gtMismatchFactor: 0.20,           // G:T wobbles ligate at ~15-25% of correct rate
  gtRiskMatchThreshold: 3,          // Minimum matches to flag as G:T risk

  // Fidelity thresholds
  defaultFidelityFallback: 0.85,    // Fallback when fidelity data unavailable
  minAcceptableFidelity: 0.90,      // Minimum fidelity before warning
  targetFidelity: 0.95,             // Target fidelity for optimization
  excellentFidelityThreshold: 0.95, // Threshold for "excellent" rating
  goodFidelityThreshold: 0.90,      // Threshold for "good" rating
  acceptableFidelityThreshold: 0.80, // Threshold for "acceptable" rating

  // Flanking sequence requirements
  optimalFlankingLength: 6,         // NEB recommended flanking length
  minFlankingLength: 4,             // Minimum acceptable flanking

  // GC content thresholds (as fractions)
  gcOptimalMin: 0.40,               // Minimum optimal GC
  gcOptimalMax: 0.60,               // Maximum optimal GC
  gcAcceptableMin: 0.30,            // Minimum acceptable GC
  gcAcceptableMax: 0.70,            // Maximum acceptable GC

  // Homology region constraints
  homologyMinLength: 15,            // Minimum homology region
  homologyMaxLength: 30,            // Maximum homology region

  // Primer pair constraints
  tmDiffWarningThreshold: 5,        // Tm difference warning threshold (°C)
  tmDiffPenaltyMultiplier: 4,       // Penalty multiplier for Tm diff > threshold

  // Scoring weights - ALIGNED with empirically calibrated DEFAULT_WEIGHTS
  // Base weights from 829-primer Döring dataset (F1=88.7%, AUC=0.930)
  // GG-specific scoring is applied as an additive layer (not replacing base)
  //
  // Key changes from previous config:
  // - terminal3DG: 0.08 → 0.18 (empirically proven critical for 3' stability)
  // - gQuadruplex: added (critical failure predictor)
  // - ggSpecific: kept as additive bonus/penalty layer
  scoringWeights: {
    // Tier 1: Critical (from calibrated weights)
    terminal3DG: 0.18,        // Was 0.08, now aligned with DEFAULT_WEIGHTS (0.20)
    gQuadruplex: 0.12,        // NEW: Critical failure predictor

    // Tier 2: Important
    hairpin: 0.10,
    homodimer: 0.10,
    composition3Prime: 0.08,  // 3' end composition

    // Tier 3: Standard
    tm: 0.10,                 // Was 0.15
    gc: 0.06,                 // Was 0.08
    gcClamp: 0.06,            // Was 0.08
    length: 0.05,             // Was 0.08
    homopolymer: 0.05,        // Was 0.06

    // GG-Specific: Applied as additive layer (not counted in base normalization)
    // This is a bonus/penalty that adjusts final score based on assembly-specific factors
    ggSpecific: 0.10,         // Was 0.15 - now additive, not normalized
  },

  // Primer quality thresholds - ALIGNED with Primer Designer (classifyQuality)
  // Previous: 85/70/55 - Now aligned with scoring.js thresholds
  excellentPrimerScore: 90,   // Was 85 - aligned with classifyQuality
  goodPrimerScore: 75,        // Was 70 - aligned with classifyQuality
  acceptablePrimerScore: 60,  // Was 55 - aligned with classifyQuality

  // Optimization limits
  maxOverhangReplacements: 5,       // Max replacements in auto-optimization
};

// ============================================================================
// OPTIMAL FLANKING SEQUENCES
// ============================================================================

/**
 * Optimal 6bp flanking sequences for Type IIS enzymes
 * Based on NEB research and enzyme cleavage efficiency studies
 *
 * Criteria for optimal flanking:
 * 1. 40-60% GC content (ideally 50%)
 * 2. No homopolymers (3+ consecutive identical bases)
 * 3. No palindromes
 * 4. Avoid sequences that form primer dimers
 * 5. Avoid sequences creating hairpins with recognition site
 */
export const OPTIMAL_FLANKING_SEQUENCES = {
  BsaI: {
    default: 'GGTGCG',     // 67% GC, balanced, no issues
    alternatives: [
      'GCGTCG',  // 67% GC
      'CGTGCG',  // 67% GC
      'TGCGAG',  // 50% GC
      'CGCTGC',  // 67% GC
      'ATGCGC',  // 50% GC
      'GCATGC',  // 50% GC
    ],
  },
  BsmBI: {
    default: 'GCTGCG',     // 67% GC
    alternatives: [
      'GCGTCG',
      'TGCGCG',
      'CGATGC',
    ],
  },
  BbsI: {
    default: 'GCGTGC',     // BbsI has 2bp spacer
    alternatives: [
      'TGCGCG',
      'CGCTGC',
      'ATGCGC',
    ],
  },
  Esp3I: {
    default: 'GCTGCG',     // Same as BsmBI (same recognition)
    alternatives: [
      'GCGTCG',
      'TGCGCG',
    ],
  },
  SapI: {
    default: 'GCTGCG',     // SapI creates 3bp overhangs
    alternatives: [
      'GCGTCG',
      'ATGCGC',
    ],
  },
};

/**
 * Optimal spacer nucleotides by enzyme
 * Based on NEB cleavage efficiency data
 */
export const OPTIMAL_SPACERS = {
  BsaI: {
    forward: 'A',   // A gives best cleavage for BsaI
    reverse: 'T',   // T (complement of A)
  },
  BsmBI: {
    forward: 'A',
    reverse: 'T',
  },
  BbsI: {
    forward: 'AA',  // BbsI needs 2bp spacer
    reverse: 'TT',
  },
  Esp3I: {
    forward: 'A',
    reverse: 'T',
  },
  SapI: {
    forward: 'A',
    reverse: 'T',
  },
};

// ============================================================================
// INPUT VALIDATION HELPERS
// ============================================================================

/**
 * Validate that input is a non-empty array
 * @param {*} arr - Input to validate
 * @param {string} name - Parameter name for error messages
 * @throws {Error} If validation fails
 */
function validateNonEmptyArray(arr, name) {
  if (!Array.isArray(arr)) {
    throw new Error(`${name} must be an array`);
  }
  if (arr.length === 0) {
    throw new Error(`${name} cannot be empty`);
  }
}

/**
 * Validate that input is a valid DNA sequence
 * @param {string} seq - Sequence to validate
 * @param {string} name - Parameter name for error messages
 * @throws {Error} If validation fails
 */
function validateDNASequence(seq, name) {
  if (typeof seq !== 'string') {
    throw new Error(`${name} must be a string`);
  }
  if (seq.length === 0) {
    throw new Error(`${name} cannot be empty`);
  }
  if (!/^[ATGCatgc]+$/.test(seq)) {
    throw new Error(`${name} contains invalid characters (only ATGC allowed)`);
  }
}

/**
 * Validate overhang sequence
 * @param {string} overhang - Overhang to validate
 * @param {number} expectedLength - Expected length (default 4)
 * @throws {Error} If validation fails
 */
function validateOverhang(overhang, expectedLength = 4) {
  validateDNASequence(overhang, 'Overhang');
  if (overhang.length !== expectedLength) {
    throw new Error(`Overhang must be exactly ${expectedLength} bases, got ${overhang.length}`);
  }
}

// ============================================================================
// G:T MISMATCH DETECTION
// ============================================================================

/**
 * Check if two bases can pair (including G:T wobbles)
 * @param {string} base1 - First base
 * @param {string} base2 - Second base
 * @param {Object} config - Configuration options
 * @returns {Object} Pairing information
 */
export function canPair(base1, base2, config = {}) {
  const { gtMismatchFactor = GG_OPTIMIZER_DEFAULTS.gtMismatchFactor } = config;

  const wc = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' };

  // Normalize to uppercase
  const b1 = (base1 || '').toUpperCase();
  const b2 = (base2 || '').toUpperCase();

  if (wc[b1] === b2) {
    return { pairs: true, isWobble: false, factor: 1.0 };
  }

  // G:T wobble pairs (accepted by T4 ligase at reduced efficiency)
  if ((b1 === 'G' && b2 === 'T') || (b1 === 'T' && b2 === 'G')) {
    return { pairs: true, isWobble: true, factor: gtMismatchFactor };
  }

  return { pairs: false, isWobble: false, factor: 0 };
}

/**
 * Calculate G:T mismatch risk for an overhang set
 * Returns pairs that could mis-ligate via G:T wobbles
 *
 * @param {string[]} overhangs - Array of overhangs to analyze
 * @param {Object} config - Configuration options
 * @returns {Array} Array of risk objects
 */
export function findGTMismatchRisks(overhangs, config = {}) {
  const {
    gtMismatchFactor = GG_OPTIMIZER_DEFAULTS.gtMismatchFactor,
    gtRiskMatchThreshold = GG_OPTIMIZER_DEFAULTS.gtRiskMatchThreshold,
  } = config;

  // Handle empty or invalid input gracefully
  if (!Array.isArray(overhangs) || overhangs.length < 2) {
    return [];
  }

  const risks = [];

  for (let i = 0; i < overhangs.length; i++) {
    for (let j = i + 1; j < overhangs.length; j++) {
      const oh1 = (overhangs[i] || '').toUpperCase();
      const oh2 = (overhangs[j] || '').toUpperCase();

      // Skip invalid overhangs
      if (!oh1 || !oh2) continue;

      const oh2rc = reverseComplement(oh2);

      // Count G:T wobbles if oh1 tried to ligate with RC of oh2
      let wobbleCount = 0;
      let matchCount = 0;
      const wobblePositions = [];

      for (let k = 0; k < oh1.length && k < oh2rc.length; k++) {
        const pair = canPair(oh1[k], oh2rc[k], { gtMismatchFactor });
        if (pair.isWobble) {
          wobbleCount++;
          wobblePositions.push(k);
        }
        if (pair.pairs) matchCount++;
      }

      // Risk: 3+ matches with at least 1 G:T wobble is dangerous
      // T4 ligase can accept these mismatches
      if (matchCount >= gtRiskMatchThreshold && wobbleCount >= 1) {
        const expectedMisLigation = Math.pow(gtMismatchFactor, wobbleCount);

        // Position-weighted risk: wobbles at 3' end (position 3) are more dangerous
        const positionWeight = wobblePositions.reduce((sum, pos) => sum + (pos + 1), 0) / wobblePositions.length;
        const positionRisk = positionWeight / oh1.length; // 0.25-1.0

        risks.push({
          overhang1: oh1,
          overhang2: oh2,
          reverseComplement2: oh2rc,
          wobblePositions,
          wobbleCount,
          matchCount,
          positionWeight,
          risk: matchCount === oh1.length ? 'critical' :
                matchCount >= gtRiskMatchThreshold && positionRisk > 0.5 ? 'high' : 'medium',
          expectedMisLigation,
          expectedMisLigationPercent: `${(expectedMisLigation * 100).toFixed(1)}%`,
          description: `${oh1} may mis-ligate with RC(${oh2})=${oh2rc} via ${wobbleCount} G:T wobble(s)`,
        });
      }
    }
  }

  return risks;
}

/**
 * Calculate enhanced fidelity including G:T mismatch penalties
 *
 * PRIMARY METHOD: Uses experimental ligation frequency matrix to calculate
 * cross-reactivity between overhangs in this specific assembly set.
 * This is more accurate than product-of-individual-fidelities.
 *
 * FALLBACK: Uses static OVERHANG_FIDELITY data when matrix unavailable.
 *
 * @param {string[]} overhangs - Array of overhangs
 * @param {string} enzyme - Enzyme name
 * @param {Object} config - Configuration options
 * @returns {Object} Enhanced fidelity data
 */
export function calculateEnhancedFidelity(overhangs, enzyme = 'BsaI', config = {}) {
  const {
    defaultFidelityFallback = GG_OPTIMIZER_DEFAULTS.defaultFidelityFallback,
    gtMismatchFactor = GG_OPTIMIZER_DEFAULTS.gtMismatchFactor,
    useMatrixData = true, // Default to matrix-based calculation (most accurate)
  } = config;

  // Handle empty or invalid input gracefully
  if (!Array.isArray(overhangs) || overhangs.length === 0) {
    return {
      baseFidelity: 0,
      baseFidelityPercent: '0.0%',
      gtAdjustedFidelity: 0,
      gtAdjustedFidelityPercent: '0.0%',
      gtPenalty: 1.0,
      gtPenaltyPercent: '0.0%',
      gtRisks: [],
      junctionFidelities: [],
      hasGTRisks: false,
      error: 'No overhangs provided',
    };
  }

  let baseFidelity = 1.0;
  let junctionFidelities = [];
  let calculationMethod = 'static';

  // PRIMARY: Use matrix-based calculation (accounts for cross-reactivity within this set)
  if (useMatrixData) {
    const matrixResult = calculateExperimentalFidelity(overhangs, enzyme);

    // Check if we got valid matrix-based results
    if (matrixResult.source === 'experimental' && matrixResult.assemblyFidelity > 0) {
      baseFidelity = matrixResult.assemblyFidelity;
      calculationMethod = 'matrix';

      // Convert junction format for consistency
      junctionFidelities = matrixResult.junctions.map(j => ({
        overhang: j.overhang,
        fidelity: j.fidelity,
        category: j.fidelity >= 0.95 ? 'excellent' :
                  j.fidelity >= 0.90 ? 'good' :
                  j.fidelity >= 0.80 ? 'medium' : 'low',
        correctFreq: j.correctFreq,
        totalFreq: j.totalFreq,
      }));
    }
  }

  // FALLBACK: Use static data if matrix calculation failed or not requested
  if (calculationMethod === 'static') {
    baseFidelity = 1.0;
    junctionFidelities = [];

    for (const oh of overhangs) {
      if (!oh) continue;

      const ohUpper = oh.toUpperCase();
      const staticData = OVERHANG_FIDELITY[ohUpper];
      const junctionFidelity = staticData?.fidelity || defaultFidelityFallback;
      const category = staticData?.category || 'unknown';

      baseFidelity *= junctionFidelity;
      junctionFidelities.push({
        overhang: ohUpper,
        fidelity: junctionFidelity,
        category,
      });
    }
  }

  // Find G:T mismatch risks (applies to both methods)
  const gtRisks = findGTMismatchRisks(overhangs, { gtMismatchFactor });

  // Adjust fidelity based on G:T risks
  // Note: For matrix calculation, cross-ligation is already accounted for,
  // but G:T wobbles may not be fully captured in the experimental data
  let gtPenalty = 1.0;
  if (calculationMethod === 'static') {
    // Only apply G:T penalty for static method (matrix already includes cross-ligation)
    gtRisks.forEach(risk => {
      gtPenalty *= (1 - risk.expectedMisLigation);
    });
  }

  const gtAdjustedFidelity = baseFidelity * gtPenalty;

  return {
    baseFidelity,
    baseFidelityPercent: `${(baseFidelity * 100).toFixed(1)}%`,
    gtAdjustedFidelity,
    gtAdjustedFidelityPercent: `${(gtAdjustedFidelity * 100).toFixed(1)}%`,
    gtPenalty,
    gtPenaltyPercent: `${((1 - gtPenalty) * 100).toFixed(1)}%`,
    gtRisks,
    junctionFidelities,
    hasGTRisks: gtRisks.length > 0,
    calculationMethod, // 'matrix' or 'static'
  };
}

// ============================================================================
// FLANKING SEQUENCE OPTIMIZATION
// ============================================================================

/**
 * Check if a sequence contains homopolymer runs
 * @param {string} seq - Sequence to check
 * @param {number} minRun - Minimum run length to flag (default 3)
 * @returns {boolean} True if homopolymer found
 */
function hasHomopolymer(seq, minRun = 3) {
  if (!seq) return false;
  const pattern = new RegExp(`(.)\\1{${minRun - 1},}`);
  return pattern.test(seq);
}

/**
 * Check if a flanking sequence could cause mispriming on the template
 *
 * Mispriming occurs when the flanking region (which shouldn't bind to template)
 * has significant complementarity to somewhere on the template, potentially
 * causing the primer to anneal at wrong locations during PCR.
 *
 * @param {string} flanking - Flanking sequence to check
 * @param {string} template - Template sequence to check against
 * @param {Object} options - Configuration options
 * @returns {Object} Mispriming risk assessment
 */
function checkFlankingMispriming(flanking, template, options = {}) {
  const {
    minMatchLength = 5,      // Minimum consecutive match to flag
    maxMismatches = 1,       // Allow 1 mismatch in longer matches
  } = options;

  if (!flanking || !template || flanking.length < minMatchLength) {
    return { risk: 'none', matches: [], penalty: 0 };
  }

  const flankingUpper = flanking.toUpperCase();
  const templateUpper = template.toUpperCase();
  const flankingRC = reverseComplement(flankingUpper);
  const matches = [];

  // Check for direct matches (flanking binds to template sense strand)
  // and reverse complement matches (flanking binds to template antisense strand)
  for (const [searchSeq, matchType] of [[flankingUpper, 'direct'], [flankingRC, 'reverse_complement']]) {
    // Sliding window search for matches
    for (let i = 0; i <= templateUpper.length - minMatchLength; i++) {
      let matchLen = 0;
      let mismatches = 0;

      for (let j = 0; j < searchSeq.length && i + j < templateUpper.length; j++) {
        if (searchSeq[j] === templateUpper[i + j]) {
          matchLen++;
        } else {
          mismatches++;
          if (mismatches > maxMismatches) break;
        }
      }

      // Flag if we found a significant match
      if (matchLen >= minMatchLength) {
        matches.push({
          type: matchType,
          position: i,
          matchLength: matchLen,
          mismatches,
          templateRegion: templateUpper.slice(i, i + searchSeq.length),
        });
      }
    }
  }

  // Determine risk level and penalty
  let risk = 'none';
  let penalty = 0;

  if (matches.length > 0) {
    const bestMatch = matches.reduce((a, b) => a.matchLength > b.matchLength ? a : b);

    if (bestMatch.matchLength >= flanking.length - 1) {
      // Near-perfect match - high risk
      risk = 'high';
      penalty = 30;
    } else if (bestMatch.matchLength >= minMatchLength + 1) {
      // Significant partial match - medium risk
      risk = 'medium';
      penalty = 15;
    } else {
      // Minimal match - low risk
      risk = 'low';
      penalty = 5;
    }
  }

  return { risk, matches, penalty, bestMatch: matches[0] };
}

/**
 * Generate candidate flanking sequences dynamically
 *
 * Creates 6bp sequences with good properties:
 * - 40-60% GC content (3-4 G/C bases)
 * - No homopolymer runs (3+ consecutive identical bases)
 * - Not palindromic
 *
 * @param {number} length - Flanking length (default 6)
 * @param {number} maxCandidates - Maximum candidates to generate
 * @returns {string[]} Array of candidate sequences
 */
function generateFlankingCandidates(length = 6, maxCandidates = 20) {
  const bases = ['A', 'T', 'G', 'C'];
  const candidates = [];

  // Generate sequences with balanced GC (3-4 out of 6)
  // Use a deterministic approach for reproducibility
  const goodPatterns = [
    'ATGCGC', 'GCATGC', 'TGCGAT', 'CGTAGC', 'GATCGC', 'CGATCG',
    'TAGCGC', 'GCTACG', 'ACGTGC', 'TGACGC', 'GCTAGC', 'CGTACG',
    'ATCGCG', 'TACGCG', 'GTACGC', 'CATGCG', 'AGCTGC', 'TCAGCG',
    'GACGTC', 'CTGACG',
  ];

  for (const pattern of goodPatterns) {
    if (candidates.length >= maxCandidates) break;

    // Validate pattern
    if (pattern.length === length &&
        !hasHomopolymer(pattern) &&
        !isPalindrome(pattern)) {
      const gc = gcContent(pattern);
      if (gc >= 0.33 && gc <= 0.67) {
        candidates.push(pattern);
      }
    }
  }

  return candidates;
}

/**
 * Check if a sequence is palindromic
 * @param {string} seq - Sequence to check
 * @returns {boolean} True if palindromic
 */
function isPalindrome(seq) {
  if (!seq) return false;
  return seq.toUpperCase() === reverseComplement(seq.toUpperCase());
}

/**
 * Calculate GC content as fraction
 * @param {string} seq - Sequence
 * @returns {number} GC content (0-1)
 */
function gcContent(seq) {
  if (!seq || seq.length === 0) return 0;
  const gc = (seq.match(/[GC]/gi) || []).length;
  return gc / seq.length;
}

/**
 * Score a flanking sequence for quality with full context
 *
 * Evaluates flanking sequence considering:
 * 1. Basic properties (length, GC, homopolymers, palindrome)
 * 2. Interaction with recognition site (combined secondary structure)
 * 3. Full primer tail secondary structure (flanking + recognition + spacer + overhang)
 * 4. Mispriming risk against template (if template provided)
 *
 * @param {string} flanking - Flanking sequence
 * @param {string} recognitionSite - Enzyme recognition site
 * @param {string} overhang - Overhang sequence
 * @param {string} homologyStart - First few bases of homology region
 * @param {Object} config - Configuration options
 * @returns {Object} Score and analysis
 */
function scoreFlanking(flanking, recognitionSite, overhang, homologyStart, config = {}) {
  const {
    optimalFlankingLength = GG_OPTIMIZER_DEFAULTS.optimalFlankingLength,
    gcAcceptableMin = GG_OPTIMIZER_DEFAULTS.gcAcceptableMin,
    gcAcceptableMax = GG_OPTIMIZER_DEFAULTS.gcAcceptableMax,
    template = null,  // NEW: template sequence for mispriming check
  } = config;

  let score = 100;
  const issues = [];

  if (!flanking) {
    return { flanking: '', score: 0, gc: 0, issues: [{ message: 'No flanking sequence' }], quality: 'poor', mispriming: null };
  }

  // Check length (should be 6bp)
  if (flanking.length < optimalFlankingLength) {
    const penalty = (optimalFlankingLength - flanking.length) * 5;
    score -= penalty;
    issues.push(`Suboptimal length: ${flanking.length}bp (recommended: ${optimalFlankingLength}bp)`);
  }

  // Check for homopolymers
  if (hasHomopolymer(flanking)) {
    score -= 15;
    issues.push('Contains homopolymer run');
  }

  // Check GC content (ideal: 40-60%)
  const gc = gcContent(flanking);
  if (gc < gcAcceptableMin || gc > gcAcceptableMax) {
    score -= 10;
    issues.push(`GC content ${(gc * 100).toFixed(0)}% outside optimal range (${gcAcceptableMin * 100}-${gcAcceptableMax * 100}%)`);
  }

  // Check for palindrome
  if (isPalindrome(flanking)) {
    score -= 10;
    issues.push('Palindromic sequence');
  }

  // Check if flanking + recognition creates hairpin potential
  const combined = flanking + (recognitionSite || '');
  if (hasHomopolymer(combined)) {
    score -= 5;
    issues.push('Creates homopolymer with recognition site');
  }

  // Check full primer prefix for secondary structure
  const primerPrefix = flanking + (recognitionSite || '') + 'A' + (overhang || '');
  let hairpinDG = 0;
  try {
    hairpinDG = calculateHairpinDG(primerPrefix);
    if (hairpinDG < -3) {
      score -= 15;
      issues.push(`Primer prefix forms stable hairpin (ΔG=${hairpinDG.toFixed(1)} kcal/mol)`);
    }
  } catch (e) {
    // Skip if calculation fails - not critical
  }

  // NEW: Check for mispriming risk against template
  let mispriming = null;
  if (template) {
    mispriming = checkFlankingMispriming(flanking, template);
    if (mispriming.penalty > 0) {
      score -= mispriming.penalty;
      if (mispriming.risk === 'high') {
        issues.push(`HIGH mispriming risk: flanking matches template at position ${mispriming.bestMatch?.position}`);
      } else if (mispriming.risk === 'medium') {
        issues.push(`Medium mispriming risk: partial match to template (${mispriming.bestMatch?.matchLength}bp)`);
      } else if (mispriming.risk === 'low') {
        issues.push(`Low mispriming risk: minor match to template`);
      }
    }
  }

  return {
    flanking,
    score: Math.max(0, score),
    gc: gc * 100,
    issues,
    quality: score >= 90 ? 'excellent' : score >= 70 ? 'good' : score >= 50 ? 'acceptable' : 'poor',
    hairpinDG,
    mispriming,
  };
}

/**
 * Select optimal flanking sequence based on primer context
 *
 * PRIORITY LOGIC:
 * 1. First try the default flanking sequence for this enzyme
 * 2. If default has no significant issues (score >= 70), use it
 * 3. If default has problems, try predefined alternatives
 * 4. If all predefined alternatives have problems, try dynamically generated candidates
 * 5. Return the best option with explanation of why it was chosen
 *
 * @param {string} enzyme - Enzyme name
 * @param {string} overhang - Overhang sequence
 * @param {string} homologyStart - First 10bp of homology region
 * @param {Object} options - Additional options
 * @param {string} options.template - Template sequence for mispriming check
 * @param {string} options.customFlanking - User-specified flanking (highest priority)
 * @param {number} options.acceptableThreshold - Minimum score to accept default (default: 70)
 * @returns {Object} Best flanking sequence with analysis
 */
export function selectOptimalFlankingSequence(enzyme, overhang, homologyStart, options = {}) {
  const {
    template = null,
    customFlanking = null,
    acceptableThreshold = 70,
  } = options;

  const enzData = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enzData) {
    console.warn(`Unknown enzyme: ${enzyme}, falling back to BsaI flanking sequences`);
  }

  const recognitionSite = enzData?.recognition || GOLDEN_GATE_ENZYMES.BsaI.recognition;
  const flankingOptions = OPTIMAL_FLANKING_SEQUENCES[enzyme] || OPTIMAL_FLANKING_SEQUENCES.BsaI;

  // PRIORITY 1: If user specified custom flanking, validate and use it
  if (customFlanking) {
    const customScore = scoreFlanking(customFlanking, recognitionSite, overhang, homologyStart, { ...options, template });
    return {
      best: customScore,
      alternatives: [],
      enzyme,
      recognitionSite,
      selectionReason: 'user_specified',
      message: 'Using user-specified flanking sequence',
    };
  }

  // PRIORITY 2: Try the default flanking first
  const defaultFlanking = flankingOptions.default || 'GGTGCG';
  const defaultScore = scoreFlanking(defaultFlanking, recognitionSite, overhang, homologyStart, { ...options, template });

  // If default is acceptable (no significant issues), use it
  if (defaultScore.score >= acceptableThreshold) {
    return {
      best: defaultScore,
      alternatives: [],
      enzyme,
      recognitionSite,
      selectionReason: 'default_acceptable',
      message: 'Using default flanking sequence (no issues detected)',
    };
  }

  // PRIORITY 3: Default has issues - try predefined alternatives
  const alternatives = (flankingOptions.alternatives || []).filter(Boolean);
  const alternativeScores = alternatives.map(flanking =>
    scoreFlanking(flanking, recognitionSite, overhang, homologyStart, { ...options, template })
  );

  // Find best alternative
  const bestAlternative = alternativeScores.reduce(
    (best, current) => (current.score > best.score ? current : best),
    { score: -1 }
  );

  // If best alternative is better than default and acceptable, use it
  if (bestAlternative.score > defaultScore.score && bestAlternative.score >= acceptableThreshold) {
    return {
      best: bestAlternative,
      alternatives: [defaultScore, ...alternativeScores.filter(s => s !== bestAlternative)].slice(0, 3),
      enzyme,
      recognitionSite,
      selectionReason: 'alternative_better',
      message: `Switched from default (score: ${defaultScore.score}) to alternative (score: ${bestAlternative.score}) due to: ${defaultScore.issues.join(', ')}`,
      defaultIssues: defaultScore.issues,
    };
  }

  // PRIORITY 4: All predefined options have issues - try dynamic generation
  const dynamicCandidates = generateFlankingCandidates(6, 20);
  const dynamicScores = dynamicCandidates.map(flanking =>
    scoreFlanking(flanking, recognitionSite, overhang, homologyStart, { ...options, template })
  );

  const bestDynamic = dynamicScores.reduce(
    (best, current) => (current.score > best.score ? current : best),
    { score: -1 }
  );

  // If dynamic candidate is better than all predefined options
  if (bestDynamic.score > defaultScore.score && bestDynamic.score > bestAlternative.score) {
    return {
      best: bestDynamic,
      alternatives: [defaultScore, bestAlternative].filter(s => s.score > 0).slice(0, 3),
      enzyme,
      recognitionSite,
      selectionReason: 'dynamic_generated',
      message: `Generated new flanking sequence (score: ${bestDynamic.score}) because all predefined options had issues`,
      defaultIssues: defaultScore.issues,
    };
  }

  // FALLBACK: Return the best of all options (even if below threshold)
  const allOptions = [defaultScore, bestAlternative, bestDynamic].filter(s => s.score >= 0);
  allOptions.sort((a, b) => b.score - a.score);

  const best = allOptions[0] || defaultScore;

  return {
    best,
    alternatives: allOptions.slice(1, 4),
    enzyme,
    recognitionSite,
    selectionReason: best === defaultScore ? 'default_fallback' : 'best_available',
    message: best.score < acceptableThreshold
      ? `Warning: Best available flanking has issues (score: ${best.score}). Consider manual review.`
      : `Selected best available option (score: ${best.score})`,
    defaultIssues: defaultScore.issues,
  };
}

// ============================================================================
// AUTOMATIC OVERHANG OPTIMIZATION
// ============================================================================

/**
 * Check if overhang matches a required pattern
 * Pattern uses '.' as wildcard, e.g., '..TG' matches any overhang ending in TG
 *
 * @param {string} overhang - Overhang to check
 * @param {string} pattern - Pattern to match (use '.' for wildcard)
 * @returns {boolean} True if matches
 */
function matchesPattern(overhang, pattern) {
  if (!pattern) return true;
  if (overhang.length !== pattern.length) return false;

  for (let i = 0; i < pattern.length; i++) {
    if (pattern[i] !== '.' && pattern[i].toUpperCase() !== overhang[i].toUpperCase()) {
      return false;
    }
  }
  return true;
}

/**
 * Find best alternative overhang considering full context
 * @param {string} currentOH - Current overhang to replace
 * @param {string[]} allOverhangs - All overhangs in the set
 * @param {number} position - Position of overhang to replace
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Options
 * @returns {Object|null} Best alternative or null
 */
function findBestAlternativeOverhang(currentOH, allOverhangs, position, enzyme, options = {}) {
  const {
    requiredPattern = null,
    defaultFidelityFallback = GG_OPTIMIZER_DEFAULTS.defaultFidelityFallback,
  } = options;

  // Get all high-fidelity overhangs
  let candidates = Object.entries(OVERHANG_FIDELITY)
    .filter(([oh, data]) => data.category === 'excellent' || data.category === 'good')
    .map(([oh, data]) => ({ overhang: oh, ...data }));

  // Filter out already used overhangs and their reverse complements
  const usedSet = new Set(
    allOverhangs.flatMap(oh => [oh.toUpperCase(), reverseComplement(oh.toUpperCase())])
  );

  candidates = candidates.filter(c => {
    const ohUpper = c.overhang.toUpperCase();
    return !usedSet.has(ohUpper) && !usedSet.has(reverseComplement(ohUpper));
  });

  // Apply pattern filter if specified
  if (requiredPattern) {
    candidates = candidates.filter(c => matchesPattern(c.overhang, requiredPattern));
  }

  // Handle edge case: no candidates available
  if (candidates.length === 0) {
    return null;
  }

  // Score each candidate in context
  const scored = candidates.map(candidate => {
    const testSet = [...allOverhangs];
    testSet[position] = candidate.overhang;

    // Calculate fidelity with this replacement
    const setFidelity = calculateEnhancedFidelity(testSet, enzyme, { defaultFidelityFallback });

    return {
      ...candidate,
      contextFidelity: setFidelity.gtAdjustedFidelity,
      gtRiskCount: setFidelity.gtRisks.length,
      // Combined score: fidelity - G:T risk penalty
      score: setFidelity.gtAdjustedFidelity - (setFidelity.gtRisks.length * 0.02),
    };
  });

  // Sort by score
  scored.sort((a, b) => b.score - a.score);

  // Return best option if it improves over current
  const best = scored[0];
  if (best) {
    const currentFidelity = OVERHANG_FIDELITY[currentOH.toUpperCase()]?.fidelity || 0.5;

    return {
      overhang: best.overhang,
      fidelity: best.fidelity,
      contextFidelity: best.contextFidelity,
      improvement: best.fidelity - currentFidelity,
      reason: `Replaced ${currentOH} (${(currentFidelity * 100).toFixed(0)}%) with ${best.overhang} (${(best.fidelity * 100).toFixed(0)}%)`,
    };
  }

  return null;
}

/**
 * Auto-optimize overhangs for maximum assembly fidelity
 *
 * Replaces problematic overhangs while respecting constraints:
 * - Required overhangs (e.g., AATG for start codon)
 * - Pattern requirements (e.g., '..TG' for start codon context)
 * - Standard fusion sites
 *
 * @param {string[]} requestedOverhangs - Initial overhang set
 * @param {Object} options - Optimization options
 * @returns {Object} Optimization result
 */
export function autoOptimizeOverhangs(requestedOverhangs, options = {}) {
  const {
    enzyme = 'BsaI',
    requiredIndices = [],     // Indices of overhangs that cannot be changed
    requiredPatterns = {},    // Map of index -> pattern (e.g., { 2: '..TG' })
    minFidelity = GG_OPTIMIZER_DEFAULTS.minAcceptableFidelity,
    maxReplacements = GG_OPTIMIZER_DEFAULTS.maxOverhangReplacements,
    defaultFidelityFallback = GG_OPTIMIZER_DEFAULTS.defaultFidelityFallback,
  } = options;

  // Handle empty or invalid input
  if (!Array.isArray(requestedOverhangs) || requestedOverhangs.length === 0) {
    return {
      original: [],
      optimized: [],
      replacements: [],
      unchangedCount: 0,
      error: 'No overhangs provided',
    };
  }

  // Normalize overhangs
  const overhangs = requestedOverhangs.map(oh => (oh || '').toUpperCase());

  // Analyze each overhang using validated static data
  const analysis = overhangs.map((oh, idx) => {
    const staticData = OVERHANG_FIDELITY[oh];
    const fidelity = staticData?.fidelity || defaultFidelityFallback;
    const category = staticData?.category || 'unknown';

    return {
      index: idx,
      overhang: oh,
      fidelity,
      category,
      isRequired: requiredIndices.includes(idx),
      requiredPattern: requiredPatterns[idx] || null,
      needsReplacement: fidelity < minFidelity || category === 'avoid' || category === 'low',
    };
  });

  // Sort by fidelity (worst first) for replacement priority
  const toReplace = analysis
    .filter(item => item.needsReplacement && !item.isRequired)
    .sort((a, b) => a.fidelity - b.fidelity);

  // Find replacements for problematic overhangs
  const optimized = [...overhangs];
  const replacements = [];

  for (const item of toReplace.slice(0, maxReplacements)) {
    const alternative = findBestAlternativeOverhang(
      item.overhang,
      optimized,
      item.index,
      enzyme,
      {
        requiredPattern: item.requiredPattern,
        defaultFidelityFallback,
      }
    );

    if (alternative && alternative.fidelity > item.fidelity) {
      optimized[item.index] = alternative.overhang;
      replacements.push({
        index: item.index,
        original: item.overhang,
        originalFidelity: item.fidelity,
        replacement: alternative.overhang,
        replacementFidelity: alternative.fidelity,
        improvement: alternative.improvement,
        reason: alternative.reason,
      });
    }
  }

  // Calculate overall fidelity before and after
  const originalFidelity = calculateEnhancedFidelity(overhangs, enzyme, { defaultFidelityFallback });
  const optimizedFidelity = calculateEnhancedFidelity(optimized, enzyme, { defaultFidelityFallback });

  return {
    original: overhangs,
    optimized,
    replacements,
    unchangedCount: overhangs.length - replacements.length,

    analysis: {
      original: analysis,
    },

    fidelityImprovement: {
      before: originalFidelity.gtAdjustedFidelity,
      beforePercent: originalFidelity.gtAdjustedFidelityPercent,
      after: optimizedFidelity.gtAdjustedFidelity,
      afterPercent: optimizedFidelity.gtAdjustedFidelityPercent,
      improvement: optimizedFidelity.gtAdjustedFidelity - originalFidelity.gtAdjustedFidelity,
      improvementPercent: `+${((optimizedFidelity.gtAdjustedFidelity - originalFidelity.gtAdjustedFidelity) * 100).toFixed(1)}%`,
    },

    gtRisks: {
      before: originalFidelity.gtRisks,
      after: optimizedFidelity.gtRisks,
      resolved: originalFidelity.gtRisks.length - optimizedFidelity.gtRisks.length,
    },
  };
}

// ============================================================================
// GOLDEN GATE PRIMER QUALITY SCORING
// ============================================================================

/**
 * Score Golden Gate specific criteria
 * @param {Object} primer - Primer object from designGoldenGatePrimers
 * @param {string} enzyme - Enzyme name
 * @param {Object} config - Configuration options
 * @returns {Object} GG-specific score
 */
function scoreGoldenGateSpecific(primer, enzyme = 'BsaI', config = {}) {
  const {
    optimalFlankingLength = GG_OPTIMIZER_DEFAULTS.optimalFlankingLength,
    homologyMinLength = GG_OPTIMIZER_DEFAULTS.homologyMinLength,
    homologyMaxLength = GG_OPTIMIZER_DEFAULTS.homologyMaxLength,
  } = config;

  let score = 1.0;
  const issues = [];

  // Safely extract structure properties
  const structure = primer?.structure || {};
  const extra = structure.extra || '';
  const homology = structure.homology || primer?.homologyRegion || '';

  // Check flanking length (should be 6bp)
  if (extra && extra.length < optimalFlankingLength) {
    const penalty = (optimalFlankingLength - extra.length) * 0.05;
    score -= penalty;
    issues.push({
      type: 'flanking_length',
      severity: 'warning',
      message: `Flanking sequence ${extra.length}bp (optimal: ${optimalFlankingLength}bp)`,
      penalty,
    });
  }

  // Check for homopolymers in flanking
  if (extra && hasHomopolymer(extra)) {
    score -= 0.05;
    issues.push({
      type: 'flanking_homopolymer',
      severity: 'minor',
      message: 'Homopolymer in flanking sequence',
      penalty: 0.05,
    });
  }

  // Check for internal recognition sites in homology
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  if (enz && homology && homology.includes(enz.recognition)) {
    score -= 0.30;
    issues.push({
      type: 'internal_site',
      severity: 'critical',
      message: `Internal ${enzyme} recognition site in homology region`,
      penalty: 0.30,
    });
  }

  // Check homology GC clamp (last 2 bases should have 1-2 G/C)
  if (homology && homology.length >= 2) {
    const last2 = homology.slice(-2);
    const gcInLast2 = (last2.match(/[GC]/gi) || []).length;
    if (gcInLast2 === 0) {
      score -= 0.10;
      issues.push({
        type: 'weak_gc_clamp',
        severity: 'warning',
        message: 'Weak 3\' GC clamp in homology region',
        penalty: 0.10,
      });
    }
  }

  // Check homology length
  const homologyLen = homology ? homology.length : 0;
  if (homologyLen > 0 && homologyLen < homologyMinLength) {
    score -= 0.15;
    issues.push({
      type: 'short_homology',
      severity: 'warning',
      message: `Short homology region: ${homologyLen}bp (minimum: ${homologyMinLength}bp)`,
      penalty: 0.15,
    });
  } else if (homologyLen > homologyMaxLength) {
    score -= 0.05;
    issues.push({
      type: 'long_homology',
      severity: 'minor',
      message: `Long homology region: ${homologyLen}bp (may increase cost)`,
      penalty: 0.05,
    });
  }

  return {
    score: Math.max(0, score),
    issues,
  };
}

/**
 * Score Golden Gate primer quality
 * Applies thermodynamic analysis to the homology region
 *
 * @param {Object} primer - Primer object from designGoldenGatePrimers
 * @param {Object} options - Scoring options
 * @returns {Object} Comprehensive quality score
 */
export function scoreGoldenGatePrimer(primer, options = {}) {
  const {
    enzyme = 'BsaI',
    weights = GG_OPTIMIZER_DEFAULTS.scoringWeights,
  } = options;

  // Handle null/undefined primer
  if (!primer) {
    return {
      composite: 0,
      quality: 'poor',
      breakdown: {},
      error: 'No primer provided',
    };
  }

  const homology = primer.structure?.homology || primer.homologyRegion || '';
  const fullSequence = primer.sequence || '';
  const homologyLength = homology.length;

  // Score the homology region (binding portion)
  const tm = primer.tm || 60;
  const gc = primer.gc || 50;

  // Basic scores
  const tmScore = scoreTm(tm);
  const gcScore = scoreGc(gc / 100);

  // Enhanced scores - length (optimal 18-25bp for homology region)
  const lengthScore = scoreLength(homologyLength, {
    optimalLow: 18,
    optimalHigh: 25,
    acceptableLow: 15,
    acceptableHigh: 30,
  });

  // GC clamp at 3' end
  const gcClampScore = homology ? scoreGcClamp(homology) : 1.0;

  // Homopolymer runs
  const homopolymerScore = homology ? scoreHomopolymer(homology) : 1.0;

  // Thermodynamic analysis
  let hairpinDG = 0;
  let homodimerDG = 0;
  let fullHairpinDG = 0;
  let terminal3DG = -8; // Default to middle of optimal range (-11 to -6)

  try {
    if (homology) {
      hairpinDG = calculateHairpinDG(homology);
      homodimerDG = calculateHomodimerDG(homology);
      // Calculate terminal 3' DG for template binding stability using NN parameters
      // This measures how well the 3' end will bind to template (not self-folding)
      const terminalResult = calculate3primeTerminalDG(homology);
      terminal3DG = terminalResult?.dG ?? -8; // Default to ideal if calculation fails
    }
    if (fullSequence) fullHairpinDG = calculateHairpinDG(fullSequence);
  } catch (e) {
    // Use neutral values if calculation fails
    terminal3DG = -8; // Default to middle of optimal range
  }

  const hairpinScore = scoreHairpin(hairpinDG);
  const homodimerScore = scoreHomodimer(homodimerDG);
  const fullHairpinScore = scoreHairpin(fullHairpinDG);
  const terminal3DGScore = scoreTerminal3DG(terminal3DG);

  // 3' composition analysis (detailed)
  // Note: score3PrimeComposition returns a NUMBER (0-1), not an object
  const composition3PrimeScore = homology
    ? score3PrimeComposition(homology, terminal3DG)
    : 1.0;

  // Build detailed 3' composition info for breakdown
  let composition3PrimeDetails = {};
  if (homology && homology.length >= 5) {
    const last5 = homology.slice(-5);
    const last2 = homology.slice(-2);
    const gcLast2 = (last2.match(/[GC]/gi) || []).length;
    const gcLast5 = (last5.match(/[GC]/gi) || []).length;
    const endsWithGC = /[GC]$/i.test(homology);
    const hasPolyAT = /[AT]{4,}/i.test(last5);

    composition3PrimeDetails = {
      last5,
      last2,
      gcLast2,
      gcLast5,
      endsWithGC,
      hasPolyAT,
      gcClampStatus: gcLast2 === 1 ? 'ideal' : gcLast2 === 2 ? 'strong' : 'weak',
    };
  }

  // G-Quadruplex analysis (critical failure predictor from calibrated weights)
  // Note: scoreGQuadruplex returns a NUMBER (0-1), not an object
  const gQuadruplexScore = homology ? scoreGQuadruplex(homology) : 1.0;

  // === GOLDEN GATE SPECIFIC CHECKS ===
  const ggSpecific = scoreGoldenGateSpecific(primer, enzyme, options);

  // === BASE SCORE CALCULATION ===
  // Weighted composite using empirically calibrated weights
  // Note: GG-specific is applied as additive layer AFTER base normalization
  const baseWeights = {
    tm: weights.tm || 0.10,
    gc: weights.gc || 0.06,
    length: weights.length || 0.05,
    gcClamp: weights.gcClamp || 0.06,
    homopolymer: weights.homopolymer || 0.05,
    terminal3DG: weights.terminal3DG || 0.18,      // Critical - empirically calibrated
    composition3Prime: weights.composition3Prime || 0.08,
    hairpin: weights.hairpin || 0.10,
    homodimer: weights.homodimer || 0.10,
    gQuadruplex: weights.gQuadruplex || 0.12,     // Critical - empirically calibrated
    fullHairpin: 0.10,                             // Account for full primer structure
  };

  // Normalize base weights to sum to 1.0
  const totalBaseWeight = Object.values(baseWeights).reduce((a, b) => a + b, 0);

  const baseComposite = (
    tmScore * baseWeights.tm +
    gcScore * baseWeights.gc +
    lengthScore * baseWeights.length +
    gcClampScore * baseWeights.gcClamp +
    homopolymerScore * baseWeights.homopolymer +
    terminal3DGScore * baseWeights.terminal3DG +
    composition3PrimeScore * baseWeights.composition3Prime +
    hairpinScore * baseWeights.hairpin +
    homodimerScore * baseWeights.homodimer +
    gQuadruplexScore * baseWeights.gQuadruplex +
    fullHairpinScore * baseWeights.fullHairpin
  ) / totalBaseWeight;

  // === GG-SPECIFIC ADDITIVE LAYER ===
  // Apply as bonus/penalty on top of base score (not normalized with it)
  // This allows GG-specific issues to adjust score without diluting critical weights
  const ggModifier = (ggSpecific.score - 0.8) * (weights.ggSpecific || 0.10);
  // ggSpecific.score ranges 0-1, baseline ~0.8 for acceptable, so modifier ranges -0.08 to +0.02

  // Final composite with GG modifier applied additively
  const composite = Math.max(0, Math.min(1, baseComposite + ggModifier));

  // === 3' COMPOSITION ANALYSIS ===
  // Detailed analysis of 3' end for design suggestions
  let threePrimeAnalysis = null;
  if (homology && homology.length >= 5) {
    const last5 = homology.slice(-5);
    const last2 = homology.slice(-2);
    const gcInLast2 = (last2.match(/[GC]/gi) || []).length;
    const gcInLast5 = (last5.match(/[GC]/gi) || []).length;
    const endsWithGC = /[GC]$/i.test(homology);

    threePrimeAnalysis = {
      last5,
      last2,
      gcInLast2,
      gcInLast5,
      endsWithGC,
      hasGCClamp: gcInLast2 >= 1,
      isStrongClamp: gcInLast2 === 2,
      terminal3DG,
    };
  }

  // Use classifyQuality for consistent tier assignment with Primer Designer
  const compositeScore = Math.round(composite * 100);
  const qualityResult = classifyQuality(compositeScore);

  return {
    composite: compositeScore,
    quality: qualityResult.tier,  // Extract string tier from object
    breakdown: {
      // Tier 1: Critical (empirically calibrated)
      terminal3DG: { score: Math.round(terminal3DGScore * 100), value: terminal3DG, weight: baseWeights.terminal3DG, tier: 'critical' },
      gQuadruplex: {
        score: Math.round(gQuadruplexScore * 100),
        weight: baseWeights.gQuadruplex,
        tier: 'critical',
      },

      // Tier 2: Important
      hairpin: { score: Math.round(hairpinScore * 100), value: hairpinDG, weight: baseWeights.hairpin, tier: 'important' },
      homodimer: { score: Math.round(homodimerScore * 100), value: homodimerDG, weight: baseWeights.homodimer, tier: 'important' },
      composition3Prime: {
        score: Math.round(composition3PrimeScore * 100),
        details: composition3PrimeDetails,
        weight: baseWeights.composition3Prime,
        tier: 'important',
      },

      // Tier 3: Standard
      tm: { score: Math.round(tmScore * 100), value: tm, weight: baseWeights.tm, tier: 'standard' },
      gc: { score: Math.round(gcScore * 100), value: gc, weight: baseWeights.gc, tier: 'standard' },
      length: { score: Math.round(lengthScore * 100), value: homologyLength, optimal: '18-25bp', weight: baseWeights.length, tier: 'standard' },
      gcClamp: { score: Math.round(gcClampScore * 100), weight: baseWeights.gcClamp, tier: 'standard' },
      homopolymer: { score: Math.round(homopolymerScore * 100), weight: baseWeights.homopolymer, tier: 'standard' },
      fullPrimerHairpin: { score: Math.round(fullHairpinScore * 100), value: fullHairpinDG, weight: baseWeights.fullHairpin, tier: 'standard' },

      // GG-Specific (additive layer)
      ggSpecific: {
        score: Math.round(ggSpecific.score * 100),
        issues: ggSpecific.issues,
        modifier: Math.round(ggModifier * 100),
        note: 'Applied as additive modifier to base score',
      },
    },
    // Include base score before GG modifier for transparency
    baseScore: Math.round(baseComposite * 100),
    ggModifier: Math.round(ggModifier * 100),
  };
}

/**
 * Score Golden Gate primer pair including heterodimer analysis
 * @param {Object} primers - Primer pair from designGoldenGatePrimers
 * @param {Object} options - Scoring options
 * @returns {Object} Pair quality score
 */
export function scoreGoldenGatePrimerPair(primers, options = {}) {
  const {
    enzyme = 'BsaI',
    tmDiffWarningThreshold = GG_OPTIMIZER_DEFAULTS.tmDiffWarningThreshold,
    tmDiffPenaltyMultiplier = GG_OPTIMIZER_DEFAULTS.tmDiffPenaltyMultiplier,
    excellentPrimerScore = GG_OPTIMIZER_DEFAULTS.excellentPrimerScore,
    goodPrimerScore = GG_OPTIMIZER_DEFAULTS.goodPrimerScore,
    acceptablePrimerScore = GG_OPTIMIZER_DEFAULTS.acceptablePrimerScore,
  } = options;

  // Handle null/undefined primers
  if (!primers || !primers.forward || !primers.reverse) {
    return {
      forward: { composite: 0, quality: 'poor', error: 'Missing primer' },
      reverse: { composite: 0, quality: 'poor', error: 'Missing primer' },
      pair: {
        score: 0,
        quality: 'poor',
        error: 'Invalid primer pair',
      },
    };
  }

  const fwdScore = scoreGoldenGatePrimer(primers.forward, { enzyme, ...options });
  const revScore = scoreGoldenGatePrimer(primers.reverse, { enzyme, ...options });

  // Heterodimer analysis
  let heterodimerDG = 0;
  let fullHeterodimerDG = 0;

  try {
    const fwdHomology = primers.forward.structure?.homology || primers.forward.homologyRegion || '';
    const revHomology = primers.reverse.homologyRegion || '';

    if (fwdHomology && revHomology) {
      heterodimerDG = calculateHeterodimerDG(fwdHomology, revHomology);
    }
    if (primers.forward.sequence && primers.reverse.sequence) {
      fullHeterodimerDG = calculateHeterodimerDG(
        primers.forward.sequence,
        primers.reverse.sequence
      );
    }
  } catch (e) {
    // Use neutral values if calculation fails
  }

  const heterodimerScore = scoreHeterodimer(fullHeterodimerDG);

  // Tm difference
  const fwdTm = primers.forward.tm || 60;
  const revTm = primers.reverse.tm || 60;
  const tmDiff = Math.abs(fwdTm - revTm);
  const tmDiffPenalty = tmDiff > tmDiffWarningThreshold
    ? Math.min(20, (tmDiff - tmDiffWarningThreshold) * tmDiffPenaltyMultiplier)
    : 0;

  // Combined pair score
  const avgPrimerScore = (fwdScore.composite + revScore.composite) / 2;
  const pairScore = Math.max(0, avgPrimerScore - tmDiffPenalty - (heterodimerDG < -8 ? 10 : 0));
  const roundedPairScore = Math.round(pairScore);

  // Use classifyQuality for consistent tier assignment with Primer Designer
  const pairQualityResult = classifyQuality(roundedPairScore);

  return {
    forward: fwdScore,
    reverse: revScore,
    pair: {
      score: roundedPairScore,
      quality: pairQualityResult.tier,  // Extract string tier from object
      heterodimer: {
        homologyRegions: heterodimerDG,
        fullPrimers: fullHeterodimerDG,
        score: Math.round(heterodimerScore * 100),
      },
      tmDifference: {
        value: tmDiff,
        penalty: tmDiffPenalty,
        acceptable: tmDiff <= tmDiffWarningThreshold,
      },
    },
  };
}

// ============================================================================
// AUTOMATIC PRIMER QUALITY OPTIMIZATION
// ============================================================================

/**
 * Calculate Tm using nearest-neighbor method (simplified Q5 approximation)
 * @param {string} seq - Primer sequence
 * @returns {number} Melting temperature in °C
 */
function calculateTmSimple(seq) {
  if (!seq || seq.length < 8) return 40;
  const gc = (seq.match(/[GC]/gi) || []).length;
  const at = seq.length - gc;
  // Wallace rule with salt correction approximation
  if (seq.length < 14) {
    return 2 * at + 4 * gc;
  }
  // Longer primers: use nearest-neighbor approximation
  return 64.9 + 41 * (gc - 16.4) / seq.length;
}

/**
 * Analyze 3' end quality of a primer
 * @param {string} seq - Primer sequence
 * @returns {Object} Analysis results
 */
function analyze3PrimeEnd(seq) {
  if (!seq || seq.length < 5) {
    return { quality: 'poor', gcClamp: 0, issues: ['Sequence too short'] };
  }

  const last2 = seq.slice(-2).toUpperCase();
  const last5 = seq.slice(-5).toUpperCase();
  const gcInLast2 = (last2.match(/[GC]/gi) || []).length;
  const gcInLast5 = (last5.match(/[GC]/gi) || []).length;

  const issues = [];
  let score = 100;

  // GC clamp check (want 1-2 G/C in last 2 bases)
  if (gcInLast2 === 0) {
    issues.push('No GC clamp (last 2 bases are AT)');
    score -= 30;
  } else if (gcInLast2 === 2 && last2 === 'GG') {
    // GG at 3' end can cause mispriming
    issues.push('Strong GG clamp may cause mispriming');
    score -= 10;
  }

  // Avoid 3+ G/C at 3' end (too strong)
  if (gcInLast5 >= 4) {
    issues.push('High GC in last 5 bases');
    score -= 15;
  }

  // Avoid homopolymers at 3' end
  if (/(.)\1{2,}$/.test(seq)) {
    issues.push('Homopolymer at 3\' end');
    score -= 25;
  }

  // Avoid ending in T (weak binding)
  if (seq.slice(-1).toUpperCase() === 'T') {
    issues.push('Ends in T (weak)');
    score -= 10;
  }

  // Determine quality tier
  let quality;
  if (score >= 90) quality = 'excellent';
  else if (score >= 70) quality = 'good';
  else if (score >= 50) quality = 'acceptable';
  else quality = 'poor';

  return {
    quality,
    score,
    gcClamp: gcInLast2,
    gcInLast5,
    last2,
    last5,
    issues,
  };
}

/**
 * Try to optimize a primer by adjusting homology length and position
 *
 * @param {string} templateSeq - Full template sequence
 * @param {boolean} isForward - True for forward primer, false for reverse
 * @param {number} currentLength - Current homology length
 * @param {Object} options - Optimization options
 * @returns {Object} Optimized homology region info
 */
function findOptimalHomologyWithQuality(templateSeq, isForward, currentLength, options = {}) {
  const {
    targetTm = 60,
    minLength = 15,
    maxLength = 30,
    maxShift = 3,  // How much to shift binding position
  } = options;

  const seq = templateSeq.toUpperCase();
  const candidates = [];

  // Try different lengths around current
  for (let len = Math.max(minLength, currentLength - 3); len <= Math.min(maxLength, currentLength + 5); len++) {
    // Try different starting positions (shift)
    for (let shift = 0; shift <= maxShift; shift++) {
      let homology;
      if (isForward) {
        // Forward primer: take from start, can shift right
        if (shift + len > seq.length) continue;
        homology = seq.slice(shift, shift + len);
      } else {
        // Reverse primer: take from end, can shift left
        if (seq.length - shift - len < 0) continue;
        homology = seq.slice(seq.length - shift - len, seq.length - shift);
      }

      if (homology.length < minLength) continue;

      // Calculate quality metrics
      const tm = calculateTmSimple(homology);
      const gc = (homology.match(/[GC]/gi) || []).length / homology.length;
      const analysis3Prime = analyze3PrimeEnd(homology);

      // Score this candidate
      let score = 100;

      // Tm scoring (optimal 58-62°C)
      const tmDiff = Math.abs(tm - targetTm);
      if (tmDiff > 5) score -= tmDiff * 3;
      else if (tmDiff > 2) score -= tmDiff * 1.5;

      // GC content (optimal 40-60%)
      if (gc < 0.35 || gc > 0.65) score -= 15;
      else if (gc < 0.40 || gc > 0.60) score -= 5;

      // 3' end quality (major factor)
      score += (analysis3Prime.score - 70) * 0.5;

      // Length penalty (prefer shorter if quality is same)
      if (len > 25) score -= (len - 25) * 2;

      // Shift penalty (prefer no shift)
      score -= shift * 5;

      candidates.push({
        homology,
        length: len,
        shift,
        tm,
        gc: gc * 100,
        analysis3Prime,
        score,
        isForward,
      });
    }
  }

  // Sort by score and return best
  candidates.sort((a, b) => b.score - a.score);

  if (candidates.length === 0) {
    // Fallback to original
    const fallbackLen = Math.min(currentLength, seq.length);
    const homology = isForward ? seq.slice(0, fallbackLen) : seq.slice(-fallbackLen);
    return {
      homology,
      length: fallbackLen,
      shift: 0,
      tm: calculateTmSimple(homology),
      gc: (homology.match(/[GC]/gi) || []).length / homology.length * 100,
      analysis3Prime: analyze3PrimeEnd(homology),
      score: 50,
      optimized: false,
    };
  }

  const best = candidates[0];
  return {
    ...best,
    optimized: best.shift > 0 || best.length !== currentLength,
    alternatives: candidates.slice(1, 4),
  };
}

/**
 * Optimize a Golden Gate primer for better quality
 *
 * @param {Object} primer - Original primer object
 * @param {string} templateSeq - Template sequence
 * @param {boolean} isForward - True for forward primer
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Options
 * @returns {Object} Optimized primer
 */
function optimizeGoldenGatePrimer(primer, templateSeq, isForward, enzyme, options = {}) {
  const structure = primer?.structure || {};
  const currentHomology = structure.homology || primer?.homologyRegion || '';
  const extra = structure.extra || 'GGTGCG';
  const recognitionSite = structure.recognitionSite || structure.bsaISite || '';
  const spacer = structure.spacer || 'A';
  const overhang = structure.overhang || '';

  // Find optimized homology
  const optimized = findOptimalHomologyWithQuality(
    templateSeq,
    isForward,
    currentHomology.length,
    options
  );

  // Build new primer sequence
  const newHomology = isForward ? optimized.homology : reverseComplement(optimized.homology);
  const newSequence = extra + recognitionSite + spacer + overhang + newHomology;

  return {
    ...primer,
    sequence: newSequence,
    length: newSequence.length,
    tm: optimized.tm,
    gc: optimized.gc,
    homologyRegion: newHomology,
    homologyLength: optimized.length,
    structure: {
      ...structure,
      homology: newHomology,
    },
    optimization: {
      original: {
        homology: currentHomology,
        length: currentHomology.length,
      },
      optimized: {
        homology: newHomology,
        length: optimized.length,
        shift: optimized.shift,
        score: optimized.score,
        analysis3Prime: optimized.analysis3Prime,
      },
      wasOptimized: optimized.optimized,
      reason: optimized.optimized
        ? `Adjusted homology for better 3' end quality (shift: ${optimized.shift}, length: ${currentHomology.length} → ${optimized.length})`
        : 'Original homology was optimal',
    },
  };
}

// ============================================================================
// OPTIMIZED PRIMER DESIGN
// ============================================================================

/**
 * Rebuild primer sequence with new flanking
 * @param {Object} primer - Original primer
 * @param {string} newFlanking - New flanking sequence
 * @param {string} enzyme - Enzyme name
 * @returns {string} New primer sequence
 */
function rebuildPrimerSequence(primer, newFlanking, enzyme) {
  const structure = primer?.structure || {};
  const recognitionSite = structure.recognitionSite || structure.bsaISite || '';
  const spacer = structure.spacer || 'A';
  const overhang = structure.overhang || '';
  const homology = structure.homology || '';

  return newFlanking + recognitionSite + spacer + overhang + homology;
}

/**
 * Design optimized Golden Gate primers with state-of-the-art enhancements
 *
 * Improvements over standard design:
 * 1. Uses optimal 6bp flanking sequences (NEB recommendation)
 * 2. Selects flanking dynamically based on template to avoid mispriming
 * 3. Prioritizes default flanking unless it causes problems
 * 4. Scores primer quality comprehensively
 *
 * @param {string} targetSeq - Target sequence to amplify
 * @param {string} leftOverhang - 5' overhang
 * @param {string} rightOverhang - 3' overhang
 * @param {Object} options - Design options
 * @param {string} options.customFlanking - User-specified flanking sequence (optional)
 * @param {boolean} options.checkMispriming - Enable mispriming check (default: true)
 * @returns {Object} Optimized primer design
 */
export function designOptimizedGoldenGatePrimers(targetSeq, leftOverhang, rightOverhang, options = {}) {
  const {
    enzyme = 'BsaI',
    customFlanking = null,
    checkMispriming = true,
    autoOptimizePrimers = true,  // NEW: Enable auto-optimization of poor quality primers
    targetTm = 60,
    ...restOptions
  } = options;

  // Input validation
  validateDNASequence(targetSeq, 'Target sequence');
  validateDNASequence(leftOverhang, 'Left overhang');
  validateDNASequence(rightOverhang, 'Right overhang');

  // Step 1: Design primers with original function (will use default 2bp flanking)
  let basePrimers = designGoldenGatePrimers(targetSeq, leftOverhang, rightOverhang, {
    enzyme,
    ...restOptions,
  });

  // Step 2: Select optimal flanking for forward primer context
  // Pass template for mispriming check if enabled
  const fwdHomologyStart = basePrimers.forward.structure?.homology?.slice(0, 10) || targetSeq.slice(0, 10);
  const fwdFlankingResult = selectOptimalFlankingSequence(
    enzyme,
    leftOverhang,
    fwdHomologyStart,
    {
      ...options,
      template: checkMispriming ? targetSeq : null,
      customFlanking: customFlanking,
    }
  );

  // Step 3: Select optimal flanking for reverse primer context
  const revHomologyStart = basePrimers.reverse.structure?.homology?.slice(0, 10) ||
                          reverseComplement(targetSeq.slice(-10));
  const revFlankingResult = selectOptimalFlankingSequence(
    enzyme,
    reverseComplement(rightOverhang),
    revHomologyStart,
    {
      ...options,
      template: checkMispriming ? targetSeq : null,
      customFlanking: customFlanking,
    }
  );

  // Step 4: Build initial optimized primers with flanking
  let fwdPrimer = {
    ...basePrimers.forward,
    structure: {
      ...basePrimers.forward.structure,
      extra: fwdFlankingResult.best.flanking,
    },
  };
  fwdPrimer.sequence = rebuildPrimerSequence(fwdPrimer, fwdFlankingResult.best.flanking, enzyme);
  fwdPrimer.length = fwdPrimer.sequence.length;

  let revPrimer = {
    ...basePrimers.reverse,
    structure: {
      ...basePrimers.reverse.structure,
      extra: revFlankingResult.best.flanking,
    },
  };
  revPrimer.sequence = rebuildPrimerSequence(revPrimer, revFlankingResult.best.flanking, enzyme);
  revPrimer.length = revPrimer.sequence.length;

  // Step 5: Score initial primers
  let initialScore = scoreGoldenGatePrimerPair({ forward: fwdPrimer, reverse: revPrimer }, { enzyme, ...options });

  // Step 6: AUTO-OPTIMIZE if primers are poor quality
  let fwdOptimization = null;
  let revOptimization = null;

  if (autoOptimizePrimers) {
    // Check forward primer quality
    const fwdQuality = initialScore.forward?.quality;
    if (fwdQuality === 'poor' || fwdQuality === 'acceptable') {
      const optimizedFwd = optimizeGoldenGatePrimer(fwdPrimer, targetSeq, true, enzyme, { targetTm });
      if (optimizedFwd.optimization?.wasOptimized) {
        fwdPrimer = optimizedFwd;
        fwdOptimization = optimizedFwd.optimization;
      }
    }

    // Check reverse primer quality
    const revQuality = initialScore.reverse?.quality;
    if (revQuality === 'poor' || revQuality === 'acceptable') {
      const optimizedRev = optimizeGoldenGatePrimer(revPrimer, targetSeq, false, enzyme, { targetTm });
      if (optimizedRev.optimization?.wasOptimized) {
        revPrimer = optimizedRev;
        revOptimization = optimizedRev.optimization;
      }
    }
  }

  // Create final optimized primer objects
  const optimizedPrimers = {
    forward: fwdPrimer,
    reverse: revPrimer,
    pcr: basePrimers.pcr,
    fidelity: basePrimers.fidelity,
    enzyme,
    warnings: basePrimers.warnings || [],
  };

  // Step 7: Re-score after optimization
  const qualityScore = scoreGoldenGatePrimerPair(optimizedPrimers, { enzyme, ...options });

  // Add optimization metadata with selection reasoning
  optimizedPrimers.optimization = {
    flanking: {
      forward: {
        original: basePrimers.forward.structure?.extra || 'GG',
        optimized: fwdFlankingResult.best.flanking,
        score: fwdFlankingResult.best.score,
        quality: fwdFlankingResult.best.quality,
        selectionReason: fwdFlankingResult.selectionReason,
        message: fwdFlankingResult.message,
        mispriming: fwdFlankingResult.best.mispriming,
      },
      reverse: {
        original: basePrimers.reverse.structure?.extra || 'GG',
        optimized: revFlankingResult.best.flanking,
        score: revFlankingResult.best.score,
        quality: revFlankingResult.best.quality,
        selectionReason: revFlankingResult.selectionReason,
        message: revFlankingResult.message,
        mispriming: revFlankingResult.best.mispriming,
      },
    },
    homology: {
      forward: fwdOptimization ? {
        wasOptimized: true,
        originalLength: fwdOptimization.originalLength,
        newLength: fwdOptimization.newLength,
        lengthChange: fwdOptimization.newLength - fwdOptimization.originalLength,
        shift: fwdOptimization.shift,
        originalTm: fwdOptimization.originalTm,
        newTm: fwdOptimization.newTm,
        originalQuality: fwdOptimization.originalQuality,
        newQuality: fwdOptimization.newQuality,
        reason: fwdOptimization.reason,
      } : { wasOptimized: false },
      reverse: revOptimization ? {
        wasOptimized: true,
        originalLength: revOptimization.originalLength,
        newLength: revOptimization.newLength,
        lengthChange: revOptimization.newLength - revOptimization.originalLength,
        shift: revOptimization.shift,
        originalTm: revOptimization.originalTm,
        newTm: revOptimization.newTm,
        originalQuality: revOptimization.originalQuality,
        newQuality: revOptimization.newQuality,
        reason: revOptimization.reason,
      } : { wasOptimized: false },
    },
    quality: qualityScore,
  };

  // Add warnings if flanking was changed from default due to issues
  if (fwdFlankingResult.selectionReason !== 'default_acceptable' &&
      fwdFlankingResult.selectionReason !== 'user_specified') {
    optimizedPrimers.warnings.push(`Forward primer: ${fwdFlankingResult.message}`);
  }
  if (revFlankingResult.selectionReason !== 'default_acceptable' &&
      revFlankingResult.selectionReason !== 'user_specified') {
    optimizedPrimers.warnings.push(`Reverse primer: ${revFlankingResult.message}`);
  }

  // Add info about homology optimization if performed
  if (fwdOptimization?.wasOptimized) {
    optimizedPrimers.warnings.push(
      `Forward primer auto-optimized: homology ${fwdOptimization.originalLength}→${fwdOptimization.newLength}bp ` +
      `(Tm: ${fwdOptimization.originalTm?.toFixed(1)}→${fwdOptimization.newTm?.toFixed(1)}°C, ` +
      `${fwdOptimization.originalQuality}→${fwdOptimization.newQuality})`
    );
  }
  if (revOptimization?.wasOptimized) {
    optimizedPrimers.warnings.push(
      `Reverse primer auto-optimized: homology ${revOptimization.originalLength}→${revOptimization.newLength}bp ` +
      `(Tm: ${revOptimization.originalTm?.toFixed(1)}→${revOptimization.newTm?.toFixed(1)}°C, ` +
      `${revOptimization.originalQuality}→${revOptimization.newQuality})`
    );
  }

  return optimizedPrimers;
}

// ============================================================================
// INTEGRATED ASSEMBLY OPTIMIZATION
// ============================================================================

/**
 * Generate recommendations based on assembly analysis
 * @param {Object} fidelity - Fidelity analysis
 * @param {Array} primerScores - Primer quality scores
 * @param {Object} config - Configuration options
 * @param {Array} primers - Optional primer objects with optimization info
 * @returns {Array} Recommendations
 */
function generateRecommendations(fidelity, primerScores, config = {}, primers = []) {
  const {
    minAcceptableFidelity = GG_OPTIMIZER_DEFAULTS.minAcceptableFidelity,
  } = config;

  const recommendations = [];

  // Count auto-optimized primers
  const optimizedPrimers = primers.filter(p =>
    p?.optimization?.homology?.forward?.wasOptimized ||
    p?.optimization?.homology?.reverse?.wasOptimized
  );

  // Fidelity recommendations
  if (fidelity && fidelity.gtAdjustedFidelity < minAcceptableFidelity) {
    recommendations.push({
      type: 'fidelity',
      severity: 'warning',
      title: 'Low Assembly Fidelity',
      message: `Assembly fidelity ${fidelity.gtAdjustedFidelityPercent} is below ${minAcceptableFidelity * 100}%. Consider reducing number of parts or using Monte Carlo optimization for overhang selection.`,
      action: 'Consider using optimizeOverhangSet() for better overhang selection',
    });
  }

  // G:T mismatch warnings
  if (fidelity?.gtRisks && fidelity.gtRisks.length > 0) {
    const criticalRisks = fidelity.gtRisks.filter(r => r.risk === 'critical' || r.risk === 'high');
    if (criticalRisks.length > 0) {
      recommendations.push({
        type: 'gt_mismatch',
        severity: 'warning',
        title: 'G:T Mismatch Risks Detected',
        message: `${criticalRisks.length} high-risk G:T mismatch pair(s) detected. T4 ligase may accept these mis-ligations.`,
        details: criticalRisks,
        action: 'Screen additional colonies or replace problematic overhangs',
      });
    } else {
      recommendations.push({
        type: 'gt_mismatch',
        severity: 'info',
        title: 'Minor G:T Mismatch Risks',
        message: `${fidelity.gtRisks.length} low-risk G:T mismatch pair(s) detected. Usually not problematic.`,
        details: fidelity.gtRisks,
      });
    }
  }

  // Primer quality recommendations
  if (primerScores && primerScores.length > 0) {
    const poorPrimers = primerScores.filter(s => s?.pair?.quality === 'poor');
    const acceptablePrimers = primerScores.filter(s => s?.pair?.quality === 'acceptable');
    const goodPrimers = primerScores.filter(s => s?.pair?.quality === 'good' || s?.pair?.quality === 'excellent');

    // Report auto-optimization results
    if (optimizedPrimers.length > 0) {
      const improvedCount = optimizedPrimers.filter(p => {
        const fwdImproved = p?.optimization?.homology?.forward?.wasOptimized &&
          p?.optimization?.homology?.forward?.newQuality !== p?.optimization?.homology?.forward?.originalQuality;
        const revImproved = p?.optimization?.homology?.reverse?.wasOptimized &&
          p?.optimization?.homology?.reverse?.newQuality !== p?.optimization?.homology?.reverse?.originalQuality;
        return fwdImproved || revImproved;
      }).length;

      recommendations.push({
        type: 'primer_optimization',
        severity: 'info',
        title: 'Auto-Optimization Applied',
        message: `${optimizedPrimers.length} primer pair(s) were auto-optimized by adjusting homology length/position.` +
          (improvedCount > 0 ? ` ${improvedCount} showed quality improvement.` : ''),
      });
    }

    // Only warn about poor primers that couldn't be improved
    if (poorPrimers.length > 0) {
      const stillPoorAfterOpt = poorPrimers.length;
      if (stillPoorAfterOpt > 0) {
        recommendations.push({
          type: 'primer_quality',
          severity: 'warning',
          title: 'Poor Primer Quality',
          message: `${stillPoorAfterOpt} primer pair(s) still scored as "poor" after auto-optimization. ` +
            `May require manual redesign or different assembly strategy.`,
          action: 'Consider splitting into smaller fragments or using different junction points',
        });
      }
    }

    if (acceptablePrimers.length > 0 && poorPrimers.length === 0) {
      recommendations.push({
        type: 'primer_quality',
        severity: 'info',
        title: 'Acceptable Primer Quality',
        message: `${acceptablePrimers.length} primer pair(s) scored as "acceptable". Should work for most applications.`,
      });
    }

    // Positive feedback for good primers
    if (goodPrimers.length > 0 && poorPrimers.length === 0 && acceptablePrimers.length === 0) {
      recommendations.push({
        type: 'primer_quality',
        severity: 'success',
        title: 'Good Primer Quality',
        message: `All ${goodPrimers.length} primer pair(s) scored as "good" or better.`,
      });
    }
  }

  return recommendations;
}

/**
 * Design fully optimized Golden Gate Assembly
 *
 * Main entry point combining all optimizations:
 * 1. Auto-optimize overhangs for maximum fidelity
 * 2. Use optimal 6bp flanking sequences
 * 3. Score primer quality comprehensively
 * 4. Detect and report G:T mismatch risks
 *
 * @param {Array} parts - Array of part objects with 'sequence' or 'seq' property (or string sequences)
 * @param {Object} options - Assembly options
 * @returns {Object} Complete optimized assembly design
 */
export function designOptimizedGoldenGateAssembly(parts, options = {}) {
  const {
    enzyme = 'BsaI',
    overhangs = null,
    autoOptimize = true,
    circular = true,               // NEW: circular vs linear assembly
    targetFidelity = GG_OPTIMIZER_DEFAULTS.targetFidelity,
    requiredOverhangIndices = [],
    requiredPatterns = {},         // NEW: pattern requirements for overhangs
  } = options;

  // Input validation
  if (!Array.isArray(parts) || parts.length === 0) {
    return {
      enzyme,
      numParts: 0,
      error: 'No parts provided',
      overhangs: { initial: [], final: [], optimization: null },
      primers: [],
      quality: { assemblyFidelity: 0, assemblyFidelityPercent: '0.0%' },
      recommendations: [{ type: 'error', severity: 'critical', message: 'No parts provided' }],
      summary: { overallQuality: 'poor' },
    };
  }

  // Step 1: Get initial overhangs - use dynamic search for 100% fidelity
  let initialOverhangs = overhangs;
  if (!initialOverhangs) {
    // For circular: numJunctions = numParts + 1 (vector flanking overhangs)
    // For linear: numJunctions = numParts (no flanking needed)
    const numJunctions = circular ? parts.length + 1 : parts.length;

    // Use dynamic search to find 100% fidelity overhang set
    // This searches the ligation matrix for combinations with zero cross-reactivity
    const optimalResult = findOptimalOverhangSet(numJunctions, enzyme, {
      minCorrectFreq: 300,
      requiredOverhangs: requiredPatterns || [],
    });

    initialOverhangs = optimalResult.overhangs;

    // Log if we couldn't achieve 100% fidelity (for debugging)
    if (!optimalResult.isPerfect) {
      console.warn(
        `Could not find 100% fidelity set for ${numJunctions} junctions. ` +
        `Best found: ${optimalResult.fidelityPercent}`
      );
    }
  }

  // Step 2: Auto-optimize overhangs if enabled
  // Note: Dynamic search already finds 100% fidelity when possible, so optimization
  // typically won't improve. But we still run it in case custom overhangs
  // were provided or for edge cases.
  let overhangOptimization = null;
  let finalOverhangs = initialOverhangs;

  if (autoOptimize && initialOverhangs.length > 0) {
    overhangOptimization = autoOptimizeOverhangs(initialOverhangs, {
      enzyme,
      requiredIndices: requiredOverhangIndices,
      requiredPatterns,
      minFidelity: targetFidelity,
    });
    finalOverhangs = overhangOptimization.optimized;
  }

  // Step 3: Design optimized primers for each part
  const primers = [];
  const errors = [];

  for (let idx = 0; idx < parts.length; idx++) {
    const part = parts[idx];
    const sequence = typeof part === 'string' ? part : (part?.sequence || part?.seq);

    // Determine overhangs for this part
    const leftOH = finalOverhangs[idx];
    const rightOH = finalOverhangs[idx + 1] || finalOverhangs[0]; // Wrap for circular

    if (!sequence || !leftOH || !rightOH) {
      errors.push({ index: idx, error: 'Missing sequence or overhang' });
      primers.push(null);
      continue;
    }

    try {
      const primer = designOptimizedGoldenGatePrimers(sequence, leftOH, rightOH, {
        enzyme,
        ...options.primerOptions,
      });
      primers.push(primer);
    } catch (e) {
      errors.push({ index: idx, error: e.message });
      primers.push(null);
    }
  }

  // Step 4: Score all primer pairs (filter out nulls)
  const validPrimers = primers.filter(p => p !== null);
  const primerScores = validPrimers.map(p => scoreGoldenGatePrimerPair(p, { enzyme }));

  // Step 5: Calculate overall assembly quality using ligation frequency matrix
  // This calculates fidelity based on actual cross-reactivity between overhangs in this specific set
  const experimentalFidelity = calculateExperimentalFidelity(finalOverhangs, enzyme);

  const overallFidelity = {
    baseFidelity: experimentalFidelity.assemblyFidelity,
    baseFidelityPercent: experimentalFidelity.assemblyFidelityPercent,
    gtAdjustedFidelity: experimentalFidelity.assemblyFidelity, // Matrix already accounts for all mis-ligations
    gtAdjustedFidelityPercent: experimentalFidelity.assemblyFidelityPercent,
    gtPenalty: 1.0,
    gtPenaltyPercent: '0.0%',
    gtRisks: [], // Matrix-based calculation handles cross-reactivity
    junctionFidelities: experimentalFidelity.junctions || [],
    hasGTRisks: false,
    source: 'experimental-matrix',
    warnings: experimentalFidelity.warnings || [],
    lowestFidelity: experimentalFidelity.lowestFidelity,
  };

  // Safe average calculation (avoid division by zero)
  const avgPrimerScore = primerScores.length > 0
    ? primerScores.reduce((sum, s) => sum + (s?.pair?.score || 0), 0) / primerScores.length
    : 0;

  // Step 6: Generate recommendations (pass primers for optimization tracking)
  const recommendations = generateRecommendations(overallFidelity, primerScores, options, primers);

  // Add errors to recommendations
  errors.forEach(err => {
    recommendations.push({
      type: 'error',
      severity: 'warning',
      title: `Part ${err.index} Error`,
      message: err.error,
    });
  });

  return {
    // Assembly metadata
    enzyme,
    numParts: parts.length,
    circular,

    // Overhangs (optimized)
    overhangs: {
      initial: initialOverhangs,
      final: finalOverhangs,
      optimization: overhangOptimization,
    },

    // Primers for each part
    primers,

    // Quality metrics
    quality: {
      assemblyFidelity: overallFidelity.gtAdjustedFidelity,
      assemblyFidelityPercent: overallFidelity.gtAdjustedFidelityPercent,
      baseFidelity: overallFidelity.baseFidelity,
      baseFidelityPercent: overallFidelity.baseFidelityPercent,
      gtPenalty: overallFidelity.gtPenalty,
      averagePrimerScore: Math.round(avgPrimerScore),
      primerScores,
      gtRisks: overallFidelity.gtRisks,
      junctionFidelities: overallFidelity.junctionFidelities,
    },

    // Recommendations
    recommendations,

    // Summary
    summary: {
      overallQuality: overallFidelity.gtAdjustedFidelity >= GG_OPTIMIZER_DEFAULTS.excellentFidelityThreshold ? 'excellent' :
                      overallFidelity.gtAdjustedFidelity >= GG_OPTIMIZER_DEFAULTS.goodFidelityThreshold ? 'good' :
                      overallFidelity.gtAdjustedFidelity >= GG_OPTIMIZER_DEFAULTS.acceptableFidelityThreshold ? 'acceptable' : 'poor',
      fidelity: overallFidelity.gtAdjustedFidelityPercent,
      primerQuality: avgPrimerScore >= GG_OPTIMIZER_DEFAULTS.excellentPrimerScore ? 'excellent' :
                     avgPrimerScore >= GG_OPTIMIZER_DEFAULTS.goodPrimerScore ? 'good' :
                     avgPrimerScore >= GG_OPTIMIZER_DEFAULTS.acceptablePrimerScore ? 'acceptable' : 'poor',
      warningCount: recommendations.filter(r => r.severity === 'warning').length,
      infoCount: recommendations.filter(r => r.severity === 'info').length,
      errorCount: errors.length,
    },
  };
}

// ============================================================================
// UI-COMPATIBLE ADAPTER
// ============================================================================

/**
 * Generate a protocol for Golden Gate assembly
 * @param {Array} parts - Array of designed parts
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Protocol object
 */
function generateProtocol(parts, enzyme) {
  const totalVolume = 20; // µL
  const dnaPerPart = 2; // µL of 75 ng/µL
  const enzymeVol = 2;
  const bufferVol = 2;
  const waterVol = totalVolume - (parts.length * dnaPerPart) - enzymeVol - bufferVol;

  return {
    title: `Golden Gate Assembly Protocol (${parts.length} parts)`,
    steps: [
      {
        step: 1,
        title: 'Prepare DNA fragments',
        details: [
          'PCR amplify each fragment using the designed primers',
          'Gel purify or PCR cleanup each product',
          'Dilute to 75 ng/µL in water or TE buffer',
        ],
      },
      {
        step: 2,
        title: 'Setup reaction',
        details: [
          `Add ${waterVol.toFixed(1)} µL nuclease-free water`,
          `Add ${bufferVol} µL T4 DNA Ligase Buffer (10X)`,
          `Add ${dnaPerPart} µL of each DNA fragment (${parts.length} fragments)`,
          `Add ${enzymeVol} µL ${enzyme}-HF + T4 DNA Ligase mix`,
        ],
      },
      {
        step: 3,
        title: 'Thermocycling',
        details: [
          '30 cycles of: 37°C for 1 min, 16°C for 1 min',
          'Final digestion: 50°C for 5 min',
          'Heat inactivation: 80°C for 5 min',
          'Hold at 4°C',
        ],
      },
      {
        step: 4,
        title: 'Transform',
        details: [
          'Transform 2-5 µL into competent cells',
          'Plate on appropriate selection media',
          'Incubate overnight at 37°C',
        ],
      },
    ],
    notes: [
      'Use equimolar amounts of each fragment for best results',
      'Expected colony count: 100-1000 per µL transformed',
      `Assembly fidelity depends on overhang selection`,
    ],
  };
}

/**
 * Adapter function that wraps designOptimizedGoldenGateAssembly
 * and returns a result structure compatible with the UI component.
 *
 * This allows the UI to use the optimized primer design while
 * maintaining backwards compatibility with the expected data structure.
 *
 * @param {Array} parts - Array of parts with { id, seq, type }
 * @param {Object} options - Options for assembly design
 * @returns {Object} Result in UI-compatible format
 */
export function designOptimizedGoldenGateAssemblyForUI(parts, options = {}) {
  const {
    enzyme = 'BsaI',
    circular = true,
    useStandardOverhangs = true,
    autoOptimize = true,
    autoDomestication = false, // If true, parts have been auto-domesticated
  } = options;

  // Call the optimized assembly function
  const optimizedResult = designOptimizedGoldenGateAssembly(parts, {
    enzyme,
    circular,
    autoOptimize,
    overhangs: useStandardOverhangs ? null : options.customOverhangs,
  });

  // Handle error case
  if (optimizedResult.error) {
    return {
      enzyme,
      parts: [],
      overhangs: [],
      circular,
      assembledSequence: '',
      assembledLength: 0,
      fidelity: {
        individual: [],
        overall: 0,
        percentage: '0.0%',
      },
      warnings: [optimizedResult.error],
      internalSiteIssues: [],
      hasInternalSites: false,
      protocol: null,
      // Include optimized data for advanced users
      _optimizedData: optimizedResult,
    };
  }

  // Get final overhangs
  const finalOverhangs = optimizedResult.overhangs.final || optimizedResult.overhangs.initial;

  // Default primer structure for error cases
  const defaultPrimerStructure = {
    extra: '',
    bsaISite: '',
    recognitionSite: '',
    spacer: '',
    overhang: '',
    homology: '',
  };

  const defaultPrimer = {
    sequence: '',
    length: 0,
    tm: 0,
    gc: 0,
    homologyRegion: '',
    homologyLength: 0,
    overhang: '',
    structure: defaultPrimerStructure,
  };

  // Get primer scores from optimized result
  const primerScores = optimizedResult.quality?.primerScores || [];

  // Transform primers to designed parts structure
  const designedParts = parts.map((part, i) => {
    const primerData = optimizedResult.primers[i];
    const leftOH = finalOverhangs[i] || '';
    const rightOH = finalOverhangs[i + 1] || '';
    const qualityScore = primerScores[i] || null;

    // Ensure forward and reverse have proper structure even when primer generation fails
    const forward = primerData?.forward ? {
      ...primerData.forward,
      structure: primerData.forward.structure || defaultPrimerStructure,
    } : { ...defaultPrimer };

    const reverse = primerData?.reverse ? {
      ...primerData.reverse,
      structure: primerData.reverse.structure || defaultPrimerStructure,
    } : { ...defaultPrimer };

    // Extract thermodynamic data from quality scores
    const fwdQuality = qualityScore?.forward || {};
    const revQuality = qualityScore?.reverse || {};
    const pairQuality = qualityScore?.pair || {};

    return {
      ...part,
      index: i + 1,
      leftOverhang: leftOH,
      rightOverhang: rightOH,
      primers: {
        forward: {
          ...forward,
          // Add quality analysis data
          qualityScore: fwdQuality.composite || null,
          qualityTier: fwdQuality.quality || null,
          breakdown: fwdQuality.breakdown || null,
        },
        reverse: {
          ...reverse,
          // Add quality analysis data
          qualityScore: revQuality.composite || null,
          qualityTier: revQuality.quality || null,
          breakdown: revQuality.breakdown || null,
        },
        warnings: primerData?.warnings || [],
        pcr: primerData?.pcr || {
          annealingTemp: Math.min(
            forward.tm || 60,
            reverse.tm || 60
          ) - 2,
          extensionTime: Math.ceil((part.seq?.length || 1000) / 1000) * 30,
          tmDifference: Math.abs((forward.tm || 60) - (reverse.tm || 60)),
        },
        // Add pair-level quality data
        pairQuality: {
          score: pairQuality.score || null,
          quality: pairQuality.quality || null,
          heterodimer: pairQuality.heterodimer || null,
          tmDifference: pairQuality.tmDifference || null,
        },
      },
    };
  });

  // Assemble the sequence
  let assembledSeq = '';
  for (const part of designedParts) {
    assembledSeq += part.seq || '';
  }

  // Build warnings array from recommendations
  const warnings = [];

  // Add G:T mismatch warnings
  if (optimizedResult.quality.gtRisks && optimizedResult.quality.gtRisks.length > 0) {
    for (const risk of optimizedResult.quality.gtRisks) {
      warnings.push(`⚠️ G:T mismatch risk between ${risk.overhang1} and ${risk.overhang2} (positions ${risk.positions?.join(', ') || 'unknown'})`);
    }
  }

  // Add fidelity warnings
  if (optimizedResult.quality.assemblyFidelity < GG_OPTIMIZER_DEFAULTS.minAcceptableFidelity) {
    warnings.push(`⚠️ Assembly fidelity (${optimizedResult.quality.assemblyFidelityPercent}) is below recommended threshold (${GG_OPTIMIZER_DEFAULTS.minAcceptableFidelity * 100}%)`);
  }

  // Add primer quality warnings
  if (optimizedResult.recommendations) {
    for (const rec of optimizedResult.recommendations) {
      if (rec.severity === 'warning' || rec.severity === 'critical') {
        warnings.push(`${rec.severity === 'critical' ? '🔴' : '⚠️'} ${rec.message}`);
      }
    }
  }

  // Check for internal restriction sites
  const internalSiteIssues = [];
  for (const part of designedParts) {
    const siteCheck = findInternalSites(part.seq || '', enzyme);
    if (siteCheck.hasSites) {
      // Check if this part was auto-domesticated (has _domesticated flag)
      const wasDomesticated = part._domesticated === true || autoDomestication;

      if (wasDomesticated) {
        // Part was auto-domesticated - show info message instead of warning
        // The site may still appear in the fragment but will be broken at assembly junction
        const issue = {
          partId: part.id,
          sites: siteCheck.sites,
          count: siteCheck.count,
          domesticated: true,
          parentPart: part._parentPart,
        };
        internalSiteIssues.push(issue);

        // Info message - site is handled by auto-domestication
        warnings.push(`✅ Part "${part.id}" - internal ${enzyme} site will be broken at assembly junction (auto-domesticated from ${part._parentPart || 'parent'})`);
      } else {
        // Part was NOT auto-domesticated - show the regular warning
        const domestication = suggestDomestication(part.seq || '', enzyme);
        const altEnzymes = findAlternativeEnzymes(part.seq || '', enzyme);

        const issue = {
          partId: part.id,
          sites: siteCheck.sites,
          count: siteCheck.count,
          domestication,
          alternativeEnzymes: altEnzymes,
        };
        internalSiteIssues.push(issue);

        // Add warning with actionable suggestions
        let suggestion = '';
        if (domestication.success && domestication.suggestions.length > 0) {
          const mut = domestication.suggestions[0].mutations[0];
          if (mut) {
            suggestion = ` Suggested fix: ${mut.originalBase}${mut.position + 1}${mut.newBase}`;
            if (mut.isSynonymous) suggestion += ' (silent mutation)';
          }
        } else if (altEnzymes.length > 0) {
          suggestion = ` Consider using ${altEnzymes[0].enzyme} instead.`;
        }

        warnings.push(`⚠️ Part "${part.id}" contains ${siteCheck.count} internal ${enzyme} site(s) - MUST be removed before assembly.${suggestion}`);
      }
    }
  }

  // Build fidelity structure matching original format
  const fidelity = {
    individual: (optimizedResult.quality.junctionFidelities || []).map((j, i) => ({
      junction: i,
      overhang: j.overhang,
      fidelity: j.fidelity,
      fidelityPercent: `${(j.fidelity * 100).toFixed(1)}%`,
    })),
    overall: optimizedResult.quality.assemblyFidelity,
    percentage: optimizedResult.quality.assemblyFidelityPercent,
    // Additional optimized data
    gtAdjusted: true,
    baseFidelity: optimizedResult.quality.baseFidelity,
    baseFidelityPercent: optimizedResult.quality.baseFidelityPercent,
  };

  return {
    // Standard fields expected by UI
    enzyme,
    parts: designedParts,
    overhangs: finalOverhangs,
    circular,
    assembledSequence: assembledSeq,
    assembledLength: assembledSeq.length,
    fidelity,
    warnings,
    internalSiteIssues,
    hasInternalSites: internalSiteIssues.length > 0,
    protocol: generateProtocol(designedParts, enzyme),

    // Additional optimized data for advanced users
    _optimizedData: {
      overhangs: optimizedResult.overhangs,
      quality: optimizedResult.quality,
      recommendations: optimizedResult.recommendations,
      summary: optimizedResult.summary,
    },
  };
}

// ============================================================================
// EXPORTS
// ============================================================================

export default {
  // Configuration
  GG_OPTIMIZER_DEFAULTS,

  // Flanking optimization
  OPTIMAL_FLANKING_SEQUENCES,
  OPTIMAL_SPACERS,
  selectOptimalFlankingSequence,

  // G:T mismatch handling
  canPair,
  findGTMismatchRisks,
  calculateEnhancedFidelity,

  // Automatic overhang optimization
  autoOptimizeOverhangs,

  // Golden Gate primer quality scoring
  scoreGoldenGatePrimer,
  scoreGoldenGatePrimerPair,

  // Optimized design functions (recommended)
  designOptimizedGoldenGatePrimers,
  designOptimizedGoldenGateAssembly,

  // UI-compatible adapter (use this for GoldenGateDesigner component)
  designOptimizedGoldenGateAssemblyForUI,
};
