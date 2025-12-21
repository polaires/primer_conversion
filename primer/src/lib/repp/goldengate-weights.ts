/**
 * Golden Gate-Specific Scoring Weights
 *
 * This module provides GG-specific scoring modifiers that work as an
 * ADDITIVE LAYER on top of the base ASSEMBLY_WEIGHTS from weightCalibration.ts.
 *
 * Philosophy:
 * - Base scoring (ASSEMBLY_WEIGHTS) captures general primer quality
 * - GG modifiers capture assembly-specific concerns (fidelity, overhangs, etc.)
 * - Final score = blend of base and GG-specific scores
 *
 * This approach:
 * 1. Leverages the calibrated base weights (F1=88.4%, AUC=0.929)
 * 2. Adds GG-specific concerns without duplicating base scoring logic
 * 3. Allows easy tuning of GG-specific factors
 */

import { ASSEMBLY_WEIGHTS, type Weights } from '../weightCalibration.js';
import { calculateCompositeScore, classifyQuality, type CompositeScoreInput } from '../scoring.js';
import { calculateFidelity, type FidelityResult, type JunctionFidelity } from './fidelity-core.js';
import { validateOverhang, type OverhangValidation } from './overhang-validation.js';

// ============================================================================
// TYPES
// ============================================================================

export interface GGContext {
  /** Enzyme being used */
  enzyme: string;

  /** All overhangs in the assembly */
  assemblyOverhangs: string[];

  /** This primer's associated overhang */
  primerOverhang: string;

  /** Fidelity result for the assembly (cached) */
  fidelity?: FidelityResult;

  /** Whether this is forward or reverse primer */
  direction: 'forward' | 'reverse';

  /** User options */
  options?: GGScoringOptions;
}

export interface GGScoringOptions {
  /** Blend ratio: 0 = 100% base, 1 = 100% GG-specific (default: 0.30) */
  ggBlendRatio?: number;

  /** Use strict quality thresholds */
  strictQuality?: boolean;

  /** Include fidelity in scoring */
  includeFidelity?: boolean;

  /** Include overhang quality in scoring */
  includeOverhangQuality?: boolean;

  /** Include flanking sequence quality in scoring */
  includeFlankingQuality?: boolean;
}

export interface GGModifierScores {
  /** Overhang fidelity contribution (0-1) */
  overhangFidelity: number;

  /** Overhang quality score (0-1) */
  overhangQuality: number;

  /** Flanking sequence quality (0-1) */
  flankingQuality: number;

  /** Recognition site safety (0-1) */
  recognitionSiteSafety: number;

  /** Overall GG modifier score (0-100) */
  score: number;

  /** Breakdown for debugging */
  breakdown: {
    [key: string]: {
      raw: number;
      weighted: number;
      weight: number;
    };
  };
}

export interface UnifiedGGScore {
  /** Final composite score (0-100) */
  composite: number;

  /** Base score from ASSEMBLY_WEIGHTS */
  baseScore: number;

  /** GG-specific modifier score */
  ggScore: number;

  /** Blend ratio used */
  blendRatio: number;

  /** Quality classification */
  quality: {
    tier: string;
    label: string;
  };

  /** GG modifier details */
  ggModifiers: GGModifierScores;

  /** Warnings and recommendations */
  warnings: string[];
}

// ============================================================================
// CONSTANTS
// ============================================================================

/**
 * Golden Gate-specific weight modifiers
 *
 * These are applied as ADDITIONS to base scoring, not replacements.
 * They capture assembly-specific concerns not in general primer scoring.
 */
export const GG_WEIGHT_MODIFIERS = {
  // GG-specific factors (not in base scoring)
  overhangFidelity: 0.35,        // Ligation fidelity from matrix data
  overhangQuality: 0.25,         // Overhang sequence quality (palindrome, homopolymer, etc.)
  flankingQuality: 0.20,         // NEB flanking sequence optimization
  recognitionSiteSafety: 0.20,   // Avoid creating new restriction sites

  // Total should sum to 1.0 for proper normalization
};

/**
 * Overrides for base weights in GG context
 *
 * Some base weights should be adjusted for GG assembly:
 * - heterodimer: Less critical because fragments are mixed after PCR
 * - tmDiff: Important for consistent annealing in one-pot assembly
 */
export const GG_WEIGHT_OVERRIDES: Record<string, number> = {
  // In GG, heterodimer between primers is less critical
  // because fragments are mixed after PCR, not during
  heterodimer: 0.04,  // Lower from 0.10 in ASSEMBLY_WEIGHTS

  // Tm difference is more important for one-pot consistency
  tmDiff: 0.06,       // Slightly higher
};

/**
 * Default GG scoring options
 */
export const DEFAULT_GG_OPTIONS: Required<GGScoringOptions> = {
  ggBlendRatio: 0.30,
  strictQuality: false,
  includeFidelity: true,
  includeOverhangQuality: true,
  includeFlankingQuality: true,
};

/**
 * Fidelity scoring thresholds
 */
const FIDELITY_THRESHOLDS = {
  excellent: 0.99,   // 99%+ fidelity → score 1.0
  good: 0.95,        // 95-99% → score 0.9
  acceptable: 0.90,  // 90-95% → score 0.7
  marginal: 0.80,    // 80-90% → score 0.4
  poor: 0.0,         // <80% → score 0.1
};

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Convert fidelity (0-1) to a score (0-1) using thresholds
 */
function fidelityToScore(fidelity: number): number {
  if (fidelity >= FIDELITY_THRESHOLDS.excellent) return 1.0;
  if (fidelity >= FIDELITY_THRESHOLDS.good) return 0.9;
  if (fidelity >= FIDELITY_THRESHOLDS.acceptable) return 0.7;
  if (fidelity >= FIDELITY_THRESHOLDS.marginal) return 0.4;
  return 0.1;
}

/**
 * Convert overhang validation score (0-100) to normalized score (0-1)
 */
function validationToScore(validation: OverhangValidation): number {
  return validation.score / 100;
}

/**
 * Score flanking sequence quality
 * Based on NEB recommendations for optimal flanking
 */
function scoreFlankingQuality(flanking: string | undefined): number {
  if (!flanking) return 0.5; // Unknown → neutral

  const seq = flanking.toUpperCase();

  let score = 1.0;

  // Check length (optimal: 6bp)
  if (seq.length < 4) score -= 0.3;
  else if (seq.length < 6) score -= 0.1;

  // Check for homopolymers
  if (/(.)\1{2,}/.test(seq)) score -= 0.2;

  // Check GC content (optimal: 40-60%)
  const gc = (seq.match(/[GC]/g) || []).length / seq.length;
  if (gc < 0.3 || gc > 0.7) score -= 0.1;

  // Check for palindrome
  const complement: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
  const rc = seq.split('').reverse().map(b => complement[b] || b).join('');
  if (seq === rc) score -= 0.15;

  return Math.max(0, score);
}

/**
 * Score recognition site safety
 * Checks if primer design might create new restriction sites
 */
function scoreRecognitionSiteSafety(
  primerSequence: string,
  enzyme: string
): number {
  // Recognition sites for common GG enzymes
  const recognitionSites: Record<string, string[]> = {
    'BsaI': ['GGTCTC', 'GAGACC'],
    'BsaI-HFv2': ['GGTCTC', 'GAGACC'],
    'BsmBI': ['CGTCTC', 'GAGACG'],
    'BsmBI-v2': ['CGTCTC', 'GAGACG'],
    'BbsI': ['GAAGAC', 'GTCTTC'],
    'BbsI-HF': ['GAAGAC', 'GTCTTC'],
    'Esp3I': ['CGTCTC', 'GAGACG'],
    'SapI': ['GCTCTTC', 'GAAGAGC'],
  };

  const sites = recognitionSites[enzyme] || recognitionSites['BsaI'];
  const seq = primerSequence.toUpperCase();

  // Check for recognition sites in primer
  for (const site of sites) {
    if (seq.includes(site)) {
      return 0.2; // Major penalty
    }
  }

  // Check for partial sites (might create site with template)
  for (const site of sites) {
    // Check if primer ends with beginning of site
    for (let i = 2; i < site.length; i++) {
      if (seq.endsWith(site.slice(0, i))) {
        return 0.7; // Minor penalty
      }
    }
  }

  return 1.0; // Safe
}

// ============================================================================
// MAIN SCORING FUNCTIONS
// ============================================================================

/**
 * Compute GG-specific modifier scores
 *
 * These scores capture Golden Gate-specific concerns:
 * - Overhang fidelity (from ligation matrix)
 * - Overhang quality (palindromes, homopolymers, etc.)
 * - Flanking sequence quality
 * - Recognition site safety
 */
export function computeGGModifiers(
  primerSequence: string,
  context: GGContext
): GGModifierScores {
  const opts = { ...DEFAULT_GG_OPTIONS, ...context.options };
  const breakdown: GGModifierScores['breakdown'] = {};

  // 1. Overhang fidelity score
  let overhangFidelityScore = 1.0;
  if (opts.includeFidelity && context.assemblyOverhangs.length > 0) {
    // Get or calculate fidelity
    const fidelity = context.fidelity ||
      calculateFidelity(context.assemblyOverhangs, context.enzyme);

    // Find this overhang's individual fidelity
    const junctionFidelity = fidelity.junctions.find(
      j => j.overhang === context.primerOverhang.toUpperCase()
    );

    if (junctionFidelity) {
      overhangFidelityScore = fidelityToScore(junctionFidelity.fidelity);
    } else {
      // Use assembly fidelity as fallback
      overhangFidelityScore = fidelityToScore(fidelity.assemblyFidelity);
    }
  }

  breakdown.overhangFidelity = {
    raw: overhangFidelityScore,
    weight: GG_WEIGHT_MODIFIERS.overhangFidelity,
    weighted: overhangFidelityScore * GG_WEIGHT_MODIFIERS.overhangFidelity,
  };

  // 2. Overhang quality score
  let overhangQualityScore = 1.0;
  if (opts.includeOverhangQuality && context.primerOverhang) {
    const validation = validateOverhang(context.primerOverhang, {
      existingOverhangs: context.assemblyOverhangs.filter(
        o => o.toUpperCase() !== context.primerOverhang.toUpperCase()
      ),
      enzyme: context.enzyme,
    });
    overhangQualityScore = validationToScore(validation);
  }

  breakdown.overhangQuality = {
    raw: overhangQualityScore,
    weight: GG_WEIGHT_MODIFIERS.overhangQuality,
    weighted: overhangQualityScore * GG_WEIGHT_MODIFIERS.overhangQuality,
  };

  // 3. Flanking sequence quality
  // Extract flanking from primer (first 6bp before recognition site)
  let flankingQualityScore = 1.0;
  if (opts.includeFlankingQuality) {
    // Assume flanking is first 6bp of primer
    const flanking = primerSequence.slice(0, 6);
    flankingQualityScore = scoreFlankingQuality(flanking);
  }

  breakdown.flankingQuality = {
    raw: flankingQualityScore,
    weight: GG_WEIGHT_MODIFIERS.flankingQuality,
    weighted: flankingQualityScore * GG_WEIGHT_MODIFIERS.flankingQuality,
  };

  // 4. Recognition site safety
  const recognitionSiteSafetyScore = scoreRecognitionSiteSafety(
    primerSequence,
    context.enzyme
  );

  breakdown.recognitionSiteSafety = {
    raw: recognitionSiteSafetyScore,
    weight: GG_WEIGHT_MODIFIERS.recognitionSiteSafety,
    weighted: recognitionSiteSafetyScore * GG_WEIGHT_MODIFIERS.recognitionSiteSafety,
  };

  // Calculate overall GG score (0-100)
  const totalWeight = Object.values(GG_WEIGHT_MODIFIERS).reduce((a, b) => a + b, 0);
  const weightedSum = Object.values(breakdown).reduce((sum, b) => sum + b.weighted, 0);
  const normalizedScore = (weightedSum / totalWeight) * 100;

  return {
    overhangFidelity: overhangFidelityScore,
    overhangQuality: overhangQualityScore,
    flankingQuality: flankingQualityScore,
    recognitionSiteSafety: recognitionSiteSafetyScore,
    score: normalizedScore,
    breakdown,
  };
}

/**
 * Compute unified Golden Gate primer score
 *
 * Combines base scoring (ASSEMBLY_WEIGHTS) with GG-specific modifiers
 * using a configurable blend ratio.
 *
 * @param primerScores - Individual feature scores from primer analysis
 * @param primerSequence - Full primer sequence
 * @param context - GG context (enzyme, overhangs, etc.)
 * @returns Unified score with breakdown
 */
export function computeUnifiedGGScore(
  primerScores: CompositeScoreInput,
  primerSequence: string,
  context: GGContext
): UnifiedGGScore {
  const opts = { ...DEFAULT_GG_OPTIONS, ...context.options };
  const warnings: string[] = [];

  // Step 1: Calculate base score using ASSEMBLY_WEIGHTS with GG overrides
  const effectiveWeights = { ...ASSEMBLY_WEIGHTS, ...GG_WEIGHT_OVERRIDES };
  const baseResult = calculateCompositeScore(primerScores, effectiveWeights);
  const baseScore = baseResult.score;

  // Step 2: Calculate GG-specific modifiers
  const ggModifiers = computeGGModifiers(primerSequence, context);
  const ggScore = ggModifiers.score;

  // Step 3: Blend scores
  const blendRatio = opts.ggBlendRatio;
  const composite = baseScore * (1 - blendRatio) + ggScore * blendRatio;

  // Step 4: Classify quality
  const qualityResult = classifyQuality(composite);

  // Step 5: Generate warnings
  if (ggModifiers.overhangFidelity < 0.7) {
    warnings.push(`Low overhang fidelity: ${(ggModifiers.overhangFidelity * 100).toFixed(0)}%`);
  }
  if (ggModifiers.overhangQuality < 0.7) {
    warnings.push('Overhang has quality issues');
  }
  if (ggModifiers.recognitionSiteSafety < 0.5) {
    warnings.push('Primer may create/contain recognition site');
  }

  // Strict quality check
  if (opts.strictQuality && qualityResult.tier === 'poor') {
    warnings.push('Primer quality is poor - consider alternatives');
  }

  return {
    composite,
    baseScore,
    ggScore,
    blendRatio,
    quality: {
      tier: qualityResult.tier,
      label: qualityResult.label,
    },
    ggModifiers,
    warnings,
  };
}

/**
 * Get effective weights for GG primer scoring
 *
 * Merges ASSEMBLY_WEIGHTS with GG-specific overrides.
 */
export function getEffectiveGGWeights(): Weights {
  return { ...ASSEMBLY_WEIGHTS, ...GG_WEIGHT_OVERRIDES };
}

/**
 * Compare two primer candidates using unified GG scoring
 */
export function compareGGPrimers(
  primer1: { scores: CompositeScoreInput; sequence: string },
  primer2: { scores: CompositeScoreInput; sequence: string },
  context: GGContext
): {
  winner: 1 | 2;
  score1: UnifiedGGScore;
  score2: UnifiedGGScore;
  marginOfVictory: number;
} {
  const score1 = computeUnifiedGGScore(primer1.scores, primer1.sequence, context);
  const score2 = computeUnifiedGGScore(primer2.scores, primer2.sequence, context);

  return {
    winner: score1.composite >= score2.composite ? 1 : 2,
    score1,
    score2,
    marginOfVictory: Math.abs(score1.composite - score2.composite),
  };
}

// ============================================================================
// EXPORTS
// ============================================================================

export { ASSEMBLY_WEIGHTS } from '../weightCalibration.js';
