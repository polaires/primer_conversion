/**
 * Piecewise Logistic Scoring Functions
 *
 * This module implements biologically-meaningful scoring functions that replace
 * linear penalties with piecewise logistic curves. These functions are based on
 * research from PrimerScore2 and empirical validation studies.
 *
 * Key principles:
 * 1. "Free zone" near optimal values (flat region with score ~1.0)
 * 2. Sharp penalty at biological thresholds
 * 3. Bounded output (interpretable 0-1 range)
 * 4. Matches observed PCR behavior
 *
 * Reference: Section 4.3 of SCORING_MANUSCRIPT_STRATEGY.md
 */

import { DEFAULT_WEIGHTS } from './weightCalibration.js';
import { SCORING_PRESETS } from './presets.js';
import { calculate3primeTerminalDG } from './tmQ5.js';

// Re-export for backward compatibility
export { SCORING_PRESETS };

// ============================================================================
// Type Definitions
// ============================================================================

export interface PiecewiseLogisticOptions {
  optimalLow?: number;
  optimalHigh?: number;
  acceptableLow?: number;
  acceptableHigh?: number;
  steepness?: number;
  floorScore?: number;
}

export interface TmScoringOptions {
  optimalLow?: number;
  optimalHigh?: number;
  acceptableLow?: number;
  acceptableHigh?: number;
  steepness?: number;
}

export interface GcScoringOptions {
  optimalLow?: number;
  optimalHigh?: number;
  acceptableLow?: number;
  acceptableHigh?: number;
  steepness?: number;
}

export interface Terminal3DGOptions {
  optimalLow?: number;
  optimalHigh?: number;
  tooLooseDecay?: number;
  tooTightDecay?: number;
}

export interface TmDiffOptions {
  freeZone?: number;
  mildZone?: number;
  moderateZone?: number;
  steepDecay?: number;
}

export interface HairpinOptions {
  threshold?: number;
  steepness?: number;
}

export interface HomodimerOptions {
  threshold?: number;
  steepness?: number;
}

export interface HeterodimerOptions {
  threshold?: number;
  steepness?: number;
}

export interface OffTargetOptions {
  singlePenalty?: number;
  multiplier?: number;
  critical?: number;
}

export interface OffTargetClassification {
  counts: {
    highRisk: number;
    mediumRisk: number;
    lowRisk: number;
  };
}

export interface LengthOptions {
  optimalLow?: number;
  optimalHigh?: number;
  acceptableLow?: number;
  acceptableHigh?: number;
  steepness?: number;
}

export interface GcClampOptions {
  idealCount?: number;
}

export interface ThreePrimeCompositionOptions {
  idealGcClamp?: number;
  optimalDGLow?: number;
  optimalDGHigh?: number;
  gcClampWeight?: number;
  terminalDGWeight?: number;
  patternWeight?: number;
}

export interface HomopolymerOptions {
  maxRun?: number;
  penaltyPerBase?: number;
}

export interface GQuadruplexOptions {
  g4MotifScore?: number;
  ggggScore?: number;
  multiGggScore?: number;
}

export interface GQuadruplexAnalysis {
  score: number;
  hasG4Motif: boolean;
  hasGGGG: boolean;
  gggCount: number;
  gggRuns: string[];
  severity: 'ok' | 'caution' | 'warning' | 'critical';
  message: string;
}

export interface AmpliconLengthOptions {
  optimalLow?: number;
  optimalHigh?: number;
  acceptableLow?: number;
  acceptableHigh?: number;
  steepness?: number;
}

export interface DistanceToROIOptions {
  optimalLow?: number;
  optimalHigh?: number;
  acceptableLow?: number;
  acceptableHigh?: number;
  steepness?: number;
}

export interface AmpliconStructureOptions {
  gcThreshold?: number;
  windowSize?: number;
  maxGcStretches?: number;
  homopolymerThreshold?: number;
}

export interface PrimerScores {
  tm: number;
  gc: number;
  length: number;
  terminal3DG: number;
  gcClamp: number;
  homopolymer: number;
  hairpin: number;
  homodimer: number;
  offTarget: number;
  gQuadruplex: number;
  threePrimeComp: number;
}

export interface PairScores {
  tmDiff: number;
  heterodimer: number;
  ampliconLength?: number;
  ampliconStructure?: number;
}

export interface CompositeScoreInput {
  tmFwd?: number;
  tmRev?: number;
  gcFwd?: number;
  gcRev?: number;
  lengthFwd?: number;
  lengthRev?: number;
  hairpinFwd?: number;
  hairpinRev?: number;
  selfDimerFwd?: number;
  selfDimerRev?: number;
  heterodimer?: number;
  tmDiff?: number;
  terminal3DG?: number;
  gcClampFwd?: number;
  gcClampRev?: number;
  threePrimeCompFwd?: number;
  threePrimeCompRev?: number;
  offTarget?: number;
  gQuadruplexFwd?: number;
  gQuadruplexRev?: number;
  homopolymerFwd?: number;
  homopolymerRev?: number;
  ampliconLength?: number;
  ampliconStructure?: number;
}

export interface ScoreBreakdownItem {
  score: number;
  weight: number;
  contribution: number;
}

export interface CompositeScoreResult {
  score: number;
  rawScore: number;
  breakdown: Record<string, ScoreBreakdownItem>;
  totalWeight: number;
}

export interface QualityTier {
  tier: 'excellent' | 'good' | 'acceptable' | 'marginal' | 'poor';
  label: string;
  description: string;
  color: string;
}

export interface Primer {
  seq: string;
  tm: number;
  gc: number;
  dg?: number;
  offTargetCount?: number;
}

export interface PrimerAnalysisResult {
  sequence: string;
  length: number;
  tm: number;
  gc: number;
  dg?: number;
  scores: PrimerScores;
  gQuadruplex: GQuadruplexAnalysis;
}

export interface PairAnalysisResult {
  tmDiff: number;
  heterodimerDG: number | null;
  ampliconLength: number | null;
  scores: PairScores;
}

export interface AnalyzePrimerPairOptions {
  heterodimerDG?: number | null;
  hairpinDGFwd?: number | null;
  hairpinDGRev?: number | null;
  homodimerDGFwd?: number | null;
  homodimerDGRev?: number | null;
  ampliconLength?: number | null;
  includeRecommendations?: boolean;
  mode?: string;
  weights?: Record<string, number> | null;
  externalWarnings?: string[] | null;
}

export interface PrimerPairAnalysis {
  forward: PrimerAnalysisResult;
  reverse: PrimerAnalysisResult | null;
  pair: PairAnalysisResult | null;
  composite: CompositeScoreResult | null;
  quality: QualityTier | null;
  warnings: string[];
  recommendations: string[];
  mode?: string;
  presetName?: string;
  criticalWarnings?: number;
  effectiveScore?: number;
}

export interface BatchDatasetEntry {
  scores: CompositeScoreInput;
  success?: boolean;
  actual?: boolean;
}

export interface FeatureDiscrimination {
  successMean: number;
  failureMean: number;
  diff: number;
}

export interface GQuadruplexStats {
  count: number;
  successCount: number;
  successRate: number | null;
}

export interface BatchAnalysisStats {
  total: number;
  success: number;
  failure: number;
  successRate: number;
  scoreDistribution: {
    success: {
      scores: number[];
      mean: number;
      std: number;
    };
    failure: {
      scores: number[];
      mean: number;
      std: number;
    };
  };
  featureDiscrimination: Record<string, FeatureDiscrimination>;
  topDiscriminativeFeatures: Array<{ feature: string } & FeatureDiscrimination>;
  gQuadruplex: {
    g4Motif: GQuadruplexStats;
    ggggRun: GQuadruplexStats;
    multiGgg: GQuadruplexStats;
    noRisk: GQuadruplexStats;
  };
}

// ============================================================================
// Core Scoring Functions
// ============================================================================

/**
 * Generic piecewise logistic function
 *
 * Creates a scoring curve with:
 * - Full score (1.0) in optimal range
 * - Linear decay in acceptable range
 * - Logistic decay outside acceptable range
 *
 * @param value - The value to score
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function piecewiseLogistic(
  value: number,
  {
    optimalLow,
    optimalHigh,
    acceptableLow,
    acceptableHigh,
    steepness = 0.5,
    floorScore = 0.0,
  }: PiecewiseLogisticOptions = {}
): number {
  // Require all range parameters
  if (
    optimalLow === undefined ||
    optimalHigh === undefined ||
    acceptableLow === undefined ||
    acceptableHigh === undefined
  ) {
    throw new Error('piecewiseLogistic requires all range parameters');
  }

  // In optimal range - perfect score
  if (value >= optimalLow && value <= optimalHigh) {
    return 1.0;
  }

  // Below optimal
  if (value < optimalLow) {
    if (value >= acceptableLow) {
      // Linear decay in acceptable range (1.0 -> 0.7)
      const ratio = (value - acceptableLow) / (optimalLow - acceptableLow);
      return 0.7 + 0.3 * ratio;
    }
    // Logistic decay below acceptable
    const excess = acceptableLow - value;
    const score = 0.7 / (1 + Math.exp(steepness * excess));
    return Math.max(floorScore, score);
  }

  // Above optimal
  if (value > optimalHigh) {
    if (value <= acceptableHigh) {
      // Linear decay in acceptable range (1.0 -> 0.7)
      const ratio = (acceptableHigh - value) / (acceptableHigh - optimalHigh);
      return 0.7 + 0.3 * ratio;
    }
    // Logistic decay above acceptable
    const excess = value - acceptableHigh;
    const score = 0.7 / (1 + Math.exp(steepness * excess));
    return Math.max(floorScore, score);
  }

  return 1.0;
}

/**
 * Tm (Melting Temperature) scoring for Sanger sequencing
 *
 * Optimal: 55-60°C (lower than RT-qPCR due to different read requirements)
 * Acceptable: 50-65°C
 * Problematic: <50°C (won't anneal) or >65°C (secondary structure risk)
 *
 * @param tm - Melting temperature in °C
 * @param options - Scoring parameters (can override defaults)
 * @returns Score from 0 to 1
 */
export function scoreTm(tm: number, options: TmScoringOptions = {}): number {
  // Guard against invalid input
  if (tm === null || tm === undefined || Number.isNaN(tm)) {
    return 0.5;  // Return neutral score for invalid input
  }

  const {
    optimalLow = 55,
    optimalHigh = 60,
    acceptableLow = 50,
    acceptableHigh = 65,
    steepness = 0.5,
  } = options;

  return piecewiseLogistic(tm, {
    optimalLow,
    optimalHigh,
    acceptableLow,
    acceptableHigh,
    steepness,
  });
}

/**
 * GC Content scoring
 *
 * Optimal: 40-60% (classic guideline)
 * Acceptable: 30-70%
 * Problematic: <30% (low stability) or >70% (secondary structure, mispriming)
 *
 * @param gc - GC content as fraction (0 to 1) or percentage (0 to 100)
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreGc(gc: number, options: GcScoringOptions = {}): number {
  // Guard against invalid input
  if (gc === null || gc === undefined || Number.isNaN(gc)) {
    return 0.5;  // Return neutral score for invalid input
  }

  // Convert to percentage if given as fraction
  const gcPct = gc <= 1 ? gc * 100 : gc;

  const {
    optimalLow = 40,
    optimalHigh = 60,
    acceptableLow = 30,
    acceptableHigh = 70,
    steepness = 0.15,  // Gentler slope for GC (wider acceptable range)
  } = options;

  return piecewiseLogistic(gcPct, {
    optimalLow,
    optimalHigh,
    acceptableLow,
    acceptableHigh,
    steepness,
  });
}

/**
 * 3' Terminal ΔG scoring
 *
 * Based on experimental finding:
 * "priming was detectable when 3'-terminal portion formed duplex
 *  more stable than -11 kcal/mol"
 *
 * Optimal range: -6 to -11 kcal/mol
 * Too loose (> -6): Won't initiate efficiently
 * Too tight (< -12): Mispriming risk (binds non-specifically)
 *
 * Note: More negative ΔG = stronger binding
 *
 * @param dG - 3' terminal ΔG in kcal/mol (negative values)
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreTerminal3DG(dG: number, options: Terminal3DGOptions = {}): number {
  const {
    optimalLow = -11,   // Most stable acceptable (more negative)
    optimalHigh = -6,   // Least stable acceptable (less negative)
    tooLooseDecay = 0.3,   // Exponential decay rate for too-loose binding
    tooTightDecay = 0.15,  // Milder decay for too-tight (still works, less specific)
  } = options;

  // In optimal range - perfect score
  if (dG >= optimalLow && dG <= optimalHigh) {
    return 1.0;
  }

  // Too loose (dG > -6, towards 0)
  // Won't initiate efficiently - exponential penalty
  if (dG > optimalHigh) {
    const excess = dG - optimalHigh;
    return Math.exp(-tooLooseDecay * excess);
  }

  // Too tight (dG < -11, very negative)
  // Still works but less specific - milder penalty
  if (dG < optimalLow) {
    const excess = optimalLow - dG;
    return Math.exp(-tooTightDecay * excess);
  }

  return 1.0;
}

/**
 * Tm Difference scoring
 *
 * Key finding from PrimerScore2: "little effect" on amplification
 * We use a generous free zone and mild penalties
 *
 * @param tmFwd - Forward primer Tm in °C
 * @param tmRev - Reverse primer Tm in °C
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreTmDiff(tmFwd: number, tmRev: number, options: TmDiffOptions = {}): number {
  const {
    freeZone = 3,      // 0-3°C: No penalty (expanded from typical 1-2°C)
    mildZone = 5,      // 3-5°C: Mild penalty
    moderateZone = 8,  // 5-8°C: Moderate penalty
    steepDecay = 0.2,  // Decay rate for >8°C
  } = options;

  const diff = Math.abs(tmFwd - tmRev);

  // 0-3°C: No penalty
  if (diff <= freeZone) {
    return 1.0;
  }

  // 3-5°C: Mild penalty (1.0 -> 0.8)
  if (diff <= mildZone) {
    return 0.9 - 0.1 * (diff - freeZone) / (mildZone - freeZone);
  }

  // 5-8°C: Moderate penalty (0.8 -> 0.5)
  if (diff <= moderateZone) {
    return 0.7 - 0.2 * (diff - mildZone) / (moderateZone - mildZone);
  }

  // >8°C: Steep exponential penalty
  return 0.5 * Math.exp(-steepDecay * (diff - moderateZone));
}

/**
 * Hairpin ΔG scoring
 *
 * Hairpin formation competes with target binding.
 * Threshold: -3 kcal/mol (more negative = more stable hairpin = worse)
 *
 * @param dG - Hairpin ΔG in kcal/mol
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreHairpin(dG: number, options: HairpinOptions = {}): number {
  const {
    threshold = -3,    // ΔG threshold for significant hairpin
    steepness = 0.8,   // Decay rate
  } = options;

  // Positive or weakly negative: no stable hairpin
  if (dG >= threshold) {
    return 1.0;
  }

  // Stable hairpin: exponential penalty
  const excess = threshold - dG;
  return Math.exp(-steepness * excess);
}

/**
 * Homodimer ΔG scoring
 *
 * Self-dimer formation competes with target binding.
 * Threshold: -6 kcal/mol (more negative = more stable dimer = worse)
 *
 * @param dG - Homodimer ΔG in kcal/mol
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreHomodimer(dG: number, options: HomodimerOptions = {}): number {
  const {
    threshold = -6,    // ΔG threshold for significant homodimer
    steepness = 0.5,   // Decay rate
  } = options;

  // Positive or weakly negative: no stable dimer
  if (dG >= threshold) {
    return 1.0;
  }

  // Stable homodimer: exponential penalty
  const excess = threshold - dG;
  return Math.exp(-steepness * excess);
}

/**
 * Heterodimer ΔG scoring
 *
 * Cross-dimer between fwd and rev primers competes with target binding.
 * Threshold: -6 kcal/mol
 *
 * @param dG - Heterodimer ΔG in kcal/mol
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreHeterodimer(dG: number, options: HeterodimerOptions = {}): number {
  const {
    threshold = -6,
    steepness = 0.5,
  } = options;

  if (dG >= threshold) {
    return 1.0;
  }

  const excess = threshold - dG;
  return Math.exp(-steepness * excess);
}

/**
 * Off-target count scoring (simple count-based)
 *
 * Based on GM1 finding that off-target is the "dominant failure factor"
 * Exponential penalty scaling for multiple off-targets
 *
 * @param count - Number of off-target binding sites
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreOffTarget(count: number, options: OffTargetOptions = {}): number {
  const {
    singlePenalty = 0.3,   // Penalty for first off-target (score = 0.7)
    multiplier = 3,        // Each additional off-target multiplies penalty
    critical = 3,          // 3+ off-targets = disqualify
  } = options;

  if (count === 0) {
    return 1.0;
  }

  if (count >= critical) {
    return 0.0;  // Disqualify
  }

  // Exponential penalty: 0.7, 0.3, 0.1 for 1, 2, 3 sites
  const penalty = singlePenalty * Math.pow(multiplier, count - 1);
  return Math.max(0, 1 - penalty);
}

/**
 * Off-target classification scoring (Type A-F based)
 *
 * Uses comprehensive off-target classification by risk level.
 * See offTargetClassification.js for detailed type definitions.
 *
 * Type A (exact/near-exact match): HIGH RISK - alternative amplification
 * Type B (antisense binding): HIGH RISK - primer-template conflict
 * Type C (partial 3' homology): MEDIUM RISK - mispriming potential
 * Type D (internal homology): LOW RISK - usually tolerable
 *
 * @param classification - Result from classifyOffTargets()
 * @returns Score from 0 to 1
 */
export function scoreOffTargetClassification(classification: OffTargetClassification): number {
  const { counts } = classification;

  // High risk sites are critical - each one severely impacts score
  if (counts.highRisk >= 3) return 0.0;  // Disqualify
  if (counts.highRisk === 2) return 0.1;
  if (counts.highRisk === 1) return 0.3;

  // Medium risk sites - moderate impact
  let score = 1.0;
  score -= counts.mediumRisk * 0.15;

  // Low risk sites - minor impact
  score -= counts.lowRisk * 0.02;

  return Math.max(0, Math.min(1, score));
}

/**
 * Primer length scoring
 *
 * Optimal: 18-24bp (standard recommendation)
 * Acceptable: 15-30bp
 *
 * @param length - Primer length in bp
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreLength(length: number, options: LengthOptions = {}): number {
  const {
    optimalLow = 18,
    optimalHigh = 24,
    acceptableLow = 15,
    acceptableHigh = 30,
    steepness = 0.3,
  } = options;

  return piecewiseLogistic(length, {
    optimalLow,
    optimalHigh,
    acceptableLow,
    acceptableHigh,
    steepness,
  });
}

/**
 * GC Clamp scoring
 *
 * Checks for G or C in the last 1-2 bases of 3' end
 * Ideal: 1 G/C in last 2 bases
 * Good: 2 G/C in last 2 bases (slightly less ideal due to mispriming risk)
 * Poor: 0 G/C in last 2 bases (weak 3' end)
 *
 * @param seq - Primer sequence
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreGcClamp(seq: string, options: GcClampOptions = {}): number {
  const {
    idealCount = 1,    // Ideal is 1 G/C in last 2 bases
  } = options;

  const last2 = seq.slice(-2).toUpperCase();
  const gcCount = (last2.match(/[GC]/g) || []).length;

  if (gcCount === idealCount) {
    return 1.0;
  }
  if (gcCount === 2) {
    return 0.85;  // Strong clamp, slight mispriming risk
  }
  if (gcCount === 0) {
    return 0.5;   // Weak 3' end
  }

  return 0.9;  // Fallback
}

/**
 * 3' End Composition scoring
 *
 * Comprehensive scoring of the 3' end region combining:
 * - GC clamp status (G/C in last 2 bases)
 * - Terminal ΔG for primer binding stability
 * - Pattern analysis (poly-A/T runs, AT-rich regions)
 *
 * This is the key feature for iterative design optimization.
 *
 * @param seq - Primer sequence
 * @param terminalDG - 3' terminal ΔG (optional, will calculate if not provided)
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function score3PrimeComposition(
  seq: string,
  terminalDG: number | null = null,
  options: ThreePrimeCompositionOptions = {}
): number {
  const {
    // GC clamp ideal is 1 G/C in last 2 bases
    idealGcClamp = 1,
    // Terminal ΔG optimal range (kcal/mol)
    optimalDGLow = -11,
    optimalDGHigh = -6,
    // Weight distribution
    gcClampWeight = 0.40,
    terminalDGWeight = 0.35,
    patternWeight = 0.25,
  } = options;

  seq = seq.toUpperCase();
  let score = 1.0;

  // Analyze 3' end
  const last5 = seq.slice(-5);
  const last2 = seq.slice(-2);
  const gcLast2 = (last2.match(/[GC]/g) || []).length;
  const gcLast5 = (last5.match(/[GC]/g) || []).length;
  const endsWithGC = /[GC]$/.test(seq);
  const hasPolyAT = /[AT]{4,}/.test(last5);
  const hasPolyA = /AAA/.test(last5);
  const hasPolyT = /TTT/.test(last5);

  // GC Clamp component (40% weight)
  let gcClampScore = 1.0;
  if (gcLast2 === idealGcClamp) {
    gcClampScore = 1.0;  // Ideal
  } else if (gcLast2 === 2) {
    gcClampScore = 0.85;  // Strong clamp, slight mispriming risk
  } else if (gcLast2 === 0) {
    gcClampScore = 0.50;  // No clamp - significant penalty
  }

  // Terminal ΔG component (35% weight)
  let dgScore = 1.0;
  if (terminalDG !== null) {
    if (terminalDG >= optimalDGLow && terminalDG <= optimalDGHigh) {
      dgScore = 1.0;  // Optimal range
    } else if (terminalDG > optimalDGHigh) {
      // Too weak (towards 0)
      const excess = terminalDG - optimalDGHigh;
      dgScore = Math.max(0.2, 1.0 - excess * 0.12);
    } else {
      // Too tight (very negative)
      const excess = optimalDGLow - terminalDG;
      dgScore = Math.max(0.5, 1.0 - excess * 0.05);  // Less severe
    }
  }

  // Pattern component (25% weight)
  let patternScore = 1.0;
  if (hasPolyAT) {
    patternScore -= 0.40;  // Poly-A/T is problematic
  }
  if (!endsWithGC) {
    patternScore -= 0.15;  // Prefer ending with G/C
  }
  if (gcLast5 <= 1) {
    patternScore -= 0.15;  // AT-rich 3' end
  }
  if (hasPolyA || hasPolyT) {
    patternScore -= 0.10;  // Triple A or T
  }
  patternScore = Math.max(0, patternScore);

  // Weighted combination
  score = gcClampScore * gcClampWeight +
          dgScore * terminalDGWeight +
          patternScore * patternWeight;

  // Normalize to ensure max is 1.0
  score = score / (gcClampWeight + terminalDGWeight + patternWeight);

  return Math.max(0, Math.min(1, Math.round(score * 1000) / 1000));
}

/**
 * Homopolymer run scoring
 *
 * Penalizes runs of 4+ identical bases (polymerase slippage risk)
 *
 * @param seq - Primer sequence
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreHomopolymer(seq: string, options: HomopolymerOptions = {}): number {
  const {
    maxRun = 3,        // Runs of 4+ are penalized
    penaltyPerBase = 0.15,  // Penalty per extra base in run
  } = options;

  // Find longest homopolymer run
  let maxFound = 1;
  let currentRun = 1;

  for (let i = 1; i < seq.length; i++) {
    if (seq[i] === seq[i - 1]) {
      currentRun++;
      maxFound = Math.max(maxFound, currentRun);
    } else {
      currentRun = 1;
    }
  }

  if (maxFound <= maxRun) {
    return 1.0;
  }

  // Penalty for each base over the limit
  const excess = maxFound - maxRun;
  return Math.max(0.3, 1 - excess * penaltyPerBase);
}

/**
 * G-Quadruplex risk scoring
 *
 * G-Quadruplexes are stable secondary structures that form in G-rich sequences.
 * They are stabilized by K+ and Mg2+ (present in PCR buffers like NEB Q5).
 * They cause polymerase arrest - Q5 especially struggles with these.
 *
 * Detection criteria:
 * - Canonical G4 motif: G3+N1-7G3+N1-7G3+N1-7G3+ (critical - disqualifies)
 * - GGGG runs (severe - high failure risk)
 * - Multiple GGG runs (moderate - potential inter-strand G4)
 *
 * @param seq - Primer sequence
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreGQuadruplex(seq: string, options: GQuadruplexOptions = {}): number {
  const {
    g4MotifScore = 0.0,      // Canonical G4 motif = complete failure
    ggggScore = 0.2,         // GGGG run = severe penalty
    multiGggScore = 0.6,     // Multiple GGG = moderate penalty
  } = options;

  const sequence = seq.toUpperCase();

  // Canonical intramolecular G-quadruplex motif
  // G{3,} followed by 1-7 any bases, repeated 4 times
  const g4Regex = /G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}/;
  const hasG4Motif = g4Regex.test(sequence);

  if (hasG4Motif) {
    return g4MotifScore;  // Critical - primer will fail
  }

  // Check for GGGG runs (very problematic for proofreading polymerases)
  const hasGGGG = /GGGG/.test(sequence);
  if (hasGGGG) {
    return ggggScore;  // Severe - may cause polymerase pausing
  }

  // Count GGG runs (for inter-strand G4 potential)
  const gggMatches = sequence.match(/GGG+/g) || [];
  if (gggMatches.length >= 2) {
    return multiGggScore;  // Moderate - potential inter-strand G4
  }

  return 1.0;  // No G-quadruplex risk
}

/**
 * Analyze G-Quadruplex risk with detailed information
 *
 * Returns both the score and diagnostic information for UI display.
 *
 * @param seq - Primer sequence
 * @returns G-Quadruplex analysis results
 */
export function analyzeGQuadruplex(seq: string): GQuadruplexAnalysis {
  const sequence = seq.toUpperCase();

  // Canonical G4 motif
  const g4Regex = /G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}/;
  const hasG4Motif = g4Regex.test(sequence);

  // GGGG runs
  const hasGGGG = /GGGG/.test(sequence);

  // Count GGG runs
  const gggMatches = sequence.match(/GGG+/g) || [];
  const gggCount = gggMatches.length;

  // Determine severity and message
  let severity: 'ok' | 'caution' | 'warning' | 'critical';
  let message: string;
  let score: number;

  if (hasG4Motif) {
    severity = 'critical';
    message = 'G-Quadruplex motif detected - primer will likely fail';
    score = 0.0;
  } else if (hasGGGG) {
    severity = 'warning';
    message = 'GGGG run detected - may cause polymerase pausing';
    score = 0.2;
  } else if (gggCount >= 2) {
    severity = 'caution';
    message = 'Multiple GGG runs - potential inter-strand G4 formation';
    score = 0.6;
  } else {
    severity = 'ok';
    message = 'No G-Quadruplex risk detected';
    score = 1.0;
  }

  return {
    score,
    hasG4Motif,
    hasGGGG,
    gggCount,
    gggRuns: gggMatches,
    severity,
    message,
  };
}

/**
 * Amplicon length scoring (Sanger-specific)
 *
 * For Sanger sequencing, optimal amplicon length depends on:
 * - Read length capability (~800-1000bp for modern sequencers)
 * - Distance to region of interest
 *
 * Optimal: 400-800bp (fits within single read)
 * Acceptable: 200-1200bp
 * Problematic: <200bp (limited coverage) or >1200bp (may need two reads)
 *
 * @param length - Amplicon length in bp
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreAmpliconLength(length: number, options: AmpliconLengthOptions = {}): number {
  const {
    optimalLow = 400,
    optimalHigh = 800,
    acceptableLow = 200,
    acceptableHigh = 1200,
    steepness = 0.005,  // Gentle slope for amplicon length
  } = options;

  return piecewiseLogistic(length, {
    optimalLow,
    optimalHigh,
    acceptableLow,
    acceptableHigh,
    steepness,
  });
}

/**
 * Distance to Region of Interest (ROI) scoring (Sanger-specific)
 *
 * For Sanger sequencing, the sequencing read starts from the primer and
 * extends ~800-1000bp. The ROI should be within the high-quality read zone.
 *
 * Quality zones for Sanger reads:
 * - 0-50bp: Low quality (near primer, noisy signal)
 * - 50-700bp: High quality read zone (optimal ROI location)
 * - 700-1000bp: Declining quality
 * - >1000bp: Unreliable
 *
 * @param distance - Distance from primer binding site to ROI start (bp)
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreDistanceToROI(distance: number, options: DistanceToROIOptions = {}): number {
  const {
    optimalLow = 100,     // Start of high-quality zone
    optimalHigh = 500,    // Sweet spot for clear reads
    acceptableLow = 50,   // Too close = noisy
    acceptableHigh = 700, // Near end of reliable reads
    steepness = 0.01,
  } = options;

  // Negative distance means ROI is before primer (sequencing wrong direction)
  if (distance < 0) {
    return 0.0;  // Cannot sequence backwards
  }

  return piecewiseLogistic(distance, {
    optimalLow,
    optimalHigh,
    acceptableLow,
    acceptableHigh,
    steepness,
  });
}

/**
 * Amplicon secondary structure scoring (Sanger-specific)
 *
 * Evaluates potential secondary structure in the amplicon region that
 * could cause sequencing read-through problems. High GC content and
 * GC-rich regions can form stable secondary structures.
 *
 * Factors considered:
 * - Overall GC content of amplicon
 * - Presence of GC-rich stretches (>70% GC in 20bp window)
 * - Potential hairpin-forming sequences
 *
 * @param sequence - Amplicon sequence
 * @param options - Scoring parameters
 * @returns Score from 0 to 1
 */
export function scoreAmpliconStructure(sequence: string, options: AmpliconStructureOptions = {}): number {
  const {
    gcThreshold = 70,        // GC% above which structure is concerning
    windowSize = 20,         // Window for GC-rich stretch detection
    maxGcStretches = 3,      // Maximum acceptable GC-rich stretches
    homopolymerThreshold = 5, // Runs of same base that cause issues
  } = options;

  if (!sequence || sequence.length < windowSize) {
    return 1.0;  // Too short to evaluate
  }

  sequence = sequence.toUpperCase();
  let score = 1.0;

  // Count GC-rich windows
  let gcRichWindows = 0;
  for (let i = 0; i <= sequence.length - windowSize; i++) {
    const window = sequence.slice(i, i + windowSize);
    const gcCount = (window.match(/[GC]/g) || []).length;
    const gcPercent = (gcCount / windowSize) * 100;
    if (gcPercent >= gcThreshold) {
      gcRichWindows++;
    }
  }

  // Penalize GC-rich stretches
  if (gcRichWindows > 0) {
    const penalty = Math.min(0.4, gcRichWindows * 0.1);
    score -= penalty;
  }

  // Check for long homopolymer runs that cause sequencing issues
  let maxRun = 1;
  let currentRun = 1;
  for (let i = 1; i < sequence.length; i++) {
    if (sequence[i] === sequence[i - 1]) {
      currentRun++;
      maxRun = Math.max(maxRun, currentRun);
    } else {
      currentRun = 1;
    }
  }

  if (maxRun >= homopolymerThreshold) {
    const runPenalty = Math.min(0.3, (maxRun - homopolymerThreshold + 1) * 0.1);
    score -= runPenalty;
  }

  // Overall GC content penalty for very high GC amplicons
  const overallGC = ((sequence.match(/[GC]/g) || []).length / sequence.length) * 100;
  if (overallGC > 65) {
    score -= Math.min(0.2, (overallGC - 65) * 0.01);
  }

  return Math.max(0, Math.min(1, score));
}

/**
 * Calculate composite score from individual feature scores
 *
 * Uses empirically-derived weights from meta-analysis
 *
 * @param scores - Object with individual feature scores (0-1)
 * @param weights - Optional custom weights
 * @returns Composite score and breakdown
 */
export function calculateCompositeScore(
  scores: CompositeScoreInput,
  weights: Record<string, number> | null = null
): CompositeScoreResult {
  // Use calibrated weights from weightCalibration.js
  // Optimized via grid search on Döring dataset (829 pairs)
  // Performance: F1=81.9%, AUC=0.848
  //
  // Key findings from calibration:
  // - offTarget is the most discriminative feature (+0.515 diff)
  // - terminal3DG is second most important (+0.487 diff)
  // - Higher Tm correlates with failure (-0.185 diff)
  const w = weights || DEFAULT_WEIGHTS;

  let totalScore = 0;
  let totalWeight = 0;
  const breakdown: Record<string, ScoreBreakdownItem> = {};

  for (const [key, score] of Object.entries(scores)) {
    const weight = w[key] || 0;
    if (weight > 0 && typeof score === 'number') {
      totalScore += score * weight;
      totalWeight += weight;
      breakdown[key] = {
        score: Math.round(score * 1000) / 1000,
        weight,
        contribution: Math.round(score * weight * 1000) / 1000,
      };
    }
  }

  // Normalize if total weight != 1
  const normalizedScore = totalWeight > 0 ? totalScore / totalWeight : 0;

  return {
    score: Math.round(normalizedScore * 100),  // 0-100 scale
    rawScore: Math.round(normalizedScore * 1000) / 1000,
    breakdown,
    totalWeight: Math.round(totalWeight * 1000) / 1000,
  };
}

/**
 * Classify quality tier based on score
 *
 * @param score - Score from 0-100
 * @returns Quality tier information
 */
export function classifyQuality(score: number): QualityTier {
  if (score >= 90) {
    return {
      tier: 'excellent',
      label: 'Excellent',
      description: 'High confidence primer pair',
      color: 'green',
    };
  }
  if (score >= 75) {
    return {
      tier: 'good',
      label: 'Good',
      description: 'Reliable primer pair',
      color: 'blue',
    };
  }
  if (score >= 60) {
    return {
      tier: 'acceptable',
      label: 'Acceptable',
      description: 'May work but consider alternatives',
      color: 'yellow',
    };
  }
  if (score >= 40) {
    return {
      tier: 'marginal',
      label: 'Marginal',
      description: 'Likely to have issues',
      color: 'orange',
    };
  }
  return {
    tier: 'poor',
    label: 'Poor',
    description: 'High failure risk - redesign recommended',
    color: 'red',
  };
}

/**
 * Build composite input object from primer scores
 *
 * SINGLE SOURCE OF TRUTH for mapping individual primer scores to the
 * composite score input format that matches DEFAULT_WEIGHTS keys.
 *
 * This function ensures consistency between:
 * - primerAnalysis.js (analyzePrimers, analyzeSinglePrimer)
 * - scoring.js (analyzePrimerPair)
 * - primers.js (calculatePairCompositeScore)
 *
 * @param fwdScores - Forward primer scores object
 * @param revScores - Reverse primer scores object (optional, null for single primer)
 * @param pairScores - Pair-level scores (optional)
 * @returns Composite input object matching DEFAULT_WEIGHTS keys
 */
export function buildCompositeInput(
  fwdScores: PrimerScores,
  revScores: PrimerScores | null = null,
  pairScores: PairScores | null = null
): CompositeScoreInput {
  // Single primer case
  if (!revScores) {
    return {
      // Tm and GC
      tmFwd: fwdScores.tm,
      gcFwd: fwdScores.gc,
      // Length
      lengthFwd: fwdScores.length,
      // Secondary structure
      hairpinFwd: fwdScores.hairpin,
      selfDimerFwd: fwdScores.homodimer,
      // 3' end stability (critical for extension)
      terminal3DG: fwdScores.terminal3DG,
      gcClampFwd: fwdScores.gcClamp,
      threePrimeCompFwd: fwdScores.threePrimeComp,
      // Specificity
      offTarget: fwdScores.offTarget,
      // G-quadruplex
      gQuadruplexFwd: fwdScores.gQuadruplex,
      // Homopolymer runs
      homopolymerFwd: fwdScores.homopolymer,
    };
  }

  // Primer pair case
  return {
    // Tm scores
    tmFwd: fwdScores.tm,
    tmRev: revScores.tm,
    // GC content scores
    gcFwd: fwdScores.gc,
    gcRev: revScores.gc,
    // Length scores
    lengthFwd: fwdScores.length,
    lengthRev: revScores.length,
    // Secondary structure - hairpin
    hairpinFwd: fwdScores.hairpin,
    hairpinRev: revScores.hairpin,
    // Secondary structure - self-dimer
    selfDimerFwd: fwdScores.homodimer,
    selfDimerRev: revScores.homodimer,
    // Cross-dimer (pair-level)
    heterodimer: pairScores?.heterodimer ?? 0.8,
    // Tm difference (pair-level)
    tmDiff: pairScores?.tmDiff ?? 1.0,
    // 3' end stability - terminal ΔG (use minimum of both primers)
    terminal3DG: Math.min(fwdScores.terminal3DG, revScores.terminal3DG),
    // 3' end stability - GC clamp
    gcClampFwd: fwdScores.gcClamp,
    gcClampRev: revScores.gcClamp,
    // 3' end composition (weighted analysis)
    threePrimeCompFwd: fwdScores.threePrimeComp,
    threePrimeCompRev: revScores.threePrimeComp,
    // Specificity - off-target (use minimum of both primers)
    offTarget: Math.min(fwdScores.offTarget, revScores.offTarget),
    // G-quadruplex risk
    gQuadruplexFwd: fwdScores.gQuadruplex,
    gQuadruplexRev: revScores.gQuadruplex,
    // Homopolymer runs
    homopolymerFwd: fwdScores.homopolymer,
    homopolymerRev: revScores.homopolymer,
    // Amplicon-related (optional, passed when available)
    ...(pairScores?.ampliconLength !== undefined && { ampliconLength: pairScores.ampliconLength }),
    ...(pairScores?.ampliconStructure !== undefined && { ampliconStructure: pairScores.ampliconStructure }),
  };
}

// SCORING_PRESETS is now imported from presets.js (single source of truth)

/**
 * Unified primer pair analysis function
 *
 * Provides comprehensive analysis for a primer pair including:
 * - Individual primer scores (Tm, GC, hairpin, homodimer, G-Quadruplex, etc.)
 * - Pair-level scores (heterodimer, Tm difference, off-target)
 * - Composite score and quality tier
 * - Detailed warnings and recommendations
 *
 * This is the canonical analysis function used by:
 * - Standalone scoring mode (score())
 * - Unified primer designer
 * - Mutagenesis module
 * - Calibration scripts
 *
 * @param fwd - Forward primer { seq, tm, gc, dg, offTargetCount }
 * @param rev - Reverse primer (optional) { seq, tm, gc, dg, offTargetCount }
 * @param options - Analysis options
 * @returns Comprehensive primer pair analysis
 */
export function analyzePrimerPair(
  fwd: Primer,
  rev: Primer | null = null,
  options: AnalyzePrimerPairOptions = {}
): PrimerPairAnalysis {
  const {
    heterodimerDG = null,  // Pre-calculated if available
    hairpinDGFwd = null,   // Pre-calculated if available
    hairpinDGRev = null,
    homodimerDGFwd = null,
    homodimerDGRev = null,
    ampliconLength = null,
    includeRecommendations = true,
    // NEW: Mode/preset support for domain-specific scoring
    mode = 'amplification',
    // NEW: Custom weights override
    weights = null,
    // NEW: Additional warnings from external analysis (e.g., mutagenesis structure check)
    externalWarnings = null,
  } = options;

  // Get preset configuration
  const preset = (SCORING_PRESETS as any)[mode] || SCORING_PRESETS.amplification;
  const effectiveWeights = weights || DEFAULT_WEIGHTS;

  const warnings: string[] = [];
  const recommendations: string[] = [];

  // Forward primer analysis
  const fwdG4 = analyzeGQuadruplex(fwd.seq);
  const fwdTerminalDG = calculate3primeTerminalDG(fwd.seq).dG;
  const fwdScores: PrimerScores = {
    tm: scoreTm(fwd.tm, preset.tmOptions),
    gc: scoreGc(fwd.gc, preset.gcOptions),
    length: scoreLength(fwd.seq.length, preset.lengthOptions),
    terminal3DG: scoreTerminal3DG(fwdTerminalDG),
    gcClamp: scoreGcClamp(fwd.seq),
    homopolymer: scoreHomopolymer(fwd.seq),
    hairpin: hairpinDGFwd !== null ? scoreHairpin(hairpinDGFwd, { threshold: preset.hairpinThreshold }) : 0.8,
    homodimer: homodimerDGFwd !== null ? scoreHomodimer(homodimerDGFwd, { threshold: preset.homodimerThreshold }) : 0.8,
    offTarget: fwd.offTargetCount !== undefined ? scoreOffTarget(fwd.offTargetCount) : 0.8,
    gQuadruplex: fwdG4.score,
    // 3' end composition score (weighted analysis of GC clamp, terminal ΔG, patterns)
    threePrimeComp: score3PrimeComposition(fwd.seq, fwdTerminalDG),
  };

  // Collect warnings from forward primer
  if (fwdG4.severity === 'critical') {
    warnings.push(`Forward: ${fwdG4.message}`);
    if (includeRecommendations) {
      recommendations.push('Redesign forward primer to avoid G-quadruplex motif');
    }
  } else if (fwdG4.severity === 'warning') {
    warnings.push(`Forward: ${fwdG4.message}`);
  }

  if (fwdScores.terminal3DG < 0.5) {
    warnings.push('Forward: Weak 3\' terminal binding');
    if (includeRecommendations) {
      recommendations.push('Consider extending forward primer for stronger 3\' end');
    }
  }

  // Critical warning for extreme length violations (forward)
  const fwdLengthLow = preset.lengthOptions?.optimalLow ?? 18;
  const fwdLengthHigh = preset.lengthOptions?.optimalHigh ?? 24;
  if (fwd.seq.length > fwdLengthHigh + 15) {
    warnings.push(`CRITICAL Forward: Primer extremely long (${fwd.seq.length}bp, optimal: ${fwdLengthLow}-${fwdLengthHigh}bp)`);
  }

  // Critical warning for severe homopolymer runs (forward)
  let fwdMaxHomopolymer = 1;
  let fwdCurrentRun = 1;
  for (let i = 1; i < fwd.seq.length; i++) {
    if (fwd.seq[i].toUpperCase() === fwd.seq[i - 1].toUpperCase()) {
      fwdCurrentRun++;
      fwdMaxHomopolymer = Math.max(fwdMaxHomopolymer, fwdCurrentRun);
    } else {
      fwdCurrentRun = 1;
    }
  }
  if (fwdMaxHomopolymer >= 6) {
    warnings.push(`CRITICAL Forward: Severe homopolymer run (${fwdMaxHomopolymer} identical bases)`);
  }

  const analysis: PrimerPairAnalysis = {
    forward: {
      sequence: fwd.seq,
      length: fwd.seq.length,
      tm: fwd.tm,
      gc: fwd.gc,
      dg: fwd.dg,
      scores: fwdScores,
      gQuadruplex: fwdG4,
    },
    reverse: null,
    pair: null,
    composite: null,
    quality: null,
    warnings,
    recommendations,
  };

  // If no reverse primer, calculate single primer composite
  if (!rev) {
    // Use unified buildCompositeInput for consistency
    const compositeInput = buildCompositeInput(fwdScores);

    analysis.composite = calculateCompositeScore(compositeInput, effectiveWeights);
    analysis.quality = classifyQuality(analysis.composite.score);
    analysis.mode = mode;
    analysis.presetName = preset.name;
    return analysis;
  }

  // Reverse primer analysis
  const revG4 = analyzeGQuadruplex(rev.seq);
  const revTerminalDG = calculate3primeTerminalDG(rev.seq).dG;
  const revScores: PrimerScores = {
    tm: scoreTm(rev.tm, preset.tmOptions),
    gc: scoreGc(rev.gc, preset.gcOptions),
    length: scoreLength(rev.seq.length, preset.lengthOptions),
    terminal3DG: scoreTerminal3DG(revTerminalDG),
    gcClamp: scoreGcClamp(rev.seq),
    homopolymer: scoreHomopolymer(rev.seq),
    hairpin: hairpinDGRev !== null ? scoreHairpin(hairpinDGRev, { threshold: preset.hairpinThreshold }) : 0.8,
    homodimer: homodimerDGRev !== null ? scoreHomodimer(homodimerDGRev, { threshold: preset.homodimerThreshold }) : 0.8,
    offTarget: rev.offTargetCount !== undefined ? scoreOffTarget(rev.offTargetCount) : 0.8,
    gQuadruplex: revG4.score,
    // 3' end composition score (weighted analysis of GC clamp, terminal ΔG, patterns)
    threePrimeComp: score3PrimeComposition(rev.seq, revTerminalDG),
  };

  // Collect warnings from reverse primer
  if (revG4.severity === 'critical') {
    warnings.push(`Reverse: ${revG4.message}`);
    if (includeRecommendations) {
      recommendations.push('Redesign reverse primer to avoid G-quadruplex motif');
    }
  } else if (revG4.severity === 'warning') {
    warnings.push(`Reverse: ${revG4.message}`);
  }

  if (revScores.terminal3DG < 0.5) {
    warnings.push('Reverse: Weak 3\' terminal binding');
    if (includeRecommendations) {
      recommendations.push('Consider extending reverse primer for stronger 3\' end');
    }
  }

  // Critical warning for extreme length violations (reverse)
  const revLengthLow = preset.lengthOptions?.optimalLow ?? 18;
  const revLengthHigh = preset.lengthOptions?.optimalHigh ?? 24;
  if (rev.seq.length > revLengthHigh + 15) {
    warnings.push(`CRITICAL Reverse: Primer extremely long (${rev.seq.length}bp, optimal: ${revLengthLow}-${revLengthHigh}bp)`);
  }

  // Critical warning for severe homopolymer runs (reverse)
  let revMaxHomopolymer = 1;
  let revCurrentRun = 1;
  for (let i = 1; i < rev.seq.length; i++) {
    if (rev.seq[i].toUpperCase() === rev.seq[i - 1].toUpperCase()) {
      revCurrentRun++;
      revMaxHomopolymer = Math.max(revMaxHomopolymer, revCurrentRun);
    } else {
      revCurrentRun = 1;
    }
  }
  if (revMaxHomopolymer >= 6) {
    warnings.push(`CRITICAL Reverse: Severe homopolymer run (${revMaxHomopolymer} identical bases)`);
  }

  analysis.reverse = {
    sequence: rev.seq,
    length: rev.seq.length,
    tm: rev.tm,
    gc: rev.gc,
    dg: rev.dg,
    scores: revScores,
    gQuadruplex: revG4,
  };

  // Pair-level analysis
  const tmDiff = Math.abs(fwd.tm - rev.tm);
  const pairScores: PairScores = {
    tmDiff: scoreTmDiff(fwd.tm, rev.tm),
    heterodimer: heterodimerDG !== null ? scoreHeterodimer(heterodimerDG, { threshold: preset.heterodimerThreshold }) : 0.8,
    ampliconLength: ampliconLength !== null ? scoreAmpliconLength(ampliconLength) : 0.9,
  };

  if (tmDiff > 5) {
    // Critical if Tm difference > 8°C (severe mismatch)
    const tmPrefix = tmDiff > 8 ? 'CRITICAL ' : '';
    warnings.push(`${tmPrefix}Tm difference: ${tmDiff.toFixed(1)}°C (>5°C may cause unequal amplification)`);
    if (includeRecommendations) {
      recommendations.push('Adjust primer lengths to match Tm values');
    }
  }

  analysis.pair = {
    tmDiff,
    heterodimerDG,
    ampliconLength,
    scores: pairScores,
  };

  // Calculate composite score using unified buildCompositeInput
  const compositeInput = buildCompositeInput(fwdScores, revScores, pairScores);

  analysis.composite = calculateCompositeScore(compositeInput, effectiveWeights);

  // Apply critical warning penalty to effective score (20 points per critical warning)
  const criticalWarnings = warnings.filter(w => w.startsWith('CRITICAL')).length;
  const effectiveScore = Math.max(0, analysis.composite.score - criticalWarnings * 20);
  analysis.quality = classifyQuality(effectiveScore);
  analysis.criticalWarnings = criticalWarnings;
  analysis.effectiveScore = effectiveScore;

  // Add mode and preset info
  analysis.mode = mode;
  analysis.presetName = preset.name;

  // Add external warnings (from domain-specific analysis like mutagenesis structure check)
  if (externalWarnings && Array.isArray(externalWarnings)) {
    warnings.push(...externalWarnings);
  }

  return analysis;
}

/**
 * Analyze a batch of primer pairs for statistics and calibration
 *
 * Provides dataset-level statistics including:
 * - Success rates by category
 * - Feature discrimination analysis (which features distinguish success/failure)
 * - G-Quadruplex prevalence and impact
 * - Score distributions
 *
 * @param dataset - Array of { scores, success } entries
 * @param weights - Weight configuration for scoring
 * @returns Batch analysis statistics
 */
export function analyzePrimerBatch(
  dataset: BatchDatasetEntry[],
  weights: Record<string, number> = DEFAULT_WEIGHTS
): BatchAnalysisStats {
  const stats: BatchAnalysisStats = {
    total: dataset.length,
    success: 0,
    failure: 0,
    successRate: 0,
    scoreDistribution: {
      success: { scores: [], mean: 0, std: 0 },
      failure: { scores: [], mean: 0, std: 0 },
    },
    featureDiscrimination: {},
    topDiscriminativeFeatures: [],
    gQuadruplex: {
      g4Motif: { count: 0, successCount: 0, successRate: 0 },
      ggggRun: { count: 0, successCount: 0, successRate: 0 },
      multiGgg: { count: 0, successCount: 0, successRate: 0 },
      noRisk: { count: 0, successCount: 0, successRate: 0 },
    },
  };

  const featureValues: Record<string, { success: number[]; failure: number[] }> = {};

  for (const entry of dataset) {
    const success = entry.success || entry.actual || false;

    if (success) {
      stats.success++;
    } else {
      stats.failure++;
    }

    // Calculate composite score
    let totalScore = 0;
    let totalWeight = 0;

    for (const [key, score] of Object.entries(entry.scores)) {
      const weight = weights[key] || 0;
      if (weight > 0 && typeof score === 'number') {
        totalScore += score * weight;
        totalWeight += weight;

        // Track feature values for discrimination analysis
        if (!featureValues[key]) {
          featureValues[key] = { success: [], failure: [] };
        }
        if (success) {
          featureValues[key].success.push(score);
        } else {
          featureValues[key].failure.push(score);
        }
      }
    }

    const compositeScore = totalWeight > 0 ? (totalScore / totalWeight) * 100 : 0;
    if (success) {
      stats.scoreDistribution.success.scores.push(compositeScore);
    } else {
      stats.scoreDistribution.failure.scores.push(compositeScore);
    }

    // G-Quadruplex analysis
    const g4Score = (entry.scores.gQuadruplexFwd ?? entry.scores.gQuadruplexRev) ?? 1.0;
    if (g4Score === 0.0) {
      stats.gQuadruplex.g4Motif.count++;
      if (success) stats.gQuadruplex.g4Motif.successCount++;
    } else if (g4Score <= 0.2) {
      stats.gQuadruplex.ggggRun.count++;
      if (success) stats.gQuadruplex.ggggRun.successCount++;
    } else if (g4Score <= 0.6) {
      stats.gQuadruplex.multiGgg.count++;
      if (success) stats.gQuadruplex.multiGgg.successCount++;
    } else {
      stats.gQuadruplex.noRisk.count++;
      if (success) stats.gQuadruplex.noRisk.successCount++;
    }
  }

  // Calculate success rate
  stats.successRate = stats.total > 0 ? stats.success / stats.total : 0;

  // Calculate score distribution statistics
  const calcStats = (arr: number[]): { mean: number; std: number } => {
    if (arr.length === 0) return { mean: 0, std: 0 };
    const mean = arr.reduce((a, b) => a + b, 0) / arr.length;
    const std = Math.sqrt(arr.reduce((sum, x) => sum + Math.pow(x - mean, 2), 0) / arr.length);
    return {
      mean: Math.round(mean * 10) / 10,
      std: Math.round(std * 10) / 10,
    };
  };

  const successStats = calcStats(stats.scoreDistribution.success.scores);
  const failureStats = calcStats(stats.scoreDistribution.failure.scores);
  stats.scoreDistribution.success.mean = successStats.mean;
  stats.scoreDistribution.success.std = successStats.std;
  stats.scoreDistribution.failure.mean = failureStats.mean;
  stats.scoreDistribution.failure.std = failureStats.std;

  // Calculate feature discrimination
  for (const [key, values] of Object.entries(featureValues)) {
    const successMean = values.success.length > 0
      ? values.success.reduce((a, b) => a + b, 0) / values.success.length
      : 0;
    const failureMean = values.failure.length > 0
      ? values.failure.reduce((a, b) => a + b, 0) / values.failure.length
      : 0;

    stats.featureDiscrimination[key] = {
      successMean: Math.round(successMean * 1000) / 1000,
      failureMean: Math.round(failureMean * 1000) / 1000,
      diff: Math.round((successMean - failureMean) * 1000) / 1000,
    };
  }

  // Sort features by discrimination power
  stats.topDiscriminativeFeatures = Object.entries(stats.featureDiscrimination)
    .sort((a, b) => Math.abs(b[1].diff) - Math.abs(a[1].diff))
    .slice(0, 10)
    .map(([key, data]) => ({ feature: key, ...data }));

  // Calculate G-Quadruplex success rates
  for (const category of Object.keys(stats.gQuadruplex) as Array<keyof typeof stats.gQuadruplex>) {
    const cat = stats.gQuadruplex[category];
    cat.successRate = cat.count > 0
      ? Math.round((cat.successCount / cat.count) * 1000) / 1000
      : null;
  }

  return stats;
}

// Re-export DEFAULT_WEIGHTS for convenience
export { DEFAULT_WEIGHTS };
