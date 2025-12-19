/**
 * Unified Assembly Core Module
 *
 * This module provides a unified API for DNA assembly methods:
 * - Gibson Assembly
 * - NEBuilder HiFi Assembly
 * - Golden Gate Assembly
 *
 * ARCHITECTURE PRINCIPLE: This module WRAPS existing primer design functions
 * without modifying them. It imports from:
 * - primers.ts (primer design core)
 * - tmQ5.ts (Q5 Tm calculation)
 * - equilibrium.ts (thermodynamic analysis)
 * - repp/goldengate.js (Golden Gate specifics)
 * - repp/assembly.js (Gibson assembly planning)
 *
 * This allows the existing amplification/mutagenesis tools to remain untouched.
 */

import { primers } from './primers.js';
import { calculateTmQ5, calculateGC, calculate3primeTerminalDG } from './tmQ5.js';
import { calculateHairpinDG, calculateHomodimerDG, calculateHeterodimerDG } from './equilibrium.js';
import { dg as foldDG } from './fold.js';
import {
  classifyQuality,
  piecewiseLogistic,
  scoreTm,
  scoreGc,
  scoreLength,
  scoreHairpin,
  scoreGcClamp,
  scoreTerminal3DG,
  scoreHomopolymer,
  score3PrimeComposition,
  calculateCompositeScore,
  analyzeGQuadruplex,
} from './scoring.js';
import { reverseComplement } from './sequenceUtils.js';
import { ANALYSIS_PRESETS } from './presets.js';
import { ASSEMBLY_WEIGHTS, OVERLAP_WEIGHTS, ANNEALING_SINGLE_WEIGHTS } from './weightCalibration.js';

// =============================================================================
// Types and Interfaces
// =============================================================================

export interface AssemblyMethod {
  name: string;
  description: string;
  overlapRange: { min: number; max: number; optimal: number };
  overlapTm?: { min: number; max: number; optimal: number };
  maxFragments: number;
  enzyme: string;
  incubationTemp: number;
  incubationTime: number;
  reference: string;
}

export interface AssemblyConfig {
  method: string;
  primerTmTarget: number;
  primerTmRange: number;
  overlapTmTarget: number;
  overlapTmRange: number;
  primerLenMin: number;
  primerLenMax: number;
  overlapLenMin: number;
  overlapLenMax: number;
  gcMin: number;
  gcMax: number;
  maxHairpinDG: number;
  maxHomodimerDG: number;
  maxHeterodimerDG: number;
  autoOptimize: boolean;
  optimizeFor: string;
  junctionSlide?: number;
  circular?: boolean;
}

export interface OverlapCandidate {
  sequence: string;
  length: number;
  offset?: number;
  tm: number;
  gc: number;
  gcPercent: string;
  terminal3DG?: number;
  hairpinDG: number;
  score?: number;
  compositeScore?: number;
  qualityTier?: string;
  isValid?: boolean;
  warnings: string[];
  scores?: Record<string, number>;
  breakdown?: {
    primary: {
      tm: { value: number; score: number; inRange: boolean };
      gc: { value: number; score: number; inRange: boolean };
      length: { value: number; score: number; inRange: boolean };
    };
    quality: Record<string, number>;
  };
}

export interface OverlapOptimization {
  optimal: OverlapCandidate;
  alternatives: OverlapCandidate[];
  allCandidates?: OverlapCandidate[];
  validCandidates?: OverlapCandidate[];
  method: string;
  scoringMethod?: string;
  optimizationSummary?: {
    totalEvaluated: number;
    validCount: number;
    bestScore: number;
    scoreDistribution: Record<string, number>;
  };
}

export interface Fragment {
  id: string;
  seq: string;
}

export interface Junction {
  junction: number;
  from: string;
  to: string;
  overlap: OverlapCandidate;
  alternatives: OverlapCandidate[];
  optimizationSummary?: {
    totalEvaluated: number;
    validCount: number;
    bestScore: number;
    scoreDistribution: Record<string, number>;
  };
}

export interface AssemblyOverlaps {
  method: string;
  fragments: { id: string; length: number }[];
  junctions: Junction[];
  quality: {
    averageScore: number;
    minimumScore: number;
    weakestJunction?: number;
    tier: string;
  };
  warnings: string[];
  circular: boolean;
  optimizationMethod: string;
}

export interface AssemblyPrimers {
  forward: {
    sequence: string;
    length: number;
    annealingRegion: string;
    annealingLength: number;
    homologyTail: string;
    tailLength: number;
    tm: number;
    tmFull: number;
    gc: number;
    gcPercent: string;
    hairpinDG: number;
    qualityScore?: number;
    qualityTier?: string;
    scoreBreakdown?: Record<string, number>;
  };
  reverse: {
    sequence: string;
    length: number;
    annealingRegion: string;
    annealingLength: number;
    homologyTail: string;
    tailLength: number;
    tm: number;
    tmFull: number;
    gc: number;
    gcPercent: string;
    hairpinDG: number;
    qualityScore?: number;
    qualityTier?: string;
    scoreBreakdown?: Record<string, number>;
  };
  pair: {
    annealingTemp: number;
    tmDifference: number;
    heterodimerDG: number;
    pairScore?: number;
    qualityTier?: string;
  };
  pcr: {
    annealingTemp: number;
    extensionTime: number;
    cycles: number;
  };
  productLength: number;
  warnings: string[];
  method: string;
  alternatives?: Array<{
    forward: AnnealingCandidate;
    reverse: AnnealingCandidate;
    pairScore: number;
    tmDiff: number;
    heterodimerDG: number;
    selectionReason: string;
  }>;
  optimization?: {
    method: string;
    candidatesEvaluated: {
      forward: number;
      reverse: number;
      pairs: number;
    };
  };
}

export interface AnnealingCandidate {
  sequence: string;
  length: number;
  offset?: number;
  tm: number;
  gc: number;
  terminal3DG: number;
  hairpinDG: number;
  scores: Record<string, number>;
  compositeScore: number;
  qualityTier: string;
  isValid: boolean;
  warnings: string[];
  isFwd: boolean;
}

export interface AssemblyDesign {
  method: string;
  methodKey: string;
  fragments: Array<{
    id: string;
    index: number;
    originalLength: number;
    primers: AssemblyPrimers;
    leftOverlap: { sequence: string; from: string };
    rightOverlap: { sequence: string; to: string };
  }>;
  junctions: Junction[];
  quality: {
    averageScore: number;
    minimumScore: number;
    weakestJunction?: number;
    tier: string;
  };
  assembly: {
    circular: boolean;
    totalLength: number;
    fragmentCount: number;
  };
  protocol: Protocol;
  cost: {
    primers: number;
    assembly: number;
    total: number;
  };
  warnings: string[];
  optimization: {
    method: string;
    autoOptimize: boolean;
    description: string;
  };
}

export interface Protocol {
  title: string;
  steps: Array<{
    step: number;
    title: string;
    details: string[];
  }>;
  tips: string[];
}

export interface SimulationResult {
  success: boolean;
  method: string;
  assembledSequence: string;
  assembledLength: number;
  gcContent: number;
  gcPercent: string;
  circular: boolean;
  circularClosure: {
    valid: boolean;
    expectedOverlap: string;
    actualOverlap: string;
    matches: boolean;
    message: string;
  } | null;
  steps: Array<{
    step: number;
    fragment: string;
    action: string;
    overlapSequence?: string;
    overlapLength?: number;
    addedLength?: number;
    position: number;
    length?: number;
    cumulativeLength: number;
  }>;
  homologyConflicts: Array<{
    fragment1: string;
    fragment1Index: number;
    fragment2: string;
    fragment2Index: number;
    homologySequence: string;
    homologyLength: number;
    position1: number;
    position2: number;
    severity: string;
    warning: string;
  }>;
  hasConflicts: boolean;
  visualization: VisualizationData;
  fragments: Array<{
    id: string;
    length: number;
    position: number;
    gcContent: number;
  }>;
}

export interface VisualizationData {
  type: string;
  totalLength: number;
  fragments: Array<{
    id: string;
    index: number;
    startPosition: number;
    endPosition: number;
    displayLength: number;
    actualLength: number;
    overlapWithNext: number;
    color: string;
    startAngle: number | null;
    endAngle: number | null;
  }>;
  junctions: Array<{
    index: number;
    from: string;
    to: string;
    position: number;
    overlapSequence?: string;
    quality?: string;
    score?: number;
  }>;
  svg: {
    width: number;
    height: number;
    radius: number | null;
    centerX: number | null;
    centerY: number | null;
  };
}

export interface JunctionValidation {
  valid: boolean;
  junctions: Array<{
    junction: number;
    from: string;
    to: string;
    hasNaturalHomology: boolean;
    naturalHomology: {
      sequence: string;
      length: number;
      hasNaturalHomology: boolean;
    } | null;
    requiresDesignedTails: boolean;
    message: string;
  }>;
  hasAllNaturalHomology: boolean;
  requiresDesignedTails: boolean;
  summary: string;
}

// =============================================================================
// Constants
// =============================================================================

/**
 * Assembly method configurations
 * Based on NEB recommendations and published protocols
 */
export const ASSEMBLY_METHODS: Record<string, AssemblyMethod> = {
  GIBSON: {
    name: 'Gibson Assembly',
    description: 'Isothermal assembly using exonuclease, polymerase, and ligase',
    overlapRange: { min: 15, max: 40, optimal: 25 },
    overlapTm: { min: 48, max: 65, optimal: 55 },
    maxFragments: 6,
    enzyme: 'T5 Exonuclease + Phusion + Taq Ligase',
    incubationTemp: 50,
    incubationTime: 60, // minutes
    reference: 'Gibson et al. 2009, Nature Methods',
  },
  NEBUILDER_HIFI: {
    name: 'NEBuilder HiFi DNA Assembly',
    description: 'NEB proprietary high-fidelity assembly mix',
    overlapRange: { min: 15, max: 35, optimal: 20 },
    overlapTm: { min: 48, max: 65, optimal: 55 },
    maxFragments: 5,
    enzyme: 'NEBuilder HiFi Master Mix',
    incubationTemp: 50,
    incubationTime: 15, // minutes for 2-3 fragments
    reference: 'NEB E5520',
  },
  GOLDEN_GATE: {
    name: 'Golden Gate Assembly',
    description: 'Type IIS restriction enzyme-based scarless assembly',
    overlapRange: { min: 4, max: 4, optimal: 4 }, // Overhang length
    maxFragments: 20,
    enzyme: 'BsaI-HFv2 + T4 DNA Ligase',
    incubationTemp: 37,
    incubationTime: 60,
    reference: 'Engler et al. 2008, PLOS ONE',
  },
};

/**
 * Default assembly configuration
 */
export const DEFAULT_ASSEMBLY_CONFIG: AssemblyConfig = {
  method: 'NEBUILDER_HIFI',
  primerTmTarget: 60,       // Target Tm for primer annealing region
  primerTmRange: 5,         // Acceptable Tm range
  overlapTmTarget: 55,      // Target Tm for overlap region
  overlapTmRange: 5,        // Acceptable overlap Tm range
  primerLenMin: 18,
  primerLenMax: 32,
  overlapLenMin: 15,
  overlapLenMax: 35,
  gcMin: 0.40,
  gcMax: 0.60,
  maxHairpinDG: -3,         // kcal/mol
  maxHomodimerDG: -6,       // kcal/mol
  maxHeterodimerDG: -6,     // kcal/mol
  // State-of-the-art optimization options
  autoOptimize: false,      // Use range-based optimization instead of target-based
  optimizeFor: 'quality',   // 'quality' | 'cost' | 'robustness'
};

// =============================================================================
// Sequence Utilities
// =============================================================================

/**
 * Calculate GC content
 */
function gcContent(seq: string): number {
  if (!seq || seq.length === 0) return 0;
  const gc = (seq.toUpperCase().match(/[GC]/g) || []).length;
  return gc / seq.length;
}

// =============================================================================
// Overlap Optimization (NEBuilder/Gibson)
// =============================================================================

/**
 * Find optimal overlap region between two fragments
 *
 * This is a key algorithm for Gibson/NEBuilder assemblies.
 * The overlap must:
 * 1. Have appropriate Tm (48-65°C for NEBuilder)
 * 2. Have balanced GC content (40-60%)
 * 3. Avoid secondary structure
 * 4. Not end with poly-T (poor polymerase extension)
 */
export function findOptimalOverlap(
  fragment1End: string,
  fragment2Start: string,
  config: Partial<AssemblyConfig> = {}
): OverlapOptimization {
  const cfg: AssemblyConfig = { ...DEFAULT_ASSEMBLY_CONFIG, ...config };
  const method = ASSEMBLY_METHODS[cfg.method] || ASSEMBLY_METHODS.NEBUILDER_HIFI;

  const minLen = method.overlapRange.min;
  const maxLen = method.overlapRange.max;
  const optimalLen = method.overlapRange.optimal;
  const targetTm = method.overlapTm!.optimal;

  const candidates: OverlapCandidate[] = [];

  // Try different overlap lengths
  for (let len = minLen; len <= maxLen; len++) {
    // The overlap is at the junction: end of frag1 = start of frag2
    // For Gibson/NEBuilder, primers add homology tails
    // The overlap sequence is from fragment2Start (what we're adding as a tail)
    const overlapSeq = fragment2Start.slice(0, len).toUpperCase();

    if (overlapSeq.length < len) continue;

    // Calculate properties
    const tm = calculateTmQ5(overlapSeq);
    const gc = gcContent(overlapSeq);

    // Check for problematic patterns
    const hasPolyT = /TTTT/.test(overlapSeq);
    const hasPolyA = /AAAA/.test(overlapSeq);
    const hasPolyG = /GGGG/.test(overlapSeq); // G-quadruplex risk
    const startsWithAT = /^[AT]/.test(overlapSeq);
    const endsWithGC = /[GC]$/.test(overlapSeq);

    // Calculate secondary structure
    let hairpinDG = 0;
    try {
      hairpinDG = foldDG(overlapSeq, 50);
    } catch (e) {
      hairpinDG = 0;
    }

    // Score the overlap
    let score = 100;

    // Tm scoring (most important)
    const tmDiff = Math.abs(tm - targetTm);
    if (tmDiff <= 2) score -= 0;
    else if (tmDiff <= 5) score -= (tmDiff - 2) * 3;
    else score -= 9 + (tmDiff - 5) * 5;

    // GC content scoring
    if (gc < cfg.gcMin) score -= (cfg.gcMin - gc) * 50;
    if (gc > cfg.gcMax) score -= (gc - cfg.gcMax) * 50;

    // Length preference (prefer optimal length)
    const lenDiff = Math.abs(len - optimalLen);
    score -= lenDiff * 0.5;

    // Penalties for problematic patterns
    if (hasPolyT) score -= 15;  // Poor 3' extension
    if (hasPolyA) score -= 10;
    if (hasPolyG) score -= 20;  // G-quadruplex risk
    if (startsWithAT) score -= 3;  // Weaker 5' end
    if (!endsWithGC) score -= 5;   // Prefer GC at 3' end

    // Secondary structure penalty
    if (hairpinDG < -2) score -= Math.abs(hairpinDG) * 3;

    // Store candidate
    const warnings: string[] = [];
    if (hasPolyT) warnings.push('Contains TTTT - may cause poor extension');
    if (hasPolyG) warnings.push('Contains GGGG - G-quadruplex risk');
    if (tm < method.overlapTm!.min) warnings.push(`Tm ${tm.toFixed(1)}°C below minimum ${method.overlapTm!.min}°C`);
    if (tm > method.overlapTm!.max) warnings.push(`Tm ${tm.toFixed(1)}°C above maximum ${method.overlapTm!.max}°C`);

    candidates.push({
      sequence: overlapSeq,
      length: len,
      tm: Math.round(tm * 10) / 10,
      gc: Math.round(gc * 1000) / 10,
      gcPercent: `${(gc * 100).toFixed(1)}%`,
      hairpinDG: Math.round(hairpinDG * 10) / 10,
      score: Math.round(score * 10) / 10,
      warnings,
    });
  }

  // Sort by score (highest first)
  candidates.sort((a, b) => (b.score || 0) - (a.score || 0));

  const best = candidates[0];

  return {
    optimal: best,
    alternatives: candidates.slice(1, 5),
    allCandidates: candidates,
    method: method.name,
  };
}

// =============================================================================
// State-of-the-Art Range-Based Optimization
// =============================================================================
// Uses unified weights from weightCalibration.js: OVERLAP_WEIGHTS, ASSEMBLY_WEIGHTS

/**
 * Score overlap Tm using range-based approach
 *
 * Unlike target-based scoring, this gives EQUAL scores to any Tm within optimal range.
 * Scoring breakdown:
 * - 1.0: Within optimal range [48-65°C for NEBuilder]
 * - 0.7-1.0: Within acceptable range (linear decay)
 * - <0.7: Outside acceptable range (steep penalty)
 */
function scoreOverlapTmRange(tm: number, method: AssemblyMethod): number {
  const optimal = method.overlapTm || { min: 48, max: 65, optimal: 55 };

  return piecewiseLogistic(tm, {
    optimalLow: optimal.min,
    optimalHigh: optimal.max,
    acceptableLow: optimal.min - 5,   // 43°C absolute minimum
    acceptableHigh: optimal.max + 3,  // 68°C absolute maximum
    steepness: 0.5,
    floorScore: 0.1,
  });
}

/**
 * Score overlap length using range-based approach
 */
function scoreOverlapLengthRange(length: number, method: AssemblyMethod): number {
  const range = method.overlapRange || { min: 15, max: 35, optimal: 20 };

  return piecewiseLogistic(length, {
    optimalLow: range.min,
    optimalHigh: range.max,
    acceptableLow: range.min - 3,     // Allow slightly shorter
    acceptableHigh: range.max + 5,    // Allow slightly longer
    steepness: 0.3,
    floorScore: 0.2,
  });
}

/**
 * Score pattern avoidance for overlap sequences
 *
 * Problematic patterns:
 * - Poly-T: Poor 3' extension by polymerase
 * - Poly-G: G-quadruplex formation risk
 * - Poly-A: AT-rich instability
 * - Palindromes: Self-annealing risk
 * - Dinucleotide repeats: Replication slippage
 */
function scorePatternAvoidance(seq: string): { score: number; patterns: Record<string, boolean> } {
  let score = 1.0;
  const patterns = {
    polyT: /TTTT/.test(seq),
    polyA: /AAAA/.test(seq),
    polyG: /GGGG/.test(seq),
    polyC: /CCCC/.test(seq),
    dinucRepeat: /(.{2})\1{3,}/.test(seq),  // e.g., ATATATAT
    lowComplexity3Prime: /(.)\1{2,}$/.test(seq), // Repeat at 3' end
  };

  // Penalties based on severity
  if (patterns.polyT) score -= 0.25;  // Most problematic for extension
  if (patterns.polyG) score -= 0.30;  // G-quadruplex is severe
  if (patterns.polyA) score -= 0.15;
  if (patterns.polyC) score -= 0.15;
  if (patterns.dinucRepeat) score -= 0.20;
  if (patterns.lowComplexity3Prime) score -= 0.15;

  return { score: Math.max(0, score), patterns };
}

/**
 * Score self-complementarity (palindrome detection)
 *
 * Palindromic sequences can self-anneal, competing with target binding.
 * This is especially problematic for overlaps where both strands meet.
 */
function scoreSelfComplementarity(seq: string): number {
  const rc = reverseComplement(seq);
  let maxMatch = 0;

  // Check for self-complementary regions (sliding window)
  for (let winLen = 6; winLen <= Math.min(12, seq.length); winLen++) {
    for (let i = 0; i <= seq.length - winLen; i++) {
      const window = seq.slice(i, i + winLen);
      // Check if this window appears in the reverse complement
      if (rc.includes(window)) {
        maxMatch = Math.max(maxMatch, winLen);
      }
    }
  }

  // Score based on longest self-complementary region
  if (maxMatch >= 10) return 0.3;  // Severe - long palindrome
  if (maxMatch >= 8) return 0.5;   // Significant palindrome
  if (maxMatch >= 6) return 0.75;  // Moderate palindrome
  return 1.0;                       // No significant self-complementarity
}

/**
 * Score cross-junction hairpin potential
 *
 * Checks if the junction region (end of frag1 + start of frag2) can form
 * a hairpin structure that would interfere with assembly.
 */
function scoreCrossJunctionHairpin(frag1End: string, frag2Start: string, overlapLen: number): number {
  // Build the junction region: last 10bp of frag1 + first 10bp of overlap
  const junctionSeq = frag1End.slice(-10) + frag2Start.slice(0, 10);

  try {
    const hairpinDG = foldDG(junctionSeq, 50);  // At assembly temperature
    if (hairpinDG < -6) return 0.3;   // Severe hairpin
    if (hairpinDG < -4) return 0.6;   // Moderate hairpin
    if (hairpinDG < -2) return 0.85;  // Weak hairpin
    return 1.0;
  } catch (e) {
    return 1.0;  // Assume no hairpin if calculation fails
  }
}

/**
 * Find optimal overlap using STATE-OF-THE-ART range-based scoring
 *
 * KEY FEATURES:
 * 1. RANGE-BASED: Both Tm=57°C and Tm=55°C score 1.0 (both within optimal range)
 * 2. POSITION FLEXIBILITY: Slides junction ±10bp to find optimal local sequence
 * 3. CROSS-JUNCTION ANALYSIS: Checks hairpin at fragment1-fragment2 boundary
 * 4. SELF-COMPLEMENTARITY: Detects palindromes that cause mispriming
 *
 * This frees up scoring capacity to differentiate by QUALITY factors:
 * - 3' end stability (GC clamp, terminal ΔG)
 * - Secondary structure avoidance
 * - Pattern avoidance (poly runs, palindromes)
 */
export function findOptimalOverlapRangeBased(
  fragment1End: string,
  fragment2Start: string,
  config: Partial<AssemblyConfig> = {}
): OverlapOptimization {
  const cfg: AssemblyConfig = { ...DEFAULT_ASSEMBLY_CONFIG, ...config };
  const method = ASSEMBLY_METHODS[cfg.method] || ASSEMBLY_METHODS.NEBUILDER_HIFI;

  const minLen = method.overlapRange.min;
  const maxLen = method.overlapRange.max;

  // Position flexibility: how far to slide the junction (bp)
  const maxSlide = cfg.junctionSlide || 10;

  const candidates: OverlapCandidate[] = [];

  // Try all lengths within the acceptable range
  for (let len = minLen; len <= maxLen; len++) {
    // POSITION FLEXIBILITY: Try sliding the junction position
    // offset=0: standard position (overlap starts at frag2 position 0)
    // offset>0: slide into frag2 (skip first N bases)
    // offset<0: would require extending into frag1 (not implemented - needs frag1 sequence)
    for (let offset = 0; offset <= maxSlide; offset++) {
      // Check if we have enough sequence
      if (offset + len > fragment2Start.length) continue;

      const overlapSeq = fragment2Start.slice(offset, offset + len).toUpperCase();
      if (overlapSeq.length < len) continue;

      // Calculate properties
      const tm = calculateTmQ5(overlapSeq);
      const gc = gcContent(overlapSeq);
      const terminal3DG = calculate3primeTerminalDG(overlapSeq)?.dG || -8;

      // Calculate secondary structure
      let hairpinDG = 0;
      try {
        hairpinDG = foldDG(overlapSeq, 50);
      } catch (e) {
        hairpinDG = 0;
      }

      // NEW: Self-complementarity scoring
      const selfCompScore = scoreSelfComplementarity(overlapSeq);

      // NEW: Cross-junction hairpin scoring
      const crossJunctionScore = scoreCrossJunctionHairpin(fragment1End, fragment2Start.slice(offset), len);

      // === RANGE-BASED SCORING ===
      // Score 1.0 for any value within optimal range
      const scores: Record<string, number> = {
        tmInRange: scoreOverlapTmRange(tm, method),
        gcInRange: scoreGc(gc, { optimalLow: 40, optimalHigh: 60 }),
        lengthInRange: scoreOverlapLengthRange(len, method),
        gcClamp: scoreGcClamp(overlapSeq),
        terminal3DG: scoreTerminal3DG(terminal3DG),
        hairpin: scoreHairpin(hairpinDG, { threshold: -3.0 }),
        patternAvoidance: scorePatternAvoidance(overlapSeq).score,
        selfComplementarity: selfCompScore,
        crossJunctionHairpin: crossJunctionScore,
        // Prefer middle of GC range (45-55%) over edges
        balancedGC: gc >= 0.45 && gc <= 0.55 ? 1.0 : 0.85,
      };

      // Use unified calculateCompositeScore with OVERLAP_WEIGHTS
      const composite = calculateCompositeScore(scores, OVERLAP_WEIGHTS);
      const compositeScore = composite.score;

      // Check if within valid ranges
      const isWithinTmRange = tm >= method.overlapTm!.min && tm <= method.overlapTm!.max;
      const isWithinGcRange = gc >= 0.30 && gc <= 0.70;
      const isWithinLengthRange = len >= minLen && len <= maxLen;
      const isValid = isWithinTmRange && isWithinGcRange && isWithinLengthRange;

      // Generate warnings
      const warnings: string[] = [];
      const { patterns } = scorePatternAvoidance(overlapSeq);
      if (patterns.polyT) warnings.push('Contains TTTT - may cause poor extension');
      if (patterns.polyG) warnings.push('Contains GGGG - G-quadruplex risk');
      if (!isWithinTmRange) warnings.push(`Tm ${tm.toFixed(1)}°C outside optimal range [${method.overlapTm!.min}-${method.overlapTm!.max}°C]`);
      if (hairpinDG < -3) warnings.push(`Stable hairpin (ΔG = ${hairpinDG.toFixed(1)} kcal/mol)`);
      if (scores.gcClamp < 0.7) warnings.push('Weak 3\' end - no GC clamp');
      if (selfCompScore < 0.6) warnings.push('Self-complementary sequence - mispriming risk');
      if (crossJunctionScore < 0.7) warnings.push('Cross-junction hairpin detected');
      if (offset > 0) warnings.push(`Junction shifted +${offset}bp into fragment`);

      // Determine quality tier
      const qualityTier = compositeScore >= 80 ? 'excellent' :
                          compositeScore >= 65 ? 'good' :
                          compositeScore >= 50 ? 'acceptable' : 'marginal';

      candidates.push({
        sequence: overlapSeq,
        length: len,
        offset,  // NEW: Track position offset
        tm: Math.round(tm * 10) / 10,
        gc: Math.round(gc * 1000) / 10,
        gcPercent: `${(gc * 100).toFixed(1)}%`,
        terminal3DG: Math.round(terminal3DG * 10) / 10,
        hairpinDG: Math.round(hairpinDG * 10) / 10,
        scores,
        compositeScore,
        qualityTier,
        isValid,
        warnings,
        // Detailed breakdown for UI
        breakdown: {
          primary: {
            tm: { value: tm, score: scores.tmInRange, inRange: isWithinTmRange },
            gc: { value: gc * 100, score: scores.gcInRange, inRange: isWithinGcRange },
            length: { value: len, score: scores.lengthInRange, inRange: isWithinLengthRange },
          },
          quality: {
            gcClamp: scores.gcClamp,
            terminal3DG: scores.terminal3DG,
            hairpin: scores.hairpin,
            patterns: scores.patternAvoidance,
            selfComplementarity: scores.selfComplementarity,
            crossJunction: scores.crossJunctionHairpin,
          },
        },
      });
    }
  }

  // Sort by composite score (quality-focused)
  candidates.sort((a, b) => (b.compositeScore || 0) - (a.compositeScore || 0));

  // Filter for valid candidates (within all ranges)
  const validCandidates = candidates.filter(c => c.isValid);

  // Use best valid candidate, or best overall if none are valid
  const best = validCandidates[0] || candidates[0];

  // Get alternatives that offer different trade-offs
  const alternatives = generateDiverseOverlapAlternatives(validCandidates, best);

  return {
    optimal: best,
    alternatives,
    allCandidates: candidates,
    validCandidates,
    method: method.name,
    scoringMethod: 'range-based',
    optimizationSummary: {
      totalEvaluated: candidates.length,
      validCount: validCandidates.length,
      bestScore: best?.compositeScore || 0,
      scoreDistribution: {
        excellent: candidates.filter(c => c.qualityTier === 'excellent').length,
        good: candidates.filter(c => c.qualityTier === 'good').length,
        acceptable: candidates.filter(c => c.qualityTier === 'acceptable').length,
        marginal: candidates.filter(c => c.qualityTier === 'marginal').length,
      },
    },
  };
}

/**
 * Generate diverse alternatives that offer different trade-offs
 *
 * Instead of just top 4 by score, find alternatives that excel in different areas:
 * - Shortest valid overlap (minimizes primer cost)
 * - Best hairpin avoidance
 * - Best 3' stability
 * - Middle Tm (most robust to temperature variations)
 */
function generateDiverseOverlapAlternatives(validCandidates: OverlapCandidate[], best: OverlapCandidate): OverlapCandidate[] {
  if (validCandidates.length <= 1) return [];

  const alternatives: OverlapCandidate[] = [];
  const seen = new Set([best?.sequence]);

  // Find shortest valid overlap
  const shortestValid = validCandidates
    .filter(c => c.sequence !== best?.sequence)
    .sort((a, b) => a.length - b.length)[0];
  if (shortestValid && !seen.has(shortestValid.sequence)) {
    alternatives.push({ ...shortestValid });
    seen.add(shortestValid.sequence);
  }

  // Find best hairpin avoidance
  const bestHairpin = validCandidates
    .filter(c => !seen.has(c.sequence))
    .sort((a, b) => (b.scores?.hairpin || 0) - (a.scores?.hairpin || 0))[0];
  if (bestHairpin && (bestHairpin.scores?.hairpin || 0) > (best?.scores?.hairpin || 0)) {
    alternatives.push({ ...bestHairpin });
    seen.add(bestHairpin.sequence);
  }

  // Find best 3' stability
  const best3Prime = validCandidates
    .filter(c => !seen.has(c.sequence))
    .sort((a, b) => ((b.scores?.gcClamp || 0) + (b.scores?.terminal3DG || 0)) - ((a.scores?.gcClamp || 0) + (a.scores?.terminal3DG || 0)))[0];
  if (best3Prime) {
    alternatives.push({ ...best3Prime });
    seen.add(best3Prime.sequence);
  }

  // Find middle Tm (most robust)
  const method = ASSEMBLY_METHODS.NEBUILDER_HIFI;
  const middleTm = (method.overlapTm!.min + method.overlapTm!.max) / 2;
  const mostRobust = validCandidates
    .filter(c => !seen.has(c.sequence))
    .sort((a, b) => Math.abs(a.tm - middleTm) - Math.abs(b.tm - middleTm))[0];
  if (mostRobust) {
    alternatives.push({ ...mostRobust });
    seen.add(mostRobust.sequence);
  }

  return alternatives.slice(0, 4);
}

/**
 * Optimize all overlaps in an assembly
 *
 * When autoOptimize is enabled, uses STATE-OF-THE-ART range-based scoring
 * that maximizes quality factors when Tm/GC/length are within acceptable ranges.
 */
export function optimizeAssemblyOverlaps(fragments: Fragment[], config: Partial<AssemblyConfig> = {}): AssemblyOverlaps {
  const cfg: AssemblyConfig = { ...DEFAULT_ASSEMBLY_CONFIG, ...config };

  if (fragments.length < 2) {
    throw new Error('Assembly requires at least 2 fragments');
  }

  const junctions: Junction[] = [];
  const useRangeBased = cfg.autoOptimize === true;

  for (let i = 0; i < fragments.length; i++) {
    const frag1 = fragments[i];
    const frag2 = fragments[(i + 1) % fragments.length]; // Circular assembly

    // Get sequences for overlap analysis
    const frag1End = frag1.seq.slice(-50).toUpperCase();
    const frag2Start = frag2.seq.slice(0, 50).toUpperCase();

    // Use range-based optimization when autoOptimize is enabled
    const overlap = useRangeBased
      ? findOptimalOverlapRangeBased(frag1End, frag2Start, cfg)
      : findOptimalOverlap(frag1End, frag2Start, cfg);

    junctions.push({
      junction: i + 1,
      from: frag1.id,
      to: frag2.id,
      overlap: overlap.optimal,
      alternatives: overlap.alternatives,
      // Include optimization metadata when using range-based
      ...(useRangeBased && overlap.optimizationSummary && {
        optimizationSummary: overlap.optimizationSummary,
      }),
    });
  }

  // Calculate overall assembly quality
  // For range-based scoring, use compositeScore; for legacy, use score
  const getScore = (j: Junction) => j.overlap.compositeScore ?? j.overlap.score ?? 0;
  const avgScore = junctions.reduce((sum, j) => sum + getScore(j), 0) / junctions.length;
  const minScore = Math.min(...junctions.map(j => getScore(j)));
  const weakestJunction = junctions.find(j => getScore(j) === minScore);

  // Collect all warnings
  const warnings: string[] = [];
  for (const junction of junctions) {
    for (const warning of junction.overlap.warnings || []) {
      warnings.push(`Junction ${junction.junction} (${junction.from} → ${junction.to}): ${warning}`);
    }
  }

  return {
    method: cfg.method,
    fragments: fragments.map(f => ({ id: f.id, length: f.seq.length })),
    junctions,
    quality: {
      averageScore: Math.round(avgScore * 10) / 10,
      minimumScore: Math.round(minScore * 10) / 10,
      weakestJunction: weakestJunction?.junction,
      tier: avgScore >= 80 ? 'excellent' : avgScore >= 60 ? 'good' : 'acceptable',
    },
    warnings,
    circular: true,
    // Indicate which optimization method was used
    optimizationMethod: useRangeBased ? 'range-based' : 'target-based',
  };
}

// =============================================================================
// Assembly Primer Design
// =============================================================================

/**
 * Design primers for one fragment in an assembly
 *
 * This function WRAPS the existing primers() function and adds:
 * - Homology tails for Gibson/NEBuilder
 * - Or Type IIS sites for Golden Gate
 */
export function designAssemblyPrimers(
  fragmentSeq: string,
  context: {
    leftOverlap?: string;
    rightOverlap?: string;
    method?: string;
  },
  options: Partial<AssemblyConfig> = {}
): AssemblyPrimers {
  const {
    leftOverlap = '',
    rightOverlap = '',
    method = 'NEBUILDER_HIFI',
  } = context;

  const cfg: AssemblyConfig = { ...DEFAULT_ASSEMBLY_CONFIG, ...options };

  // Design base primers using existing primers() function (no modification!)
  const [fwd, rev] = primers(fragmentSeq, {
    optimalTm: cfg.primerTmTarget,
    tmRange: cfg.primerTmRange,
    useCompositeScore: true,
  });

  // Get the annealing portions
  const fwdAnneal = fwd.seq;
  const revAnneal = rev.seq;

  // Add homology tails for Gibson/NEBuilder
  const fwdWithTail = leftOverlap + fwdAnneal;
  const revWithTail = reverseComplement(rightOverlap) + revAnneal;

  // Calculate full primer properties
  const fwdFullTm = calculateTmQ5(fwdWithTail);
  const revFullTm = calculateTmQ5(revWithTail);
  const fwdFullGc = gcContent(fwdWithTail);
  const revFullGc = gcContent(revWithTail);

  // Calculate thermodynamic properties for tailed primers
  const fwdHairpin = calculateHairpinDG(fwdWithTail, 55);
  const revHairpin = calculateHairpinDG(revWithTail, 55);
  const heterodimerDG = calculateHeterodimerDG(fwdWithTail, revWithTail, 55);

  // Score the assembly primers
  const warnings: string[] = [];

  if (fwdHairpin < cfg.maxHairpinDG) {
    warnings.push(`Forward primer may form hairpin (ΔG = ${fwdHairpin.toFixed(1)} kcal/mol)`);
  }
  if (revHairpin < cfg.maxHairpinDG) {
    warnings.push(`Reverse primer may form hairpin (ΔG = ${revHairpin.toFixed(1)} kcal/mol)`);
  }
  if (heterodimerDG < cfg.maxHeterodimerDG) {
    warnings.push(`Primer pair may form heterodimer (ΔG = ${heterodimerDG.toFixed(1)} kcal/mol)`);
  }

  // Calculate recommended annealing temperature
  // For 2-step PCR with tailed primers: use the annealing Tm, not full Tm
  const annealingTemp = Math.round(Math.min(fwd.tm, rev.tm));

  return {
    forward: {
      sequence: fwdWithTail,
      length: fwdWithTail.length,
      annealingRegion: fwdAnneal,
      annealingLength: fwdAnneal.length,
      homologyTail: leftOverlap,
      tailLength: leftOverlap.length,
      tm: fwd.tm,                    // Annealing Tm
      tmFull: fwdFullTm,             // Full primer Tm
      gc: Math.round(fwdFullGc * 1000) / 10,
      gcPercent: `${(fwdFullGc * 100).toFixed(1)}%`,
      hairpinDG: Math.round(fwdHairpin * 10) / 10,
    },
    reverse: {
      sequence: revWithTail,
      length: revWithTail.length,
      annealingRegion: revAnneal,
      annealingLength: revAnneal.length,
      homologyTail: reverseComplement(rightOverlap),
      tailLength: rightOverlap.length,
      tm: rev.tm,                    // Annealing Tm
      tmFull: revFullTm,             // Full primer Tm
      gc: Math.round(revFullGc * 1000) / 10,
      gcPercent: `${(revFullGc * 100).toFixed(1)}%`,
      hairpinDG: Math.round(revHairpin * 10) / 10,
    },
    pair: {
      annealingTemp,
      tmDifference: Math.round(Math.abs(fwd.tm - rev.tm) * 10) / 10,
      heterodimerDG: Math.round(heterodimerDG * 10) / 10,
    },
    pcr: {
      annealingTemp,
      extensionTime: Math.ceil(fragmentSeq.length / 1000) * 30, // 30s per kb
      cycles: 30,
    },
    productLength: fragmentSeq.length,
    warnings,
    method,
  };
}

// Uses ASSEMBLY_WEIGHTS from weightCalibration.js for annealing region scoring

/**
 * Design assembly primers with RANGE-BASED optimization
 *
 * This function explores the full design space for annealing regions,
 * scoring all candidates that fall within acceptable ranges equally,
 * then differentiating by quality factors.
 *
 * Key differences from target-based:
 * - Tm of 55°C and 60°C score equally (both within assembly range)
 * - Quality factors (GC clamp, 3' stability) become differentiators
 * - Considers full primer (annealing + overlap tail) for secondary structure
 */
export function designAssemblyPrimersOptimized(
  fragmentSeq: string,
  context: {
    leftOverlap?: string;
    rightOverlap?: string;
    method?: string;
  },
  options: Partial<AssemblyConfig> = {}
): AssemblyPrimers {
  const {
    leftOverlap = '',
    rightOverlap = '',
    method = 'NEBUILDER_HIFI',
  } = context;

  const cfg: AssemblyConfig = { ...DEFAULT_ASSEMBLY_CONFIG, ...options };
  const preset = ANALYSIS_PRESETS.assembly;

  // Define acceptable ranges for annealing region
  const annealingRanges = {
    tm: { min: preset.optimalTmRange[0], max: preset.optimalTmRange[1] },       // 48-65°C for assembly
    gc: { min: 0.40, max: 0.60 },
    length: { min: preset.optimalLengthRange[0], max: preset.optimalLengthRange[1] }, // 18-30bp
  };

  // Generate all annealing region candidates for forward primer
  const fwdCandidates = generateAnnealingCandidates(
    fragmentSeq,
    true,
    annealingRanges,
    leftOverlap
  );

  // Generate all annealing region candidates for reverse primer
  const revCandidates = generateAnnealingCandidates(
    fragmentSeq,
    false,
    annealingRanges,
    rightOverlap
  );

  // Find best pair considering both annealing quality AND full-primer interactions
  const { bestPair, pairScore, alternatives } = findBestPrimerPair(
    fwdCandidates,
    revCandidates,
    leftOverlap,
    rightOverlap,
    cfg
  );

  // Build full primers with overlap tails
  const fwdAnneal = bestPair.fwd;
  const revAnneal = bestPair.rev;

  const fwdWithTail = leftOverlap + fwdAnneal.sequence;
  const revWithTail = reverseComplement(rightOverlap) + revAnneal.sequence;

  // Calculate full primer properties
  const fwdFullTm = calculateTmQ5(fwdWithTail);
  const revFullTm = calculateTmQ5(revWithTail);
  const fwdFullGc = gcContent(fwdWithTail);
  const revFullGc = gcContent(revWithTail);

  // Calculate thermodynamic properties for tailed primers
  const fwdHairpin = calculateHairpinDG(fwdWithTail, 55);
  const revHairpin = calculateHairpinDG(revWithTail, 55);
  const heterodimerDG = calculateHeterodimerDG(fwdWithTail, revWithTail, 55);

  // Compile warnings
  const warnings: string[] = [];
  if (fwdHairpin < cfg.maxHairpinDG) {
    warnings.push(`Forward primer may form hairpin (ΔG = ${fwdHairpin.toFixed(1)} kcal/mol)`);
  }
  if (revHairpin < cfg.maxHairpinDG) {
    warnings.push(`Reverse primer may form hairpin (ΔG = ${revHairpin.toFixed(1)} kcal/mol)`);
  }
  if (heterodimerDG < cfg.maxHeterodimerDG) {
    warnings.push(`Primer pair may form heterodimer (ΔG = ${heterodimerDG.toFixed(1)} kcal/mol)`);
  }
  // Add warnings for annealing regions
  if (fwdAnneal.warnings) warnings.push(...fwdAnneal.warnings.map(w => `Forward: ${w}`));
  if (revAnneal.warnings) warnings.push(...revAnneal.warnings.map(w => `Reverse: ${w}`));

  // Calculate recommended annealing temperature
  const annealingTemp = Math.round(Math.min(fwdAnneal.tm, revAnneal.tm));

  return {
    forward: {
      sequence: fwdWithTail,
      length: fwdWithTail.length,
      annealingRegion: fwdAnneal.sequence,
      annealingLength: fwdAnneal.sequence.length,
      homologyTail: leftOverlap,
      tailLength: leftOverlap.length,
      tm: fwdAnneal.tm,
      tmFull: fwdFullTm,
      gc: Math.round(fwdFullGc * 1000) / 10,
      gcPercent: `${(fwdFullGc * 100).toFixed(1)}%`,
      hairpinDG: Math.round(fwdHairpin * 10) / 10,
      qualityScore: fwdAnneal.compositeScore,
      qualityTier: fwdAnneal.qualityTier,
      scoreBreakdown: fwdAnneal.scores,
    },
    reverse: {
      sequence: revWithTail,
      length: revWithTail.length,
      annealingRegion: revAnneal.sequence,
      annealingLength: revAnneal.sequence.length,
      homologyTail: reverseComplement(rightOverlap),
      tailLength: rightOverlap.length,
      tm: revAnneal.tm,
      tmFull: revFullTm,
      gc: Math.round(revFullGc * 1000) / 10,
      gcPercent: `${(revFullGc * 100).toFixed(1)}%`,
      hairpinDG: Math.round(revHairpin * 10) / 10,
      qualityScore: revAnneal.compositeScore,
      qualityTier: revAnneal.qualityTier,
      scoreBreakdown: revAnneal.scores,
    },
    pair: {
      annealingTemp,
      tmDifference: Math.round(Math.abs(fwdAnneal.tm - revAnneal.tm) * 10) / 10,
      heterodimerDG: Math.round(heterodimerDG * 10) / 10,
      pairScore,
      qualityTier: classifyQuality(pairScore).tier,
    },
    pcr: {
      annealingTemp,
      extensionTime: Math.ceil(fragmentSeq.length / 1000) * 30,
      cycles: 30,
    },
    productLength: fragmentSeq.length,
    warnings,
    method,
    alternatives: alternatives.slice(0, 3),
    optimization: {
      method: 'range-based',
      candidatesEvaluated: {
        forward: fwdCandidates.length,
        reverse: revCandidates.length,
        pairs: fwdCandidates.length * revCandidates.length,
      },
    },
  };
}

/**
 * Generate annealing region candidates with range-based scoring
 *
 * KEY FEATURES:
 * 1. POSITION FLEXIBILITY: Tries starting positions 0 to +5bp from fragment end
 * 2. LENGTH FLEXIBILITY: Tries all lengths within acceptable range
 * 3. SELF-COMPLEMENTARITY: Detects palindromes that interfere with PCR
 * 4. USES CENTRALIZED WEIGHTS: ANNEALING_SINGLE_WEIGHTS from weightCalibration.js
 */
function generateAnnealingCandidates(
  fragmentSeq: string,
  isFwd: boolean,
  ranges: { tm: { min: number; max: number }; gc: { min: number; max: number }; length: { min: number; max: number } },
  overlapTail: string
): AnnealingCandidate[] {
  const candidates: AnnealingCandidate[] = [];
  const seq = isFwd ? fragmentSeq : reverseComplement(fragmentSeq);

  // Position flexibility: how far to shift start position from fragment end
  const maxStartOffset = 5;

  // Try all lengths within acceptable range
  for (let len = ranges.length.min; len <= ranges.length.max; len++) {
    // POSITION FLEXIBILITY: Try different start positions
    // offset=0: standard position (starts at fragment end)
    // offset>0: shift inward (skip first N bases)
    for (let offset = 0; offset <= maxStartOffset; offset++) {
      if (offset + len > seq.length) continue;

      const annealSeq = seq.slice(offset, offset + len).toUpperCase();
      if (annealSeq.length < len) continue;

      const tm = calculateTmQ5(annealSeq);
      const gc = gcContent(annealSeq);
      const terminal3DG = calculate3primeTerminalDG(annealSeq)?.dG || -8;

      // Calculate secondary structure for annealing region alone
      let hairpinDG = 0;
      try {
        hairpinDG = foldDG(annealSeq, 55);
      } catch (e) {
        hairpinDG = 0;
      }

      // NEW: Self-complementarity scoring
      const selfCompScore = scoreSelfComplementarity(annealSeq);

      // RANGE-BASED SCORING: Score 1.0 within optimal range
      const tmScore = piecewiseLogistic(tm, {
        optimalLow: ranges.tm.min,
        optimalHigh: ranges.tm.max,
        acceptableLow: ranges.tm.min - 5,
        acceptableHigh: ranges.tm.max + 5,
        steepness: 0.5,
      });

      const gcScore = piecewiseLogistic(gc * 100, {
        optimalLow: ranges.gc.min * 100,
        optimalHigh: ranges.gc.max * 100,
        acceptableLow: 30,
        acceptableHigh: 70,
        steepness: 0.15,
      });

      const lengthScore = piecewiseLogistic(len, {
        optimalLow: ranges.length.min,
        optimalHigh: ranges.length.max,
        acceptableLow: ranges.length.min - 3,
        acceptableHigh: ranges.length.max + 5,
        steepness: 0.3,
      });

      // Build scores compatible with unified scoring system
      const homopolymerScore = scoreHomopolymer(annealSeq);
      const threePrimeComp = score3PrimeComposition(annealSeq, terminal3DG);
      const g4Analysis = analyzeGQuadruplex(annealSeq);

      const scores: Record<string, number> = {
        // Range-based primary scores
        tmInRange: tmScore,
        gcInRange: gcScore,
        lengthInRange: lengthScore,
        // Quality factors (unified scoring keys)
        gcClamp: scoreGcClamp(annealSeq),
        terminal3DG: scoreTerminal3DG(terminal3DG),
        hairpin: scoreHairpin(hairpinDG, { threshold: -3.0 }),
        homopolymer: homopolymerScore,
        threePrimeComp,
        gQuadruplex: g4Analysis.score,
        selfComplementarity: selfCompScore,
      };

      // Use centralized weights from weightCalibration.js
      const composite = calculateCompositeScore(scores, ANNEALING_SINGLE_WEIGHTS);
      const compositeScore = composite.score;

      // Check validity
      const isWithinTmRange = tm >= ranges.tm.min && tm <= ranges.tm.max;
      const isWithinGcRange = gc >= ranges.gc.min - 0.1 && gc <= ranges.gc.max + 0.1;
      const isValid = isWithinTmRange && isWithinGcRange;

      // Generate warnings
      const warnings: string[] = [];
      if (!isWithinTmRange) {
        warnings.push(`Tm ${tm.toFixed(1)}°C outside optimal range`);
      }
      if (scores.gcClamp < 0.7) {
        warnings.push('No GC clamp at 3\' end');
      }
      if (hairpinDG < -3) {
        warnings.push(`Hairpin risk (ΔG = ${hairpinDG.toFixed(1)} kcal/mol)`);
      }
      if (selfCompScore < 0.6) {
        warnings.push('Self-complementary sequence');
      }
      if (offset > 0) {
        warnings.push(`Start position shifted +${offset}bp`);
      }

      candidates.push({
        sequence: annealSeq,
        length: len,
        offset,  // NEW: Track start position offset
        tm: Math.round(tm * 10) / 10,
        gc: Math.round(gc * 1000) / 10,
        terminal3DG: Math.round(terminal3DG * 10) / 10,
        hairpinDG: Math.round(hairpinDG * 10) / 10,
        scores,
        compositeScore,
        qualityTier: compositeScore >= 75 ? 'excellent' :
                     compositeScore >= 60 ? 'good' :
                     compositeScore >= 45 ? 'acceptable' : 'marginal',
        isValid,
        warnings,
        isFwd,
      });
    }
  }

  // Sort by composite score
  candidates.sort((a, b) => b.compositeScore - a.compositeScore);

  return candidates;
}

/**
 * Find best primer pair considering pair-level interactions
 *
 * WEIGHT DISTRIBUTION (updated based on biological importance):
 * - 60% Individual annealing quality (annealing is rate-limiting for PCR)
 * - 40% Pair interactions:
 *   - 15% Tm difference (literature shows minimal effect)
 *   - 45% Heterodimer (critical in multi-fragment assembly)
 *   - 40% Full-primer hairpin (tailed primers often form stems)
 */
function findBestPrimerPair(
  fwdCandidates: AnnealingCandidate[],
  revCandidates: AnnealingCandidate[],
  leftOverlap: string,
  rightOverlap: string,
  cfg: AssemblyConfig
): {
  bestPair: { fwd: AnnealingCandidate; rev: AnnealingCandidate };
  pairScore: number;
  alternatives: Array<{
    forward: AnnealingCandidate;
    reverse: AnnealingCandidate;
    pairScore: number;
    tmDiff: number;
    heterodimerDG: number;
    selectionReason: string;
  }>;
} {
  const pairEvaluations: Array<{
    fwd: AnnealingCandidate;
    rev: AnnealingCandidate;
    pairScore: number;
    tmDiff: number;
    heterodimerDG: number;
    fwdFullHairpin: number;
    revFullHairpin: number;
    details: {
      individualScore: number;
      tmDiffScore: number;
      heterodimerScore: number;
      fullHairpinScore: number;
      weightBreakdown: {
        individual: number;
        interactions: number;
        tmDiffWeight: number;
        heterodimerWeight: number;
        hairpinWeight: number;
      };
    };
  }> = [];

  // Get top candidates to evaluate (limit to avoid combinatorial explosion)
  const topFwd = fwdCandidates.filter(c => c.isValid).slice(0, 10);
  const topRev = revCandidates.filter(c => c.isValid).slice(0, 10);

  // If no valid candidates, use best overall
  const fwdPool = topFwd.length > 0 ? topFwd : fwdCandidates.slice(0, 5);
  const revPool = topRev.length > 0 ? topRev : revCandidates.slice(0, 5);

  for (const fwd of fwdPool) {
    for (const rev of revPool) {
      // Build full primers
      const fwdFull = leftOverlap + fwd.sequence;
      const revFull = reverseComplement(rightOverlap) + rev.sequence;

      // Calculate pair-level properties
      const tmDiff = Math.abs(fwd.tm - rev.tm);
      const heterodimerDG = calculateHeterodimerDG(fwdFull, revFull, 55);
      const fwdFullHairpin = calculateHairpinDG(fwdFull, 55);
      const revFullHairpin = calculateHairpinDG(revFull, 55);

      // Score Tm difference (generous - literature shows minimal effect)
      const tmDiffScore = tmDiff <= 3 ? 1.0 :
                          tmDiff <= 5 ? 0.9 :
                          tmDiff <= 8 ? 0.7 : 0.4;

      // Score heterodimer (stricter - critical for multi-fragment assembly)
      const heterodimerScore = heterodimerDG >= -3 ? 1.0 :
                               heterodimerDG >= -5 ? 0.85 :
                               heterodimerDG >= -7 ? 0.6 :
                               heterodimerDG >= -9 ? 0.35 : 0.15;

      // Score full-primer hairpins (use relaxed threshold for tailed primers)
      const fullHairpinScore = Math.min(
        scoreHairpin(fwdFullHairpin, { threshold: -5.0 }),
        scoreHairpin(revFullHairpin, { threshold: -5.0 })
      );

      // Calculate pair composite score
      // ADJUSTED WEIGHTS: 60% individual + 40% interactions
      // This prioritizes annealing quality (rate-limiting step)
      const individualScore = (fwd.compositeScore + rev.compositeScore) / 2;

      // Pair interaction weights (within the 40%):
      // - Tm difference: 15% (literature shows minimal effect)
      // - Heterodimer: 45% (critical in multi-fragment)
      // - Full hairpin: 40% (tailed primers often form stems)
      const pairInteractionScore = (tmDiffScore * 15 + heterodimerScore * 45 + fullHairpinScore * 40);
      const pairScore = Math.round(individualScore * 0.6 + pairInteractionScore * 0.4);

      pairEvaluations.push({
        fwd,
        rev,
        pairScore,
        tmDiff,
        heterodimerDG,
        fwdFullHairpin,
        revFullHairpin,
        details: {
          individualScore,
          tmDiffScore,
          heterodimerScore,
          fullHairpinScore,
          weightBreakdown: {
            individual: 0.6,
            interactions: 0.4,
            tmDiffWeight: 0.15,
            heterodimerWeight: 0.45,
            hairpinWeight: 0.40,
          },
        },
      });
    }
  }

  // Sort by pair score
  pairEvaluations.sort((a, b) => b.pairScore - a.pairScore);

  const best = pairEvaluations[0];
  const alternatives = pairEvaluations.slice(1, 5).map(p => ({
    forward: p.fwd,
    reverse: p.rev,
    pairScore: p.pairScore,
    tmDiff: p.tmDiff,
    heterodimerDG: p.heterodimerDG,
    selectionReason: p.heterodimerDG > best.heterodimerDG ? 'Lower dimer risk' :
                     p.tmDiff < best.tmDiff ? 'Better Tm matching' :
                     Math.min(p.fwdFullHairpin, p.revFullHairpin) > Math.min(best.fwdFullHairpin, best.revFullHairpin) ? 'Better hairpin avoidance' :
                     'Alternative design',
  }));

  return {
    bestPair: { fwd: best.fwd, rev: best.rev },
    pairScore: best.pairScore,
    alternatives,
  };
}

/**
 * Design primers for a complete assembly
 */
export function designAssembly(fragments: Fragment[], options: Partial<AssemblyConfig> = {}): AssemblyDesign {
  const cfg: AssemblyConfig = { ...DEFAULT_ASSEMBLY_CONFIG, ...options };
  const method = ASSEMBLY_METHODS[cfg.method] || ASSEMBLY_METHODS.NEBUILDER_HIFI;

  if (fragments.length < 2) {
    throw new Error('Assembly requires at least 2 fragments');
  }

  if (fragments.length > method.maxFragments) {
    throw new Error(`${method.name} supports maximum ${method.maxFragments} fragments`);
  }

  // First, optimize all overlaps
  // When autoOptimize is enabled, uses range-based scoring
  const overlapPlan = optimizeAssemblyOverlaps(fragments, cfg);

  // Design primers for each fragment
  // When autoOptimize is enabled, uses range-based primer optimization
  const fragmentDesigns: AssemblyDesign['fragments'] = [];
  const useOptimizedPrimers = cfg.autoOptimize === true;

  for (let i = 0; i < fragments.length; i++) {
    const frag = fragments[i];
    const prevJunction = overlapPlan.junctions[(i - 1 + fragments.length) % fragments.length];
    const nextJunction = overlapPlan.junctions[i];

    // Get overlap sequences
    // Left overlap: comes from previous fragment's end
    // Right overlap: goes to next fragment's start
    const leftOverlap = prevJunction.overlap.sequence;
    const rightOverlap = nextJunction.overlap.sequence;

    // Design primers using either optimized (range-based) or standard (target-based) approach
    const primerDesignFn = useOptimizedPrimers ? designAssemblyPrimersOptimized : designAssemblyPrimers;
    const primerDesign = primerDesignFn(frag.seq, {
      leftOverlap: i === 0 && !cfg.circular ? '' : leftOverlap,
      rightOverlap,
      method: cfg.method,
    }, cfg);

    fragmentDesigns.push({
      id: frag.id,
      index: i + 1,
      originalLength: frag.seq.length,
      primers: primerDesign,
      leftOverlap: {
        sequence: leftOverlap,
        from: fragments[(i - 1 + fragments.length) % fragments.length].id,
      },
      rightOverlap: {
        sequence: rightOverlap,
        to: fragments[(i + 1) % fragments.length].id,
      },
    });
  }

  // Calculate total primer cost estimate
  const totalPrimerBases = fragmentDesigns.reduce((sum, f) =>
    sum + f.primers.forward.length + f.primers.reverse.length, 0
  );
  const estimatedPrimerCost = totalPrimerBases * 0.10; // ~$0.10 per base

  // Collect all warnings
  const allWarnings = [
    ...overlapPlan.warnings,
    ...fragmentDesigns.flatMap((f, i) =>
      f.primers.warnings.map(w => `Fragment ${i + 1} (${f.id}): ${w}`)
    ),
  ];

  return {
    method: method.name,
    methodKey: cfg.method,
    fragments: fragmentDesigns,
    junctions: overlapPlan.junctions,
    quality: overlapPlan.quality,
    assembly: {
      circular: cfg.circular !== false,
      totalLength: fragments.reduce((sum, f) => sum + f.seq.length, 0),
      fragmentCount: fragments.length,
    },
    protocol: generateProtocol(method, fragmentDesigns),
    cost: {
      primers: Math.round(estimatedPrimerCost * 100) / 100,
      assembly: method.name === 'NEBuilder HiFi DNA Assembly' ? 12.98 : 15.00,
      total: Math.round((estimatedPrimerCost + 12.98) * 100) / 100,
    },
    warnings: allWarnings,
    // Optimization metadata
    optimization: {
      method: useOptimizedPrimers ? 'range-based' : 'target-based',
      autoOptimize: cfg.autoOptimize || false,
      description: useOptimizedPrimers
        ? 'Range-based optimization: equal scores within optimal ranges, differentiated by quality factors'
        : 'Target-based optimization: scores penalized for deviation from target Tm/length values',
    },
  };
}

// =============================================================================
// Protocol Generation
// =============================================================================

/**
 * Generate assembly protocol
 */
function generateProtocol(method: AssemblyMethod, fragmentDesigns: AssemblyDesign['fragments']): Protocol {
  const isNEBuilder = method.name === 'NEBuilder HiFi DNA Assembly';
  const isGibson = method.name === 'Gibson Assembly';

  const steps = [
    {
      step: 1,
      title: 'Prepare DNA fragments',
      details: [
        'PCR amplify each fragment using the designed primers',
        `Use high-fidelity polymerase (Q5 or Phusion recommended)`,
        'Run on gel to verify correct size',
        'Purify PCR products (column or gel extraction)',
        'Quantify DNA concentration (Nanodrop or Qubit)',
      ],
    },
    {
      step: 2,
      title: 'Setup assembly reaction',
      details: isNEBuilder ? [
        'Prepare equimolar amounts of fragments (0.03-0.2 pmol each)',
        'For 2-3 fragments: 0.03-0.2 pmol total',
        'For 4+ fragments: 0.2-0.5 pmol total',
        'Add 10 µL NEBuilder HiFi DNA Assembly Master Mix',
        'Add water to 20 µL total volume',
      ] : [
        'Prepare equimolar amounts of fragments',
        'Mix fragments (50-100 ng each)',
        'Add 15 µL Gibson Assembly Master Mix',
        'Add water to 20 µL total volume',
      ],
    },
    {
      step: 3,
      title: 'Incubate',
      details: isNEBuilder ? [
        `Incubate at ${method.incubationTemp}°C for ${fragmentDesigns.length <= 3 ? 15 : 60} minutes`,
        '15 minutes for 2-3 fragments',
        '60 minutes for 4+ fragments',
        'Can be left up to 60 minutes without issue',
      ] : [
        `Incubate at ${method.incubationTemp}°C for ${method.incubationTime} minutes`,
      ],
    },
    {
      step: 4,
      title: 'Transform',
      details: [
        'Transform 2 µL of reaction into competent cells',
        'Use NEB 10-beta or similar high-efficiency cells',
        'Follow standard transformation protocol',
        'Plate on appropriate antibiotic selection',
      ],
    },
  ];

  return {
    title: `${method.name} Protocol`,
    steps,
    tips: [
      'Overlapping regions should have 15-25 bp homology',
      `Optimal overlap Tm is ${method.overlapTm?.optimal}°C`,
      'DpnI digest template if using plasmid-derived fragments',
      'Minimize freeze-thaw cycles of master mix',
    ],
  };
}

// =============================================================================
// Primer Export
// =============================================================================

/**
 * Export primers in various formats
 */
export function exportPrimers(assemblyDesign: AssemblyDesign, format: string = 'tsv'): string | Array<Record<string, any>> {
  const rows: Array<{ name: string; sequence: string; length: number; tm: number; notes: string }> = [];

  for (const frag of assemblyDesign.fragments) {
    rows.push({
      name: `${frag.id}_F`,
      sequence: frag.primers.forward.sequence,
      length: frag.primers.forward.length,
      tm: frag.primers.forward.tm,
      notes: `Forward primer for ${frag.id}`,
    });
    rows.push({
      name: `${frag.id}_R`,
      sequence: frag.primers.reverse.sequence,
      length: frag.primers.reverse.length,
      tm: frag.primers.reverse.tm,
      notes: `Reverse primer for ${frag.id}`,
    });
  }

  if (format === 'tsv') {
    const header = 'Name\tSequence\tLength\tTm\tNotes';
    const lines = rows.map(r => `${r.name}\t${r.sequence}\t${r.length}\t${r.tm}\t${r.notes}`);
    return [header, ...lines].join('\n');
  }

  if (format === 'csv') {
    const header = 'Name,Sequence,Length,Tm,Notes';
    const lines = rows.map(r => `${r.name},${r.sequence},${r.length},${r.tm},"${r.notes}"`);
    return [header, ...lines].join('\n');
  }

  if (format === 'json') {
    return JSON.stringify(rows, null, 2);
  }

  return rows;
}

// =============================================================================
// Assembly Simulation
// =============================================================================

/**
 * Simulate the assembly outcome - build the assembled sequence
 * and verify the assembly is valid
 */
export function simulateAssembly(
  fragments: Fragment[],
  junctions: Junction[],
  options: { circular?: boolean; method?: string } = {}
): SimulationResult {
  const { circular = true, method = 'NEBUILDER_HIFI' } = options;

  if (fragments.length < 2) {
    throw new Error('Assembly requires at least 2 fragments');
  }

  // Build the assembled sequence
  let assembledSeq = '';
  const assemblySteps: SimulationResult['steps'] = [];
  let currentPosition = 0;

  for (let i = 0; i < fragments.length; i++) {
    const frag = fragments[i];
    const junction = junctions[i];

    // For the first fragment, add the full sequence
    if (i === 0) {
      assembledSeq = frag.seq;
      assemblySteps.push({
        step: 1,
        fragment: frag.id,
        action: 'Initial fragment',
        position: 0,
        length: frag.seq.length,
        cumulativeLength: frag.seq.length,
      });
      currentPosition = frag.seq.length;
    } else {
      // For subsequent fragments, add sequence after the overlap region
      const overlapLen = junction?.overlap?.length || 20;

      // The overlap is at the start of this fragment (shared with end of previous)
      const fragWithoutOverlap = frag.seq.slice(overlapLen);

      assembledSeq += fragWithoutOverlap;

      assemblySteps.push({
        step: i + 1,
        fragment: frag.id,
        action: 'Joined via overlap',
        overlapSequence: frag.seq.slice(0, overlapLen),
        overlapLength: overlapLen,
        addedLength: fragWithoutOverlap.length,
        position: currentPosition,
        cumulativeLength: assembledSeq.length,
      });
      currentPosition = assembledSeq.length;
    }
  }

  // For circular assemblies, verify the last→first junction closes properly
  let circularClosureValid = false;
  let closureAnalysis: SimulationResult['circularClosure'] = null;

  if (circular) {
    const lastFrag = fragments[fragments.length - 1];
    const firstFrag = fragments[0];
    const lastJunction = junctions[junctions.length - 1];
    const overlapLen = lastJunction?.overlap?.length || 20;

    // The end of the assembled sequence should match the start
    const endSeq = assembledSeq.slice(-overlapLen);
    const startSeq = assembledSeq.slice(0, overlapLen);

    // Check if they would circularize correctly
    const lastFragEnd = lastFrag.seq.slice(-overlapLen);
    const firstFragStart = firstFrag.seq.slice(0, overlapLen);

    circularClosureValid = lastFragEnd === firstFragStart;

    closureAnalysis = {
      valid: circularClosureValid,
      expectedOverlap: firstFragStart,
      actualOverlap: lastFragEnd,
      matches: circularClosureValid,
      message: circularClosureValid
        ? 'Assembly circularizes correctly'
        : `Circular closure mismatch: ${lastFrag.id} end doesn't match ${firstFrag.id} start`,
    };

    // If circular, trim the redundant overlap at the end
    if (circularClosureValid) {
      // The assembled sequence already has the closure built in
      // For display purposes, we keep the linear version
    }
  }

  // Detect potential homology conflicts (internal repeats that could cause misassembly)
  const homologyConflicts = detectHomologyConflicts(fragments, junctions);

  // Calculate assembly properties
  const assemblyLength = assembledSeq.length;
  const gc = gcContent(assembledSeq);

  // Generate visualization data
  const visualization = generateVisualizationData(fragments, junctions, {
    circular,
    assembledLength: assemblyLength,
  });

  return {
    success: circular ? circularClosureValid : true,
    method,
    assembledSequence: assembledSeq,
    assembledLength: assemblyLength,
    gcContent: Math.round(gc * 1000) / 10,
    gcPercent: `${(gc * 100).toFixed(1)}%`,
    circular,
    circularClosure: closureAnalysis,
    steps: assemblySteps,
    homologyConflicts,
    hasConflicts: homologyConflicts.length > 0,
    visualization,
    fragments: fragments.map((f, i) => ({
      id: f.id,
      length: f.seq.length,
      position: assemblySteps[i]?.position || 0,
      gcContent: gcContent(f.seq) * 100,
    })),
  };
}

/**
 * Detect potential homology conflicts that could cause misassembly
 * Looks for repeated sequences between non-adjacent fragments
 */
function detectHomologyConflicts(fragments: Fragment[], junctions: Junction[]): SimulationResult['homologyConflicts'] {
  const conflicts: SimulationResult['homologyConflicts'] = [];
  const minHomologyLen = 15; // Minimum homology that could cause recombination

  for (let i = 0; i < fragments.length; i++) {
    for (let j = i + 2; j < fragments.length; j++) {
      // Skip adjacent fragments (they're supposed to have homology)
      if (j === (i + 1) % fragments.length) continue;
      if (i === 0 && j === fragments.length - 1) continue; // Circular closure

      const frag1 = fragments[i];
      const frag2 = fragments[j];

      // Check for significant shared homology
      const homology = findSharedHomology(frag1.seq, frag2.seq, minHomologyLen);

      if (homology) {
        conflicts.push({
          fragment1: frag1.id,
          fragment1Index: i,
          fragment2: frag2.id,
          fragment2Index: j,
          homologySequence: homology.sequence,
          homologyLength: homology.length,
          position1: homology.pos1,
          position2: homology.pos2,
          severity: homology.length >= 25 ? 'high' : 'medium',
          warning: `Shared ${homology.length}bp homology between ${frag1.id} and ${frag2.id} may cause misassembly`,
        });
      }
    }
  }

  return conflicts;
}

/**
 * Find the longest shared subsequence between two sequences
 */
function findSharedHomology(
  seq1: string,
  seq2: string,
  minLength: number
): { sequence: string; length: number; pos1: number; pos2: number; isReverseComplement?: boolean } | null {
  const s1 = seq1.toUpperCase();
  const s2 = seq2.toUpperCase();

  let bestMatch: { sequence: string; length: number; pos1: number; pos2: number; isReverseComplement?: boolean } | null = null;

  // Sliding window search for shared subsequences
  for (let len = Math.min(50, Math.min(s1.length, s2.length)); len >= minLength; len--) {
    for (let i = 0; i <= s1.length - len; i++) {
      const subseq = s1.slice(i, i + len);
      const pos2 = s2.indexOf(subseq);

      if (pos2 !== -1) {
        return {
          sequence: subseq,
          length: len,
          pos1: i,
          pos2: pos2,
        };
      }

      // Also check reverse complement
      const subseqRC = reverseComplement(subseq);
      const pos2RC = s2.indexOf(subseqRC);

      if (pos2RC !== -1) {
        return {
          sequence: subseq,
          length: len,
          pos1: i,
          pos2: pos2RC,
          isReverseComplement: true,
        };
      }
    }
  }

  return null;
}

/**
 * Generate visualization data for the assembly
 */
function generateVisualizationData(
  fragments: Fragment[],
  junctions: Junction[],
  options: { circular?: boolean; assembledLength?: number } = {}
): VisualizationData {
  const { circular = true, assembledLength = 0 } = options;

  // Calculate positions for each fragment in the assembled sequence
  const fragmentPositions: VisualizationData['fragments'] = [];
  let currentPos = 0;

  for (let i = 0; i < fragments.length; i++) {
    const frag = fragments[i];
    const junction = junctions[i];
    const overlapLen = junction?.overlap?.length || 20;

    const fragLength = i === 0 ? frag.seq.length : frag.seq.length - overlapLen;
    const startPos = i === 0 ? 0 : currentPos;

    fragmentPositions.push({
      id: frag.id,
      index: i,
      startPosition: startPos,
      endPosition: startPos + frag.seq.length,
      displayLength: fragLength,
      actualLength: frag.seq.length,
      overlapWithNext: overlapLen,
      color: getFragmentColor(i),
      // For circular diagram
      startAngle: circular ? (startPos / assembledLength) * 360 : null,
      endAngle: circular ? ((startPos + fragLength) / assembledLength) * 360 : null,
    });

    currentPos = startPos + fragLength;
  }

  // Generate junction visualization data
  const junctionPositions: VisualizationData['junctions'] = junctions.map((j, i) => {
    const frag = fragmentPositions[i];
    return {
      index: i,
      from: j.from || fragments[i]?.id,
      to: j.to || fragments[(i + 1) % fragments.length]?.id,
      position: frag.endPosition - (j.overlap?.length || 20),
      overlapSequence: j.overlap?.sequence,
      quality: j.overlap?.qualityTier || 'unknown',
      score: j.overlap?.score || j.overlap?.compositeScore,
    };
  });

  return {
    type: circular ? 'circular' : 'linear',
    totalLength: assembledLength,
    fragments: fragmentPositions,
    junctions: junctionPositions,
    // SVG rendering hints
    svg: {
      width: 400,
      height: circular ? 400 : 100,
      radius: circular ? 150 : null,
      centerX: circular ? 200 : null,
      centerY: circular ? 200 : null,
    },
  };
}

/**
 * Get a color for a fragment based on its index
 */
function getFragmentColor(index: number): string {
  const colors = [
    '#3b82f6', // blue
    '#10b981', // green
    '#f59e0b', // amber
    '#ef4444', // red
    '#8b5cf6', // purple
    '#06b6d4', // cyan
    '#ec4899', // pink
    '#84cc16', // lime
  ];
  return colors[index % colors.length];
}

/**
 * Verify junction compatibility - check if overlaps actually exist
 * in the provided sequences
 */
export function validateJunctions(
  fragments: Fragment[],
  options: { circular?: boolean; minOverlap?: number } = {}
): JunctionValidation {
  const { circular = true, minOverlap = 15 } = options;

  const validations: JunctionValidation['junctions'] = [];

  for (let i = 0; i < fragments.length; i++) {
    const frag1 = fragments[i];
    const frag2 = fragments[(i + 1) % fragments.length];

    // Skip last junction for linear assemblies
    if (!circular && i === fragments.length - 1) continue;

    const frag1End = frag1.seq.slice(-50).toUpperCase();
    const frag2Start = frag2.seq.slice(0, 50).toUpperCase();

    // Check if there's natural homology between the fragments
    let naturalHomology: { sequence: string; length: number; hasNaturalHomology: boolean } | null = null;
    for (let len = 35; len >= minOverlap; len--) {
      const frag1Tail = frag1End.slice(-len);
      if (frag2Start.startsWith(frag1Tail)) {
        naturalHomology = {
          sequence: frag1Tail,
          length: len,
          hasNaturalHomology: true,
        };
        break;
      }
    }

    const validation = {
      junction: i + 1,
      from: frag1.id,
      to: frag2.id,
      hasNaturalHomology: !!naturalHomology,
      naturalHomology,
      requiresDesignedTails: !naturalHomology,
      message: naturalHomology
        ? `Natural ${naturalHomology.length}bp homology found`
        : 'No natural homology - primers will add homology tails',
    };

    validations.push(validation);
  }

  const hasAllNaturalHomology = validations.every(v => v.hasNaturalHomology);
  const requiresAnyTails = validations.some(v => v.requiresDesignedTails);

  return {
    valid: true, // Assembly is always possible with designed primers
    junctions: validations,
    hasAllNaturalHomology,
    requiresDesignedTails: requiresAnyTails,
    summary: hasAllNaturalHomology
      ? 'All junctions have natural homology - simple overlap assembly possible'
      : 'Primers will add homology tails for seamless junction formation',
  };
}

// =============================================================================
// GenBank Export
// =============================================================================

/**
 * Export assembled sequence in GenBank format
 */
export function exportToGenBank(
  simulation: SimulationResult,
  options: {
    name?: string;
    description?: string;
    organism?: string;
    molType?: string;
    topology?: string;
  } = {}
): string {
  const {
    name = 'Assembly',
    description = 'DNA assembly designed with Primer Designer',
    organism = 'synthetic construct',
    molType = 'DNA',
    topology = simulation.circular ? 'circular' : 'linear',
  } = options;

  const seq = simulation.assembledSequence;
  const length = simulation.assembledLength;
  const date = new Date().toISOString().split('T')[0].replace(/-/g, '-');
  const dateFormatted = formatGenBankDate(new Date());

  let gb = '';

  // LOCUS line (standardized format)
  const locusName = name.replace(/\s+/g, '_').slice(0, 16).padEnd(16, ' ');
  gb += `LOCUS       ${locusName} ${length.toString().padStart(7, ' ')} bp    ${molType.padEnd(6, ' ')} ${topology.padEnd(8, ' ')} ${dateFormatted}\n`;

  // DEFINITION
  gb += `DEFINITION  ${description}\n`;

  // ACCESSION (placeholder)
  gb += `ACCESSION   .\n`;

  // VERSION
  gb += `VERSION     .\n`;

  // KEYWORDS
  gb += `KEYWORDS    .\n`;

  // SOURCE
  gb += `SOURCE      ${organism}\n`;
  gb += `  ORGANISM  ${organism}\n`;
  gb += `            .\n`;

  // COMMENT
  gb += `COMMENT     Assembly method: ${simulation.method || 'Unknown'}\n`;
  gb += `            Fragments: ${simulation.fragments?.length || 'Unknown'}\n`;
  gb += `            GC content: ${simulation.gcPercent || 'Unknown'}\n`;
  gb += `            Designed with Primer Designer Tool\n`;

  // FEATURES
  gb += `FEATURES             Location/Qualifiers\n`;

  // Add source feature
  gb += `     source          1..${length}\n`;
  gb += `                     /organism="${organism}"\n`;
  gb += `                     /mol_type="other DNA"\n`;

  // Add fragment features from visualization data
  if (simulation.visualization?.fragments) {
    for (const frag of simulation.visualization.fragments) {
      const start = frag.startPosition + 1; // GenBank is 1-based
      const end = frag.endPosition;

      gb += `     misc_feature    ${start}..${end}\n`;
      gb += `                     /label="${frag.id}"\n`;
      gb += `                     /note="Fragment ${frag.index + 1}: ${frag.actualLength} bp"\n`;
      gb += `                     /ApEinfo_fwdcolor="${frag.color}"\n`;
    }
  }

  // Add junction features
  if (simulation.visualization?.junctions) {
    for (const junc of simulation.visualization.junctions) {
      const pos = junc.position + 1;
      const endPos = pos + (junc.overlapSequence?.length || 20) - 1;

      gb += `     misc_feature    ${pos}..${endPos}\n`;
      gb += `                     /label="Junction ${junc.index + 1}"\n`;
      gb += `                     /note="${junc.from} -> ${junc.to}"\n`;
      if (junc.overlapSequence) {
        gb += `                     /note="Overlap: ${junc.overlapSequence}"\n`;
      }
      gb += `                     /ApEinfo_fwdcolor="#808080"\n`;
    }
  }

  // ORIGIN (sequence)
  gb += `ORIGIN\n`;

  // Format sequence in GenBank style (60 chars per line, numbered)
  const seqLower = seq.toLowerCase();
  for (let i = 0; i < seqLower.length; i += 60) {
    const lineNum = (i + 1).toString().padStart(9, ' ');
    const lineSeq = seqLower.slice(i, i + 60);

    // Split into groups of 10
    const groups: string[] = [];
    for (let j = 0; j < lineSeq.length; j += 10) {
      groups.push(lineSeq.slice(j, j + 10));
    }

    gb += `${lineNum} ${groups.join(' ')}\n`;
  }

  gb += `//\n`;

  return gb;
}

/**
 * Format date for GenBank LOCUS line
 */
function formatGenBankDate(date: Date): string {
  const months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                  'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'];
  const day = date.getDate().toString().padStart(2, '0');
  const month = months[date.getMonth()];
  const year = date.getFullYear();
  return `${day}-${month}-${year}`;
}

/**
 * Export assembly project to JSON (for save/load)
 */
export function exportProject(projectData: Record<string, any>): string {
  const project = {
    version: '1.0',
    exportedAt: new Date().toISOString(),
    tool: 'Primer Designer - Assembly Tool',
    ...projectData,
  };
  return JSON.stringify(project, null, 2);
}

/**
 * Import assembly project from JSON
 */
export function importProject(jsonString: string): Record<string, any> {
  const project = JSON.parse(jsonString);

  // Validate project structure
  if (!project.version) {
    throw new Error('Invalid project file: missing version');
  }

  if (!project.fragments || !Array.isArray(project.fragments)) {
    throw new Error('Invalid project file: missing fragments array');
  }

  // Validate each fragment has required fields
  for (const frag of project.fragments) {
    if (!frag.id || !frag.seq) {
      throw new Error('Invalid project file: fragments must have id and seq');
    }
  }

  return project;
}

/**
 * Export to FASTA format
 */
export function exportToFasta(
  simulation: SimulationResult,
  options: {
    name?: string;
    description?: string;
    lineWidth?: number;
  } = {}
): string {
  const {
    name = 'Assembly',
    description = '',
    lineWidth = 70,
  } = options;

  let fasta = `>${name}`;
  if (description) {
    fasta += ` ${description}`;
  }
  fasta += '\n';

  // Wrap sequence to specified line width
  const seq = simulation.assembledSequence;
  for (let i = 0; i < seq.length; i += lineWidth) {
    fasta += seq.slice(i, i + lineWidth) + '\n';
  }

  return fasta;
}

/**
 * Export to SnapGene DNA format (simplified)
 * Note: Full SnapGene format is proprietary; this creates a compatible subset
 */
export function exportToSnapGene(
  simulation: SimulationResult,
  options: Record<string, any> = {}
): { format: string; content: string; note: string } {
  // SnapGene can import GenBank format
  // Return GenBank with SnapGene-compatible feature colors
  return {
    format: 'genbank',
    content: exportToGenBank(simulation, {
      ...options,
      // Add SnapGene-specific color annotations
    }),
    note: 'SnapGene can import this GenBank file directly',
  };
}
