/**
 * Primer Boundary Optimizer
 *
 * Optimizes primer design at fragment junctions by adjusting the boundary
 * position between user-defined fragments. This is different from the
 * Fusion Site Optimizer which finds optimal positions from scratch.
 *
 * Use case:
 * - User provides multiple fragments with pre-defined boundaries
 * - At each junction, primers may have poor quality (low Tm, hairpins, etc.)
 * - This optimizer SHIFTS the boundary left/right to get better primer binding regions
 * - The total assembled sequence remains the same
 *
 * Example:
 *   Before (user-defined boundary, poor primer):
 *   Fragment A: ...NNNNN[AAAA]  |junction|  [GCGC]NNNNN... Fragment B
 *                    ↑ poor Tm (poly-A)    ↑ good binding
 *
 *   After (shifted boundary, better primer):
 *   Fragment A: ...NNNNN  |shifted|  [AGCT]NNNNNNNN... Fragment B
 *                           ←─5bp─→
 *                         (boundary shifted LEFT)
 *
 *   Result: Fragment A is 5bp shorter, Fragment B is 5bp longer,
 *           but primer for B now has better binding region.
 *
 * Key difference from Fusion Site Optimizer:
 * - Fusion Site Optimizer: finds optimal positions in a single sequence
 * - Primer Boundary Optimizer: adjusts user-defined boundaries for better primers
 */

import { calculateTmQ5 } from '../tmQ5.js';
import { calculateHairpinDG, calculateHomodimerDG } from '../equilibrium.js';
import { reverseComplement } from './enzymes.js';
import {
  scoreTm,
  scoreGc,
  scoreHairpin,
  scoreHomodimer,
  scoreGcClamp,
  scoreHomopolymer,
  classifyQuality,
} from '../scoring.js';
import { getOverhangFidelityExperimental } from './goldengate.js';

// ============================================================================
// TYPES
// ============================================================================

interface BoundaryOptimizerConfig {
  maxShift: number;
  minShiftThreshold: number;
  minHomologyLength: number;
  maxHomologyLength: number;
  targetTm: number;
  tmTolerance: number;
  minFragmentSize: number;
  weights: {
    primerQuality: number;
    tmMatch: number;
    balancedShift: number;
    overhangFidelity: number;
  };
  avoidPalindromes: boolean;
  avoidHomopolymers: boolean;
  minOverhangFidelity: number;
}

interface Fragment {
  id?: string;
  seq?: string;
  sequence?: string;
}

interface PrimerQualityScore {
  score: number;
  quality: string;
  tm: number;
  gc: number;
  hairpinDG: number;
  homodimerDG: number;
  gcClampScore: number;
  homopolymerScore: number;
  issues: string[];
  breakdown: {
    tm: number;
    gc: number;
    hairpin: number;
    homodimer: number;
    gcClamp: number;
    homopolymer: number;
  };
}

interface OverhangValidation {
  valid: boolean;
  issues: string[];
}

interface BoundaryOptimizationResult {
  success: boolean;
  originalPosition: number;
  optimizedPosition: number;
  shift: number;
  direction: 'left' | 'right' | 'unchanged';
  beforePrimers: {
    left: PrimerQualityScore;
    right: PrimerQualityScore;
    composite: number;
    overhang: string;
    overhangValid: boolean;
  };
  afterPrimers: {
    left: PrimerQualityScore;
    right: PrimerQualityScore;
    composite: number;
    overhang: string;
    overhangValid: boolean;
  };
  improvement: number;
  reason: string;
}

interface JunctionResult extends BoundaryOptimizationResult {
  junctionIndex: number;
  leftFragmentId: string;
  rightFragmentId: string;
  error?: string;
}

interface OptimizedFragment {
  id: string;
  seq: string;
  originalSeq: string;
  originalLength: number;
  newLength: number;
  lengthChange: number;
  error?: string;
}

interface AssemblyOptimizationResult {
  success: boolean;
  boundaries: JunctionResult[];
  optimizedFragments: OptimizedFragment[];
  summary?: {
    totalBoundaries: number;
    boundariesOptimized: number;
    boundariesUnchanged: number;
    averageImprovement: number;
    totalImprovement: number;
  };
  error?: string;
}

interface AssessmentIssue {
  junction: number;
  side: 'left' | 'right';
  fragmentId: string;
  quality: string;
  score: number;
  issues: string[];
}

interface BoundaryAssessmentResult {
  needsOptimization: boolean;
  junctionCount: number;
  averageScore: number;
  overallQuality: string;
  issues: AssessmentIssue[];
  recommendation: string;
}

interface ScoringOptions {
  targetTm?: number;
  tmTolerance?: number;
}

interface OptimizationOptions extends Partial<BoundaryOptimizerConfig> {}

// ============================================================================
// CONFIGURATION
// ============================================================================

export const BOUNDARY_OPTIMIZER_DEFAULTS: BoundaryOptimizerConfig = {
  // How far to search from the original boundary (in bp)
  maxShift: 50,

  // Minimum shift to consider (ignore tiny improvements)
  minShiftThreshold: 3,

  // Primer binding region constraints
  minHomologyLength: 15,
  maxHomologyLength: 30,
  targetTm: 60,           // °C
  tmTolerance: 5,         // ±5°C acceptable

  // Minimum fragment size after shifting
  minFragmentSize: 100,

  // Weight for different optimization criteria
  weights: {
    primerQuality: 0.40,    // Overall primer score
    tmMatch: 0.25,          // How close to target Tm
    balancedShift: 0.15,    // Prefer smaller shifts
    overhangFidelity: 0.20, // Keep good overhangs
  },

  // Overhang constraints
  avoidPalindromes: true,
  avoidHomopolymers: true,
  minOverhangFidelity: 0.90,
};

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Calculate primer quality score for a binding region
 */
function scorePrimerBindingRegion(homology: string, options: ScoringOptions = {}): PrimerQualityScore {
  const { targetTm = 60, tmTolerance = 5 } = options;

  if (!homology || homology.length < 10) {
    return {
      score: 0,
      quality: 'poor',
      tm: 0,
      gc: 0,
      hairpinDG: 0,
      homodimerDG: 0,
      gcClampScore: 0,
      homopolymerScore: 0,
      issues: ['Homology too short'],
      breakdown: {
        tm: 0,
        gc: 0,
        hairpin: 0,
        homodimer: 0,
        gcClamp: 0,
        homopolymer: 0,
      },
    };
  }

  // Calculate thermodynamic properties
  const tm = calculateTmQ5(homology);
  const gc = (homology.match(/[GC]/gi) || []).length / homology.length;
  const hairpinDG = calculateHairpinDG(homology);
  const homodimerDG = calculateHomodimerDG(homology);

  // Score each property
  const tmScore = scoreTm(tm, { targetTm } as any);  // FIXED: Type assertion for TmScoringOptions
  const gcScore = scoreGc(gc);
  const hairpinScore = scoreHairpin(hairpinDG);
  const homodimerScore = scoreHomodimer(homodimerDG);
  const gcClampScore = scoreGcClamp(homology);
  const homopolymerScore = scoreHomopolymer(homology);

  // Calculate composite score
  const composite = (
    tmScore * 0.25 +
    gcScore * 0.15 +
    hairpinScore * 0.20 +
    homodimerScore * 0.15 +
    gcClampScore * 0.15 +
    homopolymerScore * 0.10
  );

  // Identify issues
  const issues: string[] = [];
  if (Math.abs(tm - targetTm) > tmTolerance) {
    issues.push(`Tm ${tm.toFixed(1)}°C outside target range`);
  }
  if (gc < 0.35 || gc > 0.65) {
    issues.push(`GC content ${(gc * 100).toFixed(1)}% outside optimal range`);
  }
  if (hairpinDG < -3.0) {
    issues.push(`Strong hairpin (ΔG = ${hairpinDG.toFixed(1)} kcal/mol)`);
  }
  if (homodimerDG < -6.0) {
    issues.push(`Strong homodimer (ΔG = ${homodimerDG.toFixed(1)} kcal/mol)`);
  }

  // Check for homopolymers
  const homopolymerMatch = homology.match(/(.)\1{3,}/g);
  if (homopolymerMatch) {
    issues.push(`Homopolymer: ${homopolymerMatch[0]}`);
  }

  return {
    score: composite * 100,
    quality: classifyQuality(composite * 100) as unknown as string,  // FIXED: Type assertion for quality
    tm,
    gc: gc * 100,
    hairpinDG,
    homodimerDG,
    gcClampScore: gcClampScore * 100,
    homopolymerScore: homopolymerScore * 100,
    issues,
    breakdown: {
      tm: tmScore * 100,
      gc: gcScore * 100,
      hairpin: hairpinScore * 100,
      homodimer: homodimerScore * 100,
      gcClamp: gcClampScore * 100,
      homopolymer: homopolymerScore * 100,
    },
  };
}

/**
 * Get the overhang at a specific position in a sequence
 */
function getOverhangAt(sequence: string, position: number): string {
  return sequence.slice(position, position + 4).toUpperCase();
}

/**
 * Check if an overhang is acceptable (not palindrome, not homopolymer)
 */
function validateOverhang(overhang: string, options: OptimizationOptions = {}): OverhangValidation {
  const { avoidPalindromes = true, avoidHomopolymers = true } = options;

  const issues: string[] = [];
  let valid = true;

  // Check palindrome (self-complementary)
  if (avoidPalindromes) {
    const rc = reverseComplement(overhang);
    if (overhang === rc) {
      issues.push('Palindromic overhang (self-ligation risk)');
      valid = false;
    }
  }

  // Check homopolymer
  if (avoidHomopolymers) {
    if (/^(.)\1{3}$/.test(overhang)) {
      issues.push('Homopolymer overhang (poor fidelity)');
      valid = false;
    }
  }

  // Check for TNNA pattern (reduced efficiency)
  if (/^T..A$/.test(overhang)) {
    issues.push('TNNA pattern (reduced efficiency, ~70%)');
    // Don't mark invalid, just warn
  }

  return { valid, issues };
}

// ============================================================================
// MAIN OPTIMIZATION FUNCTIONS
// ============================================================================

/**
 * Analyze a junction between two fragments and find the optimal boundary position
 */
export function optimizeJunctionBoundary(
  leftSeq: string,
  rightSeq: string,
  options: OptimizationOptions = {}
): BoundaryOptimizationResult {
  const config = { ...BOUNDARY_OPTIMIZER_DEFAULTS, ...options };

  // The "full sequence" around the junction
  // leftSeq ends at position X, rightSeq starts at position X
  const fullContext = leftSeq + rightSeq;
  const originalBoundary = leftSeq.length;

  // Extract original primer binding regions
  const homologyLen = Math.min(config.maxHomologyLength, 20);

  const originalLeftHomology = leftSeq.slice(-homologyLen); // End of left fragment (reverse primer binds here)
  const originalRightHomology = rightSeq.slice(0, homologyLen); // Start of right fragment (forward primer binds here)

  const originalLeftScore = scorePrimerBindingRegion(originalLeftHomology, config);
  const originalRightScore = scorePrimerBindingRegion(originalRightHomology, config);

  const originalOverhang = getOverhangAt(fullContext, originalBoundary);
  const originalOverhangValidation = validateOverhang(originalOverhang, config);

  const originalComposite = (originalLeftScore.score + originalRightScore.score) / 2;

  // Search for better boundary positions
  let bestShift = 0;
  let bestScore = originalComposite;
  let bestResult = {
    leftScore: originalLeftScore,
    rightScore: originalRightScore,
    overhang: originalOverhang,
    overhangValidation: originalOverhangValidation,
  };

  // Search left (shift boundary into left fragment)
  for (let shift = -config.minShiftThreshold; shift >= -config.maxShift; shift--) {
    const newBoundary = originalBoundary + shift;

    // Check minimum fragment sizes
    if (newBoundary < config.minFragmentSize) continue;
    if (fullContext.length - newBoundary < config.minFragmentSize) continue;

    const newLeftHomology = fullContext.slice(Math.max(0, newBoundary - homologyLen), newBoundary);
    const newRightHomology = fullContext.slice(newBoundary, newBoundary + homologyLen);
    const newOverhang = getOverhangAt(fullContext, newBoundary);

    if (newLeftHomology.length < config.minHomologyLength) continue;
    if (newRightHomology.length < config.minHomologyLength) continue;

    const leftScore = scorePrimerBindingRegion(newLeftHomology, config);
    const rightScore = scorePrimerBindingRegion(newRightHomology, config);
    const overhangValidation = validateOverhang(newOverhang, config);

    // Skip invalid overhangs
    if (!overhangValidation.valid) continue;

    const composite = (leftScore.score + rightScore.score) / 2;

    // Apply penalty for larger shifts
    const shiftPenalty = Math.abs(shift) * 0.1;
    const adjustedScore = composite - shiftPenalty;

    if (adjustedScore > bestScore) {
      bestScore = adjustedScore;
      bestShift = shift;
      bestResult = { leftScore, rightScore, overhang: newOverhang, overhangValidation };
    }
  }

  // Search right (shift boundary into right fragment)
  for (let shift = config.minShiftThreshold; shift <= config.maxShift; shift++) {
    const newBoundary = originalBoundary + shift;

    // Check minimum fragment sizes
    if (newBoundary < config.minFragmentSize) continue;
    if (fullContext.length - newBoundary < config.minFragmentSize) continue;

    const newLeftHomology = fullContext.slice(Math.max(0, newBoundary - homologyLen), newBoundary);
    const newRightHomology = fullContext.slice(newBoundary, newBoundary + homologyLen);
    const newOverhang = getOverhangAt(fullContext, newBoundary);

    if (newLeftHomology.length < config.minHomologyLength) continue;
    if (newRightHomology.length < config.minHomologyLength) continue;

    const leftScore = scorePrimerBindingRegion(newLeftHomology, config);
    const rightScore = scorePrimerBindingRegion(newRightHomology, config);
    const overhangValidation = validateOverhang(newOverhang, config);

    // Skip invalid overhangs
    if (!overhangValidation.valid) continue;

    const composite = (leftScore.score + rightScore.score) / 2;

    // Apply penalty for larger shifts
    const shiftPenalty = Math.abs(shift) * 0.1;
    const adjustedScore = composite - shiftPenalty;

    if (adjustedScore > bestScore) {
      bestScore = adjustedScore;
      bestShift = shift;
      bestResult = { leftScore, rightScore, overhang: newOverhang, overhangValidation };
    }
  }

  // Determine direction
  let direction: 'left' | 'right' | 'unchanged' = 'unchanged';
  if (bestShift < 0) direction = 'left';
  if (bestShift > 0) direction = 'right';

  // Calculate improvement
  const improvement = Math.max(0, (bestScore - originalComposite) / 100);

  // Generate reason
  let reason = 'Original boundary is optimal';
  if (bestShift !== 0) {
    const absShift = Math.abs(bestShift);
    const leftFragment = direction === 'left' ? 'shorter' : 'longer';
    const rightFragment = direction === 'left' ? 'longer' : 'shorter';
    reason = `Shift boundary ${absShift}bp ${direction}: left fragment ${leftFragment}, right fragment ${rightFragment}. ` +
             `Primer quality improved from ${originalComposite.toFixed(1)} to ${bestScore.toFixed(1)}.`;

    // Add specific improvements
    if (bestResult.leftScore.score > originalLeftScore.score) {
      reason += ` Left primer: ${originalLeftScore.quality} → ${bestResult.leftScore.quality}.`;
    }
    if (bestResult.rightScore.score > originalRightScore.score) {
      reason += ` Right primer: ${originalRightScore.quality} → ${bestResult.rightScore.quality}.`;
    }
  }

  return {
    success: true,
    originalPosition: originalBoundary,
    optimizedPosition: originalBoundary + bestShift,
    shift: bestShift,
    direction,
    beforePrimers: {
      left: originalLeftScore,
      right: originalRightScore,
      composite: originalComposite,
      overhang: originalOverhang,
      overhangValid: originalOverhangValidation.valid,
    },
    afterPrimers: {
      left: bestResult.leftScore,
      right: bestResult.rightScore,
      composite: bestScore,
      overhang: bestResult.overhang,
      overhangValid: bestResult.overhangValidation.valid,
    },
    improvement,
    reason,
  };
}

/**
 * Optimize all boundaries in a multi-fragment assembly
 */
export function optimizeAssemblyBoundaries(
  fragments: Fragment[],
  options: OptimizationOptions = {}
): AssemblyOptimizationResult {
  const config = { ...BOUNDARY_OPTIMIZER_DEFAULTS, ...options };

  if (!Array.isArray(fragments) || fragments.length < 2) {
    return {
      success: false,
      error: 'Need at least 2 fragments to optimize boundaries',
      boundaries: [],
      optimizedFragments: [],
    };
  }

  const results: JunctionResult[] = [];
  const optimizedFragments: OptimizedFragment[] = [];

  // Running position tracker
  let cumulativeShift = 0;

  for (let i = 0; i < fragments.length - 1; i++) {
    const leftFrag = fragments[i];
    const rightFrag = fragments[i + 1];

    const leftSeq = typeof leftFrag === 'string' ? leftFrag : (leftFrag.seq || leftFrag.sequence);
    const rightSeq = typeof rightFrag === 'string' ? rightFrag : (rightFrag.seq || rightFrag.sequence);

    if (!leftSeq || !rightSeq) {
      results.push({
        junctionIndex: i,
        success: false,
        error: 'Missing sequence data',
        shift: 0,
        direction: 'unchanged',
        leftFragmentId: (leftFrag as any)?.id || `Fragment_${i + 1}`,
        rightFragmentId: (rightFrag as any)?.id || `Fragment_${i + 2}`,
      } as JunctionResult);
      continue;
    }

    // Optimize this junction
    const junctionResult = optimizeJunctionBoundary(leftSeq, rightSeq, config);
    const fullResult: JunctionResult = {
      ...junctionResult,
      junctionIndex: i,
      leftFragmentId: (leftFrag as any).id || `Fragment_${i + 1}`,
      rightFragmentId: (rightFrag as any).id || `Fragment_${i + 2}`,
    };

    results.push(fullResult);
    cumulativeShift += junctionResult.shift;
  }

  // Apply shifts to create optimized fragments
  let position = 0;
  for (let i = 0; i < fragments.length; i++) {
    const frag = fragments[i];
    const seq = typeof frag === 'string' ? frag : (frag.seq || frag.sequence);

    // Handle missing sequence
    if (!seq) {
      optimizedFragments.push({
        id: (frag as any)?.id || `Fragment_${i + 1}`,
        seq: '',
        originalSeq: '',
        originalLength: 0,
        newLength: 0,
        lengthChange: 0,
        error: 'Missing sequence data',
      });
      continue;
    }

    // Calculate adjusted start/end based on accumulated shifts
    let adjustedSeq = seq;

    // If this is not the first fragment, check if left boundary shifted
    if (i > 0 && results[i - 1]?.shift) {
      const shift = results[i - 1].shift;
      if (shift > 0) {
        // Boundary shifted right into this fragment - we gain sequence from left neighbor
        // This is handled by the previous fragment losing sequence
      } else if (shift < 0) {
        // Boundary shifted left - we lose sequence to left neighbor
        adjustedSeq = seq.slice(-shift);
      }
    }

    // If this is not the last fragment, check if right boundary shifted
    if (i < fragments.length - 1 && results[i]?.shift) {
      const shift = results[i].shift;
      if (shift < 0) {
        // Boundary shifted left - we gain sequence from right neighbor
        // We need the original sequence plus the borrowed portion
        const borrowedLen = -shift;
        const nextFrag = fragments[i + 1];
        const nextSeq = typeof nextFrag === 'string' ? nextFrag : (nextFrag.seq || nextFrag.sequence);
        adjustedSeq = adjustedSeq + (nextSeq?.slice(0, borrowedLen) || '');
      } else if (shift > 0) {
        // Boundary shifted right - we lose sequence to right neighbor
        adjustedSeq = adjustedSeq.slice(0, -shift);
      }
    }

    optimizedFragments.push({
      id: (frag as any).id || `Fragment_${i + 1}`,
      seq: adjustedSeq,
      originalSeq: seq,
      originalLength: seq.length,
      newLength: adjustedSeq.length,
      lengthChange: adjustedSeq.length - seq.length,
    });
  }

  // Calculate summary statistics
  const totalImprovement = results.reduce((sum, r) => sum + (r.improvement || 0), 0);
  const avgImprovement = totalImprovement / results.length;
  const boundariesShifted = results.filter(r => r.shift !== 0).length;

  return {
    success: true,
    boundaries: results,
    optimizedFragments,
    summary: {
      totalBoundaries: results.length,
      boundariesOptimized: boundariesShifted,
      boundariesUnchanged: results.length - boundariesShifted,
      averageImprovement: avgImprovement,
      totalImprovement,
    },
  };
}

/**
 * Quick check if boundary optimization would help
 * Useful for UI to show optimization suggestions
 */
export function assessBoundaryOptimizationPotential(
  fragments: Fragment[],
  options: OptimizationOptions = {}
): BoundaryAssessmentResult {
  const config = { ...BOUNDARY_OPTIMIZER_DEFAULTS, ...options };
  const homologyLen = 20;

  const issues: AssessmentIssue[] = [];
  let totalScore = 0;
  let junctionCount = 0;

  for (let i = 0; i < fragments.length - 1; i++) {
    const leftFrag = fragments[i];
    const rightFrag = fragments[i + 1];

    const leftSeq = typeof leftFrag === 'string' ? leftFrag : (leftFrag.seq || leftFrag.sequence);
    const rightSeq = typeof rightFrag === 'string' ? rightFrag : (rightFrag.seq || rightFrag.sequence);

    if (!leftSeq || !rightSeq) continue;

    const leftHomology = leftSeq.slice(-homologyLen);
    const rightHomology = rightSeq.slice(0, homologyLen);

    const leftScore = scorePrimerBindingRegion(leftHomology, config);
    const rightScore = scorePrimerBindingRegion(rightHomology, config);

    const avgScore = (leftScore.score + rightScore.score) / 2;
    totalScore += avgScore;
    junctionCount++;

    // Track issues
    if (leftScore.quality === 'poor' || leftScore.quality === 'acceptable') {
      issues.push({
        junction: i,
        side: 'left',
        fragmentId: (leftFrag as any).id || `Fragment_${i + 1}`,
        quality: leftScore.quality,
        score: leftScore.score,
        issues: leftScore.issues,
      });
    }
    if (rightScore.quality === 'poor' || rightScore.quality === 'acceptable') {
      issues.push({
        junction: i,
        side: 'right',
        fragmentId: (rightFrag as any).id || `Fragment_${i + 2}`,
        quality: rightScore.quality,
        score: rightScore.score,
        issues: rightScore.issues,
      });
    }
  }

  const avgScore = junctionCount > 0 ? totalScore / junctionCount : 0;
  const needsOptimization = issues.length > 0;

  return {
    needsOptimization,
    junctionCount,
    averageScore: avgScore,
    overallQuality: classifyQuality(avgScore) as unknown as string,  // FIXED: Type assertion for quality
    issues,
    recommendation: needsOptimization
      ? `${issues.length} primer binding region(s) have quality issues. Boundary optimization recommended.`
      : 'All primer binding regions have good quality. No optimization needed.',
  };
}

// ============================================================================
// EXPORTS
// ============================================================================

// Named exports for all public functions
export {
  scorePrimerBindingRegion,
  validateOverhang,
};

export default {
  optimizeJunctionBoundary,
  optimizeAssemblyBoundaries,
  assessBoundaryOptimizationPotential,
  scorePrimerBindingRegion,
  validateOverhang,
  BOUNDARY_OPTIMIZER_DEFAULTS,
};
