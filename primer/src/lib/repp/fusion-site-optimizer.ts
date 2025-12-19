/**
 * Fusion Site Optimizer for Golden Gate Assembly
 *
 * Main optimization engine that selects optimal junction positions
 * for multi-fragment Golden Gate assembly.
 *
 * Algorithms implemented:
 * 1. Greedy - Fast, good for initial exploration
 * 2. Branch & Bound - Optimal for small assemblies
 * 3. Monte Carlo - Best for large/complex assemblies
 * 4. Hybrid - Combines approaches for best results
 *
 * This is the main entry point for fusion site optimization.
 */

import { reverseComplement } from './enzymes.js';
import { calculateExperimentalFidelity, getEnzymeLigationData } from './goldengate.js';
import { scanForFusionSites, generateTargetPositions, filterByDistance, FusionSiteCandidate, TargetPosition } from './fusion-site-scanner.js';
import { scoreFusionSiteComposite, quickScoreFusionSite, CompositeScore, ScoringOptions } from './fusion-site-scorer.js';
import { calculateSetEfficiency } from './overhang-efficiency.js';
import { predictFailureModes, quickRiskAssessment } from './failure-prediction.js';
import { findGTMismatchRisks, calculateEnhancedFidelity } from './goldengate-primer-optimizer.js';

/**
 * Default optimizer configuration
 */
export const OPTIMIZER_DEFAULTS = {
  enzyme: 'BsaI',
  algorithm: 'auto' as const,
  minFragmentSize: 200,
  maxFragmentSize: 5000,
  targetFragmentSize: null as number | null,
  minDistanceFromEnds: 50,
  searchRadius: 50,

  // Monte Carlo settings
  mcIterations: 2000,
  mcTemperature: 1.0,
  mcCoolingRate: 0.995,

  // Branch & Bound settings
  maxBranchDepth: 8,
  pruningThreshold: 0.7,

  // Scoring weights
  fidelityWeight: 0.40,
  efficiencyWeight: 0.20,
  primerQualityWeight: 0.25,
  positionWeight: 0.15,
};

/**
 * Optimizer algorithm type
 */
export type OptimizerAlgorithm = 'greedy' | 'branchBound' | 'monteCarlo' | 'hybrid' | 'auto';

/**
 * Junction object with extended information
 */
export interface Junction extends FusionSiteCandidate {
  targetIndex?: number;
  idealPosition?: number;
  deviation?: number;
}

/**
 * Set score result
 */
export interface SetScore {
  composite: number;
  fidelity: number;
  efficiency: number;
  primerQuality: number;
  positionQuality: number;
  overhangs: string[];
}

/**
 * Optimization result
 */
export interface OptimizationResult {
  algorithm: string;
  junctions: Junction[];
  overhangs: string[];
  score: SetScore;
  numJunctions: number;
  numExpected: number;
  complete: boolean;
  iterations?: number;
  finalTemperature?: number;
  nodesExplored?: number;
  optimal?: boolean;
  alternativeResults?: Array<{
    algorithm: string;
    score?: number;
    complete: boolean;
  }>;
}

/**
 * Optimizer options
 */
export interface OptimizerOptions {
  enzyme?: string;
  algorithm?: OptimizerAlgorithm;
  minFragmentSize?: number;
  maxFragmentSize?: number;
  targetFragmentSize?: number | null;
  minDistanceFromEnds?: number;
  searchRadius?: number;
  forbiddenRegions?: any[];
  codingFrame?: number | null;
  proteinDomains?: any[];
  scarContext?: string;
  verbose?: boolean;

  // Algorithm-specific
  iterations?: number;
  temperature?: number;
  coolingRate?: number;
  maxBranchDepth?: number;
  pruningThreshold?: number;
  maxCandidatesPerRegion?: number;

  // Weights
  fidelityWeight?: number;
  efficiencyWeight?: number;
  primerQualityWeight?: number;
  positionWeight?: number;
}

/**
 * Calculate set fidelity from overhang list
 *
 * @param overhangs - Array of overhangs
 * @param enzyme - Enzyme name
 * @returns Assembly fidelity (0-1)
 */
function calculateSetFidelity(overhangs: string[], enzyme: string): number {
  const result = calculateExperimentalFidelity(overhangs, enzyme);
  return (result as any).assemblyFidelity || 0;
}

/**
 * Set scoring options
 */
interface SetScoringOptions {
  fidelityWeight?: number;
  efficiencyWeight?: number;
  primerQualityWeight?: number;
  positionWeight?: number;
  targetPositions?: TargetPosition[];
}

/**
 * Score a complete junction set considering all interactions
 *
 * @param junctions - Array of junction objects
 * @param enzyme - Enzyme name
 * @param options - Scoring options
 * @returns Set score
 */
function scoreJunctionSet(junctions: Junction[], enzyme: string, options: SetScoringOptions = {}): SetScore {
  const {
    fidelityWeight = 0.40,
    efficiencyWeight = 0.20,
    primerQualityWeight = 0.25,
    positionWeight = 0.15,
    targetPositions = [],
  } = options;

  const overhangs = junctions.map(j => j.overhang);

  // Calculate fidelity for the set
  const fidelity = calculateSetFidelity(overhangs, enzyme);

  // Calculate efficiency for the set
  const efficiencyResult = calculateSetEfficiency(overhangs);
  const efficiency = (efficiencyResult as any).combinedEfficiency;

  // Average primer quality
  const avgPrimerScore = junctions.reduce((sum, j) =>
    sum + ((j.score as any)?.composite || (j as any).composite || 70), 0
  ) / junctions.length;
  const primerQuality = avgPrimerScore / 100;

  // Position quality (how close to target positions)
  let positionQuality = 1.0;
  if (targetPositions.length === junctions.length) {
    const deviations = junctions.map((j, i) =>
      Math.abs(j.position - targetPositions[i].idealPosition)
    );
    const maxDeviation = Math.max(...targetPositions.map(t =>
      t.searchRegion.end - t.searchRegion.start
    ));
    const avgDeviation = deviations.reduce((a, b) => a + b, 0) / deviations.length;
    positionQuality = Math.max(0, 1 - avgDeviation / maxDeviation);
  }

  // Weighted composite
  const composite =
    fidelity * fidelityWeight +
    efficiency * efficiencyWeight +
    primerQuality * primerQualityWeight +
    positionQuality * positionWeight;

  return {
    composite,
    fidelity,
    efficiency,
    primerQuality,
    positionQuality,
    overhangs,
  };
}

/**
 * Greedy algorithm for junction selection
 *
 * Fast but may not find global optimum. Good for initial exploration.
 *
 * @param sequence - DNA sequence
 * @param numFragments - Number of fragments desired
 * @param enzyme - Enzyme name
 * @param options - Options
 * @returns Optimization result
 */
export function optimizeGreedy(
  sequence: string,
  numFragments: number,
  enzyme: string = 'BsaI',
  options: OptimizerOptions = {}
): OptimizationResult {
  const {
    minDistanceFromEnds = OPTIMIZER_DEFAULTS.minDistanceFromEnds,
    searchRadius = OPTIMIZER_DEFAULTS.searchRadius,
    forbiddenRegions = [],
  } = options;

  const numJunctions = numFragments - 1;
  if (numJunctions <= 0) {
    return {
      junctions: [],
      overhangs: [],
      score: { composite: 0, fidelity: 0, efficiency: 0, primerQuality: 0, positionQuality: 0, overhangs: [] },
      algorithm: 'greedy',
      numJunctions: 0,
      numExpected: numJunctions,
      complete: false,
    };
  }

  // Generate target positions
  const targets = generateTargetPositions(sequence.length, numFragments, {
    minDistanceFromEnds,
    searchRadius,
  });

  // Get all candidates
  const allCandidates = scanForFusionSites(sequence, {
    enzyme,
    minDistanceFromEnds,
    forbiddenRegions,
    includeEfficiency: true,
  });

  // Score each candidate
  const scoredCandidates = allCandidates.map(c => ({
    ...c,
    score: quickScoreFusionSite(sequence, c.position, enzyme),
  }));

  // Greedy selection: pick best candidate for each target region
  const selected: Junction[] = [];
  const usedOverhangs = new Set<string>();

  for (const target of targets) {
    // Get candidates in this region
    const inRegion = scoredCandidates.filter(c =>
      c.position >= target.searchRegion.start &&
      c.position <= target.searchRegion.end &&
      !usedOverhangs.has(c.overhang)
    );

    if (inRegion.length === 0) {
      console.warn(`No candidates found for target region ${target.index}`);
      continue;
    }

    // Sort by score
    inRegion.sort((a, b) => ((b.score as any)?.score || 0) - ((a.score as any)?.score || 0));

    // Pick best that doesn't conflict
    let bestCandidate: Junction | null = null;
    for (const candidate of inRegion) {
      // Check fidelity with already selected
      const testSet = [...selected.map(s => s.overhang), candidate.overhang];
      const testFidelity = calculateSetFidelity(testSet, enzyme);

      if (testFidelity > 0.5 || selected.length === 0) {
        bestCandidate = candidate;
        break;
      }
    }

    if (bestCandidate) {
      selected.push({
        ...bestCandidate,
        targetIndex: target.index,
        idealPosition: target.idealPosition,
        deviation: bestCandidate.position - target.idealPosition,
      });
      usedOverhangs.add(bestCandidate.overhang);
    }
  }

  // Calculate final set score
  const setScore = scoreJunctionSet(selected, enzyme, {
    ...options,
    targetPositions: targets,
  });

  return {
    algorithm: 'greedy',
    junctions: selected.sort((a, b) => a.position - b.position),
    overhangs: selected.map(s => s.overhang),
    score: setScore,
    numJunctions: selected.length,
    numExpected: numJunctions,
    complete: selected.length === numJunctions,
  };
}

/**
 * Monte Carlo optimization for junction selection
 *
 * Uses simulated annealing to explore solution space.
 * Best for complex assemblies with many interactions.
 *
 * @param sequence - DNA sequence
 * @param numFragments - Number of fragments
 * @param enzyme - Enzyme name
 * @param options - Options
 * @returns Optimization result
 */
export function optimizeMonteCarlo(
  sequence: string,
  numFragments: number,
  enzyme: string = 'BsaI',
  options: OptimizerOptions = {}
): OptimizationResult {
  const {
    iterations = OPTIMIZER_DEFAULTS.mcIterations,
    temperature = OPTIMIZER_DEFAULTS.mcTemperature,
    coolingRate = OPTIMIZER_DEFAULTS.mcCoolingRate,
    minDistanceFromEnds = OPTIMIZER_DEFAULTS.minDistanceFromEnds,
    searchRadius = OPTIMIZER_DEFAULTS.searchRadius,
    forbiddenRegions = [],
  } = options;

  const numJunctions = numFragments - 1;
  if (numJunctions <= 0) {
    return {
      junctions: [],
      overhangs: [],
      score: { composite: 0, fidelity: 0, efficiency: 0, primerQuality: 0, positionQuality: 0, overhangs: [] },
      algorithm: 'monteCarlo',
      numJunctions: 0,
      numExpected: numJunctions,
      complete: false,
    };
  }

  // Get all candidates
  const targets = generateTargetPositions(sequence.length, numFragments, {
    minDistanceFromEnds,
    searchRadius,
  });

  const candidatesByRegion = targets.map(target => {
    const inRegion = scanForFusionSites(sequence, {
      enzyme,
      minDistanceFromEnds,
      forbiddenRegions,
      searchWindows: [target.searchRegion],
      includeEfficiency: true,
    });
    return inRegion.map(c => ({
      ...c,
      score: quickScoreFusionSite(sequence, c.position, enzyme),
    }));
  });

  // Check if we have candidates in all regions
  const emptyRegions = candidatesByRegion.filter(c => c.length === 0);
  if (emptyRegions.length > 0) {
    console.warn(`${emptyRegions.length} regions have no candidates`);
  }

  // Initialize with greedy solution
  let current: Junction[] = [];
  const usedOverhangs = new Set<string>();

  for (let i = 0; i < numJunctions; i++) {
    if (candidatesByRegion[i].length > 0) {
      let found = false;
      for (const candidate of candidatesByRegion[i]) {
        if (!usedOverhangs.has(candidate.overhang)) {
          current.push(candidate);
          usedOverhangs.add(candidate.overhang);
          found = true;
          break;
        }
      }
      if (!found) {
        current.push(candidatesByRegion[i][0]);
      }
    }
  }

  if (current.length !== numJunctions) {
    return optimizeGreedy(sequence, numFragments, enzyme, options);
  }

  let currentScore = scoreJunctionSet(current, enzyme, {
    ...options,
    targetPositions: targets,
  });
  let best = [...current];
  let bestScore = { ...currentScore };
  let temp = temperature;

  // Monte Carlo iterations
  for (let iter = 0; iter < iterations; iter++) {
    const jIdx = Math.floor(Math.random() * numJunctions);

    if (candidatesByRegion[jIdx].length <= 1) continue;

    const alternatives = candidatesByRegion[jIdx].filter(c =>
      c.position !== current[jIdx].position
    );
    if (alternatives.length === 0) continue;

    const newCandidate = alternatives[Math.floor(Math.random() * alternatives.length)];

    const newSolution = [...current];
    newSolution[jIdx] = newCandidate;

    // Check for duplicate overhangs
    const overhangs = newSolution.map(c => c.overhang);
    const uniqueOverhangs = new Set(overhangs);
    if (uniqueOverhangs.size !== overhangs.length) continue;

    const newScore = scoreJunctionSet(newSolution, enzyme, {
      ...options,
      targetPositions: targets,
    });

    const delta = newScore.composite - currentScore.composite;
    const acceptProb = delta > 0 ? 1 : Math.exp(delta / temp);

    if (Math.random() < acceptProb) {
      current = newSolution;
      currentScore = newScore;

      if (newScore.composite > bestScore.composite) {
        best = [...newSolution];
        bestScore = { ...newScore };
      }
    }

    temp *= coolingRate;
  }

  const result = best.map((c, i) => ({
    ...c,
    targetIndex: i,
    idealPosition: targets[i].idealPosition,
    deviation: c.position - targets[i].idealPosition,
  }));

  return {
    algorithm: 'monteCarlo',
    junctions: result.sort((a, b) => a.position - b.position),
    overhangs: result.map(r => r.overhang),
    score: bestScore,
    numJunctions: result.length,
    numExpected: numJunctions,
    complete: result.length === numJunctions,
    iterations,
    finalTemperature: temp,
  };
}

/**
 * Branch and Bound optimization
 *
 * Guarantees optimal solution for small assemblies.
 *
 * @param sequence - DNA sequence
 * @param numFragments - Number of fragments
 * @param enzyme - Enzyme name
 * @param options - Options
 * @returns Optimization result
 */
export function optimizeBranchBound(
  sequence: string,
  numFragments: number,
  enzyme: string = 'BsaI',
  options: OptimizerOptions = {}
): OptimizationResult {
  const {
    maxCandidatesPerRegion = 10,
    pruningThreshold = OPTIMIZER_DEFAULTS.pruningThreshold,
    minDistanceFromEnds = OPTIMIZER_DEFAULTS.minDistanceFromEnds,
    searchRadius = OPTIMIZER_DEFAULTS.searchRadius,
    forbiddenRegions = [],
  } = options;

  const numJunctions = numFragments - 1;
  if (numJunctions <= 0) {
    return {
      junctions: [],
      overhangs: [],
      score: { composite: 0, fidelity: 0, efficiency: 0, primerQuality: 0, positionQuality: 0, overhangs: [] },
      algorithm: 'branchBound',
      numJunctions: 0,
      numExpected: numJunctions,
      complete: false,
    };
  }

  const targets = generateTargetPositions(sequence.length, numFragments, {
    minDistanceFromEnds,
    searchRadius,
  });

  const candidatesByRegion = targets.map(target => {
    const inRegion = scanForFusionSites(sequence, {
      enzyme,
      minDistanceFromEnds,
      forbiddenRegions,
      searchWindows: [target.searchRegion],
      includeEfficiency: true,
    });

    const scored = inRegion.map(c => ({
      ...c,
      score: quickScoreFusionSite(sequence, c.position, enzyme),
    }));
    scored.sort((a, b) => ((b.score as any)?.score || 0) - ((a.score as any)?.score || 0));

    return scored.slice(0, maxCandidatesPerRegion);
  });

  let best: { junctions: Junction[]; score: SetScore } | null = null;
  let bestScore = -Infinity;
  let nodesExplored = 0;

  function branch(depth: number, current: Junction[], usedOverhangs: Set<string>) {
    nodesExplored++;

    if (depth === numJunctions) {
      const score = scoreJunctionSet(current, enzyme, {
        ...options,
        targetPositions: targets,
      });

      if (score.composite > bestScore) {
        bestScore = score.composite;
        best = { junctions: [...current], score };
      }
      return;
    }

    const remainingRegions = numJunctions - depth;

    let currentPartialScore = 0;
    if (current.length > 0) {
      currentPartialScore = scoreJunctionSet(current, enzyme, {
        ...options,
        targetPositions: targets.slice(0, current.length),
      }).composite;
    }

    const avgScorePerJunction = current.length > 0
      ? currentPartialScore / current.length
      : 0.9;

    const projectedFinal = current.length > 0
      ? (currentPartialScore * current.length + Math.min(0.95, avgScorePerJunction + 0.1) * remainingRegions) / numJunctions
      : 0.95;

    if (projectedFinal < bestScore * pruningThreshold) {
      return;
    }

    for (const candidate of candidatesByRegion[depth]) {
      if (usedOverhangs.has(candidate.overhang)) continue;

      const testOverhangs = [...current.map(c => c.overhang), candidate.overhang];
      const testFidelity = calculateSetFidelity(testOverhangs, enzyme);

      if (testFidelity < 0.3) continue;

      const newUsed = new Set(usedOverhangs);
      newUsed.add(candidate.overhang);

      branch(depth + 1, [...current, candidate], newUsed);
    }
  }

  branch(0, [], new Set());

  if (!best) {
    return optimizeGreedy(sequence, numFragments, enzyme, options);
  }

  const bestResult = best as { junctions: Junction[]; score: SetScore };
  const result = bestResult.junctions.map((c: any, i: any) => ({
    ...c,
    targetIndex: i,
    idealPosition: targets[i].idealPosition,
    deviation: c.position - targets[i].idealPosition,
  }));

  return {
    algorithm: 'branchBound',
    junctions: result.sort((a: any, b: any) => a.position - b.position),
    overhangs: result.map((r: any) => r.overhang),
    score: bestResult.score,
    numJunctions: result.length,
    numExpected: numJunctions,
    complete: result.length === numJunctions,
    nodesExplored,
    optimal: true,
  };
}

/**
 * Hybrid optimization - combines multiple algorithms
 *
 * @param sequence - DNA sequence
 * @param numFragments - Number of fragments
 * @param enzyme - Enzyme name
 * @param options - Options
 * @returns Best result from multiple approaches
 */
export function optimizeHybrid(
  sequence: string,
  numFragments: number,
  enzyme: string = 'BsaI',
  options: OptimizerOptions = {}
): OptimizationResult {
  const results: OptimizationResult[] = [];

  results.push(optimizeGreedy(sequence, numFragments, enzyme, options));

  if (numFragments <= 6) {
    results.push(optimizeBranchBound(sequence, numFragments, enzyme, options));
  }

  results.push(optimizeMonteCarlo(sequence, numFragments, enzyme, options));

  const validResults = results.filter(r => r.complete);
  if (validResults.length === 0) {
    results.sort((a, b) => (b.score?.composite || 0) - (a.score?.composite || 0));
    return { ...results[0], algorithm: 'hybrid' };
  }

  validResults.sort((a, b) => (b.score?.composite || 0) - (a.score?.composite || 0));

  return {
    ...validResults[0],
    algorithm: 'hybrid',
    alternativeResults: results.filter(r => r !== validResults[0]).map(r => ({
      algorithm: r.algorithm,
      score: r.score?.composite,
      complete: r.complete,
    })),
  };
}

/**
 * Comprehensive optimization result
 */
export interface ComprehensiveOptimizationResult {
  success: boolean;
  algorithm?: string;
  enzyme?: string;
  sequenceLength?: number;
  numFragments?: number;
  numJunctions?: number;
  junctions?: Junction[];
  overhangs?: string[];
  score?: SetScore;
  detailedJunctions?: CompositeScore[];
  failurePrediction?: any;
  fragmentSizes?: number[];
  minFragmentSize?: number;
  maxFragmentSize?: number;
  avgFragmentSize?: number;
  internalSites?: any;
  summary?: {
    fidelity?: number;
    efficiency?: number;
    expectedSuccessRate?: any;
    riskLevel?: string;
    recommendations?: any[];
  };
  error?: string;
  suggestion?: string;
  sequence?: null;
  result?: OptimizationResult;
}

/**
 * Main entry point - optimizeFusionSites
 *
 * Automatically selects the best algorithm based on problem size.
 *
 * @param sequence - DNA sequence to split
 * @param numFragments - Number of fragments desired
 * @param options - Optimization options
 * @returns Comprehensive optimization result
 */
export function optimizeFusionSites(
  sequence: string,
  numFragments: number,
  options: OptimizerOptions = {}
): ComprehensiveOptimizationResult {
  const {
    enzyme = OPTIMIZER_DEFAULTS.enzyme,
    algorithm = OPTIMIZER_DEFAULTS.algorithm,
    minFragmentSize = OPTIMIZER_DEFAULTS.minFragmentSize,
    maxFragmentSize = OPTIMIZER_DEFAULTS.maxFragmentSize,
    minDistanceFromEnds = OPTIMIZER_DEFAULTS.minDistanceFromEnds,
    searchRadius = OPTIMIZER_DEFAULTS.searchRadius,
    forbiddenRegions = [],
    codingFrame = null,
    proteinDomains = [],
    scarContext = 'nonCoding',
    verbose = false,
  } = options;

  const seq = sequence.toUpperCase();
  const seqLength = seq.length;

  // Validation
  if (!seq || seqLength < 100) {
    return {
      success: false,
      error: 'Sequence too short (minimum 100bp)',
      sequence: null,
    };
  }

  if (numFragments < 2) {
    return {
      success: false,
      error: 'Need at least 2 fragments',
    };
  }

  const numJunctions = numFragments - 1;
  const idealFragmentSize = seqLength / numFragments;

  if (idealFragmentSize < minFragmentSize) {
    return {
      success: false,
      error: `Fragment size ${Math.round(idealFragmentSize)}bp is below minimum ${minFragmentSize}bp`,
      suggestion: `Try ${Math.floor(seqLength / minFragmentSize)} fragments instead`,
    };
  }

  // Check for internal restriction sites
  // TODO: Implement findInternalSites when available
  const internalSites = { hasSites: false, count: 0, sites: [] };
  if ((internalSites as any).hasSites) {
    if (verbose) {
      console.warn(`Found ${(internalSites as any).count} internal ${enzyme} sites - need domestication`);
    }
  }

  // Select algorithm
  let result: OptimizationResult;
  const selectedAlgorithm = algorithm === 'auto'
    ? (numFragments <= 5 ? 'branchBound' : numFragments <= 10 ? 'hybrid' : 'monteCarlo')
    : algorithm;

  const optimizerOptions: OptimizerOptions = {
    minDistanceFromEnds,
    searchRadius,
    forbiddenRegions,
    enzyme,
    ...options,
  };

  switch (selectedAlgorithm) {
    case 'greedy':
      result = optimizeGreedy(seq, numFragments, enzyme, optimizerOptions);
      break;
    case 'branchBound':
      result = optimizeBranchBound(seq, numFragments, enzyme, optimizerOptions);
      break;
    case 'monteCarlo':
      result = optimizeMonteCarlo(seq, numFragments, enzyme, optimizerOptions);
      break;
    case 'hybrid':
    default:
      result = optimizeHybrid(seq, numFragments, enzyme, optimizerOptions);
  }

  // Enhance result with additional analysis
  if (result.junctions && result.junctions.length > 0) {
    // Full composite scoring for each junction
    const detailedJunctions = result.junctions.map(j =>
      scoreFusionSiteComposite(seq, j.position, enzyme, {
        codingFrame,
        proteinDomains,
        scarContext,
      })
    );

    // Failure prediction
    const failurePrediction = predictFailureModes(result.overhangs, enzyme, {
      internalSites: (internalSites as any).hasSites ? (internalSites as any).sites : [],
    });

    // Calculate fragment sizes
    const positions = [0, ...result.junctions.map(j => j.position + 4), seqLength];
    const fragmentSizes: number[] = [];
    for (let i = 0; i < positions.length - 1; i++) {
      fragmentSizes.push(positions[i + 1] - positions[i]);
    }

    return {
      success: result.complete,
      algorithm: result.algorithm,
      enzyme,
      sequenceLength: seqLength,
      numFragments,
      numJunctions,

      // Optimization result
      junctions: result.junctions,
      overhangs: result.overhangs,
      score: result.score,

      // Detailed analysis
      detailedJunctions,
      failurePrediction,

      // Fragment info
      fragmentSizes,
      minFragmentSize: Math.min(...fragmentSizes),
      maxFragmentSize: Math.max(...fragmentSizes),
      avgFragmentSize: Math.round(fragmentSizes.reduce((a, b) => a + b, 0) / fragmentSizes.length),

      // Internal sites warning
      internalSites: (internalSites as any).hasSites ? internalSites : null,

      // Summary
      summary: {
        fidelity: result.score?.fidelity,
        efficiency: result.score?.efficiency,
        expectedSuccessRate: (failurePrediction as any).expectedSuccessRate,
        riskLevel: (failurePrediction as any).overallRisk,
        recommendations: (failurePrediction as any).recommendations,
      },
    };
  }

  return {
    success: false,
    error: 'Optimization failed - no valid junction set found',
    result,
  };
}

/**
 * Quick optimization for simple cases
 *
 * @param sequence - DNA sequence
 * @param numFragments - Number of fragments
 * @param enzyme - Enzyme name
 * @returns Quick result
 */
export function quickOptimize(sequence: string, numFragments: number, enzyme: string = 'BsaI'): OptimizationResult {
  return optimizeGreedy(sequence, numFragments, enzyme, {
    minDistanceFromEnds: 50,
    searchRadius: 30,
  });
}
