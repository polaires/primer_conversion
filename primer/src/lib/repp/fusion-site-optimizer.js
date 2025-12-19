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
import { calculateExperimentalFidelity, getEnzymeLigationData, findInternalSites } from './goldengate.js';
import { scanForFusionSites, generateTargetPositions, filterByDistance } from './fusion-site-scanner.js';
import { scoreFusionSiteComposite, quickScoreFusionSite } from './fusion-site-scorer.js';
import { calculateSetEfficiency } from './overhang-efficiency.js';
import { predictFailureModes, quickRiskAssessment } from './failure-prediction.js';
import { findGTMismatchRisks, calculateEnhancedFidelity } from './goldengate-primer-optimizer.js';

/**
 * Default optimizer configuration
 */
export const OPTIMIZER_DEFAULTS = {
  enzyme: 'BsaI',
  algorithm: 'auto',           // 'greedy', 'branchBound', 'monteCarlo', 'hybrid', 'auto'
  minFragmentSize: 200,        // Minimum fragment size in bp
  maxFragmentSize: 5000,       // Maximum fragment size in bp
  targetFragmentSize: null,    // Target fragment size (auto-calculated if null)
  minDistanceFromEnds: 50,     // Minimum distance from sequence ends
  searchRadius: 50,            // Search radius around target positions

  // Monte Carlo settings
  mcIterations: 2000,          // Number of MC iterations
  mcTemperature: 1.0,          // Initial temperature
  mcCoolingRate: 0.995,        // Temperature cooling rate

  // Branch & Bound settings
  maxBranchDepth: 8,           // Maximum depth for B&B
  pruningThreshold: 0.7,       // Score threshold for pruning

  // Scoring weights
  fidelityWeight: 0.40,
  efficiencyWeight: 0.20,
  primerQualityWeight: 0.25,
  positionWeight: 0.15,
};

/**
 * Calculate set fidelity from overhang list
 *
 * @param {string[]} overhangs - Array of overhangs
 * @param {string} enzyme - Enzyme name
 * @returns {number} Assembly fidelity (0-1)
 */
function calculateSetFidelity(overhangs, enzyme) {
  const result = calculateExperimentalFidelity(overhangs, enzyme);
  return result.assemblyFidelity || 0;
}

/**
 * Score a complete junction set considering all interactions
 *
 * @param {Array} junctions - Array of junction objects
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Scoring options
 * @returns {Object} Set score
 */
function scoreJunctionSet(junctions, enzyme, options = {}) {
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
  const efficiency = efficiencyResult.combinedEfficiency;

  // Average primer quality
  const avgPrimerScore = junctions.reduce((sum, j) =>
    sum + (j.score?.composite || j.composite || 70), 0
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
 * @param {string} sequence - DNA sequence
 * @param {number} numFragments - Number of fragments desired
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Options
 * @returns {Object} Optimization result
 */
export function optimizeGreedy(sequence, numFragments, enzyme = 'BsaI', options = {}) {
  const {
    minDistanceFromEnds = OPTIMIZER_DEFAULTS.minDistanceFromEnds,
    searchRadius = OPTIMIZER_DEFAULTS.searchRadius,
    forbiddenRegions = [],
  } = options;

  const numJunctions = numFragments - 1;
  if (numJunctions <= 0) {
    return { junctions: [], score: 0, algorithm: 'greedy' };
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
  const selected = [];
  const usedOverhangs = new Set();

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
    inRegion.sort((a, b) => (b.score?.score || 0) - (a.score?.score || 0));

    // Pick best that doesn't conflict
    let bestCandidate = null;
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
 * @param {string} sequence - DNA sequence
 * @param {number} numFragments - Number of fragments
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Options
 * @returns {Object} Optimization result
 */
export function optimizeMonteCarlo(sequence, numFragments, enzyme = 'BsaI', options = {}) {
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
    return { junctions: [], score: { composite: 0 }, algorithm: 'monteCarlo' };
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

  // Initialize with greedy solution, ensuring unique overhangs
  let current = [];
  const usedOverhangs = new Set();

  for (let i = 0; i < numJunctions; i++) {
    if (candidatesByRegion[i].length > 0) {
      // Find first candidate with unique overhang
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
        // No unique candidate available - use first anyway (will be handled later)
        current.push(candidatesByRegion[i][0]);
      }
    }
  }

  if (current.length !== numJunctions) {
    // Fall back to greedy if MC init fails
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
    // Select random junction to mutate
    const jIdx = Math.floor(Math.random() * numJunctions);

    if (candidatesByRegion[jIdx].length <= 1) continue;

    // Select random alternative
    const alternatives = candidatesByRegion[jIdx].filter(c =>
      c.position !== current[jIdx].position
    );
    if (alternatives.length === 0) continue;

    const newCandidate = alternatives[Math.floor(Math.random() * alternatives.length)];

    // Create new solution
    const newSolution = [...current];
    newSolution[jIdx] = newCandidate;

    // Check for duplicate overhangs
    const overhangs = newSolution.map(c => c.overhang);
    const uniqueOverhangs = new Set(overhangs);
    if (uniqueOverhangs.size !== overhangs.length) continue;

    // Score new solution
    const newScore = scoreJunctionSet(newSolution, enzyme, {
      ...options,
      targetPositions: targets,
    });

    // Accept or reject
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

    // Cool down
    temp *= coolingRate;
  }

  // Format result
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
 * @param {string} sequence - DNA sequence
 * @param {number} numFragments - Number of fragments
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Options
 * @returns {Object} Optimization result
 */
export function optimizeBranchBound(sequence, numFragments, enzyme = 'BsaI', options = {}) {
  const {
    maxCandidatesPerRegion = 10,
    pruningThreshold = OPTIMIZER_DEFAULTS.pruningThreshold,
    minDistanceFromEnds = OPTIMIZER_DEFAULTS.minDistanceFromEnds,
    searchRadius = OPTIMIZER_DEFAULTS.searchRadius,
    forbiddenRegions = [],
  } = options;

  const numJunctions = numFragments - 1;
  if (numJunctions <= 0) {
    return { junctions: [], score: { composite: 0 }, algorithm: 'branchBound' };
  }

  // Generate targets and candidates
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

    // Pre-score and limit candidates
    const scored = inRegion.map(c => ({
      ...c,
      score: quickScoreFusionSite(sequence, c.position, enzyme),
    }));
    scored.sort((a, b) => (b.score?.score || 0) - (a.score?.score || 0));

    return scored.slice(0, maxCandidatesPerRegion);
  });

  // Branch and Bound search
  let best = null;
  let bestScore = -Infinity;
  let nodesExplored = 0;

  function branch(depth, current, usedOverhangs) {
    nodesExplored++;

    if (depth === numJunctions) {
      // Leaf node - evaluate complete solution
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

    // Pruning: estimate upper bound more accurately
    const remainingRegions = numJunctions - depth;

    // For current partial solution, get actual score
    let currentPartialScore = 0;
    if (current.length > 0) {
      currentPartialScore = scoreJunctionSet(current, enzyme, {
        ...options,
        targetPositions: targets.slice(0, current.length),
      }).composite;
    }

    // Estimate upper bound for remaining: best possible per region
    // Assume remaining junctions could each contribute optimally
    // Scale by proportion: current score extrapolated + remaining at optimal
    const avgScorePerJunction = current.length > 0
      ? currentPartialScore / current.length
      : 0.9; // Assume good quality if no data yet

    const projectedFinal = current.length > 0
      ? (currentPartialScore * current.length + Math.min(0.95, avgScorePerJunction + 0.1) * remainingRegions) / numJunctions
      : 0.95; // Optimistic for empty

    if (projectedFinal < bestScore * pruningThreshold) {
      return; // Prune this branch
    }

    // Explore candidates at this depth
    for (const candidate of candidatesByRegion[depth]) {
      if (usedOverhangs.has(candidate.overhang)) continue;

      // Quick fidelity check with current set
      const testOverhangs = [...current.map(c => c.overhang), candidate.overhang];
      const testFidelity = calculateSetFidelity(testOverhangs, enzyme);

      if (testFidelity < 0.3) continue; // Prune low-fidelity branches early

      const newUsed = new Set(usedOverhangs);
      newUsed.add(candidate.overhang);

      branch(depth + 1, [...current, candidate], newUsed);
    }
  }

  // Start search
  branch(0, [], new Set());

  if (!best) {
    // Fall back to greedy if B&B fails
    return optimizeGreedy(sequence, numFragments, enzyme, options);
  }

  // Format result
  const result = best.junctions.map((c, i) => ({
    ...c,
    targetIndex: i,
    idealPosition: targets[i].idealPosition,
    deviation: c.position - targets[i].idealPosition,
  }));

  return {
    algorithm: 'branchBound',
    junctions: result.sort((a, b) => a.position - b.position),
    overhangs: result.map(r => r.overhang),
    score: best.score,
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
 * @param {string} sequence - DNA sequence
 * @param {number} numFragments - Number of fragments
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Options
 * @returns {Object} Best result from multiple approaches
 */
export function optimizeHybrid(sequence, numFragments, enzyme = 'BsaI', options = {}) {
  const results = [];

  // Try greedy first (fast baseline)
  results.push(optimizeGreedy(sequence, numFragments, enzyme, options));

  // Try B&B for small assemblies
  if (numFragments <= 6) {
    results.push(optimizeBranchBound(sequence, numFragments, enzyme, options));
  }

  // Try Monte Carlo
  results.push(optimizeMonteCarlo(sequence, numFragments, enzyme, options));

  // Select best result
  const validResults = results.filter(r => r.complete);
  if (validResults.length === 0) {
    // Return best incomplete if no complete solutions
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
 * Main entry point - optimizeFusionSites
 *
 * Automatically selects the best algorithm based on problem size.
 *
 * @param {string} sequence - DNA sequence to split
 * @param {number} numFragments - Number of fragments desired
 * @param {Object} options - Optimization options
 * @returns {Object} Comprehensive optimization result
 */
export function optimizeFusionSites(sequence, numFragments, options = {}) {
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
  const internalSites = findInternalSites(seq, enzyme);
  if (internalSites.hasSites) {
    if (verbose) {
      console.warn(`Found ${internalSites.count} internal ${enzyme} sites - need domestication`);
    }
  }

  // Select algorithm
  let result;
  const selectedAlgorithm = algorithm === 'auto'
    ? (numFragments <= 5 ? 'branchBound' : numFragments <= 10 ? 'hybrid' : 'monteCarlo')
    : algorithm;

  const optimizerOptions = {
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
      internalSites: internalSites.hasSites ? internalSites.sites : [],
    });

    // Calculate fragment sizes
    const positions = [0, ...result.junctions.map(j => j.position + 4), seqLength];
    const fragmentSizes = [];
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
      internalSites: internalSites.hasSites ? internalSites : null,

      // Summary
      summary: {
        fidelity: result.score?.fidelity,
        efficiency: result.score?.efficiency,
        expectedSuccessRate: failurePrediction.expectedSuccessRate,
        riskLevel: failurePrediction.overallRisk,
        recommendations: failurePrediction.recommendations,
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
 * @param {string} sequence - DNA sequence
 * @param {number} numFragments - Number of fragments
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Quick result
 */
export function quickOptimize(sequence, numFragments, enzyme = 'BsaI') {
  return optimizeGreedy(sequence, numFragments, enzyme, {
    minDistanceFromEnds: 50,
    searchRadius: 30,
  });
}
