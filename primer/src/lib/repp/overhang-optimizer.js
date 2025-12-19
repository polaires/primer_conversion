/**
 * Overhang Set Optimizer for Golden Gate Assembly
 *
 * JavaScript port of NEB's ggtools_optimize.pl algorithm
 * Uses Monte Carlo simulated annealing to find optimal overhang sets
 *
 * Based on: Pryor et al. (2020) PLOS ONE
 * "Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized assembly design"
 * DOI: 10.1371/journal.pone.0238592
 *
 * Original Perl implementation: Copyright (C) 2023 New England Biolabs, Inc.
 */

import { reverseComplement } from './enzymes.js';
import ligationData from './ligation-data.json';

// Physical constants for Boltzmann distribution
const BOLTZMANN_K = 1.38064852e-23;
const AVOGADRO_N = 6.02214085e-23;

/**
 * Get all valid overhangs for an enzyme (excluding palindromes)
 * @param {string} enzyme - Enzyme name
 * @returns {string[]} Array of valid overhangs
 */
function getAllOverhangs(enzyme) {
  const enzymeData = getEnzymeData(enzyme);
  if (!enzymeData) return [];

  return enzymeData.overhangs.filter(oh => {
    const rc = reverseComplement(oh);
    // Exclude palindromic overhangs (can self-ligate)
    return oh !== rc;
  });
}

/**
 * Get enzyme data from ligation data
 * @param {string} enzyme - Enzyme name (BsaI, BbsI, BsmBI, Esp3I, SapI)
 * @returns {Object|null} Enzyme data
 */
function getEnzymeData(enzyme) {
  const keyMap = {
    'BsaI': 'BsaI-HFv2',
    'BsaI-HFv2': 'BsaI-HFv2',
    'BbsI': 'BbsI-HF',
    'BbsI-HF': 'BbsI-HF',
    'BsmBI': 'BsmBI-v2',
    'BsmBI-v2': 'BsmBI-v2',
    'Esp3I': 'Esp3I',
    'SapI': 'SapI',
  };

  const dataKey = keyMap[enzyme];
  return dataKey ? ligationData.enzymes[dataKey] : null;
}

/**
 * Calculate GC count in an overhang
 * @param {string} overhang - Overhang sequence
 * @returns {number} Number of G or C bases
 */
function countGC(overhang) {
  return (overhang.match(/[GC]/gi) || []).length;
}

/**
 * Calculate AT count in an overhang
 * @param {string} overhang - Overhang sequence
 * @returns {number} Number of A or T bases
 */
function countAT(overhang) {
  return (overhang.match(/[AT]/gi) || []).length;
}

/**
 * Filter overhangs based on constraints
 * @param {string[]} overhangs - Array of overhangs to filter
 * @param {Object} constraints - Constraint options
 * @param {Object} matrix - Ligation matrix
 * @returns {string[]} Filtered overhangs
 */
function filterOverhangs(overhangs, constraints, matrix) {
  const {
    excluded = [],
    maxGC = -1,
    maxAT = -1,
    minLigationEfficiency = -1,
  } = constraints;

  const excludeSet = new Set(excluded.map(o => o.toUpperCase()));
  // Also exclude reverse complements of excluded overhangs
  excluded.forEach(o => excludeSet.add(reverseComplement(o.toUpperCase())));

  return overhangs.filter(oh => {
    const ohUpper = oh.toUpperCase();
    const rc = reverseComplement(ohUpper);

    // Skip excluded overhangs
    if (excludeSet.has(ohUpper) || excludeSet.has(rc)) return false;

    // Skip palindromic overhangs
    if (ohUpper === rc) return false;

    // Check GC content
    if (maxGC !== -1 && countGC(ohUpper) > maxGC) return false;

    // Check AT content
    if (maxAT !== -1 && countAT(ohUpper) > maxAT) return false;

    // Check minimum ligation efficiency (self-ligation with complement)
    if (minLigationEfficiency !== -1) {
      const selfLigation = (matrix[ohUpper]?.[rc] || 0) + (matrix[rc]?.[ohUpper] || 0);
      if (selfLigation < minLigationEfficiency) return false;
    }

    return true;
  });
}

/**
 * Calculate ligation fidelity for a set of overhangs
 * This is the core fidelity calculation from the NEB algorithm
 *
 * For each junction, fidelity = correct_ligation / total_possible_ligations
 * Overall fidelity = product of all junction fidelities
 *
 * @param {string[]} overhangSet - Array of overhangs in the assembly
 * @param {Object} matrix - Ligation frequency matrix
 * @returns {number} Assembly fidelity (0-1)
 */
function calculateLigationFidelity(overhangSet, matrix) {
  const n = overhangSet.length;
  let fidelity = 1.0;

  for (let i = 0; i < n; i++) {
    const oh = overhangSet[i].toUpperCase();
    const rc = reverseComplement(oh);

    // Correct ligation: oh ligates with its Watson-Crick complement
    // Both directions: oh->rc and rc->oh
    const correct = (matrix[oh]?.[rc] || 0) + (matrix[rc]?.[oh] || 0);

    // Total: sum all possible ligations with overhangs in this set
    let total = 0;

    for (let j = 0; j < n; j++) {
      const oh2 = overhangSet[j].toUpperCase();
      const rc2 = reverseComplement(oh2);

      // All four possible interactions between oh/rc and oh2/rc2
      total += matrix[oh]?.[oh2] || 0;
      total += matrix[oh]?.[rc2] || 0;
      total += matrix[rc]?.[oh2] || 0;
      total += matrix[rc]?.[rc2] || 0;
    }

    if (total > 0) {
      fidelity *= correct / total;
    }
  }

  return fidelity;
}

/**
 * Monte Carlo optimization step
 * Attempts to improve the current solution by random mutation
 *
 * @param {string[]} current - Current overhang set
 * @param {string[]} pool - Pool of available overhangs
 * @param {Object} matrix - Ligation matrix
 * @param {number} temperature - Current annealing temperature
 * @param {Set} fixedIndices - Indices of fixed overhangs that cannot be changed
 * @returns {Object} {accepted, newSet, newFidelity}
 */
function mcStep(current, pool, matrix, temperature, fixedIndices = new Set()) {
  // Track which overhangs are in the current solution (to avoid duplicates)
  const currentSet = new Set();
  current.forEach(oh => {
    currentSet.add(oh.toUpperCase());
    currentSet.add(reverseComplement(oh.toUpperCase()));
  });

  // Find variable positions (not fixed)
  const variableIndices = [];
  for (let i = 0; i < current.length; i++) {
    if (!fixedIndices.has(i)) {
      variableIndices.push(i);
    }
  }

  if (variableIndices.length === 0) {
    // Nothing to optimize
    return { accepted: false, newSet: current, newFidelity: calculateLigationFidelity(current, matrix) };
  }

  // Pick a random variable position
  const pos = variableIndices[Math.floor(Math.random() * variableIndices.length)];

  // Pick a random overhang from the pool that's not already in the solution
  const availableOverhangs = pool.filter(oh => {
    const ohUpper = oh.toUpperCase();
    const rc = reverseComplement(ohUpper);
    return !currentSet.has(ohUpper) && !currentSet.has(rc);
  });

  if (availableOverhangs.length === 0) {
    return { accepted: false, newSet: current, newFidelity: calculateLigationFidelity(current, matrix) };
  }

  const newOH = availableOverhangs[Math.floor(Math.random() * availableOverhangs.length)];

  // Create new solution
  const newSet = [...current];
  const oldOH = newSet[pos];
  newSet[pos] = newOH;

  // Calculate fidelities
  const oldFidelity = calculateLigationFidelity(current, matrix);
  const newFidelity = calculateLigationFidelity(newSet, matrix);

  // Metropolis criterion
  const delta = newFidelity - oldFidelity;

  if (delta > 0) {
    // Better solution - always accept
    return { accepted: true, newSet, newFidelity };
  } else {
    // Worse solution - accept with probability based on Boltzmann distribution
    const T = Math.pow(2, temperature);
    const prob = Math.exp(delta / (BOLTZMANN_K * T / AVOGADRO_N));

    if (Math.random() < prob) {
      return { accepted: true, newSet, newFidelity };
    } else {
      return { accepted: false, newSet: current, newFidelity: oldFidelity };
    }
  }
}

/**
 * Find appropriate annealing temperature for target acceptance ratio
 * @param {number} targetAR - Target acceptance ratio (default 0.05)
 * @param {string[]} pool - Pool of available overhangs
 * @param {number} size - Number of overhangs needed
 * @param {Object} matrix - Ligation matrix
 * @param {number} iterations - Test iterations per temperature
 * @returns {number} Appropriate temperature exponent
 */
function findAnnealingTemperature(targetAR, pool, size, matrix, iterations = 1000) {
  let exp = 0;

  // Generate initial random solution
  const shuffled = [...pool].sort(() => Math.random() - 0.5);
  let current = shuffled.slice(0, size);

  // Test acceptance ratio at current temperature
  const testAR = (temperature) => {
    let accepted = 0;
    let testSet = [...current];

    for (let i = 0; i < iterations; i++) {
      const result = mcStep(testSet, pool, matrix, temperature);
      if (result.accepted) {
        accepted++;
        testSet = result.newSet;
      }
    }

    return accepted / iterations;
  };

  let ar = testAR(exp);

  // Adjust temperature to reach target acceptance ratio
  if (ar > targetAR) {
    // Lower temperature to reduce acceptance
    let counter = 0;
    while (ar > targetAR && counter < 100) {
      exp--;
      ar = testAR(exp);
      counter++;
    }
  } else if (ar < targetAR) {
    // Raise temperature to increase acceptance
    let counter = 0;
    while (ar < targetAR && counter < 100) {
      exp++;
      ar = testAR(exp);
      counter++;
    }
  }

  return exp;
}

/**
 * Main Monte Carlo optimization algorithm
 * @param {Object} options - Optimization options
 * @returns {Object} Optimization result
 */
function runMonteCarlo(options) {
  const {
    pool,
    size,
    matrix,
    temperature,
    iterations = 1000,
    initialSet = null,
    fixedIndices = new Set(),
    onProgress = null,
  } = options;

  // Initialize solution
  let current;
  if (initialSet && initialSet.length === size) {
    current = [...initialSet];
  } else {
    // Random initialization
    const shuffled = [...pool].sort(() => Math.random() - 0.5);
    current = shuffled.slice(0, size);
  }

  let bestSet = [...current];
  let bestFidelity = calculateLigationFidelity(current, matrix);
  let currentFidelity = bestFidelity;

  let accepted = 0;
  let rejected = 0;

  for (let i = 0; i < iterations; i++) {
    const result = mcStep(current, pool, matrix, temperature, fixedIndices);

    if (result.accepted) {
      accepted++;
      current = result.newSet;
      currentFidelity = result.newFidelity;

      if (currentFidelity > bestFidelity) {
        bestSet = [...current];
        bestFidelity = currentFidelity;
      }
    } else {
      rejected++;
    }

    // Progress callback
    if (onProgress && i % 100 === 0) {
      onProgress({
        iteration: i,
        totalIterations: iterations,
        currentFidelity,
        bestFidelity,
        acceptanceRatio: accepted / (i + 1),
      });
    }
  }

  return {
    overhangs: bestSet,
    fidelity: bestFidelity,
    acceptanceRatio: accepted / iterations,
    iterations,
  };
}

/**
 * Optimize an overhang set for Golden Gate assembly
 *
 * This is the main entry point for overhang optimization.
 * Uses Monte Carlo simulated annealing to find an optimal set of overhangs
 * that maximizes assembly fidelity while respecting constraints.
 *
 * @param {Object} options - Optimization options
 * @param {number} options.numJunctions - Number of junctions (overhangs) needed
 * @param {string} [options.enzyme='BsaI'] - Enzyme to use
 * @param {string[]} [options.required=[]] - Overhangs that must be included
 * @param {string[]} [options.excluded=[]] - Overhangs that must not be used
 * @param {number} [options.maxGC=-1] - Maximum GC bases per overhang (-1 = no limit)
 * @param {number} [options.maxAT=-1] - Maximum AT bases per overhang (-1 = no limit)
 * @param {number} [options.minLigationEfficiency=-1] - Minimum self-ligation frequency (-1 = no limit)
 * @param {number} [options.iterations=10000] - Number of MC iterations
 * @param {number} [options.targetAcceptanceRatio=0.05] - Target acceptance ratio for annealing
 * @param {Function} [options.onProgress=null] - Progress callback
 * @returns {Object} Optimization result
 *
 * @example
 * // Basic optimization for 10 junctions
 * const result = optimizeOverhangSet({ numJunctions: 10 });
 *
 * @example
 * // With constraints
 * const result = optimizeOverhangSet({
 *   numJunctions: 15,
 *   required: ['GGAG', 'AATG', 'GCTT'],  // MoClo compatibility
 *   excluded: ['GATC', 'ATAT'],           // Avoid these
 *   maxGC: 3,                             // Max 3 GC per overhang
 *   enzyme: 'BsaI',
 * });
 */
export function optimizeOverhangSet(options) {
  const {
    numJunctions,
    enzyme = 'BsaI',
    required = [],
    excluded = [],
    maxGC = -1,
    maxAT = -1,
    minLigationEfficiency = -1,
    iterations = 10000,
    targetAcceptanceRatio = 0.05,
    onProgress = null,
  } = options;

  // Validate inputs
  if (!numJunctions || numJunctions < 2) {
    throw new Error('numJunctions must be at least 2');
  }

  if (numJunctions > 50) {
    throw new Error('numJunctions cannot exceed 50 (practical limit)');
  }

  // Get enzyme data
  const enzymeData = getEnzymeData(enzyme);
  if (!enzymeData) {
    throw new Error(`Unknown enzyme: ${enzyme}. Supported: BsaI, BbsI, BsmBI, Esp3I, SapI`);
  }

  const matrix = enzymeData.matrix;

  // Get all valid overhangs and apply constraints
  const allOverhangs = getAllOverhangs(enzyme);
  const filteredPool = filterOverhangs(allOverhangs, {
    excluded,
    maxGC,
    maxAT,
    minLigationEfficiency,
  }, matrix);

  // Validate required overhangs
  const requiredUpper = required.map(oh => oh.toUpperCase());
  for (const oh of requiredUpper) {
    if (!filteredPool.includes(oh)) {
      // Check if it's excluded by constraints
      if (excluded.includes(oh)) {
        throw new Error(`Required overhang ${oh} is also in excluded list`);
      }
      throw new Error(`Required overhang ${oh} is not valid or doesn't meet constraints`);
    }
  }

  // Check if we have enough overhangs
  if (filteredPool.length < numJunctions) {
    throw new Error(
      `Not enough valid overhangs (${filteredPool.length}) to create ${numJunctions} junctions. ` +
      `Try relaxing constraints (maxGC, maxAT, excluded).`
    );
  }

  // Build initial set with required overhangs
  const fixedIndices = new Set();
  const initialSet = [...requiredUpper];

  for (let i = 0; i < requiredUpper.length; i++) {
    fixedIndices.add(i);
  }

  // Fill remaining slots with random overhangs
  const remainingNeeded = numJunctions - initialSet.length;
  const usedSet = new Set(initialSet.map(oh => oh.toUpperCase()));
  initialSet.forEach(oh => usedSet.add(reverseComplement(oh.toUpperCase())));

  const availableForFill = filteredPool.filter(oh => !usedSet.has(oh.toUpperCase()));
  const shuffled = [...availableForFill].sort(() => Math.random() - 0.5);

  for (let i = 0; i < remainingNeeded && i < shuffled.length; i++) {
    initialSet.push(shuffled[i]);
  }

  if (initialSet.length < numJunctions) {
    throw new Error(
      `Could not create initial set of ${numJunctions} overhangs. ` +
      `Only ${initialSet.length} unique overhangs available after applying constraints.`
    );
  }

  // Find appropriate temperature
  const temperature = findAnnealingTemperature(
    targetAcceptanceRatio,
    filteredPool,
    numJunctions,
    matrix,
    Math.min(iterations / 10, 1000)
  );

  // Run Monte Carlo optimization
  const result = runMonteCarlo({
    pool: filteredPool,
    size: numJunctions,
    matrix,
    temperature,
    iterations,
    initialSet,
    fixedIndices,
    onProgress,
  });

  // Calculate detailed fidelity breakdown
  const junctionDetails = result.overhangs.map((oh, index) => {
    const ohUpper = oh.toUpperCase();
    const rc = reverseComplement(ohUpper);
    const correct = (matrix[ohUpper]?.[rc] || 0) + (matrix[rc]?.[ohUpper] || 0);

    let total = 0;
    for (const other of result.overhangs) {
      const oh2 = other.toUpperCase();
      const rc2 = reverseComplement(oh2);
      total += matrix[ohUpper]?.[oh2] || 0;
      total += matrix[ohUpper]?.[rc2] || 0;
      total += matrix[rc]?.[oh2] || 0;
      total += matrix[rc]?.[rc2] || 0;
    }

    const junctionFidelity = total > 0 ? correct / total : 0;

    return {
      index,
      overhang: ohUpper,
      reverseComplement: rc,
      fidelity: junctionFidelity,
      fidelityPercent: (junctionFidelity * 100).toFixed(1) + '%',
      isRequired: requiredUpper.includes(ohUpper),
      gcContent: countGC(ohUpper),
      atContent: countAT(ohUpper),
    };
  });

  // Sort by fidelity to identify weak points
  const sortedByFidelity = [...junctionDetails].sort((a, b) => a.fidelity - b.fidelity);

  return {
    enzyme,
    numJunctions,
    overhangs: result.overhangs.map(oh => oh.toUpperCase()),
    fidelity: result.fidelity,
    fidelityPercent: (result.fidelity * 100).toFixed(2) + '%',
    junctions: junctionDetails,
    weakestJunction: sortedByFidelity[0],
    strongestJunction: sortedByFidelity[sortedByFidelity.length - 1],
    constraints: {
      required: requiredUpper,
      excluded: excluded.map(oh => oh.toUpperCase()),
      maxGC,
      maxAT,
      minLigationEfficiency,
    },
    optimization: {
      iterations,
      acceptanceRatio: result.acceptanceRatio,
      temperature,
      poolSize: filteredPool.length,
    },
    source: 'optimized',
    algorithm: 'Monte Carlo simulated annealing',
    dataSource: ligationData.metadata,
  };
}

/**
 * Evaluate fidelity of a given overhang set without optimization
 * Useful for comparing custom sets against optimized ones
 *
 * @param {string[]} overhangs - Array of overhangs to evaluate
 * @param {string} [enzyme='BsaI'] - Enzyme to use
 * @returns {Object} Fidelity analysis
 */
export function evaluateOverhangSet(overhangs, enzyme = 'BsaI') {
  const enzymeData = getEnzymeData(enzyme);
  if (!enzymeData) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  const matrix = enzymeData.matrix;
  const fidelity = calculateLigationFidelity(overhangs, matrix);

  const junctionDetails = overhangs.map((oh, index) => {
    const ohUpper = oh.toUpperCase();
    const rc = reverseComplement(ohUpper);
    const correct = (matrix[ohUpper]?.[rc] || 0) + (matrix[rc]?.[ohUpper] || 0);

    let total = 0;
    for (const other of overhangs) {
      const oh2 = other.toUpperCase();
      const rc2 = reverseComplement(oh2);
      total += matrix[ohUpper]?.[oh2] || 0;
      total += matrix[ohUpper]?.[rc2] || 0;
      total += matrix[rc]?.[oh2] || 0;
      total += matrix[rc]?.[rc2] || 0;
    }

    const junctionFidelity = total > 0 ? correct / total : 0;

    return {
      index,
      overhang: ohUpper,
      reverseComplement: rc,
      fidelity: junctionFidelity,
      fidelityPercent: (junctionFidelity * 100).toFixed(1) + '%',
      gcContent: countGC(ohUpper),
      atContent: countAT(ohUpper),
    };
  });

  const sortedByFidelity = [...junctionDetails].sort((a, b) => a.fidelity - b.fidelity);

  return {
    enzyme,
    numJunctions: overhangs.length,
    overhangs: overhangs.map(oh => oh.toUpperCase()),
    fidelity,
    fidelityPercent: (fidelity * 100).toFixed(2) + '%',
    junctions: junctionDetails,
    weakestJunction: sortedByFidelity[0],
    strongestJunction: sortedByFidelity[sortedByFidelity.length - 1],
    source: 'evaluated',
  };
}

/**
 * Generate multiple optimized sets and return the best one
 * Useful for finding globally optimal solutions
 *
 * @param {Object} options - Same options as optimizeOverhangSet
 * @param {number} [options.runs=5] - Number of optimization runs
 * @returns {Object} Best optimization result
 */
export function optimizeOverhangSetMultiRun(options) {
  const { runs = 5, ...optimizeOptions } = options;

  let bestResult = null;
  const allResults = [];

  for (let i = 0; i < runs; i++) {
    const result = optimizeOverhangSet(optimizeOptions);
    allResults.push(result);

    if (!bestResult || result.fidelity > bestResult.fidelity) {
      bestResult = result;
    }
  }

  return {
    ...bestResult,
    multiRun: {
      totalRuns: runs,
      allFidelities: allResults.map(r => r.fidelity),
      averageFidelity: allResults.reduce((sum, r) => sum + r.fidelity, 0) / runs,
      bestRunIndex: allResults.indexOf(bestResult),
    },
  };
}

/**
 * Batch score multiple random overhang sets
 * Useful for understanding the fidelity distribution
 *
 * @param {Object} options - Scoring options
 * @param {number} options.numJunctions - Number of junctions
 * @param {string} [options.enzyme='BsaI'] - Enzyme to use
 * @param {number} [options.numSets=100] - Number of random sets to generate
 * @param {string[]} [options.excluded=[]] - Overhangs to exclude
 * @returns {Object} Batch scoring results
 */
export function batchScoreRandomSets(options) {
  const {
    numJunctions,
    enzyme = 'BsaI',
    numSets = 100,
    excluded = [],
  } = options;

  const enzymeData = getEnzymeData(enzyme);
  if (!enzymeData) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  const matrix = enzymeData.matrix;
  const allOverhangs = getAllOverhangs(enzyme);
  const pool = filterOverhangs(allOverhangs, { excluded }, matrix);

  const results = [];

  for (let i = 0; i < numSets; i++) {
    // Generate random set
    const shuffled = [...pool].sort(() => Math.random() - 0.5);
    const randomSet = shuffled.slice(0, numJunctions);

    // Calculate fidelity
    const fidelity = calculateLigationFidelity(randomSet, matrix);
    results.push({ overhangs: randomSet, fidelity });
  }

  // Sort by fidelity
  results.sort((a, b) => b.fidelity - a.fidelity);

  const fidelities = results.map(r => r.fidelity);
  const mean = fidelities.reduce((a, b) => a + b, 0) / fidelities.length;
  const variance = fidelities.reduce((sum, f) => sum + Math.pow(f - mean, 2), 0) / fidelities.length;
  const stdDev = Math.sqrt(variance);

  return {
    enzyme,
    numJunctions,
    numSets,
    best: results[0],
    worst: results[results.length - 1],
    statistics: {
      mean,
      stdDev,
      min: Math.min(...fidelities),
      max: Math.max(...fidelities),
      median: fidelities[Math.floor(fidelities.length / 2)],
    },
    percentiles: {
      p90: fidelities[Math.floor(fidelities.length * 0.1)],
      p75: fidelities[Math.floor(fidelities.length * 0.25)],
      p50: fidelities[Math.floor(fidelities.length * 0.5)],
      p25: fidelities[Math.floor(fidelities.length * 0.75)],
      p10: fidelities[Math.floor(fidelities.length * 0.9)],
    },
    allResults: results,
  };
}

// Export helper functions for advanced usage
export {
  calculateLigationFidelity,
  getAllOverhangs,
  filterOverhangs,
  getEnzymeData,
};
