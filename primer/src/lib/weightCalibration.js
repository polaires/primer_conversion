/**
 * Weight Calibration Module
 *
 * Implements cross-validation and grid search for optimizing scoring weights
 * based on validation datasets of primer pairs with known success/failure outcomes.
 *
 * Based on Section 2.5 of SCORING_MANUSCRIPT_STRATEGY.md
 */

import { calculateCompositeScore, classifyQuality } from './scoring.js';

/**
 * Calibrated weight configuration (Dec 2025)
 *
 * Optimized via grid search + coordinate descent on 829-entry Döring dataset
 * Performance: F1=88.4%, AUC=0.929, Cross-validated: F1=88.4% ± 2.2%
 *
 * Key findings from calibration:
 * - offTarget is DOMINANT factor (0.10 → 0.25) - confirms GM1 model
 * - terminal3DG much more important (0.05 → 0.20)
 * - G-Quadruplex validated: GGGG primers have 37.2% success vs 44% overall
 * - gQuadruplexRev > gQuadruplexFwd (0.11 vs 0.05) - reverse primer G4 more critical
 * - heterodimer less important than expected (0.10 → 0.04)
 */
export const DEFAULT_WEIGHTS = {
  // Tier 1: Critical - offTarget and terminal stability dominate
  // Optimized via grid search + coordinate descent on 829-entry Döring dataset
  // Performance: F1=88.7% ± 1.3%, AUC=0.930 ± 0.017 (5-fold CV)
  // Note: calculateCompositeScore normalizes internally, so total sum > 1 is OK
  offTarget: 0.25,       // Dominant factor - most discriminative (+0.515 diff)
  terminal3DG: 0.20,     // Very important for 3' end stability (+0.487 diff)
  gQuadruplexRev: 0.15,  // G4 in reverse primer is critical (GGGG: 37.2% vs 44% overall)

  // Tier 2: Important
  gQuadruplexFwd: 0.05,  // G4 in forward primer
  tmRev: 0.05,
  hairpinRev: 0.05,
  heterodimer: 0.06,     // Cross-dimer formation
  gcRev: 0.04,
  selfDimerFwd: 0.04,    // Self-dimer formation (optimized from 0.02)
  selfDimerRev: 0.04,
  threePrimeCompFwd: 0.04,  // 3' end composition
  threePrimeCompRev: 0.04,  // 3' end composition

  // Tier 3: Minor
  gcFwd: 0.02,
  gcClampFwd: 0.03,
  gcClampRev: 0.03,
  // Literature (PrimerScore2, ML-PCR): Tm difference has "little effect" on PCR success
  tmDiff: 0.03,
  tmFwd: 0.02,
  hairpinFwd: 0.02,
  homopolymerFwd: 0.02,
  homopolymerRev: 0.02,
  ampliconLength: 0.02,
  ampliconStructure: 0.02,
  lengthFwd: 0.01,
  lengthRev: 0.01,
  terminalBaseFwd: 0.01,
  terminalBaseRev: 0.01,
  lengthDiff: 0.01,
  distanceToROI: 0.01,
};

/**
 * Assembly-optimized weights for isothermal assembly (Gibson/NEBuilder)
 *
 * Key differences from DEFAULT_WEIGHTS:
 * - Range-based approach: Tm/GC/length within range all score equally
 * - Higher weight on quality factors (GC clamp, terminal3DG, patterns)
 * - Lower weight on off-target (assembly uses specific fragments)
 * - Higher heterodimer weight (multi-fragment reactions)
 *
 * Designed to maximize quality factors when basic constraints are met.
 */
export const ASSEMBLY_WEIGHTS = {
  // Primary constraints - within range = full score
  tmFwd: 0.08,
  tmRev: 0.08,
  gcFwd: 0.04,
  gcRev: 0.04,
  lengthFwd: 0.03,
  lengthRev: 0.03,

  // Quality factors - these differentiate when constraints met
  terminal3DG: 0.18,     // Critical for extension initiation
  gcClampFwd: 0.08,      // GC clamp for reliable priming
  gcClampRev: 0.08,
  threePrimeCompFwd: 0.06,
  threePrimeCompRev: 0.06,

  // Secondary structure
  hairpinFwd: 0.06,      // Important for tailed primers
  hairpinRev: 0.06,
  selfDimerFwd: 0.03,
  selfDimerRev: 0.03,
  heterodimer: 0.10,     // Higher for multi-fragment assemblies

  // Pattern avoidance (captured in threePrimeComp, but explicit)
  homopolymerFwd: 0.03,
  homopolymerRev: 0.03,

  // Less critical for assembly (specific fragments used)
  offTarget: 0.05,
  gQuadruplexFwd: 0.04,
  gQuadruplexRev: 0.04,

  // Pair matching
  tmDiff: 0.05,          // Important for consistent annealing
};

/**
 * Overlap-specific weights for junction optimization
 *
 * Focused on overlap sequence quality rather than primer pair dynamics.
 * Used when optimizing overlap regions between fragments.
 *
 * NOTE: heterodimer removed - overlaps don't form heterodimers with themselves.
 * Added selfComplementarity to detect palindromic sequences that cause mispriming.
 */
export const OVERLAP_WEIGHTS = {
  // Range compliance (flat score within optimal range)
  tmInRange: 0.16,
  gcInRange: 0.10,
  lengthInRange: 0.07,

  // Quality differentiation
  gcClamp: 0.14,            // 3' end G/C for stable priming
  terminal3DG: 0.14,        // 3' binding stability
  hairpin: 0.12,            // Secondary structure
  patternAvoidance: 0.12,   // Poly runs, repeats
  selfComplementarity: 0.08, // Palindromes cause self-annealing
  balancedGC: 0.05,         // Prefer 45-55% over range edges
  crossJunctionHairpin: 0.02, // Hairpin spanning junction (checked separately)
};

/**
 * Single annealing region weights (used before pair evaluation)
 *
 * These weights prioritize quality factors when Tm/GC/length are within range.
 * Moved from inline definition in assemblyCore.js for centralized management.
 */
export const ANNEALING_SINGLE_WEIGHTS = {
  // Range compliance
  tmInRange: 0.15,
  gcInRange: 0.08,
  lengthInRange: 0.07,

  // 3' end quality (critical for extension initiation)
  gcClamp: 0.16,
  terminal3DG: 0.18,

  // Secondary structure and patterns
  hairpin: 0.12,
  homopolymer: 0.06,
  threePrimeComp: 0.06,  // Reduced - overlaps with gcClamp + terminal3DG
  gQuadruplex: 0.08,     // Q5 struggles with G4
  selfComplementarity: 0.04,
};

/**
 * Weight bounds for grid search optimization
 */
export const WEIGHT_BOUNDS = {
  // Critical features: can vary more
  gQuadruplexFwd: { min: 0.05, max: 0.15, step: 0.02 },
  gQuadruplexRev: { min: 0.05, max: 0.15, step: 0.02 },
  offTarget: { min: 0.10, max: 0.30, step: 0.05 },
  terminal3DG: { min: 0.05, max: 0.20, step: 0.05 },

  // Important features: moderate variation
  tmFwd: { min: 0.02, max: 0.10, step: 0.02 },
  tmRev: { min: 0.02, max: 0.10, step: 0.02 },
  gcFwd: { min: 0.02, max: 0.08, step: 0.02 },
  gcRev: { min: 0.02, max: 0.08, step: 0.02 },
  hairpinFwd: { min: 0.02, max: 0.10, step: 0.02 },
  hairpinRev: { min: 0.02, max: 0.10, step: 0.02 },
  selfDimerFwd: { min: 0.02, max: 0.08, step: 0.02 },
  selfDimerRev: { min: 0.02, max: 0.08, step: 0.02 },
  heterodimer: { min: 0.02, max: 0.10, step: 0.02 },

  // Minor features: narrow variation
  tmDiff: { min: 0.01, max: 0.05, step: 0.01 },
  gcClampFwd: { min: 0.01, max: 0.05, step: 0.01 },
  gcClampRev: { min: 0.01, max: 0.05, step: 0.01 },
};

/**
 * Calculate performance metrics for predictions
 *
 * @param {Array} predictions - Array of {predicted: boolean, actual: boolean}
 * @returns {Object} Metrics including accuracy, precision, recall, F1, AUC
 */
export function calculateMetrics(predictions) {
  let tp = 0, tn = 0, fp = 0, fn = 0;

  for (const { predicted, actual } of predictions) {
    if (predicted && actual) tp++;
    else if (!predicted && !actual) tn++;
    else if (predicted && !actual) fp++;
    else fn++;
  }

  const accuracy = (tp + tn) / (tp + tn + fp + fn) || 0;
  const precision = tp / (tp + fp) || 0;
  const recall = tp / (tp + fn) || 0;
  const f1 = 2 * (precision * recall) / (precision + recall) || 0;
  const specificity = tn / (tn + fp) || 0;

  return {
    tp, tn, fp, fn,
    accuracy: Math.round(accuracy * 1000) / 1000,
    precision: Math.round(precision * 1000) / 1000,
    recall: Math.round(recall * 1000) / 1000,
    f1: Math.round(f1 * 1000) / 1000,
    specificity: Math.round(specificity * 1000) / 1000,
  };
}

/**
 * Calculate AUC-ROC from scored predictions
 *
 * @param {Array} scoredPredictions - Array of {score: number, actual: boolean}
 * @returns {number} AUC-ROC value (0-1)
 */
export function calculateAUC(scoredPredictions) {
  // Sort by score descending
  const sorted = [...scoredPredictions].sort((a, b) => b.score - a.score);

  let positives = 0;
  let negatives = 0;
  for (const { actual } of sorted) {
    if (actual) positives++;
    else negatives++;
  }

  if (positives === 0 || negatives === 0) return 0.5;

  // Calculate AUC using trapezoidal rule
  let auc = 0;
  let tpCount = 0;
  let fpCount = 0;
  let prevTPR = 0;
  let prevFPR = 0;

  for (const { actual } of sorted) {
    if (actual) tpCount++;
    else fpCount++;

    const tpr = tpCount / positives;
    const fpr = fpCount / negatives;

    // Trapezoidal integration
    auc += (fpr - prevFPR) * (tpr + prevTPR) / 2;

    prevTPR = tpr;
    prevFPR = fpr;
  }

  return Math.round(auc * 1000) / 1000;
}

/**
 * Score a primer pair using given weights
 *
 * @param {Object} scores - Individual feature scores for the primer pair
 * @param {Object} weights - Weight configuration
 * @returns {number} Composite score (0-100)
 */
export function scorePrimerPair(scores, weights) {
  const result = calculateCompositeScore(scores, weights);
  return result.score;
}

/**
 * Predict success/failure based on score threshold
 *
 * @param {number} score - Composite score (0-100)
 * @param {number} threshold - Score threshold for success prediction
 * @returns {boolean} Predicted success
 */
export function predictSuccess(score, threshold = 60) {
  return score >= threshold;
}

/**
 * Evaluate weights on a validation dataset
 *
 * @param {Array} dataset - Array of {scores: Object, actual: boolean} or {scores: Object, success: boolean}
 * @param {Object} weights - Weight configuration
 * @param {number} threshold - Score threshold for success prediction
 * @returns {Object} Evaluation metrics
 */
export function evaluateWeights(dataset, weights, threshold = 60) {
  const predictions = [];
  const scoredPredictions = [];

  for (const entry of dataset) {
    // Support both 'actual' and 'success' field names
    const actualSuccess = entry.actual !== undefined ? entry.actual : entry.success;
    const score = scorePrimerPair(entry.scores, weights);
    const predicted = predictSuccess(score, threshold);

    predictions.push({ predicted, actual: actualSuccess });
    scoredPredictions.push({ score, actual: actualSuccess });
  }

  const metrics = calculateMetrics(predictions);
  const auc = calculateAUC(scoredPredictions);

  return {
    ...metrics,
    auc,
    threshold,
  };
}

/**
 * Find optimal threshold for given weights
 *
 * @param {Array} dataset - Validation dataset
 * @param {Object} weights - Weight configuration
 * @returns {Object} Optimal threshold and metrics
 */
export function findOptimalThreshold(dataset, weights) {
  let bestThreshold = 60;
  let bestF1 = 0;
  let bestMetrics = null;

  // Test thresholds from 40 to 80
  for (let threshold = 40; threshold <= 80; threshold += 5) {
    const metrics = evaluateWeights(dataset, weights, threshold);

    if (metrics.f1 > bestF1) {
      bestF1 = metrics.f1;
      bestThreshold = threshold;
      bestMetrics = metrics;
    }
  }

  return {
    threshold: bestThreshold,
    metrics: bestMetrics,
  };
}

/**
 * K-fold cross-validation
 *
 * @param {Array} dataset - Full validation dataset
 * @param {Object} weights - Weight configuration
 * @param {number} k - Number of folds (default: 5)
 * @param {number} threshold - Score threshold
 * @returns {Object} Cross-validation results
 */
export function crossValidate(dataset, weights, k = 5, threshold = 60) {
  // Shuffle dataset
  const shuffled = [...dataset].sort(() => Math.random() - 0.5);
  const foldSize = Math.floor(shuffled.length / k);

  const foldMetrics = [];

  for (let i = 0; i < k; i++) {
    // Split into train and test
    const testStart = i * foldSize;
    const testEnd = i === k - 1 ? shuffled.length : (i + 1) * foldSize;
    const testSet = shuffled.slice(testStart, testEnd);

    // Evaluate on test fold
    const metrics = evaluateWeights(testSet, weights, threshold);
    foldMetrics.push(metrics);
  }

  // Average metrics across folds
  const avgMetrics = {
    accuracy: 0,
    precision: 0,
    recall: 0,
    f1: 0,
    auc: 0,
  };

  for (const metrics of foldMetrics) {
    avgMetrics.accuracy += metrics.accuracy;
    avgMetrics.precision += metrics.precision;
    avgMetrics.recall += metrics.recall;
    avgMetrics.f1 += metrics.f1;
    avgMetrics.auc += metrics.auc;
  }

  for (const key of Object.keys(avgMetrics)) {
    avgMetrics[key] = Math.round((avgMetrics[key] / k) * 1000) / 1000;
  }

  // Calculate standard deviation
  const stdMetrics = {
    accuracy: 0,
    precision: 0,
    recall: 0,
    f1: 0,
    auc: 0,
  };

  for (const metrics of foldMetrics) {
    stdMetrics.accuracy += Math.pow(metrics.accuracy - avgMetrics.accuracy, 2);
    stdMetrics.precision += Math.pow(metrics.precision - avgMetrics.precision, 2);
    stdMetrics.recall += Math.pow(metrics.recall - avgMetrics.recall, 2);
    stdMetrics.f1 += Math.pow(metrics.f1 - avgMetrics.f1, 2);
    stdMetrics.auc += Math.pow(metrics.auc - avgMetrics.auc, 2);
  }

  for (const key of Object.keys(stdMetrics)) {
    stdMetrics[key] = Math.round(Math.sqrt(stdMetrics[key] / k) * 1000) / 1000;
  }

  return {
    folds: k,
    mean: avgMetrics,
    std: stdMetrics,
    foldResults: foldMetrics,
  };
}

/**
 * Generate weight combinations for grid search
 *
 * @param {Array} features - Features to optimize
 * @param {Object} bounds - Weight bounds configuration
 * @returns {Array} Array of weight configurations
 */
export function generateWeightGrid(features, bounds = WEIGHT_BOUNDS) {
  const combinations = [{}];

  for (const feature of features) {
    const featureBounds = bounds[feature];
    if (!featureBounds) continue;

    const { min, max, step } = featureBounds;
    const newCombinations = [];

    for (const combo of combinations) {
      for (let value = min; value <= max; value += step) {
        newCombinations.push({
          ...combo,
          [feature]: Math.round(value * 100) / 100,
        });
      }
    }

    combinations.length = 0;
    combinations.push(...newCombinations);
  }

  return combinations;
}

/**
 * Grid search for optimal weights
 *
 * @param {Array} dataset - Validation dataset
 * @param {Array} featuresToOptimize - Features to include in grid search
 * @param {Object} baseWeights - Base weights for non-optimized features
 * @param {Object} options - Search options
 * @returns {Object} Best weights and metrics
 */
export function gridSearch(dataset, featuresToOptimize, baseWeights = DEFAULT_WEIGHTS, options = {}) {
  const {
    threshold = 60,
    metric = 'f1',  // Metric to optimize: 'f1', 'accuracy', 'auc'
    verbose = false,
  } = options;

  const weightCombinations = generateWeightGrid(featuresToOptimize);

  if (verbose) {
    console.log(`Testing ${weightCombinations.length} weight combinations...`);
  }

  let bestWeights = baseWeights;
  let bestScore = 0;
  let bestMetrics = null;

  for (const partialWeights of weightCombinations) {
    const testWeights = { ...baseWeights, ...partialWeights };
    const metrics = evaluateWeights(dataset, testWeights, threshold);

    if (metrics[metric] > bestScore) {
      bestScore = metrics[metric];
      bestWeights = testWeights;
      bestMetrics = metrics;
    }
  }

  return {
    weights: bestWeights,
    metrics: bestMetrics,
    optimizedMetric: metric,
    combinationsTested: weightCombinations.length,
  };
}

/**
 * Iterative optimization using coordinate descent
 *
 * @param {Array} dataset - Validation dataset
 * @param {Object} initialWeights - Starting weights
 * @param {Object} options - Optimization options
 * @returns {Object} Optimized weights and metrics
 */
export function coordinateDescent(dataset, initialWeights = DEFAULT_WEIGHTS, options = {}) {
  const {
    maxIterations = 10,
    threshold = 60,
    metric = 'f1',
    stepSize = 0.02,
    minImprovement = 0.001,
    features = Object.keys(WEIGHT_BOUNDS),
    verbose = false,
  } = options;

  let currentWeights = { ...initialWeights };
  let currentScore = evaluateWeights(dataset, currentWeights, threshold)[metric];

  for (let iter = 0; iter < maxIterations; iter++) {
    let improved = false;

    for (const feature of features) {
      const bounds = WEIGHT_BOUNDS[feature];
      if (!bounds) continue;

      // Try increasing weight
      const increasedWeights = {
        ...currentWeights,
        [feature]: Math.min(bounds.max, currentWeights[feature] + stepSize),
      };
      const increasedScore = evaluateWeights(dataset, increasedWeights, threshold)[metric];

      if (increasedScore > currentScore + minImprovement) {
        currentWeights = increasedWeights;
        currentScore = increasedScore;
        improved = true;
        if (verbose) console.log(`Iter ${iter}: ${feature} ↑ → ${metric}=${currentScore}`);
        continue;
      }

      // Try decreasing weight
      const decreasedWeights = {
        ...currentWeights,
        [feature]: Math.max(bounds.min, currentWeights[feature] - stepSize),
      };
      const decreasedScore = evaluateWeights(dataset, decreasedWeights, threshold)[metric];

      if (decreasedScore > currentScore + minImprovement) {
        currentWeights = decreasedWeights;
        currentScore = decreasedScore;
        improved = true;
        if (verbose) console.log(`Iter ${iter}: ${feature} ↓ → ${metric}=${currentScore}`);
      }
    }

    if (!improved) {
      if (verbose) console.log(`Converged after ${iter + 1} iterations`);
      break;
    }
  }

  return {
    weights: currentWeights,
    metrics: evaluateWeights(dataset, currentWeights, threshold),
    iterations: maxIterations,
  };
}

/**
 * Normalize weights to sum to 1.0
 *
 * @param {Object} weights - Weight configuration
 * @returns {Object} Normalized weights
 */
export function normalizeWeights(weights) {
  const total = Object.values(weights).reduce((a, b) => a + b, 0);
  if (total === 0) return weights;

  const normalized = {};
  for (const [key, value] of Object.entries(weights)) {
    normalized[key] = Math.round((value / total) * 1000) / 1000;
  }

  return normalized;
}

/**
 * Compare two weight configurations
 *
 * @param {Array} dataset - Validation dataset
 * @param {Object} weights1 - First weight configuration
 * @param {Object} weights2 - Second weight configuration
 * @param {string} label1 - Label for first config
 * @param {string} label2 - Label for second config
 * @returns {Object} Comparison results
 */
export function compareWeights(dataset, weights1, weights2, label1 = 'Config 1', label2 = 'Config 2') {
  const { threshold: t1, metrics: metrics1 } = findOptimalThreshold(dataset, weights1);
  const { threshold: t2, metrics: metrics2 } = findOptimalThreshold(dataset, weights2);

  return {
    [label1]: {
      threshold: t1,
      ...metrics1,
    },
    [label2]: {
      threshold: t2,
      ...metrics2,
    },
    winner: metrics1.f1 > metrics2.f1 ? label1 : label2,
    improvement: {
      accuracy: Math.round((metrics2.accuracy - metrics1.accuracy) * 1000) / 1000,
      f1: Math.round((metrics2.f1 - metrics1.f1) * 1000) / 1000,
      auc: Math.round((metrics2.auc - metrics1.auc) * 1000) / 1000,
    },
  };
}

/**
 * Generate a calibration report
 *
 * @param {Array} dataset - Validation dataset
 * @param {Object} weights - Weight configuration
 * @returns {Object} Calibration report
 */
export function generateCalibrationReport(dataset, weights) {
  const { threshold, metrics } = findOptimalThreshold(dataset, weights);
  const cv = crossValidate(dataset, weights, 5, threshold);

  // Analyze score distribution
  const scores = dataset.map(({ scores: s }) => scorePrimerPair(s, weights));
  const successScores = dataset.filter(d => d.success).map(d => scorePrimerPair(d.scores, weights));
  const failureScores = dataset.filter(d => !d.success).map(d => scorePrimerPair(d.scores, weights));

  const avgScore = (arr) => arr.length > 0 ? arr.reduce((a, b) => a + b, 0) / arr.length : 0;
  const stdScore = (arr) => {
    if (arr.length === 0) return 0;
    const avg = avgScore(arr);
    return Math.sqrt(arr.reduce((sum, x) => sum + Math.pow(x - avg, 2), 0) / arr.length);
  };

  return {
    datasetSize: dataset.length,
    successCount: successScores.length,
    failureCount: failureScores.length,
    optimalThreshold: threshold,
    metrics,
    crossValidation: cv,
    scoreDistribution: {
      all: {
        mean: Math.round(avgScore(scores) * 10) / 10,
        std: Math.round(stdScore(scores) * 10) / 10,
        min: Math.round(Math.min(...scores) * 10) / 10,
        max: Math.round(Math.max(...scores) * 10) / 10,
      },
      success: {
        mean: Math.round(avgScore(successScores) * 10) / 10,
        std: Math.round(stdScore(successScores) * 10) / 10,
      },
      failure: {
        mean: Math.round(avgScore(failureScores) * 10) / 10,
        std: Math.round(stdScore(failureScores) * 10) / 10,
      },
    },
    qualityTierBreakdown: analyzeQualityTiers(dataset, weights),
  };
}

/**
 * Analyze success rate by quality tier
 *
 * @param {Array} dataset - Validation dataset
 * @param {Object} weights - Weight configuration
 * @returns {Object} Success rates by tier
 */
function analyzeQualityTiers(dataset, weights) {
  const tiers = {
    excellent: { total: 0, success: 0 },
    good: { total: 0, success: 0 },
    acceptable: { total: 0, success: 0 },
    marginal: { total: 0, success: 0 },
    poor: { total: 0, success: 0 },
  };

  for (const { scores, success } of dataset) {
    const score = scorePrimerPair(scores, weights);
    const { tier } = classifyQuality(score);

    tiers[tier].total++;
    if (success) tiers[tier].success++;
  }

  const result = {};
  for (const [tier, data] of Object.entries(tiers)) {
    result[tier] = {
      total: data.total,
      successRate: data.total > 0
        ? Math.round((data.success / data.total) * 1000) / 1000
        : null,
    };
  }

  return result;
}
