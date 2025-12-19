/**
 * Weight Calibration Runner
 *
 * Runs weight calibration using the empirically-derived validation dataset
 * and outputs optimized weights for primer scoring.
 */

import {
  VALIDATION_DATASET,
  toCalibrationFormat,
  getDatasetStats,
  splitDataset,
} from './validationDataset.js';

import {
  DEFAULT_WEIGHTS,
  crossValidate,
  gridSearch,
  coordinateDescent,
  compareWeights,
  generateCalibrationReport,
  calculateMetrics,
  calculateAUC,
} from './weightCalibration.js';

/**
 * Run full calibration pipeline
 */
export function runCalibration(options = {}) {
  const {
    verbose = true,
    useGridSearch = true,
    useCoordinateDescent = true,
    validationFolds = 5,
  } = options;

  const results = {
    datasetStats: null,
    baselineMetrics: null,
    gridSearchResults: null,
    coordinateDescentResults: null,
    optimizedWeights: null,
    improvement: null,
  };

  // Step 1: Analyze dataset
  if (verbose) console.log('=== PRIMER SCORING WEIGHT CALIBRATION ===\n');

  const stats = getDatasetStats();
  results.datasetStats = stats;

  if (verbose) {
    console.log('Dataset Statistics:');
    console.log(`  Total entries: ${stats.total}`);
    console.log(`  Success rate: ${(stats.successRate * 100).toFixed(1)}%`);
    console.log(`  By category:`);
    for (const [cat, data] of Object.entries(stats.byCategory)) {
      const rate = ((data.success / data.total) * 100).toFixed(1);
      console.log(`    ${cat}: ${data.total} entries, ${rate}% success`);
    }
    console.log('');
  }

  // Step 2: Convert to calibration format
  const calibrationData = toCalibrationFormat(VALIDATION_DATASET);

  // Step 3: Baseline evaluation with default weights
  if (verbose) console.log('Evaluating baseline (DEFAULT_WEIGHTS)...');

  const baselineCVResult = crossValidate(calibrationData, DEFAULT_WEIGHTS, validationFolds);
  const baselineCV = baselineCVResult.mean;  // Extract mean metrics
  results.baselineMetrics = baselineCVResult;

  if (verbose) {
    console.log('Baseline Metrics:');
    console.log(`  Accuracy: ${(baselineCV.accuracy * 100).toFixed(1)}%`);
    console.log(`  Precision: ${(baselineCV.precision * 100).toFixed(1)}%`);
    console.log(`  Recall: ${(baselineCV.recall * 100).toFixed(1)}%`);
    console.log(`  F1 Score: ${(baselineCV.f1 * 100).toFixed(1)}%`);
    console.log(`  AUC-ROC: ${baselineCV.auc.toFixed(3)}`);
    console.log('');
  }

  let bestWeights = { ...DEFAULT_WEIGHTS };
  let bestMetrics = baselineCV;

  // Step 4: Grid search optimization
  if (useGridSearch) {
    if (verbose) console.log('Running grid search optimization...');

    // Focus on key features identified in literature
    const keyFeatures = ['offTarget', 'terminal3DG', 'tmFwd', 'gcFwd', 'hairpinFwd', 'heterodimer'];

    try {
      const gridResults = gridSearch(calibrationData, keyFeatures, DEFAULT_WEIGHTS, {
        threshold: 60,
        metric: 'f1',
      });

      results.gridSearchResults = gridResults;

      const gridImprovement = gridResults.metrics.f1 - baselineCV.f1;
      if (gridImprovement > 0) {
        bestWeights = gridResults.weights;
        bestMetrics = gridResults.metrics;

        if (verbose) {
          console.log('Grid Search Results:');
          console.log(`  Best F1: ${(gridResults.metrics.f1 * 100).toFixed(1)}%`);
          console.log(`  Improvement: +${(gridImprovement * 100).toFixed(1)}%`);
          console.log(`  Combinations: ${gridResults.combinationsTested}`);
          console.log('');
        }
      } else if (verbose) {
        console.log('Grid search did not improve baseline.\n');
      }
    } catch (err) {
      if (verbose) console.log(`Grid search error: ${err.message}\n`);
    }
  }

  // Step 5: Coordinate descent refinement
  if (useCoordinateDescent) {
    if (verbose) console.log('Running coordinate descent refinement...');

    try {
      const cdResults = coordinateDescent(calibrationData, bestWeights, {
        maxIterations: 50,
        minImprovement: 0.001,
        stepSize: 0.02,
        threshold: 60,
        metric: 'f1',
      });

      results.coordinateDescentResults = cdResults;

      const cdImprovement = cdResults.metrics.f1 - bestMetrics.f1;
      if (cdImprovement > 0) {
        bestWeights = cdResults.weights;
        bestMetrics = cdResults.metrics;

        if (verbose) {
          console.log('Coordinate Descent Results:');
          console.log(`  Best F1: ${(cdResults.metrics.f1 * 100).toFixed(1)}%`);
          console.log(`  Improvement: +${(cdImprovement * 100).toFixed(1)}%`);
          console.log(`  Iterations: ${cdResults.iterations}`);
          console.log('');
        }
      } else if (verbose) {
        console.log('Coordinate descent did not improve further.\n');
      }
    } catch (err) {
      if (verbose) console.log(`Coordinate descent error: ${err.message}\n`);
    }
  }

  // Step 6: Final comparison
  results.optimizedWeights = bestWeights;
  results.improvement = bestMetrics.f1 - baselineCV.f1;

  if (verbose) {
    console.log('=== CALIBRATION RESULTS ===\n');
    console.log('Optimized Weights:');
    const sortedWeights = Object.entries(bestWeights)
      .sort((a, b) => b[1] - a[1]);
    for (const [feature, weight] of sortedWeights) {
      const defaultWeight = DEFAULT_WEIGHTS[feature] || 0;
      const change = weight - defaultWeight;
      const changeStr = change === 0 ? '' : ` (${change > 0 ? '+' : ''}${change.toFixed(3)})`;
      console.log(`  ${feature}: ${weight.toFixed(3)}${changeStr}`);
    }
    console.log('');

    console.log('Performance Comparison:');
    console.log('                   Baseline    Optimized    Change');
    console.log(`  Accuracy:        ${(baselineCV.accuracy * 100).toFixed(1)}%        ${(bestMetrics.accuracy * 100).toFixed(1)}%        ${((bestMetrics.accuracy - baselineCV.accuracy) * 100).toFixed(1)}%`);
    console.log(`  Precision:       ${(baselineCV.precision * 100).toFixed(1)}%        ${(bestMetrics.precision * 100).toFixed(1)}%        ${((bestMetrics.precision - baselineCV.precision) * 100).toFixed(1)}%`);
    console.log(`  Recall:          ${(baselineCV.recall * 100).toFixed(1)}%        ${(bestMetrics.recall * 100).toFixed(1)}%        ${((bestMetrics.recall - baselineCV.recall) * 100).toFixed(1)}%`);
    console.log(`  F1 Score:        ${(baselineCV.f1 * 100).toFixed(1)}%        ${(bestMetrics.f1 * 100).toFixed(1)}%        ${((bestMetrics.f1 - baselineCV.f1) * 100).toFixed(1)}%`);
    console.log(`  AUC-ROC:         ${baselineCV.auc.toFixed(3)}        ${bestMetrics.auc.toFixed(3)}        ${(bestMetrics.auc - baselineCV.auc).toFixed(3)}`);
    console.log('');
  }

  return results;
}

/**
 * Generate weights export code
 */
export function generateWeightsExport(weights) {
  const sortedWeights = Object.entries(weights)
    .sort((a, b) => b[1] - a[1]);

  let code = '/**\n';
  code += ' * Calibrated weights for primer scoring\n';
  code += ' * Generated from empirical validation dataset\n';
  code += ` * Date: ${new Date().toISOString().split('T')[0]}\n`;
  code += ' */\n';
  code += 'export const CALIBRATED_WEIGHTS = {\n';

  for (const [feature, weight] of sortedWeights) {
    code += `  ${feature}: ${weight.toFixed(4)},\n`;
  }

  code += '};\n';

  return code;
}

/**
 * Quick validation run (for testing)
 */
export function quickValidation() {
  const calibrationData = toCalibrationFormat(VALIDATION_DATASET);
  const metrics = crossValidate(calibrationData, DEFAULT_WEIGHTS, 5);

  return {
    datasetSize: VALIDATION_DATASET.length,
    successRate: getDatasetStats().successRate,
    accuracy: metrics.accuracy,
    precision: metrics.precision,
    recall: metrics.recall,
    f1: metrics.f1,
    auc: metrics.auc,
  };
}

// Run if executed directly
if (typeof process !== 'undefined' && process.argv && process.argv[1] && process.argv[1].includes('runCalibration')) {
  console.log('\nStarting calibration...\n');
  const results = runCalibration({ verbose: true });

  console.log('\n=== EXPORTABLE WEIGHTS ===\n');
  console.log(generateWeightsExport(results.optimizedWeights));
}

export default {
  runCalibration,
  generateWeightsExport,
  quickValidation,
};
