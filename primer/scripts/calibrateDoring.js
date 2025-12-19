#!/usr/bin/env node
/**
 * Calibration Script for Döring/openPrimeR Dataset
 *
 * Uses the 829-datapoint feature_matrix.csv to calibrate primer scoring weights.
 * This replaces the synthetic 125-entry dataset with real experimental data.
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

// Import calibration functions
import {
  DEFAULT_WEIGHTS,
  evaluateWeights,
  crossValidate,
  gridSearch,
  coordinateDescent,
  findOptimalThreshold,
  calculateMetrics,
  calculateAUC,
  compareWeights,
} from '../src/lib/weightCalibration.js';

// Import scoring functions
import {
  scoreTm,
  scoreGc,
  scoreHomodimer,
  scoreHairpin,
  scoreHomopolymer as scoreHomopolymerSeq,
  scoreLength,
  scoreGcClamp as scoreGcClampSeq,
  scoreGQuadruplex,
  piecewiseLogistic,
  analyzePrimerBatch,
} from '../src/lib/scoring.js';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

/**
 * Parse CSV data into structured objects
 * Handles quoted fields containing commas properly
 */
function parseCSV(csvContent) {
  const lines = csvContent.trim().split('\n');

  // Parse CSV line respecting quoted fields
  function parseLine(line) {
    const values = [];
    let currentValue = '';
    let inQuotes = false;

    for (let i = 0; i < line.length; i++) {
      const char = line[i];

      if (char === '"') {
        inQuotes = !inQuotes;
      } else if (char === ',' && !inQuotes) {
        values.push(currentValue.replace(/"/g, '').trim());
        currentValue = '';
      } else {
        currentValue += char;
      }
    }
    // Push the last value
    values.push(currentValue.replace(/"/g, '').trim());

    return values;
  }

  const headers = parseLine(lines[0]);

  const data = [];
  for (let i = 1; i < lines.length; i++) {
    if (!lines[i].trim()) continue; // Skip empty lines

    const values = parseLine(lines[i]);
    const row = {};
    headers.forEach((h, idx) => {
      row[h] = values[idx] || '';
    });
    data.push(row);
  }

  return data;
}

/**
 * Convert Döring CSV features to scoring system format
 *
 * Maps pre-computed features from openPrimeR to our 0-1 scoring scale
 */
function convertToScores(row) {
  const scores = {};

  // Tm scoring (optimal 55-60°C for Sanger)
  const tm = parseFloat(row.melting_temp) || 0;
  scores.tmFwd = scoreTm(tm);
  scores.tmRev = scoreTm(tm); // Use same value for Fwd (no Rev data in this dataset)

  // GC content scoring (optimal 40-60%)
  const gcRatio = parseFloat(row.gc_ratio_fw) || 0;
  scores.gcFwd = scoreGc(gcRatio * 100); // Convert ratio to percentage
  scores.gcRev = scoreGc(gcRatio * 100);

  // GC clamp scoring (optimal: 1-2 G/C in last 5 bases)
  const gcClamp = parseInt(row.gc_clamp_fw) || 0;
  scores.gcClampFwd = scoreGCClampFromCount(gcClamp);
  scores.gcClampRev = scoreGCClampFromCount(gcClamp);

  // Self-dimer scoring (ΔG more negative = worse)
  const selfDimerDG = parseFloat(row.Self_Dimer_DeltaG) || 0;
  scores.selfDimerFwd = scoreHomodimer(selfDimerDG);
  scores.selfDimerRev = scoreHomodimer(selfDimerDG);

  // Hairpin/structure scoring
  const structureDG = parseFloat(row.Structure_deltaG_fw) || 0;
  scores.hairpinFwd = scoreHairpin(structureDG);
  scores.hairpinRev = scoreHairpin(structureDG);

  // Homopolymer runs
  const runs = parseInt(row.no_runs_fw) || 1;
  scores.homopolymerFwd = scoreHomopolymerFromRuns(runs);
  scores.homopolymerRev = scoreHomopolymerFromRuns(runs);

  // Primer length scoring
  const length = parseInt(row.primer_length_fw) || 20;
  scores.lengthFwd = scoreLength(length);
  scores.lengthRev = scoreLength(length);

  // Mismatch/off-target scoring based on position of 3' terminus mismatch
  const numMismatches = parseInt(row.Number_of_mismatches) || 0;
  const pos3terminus = parseInt(row.Position_3terminus) || 22; // Position from 3' end
  scores.offTarget = scoreOffTargetFromMismatches(numMismatches, pos3terminus);

  // Annealing ΔG as a proxy for binding strength
  const annealingDG = parseFloat(row.annealing_DeltaG) || -15;
  scores.terminal3DG = scoreAnnealingDG(annealingDG);

  // Heterodimer - estimate conservatively (no cross-primer data available)
  // NOTE: Do NOT use primer_efficiency here - it's the target variable!
  scores.heterodimer = 0.8; // Assume generally acceptable without specific evidence

  // G-Quadruplex scoring from primer sequence
  const primerSeq = row.Primer_Sequence || '';
  if (primerSeq) {
    scores.gQuadruplexFwd = scoreGQuadruplex(primerSeq);
    scores.gQuadruplexRev = scoreGQuadruplex(primerSeq); // Use same (no rev sequence in dataset)
  } else {
    scores.gQuadruplexFwd = 1.0;  // No sequence = assume no G4
    scores.gQuadruplexRev = 1.0;
  }

  // Set defaults for features not in the dataset
  scores.tmDiff = 0.9;  // Assume reasonable Tm matching
  scores.ampliconLength = 0.9;
  scores.ampliconStructure = 0.9;
  scores.terminalBaseFwd = 0.8;
  scores.terminalBaseRev = 0.8;
  scores.lengthDiff = 0.9;
  scores.distanceToROI = 0.9;

  return scores;
}

/**
 * Score GC clamp (number of G/C in last 5 bases of 3' end)
 * Optimal: 1-2 G/C, Acceptable: 0-3, Poor: 4-5
 */
function scoreGCClampFromCount(gcClamp) {
  if (gcClamp >= 1 && gcClamp <= 2) return 1.0;
  if (gcClamp === 0) return 0.7;
  if (gcClamp === 3) return 0.7;
  if (gcClamp === 4) return 0.4;
  return 0.2; // 5 G/C = strong 3' secondary structure risk
}

/**
 * Score off-target potential based on mismatch count and 3' position
 *
 * More mismatches = less likely to bind (good for specificity)
 * But mismatches at 3' end reduce extension efficiency
 */
function scoreOffTargetFromMismatches(numMismatches, pos3terminus) {
  // Perfect match = good specificity depends on context
  // For immunoglobulin primers, some degeneracy is expected

  // Position 3terminus indicates distance of nearest mismatch from 3' end
  // Small values (1-3) = mismatch near 3' end = bad for extension

  if (numMismatches === 0) {
    return 1.0; // Perfect match
  }

  // Mismatch position matters most
  if (pos3terminus <= 3) {
    // Mismatch in critical 3' region
    return Math.max(0.1, 0.5 - numMismatches * 0.05);
  } else if (pos3terminus <= 6) {
    // Mismatch in semi-critical region
    return Math.max(0.2, 0.7 - numMismatches * 0.03);
  } else {
    // Mismatch far from 3' end - less critical
    return Math.max(0.4, 0.9 - numMismatches * 0.02);
  }
}

/**
 * Score annealing ΔG (binding free energy)
 * More negative = stronger binding = generally better
 * But too negative can mean non-specific binding
 */
function scoreAnnealingDG(dg) {
  // Optimal: -12 to -18 kcal/mol
  // Too weak: > -8 kcal/mol
  // Too strong: < -22 kcal/mol

  return piecewiseLogistic(dg, {
    optimalLow: -18,
    optimalHigh: -12,
    acceptableLow: -22,
    acceptableHigh: -8,
    steepness: 0.3,
  });
}

/**
 * Score homopolymer runs
 */
function scoreHomopolymerFromRuns(runs) {
  if (runs <= 2) return 1.0;
  if (runs === 3) return 0.8;
  if (runs === 4) return 0.5;
  return 0.3;
}

/**
 * Load and process the Döring dataset
 */
function loadDoringDataset() {
  const csvPath = path.join(__dirname, '../paper/feature_matrix.csv');
  const csvContent = fs.readFileSync(csvPath, 'utf-8');
  const rawData = parseCSV(csvContent);

  console.log(`Loaded ${rawData.length} rows from feature_matrix.csv`);

  // Debug: Check data quality on first few rows
  console.log('\n--- Data Quality Check (first 3 rows) ---');
  for (let i = 0; i < 3 && i < rawData.length; i++) {
    const row = rawData[i];
    console.log(`Row ${i + 1}:`);
    console.log(`  melting_temp: "${row.melting_temp}" → parsed: ${parseFloat(row.melting_temp)}`);
    console.log(`  gc_ratio_fw: "${row.gc_ratio_fw}" → parsed: ${parseFloat(row.gc_ratio_fw)}`);
    console.log(`  annealing_DeltaG: "${row.annealing_DeltaG}" → parsed: ${parseFloat(row.annealing_DeltaG)}`);
    console.log(`  Experimental_Coverage: "${row.Experimental_Coverage}"`);
    const scores = convertToScores(row);
    console.log(`  Computed scores: tmFwd=${scores.tmFwd.toFixed(3)}, gcFwd=${scores.gcFwd.toFixed(3)}, offTarget=${scores.offTarget.toFixed(3)}`);
  }

  // Check a failure case too
  const failureRow = rawData.find(r => r.Experimental_Coverage === 'Unamplified');
  if (failureRow) {
    console.log(`\nFirst failure row:`);
    console.log(`  melting_temp: "${failureRow.melting_temp}" → parsed: ${parseFloat(failureRow.melting_temp)}`);
    console.log(`  gc_ratio_fw: "${failureRow.gc_ratio_fw}" → parsed: ${parseFloat(failureRow.gc_ratio_fw)}`);
    const scores = convertToScores(failureRow);
    console.log(`  Computed scores: tmFwd=${scores.tmFwd.toFixed(3)}, gcFwd=${scores.gcFwd.toFixed(3)}, offTarget=${scores.offTarget.toFixed(3)}`);
  }

  // Convert to calibration format
  const dataset = rawData.map(row => {
    const scores = convertToScores(row);
    const success = row.Experimental_Coverage === 'Amplified';

    return {
      scores,
      success,
      actual: success,
      // Keep original data for reference
      metadata: {
        primer: row.Primer,
        template: row.Template,
        efficiency: parseFloat(row.primer_efficiency) || 0,
      }
    };
  });

  // Count successes and failures
  const successes = dataset.filter(d => d.success).length;
  const failures = dataset.length - successes;

  console.log(`Dataset: ${successes} amplified (${(100*successes/dataset.length).toFixed(1)}%), ${failures} unamplified`);

  return dataset;
}

/**
 * Analyze score distributions using unified analysis function
 *
 * Uses analyzePrimerBatch from scoring.js for consistent analysis
 * across standalone, unified designer, mutagenesis, and calibration.
 */
function analyzeScoreDistributions(dataset, weights) {
  // Use unified batch analysis function
  const stats = analyzePrimerBatch(dataset, weights);

  console.log('--- Score Distribution Analysis ---');
  console.log(`Success (Amplified): mean=${stats.scoreDistribution.success.mean}, n=${stats.success}`);
  console.log(`Failure (Unamplified): mean=${stats.scoreDistribution.failure.mean}, n=${stats.failure}`);

  // Find most discriminative features
  console.log('\n--- Feature Discrimination (Success vs Failure mean) ---');
  for (const feat of stats.topDiscriminativeFeatures.slice(0, 8)) {
    console.log(`  ${feat.feature}: success=${feat.successMean.toFixed(3)}, failure=${feat.failureMean.toFixed(3)}, diff=${feat.diff.toFixed(3)}`);
  }

  // G-Quadruplex specific analysis (from unified stats)
  console.log('\n--- G-Quadruplex Analysis ---');
  const g4 = stats.gQuadruplex;
  console.log(`  Primers with canonical G4 motif: ${g4.g4Motif.count} (${g4.g4Motif.successCount} amplified, ${g4.g4Motif.successRate !== null ? (g4.g4Motif.successRate * 100).toFixed(1) : 0}% success)`);
  console.log(`  Primers with GGGG runs: ${g4.ggggRun.count} (${g4.ggggRun.successCount} amplified, ${g4.ggggRun.successRate !== null ? (g4.ggggRun.successRate * 100).toFixed(1) : 0}% success)`);
  console.log(`  Primers with multiple GGG: ${g4.multiGgg.count} (${g4.multiGgg.successCount} amplified, ${g4.multiGgg.successRate !== null ? (g4.multiGgg.successRate * 100).toFixed(1) : 0}% success)`);
  console.log(`  Primers with no G4 risk: ${g4.noRisk.count} (${g4.noRisk.successCount} amplified, ${g4.noRisk.successRate !== null ? (g4.noRisk.successRate * 100).toFixed(1) : 0}% success)`);
}

/**
 * Main calibration routine
 */
async function runCalibration() {
  console.log('='.repeat(60));
  console.log('PRIMER SCORING CALIBRATION - Döring/openPrimeR Dataset');
  console.log('='.repeat(60));
  console.log();

  // Load dataset
  const dataset = loadDoringDataset();
  console.log();

  // Analyze score distributions first
  analyzeScoreDistributions(dataset, DEFAULT_WEIGHTS);
  console.log();

  // 1. Evaluate current DEFAULT_WEIGHTS
  console.log('--- Current DEFAULT_WEIGHTS Performance ---');
  const { threshold: currentThreshold, metrics: currentMetrics } = findOptimalThreshold(dataset, DEFAULT_WEIGHTS);
  console.log(`Optimal threshold: ${currentThreshold}`);
  console.log(`Accuracy: ${(currentMetrics.accuracy * 100).toFixed(1)}%`);
  console.log(`Precision: ${(currentMetrics.precision * 100).toFixed(1)}%`);
  console.log(`Recall: ${(currentMetrics.recall * 100).toFixed(1)}%`);
  console.log(`F1: ${(currentMetrics.f1 * 100).toFixed(1)}%`);
  console.log(`AUC: ${currentMetrics.auc.toFixed(3)}`);
  console.log();

  // 2. Cross-validation with current weights
  console.log('--- 5-Fold Cross-Validation (Current Weights) ---');
  const cvCurrent = crossValidate(dataset, DEFAULT_WEIGHTS, 5, currentThreshold);
  console.log(`Mean F1: ${(cvCurrent.mean.f1 * 100).toFixed(1)}% ± ${(cvCurrent.std.f1 * 100).toFixed(1)}%`);
  console.log(`Mean AUC: ${cvCurrent.mean.auc.toFixed(3)} ± ${cvCurrent.std.auc.toFixed(3)}`);
  console.log();

  // 3. Grid search on key features
  console.log('--- Grid Search Optimization ---');
  const keyFeatures = ['offTarget', 'terminal3DG', 'gcFwd', 'selfDimerFwd', 'heterodimer', 'gQuadruplexFwd', 'gQuadruplexRev'];
  console.log(`Optimizing: ${keyFeatures.join(', ')}`);

  const gridResult = gridSearch(dataset, keyFeatures, DEFAULT_WEIGHTS, {
    threshold: currentThreshold,
    metric: 'f1',
    verbose: false,
  });

  console.log(`Combinations tested: ${gridResult.combinationsTested}`);
  console.log(`Best F1: ${(gridResult.metrics.f1 * 100).toFixed(1)}%`);
  console.log();

  // 4. Coordinate descent refinement
  console.log('--- Coordinate Descent Refinement ---');
  const cdResult = coordinateDescent(dataset, gridResult.weights, {
    maxIterations: 20,
    threshold: currentThreshold,
    metric: 'f1',
    stepSize: 0.01,
    verbose: false,
  });

  console.log(`Final F1: ${(cdResult.metrics.f1 * 100).toFixed(1)}%`);
  console.log(`Final AUC: ${cdResult.metrics.auc.toFixed(3)}`);
  console.log();

  // 5. Cross-validate optimized weights
  console.log('--- 5-Fold Cross-Validation (Optimized Weights) ---');
  const cvOptimized = crossValidate(dataset, cdResult.weights, 5, currentThreshold);
  console.log(`Mean F1: ${(cvOptimized.mean.f1 * 100).toFixed(1)}% ± ${(cvOptimized.std.f1 * 100).toFixed(1)}%`);
  console.log(`Mean AUC: ${cvOptimized.mean.auc.toFixed(3)} ± ${cvOptimized.std.auc.toFixed(3)}`);
  console.log();

  // 6. Compare weights
  console.log('--- Weight Comparison ---');
  const comparison = compareWeights(dataset, DEFAULT_WEIGHTS, cdResult.weights, 'Current', 'Optimized');
  console.log(`Winner: ${comparison.winner}`);
  console.log(`F1 improvement: ${(comparison.improvement.f1 * 100).toFixed(1)}%`);
  console.log(`AUC improvement: ${comparison.improvement.auc.toFixed(3)}`);
  console.log();

  // 7. Print optimized weights
  console.log('--- Optimized Weights ---');
  const changedWeights = {};
  for (const [key, value] of Object.entries(cdResult.weights)) {
    if (Math.abs(value - (DEFAULT_WEIGHTS[key] || 0)) > 0.005) {
      changedWeights[key] = {
        old: DEFAULT_WEIGHTS[key] || 0,
        new: value,
      };
    }
  }

  console.log('Changed weights:');
  for (const [key, { old, new: newVal }] of Object.entries(changedWeights)) {
    console.log(`  ${key}: ${old.toFixed(2)} → ${newVal.toFixed(2)}`);
  }
  console.log();

  // 8. Summary
  console.log('='.repeat(60));
  console.log('SUMMARY');
  console.log('='.repeat(60));
  console.log(`Dataset size: ${dataset.length} primer-template pairs`);
  console.log(`Current weights F1: ${(currentMetrics.f1 * 100).toFixed(1)}%`);
  console.log(`Optimized weights F1: ${(cdResult.metrics.f1 * 100).toFixed(1)}%`);
  console.log(`Improvement: ${(comparison.improvement.f1 * 100).toFixed(1)}%`);
  console.log();

  return {
    dataset,
    currentMetrics,
    optimizedMetrics: cdResult.metrics,
    optimizedWeights: cdResult.weights,
    changedWeights,
  };
}

// Run if executed directly
runCalibration().catch(console.error);
