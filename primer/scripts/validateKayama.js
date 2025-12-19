#!/usr/bin/env node
/**
 * Kayama RNN Dataset Validation Script
 *
 * Validates the Döring-calibrated primer scoring weights on an independent
 * dataset from Kayama et al. (2021) - "Prediction of PCR amplification from
 * primer and template sequences using recurrent neural network"
 *
 * This serves as a cross-domain validation:
 * - Training: Döring dataset (immunoglobulin PCR, 829 pairs)
 * - Testing: Kayama dataset (16S rRNA PCR, 72 primer sets × 31 templates)
 *
 * Key difference: Kayama provides raw sequences, so we compute all features
 * using our own scoring algorithms (unlike Döring which had pre-computed features).
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

// Import scoring functions
import {
  scoreTm,
  scoreGc,
  scoreHomodimer,
  scoreHairpin,
  scoreHomopolymer,
  scoreLength,
  scoreGcClamp,
  scoreTerminal3DG,
  piecewiseLogistic,
} from '../src/lib/scoring.js';

// Import Tm calculation
import {
  calculateTmQ5,
  calculateGC,
  calculate3primeTerminalDG,
} from '../src/lib/tmQ5.js';

// Import calibration utilities
import {
  DEFAULT_WEIGHTS,
  evaluateWeights,
  crossValidate,
  findOptimalThreshold,
  calculateMetrics,
  calculateAUC,
} from '../src/lib/weightCalibration.js';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

/**
 * Parse simple CSV (no quoted fields in Kayama data)
 */
function parseCSV(csvContent) {
  const lines = csvContent.trim().split('\n');
  const headers = lines[0].split(',');
  const data = [];

  for (let i = 1; i < lines.length; i++) {
    if (!lines[i].trim()) continue;
    const values = lines[i].split(',');
    const row = {};
    headers.forEach((h, idx) => {
      row[h.trim()] = values[idx]?.trim() || '';
    });
    data.push(row);
  }

  return data;
}

/**
 * Estimate hairpin ΔG from sequence
 * Simple heuristic based on GC content and palindromic potential
 */
function estimateHairpinDG(sequence) {
  const seq = sequence.toUpperCase();

  // Check for palindromic regions (simple approach)
  let maxPalindromeLen = 0;
  for (let i = 0; i < seq.length - 3; i++) {
    for (let len = 4; len <= Math.min(10, seq.length - i); len++) {
      const substr = seq.slice(i, i + len);
      const revComp = substr.split('').reverse().map(b => {
        const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
        return comp[b] || b;
      }).join('');

      // Check if reverse complement exists later in sequence
      if (seq.indexOf(revComp, i + len) !== -1) {
        maxPalindromeLen = Math.max(maxPalindromeLen, len);
      }
    }
  }

  // Estimate ΔG based on palindrome length and GC content
  const gc = calculateGC(seq);
  let estimatedDG = 0;

  if (maxPalindromeLen >= 4) {
    // More stable hairpin with longer stem
    estimatedDG = -1.5 * maxPalindromeLen * (0.5 + gc);
  }

  return Math.max(-15, estimatedDG); // Cap at reasonable minimum
}

/**
 * Estimate homodimer ΔG from sequence
 * Based on 3' end complementarity
 */
function estimateHomodimerDG(sequence) {
  const seq = sequence.toUpperCase();
  const last6 = seq.slice(-6);

  // Get reverse complement of last 6 bases
  const revComp = last6.split('').reverse().map(b => {
    const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
    return comp[b] || b;
  }).join('');

  // Count matches with sequence (self-complementarity at 3' end)
  let matches = 0;
  for (let i = 0; i <= seq.length - 6; i++) {
    const window = seq.slice(i, i + 6);
    let windowMatches = 0;
    for (let j = 0; j < 6; j++) {
      if (window[j] === revComp[j]) windowMatches++;
    }
    matches = Math.max(matches, windowMatches);
  }

  // Estimate ΔG based on matches
  const gc = calculateGC(seq);
  const estimatedDG = -1.2 * matches * (0.6 + gc * 0.4);

  return Math.max(-15, estimatedDG);
}

/**
 * Compute all primer scores from sequence
 */
function computePrimerScores(sequence) {
  const seq = sequence.toUpperCase().replace(/[^ATGC]/g, '');

  if (seq.length < 10) {
    return null; // Invalid primer
  }

  const scores = {};

  // Tm scoring
  const tm = calculateTmQ5(seq);
  scores.tm = tm;
  scores.tmScore = scoreTm(tm);

  // GC content
  const gc = calculateGC(seq) * 100;
  scores.gc = gc;
  scores.gcScore = scoreGc(gc);

  // Length
  scores.length = seq.length;
  scores.lengthScore = scoreLength(seq.length);

  // GC clamp
  scores.gcClampScore = scoreGcClamp(seq);

  // Homopolymer
  scores.homopolymerScore = scoreHomopolymer(seq);

  // 3' terminal ΔG
  const terminal3DG = calculate3primeTerminalDG(seq);
  scores.terminal3DG = terminal3DG.dG;
  scores.terminal3DGScore = scoreTerminal3DG(terminal3DG.dG);

  // Hairpin (estimated)
  const hairpinDG = estimateHairpinDG(seq);
  scores.hairpinDG = hairpinDG;
  scores.hairpinScore = scoreHairpin(hairpinDG);

  // Homodimer (estimated)
  const homodimerDG = estimateHomodimerDG(seq);
  scores.homodimerDG = homodimerDG;
  scores.homodimerScore = scoreHomodimer(homodimerDG);

  return scores;
}

/**
 * Convert primer pair scores to calibration format
 */
function convertToCalibrationScores(fwdScores, revScores) {
  if (!fwdScores || !revScores) return null;

  return {
    // Forward primer scores
    tmFwd: fwdScores.tmScore,
    gcFwd: fwdScores.gcScore,
    lengthFwd: fwdScores.lengthScore,
    gcClampFwd: fwdScores.gcClampScore,
    homopolymerFwd: fwdScores.homopolymerScore,
    hairpinFwd: fwdScores.hairpinScore,
    selfDimerFwd: fwdScores.homodimerScore,
    terminalBaseFwd: 0.8, // Default

    // Reverse primer scores
    tmRev: revScores.tmScore,
    gcRev: revScores.gcScore,
    lengthRev: revScores.lengthScore,
    gcClampRev: revScores.gcClampScore,
    homopolymerRev: revScores.homopolymerScore,
    hairpinRev: revScores.hairpinScore,
    selfDimerRev: revScores.homodimerScore,
    terminalBaseRev: 0.8,

    // Pair scores
    tmDiff: scoreTmDiff(fwdScores.tm, revScores.tm),
    lengthDiff: scoreLengthDiff(fwdScores.length, revScores.length),
    terminal3DG: (fwdScores.terminal3DGScore + revScores.terminal3DGScore) / 2,

    // Heterodimer - estimate from sequence complementarity
    heterodimer: 0.8, // Conservative default

    // Off-target - no BLAST data for Kayama, assume good
    offTarget: 1.0,

    // Sanger-specific defaults
    ampliconLength: 0.9,
    ampliconStructure: 0.9,
    distanceToROI: 0.9,
  };
}

/**
 * Score Tm difference
 */
function scoreTmDiff(tm1, tm2) {
  const diff = Math.abs(tm1 - tm2);
  if (diff <= 3) return 1.0;
  if (diff <= 5) return 0.8;
  if (diff <= 8) return 0.5;
  return 0.3;
}

/**
 * Score length difference
 */
function scoreLengthDiff(len1, len2) {
  const diff = Math.abs(len1 - len2);
  if (diff <= 2) return 1.0;
  if (diff <= 4) return 0.9;
  if (diff <= 6) return 0.7;
  return 0.5;
}

/**
 * Load Kayama dataset
 */
function loadKayamaDataset() {
  // Load primers
  const primersPath = path.join(__dirname, '../paper/kayama_primers.csv');
  const primersContent = fs.readFileSync(primersPath, 'utf-8');
  const primersData = parseCSV(primersContent);

  console.log(`Loaded ${primersData.length} primer entries`);

  // Group primers by set
  const primerSets = {};
  for (const row of primersData) {
    const setId = row.primer_set;
    if (!primerSets[setId]) {
      primerSets[setId] = { forward: null, reverse: null };
    }

    // Determine if forward or reverse based on name
    const name = row.primer_name.toLowerCase();
    if (name.endsWith('f') || name.includes('_f')) {
      primerSets[setId].forward = row.sequence;
    } else if (name.endsWith('r') || name.includes('_r')) {
      primerSets[setId].reverse = row.sequence;
    }
  }

  // Load PCR results
  const resultsPath = path.join(__dirname, '../paper/kayama_pcr_results.csv');
  const resultsContent = fs.readFileSync(resultsPath, 'utf-8');
  const resultsData = parseCSV(resultsContent);

  console.log(`Loaded ${resultsData.length} PCR result rows`);

  // Create dataset entries
  const dataset = [];
  let skipped = 0;

  for (const row of resultsData) {
    const setId = row.primer_set;
    const primers = primerSets[setId];

    if (!primers || !primers.forward || !primers.reverse) {
      skipped++;
      continue;
    }

    // Count successes across all templates (t1-t31)
    let successCount = 0;
    let totalTemplates = 0;

    for (let t = 1; t <= 31; t++) {
      const result = row[`t${t}`];
      if (result !== undefined && result !== '') {
        totalTemplates++;
        if (result === '1') successCount++;
      }
    }

    // Compute primer scores
    const fwdScores = computePrimerScores(primers.forward);
    const revScores = computePrimerScores(primers.reverse);

    if (!fwdScores || !revScores) {
      skipped++;
      continue;
    }

    // Convert to calibration format
    const scores = convertToCalibrationScores(fwdScores, revScores);

    // Determine success - use majority vote across templates
    // A primer pair is "successful" if it works on >50% of templates
    const successRate = totalTemplates > 0 ? successCount / totalTemplates : 0;
    const success = successRate >= 0.5;

    dataset.push({
      primerSet: setId,
      fwdSequence: primers.forward,
      revSequence: primers.reverse,
      scores,
      success,
      successRate,
      successCount,
      totalTemplates,
      metadata: {
        tmFwd: fwdScores.tm,
        tmRev: revScores.tm,
        gcFwd: fwdScores.gc,
        gcRev: revScores.gc,
      },
    });
  }

  console.log(`Created ${dataset.length} dataset entries (${skipped} skipped)`);

  return dataset;
}

/**
 * Analyze feature distributions
 */
function analyzeFeatureDistributions(dataset) {
  const successData = dataset.filter(d => d.success);
  const failureData = dataset.filter(d => !d.success);

  const avg = arr => arr.length > 0 ? arr.reduce((a, b) => a + b, 0) / arr.length : 0;

  console.log('\n--- Feature Analysis (Kayama Dataset) ---');
  console.log(`Success: ${successData.length}, Failure: ${failureData.length}`);

  const features = ['tmFwd', 'gcFwd', 'terminal3DG', 'hairpinFwd', 'selfDimerFwd'];

  for (const feature of features) {
    const successVals = successData.map(d => d.scores[feature] || 0);
    const failureVals = failureData.map(d => d.scores[feature] || 0);
    const diff = avg(successVals) - avg(failureVals);

    console.log(`  ${feature}: success=${avg(successVals).toFixed(3)}, failure=${avg(failureVals).toFixed(3)}, diff=${diff.toFixed(3)}`);
  }
}

/**
 * Load expanded dataset (individual primer-template pairs)
 * Each primer set × template combination is a separate entry
 */
function loadKayamaDatasetExpanded() {
  // Load primers
  const primersPath = path.join(__dirname, '../paper/kayama_primers.csv');
  const primersContent = fs.readFileSync(primersPath, 'utf-8');
  const primersData = parseCSV(primersContent);

  // Group primers by set
  const primerSets = {};
  for (const row of primersData) {
    const setId = row.primer_set;
    if (!primerSets[setId]) {
      primerSets[setId] = { forward: null, reverse: null };
    }

    const name = row.primer_name.toLowerCase();
    if (name.endsWith('f') || name.includes('_f')) {
      primerSets[setId].forward = row.sequence;
    } else if (name.endsWith('r') || name.includes('_r')) {
      primerSets[setId].reverse = row.sequence;
    }
  }

  // Load PCR results
  const resultsPath = path.join(__dirname, '../paper/kayama_pcr_results.csv');
  const resultsContent = fs.readFileSync(resultsPath, 'utf-8');
  const resultsData = parseCSV(resultsContent);

  // Create expanded dataset - one entry per primer-template pair
  const dataset = [];
  let skipped = 0;

  for (const row of resultsData) {
    const setId = row.primer_set;
    const primers = primerSets[setId];

    if (!primers || !primers.forward || !primers.reverse) {
      skipped++;
      continue;
    }

    // Compute primer scores once per set
    const fwdScores = computePrimerScores(primers.forward);
    const revScores = computePrimerScores(primers.reverse);

    if (!fwdScores || !revScores) {
      skipped++;
      continue;
    }

    const scores = convertToCalibrationScores(fwdScores, revScores);

    // Create entry for each template
    for (let t = 1; t <= 31; t++) {
      const result = row[`t${t}`];
      if (result === undefined || result === '') continue;

      const success = result === '1';

      dataset.push({
        primerSet: setId,
        template: t,
        scores,
        success,
        metadata: {
          tmFwd: fwdScores.tm,
          tmRev: revScores.tm,
        },
      });
    }
  }

  return dataset;
}

/**
 * Run validation
 */
function runValidation() {
  console.log('='.repeat(60));
  console.log('KAYAMA DATASET VALIDATION - Cross-Domain Test');
  console.log('='.repeat(60));
  console.log();
  console.log('Testing Döring-calibrated weights on independent Kayama dataset:');
  console.log('- Döring: Immunoglobulin PCR (829 pairs)');
  console.log('- Kayama: 16S rRNA PCR (72 primer sets × 31 templates)');
  console.log();

  // Load aggregated dataset (primer set level)
  console.log('=== AGGREGATED ANALYSIS (per primer set) ===');
  const dataset = loadKayamaDataset();

  if (dataset.length === 0) {
    console.error('ERROR: No valid dataset entries');
    return;
  }

  // Count outcomes
  const successCount = dataset.filter(d => d.success).length;
  const failureCount = dataset.length - successCount;
  console.log(`\nDataset: ${successCount} successful (${(100 * successCount / dataset.length).toFixed(1)}%), ${failureCount} failed`);

  // Analyze feature distributions
  analyzeFeatureDistributions(dataset);

  // Evaluate with DEFAULT_WEIGHTS (Döring-calibrated)
  console.log('\n--- Validation with Döring-Calibrated Weights ---');

  const { threshold, metrics } = findOptimalThreshold(dataset, DEFAULT_WEIGHTS);

  console.log(`Optimal threshold: ${threshold}`);
  console.log(`Accuracy: ${(metrics.accuracy * 100).toFixed(1)}%`);
  console.log(`Precision: ${(metrics.precision * 100).toFixed(1)}%`);
  console.log(`Recall: ${(metrics.recall * 100).toFixed(1)}%`);
  console.log(`F1 Score: ${(metrics.f1 * 100).toFixed(1)}%`);
  console.log(`AUC-ROC: ${metrics.auc.toFixed(3)}`);

  // === EXPANDED ANALYSIS ===
  console.log('\n' + '='.repeat(60));
  console.log('=== EXPANDED ANALYSIS (per primer-template pair) ===');

  const expandedDataset = loadKayamaDatasetExpanded();
  const expSuccessCount = expandedDataset.filter(d => d.success).length;
  const expFailureCount = expandedDataset.length - expSuccessCount;
  console.log(`\nDataset: ${expSuccessCount} amplified (${(100 * expSuccessCount / expandedDataset.length).toFixed(1)}%), ${expFailureCount} not amplified`);

  console.log('\n--- Validation with Döring-Calibrated Weights ---');
  const { threshold: expThreshold, metrics: expMetrics } = findOptimalThreshold(expandedDataset, DEFAULT_WEIGHTS);

  console.log(`Optimal threshold: ${expThreshold}`);
  console.log(`Accuracy: ${(expMetrics.accuracy * 100).toFixed(1)}%`);
  console.log(`Precision: ${(expMetrics.precision * 100).toFixed(1)}%`);
  console.log(`Recall: ${(expMetrics.recall * 100).toFixed(1)}%`);
  console.log(`F1 Score: ${(expMetrics.f1 * 100).toFixed(1)}%`);
  console.log(`AUC-ROC: ${expMetrics.auc.toFixed(3)}`);

  // Cross-validation on expanded dataset
  console.log('\n--- 5-Fold Cross-Validation ---');
  const cv = crossValidate(expandedDataset, DEFAULT_WEIGHTS, 5, expThreshold);
  console.log(`Mean F1: ${(cv.mean.f1 * 100).toFixed(1)}% ± ${(cv.std.f1 * 100).toFixed(1)}%`);
  console.log(`Mean AUC: ${cv.mean.auc.toFixed(3)} ± ${cv.std.auc.toFixed(3)}`);

  // Compare with baseline (random guess)
  console.log('\n--- Baseline Comparison ---');
  const baselineAccuracy = Math.max(expSuccessCount, expFailureCount) / expandedDataset.length;
  console.log(`Majority class baseline: ${(baselineAccuracy * 100).toFixed(1)}%`);
  console.log(`Our model improvement: +${((expMetrics.accuracy - baselineAccuracy) * 100).toFixed(1)}%`);

  // Summary
  console.log('\n' + '='.repeat(60));
  console.log('SUMMARY');
  console.log('='.repeat(60));
  console.log(`Training dataset: Döring (829 pairs, F1=81.9%)`);
  console.log(`Test dataset: Kayama (${expandedDataset.length} primer-template pairs)`);
  console.log(`Test F1: ${(expMetrics.f1 * 100).toFixed(1)}%`);
  console.log(`Test AUC: ${expMetrics.auc.toFixed(3)}`);
  console.log();

  // Interpretation
  if (expMetrics.f1 >= 0.70) {
    console.log('✅ GOOD: Weights generalize well to independent dataset');
  } else if (expMetrics.f1 >= 0.55) {
    console.log('⚠️  MODERATE: Some generalization, may benefit from domain-specific tuning');
  } else {
    console.log('❌ POOR: Weights may be overfit to Döring dataset');
  }

  // Key insight
  console.log('\n--- Key Insights ---');
  console.log('Note: Kayama dataset uses different PCR conditions (16S rRNA vs immunoglobulin)');
  console.log('The primer scoring captures sequence-level features, but PCR success');
  console.log('also depends on template-specific factors not in our model.');

  return {
    dataset,
    expandedDataset,
    metrics,
    expMetrics,
    threshold,
    cv,
  };
}

// Run if executed directly
runValidation();
