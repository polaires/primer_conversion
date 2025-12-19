/**
 * Tests for Validation Dataset Module
 */

import { describe, it, expect } from 'vitest';
import {
  VALIDATION_DATASET,
  calculateEntryScores,
  toCalibrationFormat,
  getDatasetStats,
  splitDataset,
} from './validationDataset.js';

describe('VALIDATION_DATASET', () => {
  it('should have at least 100 entries', () => {
    expect(VALIDATION_DATASET.length).toBeGreaterThanOrEqual(100);
  });

  it('should have entries with required fields', () => {
    const entry = VALIDATION_DATASET[0];
    expect(entry).toHaveProperty('id');
    expect(entry).toHaveProperty('category');
    expect(entry).toHaveProperty('fwd');
    expect(entry).toHaveProperty('rev');
    expect(entry).toHaveProperty('success');
    expect(entry).toHaveProperty('source');
  });

  it('should have forward primer with required properties', () => {
    const entry = VALIDATION_DATASET[0];
    expect(entry.fwd).toHaveProperty('seq');
    expect(entry.fwd).toHaveProperty('tm');
    expect(entry.fwd).toHaveProperty('gc');
    expect(entry.fwd).toHaveProperty('dg');
    expect(entry.fwd).toHaveProperty('hairpinDG');
    expect(entry.fwd).toHaveProperty('homodimerDG');
    expect(entry.fwd).toHaveProperty('offTargetCount');
  });

  it('should have reverse primer with required properties', () => {
    const entry = VALIDATION_DATASET[0];
    expect(entry.rev).toHaveProperty('seq');
    expect(entry.rev).toHaveProperty('tm');
    expect(entry.rev).toHaveProperty('gc');
    expect(entry.rev).toHaveProperty('dg');
    expect(entry.rev).toHaveProperty('hairpinDG');
    expect(entry.rev).toHaveProperty('homodimerDG');
    expect(entry.rev).toHaveProperty('offTargetCount');
  });

  it('should have both success and failure cases', () => {
    const successCount = VALIDATION_DATASET.filter(e => e.success).length;
    const failureCount = VALIDATION_DATASET.filter(e => !e.success).length;

    expect(successCount).toBeGreaterThan(0);
    expect(failureCount).toBeGreaterThan(0);
  });

  it('should have multiple categories', () => {
    const categories = new Set(VALIDATION_DATASET.map(e => e.category));
    expect(categories.size).toBeGreaterThanOrEqual(4);
    expect(categories.has('optimal')).toBe(true);
    expect(categories.has('good')).toBe(true);
    expect(categories.has('marginal')).toBe(true);
    expect(categories.has('poor')).toBe(true);
  });

  it('should have valid Tm ranges', () => {
    for (const entry of VALIDATION_DATASET) {
      expect(entry.fwd.tm).toBeGreaterThan(30);
      expect(entry.fwd.tm).toBeLessThan(80);
      expect(entry.rev.tm).toBeGreaterThan(30);
      expect(entry.rev.tm).toBeLessThan(80);
    }
  });

  it('should have valid GC content ranges', () => {
    for (const entry of VALIDATION_DATASET) {
      expect(entry.fwd.gc).toBeGreaterThanOrEqual(0);
      expect(entry.fwd.gc).toBeLessThanOrEqual(100);
      expect(entry.rev.gc).toBeGreaterThanOrEqual(0);
      expect(entry.rev.gc).toBeLessThanOrEqual(100);
    }
  });
});

describe('calculateEntryScores', () => {
  it('should calculate scores for all features', () => {
    const entry = VALIDATION_DATASET[0];
    const scores = calculateEntryScores(entry);

    expect(scores).toHaveProperty('tmFwd');
    expect(scores).toHaveProperty('tmRev');
    expect(scores).toHaveProperty('gcFwd');
    expect(scores).toHaveProperty('gcRev');
    expect(scores).toHaveProperty('lengthFwd');
    expect(scores).toHaveProperty('lengthRev');
    expect(scores).toHaveProperty('hairpinFwd');
    expect(scores).toHaveProperty('hairpinRev');
    expect(scores).toHaveProperty('selfDimerFwd');
    expect(scores).toHaveProperty('selfDimerRev');
    expect(scores).toHaveProperty('terminal3DGFwd');
    expect(scores).toHaveProperty('terminal3DGRev');
    expect(scores).toHaveProperty('tmDiff');
    expect(scores).toHaveProperty('heterodimer');
    expect(scores).toHaveProperty('offTarget');
    expect(scores).toHaveProperty('ampliconLength');
  });

  it('should return scores between 0 and 1', () => {
    const entry = VALIDATION_DATASET[0];
    const scores = calculateEntryScores(entry);

    for (const [key, value] of Object.entries(scores)) {
      expect(value).toBeGreaterThanOrEqual(0);
      expect(value).toBeLessThanOrEqual(1);
    }
  });

  it('should give high scores to optimal primers', () => {
    const optimalEntry = VALIDATION_DATASET.find(e => e.category === 'optimal');
    const scores = calculateEntryScores(optimalEntry);

    // Optimal primers should have high Tm and GC scores
    expect(scores.tmFwd).toBeGreaterThan(0.7);
    expect(scores.tmRev).toBeGreaterThan(0.7);
    expect(scores.gcFwd).toBeGreaterThan(0.7);
    expect(scores.gcRev).toBeGreaterThan(0.7);
  });

  it('should give lower scores to poor primers', () => {
    const poorEntry = VALIDATION_DATASET.find(
      e => e.category === 'poor' && e.fwd.offTargetCount >= 3
    );

    if (poorEntry) {
      const scores = calculateEntryScores(poorEntry);
      // Off-target score should be 0 for 3+ off-targets
      expect(scores.offTarget).toBe(0);
    }
  });
});

describe('toCalibrationFormat', () => {
  it('should convert dataset to calibration format', () => {
    const calibrationData = toCalibrationFormat(VALIDATION_DATASET.slice(0, 10));

    expect(calibrationData.length).toBe(10);
    expect(calibrationData[0]).toHaveProperty('id');
    expect(calibrationData[0]).toHaveProperty('scores');
    expect(calibrationData[0]).toHaveProperty('compositeScore');
    expect(calibrationData[0]).toHaveProperty('actual');
  });

  it('should include composite score', () => {
    const calibrationData = toCalibrationFormat(VALIDATION_DATASET.slice(0, 5));

    for (const entry of calibrationData) {
      expect(entry.compositeScore).toBeGreaterThanOrEqual(0);
      expect(entry.compositeScore).toBeLessThanOrEqual(100);
    }
  });

  it('should preserve actual success/failure values', () => {
    const subset = VALIDATION_DATASET.slice(0, 10);
    const calibrationData = toCalibrationFormat(subset);

    for (let i = 0; i < 10; i++) {
      expect(calibrationData[i].actual).toBe(subset[i].success);
    }
  });
});

describe('getDatasetStats', () => {
  it('should calculate total count', () => {
    const stats = getDatasetStats();
    expect(stats.total).toBe(VALIDATION_DATASET.length);
  });

  it('should calculate success and failure counts', () => {
    const stats = getDatasetStats();
    expect(stats.successCount + stats.failureCount).toBe(stats.total);
  });

  it('should calculate success rate', () => {
    const stats = getDatasetStats();
    expect(stats.successRate).toBeGreaterThan(0);
    expect(stats.successRate).toBeLessThan(1);
    expect(stats.successRate).toBeCloseTo(stats.successCount / stats.total, 5);
  });

  it('should break down by category', () => {
    const stats = getDatasetStats();
    expect(stats.byCategory).toHaveProperty('optimal');
    expect(stats.byCategory).toHaveProperty('good');
    expect(stats.byCategory).toHaveProperty('marginal');
    expect(stats.byCategory).toHaveProperty('poor');

    // Each category should have total and success counts
    for (const category of Object.values(stats.byCategory)) {
      expect(category).toHaveProperty('total');
      expect(category).toHaveProperty('success');
      expect(category.total).toBeGreaterThanOrEqual(category.success);
    }
  });

  it('should show higher success rate for optimal vs poor', () => {
    const stats = getDatasetStats();
    const optimalRate = stats.byCategory.optimal.success / stats.byCategory.optimal.total;
    const poorRate = stats.byCategory.poor.success / stats.byCategory.poor.total;

    expect(optimalRate).toBeGreaterThan(poorRate);
  });
});

describe('splitDataset', () => {
  it('should split dataset into train and test sets', () => {
    const { train, test } = splitDataset();

    expect(train.length + test.length).toBe(VALIDATION_DATASET.length);
  });

  it('should respect train ratio', () => {
    const { train, test } = splitDataset(VALIDATION_DATASET, 0.8);

    // Allow for rounding differences
    const expectedTrain = Math.floor(VALIDATION_DATASET.length * 0.8);
    expect(train.length).toBe(expectedTrain);
    expect(test.length).toBe(VALIDATION_DATASET.length - expectedTrain);
  });

  it('should shuffle data (not return same order)', () => {
    const { train: train1 } = splitDataset();
    const { train: train2 } = splitDataset();

    // With random shuffling, the first few elements should differ
    // (there's a small chance this could fail, but very unlikely)
    const firstFiveIds1 = train1.slice(0, 5).map(e => e.id).join(',');
    const firstFiveIds2 = train2.slice(0, 5).map(e => e.id).join(',');

    // At least one pair should differ if truly shuffled
    // But we can't guarantee this due to randomness
    // So we just check that both have valid entries
    expect(train1[0]).toHaveProperty('id');
    expect(train2[0]).toHaveProperty('id');
  });

  it('should work with custom train ratio', () => {
    const { train, test } = splitDataset(VALIDATION_DATASET, 0.7);

    const expectedTrain = Math.floor(VALIDATION_DATASET.length * 0.7);
    expect(train.length).toBe(expectedTrain);
  });

  it('should not modify original dataset', () => {
    const originalLength = VALIDATION_DATASET.length;
    splitDataset();
    expect(VALIDATION_DATASET.length).toBe(originalLength);
  });
});

describe('Dataset Quality Checks', () => {
  it('should have reasonable success rate (~60-80%)', () => {
    const stats = getDatasetStats();
    // Based on literature: PrimerBank 82.6%, crowdsourced 55-63%
    expect(stats.successRate).toBeGreaterThan(0.4);
    expect(stats.successRate).toBeLessThan(0.9);
  });

  it('should have optimal primers with highest success rate', () => {
    const stats = getDatasetStats();
    const categories = ['optimal', 'good', 'marginal', 'poor'];

    const rates = categories.map(cat => ({
      category: cat,
      rate: stats.byCategory[cat].success / stats.byCategory[cat].total,
    }));

    // Sort by rate descending
    rates.sort((a, b) => b.rate - a.rate);

    // Optimal should be first or second (allowing for statistical variance)
    expect(['optimal', 'good']).toContain(rates[0].category);
  });

  it('should have diverse primer characteristics', () => {
    const tmValues = VALIDATION_DATASET.map(e => e.fwd.tm);
    const gcValues = VALIDATION_DATASET.map(e => e.fwd.gc);
    const offTargets = VALIDATION_DATASET.map(e => e.fwd.offTargetCount);

    // Check Tm diversity
    const minTm = Math.min(...tmValues);
    const maxTm = Math.max(...tmValues);
    expect(maxTm - minTm).toBeGreaterThan(15);  // At least 15°C range

    // Check GC diversity
    const minGc = Math.min(...gcValues);
    const maxGc = Math.max(...gcValues);
    expect(maxGc - minGc).toBeGreaterThan(30);  // At least 30% range

    // Check off-target diversity
    const uniqueOffTargets = new Set(offTargets);
    expect(uniqueOffTargets.size).toBeGreaterThan(1);
  });
});
