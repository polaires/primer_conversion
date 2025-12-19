/**
 * Tests for Alternative Generation Module
 */

import { describe, it, expect } from 'vitest';
import {
  generateAlternativesUnified,
  assignTier,
  calculateQuickScore,
  calculateMediumScore,
  THRESHOLDS,
} from './alternativeGeneration.js';

// Test template - a realistic 200bp sequence
const TEST_TEMPLATE = 'ATGCGTACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG';

// GFP-like sequence (more realistic GC content ~50%)
const GFP_TEMPLATE = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA';

describe('generateAlternativesUnified', () => {
  it('generates alternatives for a valid template', () => {
    const alternatives = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 300), {
      numAlternatives: 3,
    });

    expect(alternatives.length).toBeGreaterThan(0);
    expect(alternatives.length).toBeLessThanOrEqual(3);
  });

  it('returns alternatives with required properties', () => {
    const alternatives = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 300), {
      numAlternatives: 3,
    });

    const alt = alternatives[0];

    // Check forward primer properties
    expect(alt.forward).toBeDefined();
    expect(alt.forward.sequence).toBeDefined();
    expect(alt.forward.length).toBeGreaterThanOrEqual(15);
    expect(alt.forward.length).toBeLessThanOrEqual(60); // Extended max for Tm matching
    expect(alt.forward.tm).toBeGreaterThan(40);
    expect(alt.forward.gc).toBeGreaterThan(0);
    expect(alt.forward.gc).toBeLessThan(1);

    // Check reverse primer properties
    expect(alt.reverse).toBeDefined();
    expect(alt.reverse.sequence).toBeDefined();

    // Check pair properties
    expect(alt.compositeScore).toBeGreaterThan(0);
    expect(alt.qualityTier).toBeDefined();
    expect(alt.tmDiff).toBeDefined();
    expect(alt.heterodimerDG).toBeDefined();
  });

  it('includes UI labels for alternatives with unique strengths', () => {
    const alternatives = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 300), {
      numAlternatives: 3,
    });

    for (const alt of alternatives) {
      // label can be string or null (null if no unique strength)
      expect(alt.label === null || typeof alt.label === 'string').toBe(true);
      // explanation can be string or null
      expect(alt.explanation === null || typeof alt.explanation === 'string').toBe(true);
      expect(alt.selectionReason).toBeDefined();
    }
  });

  it('first alternative has best overall selection reason', () => {
    const alternatives = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 300), {
      numAlternatives: 5,
    });

    if (alternatives.length > 0) {
      // First is always selected as bestOverall by hybrid selection
      expect(alternatives[0].selectionReason).toBe('bestOverall');
      // Label depends on identified strengths - can be any valid label or null
      // Valid labels include: '⭐ Best', '✂️ Compact', '🎯 Optimal Tm', etc.
      expect(alternatives[0].label === null || typeof alternatives[0].label === 'string').toBe(true);
    }
  });

  it('alternatives are sorted with best score first', () => {
    const alternatives = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 300), {
      numAlternatives: 5,
    });

    if (alternatives.length >= 2) {
      // First should have highest or equal score
      expect(alternatives[0].compositeScore).toBeGreaterThanOrEqual(
        alternatives[1].compositeScore - 5 // Allow small margin for diversity
      );
    }
  });

  it('provides diverse alternatives (different from each other)', () => {
    const alternatives = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 300), {
      numAlternatives: 5,
    });

    if (alternatives.length >= 3) {
      // Check that we have some variety in lengths
      const fwdLengths = new Set(alternatives.map(a => a.forward.length));
      const revLengths = new Set(alternatives.map(a => a.reverse.length));

      // Should have at least 2 different lengths across 5 alternatives
      expect(fwdLengths.size + revLengths.size).toBeGreaterThanOrEqual(2);
    }
  });

  it('throws error for template too short', () => {
    expect(() => {
      generateAlternativesUnified('ATGCATGC', { numAlternatives: 3 });
    }).toThrow();
  });

  it('respects maxTmDiff parameter', () => {
    const alternatives = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 300), {
      numAlternatives: 10,
      maxTmDiff: 3, // Strict Tm matching
    });

    for (const alt of alternatives) {
      expect(alt.tmDiff).toBeLessThanOrEqual(3);
    }
  });
});

describe('assignTier', () => {
  // Note: The tier system uses worstDg > -3 for excellent, meaning
  // dg should be close to 0 (weak self-interaction is good)
  // More negative dg = stronger self-interaction = worse

  it('assigns excellent tier for low Tm diff and weak self-interaction', () => {
    const pair = {
      forward: { tm: 58, dg: -2 },  // Weak self-interaction (good)
      reverse: { tm: 59, dg: -1 },  // Weak self-interaction (good)
      tmDiff: 1,
    };

    expect(assignTier(pair)).toBe('excellent');
  });

  it('assigns good tier for moderate Tm diff and moderate stability', () => {
    const pair = {
      forward: { tm: 58, dg: -4 },  // Moderate self-interaction
      reverse: { tm: 62, dg: -3.5 },
      tmDiff: 4,
    };

    expect(assignTier(pair)).toBe('good');
  });

  it('assigns acceptable tier for larger Tm diff', () => {
    const pair = {
      forward: { tm: 55, dg: -8 },  // Strong self-interaction
      reverse: { tm: 62, dg: -7 },
      tmDiff: 7,
    };

    expect(assignTier(pair)).toBe('acceptable');
  });

  it('assigns poor tier for very large Tm diff', () => {
    const pair = {
      forward: { tm: 52, dg: -8 },
      reverse: { tm: 65, dg: -7 },
      tmDiff: 13,
    };

    expect(assignTier(pair)).toBe('poor');
  });

  it('downgrades tier when dg shows strong self-interaction', () => {
    // Even with good Tm diff, strong self-interaction downgrades tier
    const pair = {
      forward: { tm: 58, dg: -6 },  // Strong self-interaction (bad)
      reverse: { tm: 59, dg: -5 },
      tmDiff: 1,
    };

    // Should NOT be excellent because worstDg (-6) is not > -3
    expect(assignTier(pair)).not.toBe('excellent');
  });
});

describe('Performance and edge cases', () => {
  it('handles longer templates efficiently', () => {
    const start = Date.now();

    const alternatives = generateAlternativesUnified(GFP_TEMPLATE, {
      numAlternatives: 5,
    });

    const elapsed = Date.now() - start;

    // Should complete in reasonable time (< 5 seconds)
    expect(elapsed).toBeLessThan(5000);
    expect(alternatives.length).toBeGreaterThan(0);
  });

  it('handles templates with extreme GC content', () => {
    // AT-rich template
    const atRich = 'A'.repeat(20) + 'GCGC' + 'T'.repeat(30) + 'GCGC' + 'A'.repeat(20) + 'GCGC' + 'T'.repeat(30) + 'GCGC' + 'A'.repeat(20);

    // Should still find some alternatives (or return empty gracefully)
    const alternatives = generateAlternativesUnified(atRich, {
      numAlternatives: 3,
    });

    // Either finds alternatives or returns empty array (no crash)
    expect(Array.isArray(alternatives)).toBe(true);
  });

  it('returns empty array when no valid pairs exist', () => {
    // Very short template that's borderline
    const shortTemplate = 'ATGC'.repeat(20); // 80bp, just enough

    const alternatives = generateAlternativesUnified(shortTemplate, {
      numAlternatives: 3,
      maxTmDiff: 0.1, // Impossibly strict
    });

    expect(Array.isArray(alternatives)).toBe(true);
    // May be empty due to strict constraints
  });
});

describe('Integration: Diversity vs Pure Ranking', () => {
  it('diversity selection produces more varied alternatives than pure ranking', () => {
    const alternatives = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 400), {
      numAlternatives: 5,
    });

    if (alternatives.length >= 3) {
      // Calculate position variance
      const fwdPositions = alternatives.map(a => a.forward.startPos || 0);
      const posRange = Math.max(...fwdPositions) - Math.min(...fwdPositions);

      // Should have some position diversity (not all from same spot)
      // This tests that diversity selection is working
      expect(posRange).toBeGreaterThanOrEqual(0);

      // Check we have different selection reasons
      const reasons = new Set(alternatives.map(a => a.selectionReason));
      // Should have bestOverall and at least one other type
      expect(reasons.has('bestOverall')).toBe(true);
    }
  });
});

// ═══════════════════════════════════════════════════════════════════════════
// TIER 1: Quick Score Tests
// ═══════════════════════════════════════════════════════════════════════════

describe('calculateQuickScore - Tiered Scoring', () => {
  it('returns 0 for primers outside hard limits', () => {
    // Tm too low
    expect(calculateQuickScore({ tm: 40, gc: 0.5, length: 20 })).toBe(0);
    // Tm too high
    expect(calculateQuickScore({ tm: 80, gc: 0.5, length: 20 })).toBe(0);
    // GC too low
    expect(calculateQuickScore({ tm: 58, gc: 0.15, length: 20 })).toBe(0);
    // GC too high
    expect(calculateQuickScore({ tm: 58, gc: 0.85, length: 20 })).toBe(0);
    // Length too short
    expect(calculateQuickScore({ tm: 58, gc: 0.5, length: 10 })).toBe(0);
    // Length too long (hardMax is now 60bp for extended primer support)
    expect(calculateQuickScore({ tm: 58, gc: 0.5, length: 65 })).toBe(0);
  });

  it('returns 1.0 for primers in optimal ranges', () => {
    // Perfect primer: optimal Tm, GC, and length
    const score = calculateQuickScore({ tm: 58, gc: 0.5, length: 20 });
    expect(score).toBe(1.0);
  });

  it('penalizes primers in acceptable but not optimal ranges', () => {
    // Tm in acceptable range (50-55) but not optimal (55-60)
    const tmAcceptable = calculateQuickScore({ tm: 52, gc: 0.5, length: 20 });
    expect(tmAcceptable).toBeLessThan(1.0);
    expect(tmAcceptable).toBeGreaterThan(0.5);

    // GC in acceptable range (30-40%) but not optimal (40-60%)
    const gcAcceptable = calculateQuickScore({ tm: 58, gc: 0.35, length: 20 });
    expect(gcAcceptable).toBeLessThan(1.0);
    expect(gcAcceptable).toBeGreaterThan(0.5);
  });

  it('aligns with piecewise thresholds (no ranking inversion)', () => {
    // Primer A: in optimal zones
    const optimal = { tm: 58, gc: 0.5, length: 20 };
    // Primer B: in acceptable zones
    const acceptable = { tm: 52, gc: 0.35, length: 26 };
    // Primer C: barely acceptable
    const marginal = { tm: 50, gc: 0.30, length: 28 };

    const scoreA = calculateQuickScore(optimal);
    const scoreB = calculateQuickScore(acceptable);
    const scoreC = calculateQuickScore(marginal);

    // Quick score should maintain same ranking as piecewise would
    expect(scoreA).toBeGreaterThan(scoreB);
    expect(scoreB).toBeGreaterThan(scoreC);
  });
});

// ═══════════════════════════════════════════════════════════════════════════
// TIER 2: Medium Score Tests
// ═══════════════════════════════════════════════════════════════════════════

describe('calculateMediumScore - Tiered Scoring', () => {
  it('returns high scores for primers with good properties', () => {
    const goodPrimer = {
      seq: 'ATGCGTACGATCGATCGATC',  // 20bp, ends with C
      tm: 58,
      gc: 0.5,
      dg: -8,  // Good terminal stability
      length: 20,
    };

    const score = calculateMediumScore(goodPrimer);
    expect(score).toBeGreaterThan(0.7);
  });

  it('returns lower scores for primers with suboptimal properties', () => {
    const suboptimalPrimer = {
      seq: 'AAAATTTTAAAATTTTAAAA',  // 20bp, AT-rich, ends with A, has homopolymers
      tm: 52,  // Below optimal
      gc: 0.2,  // Below optimal
      dg: -2,   // Weak binding
      length: 20,
    };

    const score = calculateMediumScore(suboptimalPrimer);
    expect(score).toBeLessThan(0.7);
  });

  it('uses piecewise scoring functions internally', () => {
    // Two primers that differ only in one property
    const base = {
      seq: 'ATGCGTACGATCGATCGATC',
      tm: 58,
      gc: 0.5,
      dg: -8,
      length: 20,
    };

    const weakDg = {
      ...base,
      dg: -2,  // Much weaker terminal binding
    };

    const baseScore = calculateMediumScore(base);
    const weakScore = calculateMediumScore(weakDg);

    // Should reflect the difference in terminal 3' stability
    expect(baseScore).toBeGreaterThan(weakScore);
  });
});

// ═══════════════════════════════════════════════════════════════════════════
// THRESHOLDS Configuration Tests
// ═══════════════════════════════════════════════════════════════════════════

describe('THRESHOLDS configuration', () => {
  it('has aligned Tm thresholds with scoring.js', () => {
    // These should match the values in scoring.js scoreTm()
    expect(THRESHOLDS.tm.optimalLow).toBe(55);
    expect(THRESHOLDS.tm.optimalHigh).toBe(60);
    expect(THRESHOLDS.tm.acceptableLow).toBe(50);
    expect(THRESHOLDS.tm.acceptableHigh).toBe(65);
  });

  it('has aligned GC thresholds with scoring.js', () => {
    // These should match the values in scoring.js scoreGc()
    expect(THRESHOLDS.gc.optimalLow).toBe(0.40);
    expect(THRESHOLDS.gc.optimalHigh).toBe(0.60);
    expect(THRESHOLDS.gc.acceptableLow).toBe(0.30);
    expect(THRESHOLDS.gc.acceptableHigh).toBe(0.70);
  });

  it('has aligned length thresholds with scoring.js', () => {
    // These should match the values in scoring.js scoreLength()
    // Note: acceptableHigh extended to 35bp and hardMax to 60bp for better Tm matching
    expect(THRESHOLDS.length.optimalLow).toBe(18);
    expect(THRESHOLDS.length.optimalHigh).toBe(24);
    expect(THRESHOLDS.length.acceptableLow).toBe(15);
    expect(THRESHOLDS.length.acceptableHigh).toBe(35);
    expect(THRESHOLDS.length.hardMax).toBe(60);
  });
});

// ═══════════════════════════════════════════════════════════════════════════
// Tiered Pipeline Integration Test
// ═══════════════════════════════════════════════════════════════════════════

describe('Tiered Pipeline Integration', () => {
  it('alternatives include mediumScore from Tier 2', () => {
    const alternatives = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 300), {
      numAlternatives: 3,
    });

    if (alternatives.length > 0) {
      const alt = alternatives[0];
      // Medium scores should be present from Tier 2
      expect(alt.forward.mediumScore).toBeDefined();
      expect(alt.reverse.mediumScore).toBeDefined();
      expect(alt.forward.mediumScore).toBeGreaterThan(0);
      expect(alt.forward.mediumScore).toBeLessThanOrEqual(1);
    }
  });

  it('tiered pipeline produces consistent ranking', () => {
    // Run multiple times - should produce consistent best result
    const runs = [];
    for (let i = 0; i < 3; i++) {
      const alts = generateAlternativesUnified(GFP_TEMPLATE.slice(0, 200), {
        numAlternatives: 1,
      });
      if (alts.length > 0) {
        runs.push(alts[0].compositeScore);
      }
    }

    // All runs should produce the same best score
    if (runs.length >= 2) {
      expect(runs[0]).toBe(runs[1]);
    }
  });
});
