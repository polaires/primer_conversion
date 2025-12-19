/**
 * Tests for the primers library
 * These tests verify the JavaScript port matches the Python implementation
 */

import { describe, it, expect } from 'vitest';
import { primers, score, LEN_MIN, LEN_MAX } from './primers.js';
import { gcCache, tmCache } from './tm.js';
import { dgCache } from './fold.js';
import { offTargets } from './offTargets.js';

// Helper to reverse complement
function rc(seq) {
  const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
  return seq.split('').reverse().map(c => comp[c]).join('');
}

describe('Primers', () => {
  it('should create primers without additional sequence', () => {
    const seq = 'CTACTAATAGCACACACGGGGCAATACCAGCACAAGCTAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTAAGGCACGCCAATTCGCACTGTCAGGGTCACTTGGGTGTTTTGCACTACCGTCAGGTACGCTAGTATGCGTTCTTCCTTCCAGAGGTATGTGGCTGCGTGGTCAAAAGTGCGGCATTCGTATTTGCTCCTCGTGTTTACTCTCACAAACTTGACCTGGAGATCAAGGAGATGCTTCTTGTGGAACTGGACAACGCATCAACGCAACGGATCTACGTTACAGCGT';

    const [p1, p2] = primers(seq);

    expect(p1).toBeDefined();
    expect(p2).toBeDefined();
    expect(p1.seq).toBe(p1.seq.toUpperCase());
    expect(p2.seq).toBe(p2.seq.toUpperCase());
    expect(seq.includes(p1.seq)).toBe(true);
    expect(seq.includes(rc(p2.seq))).toBe(true);
    expect(seq.startsWith(p1.seq)).toBe(true);
    expect(seq.endsWith(rc(p2.seq))).toBe(true);
    expect(p1.len).toBeGreaterThanOrEqual(LEN_MIN);
    expect(p1.len).toBeLessThanOrEqual(LEN_MAX);
    expect(p2.len).toBeGreaterThanOrEqual(LEN_MIN);
    expect(p2.len).toBeLessThanOrEqual(LEN_MAX);
    expect(p1.tm).toBeTruthy();
    expect(p2.tm).toBeTruthy();
    expect(p1.tm).toBe(p1.tmTotal);
    expect(p2.tm).toBe(p2.tmTotal);
    expect(p1.gc).toBeTruthy();
    expect(p2.gc).toBeTruthy();
    expect(p1.fwd).toBe(true);
    expect(p2.fwd).toBe(false);
    expect(p1.penalty).toBeTruthy();
    expect(p2.penalty).toBeTruthy();
  });

  it('should create primers with additional sequence on the left (fixed)', () => {
    const seq = 'CTACTAATAGCACACACGGGGCAATACCAGCACAAGCTAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTAAGGCACGCCAATTCGCACTGTCAGGGTCACTTGGGTGTTTTGCACTACCGTCAGGTACGCTAGTATGCGTTCTTCCTTCCAGAGGTATGTGGCTGCGTGGTCAAAAGTGCGGCATTCGTATTTGCTCCTCGTGTTTACTCTCACAAACTTGACCTGGAGATCAAGGAGATGCTTCTTGTGGAACTGGACAACGCATCAACGCAACGGATCTACGTTACAGCGT';
    const addFwd = 'GGTCTC';

    const [p1, p2] = primers(seq, { addFwd });

    expect(p1.seq.startsWith(addFwd)).toBe(true);
    expect((addFwd + seq).startsWith(p1.seq)).toBe(true);
    expect(p2.tm).toBe(p2.tmTotal);
    expect(p1.tm).not.toBe(p1.tmTotal);
  });

  it('should create primers with additional sequence on the right (fixed)', () => {
    const seq = 'CTACTAATAGCACACACGGGGCAATACCAGCACAAGCTAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTAAGGCACGCCAATTCGCACTGTCAGGGTCACTTGGGTGTTTTGCACTACCGTCAGGTACGCTAGTATGCGTTCTTCCTTCCAGAGGTATGTGGCTGCGTGGTCAAAAGTGCGGCATTCGTATTTGCTCCTCGTGTTTACTCTCACAAACTTGACCTGGAGATCAAGGAGATGCTTCTTGTGGAACTGGACAACGCATCAACGCAACGGATCTACGTTACAGCGT';
    const addRev = 'GGTCTC';

    const [p1, p2] = primers(seq, { addRev });

    expect(p2.seq.startsWith(addRev)).toBe(true);
    expect((addRev + rc(seq)).startsWith(p2.seq)).toBe(true);
    expect(p1.tm).toBe(p1.tmTotal);
    expect(p2.tm).not.toBe(p2.tmTotal);
  });

  it('should throw an error on an invalid input sequence', () => {
    const seq = 'AUGGCCUUUCUCGGGGGGCCUUGGUCUCUGAAAC';

    expect(() => primers(seq)).toThrow();
  });

  it('should create primers given a parent with diff-case sequence', () => {
    const [p1, p2] = primers('AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA', {
      offtargetCheck: 'ggaattacgtAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAggaccagttacagga',
    });

    expect(p1).toBeDefined();
    expect(p2).toBeDefined();
  });
});

describe('Score', () => {
  it('should score an existing pair of primers', () => {
    const [fwd, rev] = score(
      'GGTCTCAATGAGACAATA',
      'TTTCGTATGCTGACCTAG',
      '',
      'AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA'
    );

    expect(fwd.penalty).toBeGreaterThan(15);
    expect(fwd.penalty).toBeLessThan(25);
    expect(rev.penalty).toBeGreaterThan(3);
    expect(rev.penalty).toBeLessThan(10);
  });

  it('should score only a FWD primer', () => {
    const [fwd, rev] = score(
      'GGTCTCAATGAGACAATAGCACACAC',
      '',
      '',
      'AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA'
    );

    expect(fwd.penalty).toBeGreaterThan(8);
    expect(fwd.penalty).toBeLessThan(15);
    expect(rev).toBeNull();
  });

  it('should throw on short primers', () => {
    expect(() => score('ACGACTAC')).toThrow(); // too short fwd
    expect(() => score('ACGACTACGACTACGATC', 'GACTACG')).toThrow(); // too short rev
  });
});

describe('Off-targets', () => {
  it('should find and cache offtarget binding sites', () => {
    // GTGGCTAGCC is one by removed from GTGGCTAGGC in seq
    const parent = 'CTGACTCTACTTGGAAATGTGGCTAGGCCTTTGCCCACGCACCTGATCGGTCCTGTGGCTAGCCTCGTTTGCTTTTTAGGACCGGATGAACTACAGAGCATTGCAAGAATC';
    const seq = 'CTGACTCTACTTGGAAATGTGGCTAGGCCTT';

    const ot = offTargets(seq, parent);

    expect(ot[0]).toBe(0);
    expect(ot.length).toBe(seq.length);
    expect(ot.some(o => o > 0)).toBe(true);
  });
});

describe('GC Cache', () => {
  it('should calculate GC ratios for subsequences', () => {
    const seq = 'ATGCATGC';
    const gc = gcCache(seq);

    expect(gc).toBeDefined();
    expect(gc.length).toBe(seq.length);
    expect(gc[0][seq.length - 1]).toBe(0.5);
  });
});

describe('TM Cache', () => {
  it('should calculate melting temperatures for subsequences', () => {
    const seq = 'ATGCATGCATGCATGC';
    const tm = tmCache(seq);

    expect(tm).toBeDefined();
    expect(tm.length).toBe(seq.length);
    expect(typeof tm[0][seq.length - 1]).toBe('number');
    expect(tm[0][seq.length - 1]).toBeGreaterThan(0);
  });
});

describe('DG Cache', () => {
  it('should calculate free energies for subsequences', () => {
    const seq = 'ATGCATGCATGC';
    const dg = dgCache(seq);

    expect(dg).toBeDefined();
    expect(dg.length).toBe(seq.length);
    expect(typeof dg[0][seq.length - 1]).toBe('number');
  });
});

describe('Equilibrium Scoring Integration', () => {
  it('should include equilibrium metrics when scoring primer pairs', () => {
    const fwd = 'ATGCATGCATGCATGCATGC';  // 20bp primer
    const rev = 'GCATGCATGCATGCATGCAT';  // 20bp primer
    // Template: fwd binds at start, rev (RC) binds at end
    const template = fwd + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + rc(rev);

    const [p1, p2] = score(fwd, rev, template, template);

    // Check that equilibrium metrics are present
    expect(p1.scoring.equilibriumEfficiency).toBeDefined();
    expect(p1.scoring.equilibriumScore).toBeDefined();
    expect(p1.scoring.equilibriumLosses).toBeDefined();

    expect(p2.scoring.equilibriumEfficiency).toBeDefined();
    expect(p2.scoring.equilibriumScore).toBeDefined();
    expect(p2.scoring.equilibriumLosses).toBeDefined();

    // Efficiency should be between 0 and 1
    expect(p1.scoring.equilibriumEfficiency).toBeGreaterThanOrEqual(0);
    expect(p1.scoring.equilibriumEfficiency).toBeLessThanOrEqual(1);
    expect(p2.scoring.equilibriumEfficiency).toBeGreaterThanOrEqual(0);
    expect(p2.scoring.equilibriumEfficiency).toBeLessThanOrEqual(1);

    // Score should be between 0 and 100
    expect(p1.scoring.equilibriumScore).toBeGreaterThanOrEqual(0);
    expect(p1.scoring.equilibriumScore).toBeLessThanOrEqual(100);

    // Loss breakdown should have expected keys
    expect(p1.scoring.equilibriumLosses).toHaveProperty('hairpin');
    expect(p1.scoring.equilibriumLosses).toHaveProperty('homodimer');
    expect(p1.scoring.equilibriumLosses).toHaveProperty('heterodimer');
    expect(p1.scoring.equilibriumLosses).toHaveProperty('offTarget');
    expect(p1.scoring.equilibriumLosses).toHaveProperty('free');
  });

  it('should allow disabling equilibrium scoring', () => {
    const fwd = 'ATGCATGCATGCATGCATGC';
    const rev = 'GCATGCATGCATGCATGCAT';
    const template = fwd + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + rc(rev);

    const [p1, p2] = score(fwd, rev, template, template, { includeEquilibrium: false });

    // Without equilibrium, these should be null
    expect(p1.scoring.equilibriumEfficiency).toBeNull();
    expect(p1.scoring.equilibriumScore).toBeNull();
    expect(p1.scoring.equilibriumLosses).toBeNull();
  });

  it('should include equilibrium metrics in dict() output', () => {
    const fwd = 'ATGCATGCATGCATGCATGC';
    const rev = 'GCATGCATGCATGCATGCAT';
    const template = fwd + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + rc(rev);

    const [p1] = score(fwd, rev, template, template);
    const dict = p1.dict();

    expect(dict.scoring.equilibrium_efficiency).toBeDefined();
    expect(dict.scoring.equilibrium_score).toBeDefined();
    expect(dict.scoring.equilibrium_losses).toBeDefined();
  });
});

describe('Piecewise Logistic Scoring Integration', () => {
  it('should include piecewise scores when scoring primer pairs', () => {
    const fwd = 'ATGCATGCATGCATGCATGC';  // 20bp primer
    const rev = 'GCATGCATGCATGCATGCAT';  // 20bp primer
    const template = fwd + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + rc(rev);

    const [p1, p2] = score(fwd, rev, template, template);

    // Check that piecewise scores are present
    expect(p1.scoring.piecewiseScores).toBeDefined();
    expect(p1.scoring.compositeScore).toBeDefined();
    expect(p1.scoring.qualityTier).toBeDefined();

    expect(p2.scoring.piecewiseScores).toBeDefined();
    expect(p2.scoring.compositeScore).toBeDefined();
    expect(p2.scoring.qualityTier).toBeDefined();
  });

  it('should have individual piecewise scores between 0 and 1', () => {
    const fwd = 'ATGCATGCATGCATGCATGC';
    const rev = 'GCATGCATGCATGCATGCAT';
    const template = fwd + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + rc(rev);

    const [p1] = score(fwd, rev, template, template);

    // Check individual scores are in valid range
    expect(p1.scoring.piecewiseScores.tm).toBeGreaterThanOrEqual(0);
    expect(p1.scoring.piecewiseScores.tm).toBeLessThanOrEqual(1);
    expect(p1.scoring.piecewiseScores.gc).toBeGreaterThanOrEqual(0);
    expect(p1.scoring.piecewiseScores.gc).toBeLessThanOrEqual(1);
    expect(p1.scoring.piecewiseScores.length).toBeGreaterThanOrEqual(0);
    expect(p1.scoring.piecewiseScores.length).toBeLessThanOrEqual(1);
    expect(p1.scoring.piecewiseScores.hairpin).toBeGreaterThanOrEqual(0);
    expect(p1.scoring.piecewiseScores.hairpin).toBeLessThanOrEqual(1);
  });

  it('should have composite score between 0 and 100', () => {
    const fwd = 'ATGCATGCATGCATGCATGC';
    const rev = 'GCATGCATGCATGCATGCAT';
    const template = fwd + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + rc(rev);

    const [p1] = score(fwd, rev, template, template);

    expect(p1.scoring.compositeScore).toBeGreaterThanOrEqual(0);
    expect(p1.scoring.compositeScore).toBeLessThanOrEqual(100);
  });

  it('should have valid quality tier', () => {
    const fwd = 'ATGCATGCATGCATGCATGC';
    const rev = 'GCATGCATGCATGCATGCAT';
    const template = fwd + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + rc(rev);

    const [p1] = score(fwd, rev, template, template);

    const validTiers = ['excellent', 'good', 'acceptable', 'marginal', 'poor'];
    expect(validTiers).toContain(p1.scoring.qualityTier);
  });

  it('should include piecewise scores in dict() output', () => {
    const fwd = 'ATGCATGCATGCATGCATGC';
    const rev = 'GCATGCATGCATGCATGCAT';
    const template = fwd + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + rc(rev);

    const [p1] = score(fwd, rev, template, template);
    const dict = p1.dict();

    expect(dict.scoring.piecewise_scores).toBeDefined();
    expect(dict.scoring.composite_score).toBeDefined();
    expect(dict.scoring.quality_tier).toBeDefined();
  });

  it('should work for single primer scoring', () => {
    const fwd = 'ATGCATGCATGCATGCATGC';
    const template = fwd + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA';

    const [p1, p2] = score(fwd, '', template, template);

    expect(p1.scoring.piecewiseScores).toBeDefined();
    expect(p1.scoring.compositeScore).toBeDefined();
    expect(p2).toBeNull();
  });

  it('should give higher scores to well-designed primers', () => {
    // A well-designed primer: optimal length, balanced GC, varied sequence (no dimers)
    // Using non-repetitive sequence to avoid dimer issues
    const goodFwd = 'ACGTACGATCGTAGCATGC';  // 50% GC, 19bp, varied
    const goodRev = 'GCTACGATCGTAGCATGCA';  // 50% GC, 19bp, varied
    const spacer = 'CTGATCGATCGATCGATCGATCGATCGATCGATCGATCG';  // Generic spacer
    const goodTemplate = goodFwd + spacer + rc(goodRev);

    // A poorly-designed primer with multiple issues:
    // - Very low GC content
    // - Long homopolymer runs causing sequencing issues
    const badFwd = 'AAAAAAAAATATATATATAT';  // 10% GC, poly-A + poly-AT
    const badRev = 'TTTTTTTTTATATATATAT';   // 5% GC, poly-T + poly-AT
    const badTemplate = badFwd + spacer + rc(badRev);

    const [good] = score(goodFwd, goodRev, goodTemplate, goodTemplate);
    const [bad] = score(badFwd, badRev, badTemplate, badTemplate);

    // The good primer should have better individual feature scores
    // Check that GC content scoring is working correctly
    expect(good.scoring.piecewiseScores.gc).toBeGreaterThan(bad.scoring.piecewiseScores.gc);

    // Check that homopolymer scoring is working correctly
    expect(good.scoring.piecewiseScores.homopolymer).toBeGreaterThanOrEqual(bad.scoring.piecewiseScores.homopolymer);
  });

  it('should use joint Tm optimization for sequences with mismatched GC content', () => {
    // Create a sequence with GC-rich forward region (~70% GC) and AT-rich reverse region (~30% GC)
    // This is the scenario where the old algorithm would select primers based on individual
    // optimal Tm (62°C) rather than Tm MATCHING
    //
    // Forward region Tm range: ~70-83°C (high GC)
    // Reverse region Tm range: ~48-65°C (low GC)
    // Overlap exists around 65-70°C, but naive algorithm would miss it

    // Forward: ~70% GC (high Tm)
    const gcRichStart = 'GCGCATGCGCATGCGCATGCGCATGCGCATGCGCATGCGCATGCGCATGCGCATGCGC';
    // Middle: 50% GC
    const middle = 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC';
    // Reverse: ~30% GC (low Tm)
    const atRichEnd = 'ATATGCATATGCATATGCATATGCATATGCATATGCATATGCATATGCATATGCATAT';

    const seq = gcRichStart + middle + atRichEnd;

    const [p1, p2] = primers(seq);

    expect(p1).toBeDefined();
    expect(p2).toBeDefined();

    const tmDiff = Math.abs(p1.tm - p2.tm);

    // Before fix (naive top-10 by individual penalty):
    // - Forward primers selected: ~62°C optimal → 17-22bp at ~71-77°C
    // - Reverse primers selected: ~62°C optimal → 23-28bp at ~61-63°C
    // - Best pairing from those: ~10-15°C difference
    //
    // After fix (joint Tm optimization):
    // - Searches for overlapping Tm buckets
    // - Finds shorter forward primers (15-17bp at ~70-71°C)
    // - Matches with extended reverse primers (29-32bp at ~65°C)
    // - Achieves ~5-8°C Tm difference
    expect(tmDiff).toBeLessThan(10); // Joint optimization should achieve <10°C
  });

  it('should prioritize Tm matching over individual optimal Tm', () => {
    // A sequence where matching Tm requires non-optimal individual primer Tms
    // Forward: ~65% GC, Reverse: ~35% GC
    // The algorithm should sacrifice individual "optimal" Tm (62°C) for better pair matching

    const gcModerateStart = 'GCATGCATGCGCGCATGCATGCGCGCATGCATGCGCGCATGCATGCGCGCATGC';
    const middle = 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC';
    const atModerateEnd = 'ATATATGCATATATATGCATATATATGCATATATATGCATATATATGCATATAT';

    const seq = gcModerateStart + middle + atModerateEnd;

    const [p1, p2] = primers(seq);

    expect(p1).toBeDefined();
    expect(p2).toBeDefined();

    const tmDiff = Math.abs(p1.tm - p2.tm);

    // The algorithm should find Tm-matched pairs even if neither is at exactly 62°C
    // Moderate GC difference should allow for reasonable matching
    // With joint Tm optimization, should achieve tier 3 or better (<12°C diff)
    expect(tmDiff).toBeLessThan(12);
  });
});
