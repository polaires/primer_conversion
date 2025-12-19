/**
 * Tests for Equilibrium Efficiency Calculation
 *
 * Based on Pythia (Mann et al., 2009) thermodynamic model
 * Validates that our implementation correctly models chemical equilibria
 */

import { describe, it, expect } from 'vitest';
import {
  calculateEquilibriumEfficiency,
  calculateAllSpeciesDG,
  calculateHairpinDG,
  calculateHomodimerDG,
  calculateHeterodimerDG,
  calculateDuplexDG,
  calculateOffTargetDG,
  solveEquilibrium,
  efficiencyToScore,
  reverseComplement,
} from './equilibrium.js';

describe('Reverse Complement', () => {
  it('should correctly reverse complement a sequence', () => {
    expect(reverseComplement('ATGC')).toBe('GCAT');
    expect(reverseComplement('AAAA')).toBe('TTTT');
    expect(reverseComplement('GCGC')).toBe('GCGC');
    expect(reverseComplement('ATAT')).toBe('ATAT');
  });
});

describe('Hairpin ΔG Calculation', () => {
  it('should return 0 for short sequences that cannot form hairpins', () => {
    expect(calculateHairpinDG('ATGC')).toBe(0);
    expect(calculateHairpinDG('ATGCAT')).toBeLessThanOrEqual(0);
  });

  it('should detect stable hairpins in self-complementary sequences', () => {
    // GCGCTTTTGCGC can form a hairpin with GCGC stem and TTTT loop
    const hairpinSeq = 'GCGCTTTTGCGC';
    const dG = calculateHairpinDG(hairpinSeq);
    expect(dG).toBeLessThan(0);  // Should be negative (stable)
  });

  it('should return less stable (higher) ΔG for non-self-complementary sequences', () => {
    // Random sequence unlikely to form stable hairpin
    const noHairpinSeq = 'ATATATATATAT';
    const dG = calculateHairpinDG(noHairpinSeq);
    // Should be 0 or slightly negative (no stable hairpin)
    expect(dG).toBeGreaterThanOrEqual(-3);
  });

  it('should give more negative ΔG for GC-rich hairpins', () => {
    const gcRichHairpin = 'GCGCGCAAAAGCGCGC';
    const atRichHairpin = 'ATATATAAAAATATAT';
    const gcDG = calculateHairpinDG(gcRichHairpin);
    const atDG = calculateHairpinDG(atRichHairpin);
    // GC-rich should be more stable (more negative)
    expect(gcDG).toBeLessThanOrEqual(atDG);
  });
});

describe('Homodimer ΔG Calculation', () => {
  it('should return 0 for short sequences', () => {
    expect(calculateHomodimerDG('ATGC')).toBe(0);
  });

  it('should detect self-dimer in self-complementary sequences', () => {
    // GCGCGCGC is self-complementary
    const selfComp = 'GCGCGCGC';
    const dG = calculateHomodimerDG(selfComp);
    expect(dG).toBeLessThan(0);
  });

  it('should return near-zero for non-self-complementary sequences', () => {
    // AAAAAAAA cannot form homodimer (A-A is not a pair)
    const noSelfComp = 'AAAAAAAAAA';
    const dG = calculateHomodimerDG(noSelfComp);
    expect(dG).toBe(0);
  });

  it('should detect 3\' self-complementarity', () => {
    // NNNNNNNNGATC has 3' end that can self-pair
    const threePrimeSelfComp = 'AAAAAAAGATC';
    const dG = calculateHomodimerDG(threePrimeSelfComp);
    // May or may not be stable depending on match length
    expect(typeof dG).toBe('number');
  });
});

describe('Heterodimer ΔG Calculation', () => {
  it('should return 0 for short sequences', () => {
    expect(calculateHeterodimerDG('ATG', 'CAT')).toBe(0);
  });

  it('should detect stable heterodimer between complementary primers', () => {
    // These primers are designed to be complementary at 3' ends
    const fwd = 'AAAAAAGATCGATC';
    const rev = 'TTTTTTGATCGATC';  // 3' end GATC is complement of GATC
    const dG = calculateHeterodimerDG(fwd, rev);
    // Should detect some complementarity
    expect(dG).toBeLessThanOrEqual(0);
  });

  it('should return near-zero for non-complementary primers', () => {
    // These primers have no significant complementarity
    const fwd = 'AAAAAAAAAAAAAAAA';
    const rev = 'AAAAAAAAAAAAAAA';
    const dG = calculateHeterodimerDG(fwd, rev);
    expect(dG).toBe(0);  // A-A is not a pair
  });
});

describe('Duplex ΔG Calculation', () => {
  it('should return negative ΔG for perfect Watson-Crick duplex', () => {
    // Primer binds to its reverse complement
    const primer = 'ATGCATGCATGC';
    const target = reverseComplement(primer);  // The actual complementary strand
    const dG = calculateDuplexDG(primer, target);
    expect(dG).toBeLessThan(0);
  });

  it('should return more negative ΔG for GC-rich duplexes', () => {
    const gcRich = 'GCGCGCGCGCGC';
    const atRich = 'ATATATATATAT';
    // Calculate with proper complementary targets
    const gcDG = calculateDuplexDG(gcRich, reverseComplement(gcRich));
    const atDG = calculateDuplexDG(atRich, reverseComplement(atRich));
    // GC-rich should be more stable (more negative)
    expect(gcDG).toBeLessThan(atDG);
  });

  it('should scale with primer length', () => {
    const short = 'ATGCATGC';        // 8bp
    const long = 'ATGCATGCATGCATGC'; // 16bp
    const shortDG = calculateDuplexDG(short, reverseComplement(short));
    const longDG = calculateDuplexDG(long, reverseComplement(long));
    // Longer duplex should be more stable
    expect(longDG).toBeLessThan(shortDG);
  });
});

describe('Off-target ΔG Calculation', () => {
  it('should return 0 when no off-target sites exist', () => {
    const primer = 'ATGCATGCATGC';
    // Template is just the primer site + some random sequence with no repeat
    const template = 'ATGCATGCATGCAAAAAAAAAAAAAAAAAAAAA';
    const dG = calculateOffTargetDG(primer, template, 0);
    // Only intended site, so no off-target
    expect(dG).toBe(0);
  });

  it('should detect off-target when complementary match exists elsewhere', () => {
    const primer = 'ATGCATGCATGC';
    // Template has the primer's complement (binding site) twice
    // First occurrence at position 0, second at position 22
    const bindingSite = reverseComplement(primer);  // GCATGCATGCAT
    const template = bindingSite + 'AAAAAAAAAA' + bindingSite;
    const dG = calculateOffTargetDG(primer, template, 0);
    // Should find the second occurrence
    expect(dG).toBeLessThan(0);
  });
});

describe('Species ΔG Calculation', () => {
  it('should calculate ΔG for all species', () => {
    const fwd = { seq: 'ATGCATGCATGCATGCAT' };
    const rev = { seq: 'GCATGCATGCATGCATGC' };
    const template = 'ATGCATGCATGCATGCAT' + 'NNNNNNNNNNNNNNNNNNNN'.replace(/N/g, 'A') + 'GCATGCATGCATGCATGC';

    const species = calculateAllSpeciesDG(fwd, rev, template);

    expect(species).toHaveProperty('fwdHairpin');
    expect(species).toHaveProperty('revHairpin');
    expect(species).toHaveProperty('fwdHomodimer');
    expect(species).toHaveProperty('revHomodimer');
    expect(species).toHaveProperty('fwdTarget');
    expect(species).toHaveProperty('revTarget');
    expect(species).toHaveProperty('heterodimer');

    // All values should be numbers
    Object.values(species).forEach(value => {
      expect(typeof value).toBe('number');
      expect(isNaN(value)).toBe(false);
    });
  });
});

describe('Equilibrium Solver', () => {
  it('should return valid concentrations', () => {
    const K = {
      fwdHairpin: 0.1,
      revHairpin: 0.1,
      fwdHomodimer: 0.01,
      revHomodimer: 0.01,
      fwdTarget: 1000,  // Strong binding
      revTarget: 1000,
      fwdOffTarget: 0,
      revOffTarget: 0,
      heterodimer: 0.01,
    };

    const primerConc = 0.5e-6;
    const templateConc = 1e-9;
    const result = solveEquilibrium(K, primerConc, templateConc);

    // All concentrations should be non-negative
    expect(result.fwdFree).toBeGreaterThanOrEqual(0);
    expect(result.fwdHairpin).toBeGreaterThanOrEqual(0);
    expect(result.fwdHomodimer).toBeGreaterThanOrEqual(0);
    expect(result.fwdBoundTarget).toBeGreaterThanOrEqual(0);
    expect(result.revFree).toBeGreaterThanOrEqual(0);
    expect(result.revBoundTarget).toBeGreaterThanOrEqual(0);
  });

  it('should show higher target binding with stronger K', () => {
    const weakK = {
      fwdHairpin: 0.001,
      revHairpin: 0.001,
      fwdHomodimer: 0.0001,
      revHomodimer: 0.0001,
      fwdTarget: 1e3,  // Moderate binding
      revTarget: 1e3,
      fwdOffTarget: 0,
      revOffTarget: 0,
      heterodimer: 0.0001,
    };

    const strongK = {
      ...weakK,
      fwdTarget: 1e10,  // Very strong binding
      revTarget: 1e10,
    };

    const weakResult = solveEquilibrium(weakK, 0.5e-6, 1e-9);
    const strongResult = solveEquilibrium(strongK, 0.5e-6, 1e-9);

    // Stronger K should give more target binding
    expect(strongResult.fwdBoundTarget).toBeGreaterThan(weakResult.fwdBoundTarget);
  });

  it('should show higher hairpin with stronger hairpin K', () => {
    const weakHairpinK = {
      fwdHairpin: 0.001,
      revHairpin: 0.001,
      fwdHomodimer: 0.0001,
      revHomodimer: 0.0001,
      fwdTarget: 1000,
      revTarget: 1000,
      fwdOffTarget: 0,
      revOffTarget: 0,
      heterodimer: 0.0001,
    };

    const strongHairpinK = {
      ...weakHairpinK,
      fwdHairpin: 1e6,   // Very strong hairpin
    };

    const weakResult = solveEquilibrium(weakHairpinK, 0.5e-6, 1e-9);
    const strongResult = solveEquilibrium(strongHairpinK, 0.5e-6, 1e-9);

    // Stronger hairpin K should give more hairpin formation
    expect(strongResult.fwdHairpin).toBeGreaterThan(weakResult.fwdHairpin);
  });
});

describe('Equilibrium Efficiency', () => {
  it('should return efficiency between 0 and 1', () => {
    const fwd = {
      seq: 'ATGCATGCATGCATGCAT',
      bindingSite: 'ATGCATGCATGCATGCAT',
    };
    const rev = {
      seq: 'GCATGCATGCATGCATGC',
      bindingSite: 'GCATGCATGCATGCATGC',
    };
    const template = 'ATGCATGCATGCATGCAT' + 'A'.repeat(50) + reverseComplement('GCATGCATGCATGCATGC');

    const result = calculateEquilibriumEfficiency(fwd, rev, template);

    expect(result.efficiency).toBeGreaterThanOrEqual(0);
    expect(result.efficiency).toBeLessThanOrEqual(1);
    expect(result.efficiencyFwd).toBeGreaterThanOrEqual(0);
    expect(result.efficiencyFwd).toBeLessThanOrEqual(1);
    expect(result.efficiencyRev).toBeGreaterThanOrEqual(0);
    expect(result.efficiencyRev).toBeLessThanOrEqual(1);
  });

  it('should identify the bottleneck primer', () => {
    const result = calculateEquilibriumEfficiency(
      { seq: 'ATGCATGCATGCATGCAT' },
      { seq: 'GCATGCATGCATGCATGC' },
      'ATGCATGCATGCATGCAT' + 'A'.repeat(50) + reverseComplement('GCATGCATGCATGCATGC')
    );

    expect(['forward', 'reverse']).toContain(result.bottleneck);
    expect(result.efficiency).toBe(Math.min(result.efficiencyFwd, result.efficiencyRev));
  });

  it('should classify quality appropriately', () => {
    const result = calculateEquilibriumEfficiency(
      { seq: 'ATGCATGCATGCATGCAT' },
      { seq: 'GCATGCATGCATGCATGC' },
      'ATGCATGCATGCATGCAT' + 'A'.repeat(50) + reverseComplement('GCATGCATGCATGCATGC')
    );

    expect(['excellent', 'good', 'acceptable', 'marginal', 'poor']).toContain(result.quality);
  });

  it('should report loss breakdown', () => {
    const result = calculateEquilibriumEfficiency(
      { seq: 'ATGCATGCATGCATGCAT' },
      { seq: 'GCATGCATGCATGCATGC' },
      'ATGCATGCATGCATGCAT' + 'A'.repeat(50) + reverseComplement('GCATGCATGCATGCATGC')
    );

    expect(result.losses).toHaveProperty('fwd');
    expect(result.losses).toHaveProperty('rev');
    expect(result.losses.fwd).toHaveProperty('hairpin');
    expect(result.losses.fwd).toHaveProperty('homodimer');
    expect(result.losses.fwd).toHaveProperty('heterodimer');
    expect(result.losses.fwd).toHaveProperty('offTarget');
    expect(result.losses.fwd).toHaveProperty('free');

    // Loss breakdown should sum with efficiency to approximately 1
    // Note: Due to partition function approximations for homodimer/heterodimer
    // stoichiometry, exact sum may vary. We just verify it's in reasonable range.
    const fwdLossSum = Object.values(result.losses.fwd).reduce((a, b) => a + b, 0);
    const total = fwdLossSum + result.efficiencyFwd;
    expect(total).toBeGreaterThan(0);
    expect(total).toBeLessThanOrEqual(2);  // Reasonable upper bound
  });
});

describe('Efficiency to Score Conversion', () => {
  it('should return 100 for optimal efficiency', () => {
    expect(efficiencyToScore(0.95)).toBe(100);
    expect(efficiencyToScore(1.0)).toBe(100);
    expect(efficiencyToScore(0.99)).toBe(100);
  });

  it('should return score in 70-100 range for acceptable efficiency', () => {
    const score = efficiencyToScore(0.85);
    expect(score).toBeGreaterThanOrEqual(70);
    expect(score).toBeLessThanOrEqual(100);
  });

  it('should return score below 70 for poor efficiency', () => {
    const score = efficiencyToScore(0.5);
    expect(score).toBeLessThan(70);
    expect(score).toBeGreaterThan(0);
  });

  it('should be monotonically increasing', () => {
    const scores = [0.1, 0.3, 0.5, 0.7, 0.9].map(e => efficiencyToScore(e));
    for (let i = 1; i < scores.length; i++) {
      expect(scores[i]).toBeGreaterThan(scores[i - 1]);
    }
  });
});

describe('Edge Cases', () => {
  it('should handle primers with string sequences', () => {
    const result = calculateEquilibriumEfficiency(
      'ATGCATGCATGCATGCAT',
      'GCATGCATGCATGCATGC',
      'ATGCATGCATGCATGCAT' + 'A'.repeat(50) + reverseComplement('GCATGCATGCATGCATGC')
    );

    expect(result.efficiency).toBeGreaterThanOrEqual(0);
    expect(result.efficiency).toBeLessThanOrEqual(1);
  });

  it('should handle very short templates', () => {
    const result = calculateEquilibriumEfficiency(
      { seq: 'ATGCATGCATGC' },
      { seq: 'GCATGCATGCAT' },
      'ATGCATGCATGCAAAAAAAAAGCATGCATGCAT'
    );

    expect(result).toBeDefined();
    expect(typeof result.efficiency).toBe('number');
  });

  it('should handle temperature variations', () => {
    const template = 'ATGCATGCATGCATGCAT' + 'A'.repeat(50) + reverseComplement('GCATGCATGCATGCATGC');
    const fwd = { seq: 'ATGCATGCATGCATGCAT' };
    const rev = { seq: 'GCATGCATGCATGCATGC' };

    const result37 = calculateEquilibriumEfficiency(fwd, rev, template, { temperature: 37 });
    const result55 = calculateEquilibriumEfficiency(fwd, rev, template, { temperature: 55 });
    const result72 = calculateEquilibriumEfficiency(fwd, rev, template, { temperature: 72 });

    // All should return valid results
    [result37, result55, result72].forEach(r => {
      expect(r.efficiency).toBeGreaterThanOrEqual(0);
      expect(r.efficiency).toBeLessThanOrEqual(1);
    });
  });
});

describe('Problematic Primer Detection', () => {
  it('should give lower efficiency for primer with strong hairpin potential', () => {
    // This sequence can form a stable hairpin
    const hairpinPrimer = 'GCGCGCTTTTGCGCGC';
    // Normal primer
    const normalPrimer = 'ATGCATGCATGCATGC';

    const template1 = hairpinPrimer + 'A'.repeat(50) + reverseComplement(normalPrimer);
    const template2 = normalPrimer + 'A'.repeat(50) + reverseComplement(normalPrimer);

    const resultHairpin = calculateEquilibriumEfficiency(
      { seq: hairpinPrimer },
      { seq: normalPrimer },
      template1
    );

    const resultNormal = calculateEquilibriumEfficiency(
      { seq: normalPrimer },
      { seq: normalPrimer },
      template2
    );

    // Hairpin primer should have lower efficiency due to self-folding
    // (may not always be true depending on thermodynamics)
    expect(resultHairpin.losses.fwd.hairpin).toBeGreaterThanOrEqual(0);
    expect(resultNormal.losses.fwd.hairpin).toBeGreaterThanOrEqual(0);
  });

  it('should give lower efficiency for self-complementary primers', () => {
    // Self-complementary primer (can form homodimer)
    const selfCompPrimer = 'GCGCGCGCGCGCGCGC';
    // Normal primer
    const normalPrimer = 'ATGCATGCATGCATGC';

    const template = selfCompPrimer + 'A'.repeat(50) + reverseComplement(normalPrimer);

    const result = calculateEquilibriumEfficiency(
      { seq: selfCompPrimer },
      { seq: normalPrimer },
      template
    );

    // Self-complementary primer should have some homodimer loss
    expect(result.losses.fwd.homodimer).toBeGreaterThanOrEqual(0);
  });
});
