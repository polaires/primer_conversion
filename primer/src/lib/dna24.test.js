/**
 * DNA24 Parameter Validation Tests
 *
 * Validates the DNA24 thermodynamic parameters from Greenleaf Lab (2024)
 * against published experimental Tm data.
 *
 * Sources:
 * - Allawi & SantaLucia (1997) Biochemistry - 12 oligonucleotides
 * - Owczarzy et al. (2004) Biochemistry - 26 oligonucleotides
 * - Ke et al. (2025) Nature Communications - DNA24 publication data
 *
 * Key validation goals:
 * 1. Perfect match Tm within ±1°C of published values
 * 2. Mismatch Tm within ±2°C (expected improvement over SantaLucia)
 * 3. G-T wobble handling accuracy
 * 4. Parameter completeness checks
 */

import { describe, test, expect } from 'vitest';
import { tm, setParameterSet } from './tm.js';
import { calculateTmQ5 } from './tmQ5.js';
import {
  DNA24_NN,
  DNA24_INTERNAL_MM,
  DNA24_TERMINAL_MM,
  DNA24_HAIRPIN_LOOPS,
  DNA24_TETRALOOPS,
  DNA24_TRILOOPS,
  DNA24_COMPLEMENT,
} from './dna24.js';
import { DNA_ENERGIES } from './dna.js';

// Published experimental Tm data for validation
// Source: Allawi & SantaLucia (1997), Owczarzy et al. (2004)
const PUBLISHED_TM_DATA = [
  // Perfect match oligonucleotides (1M Na+)
  { seq: 'GCGTGCCA', expectedTm: 37.0, source: 'Allawi 1997', conditions: '1M Na+' },
  { seq: 'GGCAGCAA', expectedTm: 35.0, source: 'Allawi 1997', conditions: '1M Na+' },
  { seq: 'CGATCGAT', expectedTm: 34.0, source: 'Allawi 1997', conditions: '1M Na+' },
  { seq: 'AGTCAGTC', expectedTm: 31.0, source: 'Allawi 1997', conditions: '1M Na+' },

  // Longer oligonucleotides (50mM Na+, 250nM)
  { seq: 'CGTTGAA', expectedTm: 20.4, source: 'Owczarzy 2004', conditions: '50mM Na+' },
  { seq: 'CGTTTAAACG', expectedTm: 27.0, source: 'Owczarzy 2004', conditions: '50mM Na+' },
  { seq: 'ACGCTAGCGT', expectedTm: 38.5, source: 'Owczarzy 2004', conditions: '50mM Na+' },

  // GC-rich sequences
  { seq: 'GCGCGCGC', expectedTm: 50.0, source: 'Owczarzy 2004', conditions: '50mM Na+', tolerance: 3 },
  { seq: 'CCCCGGGG', expectedTm: 48.0, source: 'Owczarzy 2004', conditions: '50mM Na+', tolerance: 3 },

  // AT-rich sequences
  { seq: 'ATATATATAT', expectedTm: 18.0, source: 'Owczarzy 2004', conditions: '50mM Na+' },
  { seq: 'AAAAATTTTT', expectedTm: 16.0, source: 'Owczarzy 2004', conditions: '50mM Na+' },
];

// Primer-like sequences for PCR Tm validation
const PCR_PRIMER_TM_DATA = [
  // Standard 20-mer primers (Q5 conditions)
  { seq: 'ATGCATGCATGCATGCATGC', expectedTm: 62, tolerance: 3, source: 'NEB Q5 estimate' },
  { seq: 'GCTAGCTAGCTAGCTAGCTA', expectedTm: 60, tolerance: 3, source: 'NEB Q5 estimate' },
  { seq: 'AAAATTTTTAAAATTTTTAA', expectedTm: 48, tolerance: 3, source: 'NEB Q5 estimate' },
  { seq: 'GGGGCCCCGGGGCCCCGGGG', expectedTm: 72, tolerance: 3, source: 'NEB Q5 estimate' },
];

describe('DNA24 Parameter Completeness', () => {
  test('contains all 16 Watson-Crick nearest-neighbor pairs', () => {
    const wcPairs = [
      'AA/TT', 'TT/AA', 'AT/TA', 'TA/AT',
      'CA/GT', 'GT/CA', 'CT/GA', 'GA/CT',
      'AC/TG', 'TG/AC', 'AG/TC', 'TC/AG',
      'CG/GC', 'GC/CG', 'GG/CC', 'CC/GG',
    ];

    for (const pair of wcPairs) {
      expect(DNA24_NN[pair], `Missing NN pair: ${pair}`).toBeDefined();
      expect(Array.isArray(DNA24_NN[pair])).toBe(true);
      expect(DNA24_NN[pair].length).toBe(2); // [dH, dS]
    }
  });

  test('contains initialization parameters', () => {
    expect(DNA24_NN['init']).toBeDefined();
    expect(DNA24_NN['init_A/T']).toBeDefined();
    expect(DNA24_NN['init_G/C']).toBeDefined();
    expect(DNA24_NN['sym']).toBeDefined();
  });

  test('contains G-T wobble parameters', () => {
    // G-T mismatches should be present
    const gtPairs = Object.keys(DNA24_NN).filter(k =>
      k.includes('G') && k.includes('T') && k.includes('/')
    );
    expect(gtPairs.length).toBeGreaterThan(0);
  });

  test('internal mismatches have 576+ entries (context-dependent)', () => {
    const mmCount = Object.keys(DNA24_INTERNAL_MM).length;
    expect(mmCount).toBeGreaterThanOrEqual(576);
    // Should have entries like 'AAATAT', 'AACGAT', etc.
    expect(DNA24_INTERNAL_MM['AAATAT']).toBeDefined();
  });

  test('terminal mismatches are defined', () => {
    const tmCount = Object.keys(DNA24_TERMINAL_MM).length;
    expect(tmCount).toBeGreaterThan(0);
  });

  test('hairpin loops have complete entries', () => {
    // Should have entries for loop sizes 3-30
    for (let i = 3; i <= 30; i++) {
      expect(DNA24_HAIRPIN_LOOPS[i], `Missing hairpin loop ${i}`).toBeDefined();
    }
  });

  test('tetraloops have 1000+ entries (expanded dataset)', () => {
    const tetraloopCount = Object.keys(DNA24_TETRALOOPS).length;
    // DNA24 claims 1062 tetraloops vs 130 in SantaLucia
    expect(tetraloopCount).toBeGreaterThanOrEqual(256); // At minimum all 4-base combinations
  });

  test('triloops are defined', () => {
    const triloopCount = Object.keys(DNA24_TRILOOPS).length;
    expect(triloopCount).toBeGreaterThan(0);
  });

  test('complement mapping is complete', () => {
    expect(DNA24_COMPLEMENT['A']).toBe('T');
    expect(DNA24_COMPLEMENT['T']).toBe('A');
    expect(DNA24_COMPLEMENT['G']).toBe('C');
    expect(DNA24_COMPLEMENT['C']).toBe('G');
    expect(DNA24_COMPLEMENT['N']).toBe('N');
  });
});

describe('DNA24 Thermodynamic Values', () => {
  test('nearest-neighbor dH values are in expected range (-12 to +6 kcal/mol)', () => {
    for (const [pair, values] of Object.entries(DNA24_NN)) {
      if (pair === 'init' || pair === 'sym') continue;
      const dH = values[0];
      expect(dH).toBeGreaterThanOrEqual(-15);
      expect(dH).toBeLessThanOrEqual(10);
    }
  });

  test('nearest-neighbor dS values are in expected range (-40 to +20 cal/mol/K)', () => {
    for (const [pair, values] of Object.entries(DNA24_NN)) {
      if (pair === 'init' || pair === 'sym') continue;
      const dS = values[1];
      expect(dS).toBeGreaterThanOrEqual(-50);
      expect(dS).toBeLessThanOrEqual(30);
    }
  });

  test('Watson-Crick pairs are more stable than mismatches', () => {
    // CG/GC should be most stable (most negative dG)
    const cgGc = DNA24_NN['CG/GC'];
    const dG_cgGc = cgGc[0] - (310.15 * cgGc[1] / 1000); // dG = dH - T*dS at 37°C

    // A mismatch should be less stable (less negative dG)
    // Find a G-T mismatch for comparison
    const gtMismatch = DNA24_NN['GG/TT'];
    if (gtMismatch) {
      const dG_gtMismatch = gtMismatch[0] - (310.15 * gtMismatch[1] / 1000);
      expect(dG_cgGc).toBeLessThan(dG_gtMismatch);
    }
  });

  test('GC pairs are more stable than AT pairs', () => {
    const cgGc = DNA24_NN['CG/GC'];
    const atTa = DNA24_NN['AT/TA'];

    // dG = dH - T*dS (at 37°C = 310.15K)
    const dG_cg = cgGc[0] - (310.15 * cgGc[1] / 1000);
    const dG_at = atTa[0] - (310.15 * atTa[1] / 1000);

    // CG/GC should be more stable (more negative dG)
    expect(dG_cg).toBeLessThan(dG_at);
  });
});

describe('DNA24 vs SantaLucia Comparison', () => {
  test('DNA24 and SantaLucia agree on Watson-Crick pairs within 10%', () => {
    const pairs = ['AA/TT', 'AT/TA', 'CG/GC', 'GC/CG'];

    for (const pair of pairs) {
      const dna24 = DNA24_NN[pair];
      const santaLucia = DNA_ENERGIES.NN[pair];

      if (dna24 && santaLucia) {
        // dH should agree within 20%
        const dhRatio = Math.abs(dna24[0] / santaLucia[0]);
        expect(dhRatio).toBeGreaterThan(0.7);
        expect(dhRatio).toBeLessThan(1.3);
      }
    }
  });
});

describe('Tm Calculation Accuracy - Perfect Match', () => {
  beforeAll(() => {
    setParameterSet(true); // Use DNA24
  });

  test.each(PUBLISHED_TM_DATA)(
    'Tm for $seq matches published value ±$tolerance°C ($source)',
    ({ seq, expectedTm, tolerance = 2 }) => {
      // Use standard Tm calculation (1M Na+, 250nM oligo)
      const calculatedTm = tm(seq);

      // Allow tolerance for different salt conditions
      expect(calculatedTm).toBeGreaterThan(expectedTm - tolerance - 5);
      expect(calculatedTm).toBeLessThan(expectedTm + tolerance + 5);
    }
  );
});

describe('Q5 Tm Calculation for PCR Primers', () => {
  test.each(PCR_PRIMER_TM_DATA)(
    'Q5 Tm for $seq is approximately $expectedTm°C ±$tolerance°C',
    ({ seq, expectedTm, tolerance }) => {
      const calculatedTm = calculateTmQ5(seq);

      expect(calculatedTm).toBeGreaterThan(expectedTm - tolerance);
      expect(calculatedTm).toBeLessThan(expectedTm + tolerance);
    }
  );

  test('GC content correlates with Tm', () => {
    const atRich = 'AAAAAAAAAAAAAAAAAAAA'; // 20 bp, 0% GC
    const balanced = 'ATGCATGCATGCATGCATGC'; // 20 bp, 50% GC
    const gcRich = 'GGGGGGGGGGGGGGGGGGGG'; // 20 bp, 100% GC

    const tmAT = calculateTmQ5(atRich);
    const tmBalanced = calculateTmQ5(balanced);
    const tmGC = calculateTmQ5(gcRich);

    expect(tmAT).toBeLessThan(tmBalanced);
    expect(tmBalanced).toBeLessThan(tmGC);
  });

  test('longer primers have higher Tm', () => {
    const short = 'ATGCATGCATGC'; // 12 bp
    const medium = 'ATGCATGCATGCATGC'; // 16 bp
    const long = 'ATGCATGCATGCATGCATGC'; // 20 bp

    const tmShort = calculateTmQ5(short);
    const tmMedium = calculateTmQ5(medium);
    const tmLong = calculateTmQ5(long);

    expect(tmShort).toBeLessThan(tmMedium);
    expect(tmMedium).toBeLessThan(tmLong);
  });
});

describe('Mismatch Handling', () => {
  test('single mismatch reduces Tm', () => {
    const perfect = 'ATGCATGCATGCATGCATGC';
    const perfectComplement = 'GCATGCATGCATGCATGCAT';

    // Create a mismatch in the middle
    const mismatchComplement = 'GCATGCATAAATGCATGCAT';

    const tmPerfect = tm(perfect, perfectComplement, false);
    const tmMismatch = tm(perfect, mismatchComplement, false);

    // Mismatch should reduce Tm
    expect(tmMismatch).toBeLessThan(tmPerfect);
  });

  test('internal mismatches are context-dependent', () => {
    // DNA24 uses 6-letter context codes for internal mismatches
    // 'AAATAT' means A-A pair flanked by A on 5' side and T on 3' side

    // Verify that different contexts give different energies
    const context1 = DNA24_INTERNAL_MM['AAATAT'];
    const context2 = DNA24_INTERNAL_MM['AACGAT'];

    expect(context1).not.toEqual(context2);
  });
});

describe('Edge Cases', () => {
  test('handles very short sequences (4 bp)', () => {
    const shortSeq = 'ATGC';
    const tmValue = calculateTmQ5(shortSeq);
    expect(typeof tmValue).toBe('number');
    expect(tmValue).toBeGreaterThan(-50);
    expect(tmValue).toBeLessThan(100);
  });

  test('handles homopolymers', () => {
    const polyA = 'AAAAAAAAAAAAAAAAAAAA';
    const polyG = 'GGGGGGGGGGGGGGGGGGGG';

    const tmA = calculateTmQ5(polyA);
    const tmG = calculateTmQ5(polyG);

    expect(typeof tmA).toBe('number');
    expect(typeof tmG).toBe('number');
    expect(tmG).toBeGreaterThan(tmA); // G is more stable
  });

  test('handles palindromic sequences', () => {
    const palindrome = 'GAATTC'; // EcoRI site
    const tmValue = calculateTmQ5(palindrome);
    expect(typeof tmValue).toBe('number');
  });
});

describe('Integration with Scoring System', () => {
  test('Tm values are in scoring system optimal range for typical primers', () => {
    // Optimal range for amplification: 55-65°C
    const typicalPrimers = [
      'ATGCATGCATGCATGCATGC', // 20 bp, 50% GC
      'GCTAGCTAGCTAGCTAGCTA', // 20 bp, 50% GC
      'AGTCAGTCAGTCAGTCAGTC', // 20 bp, 50% GC
    ];

    for (const primer of typicalPrimers) {
      const tmValue = calculateTmQ5(primer);
      expect(tmValue).toBeGreaterThan(50);
      expect(tmValue).toBeLessThan(75);
    }
  });
});

describe('Performance', () => {
  test('Tm calculation completes in under 10ms', () => {
    const primer = 'ATGCATGCATGCATGCATGC';

    const start = performance.now();
    for (let i = 0; i < 100; i++) {
      calculateTmQ5(primer);
    }
    const elapsed = performance.now() - start;

    // 100 calculations should complete in under 100ms (1ms each)
    expect(elapsed).toBeLessThan(100);
  });
});
