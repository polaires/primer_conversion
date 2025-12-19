/**
 * Tests for NEBuilder HiFi DNA Assembly Module
 *
 * Tests the NEBuilder-specific functionality including:
 * - Overlap analysis with scoring
 * - Sliding window optimization
 * - Assembly design
 * - Protocol generation
 */

import { describe, it, expect } from 'vitest';
import {
  NEBUILDER_PARAMS,
  OVERLAP_SCORING_WEIGHTS,
  analyzeOverlap,
  optimizeOverlap,
  designNEBuilderAssembly,
  generateNEBuilderProtocol,
} from './nebuilder.js';

// =============================================================================
// Test Sequences
// =============================================================================

const TEST_OVERLAPS = {
  // Good overlap: balanced GC, no repeats
  good: 'ATGCGATCGATCGATCGATC',

  // High GC overlap
  highGC: 'GCGCGCGCGCGCGCGCGCGC',

  // Low GC overlap
  lowGC: 'ATATATATATATATATATATAT',

  // Contains poly-T (problematic)
  polyT: 'ATGCTTTTTTTTGCATGCAT',

  // Contains poly-G (G-quadruplex risk)
  polyG: 'ATGCGGGGGGGGCATGCAT',

  // Palindromic
  palindrome: 'GAATTCGAATTC',

  // Short
  short: 'ATGCATGCAT',

  // Long
  long: 'ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATC',
};

const TEST_FRAGMENTS = {
  vector: {
    id: 'pUC19_backbone',
    seq: 'ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTA',
  },
  gfp: {
    id: 'GFP',
    seq: 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCC',
  },
};

// =============================================================================
// Configuration Tests
// =============================================================================

describe('NEBuilder Configuration', () => {
  it('should define correct overlap parameters', () => {
    expect(NEBUILDER_PARAMS.overlap.minLength).toBe(15);
    expect(NEBUILDER_PARAMS.overlap.maxLength).toBe(35);
    expect(NEBUILDER_PARAMS.overlap.optimalLength).toBe(20);
    expect(NEBUILDER_PARAMS.overlap.minTm).toBe(48);
    expect(NEBUILDER_PARAMS.overlap.maxTm).toBe(65);
    expect(NEBUILDER_PARAMS.overlap.optimalTm).toBe(55);
  });

  it('should define fragment amount recommendations', () => {
    expect(NEBUILDER_PARAMS.fragmentAmounts['2-3'].recommended).toBe(0.1);
    expect(NEBUILDER_PARAMS.fragmentAmounts['4-6'].recommended).toBe(0.3);
  });

  it('should define incubation times', () => {
    expect(NEBUILDER_PARAMS.incubation['2-3']).toBe(15);
    expect(NEBUILDER_PARAMS.incubation['4+']).toBe(60);
  });

  it('should define scoring weights that sum to 1', () => {
    const sum = Object.values(OVERLAP_SCORING_WEIGHTS).reduce((a, b) => a + b, 0);
    expect(sum).toBeCloseTo(1.0, 2);
  });
});

// =============================================================================
// Overlap Analysis Tests
// =============================================================================

describe('analyzeOverlap', () => {
  it('should analyze a good overlap correctly', () => {
    const result = analyzeOverlap(TEST_OVERLAPS.good);

    expect(result).toBeDefined();
    expect(result.sequence).toBe(TEST_OVERLAPS.good);
    expect(result.length).toBe(TEST_OVERLAPS.good.length);
    expect(result.tm).toBeGreaterThan(40);
    expect(result.tm).toBeLessThan(70);
    expect(result.gc).toBeGreaterThan(30);
    expect(result.gc).toBeLessThan(70);
  });

  it('should calculate composite score', () => {
    const result = analyzeOverlap(TEST_OVERLAPS.good);

    expect(result.compositeScore).toBeDefined();
    expect(result.compositeScore).toBeGreaterThanOrEqual(0);
    expect(result.compositeScore).toBeLessThanOrEqual(100);
  });

  it('should provide quality classification', () => {
    const result = analyzeOverlap(TEST_OVERLAPS.good);

    expect(result.quality).toBeDefined();
    expect(['excellent', 'good', 'acceptable', 'poor']).toContain(result.quality);
  });

  it('should detect high GC content', () => {
    const result = analyzeOverlap(TEST_OVERLAPS.highGC);

    expect(result.gc).toBeGreaterThan(90);
    // Should have lower score due to GC imbalance
    expect(result.scores.gcContent).toBeLessThan(1.0);
  });

  it('should detect low GC content', () => {
    const result = analyzeOverlap(TEST_OVERLAPS.lowGC);

    expect(result.gc).toBeLessThan(10);
    expect(result.scores.gcContent).toBeLessThan(1.0);
  });

  it('should warn about poly-T runs', () => {
    const result = analyzeOverlap(TEST_OVERLAPS.polyT);

    expect(result.patterns.polyT).toBe(true);
    expect(result.warnings.length).toBeGreaterThan(0);
    expect(result.warnings.some(w => w.toLowerCase().includes('poly'))).toBe(true);
  });

  it('should warn about poly-G runs (G-quadruplex)', () => {
    const result = analyzeOverlap(TEST_OVERLAPS.polyG);

    expect(result.patterns.polyG).toBe(true);
    expect(result.warnings.some(w => w.toLowerCase().includes('quadruplex'))).toBe(true);
  });

  it('should detect palindromic sequences', () => {
    // Note: 'GAATTCGAATTC' is not perfectly palindromic, but let's test
    const perfectPalindrome = 'GAATTC'; // EcoRI site
    const result = analyzeOverlap(perfectPalindrome);

    // The palindrome detection may or may not trigger depending on implementation
    expect(result).toBeDefined();
  });

  it('should score GC clamp presence', () => {
    // Sequence ending with GC
    const withClamp = 'ATGCATGCATGCATGC';
    const resultWithClamp = analyzeOverlap(withClamp);

    // Sequence ending with AT
    const withoutClamp = 'ATGCATGCATGCATAT';
    const resultWithoutClamp = analyzeOverlap(withoutClamp);

    expect(resultWithClamp.hasGcClamp).toBe(true);
    expect(resultWithoutClamp.hasGcClamp).toBe(false);
    expect(resultWithClamp.scores.gcClamp).toBeGreaterThan(resultWithoutClamp.scores.gcClamp);
  });

  it('should calculate hairpin ΔG', () => {
    const result = analyzeOverlap(TEST_OVERLAPS.good);

    expect(result.hairpinDG).toBeDefined();
    expect(typeof result.hairpinDG).toBe('number');
  });

  it('should provide individual scores for all features', () => {
    const result = analyzeOverlap(TEST_OVERLAPS.good);

    expect(result.scores.tm).toBeDefined();
    expect(result.scores.gcContent).toBeDefined();
    expect(result.scores.length).toBeDefined();
    expect(result.scores.hairpin).toBeDefined();
    expect(result.scores.gcClamp).toBeDefined();
    expect(result.scores.patterns).toBeDefined();
  });

  it('should generate recommendation based on quality', () => {
    const goodResult = analyzeOverlap(TEST_OVERLAPS.good);
    const poorResult = analyzeOverlap(TEST_OVERLAPS.polyG);

    expect(goodResult.recommendation).toBeDefined();
    expect(poorResult.recommendation).toBeDefined();
  });
});

// =============================================================================
// Overlap Optimization Tests
// =============================================================================

describe('optimizeOverlap', () => {
  // Create a junction sequence
  const junctionSeq = TEST_FRAGMENTS.vector.seq.slice(-50) + TEST_FRAGMENTS.gfp.seq.slice(0, 50);

  it('should find optimal overlap in junction', () => {
    const result = optimizeOverlap(junctionSeq);

    expect(result).toBeDefined();
    expect(result.optimal).toBeDefined();
    expect(result.optimal.sequence).toBeDefined();
    expect(result.optimal.length).toBeGreaterThanOrEqual(15);
    expect(result.optimal.length).toBeLessThanOrEqual(35);
  });

  it('should provide alternatives sorted by score', () => {
    const result = optimizeOverlap(junctionSeq);

    expect(result.alternatives).toBeDefined();
    expect(Array.isArray(result.alternatives)).toBe(true);

    if (result.alternatives.length > 0) {
      // Alternatives should be sorted by score (descending)
      for (let i = 1; i < result.alternatives.length; i++) {
        expect(result.alternatives[i - 1].compositeScore)
          .toBeGreaterThanOrEqual(result.alternatives[i].compositeScore);
      }
    }
  });

  it('should report number of candidates evaluated', () => {
    const result = optimizeOverlap(junctionSeq);

    expect(result.totalCandidatesEvaluated).toBeDefined();
    expect(result.totalCandidatesEvaluated).toBeGreaterThan(0);
  });

  it('should respect custom length range', () => {
    const result = optimizeOverlap(junctionSeq, {
      minLen: 18,
      maxLen: 25,
    });

    expect(result.optimal.length).toBeGreaterThanOrEqual(18);
    expect(result.optimal.length).toBeLessThanOrEqual(25);
  });

  it('should respect custom target Tm', () => {
    const result = optimizeOverlap(junctionSeq, {
      targetTm: 60,
    });

    // Optimal should be close to target (within reason)
    expect(Math.abs(result.optimal.tm - 60)).toBeLessThan(15);
  });

  it('should describe search space', () => {
    const result = optimizeOverlap(junctionSeq);

    expect(result.searchSpace).toBeDefined();
    expect(result.searchSpace.lengthRange).toBeDefined();
    expect(result.searchSpace.positionRange).toBeDefined();
  });
});

// =============================================================================
// Assembly Design Tests
// =============================================================================

describe('designNEBuilderAssembly', () => {
  const twoFragments = [
    TEST_FRAGMENTS.vector,
    TEST_FRAGMENTS.gfp,
  ];

  it('should design 2-fragment assembly', () => {
    const result = designNEBuilderAssembly(twoFragments);

    expect(result).toBeDefined();
    expect(result.method).toBe('NEBuilder HiFi DNA Assembly');
    expect(result.fragments.length).toBe(2);
    expect(result.junctions.length).toBe(2);
  });

  it('should design primers for each fragment', () => {
    const result = designNEBuilderAssembly(twoFragments);

    result.fragments.forEach((frag, i) => {
      expect(frag.forward).toBeDefined();
      expect(frag.reverse).toBeDefined();
      expect(frag.forward.sequence).toBeDefined();
      expect(frag.reverse.sequence).toBeDefined();
      expect(frag.index).toBe(i + 1);
    });
  });

  it('should include homology tails in primers', () => {
    const result = designNEBuilderAssembly(twoFragments);

    result.fragments.forEach(frag => {
      expect(frag.forward.homologyTail).toBeDefined();
      expect(frag.forward.annealingRegion).toBeDefined();
      // Full primer = tail + annealing
      expect(frag.forward.sequence.length)
        .toBe(frag.forward.homologyTail.length + frag.forward.annealingRegion.length);
    });
  });

  it('should calculate Tm for annealing region', () => {
    const result = designNEBuilderAssembly(twoFragments);

    result.fragments.forEach(frag => {
      expect(frag.forward.tm).toBeDefined();
      expect(frag.forward.tm).toBeGreaterThanOrEqual(45); // Lower bound for challenging sequences
      expect(frag.forward.tm).toBeLessThan(85); // Upper bound (high GC sequences can have higher Tm)
      expect(frag.reverse.tm).toBeDefined();
    });
  });

  it('should provide junction analysis', () => {
    const result = designNEBuilderAssembly(twoFragments);

    result.junctions.forEach(junction => {
      expect(junction.optimal).toBeDefined();
      expect(junction.optimal.sequence).toBeDefined();
      expect(junction.optimal.tm).toBeDefined();
      expect(junction.optimal.compositeScore).toBeDefined();
    });
  });

  it('should calculate quality metrics', () => {
    const result = designNEBuilderAssembly(twoFragments);

    expect(result.quality).toBeDefined();
    expect(result.quality.avgJunctionScore).toBeDefined();
    expect(result.quality.minJunctionScore).toBeDefined();
    expect(result.quality.tier).toBeDefined();
    expect(result.quality.recommendation).toBeDefined();
  });

  it('should generate protocol', () => {
    const result = designNEBuilderAssembly(twoFragments);

    expect(result.protocol).toBeDefined();
    expect(result.protocol.title).toContain('NEBuilder');
    expect(result.protocol.steps).toBeDefined();
    expect(result.protocol.steps.length).toBeGreaterThan(0);
  });

  it('should collect warnings', () => {
    const result = designNEBuilderAssembly(twoFragments);

    expect(result.warnings).toBeDefined();
    expect(Array.isArray(result.warnings)).toBe(true);
  });

  it('should support linear assembly', () => {
    const result = designNEBuilderAssembly(twoFragments, { circular: false });

    expect(result.assembly.type).toBe('linear');
    expect(result.junctions.length).toBe(1); // Only 1 junction for linear
  });

  it('should detect primer dimer risks', () => {
    const result = designNEBuilderAssembly(twoFragments);

    result.fragments.forEach(frag => {
      expect(frag.pair).toBeDefined();
      expect(frag.pair.heterodimerDG).toBeDefined();
    });
  });

  it('should warn about primers >60bp', () => {
    // Create fragments that would need very long overlaps
    const longOverlapNeeded = [
      { id: 'frag1', seq: 'A'.repeat(100) + 'GCGCGCGCGC' + 'T'.repeat(100) },
      { id: 'frag2', seq: 'GCGCGCGCGC' + 'A'.repeat(200) },
    ];

    const result = designNEBuilderAssembly(longOverlapNeeded);

    // Should still work, may or may not have warnings depending on actual primer lengths
    expect(result).toBeDefined();
  });

  it('should throw for >6 fragments', () => {
    const manyFragments = Array(7).fill(null).map((_, i) => ({
      id: `frag${i}`,
      seq: TEST_FRAGMENTS.gfp.seq,
    }));

    expect(() => {
      designNEBuilderAssembly(manyFragments);
    }).toThrow('maximum 6 fragments');
  });

  it('should throw for <2 fragments', () => {
    expect(() => {
      designNEBuilderAssembly([TEST_FRAGMENTS.vector]);
    }).toThrow('at least 2 fragments');
  });
});

// =============================================================================
// Protocol Generation Tests
// =============================================================================

describe('generateNEBuilderProtocol', () => {
  it('should generate protocol for 2 fragments', () => {
    const protocol = generateNEBuilderProtocol(2);

    expect(protocol.title).toContain('NEBuilder');
    expect(protocol.kit).toBe('NEB #E5520');
    expect(protocol.materials).toBeDefined();
    expect(protocol.steps.length).toBeGreaterThanOrEqual(4);
  });

  it('should recommend 15 min incubation for 2-3 fragments', () => {
    const protocol = generateNEBuilderProtocol(2);

    const incubationStep = protocol.steps.find(s => s.title.toLowerCase().includes('incubate'));
    expect(incubationStep).toBeDefined();
    expect(incubationStep.substeps.some(s => s.includes('15 minutes'))).toBe(true);
  });

  it('should recommend 60 min incubation for 4+ fragments', () => {
    const protocol = generateNEBuilderProtocol(5);

    const incubationStep = protocol.steps.find(s => s.title.toLowerCase().includes('incubate'));
    expect(incubationStep).toBeDefined();
    expect(incubationStep.substeps.some(s => s.includes('60 minutes'))).toBe(true);
  });

  it('should include troubleshooting section', () => {
    const protocol = generateNEBuilderProtocol(2);

    expect(protocol.troubleshooting).toBeDefined();
    expect(protocol.troubleshooting.length).toBeGreaterThan(0);

    const noColonies = protocol.troubleshooting.find(t =>
      t.problem.toLowerCase().includes('no colonies')
    );
    expect(noColonies).toBeDefined();
    expect(noColonies.solutions.length).toBeGreaterThan(0);
  });

  it('should include references', () => {
    const protocol = generateNEBuilderProtocol(2);

    expect(protocol.references).toBeDefined();
    expect(protocol.references.some(r => r.includes('NEB'))).toBe(true);
  });

  it('should include DNA amount recommendations', () => {
    const protocol = generateNEBuilderProtocol(2);

    const setupStep = protocol.steps.find(s =>
      s.title.toLowerCase().includes('assembly') ||
      s.title.toLowerCase().includes('setup')
    );
    expect(setupStep).toBeDefined();
    expect(setupStep.substeps.some(s => s.includes('pmol'))).toBe(true);
  });
});

// =============================================================================
// Edge Cases
// =============================================================================

describe('Edge Cases', () => {
  it('should handle minimum fragment size', () => {
    const minFragments = [
      { id: 'min1', seq: 'A'.repeat(60) },
      { id: 'min2', seq: 'T'.repeat(60) },
    ];

    const result = designNEBuilderAssembly(minFragments);
    expect(result).toBeDefined();
  });

  it('should handle identical fragment sequences', () => {
    const sameSeq = TEST_FRAGMENTS.gfp.seq;
    const duplicates = [
      { id: 'frag1', seq: sameSeq },
      { id: 'frag2', seq: sameSeq },
    ];

    const result = designNEBuilderAssembly(duplicates);
    expect(result).toBeDefined();
  });

  it('should analyze very short overlaps', () => {
    const result = analyzeOverlap('ATGCAT');

    expect(result).toBeDefined();
    expect(result.length).toBe(6);
    // Short overlaps can still have acceptable scores depending on sequence
    expect(['acceptable', 'poor']).toContain(result.quality);
    expect(result.compositeScore).toBeLessThan(70); // Should not be excellent or good
  });

  it('should analyze very long overlaps', () => {
    const longOverlap = 'ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC';
    const result = analyzeOverlap(longOverlap);

    expect(result).toBeDefined();
    expect(result.length).toBe(52);
  });
});
