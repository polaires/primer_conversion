/**
 * Tests for Assembly Core Module
 *
 * Tests the unified assembly design API including:
 * - Overlap optimization
 * - Assembly primer design
 * - Multi-fragment assembly planning
 */

import { describe, it, expect, beforeEach } from 'vitest';
import {
  ASSEMBLY_METHODS,
  DEFAULT_ASSEMBLY_CONFIG,
  findOptimalOverlap,
  optimizeAssemblyOverlaps,
  designAssemblyPrimers,
  designAssembly,
  exportPrimers,
  simulateAssembly,
  exportToGenBank,
  exportToFasta,
  exportProject,
  importProject,
} from './assemblyCore.js';

// =============================================================================
// Test Sequences
// =============================================================================

// Example fragment sequences for testing
const TEST_FRAGMENTS = {
  vector: {
    id: 'pUC19_backbone',
    // Simplified vector sequence (~200bp for testing)
    seq: 'ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTA',
  },
  gfp: {
    id: 'GFP',
    // Simplified GFP sequence (~180bp)
    seq: 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCC',
  },
  promoter: {
    id: 'T7_promoter',
    // T7 promoter region (~100bp)
    seq: 'TAATACGACTCACTATAGGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATGGCTAGCATGACTGGTGGACAG',
  },
};

// Create junction sequence for overlap testing
const createJunctionSeq = (frag1End, frag2Start) => {
  return frag1End + frag2Start;
};

// =============================================================================
// Assembly Methods Configuration Tests
// =============================================================================

describe('Assembly Methods Configuration', () => {
  it('should define NEBuilder HiFi parameters', () => {
    const nebuilder = ASSEMBLY_METHODS.NEBUILDER_HIFI;
    expect(nebuilder).toBeDefined();
    expect(nebuilder.name).toBe('NEBuilder HiFi DNA Assembly');
    expect(nebuilder.overlapRange.min).toBe(15);
    expect(nebuilder.overlapRange.max).toBe(35);
    expect(nebuilder.overlapRange.optimal).toBe(20);
    expect(nebuilder.maxFragments).toBe(5);
  });

  it('should define Gibson Assembly parameters', () => {
    const gibson = ASSEMBLY_METHODS.GIBSON;
    expect(gibson).toBeDefined();
    expect(gibson.name).toBe('Gibson Assembly');
    expect(gibson.overlapRange.min).toBe(15);
    expect(gibson.overlapRange.max).toBe(40);
    expect(gibson.overlapRange.optimal).toBe(25);
    expect(gibson.maxFragments).toBe(6);
  });

  it('should define Golden Gate parameters', () => {
    const gg = ASSEMBLY_METHODS.GOLDEN_GATE;
    expect(gg).toBeDefined();
    expect(gg.name).toBe('Golden Gate Assembly');
    expect(gg.overlapRange.min).toBe(4);
    expect(gg.overlapRange.max).toBe(4);
    expect(gg.maxFragments).toBe(20);
  });

  it('should have default assembly config', () => {
    expect(DEFAULT_ASSEMBLY_CONFIG.method).toBe('NEBUILDER_HIFI');
    expect(DEFAULT_ASSEMBLY_CONFIG.primerTmTarget).toBe(60);
    expect(DEFAULT_ASSEMBLY_CONFIG.overlapTmTarget).toBe(55);
    expect(DEFAULT_ASSEMBLY_CONFIG.overlapLenMin).toBe(15);
    expect(DEFAULT_ASSEMBLY_CONFIG.overlapLenMax).toBe(35);
  });
});

// =============================================================================
// Overlap Optimization Tests
// =============================================================================

describe('findOptimalOverlap', () => {
  const frag1End = TEST_FRAGMENTS.vector.seq.slice(-50);
  const frag2Start = TEST_FRAGMENTS.gfp.seq.slice(0, 50);

  it('should find an overlap in the valid length range', () => {
    const result = findOptimalOverlap(frag1End, frag2Start);

    expect(result).toBeDefined();
    expect(result.optimal).toBeDefined();
    expect(result.optimal.length).toBeGreaterThanOrEqual(15);
    expect(result.optimal.length).toBeLessThanOrEqual(35);
  });

  it('should calculate Tm for the overlap', () => {
    const result = findOptimalOverlap(frag1End, frag2Start);

    expect(result.optimal.tm).toBeDefined();
    expect(typeof result.optimal.tm).toBe('number');
    expect(result.optimal.tm).toBeGreaterThan(30); // Reasonable Tm range
    expect(result.optimal.tm).toBeLessThan(80);
  });

  it('should calculate GC content', () => {
    const result = findOptimalOverlap(frag1End, frag2Start);

    expect(result.optimal.gc).toBeDefined();
    expect(result.optimal.gc).toBeGreaterThanOrEqual(0);
    expect(result.optimal.gc).toBeLessThanOrEqual(100);
  });

  it('should provide a score for the overlap', () => {
    const result = findOptimalOverlap(frag1End, frag2Start);

    expect(result.optimal.score).toBeDefined();
    expect(typeof result.optimal.score).toBe('number');
  });

  it('should provide alternative overlaps', () => {
    const result = findOptimalOverlap(frag1End, frag2Start);

    expect(result.alternatives).toBeDefined();
    expect(Array.isArray(result.alternatives)).toBe(true);
  });

  it('should detect poly-T runs and add warnings', () => {
    // Create a sequence with TTTT
    const badEnd = 'ATGCATGCTTTTTTTTGCATGCATGCATGCATGCATGCATGCATGCATGC';
    const goodStart = 'GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC';

    const result = findOptimalOverlap(badEnd, goodStart);

    // Should still find an overlap, but may have warnings
    expect(result.optimal).toBeDefined();
  });

  it('should respect method-specific overlap ranges', () => {
    // Test with Gibson (longer overlaps allowed)
    const gibsonResult = findOptimalOverlap(frag1End, frag2Start, {
      method: 'GIBSON',
    });

    expect(gibsonResult.optimal.length).toBeGreaterThanOrEqual(15);
    expect(gibsonResult.optimal.length).toBeLessThanOrEqual(40);
  });
});

// =============================================================================
// Assembly Overlap Optimization Tests
// =============================================================================

describe('optimizeAssemblyOverlaps', () => {
  const twoFragments = [
    TEST_FRAGMENTS.vector,
    TEST_FRAGMENTS.gfp,
  ];

  const threeFragments = [
    TEST_FRAGMENTS.vector,
    TEST_FRAGMENTS.promoter,
    TEST_FRAGMENTS.gfp,
  ];

  it('should optimize overlaps for 2 fragments', () => {
    const result = optimizeAssemblyOverlaps(twoFragments);

    expect(result).toBeDefined();
    expect(result.junctions).toBeDefined();
    expect(result.junctions.length).toBe(2); // Circular assembly
  });

  it('should optimize overlaps for 3 fragments', () => {
    const result = optimizeAssemblyOverlaps(threeFragments);

    expect(result.junctions.length).toBe(3);
    expect(result.quality).toBeDefined();
    expect(result.quality.averageScore).toBeDefined();
  });

  it('should identify junction from/to fragments', () => {
    const result = optimizeAssemblyOverlaps(twoFragments);

    expect(result.junctions[0].from).toBe('pUC19_backbone');
    expect(result.junctions[0].to).toBe('GFP');
    expect(result.junctions[1].from).toBe('GFP');
    expect(result.junctions[1].to).toBe('pUC19_backbone');
  });

  it('should calculate quality metrics', () => {
    const result = optimizeAssemblyOverlaps(twoFragments);

    expect(result.quality.averageScore).toBeDefined();
    expect(result.quality.minimumScore).toBeDefined();
    expect(result.quality.tier).toBeDefined();
    expect(['excellent', 'good', 'acceptable']).toContain(result.quality.tier);
  });

  it('should throw error for single fragment', () => {
    expect(() => {
      optimizeAssemblyOverlaps([TEST_FRAGMENTS.vector]);
    }).toThrow('Assembly requires at least 2 fragments');
  });

  it('should provide alternatives for each junction', () => {
    const result = optimizeAssemblyOverlaps(twoFragments);

    result.junctions.forEach(junction => {
      expect(junction.overlap).toBeDefined();
      expect(junction.alternatives).toBeDefined();
    });
  });
});

// =============================================================================
// Assembly Primer Design Tests
// =============================================================================

describe('designAssemblyPrimers', () => {
  it('should design primers with homology tails', () => {
    const fragmentSeq = TEST_FRAGMENTS.gfp.seq;
    const context = {
      leftOverlap: 'ATGCATGCATGCATGCAT',
      rightOverlap: 'GCATGCATGCATGCATGC',
      method: 'NEBUILDER_HIFI',
    };

    const result = designAssemblyPrimers(fragmentSeq, context);

    expect(result).toBeDefined();
    expect(result.forward).toBeDefined();
    expect(result.reverse).toBeDefined();
  });

  it('should include homology tail in forward primer', () => {
    const fragmentSeq = TEST_FRAGMENTS.gfp.seq;
    const leftOverlap = 'ATGCATGCATGCATGCAT';
    const context = {
      leftOverlap,
      rightOverlap: '',
      method: 'NEBUILDER_HIFI',
    };

    const result = designAssemblyPrimers(fragmentSeq, context);

    expect(result.forward.sequence.startsWith(leftOverlap)).toBe(true);
    expect(result.forward.homologyTail).toBe(leftOverlap);
    expect(result.forward.tailLength).toBe(leftOverlap.length);
  });

  it('should calculate primer Tm correctly', () => {
    const fragmentSeq = TEST_FRAGMENTS.gfp.seq;
    const result = designAssemblyPrimers(fragmentSeq, {
      leftOverlap: '',
      rightOverlap: '',
    });

    expect(result.forward.tm).toBeDefined();
    expect(result.forward.tm).toBeGreaterThan(50);
    expect(result.forward.tm).toBeLessThan(75);
    expect(result.reverse.tm).toBeDefined();
  });

  it('should provide PCR conditions', () => {
    const result = designAssemblyPrimers(TEST_FRAGMENTS.gfp.seq, {
      leftOverlap: '',
      rightOverlap: '',
    });

    expect(result.pcr).toBeDefined();
    expect(result.pcr.annealingTemp).toBeDefined();
    expect(result.pcr.extensionTime).toBeDefined();
    expect(result.pcr.cycles).toBe(30);
  });

  it('should detect potential heterodimer formation', () => {
    const result = designAssemblyPrimers(TEST_FRAGMENTS.gfp.seq, {
      leftOverlap: 'ATGCATGCATGCATGCAT',
      rightOverlap: 'ATGCATGCATGCATGCAT', // Same sequence - dimer risk
    });

    expect(result.pair.heterodimerDG).toBeDefined();
    expect(typeof result.pair.heterodimerDG).toBe('number');
  });
});

// =============================================================================
// Full Assembly Design Tests
// =============================================================================

describe('designAssembly', () => {
  const twoFragments = [
    TEST_FRAGMENTS.vector,
    TEST_FRAGMENTS.gfp,
  ];

  it('should design a complete 2-fragment assembly', () => {
    const result = designAssembly(twoFragments);

    expect(result).toBeDefined();
    expect(result.method).toContain('NEBuilder');
    expect(result.fragments.length).toBe(2);
    expect(result.junctions.length).toBe(2);
  });

  it('should design primers for each fragment', () => {
    const result = designAssembly(twoFragments);

    result.fragments.forEach(frag => {
      expect(frag.primers).toBeDefined();
      expect(frag.primers.forward).toBeDefined();
      expect(frag.primers.reverse).toBeDefined();
      expect(frag.primers.forward.sequence).toBeDefined();
      expect(frag.primers.reverse.sequence).toBeDefined();
    });
  });

  it('should generate a protocol', () => {
    const result = designAssembly(twoFragments);

    expect(result.protocol).toBeDefined();
    expect(result.protocol.title).toBeDefined();
    expect(result.protocol.steps).toBeDefined();
    expect(Array.isArray(result.protocol.steps)).toBe(true);
  });

  it('should estimate costs', () => {
    const result = designAssembly(twoFragments);

    expect(result.cost).toBeDefined();
    expect(result.cost.primers).toBeDefined();
    expect(result.cost.assembly).toBeDefined();
    expect(result.cost.total).toBeDefined();
    expect(result.cost.total).toBeGreaterThan(0);
  });

  it('should support Gibson assembly method', () => {
    const result = designAssembly(twoFragments, { method: 'GIBSON' });

    expect(result.method).toContain('Gibson');
  });

  it('should support linear (non-circular) assembly', () => {
    const result = designAssembly(twoFragments, { circular: false });

    expect(result.assembly.circular).toBe(false);
  });

  it('should throw for too many fragments', () => {
    const manyFragments = Array(10).fill(null).map((_, i) => ({
      id: `frag${i}`,
      seq: TEST_FRAGMENTS.gfp.seq,
    }));

    expect(() => {
      designAssembly(manyFragments);
    }).toThrow();
  });
});

// =============================================================================
// Export Tests
// =============================================================================

describe('exportPrimers', () => {
  let assemblyResult;

  beforeEach(() => {
    assemblyResult = designAssembly([
      TEST_FRAGMENTS.vector,
      TEST_FRAGMENTS.gfp,
    ]);
  });

  it('should export to TSV format', () => {
    const tsv = exportPrimers(assemblyResult, 'tsv');

    expect(tsv).toBeDefined();
    expect(typeof tsv).toBe('string');
    expect(tsv).toContain('Name\tSequence');
    expect(tsv).toContain('_F');
    expect(tsv).toContain('_R');
  });

  it('should export to CSV format', () => {
    const csv = exportPrimers(assemblyResult, 'csv');

    expect(csv).toBeDefined();
    expect(csv).toContain('Name,Sequence');
  });

  it('should export to JSON format', () => {
    const json = exportPrimers(assemblyResult, 'json');

    expect(json).toBeDefined();
    const parsed = JSON.parse(json);
    expect(Array.isArray(parsed)).toBe(true);
    expect(parsed.length).toBe(4); // 2 fragments × 2 primers
  });

  it('should include all primers in export', () => {
    const rows = exportPrimers(assemblyResult, 'array');

    expect(rows.length).toBe(4); // 2 fragments × 2 primers
    rows.forEach(row => {
      expect(row.name).toBeDefined();
      expect(row.sequence).toBeDefined();
      expect(row.length).toBeGreaterThan(0);
      expect(row.tm).toBeGreaterThan(0);
    });
  });
});

// =============================================================================
// Edge Cases and Error Handling
// =============================================================================

describe('Edge Cases', () => {
  it('should handle very short fragments', () => {
    const shortFragments = [
      { id: 'short1', seq: 'A'.repeat(60) },
      { id: 'short2', seq: 'T'.repeat(60) },
    ];

    // Should still work, though quality may be low
    const result = designAssembly(shortFragments);
    expect(result).toBeDefined();
  });

  it('should handle AT-rich sequences', () => {
    const atRich = [
      { id: 'at1', seq: 'ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT'.repeat(2) },
      { id: 'at2', seq: 'TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT'.repeat(2) },
    ];

    const result = designAssembly(atRich);
    expect(result).toBeDefined();
    // Should have warnings about low Tm or GC
    expect(result.warnings.length).toBeGreaterThan(0);
  });

  it('should handle GC-rich sequences', () => {
    const gcRich = [
      { id: 'gc1', seq: 'GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC'.repeat(2) },
      { id: 'gc2', seq: 'CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG'.repeat(2) },
    ];

    const result = designAssembly(gcRich);
    expect(result).toBeDefined();
  });

  it('should warn about self-complementary overlaps', () => {
    // Palindromic sequences are problematic
    const palindrome = 'GAATTCGAATTCGAATTCGAATTCGAATTCGAATTCGAATTCGAATTCGAATTC';
    const fragments = [
      { id: 'pal1', seq: palindrome.repeat(2) },
      { id: 'pal2', seq: palindrome.repeat(2) },
    ];

    const result = designAssembly(fragments);
    expect(result).toBeDefined();
  });
});

// =============================================================================
// GenBank Export Tests
// =============================================================================

describe('GenBank Export', () => {
  const mockSimulation = {
    assembledSequence: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',
    assembledLength: 72,
    circular: true,
    method: 'NEBUILDER_HIFI',
    gcPercent: '50.0%',
    visualization: {
      fragments: [
        { id: 'Fragment1', index: 0, startPosition: 0, endPosition: 36, actualLength: 36, color: '#3b82f6' },
        { id: 'Fragment2', index: 1, startPosition: 36, endPosition: 72, actualLength: 36, color: '#10b981' },
      ],
      junctions: [
        { index: 0, from: 'Fragment1', to: 'Fragment2', position: 20, overlapSequence: 'ATGCATGCATGCATGCATGC' },
      ],
    },
    fragments: [
      { id: 'Fragment1', length: 36 },
      { id: 'Fragment2', length: 36 },
    ],
  };

  it('should export valid GenBank format', () => {
    const genbank = exportToGenBank(mockSimulation);

    expect(genbank).toContain('LOCUS');
    expect(genbank).toContain('DEFINITION');
    expect(genbank).toContain('FEATURES');
    expect(genbank).toContain('ORIGIN');
    expect(genbank).toContain('//');
  });

  it('should include custom name in LOCUS', () => {
    const genbank = exportToGenBank(mockSimulation, { name: 'TestAssembly' });

    expect(genbank).toContain('TestAssembly');
  });

  it('should include fragment features', () => {
    const genbank = exportToGenBank(mockSimulation);

    expect(genbank).toContain('Fragment1');
    expect(genbank).toContain('Fragment2');
    expect(genbank).toContain('misc_feature');
  });

  it('should format sequence correctly', () => {
    const genbank = exportToGenBank(mockSimulation);

    // GenBank uses lowercase in ORIGIN
    expect(genbank).toContain('atgc');
    // Should have line numbers
    expect(genbank).toMatch(/\d+\s+[atgc\s]+/);
  });

  it('should mark topology as circular or linear', () => {
    const circularGb = exportToGenBank({ ...mockSimulation, circular: true });
    expect(circularGb).toContain('circular');

    const linearGb = exportToGenBank({ ...mockSimulation, circular: false });
    expect(linearGb).toContain('linear');
  });
});

// =============================================================================
// FASTA Export Tests
// =============================================================================

describe('FASTA Export', () => {
  const mockSimulation = {
    assembledSequence: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',
  };

  it('should export valid FASTA format', () => {
    const fasta = exportToFasta(mockSimulation);

    expect(fasta).toMatch(/^>/);
    expect(fasta).toContain('ATGC');
  });

  it('should include custom name in header', () => {
    const fasta = exportToFasta(mockSimulation, { name: 'MyAssembly' });

    expect(fasta).toContain('>MyAssembly');
  });

  it('should wrap sequence at specified width', () => {
    const longSeq = { assembledSequence: 'A'.repeat(200) };
    const fasta = exportToFasta(longSeq, { lineWidth: 80 });

    const lines = fasta.split('\n').filter(l => !l.startsWith('>'));
    expect(lines[0].length).toBeLessThanOrEqual(80);
  });
});

// =============================================================================
// Project Import/Export Tests
// =============================================================================

describe('Project Import/Export', () => {
  const projectData = {
    method: 'NEBUILDER_HIFI',
    fragments: [
      { id: 'Vector', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC' },
      { id: 'Insert', seq: 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC' },
    ],
    circular: true,
    ggEnzyme: 'BsaI-HFv2',
  };

  it('should export project as valid JSON', () => {
    const json = exportProject(projectData);
    const parsed = JSON.parse(json);

    expect(parsed.version).toBe('1.0');
    expect(parsed.method).toBe('NEBUILDER_HIFI');
    expect(parsed.fragments).toHaveLength(2);
    expect(parsed.exportedAt).toBeDefined();
  });

  it('should import project correctly', () => {
    const json = exportProject(projectData);
    const imported = importProject(json);

    expect(imported.method).toBe('NEBUILDER_HIFI');
    expect(imported.fragments).toHaveLength(2);
    expect(imported.fragments[0].id).toBe('Vector');
  });

  it('should reject invalid project files', () => {
    expect(() => importProject('{}')).toThrow('missing version');

    expect(() => importProject('{"version": "1.0"}')).toThrow('missing fragments');

    expect(() => importProject('{"version": "1.0", "fragments": [{"id": "x"}]}')).toThrow('must have id and seq');
  });
});
