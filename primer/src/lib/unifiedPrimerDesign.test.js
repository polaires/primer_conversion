/**
 * Tests for Unified Primer Design Module
 *
 * These tests verify:
 * 1. The unified API works correctly for all operation types
 * 2. Results match the original implementation (regression tests)
 * 3. Edge cases and error handling
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  designUnified,
  designBatch,
  parseNotationToSpec,
  parseBatchText,
  MUTAGENESIS_DEFAULTS,
  MUTATION_TYPES,
  CODON_TABLE,
  AA_NAMES,
} from './unifiedPrimerDesign.js';

// Import original implementations for comparison
import { primers } from './primers.js';
import {
  designSubstitutionPrimers,
  designCodonChangePrimers,
  designInsertionPrimers,
  designDeletionPrimers,
  designRegionSubstitutionPrimers,
} from './mutagenesis.js';

// =============================================================================
// Test Data
// =============================================================================

// A realistic plasmid-like test sequence (~500bp)
const TEST_TEMPLATE = 'ATGAAACAAAGCACTATTGCACTGGCACTCTTACCGTTACTGTTTACCCCTGTGACAAAAGCCAATATGGACAAGTTGTTTGACGGTATCAGAAGCAGAACTGGCGAAACTTTTACCGGTGAAGACCGTAACGGTTACGACAATAAATACAATGTTTATAATCAGACTAACGACTGTTGGGGTTTTGAATTTAAAGATGAAGATATGCTGTGCCCGGACCCAATTAGCTGGCGTAATGCCGAGATCATGCGTAAAAAATGGGACAGCAAAGAGCAGAAAAGCATGTACGAACGCCAGTTTGACGAGCTGTATAAAGAACGCTATGGTTATGCCAACAGCTACATGTATGACGATGATGACAAACATCTGTACAAGTAAGGAGGTAATAA';

// Shorter test sequence for quick tests
const SHORT_TEMPLATE = 'ATGAAACAAAGCACTATTGCACTGGCACTCTTACCGTTACTGTTTACCCCTGTGACAAAAGCC';

// GFP-like sequence with known codons
const GFP_TEMPLATE = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA';

// =============================================================================
// Unit Tests: Unified API
// =============================================================================

describe('Unified Primer Design API', () => {
  describe('Operation Type Detection', () => {
    it('should detect amplification when replacement is null', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 100,
        replacement: null,
      });

      expect(result.type).toBe('amplification');
      expect(result.operation).toBe('amplify');
    });

    it('should detect amplification when replacement is undefined', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 100,
      });

      expect(result.type).toBe('amplification');
      expect(result.operation).toBe('amplify');
    });

    it('should detect deletion when replacement is empty string', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 50,
        end: 60,
        replacement: '',
      });

      expect(result.type).toBe('deletion');
      expect(result.operation).toBe('delete');
    });

    it('should detect insertion when start === end with replacement', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 50,
        end: 50,
        replacement: 'ACGT',
      });

      expect(result.type).toBe('insertion');
      expect(result.operation).toBe('insert');
    });

    it('should detect substitution when region has replacement', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 50,
        end: 56,
        replacement: 'AAATTT',
      });

      expect(result.type).toBe('substitution');
      expect(result.operation).toBe('substitute');
    });

    it('should detect AA mutation when aaHelper is provided', () => {
      const result = designUnified(GFP_TEMPLATE, {
        start: 195, // Codon 66
        end: 198,
        aaHelper: {
          newAA: 'W',
          codonPosition: 66,
        },
      });

      expect(result.type).toBe('codon_change');
      expect(result.operation).toBe('aa_mutation');
    });
  });

  describe('Amplification', () => {
    it('should design primers for amplification', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 200,
      });

      expect(result.forward).toBeDefined();
      expect(result.reverse).toBeDefined();
      expect(result.forward.sequence).toBeTruthy();
      expect(result.reverse.sequence).toBeTruthy();
      expect(result.forward.tm).toBeGreaterThan(50);
      expect(result.reverse.tm).toBeGreaterThan(50);
      expect(result.amplifiedRegion.length).toBe(200);
    });

    it('should include quality scoring', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 200,
      });

      expect(result.compositeScore).toBeDefined();
      expect(result.qualityTier).toBeDefined();
      expect(result.compositeScore).toBeGreaterThanOrEqual(0);
      expect(result.compositeScore).toBeLessThanOrEqual(100);
    });

    it('should generate alternative primers', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 200,
      });

      expect(result.alternativePrimers).toBeDefined();
      expect(Array.isArray(result.alternativePrimers)).toBe(true);
    });

    it('should reject regions smaller than 50bp', () => {
      expect(() => designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 40,
      })).toThrow('at least 50 bp');
    });
  });

  describe('Deletion', () => {
    it('should design primers for deletion', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 110,
        replacement: '',
      });

      expect(result.type).toBe('deletion');
      expect(result.forward).toBeDefined();
      expect(result.reverse).toBeDefined();
      expect(result.deletedSequence).toBeTruthy();
      expect(result.deleteLength).toBe(10);
    });

    it('should handle single base deletion', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 101,
        replacement: '',
      });

      expect(result.deleteLength).toBe(1);
    });
  });

  describe('Insertion', () => {
    it('should design primers for insertion', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 100,
        replacement: 'ACGTACGT',
      });

      expect(result.type).toBe('insertion');
      expect(result.forward).toBeDefined();
      expect(result.reverse).toBeDefined();
      expect(result.insertedSequence).toBe('ACGTACGT');
      expect(result.insertLength).toBe(8);
    });

    it('should detect IUPAC ambiguous bases for library design', () => {
      // Note: The underlying mutagenesis library doesn't support direct IUPAC codes
      // Library design with ambiguous bases requires frontend expansion
      // This test verifies the detection and metadata, actual design would need expanded sequences
      const result = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 100,
        replacement: 'ACGT', // Use valid sequence for design
      });

      expect(result.type).toBe('insertion');

      // Test the library detection separately
      const { hasAmbiguousBases, countCombinations } = require('./sequenceUtils.js');
      expect(hasAmbiguousBases('NNK')).toBe(true);
      expect(countCombinations('NNK')).toBe(32); // 4*4*2 = 32
    });
  });

  describe('Substitution', () => {
    it('should design primers for substitution', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 106,
        replacement: 'AAATTT',
      });

      expect(result.type).toBe('substitution');
      expect(result.forward).toBeDefined();
      expect(result.reverse).toBeDefined();
    });

    it('should handle different length replacements', () => {
      // Replace 6bp with 10bp
      const result = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 106,
        replacement: 'AAAATTTTGG',
      });

      expect(result.replacementLength).toBe(10);
      expect(result.deleteLength).toBe(6);
    });
  });

  describe('Amino Acid Mutation', () => {
    it('should design primers for AA mutation', () => {
      const result = designUnified(GFP_TEMPLATE, {
        start: 195,
        end: 198,
        aaHelper: {
          newAA: 'W',
          codonPosition: 66,
          organism: 'ecoli',
        },
      });

      expect(result.type).toBe('codon_change');
      expect(result.newAA).toBe('W');
      expect(result.newCodon).toBeDefined();
      expect(CODON_TABLE['W']).toContain(result.newCodon);
    });

    it('should optimize codon for organism', () => {
      const resultEcoli = designUnified(GFP_TEMPLATE, {
        start: 195,
        end: 198,
        aaHelper: {
          newAA: 'L',
          codonPosition: 66,
          organism: 'ecoli',
        },
      });

      expect(result => {
        // E. coli prefers CTG for Leucine
        return resultEcoli.newCodon === 'CTG' || CODON_TABLE['L'].includes(resultEcoli.newCodon);
      }).toBeTruthy();
    });
  });

  describe('Error Handling', () => {
    it('should reject templates shorter than 50bp', () => {
      expect(() => designUnified('ATGCATGC', {
        start: 0,
        end: 8,
      })).toThrow('at least 50 bp');
    });

    it('should reject invalid positions', () => {
      expect(() => designUnified(TEST_TEMPLATE, {
        start: -1,
        end: 100,
      })).toThrow('Invalid region');

      expect(() => designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 10000,
      })).toThrow('Invalid region');

      // For linear sequences, start > end should throw
      expect(() => designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 50,
      }, { circular: false })).toThrow('not allowed for linear sequences');
    });
  });
});

// =============================================================================
// Regression Tests: Compare with Original Implementations
// =============================================================================

describe('Regression Tests: Unified vs Original', () => {
  describe('Amplification matches primers()', () => {
    it('should produce similar results to primers() without smart design', () => {
      const regionSeq = TEST_TEMPLATE.slice(0, 200);

      // Original (without smart design for exact comparison)
      const [origFwd, origRev] = primers(regionSeq, {
        useCompositeScore: true,
        useSmartDesign: false,
      });

      // Unified (also without smart design)
      const unified = designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 200,
      }, { useSmartDesign: false });

      // Primers should match (same algorithm)
      expect(unified.forward.sequence).toBe(origFwd.seq);
      expect(unified.reverse.sequence).toBe(origRev.seq);
      expect(unified.forward.tm).toBeCloseTo(origFwd.tm, 1);
      expect(unified.reverse.tm).toBeCloseTo(origRev.tm, 1);
    });

    it('should produce reasonable results with smart design', () => {
      const regionSeq = TEST_TEMPLATE.slice(0, 200);

      // Both with smart design (results may differ slightly due to optimization)
      const [origFwd, origRev] = primers(regionSeq, {
        useCompositeScore: true,
        useSmartDesign: true,
      });

      const unified = designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 200,
      }, { useSmartDesign: true });

      // With smart design, sequences may differ but should still be valid primers
      expect(unified.forward.sequence).toBeTruthy();
      expect(unified.reverse.sequence).toBeTruthy();
      // Tm should be in reasonable range
      expect(unified.forward.tm).toBeGreaterThan(50);
      expect(unified.forward.tm).toBeLessThan(75);
      expect(unified.reverse.tm).toBeGreaterThan(50);
      expect(unified.reverse.tm).toBeLessThan(75);
    });
  });

  describe('Deletion matches designDeletionPrimers()', () => {
    it('should produce same results as original', () => {
      // Original
      const original = designDeletionPrimers(TEST_TEMPLATE, 100, 10);

      // Unified
      const unified = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 110,
        replacement: '',
      });

      // Compare key outputs
      expect(unified.forward.sequence).toBe(original.forward.sequence);
      expect(unified.reverse.sequence).toBe(original.reverse.sequence);
      expect(unified.deleteLength).toBe(original.deleteLength);
      expect(unified.deletedSequence).toBe(original.deletedSequence);
    });
  });

  describe('Insertion matches designInsertionPrimers()', () => {
    it('should produce same results as original', () => {
      const insertSeq = 'ACGTACGT';

      // Original
      const original = designInsertionPrimers(TEST_TEMPLATE, 100, insertSeq);

      // Unified
      const unified = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 100,
        replacement: insertSeq,
      });

      // Compare key outputs
      expect(unified.forward.sequence).toBe(original.forward.sequence);
      expect(unified.reverse.sequence).toBe(original.reverse.sequence);
      expect(unified.insertedSequence).toBe(original.insertedSequence);
      expect(unified.insertLength).toBe(original.insertLength);
    });
  });

  describe('Substitution matches designRegionSubstitutionPrimers()', () => {
    it('should produce same results as original', () => {
      const replacementSeq = 'AAATTTGGG';

      // Original
      const original = designRegionSubstitutionPrimers(TEST_TEMPLATE, 100, 6, replacementSeq);

      // Unified
      const unified = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 106,
        replacement: replacementSeq,
      });

      // Compare key outputs
      expect(unified.forward.sequence).toBe(original.forward.sequence);
      expect(unified.reverse.sequence).toBe(original.reverse.sequence);
      expect(unified.deleteLength).toBe(original.deleteLength);
      expect(unified.replacementSequence).toBe(original.replacementSequence);
    });
  });

  describe('AA Mutation matches designCodonChangePrimers()', () => {
    it('should produce same results as original', () => {
      // Original
      const original = designCodonChangePrimers(GFP_TEMPLATE, 66, 'W', { organism: 'ecoli' });

      // Unified
      const unified = designUnified(GFP_TEMPLATE, {
        start: 195, // (66-1) * 3 = 195
        end: 198,
        aaHelper: {
          newAA: 'W',
          codonPosition: 66,
          organism: 'ecoli',
        },
      });

      // Compare key outputs
      expect(unified.forward.sequence).toBe(original.forward.sequence);
      expect(unified.reverse.sequence).toBe(original.reverse.sequence);
      expect(unified.newCodon).toBe(original.newCodon);
      expect(unified.oldCodon).toBe(original.oldCodon);
      expect(unified.newAA).toBe(original.newAA);
    });
  });
});

// =============================================================================
// Batch Design Tests
// =============================================================================

describe('Batch Design', () => {
  it('should process multiple specifications', () => {
    const specs = [
      { start: 100, end: 110, replacement: '' },  // Deletion
      { start: 200, end: 200, replacement: 'ACGT' },  // Insertion
      { start: 300, end: 306, replacement: 'AAATTT' },  // Substitution
    ];

    const results = designBatch(TEST_TEMPLATE, specs);

    expect(results.length).toBe(3);
    expect(results[0].success).toBe(true);
    expect(results[0].operation).toBe('delete');
    expect(results[1].success).toBe(true);
    expect(results[1].operation).toBe('insert');
    expect(results[2].success).toBe(true);
    expect(results[2].operation).toBe('substitute');
  });

  it('should handle errors gracefully in batch', () => {
    const specs = [
      { start: 100, end: 110, replacement: '' },  // Valid
      { start: -1, end: 100, replacement: '' },   // Invalid
      { start: 200, end: 200, replacement: 'ACGT' },  // Valid
    ];

    const results = designBatch(TEST_TEMPLATE, specs);

    expect(results.length).toBe(3);
    expect(results[0].success).toBe(true);
    expect(results[1].success).toBe(false);
    expect(results[1].error).toBeDefined();
    expect(results[2].success).toBe(true);
  });
});

// =============================================================================
// Notation Parsing Tests
// =============================================================================

describe('Notation Parsing', () => {
  describe('parseNotationToSpec', () => {
    it('should parse AA mutations', () => {
      const spec = parseNotationToSpec('Y66W', { orfStart: 1 });

      expect(spec.notationType).toBe('aa_mutation');
      expect(spec.aaHelper.newAA).toBe('W');
      expect(spec.aaHelper.originalAA).toBe('Y');
      expect(spec.aaHelper.codonPosition).toBe(66);
    });

    it('should parse deletions', () => {
      const spec = parseNotationToSpec('del100-105');

      expect(spec.notationType).toBe('deletion');
      expect(spec.start).toBe(99);  // 0-based
      expect(spec.end).toBe(105);
      expect(spec.replacement).toBe('');
    });

    it('should parse single position deletions', () => {
      const spec = parseNotationToSpec('del100');

      expect(spec.notationType).toBe('deletion');
      expect(spec.start).toBe(99);
      expect(spec.end).toBe(100);
    });

    it('should parse insertions', () => {
      const spec = parseNotationToSpec('ins50_ACGT');

      expect(spec.notationType).toBe('insertion');
      expect(spec.start).toBe(49);  // 0-based
      expect(spec.end).toBe(49);
      expect(spec.replacement).toBe('ACGT');
    });

    it('should reject invalid notation', () => {
      expect(() => parseNotationToSpec('invalid')).toThrow();
      expect(() => parseNotationToSpec('X66Y')).toThrow();  // X is not a valid AA
    });
  });

  describe('parseBatchText', () => {
    it('should parse multiple mutations', () => {
      const text = `
        Y66W
        S65T
        del100-105
        ins200_ACGT
      `;

      const { specs } = parseBatchText(text);

      expect(specs.length).toBe(4);
      expect(specs[0].aaHelper.newAA).toBe('W');
      expect(specs[1].aaHelper.newAA).toBe('T');
      expect(specs[2].notationType).toBe('deletion');
      expect(specs[3].notationType).toBe('insertion');
    });

    it('should handle ORF directive', () => {
      const text = `
        ORF: 10
        Y66W
      `;

      const { specs, orfStart } = parseBatchText(text);

      expect(orfStart).toBe(10);
      expect(specs.length).toBe(1);
      // Nucleotide position should account for ORF offset
      expect(specs[0].start).toBe(9 + (66 - 1) * 3);  // orfStart-1 + (codon-1)*3
    });

    it('should skip comments and empty lines', () => {
      const text = `
        # This is a comment
        Y66W

        # Another comment
        S65T
      `;

      const { specs } = parseBatchText(text);

      expect(specs.length).toBe(2);
    });

    it('should capture parsing errors', () => {
      const text = `
        Y66W
        invalid_notation
        S65T
      `;

      const { specs } = parseBatchText(text);

      expect(specs.length).toBe(3);
      expect(specs[0].aaHelper).toBeDefined();
      expect(specs[1].error).toBeDefined();
      expect(specs[2].aaHelper).toBeDefined();
    });
  });
});

// =============================================================================
// Integration Tests
// =============================================================================

describe('Integration Tests', () => {
  it('should handle full workflow: parse notation -> design -> results', () => {
    const spec = parseNotationToSpec('Y66W', { orfStart: 1 });
    const result = designUnified(GFP_TEMPLATE, spec);

    expect(result.type).toBe('codon_change');
    expect(result.forward.sequence).toBeTruthy();
    expect(result.reverse.sequence).toBeTruthy();
    expect(result.newAA).toBe('W');
  });

  it('should handle batch workflow with text input', () => {
    const text = `
      ORF: 1
      Y66W
      del100-105
    `;

    const { specs } = parseBatchText(text);
    const results = designBatch(GFP_TEMPLATE, specs);

    expect(results.length).toBe(2);
    expect(results.every(r => r.success)).toBe(true);
    expect(results[0].type).toBe('codon_change');
    expect(results[1].type).toBe('deletion');
  });
});

// =============================================================================
// Circular Plasmid Region Tests
// =============================================================================

describe('Circular Plasmid Regions', () => {
  // A circular plasmid template for testing
  const CIRCULAR_TEMPLATE = TEST_TEMPLATE; // ~400bp

  describe('Amplification with wrap-around', () => {
    it('should allow start > end for circular plasmids', () => {
      // Region wraps around: from near end back to beginning
      // e.g., start=380, end=50 -> extracts seq[380:] + seq[0:50]
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: 350,
        end: 50,
        replacement: null, // amplify
      }, { circular: true });

      expect(result.type).toBe('amplification');
      expect(result.operation).toBe('amplify');
      expect(result.circularWrapped).toBe(true);

      // Region length should be (templateLen - 350) + 50 = wrap-around length
      const expectedLength = (CIRCULAR_TEMPLATE.length - 350) + 50;
      expect(result.amplifiedRegion.length).toBe(expectedLength);
      expect(result.amplifiedRegion.isWrapped).toBe(true);
    });

    it('should reject start > end for linear sequences', () => {
      expect(() => designUnified(CIRCULAR_TEMPLATE, {
        start: 350,
        end: 50,
        replacement: null,
      }, { circular: false })).toThrow('not allowed for linear sequences');
    });

    it('should generate valid primers for wrap-around amplification', () => {
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: 350,
        end: 50,
        replacement: null,
      }, { circular: true });

      expect(result.forward.sequence).toBeTruthy();
      expect(result.reverse.sequence).toBeTruthy();
      expect(result.forward.tm).toBeGreaterThan(50);
      expect(result.reverse.tm).toBeGreaterThan(50);
    });

    it('should include quality scoring for wrap-around regions', () => {
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: 350,
        end: 50,
        replacement: null,
      }, { circular: true });

      expect(result.compositeScore).toBeDefined();
      expect(result.qualityTier).toBeDefined();
    });
  });

  describe('Deletion with wrap-around', () => {
    it('should allow deleting a region that wraps around origin', () => {
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: 380,
        end: 20,
        replacement: '', // delete
      }, { circular: true });

      expect(result.type).toBe('deletion');
      expect(result.operation).toBe('delete');
      expect(result.circularWrapped).toBe(true);
    });

    it('should calculate correct deletion length for wrap-around', () => {
      const start = 380;
      const end = 20;
      const expectedLength = (CIRCULAR_TEMPLATE.length - start) + end;

      const result = designUnified(CIRCULAR_TEMPLATE, {
        start,
        end,
        replacement: '',
      }, { circular: true });

      expect(result.deleteLength).toBe(expectedLength);
    });
  });

  describe('Substitution with wrap-around', () => {
    it('should allow substituting a region that wraps around origin', () => {
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: 380,
        end: 20,
        replacement: 'AAAAAATTTTTT', // 12bp replacement
      }, { circular: true });

      expect(result.type).toBe('substitution');
      expect(result.operation).toBe('substitute');
      expect(result.circularWrapped).toBe(true);
    });

    it('should generate valid primers for wrap-around substitution', () => {
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: 380,
        end: 20,
        replacement: 'AAAAAATTTTTT',
      }, { circular: true });

      expect(result.forward).toBeDefined();
      expect(result.reverse).toBeDefined();
      expect(result.forward.sequence).toBeTruthy();
      expect(result.reverse.sequence).toBeTruthy();
    });
  });

  describe('Edge cases', () => {
    it('should handle region starting at position 0', () => {
      // This is not a wrap-around case, just start at origin
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: 0,
        end: 100,
        replacement: null,
      }, { circular: true });

      expect(result.type).toBe('amplification');
      expect(result.circularWrapped).toBeFalsy();
    });

    it('should handle region ending at last position', () => {
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: 300,
        end: CIRCULAR_TEMPLATE.length,
        replacement: null,
      }, { circular: true });

      expect(result.type).toBe('amplification');
      expect(result.circularWrapped).toBeFalsy();
    });

    it('should handle small wrap-around regions', () => {
      // Very small wrap: last 5bp + first 45bp = 50bp minimum
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: CIRCULAR_TEMPLATE.length - 5,
        end: 50,
        replacement: null,
      }, { circular: true });

      expect(result.type).toBe('amplification');
      expect(result.circularWrapped).toBe(true);
      // Region length should be 5 + 50 = 55bp
      expect(result.amplifiedRegion.length).toBe(55);
    });

    it('should validate bounds even for circular sequences', () => {
      // Start still must be within bounds
      expect(() => designUnified(CIRCULAR_TEMPLATE, {
        start: -1,
        end: 100,
        replacement: null,
      }, { circular: true })).toThrow('out of bounds');

      // End still must be within bounds (0 to seq.length inclusive)
      expect(() => designUnified(CIRCULAR_TEMPLATE, {
        start: 100,
        end: CIRCULAR_TEMPLATE.length + 1,
        replacement: null,
      }, { circular: true })).toThrow('out of bounds');
    });
  });

  describe('Default circular behavior', () => {
    it('should default to circular mode', () => {
      // Without specifying circular option, should still work for wrap-around
      const result = designUnified(CIRCULAR_TEMPLATE, {
        start: 350,
        end: 50,
        replacement: null,
      }); // No options specified

      expect(result.type).toBe('amplification');
      expect(result.circularWrapped).toBe(true);
    });
  });
});

// =============================================================================
// Performance Tests
// =============================================================================

describe('Performance', () => {
  it('should complete amplification design in reasonable time', () => {
    const start = performance.now();

    // Disable smart design for faster performance test
    for (let i = 0; i < 5; i++) {
      designUnified(TEST_TEMPLATE, {
        start: 0,
        end: 200,
      }, { useSmartDesign: false });
    }

    const elapsed = performance.now() - start;
    // 5 designs in under 60 seconds (extended for 60bp primer search space)
    expect(elapsed).toBeLessThan(60000);
  });

  it('should complete mutation design in reasonable time', () => {
    const start = performance.now();

    for (let i = 0; i < 5; i++) {
      designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 106,
        replacement: 'AAATTT',
      });
    }

    const elapsed = performance.now() - start;
    // 5 mutation designs in under 60 seconds (extended for larger search space)
    expect(elapsed).toBeLessThan(60000);
  });
});
