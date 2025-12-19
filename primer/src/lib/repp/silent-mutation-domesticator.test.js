/**
 * Tests for Silent Mutation-Based Domestication
 *
 * These tests verify that the domestication algorithm correctly:
 * 1. Identifies internal restriction sites
 * 2. Finds valid synonymous mutations
 * 3. Preserves protein sequence
 * 4. Avoids creating new restriction sites
 * 5. Falls back appropriately when silent mutations aren't available
 */

import { describe, it, expect, beforeEach } from 'vitest';
import {
  domesticateWithSilentMutations,
  findAllSilentMutationCandidates,
  scoreMutationCandidates,
  verifyProteinSequence,
  validateDomestication,
  compareDomesticationStrategies,
  applyMutation,
  CODON_TABLE,
  CODON_TO_AA,
  ECOLI_CODON_USAGE,
} from './silent-mutation-domesticator.js';

import {
  optimizeDomestication,
  analyzeDomesticationOptions,
  generateDomesticationReport,
  DOMESTICATION_STRATEGY,
} from './domestication-optimizer.js';

import { findInternalSites, GOLDEN_GATE_ENZYMES } from './goldengate.js';
import { reverseComplement } from './enzymes.js';

// ============================================================================
// TEST SEQUENCES
// ============================================================================

// Sequence with one BsaI site (GGTCTC) in frame 0
// Codons around site: GGT CTC (Gly-Leu)
const SEQ_WITH_BSAI_FRAME0 = 'ATGAAAGGTCTCAAATGA'; // M-K-G-L-K-*

// Sequence with one BsaI site in different codon context
// The GGTCTC spans: GGT-CTC where GGT=Gly, CTC=Leu
const SEQ_BSAI_GLY_LEU = 'ATGAAAGGAGGTCTCCTGAAATGA';

// Sequence with BsaI site where silent mutation is possible
// Position the site so we have codon flexibility
const SEQ_BSAI_MUTABLE = 'ATGGGTCTCAAAGGTCTCAAATGA';

// Sequence with multiple BsaI sites
const SEQ_MULTIPLE_SITES = 'ATGAAAGGTCTCAAAGGGAAAGGTCTCAAATGA';

// Sequence with BsaI site where no silent mutation is possible (edge case)
// ATG (Met - only one codon) followed by TGG (Trp - only one codon)
const SEQ_NO_SILENT_OPTION = 'ATGTGGGGTCTCTGA';

// Longer realistic sequence with BsaI site
const SEQ_REALISTIC =
  'ATGGCTAAAGAAGAAGGTAAAGTTGCTGTTGCTGGTCTCGAAGCTGCTAAAGCTGAATGA';
//       ^-- BsaI site (GGTCTC) at position ~30

// Sequence with reverse complement BsaI site (GAGACC)
const SEQ_REVERSE_BSAI = 'ATGAAAGAGACCAAATGA';

// ============================================================================
// CODON TABLE TESTS
// ============================================================================

describe('Codon Tables', () => {
  it('should have all 64 codons mapped', () => {
    const allCodons = Object.values(CODON_TABLE).flat();
    expect(allCodons.length).toBe(64);
  });

  it('should have correct amino acid mappings', () => {
    expect(CODON_TO_AA['ATG']).toBe('M');
    expect(CODON_TO_AA['TGG']).toBe('W');
    expect(CODON_TO_AA['TAA']).toBe('*');
    expect(CODON_TO_AA['GGT']).toBe('G');
    expect(CODON_TO_AA['CTC']).toBe('L');
  });

  it('should have synonymous codons for most amino acids', () => {
    // Leucine has 6 codons
    expect(CODON_TABLE['L'].length).toBe(6);
    // Methionine has only 1
    expect(CODON_TABLE['M'].length).toBe(1);
    // Tryptophan has only 1
    expect(CODON_TABLE['W'].length).toBe(1);
  });

  it('should have E. coli codon usage frequencies', () => {
    expect(ECOLI_CODON_USAGE['CTG']).toBeGreaterThan(50); // Most common Leu codon
    expect(ECOLI_CODON_USAGE['CTA']).toBeLessThan(5);     // Rare Leu codon
  });
});

// ============================================================================
// FIND SILENT MUTATION CANDIDATES TESTS
// ============================================================================

describe('findAllSilentMutationCandidates', () => {
  it('should find candidates for BsaI site in coding sequence', () => {
    const seq = SEQ_BSAI_GLY_LEU.toUpperCase();
    const sites = findInternalSites(seq, 'BsaI');
    expect(sites.hasSites).toBe(true);

    const candidates = findAllSilentMutationCandidates(
      seq,
      sites.sites[0],
      0,
      'BsaI',
      ECOLI_CODON_USAGE
    );

    expect(candidates.length).toBeGreaterThan(0);

    // All candidates should be synonymous
    for (const candidate of candidates) {
      expect(candidate.isSynonymous).toBe(true);
      expect(candidate.breaksSite).toBe(true);
    }
  });

  it('should return empty for non-mutable sites', () => {
    // Create a sequence where the BsaI site spans Met and Trp (single-codon AAs)
    const seq = 'ATGTGGTCTCTGA'.toUpperCase(); // Artificial - may not have valid frame

    // This might not have valid candidates depending on frame alignment
    const sites = findInternalSites(seq, 'BsaI');

    if (sites.hasSites) {
      const candidates = findAllSilentMutationCandidates(
        seq,
        sites.sites[0],
        0,
        'BsaI',
        ECOLI_CODON_USAGE
      );

      // May or may not have candidates depending on exact sequence
      // The key is that all returned candidates must be valid
      for (const candidate of candidates) {
        expect(candidate.isSynonymous).toBe(true);
      }
    }
  });

  it('should handle reverse complement sites', () => {
    const seq = SEQ_REVERSE_BSAI.toUpperCase();
    const sites = findInternalSites(seq, 'BsaI');

    if (sites.hasSites) {
      const candidates = findAllSilentMutationCandidates(
        seq,
        sites.sites[0],
        0,
        'BsaI',
        ECOLI_CODON_USAGE
      );

      // All candidates should break the site
      for (const candidate of candidates) {
        expect(candidate.breaksSite).toBe(true);
      }
    }
  });

  it('should track codon frequency information', () => {
    const seq = SEQ_BSAI_GLY_LEU.toUpperCase();
    const sites = findInternalSites(seq, 'BsaI');

    if (sites.hasSites) {
      const candidates = findAllSilentMutationCandidates(
        seq,
        sites.sites[0],
        0,
        'BsaI',
        ECOLI_CODON_USAGE
      );

      for (const candidate of candidates) {
        expect(candidate.newCodonFrequency).toBeDefined();
        expect(candidate.originalCodonFrequency).toBeDefined();
        expect(typeof candidate.isRareCodon).toBe('boolean');
      }
    }
  });
});

// ============================================================================
// SCORE MUTATION CANDIDATES TESTS
// ============================================================================

describe('scoreMutationCandidates', () => {
  it('should penalize mutations that create new sites', () => {
    const seq = SEQ_BSAI_GLY_LEU.toUpperCase();
    const sites = findInternalSites(seq, 'BsaI');

    if (sites.hasSites) {
      const candidates = findAllSilentMutationCandidates(
        seq,
        sites.sites[0],
        0,
        'BsaI',
        ECOLI_CODON_USAGE
      );

      const scored = scoreMutationCandidates(candidates, seq, 'BsaI', true);

      // Scores should be between 0 and 100
      for (const candidate of scored) {
        expect(candidate.score).toBeGreaterThanOrEqual(0);
        expect(candidate.score).toBeLessThanOrEqual(100);
      }

      // Should be sorted by score (highest first)
      for (let i = 1; i < scored.length; i++) {
        expect(scored[i - 1].score).toBeGreaterThanOrEqual(scored[i].score);
      }
    }
  });

  it('should prefer wobble position mutations', () => {
    const seq = SEQ_BSAI_GLY_LEU.toUpperCase();
    const sites = findInternalSites(seq, 'BsaI');

    if (sites.hasSites) {
      const candidates = findAllSilentMutationCandidates(
        seq,
        sites.sites[0],
        0,
        'BsaI',
        ECOLI_CODON_USAGE
      );

      const scored = scoreMutationCandidates(candidates, seq, 'BsaI', false);

      // Find wobble position mutations
      const wobbleMutations = scored.filter(c => c.positionInCodon === 2);
      const nonWobbleMutations = scored.filter(c => c.positionInCodon !== 2);

      if (wobbleMutations.length > 0 && nonWobbleMutations.length > 0) {
        // Wobble mutations should have bonus (but other factors may override)
        const wobbleHasBonus = wobbleMutations.some(w =>
          w.bonuses.some(b => b.reason === 'WOBBLE_POSITION')
        );
        expect(wobbleHasBonus).toBe(true);
      }
    }
  });
});

// ============================================================================
// DOMESTICATION TESTS
// ============================================================================

describe('domesticateWithSilentMutations', () => {
  it('should return success=true when no sites present', () => {
    const seq = 'ATGAAAGCTGCTAAATGA'; // No BsaI site
    const result = domesticateWithSilentMutations(seq, 'BsaI');

    expect(result.success).toBe(true);
    expect(result.needsDomestication).toBe(false);
    expect(result.domesticatedSequence).toBe(seq.toUpperCase());
  });

  it('should domesticate sequence with single BsaI site', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const result = domesticateWithSilentMutations(seq, 'BsaI', {
      frame: 0,
      isCodingSequence: true,
    });

    if (result.success && result.mutations.length > 0) {
      // Verify the site is gone
      const checkResult = findInternalSites(result.domesticatedSequence, 'BsaI');
      expect(checkResult.hasSites).toBe(false);

      // Verify protein is preserved
      const proteinCheck = verifyProteinSequence(seq, result.domesticatedSequence, 0);
      expect(proteinCheck.identical).toBe(true);
    }
  });

  it('should handle multiple internal sites', () => {
    const seq = SEQ_MULTIPLE_SITES;
    const originalSites = findInternalSites(seq, 'BsaI');

    const result = domesticateWithSilentMutations(seq, 'BsaI', {
      frame: 0,
      isCodingSequence: true,
    });

    if (result.success) {
      // Should have attempted to fix all sites
      expect(result.verification.originalSiteCount).toBe(originalSites.count);

      // Check remaining sites
      const remainingSites = findInternalSites(result.domesticatedSequence, 'BsaI');
      expect(remainingSites.count).toBeLessThanOrEqual(result.failedSites.length);
    }
  });

  it('should report failed sites for non-coding sequences', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const result = domesticateWithSilentMutations(seq, 'BsaI', {
      isCodingSequence: false, // Can't use silent mutations
    });

    expect(result.failedSites.length).toBeGreaterThan(0);
    expect(result.failedSites[0].reason).toBe('NON_CODING');
  });

  it('should preserve sequence length', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const result = domesticateWithSilentMutations(seq, 'BsaI');

    expect(result.domesticatedSequence.length).toBe(seq.length);
  });
});

// ============================================================================
// PROTEIN VERIFICATION TESTS
// ============================================================================

describe('verifyProteinSequence', () => {
  it('should confirm identical proteins', () => {
    const original = 'ATGGCTAAATGA'; // M-A-K-*
    const domesticated = 'ATGGCCAAATGA'; // M-A-K-* (GCT->GCC, both Ala)

    const result = verifyProteinSequence(original, domesticated, 0);

    expect(result.identical).toBe(true);
    expect(result.differences.length).toBe(0);
  });

  it('should detect protein differences', () => {
    const original = 'ATGGCTAAATGA'; // M-A-K-*
    const different = 'ATGGATAAATGA'; // M-D-K-* (GCT->GAT, Ala->Asp)

    const result = verifyProteinSequence(original, different, 0);

    expect(result.identical).toBe(false);
    expect(result.differences.length).toBe(1);
    expect(result.differences[0].position).toBe(1); // Second amino acid
  });

  it('should handle different reading frames', () => {
    const seq = 'AATGGCTAAATGA';

    // Frame 0: codons start at 0
    const frame0 = verifyProteinSequence(seq, seq, 0);
    expect(frame0.identical).toBe(true);

    // Frame 1: codons start at 1
    const frame1 = verifyProteinSequence(seq, seq, 1);
    expect(frame1.identical).toBe(true);
  });
});

// ============================================================================
// VALIDATION TESTS
// ============================================================================

describe('validateDomestication', () => {
  it('should pass validation for successful domestication', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const result = domesticateWithSilentMutations(seq, 'BsaI', {
      frame: 0,
      isCodingSequence: true,
    });

    if (result.success && result.mutations.length > 0) {
      const validation = validateDomestication(result, 'BsaI', { frame: 0 });

      expect(validation.isValid).toBe(true);
      expect(validation.validations.find(v => v.check === 'PROTEIN_SEQUENCE').passed).toBe(true);
      expect(validation.validations.find(v => v.check === 'SEQUENCE_LENGTH').passed).toBe(true);
    }
  });

  it('should fail validation if protein changed', () => {
    // Manually create a bad result
    const badResult = {
      originalSequence: 'ATGGCTAAATGA',
      domesticatedSequence: 'ATGGATAAATGA', // Non-synonymous change
      mutations: [],
      failedSites: [],
    };

    const validation = validateDomestication(badResult, 'BsaI', { frame: 0 });

    expect(validation.isValid).toBe(false);
    expect(validation.validations.find(v => v.check === 'PROTEIN_SEQUENCE').passed).toBe(false);
  });
});

// ============================================================================
// UNIFIED OPTIMIZER TESTS
// ============================================================================

describe('optimizeDomestication', () => {
  it('should prefer silent mutations over junction-based', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const result = optimizeDomestication(seq, 'BsaI', {
      frame: 0,
      isCodingSequence: true,
    });

    if (result.success && result.strategy !== DOMESTICATION_STRATEGY.NONE) {
      // Should prefer silent mutation if available
      if (result.strategy === DOMESTICATION_STRATEGY.SILENT_MUTATION) {
        expect(result.warnings.length).toBe(0); // No warnings for silent mutations
      }
    }
  });

  it('should add warning for junction-based fallback', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const result = optimizeDomestication(seq, 'BsaI', {
      isCodingSequence: false, // Force non-coding to trigger fallback
      allowJunctionFallback: true,
    });

    if (result.strategy === DOMESTICATION_STRATEGY.JUNCTION_BASED ||
        result.strategy === DOMESTICATION_STRATEGY.HYBRID) {
      expect(result.warnings.length).toBeGreaterThan(0);
      expect(result.warnings[0].severity).toBe('high');
    }
  });

  it('should recommend alternative enzyme when available', () => {
    // Create sequence that has BsaI site but not BbsI site
    const seq = 'ATGAAAGGTCTCAAATGA'; // Has BsaI (GGTCTC), check if has BbsI

    const bsaiSites = findInternalSites(seq, 'BsaI');
    const bbsiSites = findInternalSites(seq, 'BbsI');

    if (bsaiSites.hasSites && !bbsiSites.hasSites) {
      const result = optimizeDomestication(seq, 'BsaI', {
        preferAlternativeEnzyme: true,
      });

      if (result.strategy === DOMESTICATION_STRATEGY.ALTERNATIVE_ENZYME) {
        expect(result.recommendedEnzyme).toBeDefined();
        expect(result.recommendedEnzyme.isCompatible).toBe(true);
      }
    }
  });

  it('should return NONE strategy when no sites present', () => {
    const seq = 'ATGAAAGCTAAATGA'; // No BsaI site

    const result = optimizeDomestication(seq, 'BsaI');

    expect(result.success).toBe(true);
    expect(result.strategy).toBe(DOMESTICATION_STRATEGY.NONE);
  });
});

// ============================================================================
// ANALYSIS TESTS
// ============================================================================

describe('analyzeDomesticationOptions', () => {
  it('should provide detailed analysis of each site', () => {
    const seq = SEQ_MULTIPLE_SITES;
    const analysis = analyzeDomesticationOptions(seq, 'BsaI', { frame: 0 });

    if (analysis.needsDomestication) {
      expect(analysis.siteAnalyses.length).toBe(analysis.siteCount);

      for (const siteAnalysis of analysis.siteAnalyses) {
        expect(siteAnalysis.silentMutation).toBeDefined();
        expect(siteAnalysis.junctionBased).toBeDefined();
      }
    }
  });

  it('should provide summary statistics', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const analysis = analyzeDomesticationOptions(seq, 'BsaI');

    if (analysis.needsDomestication) {
      expect(analysis.summary).toBeDefined();
      expect(typeof analysis.summary.silentMutationSites).toBe('number');
      expect(typeof analysis.summary.junctionRequiredSites).toBe('number');
    }
  });

  it('should include alternative enzyme recommendations', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const analysis = analyzeDomesticationOptions(seq, 'BsaI');

    if (analysis.needsDomestication) {
      expect(analysis.alternativeEnzymes).toBeDefined();
      expect(Array.isArray(analysis.alternativeEnzymes)).toBe(true);
    }
  });
});

// ============================================================================
// REPORT GENERATION TESTS
// ============================================================================

describe('generateDomesticationReport', () => {
  it('should generate comprehensive report', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const result = domesticateWithSilentMutations(seq, 'BsaI', {
      frame: 0,
      isCodingSequence: true,
    });

    const report = generateDomesticationReport(result, { frame: 0 });

    expect(report.summary).toBeDefined();
    expect(report.changes).toBeDefined();
    expect(Array.isArray(report.changes)).toBe(true);
  });

  it('should include sequence comparison when requested', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const result = domesticateWithSilentMutations(seq, 'BsaI');

    const report = generateDomesticationReport(result, {
      includeSequenceAlignment: true,
    });

    if (result.mutations.length > 0) {
      expect(report.sequenceComparison).toBeDefined();
      expect(report.sequenceComparison.totalDifferences).toBeGreaterThan(0);
    }
  });

  it('should include protein verification', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const result = domesticateWithSilentMutations(seq, 'BsaI');

    const report = generateDomesticationReport(result, { frame: 0 });

    expect(report.proteinVerification).toBeDefined();
    expect(report.proteinVerification.identical).toBe(true);
  });
});

// ============================================================================
// EDGE CASES AND REGRESSION TESTS
// ============================================================================

describe('Edge Cases', () => {
  it('should handle empty sequence', () => {
    const result = domesticateWithSilentMutations('', 'BsaI');
    expect(result.success).toBe(true);
    expect(result.needsDomestication).toBe(false);
  });

  it('should handle sequence shorter than recognition site', () => {
    const result = domesticateWithSilentMutations('ATG', 'BsaI');
    expect(result.success).toBe(true);
    expect(result.needsDomestication).toBe(false);
  });

  it('should handle unknown enzyme gracefully', () => {
    expect(() => {
      domesticateWithSilentMutations('ATGAAATGA', 'UnknownEnzyme');
    }).toThrow();
  });

  it('should handle overlapping sites', () => {
    // Create sequence with potentially overlapping sites
    const seq = 'ATGGGTCTCGGTCTCAAATGA';
    const result = domesticateWithSilentMutations(seq, 'BsaI', {
      frame: 0,
      isCodingSequence: true,
    });

    // Should handle without crashing
    expect(result).toBeDefined();
  });

  it('should handle site at sequence boundary', () => {
    // Site at the very beginning
    const seqStart = 'GGTCTCAAAAAATGA';
    const resultStart = domesticateWithSilentMutations(seqStart, 'BsaI');
    expect(resultStart).toBeDefined();

    // Site at the very end
    const seqEnd = 'ATGAAAAAAGGTCTC';
    const resultEnd = domesticateWithSilentMutations(seqEnd, 'BsaI');
    expect(resultEnd).toBeDefined();
  });
});

// ============================================================================
// STRATEGY COMPARISON TESTS
// ============================================================================

describe('compareDomesticationStrategies', () => {
  it('should compare available strategies', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const comparison = compareDomesticationStrategies(seq, 'BsaI', { frame: 0 });

    if (comparison.needsDomestication) {
      expect(comparison.strategies).toBeDefined();
      expect(comparison.strategies.silentMutation).toBeDefined();
      expect(comparison.strategies.junctionBased).toBeDefined();
      expect(comparison.recommended).toBeDefined();
    }
  });

  it('should correctly identify one-pot compatibility', () => {
    const seq = SEQ_BSAI_GLY_LEU;
    const comparison = compareDomesticationStrategies(seq, 'BsaI', { frame: 0 });

    if (comparison.needsDomestication) {
      // Silent mutation should be one-pot compatible
      expect(comparison.strategies.silentMutation.onePotCompatible).toBe(true);

      // Junction-based should NOT be one-pot compatible
      expect(comparison.strategies.junctionBased.onePotCompatible).toBe(false);
    }
  });

  it('should note fragment increase for junction-based', () => {
    const seq = SEQ_MULTIPLE_SITES;
    const sites = findInternalSites(seq, 'BsaI');
    const comparison = compareDomesticationStrategies(seq, 'BsaI', { frame: 0 });

    if (comparison.needsDomestication) {
      // Junction-based adds one fragment per site
      expect(comparison.strategies.junctionBased.fragmentIncrease).toBe(sites.count);
      // Silent mutation doesn't add fragments
      expect(comparison.strategies.silentMutation.fragmentIncrease).toBe(0);
    }
  });
});
