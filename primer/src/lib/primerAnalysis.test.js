/**
 * Tests for primerAnalysis.js - Assembly Primer Analysis
 *
 * Tests the assembly-specific primer analysis that correctly identifies
 * annealing regions and scores based on template binding, plus Golden Gate
 * enzyme site auto-detection.
 */

import { describe, it, expect } from 'vitest';
import {
  findAnnealingRegion,
  analyzeAssemblyPrimer,
  analyzeSinglePrimer,
  analyzePrimers,
  detectGoldenGateSites,
} from './primerAnalysis.js';

// Test data: A typical Golden Gate primer and its template
const GOLDEN_GATE_TEMPLATE = 'ATGAAGTCTACTGTCGCCGTCGATCTGAACGCCGTGAAGATCCAGAACCTGATGGGCTGCATGCTGGTGATCGACCTGGGCTGC';

// Golden Gate primer structure: [flanking]-[BsaI site]-[spacer]-[overhang]-[homology]
// Example: AACAGGTCTC A GGAG ATGAAGTCTACTGTCGCCGT
//          flanking(6) + BsaI(7) + spacer(1) + overhang(4) + homology(20)
const GOLDEN_GATE_FWD_PRIMER = 'AACAGGTCTCAGGAGATGAAGTCTACTGTCGCCGT';
const EXPECTED_ANNEALING_REGION = 'ATGAAGTCTACTGTCGCCGT'; // 20bp that matches template

// NEBuilder/Gibson primer: [homology tail]-[annealing region]
const NEBUILDER_FWD_PRIMER = 'GCTAGCTAGCTAGCATGAAGTCTACTGTCGCC';
const NEBUILDER_ANNEALING = 'ATGAAGTCTACTGTCGCC'; // 18bp that matches template

describe('findAnnealingRegion', () => {
  it('should find the annealing region for a Golden Gate primer', () => {
    const result = findAnnealingRegion(GOLDEN_GATE_FWD_PRIMER, GOLDEN_GATE_TEMPLATE);

    expect(result).not.toBeNull();
    expect(result.sequence).toBe(EXPECTED_ANNEALING_REGION);
    expect(result.length).toBe(20);
    expect(result.position).toBe(0); // Starts at position 0 of template
    expect(result.strand).toBe('forward');
    expect(result.tailLength).toBe(GOLDEN_GATE_FWD_PRIMER.length - 20);
    expect(result.mismatches).toBe(0);
    expect(result.identity).toBe(1.0);
  });

  it('should find annealing region for NEBuilder primer', () => {
    const result = findAnnealingRegion(NEBUILDER_FWD_PRIMER, GOLDEN_GATE_TEMPLATE);

    expect(result).not.toBeNull();
    expect(result.sequence).toBe(NEBUILDER_ANNEALING);
    expect(result.length).toBe(18);
    expect(result.tailLength).toBeGreaterThan(0);
  });

  it('should calculate Tm for the annealing region', () => {
    const result = findAnnealingRegion(GOLDEN_GATE_FWD_PRIMER, GOLDEN_GATE_TEMPLATE);

    expect(result.tm).toBeDefined();
    expect(result.tm).toBeGreaterThan(50);
    expect(result.tm).toBeLessThan(75);
    expect(result.gc).toBeDefined();
    expect(result.gcPercent).toBeDefined();
  });

  it('should return null when primer does not match template', () => {
    const nonMatchingPrimer = 'AAAAAAAAAAAAAAAAAAAA';
    const result = findAnnealingRegion(nonMatchingPrimer, GOLDEN_GATE_TEMPLATE);

    expect(result).toBeNull();
  });

  it('should find reverse primer annealing on reverse complement', () => {
    // Create a reverse primer that matches the 3' end of template (reverse complement)
    const revTemplate = 'GCAGCCCAGGTCGATCACCAGCATGCAGCCCATCAGGTTCTGGATCTTCACGGCGTTCAGATCGACGGCGACAGTAGACTTCAT';
    const revPrimer = 'AACAGGTCTCATACTGCAGCCCAGGTCGATCACC'; // tail + homology to end

    const result = findAnnealingRegion(revPrimer, GOLDEN_GATE_TEMPLATE);

    // The primer matches the reverse complement of template
    expect(result).not.toBeNull();
    expect(result.strand).toBe('reverse');
  });

  it('should support fuzzy matching with mismatches', () => {
    // Primer with 1 mismatch in annealing region
    const primerWithMismatch = 'AACAGGTCTCAGGAGATGAAGTCTACTGTCGCTGT'; // GCC -> GCT

    const result = findAnnealingRegion(primerWithMismatch, GOLDEN_GATE_TEMPLATE, {
      allowMismatch: true,
      maxMismatches: 2,
    });

    expect(result).not.toBeNull();
    expect(result.mismatches).toBeGreaterThan(0);
    expect(result.identity).toBeLessThan(1.0);
  });
});

describe('analyzeAssemblyPrimer', () => {
  it('should provide comprehensive assembly primer analysis', () => {
    const analysis = analyzeAssemblyPrimer(GOLDEN_GATE_FWD_PRIMER, GOLDEN_GATE_TEMPLATE);

    // Check structure
    expect(analysis.fullPrimer).toBeDefined();
    expect(analysis.annealingRegion).toBeDefined();
    expect(analysis.templateBinding).toBeDefined();
    expect(analysis.scores).toBeDefined();
    expect(analysis.warnings).toBeDefined();

    // Template binding should be found
    expect(analysis.templateBinding.found).toBe(true);
    expect(analysis.isAssemblyPrimer).toBe(true);

    // Full primer properties
    expect(analysis.fullPrimer.sequence).toBe(GOLDEN_GATE_FWD_PRIMER);
    expect(analysis.fullPrimer.length).toBe(GOLDEN_GATE_FWD_PRIMER.length);

    // Annealing region properties
    expect(analysis.annealingRegion.sequence).toBe(EXPECTED_ANNEALING_REGION);
    expect(analysis.annealingRegion.tm).toBeDefined();
  });

  it('should score based on annealing region properties', () => {
    const analysis = analyzeAssemblyPrimer(GOLDEN_GATE_FWD_PRIMER, GOLDEN_GATE_TEMPLATE);

    // Scores should be present
    expect(analysis.scores.annealingTm).toBeDefined();
    expect(analysis.scores.annealingGc).toBeDefined();
    expect(analysis.scores.annealingLength).toBeDefined();
    expect(analysis.scores.fullPrimerHairpin).toBeDefined();
    expect(analysis.scores.fullPrimerHomodimer).toBeDefined();

    // Scores should be between 0 and 1
    expect(analysis.scores.annealingTm).toBeGreaterThanOrEqual(0);
    expect(analysis.scores.annealingTm).toBeLessThanOrEqual(1);
  });

  it('should handle non-matching primers gracefully', () => {
    const nonMatchingPrimer = 'AAAAAAAAAAAAAAAAAAAA';
    const analysis = analyzeAssemblyPrimer(nonMatchingPrimer, GOLDEN_GATE_TEMPLATE);

    expect(analysis.templateBinding.found).toBe(false);
    expect(analysis.annealingRegion).toBeNull();
    expect(analysis.warnings.length).toBeGreaterThan(0);
    expect(analysis.warnings[0].type).toBe('templateMismatch');
  });
});

describe('analyzeSinglePrimer with assembly mode', () => {
  it('should detect assembly primer and use annealing region for scoring', () => {
    const result = analyzeSinglePrimer(GOLDEN_GATE_FWD_PRIMER, {
      mode: 'assembly',
      template: GOLDEN_GATE_TEMPLATE,
    });

    // Should have assembly-specific fields
    expect(result.isAssemblyPrimer).toBe(true);
    expect(result.annealingRegion).toBeDefined();
    expect(result.tailRegion).toBeDefined();
    expect(result.templateBinding).toBeDefined();

    // Tm should be the annealing Tm, not full primer Tm
    expect(result.tm).toBeLessThan(result.fullPrimerTm);
    expect(result.annealingRegion.tm).toBe(result.tm);
  });

  it('should still check full primer for secondary structure', () => {
    const result = analyzeSinglePrimer(GOLDEN_GATE_FWD_PRIMER, {
      mode: 'assembly',
      template: GOLDEN_GATE_TEMPLATE,
    });

    // Secondary structure should be checked on full primer
    expect(result.thermodynamics.hairpinDG).toBeDefined();
    expect(result.thermodynamics.homodimerDG).toBeDefined();
    expect(result.scores.hairpin).toBeDefined();
    expect(result.scores.homodimer).toBeDefined();
  });

  it('should fall back to standard analysis without template', () => {
    const result = analyzeSinglePrimer(GOLDEN_GATE_FWD_PRIMER, {
      mode: 'assembly',
      // No template provided
    });

    // Should not have assembly-specific fields
    expect(result.isAssemblyPrimer).toBeUndefined();
    expect(result.annealingRegion).toBeUndefined();
  });

  it('should produce different scores for assembly vs amplification mode', () => {
    const assemblyResult = analyzeSinglePrimer(GOLDEN_GATE_FWD_PRIMER, {
      mode: 'assembly',
      template: GOLDEN_GATE_TEMPLATE,
    });

    const amplificationResult = analyzeSinglePrimer(GOLDEN_GATE_FWD_PRIMER, {
      mode: 'amplification',
    });

    // Assembly mode should score the annealing region (20bp, ~60°C Tm)
    // Amplification mode scores full primer (35bp, higher Tm)
    // This should result in different Tm scores
    expect(assemblyResult.tm).not.toBe(amplificationResult.tm);

    // Assembly scoring should be more favorable since it's designed for assembly
    // (the annealing region is well-designed even if full primer looks unusual)
    expect(assemblyResult.scores.length).toBeGreaterThan(amplificationResult.scores.length);
  });
});

describe('analyzePrimers with assembly mode', () => {
  it('should analyze primer pair with template context', () => {
    // Create a simple reverse primer
    const revPrimerSeq = 'AACAGGTCTCATACTGCAGCCCAGGTCGATCACC';

    const result = analyzePrimers(
      { seq: GOLDEN_GATE_FWD_PRIMER },
      { seq: revPrimerSeq },
      { mode: 'assembly', template: GOLDEN_GATE_TEMPLATE }
    );

    expect(result.forward).toBeDefined();
    expect(result.reverse).toBeDefined();
    expect(result.pair).toBeDefined();

    // Forward should have assembly analysis if it matched
    if (result.forward.isAssemblyPrimer) {
      expect(result.forward.annealingRegion).toBeDefined();
    }
  });
});

describe('Edge cases and robustness', () => {
  it('should handle very short templates', () => {
    const shortTemplate = 'ATGAAGTCTACTGTCGCC';
    const result = findAnnealingRegion(GOLDEN_GATE_FWD_PRIMER, shortTemplate);

    // May or may not find a match depending on minimum length settings
    // Should not throw an error
    expect(() => findAnnealingRegion(GOLDEN_GATE_FWD_PRIMER, shortTemplate)).not.toThrow();
  });

  it('should handle primers that match multiple locations', () => {
    // Template with repeated sequence
    const repeatTemplate = 'ATGAAGTCTACTGTCGCCGTGATCTGAACATGAAGTCTACTGTCGCCGT';

    const result = findAnnealingRegion(GOLDEN_GATE_FWD_PRIMER, repeatTemplate);

    // Should find a match (algorithm searches from longest to shortest)
    expect(result).not.toBeNull();
    // The algorithm will find the longest annealing region, which could be at either occurrence
    expect(result.position).toBeGreaterThanOrEqual(0);
    expect(result.sequence.length).toBeGreaterThanOrEqual(15);
  });

  it('should preserve case insensitivity', () => {
    const lowercasePrimer = GOLDEN_GATE_FWD_PRIMER.toLowerCase();
    const lowercaseTemplate = GOLDEN_GATE_TEMPLATE.toLowerCase();

    const result = findAnnealingRegion(lowercasePrimer, lowercaseTemplate);

    expect(result).not.toBeNull();
    expect(result.sequence).toBe(EXPECTED_ANNEALING_REGION); // Uppercase output
  });
});

// =============================================================================
// Golden Gate Enzyme Site Detection
// =============================================================================

describe('detectGoldenGateSites', () => {
  it('should detect BsaI recognition site', () => {
    const primerWithBsaI = 'AACAGGTCTCAGGAGATGAAGTCTACT'; // Contains GGTCTC
    const result = detectGoldenGateSites(primerWithBsaI);

    expect(result.hasGoldenGateSite).toBe(true);
    expect(result.primaryEnzyme).toBe('BsaI');
    expect(result.detectedSites).toHaveLength(1);
    expect(result.detectedSites[0].enzyme).toBe('BsaI');
    expect(result.detectedSites[0].recognition).toBe('GGTCTC');
    expect(result.isLikelyGoldenGatePrimer).toBe(true);
  });

  it('should detect BsmBI recognition site', () => {
    const primerWithBsmBI = 'AACACGTCTCAGGAGATGAAGTCTACT'; // Contains CGTCTC
    const result = detectGoldenGateSites(primerWithBsmBI);

    expect(result.hasGoldenGateSite).toBe(true);
    expect(result.primaryEnzyme).toBe('BsmBI');
    expect(result.detectedSites[0].recognition).toBe('CGTCTC');
  });

  it('should detect BbsI recognition site', () => {
    const primerWithBbsI = 'AACAGAAGACNNGGAGATGAAGTCTACT'; // Contains GAAGAC
    const result = detectGoldenGateSites(primerWithBbsI);

    expect(result.hasGoldenGateSite).toBe(true);
    expect(result.primaryEnzyme).toBe('BbsI');
  });

  it('should detect SapI recognition site', () => {
    const primerWithSapI = 'AACAGCTCTTCNGGAGATGAAGTCTACT'; // Contains GCTCTTC
    const result = detectGoldenGateSites(primerWithSapI);

    expect(result.hasGoldenGateSite).toBe(true);
    expect(result.primaryEnzyme).toBe('SapI');
  });

  it('should return no detection for regular PCR primer', () => {
    const regularPrimer = 'ATGAAGTCTACTGTCGCCGT'; // No GG enzyme sites
    const result = detectGoldenGateSites(regularPrimer);

    expect(result.hasGoldenGateSite).toBe(false);
    expect(result.primaryEnzyme).toBeNull();
    expect(result.detectedSites).toHaveLength(0);
    expect(result.recommendation).toBeNull();
  });

  it('should provide recommendation when GG site detected', () => {
    const primerWithBsaI = 'AACAGGTCTCAGGAGATGAAGTCTACT';
    const result = detectGoldenGateSites(primerWithBsaI);

    expect(result.recommendation).toContain('template');
    expect(result.recommendation).toContain('Golden Gate');
  });

  it('should handle case insensitivity', () => {
    const lowercasePrimer = 'aacaggtctcaggagatgaagtctact';
    const result = detectGoldenGateSites(lowercasePrimer);

    expect(result.hasGoldenGateSite).toBe(true);
    expect(result.primaryEnzyme).toBe('BsaI');
  });

  it('should detect reverse complement sites', () => {
    // GGTCTC reverse complement is GAGACC
    const primerWithRCSite = 'AACAGAGACCAGGAGATGAAGTCTACT';
    const result = detectGoldenGateSites(primerWithRCSite);

    expect(result.hasGoldenGateSite).toBe(true);
    expect(result.detectedSites.some(s => s.orientation === 'reverse')).toBe(true);
  });
});

describe('analyzeSinglePrimer with GG detection', () => {
  it('should include GG detection info when enzyme site found without template', () => {
    const result = analyzeSinglePrimer(GOLDEN_GATE_FWD_PRIMER, {
      mode: 'amplification',
      // No template
    });

    expect(result.goldenGateDetection).not.toBeNull();
    expect(result.goldenGateDetection.hasGoldenGateSite).toBe(true);
    expect(result.goldenGateDetection.primaryEnzyme).toBe('BsaI');
  });

  it('should include warning when GG site detected without template', () => {
    const result = analyzeSinglePrimer(GOLDEN_GATE_FWD_PRIMER, {
      mode: 'amplification',
    });

    const ggWarning = result.warnings.find(w => w.type === 'goldenGateDetected');
    expect(ggWarning).toBeDefined();
    expect(ggWarning.severity).toBe('info');
    expect(ggWarning.enzyme).toBe('BsaI');
  });

  it('should not include GG warning when template is provided', () => {
    const result = analyzeSinglePrimer(GOLDEN_GATE_FWD_PRIMER, {
      mode: 'assembly',
      template: GOLDEN_GATE_TEMPLATE,
    });

    const ggWarning = result.warnings.find(w => w.type === 'goldenGateDetected');
    expect(ggWarning).toBeUndefined();
  });

  it('should not have GG detection for regular primers', () => {
    const regularPrimer = 'ATGAAGTCTACTGTCGCCGT';
    const result = analyzeSinglePrimer(regularPrimer, {
      mode: 'amplification',
    });

    expect(result.goldenGateDetection).toBeNull();
  });
});
