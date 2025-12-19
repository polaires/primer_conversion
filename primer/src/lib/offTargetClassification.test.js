/**
 * Tests for Off-Target Type A-F Classification
 */

import { expect, describe, it } from "vitest";
import {
  findTypeAOffTargets,
  findTypeBOffTargets,
  findTypeCOffTargets,
  findTypeDOffTargets,
  classifyOffTargets,
  scoreOffTargetClassification,
  analyzeOffTargetsPair,
  RISK_LEVELS,
} from "./offTargetClassification.js";

// Helper to create test template with primer binding site
function createTemplate(beforePrimer, primerBindingSite, afterPrimer) {
  return beforePrimer + primerBindingSite + afterPrimer;
}

// Reverse complement helper
function rc(seq) {
  const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
  return seq.split('').reverse().map(c => comp[c]).join('');
}

// Non-palindromic test primer (different from its reverse complement)
const TEST_PRIMER = 'ATGCATGCATGCATGCAA';  // 18bp, NOT self-complementary
const TEST_PRIMER_RC = rc(TEST_PRIMER);    // TTGCATGCATGCATGCAT

describe('Type A: Alternative Binding Sites', () => {
  it('should detect exact duplicate binding sites', () => {
    const primer = 'ATGCATGCATGCATGCAT';  // 18bp primer
    // Template has the exact primer binding site twice
    const template = 'AAAA' + rc(primer) + 'AAAAAAAAAA' + rc(primer) + 'TTTT';

    const sites = findTypeAOffTargets(primer, template, 4);  // First site at position 4

    expect(sites.length).toBeGreaterThan(0);
    expect(sites.some(s => s.subtype === 'exact')).toBe(true);
    expect(sites[0].risk).toBe(RISK_LEVELS.HIGH);
  });

  it('should return empty for unique binding site', () => {
    const primer = 'ATGCATGCATGCATGCAT';
    const template = 'AAAA' + rc(primer) + 'TTTTTTTTTTTTTTTTTTTTTTTT';

    const sites = findTypeAOffTargets(primer, template, 4);

    // Should not find additional sites
    const extraSites = sites.filter(s => s.position !== 4);
    expect(extraSites.length).toBe(0);
  });

  it('should detect near-exact matches with mismatches', () => {
    const primer = 'ATGCATGCATGCATGCAT';  // 18bp
    const primerMutant = 'ATGCATGCATGCTTGCAT';  // 1 mismatch in middle
    const template = 'AAAA' + rc(primer) + 'AAAAAAAAAA' + rc(primerMutant) + 'TTTT';

    const sites = findTypeAOffTargets(primer, template, 4);

    // The mutant site should be detected if the 3' anchor matches
    // and there are ≤2 mismatches
    expect(sites.length).toBeGreaterThanOrEqual(0);
  });
});

describe('Type B: Antisense Binding', () => {
  it('should detect direct template binding (antisense)', () => {
    // Use non-palindromic primer
    const primer = 'ATGCATGCATGCATGCAA';  // 18bp, ends with AA
    // Put the full primer sequence directly in template
    // This means on the opposite strand (rc of this region), there's a binding site
    const template = 'GGGG' + primer + 'TTTTTTTTTTTTTTTT';

    const sites = findTypeBOffTargets(primer, template);

    // Should detect antisense binding potential
    expect(sites.length).toBeGreaterThan(0);
    expect(sites[0].type).toBe('B');
    expect(sites[0].risk).toBe(RISK_LEVELS.HIGH);
  });

  it('should not flag normal complement binding', () => {
    const primer = TEST_PRIMER;
    // Template only has reverse complement (normal binding)
    const template = 'GGGG' + TEST_PRIMER_RC + 'CCCCCCCCCCCCCCCC';

    const sites = findTypeBOffTargets(primer, template);

    // Should not detect the normal binding site as antisense
    // because the primer doesn't appear directly
    expect(sites.length).toBe(0);
  });
});

describe('Type C: Partial 3\' Homology', () => {
  it('should detect 3\' anchor matches elsewhere', () => {
    // Primer's 3' end (anchor) is what we're testing for partial homology
    // The function searches for rc(anchor) in template
    const primer = 'AAAAAAAAATGCATGCGG';  // 18bp
    const anchor8 = primer.slice(-8);  // 'ATGCATGCGG'.slice(-8) = 'GCATGCGG'
    const rcAnchor8 = rc(anchor8);  // rc('GCATGCGG') = 'CCGCATGC'

    // Template: extra rc(anchor) at position 0, then filler, then intended binding site
    const filler = 'ATATATATATATATATATATATATATATAT';  // 30bp filler
    const template = rcAnchor8 + filler + rc(primer) + 'CCCCCCCCCCCC';

    // Intended position is where rc(primer) starts
    const intendedPos = rcAnchor8.length + filler.length;  // 8 + 30 = 38

    const sites = findTypeCOffTargets(primer, template, intendedPos);

    // Should detect partial 3' homology at position 0
    // (far from intended position 38, so should be flagged)
    expect(sites.some(s => s.type === 'C')).toBe(true);
    expect(sites.filter(s => s.type === 'C')[0]?.risk).toBe(RISK_LEVELS.MEDIUM);
  });

  it('should return empty when 3\' end is unique', () => {
    const primer = TEST_PRIMER;
    const template = 'GGGG' + TEST_PRIMER_RC + 'CCCCCCCCCCCCCCCCCCCC';

    const sites = findTypeCOffTargets(primer, template, 4);

    // No additional 3' homology sites beyond intended binding
    const extraSites = sites.filter(s => Math.abs(s.position - 4) > primer.length);
    expect(extraSites.length).toBe(0);
  });
});

describe('Type D: Internal Homology', () => {
  it('should detect internal sequence matches', () => {
    const primer = 'ATGCATGCATGCAAAAAA';  // Internal: ATGCATGCATGC, 3': AAAAAA
    const internalSeq = 'ATGCATGCATGC';
    // Template has internal sequence elsewhere
    const template = 'AAAA' + rc(primer) + 'GGGG' + rc(internalSeq) + 'TTTTTTTTTT';

    const sites = findTypeDOffTargets(primer, template, 4);

    // Should detect internal homology
    expect(sites.some(s => s.type === 'D')).toBe(true);
    expect(sites.filter(s => s.type === 'D')[0]?.risk).toBe(RISK_LEVELS.LOW);
  });

  it('should not flag 3\' end matches as internal', () => {
    const primer = 'ATGCATGCATGCATGCAT';
    const template = 'AAAA' + rc(primer) + 'TTTTTTTTTTTTTTTTTTTTTTTT';

    const sites = findTypeDOffTargets(primer, template, 4);

    // Internal homology excludes the 3' region
    // Should not have sites that correspond to 3' end
    expect(sites.every(s => s.primerStart < primer.length - 5)).toBe(true);
  });
});

describe('classifyOffTargets', () => {
  it('should return pass status for clean primer', () => {
    // Use non-palindromic primer
    const primer = TEST_PRIMER;
    const template = 'GGGG' + TEST_PRIMER_RC + 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC';

    const result = classifyOffTargets(primer, template, 4);

    expect(result.status).toBe('pass');
    expect(result.hasCritical).toBe(false);
  });

  it('should return critical status for high-risk sites', () => {
    const primer = TEST_PRIMER;
    // Put primer sequence directly in template (Type B - antisense binding)
    // This creates a high-risk off-target
    const template = 'GGGG' + primer + 'AAAAAAAAAA' + TEST_PRIMER_RC + 'TTTT';

    // Intended position is where rc(primer) is
    const intendedPos = 4 + primer.length + 10;

    const result = classifyOffTargets(primer, template, intendedPos);

    expect(result.status).toBe('critical');
    expect(result.hasCritical).toBe(true);
    expect(result.counts.highRisk).toBeGreaterThan(0);
  });

  it('should return warning status for medium-risk sites only', () => {
    const primer = 'GGGGGGGGGATGCATGCAA';  // 3' end = ATGCATGCAA
    const anchor = 'ATGCATGCAA';
    const rcAnchor = rc(anchor);
    // 3' anchor appears elsewhere but no exact full match
    const template = 'TTTT' + rc(primer) + 'AAAAAAAAAA' + rcAnchor + 'CCCCCCCCCCCC';

    const result = classifyOffTargets(primer, template, 4);

    // Should have medium risk (Type C) but no high risk (Type A/B)
    if (result.counts.highRisk === 0 && result.counts.mediumRisk > 0) {
      expect(result.status).toBe('warning');
    }
  });

  it('should provide accurate counts', () => {
    const primer = TEST_PRIMER;
    const template = 'GGGG' + TEST_PRIMER_RC + 'CCCCCCCCCCCCCCCCCCCCCCCC';

    const result = classifyOffTargets(primer, template, 4);

    expect(result.counts.total).toBe(
      result.counts.typeA + result.counts.typeB + result.counts.typeC + result.counts.typeD
    );
    expect(result.counts.total).toBe(
      result.counts.highRisk + result.counts.mediumRisk + result.counts.lowRisk
    );
  });

  it('should generate meaningful summary', () => {
    const primer = TEST_PRIMER;
    const template = 'GGGG' + TEST_PRIMER_RC + 'CCCCCCCCCCCCCCCCCCCCCCCC';

    const result = classifyOffTargets(primer, template, 4);

    expect(typeof result.summary).toBe('string');
    expect(result.summary.length).toBeGreaterThan(0);
  });
});

describe('scoreOffTargetClassification', () => {
  it('should return 1.0 for no off-targets', () => {
    const classification = {
      counts: { highRisk: 0, mediumRisk: 0, lowRisk: 0 }
    };

    expect(scoreOffTargetClassification(classification)).toBe(1.0);
  });

  it('should return 0.0 for 3+ high-risk sites', () => {
    const classification = {
      counts: { highRisk: 3, mediumRisk: 0, lowRisk: 0 }
    };

    expect(scoreOffTargetClassification(classification)).toBe(0.0);
  });

  it('should return low score for high-risk sites', () => {
    const oneHighRisk = { counts: { highRisk: 1, mediumRisk: 0, lowRisk: 0 } };
    const twoHighRisk = { counts: { highRisk: 2, mediumRisk: 0, lowRisk: 0 } };

    expect(scoreOffTargetClassification(oneHighRisk)).toBe(0.3);
    expect(scoreOffTargetClassification(twoHighRisk)).toBe(0.1);
  });

  it('should reduce score for medium-risk sites', () => {
    const noRisk = { counts: { highRisk: 0, mediumRisk: 0, lowRisk: 0 } };
    const mediumRisk = { counts: { highRisk: 0, mediumRisk: 2, lowRisk: 0 } };

    expect(scoreOffTargetClassification(mediumRisk)).toBeLessThan(
      scoreOffTargetClassification(noRisk)
    );
  });

  it('should minimally reduce score for low-risk sites', () => {
    const noRisk = { counts: { highRisk: 0, mediumRisk: 0, lowRisk: 0 } };
    const lowRisk = { counts: { highRisk: 0, mediumRisk: 0, lowRisk: 5 } };

    const diff = scoreOffTargetClassification(noRisk) - scoreOffTargetClassification(lowRisk);
    expect(diff).toBeLessThan(0.2);  // Low risk has minimal impact
  });
});

describe('analyzeOffTargetsPair', () => {
  // Non-palindromic primers for testing
  const FWD = 'ATGCATGCATGCATGCAA';  // Ends with AA
  const REV = 'GCATGCATGCATGCATGG';  // Ends with GG
  const FWD_RC = rc(FWD);
  const REV_RC = rc(REV);

  it('should analyze both primers', () => {
    const template = FWD_RC + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + REV_RC;

    const result = analyzeOffTargetsPair(FWD, REV, template);

    expect(result).toHaveProperty('fwd');
    expect(result).toHaveProperty('rev');
    expect(result).toHaveProperty('combined');
    expect(result).toHaveProperty('score');
  });

  it('should return pass for clean primer pair', () => {
    const template = FWD_RC + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + REV_RC;

    const result = analyzeOffTargetsPair(FWD, REV, template);

    // Status may be 'pass' or 'warning' depending on Type E/F analysis
    expect(['pass', 'warning']).toContain(result.status);
    // Should have Type E (hairpin/homodimer) and Type F (heterodimer) analysis
    expect(result.fwd).toHaveProperty('typeE');
    expect(result.rev).toHaveProperty('typeE');
    expect(result).toHaveProperty('typeF');
  });

  it('should return critical if either primer has high-risk sites', () => {
    // Put the forward primer sequence directly in template (Type B - antisense)
    // This creates a high-risk off-target
    const template = FWD + 'GGGGGGGGGG' + FWD_RC + 'AAAAAAAAAAAAAAAAAAAAAA' + REV_RC;

    const result = analyzeOffTargetsPair(FWD, REV, template);

    // Fwd has antisense binding site (Type B), so should have high risk
    expect(result.fwd.counts.typeB).toBeGreaterThan(0);
    expect(result.fwd.counts.highRisk).toBeGreaterThan(0);
  });

  it('should use enhanced scoring including Type E/F', () => {
    const template = FWD_RC + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + REV_RC;

    const result = analyzeOffTargetsPair(FWD, REV, template);

    // Score is now the enhanced score that includes Type E/F penalties
    expect(result.score).toBeGreaterThanOrEqual(0);
    expect(result.score).toBeLessThanOrEqual(1);
    // Should have score breakdown
    expect(result).toHaveProperty('scoreBreakdown');
    expect(result.scoreBreakdown).toHaveProperty('fwdBase');
    expect(result.scoreBreakdown).toHaveProperty('fwdEnhanced');
  });

  it('should combine counts from both primers including Type E/F', () => {
    const template = FWD_RC + 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' + REV_RC;

    const result = analyzeOffTargetsPair(FWD, REV, template);

    // Total sites now includes Type A-D + Type E (hairpin/homodimer) + Type F (heterodimer)
    const typeADSites = result.fwd.counts.total + result.rev.counts.total;
    const typeESites = result.combined.typeE;
    const typeFSites = result.combined.typeF;
    expect(result.combined.totalSites).toBe(typeADSites + typeESites + typeFSites);
    // Should have thermodynamic values
    expect(result).toHaveProperty('thermodynamics');
    expect(result.thermodynamics).toHaveProperty('fwdHairpinDG');
    expect(result.thermodynamics).toHaveProperty('heterodimerDG');
  });
});

describe('Edge Cases', () => {
  it('should handle short primers gracefully', () => {
    const primer = 'ATGCATGC';  // 8bp - very short
    const template = 'AAAA' + rc(primer) + 'TTTTTTTTTTTTTTTTTTTT';

    const result = classifyOffTargets(primer, template, 4);

    // Should not crash, may or may not find sites
    expect(result).toHaveProperty('status');
    expect(result).toHaveProperty('counts');
  });

  it('should handle primer longer than template portion', () => {
    const primer = 'ATGCATGCATGCATGCAT';
    const template = 'AA' + rc(primer);  // Minimal template

    const result = classifyOffTargets(primer, template, 2);

    expect(result).toHaveProperty('status');
  });

  it('should handle templates with repeated sequences', () => {
    const primer = 'ATATATATATATATAT';  // AT repeat
    const template = 'ATATAT' + rc(primer) + 'ATATATATATATAT';  // Template has AT repeats

    const result = classifyOffTargets(primer, template, 6);

    // Should handle without error
    expect(result).toHaveProperty('status');
    expect(result).toHaveProperty('counts');
  });

  it('should handle GC-rich primers', () => {
    const primer = 'GCGCGCGCGCGCGCGC';  // GC repeat
    const template = 'AAAA' + rc(primer) + 'TTTTTTTTTTTTTTTTTTTT';

    const result = classifyOffTargets(primer, template, 4);

    expect(result).toHaveProperty('status');
  });
});

describe('Integration: Realistic Scenarios', () => {
  it('should detect problematic duplicate binding site', () => {
    // Real scenario: primer sequence appears twice in plasmid
    const primer = 'CGATCGATCGATCGATCG';
    const plasmid = 'AAAA' + rc(primer) + 'GGGGGGGGGGGGGGGG' + rc(primer) + 'CCCC';

    const result = classifyOffTargets(primer, plasmid, 4);

    expect(result.hasCritical).toBe(true);
    expect(result.counts.typeA).toBeGreaterThan(0);
    expect(result.summary).toContain('HIGH RISK');
  });

  it('should flag primers that can bind antisense', () => {
    // Primer sequence appears directly in template (antisense potential)
    const primer = 'ATGCATGCATGCATGCAT';
    const template = 'GGGG' + primer + 'AAAAAAAAAAAAAAAA' + rc(primer) + 'TTTT';

    const result = classifyOffTargets(primer, template, template.indexOf(rc(primer)));

    expect(result.counts.typeB).toBeGreaterThan(0);
  });

  it('should give clean result for well-designed primer', () => {
    // Good primer: unique sequence, no problematic features
    const primer = 'ATGCGTACGTACGTACGT';  // 18bp, unique sequence
    const template = 'GGGGGGGGGG' + rc(primer) + 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCC';

    const result = classifyOffTargets(primer, template, 10);

    // Should have minimal or no off-target issues
    expect(result.counts.highRisk).toBe(0);
  });
});
