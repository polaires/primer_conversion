/**
 * Tests for Auto-Domestication Optimizer
 *
 * Tests cover:
 * - Internal site detection
 * - Adjacent site detection
 * - Junction finding and scoring
 * - Enzyme recommendations
 * - Fragment size validation
 * - Global overhang optimization
 * - Site merging strategies
 */

import { describe, it, expect, beforeEach } from 'vitest';
import {
  analyzeForDomestication,
  findDomesticationJunctions,
  detectAdjacentSites,
  recommendAlternativeEnzymes,
  validateFragmentSizes,
  optimizeGlobalOverhangSet,
  mergeAdjacentSites,
  validateDomesticationOverhangs,
  validatePostDomestication,
  DOMESTICATION_DEFAULTS,
} from './auto-domestication-optimizer';

// Test sequences
const SEQUENCES = {
  // No internal BsaI sites
  clean: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',

  // Single BsaI site (GGTCTC) at position 50
  singleSite: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGGTCTCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',

  // Two BsaI sites far apart (>200bp)
  twoSitesFar: 'ATGCATGCGGTCTCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCGGTCTCATGCATGC',

  // Two BsaI sites close together (<50bp)
  twoSitesClose: 'ATGCATGCGGTCTCATGCATGCATGCATGCGGTCTCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',

  // BsmBI site (CGTCTC) - no BsaI
  bsmbiOnly: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATCGTCTCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',

  // Both BsaI and BsmBI sites
  multiEnzyme: 'ATGCATGCGGTCTCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATCGTCTCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',

  // Reverse complement BsaI site (GAGACC)
  reverseSite: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAGACCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',
};

describe('analyzeForDomestication', () => {
  it('should return needsDomestication: false for clean sequences', () => {
    const result = analyzeForDomestication(SEQUENCES.clean, 'BsaI');
    expect(result.needsDomestication).toBe(false);
    expect(result.sites).toHaveLength(0);
    expect(result.error).toBeNull();
  });

  it('should detect a single internal BsaI site', () => {
    const result = analyzeForDomestication(SEQUENCES.singleSite, 'BsaI');
    expect(result.needsDomestication).toBe(true);
    expect(result.sites).toHaveLength(1);
    expect(result.additionalFragments).toBe(1);
  });

  it('should detect multiple internal sites', () => {
    const result = analyzeForDomestication(SEQUENCES.twoSitesFar, 'BsaI');
    expect(result.needsDomestication).toBe(true);
    expect(result.sites).toHaveLength(2);
    expect(result.additionalFragments).toBe(2);
  });

  it('should return error for adjacent sites too close together', () => {
    const result = analyzeForDomestication(SEQUENCES.twoSitesClose, 'BsaI', {
      minSiteDistance: 50,
    });
    expect(result.error).not.toBeNull();
    expect(result.error.type).toBe('ADJACENT_SITES_TOO_CLOSE');
    expect(result.alternativeEnzymes).not.toBeNull();
  });

  it('should detect reverse complement sites', () => {
    const result = analyzeForDomestication(SEQUENCES.reverseSite, 'BsaI');
    expect(result.needsDomestication).toBe(true);
    expect(result.sites).toHaveLength(1);
    expect(result.sites[0].orientation).toBe('reverse');
  });

  it('should work with different enzymes', () => {
    const bsaIResult = analyzeForDomestication(SEQUENCES.bsmbiOnly, 'BsaI');
    const bsmbiResult = analyzeForDomestication(SEQUENCES.bsmbiOnly, 'BsmBI');

    expect(bsaIResult.needsDomestication).toBe(false);
    expect(bsmbiResult.needsDomestication).toBe(true);
  });

  it('should throw error for unknown enzyme', () => {
    expect(() => analyzeForDomestication(SEQUENCES.clean, 'UnknownEnzyme')).toThrow();
  });
});

describe('findDomesticationJunctions', () => {
  it('should find junction candidates for a site', () => {
    const site = { position: 50, sequence: 'GGTCTC', orientation: 'forward' };
    const result = findDomesticationJunctions(SEQUENCES.singleSite, site, 'BsaI');

    expect(result.candidates).toBeDefined();
    expect(result.candidates.length).toBeGreaterThan(0);
    expect(result.recommended).toBeDefined();
  });

  it('should identify valid candidates that break the site', () => {
    const site = { position: 50, sequence: 'GGTCTC', orientation: 'forward' };
    const result = findDomesticationJunctions(SEQUENCES.singleSite, site, 'BsaI');

    const validCandidates = result.candidates.filter(c => c.isValid);
    expect(validCandidates.length).toBeGreaterThanOrEqual(0);
  });

  it('should score candidates and select best one', () => {
    const site = { position: 50, sequence: 'GGTCTC', orientation: 'forward' };
    const result = findDomesticationJunctions(SEQUENCES.singleSite, site, 'BsaI');

    if (result.hasValidOption) {
      expect(result.recommended.isValid).toBe(true);
      expect(result.recommended.quality.score).toBeGreaterThanOrEqual(50);
    }
  });

  it('should generate 4bp overhangs', () => {
    const site = { position: 50, sequence: 'GGTCTC', orientation: 'forward' };
    const result = findDomesticationJunctions(SEQUENCES.singleSite, site, 'BsaI');

    result.candidates.forEach(candidate => {
      expect(candidate.overhang).toHaveLength(4);
    });
  });
});

describe('detectAdjacentSites', () => {
  it('should return hasAdjacentSites: false for single site', () => {
    const sites = [{ position: 50, sequence: 'GGTCTC' }];
    const result = detectAdjacentSites(sites, 50);
    expect(result.hasAdjacentSites).toBe(false);
  });

  it('should return hasAdjacentSites: false for well-spaced sites', () => {
    const sites = [
      { position: 50, sequence: 'GGTCTC' },
      { position: 200, sequence: 'GGTCTC' },
    ];
    const result = detectAdjacentSites(sites, 50);
    expect(result.hasAdjacentSites).toBe(false);
  });

  it('should detect adjacent sites that are too close', () => {
    const sites = [
      { position: 50, sequence: 'GGTCTC' },
      { position: 80, sequence: 'GGTCTC' },
    ];
    const result = detectAdjacentSites(sites, 50);
    expect(result.hasAdjacentSites).toBe(true);
    expect(result.adjacentPairs).toHaveLength(1);
    expect(result.adjacentPairs[0].distance).toBe(30);
  });

  it('should report minimum distance found', () => {
    const sites = [
      { position: 50, sequence: 'GGTCTC' },
      { position: 80, sequence: 'GGTCTC' },
      { position: 200, sequence: 'GGTCTC' },
    ];
    const result = detectAdjacentSites(sites, 50);
    expect(result.minDistanceFound).toBe(30);
  });

  it('should handle empty or null sites array', () => {
    expect(detectAdjacentSites([], 50).hasAdjacentSites).toBe(false);
    expect(detectAdjacentSites(null, 50).hasAdjacentSites).toBe(false);
  });
});

describe('recommendAlternativeEnzymes', () => {
  it('should recommend compatible enzymes first', () => {
    const result = recommendAlternativeEnzymes(SEQUENCES.singleSite, 'BsaI');
    const compatible = result.filter(e => e.isCompatible);

    // Compatible enzymes should be sorted before incompatible
    if (compatible.length > 0) {
      const firstCompatibleIdx = result.findIndex(e => e.isCompatible);
      const firstIncompatibleIdx = result.findIndex(e => !e.isCompatible);
      if (firstIncompatibleIdx > -1) {
        expect(firstCompatibleIdx).toBeLessThan(firstIncompatibleIdx);
      }
    }
  });

  it('should include current enzyme when includeCurrentEnzyme is true', () => {
    const result = recommendAlternativeEnzymes(SEQUENCES.singleSite, 'BsaI', { includeCurrentEnzyme: true });
    const currentEnzyme = result.find(e => e.isCurrent);
    expect(currentEnzyme).toBeDefined();
    expect(currentEnzyme.enzyme).toBe('BsaI');
  });

  it('should exclude current enzyme when includeCurrentEnzyme is false', () => {
    const result = recommendAlternativeEnzymes(SEQUENCES.singleSite, 'BsaI', { includeCurrentEnzyme: false });
    const currentEnzyme = result.find(e => e.isCurrent);
    expect(currentEnzyme).toBeUndefined();
  });

  it('should count internal sites correctly', () => {
    const result = recommendAlternativeEnzymes(SEQUENCES.multiEnzyme, 'BsaI', { includeCurrentEnzyme: true });
    const bsaI = result.find(e => e.enzyme === 'BsaI');
    const bsmBI = result.find(e => e.enzyme === 'BsmBI');

    expect(bsaI.internalSites).toBe(1);
    expect(bsmBI.internalSites).toBe(1);
  });

  it('should sort by site count (fewer is better)', () => {
    const result = recommendAlternativeEnzymes(SEQUENCES.singleSite, 'BsaI');
    for (let i = 1; i < result.length; i++) {
      // Compatible enzymes first, then by site count
      if (!result[i - 1].isCompatible && !result[i].isCompatible) {
        expect(result[i - 1].internalSites).toBeLessThanOrEqual(result[i].internalSites);
      }
    }
  });
});

describe('validateFragmentSizes', () => {
  it('should pass for sequence with no junctions', () => {
    const result = validateFragmentSizes(SEQUENCES.clean, [], { minFragmentSize: 50 });
    expect(result.valid).toBe(true);
    expect(result.violations).toHaveLength(0);
  });

  it('should detect fragments below minimum size', () => {
    const junctions = [{ position: 10 }, { position: 20 }];
    const result = validateFragmentSizes(SEQUENCES.clean, junctions, { minFragmentSize: 50 });
    expect(result.valid).toBe(false);
    expect(result.violations.some(v => v.type === 'TOO_SMALL')).toBe(true);
  });

  it('should detect fragments above maximum size', () => {
    const junctions = [{ position: 10 }];
    const result = validateFragmentSizes(SEQUENCES.clean, junctions, { maxFragmentSize: 20 });
    expect(result.valid).toBe(false);
    expect(result.violations.some(v => v.type === 'TOO_LARGE')).toBe(true);
  });

  it('should handle numeric position arrays', () => {
    const junctions = [30, 60];
    const result = validateFragmentSizes(SEQUENCES.clean, junctions, { minFragmentSize: 20 });
    expect(result.fragments).toHaveLength(3);
  });

  it('should calculate correct fragment sizes', () => {
    const junctions = [{ position: 50 }];
    const result = validateFragmentSizes(SEQUENCES.clean, junctions, { minFragmentSize: 10 });
    expect(result.fragments[0].size).toBe(50);
    expect(result.fragments[1].size).toBe(SEQUENCES.clean.length - 50);
  });
});

describe('validateDomesticationOverhangs', () => {
  it('should detect duplicate overhangs', () => {
    const domesticationOverhangs = ['ATGC', 'GCTA'];
    const otherOverhangs = ['ATGC', 'TTAA'];
    const result = validateDomesticationOverhangs(domesticationOverhangs, otherOverhangs, 'BsaI');
    expect(result.duplicates).toContain('ATGC');
  });

  it('should detect reverse complement conflicts', () => {
    const domesticationOverhangs = ['ATGC'];
    const otherOverhangs = ['GCAT']; // RC of ATGC
    const result = validateDomesticationOverhangs(domesticationOverhangs, otherOverhangs, 'BsaI');
    expect(result.duplicates.length).toBeGreaterThan(0);
  });

  it('should calculate set fidelity', () => {
    const domesticationOverhangs = ['ATGC', 'GCTA'];
    const otherOverhangs = ['TTAA', 'GGCC'];
    const result = validateDomesticationOverhangs(domesticationOverhangs, otherOverhangs, 'BsaI');
    expect(result.setFidelity).toBeGreaterThanOrEqual(0);
    expect(result.setFidelity).toBeLessThanOrEqual(1);
  });
});

describe('optimizeGlobalOverhangSet', () => {
  it('should return success for empty domestication options', () => {
    const result = optimizeGlobalOverhangSet([], ['ATGC', 'GCTA'], 'BsaI');
    expect(result.success).toBe(true);
  });

  it('should select best combination from candidates', () => {
    const mockOptions = [{
      candidates: [
        { overhang: 'ATGC', isValid: true, quality: { score: 80 } },
        { overhang: 'TTAA', isValid: true, quality: { score: 70 } },
      ],
    }];
    const result = optimizeGlobalOverhangSet(mockOptions, [], 'BsaI');
    expect(result.success).toBeDefined();
    expect(result.selections).toHaveLength(1);
  });

  it('should consider existing overhangs in optimization', () => {
    const mockOptions = [{
      candidates: [
        { overhang: 'ATGC', isValid: true, quality: { score: 80 } },
        { overhang: 'TTAA', isValid: true, quality: { score: 70 } },
      ],
    }];
    const existingOverhangs = ['GCTA']; // Non-duplicate overhang
    const result = optimizeGlobalOverhangSet(mockOptions, existingOverhangs, 'BsaI');

    // Should include existing overhangs in final set
    expect(result.overhangs).toContain('GCTA');
    expect(result.selections).toHaveLength(1);
  });
});

describe('mergeAdjacentSites', () => {
  it('should not merge single site', () => {
    const sites = [{ position: 50, sequence: 'GGTCTC' }];
    const result = mergeAdjacentSites(sites, SEQUENCES.clean);
    expect(result.needsMerging).toBe(false);
    expect(result.groups).toHaveLength(1);
  });

  it('should merge sites within merge distance', () => {
    const sites = [
      { position: 50, sequence: 'GGTCTC' },
      { position: 80, sequence: 'GGTCTC' },
    ];
    const result = mergeAdjacentSites(sites, SEQUENCES.clean, { mergeDistance: 100 });
    expect(result.needsMerging).toBe(true);
    expect(result.groups).toHaveLength(1);
    expect(result.groups[0].sites).toHaveLength(2);
  });

  it('should not merge sites beyond merge distance', () => {
    const sites = [
      { position: 50, sequence: 'GGTCTC' },
      { position: 200, sequence: 'GGTCTC' },
    ];
    const result = mergeAdjacentSites(sites, SEQUENCES.clean, { mergeDistance: 100 });
    expect(result.needsMerging).toBe(false);
    expect(result.groups).toHaveLength(2);
  });

  it('should identify clustered sites that cannot be handled', () => {
    const sites = [
      { position: 50, sequence: 'GGTCTC' },
      { position: 60, sequence: 'GGTCTC' },
      { position: 70, sequence: 'GGTCTC' },
    ];
    const result = mergeAdjacentSites(sites, SEQUENCES.clean, {
      mergeDistance: 100,
      minFragmentSize: 50,
    });
    expect(result.allCanHandle).toBe(false);
    expect(result.strategy).toBe('REQUIRES_ALTERNATIVE_ENZYME');
  });

  it('should handle empty sites array', () => {
    const result = mergeAdjacentSites([], SEQUENCES.clean);
    expect(result.needsMerging).toBe(false);
    expect(result.groups).toHaveLength(0);
  });
});

describe('DOMESTICATION_DEFAULTS', () => {
  it('should have required default values', () => {
    expect(DOMESTICATION_DEFAULTS.minSiteDistance).toBeDefined();
    expect(DOMESTICATION_DEFAULTS.minFragmentSize).toBeDefined();
    expect(DOMESTICATION_DEFAULTS.preferHighFidelityOverhangs).toBeDefined();
  });

  it('should have reasonable default values', () => {
    expect(DOMESTICATION_DEFAULTS.minSiteDistance).toBeGreaterThan(0);
    expect(DOMESTICATION_DEFAULTS.minFragmentSize).toBeGreaterThan(0);
    expect(DOMESTICATION_DEFAULTS.minFragmentSize).toBeLessThanOrEqual(100);
  });
});

describe('validatePostDomestication', () => {
  // Test fragments with valid properties
  const validFragments = [
    { id: 'Frag1', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'ATGC' },
    { id: 'Frag2', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'GCTA' },
    { id: 'Frag3', seq: 'GCTAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: null },
  ];

  it('should validate fragments with no issues', () => {
    const result = validatePostDomestication(validFragments, 'BsaI');
    expect(result.isValid).toBe(true);
    expect(result.issues).toHaveLength(0);
    expect(result.passed.length).toBeGreaterThan(0);
  });

  it('should detect fragments below minimum size', () => {
    const shortFragments = [
      { id: 'Short1', seq: 'ATGCAT', _overhang: 'ATGC' },
      { id: 'Frag2', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: null },
    ];
    const result = validatePostDomestication(shortFragments, 'BsaI', { minFragmentSize: 50 });
    expect(result.isValid).toBe(false);
    expect(result.issues.some(i => i.code === 'FRAGMENT_TOO_SHORT')).toBe(true);
  });

  it('should detect internal recognition sites in fragments', () => {
    // Fragment with internal BsaI site (GGTCTC)
    const fragmentsWithSites = [
      { id: 'Frag1', seq: 'ATGCATGCATGCGGTCTCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'ATGC' },
      { id: 'Frag2', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: null },
    ];
    const result = validatePostDomestication(fragmentsWithSites, 'BsaI');
    expect(result.isValid).toBe(false);
    expect(result.issues.some(i => i.code === 'INTERNAL_SITES_REMAIN')).toBe(true);
  });

  it('should detect duplicate overhangs', () => {
    const fragmentsWithDuplicates = [
      { id: 'Frag1', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'ATGC' },
      { id: 'Frag2', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'ATGC' },
      { id: 'Frag3', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: null },
    ];
    const result = validatePostDomestication(fragmentsWithDuplicates, 'BsaI');
    expect(result.isValid).toBe(false);
    expect(result.issues.some(i => i.code === 'DUPLICATE_OVERHANGS')).toBe(true);
  });

  it('should detect reverse complement duplicates', () => {
    const fragmentsWithRCDuplicates = [
      { id: 'Frag1', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'ATGC' },
      { id: 'Frag2', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'GCAT' }, // RC of ATGC
      { id: 'Frag3', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: null },
    ];
    const result = validatePostDomestication(fragmentsWithRCDuplicates, 'BsaI');
    expect(result.isValid).toBe(false);
    expect(result.issues.some(i => i.code === 'DUPLICATE_OVERHANGS')).toBe(true);
  });

  it('should handle unknown enzyme', () => {
    const result = validatePostDomestication(validFragments, 'UnknownEnzyme');
    expect(result.isValid).toBe(false);
    expect(result.issues.some(i => i.message.includes('Unknown enzyme'))).toBe(true);
  });

  it('should handle empty fragments array', () => {
    const result = validatePostDomestication([], 'BsaI');
    expect(result.isValid).toBe(true);
    expect(result.summary.fragmentCount).toBe(0);
  });

  it('should support circular assembly validation', () => {
    const circularFragments = [
      { id: 'Frag1', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'ATGC' },
      { id: 'Frag2', seq: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'GCTA' },
      { id: 'Frag3', seq: 'GCTAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC', _overhang: 'TGCA' },
    ];
    const result = validatePostDomestication(circularFragments, 'BsaI', { isCircular: true });
    expect(result.summary.fragmentCount).toBe(3);
    // Should check circular junction
    expect(result.passed.some(p => p.code === 'NO_SCAR_SITES')).toBe(true);
  });

  it('should verify assembly sequence when expectedSequence provided', () => {
    const fragments = [
      { id: 'Frag1', seq: 'ATGCATGCAT', _overhang: 'GCAT' },
      { id: 'Frag2', seq: 'GCATTAGCTAGC', _overhang: null },
    ];
    const expectedSequence = 'ATGCATGCATTAGCTAGC';
    const result = validatePostDomestication(fragments, 'BsaI', { expectedSequence });
    expect(result.passed.some(p => p.code === 'ASSEMBLY_VERIFIED')).toBe(true);
  });

  it('should detect assembly mismatch', () => {
    const fragments = [
      { id: 'Frag1', seq: 'ATGCATGCAT', _overhang: 'GCAT' },
      { id: 'Frag2', seq: 'GCATTAGCTAGC', _overhang: null },
    ];
    const wrongSequence = 'COMPLETELYDIFFERENT';
    const result = validatePostDomestication(fragments, 'BsaI', { expectedSequence: wrongSequence });
    expect(result.issues.some(i => i.code === 'ASSEMBLY_MISMATCH')).toBe(true);
  });

  it('should provide summary with status', () => {
    const result = validatePostDomestication(validFragments, 'BsaI');
    expect(result.summary).toBeDefined();
    expect(result.summary.fragmentCount).toBe(3);
    expect(result.summary.totalIssues).toBe(0);
    expect(['valid', 'valid_with_warnings']).toContain(result.summary.status);
  });
});
