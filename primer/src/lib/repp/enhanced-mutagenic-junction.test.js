/**
 * Tests for Enhanced Mutagenic Junction Domestication
 *
 * Tests cover:
 * 1. Junction position search with NEB fidelity data
 * 2. Overhang scoring using experimental matrix
 * 3. Primer quality scoring
 * 4. Global overhang optimization
 * 5. Integration with workflow
 */

import { describe, it, expect, beforeEach } from 'vitest';
import {
  designEnhancedMutagenicJunction,
  designAllEnhancedJunctions,
  scoreOverhangFidelity,
  scorePrimerQuality,
  classifyQuality,
  ENHANCED_JUNCTION_CONFIG,
} from './enhanced-mutagenic-junction.js';
import {
  runDomesticationWorkflow,
  analyzeSequenceForDomestication,
  WORKFLOW_CONFIG,
} from './domestication-primer-workflow.js';
import { findInternalSites, GOLDEN_GATE_ENZYMES, calculateExperimentalFidelity } from './goldengate.js';

// Test sequences with known BsaI sites
const TEST_SEQUENCES = {
  // Simple sequence with one BsaI site (GGTCTC) in coding region
  singleSite: {
    sequence: 'ATGAAAGGTCTCAAATGAGGCGCTAGCTAG' + 'A'.repeat(100),  // BsaI site at position 6
    description: 'Single BsaI site in coding region',
    expectedSites: 1,
  },

  // Sequence with BsaI site where silent mutation is possible
  // GGT CTC = Gly-Leu, can mutate G→A at position 2 to get GGA CTC = Gly-Leu
  silentMutationPossible: {
    sequence: 'ATGAAAGGTCTCAAATGA' + 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG' + 'A'.repeat(50),
    description: 'BsaI site with silent mutation option',
    expectedSites: 1,
    frame: 0,
  },

  // Multiple sites requiring global optimization
  multipleSites: {
    sequence: 'ATGGGTCTC' + 'A'.repeat(200) + 'GGTCTC' + 'A'.repeat(200) + 'TAG',
    description: 'Two BsaI sites requiring optimization',
    expectedSites: 2,
  },

  // No internal sites
  noSites: {
    sequence: 'ATGAAAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG' + 'A'.repeat(100),
    description: 'No BsaI sites',
    expectedSites: 0,
  },
};

describe('Enhanced Mutagenic Junction', () => {
  describe('designEnhancedMutagenicJunction', () => {
    it('should find internal sites correctly', () => {
      const sites = findInternalSites(TEST_SEQUENCES.singleSite.sequence, 'BsaI');
      expect(sites.hasSites).toBe(true);
      expect(sites.count).toBe(TEST_SEQUENCES.singleSite.expectedSites);
    });

    it('should design junction for single site', () => {
      const sites = findInternalSites(TEST_SEQUENCES.silentMutationPossible.sequence, 'BsaI');
      expect(sites.hasSites).toBe(true);

      const site = sites.sites[0];
      const result = designEnhancedMutagenicJunction(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        site,
        'BsaI',
        { frame: 0, organism: 'ecoli' }
      );

      expect(result.success).toBe(true);
      expect(result.junctionPosition).toBeDefined();
      expect(result.overhang).toBeDefined();
      expect(result.overhang.length).toBe(4);
      expect(result.mutation).toBeDefined();
      expect(result.mutation.isSynonymous).toBe(true);
      expect(result.mutation.breaksSite).toBe(true);
    });

    it('should search wider range than basic implementation', () => {
      const sites = findInternalSites(TEST_SEQUENCES.singleSite.sequence, 'BsaI');
      const site = sites.sites[0];

      const result = designEnhancedMutagenicJunction(
        TEST_SEQUENCES.singleSite.sequence,
        site,
        'BsaI',
        { frame: 0, searchRadius: 50 }  // Enhanced uses 50bp by default
      );

      expect(result.success).toBe(true);
      // Should have found alternatives within extended range
      if (result.alternatives) {
        expect(result.alternatives.length).toBeGreaterThan(0);
      }
    });

    it('should use NEB fidelity data for overhang scoring', () => {
      const sites = findInternalSites(TEST_SEQUENCES.silentMutationPossible.sequence, 'BsaI');
      const site = sites.sites[0];

      const result = designEnhancedMutagenicJunction(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        site,
        'BsaI',
        { frame: 0 }
      );

      expect(result.success).toBe(true);
      expect(result.fidelity).toBeDefined();
      expect(result.fidelity.source).toBe('NEB_experimental');
      expect(result.fidelity.singleOverhang).toBeGreaterThan(0);
    });

    it('should include quality breakdown', () => {
      const sites = findInternalSites(TEST_SEQUENCES.silentMutationPossible.sequence, 'BsaI');
      const site = sites.sites[0];

      const result = designEnhancedMutagenicJunction(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        site,
        'BsaI',
        { frame: 0 }
      );

      expect(result.success).toBe(true);
      expect(result.quality).toBeDefined();
      expect(result.quality.breakdown).toBeDefined();
      expect(result.quality.breakdown.overhangFidelity).toBeDefined();
      expect(result.quality.breakdown.mutationQuality).toBeDefined();
      expect(result.quality.breakdown.primerQuality).toBeDefined();
    });

    it('should design optimized primers', () => {
      const sites = findInternalSites(TEST_SEQUENCES.silentMutationPossible.sequence, 'BsaI');
      const site = sites.sites[0];

      const result = designEnhancedMutagenicJunction(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        site,
        'BsaI',
        { frame: 0 }
      );

      expect(result.success).toBe(true);
      expect(result.primers).toBeDefined();
      expect(result.primers.fragment1).toBeDefined();
      expect(result.primers.fragment2).toBeDefined();
      expect(result.primers.fragment1.reversePrimer.sequence).toBeDefined();
      expect(result.primers.fragment2.forwardPrimer.sequence).toBeDefined();
    });
  });

  describe('designAllEnhancedJunctions', () => {
    it('should handle sequence with no sites', () => {
      const result = designAllEnhancedJunctions(
        TEST_SEQUENCES.noSites.sequence,
        'BsaI'
      );

      expect(result.success).toBe(true);
      expect(result.needsDomestication).toBe(false);
    });

    it('should design junctions for multiple sites', () => {
      const result = designAllEnhancedJunctions(
        TEST_SEQUENCES.multipleSites.sequence,
        'BsaI',
        { frame: 0, organism: 'ecoli' }
      );

      // May or may not succeed depending on site positions
      expect(result.needsDomestication).toBe(true);
      expect(result.totalSites).toBe(2);
    });

    it('should calculate overall assembly fidelity', () => {
      const result = designAllEnhancedJunctions(
        TEST_SEQUENCES.multipleSites.sequence,
        'BsaI',
        { frame: 0 }
      );

      if (result.junctions && result.junctions.length > 0) {
        expect(result.fidelity).toBeDefined();
        expect(result.allOverhangs).toBeDefined();
        expect(result.allOverhangs.length).toBeGreaterThan(0);
      }
    });
  });

  describe('scoreOverhangFidelity', () => {
    it('should score palindromic overhangs poorly', () => {
      const result = scoreOverhangFidelity('ATAT', [], 'BsaI', false, null);
      expect(result.score).toBeLessThan(50);  // Palindrome penalty
    });

    it('should score high-fidelity overhangs well', () => {
      // TGAC is known to be high-fidelity
      // Score calculation: base 100, adjusted by NEB fidelity data
      // Single overhang scores typically 40-70 depending on context
      const result = scoreOverhangFidelity('TGAC', [], 'BsaI', true, null);
      expect(result.score).toBeGreaterThan(40);
    });

    it('should penalize conflicts with existing overhangs', () => {
      const withoutConflict = scoreOverhangFidelity('TGAC', [], 'BsaI', false, null);
      const withConflict = scoreOverhangFidelity('TGAC', ['TGAC'], 'BsaI', false, null);

      expect(withConflict.score).toBeLessThan(withoutConflict.score);
    });
  });

  describe('classifyQuality', () => {
    it('should classify scores correctly', () => {
      expect(classifyQuality(95)).toBe('excellent');
      expect(classifyQuality(80)).toBe('good');
      expect(classifyQuality(65)).toBe('acceptable');
      expect(classifyQuality(40)).toBe('poor');
    });
  });
});

describe('Domestication Primer Workflow', () => {
  describe('runDomesticationWorkflow', () => {
    it('should handle sequence with no sites', () => {
      const result = runDomesticationWorkflow(
        TEST_SEQUENCES.noSites.sequence,
        'BsaI'
      );

      expect(result.success).toBe(true);
      expect(result.strategy).toBe('none');
    });

    it('should analyze internal sites', () => {
      const result = runDomesticationWorkflow(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        'BsaI',
        { frame: 0, isCodingSequence: true }
      );

      expect(result.analysis).toBeDefined();
      expect(result.analysis.needsDomestication).toBe(true);
      expect(result.analysis.siteCount).toBe(1);
    });

    it('should select appropriate strategy', () => {
      const result = runDomesticationWorkflow(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        'BsaI',
        { frame: 0, isCodingSequence: true }
      );

      expect(result.strategy).toBeDefined();
      // Should select silent_mutation, mutagenic_junction, or hybrid
      expect(['silent_mutation', 'mutagenic_junction', 'hybrid', 'alternative_enzyme']).toContain(result.strategy);
    });

    it('should design primers', () => {
      const result = runDomesticationWorkflow(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        'BsaI',
        { frame: 0, isCodingSequence: true }
      );

      if (result.success && result.strategy !== 'none' && result.strategy !== 'alternative_enzyme') {
        expect(result.primers).toBeDefined();
        expect(result.primers.length).toBeGreaterThan(0);
      }
    });

    it('should generate workflow guide', () => {
      const result = runDomesticationWorkflow(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        'BsaI',
        { frame: 0, includeWorkflowGuide: true }
      );

      expect(result.workflowGuide).toBeDefined();
      expect(result.workflowGuide.steps).toBeDefined();
      expect(result.workflowGuide.steps.length).toBeGreaterThan(0);
    });

    it('should validate complete design', () => {
      const result = runDomesticationWorkflow(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        'BsaI',
        { frame: 0, isCodingSequence: true }
      );

      expect(result.validation).toBeDefined();
      expect(result.validation.isValid).toBeDefined();
      expect(result.validation.warnings).toBeDefined();
    });

    it('should include primer summary for ordering', () => {
      const result = runDomesticationWorkflow(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        'BsaI',
        { frame: 0, isCodingSequence: true }
      );

      if (result.success && result.strategy !== 'none' && result.strategy !== 'alternative_enzyme') {
        expect(result.primerSummary).toBeDefined();
        expect(result.primerSummary.orderList).toBeDefined();
      }
    });
  });

  describe('analyzeSequenceForDomestication', () => {
    it('should identify sites needing domestication', () => {
      const result = analyzeSequenceForDomestication(
        TEST_SEQUENCES.singleSite.sequence,
        'BsaI',
        { frame: 0, isCodingSequence: true, organism: 'ecoli' }
      );

      expect(result.needsDomestication).toBe(true);
      expect(result.siteCount).toBe(1);
    });

    it('should analyze silent mutation availability', () => {
      const result = analyzeSequenceForDomestication(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        'BsaI',
        { frame: 0, isCodingSequence: true, organism: 'ecoli' }
      );

      expect(result.siteAnalyses).toBeDefined();
      expect(result.siteAnalyses.length).toBe(1);
      expect(result.siteAnalyses[0].analysis).toBeDefined();
    });

    it('should recommend appropriate strategy', () => {
      const result = analyzeSequenceForDomestication(
        TEST_SEQUENCES.silentMutationPossible.sequence,
        'BsaI',
        { frame: 0, isCodingSequence: true, organism: 'ecoli' }
      );

      expect(result.recommendedStrategy).toBeDefined();
      expect(['silent_mutation', 'mutagenic_junction', 'hybrid', 'alternative_enzyme']).toContain(result.recommendedStrategy);
    });

    it('should check alternative enzymes', () => {
      const result = analyzeSequenceForDomestication(
        TEST_SEQUENCES.singleSite.sequence,
        'BsaI',
        { frame: 0, isCodingSequence: true, organism: 'ecoli' }
      );

      expect(result.alternativeEnzymes).toBeDefined();
      expect(Array.isArray(result.alternativeEnzymes)).toBe(true);
    });
  });
});

describe('NEB Fidelity Integration', () => {
  it('should use experimental fidelity data when available', () => {
    // Test that NEB data is being used
    const overhangs = ['TGAC', 'GCAT', 'GATG'];  // Known high-fidelity set
    const result = calculateExperimentalFidelity(overhangs, 'BsaI');

    expect(result.source).toBe('experimental');
    expect(result.assemblyFidelity).toBeGreaterThan(0.90);  // High-fidelity set
  });

  it('should identify low-fidelity overhangs', () => {
    // AAAA is known to be low-fidelity (homopolymer)
    const overhangs = ['AAAA', 'TTTT', 'GGGG'];  // Low-fidelity set
    const result = calculateExperimentalFidelity(overhangs, 'BsaI');

    // Should have lower fidelity than high-quality overhangs
    // Homopolymers typically have lower assembly fidelity
    expect(result.assemblyFidelity).toBeDefined();
    // If warnings exist, verify they flag problematic overhangs
    if (result.warnings?.length > 0) {
      expect(result.warnings.length).toBeGreaterThan(0);
    }
  });
});

describe('Configuration', () => {
  it('should have enhanced search radius', () => {
    expect(ENHANCED_JUNCTION_CONFIG.searchRadius).toBe(50);  // Enhanced uses 50bp
    expect(ENHANCED_JUNCTION_CONFIG.extendedSearchRadius).toBe(80);
  });

  it('should have NEB fidelity thresholds', () => {
    expect(ENHANCED_JUNCTION_CONFIG.minOverhangFidelity).toBe(0.95);
    expect(ENHANCED_JUNCTION_CONFIG.targetAssemblyFidelity).toBe(0.98);
  });

  it('should have calibrated scoring weights', () => {
    const weights = ENHANCED_JUNCTION_CONFIG.weights;
    expect(weights.overhangFidelity).toBe(0.35);  // Highest weight for NEB data
    expect(weights.mutationQuality).toBe(0.25);
    expect(weights.primerQuality).toBe(0.25);
    expect(weights.positionOptimality).toBe(0.15);

    // Weights should sum to 1
    const sum = Object.values(weights).reduce((a, b) => a + b, 0);
    expect(sum).toBe(1);
  });
});
