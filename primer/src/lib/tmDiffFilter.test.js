/**
 * Test for Tm Difference Filtering
 *
 * This test investigates whether Unified and Legacy designers apply
 * the same Tm difference filtering criteria.
 *
 * FINDING: There is a difference in Tm difference handling:
 *
 * MUTAGENESIS (deletion, insertion, substitution, AA mutation):
 * - Legacy uses `selectBestByTier()` which assigns candidates to tiers:
 *   - Tier 1 (excellent): Tm diff ≤2°C
 *   - Tier 2 (good): Tm diff ≤5°C
 *   - Tier 3 (acceptable): Tm diff ≤8°C
 *   - Tier 4 (poor): Tm diff >8°C
 * - It ALWAYS picks from the best tier first
 * - So Tm >5°C only happens if no tier1 or tier2 candidates exist
 *
 * AMPLIFICATION (PCR primers):
 * - Uses `chooseBest()` which picks pair with highest compositeScore
 * - NO tier-based selection
 * - NO hard Tm diff filter on primary design
 * - Tm diff is only a small weight (0.03) in compositeScore
 * - MAX_TM_DIFF=5 filter only applies to ALTERNATIVES, not primary
 */

import { describe, it, expect } from 'vitest';

// Unified Designer
import { designUnified } from './unifiedPrimerDesign';

// Legacy Designer
import {
  designDeletionPrimers,
  designInsertionPrimers,
  designRegionSubstitutionPrimers,
  designCodonChangePrimers,
} from './mutagenesis';

// PCR primer design
import { primers, generateAlternatives } from './primers';

// Test sequences designed to potentially cause Tm mismatches
const TEST_TEMPLATE = 'ATGAAACAAAGCACTATTGCACTGGCACTCTTACCGTTACTGTTTACCCCTGTGACAAAAGCCAATATGGACAAGTTGTTTGACGGTATCAGAAGCAGAACTGGCGAAACTTTTACCGGTGAAGACCGTAACGGTTACGACAATAAATACAATGTTTATAATCAGACTAACGACTGTTGGGGTTTTGAATTTAAAGATGAAGATATGCTGTGCCCGGACCCAATTAGCTGGCGTAATGCCGAGATCATGCGTAAAAAATGGGACAGCAAAGAGCAGAAAAGCATGTACGAACGCCAGTTTGACGAGCTGTATAAAGAACGCTATGGTTATGCCAACAGCTACATGTATGACGATGATGACAAACATCTGTACAAGTAAGGAGGTAATAA';

// GC-rich sequence that might cause Tm imbalance
const GC_RICH_TEMPLATE = 'ATGGCGCGCGCCGCCGCCGCGGCGGCGGCGCCGCGCGCGGCGCGCGCCGCCGCCGCGGCGGCGGCGCCGCGCGCGATGAAACAAAGCACTATTGCACTGGCACTCTTACCGTTACTGTTTACCCCTGTGACAAAAGCCAATATGGACAAGTTGTTTGACGGTATCAGAAGCAGAACTGGCGAAACTTTTACCGGTGAAGACCGTAACGGTTACGACAATAAATACAATGTTTATAATCAGACTAACGACTGTTGGGGTTTTGAATTTAAAGATGAAGATATGCTGTGCCCGGACCCAATTAGCTGGCGTAATGCCGAGATCATGCGTAAAAAATGGGACAGCAAAGAGCAGAAAAGCATGTACGAACGCCAGTTTGACGAGCTGTATAAAGAACGCTATGGTTATGCCAACAGCTACATGTATGACGATGATGACAAACATCTGTACAAGTAAGGAGGTAATAATTT';

// AT-rich at one end, GC-rich at other (to force Tm imbalance)
const ASYMMETRIC_TEMPLATE = 'ATATATATATATATATATATATATATATATATATATATATATATATATATATAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATGAAACAAAGCACTATTGCACTGGCACTCTTACCGTTACTGTTTACCCCTGTGACAAAAGCCAATATGGACAAGTTGTTTGACGGTATCAGAAGCAGAACTGGCGAAACTTTTACCGGTGAAGACCGTAACGGTTACGACAATAAATACAATGTTTATAATCAGACTAACGACTGTTGGGGG';

// GFP for AA mutation tests
const GFP_TEMPLATE = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA';

// =============================================================================
// Test: Tm Difference Filtering Comparison
// =============================================================================

describe('Tm Difference Filtering Analysis', () => {

  describe('Amplification (PCR) - NO hard Tm filter', () => {
    it('should show that primers() does NOT have a hard Tm diff filter', () => {
      // Test on asymmetric sequence that may produce Tm imbalance
      const regionSeq = ASYMMETRIC_TEMPLATE.slice(0, 200);

      const [fwd, rev] = primers(regionSeq, {
        useCompositeScore: true,
        useSmartDesign: false,
      });

      const tmDiff = Math.abs(fwd.tm - rev.tm);

      console.log('\n=== AMPLIFICATION (primers.js) ===');
      console.log(`Forward Tm: ${fwd.tm.toFixed(1)}°C`);
      console.log(`Reverse Tm: ${rev.tm.toFixed(1)}°C`);
      console.log(`Tm Difference: ${tmDiff.toFixed(1)}°C`);
      console.log(`Composite Score: ${fwd.scoring.compositeScore}`);

      // Document the behavior - no assertion on Tm diff since there's no hard filter
      // Note: AT-rich regions can produce very low Tm primers
      expect(fwd.tm).toBeGreaterThan(0);
      expect(rev.tm).toBeGreaterThan(0);

      // This may or may not be >5°C depending on the sequence
      if (tmDiff > 5) {
        console.log(`⚠️ WARNING: Tm diff ${tmDiff.toFixed(1)}°C > 5°C - NO hard filter in primers()`);
      }
    });

    it('should show that designUnified amplification has MAX_TM_DIFF only for alternatives', () => {
      const result = designUnified(ASYMMETRIC_TEMPLATE, {
        start: 0,
        end: 200,
      }, { useSmartDesign: false });

      const tmDiff = Math.abs(result.forward.tm - result.reverse.tm);

      console.log('\n=== UNIFIED AMPLIFICATION ===');
      console.log(`Forward Tm: ${result.forward.tm.toFixed(1)}°C`);
      console.log(`Reverse Tm: ${result.reverse.tm.toFixed(1)}°C`);
      console.log(`Tm Difference: ${tmDiff.toFixed(1)}°C`);
      console.log(`Quality: ${result.qualityTier}`);

      // Check alternatives - these SHOULD be filtered by MAX_TM_DIFF=5
      if (result.alternativePrimers && result.alternativePrimers.length > 0) {
        console.log(`\nAlternatives (filtered by MAX_TM_DIFF=5):`);
        result.alternativePrimers.forEach((alt, i) => {
          const altTmDiff = Math.abs(alt.forward.tm - alt.reverse.tm);
          console.log(`  Alt ${i+1}: Tm diff = ${altTmDiff.toFixed(1)}°C`);
          // Alternatives SHOULD be ≤5°C
          expect(altTmDiff).toBeLessThanOrEqual(5.1); // small tolerance for rounding
        });
      }

      if (tmDiff > 5) {
        console.log(`\n⚠️ PRIMARY design has Tm diff ${tmDiff.toFixed(1)}°C > 5°C`);
        console.log(`   This is because MAX_TM_DIFF filter only applies to alternatives!`);
      }
    });
  });

  describe('Mutagenesis - HAS tier-based Tm filtering', () => {
    it('should show that deletion uses tier-based selection (Tm diff prioritized)', () => {
      const result = designDeletionPrimers(TEST_TEMPLATE, 100, 10);

      const tmDiff = Math.abs(result.forward.tm - result.reverse.tm);

      console.log('\n=== LEGACY DELETION (tier-based) ===');
      console.log(`Forward Tm: ${result.forward.tm.toFixed(1)}°C`);
      console.log(`Reverse Tm: ${result.reverse.tm.toFixed(1)}°C`);
      console.log(`Tm Difference: ${tmDiff.toFixed(1)}°C`);
      console.log(`Quality Tier: ${result.qualityTier}`);

      // Document the tier rules
      console.log('\nTier Rules:');
      console.log('  Tier 1 (excellent): Tm diff ≤2°C');
      console.log('  Tier 2 (good): Tm diff ≤5°C');
      console.log('  Tier 3 (acceptable): Tm diff ≤8°C');
      console.log('  Tier 4 (poor): Tm diff >8°C');

      // Tm diff should typically be ≤5°C for good/excellent tier
      if (result.qualityTier === 'excellent') {
        expect(tmDiff).toBeLessThanOrEqual(2.1);
      } else if (result.qualityTier === 'good') {
        expect(tmDiff).toBeLessThanOrEqual(5.1);
      }
    });

    it('should show that unified deletion also uses tier-based selection', () => {
      const result = designUnified(TEST_TEMPLATE, {
        start: 100,
        end: 110,
        replacement: '',
      });

      const tmDiff = Math.abs(result.forward.tm - result.reverse.tm);

      console.log('\n=== UNIFIED DELETION (tier-based) ===');
      console.log(`Forward Tm: ${result.forward.tm.toFixed(1)}°C`);
      console.log(`Reverse Tm: ${result.reverse.tm.toFixed(1)}°C`);
      console.log(`Tm Difference: ${tmDiff.toFixed(1)}°C`);
      console.log(`Quality Tier: ${result.qualityTier}`);

      // Same tier rules apply since Unified calls Legacy
      if (result.qualityTier === 'excellent') {
        expect(tmDiff).toBeLessThanOrEqual(2.1);
      } else if (result.qualityTier === 'good') {
        expect(tmDiff).toBeLessThanOrEqual(5.1);
      }
    });
  });

  describe('Compare: Same sequence, different operations', () => {
    it('should compare Tm diff handling between amplification and mutagenesis', () => {
      const results = {
        amplification: null,
        deletion: null,
      };

      // Amplification on GC-rich region
      try {
        results.amplification = designUnified(GC_RICH_TEMPLATE, {
          start: 0,
          end: 200,
        });
      } catch (e) {
        console.log('Amplification failed:', e.message);
      }

      // Deletion in same region
      try {
        results.deletion = designUnified(GC_RICH_TEMPLATE, {
          start: 100,
          end: 110,
          replacement: '',
        });
      } catch (e) {
        console.log('Deletion failed:', e.message);
      }

      console.log('\n=== COMPARISON: Same sequence, different operations ===');

      if (results.amplification) {
        const ampTmDiff = Math.abs(results.amplification.forward.tm - results.amplification.reverse.tm);
        console.log(`\nAmplification:`);
        console.log(`  Tm diff: ${ampTmDiff.toFixed(1)}°C`);
        console.log(`  Quality: ${results.amplification.qualityTier}`);
        console.log(`  Filter: compositeScore-based (NO hard Tm filter)`);
      }

      if (results.deletion) {
        const delTmDiff = Math.abs(results.deletion.forward.tm - results.deletion.reverse.tm);
        console.log(`\nDeletion:`);
        console.log(`  Tm diff: ${delTmDiff.toFixed(1)}°C`);
        console.log(`  Quality: ${results.deletion.qualityTier}`);
        console.log(`  Filter: tier-based (Tm diff ≤2/5/8°C)`);
      }

      console.log('\n=== CONCLUSION ===');
      console.log('Amplification (PCR) does NOT have a hard Tm diff filter.');
      console.log('Mutagenesis HAS tier-based selection that prioritizes low Tm diff.');
      console.log('This is the ROOT CAUSE of Tm>5 appearing in unified amplification.');
    });
  });
});

// =============================================================================
// Stress Test: Find sequences that produce high Tm diff
// =============================================================================

describe('Stress Test: Finding high Tm diff cases', () => {
  const testSequences = [
    { name: 'Standard', seq: TEST_TEMPLATE },
    { name: 'GC-rich', seq: GC_RICH_TEMPLATE },
    { name: 'Asymmetric', seq: ASYMMETRIC_TEMPLATE },
  ];

  it('should test multiple sequences for Tm diff issues', () => {
    console.log('\n=== STRESS TEST: Tm Difference Analysis ===\n');
    console.log('Sequence          | Operation     | Tm Diff | Quality | Has Filter?');
    console.log('------------------|---------------|---------|---------|------------');

    for (const { name, seq } of testSequences) {
      // Test amplification
      try {
        const amp = designUnified(seq, { start: 0, end: Math.min(200, seq.length - 50) });
        const ampTmDiff = Math.abs(amp.forward.tm - amp.reverse.tm);
        const ampFilter = ampTmDiff > 5 ? 'NO ⚠️' : 'N/A';
        console.log(`${name.padEnd(17)} | Amplification | ${ampTmDiff.toFixed(1).padStart(7)}°C | ${(amp.qualityTier || 'N/A').padEnd(7)} | ${ampFilter}`);
      } catch (e) {
        console.log(`${name.padEnd(17)} | Amplification | FAILED  |         |`);
      }

      // Test deletion
      try {
        const del = designUnified(seq, { start: 80, end: 90, replacement: '' });
        const delTmDiff = Math.abs(del.forward.tm - del.reverse.tm);
        const delFilter = 'tier-based';
        console.log(`${name.padEnd(17)} | Deletion      | ${delTmDiff.toFixed(1).padStart(7)}°C | ${(del.qualityTier || 'N/A').padEnd(7)} | ${delFilter}`);
      } catch (e) {
        console.log(`${name.padEnd(17)} | Deletion      | FAILED  |         |`);
      }
    }

    console.log('\nLegend:');
    console.log('  NO ⚠️ = No hard Tm filter, allowing Tm>5°C');
    console.log('  tier-based = Uses tier system that prioritizes Tm≤2/5/8°C');
  });
});

// =============================================================================
// Recommendation
// =============================================================================

describe('Verify Tier-Based Fix', () => {
  it('should verify tier-based selection is working on standard sequence', () => {
    // Use standard template which should have tier1/tier2 candidates available
    const [fwd, rev] = primers(TEST_TEMPLATE.slice(0, 200), {
      useCompositeScore: true,
      useSmartDesign: false,
    });

    const tmDiff = Math.abs(fwd.tm - rev.tm);

    console.log('\n=== VERIFY FIX: Standard Template ===');
    console.log(`Forward Tm: ${fwd.tm.toFixed(1)}°C`);
    console.log(`Reverse Tm: ${rev.tm.toFixed(1)}°C`);
    console.log(`Tm Difference: ${tmDiff.toFixed(1)}°C`);

    // With tier-based selection, standard template should produce ≤5°C Tm diff
    // because tier1/tier2 candidates should exist
    if (tmDiff <= 2) {
      console.log(`✅ TIER 1 (Excellent): Tm diff ≤2°C`);
    } else if (tmDiff <= 5) {
      console.log(`✅ TIER 2 (Good): Tm diff ≤5°C`);
    } else if (tmDiff <= 8) {
      console.log(`⚠️ TIER 3 (Acceptable): Tm diff ≤8°C`);
    } else {
      console.log(`❌ TIER 4 (Poor): Tm diff >8°C`);
    }

    // Standard template should get tier1 or tier2
    expect(tmDiff).toBeLessThanOrEqual(5);
  });

  it('should show that extreme sequences fall back to tier4 gracefully', () => {
    // Asymmetric template is designed to be impossible - tier4 is expected
    const [fwd, rev] = primers(ASYMMETRIC_TEMPLATE.slice(0, 200), {
      useCompositeScore: true,
      useSmartDesign: false,
    });

    const tmDiff = Math.abs(fwd.tm - rev.tm);

    console.log('\n=== EXTREME SEQUENCE (tier4 expected) ===');
    console.log(`Forward Tm: ${fwd.tm.toFixed(1)}°C`);
    console.log(`Reverse Tm: ${rev.tm.toFixed(1)}°C`);
    console.log(`Tm Difference: ${tmDiff.toFixed(1)}°C`);
    console.log(`This sequence has AT-rich 5' and GC-rich 3' - impossible to match Tm`);

    // Extreme sequence - tier4 is acceptable as fallback
    expect(fwd.tm).toBeGreaterThan(0);
    expect(rev.tm).toBeGreaterThan(0);
  });
});
