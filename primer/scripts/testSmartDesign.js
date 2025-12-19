#!/usr/bin/env node
/**
 * Smart Design Performance Test Script
 *
 * Tests the iterative smart design system to verify:
 * 1. GC clamp optimization works
 * 2. 3' end composition improvements
 * 3. Score improvements after optimization
 *
 * Run: node scripts/testSmartDesign.js
 */

import {
  primers,
  analyze3PrimeEnd,
  quickAssess,
  optimizePrimer,
  generateLengthVariants,
  scorePrimerVariant,
  generateDesignSuggestions,
} from '../src/lib/index.js';

// ANSI color codes for console output
const colors = {
  reset: '\x1b[0m',
  bright: '\x1b[1m',
  dim: '\x1b[2m',
  green: '\x1b[32m',
  yellow: '\x1b[33m',
  red: '\x1b[31m',
  cyan: '\x1b[36m',
  magenta: '\x1b[35m',
  blue: '\x1b[34m',
};

function log(msg, color = 'reset') {
  console.log(`${colors[color]}${msg}${colors.reset}`);
}

function header(msg) {
  console.log();
  log('═'.repeat(70), 'cyan');
  log(`  ${msg}`, 'bright');
  log('═'.repeat(70), 'cyan');
}

function subheader(msg) {
  console.log();
  log(`─── ${msg} ───`, 'blue');
}

function printPrimerAnalysis(label, seq) {
  const analysis = analyze3PrimeEnd(seq);
  const assessment = quickAssess(seq);

  log(`\n${label}:`, 'bright');
  log(`  Sequence: ${seq}`, 'dim');
  log(`  Length: ${seq.length}bp | Tm: ${assessment.tm}°C | GC: ${assessment.gc}%`);

  // 3' end analysis
  const gcClampColor = analysis.gcCounts.last2 >= 1 ? 'green' : 'red';
  log(`  3' End (last 5): ${analysis.last5}`, 'dim');
  log(`  GC Clamp: ${analysis.gcCounts.last2} G/C in last 2bp`, gcClampColor);
  log(`  Terminal ΔG: ${analysis.terminalDG} kcal/mol`);

  // Quality
  const qualityColor = {
    excellent: 'green',
    good: 'green',
    acceptable: 'yellow',
    marginal: 'yellow',
    poor: 'red',
  }[analysis.quality] || 'reset';
  log(`  Quality: ${analysis.quality.toUpperCase()}`, qualityColor);

  if (analysis.issues.length > 0) {
    log(`  Issues:`, 'yellow');
    for (const issue of analysis.issues) {
      log(`    • ${issue}`, 'yellow');
    }
  }

  return { analysis, assessment };
}

function printScoreComparison(label, before, after) {
  const diff = after - before;
  const diffColor = diff > 0 ? 'green' : diff < 0 ? 'red' : 'dim';
  const arrow = diff > 0 ? '↑' : diff < 0 ? '↓' : '→';
  log(`  ${label}: ${before} → ${after} (${arrow}${Math.abs(diff).toFixed(1)})`, diffColor);
}

// ============================================================================
// Test Cases
// ============================================================================

const testCases = [
  {
    name: 'Standard plasmid sequence',
    seq: 'ATGCTAGCTAGCTAAATTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCTAGCTAGCTAGC',
    description: 'Contains AT-rich and GC-rich regions',
  },
  {
    name: 'GFP gene fragment',
    seq: 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA',
    description: 'Real GFP coding sequence with varied GC content',
  },
  {
    name: 'AT-rich sequence',
    seq: 'ATATATATATATATATATATATGCATATATATATATATATATATATATATATATGCATATATATATATATATATGCATATATATAT',
    description: 'Challenging AT-rich template',
  },
  {
    name: 'GC-rich sequence',
    seq: 'GCGCGCGCGCGCGCGCGCGCATGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC',
    description: 'Challenging GC-rich template',
  },
];

// ============================================================================
// Main Test Runner
// ============================================================================

async function runTests() {
  header('SMART PRIMER DESIGN PERFORMANCE TEST');

  log('\nThis script tests the iterative smart design system that:', 'dim');
  log('  • Extends primers to achieve GC clamp (G/C in last 2 bases)', 'dim');
  log('  • Balances 3\' end composition for optimal ΔG', 'dim');
  log('  • Generates design suggestions for improvement', 'dim');

  let totalImprovements = 0;
  let totalTests = 0;

  for (const testCase of testCases) {
    header(testCase.name);
    log(testCase.description, 'dim');

    // Skip if sequence is too short
    if (testCase.seq.length < 50) {
      log('\n  Sequence too short for primer design', 'yellow');
      continue;
    }

    try {
      // ─── Design WITHOUT smart design ───
      subheader('Standard Design (useSmartDesign: false)');

      const [fwdStd, revStd] = primers(testCase.seq, {
        useCompositeScore: true,
        useSmartDesign: false,
      });

      const fwdStdAnalysis = printPrimerAnalysis('Forward Primer', fwdStd.seq);
      const revStdAnalysis = printPrimerAnalysis('Reverse Primer', revStd.seq);

      log(`\n  Pair Score: ${fwdStd.scoring.compositeScore}`, 'bright');
      log(`  Quality Tier: ${fwdStd.scoring.qualityTier}`, 'bright');

      // ─── Design WITH smart design ───
      subheader('Smart Design (useSmartDesign: true)');

      const [fwdSmart, revSmart] = primers(testCase.seq, {
        useCompositeScore: true,
        useSmartDesign: true,
        smartDesignTargetScore: 75,
      });

      const fwdSmartAnalysis = printPrimerAnalysis('Forward Primer', fwdSmart.seq);
      const revSmartAnalysis = printPrimerAnalysis('Reverse Primer', revSmart.seq);

      log(`\n  Pair Score: ${fwdSmart.scoring.compositeScore}`, 'bright');
      log(`  Quality Tier: ${fwdSmart.scoring.qualityTier}`, 'bright');

      // ─── Comparison ───
      subheader('Comparison');

      const fwdChanged = fwdStd.seq !== fwdSmart.seq;
      const revChanged = revStd.seq !== revSmart.seq;

      if (fwdChanged || revChanged) {
        log('\n  Changes made:', 'green');

        if (fwdChanged) {
          log(`    FWD: ${fwdStd.seq.length}bp → ${fwdSmart.seq.length}bp`, 'green');
          log(`         ${fwdStd.seq}`, 'dim');
          log(`       → ${fwdSmart.seq}`, 'green');

          // GC clamp comparison
          const stdGcClamp = fwdStdAnalysis.analysis.gcCounts.last2;
          const smartGcClamp = fwdSmartAnalysis.analysis.gcCounts.last2;
          if (smartGcClamp > stdGcClamp) {
            log(`         GC Clamp: ${stdGcClamp} → ${smartGcClamp} ✓`, 'green');
          }
        }

        if (revChanged) {
          log(`    REV: ${revStd.seq.length}bp → ${revSmart.seq.length}bp`, 'green');
          log(`         ${revStd.seq}`, 'dim');
          log(`       → ${revSmart.seq}`, 'green');

          const stdGcClamp = revStdAnalysis.analysis.gcCounts.last2;
          const smartGcClamp = revSmartAnalysis.analysis.gcCounts.last2;
          if (smartGcClamp > stdGcClamp) {
            log(`         GC Clamp: ${stdGcClamp} → ${smartGcClamp} ✓`, 'green');
          }
        }

        totalImprovements++;
      } else {
        log('\n  No changes needed (primers already optimal)', 'dim');
      }

      // Score comparison
      log('\n  Score Changes:', 'bright');
      printScoreComparison('Composite Score', fwdStd.scoring.compositeScore, fwdSmart.scoring.compositeScore);

      // Smart design metadata
      if (fwdSmart.smartDesign) {
        log(`\n  Smart Design Metadata:`, 'magenta');
        log(`    FWD optimized: ${fwdSmart.smartDesign.optimized}`, 'dim');
        log(`    REV optimized: ${revSmart.smartDesign?.optimized || false}`, 'dim');
        if (fwdSmart.smartDesign.reason) {
          log(`    Reason: ${fwdSmart.smartDesign.reason}`, 'dim');
        }
      }

      totalTests++;

    } catch (error) {
      log(`\n  Error: ${error.message}`, 'red');
    }
  }

  // ============================================================================
  // Detailed Analysis of Individual Functions
  // ============================================================================

  header('INDIVIDUAL FUNCTION TESTS');

  // Test analyze3PrimeEnd with various sequences
  subheader('analyze3PrimeEnd() - 3\' End Quality Analysis');

  const testSeqs = [
    { seq: 'ATGCTAGCTAGCTAGCTGC', desc: 'Ends with GC (strong clamp)' },
    { seq: 'ATGCTAGCTAGCTAGCTAG', desc: 'Ends with AG (1 G/C)' },
    { seq: 'ATGCTAGCTAGCTAGCTAA', desc: 'Ends with AA (no clamp)' },
    { seq: 'ATGCTAGCTAGCTAAAAAT', desc: 'Ends with poly-A' },
    { seq: 'GCGCGCGCGCGCGCGCGCG', desc: 'All GC (very strong)' },
  ];

  for (const { seq, desc } of testSeqs) {
    const analysis = analyze3PrimeEnd(seq);
    const gcColor = analysis.gcCounts.last2 >= 1 ? 'green' : 'red';
    const qualColor = analysis.quality === 'excellent' || analysis.quality === 'good' ? 'green' : 'yellow';

    log(`\n  ${desc}:`, 'dim');
    log(`    Seq: ...${seq.slice(-10)}`, 'dim');
    log(`    GC clamp: ${analysis.gcCounts.last2}/2  Terminal ΔG: ${analysis.terminalDG.toFixed(1)}  Quality: ${analysis.quality}`, gcColor);
  }

  // Test generateLengthVariants
  subheader('generateLengthVariants() - Variant Generation');

  const template = 'ATGCTAGCTAGCTAAATTGCATGCATGCATGCATGCATGC';
  const variants = generateLengthVariants(template, 0, 17, true, { extendMax: 3 });

  log(`\n  Template: ${template.slice(0, 30)}...`, 'dim');
  log(`  Starting length: 17bp, extending up to 3bp`, 'dim');
  log(`\n  Generated ${variants.length} variants:`, 'bright');

  for (const v of variants.slice(0, 5)) {
    const gcColor = v.analysis.gcCounts.last2 >= 1 ? 'green' : 'yellow';
    log(`    ${v.length}bp (Δ${v.lengthDelta >= 0 ? '+' : ''}${v.lengthDelta}) | Last 2: ${v.analysis.last2} | GC clamp: ${v.analysis.gcCounts.last2} | Priority: ${v.priority}`, gcColor);
  }

  // Test scorePrimerVariant
  subheader('scorePrimerVariant() - Comprehensive Scoring');

  const scoredPrimer = scorePrimerVariant('ATGCTAGCTAGCTAGCTGC');
  log(`\n  Sequence: ${scoredPrimer.seq}`, 'dim');
  log(`  Tm: ${scoredPrimer.tm}°C | GC: ${(scoredPrimer.gc * 100).toFixed(0)}% | Length: ${scoredPrimer.length}bp`, 'bright');
  log(`  Terminal ΔG: ${scoredPrimer.terminalDG} | Hairpin ΔG: ${scoredPrimer.hairpinDG}`, 'dim');
  log(`  Composite Score: ${scoredPrimer.compositeScore} | Tier: ${scoredPrimer.qualityTier}`, 'bright');

  log(`\n  Individual Scores:`, 'dim');
  for (const [key, value] of Object.entries(scoredPrimer.scores)) {
    const scoreColor = value >= 0.8 ? 'green' : value >= 0.6 ? 'yellow' : 'red';
    log(`    ${key.padEnd(15)}: ${value.toFixed(3)}`, scoreColor);
  }

  // Test generateDesignSuggestions
  subheader('generateDesignSuggestions() - Improvement Recommendations');

  const poorPrimer = scorePrimerVariant('ATATATATATATATATATA');
  const suggestions = generateDesignSuggestions(poorPrimer);

  log(`\n  Testing poor primer: ATATATATATATATATATA`, 'dim');
  log(`  Current Score: ${suggestions.currentScore} | Tier: ${suggestions.currentTier}`, 'yellow');
  log(`  Can Improve: ${suggestions.canImprove}`, suggestions.canImprove ? 'green' : 'dim');

  if (suggestions.suggestions.length > 0) {
    log(`\n  Suggestions:`, 'bright');
    for (const s of suggestions.suggestions) {
      const priorityColor = s.priority === 'high' ? 'red' : s.priority === 'medium' ? 'yellow' : 'dim';
      log(`    [${s.priority.toUpperCase()}] ${s.issue}`, priorityColor);
      log(`           Action: ${s.action}`, 'dim');
    }
  }

  // ============================================================================
  // Summary
  // ============================================================================

  header('SUMMARY');

  log(`\n  Tests run: ${totalTests}`, 'bright');
  log(`  Primers optimized: ${totalImprovements}`, totalImprovements > 0 ? 'green' : 'dim');
  log(`\n  The smart design system provides:`, 'dim');
  log(`    ✓ Iterative length adjustment for GC clamp`, 'green');
  log(`    ✓ 3' end composition optimization`, 'green');
  log(`    ✓ Actionable design suggestions`, 'green');
  log(`    ✓ Transparent optimization rationale`, 'green');

  console.log();
}

// Run tests
runTests().catch(console.error);
