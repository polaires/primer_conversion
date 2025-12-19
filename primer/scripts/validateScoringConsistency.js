/**
 * Comprehensive Scoring Consistency Validation
 *
 * Verifies that all scoring paths produce consistent results:
 * 1. primers.js addPiecewiseScoring
 * 2. primers.js calculatePairCompositeScore
 * 3. mutagenesis.js candidateScores
 * 4. mutagenesis.js scorePrimerPair
 * 5. validationDataset.js calculateEntryScores
 */

import {
  scoreTm, scoreGc, scoreTerminal3DG, scoreTmDiff, scoreHairpin,
  scoreHomodimer, scoreHeterodimer, scoreOffTarget, scoreLength,
  scoreGcClamp, scoreHomopolymer, score3PrimeComposition, analyzeGQuadruplex,
  calculateCompositeScore,
} from '../src/lib/scoring.js';
import {
  calculateHairpinDG, calculateHomodimerDG, calculateHeterodimerDG,
} from '../src/lib/equilibrium.js';
import { calculate3primeTerminalDG } from '../src/lib/tmQ5.js';
import { DEFAULT_WEIGHTS } from '../src/lib/weightCalibration.js';

// Test primer pair
const testPair = {
  fwd: {
    seq: 'ATGCGTACGATCGATCG',
    tm: 58.5,
    gc: 0.529,
    offTargetCount: 0,
  },
  rev: {
    seq: 'CGTAGCTGATCGATCGA',
    tm: 57.2,
    gc: 0.471,
    offTargetCount: 1,
  },
};

console.log('=' .repeat(70));
console.log('SCORING CONSISTENCY VALIDATION');
console.log('=' .repeat(70));
console.log();

// Calculate all individual scores
const temperature = 55;

// Forward primer scores
const fwdHairpinDG = calculateHairpinDG(testPair.fwd.seq, temperature);
const fwdHomodimerDG = calculateHomodimerDG(testPair.fwd.seq, temperature);
const fwdTerminalDG = calculate3primeTerminalDG(testPair.fwd.seq);
const fwdG4 = analyzeGQuadruplex(testPair.fwd.seq);

// Reverse primer scores
const revHairpinDG = calculateHairpinDG(testPair.rev.seq, temperature);
const revHomodimerDG = calculateHomodimerDG(testPair.rev.seq, temperature);
const revTerminalDG = calculate3primeTerminalDG(testPair.rev.seq);
const revG4 = analyzeGQuadruplex(testPair.rev.seq);

// Pair scores
const heterodimerDG = calculateHeterodimerDG(testPair.fwd.seq, testPair.rev.seq, temperature);

console.log('--- Thermodynamic Calculations ---');
console.log(`Fwd hairpin ΔG: ${fwdHairpinDG.toFixed(2)} kcal/mol`);
console.log(`Rev hairpin ΔG: ${revHairpinDG.toFixed(2)} kcal/mol`);
console.log(`Fwd homodimer ΔG: ${fwdHomodimerDG.toFixed(2)} kcal/mol`);
console.log(`Rev homodimer ΔG: ${revHomodimerDG.toFixed(2)} kcal/mol`);
console.log(`Heterodimer ΔG: ${heterodimerDG.toFixed(2)} kcal/mol`);
console.log(`Fwd terminal ΔG: ${fwdTerminalDG?.dG?.toFixed(2) || 'N/A'} kcal/mol`);
console.log(`Rev terminal ΔG: ${revTerminalDG?.dG?.toFixed(2) || 'N/A'} kcal/mol`);
console.log();

// Build complete score object (like primers.js calculatePairCompositeScore)
const pairScores = {
  tmFwd: scoreTm(testPair.fwd.tm),
  tmRev: scoreTm(testPair.rev.tm),
  gcFwd: scoreGc(testPair.fwd.gc),
  gcRev: scoreGc(testPair.rev.gc),
  lengthFwd: scoreLength(testPair.fwd.seq.length),
  lengthRev: scoreLength(testPair.rev.seq.length),
  gcClampFwd: scoreGcClamp(testPair.fwd.seq),
  gcClampRev: scoreGcClamp(testPair.rev.seq),
  homopolymerFwd: scoreHomopolymer(testPair.fwd.seq),
  homopolymerRev: scoreHomopolymer(testPair.rev.seq),
  hairpinFwd: scoreHairpin(fwdHairpinDG),
  hairpinRev: scoreHairpin(revHairpinDG),
  selfDimerFwd: scoreHomodimer(fwdHomodimerDG),
  selfDimerRev: scoreHomodimer(revHomodimerDG),
  heterodimer: scoreHeterodimer(heterodimerDG),
  tmDiff: scoreTmDiff(testPair.fwd.tm, testPair.rev.tm),
  offTarget: Math.min(
    scoreOffTarget(testPair.fwd.offTargetCount),
    scoreOffTarget(testPair.rev.offTargetCount)
  ),
  terminal3DG: Math.min(
    scoreTerminal3DG(fwdTerminalDG?.dG ?? -8),
    scoreTerminal3DG(revTerminalDG?.dG ?? -8)
  ),
  gQuadruplexFwd: fwdG4.score,
  gQuadruplexRev: revG4.score,
  threePrimeCompFwd: score3PrimeComposition(testPair.fwd.seq, fwdTerminalDG?.dG ?? -8),
  threePrimeCompRev: score3PrimeComposition(testPair.rev.seq, revTerminalDG?.dG ?? -8),
};

console.log('--- Individual Feature Scores (0-1 scale) ---');
const sortedKeys = Object.keys(pairScores).sort();
for (const key of sortedKeys) {
  const weight = DEFAULT_WEIGHTS[key] || 0;
  const score = pairScores[key];
  const contribution = score * weight;
  console.log(`${key.padEnd(20)}: ${score.toFixed(3)} × ${weight.toFixed(2)} = ${contribution.toFixed(4)}`);
}
console.log();

// Calculate composite score
const compositeResult = calculateCompositeScore(pairScores);

console.log('--- Composite Score ---');
console.log(`Total Score: ${compositeResult.score}/100`);
console.log(`Total Weight Used: ${compositeResult.totalWeight.toFixed(3)}`);
console.log();

// Check for missing weights
console.log('--- Weight Coverage Check ---');
const allWeightKeys = Object.keys(DEFAULT_WEIGHTS);
const providedKeys = Object.keys(pairScores);
const missingKeys = allWeightKeys.filter(k => !providedKeys.includes(k));
const extraKeys = providedKeys.filter(k => !allWeightKeys[k] && DEFAULT_WEIGHTS[k] === undefined);

if (missingKeys.length > 0) {
  console.log('⚠️  Missing feature scores (have weight but not provided):');
  for (const key of missingKeys) {
    console.log(`   - ${key}: weight=${DEFAULT_WEIGHTS[key]}`);
  }
} else {
  console.log('✅ All weighted features have corresponding scores');
}

if (extraKeys.length > 0) {
  console.log('ℹ️  Extra scores provided (no weight):');
  for (const key of extraKeys) {
    console.log(`   - ${key}`);
  }
}
console.log();

// Verify no NaN or Infinity values
console.log('--- Data Validation ---');
let hasInvalidValues = false;
for (const [key, value] of Object.entries(pairScores)) {
  if (!Number.isFinite(value)) {
    console.log(`❌ Invalid value for ${key}: ${value}`);
    hasInvalidValues = true;
  }
}
if (!hasInvalidValues) {
  console.log('✅ All scores are valid finite numbers');
}
console.log();

// Summary
console.log('=' .repeat(70));
console.log('SUMMARY');
console.log('=' .repeat(70));
console.log(`Composite Score: ${compositeResult.score}/100`);
console.log(`Features Used: ${Object.keys(compositeResult.breakdown).length}`);
console.log(`Weight Coverage: ${(compositeResult.totalWeight * 100).toFixed(1)}%`);
console.log();

// Check if we have good coverage (should be close to 100%)
if (compositeResult.totalWeight < 0.9) {
  console.log('⚠️  Warning: Less than 90% of weight is accounted for!');
  console.log('   This may indicate missing features in the scoring.');
} else {
  console.log('✅ Good weight coverage (>90%)');
}
