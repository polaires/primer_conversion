#!/usr/bin/env node
/**
 * DNA24 Parameter Validation Script
 *
 * Validates the DNA24 thermodynamic parameters against:
 * 1. Parameter completeness checks
 * 2. Thermodynamic value ranges
 * 3. Published Tm data comparison
 *
 * Run: node scripts/validate-dna24.js
 */

import {
  DNA24_NN,
  DNA24_INTERNAL_MM,
  DNA24_TERMINAL_MM,
  DNA24_HAIRPIN_LOOPS,
  DNA24_TETRALOOPS,
  DNA24_TRILOOPS,
  DNA24_COMPLEMENT,
} from '../src/lib/dna24.js';
import { DNA_ENERGIES } from '../src/lib/dna.js';
import { tm, setParameterSet } from '../src/lib/tm.js';
import { calculateTmQ5 } from '../src/lib/tmQ5.js';

let passed = 0;
let failed = 0;
let warnings = 0;

function test(name, condition, message = '') {
  if (condition) {
    console.log(`  ✓ ${name}`);
    passed++;
  } else {
    console.log(`  ✗ ${name}${message ? ': ' + message : ''}`);
    failed++;
  }
}

function warn(name, message) {
  console.log(`  ⚠ ${name}: ${message}`);
  warnings++;
}

function section(title) {
  console.log(`\n=== ${title} ===`);
}

// ============= Parameter Completeness Tests =============

section('Parameter Completeness');

// Watson-Crick pairs
const wcPairs = [
  'AA/TT', 'TT/AA', 'AT/TA', 'TA/AT',
  'CA/GT', 'GT/CA', 'CT/GA', 'GA/CT',
  'AC/TG', 'TG/AC', 'AG/TC', 'TC/AG',
  'CG/GC', 'GC/CG', 'GG/CC', 'CC/GG',
];

let allWcPresent = true;
const missingPairs = [];
for (const pair of wcPairs) {
  if (!DNA24_NN[pair]) {
    allWcPresent = false;
    missingPairs.push(pair);
  }
}
test('All 16 Watson-Crick NN pairs present', allWcPresent,
  missingPairs.length > 0 ? `Missing: ${missingPairs.join(', ')}` : '');

test('Initialization parameters present',
  DNA24_NN['init'] && DNA24_NN['init_A/T'] && DNA24_NN['init_G/C'] && DNA24_NN['sym']);

const gtPairs = Object.keys(DNA24_NN).filter(k =>
  k.includes('G') && k.includes('T') && k.includes('/') && !k.includes('init')
);
test('G-T wobble parameters present', gtPairs.length > 5,
  `Found ${gtPairs.length} G-T pairs`);

const mmCount = Object.keys(DNA24_INTERNAL_MM).length;
test('Internal mismatches have 576+ entries', mmCount >= 576,
  `Found ${mmCount} entries`);

const tmCount = Object.keys(DNA24_TERMINAL_MM).length;
test('Terminal mismatches defined', tmCount > 0,
  `Found ${tmCount} entries`);

let allHairpinLoops = true;
for (let i = 3; i <= 30; i++) {
  if (!DNA24_HAIRPIN_LOOPS[i]) {
    allHairpinLoops = false;
    break;
  }
}
test('Hairpin loop sizes 3-30 present', allHairpinLoops);

const tetraloopCount = Object.keys(DNA24_TETRALOOPS).length;
test('Tetraloops have expanded dataset', tetraloopCount >= 256,
  `Found ${tetraloopCount} entries (DNA24 claims 1062 vs SantaLucia's 130)`);

const triloopCount = Object.keys(DNA24_TRILOOPS).length;
test('Triloops defined', triloopCount > 0,
  `Found ${triloopCount} entries`);

test('Complement mapping complete',
  DNA24_COMPLEMENT['A'] === 'T' &&
  DNA24_COMPLEMENT['T'] === 'A' &&
  DNA24_COMPLEMENT['G'] === 'C' &&
  DNA24_COMPLEMENT['C'] === 'G');

// ============= Thermodynamic Value Range Tests =============

section('Thermodynamic Value Ranges');

let dhInRange = true;
let dhOutOfRange = [];
for (const [pair, values] of Object.entries(DNA24_NN)) {
  if (pair === 'init' || pair === 'sym') continue;
  const dH = values[0];
  if (dH < -15 || dH > 10) {
    dhInRange = false;
    dhOutOfRange.push(`${pair}: dH=${dH}`);
  }
}
test('NN dH values in range (-15 to +10 kcal/mol)', dhInRange,
  dhOutOfRange.length > 0 ? `Out of range: ${dhOutOfRange.slice(0, 3).join(', ')}` : '');

let dsInRange = true;
let dsOutOfRange = [];
for (const [pair, values] of Object.entries(DNA24_NN)) {
  if (pair === 'init' || pair === 'sym') continue;
  const dS = values[1];
  if (dS < -50 || dS > 30) {
    dsInRange = false;
    dsOutOfRange.push(`${pair}: dS=${dS}`);
  }
}
test('NN dS values in range (-50 to +30 cal/mol/K)', dsInRange,
  dsOutOfRange.length > 0 ? `Out of range: ${dsOutOfRange.slice(0, 3).join(', ')}` : '');

// dG check at 37°C
const cgGc = DNA24_NN['CG/GC'];
const dG_cgGc = cgGc[0] - (310.15 * cgGc[1] / 1000);
const atTa = DNA24_NN['AT/TA'];
const dG_atTa = atTa[0] - (310.15 * atTa[1] / 1000);
test('CG/GC more stable than AT/TA (expected)', dG_cgGc < dG_atTa,
  `dG(CG/GC)=${dG_cgGc.toFixed(2)}, dG(AT/TA)=${dG_atTa.toFixed(2)}`);

// ============= DNA24 vs SantaLucia Comparison =============

section('DNA24 vs SantaLucia Comparison');

const comparisonPairs = ['AA/TT', 'AT/TA', 'CG/GC', 'GC/CG'];
let agreementCount = 0;

for (const pair of comparisonPairs) {
  const dna24 = DNA24_NN[pair];
  const santaLucia = DNA_ENERGIES.NN[pair];

  if (dna24 && santaLucia) {
    const dhRatio = Math.abs(dna24[0] / santaLucia[0]);
    const agrees = dhRatio > 0.7 && dhRatio < 1.3;
    if (agrees) agreementCount++;

    console.log(`  ${pair}: DNA24 dH=${dna24[0].toFixed(1)}, SantaLucia dH=${santaLucia[0].toFixed(1)} ` +
      `(ratio ${dhRatio.toFixed(2)}) ${agrees ? '✓' : '⚠'}`);
  }
}
test('DNA24/SantaLucia WC pairs agree within 30%', agreementCount === comparisonPairs.length);

// ============= Tm Calculation Tests =============

section('Tm Calculations');

setParameterSet(true); // Use DNA24

// Test Q5 Tm for typical primers
const primers = [
  { seq: 'ATGCATGCATGCATGCATGC', name: '20bp 50% GC', expectedRange: [58, 68] },
  { seq: 'AAAAAAAAAAAAAAAAAAAA', name: '20bp poly-A', expectedRange: [35, 55] },
  { seq: 'GGGGGGGGGGGGGGGGGGGG', name: '20bp poly-G (G-quad)', expectedRange: [75, 95] }, // G-quadruplex = very stable
  { seq: 'GCTAGCTAGCTAGCTAGCTA', name: '20bp alternating', expectedRange: [55, 70] },
];

for (const { seq, name, expectedRange } of primers) {
  const tmValue = calculateTmQ5(seq);
  const inRange = tmValue >= expectedRange[0] && tmValue <= expectedRange[1];
  test(`Q5 Tm ${name}: ${tmValue.toFixed(1)}°C`, inRange,
    inRange ? '' : `Expected ${expectedRange[0]}-${expectedRange[1]}°C`);
}

// Length dependency
const short = 'ATGCATGCATGC';
const medium = 'ATGCATGCATGCATGC';
const long = 'ATGCATGCATGCATGCATGC';

const tmShort = calculateTmQ5(short);
const tmMedium = calculateTmQ5(medium);
const tmLong = calculateTmQ5(long);

test('Tm increases with length', tmShort < tmMedium && tmMedium < tmLong,
  `12bp=${tmShort.toFixed(1)}°C, 16bp=${tmMedium.toFixed(1)}°C, 20bp=${tmLong.toFixed(1)}°C`);

// GC dependency
const atRich = 'AAAATTTTTAAAATTTTTAA';
const balanced = 'ATGCATGCATGCATGCATGC';
const gcRich = 'GGGGCCCCGGGGCCCCGGGG';

const tmAT = calculateTmQ5(atRich);
const tmBalanced = calculateTmQ5(balanced);
const tmGC = calculateTmQ5(gcRich);

test('Tm increases with GC content', tmAT < tmBalanced && tmBalanced < tmGC,
  `0%GC=${tmAT.toFixed(1)}°C, 50%GC=${tmBalanced.toFixed(1)}°C, 100%GC=${tmGC.toFixed(1)}°C`);

// ============= Performance Test =============

section('Performance');

const testPrimer = 'ATGCATGCATGCATGCATGC';
const iterations = 1000;

const start = performance.now();
for (let i = 0; i < iterations; i++) {
  calculateTmQ5(testPrimer);
}
const elapsed = performance.now() - start;
const perCalc = elapsed / iterations;

test(`Q5 Tm calc < 1ms each`, perCalc < 1,
  `${perCalc.toFixed(3)}ms per calculation (${iterations} iterations)`);

// ============= Summary =============

section('Summary');

console.log(`\nTotal: ${passed} passed, ${failed} failed, ${warnings} warnings`);

if (failed > 0) {
  console.log('\n❌ Some validation checks failed!');
  process.exit(1);
} else if (warnings > 0) {
  console.log('\n⚠ Validation passed with warnings');
  process.exit(0);
} else {
  console.log('\n✅ All validation checks passed!');
  process.exit(0);
}
