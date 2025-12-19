#!/usr/bin/env node
/**
 * Extend pre-computed optimal overhang sets to 15, 20, 25, 30 parts
 *
 * The existing ligation-data.json has optimal sets up to 12 parts.
 * This script extends them to support larger assemblies.
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const dataPath = path.join(__dirname, '..', 'src', 'lib', 'repp', 'ligation-data.json');

// Load existing data
const data = JSON.parse(fs.readFileSync(dataPath, 'utf-8'));

// Get reverse complement
function reverseComplement(seq) {
  const comp = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' };
  return seq.split('').reverse().map(b => comp[b]).join('');
}

// Calculate assembly fidelity for a set of overhangs
function calculateFidelity(overhangs, matrix) {
  let fidelity = 1.0;

  for (const oh of overhangs) {
    const ohUpper = oh.toUpperCase();
    const rc = reverseComplement(ohUpper);

    // Correct ligation
    const correct = (matrix[ohUpper]?.[rc] || 0) + (matrix[rc]?.[ohUpper] || 0);

    // Total possible ligations in this set
    let total = 0;
    for (const oh2 of overhangs) {
      const oh2Upper = oh2.toUpperCase();
      const rc2 = reverseComplement(oh2Upper);
      total += matrix[ohUpper]?.[oh2Upper] || 0;
      total += matrix[ohUpper]?.[rc2] || 0;
      total += matrix[rc]?.[oh2Upper] || 0;
      total += matrix[rc]?.[rc2] || 0;
    }

    if (total > 0) {
      fidelity *= correct / total;
    }
  }

  return fidelity;
}

// Greedy algorithm to find good overhang sets
function findOptimalSet(numJunctions, matrix, allOverhangs) {
  console.log(`  Finding optimal set for ${numJunctions} junctions...`);

  // Filter out palindromic overhangs
  const validOverhangs = allOverhangs.filter(oh => oh !== reverseComplement(oh));

  // Score each overhang by individual fidelity
  const scored = validOverhangs.map(oh => {
    const rc = reverseComplement(oh);
    const correct = (matrix[oh]?.[rc] || 0) + (matrix[rc]?.[oh] || 0);
    let total = 0;
    for (const oh2 of allOverhangs) {
      total += matrix[oh]?.[oh2] || 0;
    }
    return { overhang: oh, individualFidelity: total > 0 ? correct / total : 0 };
  }).filter(x => x.individualFidelity > 0)
    .sort((a, b) => b.individualFidelity - a.individualFidelity);

  // Greedy selection with cross-reactivity checking
  const selected = [];
  const used = new Set();

  for (const { overhang } of scored) {
    if (selected.length >= numJunctions) break;

    const rc = reverseComplement(overhang);
    if (used.has(overhang) || used.has(rc)) continue;

    // Check cross-reactivity with already selected
    let compatible = true;
    for (const sel of selected) {
      const selRc = reverseComplement(sel);
      const cross1 = matrix[overhang]?.[selRc] || 0;
      const cross2 = matrix[sel]?.[rc] || 0;
      const correct1 = matrix[overhang]?.[rc] || 0;
      const correct2 = matrix[sel]?.[selRc] || 0;

      // Skip if cross-reactivity > 5% of correct
      if ((correct1 > 0 && cross1 > correct1 * 0.05) ||
          (correct2 > 0 && cross2 > correct2 * 0.05)) {
        compatible = false;
        break;
      }
    }

    if (compatible) {
      selected.push(overhang);
      used.add(overhang);
      used.add(rc);
    }
  }

  // Calculate final fidelity
  const fidelity = calculateFidelity(selected, matrix);

  // Find lowest junction
  let lowestFidelity = 1.0;
  let lowestOH = selected[0];

  for (const oh of selected) {
    const rc = reverseComplement(oh);
    const correct = (matrix[oh]?.[rc] || 0) + (matrix[rc]?.[oh] || 0);

    let total = 0;
    for (const oh2 of selected) {
      const oh2Rc = reverseComplement(oh2);
      total += matrix[oh]?.[oh2] || 0;
      total += matrix[oh]?.[oh2Rc] || 0;
      total += matrix[rc]?.[oh2] || 0;
      total += matrix[rc]?.[oh2Rc] || 0;
    }

    const junctionFidelity = total > 0 ? correct / total : 0;
    if (junctionFidelity < lowestFidelity) {
      lowestFidelity = junctionFidelity;
      lowestOH = oh;
    }
  }

  console.log(`    Found ${selected.length} overhangs, fidelity: ${(fidelity * 100).toFixed(2)}%`);

  return {
    overhangs: selected,
    fidelity,
    lowestJunction: { overhang: lowestOH, fidelity: lowestFidelity }
  };
}

// Process each enzyme
console.log('=== Extending Optimal Overhang Sets ===\n');

const targetSizes = [15, 20, 25, 30];

for (const [enzymeName, enzymeData] of Object.entries(data.enzymes)) {
  console.log(`Processing ${enzymeName}...`);

  const matrix = enzymeData.matrix;
  const allOverhangs = enzymeData.overhangs;

  if (!enzymeData.optimalSets) {
    enzymeData.optimalSets = {};
  }

  for (const size of targetSizes) {
    // Skip if already exists
    if (enzymeData.optimalSets[size]) {
      console.log(`  ${size} parts: already exists, skipping`);
      continue;
    }

    const result = findOptimalSet(size, matrix, allOverhangs);

    if (result.overhangs.length >= size) {
      enzymeData.optimalSets[size] = {
        overhangs: result.overhangs.slice(0, size),
        fidelity: result.fidelity,
        lowestJunction: result.lowestJunction
      };
    } else {
      console.log(`    Warning: Only found ${result.overhangs.length} compatible overhangs`);
      enzymeData.optimalSets[size] = {
        overhangs: result.overhangs,
        fidelity: result.fidelity,
        lowestJunction: result.lowestJunction,
        warning: `Only ${result.overhangs.length} compatible overhangs found`
      };
    }
  }

  console.log('');
}

// Update metadata
data.metadata.extended = new Date().toISOString();
data.metadata.optimalSetSizes = Object.keys(data.enzymes['BsaI-HFv2'].optimalSets).map(Number).sort((a, b) => a - b);

// Write updated data
fs.writeFileSync(dataPath, JSON.stringify(data, null, 2));
console.log(`\nUpdated: ${dataPath}`);

// Print summary
console.log('\n=== Summary ===');
for (const [enzymeName, enzymeData] of Object.entries(data.enzymes)) {
  console.log(`\n${enzymeName}:`);
  const sizes = Object.keys(enzymeData.optimalSets).map(Number).sort((a, b) => a - b);
  for (const size of sizes) {
    const set = enzymeData.optimalSets[size];
    console.log(`  ${size} parts: ${(set.fidelity * 100).toFixed(1)}% fidelity`);
  }
}
