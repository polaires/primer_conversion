/**
 * Extract ligation frequency data from NEB Excel files
 *
 * Based on Pryor et al. 2020 PLOS ONE paper:
 * "Enabling one-pot Golden Gate assemblies of unprecedented complexity"
 *
 * Data format: 256×256 matrix (4-base overhangs) or 64×64 (3-base for SapI)
 * where each cell contains the ligation frequency for that overhang pair.
 */

import XLSX from 'xlsx';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const paperDir = path.join(__dirname, '..', 'paper');

// Generate all possible 4-base overhangs (256 total)
function generate4BaseOverhangs() {
  const bases = ['A', 'C', 'G', 'T'];
  const overhangs = [];
  for (const b1 of bases) {
    for (const b2 of bases) {
      for (const b3 of bases) {
        for (const b4 of bases) {
          overhangs.push(b1 + b2 + b3 + b4);
        }
      }
    }
  }
  return overhangs;
}

// Generate all possible 3-base overhangs (64 total) for SapI
function generate3BaseOverhangs() {
  const bases = ['A', 'C', 'G', 'T'];
  const overhangs = [];
  for (const b1 of bases) {
    for (const b2 of bases) {
      for (const b3 of bases) {
        overhangs.push(b1 + b2 + b3);
      }
    }
  }
  return overhangs;
}

// Get reverse complement of an overhang
function reverseComplement(seq) {
  const comp = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' };
  return seq.split('').reverse().map(b => comp[b]).join('');
}

// Parse an Excel file containing ligation frequency data
function parseExcelFile(filePath, is3Base = false) {
  console.log(`\nParsing: ${path.basename(filePath)}`);

  const workbook = XLSX.read(fs.readFileSync(filePath));
  const sheetName = workbook.SheetNames[0];
  const sheet = workbook.Sheets[sheetName];

  // Convert to 2D array
  const data = XLSX.utils.sheet_to_json(sheet, { header: 1 });

  console.log(`  Sheet: ${sheetName}`);
  console.log(`  Rows: ${data.length}`);

  const expectedSize = is3Base ? 64 : 256;
  const allOverhangs = is3Base ? generate3BaseOverhangs() : generate4BaseOverhangs();
  const ohLength = is3Base ? 3 : 4;

  // Find header row (row with "Overhang" in first column and overhangs in remaining columns)
  let headerRowIdx = -1;
  let colHeaders = [];

  for (let r = 0; r < data.length; r++) {
    const row = data[r];
    if (!row) continue;

    // Check if this row contains overhang headers
    const firstCell = String(row[0] || '').trim();
    if (firstCell === 'Overhang' || firstCell === '') {
      // Check if subsequent cells are overhangs
      let ohCount = 0;
      for (let c = 1; c < row.length && c <= expectedSize + 1; c++) {
        const cell = String(row[c] || '').trim().toUpperCase();
        if (cell.length === ohLength && /^[ACGT]+$/.test(cell)) {
          ohCount++;
        }
      }
      if (ohCount >= expectedSize * 0.9) {
        headerRowIdx = r;
        // Extract column headers
        for (let c = 1; c <= expectedSize; c++) {
          const cell = String(row[c] || '').trim().toUpperCase();
          colHeaders.push(cell);
        }
        break;
      }
    }
  }

  if (headerRowIdx === -1) {
    throw new Error('Could not find header row');
  }

  console.log(`  Header row: ${headerRowIdx}`);
  console.log(`  Column headers found: ${colHeaders.length}`);
  console.log(`  Sample col headers: ${colHeaders.slice(0, 5).join(', ')}`);

  // Build full matrix indexed by overhang strings
  const matrix = {};

  // Initialize matrix for all overhangs
  for (const oh of allOverhangs) {
    matrix[oh] = {};
    for (const oh2 of allOverhangs) {
      matrix[oh][oh2] = 0;
    }
  }

  // Parse data rows (starting after header)
  let rowCount = 0;
  for (let r = headerRowIdx + 1; r < data.length; r++) {
    const row = data[r];
    if (!row || row.length < 2) continue;

    const rowHeader = String(row[0] || '').trim().toUpperCase();
    if (rowHeader.length !== ohLength || !/^[ACGT]+$/.test(rowHeader)) {
      continue;
    }

    rowCount++;

    // Parse each column value
    for (let c = 0; c < colHeaders.length; c++) {
      const colHeader = colHeaders[c];
      const value = row[c + 1]; // +1 because row[0] is the row header

      let freq = 0;
      if (typeof value === 'number') {
        freq = value;
      } else if (typeof value === 'string') {
        freq = parseFloat(value) || 0;
      }

      matrix[rowHeader][colHeader] = freq;
    }
  }

  console.log(`  Data rows parsed: ${rowCount}`);

  // Calculate statistics
  let totalEntries = 0;
  let nonZeroEntries = 0;
  let sum = 0;
  let max = 0;

  // Calculate per-overhang fidelity (correct / total)
  const overhangFidelity = {};

  for (const oh of allOverhangs) {
    const wcPartner = reverseComplement(oh);
    const correctFreq = matrix[oh][wcPartner] || 0;

    let totalFreq = 0;
    for (const oh2 of allOverhangs) {
      const freq = matrix[oh][oh2] || 0;
      totalFreq += freq;
      totalEntries++;
      if (freq > 0) {
        nonZeroEntries++;
        sum += freq;
        if (freq > max) max = freq;
      }
    }

    overhangFidelity[oh] = totalFreq > 0 ? correctFreq / totalFreq : 0;
  }

  console.log(`  Total matrix entries: ${totalEntries}`);
  console.log(`  Non-zero entries: ${nonZeroEntries}`);
  console.log(`  Max frequency: ${max}`);

  // Find best and worst overhangs
  const fidelities = Object.entries(overhangFidelity)
    .filter(([oh, f]) => f > 0)
    .sort((a, b) => b[1] - a[1]);

  console.log(`  Overhangs with data: ${fidelities.length}`);
  console.log(`  Top 5 fidelity: ${fidelities.slice(0, 5).map(([oh, f]) => `${oh}:${(f*100).toFixed(0)}%`).join(', ')}`);
  console.log(`  Bottom 5 fidelity: ${fidelities.slice(-5).map(([oh, f]) => `${oh}:${(f*100).toFixed(0)}%`).join(', ')}`);

  return {
    overhangs: allOverhangs,
    matrix,
    overhangFidelity,
    stats: {
      totalEntries,
      nonZeroEntries,
      maxFrequency: max
    }
  };
}

// Calculate assembly fidelity for a given set of overhangs
function calculateAssemblyFidelity(matrix, overhangs, allOverhangs) {
  // For each junction (overhang), calculate p(correct) = correct / total
  // Assembly fidelity = product of all p(correct)

  let assemblyFidelity = 1.0;
  const junctionFidelities = [];

  for (const oh of overhangs) {
    const wcPartner = reverseComplement(oh);
    const correctFreq = matrix[oh]?.[wcPartner] || 0;

    // Total includes only the overhangs in our set (not all 256)
    // This is the key insight from the paper - we only count mismatch potential
    // with other overhangs that are actually present in the assembly
    let totalFreq = correctFreq;

    for (const otherOh of overhangs) {
      if (otherOh === oh) continue;
      // Add mismatch frequency with other overhangs' complements
      const otherWc = reverseComplement(otherOh);
      totalFreq += matrix[oh]?.[otherWc] || 0;
    }

    const junctionFidelity = totalFreq > 0 ? correctFreq / totalFreq : 1.0;
    junctionFidelities.push({ overhang: oh, fidelity: junctionFidelity });
    assemblyFidelity *= junctionFidelity;
  }

  return {
    assemblyFidelity,
    junctionFidelities: junctionFidelities.sort((a, b) => a.fidelity - b.fidelity)
  };
}

// Find optimal overhang sets using greedy algorithm
function findOptimalSets(matrix, allOverhangs, numParts) {
  const numJunctions = numParts + 1; // Including vector junctions

  // Start with highest-fidelity overhangs
  const candidates = allOverhangs.filter(oh => {
    const wc = reverseComplement(oh);
    return (matrix[oh]?.[wc] || 0) > 0;
  });

  // Score each overhang by its individual fidelity
  const scored = candidates.map(oh => {
    const wc = reverseComplement(oh);
    const correct = matrix[oh][wc] || 0;
    let total = 0;
    for (const oh2 of allOverhangs) {
      total += matrix[oh][oh2] || 0;
    }
    return { overhang: oh, fidelity: total > 0 ? correct / total : 0 };
  }).sort((a, b) => b.fidelity - a.fidelity);

  // Greedy selection: pick overhangs that minimize cross-reactivity
  const selected = [];
  const used = new Set();

  for (const { overhang } of scored) {
    if (selected.length >= numJunctions) break;

    // Check compatibility with already selected overhangs
    let compatible = true;
    const wc = reverseComplement(overhang);

    // Skip if overhang or its complement already used
    if (used.has(overhang) || used.has(wc)) continue;

    // Check cross-reactivity with selected overhangs
    for (const sel of selected) {
      const selWc = reverseComplement(sel);
      // Check if this overhang would ligate with selected's complement
      const crossReact1 = matrix[overhang]?.[selWc] || 0;
      const crossReact2 = matrix[sel]?.[wc] || 0;
      const correct1 = matrix[overhang]?.[wc] || 0;
      const correct2 = matrix[sel]?.[selWc] || 0;

      // If cross-reactivity is >5% of correct ligation, skip
      if (crossReact1 > correct1 * 0.05 || crossReact2 > correct2 * 0.05) {
        compatible = false;
        break;
      }
    }

    if (compatible) {
      selected.push(overhang);
      used.add(overhang);
      used.add(wc);
    }
  }

  // Calculate final assembly fidelity
  const result = calculateAssemblyFidelity(matrix, selected, allOverhangs);

  return {
    overhangs: selected,
    ...result
  };
}

// Main extraction function
async function extractAllData() {
  const files = [
    { name: 'BsaI-HFv2', path: 'Ligation frequency for each overhang pair in assembly reactions with BsaI-HFv2 and T4 DNA ligase..xlsx', is3Base: false },
    { name: 'BsmBI-v2', path: 'Ligation frequency for each overhang pair in assembly reactions with BsmBI-v2 and T4 DNA ligase.xlsx', is3Base: false },
    { name: 'Esp3I', path: 'Ligation frequency for each overhang pair in assembly reactions with Esp3I and T4 DNA ligase.xlsx', is3Base: false },
    { name: 'BbsI-HF', path: 'Ligation frequency for each overhang pair in assembly reactions with BbsI-HF and T4 DNA ligase.xlsx', is3Base: false },
    { name: 'SapI', path: 'Ligation frequency for each overhang pair in assembly reactions with SapI and T4 DNA ligase.xlsx', is3Base: true },
  ];

  const results = {};

  for (const file of files) {
    const fullPath = path.join(paperDir, file.path);
    if (fs.existsSync(fullPath)) {
      try {
        results[file.name] = parseExcelFile(fullPath, file.is3Base);
      } catch (err) {
        console.error(`  Error parsing ${file.name}: ${err.message}`);
      }
    } else {
      console.error(`  File not found: ${file.path}`);
    }
  }

  return results;
}

// Generate compact data structure for the app
function generateAppData(results) {
  const output = {
    metadata: {
      source: 'Pryor et al. 2020 PLOS ONE - Data-optimized Assembly Design',
      doi: '10.1371/journal.pone.0238592',
      generated: new Date().toISOString(),
      description: 'Ligation frequency matrices for Golden Gate assembly overhang design'
    },
    enzymes: {}
  };

  for (const [enzyme, data] of Object.entries(results)) {
    if (!data || !data.matrix) continue;

    // Store compact matrix (only non-zero values)
    const compactMatrix = {};
    for (const oh1 of data.overhangs) {
      const row = {};
      for (const oh2 of data.overhangs) {
        const val = data.matrix[oh1][oh2];
        if (val > 0) {
          row[oh2] = val;
        }
      }
      if (Object.keys(row).length > 0) {
        compactMatrix[oh1] = row;
      }
    }

    // Pre-calculate optimal sets for common assembly sizes
    console.log(`\nGenerating optimal sets for ${enzyme}...`);
    const optimalSets = {};
    for (const numParts of [2, 3, 4, 5, 6, 8, 10, 12]) {
      const set = findOptimalSets(data.matrix, data.overhangs, numParts);
      optimalSets[numParts] = {
        overhangs: set.overhangs,
        fidelity: set.assemblyFidelity,
        lowestJunction: set.junctionFidelities[0]
      };
      console.log(`  ${numParts} parts: ${(set.assemblyFidelity * 100).toFixed(1)}% fidelity`);
    }

    output.enzymes[enzyme] = {
      overhangLength: enzyme === 'SapI' ? 3 : 4,
      overhangs: data.overhangs,
      matrix: compactMatrix,
      overhangFidelity: data.overhangFidelity,
      optimalSets,
      stats: data.stats
    };
  }

  return output;
}

// Run extraction
console.log('=== Extracting Ligation Frequency Data ===');
console.log('Source: Pryor et al. 2020 PLOS ONE');

extractAllData().then(results => {
  console.log('\n=== Generating App Data ===');
  const appData = generateAppData(results);

  // Write full data to JSON file
  const outputPath = path.join(__dirname, '..', 'src', 'lib', 'repp', 'ligation-data.json');
  fs.writeFileSync(outputPath, JSON.stringify(appData, null, 2));
  console.log(`\nWritten to: ${outputPath}`);

  // Print final summary
  console.log('\n=== Final Summary ===');
  for (const [enzyme, data] of Object.entries(appData.enzymes)) {
    console.log(`\n${enzyme}:`);
    console.log(`  Overhang length: ${data.overhangLength} bp`);
    console.log(`  Matrix size: ${Object.keys(data.matrix).length} overhangs with data`);

    // Show optimal sets
    console.log('  Optimal assembly fidelities:');
    for (const [parts, set] of Object.entries(data.optimalSets)) {
      console.log(`    ${parts} parts: ${(set.fidelity * 100).toFixed(1)}% (${set.overhangs.length} junctions)`);
    }
  }
});
