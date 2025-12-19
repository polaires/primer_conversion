/**
 * Validate DNA24 parameters against array melt experimental data
 *
 * Uses arr_v1_1M_n=27732.csv which contains measured dG and Tm values
 * for 27,732 hairpin sequences at 1M Na+
 */

import { createReadStream } from "fs";
import { createInterface } from "readline";
import { dg, setFoldParameterSet } from "./fold.js";
import { setParameterSet } from "./tm.js";

/**
 * Parse CSV and return array of objects
 */
async function parseCSV(filepath) {
  const results = [];
  const rl = createInterface({
    input: createReadStream(filepath),
    crlfDelay: Infinity,
  });

  let headers = null;
  for await (const line of rl) {
    if (!headers) {
      headers = line.split(",");
      continue;
    }
    const values = line.split(",");
    const row = {};
    headers.forEach((h, i) => {
      row[h] = values[i];
    });
    results.push(row);
  }
  return results;
}

/**
 * Calculate RMSE between two arrays
 */
function rmse(predicted, actual) {
  if (predicted.length !== actual.length) {
    throw new Error("Arrays must have same length");
  }
  let sumSq = 0;
  for (let i = 0; i < predicted.length; i++) {
    const diff = predicted[i] - actual[i];
    sumSq += diff * diff;
  }
  return Math.sqrt(sumSq / predicted.length);
}

/**
 * Calculate Pearson correlation coefficient
 */
function correlation(x, y) {
  const n = x.length;
  const meanX = x.reduce((a, b) => a + b, 0) / n;
  const meanY = y.reduce((a, b) => a + b, 0) / n;

  let num = 0, denX = 0, denY = 0;
  for (let i = 0; i < n; i++) {
    const dx = x[i] - meanX;
    const dy = y[i] - meanY;
    num += dx * dy;
    denX += dx * dx;
    denY += dy * dy;
  }
  return num / Math.sqrt(denX * denY);
}

/**
 * Run validation
 */
async function validate() {
  console.log("Loading experimental data...");
  const data = await parseCSV(new URL("./arr_v1_1M_n=27732.csv", import.meta.url).pathname);
  console.log(`Loaded ${data.length} sequences\n`);

  // Filter to valid entries with sequence and structure
  const hairpins = data.filter(row => {
    const seq = row.RefSeq;
    const struct = row.TargetStruct;
    const dg = parseFloat(row.dG_37);
    // Must have sequence, structure, and valid dG
    return seq && seq.length >= 6 && struct && struct.includes("(") && !isNaN(dg);
  });

  console.log(`Testing on ${Math.min(1000, hairpins.length)} hairpin sequences...\n`);

  // Test both parameter sets
  for (const useDna24 of [true, false]) {
    const paramName = useDna24 ? "DNA24" : "SantaLucia 1998";
    setFoldParameterSet(useDna24);
    setParameterSet(useDna24);

    const predicted = [];
    const actual = [];
    const errors = [];
    let skipped = 0;

    // Test full dataset
    const testSet = hairpins;

    for (const row of testSet) {
      const seq = row.RefSeq;
      const measuredDG = parseFloat(row.dG_37);

      if (!seq || isNaN(measuredDG)) {
        skipped++;
        continue;
      }

      try {
        const predictedDG = dg(seq, 37.0);

        // Skip if prediction is unreasonable
        if (predictedDG === Infinity || predictedDG === -Infinity || isNaN(predictedDG)) {
          skipped++;
          continue;
        }

        predicted.push(predictedDG);
        actual.push(measuredDG);
        errors.push(Math.abs(predictedDG - measuredDG));
      } catch (e) {
        skipped++;
      }
    }

    if (predicted.length === 0) {
      console.log(`${paramName}: No valid predictions\n`);
      continue;
    }

    const r = correlation(predicted, actual);
    const rmseVal = rmse(predicted, actual);
    const mae = errors.reduce((a, b) => a + b, 0) / errors.length;

    console.log(`=== ${paramName} ===`);
    console.log(`Sequences tested: ${predicted.length} (skipped: ${skipped})`);
    console.log(`Correlation (R): ${r.toFixed(4)}`);
    console.log(`RMSE: ${rmseVal.toFixed(3)} kcal/mol`);
    console.log(`MAE: ${mae.toFixed(3)} kcal/mol`);
    console.log();

    // Show some examples
    console.log("Sample predictions vs measured:");
    for (let i = 0; i < 5 && i < predicted.length; i++) {
      const row = testSet[i];
      console.log(`  ${row.RefSeq.substring(0, 20).padEnd(20)} Predicted: ${predicted[i].toFixed(2).padStart(6)} Measured: ${actual[i].toFixed(2).padStart(6)}`);
    }
    console.log();
  }
}

validate().catch(console.error);
