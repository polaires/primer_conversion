/**
 * Calculate the tm of a DNA sequence
 * Ported from seqfold Python library
 *
 * Updated to support DNA24 parameters (Greenleaf Lab, 2024)
 * which provide ~50% better accuracy for mismatches and G-T wobbles.
 */

import { DNA_COMPLEMENT, DNA_ENERGIES } from "./dna.js";
import { DNA24_COMPLEMENT, DNA24_ENERGIES } from "./dna24.js";

// Default to DNA24 parameters for better accuracy
let USE_DNA24 = true;

/**
 * Set which parameter set to use for Tm calculations
 * @param {boolean} useDna24 - true for DNA24 (recommended), false for SantaLucia 1998
 */
export function setParameterSet(useDna24) {
  USE_DNA24 = useDna24;
}

/**
 * Get the current energy parameters based on selected parameter set
 */
function getEnergyParams() {
  return USE_DNA24 ? DNA24_ENERGIES : DNA_ENERGIES;
}

/**
 * Get the complement map based on selected parameter set
 */
function getComplement() {
  return USE_DNA24 ? DNA24_COMPLEMENT : DNA_COMPLEMENT;
}

/**
 * Calculate the annealing temperature between seq1 and seq2.
 *
 * If seq2 is not provided, its exact complement is used.
 * In otherwords, it's assumed to be an exact match.
 *
 * This is largely influenced by Bio.SeqUtils.MeltingTemp with
 * some different defaults. Here, the reaction mixture is assumed to
 * be PCR and concentrations for Mg, Tris, K, and dNTPs are included.
 *
 * @param {string} seq1 - The seq whose tm is calculated
 * @param {string} seq2 - The seq that seq1 anneals to in 3' -> 5' direction
 * @param {boolean} pcr - Whether tm is being calculated for PCR
 * @returns {number} The estimated tm as a float
 */
export function tm(seq1, seq2 = "", pcr = true) {
  [seq1, seq2] = parseInput(seq1, seq2);
  const energies = getEnergyParams();

  // sum enthalpy and entropy. Enthalpy is first value of each tuple and
  // entropy is the second value of each tuple in:
  // SantaLucia & Hicks (2004) or DNA24 (Greenleaf Lab 2024)

  // start with initiation enthalpy and entropy
  let [dh, ds] = energies.NN["init"];

  // add in initial A/T and initial G/Cs
  const init = seq1[0] + seq1[seq1.length - 1];
  const initAt =
    (init.match(/A/g) || []).length + (init.match(/T/g) || []).length;
  const initGc =
    (init.match(/G/g) || []).length + (init.match(/C/g) || []).length;
  const [initAtH, initAtS] = energies.NN["init_A/T"];
  const [initGcH, initGcS] = energies.NN["init_G/C"];
  dh += initAt * initAtH + initGc * initGcH;
  ds += initAt * initAtS + initGc * initGcS;

  // work through each nearest neighbor pair
  for (let i = 0; i < seq1.length - 1; i++) {
    const pair = `${seq1[i]}${seq1[i + 1]}/${seq2[i]}${seq2[i + 1]}`;

    // assuming internal neighbor pair
    let pairDh = 0.0,
      pairDs = 0.0;
    if (pair in energies.NN) {
      [pairDh, pairDs] = energies.NN[pair];
    } else if (energies.INTERNAL_MM && pair in energies.INTERNAL_MM) {
      // For DNA24, INTERNAL_MM uses 6-letter context codes for folding
      // Fall back to simple mismatch penalty for Tm calculation
      [pairDh, pairDs] = [0, 0]; // Mismatch destabilizes
    }

    // overwrite if it's a terminal pair
    if (i === 0 || i === seq1.length - 2) {
      if (energies.TERMINAL_MM && pair in energies.TERMINAL_MM) {
        [pairDh, pairDs] = energies.TERMINAL_MM[pair];
      }
    }

    dh += pairDh;
    ds += pairDs;
  }

  const gc = getGc(seq1);
  return calcTm(dh, ds, pcr, gc, seq1.length);
}

/**
 * Return a TmCache where each (i, j) returns the Tm for that subspan.
 *
 * @param {string} seq1 - The seq whose tm is calculated
 * @param {string} seq2 - The seq that seq1 anneals to
 * @param {boolean} pcr - Whether for PCR
 * @returns {number[][]} 2D array where [i][j] returns tm for that range
 */
export function tmCache(seq1, seq2 = "", pcr = true) {
  [seq1, seq2] = parseInput(seq1, seq2);
  const energies = getEnergyParams();
  const n = seq1.length;

  const arrGc = gcCache(seq1);
  const arrDh = [];
  const arrDs = [];
  const arrTm = [];

  for (let i = 0; i < n; i++) {
    arrDh.push(new Array(n).fill(0.0));
    arrDs.push(new Array(n).fill(0.0));
    arrTm.push(new Array(n).fill(Infinity));
  }

  // fill in the diagonal
  for (let i = 0; i < n; i++) {
    if (i === n - 1) {
      // hackish
      arrDh[i][i] = arrDh[i - 1][i - 1];
      arrDs[i][i] = arrDs[i - 1][i - 1];
      continue;
    }

    const pair = `${seq1[i]}${seq1[i + 1]}/${seq2[i]}${seq2[i + 1]}`;
    let dhVal, dsVal;
    if (pair in energies.NN) {
      [dhVal, dsVal] = energies.NN[pair];
    } else if (energies.INTERNAL_MM && pair in energies.INTERNAL_MM) {
      [dhVal, dsVal] = energies.INTERNAL_MM[pair];
    } else {
      // Fallback for unknown pairs
      [dhVal, dsVal] = [0, 0];
    }

    arrDh[i][i] = dhVal;
    arrDs[i][i] = dsVal;
  }

  // fill in the tm array
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      arrDh[i][j] = arrDh[i][j - 1] + arrDh[j][j];
      arrDs[i][j] = arrDs[i][j - 1] + arrDs[j][j];
      arrTm[i][j] = calcTm(arrDh[i][j], arrDs[i][j], pcr, arrGc[i][j], j - i + 1);
    }
  }

  return arrTm;
}

/**
 * Return the GC ratio of each range, between i and j, in the sequence
 *
 * @param {string} seq - The sequence whose GC ratio we're querying
 * @returns {number[][]} A cache for GC ratio lookup
 */
export function gcCache(seq) {
  const n = seq.length;
  const arrGc = [];

  for (let i = 0; i < n; i++) {
    arrGc.push(new Array(n).fill(Infinity));
  }

  // fill in the diagonal
  for (let i = 0; i < n; i++) {
    if (i === n - 1) {
      // hackish
      arrGc[i][i] = arrGc[i - 1][i - 1];
      continue;
    }

    arrGc[i][i] = seq[i] === "G" || seq[i] === "C" ? 1.0 : 0.0;

    if (i === n - 2 && !arrGc[i][i]) {
      // don't ignore last pair
      arrGc[i][i] = seq[i + 1] === "G" || seq[i + 1] === "C" ? 1.0 : 0.0;
    }
  }

  // fill in the upper right of the array
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      arrGc[i][j] = arrGc[i][j - 1] + arrGc[j][j];
    }
  }

  // convert to ratios
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      arrGc[i][j] = Math.round((arrGc[i][j] / (j - i + 1)) * 10) / 10;
    }
  }

  return arrGc;
}

/**
 * Parse and prepare the input sequences. Throw if there's an issue.
 *
 * @param {string} seq1 - The main sequence whose tm is being calculated
 * @param {string} seq2 - The second sequence that seq1 is annealing to
 * @returns {[string, string]} The sequences to use for tm calculation
 */
function parseInput(seq1, seq2 = "") {
  seq1 = seq1.toUpperCase();
  if (!seq2) {
    const complement = getComplement();
    seq2 = seq1
      .split("")
      .map((c) => complement[c])
      .join("");
  }

  if (seq1.length !== seq2.length) {
    throw new Error(
      `Length mismatch between seq1 ${seq1.length} and seq2 ${seq2.length}`
    );
  }

  if (seq1.length < 2) {
    throw new Error(
      `Sequence, ${seq1.length}bp, is too short for tm calculation`
    );
  }

  return [seq1, seq2];
}

/**
 * Apply the correction formula to estimate Tm
 *
 * @param {number} dh - Accumulated enthalpy
 * @param {number} ds - Accumulated entropy
 * @param {boolean} pcr - Whether this is for PCR or not
 * @param {number} gc - The GC% of the sequence
 * @param {number} seqLen - The length of the sequence
 * @returns {number} The estimated tm
 */
function calcTm(dh, ds, pcr, gc, seqLen) {
  // adjust salt based on mode
  let seq1Conc, seq2Conc, Na, K, Tris, Mg, dNTPs;
  if (pcr) {
    seq1Conc = 250.0;
    seq2Conc = 0.0;
    Na = 0;
    K = 50;
    Tris = 2; // see Thermo
    Mg = 1.5; // see NEB
    dNTPs = 0.2; // see NEB
  } else {
    seq1Conc = 25.0;
    seq2Conc = 25.0;
    Na = 50;
    K = 0;
    Tris = 0;
    Mg = 0;
    dNTPs = 0;
  }

  // salt correction for deltaS
  // copied from Bio.SeqUtils' use of a decision tree by:
  // Owczarzy et al. (2008), Biochemistry 47: 5336-5353
  const Mon = Na + K + Tris / 2.0; // monovalent ions
  let mg = Mg * 1e-3; // Lowercase ions (mg, mon, dntps) are molar
  const mon = Mon * 1e-3;

  // coefficients to a multi-variate from the paper
  let a = 3.92,
    b = -0.911,
    c = 6.26,
    d = 1.42,
    e = -48.2,
    f = 52.5,
    g = 8.31;

  if (dNTPs > 0) {
    const dntps = dNTPs * 1e-3;
    const ka = 3e4; // Dissociation constant for Mg:dNTP
    // Free Mg2+ calculation:
    mg =
      (-(ka * dntps - ka * mg + 1.0) +
        Math.sqrt(
          Math.pow(ka * dntps - ka * mg + 1.0, 2) + 4.0 * ka * mg
        )) /
      (2.0 * ka);
  }

  let corr;
  if (Mon > 0) {
    const R = Math.sqrt(mg) / mon;
    if (R < 0.22) {
      return (
        ((4.29 * gc) / 100 - 3.95) * 1e-5 * Math.log(mon) +
        9.4e-6 * Math.pow(Math.log(mon), 2)
      );
    } else if (R < 6.0) {
      a = 3.92 * (0.843 - 0.352 * Math.sqrt(mon) * Math.log(mon));
      d =
        1.42 *
        (1.279 - 4.03e-3 * Math.log(mon) - 8.03e-3 * Math.pow(Math.log(mon), 2));
      g =
        8.31 *
        (0.486 - 0.258 * Math.log(mon) + 5.25e-3 * Math.pow(Math.log(mon), 3));
    }
  }

  corr =
    (a +
      b * Math.log(mg) +
      (gc / 100) * (c + d * Math.log(mg)) +
      (1 / (2.0 * (seqLen - 1))) *
        (e + f * Math.log(mg) + g * Math.pow(Math.log(mg), 2))) *
    1e-5;

  // tm with concentration consideration
  const k = (seq1Conc - seq2Conc / 2.0) * 1e-9;
  const R_const = 1.9872;
  let est = (dh * 1000.0) / (ds + R_const * Math.log(k)) - 273.15;

  // add in salt correction
  est = 1 / (1 / (est + 273.15) + corr) - 273.1;

  return Math.round(est * 10) / 10;
}

/**
 * Return the GC ratio of a sequence.
 *
 * @param {string} seq - The sequence
 * @returns {number} The GC ratio
 */
function getGc(seq) {
  const gcCount =
    (seq.match(/G/g) || []).length + (seq.match(/C/g) || []).length;
  return gcCount / seq.length;
}
