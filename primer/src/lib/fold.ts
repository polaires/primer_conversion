/**
 * Predict nucleic acid secondary structure
 * Ported from seqfold Python library
 *
 * Updated to support DNA24 parameters (Greenleaf Lab, 2024)
 * which provide improved accuracy for hairpins and mismatches.
 */

import { DNA_ENERGIES } from "./dna";
import { DNA24_ENERGIES } from "./dna24.js";

// Structure representing a folding structure with energy, description, and base pairs
export interface FoldStructure {
  e: number;          // Energy value
  desc: string;       // Description of the structure type
  ij: number[][];     // Array of [i, j] base pair indices
}

// Energy parameter set containing thermodynamic parameters
export interface EnergyParameters {
  COMPLEMENT: Record<string, string>;
  NN: Record<string, [number, number]>;
  DE: Record<string, [number, number]>;
  TERMINAL_MM?: Record<string, [number, number]>;
  HAIRPIN_MM?: Record<string, [number, number]>;
  INTERNAL_MM?: Record<string, [number, number]>;
  TRI_TETRA_LOOPS?: Record<string, [number, number]>;
  HAIRPIN_LOOPS: Record<number, [number, number]>;
  BULGE_LOOPS: Record<number, [number, number]>;
  INTERNAL_LOOPS: Record<number, [number, number]>;
  MULTIBRANCH: [number, number, number, number];
}

// Default to DNA24 parameters for better accuracy
let USE_DNA24 = true;

/**
 * Set which parameter set to use for folding
 * @param useDna24 - true for DNA24 (recommended), false for SantaLucia 1998
 */
export function setFoldParameterSet(useDna24: boolean): void {
  USE_DNA24 = useDna24;
}

/**
 * Get the current energy parameters based on selected parameter set
 */
function getEnergyParams(): EnergyParameters {
  return USE_DNA24 ? DNA24_ENERGIES : DNA_ENERGIES;
}

/**
 * Convert a pair string to DNA24 terminal mismatch format
 * "AA/TT" -> "AATT" or "AC/TG" -> "ACTG"
 * @param pair - Pair in "XX/YY" format
 * @returns 4-character code for DNA24 lookup
 */
function pairToTerminalCode(pair: string): string {
  // pair format: "XY/ZW" -> "XYZW" (removing the slash)
  return pair.replace("/", "");
}

/**
 * Get 6-character context code for DNA24 internal mismatch lookup
 * The context is: seq[i], seq[i+1], seq[i+2], seq[j-1], seq[j], seq[j+1]
 * where i,j are the outer closing base pair positions
 * @param seq - The full sequence
 * @param i - 5' position of closing pair
 * @param j - 3' position of closing pair
 * @returns 6-character context code or null if invalid
 */
function getInternalMismatchContext(seq: string, i: number, j: number): string | null {
  // For internal mismatches, we need context around the mismatch
  // DNA24 uses: pos i, i+1, i+2 and j-1, j, j+1 (relative to the closing pair)
  if (i + 2 >= seq.length || j - 1 < 0) {
    return null;
  }
  return seq[i] + seq[i + 1] + seq[i + 2] + seq[j - 1] + seq[j] + (j + 1 < seq.length ? seq[j + 1] : "");
}

/**
 * Look up terminal mismatch energy with format handling for DNA24 vs legacy
 * @param pair - Pair in "XX/YY" format
 * @param emap - Energy parameters
 * @returns [dH, dS] or null if not found
 */
function lookupTerminalMM(pair: string, emap: EnergyParameters): [number, number] | null {
  if (USE_DNA24) {
    // DNA24 uses 4-character codes without slash
    const code = pairToTerminalCode(pair);
    if (emap.TERMINAL_MM && code in emap.TERMINAL_MM) {
      return emap.TERMINAL_MM[code];
    }
  } else {
    // Legacy format uses slash notation
    if (emap.TERMINAL_MM && pair in emap.TERMINAL_MM) {
      return emap.TERMINAL_MM[pair];
    }
  }
  return null;
}

/**
 * Look up hairpin mismatch energy (DNA24 has separate HAIRPIN_MM)
 * @param pair - Pair in "XX/YY" format
 * @param emap - Energy parameters
 * @returns [dH, dS] or null if not found
 */
function lookupHairpinMM(pair: string, emap: EnergyParameters): [number, number] | null {
  if (USE_DNA24) {
    // DNA24 has dedicated hairpin mismatch parameters
    const code = pairToTerminalCode(pair);
    if (emap.HAIRPIN_MM && code in emap.HAIRPIN_MM) {
      return emap.HAIRPIN_MM[code];
    }
    // Fall back to terminal MM if hairpin MM not found
    if (emap.TERMINAL_MM && code in emap.TERMINAL_MM) {
      return emap.TERMINAL_MM[code];
    }
  } else {
    // Legacy format uses terminal MM for hairpins too
    if (emap.TERMINAL_MM && pair in emap.TERMINAL_MM) {
      return emap.TERMINAL_MM[pair];
    }
  }
  return null;
}

/**
 * Look up internal mismatch energy with format handling for DNA24 vs legacy
 * @param pair - Pair in "XX/YY" format
 * @param emap - Energy parameters
 * @param seq - Full sequence (needed for DNA24 context)
 * @param i - Position i (for context lookup)
 * @param j - Position j (for context lookup)
 * @returns [dH, dS] or null if not found
 */
function lookupInternalMM(
  pair: string,
  emap: EnergyParameters,
  seq: string | null = null,
  i: number = -1,
  j: number = -1
): [number, number] | null {
  if (USE_DNA24 && seq !== null && i >= 0 && j >= 0) {
    // DNA24 uses 6-character context codes
    const context = getInternalMismatchContext(seq, i, j);
    if (context && emap.INTERNAL_MM && context in emap.INTERNAL_MM) {
      return emap.INTERNAL_MM[context];
    }
    // Fall back to simple 4-char lookup
    const code = pairToTerminalCode(pair);
    if (emap.TERMINAL_MM && code in emap.TERMINAL_MM) {
      return emap.TERMINAL_MM[code];
    }
  } else if (!USE_DNA24) {
    // Legacy format uses slash notation
    if (emap.INTERNAL_MM && pair in emap.INTERNAL_MM) {
      return emap.INTERNAL_MM[pair];
    }
  }
  return null;
}

const STRUCT_DEFAULT: FoldStructure = { e: -Infinity, desc: "", ij: [] };
const STRUCT_NULL: FoldStructure = { e: Infinity, desc: "", ij: [] };

/**
 * Create a structure object
 * @param e - Energy
 * @param desc - Description
 * @param ij - Array of [i, j] pairs
 * @returns Structure object
 */
function createStruct(e: number = -Infinity, desc: string = "", ij: number[][] = []): FoldStructure {
  return { e, desc, ij: ij.slice() };
}

/**
 * Fold the DNA sequence and return the lowest free energy score.
 * Based on Zuker and Stiegler, 1981
 *
 * @param seq - The sequence to fold
 * @param temp - The temperature in Celsius
 * @returns A list of structures
 */
export function fold(seq: string, temp: number = 37.0): FoldStructure[] {
  const [vCache, wCache] = createCache(seq, temp);
  const n = seq.length;

  // get the minimum free energy structure out of the cache
  return traceback(0, n - 1, vCache, wCache);
}

/**
 * Fold the sequence and return just the delta G of the structure
 *
 * @param seq - The sequence to fold
 * @param temp - The temperature to fold at
 * @returns The minimum free energy of the folded sequence
 */
export function dg(seq: string, temp: number = 37.0): number {
  const structs = fold(seq, temp);
  const dgSum = structs.reduce((sum, s) => sum + s.e, 0);
  // Handle Infinity (no valid structure found) - return 0 (no secondary structure)
  if (!Number.isFinite(dgSum)) {
    return 0;
  }
  return Math.round(dgSum * 100) / 100;
}

/**
 * Fold a nucleic acid sequence and return the estimated dg of each (i,j) pairing.
 *
 * @param seq - The nucleic acid sequence to fold
 * @param temp - The temperature to fold at
 * @returns A 2D matrix where each (i, j) pairing corresponds to the MFE
 */
export function dgCache(seq: string, temp: number = 37.0): number[][] {
  const [, wCache] = createCache(seq, temp);

  const cache: number[][] = [];
  for (const row of wCache) {
    cache.push(row.map((s) => s.e));
  }

  return cache;
}

/**
 * Create caches for the w_cache and v_cache
 *
 * @param seq - The sequence to fold
 * @param temp - The temperature to fold at
 * @returns The v_cache and the w_cache
 */
function createCache(seq: string, temp: number = 37.0): [FoldStructure[][], FoldStructure[][]] {
  seq = seq.toUpperCase();
  temp = temp + 273.15; // kelvin

  const emap = getEnergyParams();

  const n = seq.length;
  const vCache: FoldStructure[][] = [];
  const wCache: FoldStructure[][] = [];

  for (let i = 0; i < n; i++) {
    vCache.push(new Array(n).fill(null).map(() => ({ ...STRUCT_DEFAULT })));
    wCache.push(new Array(n).fill(null).map(() => ({ ...STRUCT_DEFAULT })));
  }

  // fill the cache
  computeW(seq, 0, n - 1, temp, vCache, wCache, emap);

  return [vCache, wCache];
}

/**
 * Find and return the lowest free energy structure in Sij subsequence
 * Figure 2B in Zuker and Stiegler, 1981
 *
 * @param seq - The sequence being folded
 * @param i - The start index
 * @param j - The end index (inclusive)
 * @param temp - The temperature in Kelvin
 * @param vCache - Free energy cache for if i and j bp
 * @param wCache - Free energy cache for lowest energy structure
 * @param emap - Energy map
 * @returns The free energy for the subsequence
 */
function computeW(
  seq: string,
  i: number,
  j: number,
  temp: number,
  vCache: FoldStructure[][],
  wCache: FoldStructure[][],
  emap: EnergyParameters
): FoldStructure {
  if (wCache[i][j].e !== -Infinity) {
    return wCache[i][j];
  }

  if (j - i < 4) {
    wCache[i][j] = { ...STRUCT_NULL };
    return wCache[i][j];
  }

  const w1 = computeW(seq, i + 1, j, temp, vCache, wCache, emap);
  const w2 = computeW(seq, i, j - 1, temp, vCache, wCache, emap);
  const w3 = computeV(seq, i, j, temp, vCache, wCache, emap);

  let w4: FoldStructure = { ...STRUCT_NULL };
  for (let k = i + 1; k < j - 1; k++) {
    const w4Test = multiBranch(seq, i, k, j, temp, vCache, wCache, emap, false);

    if (w4Test && w4Test.e < w4.e) {
      w4 = w4Test;
    }
  }

  const w = minStruct(w1, w2, w3, w4);
  wCache[i][j] = w;
  return w;
}

/**
 * Find, store and return the minimum free energy of the structure between i and j
 *
 * @param seq - The sequence being folded
 * @param i - The start index
 * @param j - The end index (inclusive)
 * @param temp - The temperature in Kelvin
 * @param vCache - Free energy cache for if i and j bp
 * @param wCache - Free energy cache for lowest energy structure
 * @param emap - Energy map
 * @returns The minimum energy folding structure
 */
function computeV(
  seq: string,
  i: number,
  j: number,
  temp: number,
  vCache: FoldStructure[][],
  wCache: FoldStructure[][],
  emap: EnergyParameters
): FoldStructure {
  if (vCache[i][j].e !== -Infinity) {
    return vCache[i][j];
  }

  // the ends must basepair for V(i,j)
  if (emap.COMPLEMENT[seq[i]] !== seq[j]) {
    vCache[i][j] = { ...STRUCT_NULL };
    return vCache[i][j];
  }

  // if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
  let isolatedOuter = true;
  if (i && j < seq.length - 1) {
    isolatedOuter = emap.COMPLEMENT[seq[i - 1]] !== seq[j + 1];
  }
  const isolatedInner = emap.COMPLEMENT[seq[i + 1]] !== seq[j - 1];

  if (isolatedOuter && isolatedInner) {
    vCache[i][j] = createStruct(1600);
    return vCache[i][j];
  }

  // E1 = FH(i, j); hairpin
  const pair = getPair(seq, i, i + 1, j, j - 1);
  const e1 = createStruct(hairpin(seq, i, j, temp, emap), "HAIRPIN:" + pair);

  if (j - i === 4) {
    // small hairpin; 4bp
    vCache[i][j] = e1;
    wCache[i][j] = e1;
    return vCache[i][j];
  }

  // E2 = min{FL(i, j, i', j') + V(i', j')}, i<i'<j'<j
  const n = seq.length;
  let e2 = createStruct(Infinity);

  for (let i1 = i + 1; i1 < j - 4; i1++) {
    for (let j1 = i1 + 4; j1 < j; j1++) {
      // i1 and j1 must match
      if (emap.COMPLEMENT[seq[i1]] !== seq[j1]) {
        continue;
      }

      const pairStr = getPair(seq, i, i1, j, j1);
      const pairLeft = getPair(seq, i, i + 1, j, j - 1);
      const pairRight = getPair(seq, i1 - 1, i1, j1 + 1, j1);
      const pairInner = pairLeft in emap.NN || pairRight in emap.NN;

      const stack = i1 === i + 1 && j1 === j - 1;
      const bulgeLeft = i1 > i + 1;
      const bulgeRight = j1 < j - 1;

      let e2Test = Infinity,
        e2TestType = "";

      if (stack) {
        // it's a neighboring/stacking pair in a helix
        e2Test = stackEnergy(seq, i, i1, j, j1, temp, emap);
        e2TestType = `STACK:${pairStr}`;

        if ((i > 0 && j === n - 1) || (i === 0 && j < n - 1)) {
          e2TestType = `STACK_DE:${pairStr}`;
        }
      } else if (bulgeLeft && bulgeRight && !pairInner) {
        // it's an interior loop
        e2Test = internalLoop(seq, i, i1, j, j1, temp, emap);
        e2TestType = `INTERIOR_LOOP:${i1 - i}/${j - j1}`;

        if (i1 - i === 2 && j - j1 === 2) {
          const loopLeft = seq.slice(i, i1 + 1);
          const loopRight = seq.slice(j1, j + 1);
          e2TestType = `STACK:${loopLeft}/${loopRight
            .split("")
            .reverse()
            .join("")}`;
        }
      } else if (bulgeLeft && !bulgeRight) {
        // it's a bulge on the left side
        e2Test = bulge(seq, i, i1, j, j1, temp, emap);
        e2TestType = `BULGE:${i1 - i}`;
      } else if (!bulgeLeft && bulgeRight) {
        // it's a bulge on the right side
        e2Test = bulge(seq, i, i1, j, j1, temp, emap);
        e2TestType = `BULGE:${j - j1}`;
      } else {
        // it's basically a hairpin, only outside bp match
        continue;
      }

      // add V(i', j')
      e2Test += computeV(seq, i1, j1, temp, vCache, wCache, emap).e;
      if (e2Test !== -Infinity && e2Test < e2.e) {
        e2 = createStruct(e2Test, e2TestType, [[i1, j1]]);
      }
    }
  }

  // E3 = min{W(i+1,i') + W(i'+1,j-1)}, i+1<i'<j-2
  let e3: FoldStructure = { ...STRUCT_NULL };
  if (!isolatedOuter || !i || j === seq.length - 1) {
    for (let k = i + 1; k < j - 1; k++) {
      const e3Test = multiBranch(seq, i, k, j, temp, vCache, wCache, emap, true);

      if (e3Test && e3Test.e < e3.e) {
        e3 = e3Test;
      }
    }
  }

  const e = minStruct(e1, e2, e3);
  vCache[i][j] = e;
  return e;
}

/**
 * Return a stack representation, a key for the NN maps
 */
function getPair(s: string, i: number, i1: number, j: number, j1: number): string {
  return (
    (i >= 0 ? s[i] : ".") +
    (i1 >= 0 ? s[i1] : ".") +
    "/" +
    (j >= 0 ? s[j] : ".") +
    (j1 >= 0 ? s[j1] : ".")
  );
}

/**
 * Return the struct with the lowest free energy that isn't -inf (undef)
 */
function minStruct(...structs: FoldStructure[]): FoldStructure {
  let s: FoldStructure = { ...STRUCT_NULL };
  for (const struct of structs) {
    if (struct.e !== -Infinity && struct.e < s.e) {
      s = struct;
    }
  }
  return s;
}

/**
 * Find the free energy given delta h, s and temp
 */
function calcDG(dH: number, dS: number, temp: number): number {
  return dH - temp * (dS / 1000.0);
}

/**
 * Jacobson-Stockmayer extrapolation formula
 */
function jsExtrapolation(queryLen: number, knownLen: number, dgX: number, temp: number): number {
  const gasConstant = 1.9872e-3;
  return dgX + 2.44 * gasConstant * temp * Math.log(queryLen / knownLen);
}

/**
 * Get the free energy for a stack
 */
function stackEnergy(
  seq: string,
  i: number,
  i1: number,
  j: number,
  j1: number,
  temp: number,
  emap: EnergyParameters
): number {
  if ([i, i1, j, j1].some((x) => x >= seq.length)) {
    return 0.0;
  }

  const pair = getPair(seq, i, i1, j, j1);

  if ([i, i1, j, j1].some((x) => x === -1)) {
    // it's a dangling end
    if (pair in emap.DE) {
      const [dH, dS] = emap.DE[pair];
      return calcDG(dH, dS, temp);
    }
    return 0.0;
  }

  if (i > 0 && j < seq.length - 1) {
    // it's internal
    let dH = 0, dS = 0;
    if (pair in emap.NN) {
      [dH, dS] = emap.NN[pair];
    } else {
      const mm = lookupInternalMM(pair, emap, seq, i, j);
      if (mm) [dH, dS] = mm;
    }
    return calcDG(dH, dS, temp);
  }

  if (i === 0 && j === seq.length - 1) {
    // it's terminal
    let dH = 0, dS = 0;
    if (pair in emap.NN) {
      [dH, dS] = emap.NN[pair];
    } else {
      const mm = lookupTerminalMM(pair, emap);
      if (mm) [dH, dS] = mm;
    }
    return calcDG(dH, dS, temp);
  }

  if (i > 0 && j === seq.length - 1) {
    // it's dangling on left
    let dH = 0, dS = 0;
    if (pair in emap.NN) {
      [dH, dS] = emap.NN[pair];
    } else {
      const mm = lookupTerminalMM(pair, emap);
      if (mm) [dH, dS] = mm;
    }
    let dgVal = calcDG(dH, dS, temp);

    const pairDe = seq[i - 1] + seq[i] + "/." + seq[j];
    if (pairDe in emap.DE) {
      const [dH2, dS2] = emap.DE[pairDe];
      dgVal += calcDG(dH2, dS2, temp);
    }
    return dgVal;
  }

  if (i === 0 && j < seq.length - 1) {
    // it's dangling on right
    let dH = 0, dS = 0;
    if (pair in emap.NN) {
      [dH, dS] = emap.NN[pair];
    } else {
      const mm = lookupTerminalMM(pair, emap);
      if (mm) [dH, dS] = mm;
    }
    let dgVal = calcDG(dH, dS, temp);

    const pairDe = "." + seq[i] + "/" + seq[j + 1] + seq[j];
    if (pairDe in emap.DE) {
      const [dH2, dS2] = emap.DE[pairDe];
      dgVal += calcDG(dH2, dS2, temp);
    }
    return dgVal;
  }

  return 0;
}

/**
 * Calculate the free energy of a hairpin
 */
function hairpin(seq: string, i: number, j: number, temp: number, emap: EnergyParameters): number {
  if (j - i < 4) {
    return Infinity;
  }

  const hairpinSeq = seq.slice(i, j + 1);
  const hairpinLen = hairpinSeq.length - 2;
  const pair = getPair(seq, i, i + 1, j, j - 1);

  if (emap.COMPLEMENT[hairpinSeq[0]] !== hairpinSeq[hairpinSeq.length - 1]) {
    throw new Error("not known terminal pair");
  }

  let dgVal = 0.0;

  if (emap.TRI_TETRA_LOOPS && hairpinSeq in emap.TRI_TETRA_LOOPS) {
    const [dH, dS] = emap.TRI_TETRA_LOOPS[hairpinSeq];
    dgVal = calcDG(dH, dS, temp);
  }

  // add penalty based on size
  if (hairpinLen in emap.HAIRPIN_LOOPS) {
    const [dH, dS] = emap.HAIRPIN_LOOPS[hairpinLen];
    dgVal += calcDG(dH, dS, temp);
  } else {
    // it's too large, extrapolate
    const [dH, dS] = emap.HAIRPIN_LOOPS[30];
    const dgInc = calcDG(dH, dS, temp);
    dgVal += jsExtrapolation(hairpinLen, 30, dgInc, temp);
  }

  // add penalty for a terminal mismatch (use hairpin-specific params for DNA24)
  if (hairpinLen > 3) {
    const mm = lookupHairpinMM(pair, emap);
    if (mm) {
      const [dH, dS] = mm;
      dgVal += calcDG(dH, dS, temp);
    }
  }

  // add penalty if length 3 and AT closing
  if (
    hairpinLen === 3 &&
    (hairpinSeq[0] === "A" || hairpinSeq[hairpinSeq.length - 1] === "A")
  ) {
    dgVal += 0.5;
  }

  return dgVal;
}

/**
 * Calculate the free energy associated with a bulge
 */
function bulge(
  seq: string,
  i: number,
  i1: number,
  j: number,
  j1: number,
  temp: number,
  emap: EnergyParameters
): number {
  const loopLen = Math.max(i1 - i - 1, j - j1 - 1);
  if (loopLen <= 0) {
    throw new Error("Invalid bulge");
  }

  let dgVal: number;

  // add penalty based on size
  if (loopLen in emap.BULGE_LOOPS) {
    const [dH, dS] = emap.BULGE_LOOPS[loopLen];
    dgVal = calcDG(dH, dS, temp);
  } else {
    // it's too large, extrapolate
    const [dH, dS] = emap.BULGE_LOOPS[30];
    dgVal = calcDG(dH, dS, temp);
    dgVal = jsExtrapolation(loopLen, 30, dgVal, temp);
  }

  if (loopLen === 1) {
    // if len 1, include the delta G of intervening NN
    dgVal += stackEnergy(seq, i, i1, j, j1, temp, emap);
  }

  // penalize AT terminal bonds
  if ([i, i1, j, j1].some((k) => seq[k] === "A")) {
    dgVal += 0.5;
  }

  return dgVal;
}

/**
 * Calculate the free energy of an internal loop
 */
function internalLoop(
  seq: string,
  i: number,
  i1: number,
  j: number,
  j1: number,
  temp: number,
  emap: EnergyParameters
): number {
  const loopLeft = i1 - i - 1;
  const loopRight = j - j1 - 1;
  const loopLen = loopLeft + loopRight;

  if (loopLeft < 1 || loopRight < 1) {
    throw new Error("Invalid internal loop");
  }

  // single bp mismatch, sum up the two single mismatch pairs
  if (loopLeft === 1 && loopRight === 1) {
    const mmLeft = stackEnergy(seq, i, i1, j, j1, temp, emap);
    const mmRight = stackEnergy(seq, i1 - 1, i1, j1 + 1, j1, temp, emap);
    return mmLeft + mmRight;
  }

  let dgVal: number;

  // apply a penalty based on loop size
  if (loopLen in emap.INTERNAL_LOOPS) {
    const [dH, dS] = emap.INTERNAL_LOOPS[loopLen];
    dgVal = calcDG(dH, dS, temp);
  } else {
    // it's too large an internal loop, extrapolate
    const [dH, dS] = emap.INTERNAL_LOOPS[30];
    dgVal = calcDG(dH, dS, temp);
    dgVal = jsExtrapolation(loopLen, 30, dgVal, temp);
  }

  // apply an asymmetry penalty
  const loopAsymmetry = Math.abs(loopLeft - loopRight);
  dgVal += 0.3 * loopAsymmetry;

  // apply penalty based on the mismatching pairs on either side of the loop
  const pairLeftMm = getPair(seq, i, i + 1, j, j - 1);
  const mmLeft = lookupTerminalMM(pairLeftMm, emap);
  if (mmLeft) {
    const [dH, dS] = mmLeft;
    dgVal += calcDG(dH, dS, temp);
  }

  const pairRightMm = getPair(seq, i1 - 1, i1, j1 + 1, j1);
  const mmRight = lookupTerminalMM(pairRightMm, emap);
  if (mmRight) {
    const [dH, dS] = mmRight;
    dgVal += calcDG(dH, dS, temp);
  }

  return dgVal;
}

/**
 * Calculate a multi-branch energy penalty using a linear formula
 */
function multiBranch(
  seq: string,
  i: number,
  k: number,
  j: number,
  temp: number,
  vCache: FoldStructure[][],
  wCache: FoldStructure[][],
  emap: EnergyParameters,
  helix: boolean = false
): FoldStructure {
  let left: FoldStructure, right: FoldStructure;

  if (helix) {
    left = computeW(seq, i + 1, k, temp, vCache, wCache, emap);
    right = computeW(seq, k + 1, j - 1, temp, vCache, wCache, emap);
  } else {
    left = computeW(seq, i, k, temp, vCache, wCache, emap);
    right = computeW(seq, k + 1, j, temp, vCache, wCache, emap);
  }

  if (!left || left.e === Infinity || !right || right.e === Infinity) {
    return { ...STRUCT_NULL };
  }

  // gather all branches of this multi-branch structure
  const branches: number[][] = [];

  function addBranch(s: FoldStructure): void {
    if (!s || s.e === Infinity || !s.ij || s.ij.length === 0) {
      return;
    }
    if (s.ij.length === 1) {
      branches.push(s.ij[0]);
      return;
    }
    for (const [i1, j1] of s.ij) {
      addBranch(computeW(seq, i1, j1, temp, vCache, wCache, emap));
    }
  }

  addBranch(left);
  addBranch(right);

  // this isn't multi-branched
  if (branches.length < 2) {
    return { ...STRUCT_NULL };
  }

  // if there's a helix, i,j counts as well
  if (helix) {
    branches.push([i, j]);
  }

  // count up unpaired bp and asymmetry
  const branchesCount = branches.length;
  let unpaired = 0;
  let eSum = 0.0;

  for (let index = 0; index < branches.length; index++) {
    const [i2, j2] = branches[index];
    const [, j1] = branches[(index - 1 + branches.length) % branches.length];
    const [i3, j3] = branches[(index + 1) % branches.length];

    let unpairedLeft = 0;
    let unpairedRight = 0;
    let de = 0.0;

    if (index === branches.length - 1 && !helix) {
      // pass
    } else if (i3 === i && j3 === j) {
      unpairedLeft = i2 - j1 - 1;
      unpairedRight = j3 - j2 - 1;

      if (unpairedLeft && unpairedRight) {
        de = stackEnergy(seq, i2 - 1, i2, j2 + 1, j2, temp, emap);
      } else if (unpairedRight) {
        de = stackEnergy(seq, -1, i2, j2 + 1, j2, temp, emap);
        if (unpairedRight === 1) {
          de = Math.min(stackEnergy(seq, i3, -1, j3, j3 - 1, temp, emap), de);
        }
      }
    } else if (i2 === i && j2 === j) {
      unpairedLeft = j2 - j1 - 1;
      unpairedRight = i3 - i2 - 1;

      if (unpairedLeft && unpairedRight) {
        de = stackEnergy(seq, i2 - 1, i2, j2 + 1, j2, temp, emap);
      } else if (unpairedRight) {
        de = stackEnergy(seq, i2, i2 + 1, j2, -1, temp, emap);
        if (unpairedRight === 1) {
          de = Math.min(stackEnergy(seq, i3 - 1, i3, -1, j3, temp, emap), de);
        }
      }
    } else {
      unpairedLeft = i2 - j1 - 1;
      unpairedRight = i3 - j2 - 1;

      if (unpairedLeft && unpairedRight) {
        de = stackEnergy(seq, i2 - 1, i2, j2 + 1, j2, temp, emap);
      } else if (unpairedRight) {
        de = stackEnergy(seq, -1, i2, j2 + 1, j2, temp, emap);
        if (unpairedRight === 1) {
          de = Math.min(stackEnergy(seq, i2 - 1, i2, j2 + 1, j2, temp, emap), de);
        }
      }
    }

    eSum += de;
    unpaired += Math.max(0, unpairedRight);

    if (!(i2 === i && j2 === j)) {
      // add energy
      eSum += computeW(seq, i2, j2, temp, vCache, wCache, emap).e;
    }
  }

  // penalty for unmatched bp and multi-branch
  const [a, b, c, d] = emap.MULTIBRANCH;
  let eMultibranch = a + b * branches.length + c * unpaired;

  if (unpaired === 0) {
    eMultibranch = a + d;
  }

  // energy of min-energy neighbors
  const e = eMultibranch + eSum;

  // pointer to next structures
  const branchesResult = helix ? branches.slice(0, -1) : branches;

  return createStruct(
    e,
    `BIFURCATION:${unpaired}n/${branchesCount}h`,
    branchesResult
  );
}

/**
 * Traceback through the V(i,j) and W(i,j) caches to find the structure
 */
function traceback(i: number, j: number, vCache: FoldStructure[][], wCache: FoldStructure[][]): FoldStructure[] {
  // move i,j down-left to start coordinates
  let sW = wCache[i][j];
  if (!sW.desc.includes("HAIRPIN")) {
    while (
      i + 1 < wCache.length &&
      wCache[i + 1][j].e === sW.e &&
      wCache[i + 1][j].desc === sW.desc
    ) {
      i++;
    }
    while (
      j - 1 >= 0 &&
      wCache[i][j - 1].e === sW.e &&
      wCache[i][j - 1].desc === sW.desc
    ) {
      j--;
    }
  }

  const structs: FoldStructure[] = [];

  while (true) {
    let sV = vCache[i][j];

    // multibrach structures are only in the w_cache.
    if (sW.ij.length > 1) {
      sV = sW;
    }

    structs.push(createStruct(sV.e, sV.desc, [[i, j]]));

    // it's a multibranch
    if (sV.ij.length > 1) {
      let eSum = 0.0;
      let resultStructs = trackbackEnergy(structs);
      const branchStructs: FoldStructure[] = [];

      for (const [i1, j1] of sV.ij) {
        const tb = traceback(i1, j1, vCache, wCache);
        if (tb && tb.length > 0 && tb[0].ij && tb[0].ij.length > 0) {
          const [i2, j2] = tb[0].ij[0];
          eSum += wCache[i2][j2].e;
          branchStructs.push(...tb);
        }
      }

      const last = resultStructs[resultStructs.length - 1];
      resultStructs[resultStructs.length - 1] = createStruct(
        Math.round((last.e - eSum) * 10) / 10,
        last.desc,
        last.ij
      );

      return resultStructs.concat(branchStructs);
    }

    // it's a stack, bulge, etc
    // there's another single structure beyond this
    if (sV.ij.length === 1) {
      [i, j] = sV.ij[0];
      sW = wCache[i][j];
      continue;
    }

    // it's a hairpin, end of structure
    return trackbackEnergy(structs);
  }
}

/**
 * Add energy to each structure, based on how W(i,j) differs from the one after
 */
function trackbackEnergy(structs: FoldStructure[]): FoldStructure[] {
  const structsE: FoldStructure[] = [];
  for (let index = 0; index < structs.length; index++) {
    const struct = structs[index];
    const eNext = index === structs.length - 1 ? 0.0 : structs[index + 1].e;
    const eCorrected = Math.round((struct.e - eNext) * 10) / 10;
    structsE.push(createStruct(eCorrected, struct.desc, struct.ij));
  }
  return structsE;
}

/**
 * Fold the sequence and return a consolidated structure result
 * This is an alias/wrapper for fold() that returns a single structure object
 * with consolidated energy, base pairs, and description.
 *
 * @param seq - The sequence to fold
 * @param temp - The temperature in Celsius
 * @returns Consolidated structure: { e: number, ij: [[i,j],...], desc: string }
 */
export function foldSequence(seq: string, temp: number = 37.0): FoldStructure {
  if (!seq || seq.length < 6) {
    return { e: 0, ij: [], desc: '' };
  }

  try {
    const structures = fold(seq, temp);
    if (!structures || structures.length === 0) {
      return { e: 0, ij: [], desc: '' };
    }

    // Consolidate all structures into one result
    let totalEnergy = 0;
    const allPairs: number[][] = [];
    const descriptions: string[] = [];

    for (const struct of structures) {
      totalEnergy += struct.e;
      if (struct.desc) {
        descriptions.push(struct.desc);
      }
      if (struct.ij) {
        for (const [i, j] of struct.ij) {
          // Only add if this is a valid base pair (not a pointer to inner structure)
          if (i < seq.length && j < seq.length && i !== j) {
            allPairs.push([i, j]);
          }
        }
      }
    }

    // Round total energy
    totalEnergy = Math.round(totalEnergy * 100) / 100;

    // Get primary structure type
    const primaryDesc = descriptions[0]?.split(':')[0] || '';

    return {
      e: totalEnergy,
      ij: allPairs,
      desc: primaryDesc
    };
  } catch (e) {
    console.error('Fold computation failed:', e);
    return { e: 0, ij: [], desc: '' };
  }
}

/**
 * Convert base pair indices to dot-bracket notation for fornac visualization
 * @param seq - The nucleotide sequence
 * @param basePairs - Array of [i, j] base pair indices
 * @returns Dot-bracket notation string (e.g., "((.......))")
 */
export function basePairsToDotBracket(seq: string, basePairs: number[][]): string {
  if (!seq) return '';

  const structure = Array(seq.length).fill('.');

  if (basePairs && basePairs.length > 0) {
    // Sort pairs to ensure proper nesting (opening comes before closing)
    const sortedPairs = [...basePairs].sort((a, b) => a[0] - b[0]);

    for (const [i, j] of sortedPairs) {
      if (i >= 0 && i < seq.length && j >= 0 && j < seq.length && i < j) {
        structure[i] = '(';
        structure[j] = ')';
      }
    }
  }

  return structure.join('');
}
