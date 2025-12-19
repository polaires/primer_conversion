/**
 * Multi-Polymerase Tm Calculator Library
 *
 * Implements different Tm calculation algorithms for various DNA polymerases.
 * Each polymerase has unique buffer formulations and thermodynamic requirements.
 *
 * Supported Polymerases:
 * - NEB Q5: High-fidelity, SantaLucia + Owczarzy + Q5 empirical correction
 * - Thermo Phusion: High-fidelity, Breslauer + Schildkraut salt correction
 * - Taq: Standard polymerase, basic nearest-neighbor calculation
 * - NEB OneTaq: Hot-start Taq blend, optimized for broader Tm range
 * - Deep Vent: Thermococcus polymerase, high thermostability
 *
 * References:
 * - SantaLucia J Jr. (1998). PNAS 95(4):1460-5 [Unified NN parameters]
 * - Breslauer KJ, et al. (1986). PNAS 83(11):3746-50 [Original NN parameters]
 * - Owczarzy R, et al. (2008). Biochemistry 47(19):5336-53 [Mg2+ correction]
 * - Schildkraut C, Lifson S. (1965). Biopolymers 3(2):195-208 [Na+ correction]
 */

import {
  NN_PARAMS,
  OWCZARZY,
  R_GAS_CONSTANT as R,
} from './thermoConstants.js';

// Breslauer (1986) nearest-neighbor parameters (kcal/mol and cal/mol·K)
// Used by Phusion and some other polymerases
const BRESLAUER_PARAMS = {
  'AA': [-9.1, -24.0], 'TT': [-9.1, -24.0],
  'AT': [-8.6, -23.9],
  'TA': [-6.0, -16.9],
  'CA': [-5.8, -12.9], 'TG': [-5.8, -12.9],
  'GT': [-6.5, -17.3], 'AC': [-6.5, -17.3],
  'CT': [-7.8, -20.8], 'AG': [-7.8, -20.8],
  'GA': [-5.6, -13.5], 'TC': [-5.6, -13.5],
  'CG': [-11.9, -27.8],
  'GC': [-11.1, -26.7],
  'GG': [-11.0, -26.6], 'CC': [-11.0, -26.6],
};

/**
 * Calculate GC content of a sequence
 */
export function calculateGC(seq) {
  if (!seq || typeof seq !== 'string' || seq.length === 0) return 0;
  const upper = seq.toUpperCase();
  let gc = 0, valid = 0;
  for (const c of upper) {
    if (c === 'G' || c === 'C') { gc++; valid++; }
    else if (c === 'A' || c === 'T') { valid++; }
  }
  return valid === 0 ? 0 : gc / valid;
}

/**
 * Polymerase definitions with their specific properties
 */
export const POLYMERASES = {
  q5: {
    id: 'q5',
    name: 'NEB Q5',
    fullName: 'Q5 High-Fidelity DNA Polymerase',
    manufacturer: 'New England Biolabs',
    type: 'High-Fidelity',
    fidelity: '280× Taq',
    processivity: 'High',
    extension: '20-30 sec/kb',
    hotStart: true,
    proofReading: true,
    method: 'SantaLucia 1998 + Owczarzy Mg²⁺ + Q5 empirical',
    defaults: {
      primerConc: 500,   // nM
      mgConc: 2.0,       // mM
      annealOffset: 1,   // °C added to lower Tm
      maxAnneal: 72,     // Maximum annealing temp
    },
    protocol: {
      initialDenature: { temp: 98, time: '30 sec' },
      denature: { temp: 98, time: '10 sec' },
      extension: { temp: 72, time: '20-30 sec/kb' },
      finalExtension: { temp: 72, time: '2 min' },
      cycles: '25-35',
    },
    notes: 'Best for cloning, mutagenesis, and applications requiring high fidelity. Uses 2-step or 3-step PCR.',
    color: '#0d9488', // Teal - distinctive for NEB Q5
  },

  phusion: {
    id: 'phusion',
    name: 'Phusion',
    fullName: 'Phusion High-Fidelity DNA Polymerase',
    manufacturer: 'Thermo Fisher',
    type: 'High-Fidelity',
    fidelity: '52× Taq',
    processivity: 'High',
    extension: '15-30 sec/kb',
    hotStart: false,
    proofReading: true,
    method: 'Breslauer 1986 + Schildkraut salt correction',
    defaults: {
      primerConc: 500,
      mgConc: 1.5,
      annealOffset: 3,    // Phusion anneals 3°C above lower Tm
      maxAnneal: 72,
    },
    protocol: {
      initialDenature: { temp: 98, time: '30 sec' },
      denature: { temp: 98, time: '10 sec' },
      extension: { temp: 72, time: '15-30 sec/kb' },
      finalExtension: { temp: 72, time: '5-10 min' },
      cycles: '25-35',
    },
    notes: 'Fast extension rate, suitable for long and difficult templates. Generates blunt-ended products.',
    color: '#3498db',
  },

  taq: {
    id: 'taq',
    name: 'Taq',
    fullName: 'Taq DNA Polymerase',
    manufacturer: 'Various',
    type: 'Standard',
    fidelity: '1× (baseline)',
    processivity: 'Medium',
    extension: '1 min/kb',
    hotStart: false,
    proofReading: false,
    method: 'Basic nearest-neighbor (SantaLucia) + Na⁺ correction',
    defaults: {
      primerConc: 200,
      mgConc: 1.5,
      annealOffset: -5,   // Taq anneals 5°C below average Tm
      maxAnneal: 65,
    },
    protocol: {
      initialDenature: { temp: 95, time: '2-3 min' },
      denature: { temp: 95, time: '30 sec' },
      extension: { temp: 72, time: '1 min/kb' },
      finalExtension: { temp: 72, time: '5-10 min' },
      cycles: '25-35',
    },
    notes: 'Classic workhorse polymerase. Adds A-overhangs suitable for TA cloning. Lower cost, lower fidelity.',
    color: '#27ae60',
  },

  onetaq: {
    id: 'onetaq',
    name: 'OneTaq',
    fullName: 'OneTaq DNA Polymerase',
    manufacturer: 'New England Biolabs',
    type: 'Hot-Start Taq Blend',
    fidelity: '~2× Taq',
    processivity: 'Medium-High',
    extension: '1 min/kb',
    hotStart: true,
    proofReading: false,
    method: 'SantaLucia + Owczarzy Mg²⁺ (simplified)',
    defaults: {
      primerConc: 200,
      mgConc: 1.8,
      annealOffset: 0,    // Use calculated Tm directly
      maxAnneal: 68,
    },
    protocol: {
      initialDenature: { temp: 94, time: '30 sec' },
      denature: { temp: 94, time: '30 sec' },
      extension: { temp: 68, time: '1 min/kb' },
      finalExtension: { temp: 68, time: '5 min' },
      cycles: '25-35',
    },
    notes: 'Hot-start Taq blend optimized for routine PCR. Good balance of convenience and performance.',
    color: '#9b59b6',
  },

  deepvent: {
    id: 'deepvent',
    name: 'Deep Vent',
    fullName: 'Deep Vent DNA Polymerase',
    manufacturer: 'New England Biolabs',
    type: 'Archaeal High-Fidelity',
    fidelity: '4-5× Taq',
    processivity: 'Low',
    extension: '2 min/kb',
    hotStart: false,
    proofReading: true,
    method: 'SantaLucia + Owczarzy Mg²⁺',
    defaults: {
      primerConc: 500,
      mgConc: 2.0,
      annealOffset: 0,
      maxAnneal: 72,
    },
    protocol: {
      initialDenature: { temp: 95, time: '2-3 min' },
      denature: { temp: 95, time: '30 sec' },
      extension: { temp: 72, time: '2 min/kb' },
      finalExtension: { temp: 72, time: '10 min' },
      cycles: '25-35',
    },
    notes: 'From Thermococcus, extremely thermostable (half-life 23h at 95°C). Generates blunt ends.',
    color: '#e67e22',
  },

  kapa: {
    id: 'kapa',
    name: 'KAPA HiFi',
    fullName: 'KAPA HiFi DNA Polymerase',
    manufacturer: 'Roche (KAPA Biosystems)',
    type: 'High-Fidelity',
    fidelity: '100× Taq',
    processivity: 'Very High',
    extension: '15-60 sec/kb',
    hotStart: true,
    proofReading: true,
    method: 'SantaLucia + Owczarzy Mg²⁺',
    defaults: {
      primerConc: 300,
      mgConc: 2.5,
      annealOffset: 0,
      maxAnneal: 72,
    },
    protocol: {
      initialDenature: { temp: 95, time: '3 min' },
      denature: { temp: 98, time: '20 sec' },
      extension: { temp: 72, time: '15-60 sec/kb' },
      finalExtension: { temp: 72, time: '1 min/kb' },
      cycles: '25-35',
    },
    notes: 'Engineered B-family polymerase. High speed, high fidelity. Popular for NGS library prep.',
    color: '#1abc9c',
  },
};

/**
 * ===================================================================
 * TM CALCULATION METHODS
 * ===================================================================
 */

/**
 * Calculate Tm using NEB Q5 algorithm (SantaLucia + Owczarzy + Q5 empirical)
 */
export function calculateTmQ5(sequence, options = {}) {
  const { primerConc = 500, mgConc = 2.0 } = options;
  const seq = sequence.toUpperCase().replace(/[^ATGC]/g, '');

  if (seq.length < 2) throw new Error('Sequence must be at least 2 nucleotides');

  // Step 1: Nearest-neighbor (SantaLucia)
  let dH = 0.2, dS = -5.7;
  for (let i = 0; i < seq.length - 1; i++) {
    const dinuc = seq[i] + seq[i + 1];
    const params = NN_PARAMS[dinuc];
    if (params) { dH += params[0]; dS += params[1]; }
  }

  // Terminal AT penalties
  if (seq[0] === 'A' || seq[0] === 'T') { dH += 2.2; dS += 6.9; }
  if (seq[seq.length - 1] === 'A' || seq[seq.length - 1] === 'T') { dH += 2.2; dS += 6.9; }

  const Ct = primerConc * 1e-9;
  const Tm_1M = (dH * 1000) / (dS + R * Math.log(Ct / 4)) - 273.15;

  // Step 2: Owczarzy Mg²⁺ correction
  const N = seq.length;
  const fGC = calculateGC(seq);
  const Mg = mgConc * 1e-3;
  const lnMg = Math.log(Mg);

  const correction =
    (OWCZARZY.a + OWCZARZY.b * lnMg + fGC * (OWCZARZY.c + OWCZARZY.d * lnMg)) +
    (1 / (2 * (N - 1))) * (OWCZARZY.e + OWCZARZY.f * lnMg + OWCZARZY.g * lnMg * lnMg);

  const invTm = (1 / (Tm_1M + 273.15)) + correction;
  const Tm_Mg = (1 / invTm) - 273.15;

  // Step 3: Q5 empirical correction
  const Tm_NEB = 18.25 + 0.949 * Tm_Mg + 8.67 * fGC - 5.25 * Math.log(N) + 0.12 * N - 77.05 / N;

  return Math.round(Tm_NEB);
}

/**
 * Calculate Tm using Breslauer method + Schildkraut salt correction (Phusion)
 */
export function calculateTmPhusion(sequence, options = {}) {
  const { primerConc = 500, naConc = 50 } = options; // naConc in mM
  const seq = sequence.toUpperCase().replace(/[^ATGC]/g, '');

  if (seq.length < 2) throw new Error('Sequence must be at least 2 nucleotides');

  // Use Breslauer parameters
  let dH = 0, dS = 0;
  for (let i = 0; i < seq.length - 1; i++) {
    const dinuc = seq[i] + seq[i + 1];
    const params = BRESLAUER_PARAMS[dinuc];
    if (params) { dH += params[0]; dS += params[1]; }
  }

  // Initiation parameter for Breslauer
  dH += -2.0;
  dS += -10.0;

  // Calculate Tm at standard conditions
  const Ct = primerConc * 1e-9;
  const Tm_basic = (dH * 1000) / (dS + R * Math.log(Ct / 4)) - 273.15;

  // Schildkraut salt correction: Tm = Tm_basic + 16.6 * log10([Na+])
  const Na = naConc * 1e-3; // Convert mM to M
  const Tm_corrected = Tm_basic + 16.6 * Math.log10(Na);

  return Math.round(Tm_corrected);
}

/**
 * Calculate Tm using basic nearest-neighbor (standard Taq)
 * Uses simplified calculation suitable for Taq polymerase
 */
export function calculateTmTaq(sequence, options = {}) {
  const { primerConc = 200, naConc = 50 } = options;
  const seq = sequence.toUpperCase().replace(/[^ATGC]/g, '');

  if (seq.length < 2) throw new Error('Sequence must be at least 2 nucleotides');

  // Use SantaLucia parameters
  let dH = 0.2, dS = -5.7;
  for (let i = 0; i < seq.length - 1; i++) {
    const dinuc = seq[i] + seq[i + 1];
    const params = NN_PARAMS[dinuc];
    if (params) { dH += params[0]; dS += params[1]; }
  }

  // Terminal penalties
  if (seq[0] === 'A' || seq[0] === 'T') { dH += 2.2; dS += 6.9; }
  if (seq[seq.length - 1] === 'A' || seq[seq.length - 1] === 'T') { dH += 2.2; dS += 6.9; }

  const Ct = primerConc * 1e-9;
  const Tm_basic = (dH * 1000) / (dS + R * Math.log(Ct / 4)) - 273.15;

  // Simple Na+ correction
  const Na = naConc * 1e-3;
  const Tm_corrected = Tm_basic + 16.6 * Math.log10(Na);

  return Math.round(Tm_corrected);
}

/**
 * Calculate Tm using SantaLucia + Owczarzy (OneTaq, Deep Vent, KAPA)
 */
export function calculateTmOwczarzy(sequence, options = {}) {
  const { primerConc = 500, mgConc = 2.0 } = options;
  const seq = sequence.toUpperCase().replace(/[^ATGC]/g, '');

  if (seq.length < 2) throw new Error('Sequence must be at least 2 nucleotides');

  // Step 1: Nearest-neighbor (SantaLucia)
  let dH = 0.2, dS = -5.7;
  for (let i = 0; i < seq.length - 1; i++) {
    const dinuc = seq[i] + seq[i + 1];
    const params = NN_PARAMS[dinuc];
    if (params) { dH += params[0]; dS += params[1]; }
  }

  if (seq[0] === 'A' || seq[0] === 'T') { dH += 2.2; dS += 6.9; }
  if (seq[seq.length - 1] === 'A' || seq[seq.length - 1] === 'T') { dH += 2.2; dS += 6.9; }

  const Ct = primerConc * 1e-9;
  const Tm_1M = (dH * 1000) / (dS + R * Math.log(Ct / 4)) - 273.15;

  // Step 2: Owczarzy Mg²⁺ correction
  const N = seq.length;
  const fGC = calculateGC(seq);
  const Mg = mgConc * 1e-3;
  const lnMg = Math.log(Mg);

  const correction =
    (OWCZARZY.a + OWCZARZY.b * lnMg + fGC * (OWCZARZY.c + OWCZARZY.d * lnMg)) +
    (1 / (2 * (N - 1))) * (OWCZARZY.e + OWCZARZY.f * lnMg + OWCZARZY.g * lnMg * lnMg);

  const invTm = (1 / (Tm_1M + 273.15)) + correction;
  const Tm_Mg = (1 / invTm) - 273.15;

  return Math.round(Tm_Mg);
}

/**
 * ===================================================================
 * UNIFIED CALCULATION INTERFACE
 * ===================================================================
 */

/**
 * Calculate Tm for any polymerase
 *
 * @param {string} sequence - DNA primer sequence
 * @param {string} polymeraseId - Polymerase ID (q5, phusion, taq, onetaq, deepvent, kapa)
 * @param {Object} options - Optional concentration overrides
 * @returns {number} Melting temperature in °C
 */
export function calculateTm(sequence, polymeraseId = 'q5', options = {}) {
  const polymerase = POLYMERASES[polymeraseId];
  if (!polymerase) {
    throw new Error(`Unknown polymerase: ${polymeraseId}`);
  }

  const mergedOptions = { ...polymerase.defaults, ...options };

  switch (polymeraseId) {
    case 'q5':
      return calculateTmQ5(sequence, mergedOptions);
    case 'phusion':
      return calculateTmPhusion(sequence, mergedOptions);
    case 'taq':
      return calculateTmTaq(sequence, mergedOptions);
    case 'onetaq':
    case 'deepvent':
    case 'kapa':
      return calculateTmOwczarzy(sequence, mergedOptions);
    default:
      return calculateTmOwczarzy(sequence, mergedOptions);
  }
}

/**
 * Calculate annealing temperature for a primer pair
 *
 * @param {string} primer1 - Forward primer sequence
 * @param {string} primer2 - Reverse primer sequence
 * @param {string} polymeraseId - Polymerase ID
 * @param {Object} options - Optional overrides
 * @returns {Object} Annealing temperature result
 */
export function calculateAnnealing(primer1, primer2, polymeraseId = 'q5', options = {}) {
  const polymerase = POLYMERASES[polymeraseId];
  if (!polymerase) {
    throw new Error(`Unknown polymerase: ${polymeraseId}`);
  }

  const mergedOptions = { ...polymerase.defaults, ...options };
  const { annealOffset, maxAnneal } = mergedOptions;

  const tm1 = calculateTm(primer1, polymeraseId, mergedOptions);
  const tm2 = calculateTm(primer2, polymeraseId, mergedOptions);

  const tmLower = Math.min(tm1, tm2);
  const tmHigher = Math.max(tm1, tm2);

  // Calculate annealing temp based on polymerase type
  let annealingTemp;
  if (polymeraseId === 'taq') {
    // Taq: Ta = (Tm1 + Tm2) / 2 + offset (usually -5)
    annealingTemp = Math.round((tm1 + tm2) / 2) + annealOffset;
  } else {
    // High-fidelity: Ta = Tm_lower + offset
    annealingTemp = tmLower + annealOffset;
  }

  // Cap at maximum
  const isCapped = annealingTemp > maxAnneal;
  annealingTemp = Math.min(annealingTemp, maxAnneal);

  return {
    tm1,
    tm2,
    tmLower,
    tmHigher,
    tmDifference: Math.abs(tm1 - tm2),
    annealingTemp,
    isCapped,
    polymerase: polymerase.name,
    method: polymerase.method,
  };
}

/**
 * Get detailed Tm calculation for display
 */
export function calculateTmDetailed(sequence, polymeraseId = 'q5', options = {}) {
  const polymerase = POLYMERASES[polymeraseId];
  if (!polymerase) {
    throw new Error(`Unknown polymerase: ${polymeraseId}`);
  }

  const mergedOptions = { ...polymerase.defaults, ...options };
  const seq = sequence.toUpperCase().replace(/[^ATGC]/g, '');

  if (seq.length < 2) {
    throw new Error('Sequence must be at least 2 nucleotides');
  }

  const tm = calculateTm(sequence, polymeraseId, mergedOptions);
  const gc = calculateGC(seq);

  return {
    sequence: seq,
    length: seq.length,
    gcContent: gc,
    gcPercent: (gc * 100).toFixed(1) + '%',
    tm,
    polymerase: polymerase.name,
    method: polymerase.method,
    options: mergedOptions,
  };
}

/**
 * Get list of available polymerases for UI
 */
export function getPolymeraseList() {
  return Object.values(POLYMERASES).map(p => ({
    id: p.id,
    name: p.name,
    fullName: p.fullName,
    manufacturer: p.manufacturer,
    type: p.type,
    color: p.color,
  }));
}

/**
 * Get polymerase info by ID
 */
export function getPolymeraseInfo(polymeraseId) {
  return POLYMERASES[polymeraseId] || null;
}

export default {
  calculateTm,
  calculateAnnealing,
  calculateTmDetailed,
  calculateGC,
  getPolymeraseList,
  getPolymeraseInfo,
  POLYMERASES,
};
