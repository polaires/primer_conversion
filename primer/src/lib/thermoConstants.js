/**
 * Centralized Thermodynamic Constants
 *
 * This module provides a single source of truth for all thermodynamic
 * parameters used across the primer design system.
 *
 * References:
 * - SantaLucia J Jr. (1998). PNAS 95(4):1460-5
 * - SantaLucia J Jr. & Hicks D (2004) Annu Rev Biophys Biomol Struct 33:415-40
 * - Owczarzy R, et al. (2008). Biochemistry 47(19):5336-53
 * - Allawi HT, SantaLucia J Jr. (1997-1999) - Mismatch parameters
 * - Bommarito S et al. (2000) - Terminal/dangling effects
 */

// =============================================================================
// Watson-Crick Nearest-Neighbor Parameters (SantaLucia & Hicks 2004)
// =============================================================================

/**
 * Nearest-neighbor parameters for perfectly matched DNA duplexes
 * Used by NEB Q5 Tm Calculator
 *
 * Array format: [ΔH (kcal/mol), ΔS (cal/mol·K)]
 */
export const NN_PARAMS_ARRAY = {
  'AA': [-7.6, -21.3],
  'TT': [-7.6, -21.3],
  'AT': [-7.2, -20.4],
  'TA': [-7.2, -21.3],
  'CA': [-8.5, -22.7],
  'TG': [-8.5, -22.7],
  'GT': [-8.4, -22.4],
  'AC': [-8.4, -22.4],
  'CT': [-7.8, -21.0],
  'AG': [-7.8, -21.0],
  'GA': [-8.2, -22.2],
  'TC': [-8.2, -22.2],
  'CG': [-10.6, -27.2],
  'GC': [-9.8, -24.4],
  'GG': [-8.0, -19.9],
  'CC': [-8.0, -19.9],
};

/**
 * Same parameters in object format for modules that prefer named properties
 * Object format: { dH: kcal/mol, dS: cal/(mol·K) }
 */
export const NN_PARAMS_OBJECT = {
  'AA': { dH: -7.6, dS: -21.3 },
  'TT': { dH: -7.6, dS: -21.3 },
  'AT': { dH: -7.2, dS: -20.4 },
  'TA': { dH: -7.2, dS: -21.3 },
  'CA': { dH: -8.5, dS: -22.7 },
  'TG': { dH: -8.5, dS: -22.7 },
  'GT': { dH: -8.4, dS: -22.4 },
  'AC': { dH: -8.4, dS: -22.4 },
  'CT': { dH: -7.8, dS: -21.0 },
  'AG': { dH: -7.8, dS: -21.0 },
  'GA': { dH: -8.2, dS: -22.2 },
  'TC': { dH: -8.2, dS: -22.2 },
  'CG': { dH: -10.6, dS: -27.2 },
  'GC': { dH: -9.8, dS: -24.4 },
  'GG': { dH: -8.0, dS: -19.9 },
  'CC': { dH: -8.0, dS: -19.9 },
};

// =============================================================================
// Initiation Parameters (SantaLucia & Hicks 2004)
// =============================================================================

/**
 * Initiation and terminal penalty parameters
 * Used for Tm calculations
 */
export const INIT_PARAMS = {
  // Terminal AT penalty (applies to each end with A or T)
  terminalAT: { dH: 2.3, dS: 4.1 },
  // Initiation with GC (less common, typically absorbed into NN)
  initGC: { dH: 0.1, dS: -2.8 },
  // Initiation with AT
  initAT: { dH: 2.3, dS: 4.1 },
};

// =============================================================================
// Owczarzy 2008 Mg²⁺ Salt Correction Coefficients
// =============================================================================

/**
 * Coefficients for Owczarzy Mg²⁺ correction (Biochemistry 2008)
 * Used in Q5 Tm calculations
 */
export const OWCZARZY_COEFFICIENTS = {
  a: 3.92e-5,
  b: -9.11e-6,
  c: 6.26e-5,
  d: 1.42e-5,
  e: -4.82e-4,
  f: 5.25e-4,
  g: 8.31e-5,
};

// =============================================================================
// Gas Constant
// =============================================================================

/**
 * Gas constant in cal/(mol·K)
 */
export const R_GAS_CONSTANT = 1.987;

// =============================================================================
// Mismatch Parameters (Allawi & SantaLucia publications)
// =============================================================================

/**
 * Single mismatch parameters
 * G·T mismatches (Allawi & SantaLucia 1997)
 * G·A mismatches (Allawi & SantaLucia 1998)
 * C·T mismatches (Allawi & SantaLucia 1998)
 * A·C mismatches (Allawi & SantaLucia 1998)
 */
export const NN_MISMATCH = {
  // G·T mismatches (wobble pairs - relatively stable)
  'GT/CA': { dH: -0.8, dS: -1.4 },
  'TG/AC': { dH: -1.0, dS: -2.3 },
  'GT/TA': { dH: -1.6, dS: -4.0 },
  'TG/AT': { dH: -1.6, dS: -4.0 },
  'GT/GA': { dH: -2.8, dS: -7.0 },
  'TG/AG': { dH: -2.8, dS: -7.0 },
  'GT/GG': { dH: -3.2, dS: -8.4 },
  'GG/TG': { dH: -3.2, dS: -8.4 },
  'GT/CG': { dH: -4.0, dS: -10.4 },
  'CG/GT': { dH: -4.0, dS: -10.4 },

  // G·A mismatches
  'GA/CA': { dH: -0.7, dS: -0.8 },
  'AG/AC': { dH: -0.7, dS: -0.8 },
  'GA/TA': { dH: -0.5, dS: 0.0 },
  'AG/AT': { dH: -0.5, dS: 0.0 },
  'GA/GA': { dH: 0.5, dS: 3.0 },
  'AA/GA': { dH: 0.7, dS: 3.0 },
  'GA/AA': { dH: 0.7, dS: 3.0 },

  // C·T mismatches
  'CT/GA': { dH: 0.2, dS: 0.7 },
  'TC/AG': { dH: 0.2, dS: 0.7 },
  'CT/CA': { dH: 0.7, dS: 1.0 },
  'TC/AC': { dH: 0.7, dS: 1.0 },
  'CT/TA': { dH: 1.2, dS: 2.0 },
  'TC/AT': { dH: 1.2, dS: 2.0 },

  // A·C mismatches
  'AC/CA': { dH: 2.3, dS: 4.6 },
  'CA/AC': { dH: 2.3, dS: 4.6 },
  'AC/TA': { dH: 5.3, dS: 14.6 },
  'CA/AT': { dH: 5.3, dS: 14.6 },

  // A·A mismatches (Peyret 1999)
  'AA/TA': { dH: 1.2, dS: 1.7 },
  'AA/AA': { dH: 4.7, dS: 12.9 },
  'AA/CA': { dH: 0.6, dS: -0.6 },
  'AA/GA': { dH: -0.9, dS: -4.2 },

  // T·T mismatches
  'TT/AT': { dH: -0.2, dS: -1.5 },
  'TT/TT': { dH: -1.0, dS: -4.4 },

  // C·C mismatches
  'CC/GC': { dH: 0.6, dS: -0.6 },
  'CC/CC': { dH: 3.6, dS: 8.9 },

  // G·G mismatches
  'GG/CG': { dH: -3.1, dS: -9.5 },
  'GG/GG': { dH: -1.4, dS: -6.2 },
};

// =============================================================================
// Terminal Mismatch Corrections (Bommarito 2000)
// =============================================================================

/**
 * Terminal mismatches are less destabilizing than internal ones
 */
export const TERMINAL_CORRECTIONS = {
  '5prime': { dH: 0.4, dS: 0.5 },   // 5' mismatches less destabilizing
  '3prime': { dH: -0.6, dS: -1.0 }, // 3' mismatches more destabilizing (critical!)
};

// =============================================================================
// Dangling End Parameters (Bommarito 2000)
// =============================================================================

/**
 * Parameters for single unpaired nucleotides at the ends
 */
export const DANGLING_END_CORRECTIONS = {
  // 5' dangling ends (less impact on stability)
  '5_A': { dH: 0.2, dS: 2.3 },
  '5_T': { dH: -6.9, dS: -20.0 },
  '5_G': { dH: -3.9, dS: -10.9 },
  '5_C': { dH: -4.4, dS: -12.6 },
  // 3' dangling ends (more impact on stability)
  '3_A': { dH: -0.7, dS: -0.8 },
  '3_T': { dH: -0.5, dS: -1.1 },
  '3_G': { dH: -5.9, dS: -16.5 },
  '3_C': { dH: -2.1, dS: -3.9 },
};

// =============================================================================
// Consecutive Mismatch Parameters (Peyret 1999)
// =============================================================================

/**
 * Tandem mismatches have different thermodynamics than single mismatches
 */
export const CONSECUTIVE_MISMATCH_CORRECTION = {
  // Penalty for consecutive mismatches (per additional mismatch)
  additionalPenalty: { dH: 0.5, dS: 1.5 },
  // Maximum destabilization cap
  maxPenalty: { dH: 5.0, dS: 15.0 },
};

// =============================================================================
// Backward Compatibility Aliases
// =============================================================================

// Alias for modules expecting 'NN_PARAMS' in array format (tmQ5.js style)
export const NN_PARAMS = NN_PARAMS_ARRAY;

// Alias for modules expecting 'NN_MATCHED' in object format (mutagenesis.js style)
export const NN_MATCHED = NN_PARAMS_OBJECT;

// Alias for Owczarzy coefficients
export const OWCZARZY = OWCZARZY_COEFFICIENTS;
