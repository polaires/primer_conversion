/**
 * Sequencing Primer Design Module
 *
 * Designs primers optimized for Sanger sequencing based on peer-reviewed research:
 * - SantaLucia & Hicks 2004 nearest-neighbor thermodynamics
 * - Vallone & Butler 2004 - AutoDimer algorithm
 * - Unified scoring from primerAnalysis.js (calibrated on Döring dataset)
 *
 * Key differences from PCR primers:
 * - No salt correction (primers shipped dry for commercial sequencing)
 * - Tm range 50-65°C (optimal 55-60°C)
 * - Emphasis on avoiding secondary structures
 * - Position 50-600 bp upstream of target for readable signal
 *
 * This module follows the unified primer design patterns:
 * - Uses analyzeSinglePrimer() from primerAnalysis.js for scoring
 * - Uses equilibrium.js for hairpin/dimer calculations (no duplication)
 * - Supports circular sequences for plasmid sequencing
 * - Generates alternative primers at each position
 */

import { calculateGC, calculate3primeTerminalDG } from './tmQ5.js';
import { offTargets } from './offTargets.js';
import {
  scoreTm,
  scoreGc,
  scoreTerminal3DG,
  scoreHairpin,
  scoreHomodimer,
  scoreOffTarget,
  scoreLength,
  scoreGcClamp,
  scoreHomopolymer,
  calculateCompositeScore,
  classifyQuality,
  analyzeGQuadruplex,
  score3PrimeComposition,
} from './scoring.js';
import { calculateHairpinDG, calculateHomodimerDG } from './equilibrium.js';
import { ANALYSIS_PRESETS } from './presets.js';
import { NN_PARAMS_OBJECT as NN_PARAMS, INIT_PARAMS } from './thermoConstants.js';
import { reverseComplement } from './sequenceUtils.js';
import { analyzeSinglePrimer } from './primerAnalysis.js';

// =============================================================================
// Default Parameters for Sequencing Primers
// =============================================================================

export const SEQUENCING_DEFAULTS = {
  // Primer parameters (consensus: 18-24 bp optimal)
  minPrimerLength: 18,
  maxPrimerLength: 24,
  optimalPrimerLength: 20,

  // Tm parameters for cycle sequencing
  // Target range: 50-65°C, ideally 55-60°C for BigDye chemistry
  // Note: Using SantaLucia 1998 NN without salt correction
  minTm: 50,
  maxTm: 65,
  optimalTm: 58,

  // GC parameters (40-60%, optimal 45-55%)
  minGC: 0.40,
  maxGC: 0.60,
  optimalGC: 0.50,

  // GC Clamp parameters
  gcClampRequired: true,          // Require G/C in last 5 bases
  minGCinLast5: 1,                // Minimum G/C in last 5 bases
  maxGCinLast5: 2,                // Maximum G/C in last 5 bases (avoid >3)
  maxConsecutiveGCat3prime: 3,    // Max consecutive G/C at 3' end

  // Homopolymer and repeat limits
  maxPolyN: 4,                    // Max consecutive identical bases
  maxDinucleotideRepeat: 4,       // Max dinucleotide repeat units (e.g., ATATATAT = 4)

  // Secondary structure thresholds (ΔG in kcal/mol)
  maxHairpinDG: -2.0,             // 3' end hairpin threshold
  maxInternalHairpinDG: -3.0,     // Internal hairpin threshold
  maxSelfDimerDG: -5.0,           // 3' end self-dimer threshold
  maxInternalDimerDG: -6.0,       // Internal self-dimer threshold

  // Tm balance (5' vs 3' half)
  maxTmImbalance: 5,              // Max Tm difference between halves

  // Off-target checking
  checkOffTargets: true,          // Enable off-target binding site detection
  maxOffTargets: 0,               // Maximum allowed off-target sites (0 = ideal)

  // Read parameters
  readStartOffset: 50,            // First 15-40 bp are poor quality
  optimalReadLength: 500,         // High-quality reads: 50-500 bp after primer
  minReadLength: 400,
  maxReadLength: 700,

  // Primer spacing for walking
  primerSpacing: 350,             // Spacing to ensure overlap with 500bp reads
  spacingTolerance: 50,
  minOverlap: 150,                // At least 150bp overlap between reads

  // Coverage strategy
  bidirectionalCoverage: true,
  ensureEdgeCoverage: true,

  // Circular sequence support (for plasmid sequencing)
  circular: false,                // Set to true for plasmids/circular templates

  // Alternative primer generation
  generateAlternatives: true,     // Generate alternative primers at each position
  maxAlternatives: 5,             // Maximum alternatives per position

  // Scoring weights (based on BatchPrimer3)
  weights: {
    length: 0.5,
    tm: 2.0,
    gc: 1.5,
    gcClamp: 3.0,
    polyN: 5.0,
    dinucRepeat: 4.0,
    hairpin: 4.0,
    selfDimer: 5.0,
    tmBalance: 1.0,
    offTarget: 15.0,              // Strong penalty for off-target binding
  },
};

// =============================================================================
// Thermodynamic Calculations
// =============================================================================

/**
 * Calculate Tm using SantaLucia 1998 nearest-neighbor method
 * WITHOUT salt correction (for commercial sequencing - primers shipped dry)
 *
 * Tm = ΔH / (ΔS + R × ln(C/4)) - 273.15
 *
 * @param {string} seq - Primer sequence
 * @param {number} primerConc - Primer concentration in M (default 0.25 µM)
 * @returns {number} Melting temperature in °C
 */
export function calculateTmNN(seq, primerConc = 0.25e-6) {
  const sequence = seq.toUpperCase();
  const R = 1.987; // Gas constant cal/(mol·K)

  let dH = 0;
  let dS = 0;

  // Sum nearest-neighbor contributions
  for (let i = 0; i < sequence.length - 1; i++) {
    const dinuc = sequence.slice(i, i + 2);
    const params = NN_PARAMS[dinuc];
    if (params) {
      dH += params.dH;
      dS += params.dS;
    }
  }

  // Add initiation parameters
  const firstBase = sequence[0];
  const lastBase = sequence[sequence.length - 1];

  // Terminal AT penalty
  if (firstBase === 'A' || firstBase === 'T') {
    dH += INIT_PARAMS.terminalAT.dH;
    dS += INIT_PARAMS.terminalAT.dS;
  }
  if (lastBase === 'A' || lastBase === 'T') {
    dH += INIT_PARAMS.terminalAT.dH;
    dS += INIT_PARAMS.terminalAT.dS;
  }

  // Symmetry correction for self-complementary sequences
  if (isSymmetric(sequence)) {
    dS += -1.4; // Symmetry correction
  }

  // Calculate Tm (convert dH to cal/mol for consistency with dS)
  // Tm = ΔH / (ΔS + R × ln(Ct/4)) - 273.15
  const Ct = primerConc;
  const tm = (dH * 1000) / (dS + R * Math.log(Ct / 4)) - 273.15;

  return Math.round(tm * 10) / 10;
}

/**
 * Check if sequence is self-complementary (palindromic)
 */
function isSymmetric(seq) {
  const rc = reverseComplement(seq);
  return seq === rc;
}

/**
 * Calculate ΔG for a sequence at given temperature
 * Used for hairpin and dimer stability assessment
 */
function calculateDeltaG(seq, tempC = 37) {
  const sequence = seq.toUpperCase();
  const tempK = tempC + 273.15;

  let dH = 0;
  let dS = 0;

  for (let i = 0; i < sequence.length - 1; i++) {
    const dinuc = sequence.slice(i, i + 2);
    const params = NN_PARAMS[dinuc];
    if (params) {
      dH += params.dH;
      dS += params.dS;
    }
  }

  // ΔG = ΔH - T × ΔS (convert dS to kcal)
  return dH - (tempK * dS / 1000);
}

// =============================================================================
// Circular Sequence Helpers (aligned with unifiedPrimerDesign.js pattern)
// =============================================================================

/**
 * Prepare a working sequence for circular designs by doubling the template.
 * This allows standard algorithms to work across the origin.
 *
 * @param {string} seq - The template sequence
 * @param {boolean} circular - Whether sequence is circular
 * @returns {Object} { workingSeq, originalLength }
 */
function prepareCircularWorkingSequence(seq, circular) {
  if (!circular) {
    return {
      workingSeq: seq,
      originalLength: seq.length,
      isCircular: false,
    };
  }

  // Double the sequence to handle positions that wrap around origin
  return {
    workingSeq: seq + seq,
    originalLength: seq.length,
    isCircular: true,
  };
}

/**
 * Adjust primer coordinates back to circular range after designing on doubled sequence.
 *
 * @param {Object} primer - Primer object with position
 * @param {number} seqLength - Original sequence length
 * @returns {Object} Primer with adjusted coordinates
 */
function adjustCircularCoordinates(primer, seqLength) {
  const result = { ...primer };
  if (result.position !== undefined && result.position >= seqLength) {
    result.position = result.position % seqLength;
    result.wrapsOrigin = true;
  }
  return result;
}

/**
 * Calculate position considering circular wrap-around.
 *
 * @param {number} position - Position in working sequence
 * @param {number} seqLength - Original sequence length
 * @param {boolean} circular - Whether sequence is circular
 * @returns {number} Adjusted position
 */
function circularPosition(position, seqLength, circular) {
  if (!circular) return position;
  return ((position % seqLength) + seqLength) % seqLength;
}

// =============================================================================
// Secondary Structure Detection (uses equilibrium.js - no duplication)
// =============================================================================
// NOTE: Hairpin and self-dimer detection now uses the shared functions from
// equilibrium.js (calculateHairpinDG, calculateHomodimerDG) to avoid code
// duplication and ensure consistent thermodynamic calculations.
//
// Legacy detectHairpins() and detectSelfDimer() functions have been removed.
// Use the equilibrium.js functions directly or via analyzeSinglePrimer().
// =============================================================================

// =============================================================================
// Sequence Composition Analysis
// =============================================================================

/**
 * Detect homopolymer runs (consecutive identical bases)
 */
function detectHomopolymers(seq) {
  const matches = seq.toUpperCase().match(/(.)\1{2,}/g) || [];
  return matches.map(m => ({
    base: m[0],
    length: m.length,
    position: seq.toUpperCase().indexOf(m),
  })).sort((a, b) => b.length - a.length);
}

/**
 * Detect dinucleotide repeats (e.g., ATATAT, GCGCGC)
 */
function detectDinucleotideRepeats(seq) {
  const sequence = seq.toUpperCase();
  const repeats = [];

  // Check each dinucleotide pattern
  const patterns = ['AT', 'TA', 'GC', 'CG', 'AC', 'CA', 'GT', 'TG', 'AG', 'GA', 'CT', 'TC'];

  for (const pattern of patterns) {
    const regex = new RegExp(`(${pattern}){2,}`, 'g');
    let match;
    while ((match = regex.exec(sequence)) !== null) {
      repeats.push({
        pattern,
        sequence: match[0],
        units: match[0].length / 2,
        position: match.index,
      });
    }
  }

  return repeats.sort((a, b) => b.units - a.units);
}

/**
 * Analyze GC clamp at 3' end
 */
function analyzeGCClamp(seq) {
  const sequence = seq.toUpperCase();
  const last5 = sequence.slice(-5);
  const last3 = sequence.slice(-3);

  let gcInLast5 = 0;
  let consecutiveGCat3 = 0;
  let currentStreak = 0;

  for (const base of last5) {
    if (base === 'G' || base === 'C') gcInLast5++;
  }

  // Count consecutive G/C at 3' end
  for (let i = sequence.length - 1; i >= 0; i--) {
    if (sequence[i] === 'G' || sequence[i] === 'C') {
      currentStreak++;
    } else {
      break;
    }
  }
  consecutiveGCat3 = currentStreak;

  return {
    gcInLast5,
    consecutiveGCat3,
    last5,
    hasProperClamp: gcInLast5 >= 1 && gcInLast5 <= 2 && consecutiveGCat3 <= 3,
    lastBase: sequence[sequence.length - 1],
    hasGCatEnd: sequence[sequence.length - 1] === 'G' || sequence[sequence.length - 1] === 'C',
  };
}

/**
 * Calculate Tm balance between 5' and 3' halves
 */
function calculateTmBalance(seq) {
  const sequence = seq.toUpperCase();
  const midpoint = Math.floor(sequence.length / 2);

  const fiveHalf = sequence.slice(0, midpoint);
  const threeHalf = sequence.slice(midpoint);

  const tm5 = calculateTmNN(fiveHalf);
  const tm3 = calculateTmNN(threeHalf);

  return {
    tm5prime: tm5,
    tm3prime: tm3,
    difference: Math.abs(tm5 - tm3),
    isBalanced: Math.abs(tm5 - tm3) <= 5,
  };
}

// =============================================================================
// Unified Scoring Function (uses primerAnalysis.js)
// =============================================================================

/**
 * Calculate comprehensive quality score for a sequencing primer
 *
 * This function now uses the unified analyzeSinglePrimer() from primerAnalysis.js
 * to ensure consistent scoring across all primer design modes while maintaining
 * sequencing-specific parameters (no salt Tm, stricter structure thresholds).
 *
 * @param {string} sequence - Primer sequence
 * @param {Object} options - Scoring options
 * @param {string} template - Template sequence for off-target checking (optional)
 * @returns {Object} Quality assessment with detailed analysis (backward compatible)
 */
export function scorePrimerQuality(sequence, options = {}, template = null) {
  const opts = { ...SEQUENCING_DEFAULTS, ...options };
  const weights = { ...SEQUENCING_DEFAULTS.weights, ...options.weights };
  const seq = sequence.toUpperCase();

  // Calculate Tm using sequencing-specific method (no salt correction)
  const tm = calculateTmNN(seq);
  const gc = calculateGC(seq);

  // Off-target analysis
  let offTargetCount = 0;
  if (opts.checkOffTargets && template) {
    const otCache = offTargets(seq, template);
    offTargetCount = otCache[seq.length - 1] || 0;
  }

  // Use unified analysis from primerAnalysis.js
  // This provides: G-Quadruplex, 3' composition, hairpin/dimer via equilibrium.js
  const unifiedAnalysis = analyzeSinglePrimer(
    { seq, tm, gc, offTargetCount },
    {
      mode: 'sequencing',
      template,
      temperature: 55,
      include3PrimeAnalysis: true,
      includeGQuadruplex: true,
    }
  );

  // Get thermodynamic values from equilibrium.js (via unified analysis)
  const hairpinDG = unifiedAnalysis.thermodynamics?.hairpinDG ?? calculateHairpinDG(seq, 55);
  const homodimerDG = unifiedAnalysis.thermodynamics?.homodimerDG ?? calculateHomodimerDG(seq, 55);
  const terminal3DG = unifiedAnalysis.thermodynamics?.terminal3DG ?? calculate3primeTerminalDG(seq).dG;

  // Composition analysis (still needed for detailed output)
  const homopolymers = detectHomopolymers(seq);
  const dinucRepeats = detectDinucleotideRepeats(seq);
  const gcClamp = analyzeGCClamp(seq);
  const tmBalance = calculateTmBalance(seq);

  // Build issues and warnings lists from unified analysis
  const issues = [];
  const warnings = [];
  let penalty = 0;

  // === Length scoring ===
  if (seq.length < opts.minPrimerLength) {
    issues.push(`Too short (${seq.length} bp, min ${opts.minPrimerLength})`);
    penalty += 50;
  } else if (seq.length > opts.maxPrimerLength) {
    issues.push(`Too long (${seq.length} bp, max ${opts.maxPrimerLength})`);
    penalty += (seq.length - opts.maxPrimerLength) * weights.length * 5;
  } else {
    penalty += Math.abs(seq.length - opts.optimalPrimerLength) * weights.length;
  }

  // === Tm scoring ===
  if (tm < opts.minTm) {
    issues.push(`Tm too low (${tm}°C, min ${opts.minTm}°C)`);
    penalty += (opts.minTm - tm) * weights.tm;
  } else if (tm > opts.maxTm) {
    issues.push(`Tm too high (${tm}°C, max ${opts.maxTm}°C)`);
    penalty += (tm - opts.maxTm) * weights.tm;
  } else {
    penalty += Math.abs(tm - opts.optimalTm) * weights.tm * 0.2;
  }

  // === GC content scoring ===
  if (gc < opts.minGC) {
    issues.push(`GC too low (${(gc * 100).toFixed(0)}%, min ${opts.minGC * 100}%)`);
    penalty += (opts.minGC - gc) * 100 * weights.gc;
  } else if (gc > opts.maxGC) {
    issues.push(`GC too high (${(gc * 100).toFixed(0)}%, max ${opts.maxGC * 100}%)`);
    penalty += (gc - opts.maxGC) * 100 * weights.gc;
  } else {
    penalty += Math.abs(gc - opts.optimalGC) * 100 * weights.gc * 0.1;
  }

  // === GC Clamp scoring ===
  if (opts.gcClampRequired) {
    if (gcClamp.gcInLast5 < opts.minGCinLast5) {
      warnings.push(`Weak 3' end (${gcClamp.gcInLast5} G/C in last 5 bases)`);
      penalty += (opts.minGCinLast5 - gcClamp.gcInLast5) * weights.gcClamp;
    } else if (gcClamp.gcInLast5 > opts.maxGCinLast5) {
      warnings.push(`Strong 3' end (${gcClamp.gcInLast5} G/C in last 5 bases)`);
      penalty += (gcClamp.gcInLast5 - opts.maxGCinLast5) * weights.gcClamp;
    }
    if (gcClamp.consecutiveGCat3 > opts.maxConsecutiveGCat3prime) {
      warnings.push(`${gcClamp.consecutiveGCat3} consecutive G/C at 3' end`);
      penalty += (gcClamp.consecutiveGCat3 - opts.maxConsecutiveGCat3prime) * weights.gcClamp * 2;
    }
  }

  // === Homopolymer scoring ===
  const maxHomopolymer = homopolymers[0];
  if (maxHomopolymer && maxHomopolymer.length > opts.maxPolyN) {
    issues.push(`Poly-${maxHomopolymer.base} run of ${maxHomopolymer.length} bases`);
    penalty += (maxHomopolymer.length - opts.maxPolyN) * weights.polyN;
  } else if (maxHomopolymer && maxHomopolymer.length === opts.maxPolyN) {
    warnings.push(`Poly-${maxHomopolymer.base} run of ${maxHomopolymer.length} bases`);
    penalty += weights.polyN * 0.5;
  }

  // === Dinucleotide repeat scoring ===
  const maxRepeat = dinucRepeats[0];
  if (maxRepeat && maxRepeat.units > opts.maxDinucleotideRepeat) {
    warnings.push(`${maxRepeat.pattern} repeat (${maxRepeat.units} units)`);
    penalty += (maxRepeat.units - opts.maxDinucleotideRepeat) * weights.dinucRepeat;
  }

  // === Hairpin scoring (using equilibrium.js values) ===
  if (hairpinDG < opts.maxHairpinDG) {
    issues.push(`Hairpin structure (ΔG = ${hairpinDG.toFixed(1)} kcal/mol)`);
    penalty += Math.abs(hairpinDG - opts.maxHairpinDG) * weights.hairpin;
  } else if (hairpinDG < opts.maxInternalHairpinDG) {
    warnings.push(`Potential hairpin (ΔG = ${hairpinDG.toFixed(1)} kcal/mol)`);
    penalty += Math.abs(hairpinDG - opts.maxInternalHairpinDG) * weights.hairpin * 0.5;
  }

  // === Self-dimer scoring (using equilibrium.js values) ===
  if (homodimerDG < opts.maxSelfDimerDG) {
    issues.push(`Self-dimer (ΔG = ${homodimerDG.toFixed(1)} kcal/mol)`);
    penalty += Math.abs(homodimerDG - opts.maxSelfDimerDG) * weights.selfDimer;
  } else if (homodimerDG < opts.maxInternalDimerDG) {
    warnings.push(`Self-dimer potential (ΔG = ${homodimerDG.toFixed(1)} kcal/mol)`);
    penalty += Math.abs(homodimerDG - opts.maxInternalDimerDG) * weights.selfDimer * 0.5;
  }

  // === G-Quadruplex scoring (NEW - from unified analysis) ===
  const gQuadruplex = unifiedAnalysis.gQuadruplex;
  if (gQuadruplex?.severity === 'critical') {
    issues.push(gQuadruplex.message);
    penalty += 20;
  } else if (gQuadruplex?.severity === 'warning') {
    warnings.push(gQuadruplex.message);
    penalty += 10;
  }

  // === 3' Composition scoring (NEW - from unified analysis) ===
  if (unifiedAnalysis.analysis3Prime?.quality === 'poor') {
    warnings.push(`3' end quality issues: ${unifiedAnalysis.analysis3Prime.issues.join(', ')}`);
    penalty += 5;
  }

  // === Tm balance scoring ===
  if (!tmBalance.isBalanced) {
    warnings.push(`Tm imbalance (5': ${tmBalance.tm5prime}°C, 3': ${tmBalance.tm3prime}°C)`);
    penalty += (tmBalance.difference - opts.maxTmImbalance) * weights.tmBalance;
  }

  // === Off-target scoring ===
  if (opts.checkOffTargets && offTargetCount > opts.maxOffTargets) {
    if (offTargetCount >= 3) {
      issues.push(`${offTargetCount} off-target binding sites`);
      penalty += offTargetCount * weights.offTarget;
    } else {
      warnings.push(`${offTargetCount} off-target binding site${offTargetCount > 1 ? 's' : ''}`);
      penalty += offTargetCount * weights.offTarget * 0.5;
    }
  }

  // Add warnings from unified analysis that we haven't added yet
  for (const w of unifiedAnalysis.warnings || []) {
    const msg = w.message;
    if (!warnings.includes(msg) && !issues.some(i => i.includes(msg.slice(0, 20)))) {
      if (w.severity === 'critical') {
        issues.push(msg);
      } else {
        warnings.push(msg);
      }
    }
  }

  // Get piecewise scores from unified analysis
  const piecewiseScores = unifiedAnalysis.scores || {
    tm: scoreTm(tm, ANALYSIS_PRESETS.sequencing.tmOptions),
    gc: scoreGc(gc, ANALYSIS_PRESETS.sequencing.gcOptions),
    length: scoreLength(seq.length, ANALYSIS_PRESETS.sequencing.lengthOptions),
    gcClamp: scoreGcClamp(seq),
    homopolymer: scoreHomopolymer(seq),
    hairpin: scoreHairpin(hairpinDG, { threshold: ANALYSIS_PRESETS.sequencing.hairpinThreshold }),
    homodimer: scoreHomodimer(homodimerDG, { threshold: ANALYSIS_PRESETS.sequencing.homodimerThreshold }),
    offTarget: scoreOffTarget(offTargetCount),
    terminal3DG: scoreTerminal3DG(terminal3DG),
    gQuadruplex: gQuadruplex?.score ?? 1.0,
    threePrimeComp: score3PrimeComposition(seq, terminal3DG),
  };

  // Calculate composite score using calibrated weights
  const compositeResult = calculateCompositeScore({
    tmFwd: piecewiseScores.tm,
    gcFwd: piecewiseScores.gc,
    lengthFwd: piecewiseScores.length,
    gcClampFwd: piecewiseScores.gcClamp,
    homopolymerFwd: piecewiseScores.homopolymer,
    hairpinFwd: piecewiseScores.hairpin,
    selfDimerFwd: piecewiseScores.homodimer,
    offTarget: piecewiseScores.offTarget,
    terminal3DG: piecewiseScores.terminal3DG,
    gQuadruplexFwd: piecewiseScores.gQuadruplex,
    threePrimeCompFwd: piecewiseScores.threePrimeComp,
  });

  const compositeScore = compositeResult.score;
  const qualityFromComposite = classifyQuality(compositeScore);

  // Determine overall quality
  let quality;
  if (issues.length === 0 && penalty < 5) {
    quality = 'excellent';
  } else if (issues.length === 0 && penalty < 15) {
    quality = 'good';
  } else if (issues.length === 0 && penalty < 30) {
    quality = 'acceptable';
  } else {
    quality = 'poor';
  }

  // Use calibrated quality tier when no critical issues
  const finalQualityTier = issues.length === 0 ? qualityFromComposite.tier : quality;

  return {
    sequence: seq,
    length: seq.length,
    tm,
    gc,
    gcPercent: (gc * 100).toFixed(1) + '%',

    // Detailed analysis
    gcClamp,
    hairpinAnalysis: {
      hasHairpin: hairpinDG < -1.0,
      worstDG: hairpinDG,
    },
    selfDimerAnalysis: {
      hasDimer: homodimerDG < -3.0,
      worstDG: homodimerDG,
    },
    tmBalance,
    homopolymers: maxHomopolymer || null,
    dinucRepeats: maxRepeat || null,
    offTargetCount,

    // NEW: G-Quadruplex and 3' composition analysis
    gQuadruplex,
    analysis3Prime: unifiedAnalysis.analysis3Prime,

    // Legacy compatibility
    hasGCClamp: gcClamp.hasGCatEnd,
    hasStrongGCClamp: gcClamp.gcInLast5 >= 2,

    // Quality assessment (legacy)
    quality,
    penalty: Math.round(penalty * 10) / 10,
    issues,
    warnings,

    // Calibrated scoring (unified with primerAnalysis.js)
    piecewiseScores,
    compositeScore,
    qualityTier: finalQualityTier,
    scoring: {
      penalty: Math.round(penalty * 10) / 10,
      piecewiseScores,
      compositeScore,
      qualityTier: finalQualityTier,
    },

    // Thermodynamics from equilibrium.js
    thermodynamics: {
      hairpinDG,
      homodimerDG,
      terminal3DG,
    },
  };
}

// =============================================================================
// Primer Design Functions
// =============================================================================

/**
 * Design sequencing primers for a DNA template
 *
 * Uses a bidirectional walking primer strategy based on:
 * - Geneious primer spacing algorithm (450-500 bp ± 50)
 * - ABI 3730/BigDye v3.1 read specifications (50-500 bp usable)
 * - Edge coverage with complementary directional primers
 *
 * Features (aligned with unified primer design patterns):
 * - Circular sequence support for plasmid sequencing
 * - Alternative primer generation at each position
 * - Unified scoring via analyzeSinglePrimer()
 *
 * @param {string} template - Template sequence
 * @param {Object} options - Design options
 * @param {boolean} options.circular - Whether template is circular (plasmid)
 * @param {boolean} options.generateAlternatives - Generate alternative primers
 * @param {number} options.maxAlternatives - Max alternatives per position
 * @returns {Object} Designed primers with coverage map and alternatives
 */
export function designSequencingPrimers(template, options = {}) {
  const opts = { ...SEQUENCING_DEFAULTS, ...options };
  const originalSeq = template.toUpperCase();

  // Handle circular sequences by doubling the template
  const { workingSeq, originalLength, isCircular } = prepareCircularWorkingSequence(
    originalSeq,
    opts.circular
  );
  const seq = isCircular ? workingSeq : originalSeq;

  if (seq.length < 100) {
    throw new Error('Template too short for sequencing primer design (min 100 bp)');
  }

  const primers = [];
  const targetPositions = [];

  // Calculate effective coverage region per primer
  const effectiveReadStart = opts.readStartOffset;
  const effectiveReadLength = opts.optimalReadLength;
  const stepSize = opts.primerSpacing || (effectiveReadLength - opts.minOverlap);

  // === EDGE COVERAGE: Beginning of sequence ===
  if (opts.ensureEdgeCoverage) {
    const idealRevPos = effectiveReadStart * 2 + opts.optimalPrimerLength;
    const revEdgePos = Math.min(idealRevPos, seq.length - opts.minPrimerLength);
    if (revEdgePos > opts.minPrimerLength) {
      targetPositions.push({
        position: revEdgePos,
        direction: 'reverse',
        priority: 'edge',
        readStart: 0,
        readEnd: Math.max(revEdgePos - opts.readStartOffset, 50),
      });
    }
  }

  // === FORWARD PRIMERS: Walk from start ===
  let currentPos = 0;
  while (currentPos < seq.length) {
    const readStart = currentPos + effectiveReadStart;
    const readEnd = Math.min(readStart + effectiveReadLength, seq.length);

    if (readStart < seq.length && readEnd > readStart) {
      targetPositions.push({
        position: currentPos,
        direction: 'forward',
        readStart,
        readEnd,
      });
    }

    currentPos += stepSize;
  }

  // === REVERSE PRIMERS: Walk from end ===
  if (opts.bidirectionalCoverage || seq.length > 800) {
    currentPos = seq.length;
    while (currentPos > effectiveReadStart) {
      const readEnd = currentPos - effectiveReadStart;
      const readStart = Math.max(readEnd - effectiveReadLength, 0);

      if (readEnd > 0 && readEnd > readStart) {
        targetPositions.push({
          position: currentPos,
          direction: 'reverse',
          readStart,
          readEnd,
        });
      }

      currentPos -= stepSize;
    }
  }

  // === EDGE COVERAGE: End of sequence ===
  if (opts.ensureEdgeCoverage && seq.length > effectiveReadLength) {
    const idealFwdPos = Math.max(0, seq.length - effectiveReadLength - effectiveReadStart - opts.optimalPrimerLength);
    targetPositions.push({
      position: idealFwdPos,
      direction: 'forward',
      priority: 'edge',
      readStart: idealFwdPos + opts.optimalPrimerLength + effectiveReadStart,
      readEnd: seq.length,
    });

    const revEndPos = seq.length;
    targetPositions.push({
      position: revEndPos,
      direction: 'reverse',
      priority: 'edge',
      readStart: Math.max(0, seq.length - effectiveReadStart - effectiveReadLength),
      readEnd: seq.length - effectiveReadStart,
    });
  }

  // === DESIGN OPTIMAL PRIMERS ===
  const baseSearchWindow = opts.spacingTolerance || 50;

  for (const target of targetPositions) {
    const candidates = [];
    const searchWindow = target.priority === 'edge' ? baseSearchWindow * 2 : baseSearchWindow;

    if (target.direction === 'forward') {
      const windowStart = Math.max(0, target.position - searchWindow);
      const windowEnd = Math.min(seq.length - opts.minPrimerLength, target.position + searchWindow);

      for (let pos = windowStart; pos <= windowEnd; pos++) {
        for (let len = opts.minPrimerLength; len <= opts.maxPrimerLength; len++) {
          const endPos = pos + len;
          if (endPos > seq.length) continue;

          const primerSeq = seq.slice(pos, endPos);
          if (primerSeq.length < opts.minPrimerLength) continue;

          const quality = scorePrimerQuality(primerSeq, opts, seq);
          const distancePenalty = Math.abs(pos - target.position) * 0.01;

          candidates.push({
            ...quality,
            penalty: quality.penalty + distancePenalty,
            position: pos,
            direction: 'forward',
            name: `Seq_F${Math.floor(pos / 100) * 100}`,
            readStart: pos + len + opts.readStartOffset,
            readEnd: Math.min(pos + len + opts.readStartOffset + effectiveReadLength, seq.length),
            priority: target.priority,
          });
        }
      }
    } else {
      const windowStart = Math.max(opts.minPrimerLength, target.position - searchWindow);
      const windowEnd = Math.min(seq.length, target.position + searchWindow);

      for (let pos = windowStart; pos <= windowEnd; pos++) {
        for (let len = opts.minPrimerLength; len <= opts.maxPrimerLength; len++) {
          const startPos = pos - len;
          if (startPos < 0) continue;

          const primerSeq = reverseComplement(seq.slice(startPos, pos));
          if (primerSeq.length < opts.minPrimerLength) continue;

          const quality = scorePrimerQuality(primerSeq, opts, seq);
          const distancePenalty = Math.abs(pos - target.position) * 0.01;

          candidates.push({
            ...quality,
            penalty: quality.penalty + distancePenalty,
            position: startPos,
            direction: 'reverse',
            name: `Seq_R${Math.floor(startPos / 100) * 100}`,
            readStart: Math.max(startPos - opts.readStartOffset - effectiveReadLength, 0),
            readEnd: Math.max(startPos - opts.readStartOffset, 0),
            priority: target.priority,
          });
        }
      }
    }

    // Select best candidate and generate alternatives
    if (candidates.length > 0) {
      candidates.sort((a, b) => a.penalty - b.penalty);
      const best = candidates[0];

      // Adjust position for circular sequences
      if (isCircular && best.position >= originalLength) {
        best.position = best.position % originalLength;
        best.wrapsOrigin = true;
      }

      const minDistance = target.priority === 'edge' ? 100 : stepSize / 2;
      const overlaps = primers.some(p =>
        p.direction === best.direction &&
        Math.abs(p.position - best.position) < minDistance
      );

      const qualityOk = target.priority === 'edge'
        ? best.penalty < 50
        : best.quality !== 'poor';

      if (!overlaps && qualityOk) {
        // Generate alternatives if enabled (aligned with unified pattern)
        if (opts.generateAlternatives && candidates.length > 1) {
          const maxAlts = opts.maxAlternatives || 5;
          best.alternatives = candidates
            .slice(1, maxAlts + 1)
            // Use compositeScore (calibrated) instead of penalty-based quality
            .filter(alt => alt.compositeScore >= 50)
            .map(alt => {
              // Adjust position for circular
              if (isCircular && alt.position >= originalLength) {
                alt.position = alt.position % originalLength;
                alt.wrapsOrigin = true;
              }
              return {
                sequence: alt.sequence,
                length: alt.length,
                tm: alt.tm,
                gc: alt.gc,
                gcPercent: alt.gcPercent,
                position: alt.position,
                direction: alt.direction,
                compositeScore: alt.compositeScore,
                qualityTier: alt.qualityTier,
                penalty: alt.penalty,
                scoreDelta: alt.compositeScore - best.compositeScore,
              };
            });
        }
        primers.push(best);
      }
    }
  }

  // Calculate initial coverage
  let coverage = calculateCoverage(primers, seq.length, opts);

  // === RESCUE PRIMERS ===
  if (coverage.gaps.length > 0) {
    const rescuePrimers = [];

    for (const gap of coverage.gaps) {
      if (gap.length < 50) continue;

      // Forward rescue primer
      const targetFwdPos = Math.max(0, gap.start - opts.readStartOffset - opts.optimalPrimerLength - 50);
      const fwdCandidates = [];
      const rescueWindow = 150;

      for (let pos = Math.max(0, targetFwdPos - rescueWindow); pos <= Math.min(seq.length - opts.minPrimerLength, targetFwdPos + rescueWindow); pos++) {
        for (let len = opts.minPrimerLength; len <= opts.maxPrimerLength; len++) {
          if (pos + len > seq.length) continue;

          const primerSeq = seq.slice(pos, pos + len);
          const quality = scorePrimerQuality(primerSeq, opts, seq);

          const adjustedPenalty = quality.tm > opts.maxTm && quality.tm <= 75
            ? quality.penalty - (quality.tm - opts.maxTm) * 1.5
            : quality.penalty;

          fwdCandidates.push({
            ...quality,
            penalty: adjustedPenalty,
            position: pos,
            direction: 'forward',
            readStart: pos + len + opts.readStartOffset,
            readEnd: Math.min(pos + len + opts.readStartOffset + effectiveReadLength, seq.length),
            isRescue: true,
          });
        }
      }

      if (fwdCandidates.length > 0) {
        fwdCandidates.sort((a, b) => a.penalty - b.penalty);
        const best = fwdCandidates[0];

        // Adjust position for circular sequences
        if (isCircular && best.position >= originalLength) {
          best.position = best.position % originalLength;
          best.wrapsOrigin = true;
        }

        const overlaps = primers.concat(rescuePrimers).some(p =>
          p.direction === 'forward' &&
          Math.abs(p.position - best.position) < 100
        );

        if (!overlaps && best.penalty < 50) {  // Increased threshold for rescue primers
          best.name = `Seq_F_rescue`;
          best.warnings.push('Rescue primer with adjusted quality threshold');
          rescuePrimers.push(best);
        }
      }

      // Reverse rescue primer
      const targetRevPos = gap.end + opts.readStartOffset + opts.optimalPrimerLength;
      const revCandidates = [];

      for (let pos = Math.max(opts.minPrimerLength, targetRevPos - rescueWindow); pos <= Math.min(seq.length, targetRevPos + rescueWindow); pos++) {
        for (let len = opts.minPrimerLength; len <= opts.maxPrimerLength; len++) {
          const startPos = pos - len;
          if (startPos < 0) continue;

          const primerSeq = reverseComplement(seq.slice(startPos, pos));
          if (primerSeq.length < opts.minPrimerLength) continue;

          const quality = scorePrimerQuality(primerSeq, opts, seq);

          const adjustedPenalty = quality.tm > opts.maxTm && quality.tm <= 75
            ? quality.penalty - (quality.tm - opts.maxTm) * 1.5
            : quality.penalty;

          const readEnd = startPos - opts.readStartOffset;
          const readStart = Math.max(0, readEnd - effectiveReadLength);

          if (readStart <= gap.start && readEnd >= gap.end) {
            revCandidates.push({
              ...quality,
              penalty: adjustedPenalty,
              position: startPos,
              direction: 'reverse',
              readStart,
              readEnd,
              isRescue: true,
            });
          }
        }
      }

      if (revCandidates.length > 0) {
        revCandidates.sort((a, b) => a.penalty - b.penalty);
        const best = revCandidates[0];

        // Adjust position for circular sequences
        if (isCircular && best.position >= originalLength) {
          best.position = best.position % originalLength;
          best.wrapsOrigin = true;
        }

        const overlaps = primers.concat(rescuePrimers).some(p =>
          p.direction === 'reverse' &&
          Math.abs(p.position - best.position) < 100
        );

        if (!overlaps && best.penalty < 50) {  // Increased threshold for rescue primers
          best.name = `Seq_R_rescue`;
          best.warnings.push('Rescue primer with adjusted quality threshold');
          rescuePrimers.push(best);
        }
      }
    }

    primers.push(...rescuePrimers);
    coverage = calculateCoverage(primers, seq.length, opts);
  }

  // Sort and rename primers
  primers.sort((a, b) => {
    if (a.direction !== b.direction) {
      return a.direction === 'forward' ? -1 : 1;
    }
    return a.position - b.position;
  });

  let fwdCount = 1;
  let revCount = 1;
  for (const primer of primers) {
    if (primer.direction === 'forward') {
      primer.name = primer.isRescue ? `Seq_F${fwdCount++}*` : `Seq_F${fwdCount++}`;
    } else {
      primer.name = primer.isRescue ? `Seq_R${revCount++}*` : `Seq_R${revCount++}`;
    }
  }

  // Recalculate coverage using original length for circular sequences
  const coverageLength = isCircular ? originalLength : seq.length;
  coverage = calculateCoverage(primers, coverageLength, opts);

  // Calculate Tm compatibility analysis for the primer set
  const tmCompatibility = analyzeTmCompatibility(primers, opts);

  return {
    template: originalSeq,
    templateLength: originalLength,
    primers,
    primerCount: primers.length,
    forwardPrimers: primers.filter(p => p.direction === 'forward'),
    reversePrimers: primers.filter(p => p.direction === 'reverse'),
    coverage,
    estimatedReads: primers.length,
    estimatedCost: primers.length * 8,
    hasRescuePrimers: primers.some(p => p.isRescue),
    hasAlternatives: primers.some(p => p.alternatives && p.alternatives.length > 0),
    // Circular sequence info
    isCircular,
    circularWrapped: primers.some(p => p.wrapsOrigin),
    // Tm compatibility analysis
    tmCompatibility,
  };
}

/**
 * Analyze Tm compatibility across a primer set
 * Good sequencing primer sets should have similar Tms for consistent reaction conditions
 *
 * @param {Array} primers - Array of primer objects
 * @param {Object} opts - Options including optimal Tm
 * @returns {Object} Tm compatibility analysis
 */
function analyzeTmCompatibility(primers, opts) {
  if (primers.length === 0) {
    return { compatible: true, range: 0, mean: 0, stdDev: 0, outliers: [] };
  }

  const tms = primers.map(p => p.tm);
  const minTm = Math.min(...tms);
  const maxTm = Math.max(...tms);
  const range = maxTm - minTm;
  const mean = tms.reduce((a, b) => a + b, 0) / tms.length;

  // Calculate standard deviation
  const squaredDiffs = tms.map(tm => Math.pow(tm - mean, 2));
  const avgSquaredDiff = squaredDiffs.reduce((a, b) => a + b, 0) / tms.length;
  const stdDev = Math.sqrt(avgSquaredDiff);

  // Identify outliers (more than 5°C from mean or outside acceptable range)
  const outliers = primers.filter(p =>
    Math.abs(p.tm - mean) > 5 ||
    p.tm < opts.minTm ||
    p.tm > opts.maxTm
  ).map(p => ({
    name: p.name,
    tm: p.tm,
    deviation: Math.round((p.tm - mean) * 10) / 10,
  }));

  // Compatible if range is ≤5°C and no outliers outside acceptable range
  const compatible = range <= 5 && outliers.length === 0;

  return {
    compatible,
    range: Math.round(range * 10) / 10,
    minTm: Math.round(minTm * 10) / 10,
    maxTm: Math.round(maxTm * 10) / 10,
    mean: Math.round(mean * 10) / 10,
    stdDev: Math.round(stdDev * 10) / 10,
    outliers,
    recommendation: compatible
      ? 'All primers have compatible Tm values'
      : range > 5
        ? `Tm range of ${range.toFixed(1)}°C may affect reaction consistency`
        : `${outliers.length} primer(s) outside optimal Tm range`,
  };
}

/**
 * Design a single sequencing primer at a specific position
 *
 * @param {string} template - Template sequence
 * @param {number} position - Target position
 * @param {string} direction - 'forward' or 'reverse'
 * @param {Object} options - Design options
 * @returns {Object} Best primer with optional alternatives
 */
export function designPrimerAtPosition(template, position, direction = 'forward', options = {}) {
  const opts = { ...SEQUENCING_DEFAULTS, ...options };
  const seq = template.toUpperCase();

  const candidates = [];

  for (let len = opts.minPrimerLength; len <= opts.maxPrimerLength; len++) {
    let primerSeq;
    let startPos;

    if (direction === 'forward') {
      startPos = Math.max(0, position - len);
      primerSeq = seq.slice(startPos, startPos + len);
    } else {
      startPos = Math.min(seq.length - len, position);
      primerSeq = reverseComplement(seq.slice(startPos, startPos + len));
    }

    if (primerSeq.length < opts.minPrimerLength) continue;

    const quality = scorePrimerQuality(primerSeq, opts, seq);

    candidates.push({
      ...quality,
      position: startPos,
      direction,
      targetPosition: position,
    });
  }

  if (candidates.length === 0) {
    throw new Error('Could not design primer at specified position');
  }

  candidates.sort((a, b) => a.penalty - b.penalty);
  const best = candidates[0];

  // Add alternatives if requested
  if (opts.generateAlternatives && candidates.length > 1) {
    const maxAlts = opts.maxAlternatives || 5;
    best.alternatives = candidates
      .slice(1, maxAlts + 1)
      // Use compositeScore (calibrated) instead of penalty-based quality
      .filter(alt => alt.compositeScore >= 50)
      .map(alt => ({
        sequence: alt.sequence,
        length: alt.length,
        tm: alt.tm,
        gc: alt.gc,
        gcPercent: alt.gcPercent,
        position: alt.position,
        compositeScore: alt.compositeScore,
        qualityTier: alt.qualityTier,
        scoreDelta: alt.compositeScore - best.compositeScore,
      }));
  }

  return best;
}

// =============================================================================
// Coverage Calculation
// =============================================================================

function calculateCoverage(primers, templateLength, opts) {
  const coverageArray = new Array(templateLength).fill(0);

  for (const primer of primers) {
    const start = primer.readStart ?? (primer.position + opts.readStartOffset);
    const end = primer.readEnd ?? Math.min(start + opts.optimalReadLength, templateLength);

    for (let i = Math.max(0, start); i < Math.min(end, templateLength); i++) {
      coverageArray[i]++;
    }
  }

  const covered = coverageArray.filter(c => c > 0).length;
  const doubleCovered = coverageArray.filter(c => c >= 2).length;

  return {
    totalBases: templateLength,
    coveredBases: covered,
    coveragePercent: ((covered / templateLength) * 100).toFixed(1) + '%',
    doubleCoveredBases: doubleCovered,
    doubleCoveragePercent: ((doubleCovered / templateLength) * 100).toFixed(1) + '%',
    gaps: findGaps(coverageArray),
  };
}

function findGaps(coverageArray) {
  const gaps = [];
  let gapStart = null;

  for (let i = 0; i < coverageArray.length; i++) {
    if (coverageArray[i] === 0 && gapStart === null) {
      gapStart = i;
    } else if (coverageArray[i] > 0 && gapStart !== null) {
      gaps.push({ start: gapStart, end: i - 1, length: i - gapStart });
      gapStart = null;
    }
  }

  if (gapStart !== null) {
    gaps.push({ start: gapStart, end: coverageArray.length - 1, length: coverageArray.length - gapStart });
  }

  return gaps;
}

// =============================================================================
// Backward Compatibility Exports
// =============================================================================

/**
 * @deprecated Use calculateHairpinDG from equilibrium.js instead
 * Wrapper for backward compatibility - returns simplified hairpin analysis
 */
export function detectHairpins(seq) {
  const hairpinDG = calculateHairpinDG(seq, 55);
  const hasHairpin = hairpinDG < -1.0;
  const hasProblem = hairpinDG < -2.0;

  return {
    hairpins: hasHairpin ? [{ dG: hairpinDG, is3prime: true }] : [],
    worstHairpin: hasHairpin ? { dG: hairpinDG } : null,
    worst3primeHairpin: hasHairpin ? { dG: hairpinDG } : null,
    hasProblem,
    // Note: This is a simplified wrapper. For full analysis, use analyzeSinglePrimer()
    _deprecated: true,
    _useInstead: 'calculateHairpinDG from equilibrium.js or analyzeSinglePrimer from primerAnalysis.js',
  };
}

/**
 * @deprecated Use calculateHomodimerDG from equilibrium.js instead
 * Wrapper for backward compatibility - returns simplified self-dimer analysis
 */
export function detectSelfDimer(seq) {
  const homodimerDG = calculateHomodimerDG(seq, 55);
  const hasDimer = homodimerDG < -3.0;
  const hasProblem = homodimerDG < -5.0;

  return {
    dimers: hasDimer ? [{ dG: homodimerDG, is3prime: true }] : [],
    worstDimer: hasDimer ? { dG: homodimerDG } : null,
    worst3primeDimer: hasDimer ? { dG: homodimerDG } : null,
    hasProblem,
    // Note: This is a simplified wrapper. For full analysis, use analyzeSinglePrimer()
    _deprecated: true,
    _useInstead: 'calculateHomodimerDG from equilibrium.js or analyzeSinglePrimer from primerAnalysis.js',
  };
}

// Re-export from other modules for convenience
export { reverseComplement };
export { calculateHairpinDG, calculateHomodimerDG } from './equilibrium.js';
export { analyzeSinglePrimer } from './primerAnalysis.js';
export { ANALYSIS_PRESETS } from './presets.js';
