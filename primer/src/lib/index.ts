/**
 * Primers library - PCR primer design tool
 * JavaScript port of https://github.com/Lattice-Automation/primers
 *
 * Now supports DNA24 thermodynamic parameters (Greenleaf Lab, 2024)
 * for improved accuracy in secondary structure and mismatch predictions.
 */

export { primers, create, score, generateAlternatives, LEN_MIN, LEN_MAX, LEN_MAX_EXTENDED } from "./primers.js";
export { gcCache, tmCache, tm, setParameterSet } from "./tm.js";
export {
  calculateTmQ5,
  calculateAnnealingQ5,
  calculateTmQ5Detailed,
  calculateGC,
  batchCalculateTmQ5,
  Q5_DEFAULTS,
} from "./tmQ5.js";
export { dgCache, dg, fold, setFoldParameterSet } from "./fold.js";
export { offTargets } from "./offTargets.js";

// DNA24 parameter configuration
import { setParameterSet } from "./tm.js";
import { setFoldParameterSet } from "./fold.js";

/**
 * Configure which thermodynamic parameter set to use for calculations.
 *
 * DNA24 (default): Greenleaf Lab 2024 parameters from high-throughput
 * array melt experiments. Provides ~35% better accuracy for secondary
 * structure predictions and ~50% better mismatch predictions.
 *
 * SantaLucia 1998: Classic parameters widely used in primer design tools.
 *
 * @param useDna24 - true for DNA24 (recommended), false for SantaLucia 1998
 *
 * @example
 * // Use DNA24 parameters (default, recommended)
 * setThermodynamicParams(true);
 *
 * @example
 * // Use legacy SantaLucia 1998 parameters
 * setThermodynamicParams(false);
 */
export function setThermodynamicParams(useDna24: boolean = true): void {
  setParameterSet(useDna24);
  setFoldParameterSet(useDna24);
}

// Sequencing primer design
export {
  designSequencingPrimers,
  designPrimerAtPosition,
  scorePrimerQuality,
  SEQUENCING_DEFAULTS,
} from "./sequencing.js";

// Smart primer design with iterative optimization
export {
  analyze3PrimeEnd,
  optimizePrimer,
  optimizePrimerPair,
  scorePrimerVariant,
  generateDesignSuggestions,
  quickAssess,
  generateLengthVariants,
  score3PrimeComposition,
  SMART_DESIGN_CONSTRAINTS,
} from "./smartPrimers.js";

// Scoring functions
export {
  scoreTm,
  scoreGc,
  scoreTerminal3DG,
  scoreTmDiff,
  scoreHairpin,
  scoreHomodimer,
  scoreHeterodimer,
  scoreOffTarget,
  scoreLength,
  scoreGcClamp,
  scoreHomopolymer,
  score3PrimeComposition as score3PrimeComp,
  calculateCompositeScore,
  classifyQuality,
  piecewiseLogistic,
  DEFAULT_WEIGHTS,
  SCORING_PRESETS,
  analyzePrimerPair as analyzeWithScoring,  // Canonical scoring analysis
} from "./scoring.js";

// Unified primer analysis (single entry point for all design modes)
export {
  analyzePrimers,
  analyzeSinglePrimer,
  analyzePrimerPairForMutagenesis,
  quickPrimerCheck,
  calculateTm as calculateTmUnified,
  ANALYSIS_PRESETS,
} from "./primerAnalysis.js";

// Unified presets (single source of truth for all modes)
export { getPreset } from "./presets.js";

// Site-directed mutagenesis (base changer)
export {
  designSubstitutionPrimers,
  designCodonChangePrimers,
  designInsertionPrimers,
  designDeletionPrimers,
  designRegionSubstitutionPrimers,
  designPrimersFromNotation,
  parseMutationNotation,
  MUTATION_TYPES,
  MUTAGENESIS_DEFAULTS,
  CODON_TABLE,
  CODON_TO_AA,
  AA_NAMES,
  // Enhanced features
  designPrimers,
  analyzePrimerPair,
  batchDesignMutations,
  compareTmMethods,
  calculateMismatchedTm,
  checkMutantSecondaryStructure,
  checkHeterodimer,
  checkPrimerSpecificity,
  selectOptimalCodon,
  checkGQuadruplexRisk,
  score3primeTerminalBase,
  scorePrimerPair,
  // Thermodynamic parameters
  NN_MATCHED,
  NN_MISMATCH,
  DANGLING_END_CORRECTIONS,
  TERMINAL_CORRECTIONS,
  CONSECUTIVE_MISMATCH_CORRECTION,
  // Evidence-based dimer thresholds (IDT, Premier Biosoft references)
  DIMER_THRESHOLDS,
  getThreshold,
  classifyDimerSeverity,
} from "./mutagenesis.js";
