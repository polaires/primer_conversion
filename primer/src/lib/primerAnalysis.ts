/**
 * Unified Primer Analysis Module
 *
 * Provides a single entry point for all primer analysis across design modes.
 * This module coordinates between:
 * - scoring.js (piecewise scoring functions)
 * - equilibrium.js (thermodynamic calculations)
 * - offTargetClassification.js (off-target analysis)
 * - smartPrimers.js (3' end analysis)
 *
 * Key features:
 * - Domain-specific presets (amplification, mutagenesis, sequencing, assembly)
 * - Unified analysis output format
 * - Backward-compatible with existing modules
 * - Extensible for future analysis types
 */

import {
  scoreTm,
  scoreGc,
  scoreLength,
  scoreGcClamp,
  scoreHomopolymer,
  scoreHairpin,
  scoreHomodimer,
  scoreHeterodimer,
  scoreTerminal3DG,
  scoreTmDiff,
  scoreOffTarget,
  scoreGQuadruplex,
  analyzeGQuadruplex,
  score3PrimeComposition,
  calculateCompositeScore,
  classifyQuality,
  buildCompositeInput,
  DEFAULT_WEIGHTS,
  type GQuadruplexAnalysis,
  type CompositeScoreResult,
  type ScoreBreakdownItem,
} from './scoring.js';

import {
  calculateHairpinDG,
  calculateHomodimerDG,
  calculateHeterodimerDG,
} from './equilibrium.js';

import { calculateTmQ5, calculateGC, calculate3primeTerminalDG } from './tmQ5.js';
import { tm } from './tm.js';
import { analyze3PrimeEnd, classify3PrimeStructureSeverity } from './smartPrimers.js';
import { ANALYSIS_PRESETS } from './presets.js';
import { reverseComplement } from './sequenceUtils.js';

// Re-export for backward compatibility
export { ANALYSIS_PRESETS };

// =============================================================================
// Type Definitions
// =============================================================================

export type AnalysisMode = 'amplification' | 'mutagenesis' | 'sequencing' | 'assembly' | 'goldengate';
export type WarningSeverity = 'ok' | 'info' | 'caution' | 'warning' | 'critical';
export type WarningType =
  | 'templateMismatch'
  | 'annealingTm'
  | 'mismatch'
  | 'hairpin'
  | 'homodimer'
  | 'heterodimer'
  | 'gQuadruplex'
  | 'terminal3DG'
  | '3primeEnd'
  | 'length'
  | 'homopolymer'
  | 'tmDiff'
  | 'goldenGateDetected';

export type Strand = 'forward' | 'reverse';
export type QualityTier = 'excellent' | 'good' | 'acceptable' | 'marginal' | 'poor';

export interface GoldenGateEnzyme {
  recognition: string;
  name: string;
  alias: string | null;
}

export interface GoldenGateSite {
  enzyme: string;
  recognition: string;
  position: number;
  orientation: 'forward' | 'reverse';
  alias: string | null;
}

export interface GoldenGateDetectionResult {
  hasGoldenGateSite: boolean;
  detectedSites: GoldenGateSite[];
  primaryEnzyme: string | null;
  isLikelyGoldenGatePrimer: boolean;
  recommendation: string | null;
}

export interface AnnealingRegion {
  sequence: string;
  length: number;
  position: number;
  strand: Strand;
  tailLength: number;
  tailSequence: string;
  mismatches: number;
  identity: number;
  tm: number | null;
  gc: number | null;
  gcPercent: string;
  withinOptimalTm: boolean;
  terminal3DG?: number;
}

export interface FuzzyMatch {
  position: number;
  strand: Strand;
  mismatches: number;
}

export interface TemplateBinding {
  found: boolean;
  warnings: string[];
}

export interface PrimerWarning {
  type: WarningType;
  severity: WarningSeverity;
  message: string;
  primer?: 'forward' | 'reverse' | 'pair';
  value?: number;
  dG?: number;
  issues?: string[];
  enzyme?: string;
  sites?: GoldenGateSite[];
}

export interface PrimerScores {
  tm: number;
  gc: number;
  length: number;
  gcClamp: number;
  homopolymer: number;
  hairpin: number;
  homodimer: number;
  terminal3DG: number;
  offTarget: number;
  gQuadruplex: number;
  threePrimeComp: number;
}

export interface AssemblyPrimerScores {
  annealingTm: number;
  annealingGc: number;
  annealingLength: number;
  annealingTerminal3DG: number;
  annealingGcClamp: number;
  fullPrimerHairpin: number;
  fullPrimerHomodimer: number;
  homopolymer: number;
  gQuadruplex?: number;
}

export interface Thermodynamics {
  hairpinDG: number;
  homodimerDG: number;
  terminal3DG: number;
}

export interface ThreePrimeAnalysis {
  quality: QualityTier;
  issues: string[];
  gcClamp: number;
  terminalDG: number;
  hasProblematicPattern: boolean;
}

export interface AnnealingRegionInfo {
  sequence: string;
  length: number;
  tm: number;
  gc: number;
  gcPercent: string;
  position: number;
  strand: Strand;
  terminal3DG: number;
}

export interface TailRegionInfo {
  sequence: string;
  length: number;
}

export interface FullPrimerInfo {
  sequence: string;
  length: number;
  tm: number;
  gc: number;
  gcPercent: string;
  hairpinDG: number;
  homodimerDG: number;
}

export interface AssemblyPrimerAnalysis {
  fullPrimer: FullPrimerInfo;
  annealingRegion: AnnealingRegion | null;
  templateBinding: TemplateBinding;
  isAssemblyPrimer: boolean;
  scores: Partial<AssemblyPrimerScores>;
  warnings: PrimerWarning[];
  gQuadruplex?: GQuadruplexAnalysis;
}

export interface SinglePrimerAnalysis {
  sequence: string;
  length: number;
  tm: number;
  fullPrimerTm?: number;
  gc: number;
  gcPercent: string;
  thermodynamics: Thermodynamics;
  scores: PrimerScores;
  gQuadruplex: GQuadruplexAnalysis | null;
  analysis3Prime: ThreePrimeAnalysis | null;
  warnings: PrimerWarning[];
  goldenGateDetection?: GoldenGateDetectionResult | null;
  isAssemblyPrimer?: boolean;
  annealingRegion?: AnnealingRegionInfo;
  tailRegion?: TailRegionInfo | null;
  templateBinding?: TemplateBinding;
}

export interface PairScores {
  tmDiff: number;
  heterodimer: number;
}

export interface PairInfo {
  tmDiff: number;
  heterodimerDG: number;
  ampliconLength: number | null;
  scores: PairScores;
}

export interface QualityClassification {
  tier: QualityTier;
  description: string;
  color: string;
}

export interface PrimerPairAnalysis {
  mode: AnalysisMode;
  preset: string;
  forward: SinglePrimerAnalysis;
  reverse: SinglePrimerAnalysis | null;
  pair: PairInfo | null;
  composite: CompositeScoreResult;
  quality: QualityClassification;
  effectiveScore?: number;
  criticalWarnings?: number;
  warnings: PrimerWarning[];
  recommendations: string[];
  isAcceptable?: boolean;
}

export interface QuickCheckResult {
  sequence: string;
  length: number;
  gc: number;
  gcClamp: number;
  quality: QualityTier;
  issues: string[];
  passesQuickCheck: boolean;
}

export interface PrimerInput {
  seq: string;
  tm?: number;
  gc?: number;
  offTargetCount?: number;
}

export interface FindAnnealingRegionOptions {
  minLength?: number;
  maxLength?: number;
  allowMismatch?: boolean;
  maxMismatches?: number;
}

export interface AnalyzeAssemblyPrimerOptions {
  targetTm?: number;
  temperature?: number;
  minLength?: number;
  maxLength?: number;
  allowMismatch?: boolean;
  maxMismatches?: number;
}

export interface AnalyzeSinglePrimerOptions {
  mode?: AnalysisMode;
  template?: string | null;
  temperature?: number;
  include3PrimeAnalysis?: boolean;
  includeGQuadruplex?: boolean;
}

export interface AnalyzePrimersOptions {
  mode?: AnalysisMode;
  template?: string | null;
  temperature?: number;
  ampliconLength?: number | null;
  include3PrimeAnalysis?: boolean;
  includeGQuadruplex?: boolean;
  includeRecommendations?: boolean;
  heterodimerDG?: number | null;
  structureWarnings?: PrimerWarning[] | null;
}

export interface QuickPrimerCheckOptions {
  mode?: AnalysisMode;
}

export interface MutagenesisAnalysisResult {
  forward: {
    sequence: string;
    length: number;
    tm: number;
    gc: number;
    gcPercent: string;
    hairpinDG: number;
    selfDimerDG: number;
    foldDG: number;
  };
  reverse: {
    sequence: string;
    length: number;
    tm: number;
    gc: number;
    gcPercent: string;
    hairpinDG: number;
    selfDimerDG: number;
    foldDG: number;
  };
  pair: {
    tmDifference: number;
    annealingTemp: number;
    heterodimerDG: number;
  };
  offTargets: null;
  warnings: PrimerWarning[];
  quality: QualityTier;
  isAcceptable: boolean;
  compositeScore: number;
  gQuadruplex: {
    forward: GQuadruplexAnalysis | null;
    reverse: GQuadruplexAnalysis | null;
  };
}

// =============================================================================
// Assembly Primer Annealing Region Detection
// =============================================================================
//
// NAMING CONVENTIONS:
// - "annealing region" = portion of primer that binds to template (3' end)
// - "tail region" = non-binding portion (5' end, contains enzyme sites, overhangs)
// - "full primer" = complete primer sequence
//
// In Golden Gate context:
//   Primer structure: [flanking]-[enzyme]-[spacer]-[overhang]-[homology]
//   - "homology" in GG === "annealing region" here
//   - Everything before homology === "tail region" here
//
// In Gibson/NEBuilder context:
//   Primer structure: [homology tail]-[annealing region]
//   - The 5' homology tail provides overlap between fragments
//   - The 3' annealing region binds template for PCR amplification
//
// SCORING APPROACH:
// - Tm, GC, length scores based on annealing region (determines PCR efficiency)
// - Secondary structure scores on full primer (hairpin/dimer affects all of it)
//
// RELATED FILES:
// - goldengate-primer-optimizer.js: GG-specific scoring with enzyme site validation
// - nebuilder.js: NEBuilder primer design with overlap Tm optimization
// - assemblyCore.js: Shared assembly primer design logic
// =============================================================================

/**
 * Constants for assembly primer analysis
 * Based on NEB recommendations and peer-reviewed literature
 */
const ASSEMBLY_CONSTANTS = {
  MIN_ANNEALING_LENGTH: 15,  // Minimum binding region (NEB recommendation)
  MAX_ANNEALING_LENGTH: 35,  // Maximum to search for
  TARGET_ANNEALING_TM: 60,   // Target Tm for Q5 (NEB recommends 60-72°C)
  MIN_ANNEALING_TM: 50,      // Minimum acceptable Tm
  MAX_ANNEALING_TM: 72,      // Maximum acceptable Tm
  MIN_MATCH_IDENTITY: 0.9,   // 90% identity for fuzzy matching
};

/**
 * Golden Gate enzyme recognition sites for auto-detection
 * These are Type IIS restriction enzymes commonly used in Golden Gate assembly
 */
const GOLDEN_GATE_ENZYMES: Record<string, GoldenGateEnzyme> = {
  BsaI:  { recognition: 'GGTCTC', name: 'BsaI',  alias: null },
  BsmBI: { recognition: 'CGTCTC', name: 'BsmBI', alias: 'Esp3I' },
  BbsI:  { recognition: 'GAAGAC', name: 'BbsI',  alias: null },
  SapI:  { recognition: 'GCTCTTC', name: 'SapI', alias: 'BspQI' },
};

/**
 * Detect Golden Gate enzyme recognition sites in a primer sequence.
 *
 * This function scans a primer for common Type IIS restriction enzyme sites
 * used in Golden Gate assembly (BsaI, BsmBI, BbsI, SapI, etc.).
 *
 * @param primerSeq - Primer sequence to analyze
 * @returns Detection result with found sites and recommendations
 */
export function detectGoldenGateSites(primerSeq: string): GoldenGateDetectionResult {
  const primer = primerSeq.toUpperCase();
  const detectedSites: GoldenGateSite[] = [];

  for (const [enzymeName, enzyme] of Object.entries(GOLDEN_GATE_ENZYMES)) {
    const recognition = enzyme.recognition;
    const rcRecognition = reverseComplement(recognition);

    // Check forward orientation
    let pos = primer.indexOf(recognition);
    if (pos !== -1) {
      detectedSites.push({
        enzyme: enzymeName,
        recognition,
        position: pos,
        orientation: 'forward',
        alias: enzyme.alias,
      });
    }

    // Check reverse complement (for reverse primers)
    pos = primer.indexOf(rcRecognition);
    if (pos !== -1 && rcRecognition !== recognition) {
      detectedSites.push({
        enzyme: enzymeName,
        recognition: rcRecognition,
        position: pos,
        orientation: 'reverse',
        alias: enzyme.alias,
      });
    }
  }

  const hasGoldenGateSite = detectedSites.length > 0;
  const primaryEnzyme = detectedSites.length > 0 ? detectedSites[0].enzyme : null;

  return {
    hasGoldenGateSite,
    detectedSites,
    primaryEnzyme,
    isLikelyGoldenGatePrimer: hasGoldenGateSite && primer.length > 25,
    recommendation: hasGoldenGateSite
      ? 'Golden Gate enzyme site detected. Provide a template sequence for accurate annealing region analysis.'
      : null,
  };
}

/**
 * Find the annealing region of a primer on a template sequence.
 *
 * For assembly primers (Golden Gate, Gibson, NEBuilder), the primer structure is:
 *   [5' tail/flanking]-[enzyme site]-[spacer]-[overhang]-[annealing region]
 *
 * Only the 3' annealing region actually binds to the template during PCR.
 * This function identifies that region by finding where the primer matches the template.
 *
 * @param primerSeq - Full primer sequence
 * @param templateSeq - Template sequence to search
 * @param options - Search options
 * @returns Annealing region info or null if not found
 */
export function findAnnealingRegion(
  primerSeq: string,
  templateSeq: string,
  options: FindAnnealingRegionOptions = {}
): AnnealingRegion | null {
  const {
    minLength = ASSEMBLY_CONSTANTS.MIN_ANNEALING_LENGTH,
    maxLength = ASSEMBLY_CONSTANTS.MAX_ANNEALING_LENGTH,
    allowMismatch = true,
    maxMismatches = 2,
  } = options;

  const primer = primerSeq.toUpperCase();
  const template = templateSeq.toUpperCase();
  const templateRC = reverseComplement(template);

  // Strategy: Search from 3' end of primer (longest to shortest)
  // The 3' end is what actually anneals to the template

  let bestMatch: AnnealingRegion | null = null;

  for (let len = Math.min(maxLength, primer.length); len >= minLength; len--) {
    const annealingCandidate = primer.slice(-len);

    // Search on forward strand
    let matchPos = template.indexOf(annealingCandidate);
    if (matchPos !== -1) {
      bestMatch = {
        sequence: annealingCandidate,
        length: len,
        position: matchPos,
        strand: 'forward',
        tailLength: primer.length - len,
        tailSequence: primer.slice(0, primer.length - len),
        mismatches: 0,
        identity: 1.0,
        tm: null,
        gc: null,
        gcPercent: '',
        withinOptimalTm: false,
      };
      break;
    }

    // Search on reverse complement strand (for reverse primers)
    matchPos = templateRC.indexOf(annealingCandidate);
    if (matchPos !== -1) {
      bestMatch = {
        sequence: annealingCandidate,
        length: len,
        position: template.length - matchPos - len,
        strand: 'reverse',
        tailLength: primer.length - len,
        tailSequence: primer.slice(0, primer.length - len),
        mismatches: 0,
        identity: 1.0,
        tm: null,
        gc: null,
        gcPercent: '',
        withinOptimalTm: false,
      };
      break;
    }

    // Try fuzzy matching if allowed (for primers with minor mismatches)
    if (allowMismatch && len >= 18) {
      const fuzzyMatch = findFuzzyMatch(annealingCandidate, template, templateRC, maxMismatches);
      if (fuzzyMatch) {
        bestMatch = {
          sequence: annealingCandidate,
          length: len,
          position: fuzzyMatch.position,
          strand: fuzzyMatch.strand,
          tailLength: primer.length - len,
          tailSequence: primer.slice(0, primer.length - len),
          mismatches: fuzzyMatch.mismatches,
          identity: 1 - (fuzzyMatch.mismatches / len),
          tm: null,
          gc: null,
          gcPercent: '',
          withinOptimalTm: false,
        };
        break;
      }
    }
  }

  if (!bestMatch) {
    return null;
  }

  // Calculate Tm and GC for the annealing region
  try {
    bestMatch.tm = calculateTmQ5(bestMatch.sequence);
    bestMatch.gc = calculateGC(bestMatch.sequence);
    bestMatch.gcPercent = `${Math.round(bestMatch.gc * 100)}%`;
    bestMatch.withinOptimalTm = bestMatch.tm >= ASSEMBLY_CONSTANTS.MIN_ANNEALING_TM &&
                                bestMatch.tm <= ASSEMBLY_CONSTANTS.MAX_ANNEALING_TM;
  } catch (e) {
    bestMatch.tm = null;
    bestMatch.gc = null;
    bestMatch.withinOptimalTm = false;
  }

  return bestMatch;
}

/**
 * Find a fuzzy match allowing some mismatches
 * @private
 */
function findFuzzyMatch(
  query: string,
  template: string,
  templateRC: string,
  maxMismatches: number
): FuzzyMatch | null {
  const queryLen = query.length;

  // Search forward strand
  for (let i = 0; i <= template.length - queryLen; i++) {
    const segment = template.slice(i, i + queryLen);
    const mismatches = countMismatches(query, segment);
    if (mismatches <= maxMismatches) {
      return { position: i, strand: 'forward', mismatches };
    }
  }

  // Search reverse strand
  for (let i = 0; i <= templateRC.length - queryLen; i++) {
    const segment = templateRC.slice(i, i + queryLen);
    const mismatches = countMismatches(query, segment);
    if (mismatches <= maxMismatches) {
      return {
        position: template.length - i - queryLen,
        strand: 'reverse',
        mismatches
      };
    }
  }

  return null;
}

/**
 * Count mismatches between two sequences
 * @private
 */
function countMismatches(seq1: string, seq2: string): number {
  let mismatches = 0;
  for (let i = 0; i < seq1.length; i++) {
    if (seq1[i] !== seq2[i]) mismatches++;
  }
  return mismatches;
}

/**
 * Analyze an assembly primer with template context
 *
 * This function provides specialized analysis for Golden Gate, Gibson, and NEBuilder
 * assembly primers. It:
 * 1. Identifies the annealing region (portion that binds template)
 * 2. Scores the annealing region for PCR efficiency
 * 3. Checks the full primer for secondary structure issues
 * 4. Validates template binding
 *
 * @param primerSeq - Full primer sequence
 * @param templateSeq - Template sequence
 * @param options - Analysis options
 * @returns Comprehensive assembly primer analysis
 */
export function analyzeAssemblyPrimer(
  primerSeq: string,
  templateSeq: string,
  options: AnalyzeAssemblyPrimerOptions = {}
): AssemblyPrimerAnalysis {
  const {
    targetTm = ASSEMBLY_CONSTANTS.TARGET_ANNEALING_TM,
    temperature = 55,
  } = options;

  const primer = primerSeq.toUpperCase();
  const template = templateSeq.toUpperCase();

  // Find the annealing region
  const annealingRegion = findAnnealingRegion(primer, template, options);

  // Calculate full primer properties (for secondary structure)
  const fullPrimerTm = calculateTmQ5(primer);
  const fullPrimerGc = calculateGC(primer);
  const fullHairpinDG = calculateHairpinDG(primer, temperature);
  const fullHomodimerDG = calculateHomodimerDG(primer, temperature);

  // Build analysis result
  const analysis: AssemblyPrimerAnalysis = {
    fullPrimer: {
      sequence: primer,
      length: primer.length,
      tm: Math.round(fullPrimerTm * 10) / 10,
      gc: Math.round(fullPrimerGc * 100) / 100,
      gcPercent: `${Math.round(fullPrimerGc * 100)}%`,
      hairpinDG: Math.round(fullHairpinDG * 100) / 100,
      homodimerDG: Math.round(fullHomodimerDG * 100) / 100,
    },
    annealingRegion: null,
    templateBinding: {
      found: false,
      warnings: [],
    },
    isAssemblyPrimer: false,
    scores: {},
    warnings: [],
  };

  if (!annealingRegion) {
    // Primer doesn't match template - might be wrong template or fully synthetic
    analysis.templateBinding.warnings.push(
      'Primer does not match provided template - cannot identify annealing region'
    );
    analysis.warnings.push({
      type: 'templateMismatch',
      severity: 'warning',
      message: 'Primer sequence not found in template. Using full primer for analysis.',
    });
    return analysis;
  }

  // We found the annealing region
  analysis.annealingRegion = annealingRegion;
  analysis.templateBinding.found = true;
  analysis.isAssemblyPrimer = annealingRegion.tailLength > 0;

  // Calculate annealing region 3' terminal ΔG
  const annealingTerminalDG = calculate3primeTerminalDG(annealingRegion.sequence).dG;

  // Add annealing-specific properties
  analysis.annealingRegion.terminal3DG = Math.round(annealingTerminalDG * 100) / 100;

  // Generate scores based on annealing region (what matters for PCR)
  const preset = ANALYSIS_PRESETS.assembly;

  analysis.scores = {
    // Score annealing region properties (these determine PCR success)
    annealingTm: scoreTm(annealingRegion.tm!, {
      optimalLow: preset.optimalTmRange[0],
      optimalHigh: preset.optimalTmRange[1],
    }),
    annealingGc: scoreGc(annealingRegion.gc!, {
      optimalLow: preset.optimalGcRange[0],
      optimalHigh: preset.optimalGcRange[1],
    }),
    annealingLength: scoreLength(annealingRegion.length, {
      optimalLow: preset.optimalLengthRange[0],
      optimalHigh: preset.optimalLengthRange[1],
    }),
    annealingTerminal3DG: scoreTerminal3DG(annealingTerminalDG),
    annealingGcClamp: scoreGcClamp(annealingRegion.sequence),

    // Score full primer secondary structure (affects primer-primer interactions)
    fullPrimerHairpin: scoreHairpin(fullHairpinDG, { threshold: preset.thresholds.hairpin }),
    fullPrimerHomodimer: scoreHomodimer(fullHomodimerDG, { threshold: preset.thresholds.homodimer }),

    // Homopolymer and G-quadruplex checked on full primer
    homopolymer: scoreHomopolymer(primer),
  };

  // Check for G-quadruplex in full primer
  const gQuadruplexAnalysis = analyzeGQuadruplex(primer);
  analysis.scores.gQuadruplex = gQuadruplexAnalysis?.score ?? 1.0;
  analysis.gQuadruplex = gQuadruplexAnalysis;

  // Generate warnings
  if (!annealingRegion.withinOptimalTm) {
    const severity: WarningSeverity = annealingRegion.tm! < ASSEMBLY_CONSTANTS.MIN_ANNEALING_TM ? 'warning' : 'info';
    analysis.warnings.push({
      type: 'annealingTm',
      severity,
      value: annealingRegion.tm!,
      message: `Annealing region Tm (${annealingRegion.tm!.toFixed(1)}°C) outside optimal range (${ASSEMBLY_CONSTANTS.MIN_ANNEALING_TM}-${ASSEMBLY_CONSTANTS.MAX_ANNEALING_TM}°C)`,
    });
  }

  if (annealingRegion.mismatches > 0) {
    analysis.warnings.push({
      type: 'mismatch',
      severity: annealingRegion.mismatches > 1 ? 'warning' : 'info',
      value: annealingRegion.mismatches,
      message: `Primer has ${annealingRegion.mismatches} mismatch(es) with template`,
    });
  }

  if (fullHairpinDG < preset.thresholds.hairpin) {
    analysis.warnings.push({
      type: 'hairpin',
      severity: fullHairpinDG < -6.0 ? 'critical' : 'warning',
      dG: fullHairpinDG,
      message: `Full primer may form stable hairpin (ΔG = ${fullHairpinDG.toFixed(1)} kcal/mol)`,
    });
  }

  if (fullHomodimerDG < preset.thresholds.homodimer) {
    analysis.warnings.push({
      type: 'homodimer',
      severity: fullHomodimerDG < -9.0 ? 'critical' : 'warning',
      dG: fullHomodimerDG,
      message: `Full primer may form stable self-dimer (ΔG = ${fullHomodimerDG.toFixed(1)} kcal/mol)`,
    });
  }

  if (gQuadruplexAnalysis?.severity === 'critical' || gQuadruplexAnalysis?.severity === 'warning') {
    analysis.warnings.push({
      type: 'gQuadruplex',
      severity: gQuadruplexAnalysis.severity,
      message: gQuadruplexAnalysis.message,
    });
  }

  if (analysis.scores.annealingTerminal3DG! < 0.5) {
    analysis.warnings.push({
      type: 'terminal3DG',
      severity: 'warning',
      dG: annealingTerminalDG,
      message: `Weak 3' terminal binding in annealing region (ΔG = ${annealingTerminalDG.toFixed(1)} kcal/mol)`,
    });
  }

  return analysis;
}

// ANALYSIS_PRESETS is now imported from presets.js (single source of truth)

/**
 * Calculate Tm using the appropriate method
 *
 * @param seq - Primer sequence
 * @param method - 'santaLucia' | 'q5' | 'auto'
 * @returns Melting temperature in °C
 */
export function calculateTm(seq: string, method: 'santaLucia' | 'q5' | 'auto' = 'auto'): number {
  if (method === 'q5') {
    return calculateTmQ5(seq);
  }
  if (method === 'santaLucia') {
    return tm(seq, '', true);
  }
  // Auto: use Q5 for longer primers (mutagenesis), SantaLucia for shorter
  return seq.length > 30 ? calculateTmQ5(seq) : tm(seq, '', true);
}

/**
 * Analyze a single primer comprehensively
 *
 * For assembly modes (Golden Gate, Gibson, NEBuilder) with a template provided,
 * this function identifies the annealing region and scores based on that portion,
 * while still checking the full primer for secondary structure issues.
 *
 * @param primer - Primer object or sequence string
 * @param options - Analysis options
 * @returns Comprehensive single primer analysis
 */
export function analyzeSinglePrimer(
  primer: PrimerInput | string,
  options: AnalyzeSinglePrimerOptions = {}
): SinglePrimerAnalysis {
  const {
    mode = 'amplification',
    template = null,
    temperature = 55,
    include3PrimeAnalysis = true,
    includeGQuadruplex = true,
  } = options;

  const preset = ANALYSIS_PRESETS[mode] || ANALYSIS_PRESETS.amplification;
  const seq = typeof primer === 'string' ? primer : primer.seq;
  const seqUpper = seq.toUpperCase();
  const offTargetCount = typeof primer === 'object' ? (primer.offTargetCount ?? 0) : 0;

  // Check if we're in an assembly-related mode (for annealing region detection)
  const isAssemblyMode = mode === 'assembly' || mode === 'goldengate';

  // ==========================================================================
  // ASSEMBLY MODE WITH TEMPLATE: Use specialized analysis
  // ==========================================================================
  // For assembly primers (Golden Gate, Gibson, NEBuilder), the full primer
  // includes tails/enzyme sites that don't anneal to the template.
  // We need to find the annealing region and score based on that.
  // ==========================================================================
  if (isAssemblyMode && template) {
    const assemblyAnalysis = analyzeAssemblyPrimer(seqUpper, template, { temperature });

    if (assemblyAnalysis.templateBinding.found && assemblyAnalysis.annealingRegion) {
      const annealingRegion = assemblyAnalysis.annealingRegion;

      // Use annealing region properties for scoring (what matters for PCR)
      const annealingTm = annealingRegion.tm!;
      const annealingGc = annealingRegion.gc!;
      const annealingTerminalDG = annealingRegion.terminal3DG!;

      // Full primer properties (for secondary structure)
      const fullPrimer = assemblyAnalysis.fullPrimer;

      // 3' end analysis on annealing region (since that's what extends)
      const analysis3Prime = include3PrimeAnalysis ? analyze3PrimeEnd(annealingRegion.sequence) : null;

      // Calculate scores using annealing region for Tm/GC/length, full primer for structure
      const scores: PrimerScores = {
        // Annealing region scores (determine PCR success)
        tm: scoreTm(annealingTm, {
          optimalLow: preset.optimalTmRange[0],
          optimalHigh: preset.optimalTmRange[1],
        }),
        gc: scoreGc(annealingGc, {
          optimalLow: preset.optimalGcRange[0],
          optimalHigh: preset.optimalGcRange[1],
        }),
        length: scoreLength(annealingRegion.length, {
          optimalLow: preset.optimalLengthRange[0],
          optimalHigh: preset.optimalLengthRange[1],
        }),
        gcClamp: scoreGcClamp(annealingRegion.sequence),
        terminal3DG: scoreTerminal3DG(annealingTerminalDG),
        threePrimeComp: score3PrimeComposition(annealingRegion.sequence, annealingTerminalDG),

        // Full primer scores (secondary structure affects entire primer)
        hairpin: assemblyAnalysis.scores.fullPrimerHairpin!,
        homodimer: assemblyAnalysis.scores.fullPrimerHomodimer!,
        homopolymer: assemblyAnalysis.scores.homopolymer!,
        gQuadruplex: assemblyAnalysis.scores.gQuadruplex ?? 1.0,
        offTarget: scoreOffTarget(offTargetCount),
      };

      // Collect warnings from assembly analysis
      const warnings = [...assemblyAnalysis.warnings];

      // Add 3' end warning if applicable
      if (analysis3Prime?.quality === 'poor') {
        warnings.push({
          type: '3primeEnd',
          severity: 'warning',
          issues: analysis3Prime.issues,
          message: `Annealing region 3' end quality issues: ${analysis3Prime.issues.join(', ')}`,
        });
      }

      return {
        sequence: seqUpper,
        length: seqUpper.length,
        // Report both full primer and annealing Tm for clarity
        tm: Math.round(annealingTm * 10) / 10,
        fullPrimerTm: Math.round(fullPrimer.tm * 10) / 10,
        gc: Math.round(annealingGc * 100) / 100,
        gcPercent: annealingRegion.gcPercent,
        thermodynamics: {
          hairpinDG: fullPrimer.hairpinDG,
          homodimerDG: fullPrimer.homodimerDG,
          terminal3DG: annealingTerminalDG,
        },
        scores,
        gQuadruplex: assemblyAnalysis.gQuadruplex ?? null,
        analysis3Prime,
        warnings,
        // Assembly-specific fields
        isAssemblyPrimer: true,
        annealingRegion: {
          sequence: annealingRegion.sequence,
          length: annealingRegion.length,
          tm: Math.round(annealingTm * 10) / 10,
          gc: Math.round(annealingGc * 100) / 100,
          gcPercent: annealingRegion.gcPercent,
          position: annealingRegion.position,
          strand: annealingRegion.strand,
          terminal3DG: annealingTerminalDG,
        },
        tailRegion: annealingRegion.tailLength > 0 ? {
          sequence: annealingRegion.tailSequence,
          length: annealingRegion.tailLength,
        } : null,
        templateBinding: assemblyAnalysis.templateBinding,
      };
    }
    // If template binding not found, fall through to standard analysis
    // but include the warning
  }

  // ==========================================================================
  // STANDARD ANALYSIS (non-assembly modes or no template)
  // ==========================================================================

  // Calculate basic properties if not provided
  const primerTm = typeof primer === 'object' && primer.tm !== undefined
    ? primer.tm
    : calculateTm(seqUpper, preset.tmMethod);
  const primerGc = typeof primer === 'object' && primer.gc !== undefined
    ? primer.gc
    : calculateGC(seqUpper);

  // Calculate thermodynamic values
  const hairpinDG = calculateHairpinDG(seqUpper, temperature);
  const homodimerDG = calculateHomodimerDG(seqUpper, temperature);

  // Calculate 3' terminal ΔG (actual binding strength, not self-folding)
  const terminalDG = calculate3primeTerminalDG(seqUpper).dG;

  // 3' end analysis
  const analysis3Prime = include3PrimeAnalysis ? analyze3PrimeEnd(seqUpper) : null;

  // G-quadruplex analysis
  const gQuadruplexAnalysis = includeGQuadruplex ? analyzeGQuadruplex(seqUpper) : null;

  // Calculate individual scores
  const scores: PrimerScores = {
    tm: scoreTm(primerTm, {
      optimalLow: preset.optimalTmRange[0],
      optimalHigh: preset.optimalTmRange[1],
    }),
    gc: scoreGc(primerGc, {
      optimalLow: preset.optimalGcRange[0],
      optimalHigh: preset.optimalGcRange[1],
    }),
    length: scoreLength(seqUpper.length, {
      optimalLow: preset.optimalLengthRange[0],
      optimalHigh: preset.optimalLengthRange[1],
    }),
    gcClamp: scoreGcClamp(seqUpper),
    homopolymer: scoreHomopolymer(seqUpper),
    hairpin: scoreHairpin(hairpinDG, { threshold: preset.thresholds.hairpin }),
    homodimer: scoreHomodimer(homodimerDG, { threshold: preset.thresholds.homodimer }),
    terminal3DG: scoreTerminal3DG(terminalDG),
    offTarget: scoreOffTarget(offTargetCount),
    gQuadruplex: gQuadruplexAnalysis?.score ?? 1.0,
    // 3' end composition score (weighted analysis of GC clamp, terminal ΔG, patterns)
    threePrimeComp: score3PrimeComposition(seqUpper, terminalDG),
  };

  // Collect warnings
  const warnings: PrimerWarning[] = [];

  if (gQuadruplexAnalysis?.severity === 'critical') {
    warnings.push({
      type: 'gQuadruplex',
      severity: 'critical',
      message: gQuadruplexAnalysis.message,
    });
  } else if (gQuadruplexAnalysis?.severity === 'warning') {
    warnings.push({
      type: 'gQuadruplex',
      severity: 'warning',
      message: gQuadruplexAnalysis.message,
    });
  }

  if (hairpinDG < preset.thresholds.hairpin) {
    const severity: WarningSeverity = hairpinDG < (preset.thresholds.hairpinCritical ?? -6.0) ? 'critical' : 'warning';
    warnings.push({
      type: 'hairpin',
      severity,
      dG: hairpinDG,
      message: `Stable hairpin detected (ΔG = ${hairpinDG.toFixed(1)} kcal/mol)`,
    });
  }

  if (homodimerDG < preset.thresholds.homodimer) {
    const severity: WarningSeverity = homodimerDG < (preset.thresholds.homodimerCritical ?? -9.0) ? 'critical' : 'warning';
    warnings.push({
      type: 'homodimer',
      severity,
      dG: homodimerDG,
      message: `Stable self-dimer detected (ΔG = ${homodimerDG.toFixed(1)} kcal/mol)`,
    });
  }

  if (scores.terminal3DG < 0.5) {
    warnings.push({
      type: 'terminal3DG',
      severity: 'warning',
      dG: terminalDG,
      message: `Weak 3' terminal binding (ΔG = ${terminalDG.toFixed(1)} kcal/mol)`,
    });
  }

  if (analysis3Prime?.quality === 'poor') {
    warnings.push({
      type: '3primeEnd',
      severity: 'warning',
      issues: analysis3Prime.issues,
      message: `3' end quality issues: ${analysis3Prime.issues.join(', ')}`,
    });
  }

  // Warning for extreme length violations
  const lengthLow = preset.optimalLengthRange[0];
  const lengthHigh = preset.optimalLengthRange[1];
  if (seqUpper.length > lengthHigh + 15) {
    // Primer is extremely long (e.g., 50bp for PCR optimal 18-24)
    warnings.push({
      type: 'length',
      severity: 'critical',
      value: seqUpper.length,
      message: `Primer extremely long: ${seqUpper.length}bp (optimal: ${lengthLow}-${lengthHigh}bp)`,
    });
  } else if (seqUpper.length > lengthHigh + 6 || seqUpper.length < lengthLow - 3) {
    // Primer is significantly outside optimal range
    warnings.push({
      type: 'length',
      severity: 'warning',
      value: seqUpper.length,
      message: `Primer length outside optimal range: ${seqUpper.length}bp (optimal: ${lengthLow}-${lengthHigh}bp)`,
    });
  }

  // Warning for severe homopolymer runs (5+ identical bases)
  let maxHomopolymerRun = 1;
  let currentRun = 1;
  for (let i = 1; i < seqUpper.length; i++) {
    if (seqUpper[i] === seqUpper[i - 1]) {
      currentRun++;
      maxHomopolymerRun = Math.max(maxHomopolymerRun, currentRun);
    } else {
      currentRun = 1;
    }
  }
  if (maxHomopolymerRun >= 6) {
    warnings.push({
      type: 'homopolymer',
      severity: 'critical',
      value: maxHomopolymerRun,
      message: `Severe homopolymer run: ${maxHomopolymerRun} identical bases (polymerase slippage risk)`,
    });
  } else if (maxHomopolymerRun >= 5) {
    warnings.push({
      type: 'homopolymer',
      severity: 'warning',
      value: maxHomopolymerRun,
      message: `Homopolymer run: ${maxHomopolymerRun} identical bases (potential polymerase slippage)`,
    });
  }

  // Auto-detect Golden Gate enzyme sites and suggest template/mode switch if found
  const ggDetection = detectGoldenGateSites(seqUpper);
  if (ggDetection.hasGoldenGateSite) {
    // Check if user is NOT in goldengate mode (even assembly mode should warn)
    if (mode !== 'goldengate') {
      warnings.push({
        type: 'goldenGateDetected',
        severity: 'warning',
        enzyme: ggDetection.primaryEnzyme!,
        sites: ggDetection.detectedSites,
        message: `${ggDetection.primaryEnzyme} site detected but using ${mode} mode. Switch to Golden Gate mode for accurate scoring.`,
      });
    } else if (!template) {
      // In goldengate mode but no template
      warnings.push({
        type: 'goldenGateDetected',
        severity: 'info',
        enzyme: ggDetection.primaryEnzyme!,
        sites: ggDetection.detectedSites,
        message: `${ggDetection.primaryEnzyme} site detected. Add a template sequence for accurate annealing region Tm analysis.`,
      });
    }
  }

  return {
    sequence: seqUpper,
    length: seqUpper.length,
    tm: Math.round(primerTm * 10) / 10,
    gc: Math.round(primerGc * 100) / 100,
    gcPercent: `${Math.round(primerGc * 100)}%`,
    thermodynamics: {
      hairpinDG: Math.round(hairpinDG * 100) / 100,
      homodimerDG: Math.round(homodimerDG * 100) / 100,
      terminal3DG: Math.round(terminalDG * 100) / 100,
    },
    scores,
    gQuadruplex: gQuadruplexAnalysis,
    analysis3Prime,
    warnings,
    // Include GG detection info for UI to display
    goldenGateDetection: ggDetection.hasGoldenGateSite ? ggDetection : null,
  };
}

/**
 * Unified primer pair analysis function
 *
 * This is the main entry point for comprehensive primer pair analysis.
 * It provides consistent analysis across all design modes while allowing
 * domain-specific customization through presets.
 *
 * @param fwd - Forward primer (object with seq/tm/gc/dg or sequence string)
 * @param rev - Reverse primer (optional)
 * @param options - Analysis options
 * @returns Comprehensive primer pair analysis
 */
export function analyzePrimers(
  fwd: PrimerInput | string,
  rev: PrimerInput | string | null = null,
  options: AnalyzePrimersOptions = {}
): PrimerPairAnalysis {
  const {
    mode = 'amplification',
    template = null,
    temperature = 55,
    ampliconLength = null,
    include3PrimeAnalysis = true,
    includeGQuadruplex = true,
    includeRecommendations = true,
    // Allow passing pre-calculated thermodynamic values
    heterodimerDG = null,
    // Domain-specific warnings (from mutagenesis checkMutantSecondaryStructure)
    structureWarnings = null,
  } = options;

  const preset = ANALYSIS_PRESETS[mode] || ANALYSIS_PRESETS.amplification;

  // Analyze forward primer
  const fwdAnalysis = analyzeSinglePrimer(fwd, {
    mode,
    template,
    temperature,
    include3PrimeAnalysis,
    includeGQuadruplex,
  });

  // If no reverse primer, return single primer analysis
  if (!rev) {
    const singleScores = {
      tmFwd: fwdAnalysis.scores.tm,
      gcFwd: fwdAnalysis.scores.gc,
      lengthFwd: fwdAnalysis.scores.length,
      hairpinFwd: fwdAnalysis.scores.hairpin,
      selfDimerFwd: fwdAnalysis.scores.homodimer,
      terminal3DG: fwdAnalysis.scores.terminal3DG,
      offTarget: fwdAnalysis.scores.offTarget,
      gQuadruplexFwd: fwdAnalysis.scores.gQuadruplex,
      // Include homopolymer and GC clamp scores
      homopolymerFwd: fwdAnalysis.scores.homopolymer,
      gcClampFwd: fwdAnalysis.scores.gcClamp,
      // Include 3' end composition score
      threePrimeCompFwd: fwdAnalysis.scores.threePrimeComp,
    };

    const composite = calculateCompositeScore(singleScores, preset.weights);
    const quality = classifyQuality(composite.score);

    return {
      mode,
      preset: preset.name,
      forward: fwdAnalysis,
      reverse: null,
      pair: null,
      composite,
      quality,
      warnings: fwdAnalysis.warnings,
      recommendations: includeRecommendations ? generateRecommendations(fwdAnalysis, null, preset) : [],
    };
  }

  // Analyze reverse primer
  const revAnalysis = analyzeSinglePrimer(rev, {
    mode,
    template,
    temperature,
    include3PrimeAnalysis,
    includeGQuadruplex,
  });

  // Calculate pair-level thermodynamics
  const fwdSeq = typeof fwd === 'string' ? fwd.toUpperCase() : fwd.seq.toUpperCase();
  const revSeq = typeof rev === 'string' ? rev.toUpperCase() : rev.seq.toUpperCase();
  const calculatedHeterodimerDG = heterodimerDG ?? calculateHeterodimerDG(fwdSeq, revSeq, temperature);

  // Pair-level scores
  const tmDiff = Math.abs(fwdAnalysis.tm - revAnalysis.tm);
  const pairScores: PairScores = {
    tmDiff: scoreTmDiff(fwdAnalysis.tm, revAnalysis.tm),
    heterodimer: scoreHeterodimer(calculatedHeterodimerDG, { threshold: preset.thresholds.heterodimer }),
  };

  // Collect all warnings
  const allWarnings: PrimerWarning[] = [
    ...fwdAnalysis.warnings.map(w => ({ ...w, primer: 'forward' as const })),
    ...revAnalysis.warnings.map(w => ({ ...w, primer: 'reverse' as const })),
  ];

  // Add pair-level warnings
  if (tmDiff > 5) {
    allWarnings.push({
      type: 'tmDiff',
      primer: 'pair',
      severity: tmDiff > 8 ? 'critical' : 'warning',
      value: tmDiff,
      message: `Tm difference: ${tmDiff.toFixed(1)}°C (>5°C may cause unequal amplification)`,
    });
  }

  if (calculatedHeterodimerDG < preset.thresholds.heterodimer) {
    const severity: WarningSeverity = calculatedHeterodimerDG < (preset.thresholds.heterodimerCritical ?? -9.0) ? 'critical' : 'warning';
    allWarnings.push({
      type: 'heterodimer',
      primer: 'pair',
      severity,
      dG: calculatedHeterodimerDG,
      message: `Stable heterodimer detected (ΔG = ${calculatedHeterodimerDG.toFixed(1)} kcal/mol)`,
    });
  }

  // Add domain-specific warnings
  if (structureWarnings && Array.isArray(structureWarnings)) {
    allWarnings.push(...structureWarnings);
  }

  // Calculate composite score using unified buildCompositeInput
  const compositeInput = buildCompositeInput(fwdAnalysis.scores, revAnalysis.scores, pairScores);

  const composite = calculateCompositeScore(compositeInput, preset.weights);
  const quality = classifyQuality(composite.score);

  // Determine overall quality based on warnings
  let effectiveQuality = quality;
  const criticalWarnings = allWarnings.filter(w => w.severity === 'critical').length;
  const effectiveScore = Math.max(0, composite.score - criticalWarnings * 20);
  if (criticalWarnings > 0 && quality.tier !== 'poor') {
    effectiveQuality = classifyQuality(effectiveScore);
  }

  return {
    mode,
    preset: preset.name,
    forward: fwdAnalysis,
    reverse: revAnalysis,
    pair: {
      tmDiff: Math.round(tmDiff * 10) / 10,
      heterodimerDG: Math.round(calculatedHeterodimerDG * 100) / 100,
      ampliconLength,
      scores: pairScores,
    },
    composite,
    quality: effectiveQuality,
    // Effective score after critical warning penalties (-20 per critical warning)
    effectiveScore,
    criticalWarnings,
    warnings: allWarnings,
    recommendations: includeRecommendations ? generateRecommendations(fwdAnalysis, revAnalysis, preset) : [],
    // For backward compatibility
    isAcceptable: effectiveQuality.tier !== 'poor',
  };
}

/**
 * Generate recommendations based on analysis results
 *
 * @param fwdAnalysis - Forward primer analysis
 * @param revAnalysis - Reverse primer analysis (optional)
 * @param preset - Analysis preset
 * @returns Array of recommendation strings
 */
function generateRecommendations(
  fwdAnalysis: SinglePrimerAnalysis,
  revAnalysis: SinglePrimerAnalysis | null,
  preset: any
): string[] {
  const recommendations: string[] = [];

  // Forward primer recommendations
  if (fwdAnalysis.scores.gQuadruplex < 0.5) {
    recommendations.push('Redesign forward primer to avoid G-quadruplex motif');
  }
  if (fwdAnalysis.scores.terminal3DG < 0.5) {
    recommendations.push('Consider extending forward primer for stronger 3\' end');
  }
  if (fwdAnalysis.scores.hairpin < 0.5) {
    recommendations.push('Forward primer may form stable hairpin - consider adjusting length/position');
  }

  // Reverse primer recommendations
  if (revAnalysis) {
    if (revAnalysis.scores.gQuadruplex < 0.5) {
      recommendations.push('Redesign reverse primer to avoid G-quadruplex motif');
    }
    if (revAnalysis.scores.terminal3DG < 0.5) {
      recommendations.push('Consider extending reverse primer for stronger 3\' end');
    }
    if (revAnalysis.scores.hairpin < 0.5) {
      recommendations.push('Reverse primer may form stable hairpin - consider adjusting length/position');
    }

    // Pair recommendations
    const tmDiff = Math.abs(fwdAnalysis.tm - revAnalysis.tm);
    if (tmDiff > 5) {
      recommendations.push('Adjust primer lengths to match Tm values');
    }
  }

  return recommendations;
}

/**
 * Adapter function for mutagenesis module
 *
 * Transforms the unified analysis output to match the existing mutagenesis
 * analyzePrimerPair() return format for backward compatibility.
 *
 * @param fwdPrimer - Forward primer sequence
 * @param revPrimer - Reverse primer sequence
 * @param template - Template sequence (optional)
 * @param options - Analysis options
 * @returns Analysis in mutagenesis format
 */
export function analyzePrimerPairForMutagenesis(
  fwdPrimer: string,
  revPrimer: string,
  template: string | null = null,
  options: Partial<AnalyzePrimersOptions> = {}
): MutagenesisAnalysisResult {
  const result = analyzePrimers(
    { seq: fwdPrimer },
    { seq: revPrimer },
    {
      mode: 'mutagenesis',
      template,
      ...options,
    }
  );

  // Transform to existing mutagenesis format
  return {
    forward: {
      sequence: result.forward.sequence,
      length: result.forward.length,
      tm: result.forward.tm,
      gc: result.forward.gc,
      gcPercent: result.forward.gcPercent,
      hairpinDG: result.forward.thermodynamics.hairpinDG,
      selfDimerDG: result.forward.thermodynamics.homodimerDG,
      foldDG: result.forward.thermodynamics.hairpinDG,  // Backward compatibility
    },
    reverse: {
      sequence: result.reverse!.sequence,
      length: result.reverse!.length,
      tm: result.reverse!.tm,
      gc: result.reverse!.gc,
      gcPercent: result.reverse!.gcPercent,
      hairpinDG: result.reverse!.thermodynamics.hairpinDG,
      selfDimerDG: result.reverse!.thermodynamics.homodimerDG,
      foldDG: result.reverse!.thermodynamics.hairpinDG,  // Backward compatibility
    },
    pair: {
      tmDifference: result.pair!.tmDiff,
      annealingTemp: Math.min(result.forward.tm, result.reverse!.tm) - 5,
      heterodimerDG: result.pair!.heterodimerDG,
    },
    offTargets: null,  // Calculated separately in mutagenesis
    warnings: result.warnings,
    quality: result.quality.tier,
    isAcceptable: result.isAcceptable!,
    // Extended analysis (new)
    compositeScore: result.composite.score,
    gQuadruplex: {
      forward: result.forward.gQuadruplex,
      reverse: result.reverse!.gQuadruplex,
    },
  };
}

/**
 * Quick primer quality check
 *
 * Performs a fast quality assessment without full thermodynamic calculations.
 * Useful for filtering candidates before detailed analysis.
 *
 * @param seq - Primer sequence
 * @param options - Quick check options
 * @returns Quick quality assessment
 */
export function quickPrimerCheck(seq: string, options: QuickPrimerCheckOptions = {}): QuickCheckResult {
  const {
    mode = 'amplification',
  } = options;

  const preset = ANALYSIS_PRESETS[mode] || ANALYSIS_PRESETS.amplification;
  const seqUpper = seq.toUpperCase();
  const issues: string[] = [];

  // Length check
  if (seqUpper.length < preset.optimalLengthRange[0]) {
    issues.push(`Length too short: ${seqUpper.length}bp`);
  } else if (seqUpper.length > preset.optimalLengthRange[1] + 10) {
    issues.push(`Length too long: ${seqUpper.length}bp`);
  }

  // GC content check
  const gc = calculateGC(seqUpper);
  if (gc < 0.30) {
    issues.push(`GC content too low: ${Math.round(gc * 100)}%`);
  } else if (gc > 0.75) {
    issues.push(`GC content too high: ${Math.round(gc * 100)}%`);
  }

  // GC clamp check
  const last2 = seqUpper.slice(-2);
  const gcLast2 = (last2.match(/[GC]/g) || []).length;
  if (gcLast2 === 0) {
    issues.push('No GC clamp (0 G/C in last 2 bases)');
  }

  // Homopolymer check
  if (/(.)\1{3,}/.test(seqUpper)) {
    issues.push('Homopolymer run (4+ identical bases)');
  }

  // G-quadruplex quick check
  if (/G{4,}/.test(seqUpper)) {
    issues.push('GGGG run detected - G-quadruplex risk');
  }

  const qualityLevel: QualityTier =
    issues.length === 0 ? 'excellent' :
    issues.length === 1 ? 'good' :
    issues.length === 2 ? 'acceptable' :
    issues.length === 3 ? 'marginal' : 'poor';

  return {
    sequence: seqUpper,
    length: seqUpper.length,
    gc: Math.round(gc * 100),
    gcClamp: gcLast2,
    quality: qualityLevel,
    issues,
    passesQuickCheck: issues.length <= 2,
  };
}

// Re-export utilities for convenience
export {
  ANALYSIS_PRESETS as presets,
  calculateCompositeScore,
  classifyQuality,
  DEFAULT_WEIGHTS,
  // Assembly primer analysis (for direct use by assembly modules)
  // Note: findAnnealingRegion, analyzeAssemblyPrimer, detectGoldenGateSites
  // are already exported at their function definitions
};
