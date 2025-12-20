/**
 * Enhanced Site-Directed Mutagenesis Module (Combined with PCR Primer Design)
 *
 * Combines NEB Base Changer algorithm with Q5 SDM kit design principles and
 * integrates PCR primer design functionality.
 *
 * Key features:
 * - Mismatch Tm calculation using Allawi/SantaLucia parameters with full nearest-neighbor thermodynamics
 * - Back-to-back primer design (non-overlapping for Q5 SDM)
 * - QuikChange-style overlapping primers (traditional)
 * - Optimal codon selection (minimum changes + codon usage)
 * - Off-target binding site detection
 * - Secondary structure analysis for mutant primers (using Zuker algorithm)
 * - PCR primer design integration
 * - Dangling end effects
 * - Consecutive mismatch handling
 * - Terminal vs internal mismatch differentiation
 *
 * References:
 * - SantaLucia J Jr. (1998) PNAS 95:1460-5
 * - SantaLucia J Jr. & Hicks D (2004) Annu Rev Biophys Biomol Struct 33:415-40
 * - Allawi HT, SantaLucia J Jr. (1997) Biochemistry 36:10581-94 (G·T mismatches)
 * - Allawi HT, SantaLucia J Jr. (1998) Biochemistry 37:2170-9 (G·A, C·T, A·C)
 * - Allawi HT, SantaLucia J Jr. (1998) Nucleic Acids Res 26:2694-701 (internal loops)
 * - Peyret N et al. (1999) Biochemistry 38:3468-77 (tandem mismatches)
 * - Bommarito S et al. (2000) Nucleic Acids Res 28:1929-34 (terminal/dangling effects)
 * - Owczarzy R et al. (2008) Biochemistry 47:5336-53 (salt corrections)
 * - Q5 Site-Directed Mutagenesis Kit (NEB #E0554)
 * - Vallone PM, Butler JM (2004) AutoDimer algorithm
 */

import { calculateTmQ5, calculateAnnealingQ5, calculateGC, findOptimalLengthForTm, getValidPrimerCandidates, findBestTmMatchedPair, calculate3primeTerminalDG } from './tmQ5.js';
import { dg, fold } from './fold.js';
import { offTargets } from './offTargets.js';
import { DNA_INTERNAL_MM, DNA_TERMINAL_MM, DNA_DE, DNA_NN } from './dna';
import { tm as calculateTmGeneral, tmCache, gcCache } from './tm.js';
// Calibrated scoring system (Döring dataset, F1=81.9%, AUC=0.848)
import {
  scoreTm, scoreGc, scoreTerminal3DG, scoreHairpin, scoreHomodimer,
  scoreHeterodimer, scoreTmDiff, scoreOffTarget, scoreLength, scoreGcClamp,
  scoreHomopolymer, scoreGQuadruplex, analyzeGQuadruplex,
  score3PrimeComposition,  // 3' end composition scoring (8% total weight)
  calculateCompositeScore, classifyQuality,
} from './scoring.js';
import { calculateHairpinDG, calculateHomodimerDG, calculateHeterodimerDG } from './equilibrium.js';
import { classify3PrimeStructureSeverity } from './smartPrimers.js';
// Unified analysis for consistent scoring across modules
import { analyzePrimers as unifiedAnalyzePrimers } from './primerAnalysis.js';
// Centralized thermodynamic constants
import {
  NN_MATCHED,
  NN_MISMATCH,
  TERMINAL_CORRECTIONS,
  DANGLING_END_CORRECTIONS,
  CONSECUTIVE_MISMATCH_CORRECTION,
} from './thermoConstants.js';
import { reverseComplement } from './sequenceUtils.js';

// =============================================================================
// Type Definitions
// =============================================================================

export interface ThresholdRange {
  ideal: number;
  warning: number;
  critical: number;
}

export interface DimerThresholds {
  hairpin: {
    threePrime: ThresholdRange;
    internal: ThresholdRange;
  };
  selfDimer: {
    threePrime: ThresholdRange;
    internal: ThresholdRange;
  };
  heterodimer: {
    threePrime: ThresholdRange;
    internal: ThresholdRange;
  };
  threePrimeExtensible: {
    minConsecutiveBp: number;
    minOverlapBp: number;
  };
}

export interface MutagenesisDefaults {
  minAnnealingLength: number;
  maxAnnealingLength: number;
  minPrimerLength: number;
  maxPrimerLength: number;
  minFlankingLength: number;
  maxFlankingLength: number;
  minTm: number;
  maxTm: number;
  minGC: number;
  maxGC: number;
  maxPolyN: number;
  minDg: number;
  gcClampRequired: boolean;
  checkOffTargets: boolean;
  maxOffTargets: number;
  strategy: 'back-to-back' | 'overlapping';
  organism: 'ecoli' | 'human' | null;
  circular: boolean;
  confineTo5Tails: boolean;
}

export interface MutationTypes {
  SUBSTITUTION: 'substitution';
  INSERTION: 'insertion';
  DELETION: 'deletion';
  CODON_CHANGE: 'codon_change';
}

export interface Mismatch {
  position: number;
  primerBase: string;
  templateBase: string;
  type: string;
  isTerminal: boolean;
  is3prime: boolean;
}

export interface MismatchedTmResult {
  tm: number | null;
  willNotBind: boolean;
  mismatchFraction: number;
  dH: number;
  dS: number;
  mismatches: Mismatch[];
  mismatchCount: number;
  consecutiveMismatchCount: number;
  maxConsecutiveMismatches: number;
  hasCritical3primeMismatch: boolean;
  hasTerminalMismatch: boolean;
  has5primeMismatch: boolean;
  has3primeMismatch: boolean;
  thermodynamics: {
    dH: number;
    dS: number;
    saltCorrection: number;
  };
}

export interface FoldStructure {
  e: number;
  desc?: string;
  ij?: [number, number][];
}

export interface SecondaryStructureWarning {
  type: string;
  dG?: number;
  severity: 'critical' | 'warning' | 'info';
  involves3Prime?: boolean;
  message: string;
  tooltip?: string;
  details?: any;
  severityLevel?: string;
}

export interface SecondaryStructureCheck {
  hairpinDG: number;
  selfDimerDG: number;
  foldDG: number;
  structures: FoldStructure[];
  warnings: SecondaryStructureWarning[];
  isAcceptable: boolean;
  changeFromOriginal: {
    hairpinChange: number;
    selfDimerChange: number;
    foldChange: number;
    madeWorse: boolean;
  } | null;
}

export interface AlignedBase {
  pos1: number;
  pos2: number;
  base1: string;
  base2: string;
}

export interface HeterodimerAlignment {
  offset: number;
  alignedBases: AlignedBase[];
  complementaryRegionLength: number;
  maxConsecutive: number;
  involves3Prime1: boolean;
  involves3Prime2: boolean;
  involves5Prime1: boolean;
  involves5Prime2: boolean;
}

export interface HeterodimerResult {
  heterodimerDG: number;
  alignment: HeterodimerAlignment | null;
  warnings: SecondaryStructureWarning[];
  isAcceptable: boolean;
  isExpectedOverlap?: boolean;
  overlapLength?: number;
  overlapPercent?: number;
  severity: 'safe' | 'warning' | 'critical';
  dimerType: string;
  involves3Prime1?: boolean;
  involves3Prime2?: boolean;
  maxConsecutive?: number;
  message?: string;
}

export interface CodonCandidate {
  codon: string;
  changes: number;
  positions: number[];
  usage: number;
  score: number;
}

export interface CodonSelection {
  selectedCodon: string;
  nucleotideChanges: number;
  codonUsage: number;
  allCandidates: CodonCandidate[];
  changedPositions: number[];
}

export interface BindingSite {
  position: number;
  strand: 'forward' | 'reverse';
  mismatches: number;
  sequence: string;
}

export interface SpecificityResult {
  offTargetCount: number;
  forwardSites: number;
  reverseSites: number;
  bindingSites: BindingSite[];
  isSpecific: boolean;
}

export interface PrimerInfo {
  sequence: string;
  length: number;
  tm: number;
  gc: number;
  gcPercent: string;
  hasGCClamp: boolean;
  start: number;
  end: number;
  dg?: number;
  hairpinDG?: number;
  selfDimerDG?: number;
  foldDG?: number;
  tmMismatch?: number | null;
  mismatchCount?: number;
  willNotBind?: boolean;
  mismatchFraction?: number;
  offTargets?: number;
  terminalDG?: { dG: number; sequence: string };
  terminalBase?: { base: string; penalty: number; classification: string; hasGCClamp: boolean };
  gQuadruplex?: { hasG4Motif: boolean; hasGGGG: boolean; score: number };
  isRescue?: boolean;
  isDeletion?: boolean;
  scores?: any;
}

export interface CandidatePair {
  forward: PrimerInfo;
  reverse: PrimerInfo;
  design: 'back-to-back' | 'overlapping';
  tmDiff?: number;
  penalty: number;
  leftFlank?: number;
  rightFlank?: number;
  splitPointOffset?: number;
  contextIssues?: any[];
  circularWrapped?: boolean;
  warnings?: SecondaryStructureWarning[];
  compositeScore?: number;
  piecewiseScores?: Record<string, number>;
  qualityTier?: string;
  tierQuality?: string;
}

export interface ProtocolStep {
  name: string;
  temp?: string;
  time?: string;
  substeps?: ProtocolStep[];
}

export interface Protocol {
  name: string;
  steps: ProtocolStep[];
  notes: string[];
}

export interface MutagenesisResult {
  type: string;
  originalSequence: string;
  mutatedSequence: string;
  position: number;
  change: string;
  description: string;
  forward: PrimerInfo;
  reverse: PrimerInfo;
  design: string;
  annealingTemp: number;
  tmDifference: number;
  quality: string;
  penalty: number;
  compositeScore?: number;
  piecewiseScores?: Record<string, number>;
  qualityTier?: string;
  structureCheck: SecondaryStructureCheck;
  splitPointOffset?: number;
  alternateDesigns?: CandidatePair[];
  crossStrategyAlternates?: CandidatePair[];
  protocol: Protocol;
  warnings?: SecondaryStructureWarning[];
  // Type-specific fields
  oldBase?: string;
  newBase?: string;
  nucleotidePosition?: number;
  oldCodon?: string;
  newCodon?: string;
  oldAA?: string;
  newAA?: string;
  codonChanges?: number;
  codonUsage?: number;
  alternativeCodons?: CodonCandidate[];
  insertedSequence?: string;
  insertLength?: number;
  deletedSequence?: string;
  deleteLength?: number;
  replacementSequence?: string;
  replacementLength?: number;
  analysis?: PrimerPairAnalysis;
}

export interface ParsedMutation {
  type: string;
  oldAA?: string;
  position: number;
  newAA?: string;
  notation: string;
  oldBase?: string;
  newBase?: string;
  length?: number;
  sequence?: string;
}

export interface GQuadruplexRisk {
  hasG4Motif: boolean;
  hasGGGG: boolean;
  gggCount: number;
  gggRuns: number[];
  severity: string;
  message: string;
  score: number;
}

export interface TerminalBaseScore {
  base: string;
  penalty: number;
  classification: string;
  hasGCClamp: boolean;
}

export interface PrimerScoreBreakdown {
  tmDiff: { value: number; penalty: number };
  tmDistance: { penalty: number };
  length: { fwd: number; rev: number; penalty: number };
  terminalDG: {
    fwd: { dG: number; sequence: string };
    rev: { dG: number; sequence: string };
    penalty: number;
  };
  terminalBase: {
    fwd: TerminalBaseScore;
    rev: TerminalBaseScore;
    penalty: number;
  };
  gQuadruplex: {
    fwd: GQuadruplexRisk;
    rev: GQuadruplexRisk;
    penalty: number;
  };
}

export interface PrimerScoreResult {
  score: number;
  notes: string[];
  breakdown: PrimerScoreBreakdown;
  fwd: PrimerInfo;
  rev: PrimerInfo;
  summary: {
    tmDiff: number;
    avgTm: number;
    totalLength: number;
    hasG4Risk: boolean;
    hasGCClamp: boolean;
  };
  compositeScore: number | null;
  qualityTier: string | null;
  piecewiseScores: Record<string, number> | null;
}

export interface PrimerPairAnalysis {
  forward: PrimerInfo;
  reverse: PrimerInfo;
  pair: {
    tmDifference: number;
    annealingTemp: number;
    heterodimerDG: number;
    scores?: any;
  };
  offTargets: {
    forward: SpecificityResult;
    reverse: SpecificityResult;
  } | null;
  warnings: SecondaryStructureWarning[];
  quality: string;
  isAcceptable: boolean;
  compositeScore?: number;
  qualityTier?: string;
  gQuadruplex?: {
    forward: GQuadruplexRisk;
    reverse: GQuadruplexRisk;
  };
}

export interface BatchResult {
  mutation: string;
  success: boolean;
  result?: MutagenesisResult;
  error?: string;
}

export interface BatchDesignResult {
  results: BatchResult[];
  successCount: number;
  failureCount: number;
}

export interface TmComparisonResult {
  sequence: string;
  length: number;
  gcContent: number;
  methods: {
    q5: number;
    general: number;
    mismatch?: number | null;
    simple: number;
  };
  mismatchDetails?: MismatchedTmResult;
}

// =============================================================================
// Evidence-Based Secondary Structure Thresholds
// =============================================================================
/**
 * Industry-standard thresholds for primer secondary structure analysis.
 *
 * IMPORTANT: 3' end involvement is MORE CRITICAL than internal structures because:
 * - 3' end is where DNA polymerase extends from
 * - Even weak 3' structures can cause primer-dimer artifacts
 * - Hot start polymerases do NOT prevent 3' extensible dimers
 *
 * References:
 * -----------
 * 1. IDT (Integrated DNA Technologies) OligoAnalyzer Guidelines:
 *    https://www.idtdna.com/pages/support/faqs/how-can-i-check-my-pcr-primers-using-the-oligoanalyzer-program-to-ensure-there-are-no-significant-primer-design-issues-
 *    - Self-dimers/heterodimers: 3' end ≥ -6 kcal/mol, internal ≥ -9 kcal/mol
 *    - "IDT has found that primers with dimers at -9 kcal/mol can cause problems"
 *
 * 2. Premier Biosoft Primer Design Technical Notes:
 *    http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
 *    - Hairpin: 3' end ≥ -2 kcal/mol, internal ≥ -3 kcal/mol
 *    - Self-dimer: 3' end ≥ -5 kcal/mol, internal ≥ -6 kcal/mol
 *    - Cross-dimer: 3' end ≥ -5 kcal/mol, internal ≥ -6 kcal/mol
 *
 * 3. Hot Start Limitations (Nucleosides, Nucleotides & Nucleic Acids, 2020):
 *    https://pubmed.ncbi.nlm.nih.gov/32799617/
 *    - "Even hot-start DNA polymerases do not prevent PD formation if primers
 *       have stable 3'-overlapping"
 *    - ">2 overlapping nucleotides at 3' cause considerable primer dimer accumulation"
 *
 * 4. Commercial Hot Start Variability (BioTechniques, 2013):
 *    https://www.tandfonline.com/doi/full/10.2144/000114481
 *    - 7-12 of 17 tested hot start polymerases showed activity before thermal activation
 *    - Hot start is not a reliable solution for poor primer design
 */
export const DIMER_THRESHOLDS: DimerThresholds = {
  // Hairpin thresholds (self-folding of single primer)
  hairpin: {
    threePrime: {
      ideal: -2.0,      // Premier Biosoft: "3' end hairpin with ΔG of -2 kcal/mol tolerated"
      warning: -3.0,    // Approaching problematic
      critical: -4.0,   // Will likely cause issues
    },
    internal: {
      ideal: -3.0,      // Premier Biosoft: "internal hairpin with ΔG of -3 kcal/mol tolerated"
      warning: -5.0,    // Monitor for issues
      critical: -6.0,   // IDT: significant problems at this level
    },
  },

  // Self-dimer thresholds (primer binding to itself)
  selfDimer: {
    threePrime: {
      ideal: -5.0,      // Premier Biosoft recommendation
      warning: -6.0,    // IDT 3' threshold
      critical: -8.0,   // Severe - will form stable dimers
    },
    internal: {
      ideal: -6.0,      // Premier Biosoft recommendation
      warning: -8.0,    // Approaching IDT critical
      critical: -9.0,   // IDT: "significant problems at -9 kcal/mol"
    },
  },

  // Heterodimer/Cross-dimer thresholds (two primers binding)
  heterodimer: {
    threePrime: {
      ideal: -5.0,      // Premier Biosoft recommendation
      warning: -6.0,    // IDT 3' threshold
      critical: -8.0,   // Severe binding at 3' ends
    },
    internal: {
      ideal: -6.0,      // Premier Biosoft recommendation
      warning: -8.0,    // Approaching critical
      critical: -9.0,   // IDT threshold
    },
  },

  // 3' Extensible Dimer - ALWAYS CRITICAL regardless of ΔG
  // This is when both primers' 3' ends are involved with sufficient overlap
  threePrimeExtensible: {
    minConsecutiveBp: 4,  // 4+ consecutive base pairs = extensible
    minOverlapBp: 3,      // Literature: >2 bp overlap causes problems
    // These are ALWAYS critical - polymerase can extend from hybridized 3' ends
  },
};

/**
 * Get appropriate threshold based on position involvement
 * @param structureType - 'hairpin', 'selfDimer', or 'heterodimer'
 * @param involves3Prime - Whether 3' end is involved
 * @returns Threshold object with ideal, warning, critical values
 */
export function getThreshold(structureType: string, involves3Prime: boolean = false): ThresholdRange {
  const thresholds = DIMER_THRESHOLDS[structureType as keyof typeof DIMER_THRESHOLDS];
  if (!thresholds || typeof thresholds === 'object' && 'minConsecutiveBp' in thresholds) {
    console.warn(`Unknown structure type: ${structureType}`);
    return { ideal: -3.0, warning: -6.0, critical: -9.0 };
  }
  return involves3Prime ? thresholds.threePrime : thresholds.internal;
}

/**
 * Classify severity based on evidence-based thresholds
 * @param deltaG - The ΔG value in kcal/mol
 * @param structureType - 'hairpin', 'selfDimer', or 'heterodimer'
 * @param involves3Prime - Whether 3' end is involved
 * @returns 'pass', 'warning', or 'fail'
 */
export function classifyDimerSeverity(deltaG: number, structureType: string, involves3Prime: boolean = false): 'pass' | 'warning' | 'fail' {
  const threshold = getThreshold(structureType, involves3Prime);

  if (deltaG >= threshold.ideal) return 'pass';
  if (deltaG >= threshold.warning) return 'warning';
  return 'fail';
}

// =============================================================================
// Codon Tables and Usage
// =============================================================================

const CODON_TABLE: Record<string, string[]> = {
  'F': ['TTT', 'TTC'],
  'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
  'I': ['ATT', 'ATC', 'ATA'],
  'M': ['ATG'],
  'V': ['GTT', 'GTC', 'GTA', 'GTG'],
  'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
  'P': ['CCT', 'CCC', 'CCA', 'CCG'],
  'T': ['ACT', 'ACC', 'ACA', 'ACG'],
  'A': ['GCT', 'GCC', 'GCA', 'GCG'],
  'Y': ['TAT', 'TAC'],
  'H': ['CAT', 'CAC'],
  'Q': ['CAA', 'CAG'],
  'N': ['AAT', 'AAC'],
  'K': ['AAA', 'AAG'],
  'D': ['GAT', 'GAC'],
  'E': ['GAA', 'GAG'],
  'C': ['TGT', 'TGC'],
  'W': ['TGG'],
  'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
  'G': ['GGT', 'GGC', 'GGA', 'GGG'],
  '*': ['TAA', 'TAG', 'TGA'],
};

const CODON_TO_AA: Record<string, string> = {};
for (const [aa, codons] of Object.entries(CODON_TABLE)) {
  for (const codon of codons) {
    CODON_TO_AA[codon] = aa;
  }
}

const AA_NAMES: Record<string, string> = {
  'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
  'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
  'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
  'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
  '*': 'Stop',
};

/**
 * E. coli codon usage table (fraction per amino acid)
 * From Kazusa codon database
 */
const CODON_USAGE_ECOLI: Record<string, number> = {
  'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13,
  'CTT': 0.12, 'CTC': 0.10, 'CTA': 0.04, 'CTG': 0.47,
  'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1.00,
  'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.17, 'GTG': 0.35,
  'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14,
  'AGT': 0.16, 'AGC': 0.25, 'CCT': 0.18, 'CCC': 0.13,
  'CCA': 0.20, 'CCG': 0.49, 'ACT': 0.19, 'ACC': 0.40,
  'ACA': 0.17, 'ACG': 0.25, 'GCT': 0.18, 'GCC': 0.26,
  'GCA': 0.23, 'GCG': 0.33, 'TAT': 0.59, 'TAC': 0.41,
  'CAT': 0.57, 'CAC': 0.43, 'CAA': 0.34, 'CAG': 0.66,
  'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26,
  'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32,
  'TGT': 0.46, 'TGC': 0.54, 'TGG': 1.00, 'CGT': 0.36,
  'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11, 'AGA': 0.07,
  'AGG': 0.04, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13,
  'GGG': 0.15, 'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.30,
};

/**
 * Human codon usage table
 */
const CODON_USAGE_HUMAN: Record<string, number> = {
  'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13,
  'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.41,
  'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16, 'ATG': 1.00,
  'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47,
  'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06,
  'AGT': 0.15, 'AGC': 0.24, 'CCT': 0.28, 'CCC': 0.33,
  'CCA': 0.27, 'CCG': 0.11, 'ACT': 0.24, 'ACC': 0.36,
  'ACA': 0.28, 'ACG': 0.12, 'GCT': 0.26, 'GCC': 0.40,
  'GCA': 0.23, 'GCG': 0.11, 'TAT': 0.43, 'TAC': 0.57,
  'CAT': 0.41, 'CAC': 0.59, 'CAA': 0.25, 'CAG': 0.75,
  'AAT': 0.46, 'AAC': 0.54, 'AAA': 0.42, 'AAG': 0.58,
  'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58,
  'TGT': 0.45, 'TGC': 0.55, 'TGG': 1.00, 'CGT': 0.08,
  'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21, 'AGA': 0.20,
  'AGG': 0.20, 'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25,
  'GGG': 0.25, 'TAA': 0.28, 'TAG': 0.20, 'TGA': 0.52,
};

// =============================================================================
// Default Parameters
// =============================================================================

export const MUTAGENESIS_DEFAULTS: MutagenesisDefaults = {
  // Primer annealing parameters (NEB BaseChanger style)
  // The annealing portion is the 3' region that binds to template
  // NEB recommends minimum 15bp total primer length for specificity
  minAnnealingLength: 15,  // NEB minimum for specificity
  maxAnnealingLength: 35,  // Maximum annealing region

  // Legacy parameters for total primer length (including 5' tail)
  minPrimerLength: 15,
  maxPrimerLength: 60,

  // Flanking length for overlapping primer design
  minFlankingLength: 10,   // Minimum bp on each side of mutation
  maxFlankingLength: 25,   // Maximum bp on each side of mutation

  // Tm parameters (NEB BaseChanger style)
  // Primer is extended until minTm is reached
  // Higher minTm = longer primers = more specificity
  minTm: 55,   // NEB default: extend until at least 55°C
  maxTm: 72,   // Cap at 72°C

  // GC parameters
  minGC: 0.40,
  maxGC: 0.60,

  // Quality
  maxPolyN: 4,
  minDg: -5.0,
  gcClampRequired: true,

  // Off-target
  checkOffTargets: false, // Disabled by default for performance - runs on top candidates only
  maxOffTargets: 2,

  // Design strategy
  strategy: 'back-to-back', // 'back-to-back' (Q5 SDM) or 'overlapping' (QuikChange)

  // Codon optimization
  organism: 'ecoli', // 'ecoli', 'human', or null for no preference

  // Sequence topology
  circular: true, // Default to circular (most templates are plasmids)

  // 5' Tail confinement - keep mutations in non-annealing 5' overhang
  // When enabled, the 3' annealing region is 100% template-complementary
  // This improves PCR efficiency but may not work for all mutation types
  confineTo5Tails: false,
};

export const MUTATION_TYPES: MutationTypes = {
  SUBSTITUTION: 'substitution',
  INSERTION: 'insertion',
  DELETION: 'deletion',
  CODON_CHANGE: 'codon_change',
};

// Gas constant
const R = 1.987; // cal/(mol·K)

// =============================================================================
// Mismatch Tm Calculation (Novel Feature)
// =============================================================================

/**
 * Get the complementary base
 */
function getComplement(base: string): string {
  const comp: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
  return comp[base] || base;
}

/**
 * Check if two bases form a Watson-Crick pair
 */
function isWatsonCrickPair(base1: string, base2: string): boolean {
  return (base1 === 'A' && base2 === 'T') ||
         (base1 === 'T' && base2 === 'A') ||
         (base1 === 'G' && base2 === 'C') ||
         (base1 === 'C' && base2 === 'G');
}

/**
 * Get mismatch type key for looking up parameters
 * Format: 'primer_dinuc/template_strand_dinuc'
 */
function getMismatchKey(primerDinuc: string, templateStrandDinuc: string): string | null {
  // Try both orientations (tables may have either orientation)
  const key1 = `${primerDinuc}/${templateStrandDinuc}`;
  const key2 = `${templateStrandDinuc}/${primerDinuc}`;

  if (NN_MISMATCH[key1]) return key1;
  if (NN_MISMATCH[key2]) return key2;

  // Fallback: no matching key found
  return null;
}

/**
 * Calculate Tm for a primer with mismatches to template
 *
 * This is the enhanced NEB Base Changer algorithm - accounts for destabilization
 * caused by mismatched bases using Allawi/SantaLucia parameters.
 *
 * Enhanced features (vs original NEBaseChanger):
 * - Consecutive mismatch handling
 * - Dangling end effects
 * - Terminal vs internal mismatch differentiation
 * - Full Owczarzy Mg²⁺ salt correction
 *
 * @param primerSeq - Primer sequence (5'→3')
 * @param templateSeq - Coding/sense strand sequence (5'→3') - will be converted
 *                               to template strand internally for NN calculation
 * @param options - Calculation options
 * @returns Tm calculation result with details
 */
export function calculateMismatchedTm(primerSeq: string, templateSeq: string, options: {
  primerConc?: number;
  naConc?: number;
  mgConc?: number;
  useDanglingEnds?: boolean;
} = {}): MismatchedTmResult {
  const {
    primerConc = 500, // nM
    naConc = 50,      // mM Na+ (for salt correction)
    mgConc = 2.0,     // mM Mg2+
    useDanglingEnds = true,
  } = options;

  const primer = primerSeq.toUpperCase();
  const templateCoding = templateSeq.toUpperCase();

  if (primer.length !== templateCoding.length) {
    throw new Error('Primer and template must be same length for mismatch Tm calculation');
  }

  // Convert coding strand to template strand (complement) for proper NN calculation
  // The primer binds to the antisense/template strand, not the coding strand
  // This is consistent with how tm.js handles it - seq2 is the complement of seq1
  const template = templateCoding.split('').map(b => getComplement(b)).join('');

  // Initialize with duplex initiation (SantaLucia 1998)
  let dH = 0.2;   // kcal/mol
  let dS = -5.7;  // cal/(mol·K)

  const mismatches: Mismatch[] = [];
  let consecutiveMismatchCount = 0;
  let maxConsecutive = 0;
  let currentConsecutive = 0;

  // First pass: identify all mismatches and their positions
  for (let i = 0; i < primer.length; i++) {
    if (!isWatsonCrickPair(primer[i], template[i])) {
      mismatches.push({
        position: i,
        primerBase: primer[i],
        templateBase: template[i],
        type: `${primer[i]}·${template[i]}`,
        isTerminal: i === 0 || i === primer.length - 1,
        is3prime: i >= primer.length - 3,
      });
      currentConsecutive++;
      maxConsecutive = Math.max(maxConsecutive, currentConsecutive);
    } else {
      if (currentConsecutive > 1) {
        consecutiveMismatchCount += currentConsecutive - 1;
      }
      currentConsecutive = 0;
    }
  }
  if (currentConsecutive > 1) {
    consecutiveMismatchCount += currentConsecutive - 1;
  }

  // Calculate nearest-neighbor contributions
  for (let i = 0; i < primer.length - 1; i++) {
    const primerDinuc = primer[i] + primer[i + 1];
    const templateDinuc = template[i] + template[i + 1];

    // Check for mismatches in this dinucleotide
    const mismatch1 = !isWatsonCrickPair(primer[i], template[i]);
    const mismatch2 = !isWatsonCrickPair(primer[i + 1], template[i + 1]);
    const isTerminalPair = (i === 0 || i === primer.length - 2);

    if (!mismatch1 && !mismatch2) {
      // Perfect match - use standard NN params
      const params = NN_MATCHED[primerDinuc];
      if (params) {
        dH += params.dH;
        dS += params.dS;
      }
    } else {
      // Contains mismatch - look up mismatch parameters
      // First try our local mismatch table (format: 'primer_dinuc/template_dinuc')
      const key = getMismatchKey(primerDinuc, templateDinuc);

      if (key && NN_MISMATCH[key]) {
        const params = NN_MISMATCH[key];
        dH += params.dH;
        dS += params.dS;
      } else {
        // Try the comprehensive DNA_INTERNAL_MM or DNA_TERMINAL_MM from dna.js
        const nnKey = `${primerDinuc}/${templateDinuc}`;
        const nnKeyRev = `${templateDinuc}/${primerDinuc}`;

        if (isTerminalPair && DNA_TERMINAL_MM[nnKey]) {
          const [mmDh, mmDs] = DNA_TERMINAL_MM[nnKey];
          dH += mmDh;
          dS += mmDs;
        } else if (isTerminalPair && DNA_TERMINAL_MM[nnKeyRev]) {
          const [mmDh, mmDs] = DNA_TERMINAL_MM[nnKeyRev];
          dH += mmDh;
          dS += mmDs;
        } else if (DNA_INTERNAL_MM[nnKey]) {
          const [mmDh, mmDs] = DNA_INTERNAL_MM[nnKey];
          dH += mmDh;
          dS += mmDs;
        } else if (DNA_INTERNAL_MM[nnKeyRev]) {
          const [mmDh, mmDs] = DNA_INTERNAL_MM[nnKeyRev];
          dH += mmDh;
          dS += mmDs;
        } else {
          // Fallback: use average mismatch penalty
          dH += 1.0;  // Destabilizing (positive)
          dS += 2.5;
        }
      }
    }
  }

  // Terminal AT penalties
  if (primer[0] === 'A' || primer[0] === 'T') {
    dH += 2.2;
    dS += 6.9;
  }
  if (primer[primer.length - 1] === 'A' || primer[primer.length - 1] === 'T') {
    dH += 2.2;
    dS += 6.9;
  }

  // Terminal mismatch corrections (Bommarito 2000)
  const terminalMismatches = mismatches.filter(m => m.isTerminal);
  for (const mm of terminalMismatches) {
    if (mm.position === 0) {
      // 5' terminal mismatch - less destabilizing
      dH += TERMINAL_CORRECTIONS['5prime'].dH;
      dS += TERMINAL_CORRECTIONS['5prime'].dS;
    }
    if (mm.position === primer.length - 1) {
      // 3' terminal mismatch - more destabilizing (critical for extension!)
      dH += TERMINAL_CORRECTIONS['3prime'].dH;
      dS += TERMINAL_CORRECTIONS['3prime'].dS;
    }
  }

  // Dangling end effects (Bommarito 2000)
  // Apply if using back-to-back design where primers have 5' overhangs
  if (useDanglingEnds && mismatches.length > 0) {
    // 5' dangling end effect
    if (mismatches.some(m => m.position === 0)) {
      const base = primer[0];
      const deKey = `5_${base}`;
      if (DANGLING_END_CORRECTIONS[deKey]) {
        dH += DANGLING_END_CORRECTIONS[deKey].dH;
        dS += DANGLING_END_CORRECTIONS[deKey].dS;
      }
    }
    // 3' dangling end effect (more critical)
    if (mismatches.some(m => m.position === primer.length - 1)) {
      const base = primer[primer.length - 1];
      const deKey = `3_${base}`;
      if (DANGLING_END_CORRECTIONS[deKey]) {
        dH += DANGLING_END_CORRECTIONS[deKey].dH;
        dS += DANGLING_END_CORRECTIONS[deKey].dS;
      }
    }
  }

  // Consecutive mismatch penalty (Peyret 1999)
  // Tandem mismatches are more destabilizing than isolated ones
  if (consecutiveMismatchCount > 0) {
    const penaltyDh = Math.min(
      consecutiveMismatchCount * CONSECUTIVE_MISMATCH_CORRECTION.additionalPenalty.dH,
      CONSECUTIVE_MISMATCH_CORRECTION.maxPenalty.dH
    );
    const penaltyDs = Math.min(
      consecutiveMismatchCount * CONSECUTIVE_MISMATCH_CORRECTION.additionalPenalty.dS,
      CONSECUTIVE_MISMATCH_CORRECTION.maxPenalty.dS
    );
    dH += penaltyDh;
    dS += penaltyDs;
  }

  // Salt correction (Owczarzy 2008)
  // Use comprehensive Mg²⁺ correction when Mg is present
  let saltCorr: number;
  if (mgConc > 0) {
    // Owczarzy Mg²⁺ correction coefficients
    const a = 3.92e-5;
    const b = -9.11e-6;
    const c = 6.26e-5;
    const d = 1.42e-5;
    const e = -4.82e-4;
    const f = 5.25e-4;
    const g = 8.31e-5;

    const N = primer.length;
    const fGC = calculateGC(primer);
    const Mg = mgConc * 1e-3;
    const lnMg = Math.log(Mg);

    // This correction will be applied to the final Tm
    saltCorr = (a + b * lnMg + fGC * (c + d * lnMg)) +
               (1 / (2 * (N - 1))) * (e + f * lnMg + g * lnMg * lnMg);
  } else {
    // Simplified Na⁺ correction
    saltCorr = 0;
    dS += 0.368 * (primer.length - 1) * Math.log(naConc * 1e-3);
  }

  // Calculate Tm at 1M
  const Ct = primerConc * 1e-9;
  let Tm = (dH * 1000) / (dS + R * Math.log(Ct / 4)) - 273.15;

  // Apply Mg²⁺ correction if applicable
  if (mgConc > 0 && saltCorr !== 0) {
    const invTm = (1 / (Tm + 273.15)) + saltCorr;
    Tm = (1 / invTm) - 273.15;
  }

  // Identify critical mismatches
  const hasCritical3primeMismatch = mismatches.some(m => m.is3prime);
  const hasTerminalMismatch = mismatches.some(m => m.isTerminal);

  // Sanity check: if too many mismatches or Tm is unreasonable, primer won't bind
  // More than 50% mismatches means the primer essentially won't hybridize
  const mismatchFraction = mismatches.length / primer.length;
  const tooManyMismatches = mismatchFraction > 0.5;
  const unreasonableTm = Tm < -50 || !isFinite(Tm);
  const willNotBind = tooManyMismatches || unreasonableTm;

  return {
    tm: willNotBind ? null : Math.round(Tm * 10) / 10,
    willNotBind,
    mismatchFraction: Math.round(mismatchFraction * 100),
    dH,
    dS,
    mismatches,
    mismatchCount: mismatches.length,
    consecutiveMismatchCount,
    maxConsecutiveMismatches: maxConsecutive,
    hasCritical3primeMismatch,
    hasTerminalMismatch,
    has5primeMismatch: mismatches.some(m => m.position === 0),
    has3primeMismatch: mismatches.some(m => m.position === primer.length - 1),
    thermodynamics: {
      dH,
      dS,
      saltCorrection: saltCorr,
    },
  };
}

// =============================================================================
// Secondary Structure Analysis for Mutant Primers
// =============================================================================

/**
 * Check if mutation introduces problematic secondary structures
 * Based on AutoDimer algorithm (Vallone & Butler 2004)
 *
 * Enhanced features:
 * - Uses full Zuker algorithm from fold.js
 * - Detects 3' end involvement
 * - Compares with original primer (if provided)
 * - Heterodimer detection for primer pairs
 * - Evidence-based thresholds from IDT and Premier Biosoft
 */
export function checkMutantSecondaryStructure(mutantPrimer: string, originalPrimer: string | null = null, options: {
  temperature?: number;
  hairpinThreshold?: number;
  criticalHairpinThreshold?: number;
  selfDimerThreshold?: number;
  criticalSelfDimerThreshold?: number;
} = {}): SecondaryStructureCheck {
  // Use evidence-based thresholds from DIMER_THRESHOLDS (see references above)
  // Users can override if needed, but defaults are literature-supported
  const {
    temperature = 37,
    // Hairpin thresholds (internal - 3' involvement checked separately)
    hairpinThreshold = DIMER_THRESHOLDS.hairpin.internal.ideal,           // -3.0
    criticalHairpinThreshold = DIMER_THRESHOLDS.hairpin.internal.critical, // -6.0
    // Self-dimer thresholds (internal)
    selfDimerThreshold = DIMER_THRESHOLDS.selfDimer.internal.ideal,        // -6.0
    criticalSelfDimerThreshold = DIMER_THRESHOLDS.selfDimer.internal.critical, // -9.0
  } = options;

  const seq = mutantPrimer.toUpperCase();
  const warnings: SecondaryStructureWarning[] = [];

  // Use full Zuker algorithm from fold.js for ACCURATE hairpin analysis
  // This replaces the simplified estimateHairpinDG() for primary values
  let foldDg: number;
  let structures: FoldStructure[];
  try {
    structures = fold(seq, temperature) as any;
    foldDg = structures.reduce((sum, s) => sum + s.e, 0);
    foldDg = Math.round(foldDg * 100) / 100;
  } catch (e) {
    foldDg = dg(seq, temperature);
    structures = [];
  }

  // Calculate self-dimer ΔG using alignment-based algorithm
  const selfDimerDg = estimateSelfDimerDG(seq);

  // Use foldDg as the PRIMARY hairpin value (accurate Zuker algorithm)
  // This ensures consistency between Quick Status, Analysis Details, and Visualizations
  const hairpinDg = foldDg;

  // Check if 3' end is involved in hairpin structure
  const basePairs = structures.flatMap(s => s.ij || []);
  const hairpinInvolves3Prime = basePairs.some(([i, j]) =>
    i >= seq.length - 5 || j >= seq.length - 5
  );

  // Use appropriate threshold based on 3' involvement
  const effectiveHairpinThreshold = hairpinInvolves3Prime
    ? DIMER_THRESHOLDS.hairpin.threePrime.ideal      // -2.0 (stricter)
    : hairpinThreshold;                               // -3.0 (default)
  const effectiveCriticalHairpin = hairpinInvolves3Prime
    ? DIMER_THRESHOLDS.hairpin.threePrime.critical   // -4.0 (stricter)
    : criticalHairpinThreshold;                       // -6.0 (default)

  // Analyze hairpin structures using accurate foldDg and position-aware thresholds
  if (foldDg < effectiveHairpinThreshold) {
    warnings.push({
      type: 'hairpin',
      dG: foldDg,
      severity: foldDg < effectiveCriticalHairpin ? 'critical' : 'warning',
      involves3Prime: hairpinInvolves3Prime,
      message: hairpinInvolves3Prime
        ? `Stable hairpin at 3' end (ΔG = ${foldDg.toFixed(1)} kcal/mol) - may block extension`
        : `Stable hairpin detected (ΔG = ${foldDg.toFixed(1)} kcal/mol)`,
      details: structures.filter(s => s.desc && s.desc.includes('HAIRPIN')),
    });
  }

  // Analyze self-dimer structures
  // Note: Self-dimer 3' involvement requires alignment analysis (done in estimateSelfDimerDG)
  if (selfDimerDg < selfDimerThreshold) {
    warnings.push({
      type: 'self-dimer',
      dG: selfDimerDg,
      severity: selfDimerDg < criticalSelfDimerThreshold ? 'critical' : 'warning',
      message: `Stable self-dimer detected (ΔG = ${selfDimerDg.toFixed(1)} kcal/mol)`,
    });
  }

  // Use unified severity classification for 3' structure (basePairs already computed above)
  const threePrimeSeverity = classify3PrimeStructureSeverity({
    energy: foldDg,
    basePairs,
    seqLength: seq.length,
  });

  // Only add warning if severity warrants it
  if (threePrimeSeverity.shouldWarn) {
    warnings.push({
      type: '3prime-structure',
      dG: foldDg,
      severity: threePrimeSeverity.level === 'critical' ? 'critical' : 'warning',
      severityLevel: threePrimeSeverity.level,
      message: threePrimeSeverity.message || '',
      tooltip: threePrimeSeverity.tooltip,
      details: {
        structureInvolves3prime: threePrimeSeverity.level !== 'none',
        severityInfo: threePrimeSeverity,
      },
    });
  }

  // Compare with original primer if provided
  let changeAnalysis: { hairpinChange: number; selfDimerChange: number; foldChange: number; madeWorse: boolean } | null = null;
  if (originalPrimer) {
    const origSeq = originalPrimer.toUpperCase();
    const origFold = dg(origSeq, temperature);
    const origSelfDimer = estimateSelfDimerDG(origSeq);

    changeAnalysis = {
      hairpinChange: foldDg - origFold,
      selfDimerChange: selfDimerDg - origSelfDimer,
      foldChange: foldDg - origFold,
      madeWorse: (foldDg < origFold - 1.0) || (selfDimerDg < origSelfDimer - 1.0),
    };

    if (changeAnalysis.madeWorse) {
      warnings.push({
        type: 'stability-change',
        severity: 'warning',
        message: `Mutation significantly destabilizes secondary structure`,
        details: changeAnalysis,
      });
    }
  }

  return {
    // hairpinDG now uses accurate Zuker algorithm (same as foldDG)
    // This ensures consistency across all displays
    hairpinDG: hairpinDg,
    selfDimerDG: selfDimerDg,
    foldDG: foldDg,
    structures,
    warnings,
    isAcceptable: warnings.filter(w => w.severity === 'critical').length === 0,
    changeFromOriginal: changeAnalysis,
  };
}

/**
 * Check if 3' end is involved in any secondary structure
 */
function check3primeInvolvement(seq: string, structures: FoldStructure[]): boolean {
  const n = seq.length;
  const threePrimePositions = [n - 1, n - 2, n - 3];

  for (const struct of structures) {
    if (struct.ij && struct.ij.length > 0) {
      for (const [i, j] of struct.ij) {
        if (threePrimePositions.includes(i) || threePrimePositions.includes(j)) {
          return true;
        }
      }
    }
  }
  return false;
}

/**
 * Check for heterodimer formation between two primers
 * Important for primer pair design
 *
 * For overlapping (QuikChange-style) mutagenesis primers, the forward and reverse
 * primers are intentionally designed as reverse complements of each other. In this
 * case, heterodimer formation is EXPECTED and should not be flagged as a warning.
 *
 * ENHANCED: Now includes 3' involvement analysis for accurate severity classification.
 * 3' extensible dimers (both 3' ends involved with ≥4 consecutive bp) are CRITICAL
 * as they cause primer extension artifacts. This matches the CrossDimerZipper visualization.
 *
 * Evidence-based thresholds from:
 * - IDT OligoAnalyzer: 3' end ≥ -6 kcal/mol, internal ≥ -9 kcal/mol
 * - Premier Biosoft: 3' end ≥ -5 kcal/mol, internal ≥ -6 kcal/mol
 * - Literature: Hot start does NOT prevent 3' extensible dimers (PMID: 32799617)
 *
 * @param primer1 - Forward primer sequence
 * @param primer2 - Reverse primer sequence
 * @param options - Configuration options
 */
export function checkHeterodimer(primer1: string, primer2: string, options: {
  threshold?: number;
  criticalThreshold?: number;
  threePrimeThreshold?: number;
  designType?: 'overlapping' | 'back-to-back' | null;
} = {}): HeterodimerResult {
  // Use evidence-based thresholds from DIMER_THRESHOLDS
  const {
    threshold = DIMER_THRESHOLDS.heterodimer.internal.ideal,           // -6.0
    criticalThreshold = DIMER_THRESHOLDS.heterodimer.internal.critical, // -9.0
    threePrimeThreshold = DIMER_THRESHOLDS.heterodimer.threePrime.warning, // -6.0 (stricter)
    designType = null,
  } = options;

  const seq1 = primer1.toUpperCase();
  const seq2 = primer2.toUpperCase();
  const rc2 = reverseComplement(seq2);

  // For overlapping design, check if primers are exact complements (expected behavior)
  const isOverlappingDesign = designType === 'overlapping';
  const isExactComplement = seq1 === rc2;

  let minDG = 0;
  let bestAlignment: HeterodimerAlignment | null = null;
  let bestMaxConsecutive = 0;

  // Check all possible alignments with enhanced tracking
  for (let offset = -seq1.length + 4; offset <= seq2.length - 4; offset++) {
    let dG = 0;
    let consecutive = 0;
    let maxConsecutive = 0;
    let alignedBases: AlignedBase[] = [];

    for (let i = 0; i < seq1.length; i++) {
      const j = i + offset;
      if (j >= 0 && j < rc2.length) {
        const pos2InSeq2 = seq2.length - 1 - j;
        if (isWatsonCrickPair(seq1[i], seq2[pos2InSeq2])) {
          consecutive++;
          maxConsecutive = Math.max(maxConsecutive, consecutive);
          alignedBases.push({
            pos1: i,
            pos2: pos2InSeq2,
            base1: seq1[i],
            base2: seq2[pos2InSeq2]
          });
          if (consecutive >= 3) {
            const dinuc = seq1[i - 1] + seq1[i];
            if (NN_MATCHED[dinuc]) {
              dG += NN_MATCHED[dinuc].dH - 310.15 * (NN_MATCHED[dinuc].dS / 1000);
            }
          }
        } else {
          consecutive = 0;
        }
      }
    }

    // Weight by both energy and consecutive matches for best alignment selection
    const weightedScore = Math.abs(dG) + (maxConsecutive * 2);
    const currentBestScore = Math.abs(minDG) + (bestMaxConsecutive * 2);

    if (weightedScore > currentBestScore || (dG < minDG && weightedScore >= currentBestScore)) {
      minDG = dG;
      bestMaxConsecutive = maxConsecutive;
      bestAlignment = {
        offset,
        alignedBases,
        complementaryRegionLength: alignedBases.length,
        maxConsecutive,
        // Track 3' involvement (last 5 bases of each primer)
        involves3Prime1: alignedBases.some(p => p.pos1 >= seq1.length - 5),
        involves3Prime2: alignedBases.some(p => p.pos2 >= seq2.length - 5),
        involves5Prime1: alignedBases.some(p => p.pos1 < 5),
        involves5Prime2: alignedBases.some(p => p.pos2 < 5),
      };
    }
  }

  const warnings: SecondaryStructureWarning[] = [];

  // For overlapping design with exact complement primers, heterodimer is expected
  if (isOverlappingDesign && isExactComplement) {
    // No warnings - this is correct QuikChange primer design
    return {
      heterodimerDG: minDG,
      alignment: bestAlignment,
      warnings: [],
      isAcceptable: true,
      isExpectedOverlap: true,
      overlapLength: seq1.length,
      message: `Full primer overlap expected for QuikChange design (${seq1.length} bp)`,
      severity: 'safe',
      dimerType: 'expected_overlap',
    };
  }

  // For overlapping design with partial complement, validate the overlap region
  if (isOverlappingDesign) {
    // Calculate the overlap percentage
    const overlapPercent = bestAlignment ?
      (bestAlignment.complementaryRegionLength / Math.min(seq1.length, seq2.length)) * 100 : 0;

    if (overlapPercent < 90) {
      warnings.push({
        type: 'incomplete_overlap',
        severity: 'warning',
        message: `Overlapping primers have only ${overlapPercent.toFixed(0)}% complementarity (expected ~100%)`,
        details: bestAlignment,
      });
    }

    return {
      heterodimerDG: minDG,
      alignment: bestAlignment,
      warnings,
      isAcceptable: overlapPercent >= 90,
      isExpectedOverlap: overlapPercent >= 90,
      overlapLength: bestAlignment?.complementaryRegionLength || 0,
      overlapPercent,
      severity: 'safe',
      dimerType: 'ligation_junction',
    };
  }

  // ==========================================================================
  // UNIFIED SEVERITY CLASSIFICATION for non-overlapping designs
  // This matches the CrossDimerZipper visualization logic exactly
  // ==========================================================================
  let severity: 'safe' | 'warning' | 'critical' = 'safe';
  let dimerType = 'none';

  const involves3Prime1 = bestAlignment?.involves3Prime1 || false;
  const involves3Prime2 = bestAlignment?.involves3Prime2 || false;
  const maxConsecutive = bestAlignment?.maxConsecutive || 0;
  const only5PrimeOverlap = (bestAlignment?.involves5Prime1 || bestAlignment?.involves5Prime2) &&
                            !involves3Prime1 && !involves3Prime2;

  // 3' Extensible Dimer: MOST CRITICAL - both 3' ends involved with strong binding
  // This causes primer extension artifacts where primers extend off each other
  if (involves3Prime1 && involves3Prime2 && maxConsecutive >= 4) {
    severity = 'critical';
    dimerType = '3prime_extensible';
    warnings.push({
      type: 'heterodimer_3prime_extensible',
      dG: minDG,
      severity: 'critical',
      message: `CRITICAL: 3' Extensible Dimer - both primer 3' ends involved with ${maxConsecutive} consecutive bp`,
      tooltip: "Both primers' 3' ends are hybridized together. DNA polymerase can extend these, creating primer-dimer artifacts that compete with your target amplification.",
      details: bestAlignment,
    });
  }
  // One 3' end involved with moderate binding
  else if ((involves3Prime1 || involves3Prime2) && maxConsecutive >= 3) {
    severity = 'warning';
    dimerType = 'internal_with_3prime';
    warnings.push({
      type: 'heterodimer_3prime',
      dG: minDG,
      severity: 'warning',
      message: `Warning: Heterodimer involves ${involves3Prime1 ? 'forward' : 'reverse'} primer 3' end (${maxConsecutive} consecutive bp)`,
      tooltip: "One primer's 3' end is involved in the dimer. This may reduce PCR efficiency.",
      details: bestAlignment,
    });
  }
  // Internal binding only (5' ends or middle) - less concerning
  else if (maxConsecutive >= 4) {
    severity = 'warning';
    dimerType = 'internal';
    if (minDG < threshold) {
      warnings.push({
        type: 'heterodimer_internal',
        dG: minDG,
        severity: 'warning',
        message: `Internal heterodimer (ΔG = ${minDG.toFixed(1)} kcal/mol, ${maxConsecutive} consecutive bp)`,
        tooltip: "Primers bind internally, not at 3' ends. This is less severe but may reduce efficiency at high primer concentrations.",
        details: bestAlignment,
      });
    }
  }
  // Threshold-based fallback for weak but notable interactions
  else if (minDG < threshold) {
    severity = minDG < criticalThreshold ? 'critical' : 'warning';
    dimerType = 'minor';
    warnings.push({
      type: 'heterodimer',
      dG: minDG,
      severity: severity,
      message: `Primers form heterodimer (ΔG = ${minDG.toFixed(1)} kcal/mol)`,
      details: bestAlignment,
    });
  }

  // Determine acceptability based on severity
  const isAcceptable = severity !== 'critical';

  return {
    heterodimerDG: minDG,
    alignment: bestAlignment,
    warnings,
    isAcceptable,
    severity,
    dimerType,
    // Expose 3' involvement for displays that want to show this detail
    involves3Prime1,
    involves3Prime2,
    maxConsecutive,
  };
}

/**
 * Estimate hairpin ΔG using simplified algorithm
 */
function estimateHairpinDG(seq: string): number {
  let minDG = 0;

  // Check for potential hairpin loops (min 4nt loop)
  for (let loopStart = 4; loopStart < seq.length - 4; loopStart++) {
    for (let loopLen = 4; loopLen <= 8 && loopStart + loopLen < seq.length; loopLen++) {
      let stemLen = 0;
      let dG = 0;

      // Loop penalty
      dG += 4.0 + 1.75 * Math.log(loopLen);

      // Check stem formation
      for (let i = 0; i < Math.min(loopStart, seq.length - loopStart - loopLen); i++) {
        const base5 = seq[loopStart - 1 - i];
        const base3 = seq[loopStart + loopLen + i];

        if (isWatsonCrickPair(base5, base3)) {
          stemLen++;
          const dinuc = base5 + seq[loopStart - i] || base5;
          if (NN_MATCHED[dinuc]) {
            dG += NN_MATCHED[dinuc].dH - 310.15 * (NN_MATCHED[dinuc].dS / 1000);
          }
        } else {
          break;
        }
      }

      if (stemLen >= 3 && dG < minDG) {
        minDG = dG;
      }
    }
  }

  return minDG;
}

/**
 * Estimate self-dimer ΔG
 */
function estimateSelfDimerDG(seq: string): number {
  const rc = reverseComplement(seq);
  let minDG = 0;

  // Check all possible alignments
  for (let offset = -seq.length + 4; offset <= seq.length - 4; offset++) {
    let dG = 0;
    let matchLen = 0;

    for (let i = 0; i < seq.length; i++) {
      const j = i + offset;
      if (j >= 0 && j < seq.length) {
        if (seq[i] === rc[j] && isWatsonCrickPair(seq[i], seq[seq.length - 1 - j])) {
          matchLen++;
          if (matchLen >= 3) {
            const dinuc = seq[i - 1] + seq[i];
            if (NN_MATCHED[dinuc]) {
              dG += NN_MATCHED[dinuc].dH - 310.15 * (NN_MATCHED[dinuc].dS / 1000);
            }
          }
        } else {
          matchLen = 0;
        }
      }
    }

    if (dG < minDG) {
      minDG = dG;
    }
  }

  return minDG;
}

// =============================================================================
// Optimal Codon Selection
// =============================================================================

/**
 * Select optimal codon for amino acid change using minimum parsimony
 * with optional codon usage optimization
 *
 * @param originalCodon - Current codon
 * @param targetAA - Target amino acid (one-letter code)
 * @param options - Selection options
 * @returns Selected codon with details
 */
export function selectOptimalCodon(originalCodon: string, targetAA: string, options: {
  organism?: 'ecoli' | 'human' | null;
} = {}): CodonSelection {
  const { organism = 'ecoli' } = options;

  const targetCodons = CODON_TABLE[targetAA.toUpperCase()];
  if (!targetCodons) {
    throw new Error(`Invalid amino acid: ${targetAA}`);
  }

  const usageTable = organism === 'ecoli' ? CODON_USAGE_ECOLI :
                     organism === 'human' ? CODON_USAGE_HUMAN : null;

  let bestCodon: string | null = null;
  let minChanges = 4;
  let bestUsage = 0;

  const candidates: CodonCandidate[] = [];

  for (const codon of targetCodons) {
    // Count nucleotide changes (Hamming distance)
    let changes = 0;
    const positions: number[] = [];
    for (let i = 0; i < 3; i++) {
      if (codon[i] !== originalCodon[i]) {
        changes++;
        positions.push(i + 1);
      }
    }

    const usage = usageTable ? (usageTable[codon] || 0) : 0.5;

    candidates.push({
      codon,
      changes,
      positions,
      usage,
      score: changes * 10 - usage * 5, // Lower is better
    });

    // Selection criteria: minimum changes first, then usage
    if (changes < minChanges || (changes === minChanges && usage > bestUsage)) {
      minChanges = changes;
      bestCodon = codon;
      bestUsage = usage;
    }
  }

  candidates.sort((a, b) => a.score - b.score);

  return {
    selectedCodon: bestCodon!,
    nucleotideChanges: minChanges,
    codonUsage: bestUsage,
    allCandidates: candidates,
    changedPositions: candidates.find(c => c.codon === bestCodon)?.positions || [],
  };
}

// =============================================================================
// Off-Target Checking
// =============================================================================

/**
 * Check primer for off-target binding sites
 */
export function checkPrimerSpecificity(primer: string, fullTemplate: string, maxMismatches: number = 3): SpecificityResult {
  const seq = primer.toUpperCase();
  const template = fullTemplate.toUpperCase();
  const rc = reverseComplement(seq);

  const bindingSites: BindingSite[] = [];

  // Use existing offTargets function for efficient checking
  const otForward = offTargets(seq, template);
  const otReverse = offTargets(rc, template);

  // Count binding sites
  let forwardSites = 0;
  let reverseSites = 0;

  for (let i = 0; i < template.length; i++) {
    if (otForward[i] > 0) forwardSites++;
    if (otReverse[i] > 0) reverseSites++;
  }

  // Also do detailed sliding window for report
  for (let i = 0; i <= template.length - seq.length; i++) {
    const region = template.slice(i, i + seq.length);

    // Check forward
    let mismatches = 0;
    for (let j = 0; j < seq.length; j++) {
      if (seq[j] !== region[j]) mismatches++;
    }
    if (mismatches > 0 && mismatches <= maxMismatches) {
      bindingSites.push({
        position: i,
        strand: 'forward',
        mismatches,
        sequence: region,
      });
    }

    // Check reverse
    mismatches = 0;
    for (let j = 0; j < rc.length; j++) {
      if (rc[j] !== region[j]) mismatches++;
    }
    if (mismatches > 0 && mismatches <= maxMismatches) {
      bindingSites.push({
        position: i,
        strand: 'reverse',
        mismatches,
        sequence: region,
      });
    }
  }

  return {
    offTargetCount: bindingSites.length,
    forwardSites,
    reverseSites,
    bindingSites: bindingSites.slice(0, 10), // Limit report
    isSpecific: bindingSites.length <= 2,
  };
}

// =============================================================================
// Back-to-Back Primer Design (Q5 SDM Style)
// =============================================================================

/**
 * Design back-to-back (non-overlapping) primers for Q5 SDM
 *
 * Key principle: 5' ends of primers are adjacent (back-to-back)
 * - Forward primer contains mutation, 3' end binds downstream
 * - Reverse primer's 5' end is immediately upstream of forward's 5' end
 * - This allows exponential amplification
 *
 * @param template - Template sequence
 * @param mutPosition - Mutation position
 * @param mutatedSeq - Mutated template sequence
 * @param options - Design options
 */
function designBackToBackPrimers(
  template: string,
  mutPosition: number,
  mutatedSeq: string,
  mutLength: number,
  options: Partial<MutagenesisDefaults> & { exhaustiveSearch?: boolean },
  deletionLength: number = 0
): { candidates: CandidatePair[]; diagnostics: any } {
  const opts = { ...MUTAGENESIS_DEFAULTS, ...options };
  const seq = template.toUpperCase();
  const mutSeq = mutatedSeq.toUpperCase();

  const candidates: CandidatePair[] = [];
  const isDeletion = deletionLength > 0 && mutLength === 0;
  const isInsertion = mutLength > 0 && deletionLength === 0;

  // Diagnostic tracking for failure analysis
  const diagnostics: any = {
    fwdRegionTooShort: 0,
    revRegionTooShort: 0,
    noFwdCandidatesInTmWindow: 0,
    noRevCandidatesInTmWindow: 0,
    fwdTmTooLow: { count: 0, lowestTm: Infinity, sequence: null },
    fwdTmTooHigh: { count: 0, highestTm: -Infinity, sequence: null },
    revTmTooLow: { count: 0, lowestTm: Infinity, sequence: null },
    revTmTooHigh: { count: 0, highestTm: -Infinity, sequence: null },
    fwdGcExtreme: { count: 0, examples: [] },
    revGcExtreme: { count: 0, examples: [] },
    positionsExplored: 0,
    pairsEvaluated: 0,
  };

  // For back-to-back design:
  // Forward primer: contains mutation, extends 3' (downstream)
  // Reverse primer: 5' end adjacent to forward's 5' end, extends 3' (upstream on - strand)
  //
  // JOINT TM OPTIMIZATION ALGORITHM:
  // Instead of greedily stopping at the first valid length, we:
  // 1. Generate ALL valid candidates for forward primer (within Tm window)
  // 2. Generate ALL valid candidates for reverse primer (within Tm window)
  // 3. Find the best-matched pair that minimizes Tm difference

  const minTm = opts.minTm || 55;  // NEB default: 55°C minimum
  const maxTotalTm = opts.maxTm || 72;  // Cap total primer Tm
  const minAnnealing = opts.minAnnealingLength || 15;  // NEB minimum for specificity
  const maxAnnealing = opts.maxAnnealingLength || 35;

  // NEB BaseChanger approach for 5' flanking:
  // - Substitutions: mutation ~10nt from 5' end (can be reduced for GC-rich sequences)
  // - Insertions: mutation at 5' end (minimal 5' context)
  // - Deletions: forward primer spans junction, need binding on both sides
  //   but total length should be controlled by Tm (not forced to be too long)
  let max5primeFlank: number;
  let min5primeFlank: number;  // Minimum 5' context for binding stability
  if (isDeletion) {
    // For deletions, the forward primer spans the junction
    // 5' side binds before deletion, 3' side binds after deletion
    // We want adequate binding on both sides but Tm should control total length
    max5primeFlank = 15;  // Maximum 5' binding
    min5primeFlank = 5;   // Reduced from 10 - let Tm control total length
  } else if (isInsertion) {
    max5primeFlank = 5;   // Insertions: mutation near 5' end
    min5primeFlank = 0;
  } else {
    // Substitutions: start small and increase if needed (prefer shorter primers)
    max5primeFlank = 10;  // Maximum ~10nt from 5' end (NEB default)
    min5primeFlank = 3;   // Minimum 3nt for basic binding stability at 5' end
  }

  const minFwd5prime = Math.max(0, mutPosition - max5primeFlank);
  const maxFwd5prime = Math.max(0, mutPosition - min5primeFlank);

  // Iterate over 5' flanking positions
  for (let fwd5prime = maxFwd5prime; fwd5prime >= minFwd5prime; fwd5prime--) {
    diagnostics.positionsExplored++;

    // For deletions, ensure minimum binding on 5' side of junction
    if (isDeletion) {
      const binding5prime = mutPosition - fwd5prime;
      if (binding5prime < 5) continue;
    }

    // Calculate the annealing region (3' portion after mutation)
    const annealingStart = mutPosition + mutLength;
    const annealingRegion = mutSeq.slice(annealingStart);

    if (annealingRegion.length < minAnnealing) {
      diagnostics.fwdRegionTooShort++;
      continue;
    }

    // Reverse primer: 5' end is at fwd5prime - 1, extends upstream
    const rev5prime = fwd5prime - 1;
    if (rev5prime < 0) continue;

    const revRegion = seq.slice(0, rev5prime + 1);

    // Allow shorter reverse primers (min 10bp) for extreme Tm matching cases
    // GC-rich primers can maintain adequate Tm even at shorter lengths
    const adaptiveMinAnnealing = Math.max(10, Math.min(minAnnealing, revRegion.length));

    if (revRegion.length < 10) {
      diagnostics.revRegionTooShort++;
      continue;
    }

    // Reverse the region for candidate generation (3' to 5' direction)
    const revRegionFlipped = revRegion.split('').reverse().join('');

    // =================================================================
    // JOINT TM OPTIMIZATION: Generate ALL valid candidates for both primers
    // =================================================================

    // Get ALL valid reverse primer candidates within Tm window
    // Use adaptive minimum length for regions with limited space
    let revCandidates = getValidPrimerCandidates(revRegionFlipped, {
      minTm,
      maxTm: maxTotalTm,
      minLength: adaptiveMinAnnealing,
      maxLength: Math.min(maxAnnealing, revRegion.length),
    });

    // RESCUE MODE: If no candidates in normal range, try with relaxed soft cap
    // This handles GC-rich regions where even minimum length exceeds 72°C
    // Hard cap at 76°C prevents truly problematic primers
    let revUsedRescue = false;
    if (revCandidates.length === 0) {
      revCandidates = getValidPrimerCandidates(revRegionFlipped, {
        minTm,
        maxTm: maxTotalTm,
        absoluteMaxTm: 76,  // Hard cap - Q5 limit
        minLength: adaptiveMinAnnealing,
        maxLength: Math.min(maxAnnealing, revRegion.length),
        rescueMode: true,
      });
      if (revCandidates.length > 0) {
        revUsedRescue = true;
        diagnostics.usedRescueMode = true;
      }
    }

    if (revCandidates.length === 0) {
      // Track WHY no reverse candidates - check Tm of available lengths
      diagnostics.noRevCandidatesInTmWindow++;
      for (let len = minAnnealing; len <= Math.min(maxAnnealing, revRegion.length); len++) {
        const testSeq = revRegionFlipped.slice(0, len);
        const testTm = calculateTmQ5(testSeq);
        const testGc = calculateGC(testSeq);
        if (testTm < minTm) {
          diagnostics.revTmTooLow.count++;
          if (testTm < diagnostics.revTmTooLow.lowestTm) {
            diagnostics.revTmTooLow.lowestTm = testTm;
            diagnostics.revTmTooLow.sequence = testSeq;
            diagnostics.revTmTooLow.gc = testGc;
          }
        } else if (testTm > 76) {  // Check against hard cap
          diagnostics.revTmTooHigh.count++;
          if (testTm > diagnostics.revTmTooHigh.highestTm) {
            diagnostics.revTmTooHigh.highestTm = testTm;
            diagnostics.revTmTooHigh.sequence = testSeq;
            diagnostics.revTmTooHigh.gc = testGc;
          }
        }
        if (testGc < 0.3 || testGc > 0.7) {
          diagnostics.revGcExtreme.count++;
        }
      }
      continue;
    }

    // For forward primer, generate candidates by varying the annealing length
    const fwdCandidates: any[] = [];
    let fwdTmIssues = { tooLow: 0, tooHigh: 0 };

    // Calculate original template position for 3' annealing when using confineTo5Tails
    // For insertions: original position = mutPosition (insertion doesn't consume template bases)
    // For deletions: original position = mutPosition + deletionLength (skip deleted region)
    // For substitutions: original position = mutPosition + mutLength (same as mutSeq)
    const originalAnnealingStart = isDeletion
      ? mutPosition + deletionLength
      : isInsertion
        ? mutPosition
        : mutPosition + mutLength;

    for (let annealingLen = minAnnealing; annealingLen <= Math.min(maxAnnealing, annealingRegion.length); annealingLen++) {
      const fwd3prime = annealingStart + annealingLen;
      if (fwd3prime > mutSeq.length) continue;

      let fwdSeq: string;
      if (opts.confineTo5Tails) {
        // Confine mutations to 5' tail: 3' annealing region uses ORIGINAL template
        // Structure: [5' tail with mutation from mutSeq] + [3' annealing from original seq]
        const tailEnd = mutPosition + mutLength;  // End of mutation in mutSeq
        const tail5prime = mutSeq.slice(fwd5prime, tailEnd);  // 5' portion with mutation
        const anneal3prime = seq.slice(originalAnnealingStart, originalAnnealingStart + annealingLen);  // From original

        // Check if we have enough original template for annealing
        if (anneal3prime.length < annealingLen) continue;

        fwdSeq = tail5prime + anneal3prime;
      } else {
        // Standard design: entire primer from mutated sequence
        fwdSeq = mutSeq.slice(fwd5prime, fwd3prime);
      }

      // Skip if total primer is too short or too long
      if (fwdSeq.length < opts.minPrimerLength || fwdSeq.length > opts.maxPrimerLength) continue;

      const fwdTm = calculateTmQ5(fwdSeq);
      const fwdGc = calculateGC(fwdSeq);

      // Hard ceiling: physically problematic above 76°C (Q5 limit)
      if (fwdTm > 76) {
        fwdTmIssues.tooHigh++;
        continue;
      }

      // Include candidates up to hard cap - let penalty system decide
      // Candidates above normal max (72°C) are marked as rescue candidates
      if (fwdTm >= minTm) {
        fwdCandidates.push({
          length: fwdSeq.length,
          annealingLength: annealingLen,
          tm: fwdTm,
          gc: fwdGc,
          gcPercent: (fwdGc * 100).toFixed(1) + '%',
          sequence: fwdSeq,
          hasGCClamp: fwdSeq[fwdSeq.length - 1] === 'G' || fwdSeq[fwdSeq.length - 1] === 'C',
          start: fwd5prime,
          end: fwd3prime,
          isRescue: fwdTm > maxTotalTm,  // Mark high-Tm candidates
          confinedTo5Tails: opts.confineTo5Tails,  // Track if mutation confined to 5' tail
        });
      } else {
        // Track why Tm failed (too low)
        fwdTmIssues.tooLow++;
        if (fwdTm < diagnostics.fwdTmTooLow.lowestTm) {
          diagnostics.fwdTmTooLow.lowestTm = fwdTm;
          diagnostics.fwdTmTooLow.sequence = fwdSeq;
          diagnostics.fwdTmTooLow.gc = fwdGc;
        }
        if (fwdGc < 0.3 || fwdGc > 0.7) {
          diagnostics.fwdGcExtreme.count++;
        }
      }
    }

    // Track if any forward candidates are using rescue mode (Tm > 72°C)
    const fwdUsedRescue = fwdCandidates.some(c => c.isRescue);
    if (fwdUsedRescue) {
      diagnostics.usedRescueMode = true;
    }

    if (fwdCandidates.length === 0) {
      diagnostics.noFwdCandidatesInTmWindow++;
      diagnostics.fwdTmTooLow.count += fwdTmIssues.tooLow;
      diagnostics.fwdTmTooHigh.count += fwdTmIssues.tooHigh;
      continue;
    }

    // =================================================================
    // JOINT OPTIMIZATION: Find best Tm-matched pairs from Cartesian product
    // =================================================================

    for (const fwd of fwdCandidates) {
      for (const rev of revCandidates) {
        diagnostics.pairsEvaluated++;

        // Convert reverse candidate to actual reverse complement sequence
        const rev3prime = rev5prime - rev.length + 1;
        if (rev3prime < 0) continue;

        const revSeq = reverseComplement(seq.slice(rev3prime, rev5prime + 1));
        const revTm = calculateTmQ5(revSeq);
        const revGc = calculateGC(revSeq);
        const revHasClamp = revSeq[revSeq.length - 1] === 'G' || revSeq[revSeq.length - 1] === 'C';

        // =================================================================
        // SCORING: Quadratic/Exponential "Soft Walls" Penalty System
        // =================================================================
        // Philosophy: Small errors are cheap, large errors are prohibitively
        // expensive. This prevents good GC from "paying off" bad Tm mismatch.
        // =================================================================
        let penalty = 0;
        const tmDiff = Math.abs(fwd.tm - revTm);

        // --- 1. Tm Difference (Critical Constraint) ---
        // Quadratic penalty with 1°C "free zone"
        // 0-1°C → 0 pts, 2°C → 2 pts, 3°C → 8 pts, 5°C → 32 pts
        const effectiveTmDiff = Math.max(0, tmDiff - 1.0);
        penalty += 2.0 * Math.pow(effectiveTmDiff, 2);

        // --- 2. Minimum Tm Floor (Hard Wall) ---
        // Being under minTm is a critical failure - quadratic penalty
        // 1° under → 10 pts, 3° under → 90 pts, 5° under → 250 pts
        if (fwd.tm < minTm) {
          penalty += 10.0 * Math.pow(minTm - fwd.tm, 2);
        }
        if (revTm < minTm) {
          penalty += 10.0 * Math.pow(minTm - revTm, 2);
        }

        // --- 3. Maximum Tm Ceiling (Soft Wall) ---
        // High Tm is less dangerous than low Tm - keep linear but steep
        if (fwd.tm > maxTotalTm) penalty += (fwd.tm - maxTotalTm) * 8;
        if (revTm > maxTotalTm) penalty += (revTm - maxTotalTm) * 8;

        // --- 4. Length Penalty (Safe Zone: 15-28bp) ---
        // Short primers with adequate Tm are desirable (high specificity)
        // Only penalize primers exceeding optimal max length
        const optimalMaxLen = 28;  // SDM often needs longer primers for mutation context
        if (fwd.length > optimalMaxLen) {
          penalty += (fwd.length - optimalMaxLen) * 1.5;
        }
        if (rev.length > optimalMaxLen) {
          penalty += (rev.length - optimalMaxLen) * 1.5;
        }

        // --- 5. GC Content (Asymmetric Parabolic Bowl) ---
        // Centered on 50%, but high GC (>60%) is more problematic
        // (secondary structure, mispriming risk)
        const targetGC = 0.50;
        const fwdGcDiff = fwd.gc - targetGC;
        const revGcDiff = revGc - targetGC;

        // High GC: steeper penalty (150x), Low GC: milder (80x)
        // 10% deviation → 1.5 or 0.8 pts, 20% deviation → 6 or 3.2 pts
        penalty += (fwdGcDiff > 0 ? 150 : 80) * Math.pow(fwdGcDiff, 2);
        penalty += (revGcDiff > 0 ? 150 : 80) * Math.pow(revGcDiff, 2);

        // --- 6. GC Clamp (Now handled by score3primeTerminalBase in Top 10 phase) ---
        // Removed redundant check - terminal base scoring handles this better

        // Template region for mismatch Tm calculation
        let templateRegion: string | null = null;
        if (!isDeletion) {
          const templateStart = fwd5prime;
          const templateEnd = fwd.end - mutLength + (isInsertion ? 0 : mutLength);
          if (templateEnd <= seq.length && templateEnd - templateStart === fwd.length) {
            templateRegion = seq.slice(templateStart, templateEnd);
          }
        }

        // Track if rescue mode was used for either primer
        const fwdIsRescue = fwd.isRescue || false;
        const revIsRescue = rev.isRescue || false;
        const usedRescueMode = fwdIsRescue || revIsRescue;

        // Calculate dG for both primers (hairpin/secondary structure stability)
        const fwdDg = dg(fwd.sequence, 37);
        const revDg = dg(revSeq, 37);

        const candidateEntry: CandidatePair = {
          forward: {
            sequence: fwd.sequence,
            length: fwd.length,
            tm: fwd.tm,
            gc: fwd.gc,
            gcPercent: fwd.gcPercent,
            dg: fwdDg,
            hasGCClamp: fwd.hasGCClamp,
            start: fwd.start,
            end: fwd.end,
            _templateRegion: templateRegion,
            _isDeletion: isDeletion,
            isRescue: fwdIsRescue,
          } as any,
          reverse: {
            sequence: revSeq,
            length: revSeq.length,
            tm: revTm,
            gc: revGc,
            gcPercent: (revGc * 100).toFixed(1) + '%',
            dg: revDg,
            hasGCClamp: revHasClamp,
            start: rev3prime,
            end: rev5prime + 1,
            isRescue: revIsRescue,
          },
          design: 'back-to-back',
          tmDiff,
          penalty,
        };

        // Add warning if rescue mode was used (GC-rich region required elevated Tm)
        if (usedRescueMode) {
          candidateEntry.warnings = candidateEntry.warnings || [];
          candidateEntry.warnings.push({
            type: 'rescue_mode',
            severity: 'info',
            message: `GC-rich region: ${fwdIsRescue ? 'forward' : ''}${fwdIsRescue && revIsRescue ? ' and ' : ''}${revIsRescue ? 'reverse' : ''} primer Tm exceeds normal range (>72°C). Consider gradient PCR optimization.`,
          });
        }

        candidates.push(candidateEntry);
      }
    }
  }

  return { candidates, diagnostics };
}

// =============================================================================
// Overlapping Primer Design (QuikChange Style)
// =============================================================================

/**
 * Design overlapping primers (traditional QuikChange)
 */
function designOverlappingPrimers(
  template: string,
  mutPosition: number,
  mutatedSeq: string,
  mutLength: number,
  options: Partial<MutagenesisDefaults>
): CandidatePair[] {
  const opts = { ...MUTAGENESIS_DEFAULTS, ...options };
  const mutSeq = mutatedSeq.toUpperCase();

  const candidates: CandidatePair[] = [];

  for (let leftFlank = opts.minFlankingLength; leftFlank <= opts.maxFlankingLength; leftFlank++) {
    for (let rightFlank = opts.minFlankingLength; rightFlank <= opts.maxFlankingLength; rightFlank++) {
      const fwdStart = Math.max(0, mutPosition - leftFlank);
      const fwdEnd = Math.min(mutSeq.length, mutPosition + mutLength + rightFlank);

      const fwdSeq = mutSeq.slice(fwdStart, fwdEnd);
      const revSeq = reverseComplement(fwdSeq);

      if (fwdSeq.length < opts.minPrimerLength || fwdSeq.length > opts.maxPrimerLength) continue;

      const fwdTm = calculateTmQ5(fwdSeq);
      const revTm = calculateTmQ5(revSeq);
      const fwdGc = calculateGC(fwdSeq);
      const fwdDg = dg(fwdSeq, 37);

      let penalty = 0;

      // NEB-style: penalize below minTm, mild penalty above maxTm
      if (fwdTm < opts.minTm) penalty += (opts.minTm - fwdTm) * 10;
      if (fwdTm > opts.maxTm) penalty += (fwdTm - opts.maxTm) * 3;

      if (fwdGc < opts.minGC) penalty += (opts.minGC - fwdGc) * 20;
      if (fwdGc > opts.maxGC) penalty += (fwdGc - opts.maxGC) * 20;

      // Prefer shorter primers (more specific) when Tm is adequate
      const minLen = opts.minAnnealingLength || 15;
      penalty += Math.max(0, fwdSeq.length - minLen) * 0.2;

      if (fwdDg < opts.minDg) penalty += Math.abs(fwdDg - opts.minDg) * 2;

      const fwdHasClamp = fwdSeq[fwdSeq.length - 1] === 'G' || fwdSeq[fwdSeq.length - 1] === 'C';
      const revHasClamp = revSeq[revSeq.length - 1] === 'G' || revSeq[revSeq.length - 1] === 'C';
      if (opts.gcClampRequired && !fwdHasClamp) penalty += 5;

      candidates.push({
        forward: {
          sequence: fwdSeq,
          length: fwdSeq.length,
          tm: fwdTm,
          gc: fwdGc,
          gcPercent: (fwdGc * 100).toFixed(1) + '%',
          dg: fwdDg,
          hasGCClamp: fwdHasClamp,
          start: fwdStart,
          end: fwdEnd,
        },
        reverse: {
          sequence: revSeq,
          length: revSeq.length,
          tm: revTm,
          gc: calculateGC(revSeq),
          gcPercent: (calculateGC(revSeq) * 100).toFixed(1) + '%',
          dg: dg(revSeq, 37),
          hasGCClamp: revHasClamp,
          // Overlapping design: reverse primer has same position as forward (they're complementary)
          start: fwdStart,
          end: fwdEnd,
        },
        design: 'overlapping',
        leftFlank,
        rightFlank,
        penalty,
      });
    }
  }

  return candidates;
}

// =============================================================================
// Main Design Functions
// =============================================================================

// [Due to length, I'll continue this in the next message. For now, I'll add the helper functions and then continue with the main design function]

/**
 * Assign candidates to quality tiers based on hard thresholds.
 * This is more interpretable than pure Pareto sorting.
 *
 * Tier 1 (Excellent): Tm diff ≤2°C AND structure ΔG > -3
 * Tier 2 (Good): Tm diff ≤5°C AND structure ΔG > -5
 * Tier 3 (Acceptable): Tm diff ≤8°C
 * Tier 4 (Poor): Everything else
 */
function assignTiers(candidates: CandidatePair[], options: { minScoreForTier1?: number } = {}): {
  tier1: CandidatePair[];
  tier2: CandidatePair[];
  tier3: CandidatePair[];
  tier4: CandidatePair[];
} {
  const {
    // Minimum composite score for tier 1 (0-100 scale)
    // Score >= 70 is required for tier 1 (consistent with primers.js)
    minScoreForTier1 = 70,
  } = options;

  const tiers = {
    tier1: [] as CandidatePair[], // Excellent
    tier2: [] as CandidatePair[], // Good
    tier3: [] as CandidatePair[], // Acceptable
    tier4: [] as CandidatePair[], // Poor
  };

  for (const c of candidates) {
    const tmDiff = Math.abs(c.forward.tm - c.reverse.tm);
    const fwdDg = c.forward.dg || 0;
    const revDg = c.reverse.dg || 0;
    const worstDg = Math.min(fwdDg, revDg);

    // Check composite score threshold - low score candidates should not be tier 1
    // compositeScore is 0-100 (higher = better), so we check >= threshold
    const hasCompositeScore = typeof c.compositeScore === 'number';
    const scoreTooLow = hasCompositeScore && (c.compositeScore ?? 0) < minScoreForTier1;

    if (tmDiff <= 2 && worstDg > -3 && !scoreTooLow) {
      tiers.tier1.push(c);
    } else if (tmDiff <= 5 && worstDg > -5) {
      // Candidates with low score but good Tm/dG go to tier 2
      tiers.tier2.push(c);
    } else if (tmDiff <= 8) {
      tiers.tier3.push(c);
    } else {
      tiers.tier4.push(c);
    }
  }

  // Sort each tier by composite score (descending - higher is better)
  // Fall back to penalty (ascending) if no composite score
  for (const tier of Object.values(tiers)) {
    tier.sort((a, b) => {
      if (typeof a.compositeScore === 'number' && typeof b.compositeScore === 'number') {
        return b.compositeScore - a.compositeScore; // Higher score first
      }
      return a.penalty - b.penalty; // Lower penalty first (fallback)
    });
  }

  return tiers;
}

/**
 * Select best candidate using tiered selection.
 * Picks from the best non-empty tier, sorted by penalty within tier.
 */
function selectBestByTier(candidates: CandidatePair[]): (CandidatePair & { tierQuality?: string }) | null {
  const tiers = assignTiers(candidates);

  if (tiers.tier1.length > 0) return { ...tiers.tier1[0], tierQuality: 'excellent' };
  if (tiers.tier2.length > 0) return { ...tiers.tier2[0], tierQuality: 'good' };
  if (tiers.tier3.length > 0) return { ...tiers.tier3[0], tierQuality: 'acceptable' };
  if (tiers.tier4.length > 0) return { ...tiers.tier4[0], tierQuality: 'poor' };

  return null;
}

/**
 * Check for problematic sequence patterns at a position
 * Used by sliding window to avoid bad junction contexts
 */
function checkSequenceContext(seq: string, pos: number): { hasProblems: boolean; issues: any[] } {
  const issues: any[] = [];
  const windowSize = 6;
  const start = Math.max(0, pos - windowSize);
  const end = Math.min(seq.length, pos + windowSize);
  const context = seq.slice(start, end).toUpperCase();

  // Check for homopolymer runs (GGGG, CCCC, AAAA, TTTT)
  if (/(.)\1{3,}/.test(context)) {
    const match = context.match(/(.)\1{3,}/);
    if (match) issues.push({ type: 'homopolymer', pattern: match[0], severity: 'warning' });
  }

  // Check for dinucleotide repeats (ATAT, GCGC, etc.)
  if (/(..)\1{2,}/.test(context)) {
    const match = context.match(/(..)\1{2,}/);
    if (match) issues.push({ type: 'dinuc_repeat', pattern: match[0], severity: 'warning' });
  }

  // Check for high GC content in local context (potential hairpin)
  const localGC = (context.match(/[GC]/g) || []).length / context.length;
  if (localGC > 0.8) {
    issues.push({ type: 'high_gc', gc: localGC, severity: 'caution' });
  }

  // Check for palindromic sequences (potential hairpin)
  const rc = reverseComplement(context);
  if (context.includes(rc.slice(2, 6)) || rc.includes(context.slice(2, 6))) {
    issues.push({ type: 'palindrome', severity: 'warning' });
  }

  return {
    hasProblems: issues.some(i => i.severity === 'warning'),
    issues,
  };
}

// =============================================================================
// Public API Functions for Mutagenesis Primer Design
// =============================================================================

/**
 * Parse mutation notation (e.g., "A123G", "K45R", "del100-110")
 */
export function parseMutationNotation(notation: string): any {
  const upper = notation.trim().toUpperCase();

  // Point mutation: A123G (nucleotide) or K45R (amino acid)
  const pointMatch = upper.match(/^([A-Z])(\d+)([A-Z])$/);
  if (pointMatch) {
    const [, original, posStr, replacement] = pointMatch;
    const position = parseInt(posStr, 10);
    const isAA = 'ACDEFGHIKLMNPQRSTVWY'.includes(original);
    return {
      type: isAA ? 'substitution' : 'point',
      position: position - 1,  // Convert to 0-indexed
      original,
      replacement,
      notation: upper
    };
  }

  // Deletion: del100-110 or Δ100-110
  const delMatch = upper.match(/^(?:DEL|Δ)(\d+)(?:-(\d+))?$/);
  if (delMatch) {
    const start = parseInt(delMatch[1], 10) - 1;
    const end = delMatch[2] ? parseInt(delMatch[2], 10) - 1 : start;
    return {
      type: 'deletion',
      position: start,
      length: end - start + 1,
      notation: upper
    };
  }

  // Insertion: ins100_ACGT
  const insMatch = upper.match(/^INS(\d+)_([ACGT]+)$/);
  if (insMatch) {
    return {
      type: 'insertion',
      position: parseInt(insMatch[1], 10) - 1,
      insertion: insMatch[2],
      notation: upper
    };
  }

  return null;
}

/**
 * Design primers for deleting a region from template
 * @param template - Template sequence
 * @param start - Start position of deletion (0-indexed)
 * @param length - Number of bases to delete
 * @param options - Design options
 */
export function designDeletionPrimers(
  template: string,
  start: number,
  length: number,
  options: Partial<MutagenesisDefaults> = {}
): any {
  const seq = template.toUpperCase();
  const opts = { ...MUTAGENESIS_DEFAULTS, ...options };

  // Create mutated sequence by removing the deletion region
  const mutatedSeq = seq.slice(0, start) + seq.slice(start + length);

  // Design primers using back-to-back approach
  const result = designBackToBackPrimers(seq, start, mutatedSeq, 0, opts, length);

  const best = selectBestByTier(result.candidates);

  if (!best) {
    return {
      forward: { sequence: '', tm: 0, gc: 0 },
      reverse: { sequence: '', tm: 0, gc: 0 },
      error: 'No suitable primers found for deletion'
    };
  }

  return {
    forward: best.forward,
    reverse: best.reverse,
    mutatedSequence: mutatedSeq,
    quality: best.tierQuality || 'unknown',
    qualityTier: best.tierQuality || 'unknown',
    compositeScore: best.compositeScore ?? 75,
    effectiveScore: best.compositeScore ?? 75,
    alternateDesigns: result.candidates.filter(c => c !== best).slice(0, 10),
  };
}

/**
 * Design primers for inserting a sequence into template
 * @param template - Template sequence
 * @param position - Position to insert at (0-indexed)
 * @param insertion - Sequence to insert
 * @param options - Design options
 */
export function designInsertionPrimers(
  template: string,
  position: number,
  insertion: string,
  options: Partial<MutagenesisDefaults> = {}
): any {
  const seq = template.toUpperCase();
  const insertSeq = insertion.toUpperCase();
  const opts = { ...MUTAGENESIS_DEFAULTS, ...options };

  // Create mutated sequence by inserting the new sequence
  const mutatedSeq = seq.slice(0, position) + insertSeq + seq.slice(position);

  // Design primers
  const result = designBackToBackPrimers(seq, position, mutatedSeq, insertSeq.length, opts, 0);

  const best = selectBestByTier(result.candidates);

  if (!best) {
    return {
      forward: { sequence: '', tm: 0, gc: 0 },
      reverse: { sequence: '', tm: 0, gc: 0 },
      protocol: { name: 'Insertion', steps: [], notes: [] },
      error: 'No suitable primers found for insertion'
    };
  }

  return {
    forward: best.forward,
    reverse: best.reverse,
    protocol: generateProtocol('insertion', best, opts),
    mutatedSequence: mutatedSeq,
    quality: best.tierQuality || 'unknown',
    qualityTier: best.tierQuality || 'unknown',
    compositeScore: best.compositeScore ?? 75,
    effectiveScore: best.compositeScore ?? 75,
    alternateDesigns: result.candidates.filter(c => c !== best).slice(0, 10),
  };
}

/**
 * Design primers for substituting a region in template
 * @param template - Template sequence
 * @param start - Start position of substitution (0-indexed)
 * @param length - Length of region to replace
 * @param replacement - Replacement sequence
 * @param options - Design options
 */
export function designRegionSubstitutionPrimers(
  template: string,
  start: number,
  length: number,
  replacement: string,
  options: Partial<MutagenesisDefaults> = {}
): any {
  const seq = template.toUpperCase();
  const replaceSeq = replacement.toUpperCase();
  const opts = { ...MUTAGENESIS_DEFAULTS, ...options };

  // Create mutated sequence
  const mutatedSeq = seq.slice(0, start) + replaceSeq + seq.slice(start + length);

  // Determine if this is more like a deletion or insertion
  const netChange = replaceSeq.length - length;

  const result = designBackToBackPrimers(
    seq,
    start,
    mutatedSeq,
    replaceSeq.length,
    opts,
    length
  );

  const best = selectBestByTier(result.candidates);

  if (!best) {
    return {
      forward: { sequence: '', tm: 0, gc: 0 },
      reverse: { sequence: '', tm: 0, gc: 0 },
      protocol: { name: 'Substitution', steps: [], notes: [] },
      error: 'No suitable primers found for substitution'
    };
  }

  return {
    forward: best.forward,
    reverse: best.reverse,
    protocol: generateProtocol('substitution', best, opts),
    mutatedSequence: mutatedSeq,
    quality: best.tierQuality || 'unknown',
    qualityTier: best.tierQuality || 'unknown',
    compositeScore: best.compositeScore ?? 75,
    effectiveScore: best.compositeScore ?? 75,
    alternateDesigns: result.candidates.filter(c => c !== best).slice(0, 10),
  };
}

/**
 * Design primers for changing a codon to encode a different amino acid
 * @param template - Template sequence
 * @param codonPosition - Position of the first base of the codon (0-indexed)
 * @param targetAA - Target amino acid (single letter code)
 * @param options - Design options
 */
export function designCodonChangePrimers(
  template: string,
  codonPosition: number,
  targetAA: string,
  options: Partial<MutagenesisDefaults> = {}
): any {
  const seq = template.toUpperCase();
  const opts = { ...MUTAGENESIS_DEFAULTS, ...options };

  // Get original codon
  const originalCodon = seq.slice(codonPosition, codonPosition + 3);

  // Select optimal codon for target amino acid
  const codonSelection = selectOptimalCodon(originalCodon, targetAA.toUpperCase(), opts);

  if (!codonSelection.selectedCodon) {
    return {
      forward: { sequence: '', tm: 0, gc: 0 },
      reverse: { sequence: '', tm: 0, gc: 0 },
      protocol: { name: 'Codon Change', steps: [], notes: [] },
      error: `No valid codon found for amino acid ${targetAA}`
    };
  }

  // Create mutated sequence with new codon
  const mutatedSeq = seq.slice(0, codonPosition) + codonSelection.selectedCodon + seq.slice(codonPosition + 3);

  // Design primers
  const result = designBackToBackPrimers(seq, codonPosition, mutatedSeq, 3, opts, 3);

  const best = selectBestByTier(result.candidates);

  if (!best) {
    return {
      forward: { sequence: '', tm: 0, gc: 0 },
      reverse: { sequence: '', tm: 0, gc: 0 },
      protocol: { name: 'Codon Change', steps: [], notes: [] },
      error: 'No suitable primers found for codon change'
    };
  }

  return {
    forward: best.forward,
    reverse: best.reverse,
    protocol: generateProtocol('codon_change', best, opts),
    mutatedSequence: mutatedSeq,
    codonChange: {
      original: originalCodon,
      new: codonSelection.selectedCodon,
      targetAA,
      changes: codonSelection.changedPositions
    },
    quality: best.tierQuality || 'unknown',
    qualityTier: best.tierQuality || 'unknown',
    compositeScore: best.compositeScore ?? 75,
    effectiveScore: best.compositeScore ?? 75,
    alternateDesigns: result.candidates.filter(c => c !== best).slice(0, 10),
  };
}

/**
 * Generate protocol steps for mutagenesis
 */
function generateProtocol(
  mutationType: string,
  candidate: CandidatePair,
  _options: MutagenesisDefaults
): any {
  const steps = [
    {
      name: 'Prepare PCR master mix',
      temp: 'RT',
      time: '10 min'
    },
    {
      name: 'Set up PCR reaction',
      temp: 'RT',
      time: '5 min'
    },
    {
      name: 'Run PCR cycling',
      temp: `${Math.min(candidate.forward.tm, candidate.reverse.tm).toFixed(1)}°C annealing`,
      time: '2-3 hours'
    },
    {
      name: 'DpnI digest template',
      temp: '37°C',
      time: '1 hour'
    },
    {
      name: 'Transform competent cells',
      temp: '37°C',
      time: 'overnight'
    }
  ];

  const notes = [
    `Mutation type: ${mutationType}`,
    `Tm difference: ${Math.abs(candidate.forward.tm - candidate.reverse.tm).toFixed(1)}°C`,
    `Forward primer: ${candidate.forward.sequence}`,
    `Reverse primer: ${candidate.reverse.sequence}`
  ];

  return { name: `${mutationType} Mutagenesis Protocol`, steps, notes };
}

// [Continue with the remaining functions - I'll split this into multiple messages due to length constraints]
