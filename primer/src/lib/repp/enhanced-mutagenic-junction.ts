/**
 * Enhanced Mutagenic Junction Domestication for Golden Gate Assembly
 *
 * State-of-the-art domestication using NEB experimental fidelity data.
 *
 * Key improvements over basic implementation:
 * 1. Uses NEB's 256×256 ligation frequency matrix (Pryor et al. 2020)
 * 2. Wider search range (±50bp) for optimal junction positions
 * 3. Considers global overhang set compatibility
 * 4. Comprehensive primer quality scoring using thermodynamic analysis
 * 5. Reading frame-aware junction placement
 * 6. Integration-ready with Unified Primer Designer
 *
 * References:
 * - Pryor et al. (2020) PLOS ONE - Data-optimized assembly design
 * - Potapov et al. (2018) ACS Synth Biol - Overhang fidelity profiling
 */

import {
  GOLDEN_GATE_ENZYMES,
  findInternalSites,
  calculateExperimentalFidelity,
  getEnzymeLigationData,
  getOverhangFidelityExperimental,
} from './goldengate.js';
import { reverseComplement } from './enzymes.js';
import {
  CODON_TO_AA,
  CODON_TABLE,
  ECOLI_CODON_USAGE,
  YEAST_CODON_USAGE,
} from './silent-mutation-domesticator.js';
import { calculateHairpinDG, calculateHomodimerDG, calculateHeterodimerDG } from '../equilibrium.js';
import { calculateTmQ5, calculate3primeTerminalDG } from '../tmQ5.js';
import {
  scoreTm,
  scoreGc,
  scoreHairpin,
  scoreHomodimer,
  scoreLength,
  scoreGcClamp,
  scoreHomopolymer,
  scoreTerminal3DG,
  scoreGQuadruplex,
} from '../scoring.js';
import { OPTIMAL_FLANKING_SEQUENCES, OPTIMAL_SPACERS } from './goldengate-primer-optimizer.js';

// ============================================================================
// ENHANCED CONFIGURATION
// ============================================================================

export const ENHANCED_JUNCTION_CONFIG = {
  // Junction search parameters (ENHANCED)
  searchRadius: 50,
  extendedSearchRadius: 80,

  // Primer design parameters
  homologyLength: {
    optimal: 20,
    min: 15,
    max: 30,
  },
  targetTm: 62,
  tmTolerance: 5,
  useAdaptiveHomology: true,

  // Overhang selection (ENHANCED with NEB data)
  preferHighFidelityOverhangs: true,
  minOverhangFidelity: 0.95,
  targetAssemblyFidelity: 0.98,

  // Mutation preferences
  preferWobblePosition: true,
  avoidRareCodons: true,
  rareCodoThreshold: 10.0,

  // Scoring weights (CALIBRATED)
  weights: {
    overhangFidelity: 0.35,
    mutationQuality: 0.25,
    primerQuality: 0.25,
    positionOptimality: 0.15,
  },

  // Quality thresholds
  qualityThresholds: {
    excellent: 90,
    good: 75,
    acceptable: 60,
    poor: 0,
  },
} as const;

// ============================================================================
// TYPES
// ============================================================================

interface InternalSite {
  position: number;
  sequence: string;
  orientation: 'forward' | 'reverse';
}

interface CodonUsageTable {
  [codon: string]: number;
}

interface PrimerCoverage {
  fragment1Primer: Array<{ position: number; indexInSite: number }>;
  fragment2Primer: Array<{ position: number; indexInSite: number }>;
  totalCoverage: number;
  fragment1HomologyRange: { start: number; end: number };
  fragment2HomologyRange: { start: number; end: number };
}

interface MutationScoreBreakdown {
  wobbleBonus?: number;
  middlePositionPenalty?: number;
  rareCodonPenalty?: number;
  frequencyPreservedBonus?: number;
  positionBonus?: number;
}

interface Mutation {
  sequencePosition: number;
  positionInSite: number;
  positionInCodon: number;
  originalBase: string;
  newBase: string;
  originalCodon: string;
  newCodon: string;
  aminoAcid: string;
  codonFrequency: number;
  originalCodonFrequency: number;
  inFragment1Primer: boolean;
  inFragment2Primer: boolean;
  score: number;
  scoreBreakdown: MutationScoreBreakdown;
  isSynonymous: boolean;
  breaksSite: boolean;
}

interface JunctionCandidate {
  junctionPosition: number;
  overhang: string;
  coverage: PrimerCoverage;
  mutations: Mutation[];
  bestMutation: Mutation;
  positionOptimality: number;
  overhangFidelity: number | null;
  primerQuality: PrimerQualityScore | null;
  overallScore: number | null;
  scoring?: {
    overhangFidelityScore: number;
    mutationScore: number;
    primerScore: number;
    positionScore: number;
  };
  assemblyFidelityWithExisting?: number;
}

interface HomologyAnalysis {
  score: number;
  tm: number;
  gc: number;
  length: number;
  gcClamp: boolean;
  issues: string[];
}

interface PrimerQualityScore {
  score: number;
  fragment1: HomologyAnalysis;
  fragment2: HomologyAnalysis;
  tmDifference: number;
}

interface OverhangFidelityScore {
  score: number;
  singleFidelity: number;
  assemblyFidelity: number;
}

interface PrimerComponent {
  flanking: string;
  recognition: string;
  spacer: string;
  overhang: string;
  homology: string;
}

interface PrimerInfo {
  fullSequence?: string;
  sequence: string;
  components: PrimerComponent;
  hasMutation: boolean;
  mutationPosition: number | null;
  tm?: number;
  length: number;
  homologyLength: number;
}

interface PrimerDesign {
  fragment1: {
    reversePrimer: PrimerInfo;
    note: string;
  };
  fragment2: {
    forwardPrimer: PrimerInfo;
    note: string;
  };
  mutation: {
    applied: boolean;
    inFragment1Primer: boolean;
    inFragment2Primer: boolean;
    change: string;
    codonChange: string;
  };
  tmDifference: number;
  optimizedFlanking: string;
}

interface ValidationWarning {
  code: string;
  message: string;
  details?: any;
}

interface JunctionValidation {
  isValid: boolean;
  warnings: ValidationWarning[];
  errors: ValidationWarning[];
  overallQuality: string;
}

interface QualityBreakdown {
  overhangFidelity: number;
  mutationQuality: number;
  primerQuality: number;
  positionOptimality: number;
}

interface FidelityInfo {
  singleOverhang: number;
  withExisting: number;
  source: 'NEB_experimental' | 'calculated';
}

interface AlternativeJunction {
  junctionPosition: number;
  overhang: string;
  overallScore: number;
  overhangFidelity: number;
}

interface EnhancedJunctionResult {
  success: boolean;
  site: InternalSite;
  junctionPosition?: number;
  overhang?: string;
  quality?: {
    overall: number;
    tier: string;
    breakdown: QualityBreakdown;
  };
  fidelity?: FidelityInfo;
  mutation?: Mutation;
  alternativeMutations?: Mutation[];
  primers?: PrimerDesign;
  coverage?: PrimerCoverage;
  validation?: JunctionValidation;
  alternatives?: AlternativeJunction[];
  message?: string;
  error?: string;
  searchedRange?: {
    start: number;
    end: number;
  };
}

interface EnhancedJunctionOptions {
  frame?: number;
  organism?: string;
  existingOverhangs?: string[];
  searchRadius?: number;
  targetTm?: number;
  requireHighFidelity?: boolean;
}

interface BatchJunctionResult {
  success: boolean;
  needsDomestication: boolean;
  totalSites: number;
  sites: InternalSite[];
  junctions: EnhancedJunctionResult[];
  failedSites: Array<{
    site: InternalSite;
    error: string;
    message: string;
  }>;
  allOverhangs: string[];
  fidelity: {
    assembly: number;
    source: string;
  };
  additionalFragments: number;
  optimized: boolean;
  message: string;
}

type QualityTier = 'excellent' | 'good' | 'acceptable' | 'poor';

// ============================================================================
// MAIN ENHANCED FUNCTION
// ============================================================================

/**
 * Design mutagenic junction with state-of-the-art optimization
 */
export function designEnhancedMutagenicJunction(
  sequence: string,
  site: InternalSite,
  enzyme: string = 'BsaI',
  options: EnhancedJunctionOptions = {}
): EnhancedJunctionResult {
  const {
    frame = 0,
    organism = 'ecoli',
    existingOverhangs = [],
    searchRadius = ENHANCED_JUNCTION_CONFIG.searchRadius,
    targetTm = ENHANCED_JUNCTION_CONFIG.targetTm,
    requireHighFidelity = true,
  } = options;

  const seq = sequence.toUpperCase();
  const enz = GOLDEN_GATE_ENZYMES[enzyme];

  if (!enz) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  const recognition = enz.recognition;
  const overhangLen = enz.overhangLength || 4;
  const codonUsage = organism === 'yeast' ? YEAST_CODON_USAGE : ECOLI_CODON_USAGE;

  const siteStart = site.position;
  const siteEnd = siteStart + recognition.length;

  // Get NEB ligation data for accurate fidelity calculation
  const enzymeData = getEnzymeLigationData(enzyme);
  const hasExperimentalData = !!enzymeData;

  // Phase 1: Enumerate all candidate junction positions
  const candidates = enumerateJunctionCandidates(
    seq,
    site,
    enzyme,
    {
      searchRadius,
      frame,
      codonUsage,
      existingOverhangs,
      hasExperimentalData,
    }
  );

  if (candidates.length === 0) {
    // Try extended search
    const extendedCandidates = enumerateJunctionCandidates(
      seq,
      site,
      enzyme,
      {
        searchRadius: ENHANCED_JUNCTION_CONFIG.extendedSearchRadius,
        frame,
        codonUsage,
        existingOverhangs,
        hasExperimentalData,
      }
    );

    if (extendedCandidates.length === 0) {
      return {
        success: false,
        site,
        error: 'NO_VALID_JUNCTION',
        message: `No valid junction position found for site at ${siteStart}`,
        searchedRange: {
          start: Math.max(0, siteStart - ENHANCED_JUNCTION_CONFIG.extendedSearchRadius),
          end: Math.min(seq.length, siteEnd + ENHANCED_JUNCTION_CONFIG.extendedSearchRadius),
        },
      };
    }

    candidates.push(...extendedCandidates);
  }

  // Phase 2: Score all candidates using enhanced scoring
  const scoredCandidates = scoreJunctionCandidates(
    candidates,
    seq,
    enzyme,
    {
      existingOverhangs,
      targetTm,
      hasExperimentalData,
      enzymeData,
    }
  );

  // Phase 3: Select best candidate
  const best = scoredCandidates[0];

  // Phase 4: Design optimized primers for the junction
  const primerDesign = designOptimizedMutagenicPrimers(
    seq,
    best,
    enzyme,
    {
      targetTm,
      organism,
    }
  );

  // Phase 5: Validate and generate result
  const validationResult = validateJunctionDesign(
    best,
    primerDesign,
    existingOverhangs,
    enzyme
  );

  return {
    success: true,
    site,

    // Junction details
    junctionPosition: best.junctionPosition,
    overhang: best.overhang,

    // Quality metrics (ENHANCED)
    quality: {
      overall: best.overallScore!,
      tier: classifyQuality(best.overallScore!),
      breakdown: {
        overhangFidelity: best.scoring!.overhangFidelityScore,
        mutationQuality: best.scoring!.mutationScore,
        primerQuality: best.scoring!.primerScore,
        positionOptimality: best.scoring!.positionScore,
      },
    },

    // NEB fidelity data (ENHANCED)
    fidelity: {
      singleOverhang: best.overhangFidelity!,
      withExisting: best.assemblyFidelityWithExisting!,
      source: hasExperimentalData ? 'NEB_experimental' : 'calculated',
    },

    // Mutation details
    mutation: best.bestMutation,
    alternativeMutations: best.mutations.slice(1, 5),

    // Primer details (ENHANCED)
    primers: primerDesign,

    // Coverage analysis
    coverage: best.coverage,

    // Validation
    validation: validationResult,

    // Alternative options
    alternatives: scoredCandidates.slice(1, 10).map(c => ({
      junctionPosition: c.junctionPosition,
      overhang: c.overhang,
      overallScore: c.overallScore!,
      overhangFidelity: c.overhangFidelity!,
    })),

    message: generateJunctionMessage(best, primerDesign),
  };
}

// ============================================================================
// CANDIDATE ENUMERATION (ENHANCED)
// ============================================================================

/**
 * Enumerate all candidate junction positions with enhanced search
 */
function enumerateJunctionCandidates(
  sequence: string,
  site: InternalSite,
  enzyme: string,
  options: {
    searchRadius: number;
    frame: number;
    codonUsage: CodonUsageTable;
    existingOverhangs: string[];
    hasExperimentalData: boolean;
  }
): JunctionCandidate[] {
  const {
    searchRadius,
    frame,
    codonUsage,
    existingOverhangs,
    hasExperimentalData,
  } = options;

  const candidates: JunctionCandidate[] = [];
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  const recognition = site.sequence;
  const siteStart = site.position;
  const siteEnd = siteStart + recognition.length;
  const overhangLen = enz.overhangLength || 4;

  // Enhanced search range
  const searchStart = Math.max(0, siteStart - searchRadius);
  const searchEnd = Math.min(sequence.length - overhangLen, siteEnd + searchRadius);

  for (let junctionPos = searchStart; junctionPos <= searchEnd; junctionPos++) {
    // Extract overhang at this position
    const overhang = sequence.slice(junctionPos, junctionPos + overhangLen);

    // Skip palindromic overhangs (can self-ligate)
    if (isPalindrome(overhang)) continue;

    // Calculate primer coverage of the restriction site
    const coverage = calculatePrimerCoverage(
      junctionPos,
      siteStart,
      siteEnd,
      sequence.length,
      ENHANCED_JUNCTION_CONFIG.homologyLength.optimal
    );

    // Need at least one base of site covered by primers
    if (coverage.totalCoverage === 0) continue;

    // Find all possible silent mutations at covered positions
    const mutations = findSilentMutationsForCoverage(
      sequence,
      site,
      coverage,
      frame,
      enzyme,
      codonUsage
    );

    // Need at least one valid mutation
    if (mutations.length === 0) continue;

    // Calculate position optimality (prefer keeping codons intact)
    const positionOptimality = calculatePositionOptimality(
      junctionPos,
      siteStart,
      siteEnd,
      frame,
      overhangLen
    );

    candidates.push({
      junctionPosition: junctionPos,
      overhang,
      coverage,
      mutations,
      bestMutation: mutations[0],
      positionOptimality,
      // Will be scored in next phase
      overhangFidelity: null,
      primerQuality: null,
      overallScore: null,
    });
  }

  return candidates;
}

/**
 * Calculate which parts of the restriction site are covered by primers
 */
function calculatePrimerCoverage(
  junctionPos: number,
  siteStart: number,
  siteEnd: number,
  seqLen: number,
  homologyLen: number
): PrimerCoverage {
  const overhangLen = 4;

  // Fragment 1 (upstream) reverse primer coverage
  const frag1HomologyEnd = junctionPos + overhangLen;
  const frag1HomologyStart = Math.max(0, frag1HomologyEnd - homologyLen);

  // Fragment 2 (downstream) forward primer coverage
  const frag2HomologyStart = junctionPos;
  const frag2HomologyEnd = Math.min(seqLen, junctionPos + homologyLen);

  // Which site bases are covered by which primer
  const frag1Coverage: Array<{ position: number; indexInSite: number }> = [];
  const frag2Coverage: Array<{ position: number; indexInSite: number }> = [];

  for (let pos = siteStart; pos < siteEnd; pos++) {
    const indexInSite = pos - siteStart;

    if (pos >= frag1HomologyStart && pos < frag1HomologyEnd) {
      frag1Coverage.push({ position: pos, indexInSite });
    }
    if (pos >= frag2HomologyStart && pos < frag2HomologyEnd) {
      frag2Coverage.push({ position: pos, indexInSite });
    }
  }

  // Unique positions covered
  const allPositions = new Set([
    ...frag1Coverage.map(c => c.position),
    ...frag2Coverage.map(c => c.position),
  ]);

  return {
    fragment1Primer: frag1Coverage,
    fragment2Primer: frag2Coverage,
    totalCoverage: allPositions.size,
    fragment1HomologyRange: { start: frag1HomologyStart, end: frag1HomologyEnd },
    fragment2HomologyRange: { start: frag2HomologyStart, end: frag2HomologyEnd },
  };
}

/**
 * Calculate position optimality (prefer keeping codons intact)
 */
function calculatePositionOptimality(
  junctionPos: number,
  siteStart: number,
  siteEnd: number,
  frame: number,
  overhangLen: number
): number {
  let score = 100;

  // Prefer junctions at codon boundaries
  const adjustedPos = (junctionPos - frame);
  const posInCodon = ((adjustedPos % 3) + 3) % 3;

  if (posInCodon === 0) {
    // Junction at codon boundary - ideal
    score += 10;
  } else if (posInCodon === 2) {
    // Junction after wobble position - acceptable
    score += 5;
  }
  // posInCodon === 1 is worst (splits codon in middle)

  // Prefer junctions closer to the site (less disruption)
  const siteCenter = (siteStart + siteEnd) / 2;
  const distanceFromSite = Math.abs(junctionPos + overhangLen/2 - siteCenter);
  const maxDistance = 50;
  const distancePenalty = Math.min(20, (distanceFromSite / maxDistance) * 20);
  score -= distancePenalty;

  return Math.max(0, Math.min(100, score));
}

// ============================================================================
// SILENT MUTATION FINDING (ENHANCED)
// ============================================================================

/**
 * Find silent mutations for positions covered by primers
 */
function findSilentMutationsForCoverage(
  sequence: string,
  site: InternalSite,
  coverage: PrimerCoverage,
  frame: number,
  enzyme: string,
  codonUsage: CodonUsageTable
): Mutation[] {
  const mutations: Mutation[] = [];
  const recognition = site.sequence;
  const siteStart = site.position;
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  const enzRecognition = enz.recognition;
  const enzRecognitionRC = reverseComplement(enzRecognition);

  // Combine all covered positions
  const allCovered = new Map<number, { position: number; indexInSite: number; inFrag1: boolean; inFrag2: boolean }>();

  for (const base of coverage.fragment1Primer) {
    if (!allCovered.has(base.position)) {
      allCovered.set(base.position, {
        ...base,
        inFrag1: true,
        inFrag2: false,
      });
    } else {
      allCovered.get(base.position)!.inFrag1 = true;
    }
  }

  for (const base of coverage.fragment2Primer) {
    if (!allCovered.has(base.position)) {
      allCovered.set(base.position, {
        ...base,
        inFrag1: false,
        inFrag2: true,
      });
    } else {
      allCovered.get(base.position)!.inFrag2 = true;
    }
  }

  // For each covered position, try all possible mutations
  for (const [seqPos, baseInfo] of allCovered) {
    const originalBase = sequence[seqPos];
    const posInSite = baseInfo.indexInSite;

    // Determine codon context
    const adjustedPos = seqPos - frame;
    if (adjustedPos < 0) continue;

    const codonIndex = Math.floor(adjustedPos / 3);
    const codonStart = codonIndex * 3 + frame;
    const posInCodon = seqPos - codonStart;

    if (codonStart < 0 || codonStart + 3 > sequence.length) continue;
    if (posInCodon < 0 || posInCodon > 2) continue;

    const originalCodon = sequence.slice(codonStart, codonStart + 3);
    const originalAA = CODON_TO_AA[originalCodon];

    if (!originalAA) continue;

    // Try each alternative base
    for (const newBase of ['A', 'T', 'G', 'C']) {
      if (newBase === originalBase) continue;

      // Create new codon
      const newCodon =
        originalCodon.slice(0, posInCodon) +
        newBase +
        originalCodon.slice(posInCodon + 1);

      const newAA = CODON_TO_AA[newCodon];

      // Must be synonymous
      if (newAA !== originalAA) continue;

      // Check if mutation breaks the site
      const newSiteSeq =
        recognition.slice(0, posInSite) +
        newBase +
        recognition.slice(posInSite + 1);

      const breaksForward = newSiteSeq !== enzRecognition;
      const breaksReverse = newSiteSeq !== enzRecognitionRC;
      const siteIsBroken = site.orientation === 'forward' ? breaksForward : breaksReverse;

      if (!siteIsBroken) continue;

      // Score this mutation (ENHANCED scoring)
      const mutationScore = scoreMutation(
        originalCodon,
        newCodon,
        posInCodon,
        posInSite,
        recognition.length,
        codonUsage
      );

      mutations.push({
        sequencePosition: seqPos,
        positionInSite: posInSite,
        positionInCodon: posInCodon,
        originalBase,
        newBase,
        originalCodon,
        newCodon,
        aminoAcid: originalAA,
        codonFrequency: codonUsage[newCodon] || 0,
        originalCodonFrequency: codonUsage[originalCodon] || 0,
        inFragment1Primer: baseInfo.inFrag1,
        inFragment2Primer: baseInfo.inFrag2,
        score: mutationScore.score,
        scoreBreakdown: mutationScore.breakdown,
        isSynonymous: true,
        breaksSite: true,
      });
    }
  }

  // Sort by score
  mutations.sort((a, b) => b.score - a.score);

  return mutations;
}

/**
 * Score a silent mutation (ENHANCED)
 */
function scoreMutation(
  originalCodon: string,
  newCodon: string,
  posInCodon: number,
  posInSite: number,
  siteLen: number,
  codonUsage: CodonUsageTable
): { score: number; breakdown: MutationScoreBreakdown } {
  const breakdown: MutationScoreBreakdown = {};
  let score = 80;

  const newFreq = codonUsage[newCodon] || 0;
  const originalFreq = codonUsage[originalCodon] || 0;

  // Wobble position preference (3rd position safer)
  if (posInCodon === 2) {
    score += 15;
    breakdown.wobbleBonus = 15;
  } else if (posInCodon === 1) {
    score -= 5;
    breakdown.middlePositionPenalty = -5;
  }

  // Codon frequency scoring
  if (newFreq < ENHANCED_JUNCTION_CONFIG.rareCodoThreshold) {
    const penalty = Math.min(25, (ENHANCED_JUNCTION_CONFIG.rareCodoThreshold - newFreq) * 2);
    score -= penalty;
    breakdown.rareCodonPenalty = -penalty;
  } else if (newFreq >= originalFreq * 0.8) {
    score += 10;
    breakdown.frequencyPreservedBonus = 10;
  }

  // Position in site preference (middle is more robust)
  const distFromMiddle = Math.abs(posInSite - (siteLen - 1) / 2);
  const positionBonus = (1 - distFromMiddle / (siteLen / 2)) * 5;
  score += positionBonus;
  breakdown.positionBonus = positionBonus;

  return {
    score: Math.max(0, Math.min(100, score)),
    breakdown,
  };
}

// ============================================================================
// CANDIDATE SCORING (ENHANCED with NEB data)
// ============================================================================

/**
 * Score all junction candidates using NEB experimental fidelity data
 */
function scoreJunctionCandidates(
  candidates: JunctionCandidate[],
  sequence: string,
  enzyme: string,
  options: {
    existingOverhangs: string[];
    targetTm: number;
    hasExperimentalData: boolean;
    enzymeData: any;
  }
): JunctionCandidate[] {
  const {
    existingOverhangs,
    targetTm,
    hasExperimentalData,
    enzymeData,
  } = options;

  const weights = ENHANCED_JUNCTION_CONFIG.weights;

  for (const candidate of candidates) {
    // 1. Score overhang fidelity using NEB data (CRITICAL)
    const fidelityResult = scoreOverhangFidelity(
      candidate.overhang,
      existingOverhangs,
      enzyme,
      hasExperimentalData,
      enzymeData
    );
    candidate.overhangFidelity = fidelityResult.singleFidelity;
    candidate.assemblyFidelityWithExisting = fidelityResult.assemblyFidelity;

    // 2. Score mutation quality
    const mutationScore = candidate.bestMutation.score;

    // 3. Score primer quality (thermodynamics)
    const primerScore = scorePrimerQuality(
      sequence,
      candidate.junctionPosition,
      candidate.coverage,
      candidate.bestMutation,
      targetTm
    );
    candidate.primerQuality = primerScore;

    // 4. Position optimality (already calculated)
    const positionScore = candidate.positionOptimality;

    // Store scoring breakdown
    candidate.scoring = {
      overhangFidelityScore: fidelityResult.score,
      mutationScore,
      primerScore: primerScore.score,
      positionScore,
    };

    // Calculate weighted overall score
    candidate.overallScore = Math.round(
      weights.overhangFidelity * fidelityResult.score +
      weights.mutationQuality * mutationScore +
      weights.primerQuality * primerScore.score +
      weights.positionOptimality * positionScore
    );
  }

  // Sort by overall score
  candidates.sort((a, b) => b.overallScore! - a.overallScore!);

  return candidates;
}

/**
 * Score overhang fidelity using NEB experimental data
 */
export function scoreOverhangFidelity(
  overhang: string,
  existingOverhangs: string[],
  enzyme: string,
  hasExperimentalData: boolean,
  enzymeData: any
): OverhangFidelityScore {
  const ohUpper = overhang.toUpperCase();

  // Get single overhang fidelity from NEB data
  let singleFidelity = 0.85; // default
  let score = 70;

  if (hasExperimentalData) {
    const fidelityData = getOverhangFidelityExperimental(ohUpper, enzyme);
    singleFidelity = fidelityData.fidelity;

    // Convert fidelity to score (0.95+ = 100, 0.80 = 60)
    if (singleFidelity >= 0.99) {
      score = 100;
    } else if (singleFidelity >= 0.95) {
      score = 90 + (singleFidelity - 0.95) * 200;
    } else if (singleFidelity >= 0.90) {
      score = 80 + (singleFidelity - 0.90) * 200;
    } else if (singleFidelity >= 0.80) {
      score = 60 + (singleFidelity - 0.80) * 200;
    } else {
      score = Math.max(0, singleFidelity * 75);
    }
  } else {
    // Fallback: use heuristic scoring
    score = scoreOverhangHeuristic(ohUpper);
    singleFidelity = score / 100;
  }

  // Calculate assembly fidelity with existing overhangs
  let assemblyFidelity = singleFidelity;

  if (existingOverhangs.length > 0 && hasExperimentalData) {
    const allOverhangs = [...existingOverhangs, ohUpper];
    const fidelityResult = calculateExperimentalFidelity(allOverhangs, enzyme);
    assemblyFidelity = fidelityResult.assemblyFidelity;

    // Penalize if new overhang reduces assembly fidelity significantly
    if (assemblyFidelity < 0.90) {
      score -= 20;
    } else if (assemblyFidelity < 0.95) {
      score -= 10;
    }
  }

  // Check for conflicts with existing overhangs
  for (const existing of existingOverhangs) {
    if (ohUpper === existing.toUpperCase() ||
        ohUpper === reverseComplement(existing.toUpperCase())) {
      score -= 50; // Major penalty for duplicate
      break;
    }
  }

  return {
    score: Math.max(0, Math.min(100, score)),
    singleFidelity,
    assemblyFidelity,
  };
}

/**
 * Heuristic overhang scoring (fallback when no NEB data)
 */
function scoreOverhangHeuristic(overhang: string): number {
  let score = 70;

  // Palindrome check
  if (isPalindrome(overhang)) {
    score -= 40;
  }

  // Homopolymer check
  if (/^(.)\1{3}$/.test(overhang)) {
    score -= 30;
  } else if (/(.)\1{2}/.test(overhang)) {
    score -= 15;
  }

  // GC content
  const gc = (overhang.match(/[GC]/g) || []).length / overhang.length;
  if (gc === 0.5) {
    score += 15;
  } else if (gc < 0.25 || gc > 0.75) {
    score -= 15;
  }

  // TNNA pattern bonus (high efficiency)
  if (/^T..A$/.test(overhang)) {
    score += 10;
  }

  return Math.max(0, Math.min(100, score));
}

/**
 * Score primer quality using thermodynamic analysis
 */
export function scorePrimerQuality(
  sequence: string,
  junctionPos: number,
  coverage: PrimerCoverage,
  mutation: Mutation,
  targetTm: number
): PrimerQualityScore {
  const overhangLen = 4;
  const homologyLen = ENHANCED_JUNCTION_CONFIG.homologyLength.optimal;

  // Get homology regions
  const frag1HomologyStart = Math.max(0, junctionPos + overhangLen - homologyLen);
  const frag1HomologyEnd = junctionPos + overhangLen;
  let frag1Homology = sequence.slice(frag1HomologyStart, frag1HomologyEnd);

  const frag2HomologyStart = junctionPos;
  const frag2HomologyEnd = Math.min(sequence.length, junctionPos + homologyLen);
  let frag2Homology = sequence.slice(frag2HomologyStart, frag2HomologyEnd);

  // Apply mutation to homology
  const mutPos = mutation.sequencePosition;

  if (mutation.inFragment1Primer && mutPos >= frag1HomologyStart && mutPos < frag1HomologyEnd) {
    const relPos = mutPos - frag1HomologyStart;
    frag1Homology = frag1Homology.slice(0, relPos) + mutation.newBase + frag1Homology.slice(relPos + 1);
  }

  if (mutation.inFragment2Primer && mutPos >= frag2HomologyStart && mutPos < frag2HomologyEnd) {
    const relPos = mutPos - frag2HomologyStart;
    frag2Homology = frag2Homology.slice(0, relPos) + mutation.newBase + frag2Homology.slice(relPos + 1);
  }

  // Score each homology region
  const frag1Score = analyzeHomologyRegion(frag1Homology, false, targetTm);
  const frag2Score = analyzeHomologyRegion(frag2Homology, true, targetTm);

  // Combined score
  const combinedScore = (frag1Score.score + frag2Score.score) / 2;

  // Tm difference penalty
  const tmDiff = Math.abs(frag1Score.tm - frag2Score.tm);
  let tmPenalty = 0;
  if (tmDiff > 5) {
    tmPenalty = Math.min(20, (tmDiff - 5) * 4);
  }

  return {
    score: Math.max(0, combinedScore - tmPenalty),
    fragment1: frag1Score,
    fragment2: frag2Score,
    tmDifference: tmDiff,
  };
}

/**
 * Analyze a homology region for primer quality
 */
function analyzeHomologyRegion(region: string, isForward: boolean, targetTm: number): HomologyAnalysis {
  if (!region || region.length < 10) {
    return { score: 30, tm: 0, gc: 0, issues: ['Too short'], length: region?.length || 0, gcClamp: false };
  }

  const seq = region.toUpperCase();
  const issues: string[] = [];

  // Calculate Tm
  let tm = 62;
  try {
    tm = calculateTmQ5(seq);
  } catch (e) {
    // Fallback
    const gc = (seq.match(/[GC]/g) || []).length;
    tm = 64.9 + 41 * (gc - 16.4) / seq.length;
  }

  // GC content
  const gcCount = (seq.match(/[GC]/g) || []).length;
  const gc = (gcCount / seq.length) * 100;

  // Score Tm (target ±5°C is ideal)
  const tmDiff = Math.abs(tm - targetTm);
  let tmScore = 100;
  if (tmDiff > 10) {
    tmScore = 50;
    issues.push(`Tm ${tm.toFixed(1)}°C far from target`);
  } else if (tmDiff > 5) {
    tmScore = 75;
    issues.push(`Tm ${tm.toFixed(1)}°C outside ideal range`);
  }

  // Score GC (40-60% ideal)
  let gcScore = 100;
  if (gc < 30 || gc > 70) {
    gcScore = 50;
    issues.push(`GC ${gc.toFixed(0)}% outside acceptable range`);
  } else if (gc < 40 || gc > 60) {
    gcScore = 75;
  }

  // Check for hairpin
  let hairpinScore = 100;
  try {
    const hairpinDG = calculateHairpinDG(seq);
    if (hairpinDG < -5) {
      hairpinScore = 50;
      issues.push('Strong hairpin potential');
    } else if (hairpinDG < -3) {
      hairpinScore = 75;
    }
  } catch (e) {}

  // Check GC clamp
  const last2 = seq.slice(-2);
  const gcClamp = (last2.match(/[GC]/g) || []).length;
  let clampScore = 100;
  if (gcClamp === 0) {
    clampScore = 60;
    issues.push('No GC clamp');
  } else if (gcClamp === 2) {
    clampScore = 100;
  } else {
    clampScore = 80;
  }

  // Check G-quadruplex
  let gQuadScore = 100;
  if (/GGG.{1,7}GGG.{1,7}GGG.{1,7}GGG/.test(seq)) {
    gQuadScore = 40;
    issues.push('G-quadruplex risk');
  }

  // Weighted composite
  const score = (
    tmScore * 0.30 +
    gcScore * 0.20 +
    hairpinScore * 0.20 +
    clampScore * 0.15 +
    gQuadScore * 0.15
  );

  return {
    score: Math.round(score),
    tm,
    gc,
    length: seq.length,
    gcClamp: gcClamp > 0,
    issues,
  };
}

// ============================================================================
// PRIMER DESIGN (ENHANCED)
// ============================================================================

/**
 * Calculate optimal homology length to hit target Tm
 */
function calculateAdaptiveHomologyLength(
  sequence: string,
  startPos: number,
  direction: 'forward' | 'reverse',
  targetTm: number,
  config: { min: number; max: number; optimal: number }
): number {
  const { min, max, optimal } = config;

  let bestLen = optimal;
  let bestTmDiff = Infinity;

  for (let len = min; len <= max; len++) {
    let region;
    if (direction === 'forward') {
      region = sequence.slice(startPos, Math.min(sequence.length, startPos + len));
    } else {
      region = sequence.slice(Math.max(0, startPos - len), startPos);
    }

    if (region.length < min) continue;

    let tm = 62;
    try {
      tm = calculateTmQ5(region);
    } catch (e) {
      // Fallback
      const gc = (region.match(/[GC]/g) || []).length;
      tm = 64.9 + 41 * (gc - 16.4) / region.length;
    }

    const tmDiff = Math.abs(tm - targetTm);
    if (tmDiff < bestTmDiff) {
      bestTmDiff = tmDiff;
      bestLen = len;
    }

    // Stop early if we hit target within tolerance
    if (tmDiff <= ENHANCED_JUNCTION_CONFIG.tmTolerance) {
      break;
    }
  }

  return bestLen;
}

/**
 * Design optimized mutagenic primers
 */
function designOptimizedMutagenicPrimers(
  sequence: string,
  candidate: JunctionCandidate,
  enzyme: string,
  options: { targetTm: number; organism: string }
): PrimerDesign {
  const { targetTm, organism } = options;
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  const recognition = enz.recognition;
  const overhangLen = candidate.overhang.length;
  const junctionPos = candidate.junctionPosition;
  const mutation = candidate.bestMutation;

  // Get optimal flanking sequence for enzyme
  const flankingConfig = OPTIMAL_FLANKING_SEQUENCES[enzyme] || OPTIMAL_FLANKING_SEQUENCES.BsaI;
  const flanking = flankingConfig.default || 'GGTGCG';

  // Get optimal spacer
  const spacerConfig = OPTIMAL_SPACERS[enzyme] || OPTIMAL_SPACERS.BsaI;
  const spacer = (spacerConfig as any).default || 'A';  // FIXED: Type assertion for default property

  // Calculate adaptive homology length based on Tm
  const homologyConfig = ENHANCED_JUNCTION_CONFIG.homologyLength;
  let frag1HomologyLen = homologyConfig.optimal;
  let frag2HomologyLen = homologyConfig.optimal;

  if (ENHANCED_JUNCTION_CONFIG.useAdaptiveHomology) {
    frag1HomologyLen = calculateAdaptiveHomologyLength(
      sequence,
      junctionPos + overhangLen,
      'reverse',
      targetTm,
      homologyConfig
    ) as any;  // FIXED: Type assertion for number assignment
    frag2HomologyLen = calculateAdaptiveHomologyLength(
      sequence,
      junctionPos,
      'forward',
      targetTm,
      homologyConfig
    ) as any;  // FIXED: Type assertion for number assignment
  }

  // Fragment 1 reverse primer homology
  const frag1HomologyStart = Math.max(0, junctionPos + overhangLen - frag1HomologyLen);
  const frag1HomologyEnd = junctionPos + overhangLen;
  let frag1Homology = sequence.slice(frag1HomologyStart, frag1HomologyEnd);

  // Fragment 2 forward primer homology
  const frag2HomologyStart = junctionPos;
  const frag2HomologyEnd = Math.min(sequence.length, junctionPos + frag2HomologyLen);
  let frag2Homology = sequence.slice(frag2HomologyStart, frag2HomologyEnd);

  // Apply mutation
  const mutPos = mutation.sequencePosition;

  if (mutation.inFragment1Primer && mutPos >= frag1HomologyStart && mutPos < frag1HomologyEnd) {
    const relPos = mutPos - frag1HomologyStart;
    frag1Homology = frag1Homology.slice(0, relPos) + mutation.newBase + frag1Homology.slice(relPos + 1);
  }

  if (mutation.inFragment2Primer && mutPos >= frag2HomologyStart && mutPos < frag2HomologyEnd) {
    const relPos = mutPos - frag2HomologyStart;
    frag2Homology = frag2Homology.slice(0, relPos) + mutation.newBase + frag2Homology.slice(relPos + 1);
  }

  // Build full primer sequences
  const fragment2ForwardPrimer: PrimerInfo = {
    sequence: flanking + recognition + spacer + candidate.overhang + frag2Homology,
    components: {
      flanking,
      recognition,
      spacer,
      overhang: candidate.overhang,
      homology: frag2Homology,
    },
    hasMutation: mutation.inFragment2Primer,
    mutationPosition: mutation.inFragment2Primer ?
      (flanking.length + recognition.length + spacer.length + overhangLen + (mutPos - frag2HomologyStart)) : null,
    length: 0,
    homologyLength: frag2Homology.length,
  };
  fragment2ForwardPrimer.length = fragment2ForwardPrimer.sequence.length;

  const fragment1ReversePrimer: PrimerInfo = {
    sequence: flanking + recognition + spacer +
              reverseComplement(candidate.overhang) + reverseComplement(frag1Homology),
    components: {
      flanking,
      recognition,
      spacer,
      overhang: reverseComplement(candidate.overhang),
      homology: reverseComplement(frag1Homology),
    },
    hasMutation: mutation.inFragment1Primer,
    mutationPosition: mutation.inFragment1Primer ?
      (flanking.length + recognition.length + spacer.length + overhangLen +
       (frag1HomologyEnd - 1 - mutPos)) : null,
    length: 0,
    homologyLength: frag1Homology.length,
  };
  fragment1ReversePrimer.length = fragment1ReversePrimer.sequence.length;

  // Calculate primer properties
  let fwdTm = 0, revTm = 0;
  try {
    fwdTm = calculateTmQ5(fragment2ForwardPrimer.components.homology);
    revTm = calculateTmQ5(fragment1ReversePrimer.components.homology);
  } catch (e) {}

  fragment2ForwardPrimer.tm = fwdTm;
  fragment1ReversePrimer.tm = revTm;

  return {
    fragment1: {
      reversePrimer: fragment1ReversePrimer,
      note: 'Reverse primer for Fragment 1 (upstream of junction)',
    },
    fragment2: {
      forwardPrimer: fragment2ForwardPrimer,
      note: 'Forward primer for Fragment 2 (downstream of junction)',
    },
    mutation: {
      applied: true,
      inFragment1Primer: mutation.inFragment1Primer,
      inFragment2Primer: mutation.inFragment2Primer,
      change: `${mutation.originalBase}→${mutation.newBase} at position ${mutPos}`,
      codonChange: `${mutation.originalCodon}→${mutation.newCodon} (${mutation.aminoAcid})`,
    },
    tmDifference: Math.abs(fwdTm - revTm),
    optimizedFlanking: flanking,
  };
}

// ============================================================================
// VALIDATION AND UTILITIES
// ============================================================================

/**
 * Validate the junction design
 */
function validateJunctionDesign(
  candidate: JunctionCandidate,
  primerDesign: PrimerDesign,
  existingOverhangs: string[],
  enzyme: string
): JunctionValidation {
  const warnings: ValidationWarning[] = [];
  const errors: ValidationWarning[] = [];

  // Check overhang fidelity
  if (candidate.overhangFidelity! < 0.90) {
    warnings.push({
      code: 'LOW_OVERHANG_FIDELITY',
      message: `Overhang fidelity is ${(candidate.overhangFidelity! * 100).toFixed(1)}% (target: >95%)`,
    });
  }

  // Check assembly fidelity
  if (candidate.assemblyFidelityWithExisting! < 0.90) {
    warnings.push({
      code: 'LOW_ASSEMBLY_FIDELITY',
      message: `Assembly fidelity with existing overhangs is ${(candidate.assemblyFidelityWithExisting! * 100).toFixed(1)}%`,
    });
  }

  // Check Tm difference
  if (primerDesign.tmDifference > 5) {
    warnings.push({
      code: 'HIGH_TM_DIFF',
      message: `Primer Tm difference is ${primerDesign.tmDifference.toFixed(1)}°C (target: <5°C)`,
    });
  }

  // Check primer quality
  if (candidate.primerQuality!.score < 60) {
    warnings.push({
      code: 'LOW_PRIMER_QUALITY',
      message: 'Primer quality is below optimal',
      details: candidate.primerQuality,
    });
  }

  // Check mutation score
  if (candidate.bestMutation.score < 60) {
    warnings.push({
      code: 'LOW_MUTATION_QUALITY',
      message: 'Mutation uses rare codon or non-optimal position',
      details: candidate.bestMutation.scoreBreakdown,
    });
  }

  return {
    isValid: errors.length === 0,
    warnings,
    errors,
    overallQuality: classifyQuality(candidate.overallScore!),
  };
}

/**
 * Check if overhang is palindromic
 */
function isPalindrome(overhang: string): boolean {
  return overhang.toUpperCase() === reverseComplement(overhang.toUpperCase());
}

/**
 * Classify quality tier
 */
export function classifyQuality(score: number): QualityTier {
  const thresholds = ENHANCED_JUNCTION_CONFIG.qualityThresholds;
  if (score >= thresholds.excellent) return 'excellent';
  if (score >= thresholds.good) return 'good';
  if (score >= thresholds.acceptable) return 'acceptable';
  return 'poor';
}

/**
 * Generate human-readable message
 */
function generateJunctionMessage(candidate: JunctionCandidate, primerDesign: PrimerDesign): string {
  const quality = classifyQuality(candidate.overallScore!);
  const mut = candidate.bestMutation;

  return `${quality.charAt(0).toUpperCase() + quality.slice(1)} junction at position ${candidate.junctionPosition} ` +
         `with overhang ${candidate.overhang} (${(candidate.overhangFidelity! * 100).toFixed(1)}% fidelity). ` +
         `Silent mutation: ${mut.originalBase}→${mut.newBase} (${mut.originalCodon}→${mut.newCodon}, ${mut.aminoAcid}).`;
}

// ============================================================================
// BATCH PROCESSING
// ============================================================================

/**
 * Design enhanced mutagenic junctions for all internal sites
 * with global overhang optimization
 */
export function designAllEnhancedJunctions(
  sequence: string,
  enzyme: string = 'BsaI',
  options: {
    frame?: number;
    organism?: string;
    existingOverhangs?: string[];
    optimizeGlobalFidelity?: boolean;
  } = {}
): BatchJunctionResult {
  const {
    frame = 0,
    organism = 'ecoli',
    existingOverhangs = [],
    optimizeGlobalFidelity = true,
  } = options;

  const internalSites = findInternalSites(sequence, enzyme);

  if (!internalSites.hasSites) {
    return {
      success: true,
      needsDomestication: false,
      totalSites: 0,
      sites: [],
      junctions: [],
      failedSites: [],
      allOverhangs: existingOverhangs,
      fidelity: { assembly: 1.0, source: 'NEB_experimental' },
      additionalFragments: 0,
      optimized: false,
      message: `No internal ${enzyme} sites found`,
    };
  }

  const junctions: EnhancedJunctionResult[] = [];
  const failedSites: Array<{ site: InternalSite; error: string; message: string }> = [];
  const allOverhangs = [...existingOverhangs];

  // Design junction for each site
  for (const site of internalSites.sites) {
    const result = designEnhancedMutagenicJunction(sequence, site, enzyme, {
      frame,
      organism,
      existingOverhangs: allOverhangs,
    });

    if (result.success) {
      junctions.push(result);
      allOverhangs.push(result.overhang!);
    } else {
      failedSites.push({
        site,
        error: result.error || 'UNKNOWN',
        message: result.message || 'Failed to design junction',
      });
    }
  }

  // Calculate overall assembly fidelity
  let fidelityResult: any = null;
  if (junctions.length > 0) {
    try {
      fidelityResult = calculateExperimentalFidelity(allOverhangs, enzyme);
    } catch (e) {}
  }

  // Optionally optimize global overhang set
  let optimizedResult: any = null;
  if (optimizeGlobalFidelity && junctions.length > 1 && fidelityResult?.assemblyFidelity < 0.95) {
    optimizedResult = optimizeGlobalOverhangSelection(
      sequence,
      internalSites.sites,
      enzyme,
      { frame, organism, existingOverhangs }
    );
  }

  return {
    success: failedSites.length === 0,
    needsDomestication: true,
    totalSites: internalSites.count,
    sites: internalSites.sites,
    junctions: optimizedResult?.junctions || junctions,
    failedSites,
    allOverhangs: optimizedResult?.allOverhangs || allOverhangs,
    fidelity: {
      assembly: optimizedResult?.assemblyFidelity || fidelityResult?.assemblyFidelity,
      source: 'NEB_experimental',
    },
    additionalFragments: junctions.length,
    optimized: !!optimizedResult,
    message: generateBatchMessage(junctions, failedSites, internalSites.count),
  };
}

/**
 * Optimize global overhang selection for maximum assembly fidelity
 */
function optimizeGlobalOverhangSelection(
  sequence: string,
  sites: InternalSite[],
  enzyme: string,
  options: { frame: number; organism: string; existingOverhangs: string[] }
): { junctions: EnhancedJunctionResult[]; allOverhangs: string[]; assemblyFidelity: number } | null {
  const { frame, organism, existingOverhangs } = options;

  // Get all valid junction positions for each site
  const siteCandidates = sites.map(site => {
    const candidates = enumerateJunctionCandidates(
      sequence,
      site,
      enzyme,
      {
        searchRadius: ENHANCED_JUNCTION_CONFIG.searchRadius,
        frame,
        codonUsage: organism === 'yeast' ? YEAST_CODON_USAGE : ECOLI_CODON_USAGE,
        existingOverhangs,
        hasExperimentalData: true,
      }
    );
    return { site, candidates };
  });

  // Try to find combination with highest overall fidelity
  // For small number of sites, try exhaustive search
  if (sites.length <= 4) {
    return exhaustiveOverhangSearch(siteCandidates, sequence, enzyme, options);
  }

  // For larger numbers, use greedy approach with backtracking
  return greedyOverhangSearch(siteCandidates, sequence, enzyme, options);
}

/**
 * Exhaustive search for optimal overhang combination (small assemblies)
 */
function exhaustiveOverhangSearch(
  siteCandidates: Array<{ site: InternalSite; candidates: JunctionCandidate[] }>,
  sequence: string,
  enzyme: string,
  options: { existingOverhangs: string[]; frame: number; organism: string }
): { junctions: EnhancedJunctionResult[]; allOverhangs: string[]; assemblyFidelity: number } | null {
  const { existingOverhangs, frame, organism } = options;

  // Limit candidates per site for tractability
  const maxCandidatesPerSite = 10;
  const limitedCandidates = siteCandidates.map(sc => ({
    site: sc.site,
    candidates: sc.candidates.slice(0, maxCandidatesPerSite),
  }));

  let bestFidelity = 0;
  let bestCombination: { overhangs: string[]; junctions: Array<{ site: InternalSite; candidate: JunctionCandidate }> } | null = null;

  // Recursive enumeration
  function enumerate(
    siteIndex: number,
    currentOverhangs: string[],
    currentJunctions: Array<{ site: InternalSite; candidate: JunctionCandidate }>
  ): void {
    if (siteIndex >= limitedCandidates.length) {
      // Evaluate this combination
      const allOverhangs = [...existingOverhangs, ...currentOverhangs];
      try {
        const fidelity = calculateExperimentalFidelity(allOverhangs, enzyme);
        if (fidelity.assemblyFidelity > bestFidelity) {
          bestFidelity = fidelity.assemblyFidelity;
          bestCombination = {
            overhangs: [...currentOverhangs],
            junctions: [...currentJunctions],
          };
        }
      } catch (e) {}
      return;
    }

    const { candidates, site } = limitedCandidates[siteIndex];

    for (const candidate of candidates) {
      // Skip if overhang conflicts with existing
      const hasConflict = currentOverhangs.some(oh =>
        oh === candidate.overhang ||
        oh === reverseComplement(candidate.overhang)
      );
      if (hasConflict) continue;

      enumerate(
        siteIndex + 1,
        [...currentOverhangs, candidate.overhang],
        [...currentJunctions, { site, candidate }]
      );
    }
  }

  enumerate(0, [], []);

  if (!bestCombination) {
    return null;
  }

  // Build full junction designs for best combination
  const combo = bestCombination as { overhangs: string[]; junctions: Array<{ site: InternalSite; candidate: JunctionCandidate }> };
  const junctions = combo.junctions.map(({ site, candidate }: { site: InternalSite; candidate: JunctionCandidate }) => {
    return designEnhancedMutagenicJunction(sequence, site, enzyme, {
      frame,
      organism,
      existingOverhangs: [...existingOverhangs, ...combo.overhangs.filter(oh => oh !== candidate.overhang)],
    });
  });

  return {
    junctions,
    allOverhangs: [...existingOverhangs, ...combo.overhangs],
    assemblyFidelity: bestFidelity,
  };
}

/**
 * Greedy search with backtracking (larger assemblies)
 */
function greedyOverhangSearch(
  siteCandidates: Array<{ site: InternalSite; candidates: JunctionCandidate[] }>,
  sequence: string,
  enzyme: string,
  options: { existingOverhangs: string[]; frame: number; organism: string }
): { junctions: EnhancedJunctionResult[]; allOverhangs: string[]; assemblyFidelity: number } | null {
  const { existingOverhangs, frame, organism } = options;

  const selectedOverhangs: string[] = [];
  const selectedJunctions: EnhancedJunctionResult[] = [];

  for (const { site, candidates } of siteCandidates) {
    // Score each candidate considering already selected overhangs
    const allExisting = [...existingOverhangs, ...selectedOverhangs];

    let bestCandidate: JunctionCandidate | null = null;
    let bestFidelity = 0;

    for (const candidate of candidates.slice(0, 20)) {
      // Skip conflicts
      const hasConflict = allExisting.some(oh =>
        oh === candidate.overhang ||
        oh === reverseComplement(candidate.overhang)
      );
      if (hasConflict) continue;

      // Calculate fidelity with this overhang
      const testOverhangs = [...allExisting, candidate.overhang];
      try {
        const fidelity = calculateExperimentalFidelity(testOverhangs, enzyme);
        if (fidelity.assemblyFidelity > bestFidelity) {
          bestFidelity = fidelity.assemblyFidelity;
          bestCandidate = candidate;
        }
      } catch (e) {}
    }

    if (bestCandidate) {
      selectedOverhangs.push(bestCandidate.overhang);

      const junction = designEnhancedMutagenicJunction(sequence, site, enzyme, {
        frame,
        organism,
        existingOverhangs: allExisting,
      });
      selectedJunctions.push(junction);
    }
  }

  const allOverhangs = [...existingOverhangs, ...selectedOverhangs];
  let finalFidelity = 0;
  try {
    finalFidelity = calculateExperimentalFidelity(allOverhangs, enzyme).assemblyFidelity;
  } catch (e) {}

  return {
    junctions: selectedJunctions,
    allOverhangs,
    assemblyFidelity: finalFidelity,
  };
}

/**
 * Generate batch message
 */
function generateBatchMessage(
  junctions: EnhancedJunctionResult[],
  failedSites: Array<{ site: InternalSite; error: string; message: string }>,
  totalSites: number
): string {
  const parts: string[] = [];

  if (junctions.length > 0) {
    const avgQuality = junctions.reduce((sum, j) => sum + (j.quality?.overall || 0), 0) / junctions.length;
    parts.push(`Designed ${junctions.length} enhanced junction(s) (avg quality: ${avgQuality.toFixed(0)})`);
  }

  if (failedSites.length > 0) {
    parts.push(`${failedSites.length} site(s) could not be domesticated`);
  }

  return parts.join('. ') + '.';
}
