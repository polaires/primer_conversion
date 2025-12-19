/**
 * Silent Mutation-Based Domestication for Golden Gate Assembly
 *
 * This module implements a scientifically correct domestication strategy that uses
 * silent (synonymous) mutations to remove internal restriction sites while:
 * 1. Preserving the exact amino acid sequence
 * 2. Maintaining assembly efficiency (no additional fragments)
 * 3. Not creating new problematic restriction sites
 *
 * The junction-based approach is provided as a fallback for non-coding regions
 * or when silent mutations are not possible.
 *
 * KEY INSIGHT: In one-pot Golden Gate assembly, the enzyme is present throughout
 * the reaction. If the assembled product recreates an internal site, it will be
 * cut again, creating an equilibrium that reduces efficiency. Silent mutations
 * permanently remove the site from the sequence.
 */

import { GOLDEN_GATE_ENZYMES, findInternalSites, calculateExperimentalFidelity } from './goldengate.js';
import { reverseComplement } from './enzymes.js';

// ============================================================================
// TYPES
// ============================================================================

interface CodonTable {
  [aminoAcid: string]: string[];
}

interface CodonToAA {
  [codon: string]: string;
}

interface CodonUsageTable {
  [codon: string]: number;
}

interface DomesticationConfig {
  weights: {
    siteBreaking: number;
    noNewSites: number;
    codonFrequency: number;
    minimalChange: number;
    positionPreference: number;
  };
  minCodonUsage: number;
  rareCodoThreshold: number;
  allowJunctionFallback: boolean;
  warnOnJunctionFallback: boolean;
  newSiteDetectionWindow: number;
  requireFrameValidation: boolean;
  checkAllGoldenGateEnzymes: boolean;
}

interface InternalSite {
  position: number;
  sequence: string;
  orientation: 'forward' | 'reverse';
}

interface MutationCandidate {
  sequencePosition: number;
  positionInSite: number;
  positionInCodon: number;
  codonStart: number;
  originalBase: string;
  newBase: string;
  originalCodon: string;
  newCodon: string;
  aminoAcid: string;
  originalSiteSequence: string;
  newSiteSequence: string;
  siteOrientation: string;
  originalCodonFrequency: number;
  newCodonFrequency: number;
  frequencyRatio: number;
  isSynonymous: boolean;
  breaksSite: boolean;
  isRareCodon: boolean;
}

interface Penalty {
  reason: string;
  penalty: number;
  details: string | NewSitesResult;
}

interface Bonus {
  reason: string;
  bonus: number;
  details: string;
}

interface ScoredMutation extends MutationCandidate {
  score: number;
  penalties: Penalty[];
  bonuses: Bonus[];
  newSitesCreated: NewSitesResult;
}

interface NewSite {
  enzyme: string;
  orientation: string;
  recognition: string;
  addedCount: number;
}

interface NewSitesResult {
  createsNewSite: boolean;
  newSites: NewSite[];
  enzymesChecked: string[];
}

interface SiteMutation {
  site: InternalSite;
  mutation: ScoredMutation;
  alternatives: ScoredMutation[];
}

interface FailedSite {
  site: InternalSite;
  reason: string;
  message: string;
}

interface FallbackInfo {
  required: boolean;
  sites: FailedSite[];
  warning: string;
  message: string;
}

interface VerificationResult {
  originalSiteCount: number;
  remainingSiteCount: number;
  mutationsApplied: number;
  sitesRequiringFallback: number;
}

interface DomesticationResult {
  success: boolean;
  needsDomestication: boolean;
  method: string;
  originalSequence: string;
  domesticatedSequence: string;
  mutations: SiteMutation[];
  sites?: InternalSite[];
  failedSites?: FailedSite[];
  fallbackInfo?: FallbackInfo | null;
  verification?: VerificationResult;
  message: string;
}

interface DomesticationOptions {
  frame?: number;
  isCodingSequence?: boolean;
  organism?: 'ecoli' | 'yeast';
  allowJunctionFallback?: boolean;
  existingOverhangs?: string[];
  checkAllEnzymes?: boolean;
}

interface ProteinDifference {
  position: number;
  original: string;
  domesticated: string;
}

interface ProteinVerificationResult {
  identical: boolean;
  originalProtein: string;
  domesticatedProtein: string;
  originalLength: number;
  domesticatedLength: number;
  differences: ProteinDifference[];
  message: string;
}

interface ValidationCheck {
  check: string;
  passed: boolean;
  message: string;
  details?: any;
}

interface ValidationResult {
  isValid: boolean;
  validations: ValidationCheck[];
  summary: {
    totalChecks: number;
    passed: number;
    failed: number;
  };
}

interface ValidationOptions {
  frame?: number;
  checkAllEnzymes?: boolean;
}

interface DomesticationStrategy {
  name: string;
  available: boolean;
  fragmentIncrease: number;
  sequenceChange: boolean;
  proteinChange: boolean;
  onePotCompatible: boolean;
  mutations?: SiteMutation[];
  pros: string[];
  cons: string[];
}

interface StrategyComparison {
  needsDomestication: boolean;
  internalSites?: number;
  strategies?: {
    silentMutation: DomesticationStrategy;
    junctionBased: DomesticationStrategy;
  };
  recommended?: string;
  reason?: string;
  message?: string;
}

interface ComparisonOptions {
  frame?: number;
  existingOverhangs?: string[];
}

// ============================================================================
// CODON TABLES
// ============================================================================

/**
 * Standard genetic code: Amino acid -> list of codons
 */
export const CODON_TABLE: CodonTable = {
  F: ['TTT', 'TTC'],
  L: ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
  I: ['ATT', 'ATC', 'ATA'],
  M: ['ATG'],
  V: ['GTT', 'GTC', 'GTA', 'GTG'],
  S: ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
  P: ['CCT', 'CCC', 'CCA', 'CCG'],
  T: ['ACT', 'ACC', 'ACA', 'ACG'],
  A: ['GCT', 'GCC', 'GCA', 'GCG'],
  Y: ['TAT', 'TAC'],
  H: ['CAT', 'CAC'],
  Q: ['CAA', 'CAG'],
  N: ['AAT', 'AAC'],
  K: ['AAA', 'AAG'],
  D: ['GAT', 'GAC'],
  E: ['GAA', 'GAG'],
  C: ['TGT', 'TGC'],
  W: ['TGG'],
  R: ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
  G: ['GGT', 'GGC', 'GGA', 'GGG'],
  '*': ['TAA', 'TAG', 'TGA'],
};

/**
 * Reverse lookup: codon -> amino acid
 */
export const CODON_TO_AA: CodonToAA = {};
for (const [aa, codons] of Object.entries(CODON_TABLE)) {
  for (const codon of codons) {
    CODON_TO_AA[codon] = aa;
  }
}

/**
 * E. coli codon usage frequency (per 1000 codons)
 * Used for scoring mutations - prefer common codons
 * Source: Kazusa codon usage database
 */
export const ECOLI_CODON_USAGE: CodonUsageTable = {
  TTT: 22.0, TTC: 16.6, TTA: 13.8, TTG: 13.7,
  CTT: 11.0, CTC: 11.0, CTA: 3.9,  CTG: 52.6,
  ATT: 30.3, ATC: 25.1, ATA: 4.5,  ATG: 27.8,
  GTT: 18.4, GTC: 15.3, GTA: 10.9, GTG: 26.3,
  TCT: 8.6,  TCC: 8.8,  TCA: 7.2,  TCG: 8.9,
  CCT: 7.0,  CCC: 5.5,  CCA: 8.4,  CCG: 23.1,
  ACT: 9.0,  ACC: 23.4, ACA: 7.1,  ACG: 14.4,
  GCT: 15.5, GCC: 25.5, GCA: 20.3, GCG: 33.6,
  TAT: 16.2, TAC: 12.2, TAA: 2.0,  TAG: 0.2,
  CAT: 12.9, CAC: 9.7,  CAA: 15.3, CAG: 28.8,
  AAT: 17.7, AAC: 21.7, AAA: 33.6, AAG: 10.3,
  GAT: 32.1, GAC: 19.1, GAA: 39.6, GAG: 17.8,
  TGT: 5.2,  TGC: 6.5,  TGA: 1.0,  TGG: 15.2,
  CGT: 20.9, CGC: 22.0, CGA: 3.6,  CGG: 5.4,
  AGT: 8.8,  AGC: 16.0, AGA: 2.1,  AGG: 1.2,
  GGT: 24.7, GGC: 29.6, GGA: 7.9,  GGG: 11.0,
};

/**
 * Yeast (S. cerevisiae) codon usage frequency
 */
export const YEAST_CODON_USAGE: CodonUsageTable = {
  TTT: 26.1, TTC: 18.4, TTA: 26.2, TTG: 27.2,
  CTT: 12.3, CTC: 5.4,  CTA: 13.5, CTG: 10.5,
  ATT: 30.1, ATC: 17.2, ATA: 17.8, ATG: 20.9,
  GTT: 22.1, GTC: 11.8, GTA: 11.8, GTG: 10.8,
  TCT: 23.5, TCC: 14.2, TCA: 18.7, TCG: 8.6,
  CCT: 13.5, CCC: 6.8,  CCA: 18.3, CCG: 5.3,
  ACT: 20.3, ACC: 12.7, ACA: 17.8, ACG: 8.0,
  GCT: 21.2, GCC: 12.6, GCA: 16.2, GCG: 6.2,
  TAT: 18.8, TAC: 14.8, TAA: 1.1,  TAG: 0.5,
  CAT: 13.6, CAC: 7.8,  CAA: 27.3, CAG: 12.1,
  AAT: 35.7, AAC: 24.8, AAA: 41.9, AAG: 30.8,
  GAT: 37.6, GAC: 20.2, GAA: 45.6, GAG: 19.2,
  TGT: 8.1,  TGC: 4.8,  TGA: 0.7,  TGG: 10.4,
  CGT: 6.4,  CGC: 2.6,  CGA: 3.0,  CGG: 1.7,
  AGT: 14.2, AGC: 9.8,  AGA: 21.3, AGG: 9.2,
  GGT: 23.9, GGC: 9.8,  GGA: 10.9, GGG: 6.0,
};

// ============================================================================
// CONFIGURATION
// ============================================================================

export const DOMESTICATION_CONFIG: DomesticationConfig = {
  // Mutation selection weights
  weights: {
    siteBreaking: 1.0,        // Must break site (binary - required)
    noNewSites: 0.95,         // Avoid creating new restriction sites
    codonFrequency: 0.3,      // Prefer common codons (lower weight - not priority)
    minimalChange: 0.2,       // Prefer single nucleotide changes
    positionPreference: 0.1,  // Slight preference for middle positions
  },

  // Thresholds
  minCodonUsage: 5.0,         // Minimum acceptable codon frequency (per 1000)
  rareCodoThreshold: 10.0,    // Below this is considered "rare"

  // Fallback settings
  allowJunctionFallback: true,
  warnOnJunctionFallback: true,

  // Safety settings (NEW)
  newSiteDetectionWindow: 30, // Expanded from 20bp to catch edge cases
  requireFrameValidation: true,
  checkAllGoldenGateEnzymes: true,
};

// ============================================================================
// MAIN DOMESTICATION FUNCTIONS
// ============================================================================

/**
 * Main entry point: Domesticate a sequence using silent mutations
 */
export function domesticateWithSilentMutations(
  sequence: string,
  enzyme: string = 'BsaI',
  options: DomesticationOptions = {}
): DomesticationResult {
  const {
    frame = 0,                    // Reading frame (0, 1, or 2)
    isCodingSequence = true,      // Is this a coding sequence?
    organism = 'ecoli',           // For codon usage tables
    allowJunctionFallback = true, // Fall back to junction-based if needed
    existingOverhangs = [],       // Existing overhangs in assembly (for fidelity calc)
    checkAllEnzymes = false,      // Check for sites of all GG enzymes
  } = options;

  const seq = sequence.toUpperCase();
  const internalSites = findInternalSites(seq, enzyme);

  // No domestication needed
  if (!internalSites.hasSites) {
    return {
      success: true,
      needsDomestication: false,
      method: 'none',
      originalSequence: sequence,
      domesticatedSequence: sequence,
      mutations: [],
      sites: [],
      message: `No internal ${enzyme} sites found - sequence is already compatible`,
    };
  }

  const codonUsage = organism === 'yeast' ? YEAST_CODON_USAGE : ECOLI_CODON_USAGE;

  // For each internal site, find the best silent mutation
  const mutations: SiteMutation[] = [];
  const failedSites: FailedSite[] = [];
  let domesticatedSeq = seq;

  for (const site of internalSites.sites) {
    if (isCodingSequence) {
      // Find all possible silent mutations for this site
      const candidates = findAllSilentMutationCandidates(
        domesticatedSeq,
        site,
        frame,
        enzyme,
        codonUsage
      );

      // Score and rank candidates
      const scoredCandidates = scoreMutationCandidates(
        candidates,
        domesticatedSeq,
        enzyme,
        checkAllEnzymes
      );

      if (scoredCandidates.length > 0) {
        // Apply the best mutation
        const bestMutation = scoredCandidates[0];
        domesticatedSeq = applyMutation(domesticatedSeq, bestMutation);

        mutations.push({
          site,
          mutation: bestMutation,
          alternatives: scoredCandidates.slice(1, 5), // Keep top 5 alternatives
        });
      } else {
        failedSites.push({
          site,
          reason: 'NO_SYNONYMOUS_MUTATION',
          message: `No synonymous mutation found that breaks the ${enzyme} site at position ${site.position}`,
        });
      }
    } else {
      // Non-coding sequence - can't use silent mutations
      failedSites.push({
        site,
        reason: 'NON_CODING',
        message: `Site at position ${site.position} is in non-coding region - silent mutations not applicable`,
      });
    }
  }

  // Verify domestication success
  const verifyResult = findInternalSites(domesticatedSeq, enzyme);

  // Handle failed sites with junction-based fallback
  let fallbackInfo: FallbackInfo | null = null;
  if (failedSites.length > 0 && allowJunctionFallback) {
    fallbackInfo = {
      required: true,
      sites: failedSites,
      warning: 'JUNCTION_FALLBACK_REQUIRED',
      message: `${failedSites.length} site(s) require junction-based domestication. ` +
               'Note: Junction-based domestication adds fragments and recreates sites ' +
               'in one-pot reactions. Consider sequential assembly (digest → heat-kill → ligate) ' +
               'or using a different enzyme.',
    };
  }

  return {
    success: verifyResult.count === 0 || (failedSites.length > 0 && allowJunctionFallback),
    needsDomestication: true,
    method: failedSites.length === 0 ? 'silent_mutation' : 'hybrid',
    originalSequence: sequence,
    domesticatedSequence: domesticatedSeq,
    mutations,
    failedSites,
    fallbackInfo,
    verification: {
      originalSiteCount: internalSites.count,
      remainingSiteCount: verifyResult.count,
      mutationsApplied: mutations.length,
      sitesRequiringFallback: failedSites.length,
    },
    message: generateSummaryMessage(mutations, failedSites, enzyme),
  };
}

/**
 * Find ALL possible silent mutations that could break a restriction site
 *
 * For a 6bp recognition site like GGTCTC, each position overlaps with 1-2 codons.
 * We enumerate all possible single-nucleotide changes that:
 * 1. Break the recognition sequence
 * 2. Preserve the amino acid (synonymous)
 */
export function findAllSilentMutationCandidates(
  sequence: string,
  site: InternalSite,
  frame: number,
  enzyme: string,
  codonUsage: CodonUsageTable
): MutationCandidate[] {
  const candidates: MutationCandidate[] = [];
  const recognition = site.sequence;
  const siteStart = site.position;
  const enz = GOLDEN_GATE_ENZYMES[enzyme];

  if (!enz) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  const enzRecognition = enz.recognition;
  const enzRecognitionRC = reverseComplement(enzRecognition);

  // For each position in the recognition sequence
  for (let posInSite = 0; posInSite < recognition.length; posInSite++) {
    const seqPos = siteStart + posInSite;
    const originalBase = recognition[posInSite];

    // Determine which codon this position is in
    const adjustedPos = seqPos - frame;
    if (adjustedPos < 0) continue;

    const codonIndex = Math.floor(adjustedPos / 3);
    const codonStart = codonIndex * 3 + frame;
    const posInCodon = (seqPos - codonStart);

    // Validate codon boundaries
    if (codonStart < 0 || codonStart + 3 > sequence.length) continue;
    if (posInCodon < 0 || posInCodon > 2) continue;

    const originalCodon = sequence.slice(codonStart, codonStart + 3);
    const originalAA = CODON_TO_AA[originalCodon];

    if (!originalAA) continue; // Invalid codon

    // Try each alternative base
    for (const newBase of ['A', 'T', 'G', 'C']) {
      if (newBase === originalBase) continue;

      // Create new codon with this mutation
      const newCodon =
        originalCodon.slice(0, posInCodon) +
        newBase +
        originalCodon.slice(posInCodon + 1);

      const newAA = CODON_TO_AA[newCodon];

      // Must be synonymous (same amino acid)
      if (newAA !== originalAA) continue;

      // Check if this mutation breaks the recognition site
      const newSiteSeq =
        recognition.slice(0, posInSite) +
        newBase +
        recognition.slice(posInSite + 1);

      // Site is broken if it no longer matches recognition or its RC
      const breaksForward = newSiteSeq !== enzRecognition;
      const breaksReverse = newSiteSeq !== enzRecognitionRC;
      const siteIsBroken = site.orientation === 'forward' ? breaksForward : breaksReverse;

      if (!siteIsBroken) continue;

      // This is a valid candidate!
      const codonFrequency = codonUsage[newCodon] || 0;
      const originalFrequency = codonUsage[originalCodon] || 0;

      candidates.push({
        // Position information
        sequencePosition: seqPos,
        positionInSite: posInSite,
        positionInCodon: posInCodon,
        codonStart: codonStart,

        // Mutation details
        originalBase,
        newBase,
        originalCodon,
        newCodon,
        aminoAcid: originalAA,

        // Site breaking info
        originalSiteSequence: recognition,
        newSiteSequence: newSiteSeq,
        siteOrientation: site.orientation,

        // Codon usage info
        originalCodonFrequency: originalFrequency,
        newCodonFrequency: codonFrequency,
        frequencyRatio: originalFrequency > 0 ? codonFrequency / originalFrequency : 0,

        // Flags
        isSynonymous: true,
        breaksSite: true,
        isRareCodon: codonFrequency < DOMESTICATION_CONFIG.rareCodoThreshold,
      });
    }
  }

  return candidates;
}

/**
 * Score and rank mutation candidates
 *
 * Scoring priorities (in order):
 * 1. Must break the site (all candidates do, so this is filtered earlier)
 * 2. Must not create new restriction sites (critical)
 * 3. Codon frequency (prefer common codons, but not priority)
 * 4. Position preference (middle of site slightly better)
 */
export function scoreMutationCandidates(
  candidates: MutationCandidate[],
  sequence: string,
  enzyme: string,
  checkAllEnzymes: boolean = false
): ScoredMutation[] {
  const weights = DOMESTICATION_CONFIG.weights;
  const scoredCandidates: ScoredMutation[] = [];

  for (const candidate of candidates) {
    let score = 100; // Start with perfect score
    const penalties: Penalty[] = [];
    const bonuses: Bonus[] = [];

    // Apply the mutation temporarily to check for new sites
    const mutatedSeq = applyMutation(sequence, candidate);

    // CRITICAL: Check if mutation creates new restriction sites
    const newSitesCreated = checkForNewSites(
      sequence,
      mutatedSeq,
      candidate.sequencePosition,
      enzyme,
      checkAllEnzymes
    );

    if (newSitesCreated.createsNewSite) {
      score -= 80 * weights.noNewSites;
      penalties.push({
        reason: 'CREATES_NEW_SITE',
        penalty: 80,
        details: newSitesCreated,
      });
    }

    // Codon frequency scoring (lower priority)
    if (candidate.isRareCodon) {
      const penaltyAmount = 20 * weights.codonFrequency;
      score -= penaltyAmount;
      penalties.push({
        reason: 'RARE_CODON',
        penalty: penaltyAmount,
        details: `${candidate.newCodon} frequency: ${candidate.newCodonFrequency.toFixed(1)}/1000`,
      });
    } else if (candidate.frequencyRatio >= 0.8) {
      // Bonus for maintaining similar codon frequency
      const bonusAmount = 10 * weights.codonFrequency;
      score += bonusAmount;
      bonuses.push({
        reason: 'GOOD_CODON_FREQUENCY',
        bonus: bonusAmount,
        details: `${candidate.newCodon} frequency: ${candidate.newCodonFrequency.toFixed(1)}/1000`,
      });
    }

    // Position preference: middle of recognition site is slightly better
    // (more robust to nearby sequence variations)
    const recognitionLength = candidate.originalSiteSequence.length;
    const distFromMiddle = Math.abs(candidate.positionInSite - (recognitionLength - 1) / 2);
    const positionScore = (1 - distFromMiddle / (recognitionLength / 2)) * 5 * weights.positionPreference;
    score += positionScore;

    // Wobble position preference (3rd position mutations often safer)
    if (candidate.positionInCodon === 2) {
      score += 5 * weights.minimalChange;
      bonuses.push({
        reason: 'WOBBLE_POSITION',
        bonus: 5,
        details: 'Mutation at wobble position (3rd codon position)',
      });
    }

    scoredCandidates.push({
      ...candidate,
      score: Math.max(0, Math.min(100, score)),
      penalties,
      bonuses,
      newSitesCreated,
    });
  }

  // Sort by score (highest first)
  scoredCandidates.sort((a, b) => b.score - a.score);

  return scoredCandidates;
}

/**
 * Check if a mutation creates new restriction sites
 */
export function checkForNewSites(
  originalSeq: string,
  mutatedSeq: string,
  mutationPos: number,
  currentEnzyme: string,
  checkAll: boolean = false
): NewSitesResult {
  const enzymesToCheck = checkAll
    ? Object.keys(GOLDEN_GATE_ENZYMES)
    : [currentEnzyme];

  const newSites: NewSite[] = [];
  // EXPANDED: Use configurable window size (default 30bp, was 20bp)
  // This catches edge cases where new sites could be created just outside the old window
  const windowSize = DOMESTICATION_CONFIG.newSiteDetectionWindow || 30;

  const windowStart = Math.max(0, mutationPos - windowSize);
  const windowEnd = Math.min(mutatedSeq.length, mutationPos + windowSize);

  for (const enzName of enzymesToCheck) {
    const enz = GOLDEN_GATE_ENZYMES[enzName];
    if (!enz) continue;

    const recognition = enz.recognition;
    const recognitionRC = reverseComplement(recognition);

    // Check the window region for new sites
    const originalWindow = originalSeq.slice(windowStart, windowEnd);
    const mutatedWindow = mutatedSeq.slice(windowStart, windowEnd);

    // Count sites in original vs mutated
    const originalForward = countOccurrences(originalWindow, recognition);
    const mutatedForward = countOccurrences(mutatedWindow, recognition);
    const originalReverse = countOccurrences(originalWindow, recognitionRC);
    const mutatedReverse = countOccurrences(mutatedWindow, recognitionRC);

    if (mutatedForward > originalForward) {
      newSites.push({
        enzyme: enzName,
        orientation: 'forward',
        recognition,
        addedCount: mutatedForward - originalForward,
      });
    }

    if (mutatedReverse > originalReverse) {
      newSites.push({
        enzyme: enzName,
        orientation: 'reverse',
        recognition: recognitionRC,
        addedCount: mutatedReverse - originalReverse,
      });
    }
  }

  return {
    createsNewSite: newSites.length > 0,
    newSites,
    enzymesChecked: enzymesToCheck,
  };
}

/**
 * Apply a mutation to a sequence
 */
export function applyMutation(sequence: string, mutation: MutationCandidate | ScoredMutation): string {
  const pos = mutation.sequencePosition;
  return sequence.slice(0, pos) + mutation.newBase + sequence.slice(pos + 1);
}

/**
 * Apply multiple mutations to a sequence
 */
export function applyMutations(sequence: string, mutations: SiteMutation[]): string {
  // Sort by position descending to avoid position shifting issues
  const sorted = [...mutations].sort((a, b) => b.mutation.sequencePosition - a.mutation.sequencePosition);

  let result = sequence;
  for (const mutation of sorted) {
    result = applyMutation(result, mutation.mutation);
  }
  return result;
}

// ============================================================================
// VALIDATION & VERIFICATION
// ============================================================================

/**
 * Verify that domesticated sequence produces identical protein
 */
export function verifyProteinSequence(
  originalSeq: string,
  domesticatedSeq: string,
  frame: number = 0
): ProteinVerificationResult {
  const translateSequence = (seq: string, fr: number): string => {
    const protein: string[] = [];
    for (let i = fr; i + 3 <= seq.length; i += 3) {
      const codon = seq.slice(i, i + 3).toUpperCase();
      const aa = CODON_TO_AA[codon] || '?';
      protein.push(aa);
      if (aa === '*') break; // Stop at stop codon
    }
    return protein.join('');
  };

  const originalProtein = translateSequence(originalSeq, frame);
  const domesticatedProtein = translateSequence(domesticatedSeq, frame);

  const identical = originalProtein === domesticatedProtein;

  // Find differences if any
  const differences: ProteinDifference[] = [];
  if (!identical) {
    for (let i = 0; i < Math.max(originalProtein.length, domesticatedProtein.length); i++) {
      if (originalProtein[i] !== domesticatedProtein[i]) {
        differences.push({
          position: i,
          original: originalProtein[i] || '-',
          domesticated: domesticatedProtein[i] || '-',
        });
      }
    }
  }

  return {
    identical,
    originalProtein,
    domesticatedProtein,
    originalLength: originalProtein.length,
    domesticatedLength: domesticatedProtein.length,
    differences,
    message: identical
      ? 'Protein sequences are identical'
      : `Protein sequences differ at ${differences.length} position(s)`,
  };
}

/**
 * Comprehensive validation of domestication result
 */
export function validateDomestication(
  domesticationResult: DomesticationResult,
  enzyme: string = 'BsaI',
  options: ValidationOptions = {}
): ValidationResult {
  const {
    frame = 0,
    checkAllEnzymes = false,
  } = options;

  const validations: ValidationCheck[] = [];
  let isValid = true;

  // 1. Verify protein sequence preserved
  const proteinCheck = verifyProteinSequence(
    domesticationResult.originalSequence,
    domesticationResult.domesticatedSequence,
    frame
  );

  validations.push({
    check: 'PROTEIN_SEQUENCE',
    passed: proteinCheck.identical,
    message: proteinCheck.message,
    details: proteinCheck,
  });

  if (!proteinCheck.identical) {
    isValid = false;
  }

  // 2. Verify no internal sites remain
  const siteCheck = findInternalSites(domesticationResult.domesticatedSequence, enzyme);

  validations.push({
    check: 'NO_INTERNAL_SITES',
    passed: !siteCheck.hasSites,
    message: siteCheck.hasSites
      ? `${siteCheck.count} internal ${enzyme} site(s) remain`
      : `No internal ${enzyme} sites in domesticated sequence`,
    details: siteCheck,
  });

  if (siteCheck.hasSites && domesticationResult.failedSites && domesticationResult.failedSites.length === 0) {
    isValid = false;
  }

  // 3. Check for new sites created by mutations
  if (checkAllEnzymes) {
    for (const enzName of Object.keys(GOLDEN_GATE_ENZYMES)) {
      if (enzName === enzyme) continue;

      const originalSites = findInternalSites(domesticationResult.originalSequence, enzName);
      const newSites = findInternalSites(domesticationResult.domesticatedSequence, enzName);

      if (newSites.count > originalSites.count) {
        validations.push({
          check: `NEW_SITES_${enzName}`,
          passed: false,
          message: `Domestication created ${newSites.count - originalSites.count} new ${enzName} site(s)`,
          details: { originalCount: originalSites.count, newCount: newSites.count },
        });
        // This is a warning, not a failure
      }
    }
  }

  // 4. Verify sequence length unchanged
  const lengthCheck = domesticationResult.originalSequence.length ===
                      domesticationResult.domesticatedSequence.length;

  validations.push({
    check: 'SEQUENCE_LENGTH',
    passed: lengthCheck,
    message: lengthCheck
      ? 'Sequence length unchanged'
      : `Sequence length changed from ${domesticationResult.originalSequence.length} to ${domesticationResult.domesticatedSequence.length}`,
  });

  if (!lengthCheck) {
    isValid = false;
  }

  return {
    isValid,
    validations,
    summary: {
      totalChecks: validations.length,
      passed: validations.filter(v => v.passed).length,
      failed: validations.filter(v => !v.passed).length,
    },
  };
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Count occurrences of a substring in a string
 */
function countOccurrences(str: string, substr: string): number {
  let count = 0;
  let pos = str.indexOf(substr);
  while (pos !== -1) {
    count++;
    pos = str.indexOf(substr, pos + 1);
  }
  return count;
}

/**
 * Generate summary message for domestication result
 */
function generateSummaryMessage(mutations: SiteMutation[], failedSites: FailedSite[], enzyme: string): string {
  const parts: string[] = [];

  if (mutations.length > 0) {
    parts.push(`Applied ${mutations.length} silent mutation(s) to remove internal ${enzyme} sites`);
  }

  if (failedSites.length > 0) {
    parts.push(`${failedSites.length} site(s) could not be removed with silent mutations`);
  }

  if (mutations.length === 0 && failedSites.length === 0) {
    return 'No domestication performed';
  }

  return parts.join('. ') + '.';
}

// ============================================================================
// COMPARISON & ANALYSIS
// ============================================================================

/**
 * Compare domestication strategies for a given sequence
 */
export function compareDomesticationStrategies(
  sequence: string,
  enzyme: string = 'BsaI',
  options: ComparisonOptions = {}
): StrategyComparison {
  const { frame = 0, existingOverhangs = [] } = options;

  const internalSites = findInternalSites(sequence, enzyme);

  if (!internalSites.hasSites) {
    return {
      needsDomestication: false,
      message: 'No internal sites found - no domestication needed',
    };
  }

  // Strategy 1: Silent mutations
  const silentResult = domesticateWithSilentMutations(sequence, enzyme, {
    frame,
    isCodingSequence: true,
    allowJunctionFallback: false,
  });

  // Strategy 2: Junction-based (for comparison)
  // Import from auto-domestication-optimizer would be needed here
  // For now, we'll estimate the impact

  const strategies = {
    silentMutation: {
      name: 'Silent Mutation',
      available: (silentResult.failedSites?.length || 0) === 0,
      fragmentIncrease: 0,
      sequenceChange: true,
      proteinChange: false,
      onePotCompatible: true,
      mutations: silentResult.mutations,
      pros: [
        'No additional fragments needed',
        'Works in one-pot reactions',
        'Maintains assembly efficiency',
      ],
      cons: [
        'Slightly changes DNA sequence',
        'May affect codon usage',
        'Only works for coding sequences',
      ],
    },
    junctionBased: {
      name: 'Junction-Based',
      available: true,
      fragmentIncrease: internalSites.count,
      sequenceChange: false,
      proteinChange: false,
      onePotCompatible: false, // THE KEY ISSUE
      pros: [
        'Preserves exact DNA sequence',
        'Works for non-coding regions',
      ],
      cons: [
        'Adds additional fragments (reduced efficiency)',
        'NOT compatible with one-pot reactions',
        'Requires sequential assembly protocol',
        'Higher cost (more primers)',
      ],
    },
  };

  // Recommendation
  let recommended = 'silentMutation';
  let reason = 'Silent mutations maintain assembly efficiency and work in one-pot reactions';

  if (!strategies.silentMutation.available) {
    recommended = 'junctionBased';
    reason = 'Silent mutations not available for all sites - junction-based required';
  }

  return {
    needsDomestication: true,
    internalSites: internalSites.count,
    strategies,
    recommended,
    reason,
  };
}

// ============================================================================
// EXPORTS
// ============================================================================

export {
  domesticateWithSilentMutations as domesticate,
};
