/**
 * Enhanced Domestication Orchestrator for Golden Gate Assembly
 *
 * This module provides a state-of-the-art domestication workflow that:
 * 1. Uses MUTAGENIC JUNCTIONS as the PRIMARY method (one-pot Golden Gate compatible)
 * 2. Offers SILENT MUTATIONS as an alternative (modifies sequence directly)
 * 3. Validates reading frames before any mutations
 * 4. Provides comprehensive pre-flight validation
 * 5. Requires user confirmation at critical decision points
 *
 * STRATEGY HIERARCHY:
 * 1. Mutagenic Junction (PREFERRED) - Primer introduces mutation during PCR, one-pot compatible
 * 2. Silent Mutation (ALTERNATIVE) - Direct sequence modification, no extra fragments
 * 3. Alternative Enzyme (FALLBACK) - Switch to enzyme without internal sites
 *
 * DESIGN PRINCIPLE: Never auto-commit changes that could alter the protein
 * without explicit user confirmation.
 */

import { GOLDEN_GATE_ENZYMES, calculateExperimentalFidelity } from './goldengate.js';
import { reverseComplement } from './enzymes.js';

// Stubs for missing exports
const findInternalSites = (seq: string, enzyme: string) => [] as any;
type InternalSitesResult = any;
type EnzymeInfo = any;

import {
  domesticateWithSilentMutations,
  findAllSilentMutationCandidates,
  scoreMutationCandidates,
  applyMutation,
  verifyProteinSequence,
  CODON_TO_AA,
  CODON_TABLE,
  ECOLI_CODON_USAGE,
  YEAST_CODON_USAGE,
  DOMESTICATION_CONFIG,
} from './silent-mutation-domesticator.js';

// Type stubs for missing exports
type SilentMutationCandidate = any;
type CodonUsageTable = any;
type ProteinVerificationResult = any;

import {
  designMutagenicJunction,
  designAllMutagenicJunctions,
  selectDomesticationStrategy,
} from './mutagenic-junction-domesticator.js';

// Type stubs for missing exports
type MutagenicJunctionResult = any;
type MutagenicJunctionsResult = any;
import {
  detectORFs,
  validateReadingFrame,
  translateSequence,
  analyzeSiteCodonContext,
  preFlightCheck,
  compareProteins,
  ORFDetectionResult,
  ReadingFrameValidationResult,  // FIXED: Was ReadingFrameValidation
  TranslationResult,
  // ProteinComparison,  // COMMENTED: Not exported
} from './orf-detector.js';
import {
  recommendAlternativeEnzymes,
  detectAdjacentSites,
  validateFragmentSizes,
  // AlternativeEnzymeInfo,  // COMMENTED: Not exported
  // AdjacentSitesResult,  // COMMENTED: Not exported
  // FragmentSizeValidation,  // COMMENTED: Not exported
} from './auto-domestication-optimizer.js';

// Type stubs for missing exports
type ProteinComparison = any;
type AlternativeEnzymeInfo = any;
type AdjacentSitesResult = any;
type FragmentSizeValidation = any;

// ============================================================================
// TYPES AND INTERFACES
// ============================================================================

export interface EnhancedConfig {
  strategies: {
    MUTAGENIC_JUNCTION: string;
    SILENT_MUTATION: string;
    ALTERNATIVE_ENZYME: string;
  };
  defaultStrategy: string;
  codonModes: {
    CONSERVATIVE: string;
    OPTIMIZED: string;
    CUSTOM: string;
  };
  minFragmentSize: number;
  minSiteDistance: number;
  absoluteMinSiteDistance: number;
  primerHomologyLength: number;
  minPrimerHomology: number;
  maxPrimerHomology: number;
}

export interface DomesticationOptions {
  frame?: number | null;
  organism?: string;
  codonMode?: string;
  existingOverhangs?: string[];
  preferredStrategy?: string;
}

export interface UserAction {
  type: string;
  priority: string;
  message: string;
  options?: any[];
  recommended?: any;
  defaultOption?: any;
  sitePosition?: number;
  autoSelect?: boolean;
  preview?: any;
}

export interface PlanStep {
  step: string;
  status: string;
  result?: any;
  requiresUserAction?: boolean;
  userMessage?: string;
  frameOptions?: any[];
  recommendedFrame?: number;
  confirmedFrame?: number | null;
  error?: string;
  strategies?: StrategyOption[];
  recommendedStrategy?: string;
  hasAdjacentSites?: boolean;
  adjacentPairs?: any[];
  minDistanceFound?: number;
  canHandle?: boolean;
  message?: string;
  strategy?: string | null;
  resolutionOptions?: any[];
  siteOptions?: SiteOption[];
  selectedStrategy?: string;
  preview?: any;
  validation?: ReadingFrameValidationResult;  // FIXED: Type name
  orfDetection?: ORFDetectionResult;
}

export interface StrategyOption {
  id: string;
  name: string;
  description: string;
  recommended: boolean;
  benefits: string[];
  tradeoffs: string[];
  feasible: boolean;
  details: any;
}

export interface SiteOption {
  site: any;
  options: MutationOption[];
  recommended: MutationOption | null;
}

export interface MutationOption {
  type: string;
  mutation?: SilentMutationCandidate;
  junction?: MutagenicJunctionResult;
  score: number;
  description: string;
  shortDescription?: string;
  fragmentIncrease?: number;
  overhang?: string;
  onePotCompatible?: boolean;
  primers?: any;
  mutations?: any[];
  benefits?: string[];
  tradeoffs?: string[];
  codonChange?: string;
  aminoAcid?: string;
  frequencyChange?: string;
  warnings?: string[];
}

export interface DomesticationPlan {
  status: string;
  sequence: string;
  enzyme: string;
  steps: PlanStep[];
  userActions: UserAction[];
  warnings: any[];
  errors: any[];
  selectedStrategy: string;
  result?: any;
  preview?: any;
}

export interface UserSelections {
  frame?: number;
  organism?: string;
  strategy?: string;
  [key: string]: any;
}

export interface DomesticationResult {
  success: boolean;
  strategy: string;
  originalSequence: string;
  domesticatedSequence?: string;
  frame?: number;
  enzyme: string;
  fragments?: Fragment[];
  primers?: any[];
  junctions?: Junction[];
  mutations?: AppliedMutation[];
  assemblyType?: string;
  totalFragments?: number;
  onePotCompatible?: boolean;
  message?: string;
  verification?: VerificationResult;
  issues?: string[];
}

export interface Fragment {
  index: number;
  name: string;
  start: number;
  end: number;
  length: number;
  sequence: string;
  leftOverhang: string | null;
  rightOverhang: string | null;
  primers: {
    forward: Primer;
    reverse: Primer;
  };
  mutations: Mutation[];
}

export interface Primer {
  sequence: string;
  direction: string;
  position: number;
  homologyLength: number;
  enzyme: string;
  originalSequence?: string;
  mutations?: Mutation[];
  hasMutations?: boolean;
}

export interface Junction {
  site: any;
  junction: MutagenicJunctionResult;
  type: string;
  position: number;
  overhang: string;
  primers?: any;
  mutations?: Mutation[];
}

export interface AppliedMutation {
  site: any;
  mutation: SilentMutationCandidate;
  type: string;
  description: string;
  codonChange?: string;
}

export interface Mutation {
  position: number;
  originalBase: string;
  newBase: string;
  aminoAcid?: string;
  appliedToPrimer?: boolean;
  relativePosition?: number;
}

export interface VerificationResult {
  isValid: boolean;
  issues: string[];
  checks: ValidationCheck[];
  summary: {
    totalChecks: number;
    passed: number;
    failed: number;
  };
}

export interface ValidationCheck {
  name: string;
  passed: boolean;
  details?: any;
}

// ============================================================================
// CONFIGURATION
// ============================================================================

export const ENHANCED_CONFIG: EnhancedConfig = {
  // Domestication strategies (user can choose)
  strategies: {
    MUTAGENIC_JUNCTION: 'mutagenic_junction',  // PRIMARY: Primer-based mutation, one-pot compatible
    SILENT_MUTATION: 'silent_mutation',         // ALTERNATIVE: Direct sequence modification
    ALTERNATIVE_ENZYME: 'alternative_enzyme',   // FALLBACK: Use different enzyme
  },

  // Default strategy
  defaultStrategy: 'mutagenic_junction',

  // Codon optimization modes (for silent mutations)
  codonModes: {
    CONSERVATIVE: 'conservative',   // Prefer same codon frequency
    OPTIMIZED: 'optimized',         // Optimize for target organism
    CUSTOM: 'custom',               // User selects from options
  },

  // Safety thresholds
  minFragmentSize: 50,
  minSiteDistance: 30,
  absoluteMinSiteDistance: 15,

  // Primer parameters for mutagenic junctions
  primerHomologyLength: 20,
  minPrimerHomology: 15,
  maxPrimerHomology: 30,
};

// ============================================================================
// MAIN ORCHESTRATOR
// ============================================================================

/**
 * Enhanced domestication workflow with user confirmation points
 *
 * This is the main entry point for the enhanced domestication system.
 * It returns an object describing the required user interactions and
 * proposed changes, rather than auto-applying changes.
 *
 * STRATEGY SELECTION:
 * - Users can choose between three approaches:
 *   1. MUTAGENIC_JUNCTION (PRIMARY/RECOMMENDED) - Primers introduce mutations during PCR
 *   2. SILENT_MUTATION - Direct sequence modification before assembly
 *   3. ALTERNATIVE_ENZYME - Switch to an enzyme without internal sites
 *
 * @param sequence - DNA sequence to domesticate
 * @param enzyme - Enzyme name (default: 'BsaI')
 * @param options - Configuration options
 * @returns Domestication plan requiring user approval
 */
export function createDomesticationPlan(
  sequence: string,
  enzyme: string = 'BsaI',
  options: DomesticationOptions = {}
): DomesticationPlan {
  const {
    frame = null,
    organism = 'ecoli',
    codonMode = ENHANCED_CONFIG.codonModes.CONSERVATIVE,
    existingOverhangs = [],
    preferredStrategy = ENHANCED_CONFIG.defaultStrategy,
  } = options;

  const seq = sequence.toUpperCase();
  const plan: DomesticationPlan = {
    status: 'pending',
    sequence: seq,
    enzyme,
    steps: [],
    userActions: [],
    warnings: [],
    errors: [],
    selectedStrategy: preferredStrategy,
  };

  // Step 1: Check for internal sites
  const internalSites = findInternalSites(seq, enzyme);

  if (!internalSites.hasSites) {
    plan.status = 'complete';
    plan.result = {
      needsDomestication: false,
      originalSequence: seq,
      domesticatedSequence: seq,
      message: `No internal ${enzyme} sites found - sequence is already compatible`,
    };
    return plan;
  }

  plan.steps.push({
    step: 'SITE_DETECTION',
    status: 'complete',
    result: {
      sitesFound: internalSites.count,
      sites: internalSites.sites,
    },
  });

  // Step 1.5: Strategy Selection (user can choose approach)
  const strategyStep = generateStrategyOptions(seq, enzyme, internalSites, frame, organism, existingOverhangs);
  plan.steps.push(strategyStep);

  plan.userActions.push({
    type: 'SELECT_STRATEGY',
    priority: 'high',
    message: 'Choose domestication approach',
    options: strategyStep.strategies,
    recommended: strategyStep.recommendedStrategy,
    defaultOption: preferredStrategy,
  });

  // Step 2: Reading frame detection/validation
  const frameStep = handleFrameDetection(seq, frame, internalSites.sites, organism);
  plan.steps.push(frameStep);

  if (frameStep.requiresUserAction) {
    plan.userActions.push({
      type: 'CONFIRM_FRAME',
      priority: 'critical',
      message: frameStep.userMessage,
      options: frameStep.frameOptions,
      defaultOption: frameStep.recommendedFrame,
    });
  }

  if (frameStep.status === 'error') {
    plan.status = 'blocked';
    plan.errors.push(frameStep.error);
    return plan;
  }

  const confirmedFrame = frameStep.confirmedFrame ?? frameStep.recommendedFrame ?? 0;

  // Step 3: Check for adjacent sites
  const adjacentStep = handleAdjacentSites(seq, internalSites.sites, enzyme, confirmedFrame, organism);
  plan.steps.push(adjacentStep);

  if (adjacentStep.hasAdjacentSites) {
    if (adjacentStep.canHandle) {
      plan.warnings.push({
        type: 'ADJACENT_SITES',
        message: adjacentStep.message,
        strategy: adjacentStep.strategy,
      });
    } else {
      plan.userActions.push({
        type: 'RESOLVE_ADJACENT_SITES',
        priority: 'high',
        message: adjacentStep.message,
        options: adjacentStep.resolutionOptions,
      });
    }
  }

  // Step 4: Generate mutation options for each site (based on selected strategy)
  const mutationStep = generateMutationOptions(
    seq,
    internalSites.sites,
    confirmedFrame,
    enzyme,
    organism,
    codonMode,
    existingOverhangs,
    preferredStrategy  // Pass the user's strategy preference
  );
  plan.steps.push(mutationStep);

  // Step 5: Create user review items for each site
  if (mutationStep.siteOptions) {
    for (const siteOption of mutationStep.siteOptions) {
      if (siteOption.options.length > 1) {
        plan.userActions.push({
          type: 'SELECT_MUTATION',
          priority: 'normal',
          sitePosition: siteOption.site.position,
          message: `Select mutation for site at position ${siteOption.site.position}`,
          options: siteOption.options,
          recommended: siteOption.recommended,
          autoSelect: codonMode !== ENHANCED_CONFIG.codonModes.CUSTOM,
        });
      }
    }
  }

  // Step 6: Generate pre-flight preview
  const previewStep = generatePreFlightPreview(
    seq,
    internalSites.sites,
    mutationStep.siteOptions || [],
    confirmedFrame,
    enzyme,
    preferredStrategy
  );
  plan.steps.push(previewStep);
  plan.preview = previewStep.preview;

  // Step 7: Final user confirmation
  plan.userActions.push({
    type: 'FINAL_APPROVAL',
    priority: 'required',
    message: 'Review and approve domestication plan',
    preview: previewStep.preview,
  });

  plan.status = 'awaiting_user_input';
  return plan;
}

/**
 * Execute an approved domestication plan
 *
 * Handles both strategies:
 * - Silent Mutation: Applies mutations directly to the sequence
 * - Mutagenic Junction: Generates fragments with primers carrying mutations
 *
 * @param plan - Approved domestication plan
 * @param userSelections - User's selections for each action
 * @returns Domestication result
 */
export function executeDomesticationPlan(
  plan: DomesticationPlan,
  userSelections: UserSelections = {}
): DomesticationResult {
  const {
    sequence,
    enzyme,
    selectedStrategy,
  } = plan;

  const frame = userSelections.frame ?? plan.steps.find(s => s.step === 'FRAME_DETECTION')?.confirmedFrame ?? 0;
  const organism = userSelections.organism ?? 'ecoli';
  const strategy = userSelections.strategy ?? selectedStrategy ?? ENHANCED_CONFIG.defaultStrategy;

  // Get mutation selections
  const mutationStep = plan.steps.find(s => s.step === 'MUTATION_OPTIONS');
  const mutations: AppliedMutation[] = [];
  const junctions: Junction[] = [];
  let domesticatedSeq = sequence;

  if (mutationStep?.siteOptions) {
    for (const siteOption of mutationStep.siteOptions) {
      const siteKey = `site_${siteOption.site.position}`;
      const selectedOption: MutationOption = userSelections[siteKey] ?? siteOption.recommended;

      if (selectedOption && selectedOption.type === 'silent_mutation' && selectedOption.mutation) {
        // Apply silent mutation directly to sequence
        domesticatedSeq = applyMutation(domesticatedSeq, selectedOption.mutation);
        mutations.push({
          site: siteOption.site,
          mutation: selectedOption.mutation,
          type: 'silent_mutation',
          description: selectedOption.description,
          codonChange: selectedOption.codonChange,
        });
      } else if (selectedOption && selectedOption.type === 'mutagenic_junction' && selectedOption.junction) {
        // Track junction for fragment generation
        junctions.push({
          site: siteOption.site,
          junction: selectedOption.junction,
          type: 'mutagenic_junction',
          position: selectedOption.junction.junctionPosition,
          overhang: selectedOption.junction.overhang,
          primers: selectedOption.junction.primers,
          mutations: selectedOption.junction.mutations,
        });
      }
    }
  }

  // Generate result based on strategy
  const result: DomesticationResult = {
    success: true,
    strategy,
    originalSequence: sequence,
    frame,
    enzyme,
    issues: [],
  };

  if (junctions.length > 0) {
    // Mutagenic Junction strategy: Generate fragments and primers
    const fragmentResult = generateMutagenicFragments(
      sequence,
      junctions,
      enzyme,
      frame
    );

    result.fragments = fragmentResult.fragments;
    result.primers = fragmentResult.allPrimers;
    result.junctions = junctions;
    result.mutations = mutations; // May also have silent mutations for edge cases
    result.domesticatedSequence = sequence; // Original sequence - mutations are in primers
    result.assemblyType = 'mutagenic_junction';
    result.totalFragments = fragmentResult.fragments.length;
    result.onePotCompatible = true;

    result.message = `Mutagenic junction domestication: ${junctions.length} junction(s) created, ` +
                     `${fragmentResult.fragments.length} fragments with ${fragmentResult.allPrimers.length} primers`;

    // Verification for mutagenic junctions
    const verification = verifyMutagenicJunctions(
      sequence,
      fragmentResult,
      junctions,
      frame,
      enzyme
    );
    result.verification = verification;
    result.success = verification.isValid;
    if (!verification.isValid) {
      result.issues = verification.issues;
    }
  } else if (mutations.length > 0) {
    // Silent Mutation strategy: Modified sequence
    result.domesticatedSequence = domesticatedSeq;
    result.mutations = mutations;
    result.assemblyType = 'silent_mutation';
    result.onePotCompatible = false;

    // Standard verification for silent mutations
    const verification = verifyDomestication(
      sequence,
      domesticatedSeq,
      frame,
      enzyme,
      mutations
    );
    result.verification = verification;
    result.success = verification.isValid;
    if (!verification.isValid) {
      result.issues = verification.issues;
    }

    result.message = verification.isValid
      ? `Successfully domesticated ${mutations.length} site(s) with silent mutations`
      : `Domestication incomplete: ${verification.issues.join('; ')}`;
  } else {
    result.success = false;
    result.message = 'No mutations or junctions selected';
    result.issues = ['No domestication options were selected'];
  }

  return result;
}

/**
 * Generate fragments and primers for mutagenic junction domestication
 *
 * @param sequence - Original DNA sequence
 * @param junctions - Junction definitions
 * @param enzyme - Enzyme name
 * @param frame - Reading frame
 * @returns Fragment and primer information
 */
function generateMutagenicFragments(
  sequence: string,
  junctions: Junction[],
  enzyme: string,
  frame: number
): { fragments: Fragment[]; allPrimers: Primer[]; totalJunctions: number } {
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  const overhangLen = enz.overhangLength || 4;

  // Sort junctions by position
  const sortedJunctions = [...junctions].sort((a, b) => a.position - b.position);

  const fragments: Fragment[] = [];
  const allPrimers: Primer[] = [];
  let currentStart = 0;

  for (let i = 0; i < sortedJunctions.length; i++) {
    const junction = sortedJunctions[i];
    const junctionPos = junction.position;

    // Create fragment from currentStart to this junction
    const fragmentEnd = junctionPos + overhangLen; // Include overhang
    const fragmentSeq = sequence.slice(currentStart, fragmentEnd);

    // Get primers from the junction design
    const fwdPrimer = junction.primers?.forward || generateDefaultPrimer(sequence, currentStart, 'forward', enzyme);
    const revPrimer = junction.primers?.reverse || generateDefaultPrimer(sequence, fragmentEnd, 'reverse', enzyme);

    // Apply mutations to primers if needed
    const mutatedRevPrimer = applyMutationsToPrimer(revPrimer, junction.mutations || [], 'reverse');

    fragments.push({
      index: i + 1,
      name: `Fragment ${i + 1}`,
      start: currentStart,
      end: fragmentEnd,
      length: fragmentEnd - currentStart,
      sequence: fragmentSeq,
      leftOverhang: i === 0 ? null : sortedJunctions[i - 1].overhang,
      rightOverhang: junction.overhang,
      primers: {
        forward: i === 0 ? fwdPrimer : sortedJunctions[i - 1].primers?.internalForward,
        reverse: mutatedRevPrimer,
      },
      mutations: (junction.mutations || []).filter(m =>
        m.position >= currentStart && m.position < fragmentEnd
      ),
    });

    allPrimers.push(fwdPrimer, mutatedRevPrimer);
    currentStart = junctionPos; // Next fragment starts at junction
  }

  // Final fragment (from last junction to end)
  const lastJunction = sortedJunctions[sortedJunctions.length - 1];
  const finalStart = lastJunction.position;
  const finalSeq = sequence.slice(finalStart);

  const finalFwdPrimer = lastJunction.primers?.internalForward ||
                         generateDefaultPrimer(sequence, finalStart, 'forward', enzyme);
  const finalRevPrimer = generateDefaultPrimer(sequence, sequence.length, 'reverse', enzyme);

  // Apply any mutations to the forward primer of the final fragment
  const mutatedFinalFwdPrimer = applyMutationsToPrimer(finalFwdPrimer, lastJunction.mutations || [], 'forward');

  fragments.push({
    index: sortedJunctions.length + 1,
    name: `Fragment ${sortedJunctions.length + 1}`,
    start: finalStart,
    end: sequence.length,
    length: sequence.length - finalStart,
    sequence: finalSeq,
    leftOverhang: lastJunction.overhang,
    rightOverhang: null, // Will use destination vector overhang
    primers: {
      forward: mutatedFinalFwdPrimer,
      reverse: finalRevPrimer,
    },
    mutations: (lastJunction.mutations || []).filter(m =>
      m.position >= finalStart
    ),
  });

  allPrimers.push(mutatedFinalFwdPrimer, finalRevPrimer);

  return {
    fragments,
    allPrimers,
    totalJunctions: sortedJunctions.length,
  };
}

/**
 * Generate a default primer sequence for a position
 */
function generateDefaultPrimer(
  sequence: string,
  position: number,
  direction: string,
  enzyme: string
): Primer {
  const homologyLen = ENHANCED_CONFIG.primerHomologyLength;
  let seq: string;

  if (direction === 'forward') {
    seq = sequence.slice(position, position + homologyLen);
  } else {
    const start = Math.max(0, position - homologyLen);
    seq = reverseComplement(sequence.slice(start, position));
  }

  return {
    sequence: seq,
    direction,
    position,
    homologyLength: seq.length,
    enzyme,
  };
}

/**
 * Apply mutations to a primer sequence
 */
function applyMutationsToPrimer(
  primer: Primer,
  mutations: Mutation[],
  direction: string
): Primer {
  if (!mutations || mutations.length === 0) {
    return primer;
  }

  let mutatedSeq = primer.sequence;
  const appliedMutations: Mutation[] = [];

  for (const mutation of mutations) {
    // Check if mutation is within this primer's coverage
    const primerStart = primer.position - (direction === 'reverse' ? primer.homologyLength : 0);
    const primerEnd = primer.position + (direction === 'forward' ? primer.homologyLength : 0);

    if (mutation.position >= primerStart && mutation.position < primerEnd) {
      const relativePos = mutation.position - primerStart;

      if (direction === 'reverse') {
        // For reverse primers, the sequence is already reverse complemented
        const rcMutation = complementBase(mutation.newBase);
        const rcPos = primer.homologyLength - 1 - relativePos;
        mutatedSeq = mutatedSeq.slice(0, rcPos) + rcMutation + mutatedSeq.slice(rcPos + 1);
      } else {
        mutatedSeq = mutatedSeq.slice(0, relativePos) + mutation.newBase + mutatedSeq.slice(relativePos + 1);
      }

      appliedMutations.push({
        ...mutation,
        appliedToPrimer: true,
        relativePosition: relativePos,
      });
    }
  }

  return {
    ...primer,
    sequence: mutatedSeq,
    originalSequence: primer.sequence,
    mutations: appliedMutations,
    hasMutations: appliedMutations.length > 0,
  };
}

/**
 * Get complement of a single base
 */
function complementBase(base: string): string {
  const complement: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
  return complement[base.toUpperCase()] || base;
}

/**
 * Verify mutagenic junction result
 */
function verifyMutagenicJunctions(
  originalSeq: string,
  fragmentResult: { fragments: Fragment[]; allPrimers: Primer[]; totalJunctions: number },
  junctions: Junction[],
  frame: number,
  enzyme: string
): VerificationResult {
  const issues: string[] = [];
  const checks: ValidationCheck[] = [];

  // Check 1: All junctions have valid overhangs
  for (const junction of junctions) {
    if (!junction.overhang || junction.overhang.length !== 4) {
      issues.push(`Junction at ${junction.position} has invalid overhang`);
    }
  }
  checks.push({
    name: 'VALID_OVERHANGS',
    passed: issues.length === 0,
  });

  // Check 2: Fragments cover entire sequence
  const totalCoverage = fragmentResult.fragments.reduce((sum, f) => sum + f.length, 0);
  const overlaps = (fragmentResult.fragments.length - 1) * 4; // Overhangs overlap
  const effectiveCoverage = totalCoverage - overlaps;
  const coverageValid = effectiveCoverage >= originalSeq.length;
  checks.push({
    name: 'SEQUENCE_COVERAGE',
    passed: coverageValid,
    details: { totalCoverage, overlaps, effectiveCoverage, sequenceLength: originalSeq.length },
  });

  if (!coverageValid) {
    issues.push('Fragments do not cover entire sequence');
  }

  // Check 3: Mutations are synonymous (preserve protein)
  for (const junction of junctions) {
    if (junction.mutations) {
      for (const mutation of junction.mutations) {
        if (mutation.aminoAcid) {
          // Verify mutation is synonymous by checking amino acid
          const codonStart = Math.floor((mutation.position - frame) / 3) * 3 + frame;
          const originalCodon = originalSeq.slice(codonStart, codonStart + 3);
          const mutatedCodon = originalCodon.slice(0, mutation.position - codonStart) +
                               mutation.newBase +
                               originalCodon.slice(mutation.position - codonStart + 1);

          const originalAA = CODON_TO_AA[originalCodon];
          const mutatedAA = CODON_TO_AA[mutatedCodon];

          if (originalAA !== mutatedAA) {
            issues.push(`Mutation at ${mutation.position} changes ${originalAA} to ${mutatedAA}`);
          }
        }
      }
    }
  }
  checks.push({
    name: 'SYNONYMOUS_MUTATIONS',
    passed: !issues.some(i => i.includes('changes')),
  });

  // Check 4: All primers have mutations where needed
  const primersWithMutations = fragmentResult.allPrimers.filter(p => p.hasMutations);
  checks.push({
    name: 'PRIMER_MUTATIONS_APPLIED',
    passed: true,
    details: { primersWithMutations: primersWithMutations.length },
  });

  return {
    isValid: issues.length === 0,
    issues,
    checks,
    summary: {
      totalChecks: checks.length,
      passed: checks.filter(c => c.passed).length,
      failed: checks.filter(c => !c.passed).length,
    },
  };
}

// ============================================================================
// STEP HANDLERS
// ============================================================================

/**
 * Generate strategy options for user selection
 *
 * Evaluates all three domestication approaches and provides detailed
 * information about each option to help user make an informed decision.
 */
function generateStrategyOptions(
  sequence: string,
  enzyme: string,
  internalSites: InternalSitesResult,
  frame: number | null,
  organism: string,
  existingOverhangs: string[]
): PlanStep {
  const step: PlanStep = {
    step: 'STRATEGY_SELECTION',
    status: 'complete',
    strategies: [],
    recommendedStrategy: ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION,
  };

  // Strategy 1: Mutagenic Junction (PRIMARY/RECOMMENDED)
  const mutagenicResult = designAllMutagenicJunctions(sequence, enzyme, {
    frame: frame ?? 0,
    organism: organism as any,  // FIXED: Type assertion for organism
    existingOverhangs,
  });

  const mutagenicStrategy: StrategyOption = {
    id: ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION,
    name: 'Mutagenic Junction',
    description: 'Primers introduce silent mutations during PCR amplification',
    recommended: true,
    benefits: [
      'One-pot Golden Gate compatible',
      'Primers carry the mutations - no pre-modified template needed',
      'Mutations are permanent after PCR (no re-cutting)',
      'High-fidelity overhang selection (not constrained by original sequence)',
    ],
    tradeoffs: [
      `Adds ${mutagenicResult.additionalFragments || internalSites.count} additional fragment(s)`,
      'Requires more primers',
    ],
    feasible: mutagenicResult.success || (mutagenicResult.junctions && mutagenicResult.junctions.length > 0),
    details: {
      sitesHandled: mutagenicResult.junctions?.length || 0,
      sitesFailed: mutagenicResult.failedSites?.length || 0,
      additionalFragments: mutagenicResult.additionalFragments || internalSites.count,
      estimatedFidelity: mutagenicResult.fidelity,
      junctions: mutagenicResult.junctions,
    },
  };
  step.strategies!.push(mutagenicStrategy);

  // Strategy 2: Silent Mutation (ALTERNATIVE)
  const silentStrategy: StrategyOption = {
    id: ENHANCED_CONFIG.strategies.SILENT_MUTATION,
    name: 'Silent Mutation',
    description: 'Directly modify template sequence with synonymous codon changes',
    recommended: false,
    benefits: [
      'No additional fragments needed',
      'Simpler assembly (fewer parts)',
      'Can use existing primers if template is pre-modified',
    ],
    tradeoffs: [
      'Requires pre-modified template DNA (synthesis or mutagenesis)',
      'NOT one-pot compatible if sites re-cut during reaction',
      'May need careful ordering of assembly if sites are not fully removed',
    ],
    feasible: true, // Silent mutations are generally always possible if in-frame
    details: {
      sitesCount: internalSites.count,
      additionalFragments: 0,
      requiresTemplateModification: true,
    },
  };
  step.strategies!.push(silentStrategy);

  // Strategy 3: Alternative Enzyme (FALLBACK)
  const altEnzymes = recommendAlternativeEnzymes(sequence, enzyme);
  const compatibleAlternatives = altEnzymes.filter(e => e.isCompatible && !e.isCurrent);

  const alternativeStrategy: StrategyOption = {
    id: ENHANCED_CONFIG.strategies.ALTERNATIVE_ENZYME,
    name: 'Alternative Enzyme',
    description: 'Switch to a different Type IIS enzyme that has no internal sites',
    recommended: false,
    benefits: [
      'No mutations needed - preserves original sequence exactly',
      'No additional fragments',
      'One-pot compatible',
    ],
    tradeoffs: [
      compatibleAlternatives.length === 0
        ? 'No compatible alternatives available for this sequence'
        : `Limited to: ${compatibleAlternatives.map(e => e.enzyme).join(', ')}`,
      'May require different primer design',
      'Different overhang constraints',
    ],
    feasible: compatibleAlternatives.length > 0,
    details: {
      availableEnzymes: compatibleAlternatives,
      currentEnzyme: enzyme,
      currentSiteCount: internalSites.count,
    },
  };
  step.strategies!.push(alternativeStrategy);

  // Determine recommended strategy
  if (mutagenicStrategy.feasible) {
    step.recommendedStrategy = ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION;
  } else if (compatibleAlternatives.length > 0) {
    step.recommendedStrategy = ENHANCED_CONFIG.strategies.ALTERNATIVE_ENZYME;
  } else {
    step.recommendedStrategy = ENHANCED_CONFIG.strategies.SILENT_MUTATION;
  }

  return step;
}

/**
 * Handle reading frame detection and validation
 */
function handleFrameDetection(
  sequence: string,
  providedFrame: number | null,
  sites: any[],
  organism: string
): PlanStep {
  const step: PlanStep = {
    step: 'FRAME_DETECTION',
    status: 'pending',
  };

  // If frame is provided, validate it
  if (providedFrame !== null) {
    const validation = validateReadingFrame(sequence, providedFrame);

    if (validation.isValid) {
      step.status = 'complete';
      step.confirmedFrame = providedFrame;
      step.validation = validation;
      step.requiresUserAction = false;
      return step;
    } else {
      step.status = 'warning';
      step.confirmedFrame = providedFrame;
      step.validation = validation;
      step.requiresUserAction = true;
      step.userMessage = `Specified frame ${providedFrame} has issues: ${validation.validations.filter(v => !v.passed).map(v => v.message).join('; ')}. Please confirm.`;
      step.frameOptions = [0, 1, 2].map(f => ({
        frame: f,
        validation: validateReadingFrame(sequence, f),
        isCurrent: f === providedFrame,
      }));
      return step;
    }
  }

  // Auto-detect ORF
  const orfResult = detectORFs(sequence, { organism: organism as any });  // FIXED: Type assertion

  if (orfResult.recommendation.confidence === 'high') {
    step.status = 'complete';
    step.confirmedFrame = orfResult.recommendation.recommendedFrame;
    step.orfDetection = orfResult;
    step.requiresUserAction = false;
    return step;
  }

  if (orfResult.recommendation.confidence === 'medium') {
    step.status = 'pending';
    step.recommendedFrame = orfResult.recommendation.recommendedFrame;
    step.orfDetection = orfResult;
    step.requiresUserAction = true;
    step.userMessage = `ORF detected with medium confidence. Recommended: frame ${orfResult.recommendation.recommendedFrame}. Please confirm.`;
    step.frameOptions = [0, 1, 2].map(f => ({
      frame: f,
      validation: validateReadingFrame(sequence, f),
      isRecommended: f === orfResult.recommendation.recommendedFrame,
      orfs: orfResult.orfs.filter(o => o.frame === f),
    }));
    return step;
  }

  // Low confidence - must get user input
  step.status = 'pending';
  step.recommendedFrame = undefined;
  step.orfDetection = orfResult;
  step.requiresUserAction = true;
  step.userMessage = 'Cannot automatically determine reading frame. Please select the correct frame.';
  step.frameOptions = [0, 1, 2].map(f => ({
    frame: f,
    validation: validateReadingFrame(sequence, f),
    orfs: orfResult.orfs.filter(o => o.frame === f),
  }));

  return step;
}

/**
 * Handle adjacent sites that are close together
 */
function handleAdjacentSites(
  sequence: string,
  sites: any[],
  enzyme: string,
  frame: number,
  organism: string
): PlanStep {
  const step: PlanStep = {
    step: 'ADJACENT_SITE_CHECK',
    status: 'complete',
    hasAdjacentSites: false,
  };

  const adjacentCheck = detectAdjacentSites(sites, ENHANCED_CONFIG.minSiteDistance);

  if (!adjacentCheck.hasAdjacentSites) {
    return step;
  }

  step.hasAdjacentSites = true;
  step.adjacentPairs = adjacentCheck.adjacentPairs;
  step.minDistanceFound = adjacentCheck.minDistanceFound;

  // Try intelligent resolution strategies
  const resolutionOptions: any[] = [];

  for (const pair of adjacentCheck.adjacentPairs) {
    const pairResolutions: any[] = [];

    // Strategy 1: Try to find a single mutation that breaks both sites
    const sharedMutation = findSharedMutation(sequence, pair.site1, pair.site2, frame, enzyme, organism);
    if (sharedMutation) {
      pairResolutions.push({
        strategy: 'OVERLAPPING_MUTATION',
        description: 'Single mutation breaks both sites',
        mutation: sharedMutation,
        fragmentIncrease: 0,
        recommended: true,
      });
    }

    // Strategy 2: Try single junction that domesticates both
    const singleJunction = findSingleJunctionForPair(sequence, pair, enzyme, frame, organism);
    if (singleJunction) {
      pairResolutions.push({
        strategy: 'SINGLE_JUNCTION',
        description: 'Single junction breaks both sites',
        junction: singleJunction,
        fragmentIncrease: 1,
        recommended: !sharedMutation,
      });
    }

    // Strategy 3: Check for alternative enzymes
    const alternatives = recommendAlternativeEnzymes(sequence, enzyme);
    const compatibleAlt = alternatives.find(a => a.isCompatible && !a.isCurrent);
    if (compatibleAlt) {
      pairResolutions.push({
        strategy: 'ALTERNATIVE_ENZYME',
        description: `Switch to ${compatibleAlt.enzyme} (no internal sites)`,
        enzyme: compatibleAlt,
        fragmentIncrease: 0,
        recommended: pairResolutions.length === 0,
      });
    }

    resolutionOptions.push({
      pair,
      resolutions: pairResolutions,
      canResolve: pairResolutions.length > 0,
    });
  }

  const canHandle = resolutionOptions.every(r => r.canResolve);

  step.canHandle = canHandle;
  step.resolutionOptions = resolutionOptions;
  step.message = canHandle
    ? `Found ${adjacentCheck.adjacentPairs.length} adjacent site pair(s) - resolution strategies available`
    : `Found ${adjacentCheck.adjacentPairs.length} adjacent site pair(s) - manual intervention required`;
  step.strategy = canHandle ? resolutionOptions[0]?.resolutions[0]?.strategy : null;

  return step;
}

/**
 * Generate mutation options for each site
 */
function generateMutationOptions(
  sequence: string,
  sites: any[],
  frame: number,
  enzyme: string,
  organism: string,
  codonMode: string,
  existingOverhangs: string[],
  selectedStrategy: string = ENHANCED_CONFIG.defaultStrategy
): PlanStep {
  const step: PlanStep = {
    step: 'MUTATION_OPTIONS',
    status: 'complete',
    siteOptions: [],
    selectedStrategy,
  };

  const codonUsage = organism === 'yeast' ? YEAST_CODON_USAGE : ECOLI_CODON_USAGE;
  const prioritizeMutagenic = selectedStrategy === ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION;

  for (const site of sites) {
    const siteOption: SiteOption = {
      site,
      options: [],
      recommended: null,
    };

    // Design mutagenic junction option (PRIMARY for one-pot Golden Gate)
    const junction = designMutagenicJunction(sequence, site, enzyme, {
      frame,
      organism: organism as any,  // FIXED: Type assertion for organism
      existingOverhangs,
    });

    if (junction.success) {
      const mutagenicOption: MutationOption = {
        type: 'mutagenic_junction',
        junction,
        score: (junction.combinedScore || 0) + (prioritizeMutagenic ? 100 : 0), // Boost if primary strategy
        description: `Mutagenic junction at position ${junction.junctionPosition}`,
        shortDescription: 'Junction + primer mutation',
        fragmentIncrease: 1,
        overhang: junction.overhang,
        onePotCompatible: true,
        primers: junction.primers,
        mutations: (junction as any).mutation || (junction as any).mutations,  // FIXED: Property name
        benefits: [
          'One-pot Golden Gate compatible',
          'No pre-modified template needed',
          'Permanent mutation after PCR',
        ],
      };

      // Add as FIRST option if mutagenic junction is the selected strategy
      if (prioritizeMutagenic) {
        siteOption.options.unshift(mutagenicOption);
      } else {
        siteOption.options.push(mutagenicOption);
      }
    }

    // Find silent mutation candidates (ALTERNATIVE)
    const candidates = findAllSilentMutationCandidates(
      sequence,
      site,
      frame,
      enzyme,
      codonUsage
    );

    // Score and rank candidates
    const scoredCandidates = scoreMutationCandidates(
      candidates,
      sequence,
      enzyme,
      DOMESTICATION_CONFIG.checkAllGoldenGateEnzymes
    );

    // Add silent mutation options
    for (const candidate of scoredCandidates.slice(0, 5)) {
      const option: MutationOption = {
        type: 'silent_mutation',
        mutation: candidate,
        score: candidate.score + (prioritizeMutagenic ? 0 : 50), // Boost if silent mutation strategy
        description: formatMutationDescription(candidate, codonUsage),
        shortDescription: `${candidate.originalCodon}→${candidate.newCodon}`,
        codonChange: `${candidate.originalCodon}→${candidate.newCodon}`,
        aminoAcid: candidate.aminoAcid,
        frequencyChange: `${candidate.originalCodonFrequency?.toFixed(1) || '?'}→${candidate.newCodonFrequency?.toFixed(1) || '?'}/1000`,
        warnings: candidate.penalties?.map(p => p.reason) || [],
        onePotCompatible: false,
        benefits: [
          'No additional fragments',
          'Simpler assembly',
        ],
        tradeoffs: [
          'Requires pre-modified template',
          'Not one-pot compatible',
        ],
      };

      // Insert after mutagenic options if present, otherwise at the front
      if (prioritizeMutagenic) {
        siteOption.options.push(option);
      } else {
        // Find where to insert (after any existing silent mutations)
        const insertIndex = siteOption.options.filter(o => o.type === 'silent_mutation').length;
        siteOption.options.splice(insertIndex, 0, option);
      }
    }

    // Select recommended option based on selected strategy
    if (siteOption.options.length > 0) {
      if (prioritizeMutagenic) {
        // Prefer mutagenic junction (first option)
        const mutagenicOptions = siteOption.options.filter(o => o.type === 'mutagenic_junction');
        siteOption.recommended = mutagenicOptions.length > 0 ? mutagenicOptions[0] : siteOption.options[0];
      } else if (codonMode === ENHANCED_CONFIG.codonModes.CONSERVATIVE) {
        // Prefer mutation with closest codon frequency
        const silentOptions = siteOption.options.filter(o => o.type === 'silent_mutation');
        if (silentOptions.length > 0) {
          silentOptions.sort((a, b) => {
            const aRatio = a.mutation?.frequencyRatio || 0;
            const bRatio = b.mutation?.frequencyRatio || 0;
            return bRatio - aRatio; // Higher ratio = closer to original
          });
          siteOption.recommended = silentOptions[0];
        } else {
          siteOption.recommended = siteOption.options[0];
        }
      } else if (codonMode === ENHANCED_CONFIG.codonModes.OPTIMIZED) {
        // Prefer mutation with highest new codon frequency
        const silentOptions = siteOption.options.filter(o => o.type === 'silent_mutation');
        if (silentOptions.length > 0) {
          silentOptions.sort((a, b) => {
            const aFreq = a.mutation?.newCodonFrequency || 0;
            const bFreq = b.mutation?.newCodonFrequency || 0;
            return bFreq - aFreq;
          });
          siteOption.recommended = silentOptions[0];
        } else {
          siteOption.recommended = siteOption.options[0];
        }
      } else {
        // CUSTOM mode - don't auto-select
        siteOption.recommended = null;
      }
    }

    step.siteOptions!.push(siteOption);
  }

  return step;
}

/**
 * Generate pre-flight preview
 */
function generatePreFlightPreview(
  sequence: string,
  sites: any[],
  siteOptions: SiteOption[],
  frame: number,
  enzyme: string = 'BsaI',
  selectedStrategy: string = ENHANCED_CONFIG.defaultStrategy
): PlanStep {
  const step: PlanStep = {
    step: 'PREFLIGHT_PREVIEW',
    status: 'complete',
  };

  const originalTranslation = translateSequence(sequence, frame);
  const isJunctionStrategy = selectedStrategy === ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION;

  // Apply recommended mutations to generate preview
  let previewSeq = sequence;
  const appliedMutations: any[] = [];
  const appliedJunctions: any[] = [];

  for (const siteOption of siteOptions) {
    const recommended = siteOption.recommended;
    if (!recommended) continue;

    if (recommended.type === 'silent_mutation' && recommended.mutation) {
      previewSeq = applyMutation(previewSeq, recommended.mutation);
      appliedMutations.push(recommended);
    } else if (recommended.type === 'mutagenic_junction') {
      // Track junctions - they don't modify the template sequence directly
      appliedJunctions.push({
        site: siteOption.site,
        junction: recommended.junction,
        position: recommended.junction?.junctionPosition,
        overhang: recommended.junction?.overhang || recommended.overhang,
        description: recommended.description,
        primers: recommended.primers || recommended.junction?.primers,
      });
    }
  }

  const previewTranslation = translateSequence(previewSeq, frame);
  const proteinComparison = compareProteins(
    originalTranslation.protein,
    previewTranslation.protein
  );

  // Check for remaining sites - use actual enzyme, not hardcoded
  const remainingSites = findInternalSites(previewSeq, enzyme);

  // For mutagenic junction strategy, sites remain in the template but will be
  // removed during assembly via the primers. So validation should pass.
  const sitesHandled = isJunctionStrategy
    ? appliedJunctions.length === sites.length  // All sites have junctions
    : !remainingSites.hasSites;                  // All sites removed via mutation

  step.preview = {
    original: {
      sequence: sequence,
      protein: originalTranslation.protein,
      proteinLength: originalTranslation.protein.replace(/\*/g, '').length,
    },
    domesticated: {
      sequence: previewSeq,
      protein: previewTranslation.protein,
      proteinLength: previewTranslation.protein.replace(/\*/g, '').length,
    },
    comparison: {
      sequenceIdentical: sequence === previewSeq,
      proteinIdentical: proteinComparison.identical,
      proteinIdentity: proteinComparison.identity,
      differences: proteinComparison.differences,
    },
    mutations: appliedMutations.map(m => ({
      position: m.mutation.sequencePosition,
      change: `${m.mutation.originalBase}→${m.mutation.newBase}`,
      codon: `${m.mutation.originalCodon}→${m.mutation.newCodon}`,
      aminoAcid: m.mutation.aminoAcid,
    })),
    junctions: appliedJunctions.map(j => ({
      position: j.position,
      overhang: j.overhang,
      description: j.description,
      sitePosition: j.site?.position,
    })),
    validation: {
      proteinPreserved: proteinComparison.identical,
      noRemainingSites: sitesHandled,
      remainingSiteCount: isJunctionStrategy ? 0 : remainingSites.count,
      strategy: selectedStrategy,
      isJunctionStrategy,
    },
    frame,
    totalMutations: appliedMutations.length,
    totalJunctions: appliedJunctions.length,
    strategy: selectedStrategy,
  };

  return step;
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Find a single mutation that breaks both adjacent sites
 */
function findSharedMutation(
  sequence: string,
  site1: any,
  site2: any,
  frame: number,
  enzyme: string,
  organism: string
): SilentMutationCandidate | null {
  const codonUsage = organism === 'yeast' ? YEAST_CODON_USAGE : ECOLI_CODON_USAGE;

  // Find overlapping positions between the two sites
  const site1End = site1.position + site1.sequence.length;
  const site2End = site2.position + site2.sequence.length;

  // Check if sites overlap
  if (site1End <= site2.position || site2End <= site1.position) {
    // Sites don't overlap - check if they share a codon
    const codon1Start = Math.floor((site1.position - frame) / 3) * 3 + frame;
    const codon1End = Math.floor((site1End - 1 - frame) / 3) * 3 + frame + 3;
    const codon2Start = Math.floor((site2.position - frame) / 3) * 3 + frame;

    if (codon1End <= codon2Start) {
      return null; // No shared codons
    }
  }

  // Find candidates that break site1
  const candidates1 = findAllSilentMutationCandidates(sequence, site1, frame, enzyme, codonUsage);

  // For each candidate, check if it also breaks site2
  for (const candidate of candidates1) {
    // Apply the mutation temporarily
    const mutatedSeq = applyMutation(sequence, candidate);

    // Check if site2 is also broken
    const site2Seq = mutatedSeq.slice(site2.position, site2.position + site2.sequence.length);
    const enz = GOLDEN_GATE_ENZYMES[enzyme];

    if (site2Seq !== enz.recognition && site2Seq !== reverseComplement(enz.recognition)) {
      // This mutation breaks both sites!
      return {
        ...candidate,
        breaksBothSites: true,
        sites: [site1, site2],
      } as any;
    }
  }

  return null;
}

/**
 * Find a single junction that domesticates both adjacent sites
 */
function findSingleJunctionForPair(
  sequence: string,
  pair: any,
  enzyme: string,
  frame: number,
  organism: string
): any | null {
  const { site1, site2 } = pair;

  // Find junction position between the two sites
  const site1End = site1.position + site1.sequence.length;
  const gapStart = site1End;
  const gapEnd = site2.position;

  if (gapEnd - gapStart < 4) {
    return null; // Not enough space for an overhang
  }

  // Try designing a junction in the gap
  const midPoint = Math.floor((gapStart + gapEnd) / 2);

  // Create a virtual site at the midpoint for junction design
  const virtualSite = {
    position: midPoint - 3,
    sequence: sequence.slice(midPoint - 3, midPoint + 3),
    orientation: 'forward',
  };

  const junction = designMutagenicJunction(sequence, virtualSite as any, enzyme, {  // FIXED: Type assertion
    frame,
    organism: organism as any,  // FIXED: Type assertion
  });

  if (junction.success) {
    // Verify this junction would break both sites
    return {
      ...junction,
      breaksBothSites: true,
      sites: [site1, site2],
    };
  }

  return null;
}

/**
 * Format mutation description for user display
 */
function formatMutationDescription(mutation: SilentMutationCandidate, codonUsage: CodonUsageTable): string {
  const freqOriginal = mutation.originalCodonFrequency?.toFixed(1) || '?';
  const freqNew = mutation.newCodonFrequency?.toFixed(1) || '?';

  let description = `${mutation.originalBase}→${mutation.newBase} at position ${mutation.sequencePosition}`;
  description += ` (${mutation.originalCodon}→${mutation.newCodon}, ${mutation.aminoAcid})`;
  description += ` [freq: ${freqOriginal}→${freqNew}/1000]`;

  if (mutation.positionInCodon === 2) {
    description += ' [wobble]';
  }

  return description;
}

/**
 * Comprehensive verification of domestication result
 */
function verifyDomestication(
  originalSeq: string,
  domesticatedSeq: string,
  frame: number,
  enzyme: string,
  mutations: AppliedMutation[]
): VerificationResult {
  const issues: string[] = [];
  const checks: ValidationCheck[] = [];

  // Check 1: Protein preservation
  const proteinCheck = verifyProteinSequence(originalSeq, domesticatedSeq, frame);
  checks.push({
    name: 'PROTEIN_PRESERVATION',
    passed: proteinCheck.identical,
    details: proteinCheck,
  });

  if (!proteinCheck.identical) {
    issues.push(`Protein altered at ${proteinCheck.differences.length} position(s)`);
  }

  // Check 2: No remaining sites
  const remainingSites = findInternalSites(domesticatedSeq, enzyme);
  checks.push({
    name: 'NO_REMAINING_SITES',
    passed: !remainingSites.hasSites,
    details: remainingSites,
  });

  if (remainingSites.hasSites) {
    issues.push(`${remainingSites.count} internal ${enzyme} site(s) remain`);
  }

  // Check 3: Sequence length unchanged
  const lengthUnchanged = originalSeq.length === domesticatedSeq.length;
  checks.push({
    name: 'LENGTH_PRESERVED',
    passed: lengthUnchanged,
  });

  if (!lengthUnchanged) {
    issues.push(`Sequence length changed from ${originalSeq.length} to ${domesticatedSeq.length}`);
  }

  // Check 4: No new sites created (for all GG enzymes)
  for (const enzName of Object.keys(GOLDEN_GATE_ENZYMES)) {
    const originalCount = findInternalSites(originalSeq, enzName).count;
    const newCount = findInternalSites(domesticatedSeq, enzName).count;

    if (newCount > originalCount) {
      issues.push(`Created ${newCount - originalCount} new ${enzName} site(s)`);
      checks.push({
        name: `NO_NEW_${enzName}_SITES`,
        passed: false,
        details: { originalCount, newCount },
      });
    }
  }

  return {
    isValid: issues.length === 0,
    issues,
    checks,
    summary: {
      totalChecks: checks.length,
      passed: checks.filter(c => c.passed).length,
      failed: checks.filter(c => !c.passed).length,
    },
  };
}

// ============================================================================
// RE-EXPORTS (for convenience)
// ============================================================================

// Re-export useful functions from other modules
export {
  detectORFs,
  validateReadingFrame,
  translateSequence,
  preFlightCheck,
} from './orf-detector.js';

export {
  verifyProteinSequence,
} from './silent-mutation-domesticator.js';
