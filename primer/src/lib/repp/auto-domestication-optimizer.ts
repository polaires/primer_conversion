/**
 * Auto-Domestication Optimizer for Golden Gate Assembly
 *
 * Automatically detects internal restriction sites and integrates domestication
 * into the Fusion Site Optimizer workflow. Instead of mutating the sequence,
 * this approach places junction positions within internal sites, using the
 * assembly itself to break and remove problematic sites.
 *
 * Key concept: Each internal site becomes an additional fragment boundary,
 * with primers designed to create overhangs that don't recreate the site.
 */

import { GOLDEN_GATE_ENZYMES, calculateExperimentalFidelity, getEnzymeLigationData, findInternalSites } from './goldengate.js';
import { reverseComplement } from './enzymes.js';
import {
  validateOverhang,
  checkSiteRecreation,
} from './overhang-validation.js';

/**
 * Configuration for auto-domestication
 */
export const DOMESTICATION_DEFAULTS = {
  // Prefer junction positions that maximize overhang quality
  preferHighFidelityOverhangs: true,
  // Minimum distance from site edges for junction
  minDistanceFromSiteEdge: 0,
  // Consider biological context (coding regions)
  respectCodingFrames: true,
  // Allow manual override of domestication junctions
  allowManualOverride: true,
  // Minimum distance between adjacent sites for auto-domestication to work
  minSiteDistance: 50,
  // Default minimum fragment size for validation
  minFragmentSize: 50,
} as const;

// ============================================================================
// TYPES
// ============================================================================

interface InternalSite {
  position: number;
  sequence: string;
  orientation: 'forward' | 'reverse';
  index?: number;
}

interface AdjacentSitePair {
  site1: InternalSite & { index: number };
  site2: InternalSite & { index: number };
  distance: number;
  requiredDistance: number;
}

interface AdjacentSitesAnalysis {
  hasAdjacentSites: boolean;
  adjacentPairs: AdjacentSitePair[];
  minDistanceFound: number;
  totalSites: number;
}

interface AlternativeEnzyme {
  enzyme: string;
  fullName: string;
  recognition: string;
  overhangLength: number;
  internalSites: number;
  forwardSites: number;
  reverseSites: number;
  isCompatible: boolean;
  hasLigationData: boolean;
  isCurrent: boolean;
}

interface Fragment {
  start: number;
  end: number;
  size: number;
  index: number;
}

interface FragmentViolation {
  type: 'TOO_SMALL' | 'TOO_LARGE';
  fragmentIndex: number;
  start: number;
  end: number;
  size: number;
  minRequired?: number;
  maxAllowed?: number;
  message: string;
}

interface FragmentSizeValidation {
  valid: boolean;
  fragments: Fragment[];
  violations: FragmentViolation[];
  minFragmentSize: number;
  maxFragmentSize: number;
}

interface OverhangQuality {
  score: number;
  issues: string[];
  benefits: string[];
  gc: number;
  isPalindrome: boolean;
}

interface JunctionCandidate {
  position: number;
  overhang: string;
  quality: OverhangQuality;
  breaksForward: boolean;
  breaksReverse: boolean;
  recreatesSite: boolean;
  distanceFromSiteCenter: number;
  isValid: boolean;
}

interface DomesticationJunctionOptions {
  site: InternalSite;
  sitePosition: number;
  siteSequence: string;
  orientation: 'forward' | 'reverse';
  candidates: JunctionCandidate[];
  recommended: JunctionCandidate | undefined;
  hasValidOption: boolean;
}

interface DomesticationError {
  type: string;
  message: string;
  details?: any;
  minDistanceFound?: number;
}

interface DomesticationAnalysis {
  needsDomestication: boolean;
  sites: InternalSite[];
  domesticationOptions?: DomesticationJunctionOptions[];
  additionalFragments: number;
  error: DomesticationError | null;
  alternativeEnzymes: AlternativeEnzyme[] | null;
  recommendedEnzyme?: AlternativeEnzyme;
  message: string;
}

interface RequiredJunction {
  position: number;
  overhang: string;
  isRequired: boolean;
  reason: string;
  site: InternalSite;
  quality: OverhangQuality;
}

interface ExcludeZone {
  start: number;
  end: number;
}

interface JunctionTarget {
  idealPosition: number;
  searchRadius: number;
  index: number;
}

interface DomesticationOptimizationResult {
  needsDomestication: boolean;
  internalSites: InternalSite[];
  domesticationJunctions: RequiredJunction[];
  additionalFragmentsNeeded: number;
  totalFragments: number;
  totalJunctions: number;
  remainingJunctionsToOptimize: number;
  targetPositions: JunctionTarget[];
  message: string;
  constraints: {
    excludeZones: ExcludeZone[];
    [key: string]: any;
  };
}

interface ValidationIssue {
  type: string;
  code: string;
  message: string;
  details?: any;
  recommendation?: string;
}

interface ValidationCheck {
  code: string;
  message: string;
}

interface PostDomesticationValidation {
  isValid: boolean;
  issues: ValidationIssue[];
  warnings: ValidationIssue[];
  passed: ValidationCheck[];
  summary: {
    fragmentCount: number;
    totalIssues: number;
    totalWarnings: number;
    totalPassed: number;
    status: 'valid' | 'valid_with_warnings' | 'invalid';
  };
}

interface FragmentWithOverhang {
  id?: string;
  seq?: string;
  _overhang?: string;
}

interface SiteGroup {
  sites: InternalSite[];
  type: 'single' | 'adjacent' | 'clustered';
  span?: number;
  canHandle: boolean;
  strategy: string;
  reason?: string;
  fragmentCount?: number;
}

interface MergedSitesResult {
  groups: SiteGroup[];
  needsMerging: boolean;
  allCanHandle?: boolean;
  totalSites?: number;
  groupCount?: number;
  strategy?: string;
}

interface DomesticationSelection {
  overhang: string;
  position: number;
  quality: OverhangQuality;
}

interface OverhangValidation {
  isValid: boolean;
  duplicates: string[];
  setFidelity: number;
  problematicPairs: any[];
  recommendation: string;
}

interface GlobalOverhangOptimizationResult {
  success: boolean;
  selections: DomesticationSelection[];
  overhangs: string[];
  fidelity: number;
  validation: OverhangValidation;
  searchMethod: 'exhaustive' | 'greedy';
  combinationsSearched: number;
  error?: string;
  message?: string;
}

// ============================================================================
// FUNCTIONS
// ============================================================================

/**
 * Detect adjacent internal sites that are too close for auto-domestication
 * Sites closer than minDistance bp cannot be reliably split into separate fragments
 */
export function detectAdjacentSites(
  sites: InternalSite[] | null | undefined,
  minDistance: number = DOMESTICATION_DEFAULTS.minSiteDistance
): AdjacentSitesAnalysis {
  if (!sites || sites.length < 2) {
    return { hasAdjacentSites: false, adjacentPairs: [], minDistanceFound: Infinity, totalSites: sites?.length || 0 };
  }

  // Sort sites by position
  const sortedSites = [...sites].sort((a, b) => a.position - b.position);
  const adjacentPairs: AdjacentSitePair[] = [];
  let minDistanceFound = Infinity;

  for (let i = 0; i < sortedSites.length - 1; i++) {
    const site1 = sortedSites[i];
    const site2 = sortedSites[i + 1];
    const distance = site2.position - site1.position;

    if (distance < minDistanceFound) {
      minDistanceFound = distance;
    }

    if (distance < minDistance) {
      adjacentPairs.push({
        site1: { ...site1, index: i },
        site2: { ...site2, index: i + 1 },
        distance,
        requiredDistance: minDistance,
      });
    }
  }

  return {
    hasAdjacentSites: adjacentPairs.length > 0,
    adjacentPairs,
    minDistanceFound,
    totalSites: sites.length,
  };
}

/**
 * Find alternative enzymes that don't have internal sites in the sequence
 * Returns enzymes ranked by compatibility (fewer sites = better)
 * Includes current enzyme in results for comparison
 */
export function recommendAlternativeEnzymes(
  sequence: string,
  currentEnzyme: string = 'BsaI',
  options: { includeCurrentEnzyme?: boolean } = {}
): AlternativeEnzyme[] {
  const { includeCurrentEnzyme = true } = options;
  const alternatives: AlternativeEnzyme[] = [];
  const seq = sequence.toUpperCase();

  for (const [enzymeName, enzymeData] of Object.entries(GOLDEN_GATE_ENZYMES)) {
    const recognition = enzymeData.recognition;
    const recognitionRC = reverseComplement(recognition);

    // Count forward sites
    let forwardCount = 0;
    let pos = seq.indexOf(recognition);
    while (pos !== -1) {
      forwardCount++;
      pos = seq.indexOf(recognition, pos + 1);
    }

    // Count reverse sites
    let reverseCount = 0;
    pos = seq.indexOf(recognitionRC);
    while (pos !== -1) {
      reverseCount++;
      pos = seq.indexOf(recognitionRC, pos + 1);
    }

    const totalSites = forwardCount + reverseCount;
    const isCurrent = enzymeName === currentEnzyme;

    // Skip current enzyme if not included
    if (isCurrent && !includeCurrentEnzyme) continue;

    alternatives.push({
      enzyme: enzymeName,
      fullName: enzymeData.fullName || enzymeName,
      recognition: enzymeData.recognition,
      overhangLength: enzymeData.overhangLength || 4,
      internalSites: totalSites,
      forwardSites: forwardCount,
      reverseSites: reverseCount,
      isCompatible: totalSites === 0,
      hasLigationData: !!enzymeData.dataKey,
      isCurrent,
    });
  }

  // Sort by: compatible first, then by site count (ascending), then by ligation data availability
  // Current enzyme is sorted normally but marked for UI display
  alternatives.sort((a, b) => {
    if (a.isCompatible !== b.isCompatible) return Number(b.isCompatible) - Number(a.isCompatible);
    if (a.internalSites !== b.internalSites) return a.internalSites - b.internalSites;
    if (a.hasLigationData !== b.hasLigationData) return Number(b.hasLigationData) - Number(a.hasLigationData);
    return 0;
  });

  return alternatives;
}

/**
 * Validate that domestication-created fragments meet size constraints
 */
export function validateFragmentSizes(
  sequence: string,
  junctions: (number | { position: number })[] | null | undefined,
  constraints: { minFragmentSize?: number; maxFragmentSize?: number } = {}
): FragmentSizeValidation {
  const {
    minFragmentSize = DOMESTICATION_DEFAULTS.minFragmentSize,
    maxFragmentSize = 10000,
  } = constraints;

  if (!junctions || junctions.length === 0) {
    return {
      valid: true,
      fragments: [{ start: 0, end: sequence.length, size: sequence.length, index: 0 }],
      violations: [],
      minFragmentSize,
      maxFragmentSize,
    };
  }

  // Sort junction positions
  const positions = [
    0,
    ...junctions.map(j => typeof j === 'number' ? j : j.position).sort((a, b) => a - b),
    sequence.length
  ];
  const fragments: Fragment[] = [];
  const violations: FragmentViolation[] = [];

  for (let i = 0; i < positions.length - 1; i++) {
    const start = positions[i];
    const end = positions[i + 1];
    const size = end - start;

    fragments.push({ start, end, size, index: i });

    if (size < minFragmentSize) {
      violations.push({
        type: 'TOO_SMALL',
        fragmentIndex: i,
        start,
        end,
        size,
        minRequired: minFragmentSize,
        message: `Fragment ${i + 1} (${start}-${end}) is ${size}bp, below minimum ${minFragmentSize}bp`,
      });
    }

    if (size > maxFragmentSize) {
      violations.push({
        type: 'TOO_LARGE',
        fragmentIndex: i,
        start,
        end,
        size,
        maxAllowed: maxFragmentSize,
        message: `Fragment ${i + 1} (${start}-${end}) is ${size}bp, above maximum ${maxFragmentSize}bp`,
      });
    }
  }

  return {
    valid: violations.length === 0,
    fragments,
    violations,
    minFragmentSize,
    maxFragmentSize,
  };
}

/**
 * Analyze a sequence for internal sites that need domestication
 */
export function analyzeForDomestication(
  sequence: string,
  enzyme: string = 'BsaI',
  options: { minFragmentSize?: number; minSiteDistance?: number } = {}
): DomesticationAnalysis {
  const {
    minFragmentSize = DOMESTICATION_DEFAULTS.minFragmentSize,
    minSiteDistance = DOMESTICATION_DEFAULTS.minSiteDistance,
  } = options;

  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enz) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  const internalSites = findInternalSites(sequence, enzyme);

  if (!internalSites.hasSites) {
    return {
      needsDomestication: false,
      sites: [],
      additionalFragments: 0,
      message: 'Sequence is compatible - no internal sites found',
      error: null,
      alternativeEnzymes: null,
    };
  }

  // Check for adjacent sites that are too close together
  const adjacentCheck = detectAdjacentSites(internalSites.sites, minSiteDistance);

  if (adjacentCheck.hasAdjacentSites) {
    // Get alternative enzyme recommendations
    const alternatives = recommendAlternativeEnzymes(sequence, enzyme);
    const compatibleAlternatives = alternatives.filter(a => a.isCompatible);
    const bestAlternative = compatibleAlternatives[0] || alternatives[0];

    return {
      needsDomestication: true,
      sites: internalSites.sites,
      domesticationOptions: [],
      additionalFragments: internalSites.count,
      error: {
        type: 'ADJACENT_SITES_TOO_CLOSE',
        message: `Cannot auto-domesticate: ${adjacentCheck.adjacentPairs.length} pair(s) of internal ${enzyme} sites are too close together (minimum distance: ${minSiteDistance}bp).`,
        details: adjacentCheck.adjacentPairs.map(pair => ({
          site1Position: pair.site1.position,
          site2Position: pair.site2.position,
          distance: pair.distance,
          requiredDistance: minSiteDistance,
        })),
        minDistanceFound: adjacentCheck.minDistanceFound,
      },
      alternativeEnzymes: alternatives,
      recommendedEnzyme: bestAlternative,
      message: compatibleAlternatives.length > 0
        ? `Sites too close for auto-domestication. Recommended: switch to ${bestAlternative.enzyme} (${bestAlternative.fullName}) which has no internal sites.`
        : `Sites too close for auto-domestication. Consider using ${bestAlternative.enzyme} (${bestAlternative.internalSites} site${bestAlternative.internalSites !== 1 ? 's' : ''}).`,
    };
  }

  // For each internal site, find optimal junction positions
  const domesticationOptions = internalSites.sites.map((site: any) => {
    return findDomesticationJunctions(sequence, site, enzyme);
  });

  // Validate that resulting fragments meet size constraints
  const junctionPositions = domesticationOptions
    .filter((opt: any) => opt.hasValidOption && opt.recommended)
    .map((opt: any) => opt.recommended!.position);

  const fragmentValidation = validateFragmentSizes(sequence, junctionPositions, { minFragmentSize });

  if (!fragmentValidation.valid) {
    // Get alternative enzyme recommendations
    const alternatives = recommendAlternativeEnzymes(sequence, enzyme);
    const compatibleAlternatives = alternatives.filter(a => a.isCompatible);
    const bestAlternative = compatibleAlternatives[0] || alternatives[0];

    return {
      needsDomestication: true,
      sites: internalSites.sites,
      domesticationOptions,
      additionalFragments: internalSites.count,
      error: {
        type: 'FRAGMENT_SIZE_VIOLATION',
        message: `Auto-domestication would create fragments smaller than ${minFragmentSize}bp minimum.`,
        details: {
          violations: fragmentValidation.violations,
          fragments: fragmentValidation.fragments,
        },
      },
      alternativeEnzymes: alternatives,
      recommendedEnzyme: bestAlternative,
      message: compatibleAlternatives.length > 0
        ? `Domestication creates undersized fragments. Recommended: switch to ${bestAlternative.enzyme} (${bestAlternative.fullName}).`
        : `Domestication creates undersized fragments. ${bestAlternative.enzyme} has fewest sites (${bestAlternative.internalSites}).`,
    };
  }

  // Check if any domestication options lack valid junctions
  const invalidOptions = domesticationOptions.filter((opt: any) => !opt.hasValidOption);
  if (invalidOptions.length > 0) {
    const alternatives = recommendAlternativeEnzymes(sequence, enzyme);
    const compatibleAlternatives = alternatives.filter(a => a.isCompatible);
    const bestAlternative = compatibleAlternatives[0] || alternatives[0];

    return {
      needsDomestication: true,
      sites: internalSites.sites,
      domesticationOptions,
      additionalFragments: internalSites.count,
      error: {
        type: 'NO_VALID_JUNCTION',
        message: `${invalidOptions.length} site(s) have no valid junction options (all candidates would recreate the site or have low quality).`,
        details: invalidOptions.map((opt: any) => ({
          position: opt.sitePosition,
          sequence: opt.siteSequence,
          candidateCount: opt.candidates?.length || 0,
        })),
      },
      alternativeEnzymes: alternatives,
      recommendedEnzyme: bestAlternative,
      message: compatibleAlternatives.length > 0
        ? `Some sites cannot be domesticated. Recommended: switch to ${bestAlternative.enzyme} (${bestAlternative.fullName}).`
        : `Some sites cannot be domesticated. Consider manual sequence modification.`,
    };
  }

  return {
    needsDomestication: true,
    sites: internalSites.sites,
    domesticationOptions,
    additionalFragments: internalSites.count,
    error: null,
    alternativeEnzymes: null,
    message: `Found ${internalSites.count} internal ${enzyme} site(s). ` +
             `Assembly will require ${internalSites.count} additional fragment(s) for domestication.`,
  };
}

/**
 * Find optimal junction positions within an internal site to break it
 *
 * The key is to place the junction so that:
 * 1. The recognition sequence is split (broken)
 * 2. The resulting 4bp overhang is high-quality (good fidelity)
 * 3. The assembled product doesn't recreate the site
 */
export function findDomesticationJunctions(
  sequence: string,
  site: InternalSite,
  enzyme: string = 'BsaI'
): DomesticationJunctionOptions {
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enz) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }
  const recognition = enz.recognition; // e.g., 'GGTCTC' for BsaI
  const cutOffset = enz.cutOffset || 1; // Distance from recognition to cut site
  const overhangLen = enz.overhangLength || 4; // Usually 4bp for Type IIS enzymes

  const siteStart = site.position;
  const siteEnd = siteStart + recognition.length;
  const seq = sequence.toUpperCase();

  // Find all possible junction positions within and around the site
  const candidates: JunctionCandidate[] = [];

  // Calculate valid junction range
  const minJunctionPos = Math.max(0, siteStart - overhangLen + 1);
  const maxJunctionPos = Math.min(seq.length - overhangLen, siteEnd - 1);

  for (let pos = minJunctionPos; pos <= maxJunctionPos; pos++) {
    const overhang = seq.slice(pos, pos + overhangLen);

    // Check if this junction position actually breaks the site
    const breaksForward = doesJunctionBreakSite(pos, overhangLen, siteStart, siteEnd, recognition, false);
    const breaksReverse = doesJunctionBreakSite(pos, overhangLen, siteStart, siteEnd, recognition, true);

    if (breaksForward || breaksReverse) {
      // Score this overhang
      const quality = scoreOverhangForDomestication(overhang, enzyme, seq, pos);

      // Check if assembly would recreate the site
      const recreatesSite = wouldRecreateRecognitionSite(seq, pos, overhangLen, recognition, enzyme);

      candidates.push({
        position: pos,
        overhang,
        quality,
        breaksForward,
        breaksReverse,
        recreatesSite,
        distanceFromSiteCenter: Math.abs(pos + overhangLen/2 - (siteStart + recognition.length/2)),
        isValid: !recreatesSite && quality.score >= 50,
      });
    }
  }

  // Sort by quality (valid first, then by score)
  candidates.sort((a, b) => {
    if (a.isValid !== b.isValid) return Number(b.isValid) - Number(a.isValid);
    return b.quality.score - a.quality.score;
  });

  // Get the best valid candidate
  const bestCandidate = candidates.find(c => c.isValid) || candidates[0];

  return {
    site,
    sitePosition: siteStart,
    siteSequence: site.sequence,
    orientation: site.orientation,
    candidates,
    recommended: bestCandidate,
    hasValidOption: candidates.some(c => c.isValid),
  };
}

/**
 * Check if a junction at given position breaks a recognition site
 */
function doesJunctionBreakSite(
  junctionPos: number,
  overhangLen: number,
  siteStart: number,
  siteEnd: number,
  recognition: string,
  isReverse: boolean
): boolean {
  // The overhang spans [junctionPos, junctionPos + overhangLen)
  // The site spans [siteStart, siteEnd)
  // Junction breaks site if overhang overlaps with site interior

  const overhangStart = junctionPos;
  const overhangEnd = junctionPos + overhangLen;

  // Check if overhang overlaps with the recognition site
  const overlapStart = Math.max(overhangStart, siteStart);
  const overlapEnd = Math.min(overhangEnd, siteEnd);

  // Must have actual overlap (at least 1bp inside the site)
  if (overlapStart >= overlapEnd) {
    return false;
  }

  // For a site to be truly broken, the overhang should split the recognition
  // sequence, not just touch its edges. Check that overlap is in the interior.
  const isInteriorOverlap = overlapStart > siteStart || overlapEnd < siteEnd;

  return isInteriorOverlap;
}

/**
 * Score an overhang for use in domestication
 *
 * Uses the unified overhang-validation module for consistent scoring.
 */
function scoreOverhangForDomestication(
  overhang: string,
  enzyme: string,
  sequence: string,
  position: number
): OverhangQuality {
  // Use the unified validation module
  const validation = validateOverhang(overhang, {
    enzyme,
  });

  // Convert unified validation to local OverhangQuality format
  const issues: string[] = validation.issues
    .filter(i => i.severity === 'error' || i.severity === 'warning')
    .map(i => i.message);

  const benefits: string[] = [];

  // Check for TNNA pattern benefit (not flagged as issue in unified module)
  if (/^T..A$/.test(overhang.toUpperCase())) {
    benefits.push('TNNA pattern - high ligation efficiency');
  }

  // gcContent is 0-1, convert to percent for comparison
  const gcPercent = validation.details.gcContent * 100;

  // Check for optimal GC content (45-55%)
  if (gcPercent >= 45 && gcPercent <= 55) {
    benefits.push('Optimal GC content');
  }

  return {
    score: validation.score,
    issues,
    benefits,
    gc: gcPercent,
    isPalindrome: validation.details.isPalindrome,
  };
}

/**
 * Check if assembling with this junction would recreate the recognition site
 *
 * Uses the unified checkSiteRecreation function from overhang-validation.ts
 * for comprehensive site recreation analysis.
 */
function wouldRecreateRecognitionSite(
  sequence: string,
  junctionPos: number,
  overhangLen: number,
  recognition: string,
  enzyme: string
): boolean {
  const seq = sequence.toUpperCase();

  // Get the overhang
  const overhang = seq.slice(junctionPos, junctionPos + overhangLen);

  // Get flanking sequences for comprehensive analysis
  const flankingLen = 10; // Check 10bp on each side
  const upstreamSeq = seq.slice(Math.max(0, junctionPos - flankingLen), junctionPos);
  const downstreamSeq = seq.slice(junctionPos + overhangLen, junctionPos + overhangLen + flankingLen);

  // Use the unified site recreation check
  const result = checkSiteRecreation(upstreamSeq, downstreamSeq, overhang, enzyme);

  // Return true if there's any high risk of site recreation
  return result.recreatesSite || result.risk === 'high';
}

/**
 * Integrate domestication junctions into fusion site optimization
 */
export function optimizeWithDomestication(
  sequence: string,
  requestedFragments: number,
  enzyme: string = 'BsaI',
  options: {
    weights?: any;
    constraints?: any;
    bioContext?: any;
  } = {}
): DomesticationOptimizationResult {
  const {
    weights = {},
    constraints = {},
    bioContext = {},
  } = options;

  // Step 1: Analyze for domestication needs
  const domesticationAnalysis = analyzeForDomestication(sequence, enzyme);

  // Step 2: Calculate total fragments needed
  const domesticationJunctions = domesticationAnalysis.needsDomestication
    ? (domesticationAnalysis.domesticationOptions || []).filter(opt => opt.hasValidOption)
    : [];

  const additionalFragments = domesticationJunctions.length;
  const totalFragments = requestedFragments + additionalFragments;
  const totalJunctions = totalFragments - 1;

  // Step 3: Create required junction positions from domestication
  const requiredJunctions: RequiredJunction[] = domesticationJunctions.map(opt => ({
    position: opt.recommended!.position,
    overhang: opt.recommended!.overhang,
    isRequired: true,
    reason: 'domestication',
    site: opt.site,
    quality: opt.recommended!.quality,
  }));

  // Step 4: Calculate remaining junctions needed
  const remainingJunctions = Math.max(0, (requestedFragments - 1) - requiredJunctions.length);

  // Step 5: Generate target positions for remaining junctions
  const excludeZones: ExcludeZone[] = requiredJunctions.map(j => ({
    start: j.position - 200,
    end: j.position + 200,
  }));

  const targetPositions = generateUserJunctionTargets(
    sequence.length,
    remainingJunctions,
    excludeZones,
    constraints
  );

  return {
    needsDomestication: domesticationAnalysis.needsDomestication,
    internalSites: domesticationAnalysis.sites,
    domesticationJunctions: requiredJunctions,
    additionalFragmentsNeeded: additionalFragments,
    totalFragments,
    totalJunctions: requiredJunctions.length + remainingJunctions,
    remainingJunctionsToOptimize: remainingJunctions,
    targetPositions,
    message: domesticationAnalysis.needsDomestication
      ? `Auto-domestication: Adding ${additionalFragments} junction(s) to remove internal ${enzyme} sites. ` +
        `Total assembly: ${totalFragments} fragments.`
      : `No internal sites found. Proceeding with ${requestedFragments} fragments.`,
    constraints: {
      ...constraints,
      excludeZones,
    },
  };
}

/**
 * Generate target positions for user-requested junctions,
 * avoiding exclusion zones (domestication junctions)
 */
function generateUserJunctionTargets(
  seqLength: number,
  numJunctions: number,
  excludeZones: ExcludeZone[],
  constraints: any = {}
): JunctionTarget[] {
  if (numJunctions === 0) return [];

  const {
    minFragmentSize = 200,
    maxFragmentSize = 5000,
    minDistanceFromEnds = 50,
  } = constraints;

  // Calculate ideal spacing
  const usableLength = seqLength - (2 * minDistanceFromEnds);
  const idealSpacing = usableLength / (numJunctions + 1);

  const targets: JunctionTarget[] = [];

  for (let i = 1; i <= numJunctions; i++) {
    const idealPos = minDistanceFromEnds + Math.round(idealSpacing * i);

    // Check if in exclusion zone
    const inExclusionZone = excludeZones.some(
      zone => idealPos >= zone.start && idealPos <= zone.end
    );

    // If in exclusion zone, find nearest valid position
    let adjustedPos = idealPos;
    if (inExclusionZone) {
      // Find nearest edge of exclusion zone
      const nearestZone = excludeZones.find(
        zone => idealPos >= zone.start && idealPos <= zone.end
      );
      if (nearestZone) {
        const distToStart = Math.abs(idealPos - nearestZone.start);
        const distToEnd = Math.abs(idealPos - nearestZone.end);
        adjustedPos = distToStart < distToEnd ? nearestZone.start - 50 : nearestZone.end + 50;
      }
    }

    targets.push({
      idealPosition: adjustedPos,
      searchRadius: 50,
      index: i - 1,
    });
  }

  return targets;
}

/**
 * Get a summary of domestication requirements for UI display
 */
export function getDomesticationSummary(sequence: string, enzyme: string = 'BsaI') {
  const analysis = analyzeForDomestication(sequence, enzyme);

  if (!analysis.needsDomestication) {
    return {
      status: 'compatible',
      icon: '✓',
      title: 'No Internal Sites',
      description: `Sequence is compatible with ${enzyme} assembly`,
      sites: [],
      action: null,
    };
  }

  const siteDescriptions = analysis.sites.map((site, i) => {
    const option = analysis.domesticationOptions ? analysis.domesticationOptions[i] : undefined;
    return {
      position: site.position,
      sequence: site.sequence,
      orientation: site.orientation,
      recommendedJunction: option?.recommended,
      hasValidOption: option?.hasValidOption,
    };
  });

  const allValid = siteDescriptions.every(s => s.hasValidOption);

  return {
    status: allValid ? 'auto-fixable' : 'needs-attention',
    icon: allValid ? '⚡' : '⚠',
    title: `${analysis.sites.length} Internal Site${analysis.sites.length > 1 ? 's' : ''} Found`,
    description: allValid
      ? `Will add ${analysis.additionalFragments} junction(s) to auto-domesticate`
      : `Some sites may require manual adjustment or alternative enzyme`,
    sites: siteDescriptions,
    additionalFragments: analysis.additionalFragments,
    action: allValid ? 'auto-domesticate' : 'review-required',
  };
}

/**
 * Validate that a set of overhangs (including domestication overhangs)
 * won't have cross-ligation issues
 */
export function validateDomesticationOverhangs(
  domesticationOverhangs: string[],
  otherOverhangs: string[],
  enzyme: string
): OverhangValidation {
  const allOverhangs = [...domesticationOverhangs, ...otherOverhangs];

  // Check for duplicates
  const seen = new Set<string>();
  const duplicates: string[] = [];
  for (const oh of allOverhangs) {
    if (seen.has(oh)) {
      duplicates.push(oh);
    }
    seen.add(oh);
    // Also check reverse complement
    const rc = reverseComplement(oh);
    if (seen.has(rc)) {
      duplicates.push(`${oh}/${rc}`);
    }
  }

  // Calculate set fidelity
  const fidelityResult = calculateExperimentalFidelity(allOverhangs, enzyme);

  return {
    isValid: duplicates.length === 0 && fidelityResult.assemblyFidelity >= 0.85,
    duplicates,
    setFidelity: fidelityResult.assemblyFidelity,
    problematicPairs: (fidelityResult as any).problematicPairs || [],
    recommendation: fidelityResult.assemblyFidelity >= 0.95
      ? 'Excellent overhang set'
      : fidelityResult.assemblyFidelity >= 0.85
      ? 'Good overhang set'
      : 'Consider adjusting junction positions for better fidelity',
  };
}

/**
 * Optimize the global overhang set by selecting the best combination of
 * domestication junction candidates that maximizes overall assembly fidelity.
 */
export function optimizeGlobalOverhangSet(
  domesticationOptions: DomesticationJunctionOptions[] | null,
  existingOverhangs: string[] = [],
  enzyme: string = 'BsaI',
  options: { maxCombinations?: number; minFidelity?: number } = {}
): GlobalOverhangOptimizationResult {
  const {
    maxCombinations = 10000,
    minFidelity = 0.85,
  } = options;

  // If no domestication needed, just validate existing overhangs
  if (!domesticationOptions || domesticationOptions.length === 0) {
    const validation = validateDomesticationOverhangs([], existingOverhangs, enzyme);
    return {
      success: true,
      selections: [],
      overhangs: existingOverhangs,
      fidelity: validation.setFidelity,
      validation,
      searchMethod: 'exhaustive',
      combinationsSearched: 0,
    };
  }

  // Get valid candidates for each site
  const siteCandidates = domesticationOptions.map(opt => {
    const valid = (opt.candidates || []).filter(c => c.isValid);
    // If no valid candidates, include best invalid as fallback
    return valid.length > 0 ? valid : (opt.candidates || []).slice(0, 3);
  });

  // Check if any site has no candidates at all
  const emptySites = siteCandidates.filter(c => c.length === 0);
  if (emptySites.length > 0) {
    return {
      success: false,
      selections: [],
      overhangs: [],
      fidelity: 0,
      validation: { isValid: false, duplicates: [], setFidelity: 0, problematicPairs: [], recommendation: '' },
      searchMethod: 'exhaustive',
      combinationsSearched: 0,
      error: 'NO_CANDIDATES',
      message: `${emptySites.length} site(s) have no junction candidates`,
    };
  }

  // Calculate total combinations
  const totalCombinations = siteCandidates.reduce((acc, c) => acc * c.length, 1);

  let bestSelection: DomesticationSelection[] | null = null;
  let bestFidelity = 0;
  let bestValidation: OverhangValidation | null = null;

  if (totalCombinations <= maxCombinations) {
    // Exhaustive search for small search spaces
    const result = searchAllCombinations(siteCandidates, existingOverhangs, enzyme);
    bestSelection = result.selection;
    bestFidelity = result.fidelity;
    bestValidation = result.validation;
  } else {
    // Greedy optimization with local search for large spaces
    const result = greedyOptimizeOverhangs(siteCandidates, existingOverhangs, enzyme);
    bestSelection = result.selection;
    bestFidelity = result.fidelity;
    bestValidation = result.validation;
  }

  // Build result
  const selectedOverhangs = bestSelection.map(s => s.overhang);
  const allOverhangs = [...selectedOverhangs, ...existingOverhangs];

  return {
    success: bestFidelity >= minFidelity,
    selections: bestSelection,
    overhangs: allOverhangs,
    fidelity: bestFidelity,
    validation: bestValidation!,
    searchMethod: totalCombinations <= maxCombinations ? 'exhaustive' : 'greedy',
    combinationsSearched: Math.min(totalCombinations, maxCombinations),
  };
}

/**
 * Exhaustive search through all combinations (for small search spaces)
 */
function searchAllCombinations(
  siteCandidates: JunctionCandidate[][],
  existingOverhangs: string[],
  enzyme: string
): { selection: DomesticationSelection[]; fidelity: number; validation: OverhangValidation } {
  const indices = siteCandidates.map(() => 0);
  let bestSelection: DomesticationSelection[] | null = null;
  let bestFidelity = 0;
  let bestValidation: OverhangValidation | null = null;

  const evaluate = () => {
    const selection: DomesticationSelection[] = indices.map((idx, siteIdx) => {
      const candidate = siteCandidates[siteIdx][idx];
      return {
        overhang: candidate.overhang,
        position: candidate.position,
        quality: candidate.quality,
      };
    });
    const overhangs = selection.map(s => s.overhang);
    const validation = validateDomesticationOverhangs(overhangs, existingOverhangs, enzyme);

    if (validation.setFidelity > bestFidelity) {
      bestFidelity = validation.setFidelity;
      bestSelection = selection;
      bestValidation = validation;
    }
  };

  // Iterate through all combinations
  evaluate();
  while (true) {
    // Increment indices (like counting in mixed-radix)
    let pos = 0;
    while (pos < indices.length) {
      indices[pos]++;
      if (indices[pos] < siteCandidates[pos].length) {
        break;
      }
      indices[pos] = 0;
      pos++;
    }
    if (pos >= indices.length) break; // All combinations exhausted

    evaluate();
  }

  return { selection: bestSelection!, fidelity: bestFidelity, validation: bestValidation! };
}

/**
 * Greedy optimization with local improvements (for large search spaces)
 */
function greedyOptimizeOverhangs(
  siteCandidates: JunctionCandidate[][],
  existingOverhangs: string[],
  enzyme: string
): { selection: DomesticationSelection[]; fidelity: number; validation: OverhangValidation } {
  // Start with highest-scoring candidate for each site
  let selection: DomesticationSelection[] = siteCandidates.map(candidates => {
    const best = candidates[0];
    return {
      overhang: best.overhang,
      position: best.position,
      quality: best.quality,
    };
  });
  let overhangs = selection.map(s => s.overhang);
  let validation = validateDomesticationOverhangs(overhangs, existingOverhangs, enzyme);
  let bestFidelity = validation.setFidelity;

  // Iteratively improve each position
  let improved = true;
  let iterations = 0;
  const maxIterations = 100;

  while (improved && iterations < maxIterations) {
    improved = false;
    iterations++;

    for (let siteIdx = 0; siteIdx < siteCandidates.length; siteIdx++) {
      const currentCandidate = selection[siteIdx];

      // Try each alternative candidate for this site
      for (const candidate of siteCandidates[siteIdx]) {
        if (candidate.overhang === currentCandidate.overhang &&
            candidate.position === currentCandidate.position) continue;

        // Create new selection with this candidate
        const newSelection = [...selection];
        newSelection[siteIdx] = {
          overhang: candidate.overhang,
          position: candidate.position,
          quality: candidate.quality,
        };

        const newOverhangs = newSelection.map(s => s.overhang);
        const newValidation = validateDomesticationOverhangs(newOverhangs, existingOverhangs, enzyme);

        if (newValidation.setFidelity > bestFidelity) {
          selection = newSelection;
          overhangs = newOverhangs;
          validation = newValidation;
          bestFidelity = newValidation.setFidelity;
          improved = true;
        }
      }
    }
  }

  return { selection, fidelity: bestFidelity, validation };
}

/**
 * Merge adjacent internal sites that are close together into a single
 * domestication group.
 */
export function mergeAdjacentSites(
  sites: InternalSite[] | null,
  sequence: string,
  options: { mergeDistance?: number; minFragmentSize?: number } = {}
): MergedSitesResult {
  const {
    mergeDistance = 100,
    minFragmentSize = 50,
  } = options;

  if (!sites || sites.length === 0) {
    return { groups: [], needsMerging: false };
  }

  if (sites.length === 1) {
    return {
      groups: [{ sites: [sites[0]], type: 'single', canHandle: true, strategy: 'Place junction within site to break it' }],
      needsMerging: false,
    };
  }

  // Sort sites by position
  const sortedSites = [...sites].sort((a, b) => a.position - b.position);
  const groups: SiteGroup[] = [];
  let currentGroup = [sortedSites[0]];

  for (let i = 1; i < sortedSites.length; i++) {
    const prevSite = sortedSites[i - 1];
    const currSite = sortedSites[i];
    const distance = currSite.position - prevSite.position;

    if (distance < mergeDistance) {
      // Sites are close, add to current group
      currentGroup.push(currSite);
    } else {
      // Sites are far apart, start new group
      groups.push(analyzeGroup(currentGroup, sequence, minFragmentSize));
      currentGroup = [currSite];
    }
  }

  // Don't forget the last group
  groups.push(analyzeGroup(currentGroup, sequence, minFragmentSize));

  const needsMerging = groups.some(g => g.sites.length > 1);
  const allCanHandle = groups.every(g => g.canHandle);

  return {
    groups,
    needsMerging,
    allCanHandle,
    totalSites: sites.length,
    groupCount: groups.length,
    strategy: allCanHandle
      ? (needsMerging ? 'MERGE_ADJACENT' : 'STANDARD')
      : 'REQUIRES_ALTERNATIVE_ENZYME',
  };
}

/**
 * Analyze a group of adjacent sites to determine handling strategy
 */
function analyzeGroup(sites: InternalSite[], sequence: string, minFragmentSize: number): SiteGroup {
  if (sites.length === 1) {
    return {
      sites,
      type: 'single',
      canHandle: true,
      strategy: 'Place junction within site to break it',
    };
  }

  // Multiple adjacent sites
  const firstPos = sites[0].position;
  const lastPos = sites[sites.length - 1].position;
  const span = lastPos - firstPos + 6; // +6 for recognition length

  // Check if the span is too small for multiple fragments
  const minRequiredSpace = (sites.length) * minFragmentSize;

  if (span < minRequiredSpace) {
    return {
      sites,
      type: 'clustered',
      span,
      canHandle: false,
      strategy: 'Sites too close - consider alternative enzyme',
      reason: `${sites.length} sites within ${span}bp span require at least ${minRequiredSpace}bp`,
    };
  }

  // Sites can potentially be handled with careful junction placement
  return {
    sites,
    type: 'adjacent',
    span,
    canHandle: true,
    strategy: 'Place junctions to create fragments of at least ' + minFragmentSize + 'bp',
    fragmentCount: sites.length + 1,
  };
}

/**
 * Post-domestication validation - verify that the domesticated assembly is valid
 */
export function validatePostDomestication(
  fragments: FragmentWithOverhang[],
  enzyme: string = 'BsaI',
  options: { minFragmentSize?: number; expectedSequence?: string | null; isCircular?: boolean } = {}
): PostDomesticationValidation {
  const {
    minFragmentSize = DOMESTICATION_DEFAULTS.minFragmentSize,
    expectedSequence = null,
    isCircular = false,
  } = options;

  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enz) {
    return {
      isValid: false,
      issues: [{ type: 'error', code: 'UNKNOWN_ENZYME', message: `Unknown enzyme: ${enzyme}` }],
      warnings: [],
      passed: [],
      summary: { fragmentCount: 0, totalIssues: 1, totalWarnings: 0, totalPassed: 0, status: 'invalid' },
    };
  }

  const issues: ValidationIssue[] = [];
  const warnings: ValidationIssue[] = [];
  const passed: ValidationCheck[] = [];
  const recognition = enz.recognition;
  const recognitionRC = reverseComplement(recognition);
  const overhangLen = enz.overhangLength || 4;

  // Validation 1: Check fragment sizes
  const fragmentSizeIssues = [];
  for (let i = 0; i < fragments.length; i++) {
    const frag = fragments[i];
    const fragLen = frag.seq?.length || 0;
    if (fragLen < minFragmentSize) {
      fragmentSizeIssues.push({
        fragmentIndex: i,
        fragmentId: frag.id || `Fragment ${i + 1}`,
        length: fragLen,
        minRequired: minFragmentSize,
      });
    }
  }

  if (fragmentSizeIssues.length > 0) {
    issues.push({
      type: 'error',
      code: 'FRAGMENT_TOO_SHORT',
      message: `${fragmentSizeIssues.length} fragment(s) below minimum size (${minFragmentSize}bp)`,
      details: fragmentSizeIssues,
    });
  } else {
    passed.push({
      code: 'FRAGMENT_SIZES',
      message: `All ${fragments.length} fragments meet minimum size requirement`,
    });
  }

  // Validation 2: Check for recognition sites in fragments
  const sitesInFragments = [];
  for (let i = 0; i < fragments.length; i++) {
    const frag = fragments[i];
    const seq = (frag.seq || '').toUpperCase();

    // Check for internal sites (not at ends where primers will add new sites)
    const searchRegion = seq.slice(overhangLen, seq.length - overhangLen);
    if (searchRegion.includes(recognition) || searchRegion.includes(recognitionRC)) {
      sitesInFragments.push({
        fragmentIndex: i,
        fragmentId: frag.id || `Fragment ${i + 1}`,
        hasForward: searchRegion.includes(recognition),
        hasReverse: searchRegion.includes(recognitionRC),
      });
    }
  }

  if (sitesInFragments.length > 0) {
    issues.push({
      type: 'error',
      code: 'INTERNAL_SITES_REMAIN',
      message: `${sitesInFragments.length} fragment(s) still contain internal ${enzyme} sites`,
      details: sitesInFragments,
    });
  } else {
    passed.push({
      code: 'NO_INTERNAL_SITES',
      message: `No internal ${enzyme} sites in any fragment`,
    });
  }

  // Validation 3: Check overhang uniqueness and fidelity
  const overhangs = fragments
    .filter(f => f._overhang)
    .map(f => f._overhang!.toUpperCase());

  if (overhangs.length > 0) {
    // Check for duplicates
    const seen = new Set<string>();
    const duplicates: string[] = [];
    for (const oh of overhangs) {
      const rc = reverseComplement(oh);
      if (seen.has(oh)) {
        duplicates.push(oh);
      }
      if (seen.has(rc) && oh !== rc) {
        duplicates.push(`${oh} (reverse complement of existing overhang)`);
      }
      seen.add(oh);
    }

    if (duplicates.length > 0) {
      issues.push({
        type: 'error',
        code: 'DUPLICATE_OVERHANGS',
        message: `Duplicate overhangs detected: ${duplicates.join(', ')}`,
        details: duplicates,
      });
    } else {
      passed.push({
        code: 'UNIQUE_OVERHANGS',
        message: `All ${overhangs.length} overhangs are unique`,
      });
    }

    // Check fidelity if we have ligation data
    try {
      const fidelityResult = calculateExperimentalFidelity(overhangs, enzyme);
      if (fidelityResult.assemblyFidelity < 0.85) {
        warnings.push({
          type: 'warning',
          code: 'LOW_FIDELITY',
          message: `Assembly fidelity is ${(fidelityResult.assemblyFidelity * 100).toFixed(1)}% (recommended: >85%)`,
          details: {
            fidelity: fidelityResult.assemblyFidelity,
            problematicPairs: (fidelityResult as any).problematicPairs || [],
          },
        });
      } else if (fidelityResult.assemblyFidelity < 0.95) {
        passed.push({
          code: 'FIDELITY_GOOD',
          message: `Assembly fidelity: ${(fidelityResult.assemblyFidelity * 100).toFixed(1)}% (good)`,
        });
      } else {
        passed.push({
          code: 'FIDELITY_EXCELLENT',
          message: `Assembly fidelity: ${(fidelityResult.assemblyFidelity * 100).toFixed(1)}% (excellent)`,
        });
      }
    } catch (e) {
      // Fidelity calculation may not be available for all enzymes
      warnings.push({
        type: 'warning',
        code: 'FIDELITY_UNKNOWN',
        message: `Could not calculate fidelity for ${enzyme}`,
      });
    }
  }

  // Validation 4: Simulate junction scars for site recreation
  const scarIssues = [];
  for (let i = 0; i < fragments.length - 1; i++) {
    const upstreamFrag = fragments[i];
    const downstreamFrag = fragments[i + 1];
    const overhang = upstreamFrag._overhang;

    if (!overhang) continue;

    // Get the context around the junction
    const upstreamSeq = (upstreamFrag.seq || '').toUpperCase();
    const downstreamSeq = (downstreamFrag.seq || '').toUpperCase();

    // The scar region includes bases before and after the overhang
    const contextSize = recognition.length;
    const upstreamContext = upstreamSeq.slice(-contextSize);
    const overhangSeq = overhang.toUpperCase();
    const downstreamContext = downstreamSeq.slice(0, contextSize);

    // Simulate the assembled scar
    const scarRegion = upstreamContext + downstreamContext;

    if (scarRegion.includes(recognition) || scarRegion.includes(recognitionRC)) {
      scarIssues.push({
        junctionIndex: i,
        upstream: upstreamFrag.id || `Fragment ${i + 1}`,
        downstream: downstreamFrag.id || `Fragment ${i + 2}`,
        overhang: overhang,
        scarRegion: scarRegion,
        recreatedSite: scarRegion.includes(recognition) ? 'forward' : 'reverse',
      });
    }
  }

  // Check circular assembly junction (last to first)
  if (isCircular && fragments.length > 1) {
    const lastFrag = fragments[fragments.length - 1];
    const firstFrag = fragments[0];
    const lastOverhang = lastFrag._overhang;

    if (lastOverhang) {
      const lastSeq = (lastFrag.seq || '').toUpperCase();
      const firstSeq = (firstFrag.seq || '').toUpperCase();
      const contextSize = recognition.length;
      const upstreamContext = lastSeq.slice(-contextSize);
      const downstreamContext = firstSeq.slice(0, contextSize);
      const scarRegion = upstreamContext + downstreamContext;

      if (scarRegion.includes(recognition) || scarRegion.includes(recognitionRC)) {
        scarIssues.push({
          junctionIndex: 'circular' as any,
          upstream: lastFrag.id || `Fragment ${fragments.length}`,
          downstream: firstFrag.id || 'Fragment 1',
          overhang: lastOverhang,
          scarRegion: scarRegion,
          recreatedSite: scarRegion.includes(recognition) ? 'forward' : 'reverse',
        });
      }
    }
  }

  if (scarIssues.length > 0) {
    issues.push({
      type: 'error',
      code: 'SITE_RECREATED_AT_SCAR',
      message: `${scarIssues.length} junction scar(s) recreate ${enzyme} recognition site`,
      details: scarIssues,
      recommendation: 'Adjust junction positions to avoid recreating the recognition site at the assembly scar',
    });
  } else if (fragments.length > 1) {
    passed.push({
      code: 'NO_SCAR_SITES',
      message: `No ${enzyme} sites recreated at ${fragments.length - 1} junction scar(s)`,
    });
  }

  // Validation 5: Assembly simulation (if expected sequence provided)
  if (expectedSequence) {
    // Simulate assembly by concatenating fragments at overhangs
    let assembled = '';
    for (let i = 0; i < fragments.length; i++) {
      const frag = fragments[i];
      const seq = frag.seq || '';

      if (i === 0) {
        // First fragment: include everything up to overhang (inclusive)
        assembled = seq;
      } else {
        // Subsequent fragments: skip the overhang region (already included from previous)
        const prevOverhang = fragments[i - 1]._overhang || '';
        // Fragment starts with the overhang, skip it
        assembled += seq.slice(prevOverhang.length);
      }
    }

    const expectedUpper = expectedSequence.toUpperCase();
    const assembledUpper = assembled.toUpperCase();

    if (assembledUpper === expectedUpper) {
      passed.push({
        code: 'ASSEMBLY_VERIFIED',
        message: 'Assembled sequence matches expected sequence',
      });
    } else {
      const lengthDiff = assembled.length - expectedSequence.length;
      issues.push({
        type: 'error',
        code: 'ASSEMBLY_MISMATCH',
        message: `Assembled sequence (${assembled.length}bp) doesn't match expected (${expectedSequence.length}bp)`,
        details: { lengthDifference: lengthDiff },
      });
    }
  }

  const isValid = issues.length === 0;

  return {
    isValid,
    issues,
    warnings,
    passed,
    summary: {
      fragmentCount: fragments.length,
      totalIssues: issues.length,
      totalWarnings: warnings.length,
      totalPassed: passed.length,
      status: isValid
        ? (warnings.length > 0 ? 'valid_with_warnings' : 'valid')
        : 'invalid',
    },
  };
}
