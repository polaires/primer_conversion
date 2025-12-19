/**
 * Fusion Site Scanner for Golden Gate Assembly
 *
 * Scans DNA sequences to find optimal candidate positions for fusion sites.
 * Implements the candidate generator component of the Fusion Site Optimizer.
 *
 * Key features:
 * 1. Identifies all valid 4bp overhang positions
 * 2. Applies hard filters (palindromes, homopolymers)
 * 3. Respects search windows and forbidden regions
 * 4. Pre-computes useful metrics for scoring
 */

import { reverseComplement } from './enzymes.js';
import { GOLDEN_GATE_ENZYMES, getEnzymeLigationData, getOverhangFidelityExperimental } from './goldengate.js';
import { isPalindrome, isHomopolymer, calculateEfficiency, countGC } from './overhang-efficiency.js';
import { checkSiteCreation } from './site-creation-check.js';

/**
 * Search window region
 */
export interface SearchWindow {
  start: number;
  end: number;
}

/**
 * Forbidden region
 */
export interface ForbiddenRegion {
  start: number;
  end: number;
  reason?: string;
}

/**
 * Fusion site candidate
 */
export interface FusionSiteCandidate {
  position: number;
  overhang: string;
  reverseComplement: string;
  baseFidelity: number;
  gcContent: number;
  gcCount: number;
  isTNNA: boolean;
  isHighGC: boolean;
  isLowGC: boolean;
  upstreamContext: string;
  downstreamContext: string;
  efficiency?: number;
  efficiencyWarnings?: string[];
  isOptimalEfficiency?: boolean;
  isAcceptableEfficiency?: boolean;
  siteCreationRisk?: boolean;
  siteCreationSeverity?: string;
  siteCreationWarnings?: string[];
  combinedScore?: number;
  score?: any;
}

/**
 * Scanning options
 */
export interface ScanOptions {
  enzyme?: string;
  searchWindows?: SearchWindow[] | null;
  forbiddenRegions?: ForbiddenRegion[];
  minDistanceFromEnds?: number;
  minDistanceBetween?: number;
  includeEfficiency?: boolean;
  includeSiteCreation?: boolean;
  maxCandidates?: number;
}

/**
 * Check if a position is within any of the allowed windows
 * @param pos - Position to check
 * @param windows - Array of search windows
 * @returns True if in a window
 */
function isInWindow(pos: number, windows: SearchWindow[] | null): boolean {
  if (!windows || windows.length === 0) return true;
  return windows.some(w => pos >= w.start && pos <= w.end);
}

/**
 * Check if a position is in a forbidden region
 * @param pos - Position to check
 * @param regions - Array of forbidden regions
 * @returns True if in forbidden region
 */
function isInForbiddenRegion(pos: number, regions: ForbiddenRegion[]): boolean {
  if (!regions || regions.length === 0) return false;
  return regions.some(r => pos >= r.start && pos < r.end);
}

/**
 * Get base overhang fidelity from experimental matrix
 * Uses centralized function from goldengate.js
 * @param overhang - 4bp overhang sequence
 * @param enzyme - Enzyme name
 * @returns Base fidelity (0-1)
 */
function getOverhangFidelity(overhang: string, enzyme: string): number {
  const result = getOverhangFidelityExperimental(overhang, enzyme);
  return result.fidelity;
}

/**
 * Scan sequence for all potential fusion sites
 *
 * This is the core scanning function that identifies candidate junction
 * positions in a DNA sequence.
 *
 * @param sequence - DNA sequence to scan
 * @param options - Scanning options
 * @returns Array of candidate fusion sites
 */
export function scanForFusionSites(sequence: string, options: ScanOptions = {}): FusionSiteCandidate[] {
  const {
    enzyme = 'BsaI',
    searchWindows = null,
    forbiddenRegions = [],
    minDistanceFromEnds = 50,
    minDistanceBetween = 100,
    includeEfficiency = true,
    includeSiteCreation = false,
    maxCandidates = 1000,
  } = options;

  if (!sequence || sequence.length < minDistanceFromEnds * 2 + 4) {
    return [];
  }

  const candidates: FusionSiteCandidate[] = [];
  const seq = sequence.toUpperCase();

  // Determine enzyme overhang length
  const enzData = (GOLDEN_GATE_ENZYMES as any)[enzyme];
  const overhangLength = enzData?.overhangLength || 4;

  // Scan through valid positions
  const startPos = minDistanceFromEnds;
  const endPos = seq.length - minDistanceFromEnds - overhangLength;

  for (let pos = startPos; pos <= endPos; pos++) {
    // Check if position is in allowed window
    if (searchWindows && !isInWindow(pos, searchWindows)) continue;
    if (isInForbiddenRegion(pos, forbiddenRegions)) continue;

    const overhang = seq.slice(pos, pos + overhangLength);

    // Validate overhang contains only valid bases
    if (!/^[ATGC]+$/i.test(overhang)) continue;

    // Hard filters - exclude completely
    if (isPalindrome(overhang)) continue;
    if (isHomopolymer(overhang)) continue;

    // Get base fidelity from matrix
    const baseFidelity = getOverhangFidelity(overhang, enzyme);
    if (baseFidelity === 0) continue; // No ligation data

    // Build candidate object
    const candidate: FusionSiteCandidate = {
      position: pos,
      overhang,
      reverseComplement: reverseComplement(overhang),
      baseFidelity,

      // Pre-compute GC content
      gcContent: countGC(overhang) / overhangLength,
      gcCount: countGC(overhang),

      // Pre-compute pattern flags
      isTNNA: /^T..A$/i.test(overhang),
      isHighGC: countGC(overhang) === overhangLength,
      isLowGC: countGC(overhang) === 0,

      // Context information
      upstreamContext: seq.slice(Math.max(0, pos - 10), pos),
      downstreamContext: seq.slice(pos + overhangLength, Math.min(seq.length, pos + overhangLength + 10)),
    };

    // Optional: Calculate efficiency metrics
    if (includeEfficiency) {
      const efficiency = calculateEfficiency(overhang);
      candidate.efficiency = (efficiency as any).efficiency;
      candidate.efficiencyWarnings = (efficiency as any).warnings;
      candidate.isOptimalEfficiency = (efficiency as any).isOptimal;
      candidate.isAcceptableEfficiency = (efficiency as any).isAcceptable;
    }

    // Optional: Check for site creation (expensive)
    if (includeSiteCreation) {
      const siteCheck = checkSiteCreation(seq, pos, enzyme);
      candidate.siteCreationRisk = (siteCheck as any).hasRisk;
      candidate.siteCreationSeverity = (siteCheck as any).severity;
      if ((siteCheck as any).hasRisk) {
        candidate.siteCreationWarnings = (siteCheck as any).risks.map((r: any) => r.message);
      }
    }

    candidates.push(candidate);

    // Limit total candidates
    if (candidates.length >= maxCandidates) break;
  }

  return candidates;
}

/**
 * Ranking options
 */
interface RankOptions extends ScanOptions {
  topN?: number;
  sortBy?: 'fidelity' | 'efficiency' | 'combined';
}

/**
 * Statistics for scanned candidates
 */
export interface ScanStatistics {
  totalCandidates: number;
  averageFidelity: number;
  averageEfficiency: number;
  tnnaCount: number;
  highGCCount: number;
  lowGCCount: number;
  optimalCount: number;
  acceptableCount: number;
}

/**
 * Scan and rank result
 */
export interface ScanAndRankResult {
  candidates: FusionSiteCandidate[];
  allCandidates: FusionSiteCandidate[];
  stats: ScanStatistics;
  sequenceLength: number;
  coverage: number;
}

/**
 * Scan and rank fusion sites by quality
 *
 * @param sequence - DNA sequence to scan
 * @param options - Scanning options
 * @returns Ranked candidates and statistics
 */
export function scanAndRankFusionSites(sequence: string, options: RankOptions = {}): ScanAndRankResult {
  const {
    topN = 100,
    sortBy = 'fidelity',
    ...scanOptions
  } = options;

  // Get all candidates
  const candidates = scanForFusionSites(sequence, {
    ...scanOptions,
    includeEfficiency: true,
  });

  // Calculate combined score for each candidate
  const scored = candidates.map(c => ({
    ...c,
    combinedScore: c.baseFidelity * (c.efficiency || 1.0),
  }));

  // Sort based on criteria
  let sorted: FusionSiteCandidate[];
  switch (sortBy) {
    case 'efficiency':
      sorted = scored.sort((a, b) => (b.efficiency || 0) - (a.efficiency || 0));
      break;
    case 'combined':
      sorted = scored.sort((a, b) => (b.combinedScore || 0) - (a.combinedScore || 0));
      break;
    case 'fidelity':
    default:
      sorted = scored.sort((a, b) => b.baseFidelity - a.baseFidelity);
  }

  // Statistics
  const stats: ScanStatistics = {
    totalCandidates: candidates.length,
    averageFidelity: candidates.reduce((s, c) => s + c.baseFidelity, 0) / candidates.length || 0,
    averageEfficiency: candidates.reduce((s, c) => s + (c.efficiency || 1), 0) / candidates.length || 0,
    tnnaCount: candidates.filter(c => c.isTNNA).length,
    highGCCount: candidates.filter(c => c.isHighGC).length,
    lowGCCount: candidates.filter(c => c.isLowGC).length,
    optimalCount: candidates.filter(c => c.isOptimalEfficiency).length,
    acceptableCount: candidates.filter(c => c.isAcceptableEfficiency && !c.isOptimalEfficiency).length,
  };

  return {
    candidates: sorted.slice(0, topN),
    allCandidates: sorted,
    stats,
    sequenceLength: sequence.length,
    coverage: candidates.length / (sequence.length - 100), // Approximate coverage
  };
}

/**
 * Target region with weight
 */
export interface TargetRegion {
  start: number;
  end: number;
  weight?: number;
}

/**
 * Target region result
 */
interface TargetRegionResult {
  regionIndex: number;
  region: TargetRegion;
  candidates: FusionSiteCandidate[];
  totalFound: number;
  bestCandidate: FusionSiteCandidate | null;
}

/**
 * Target regions scan result
 */
export interface TargetRegionsScanResult {
  regions: TargetRegionResult[];
  totalRegions: number;
  regionsWithCandidates: number;
  emptyRegions: number[];
  allCandidates: FusionSiteCandidate[];
}

/**
 * Target scanning options
 */
interface TargetScanOptions extends ScanOptions {
  candidatesPerRegion?: number;
}

/**
 * Find candidates in specific target regions
 *
 * @param sequence - DNA sequence
 * @param targetRegions - Array of target regions
 * @param options - Scanning options
 * @returns Candidates organized by region
 */
export function scanTargetRegions(
  sequence: string,
  targetRegions: TargetRegion[],
  options: TargetScanOptions = {}
): TargetRegionsScanResult {
  const {
    candidatesPerRegion = 10,
    enzyme = 'BsaI',
    ...scanOptions
  } = options;

  const results: TargetRegionResult[] = targetRegions.map((region, index) => {
    const regionCandidates = scanForFusionSites(sequence, {
      ...scanOptions,
      enzyme,
      searchWindows: [{ start: region.start, end: region.end }],
      includeEfficiency: true,
    });

    // Sort by combined quality
    const sorted = regionCandidates
      .map(c => ({
        ...c,
        combinedScore: c.baseFidelity * (c.efficiency || 1.0),
      }))
      .sort((a, b) => (b.combinedScore || 0) - (a.combinedScore || 0));

    return {
      regionIndex: index,
      region,
      candidates: sorted.slice(0, candidatesPerRegion),
      totalFound: regionCandidates.length,
      bestCandidate: sorted[0] || null,
    };
  });

  // Check for regions with no candidates
  const emptyRegions = results.filter(r => r.totalFound === 0);

  return {
    regions: results,
    totalRegions: targetRegions.length,
    regionsWithCandidates: results.filter(r => r.totalFound > 0).length,
    emptyRegions: emptyRegions.map(r => r.regionIndex),
    allCandidates: results.flatMap(r => r.candidates),
  };
}

/**
 * Target position
 */
export interface TargetPosition {
  index: number;
  idealPosition: number;
  searchRegion: SearchWindow;
  expectedFragmentSize: number;
}

/**
 * Target position options
 */
interface TargetPositionOptions {
  minDistanceFromEnds?: number;
  searchRadius?: number;
}

/**
 * Generate evenly-spaced junction targets for a sequence
 *
 * @param sequenceLength - Length of the sequence
 * @param numFragments - Number of fragments desired
 * @param options - Options
 * @returns Target positions and search regions
 */
export function generateTargetPositions(
  sequenceLength: number,
  numFragments: number,
  options: TargetPositionOptions = {}
): TargetPosition[] {
  const {
    minDistanceFromEnds = 50,
    searchRadius = 50,
  } = options;

  const numJunctions = numFragments - 1;
  if (numJunctions <= 0) {
    return [];
  }

  // Calculate ideal fragment size
  const usableLength = sequenceLength - 2 * minDistanceFromEnds;
  const idealFragmentSize = usableLength / numFragments;

  const targets: TargetPosition[] = [];
  for (let i = 1; i <= numJunctions; i++) {
    const idealPosition = minDistanceFromEnds + Math.round(idealFragmentSize * i);

    targets.push({
      index: i - 1,
      idealPosition,
      searchRegion: {
        start: Math.max(minDistanceFromEnds, idealPosition - searchRadius),
        end: Math.min(sequenceLength - minDistanceFromEnds - 4, idealPosition + searchRadius),
      },
      expectedFragmentSize: idealFragmentSize,
    });
  }

  return targets;
}

/**
 * Filter candidates to ensure minimum distance between selections
 *
 * @param candidates - Array of candidate objects
 * @param minDistance - Minimum distance between selections
 * @param numToSelect - Number of candidates to select
 * @returns Filtered candidates
 */
export function filterByDistance(
  candidates: FusionSiteCandidate[],
  minDistance: number,
  numToSelect: number
): FusionSiteCandidate[] {
  if (candidates.length === 0) return [];

  // Sort by combined score (best first)
  const sorted = [...candidates].sort((a, b) =>
    ((b.combinedScore || b.baseFidelity) - (a.combinedScore || a.baseFidelity))
  );

  const selected: FusionSiteCandidate[] = [];

  for (const candidate of sorted) {
    // Check distance from all already-selected candidates
    const farEnough = selected.every(sel =>
      Math.abs(sel.position - candidate.position) >= minDistance
    );

    if (farEnough) {
      selected.push(candidate);
      if (selected.length >= numToSelect) break;
    }
  }

  // Sort by position for output
  return selected.sort((a, b) => a.position - b.position);
}

/**
 * Region coverage info
 */
interface RegionCoverage extends TargetPosition {
  candidatesFound: number;
  hasSufficientCandidates: boolean;
}

/**
 * Feasibility assessment result
 */
export interface FeasibilityAssessment {
  feasible: boolean;
  sequenceLength?: number;
  numFragments?: number;
  numJunctions?: number;
  totalCandidates?: number;
  regionCoverage?: RegionCoverage[];
  insufficientRegions?: number[];
  recommendation?: string;
  reason?: string;
}

/**
 * Feasibility options
 */
interface FeasibilityOptions {
  enzyme?: string;
  minFragmentSize?: number;
  minCandidatesPerRegion?: number;
}

/**
 * Quick scan to determine if a sequence is suitable for fusion site optimization
 *
 * @param sequence - DNA sequence
 * @param numFragments - Desired number of fragments
 * @param options - Options
 * @returns Feasibility assessment
 */
export function assessFeasibility(
  sequence: string,
  numFragments: number,
  options: FeasibilityOptions = {}
): FeasibilityAssessment {
  const {
    enzyme = 'BsaI',
    minFragmentSize = 200,
    minCandidatesPerRegion = 3,
  } = options;

  const numJunctions = numFragments - 1;
  const seqLength = sequence.length;

  // Basic checks
  if (seqLength < minFragmentSize * numFragments) {
    return {
      feasible: false,
      reason: `Sequence too short (${seqLength}bp) for ${numFragments} fragments of minimum ${minFragmentSize}bp`,
    };
  }

  // Quick scan
  const candidates = scanForFusionSites(sequence, {
    enzyme,
    includeEfficiency: true,
    maxCandidates: 2000,
  });

  // Generate target regions
  const targets = generateTargetPositions(seqLength, numFragments);

  // Check each region
  const regionCoverage: RegionCoverage[] = targets.map(target => {
    const inRegion = candidates.filter(c =>
      c.position >= target.searchRegion.start &&
      c.position <= target.searchRegion.end
    );
    return {
      ...target,
      candidatesFound: inRegion.length,
      hasSufficientCandidates: inRegion.length >= minCandidatesPerRegion,
    };
  });

  const insufficientRegions = regionCoverage.filter(r => !r.hasSufficientCandidates);

  return {
    feasible: insufficientRegions.length === 0,
    sequenceLength: seqLength,
    numFragments,
    numJunctions,
    totalCandidates: candidates.length,
    regionCoverage,
    insufficientRegions: insufficientRegions.map(r => r.index),
    recommendation: insufficientRegions.length === 0
      ? 'Sufficient candidates found - proceed with optimization'
      : `Regions ${insufficientRegions.map(r => r.index).join(', ')} have insufficient candidates`,
  };
}

export {
  isInWindow,
  isInForbiddenRegion,
  getOverhangFidelity,
};
