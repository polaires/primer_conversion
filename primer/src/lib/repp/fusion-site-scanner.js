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
 * Check if a position is within any of the allowed windows
 * @param {number} pos - Position to check
 * @param {Array} windows - Array of {start, end} objects
 * @returns {boolean} True if in a window
 */
function isInWindow(pos, windows) {
  if (!windows || windows.length === 0) return true;
  return windows.some(w => pos >= w.start && pos <= w.end);
}

/**
 * Check if a position is in a forbidden region
 * @param {number} pos - Position to check
 * @param {Array} regions - Array of {start, end} objects
 * @returns {boolean} True if in forbidden region
 */
function isInForbiddenRegion(pos, regions) {
  if (!regions || regions.length === 0) return false;
  return regions.some(r => pos >= r.start && pos < r.end);
}

/**
 * Get base overhang fidelity from experimental matrix
 * Uses centralized function from goldengate.js
 * @param {string} overhang - 4bp overhang sequence
 * @param {string} enzyme - Enzyme name
 * @returns {number} Base fidelity (0-1)
 */
function getOverhangFidelity(overhang, enzyme) {
  const result = getOverhangFidelityExperimental(overhang, enzyme);
  return result.fidelity;
}

/**
 * Scan sequence for all potential fusion sites
 *
 * This is the core scanning function that identifies candidate junction
 * positions in a DNA sequence.
 *
 * @param {string} sequence - DNA sequence to scan
 * @param {Object} options - Scanning options
 * @returns {Array} Array of candidate fusion sites
 */
export function scanForFusionSites(sequence, options = {}) {
  const {
    enzyme = 'BsaI',
    searchWindows = null,          // [{start, end}, ...] - allowed regions
    forbiddenRegions = [],         // [{start, end, reason}, ...] - excluded regions
    minDistanceFromEnds = 50,      // Minimum distance from sequence ends
    minDistanceBetween = 100,      // Minimum distance between candidates (for pre-filtering)
    includeEfficiency = true,      // Calculate efficiency metrics
    includeSiteCreation = false,   // Check for site creation (expensive)
    maxCandidates = 1000,          // Maximum candidates to return
  } = options;

  if (!sequence || sequence.length < minDistanceFromEnds * 2 + 4) {
    return [];
  }

  const candidates = [];
  const seq = sequence.toUpperCase();

  // Determine enzyme overhang length
  const enzData = GOLDEN_GATE_ENZYMES[enzyme];
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

    // ═══════════════════════════════════════════════════════════════════
    // HARD FILTERS - exclude completely
    // ═══════════════════════════════════════════════════════════════════

    // Exclude palindromic overhangs (can self-ligate)
    if (isPalindrome(overhang)) continue;

    // Exclude homopolymers (very poor fidelity)
    if (isHomopolymer(overhang)) continue;

    // Get base fidelity from matrix
    const baseFidelity = getOverhangFidelity(overhang, enzyme);
    if (baseFidelity === 0) continue; // No ligation data

    // ═══════════════════════════════════════════════════════════════════
    // Build candidate object
    // ═══════════════════════════════════════════════════════════════════

    const candidate = {
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
      candidate.efficiency = efficiency.efficiency;
      candidate.efficiencyWarnings = efficiency.warnings;
      candidate.isOptimalEfficiency = efficiency.isOptimal;
      candidate.isAcceptableEfficiency = efficiency.isAcceptable;
    }

    // Optional: Check for site creation (expensive)
    if (includeSiteCreation) {
      const siteCheck = checkSiteCreation(seq, pos, enzyme);
      candidate.siteCreationRisk = siteCheck.hasRisk;
      candidate.siteCreationSeverity = siteCheck.severity;
      if (siteCheck.hasRisk) {
        candidate.siteCreationWarnings = siteCheck.risks.map(r => r.message);
      }
    }

    candidates.push(candidate);

    // Limit total candidates
    if (candidates.length >= maxCandidates) break;
  }

  return candidates;
}

/**
 * Scan and rank fusion sites by quality
 *
 * @param {string} sequence - DNA sequence to scan
 * @param {Object} options - Scanning options
 * @returns {Object} Ranked candidates and statistics
 */
export function scanAndRankFusionSites(sequence, options = {}) {
  const {
    topN = 100,  // Number of top candidates to return
    sortBy = 'fidelity',  // 'fidelity', 'efficiency', 'combined'
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
  let sorted;
  switch (sortBy) {
    case 'efficiency':
      sorted = scored.sort((a, b) => (b.efficiency || 0) - (a.efficiency || 0));
      break;
    case 'combined':
      sorted = scored.sort((a, b) => b.combinedScore - a.combinedScore);
      break;
    case 'fidelity':
    default:
      sorted = scored.sort((a, b) => b.baseFidelity - a.baseFidelity);
  }

  // Statistics
  const stats = {
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
 * Find candidates in specific target regions
 *
 * @param {string} sequence - DNA sequence
 * @param {Array} targetRegions - Array of {start, end, weight} objects
 * @param {Object} options - Scanning options
 * @returns {Object} Candidates organized by region
 */
export function scanTargetRegions(sequence, targetRegions, options = {}) {
  const {
    candidatesPerRegion = 10,
    enzyme = 'BsaI',
    ...scanOptions
  } = options;

  const results = targetRegions.map((region, index) => {
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
      .sort((a, b) => b.combinedScore - a.combinedScore);

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
 * Generate evenly-spaced junction targets for a sequence
 *
 * @param {number} sequenceLength - Length of the sequence
 * @param {number} numFragments - Number of fragments desired
 * @param {Object} options - Options
 * @returns {Array} Target positions and search regions
 */
export function generateTargetPositions(sequenceLength, numFragments, options = {}) {
  const {
    minDistanceFromEnds = 50,
    searchRadius = 50,  // How far from ideal to search
  } = options;

  const numJunctions = numFragments - 1;
  if (numJunctions <= 0) {
    return [];
  }

  // Calculate ideal fragment size
  const usableLength = sequenceLength - 2 * minDistanceFromEnds;
  const idealFragmentSize = usableLength / numFragments;

  const targets = [];
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
 * @param {Array} candidates - Array of candidate objects
 * @param {number} minDistance - Minimum distance between selections
 * @param {number} numToSelect - Number of candidates to select
 * @returns {Array} Filtered candidates
 */
export function filterByDistance(candidates, minDistance, numToSelect) {
  if (candidates.length === 0) return [];

  // Sort by combined score (best first)
  const sorted = [...candidates].sort((a, b) =>
    (b.combinedScore || b.baseFidelity) - (a.combinedScore || a.baseFidelity)
  );

  const selected = [];

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
 * Quick scan to determine if a sequence is suitable for fusion site optimization
 *
 * @param {string} sequence - DNA sequence
 * @param {number} numFragments - Desired number of fragments
 * @param {Object} options - Options
 * @returns {Object} Feasibility assessment
 */
export function assessFeasibility(sequence, numFragments, options = {}) {
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
  const regionCoverage = targets.map(target => {
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
