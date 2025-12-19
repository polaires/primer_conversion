/**
 * Diversity Selection Module
 *
 * Implements geometric mean-based diversity selection for primer candidates.
 * This approach balances quality scores with positional/structural diversity
 * without requiring a tunable weight parameter.
 *
 * Key insight: Geometric mean √(quality × diversity) naturally penalizes
 * extreme imbalance between the two objectives.
 */

/**
 * Calculate normalized distance between two primer candidates
 * Uses position, length, and Tm as diversity dimensions
 *
 * @param {Object} a - First candidate
 * @param {Object} b - Second candidate
 * @param {Object} ranges - Normalization ranges
 * @returns {number} Normalized distance (0-1)
 */
export function calculateDistance(a, b, ranges) {
  const {
    maxPositionDiff = 100,  // Typical template region size
    maxLengthDiff = 10,     // LEN_MAX - LEN_MIN roughly
    maxTmDiff = 15,         // Typical Tm range
  } = ranges;

  const posDiff = Math.abs((a.startPos || 0) - (b.startPos || 0)) / maxPositionDiff;
  const lenDiff = Math.abs(a.length - b.length) / maxLengthDiff;
  const tmDiff = Math.abs(a.tm - b.tm) / maxTmDiff;

  // Euclidean distance in normalized space, capped at 1.0
  const distance = Math.sqrt(posDiff * posDiff + lenDiff * lenDiff + tmDiff * tmDiff);
  return Math.min(1.0, distance / Math.sqrt(3)); // Normalize by max possible distance
}

/**
 * Select diverse candidates using geometric mean of quality and diversity
 *
 * Algorithm:
 * 1. First pick: best by score (no diversity yet)
 * 2. Subsequent picks: maximize √(quality × minDistanceToSelected)
 *
 * This naturally balances quality and diversity without a weight parameter.
 *
 * @param {Array} candidates - Array of candidates with score and position info
 * @param {number} numToSelect - Number of candidates to select
 * @param {Object} options - Selection options
 * @returns {Array} Selected diverse candidates
 */
export function selectDiverseCandidates(candidates, numToSelect, options = {}) {
  const {
    scoreKey = 'score',           // Property name for score
    minScore = 0.0,               // Minimum acceptable score
    ranges = {},                  // Distance normalization ranges
  } = options;

  if (candidates.length === 0) return [];

  // Sort by score descending
  const sorted = [...candidates].sort((a, b) => b[scoreKey] - a[scoreKey]);

  // Filter by minimum score FIRST (before any early returns)
  const viable = sorted.filter(c => c[scoreKey] >= minScore);
  if (viable.length === 0) return sorted.slice(0, numToSelect);
  if (viable.length <= numToSelect) return viable;

  const selected = [];

  // First pick: best by score (no diversity consideration yet)
  selected.push(viable[0]);
  const remaining = viable.slice(1);

  // Subsequent picks using geometric mean
  while (selected.length < numToSelect && remaining.length > 0) {
    let bestIdx = 0;
    let bestValue = -Infinity;

    for (let i = 0; i < remaining.length; i++) {
      const candidate = remaining[i];

      // Quality component: normalized to 0.1-1.0 range
      // Use 0.1 floor to ensure lowest-scoring viable candidates still have non-zero quality
      const maxScore = viable[0][scoreKey];
      const minViableScore = viable[viable.length - 1][scoreKey];
      const scoreRange = maxScore - minViableScore || 1;
      const rawQuality = (candidate[scoreKey] - minViableScore) / scoreRange;
      // Map 0-1 to 0.1-1.0 so even lowest score has meaningful quality
      const quality = 0.1 + rawQuality * 0.9;

      // Diversity component: minimum distance to any already-selected candidate
      const minDistance = Math.min(
        ...selected.map(s => calculateDistance(candidate, s, ranges))
      );

      // Add small epsilon to avoid zero distance issues
      const adjustedDistance = minDistance + 0.01;

      // Geometric mean: √(quality × diversity)
      // This naturally balances the two without a weight parameter
      const value = Math.sqrt(quality * adjustedDistance);

      if (value > bestValue) {
        bestValue = value;
        bestIdx = i;
      }
    }

    selected.push(remaining.splice(bestIdx, 1)[0]);
  }

  return selected;
}

/**
 * Select Pareto-optimal candidates across multiple objectives
 *
 * A candidate is Pareto-optimal if no other candidate is better in ALL objectives.
 *
 * @param {Array} candidates - Array of candidates
 * @param {Array} objectives - Array of objective functions (higher = better)
 * @returns {Array} Pareto-optimal candidates
 */
export function selectParetoOptimal(candidates, objectives) {
  if (candidates.length === 0) return [];

  const dominated = new Set();

  for (let i = 0; i < candidates.length; i++) {
    if (dominated.has(i)) continue;

    for (let j = 0; j < candidates.length; j++) {
      if (i === j || dominated.has(j)) continue;

      // Check if candidate i dominates candidate j
      let dominatesAll = true;
      let betterInSome = false;

      for (const obj of objectives) {
        const valI = obj(candidates[i]);
        const valJ = obj(candidates[j]);

        if (valI < valJ) dominatesAll = false;
        if (valI > valJ) betterInSome = true;
      }

      if (dominatesAll && betterInSome) {
        dominated.add(j);
      }
    }
  }

  return candidates.filter((_, idx) => !dominated.has(idx));
}

/**
 * Hybrid selection: top N by score + Pareto-diverse alternatives
 *
 * Guarantees user expectations (top scores) while providing genuinely different options.
 *
 * @param {Array} candidates - Array of candidates
 * @param {number} numToSelect - Total number to select
 * @param {Object} options - Selection options
 * @returns {Array} Selected candidates with labels
 */
export function selectHybrid(candidates, numToSelect, options = {}) {
  const {
    scoreKey = 'score',
    numByScore = Math.ceil(numToSelect / 2),
    paretoObjectives = null,
    ranges = {},
  } = options;

  if (candidates.length === 0) return [];
  if (candidates.length <= numToSelect) {
    return candidates.map((c, i) => ({
      ...c,
      selectionReason: i === 0 ? 'bestOverall' : 'included',
    }));
  }

  // Sort by score
  const sorted = [...candidates].sort((a, b) => b[scoreKey] - a[scoreKey]);

  const selected = [];

  // Add top N by score
  const topByScore = sorted.slice(0, numByScore);
  topByScore.forEach((c, i) => {
    selected.push({
      ...c,
      selectionReason: i === 0 ? 'bestOverall' : 'topByScore',
    });
  });

  // Remaining slots filled with diverse Pareto-optimal candidates
  const numRemaining = numToSelect - selected.length;
  if (numRemaining > 0) {
    const remaining = sorted.filter(c =>
      !selected.some(s => s === c || (s.id && s.id === c.id))
    );

    if (paretoObjectives) {
      // Use Pareto selection
      const paretoSet = selectParetoOptimal(remaining, paretoObjectives);

      // Select most diverse from Pareto front
      const diverse = selectDiverseCandidates(paretoSet, numRemaining, {
        scoreKey,
        ranges,
      });

      diverse.forEach(c => {
        selected.push({
          ...c,
          selectionReason: 'paretoOptimal',
        });
      });
    } else {
      // Fallback: just use diversity selection
      const diverse = selectDiverseCandidates(remaining, numRemaining, {
        scoreKey,
        ranges,
      });

      diverse.forEach(c => {
        selected.push({
          ...c,
          selectionReason: 'diverse',
        });
      });
    }
  }

  return selected;
}

/**
 * Identify the unique strength of a candidate relative to others
 * Only assigns a strength if this candidate is THE best (not tied)
 *
 * @param {Object} candidate - The candidate to analyze
 * @param {Array} allCandidates - All candidates for comparison
 * @param {number} index - Index of this candidate in the list
 * @returns {Array} Array of strength identifiers (usually 0-1 items)
 */
export function identifyStrengths(candidate, allCandidates, index = -1) {
  const strengths = [];

  if (allCandidates.length === 0) return strengths;

  const candidateScore = candidate.score || candidate.compositeScore || 0;
  const candidateTmDiff = Math.abs(candidate.tmDiff || (candidate.forward?.tm - candidate.reverse?.tm) || 0);
  const candidateDG = candidate.heterodimerDG || candidate.heterodimer?.dg || -10;
  const candidateLen = candidate.ampliconLength || 0;
  const candidatePrimerLen = (candidate.forward?.length || 0) + (candidate.reverse?.length || 0);

  // Calculate stats for all candidates
  const scores = allCandidates.map(c => c.score || c.compositeScore || 0);
  const tmDiffs = allCandidates.map(c => Math.abs(c.tmDiff || (c.forward?.tm - c.reverse?.tm) || 0));
  const dimerDGs = allCandidates.map(c => c.heterodimerDG || c.heterodimer?.dg || -10);
  const ampliconLens = allCandidates.map(c => c.ampliconLength || 0).filter(l => l > 0);
  const primerLens = allCandidates.map(c => (c.forward?.length || 0) + (c.reverse?.length || 0));

  // Best overall score - only if uniquely best (not tied) and gap >= 1
  const maxScore = Math.max(...scores);
  const countAtMax = scores.filter(s => s >= maxScore).length;
  const secondBestScore = Math.max(...scores.filter(s => s < maxScore), 0);
  if (candidateScore >= maxScore && countAtMax === 1 && maxScore - secondBestScore >= 1) {
    strengths.push('bestOverall');
    return strengths; // bestOverall is exclusive
  }

  // Best Tm matching - only if uniquely best and significantly better than others
  const minTmDiff = Math.min(...tmDiffs);
  const countAtMin = tmDiffs.filter(t => t <= minTmDiff + 0.05).length;
  const secondMinTmDiff = Math.min(...tmDiffs.filter(t => t > minTmDiff + 0.05), 99);
  if (candidateTmDiff <= minTmDiff + 0.05 && countAtMin === 1 && secondMinTmDiff - minTmDiff >= 0.3) {
    strengths.push('bestTmMatch');
    return strengths;
  }

  // Safest dimer - only if uniquely safer
  const safestDG = Math.max(...dimerDGs);
  const countAtSafest = dimerDGs.filter(d => d >= safestDG - 0.1).length;
  const secondSafestDG = Math.max(...dimerDGs.filter(d => d < safestDG - 0.1), -20);
  if (candidateDG >= safestDG - 0.1 && countAtSafest === 1 && safestDG - secondSafestDG >= 0.5) {
    strengths.push('safestDimer');
    return strengths;
  }

  // Shortest amplicon - only if uniquely shortest
  if (ampliconLens.length > 0 && candidateLen > 0) {
    const minLen = Math.min(...ampliconLens);
    const countAtMin = ampliconLens.filter(l => l <= minLen).length;
    const secondMinLen = Math.min(...ampliconLens.filter(l => l > minLen), 9999);
    if (candidateLen <= minLen && countAtMin === 1 && secondMinLen - minLen >= 5) {
      strengths.push('shortestAmplicon');
      return strengths;
    }
  }

  // Shortest primers - compact design
  const minPrimerLen = Math.min(...primerLens.filter(l => l > 0));
  if (candidatePrimerLen > 0 && candidatePrimerLen <= minPrimerLen && primerLens.filter(l => l === minPrimerLen).length === 1) {
    strengths.push('compact');
    return strengths;
  }

  // Longest primers - maximum specificity
  const maxPrimerLen = Math.max(...primerLens);
  if (candidatePrimerLen >= maxPrimerLen && primerLens.filter(l => l === maxPrimerLen).length === 1) {
    strengths.push('specific');
    return strengths;
  }

  return strengths;
}

/**
 * Identify quality badges for a primer pair
 * These are non-exclusive features that indicate good properties
 * Uses SVG path data for professional icons
 *
 * @param {Object} pair - Primer pair with forward/reverse primers
 * @returns {Array} Array of badge objects with id, svgPath, label, and tooltip
 */
export function identifyBadges(pair) {
  const badges = [];

  if (!pair) return badges;

  const fwd = pair.forward || {};
  const rev = pair.reverse || {};
  const heterodimerDG = pair.heterodimerDG ?? pair.heterodimer?.dg ?? -10;
  const tmDiff = Math.abs(pair.tmDiff ?? (fwd.tm - rev.tm) ?? 5);

  // GC Clamp: Both primers end with G or C
  const fwdGC = fwd.hasGCClamp ?? /[GC]$/i.test(fwd.sequence || '');
  const revGC = rev.hasGCClamp ?? /[GC]$/i.test(rev.sequence || '');
  if (fwdGC && revGC) {
    badges.push({
      id: 'gcClamp',
      // Lock icon path
      svgPath: 'M12 17a2 2 0 100-4 2 2 0 000 4zm6-7V8a6 6 0 10-12 0v2a2 2 0 00-2 2v6a2 2 0 002 2h12a2 2 0 002-2v-6a2 2 0 00-2-2zm-8-2a4 4 0 118 0v2H8V8z',
      label: 'GC',
      tooltip: 'Both primers have GC clamp - stable 3\' binding',
    });
  }

  // Low Tm difference: < 1°C
  if (tmDiff <= 1.0) {
    badges.push({
      id: 'lowTmDiff',
      // Target/crosshair icon path
      svgPath: 'M12 2a10 10 0 100 20 10 10 0 000-20zm0 18a8 8 0 110-16 8 8 0 010 16zm0-14a6 6 0 100 12 6 6 0 000-12zm0 10a4 4 0 110-8 4 4 0 010 8zm0-6a2 2 0 100 4 2 2 0 000-4z',
      label: 'Tm',
      tooltip: `ΔTm ${tmDiff.toFixed(1)}°C - excellent annealing consistency`,
    });
  }

  // Safe heterodimer: ΔG > -6 (threshold for stable dimer)
  if (heterodimerDG > -6) {
    badges.push({
      id: 'safeDimer',
      // Shield check icon path
      svgPath: 'M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10zm-1.5-5.5l-3-3 1.5-1.5 1.5 1.5 4-4 1.5 1.5-5.5 5.5z',
      label: 'Safe',
      tooltip: `Heterodimer ΔG ${heterodimerDG.toFixed(1)} - no stable primer-dimer`,
    });
  }

  // Optimal primer lengths: both 18-24bp
  const fwdLen = fwd.length || fwd.sequence?.length || 0;
  const revLen = rev.length || rev.sequence?.length || 0;
  if (fwdLen >= 18 && fwdLen <= 24 && revLen >= 18 && revLen <= 24) {
    badges.push({
      id: 'optimalLength',
      // Ruler icon path
      svgPath: 'M3 5v14a2 2 0 002 2h14a2 2 0 002-2V5a2 2 0 00-2-2H5a2 2 0 00-2 2zm16 0v2h-2V5h2zm-4 0v4h-2V5h2zm-4 0v2h-2V5h2zm-4 0v4H5V5h2zm-2 6h2v2H5v-2zm4 0h2v2H9v-2zm4 0h2v2h-2v-2zm4 0h2v2h-2v-2z',
      label: 'Len',
      tooltip: `Optimal lengths (${fwdLen}+${revLen}bp) - balanced specificity/cost`,
    });
  }

  // Good 3' terminal stability: dg between -11 and -6 (stable but not too sticky)
  const fwdDg = fwd.dg ?? fwd.terminal3DG ?? -8;
  const revDg = rev.dg ?? rev.terminal3DG ?? -8;
  if (fwdDg >= -11 && fwdDg <= -6 && revDg >= -11 && revDg <= -6) {
    badges.push({
      id: 'stable3Prime',
      // Bolt/lightning icon path
      svgPath: 'M13 2L3 14h9l-1 8 10-12h-9l1-8z',
      label: "3'",
      tooltip: `Good 3' stability (${fwdDg.toFixed(1)}/${revDg.toFixed(1)}) - efficient priming`,
    });
  }

  return badges;
}

/**
 * Label SVG icon paths for professional rendering
 */
export const LABEL_SVG_PATHS = {
  bestOverall: 'M12 2l3.09 6.26L22 9.27l-5 4.87 1.18 6.88L12 17.77l-6.18 3.25L7 14.14 2 9.27l6.91-1.01L12 2z', // star
  bestTmMatch: 'M12 2a10 10 0 100 20 10 10 0 000-20zm0 18a8 8 0 110-16 8 8 0 010 16zm0-14a6 6 0 100 12 6 6 0 000-12zm0 10a4 4 0 110-8 4 4 0 010 8zm0-6a2 2 0 100 4 2 2 0 000-4z', // target
  safestDimer: 'M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10zm-1.5-5.5l-3-3 1.5-1.5 1.5 1.5 4-4 1.5 1.5-5.5 5.5z', // shield-check
  shortestAmplicon: 'M4 14h4v4H6v-2H4v-2zm0-4h2V8h2V6H4v4zm12 6h-2v2h4v-4h-2v2zm-2-6V8h2v2h2V6h-4v4z', // minimize
  compact: 'M4 14h4v4H6v-2H4v-2zm0-4h2V8h2V6H4v4zm12 6h-2v2h4v-4h-2v2zm-2-6V8h2v2h2V6h-4v4z', // minimize
  specific: 'M21 11V3h-8l3.29 3.29-10 10L3 13v8h8l-3.29-3.29 10-10L21 11z', // maximize
};

/**
 * Generate a label object with SVG icon and text based on strengths
 *
 * @param {Array} strengths - Array of strength identifiers
 * @param {string} tier - Quality tier (excellent, good, etc.) - not shown in label
 * @returns {Object|null} Label object with {key, text, svgPath}, or null if no unique strength
 */
export function generateLabel(strengths, tier = '') {
  // Return structured label with icon path for professional rendering
  if (strengths.includes('bestOverall')) {
    return { key: 'bestOverall', text: 'Best', svgPath: LABEL_SVG_PATHS.bestOverall };
  }
  if (strengths.includes('bestTmMatch')) {
    return { key: 'bestTmMatch', text: 'Tm', svgPath: LABEL_SVG_PATHS.bestTmMatch };
  }
  if (strengths.includes('safestDimer')) {
    return { key: 'safestDimer', text: 'Safe', svgPath: LABEL_SVG_PATHS.safestDimer };
  }
  if (strengths.includes('shortestAmplicon')) {
    return { key: 'shortestAmplicon', text: 'Short', svgPath: LABEL_SVG_PATHS.shortestAmplicon };
  }
  if (strengths.includes('compact')) {
    return { key: 'compact', text: 'Compact', svgPath: LABEL_SVG_PATHS.compact };
  }
  if (strengths.includes('specific')) {
    return { key: 'specific', text: 'Specific', svgPath: LABEL_SVG_PATHS.specific };
  }

  return null; // No label if no unique strength
}

/**
 * Generate an explanation for why this alternative was selected
 * Only returns a value if it adds insight beyond what's visible in the card
 *
 * @param {Object} candidate - The candidate
 * @param {Array} strengths - Its identified strengths
 * @param {Object} reference - Reference candidate for comparison
 * @returns {string|null} Explanation text or null if nothing useful to add
 */
export function generateExplanation(candidate, strengths, reference = null) {
  // Only provide explanations for labeled alternatives - explain the practical benefit
  if (strengths.includes('bestOverall')) {
    return 'Top composite score across all metrics';
  }

  if (strengths.includes('bestTmMatch')) {
    return 'Most consistent annealing temperature';
  }

  if (strengths.includes('safestDimer')) {
    return 'Lowest risk of primer-primer binding';
  }

  if (strengths.includes('shortestAmplicon')) {
    return 'Faster PCR cycles, easier gel analysis';
  }

  if (strengths.includes('compact')) {
    return 'Shortest primers - lower synthesis cost';
  }

  if (strengths.includes('specific')) {
    return 'Longest primers - maximum specificity';
  }

  // No explanation for unlabeled alternatives - visible metrics tell the story
  return null;
}
