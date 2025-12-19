/**
 * Alternative Generation Module for Unified Designer
 *
 * Implements state-of-the-art primer alternative generation with:
 * 1. Tiered scoring architecture (quick → medium → full) to balance speed and accuracy
 * 2. Unified piecewise scoring throughout pipeline (no ranking inversion)
 * 3. Geometric mean diversity selection (parameter-free)
 * 4. Hybrid selection: top by score + Pareto-diverse alternatives
 * 5. UI-friendly labels explaining each alternative's strengths
 *
 * Tiered Scoring Architecture:
 * ┌─────────────────────────────────────────────────────────────────┐
 * │ TIER 1: Quick Screening (~1ms for 1000 candidates)             │
 * │ - Hard rejection: Tm, GC, length outside acceptable bounds     │
 * │ - Soft penalties aligned with piecewise optimal zones          │
 * │ - Purpose: Eliminate obviously bad candidates fast             │
 * ├─────────────────────────────────────────────────────────────────┤
 * │ TIER 2: Medium Scoring (~20ms for 100 candidates)              │
 * │ - Individual primer piecewise scores (Tm, GC, length, etc.)    │
 * │ - Applied after diversity selection narrows the pool           │
 * │ - Purpose: Rank primers for pair evaluation                    │
 * ├─────────────────────────────────────────────────────────────────┤
 * │ TIER 3: Full Pair Scoring (~50ms for promising pairs)          │
 * │ - Heterodimer ΔG calculation (expensive O(n²))                 │
 * │ - Full composite score with all weights                        │
 * │ - Purpose: Final ranking of primer pairs                       │
 * └─────────────────────────────────────────────────────────────────┘
 */

import {
  selectDiverseCandidates,
  selectHybrid,
  selectParetoOptimal,
  identifyStrengths,
  generateLabel,
  generateExplanation,
} from './diversitySelection.js';
import {
  scoreTm,
  scoreGc,
  scoreTerminal3DG,
  scoreTmDiff,
  scoreHairpin,
  scoreHomodimer,
  scoreHeterodimer,
  scoreGcClamp,
  scoreLength,
  scoreHomopolymer,
  score3PrimeComposition,
  calculateCompositeScore,
  classifyQuality,
} from './scoring.js';
import { calculateHeterodimerDG } from './equilibrium.js';
import { reverseComplement } from './sequenceUtils.js';
import { gcCache, tmCache } from './tm.js';
import { dgCache } from './fold.js';
import { calculate3primeTerminalDG } from './tmQ5.js';

// Import constants from primers.js for consistency
import { LEN_MIN, LEN_MAX, LEN_MAX_EXTENDED } from './primers.js';

// Thresholds aligned with piecewise scoring (from scoring.js)
const THRESHOLDS = {
  tm: {
    optimalLow: 55,
    optimalHigh: 60,
    acceptableLow: 50,
    acceptableHigh: 65,
    hardMin: 45,
    hardMax: 72,
  },
  gc: {
    optimalLow: 0.40,
    optimalHigh: 0.60,
    acceptableLow: 0.30,
    acceptableHigh: 0.70,
    hardMin: 0.20,
    hardMax: 0.80,
  },
  length: {
    optimalLow: 18,
    optimalHigh: 24,
    acceptableLow: 15,
    acceptableHigh: 35,
    hardMin: 12,
    hardMax: 60, // Extended to support longer primers for Tm matching
  },
};

/**
 * TIER 1: Quick score function for initial screening
 *
 * Designed to be fast (no complex calculations) while maintaining
 * ranking alignment with piecewise scoring to avoid inversion.
 *
 * Key insight: The quick score doesn't need to be perfectly accurate -
 * it just needs to NOT exclude good candidates. False positives are
 * acceptable; false negatives are not.
 *
 * @param {Object} primer - Primer with tm, gc, length properties
 * @returns {number} Quick score 0-1
 */
function calculateQuickScore(primer) {
  const { tm, gc, length } = THRESHOLDS;

  // Hard rejections - must match piecewise "floor score" conditions
  if (primer.tm < tm.hardMin || primer.tm > tm.hardMax) return 0;
  if (primer.gc < gc.hardMin || primer.gc > gc.hardMax) return 0;
  if (primer.length < length.hardMin || primer.length > length.hardMax) return 0;

  let score = 1.0;

  // Tm scoring - aligned with piecewise zones
  if (primer.tm >= tm.optimalLow && primer.tm <= tm.optimalHigh) {
    // In optimal zone - full score
  } else if (primer.tm >= tm.acceptableLow && primer.tm <= tm.acceptableHigh) {
    // In acceptable zone - linear decay matching piecewise (1.0 → 0.7)
    const distFromOptimal = primer.tm < tm.optimalLow
      ? tm.optimalLow - primer.tm
      : primer.tm - tm.optimalHigh;
    const acceptableRange = primer.tm < tm.optimalLow
      ? tm.optimalLow - tm.acceptableLow
      : tm.acceptableHigh - tm.optimalHigh;
    score -= 0.3 * (distFromOptimal / acceptableRange);
  } else {
    // Below acceptable - steep penalty
    score -= 0.5;
  }

  // GC scoring - aligned with piecewise zones
  if (primer.gc >= gc.optimalLow && primer.gc <= gc.optimalHigh) {
    // In optimal zone - full score
  } else if (primer.gc >= gc.acceptableLow && primer.gc <= gc.acceptableHigh) {
    // In acceptable zone - linear decay
    const distFromOptimal = primer.gc < gc.optimalLow
      ? gc.optimalLow - primer.gc
      : primer.gc - gc.optimalHigh;
    const acceptableRange = primer.gc < gc.optimalLow
      ? gc.optimalLow - gc.acceptableLow
      : gc.acceptableHigh - gc.optimalHigh;
    score -= 0.2 * (distFromOptimal / acceptableRange);
  } else {
    // Outside acceptable - penalty
    score -= 0.4;
  }

  // Length scoring - aligned with piecewise zones
  if (primer.length >= length.optimalLow && primer.length <= length.optimalHigh) {
    // In optimal zone - full score
  } else if (primer.length >= length.acceptableLow && primer.length <= length.acceptableHigh) {
    // In acceptable zone - mild penalty
    score -= 0.1;
  } else {
    // Outside acceptable - penalty
    score -= 0.3;
  }

  // Ensure primers that pass all hard limits get at least a minimum positive score
  // This prevents cumulative soft penalties from excluding edge-case primers
  // that are still valid (e.g., very GC-rich or AT-rich regions)
  const MIN_PASSING_SCORE = 0.01;
  const clampedScore = Math.max(0, Math.min(1, score));

  // If all hard limits pass but score is 0, give it the minimum passing score
  return clampedScore === 0 ? MIN_PASSING_SCORE : clampedScore;
}

/**
 * TIER 2: Medium scoring for individual primers
 *
 * Uses actual piecewise scoring functions but only for individual primers,
 * not pair-level calculations. This provides accurate ranking while
 * deferring expensive heterodimer calculations.
 *
 * @param {Object} primer - Primer object with seq, tm, gc, dg
 * @returns {number} Medium score 0-1 (average of piecewise scores)
 */
function calculateMediumScore(primer) {
  const scores = calculatePrimerPiecewiseScores(primer);

  // Weighted average of individual primer scores
  // Weights reflect relative importance for individual primer quality
  const weights = {
    tm: 0.25,
    gc: 0.20,
    terminal3DG: 0.20,
    length: 0.10,
    gcClamp: 0.15,
    homopolymer: 0.10,
  };

  let totalScore = 0;
  let totalWeight = 0;

  totalScore += scores.tm * weights.tm;
  totalScore += scores.gc * weights.gc;
  totalScore += scores.terminal3DG * weights.terminal3DG;
  totalScore += scores.length * weights.length;
  totalScore += scores.gcClamp * weights.gcClamp;
  totalScore += scores.homopolymer * weights.homopolymer;

  totalWeight = Object.values(weights).reduce((a, b) => a + b, 0);

  return totalScore / totalWeight;
}

/**
 * Calculate piecewise scores for a single primer
 *
 * @param {Object} primer - Primer object with seq, tm, gc, dg
 * @returns {Object} Piecewise scores
 */
function calculatePrimerPiecewiseScores(primer) {
  const terminal3DGValue = calculate3primeTerminalDG(primer.seq).dG;
  return {
    tm: scoreTm(primer.tm),
    gc: scoreGc(primer.gc),
    terminal3DG: scoreTerminal3DG(terminal3DGValue),
    length: scoreLength(primer.length),
    gcClamp: scoreGcClamp(primer.seq),
    homopolymer: scoreHomopolymer(primer.seq),
    threePrimeComp: score3PrimeComposition(primer.seq, terminal3DGValue),
  };
}

/**
 * Calculate pair-level piecewise scores
 *
 * @param {Object} fwd - Forward primer
 * @param {Object} rev - Reverse primer
 * @param {number} annealingTemp - Annealing temperature for dimer calculations
 * @returns {Object} Pair scores and composite
 */
function calculatePairScores(fwd, rev, annealingTemp = 55) {
  const fwdScores = calculatePrimerPiecewiseScores(fwd);
  const revScores = calculatePrimerPiecewiseScores(rev);

  // Pair-level calculations
  const tmDiff = Math.abs(fwd.tm - rev.tm);
  const heterodimerDG = calculateHeterodimerDG(fwd.seq, rev.seq, annealingTemp);

  const pairScores = {
    tmDiff: scoreTmDiff(fwd.tm, rev.tm),
    heterodimer: scoreHeterodimer(heterodimerDG),
  };

  // Combine into scoring object for composite calculation
  const allScores = {
    tmFwd: fwdScores.tm,
    tmRev: revScores.tm,
    gcFwd: fwdScores.gc,
    gcRev: revScores.gc,
    terminal3DG: (fwdScores.terminal3DG + revScores.terminal3DG) / 2,
    lengthFwd: fwdScores.length,
    lengthRev: revScores.length,
    gcClampFwd: fwdScores.gcClamp,
    gcClampRev: revScores.gcClamp,
    homopolymerFwd: fwdScores.homopolymer,
    homopolymerRev: revScores.homopolymer,
    threePrimeCompFwd: fwdScores.threePrimeComp,
    threePrimeCompRev: revScores.threePrimeComp,
    tmDiff: pairScores.tmDiff,
    heterodimer: pairScores.heterodimer,
  };

  const composite = calculateCompositeScore(allScores);
  const quality = classifyQuality(composite.score);

  return {
    scores: allScores,
    compositeScore: composite.score,
    qualityTier: quality.tier,
    tmDiff,
    heterodimerDG,
    fwdScores,
    revScores,
  };
}

/**
 * Generate primer candidates from a sequence region
 *
 * @param {string} seq - Template sequence
 * @param {boolean} isForward - Direction
 * @param {Object} options - Generation options
 * @returns {Array} Array of primer candidates with position info
 */
function generateCandidates(seq, isForward, options = {}) {
  const {
    minLength = LEN_MIN,
    maxLength = LEN_MAX_EXTENDED, // Default to extended length for better Tm matching
    addLen = 0,
  } = options;

  // Pre-compute caches
  const gc = gcCache(seq);
  const tm = tmCache(seq);
  const dg = dgCache(seq);

  const candidates = [];

  // Generate candidates at each valid position
  const maxStart = addLen - 0 + 1; // Simplified from range calculation
  for (let start = 0; start < Math.min(maxStart, seq.length - minLength); start++) {
    for (let len = minLength; len <= maxLength; len++) {
      const end = start + len;
      if (end > seq.length) continue;

      const primerSeq = seq.slice(start, end);
      const primerTm = tm[start + addLen]?.[end - 1] || tm[start]?.[end - 1];
      const primerGc = gc[start]?.[end - 1];
      const primerDg = dg[start]?.[end - 1];

      if (primerTm === undefined || primerGc === undefined) continue;

      const candidate = {
        seq: primerSeq,
        length: len,
        tm: primerTm,
        gc: primerGc,
        dg: primerDg || 0,
        startPos: start,
        endPos: end,
        isForward,
      };

      // Calculate quick score for initial filtering
      candidate.quickScore = calculateQuickScore(candidate);

      if (candidate.quickScore > 0) {
        candidates.push(candidate);
      }
    }
  }

  return candidates;
}

/**
 * Generate alternatives with diversity selection
 *
 * Main entry point for the unified designer.
 * Uses a three-tier scoring architecture for optimal speed/accuracy balance.
 *
 * @param {string} template - Template sequence for the region to amplify
 * @param {Object} options - Generation options
 * @returns {Array} Array of primer pair alternatives with labels
 */
export function generateAlternativesUnified(template, options = {}) {
  const {
    numAlternatives = 8,
    maxCandidatesPerDirection = 50,
    topCandidatesForPairing = 20,  // Tier 2 reduces to this before pairing
    annealingTemp = 55,
    maxTmDiff = 8, // Allow up to tier 3
    earlyRejectTmDiff = 10, // Quick rejection threshold (looser than final)
  } = options;

  const seq = template.toUpperCase();

  if (seq.length < LEN_MAX * 2) {
    throw new Error(`Template too short for alternative generation: ${seq.length}bp`);
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // TIER 1: Quick Screening - Generate and filter candidates (~1ms)
  // ═══════════════════════════════════════════════════════════════════════════

  // Generate forward candidates from 5' end
  const fwdRegion = seq.slice(0, LEN_MAX + 10);
  const fwdCandidates = generateCandidates(fwdRegion, true, options);

  // Generate reverse candidates from 3' end (reverse complement for proper orientation)
  // This matches primers.js which reverse complements the 3' region before primer generation
  const revRegion = reverseComplement(seq.slice(-LEN_MAX - 10));
  const revCandidates = generateCandidates(revRegion, false, options);

  if (fwdCandidates.length === 0 || revCandidates.length === 0) {
    return [];
  }

  // Apply diversity-aware selection using quick scores
  // This prevents all candidates clustering in one "sweet spot"
  const diverseRanges = {
    maxPositionDiff: 20,
    maxLengthDiff: LEN_MAX - LEN_MIN,
    maxTmDiff: 15,
  };

  const diverseFwd = selectDiverseCandidates(
    fwdCandidates,
    maxCandidatesPerDirection,
    { scoreKey: 'quickScore', ranges: diverseRanges }
  );

  const diverseRev = selectDiverseCandidates(
    revCandidates,
    maxCandidatesPerDirection,
    { scoreKey: 'quickScore', ranges: diverseRanges }
  );

  // ═══════════════════════════════════════════════════════════════════════════
  // TIER 2: Medium Scoring - Apply piecewise scores to individuals (~20ms)
  // ═══════════════════════════════════════════════════════════════════════════

  // Calculate medium scores for diverse candidates
  const scoredFwd = diverseFwd.map(p => ({
    ...p,
    mediumScore: calculateMediumScore(p),
    piecewiseScores: calculatePrimerPiecewiseScores(p),
  }));

  const scoredRev = diverseRev.map(p => ({
    ...p,
    mediumScore: calculateMediumScore(p),
    piecewiseScores: calculatePrimerPiecewiseScores(p),
  }));

  // Re-select using medium scores for better accuracy
  const topFwd = selectDiverseCandidates(
    scoredFwd,
    topCandidatesForPairing,
    { scoreKey: 'mediumScore', ranges: diverseRanges }
  );

  const topRev = selectDiverseCandidates(
    scoredRev,
    topCandidatesForPairing,
    { scoreKey: 'mediumScore', ranges: diverseRanges }
  );

  // ═══════════════════════════════════════════════════════════════════════════
  // TIER 3: Full Pair Scoring - Heterodimer and composite (~50ms)
  // ═══════════════════════════════════════════════════════════════════════════

  const pairs = [];

  for (const fwd of topFwd) {
    for (const rev of topRev) {
      // Stage 3a: Quick pair rejection (Tm diff check - very fast)
      const quickTmDiff = Math.abs(fwd.tm - rev.tm);
      if (quickTmDiff > earlyRejectTmDiff) continue;

      // Stage 3b: Full scoring including heterodimer (expensive)
      const pairScoring = calculatePairScores(fwd, rev, annealingTemp);

      // Final Tm diff check with accurate value
      if (pairScoring.tmDiff > maxTmDiff) continue;

      // Early rejection for severe heterodimer (ΔG < -12 kcal/mol)
      // At this level, heterodimer formation will likely cause PCR failure
      // Score would be ~0.05 (95% penalty), so reject early to save computation
      if (pairScoring.heterodimerDG < -12) continue;

      // Calculate amplicon length
      const ampliconLength = seq.length - fwd.startPos - (revRegion.length - rev.endPos);

      pairs.push({
        forward: {
          sequence: fwd.seq,
          length: fwd.length,
          tm: Math.round(fwd.tm * 10) / 10,
          gc: fwd.gc,
          dg: fwd.dg,
          gcPercent: `${(fwd.gc * 100).toFixed(1)}%`,
          hasGCClamp: /[GC]$/.test(fwd.seq),
          startPos: fwd.startPos,
          mediumScore: Math.round(fwd.mediumScore * 100) / 100,
        },
        reverse: {
          sequence: rev.seq,
          length: rev.length,
          tm: Math.round(rev.tm * 10) / 10,
          gc: rev.gc,
          dg: rev.dg,
          gcPercent: `${(rev.gc * 100).toFixed(1)}%`,
          hasGCClamp: /[GC]$/.test(rev.seq),
          startPos: rev.startPos,
          mediumScore: Math.round(rev.mediumScore * 100) / 100,
        },
        compositeScore: pairScoring.compositeScore,
        qualityTier: pairScoring.qualityTier,
        tmDiff: Math.round(pairScoring.tmDiff * 10) / 10,
        heterodimerDG: Math.round(pairScoring.heterodimerDG * 10) / 10,
        ampliconLength,
        scores: pairScoring.scores,
        // For diversity selection
        score: pairScoring.compositeScore,
        startPos: fwd.startPos,
        length: fwd.length + rev.length,
        tm: (fwd.tm + rev.tm) / 2,
      });
    }
  }

  if (pairs.length === 0) {
    return [];
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // FINAL: Hybrid Selection + Labeling
  // ═══════════════════════════════════════════════════════════════════════════

  const paretoObjectives = [
    p => p.compositeScore,           // Higher score better
    p => -p.tmDiff,                  // Lower Tm diff better
    p => p.heterodimerDG,            // Higher (less negative) better
    p => -p.ampliconLength,          // Shorter amplicon as option
  ];

  const diversePairRanges = {
    maxPositionDiff: seq.length,
    maxLengthDiff: 20,
    maxTmDiff: 15,
  };

  const selectedPairs = selectHybrid(pairs, numAlternatives, {
    scoreKey: 'compositeScore',
    numByScore: Math.ceil(numAlternatives / 2),
    paretoObjectives,
    ranges: diversePairRanges,
  });

  // Add labels and explanations
  const reference = selectedPairs[0];

  return selectedPairs.map((pair, index) => {
    const strengths = identifyStrengths(pair, pairs);
    const label = generateLabel(strengths, pair.qualityTier);
    const explanation = index === 0
      ? 'Highest composite score across all metrics'
      : generateExplanation(pair, strengths, reference);

    return {
      forward: pair.forward,
      reverse: pair.reverse,
      compositeScore: pair.compositeScore,
      qualityTier: pair.qualityTier,
      tmDiff: pair.tmDiff,
      heterodimerDG: pair.heterodimerDG,
      ampliconLength: pair.ampliconLength,
      // UI labeling
      label,
      explanation,
      strengths,
      selectionReason: pair.selectionReason,
      // Keep scores for debugging
      scores: pair.scores,
    };
  });
}

/**
 * Tier assignment for backward compatibility
 *
 * @param {Object} pair - Primer pair
 * @returns {string} Tier label
 */
export function assignTier(pair) {
  const tmDiff = pair.tmDiff || Math.abs(pair.forward.tm - pair.reverse.tm);
  const worstDg = Math.min(pair.forward.dg || 0, pair.reverse.dg || 0);

  if (tmDiff <= 2 && worstDg > -3) return 'excellent';
  if (tmDiff <= 5 && worstDg > -5) return 'good';
  if (tmDiff <= 8) return 'acceptable';
  return 'poor';
}

// ═══════════════════════════════════════════════════════════════════════════════
// ENHANCED ALTERNATIVE GENERATION WITH SLIDING WINDOW + AMPLICON SCORING
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Generate alternatives with sliding window and amplicon-level scoring
 *
 * This enhanced version provides:
 * 1. Position-flexible primer generation (not just 5'/3' ends)
 * 2. Amplicon-level scoring (secondary structure, repeats, GC distribution)
 * 3. Application-specific presets (qPCR, colony PCR, etc.)
 * 4. Target amplicon length optimization
 *
 * @param {string} template - Full template sequence
 * @param {Object} options - Generation options
 * @returns {Array} Array of primer pair alternatives with comprehensive scoring
 */
export function generateAlternativesEnhanced(template, options = {}) {
  const {
    // Basic options
    numAlternatives = 10,
    annealingTemp = 55,
    maxTmDiff = 8,
    earlyRejectTmDiff = 10,

    // Sliding window options
    useSlidingWindow = true,
    fwdSearchStart = 0,
    fwdSearchEnd = null,
    revSearchStart = null,
    revSearchEnd = null,
    windowStep = 5,

    // Amplicon options
    minAmpliconLength = 100,
    maxAmpliconLength = 2000,
    targetAmpliconLength = null,
    application = 'pcr',  // 'pcr', 'qpcr', 'sequencing'

    // Scoring options
    includeAmpliconAnalysis = true,
    ampliconScoreWeight = 0.15,  // Weight of amplicon score in final ranking

    // Candidate limits
    maxCandidatesPerDirection = 80,
    topCandidatesForPairing = 25,
  } = options;

  const seq = template.toUpperCase();
  const seqLen = seq.length;

  if (seqLen < minAmpliconLength + LEN_MIN * 2) {
    throw new Error(`Template too short: ${seqLen}bp (need at least ${minAmpliconLength + LEN_MIN * 2}bp)`);
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // STAGE 1: Generate candidates using sliding window or fixed-end
  // ═══════════════════════════════════════════════════════════════════════════

  let fwdCandidates, revCandidates;

  if (useSlidingWindow) {
    const slidingResult = generateSlidingWindowCandidates(seq, {
      fwdStart: fwdSearchStart,
      fwdEnd: fwdSearchEnd,
      revStart: revSearchStart,
      revEnd: revSearchEnd,
      minAmpliconLength,
      maxAmpliconLength,
      targetAmpliconLength,
      windowStep,
      maxCandidatesPerRegion: maxCandidatesPerDirection,
    });
    fwdCandidates = slidingResult.forward;
    revCandidates = slidingResult.reverse;
  } else {
    // Fall back to fixed-end generation
    const fwdRegion = seq.slice(0, LEN_MAX + 10);
    fwdCandidates = generateCandidates(fwdRegion, true, options);

    const revRegion = reverseComplement(seq.slice(-LEN_MAX - 10));
    revCandidates = generateCandidates(revRegion, false, options);
  }

  if (fwdCandidates.length === 0 || revCandidates.length === 0) {
    return [];
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // STAGE 2: Medium scoring for individuals
  // ═══════════════════════════════════════════════════════════════════════════

  const scoredFwd = fwdCandidates.map(p => ({
    ...p,
    mediumScore: calculateMediumScore(p),
    piecewiseScores: calculatePrimerPiecewiseScores(p),
  }));

  const scoredRev = revCandidates.map(p => ({
    ...p,
    mediumScore: calculateMediumScore(p),
    piecewiseScores: calculatePrimerPiecewiseScores(p),
  }));

  // Diversity-aware selection
  const diverseRanges = {
    maxPositionDiff: seqLen / 4,
    maxLengthDiff: LEN_MAX - LEN_MIN,
    maxTmDiff: 15,
  };

  const topFwd = selectDiverseCandidates(
    scoredFwd,
    topCandidatesForPairing,
    { scoreKey: 'mediumScore', ranges: diverseRanges }
  );

  const topRev = selectDiverseCandidates(
    scoredRev,
    topCandidatesForPairing,
    { scoreKey: 'mediumScore', ranges: diverseRanges }
  );

  // ═══════════════════════════════════════════════════════════════════════════
  // STAGE 3: Full pair scoring with amplicon analysis
  // ═══════════════════════════════════════════════════════════════════════════

  const pairs = [];

  for (const fwd of topFwd) {
    for (const rev of topRev) {
      // Quick rejection
      const quickTmDiff = Math.abs(fwd.tm - rev.tm);
      if (quickTmDiff > earlyRejectTmDiff) continue;

      // Calculate amplicon boundaries
      const fwdEnd = fwd.templatePos !== undefined ? fwd.templatePos + fwd.length : fwd.endPos;
      const revStart = rev.templatePos !== undefined ? rev.templatePos : seqLen - LEN_MAX - 10 + rev.startPos;

      // Calculate amplicon length
      let ampliconLength;
      if (useSlidingWindow && fwd.templatePos !== undefined && rev.templateEndPos !== undefined) {
        ampliconLength = rev.templateEndPos - fwd.templatePos;
      } else {
        ampliconLength = seq.length - fwd.startPos - (LEN_MAX + 10 - rev.endPos);
      }

      // Amplicon length constraints
      if (ampliconLength < minAmpliconLength || ampliconLength > maxAmpliconLength) continue;

      // Full pair scoring
      const pairScoring = calculatePairScores(fwd, rev, annealingTemp);

      if (pairScoring.tmDiff > maxTmDiff) continue;
      if (pairScoring.heterodimerDG < -12) continue;

      // Amplicon-level analysis
      let ampliconAnalysis = null;
      let adjustedScore = pairScoring.compositeScore;

      if (includeAmpliconAnalysis) {
        const ampliconStart = fwd.templatePos ?? fwd.startPos;
        const ampliconEnd = rev.templateEndPos ?? (seqLen - (LEN_MAX + 10 - rev.endPos));
        const ampliconSeq = seq.slice(ampliconStart, ampliconEnd);

        if (ampliconSeq.length > 0) {
          ampliconAnalysis = scoreAmplicon(ampliconSeq, {
            targetLength: targetAmpliconLength,
            application,
            optimalLengthRange: application === 'qpcr' ? [70, 150] : [200, 800],
          });

          // Blend amplicon score into composite
          adjustedScore = pairScoring.compositeScore * (1 - ampliconScoreWeight) +
                          (ampliconAnalysis.score * 100) * ampliconScoreWeight;
        }
      }

      pairs.push({
        forward: {
          sequence: fwd.seq,
          length: fwd.length,
          tm: Math.round(fwd.tm * 10) / 10,
          gc: fwd.gc,
          dg: fwd.dg,
          gcPercent: `${(fwd.gc * 100).toFixed(1)}%`,
          hasGCClamp: /[GC]$/.test(fwd.seq),
          startPos: fwd.startPos,
          templatePos: fwd.templatePos,
          mediumScore: Math.round(fwd.mediumScore * 100) / 100,
        },
        reverse: {
          sequence: rev.seq,
          length: rev.length,
          tm: Math.round(rev.tm * 10) / 10,
          gc: rev.gc,
          dg: rev.dg,
          gcPercent: `${(rev.gc * 100).toFixed(1)}%`,
          hasGCClamp: /[GC]$/.test(rev.seq),
          startPos: rev.startPos,
          templatePos: rev.templatePos,
          mediumScore: Math.round(rev.mediumScore * 100) / 100,
        },
        compositeScore: Math.round(adjustedScore * 10) / 10,
        primerScore: pairScoring.compositeScore,
        ampliconScore: ampliconAnalysis ? Math.round(ampliconAnalysis.score * 100) : null,
        qualityTier: pairScoring.qualityTier,
        tmDiff: Math.round(pairScoring.tmDiff * 10) / 10,
        heterodimerDG: Math.round(pairScoring.heterodimerDG * 10) / 10,
        ampliconLength,
        ampliconAnalysis,
        scores: pairScoring.scores,
        // For diversity selection
        score: adjustedScore,
        startPos: fwd.templatePos ?? fwd.startPos,
        length: fwd.length + rev.length,
        tm: (fwd.tm + rev.tm) / 2,
      });
    }
  }

  if (pairs.length === 0) {
    return [];
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // STAGE 4: Hybrid selection with amplicon considerations
  // ═══════════════════════════════════════════════════════════════════════════

  const paretoObjectives = [
    p => p.compositeScore,           // Higher score better
    p => -p.tmDiff,                  // Lower Tm diff better
    p => p.heterodimerDG,            // Higher (less negative) better
    p => p.ampliconScore ?? 80,      // Higher amplicon score better
  ];

  // Add amplicon length objective based on target
  if (targetAmpliconLength !== null) {
    paretoObjectives.push(p => -Math.abs(p.ampliconLength - targetAmpliconLength));
  } else {
    paretoObjectives.push(p => -p.ampliconLength);  // Shorter as option
  }

  const diversePairRanges = {
    maxPositionDiff: seqLen,
    maxLengthDiff: 20,
    maxTmDiff: 15,
  };

  const selectedPairs = selectHybrid(pairs, numAlternatives, {
    scoreKey: 'compositeScore',
    numByScore: Math.ceil(numAlternatives / 2),
    paretoObjectives,
    ranges: diversePairRanges,
  });

  // Add labels and explanations
  const reference = selectedPairs[0];

  return selectedPairs.map((pair, index) => {
    const strengths = identifyStrengthsEnhanced(pair, pairs);
    const label = generateLabelEnhanced(strengths, pair);
    const explanation = index === 0
      ? 'Best overall score combining primer quality and amplicon characteristics'
      : generateExplanationEnhanced(pair, strengths, reference);

    return {
      forward: pair.forward,
      reverse: pair.reverse,
      compositeScore: pair.compositeScore,
      primerScore: pair.primerScore,
      ampliconScore: pair.ampliconScore,
      qualityTier: pair.qualityTier,
      tmDiff: pair.tmDiff,
      heterodimerDG: pair.heterodimerDG,
      ampliconLength: pair.ampliconLength,
      ampliconAnalysis: pair.ampliconAnalysis,
      // UI labeling
      label,
      explanation,
      strengths,
      selectionReason: pair.selectionReason,
      scores: pair.scores,
    };
  });
}

/**
 * Identify strengths including amplicon characteristics
 */
function identifyStrengthsEnhanced(candidate, allCandidates) {
  const strengths = identifyStrengths(candidate, allCandidates);

  // Add amplicon-specific strengths
  if (candidate.ampliconAnalysis && allCandidates.length > 1) {
    const ampliconScores = allCandidates
      .filter(c => c.ampliconScore !== null)
      .map(c => c.ampliconScore);

    if (ampliconScores.length > 0) {
      const maxAmpliconScore = Math.max(...ampliconScores);
      if (candidate.ampliconScore >= maxAmpliconScore - 1 && !strengths.includes('bestOverall')) {
        strengths.push('cleanAmplicon');
      }
    }

    // Check for optimal amplicon length match
    if (candidate.ampliconAnalysis.lengthScore >= 0.95) {
      strengths.push('optimalSize');
    }
  }

  return strengths;
}

/**
 * SVG icon paths for enhanced labels
 */
const ENHANCED_LABEL_ICONS = {
  bestOverall: 'M12 2l3.09 6.26L22 9.27l-5 4.87 1.18 6.88L12 17.77l-6.18 3.25L7 14.14 2 9.27l6.91-1.01L12 2z', // star
  bestTmMatch: 'M12 2a10 10 0 100 20 10 10 0 000-20zm0 18a8 8 0 110-16 8 8 0 010 16zm0-14a6 6 0 100 12 6 6 0 000-12zm0 10a4 4 0 110-8 4 4 0 010 8zm0-6a2 2 0 100 4 2 2 0 000-4z', // target
  safestDimer: 'M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10zm-1.5-5.5l-3-3 1.5-1.5 1.5 1.5 4-4 1.5 1.5-5.5 5.5z', // shield-check
  shortestAmplicon: 'M4 14h4v4H6v-2H4v-2zm0-4h2V8h2V6H4v4zm12 6h-2v2h4v-4h-2v2zm-2-6V8h2v2h2V6h-4v4z', // minimize
  cleanAmplicon: 'M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41L9 16.17z', // check
  optimalSize: 'M3 5v14a2 2 0 002 2h14a2 2 0 002-2V5a2 2 0 00-2-2H5a2 2 0 00-2 2zm16 0v2h-2V5h2zm-4 0v4h-2V5h2zm-4 0v2h-2V5h2zm-4 0v4H5V5h2z', // ruler
  compact: 'M4 14h4v4H6v-2H4v-2zm0-4h2V8h2V6H4v4zm12 6h-2v2h4v-4h-2v2zm-2-6V8h2v2h2V6h-4v4z', // minimize
  specific: 'M21 11V3h-8l3.29 3.29-10 10L3 13v8h8l-3.29-3.29 10-10L21 11z', // maximize
};

/**
 * Generate label including amplicon considerations
 * Returns structured object with SVG icon for professional rendering
 */
function generateLabelEnhanced(strengths, pair) {
  // Standard labels with SVG icons
  if (strengths.includes('bestOverall')) {
    return { key: 'bestOverall', text: 'Best', svgPath: ENHANCED_LABEL_ICONS.bestOverall };
  }
  if (strengths.includes('bestTmMatch')) {
    return { key: 'bestTmMatch', text: 'Tm', svgPath: ENHANCED_LABEL_ICONS.bestTmMatch };
  }
  if (strengths.includes('safestDimer')) {
    return { key: 'safestDimer', text: 'Safe', svgPath: ENHANCED_LABEL_ICONS.safestDimer };
  }
  if (strengths.includes('shortestAmplicon')) {
    return { key: 'shortestAmplicon', text: 'Short', svgPath: ENHANCED_LABEL_ICONS.shortestAmplicon };
  }

  // Amplicon-aware labels
  if (strengths.includes('cleanAmplicon')) {
    return { key: 'cleanAmplicon', text: 'Clean', svgPath: ENHANCED_LABEL_ICONS.cleanAmplicon };
  }
  if (strengths.includes('optimalSize')) {
    return { key: 'optimalSize', text: 'Optimal', svgPath: ENHANCED_LABEL_ICONS.optimalSize };
  }
  if (strengths.includes('compact')) {
    return { key: 'compact', text: 'Compact', svgPath: ENHANCED_LABEL_ICONS.compact };
  }
  if (strengths.includes('specific')) {
    return { key: 'specific', text: 'Specific', svgPath: ENHANCED_LABEL_ICONS.specific };
  }

  return null;
}

/**
 * Generate explanation including amplicon characteristics
 */
function generateExplanationEnhanced(candidate, strengths, reference) {
  if (strengths.includes('cleanAmplicon')) {
    return 'Amplicon with minimal secondary structure risk';
  }
  if (strengths.includes('optimalSize')) {
    return 'Amplicon length in optimal range for application';
  }

  // Amplicon-specific comparison
  if (candidate.ampliconAnalysis && reference.ampliconAnalysis) {
    if (candidate.ampliconScore > reference.ampliconScore + 5) {
      return `Better amplicon quality (${candidate.ampliconScore} vs ${reference.ampliconScore})`;
    }
  }

  // Fall back to standard explanation
  return generateExplanation(candidate, strengths, reference);
}

// Export scoring functions for testing
export {
  calculateQuickScore,
  calculateMediumScore,
  calculatePrimerPiecewiseScores,
  calculatePairScores,
  THRESHOLDS,
};
