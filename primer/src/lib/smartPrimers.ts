/**
 * Smart Primer Design Module
 *
 * Implements iterative recalibration for optimal primer design.
 * Key features:
 * - GC clamp optimization (extend length to achieve G/C at 3' end)
 * - 3' end composition balancing (target optimal ΔG range)
 * - Multi-variant evaluation (generate and compare alternatives)
 *
 * Based on research showing that 3' end stability is critical for primer function:
 * - Priming requires 3' terminal ΔG < -11 kcal/mol for stable duplex
 * - GC clamp (G/C in last 2 bases) improves extension efficiency
 * - AT-rich 3' ends have higher mispriming risk
 */

import { tm } from './tm.js';
import { calculate3primeTerminalDG } from './tmQ5.js';
import {
  scoreGcClamp,
  scoreTerminal3DG,
  scoreLength,
  scoreGc,
  scoreTm,
  scoreHomopolymer,
  calculateCompositeScore,
  classifyQuality,
} from './scoring.js';
import {
  calculateHairpinDG,
  calculateHomodimerDG,
} from './equilibrium.js';

// =============================================================================
// Types and Interfaces
// =============================================================================

export interface SmartDesignConstraints {
  minLength: number;
  maxLength: number;
  optimalLengthRange: [number, number];
  gcClampTarget: number;
  gcClampMinimum: number;
  terminal3DG: {
    optimal: [number, number];
    acceptable: [number, number];
    minimum: number;
    maximum: number;
  };
  minAcceptableScore: number;
  targetScore: number;
  maxIterations: number;
  maxVariantsPerPosition: number;
}

export interface Analysis3Prime {
  sequence: string;
  last5: string;
  last3: string;
  last2: string;
  lastBase: string;
  gcCounts: {
    last2: number;
    last3: number;
    last5: number;
  };
  terminalDG: number;
  hasGcClamp: boolean;
  endsWithGC: boolean;
  patterns: {
    hasPolyA: boolean;
    hasPolyT: boolean;
    hasPolyAT: boolean;
  };
  quality: string;
  issues: string[];
  canImprove: boolean;
}

export interface StructureSeverity {
  level: string;
  label: string;
  shortLabel: string;
  message: string | null;
  tooltip: string;
  color: string;
  icon: string;
  shouldWarn: boolean;
  actionRequired?: boolean;
}

export interface StructureBadge {
  text: string;
  className: string;
  bgColor: string;
  textColor: string;
  borderColor: string;
  severity: string;
  tooltip: string;
  dG: number;
}

export interface LengthVariant {
  seq: string;
  length: number;
  lengthDelta: number;
  startPos: number;
  analysis: Analysis3Prime;
  priority: number;
}

export interface OptimizedPrimer {
  seq: string;
  length: number;
  originalSeq?: string;
  originalLength?: number;
  optimized: boolean;
  lengthDelta?: number;
  reason: string;
  analysis: Analysis3Prime;
  originalAnalysis?: Analysis3Prime;
  variants: LengthVariant[];
}

export interface PrimerScore {
  seq: string;
  length: number;
  tm: number;
  gc: number;
  terminalDG: number;
  hairpinDG: number;
  homodimerDG: number;
  scores: Record<string, number>;
  compositeScore: number;
  qualityTier: string;
  analysis3Prime: Analysis3Prime;
}

export interface DesignSuggestions {
  currentScore: number;
  currentTier: string;
  canImprove: boolean;
  suggestions: Array<{
    priority: string;
    issue: string;
    action: string;
    expectedImprovement: string;
  }>;
  summary: string;
}

export interface OptimizedPair {
  forward: OptimizedPrimer & {
    scoring: PrimerScore;
    suggestions: DesignSuggestions;
  };
  reverse: OptimizedPrimer & {
    scoring: PrimerScore;
    suggestions: DesignSuggestions;
  };
  pairScore: number;
  pairQuality: { tier: string; color: string; label: string; description: string };
  tmDiff: number;
  optimizationSummary: {
    fwdOptimized: boolean;
    revOptimized: boolean;
    fwdReason: string;
    revReason: string;
    improvements: string[];
  };
  meetsTargetScore: boolean;
}

export interface QuickAssessment {
  seq: string;
  length: number;
  tm: number;
  gc: number;
  gcClamp: number;
  terminal3DG: number;
  quality: string;
  issues: string[];
  needsOptimization: boolean;
  canImproveWith: string | null;
}

// =============================================================================
// Constants
// =============================================================================

/**
 * Design constraints for smart optimization
 */
export const SMART_DESIGN_CONSTRAINTS: SmartDesignConstraints = {
  // Length bounds
  minLength: 15,
  maxLength: 32,
  optimalLengthRange: [18, 24],

  // GC clamp requirements
  gcClampTarget: 1,        // Ideal: 1 G/C in last 2 bases
  gcClampMinimum: 1,       // Minimum acceptable

  // 3' terminal ΔG requirements (kcal/mol)
  terminal3DG: {
    optimal: [-11, -6],    // Optimal range
    acceptable: [-14, -4], // Acceptable range
    minimum: -15,          // Below this is too tight
    maximum: -3,           // Above this is too loose
  },

  // Score thresholds
  minAcceptableScore: 60,
  targetScore: 75,

  // Iteration limits
  maxIterations: 5,
  maxVariantsPerPosition: 7,
};

// =============================================================================
// Analysis Functions
// =============================================================================

/**
 * Analyze 3' end composition quality
 */
export function analyze3PrimeEnd(seq: string): Analysis3Prime {
  seq = seq.toUpperCase();
  const last5 = seq.slice(-5);
  const last3 = seq.slice(-3);
  const last2 = seq.slice(-2);
  const lastBase = seq.slice(-1);

  // Count G/C in different windows
  const gcLast2 = (last2.match(/[GC]/g) || []).length;
  const gcLast3 = (last3.match(/[GC]/g) || []).length;
  const gcLast5 = (last5.match(/[GC]/g) || []).length;

  // Check for problematic patterns
  const hasPolyA = /AAA/.test(last5);
  const hasPolyT = /TTT/.test(last5);
  const hasPolyAT = /[AT]{4,}/.test(last5);
  const endsWithAT = /[AT]$/.test(seq);
  const endsWithGC = /[GC]$/.test(seq);

  // Calculate 3' terminal ΔG (actual binding strength of last 5 bases)
  const terminalDG = calculate3primeTerminalDG(seq).dG;

  // Determine quality tier
  let quality = 'excellent';
  const issues: string[] = [];

  if (gcLast2 === 0) {
    quality = 'poor';
    issues.push('No GC clamp (0 G/C in last 2 bases)');
  } else if (gcLast2 === 2 && !hasPolyAT) {
    quality = 'good';  // Strong clamp but slight mispriming risk
  }

  if (hasPolyAT) {
    quality = quality === 'excellent' ? 'acceptable' : 'poor';
    issues.push('Poly-A/T run in 3\' end');
  }

  if (terminalDG > -6) {
    quality = 'poor';
    issues.push(`3' terminal ΔG too weak: ${terminalDG.toFixed(1)} kcal/mol`);
  } else if (terminalDG < -14) {
    quality = quality === 'excellent' ? 'acceptable' : quality;
    issues.push(`3' terminal ΔG very strong: ${terminalDG.toFixed(1)} kcal/mol (mispriming risk)`);
  }

  if (gcLast5 <= 1 && quality !== 'poor') {
    quality = 'acceptable';
    issues.push('AT-rich 3\' end (≤1 G/C in last 5 bases)');
  }

  return {
    sequence: seq,
    last5,
    last3,
    last2,
    lastBase,
    gcCounts: {
      last2: gcLast2,
      last3: gcLast3,
      last5: gcLast5,
    },
    terminalDG: Math.round(terminalDG * 100) / 100,
    hasGcClamp: gcLast2 >= 1,
    endsWithGC,
    patterns: {
      hasPolyA,
      hasPolyT,
      hasPolyAT,
    },
    quality,
    issues,
    canImprove: quality !== 'excellent' && issues.length > 0,
  };
}

/**
 * Classify 3' structure severity for user guidance
 *
 * This function provides tiered severity levels that align with actual PCR impact,
 * helping users understand when they MUST act vs when it's informational.
 *
 * Severity Levels:
 * - critical: Last 3 bases trapped + ΔG ≤ -4 → Extension will fail
 * - warning: Last 5 bases involved + ΔG ≤ -3 → May reduce efficiency
 * - info: Last 10 bases involved + ΔG > -2 → Usually not problematic
 * - none: No 3' involvement → Good to go
 */
export function classify3PrimeStructureSeverity({
  energy,
  basePairs = [],
  seqLength,
  structure = null,
}: {
  energy: number;
  basePairs?: Array<[number, number]>;
  seqLength: number;
  structure?: { has3PrimeInStem?: boolean; has3PrimeInLoop?: boolean } | null;
}): StructureSeverity {
  // Find which 3' positions are involved in base pairs
  const threePrimePositions = {
    last3: new Set([seqLength - 1, seqLength - 2, seqLength - 3]),
    last5: new Set([seqLength - 1, seqLength - 2, seqLength - 3, seqLength - 4, seqLength - 5]),
    last10: new Set(Array.from({ length: 10 }, (_, i) => seqLength - 1 - i)),
  };

  let involvesLast3 = false;
  let involvesLast5 = false;
  let involvesLast10 = false;

  for (const [i, j] of basePairs) {
    if (threePrimePositions.last3.has(i) || threePrimePositions.last3.has(j)) {
      involvesLast3 = true;
    }
    if (threePrimePositions.last5.has(i) || threePrimePositions.last5.has(j)) {
      involvesLast5 = true;
    }
    if (threePrimePositions.last10.has(i) || threePrimePositions.last10.has(j)) {
      involvesLast10 = true;
    }
  }

  // Also check structure info if provided
  if (structure) {
    if (structure.has3PrimeInStem || structure.has3PrimeInLoop) {
      involvesLast10 = true;
    }
  }

  // No 3' involvement at all
  if (!involvesLast10) {
    return {
      level: 'none',
      label: 'No 3\' Structure',
      shortLabel: '✓ Free',
      message: null,
      tooltip: 'The 3\' end is free and available for template binding. This is ideal for primer function.',
      color: 'success',
      icon: '✓',
      shouldWarn: false,
    };
  }

  // CRITICAL: Last 3 bases trapped + stable structure
  if (involvesLast3 && energy <= -4) {
    return {
      level: 'critical',
      label: '3\' End Blocked',
      shortLabel: '✗ Blocked',
      message: '3\' end is trapped in a stable structure - extension will fail',
      tooltip: `Critical Issue: The last 3 bases of the primer are base-paired in a stable hairpin (ΔG = ${energy.toFixed(1)} kcal/mol).

Why this matters:
• DNA polymerase REQUIRES a free 3' OH to begin extension
• This structure is stable enough to persist at annealing temperature
• PCR will likely fail or have very low efficiency

Action Required: Redesign the primer by adjusting length or position to free the 3' end.`,
      color: 'danger',
      icon: '🔴',
      shouldWarn: true,
      actionRequired: true,
    };
  }

  // WARNING: Last 5 bases involved + moderate structure
  if (involvesLast5 && energy <= -3) {
    return {
      level: 'warning',
      label: '3\' Structure Risk',
      shortLabel: '⚠ Risk',
      message: '3\' end partially involved in structure - may reduce efficiency',
      tooltip: `Warning: Bases near the 3' end are involved in secondary structure (ΔG = ${energy.toFixed(1)} kcal/mol).

Why this matters:
• The structure may compete with template binding
• PCR may work but with reduced efficiency
• Higher primer concentrations or DMSO may help

Recommendation: Consider optimizing if you experience low yields. Try:
• Adjust primer length by 1-2 bases
• Use touchdown PCR protocol
• Add 2-5% DMSO to reduce secondary structure`,
      color: 'warning',
      icon: '🟠',
      shouldWarn: true,
      actionRequired: false,
    };
  }

  // INFO: Structure present but weak or not involving critical bases
  if (involvesLast10 && energy > -2) {
    return {
      level: 'info',
      label: 'Weak 3\' Structure',
      shortLabel: 'ℹ Weak',
      message: 'Weak structure near 3\' end - typically not problematic',
      tooltip: `Information: A weak secondary structure is detected near the 3' end (ΔG = ${energy.toFixed(1)} kcal/mol).

Why this is usually OK:
• The structure is too weak to persist at annealing temperature
• It will "melt" when the primer binds to the template
• Normal PCR conditions should work fine

No action needed unless you experience unexpected issues.`,
      color: 'info',
      icon: '🟡',
      shouldWarn: false,
      actionRequired: false,
    };
  }

  // MODERATE: Structure in last 10 but not last 5, or moderate energy
  if (involvesLast10) {
    const isModerateSeverity = energy <= -3;
    return {
      level: isModerateSeverity ? 'moderate' : 'low',
      label: isModerateSeverity ? 'Moderate 3\' Structure' : 'Minor 3\' Structure',
      shortLabel: isModerateSeverity ? '△ Moderate' : 'ℹ Minor',
      message: isModerateSeverity
        ? 'Structure near 3\' region - monitor PCR efficiency'
        : 'Minor structure detected - usually fine',
      tooltip: `${isModerateSeverity ? 'Moderate' : 'Minor'} structure detected near the 3' end (ΔG = ${energy.toFixed(1)} kcal/mol).

The structure involves bases 5-10 from the 3' end but not the critical terminal bases.

${isModerateSeverity
  ? 'Monitor PCR efficiency. If yields are low, consider optimization.'
  : 'This typically does not affect PCR performance. Proceed with normal protocols.'}`,
      color: isModerateSeverity ? 'warning' : 'info',
      icon: isModerateSeverity ? '🟡' : 'ℹ️',
      shouldWarn: isModerateSeverity,
      actionRequired: false,
    };
  }

  // Default: shouldn't reach here, but be safe
  return {
    level: 'none',
    label: 'No Significant Structure',
    shortLabel: '✓ OK',
    message: null,
    tooltip: 'No significant secondary structure detected at the 3\' end.',
    color: 'success',
    icon: '✓',
    shouldWarn: false,
  };
}

/**
 * Get 3' structure badge for alternative designs
 *
 * Returns a compact badge object for displaying in alternative design cards.
 */
export function get3PrimeStructureBadge({
  energy,
  basePairs = [],
  seqLength,
  structure = null,
}: {
  energy: number;
  basePairs?: Array<[number, number]>;
  seqLength: number;
  structure?: { has3PrimeInStem?: boolean; has3PrimeInLoop?: boolean } | null;
}): StructureBadge | null {
  const severity = classify3PrimeStructureSeverity({ energy, basePairs, seqLength, structure });

  if (severity.level === 'none') {
    return null; // No badge for good primers
  }

  const badgeConfig: Record<string, {
    text: string;
    className: string;
    bgColor: string;
    textColor: string;
    borderColor: string;
  }> = {
    critical: {
      text: "3' Blocked",
      className: 'structure-badge-critical',
      bgColor: '#fee2e2',
      textColor: '#dc2626',
      borderColor: '#fecaca',
    },
    warning: {
      text: "3' Risk",
      className: 'structure-badge-warning',
      bgColor: '#fef3c7',
      textColor: '#d97706',
      borderColor: '#fde68a',
    },
    moderate: {
      text: "3' Moderate",
      className: 'structure-badge-moderate',
      bgColor: '#fef9c3',
      textColor: '#ca8a04',
      borderColor: '#fef08a',
    },
    low: {
      text: "3' Minor",
      className: 'structure-badge-low',
      bgColor: '#e0f2fe',
      textColor: '#0284c7',
      borderColor: '#bae6fd',
    },
    info: {
      text: "3' Weak",
      className: 'structure-badge-info',
      bgColor: '#f0f9ff',
      textColor: '#0369a1',
      borderColor: '#e0f2fe',
    },
  };

  const config = badgeConfig[severity.level] || badgeConfig.info;

  return {
    ...config,
    severity: severity.level,
    tooltip: severity.tooltip,
    dG: energy,
  };
}

/**
 * Generate length variants for a primer at the same binding position
 *
 * Extends or contracts the primer to achieve better 3' end composition.
 */
export function generateLengthVariants(
  fullSeq: string,
  startPos: number,
  currentLength: number,
  isFwd: boolean,
  options: {
    minLength?: number;
    maxLength?: number;
    extendMax?: number;
    contractMax?: number;
  } = {}
): LengthVariant[] {
  const {
    minLength = SMART_DESIGN_CONSTRAINTS.minLength,
    maxLength = SMART_DESIGN_CONSTRAINTS.maxLength,
    extendMax = 3,
    contractMax = 2,
  } = options;

  const variants: LengthVariant[] = [];

  // Calculate the range of lengths to try
  const minLen = Math.max(minLength, currentLength - contractMax);
  const maxLen = Math.min(maxLength, currentLength + extendMax);

  for (let len = minLen; len <= maxLen; len++) {
    let seq: string;

    if (isFwd) {
      // Forward primer: extend/contract from 3' end (into sequence)
      const endPos = startPos + len;
      if (endPos > fullSeq.length) continue;
      seq = fullSeq.slice(startPos, endPos);
    } else {
      // Reverse primer: extend/contract from 3' end (away from sequence start)
      // The sequence is already reverse complemented in the caller
      const endPos = startPos + len;
      if (endPos > fullSeq.length) continue;
      seq = fullSeq.slice(startPos, endPos);
    }

    if (seq.length !== len) continue;

    const analysis = analyze3PrimeEnd(seq);
    const lengthDelta = len - currentLength;

    variants.push({
      seq,
      length: len,
      lengthDelta,
      startPos,
      analysis,
      priority: calculateVariantPriority(analysis, lengthDelta),
    });
  }

  // Sort by priority (higher is better)
  variants.sort((a, b) => b.priority - a.priority);

  return variants;
}

/**
 * Calculate priority score for a variant
 */
function calculateVariantPriority(analysis: Analysis3Prime, lengthDelta: number): number {
  let priority = 0;

  // Strong preference for GC clamp
  if (analysis.gcCounts.last2 === 1) {
    priority += 100;  // Ideal
  } else if (analysis.gcCounts.last2 === 2) {
    priority += 80;   // Strong but slight mispriming risk
  } else {
    priority -= 50;   // No clamp
  }

  // Prefer ending with G/C
  if (analysis.endsWithGC) {
    priority += 30;
  }

  // Good terminal ΔG
  const { terminal3DG } = SMART_DESIGN_CONSTRAINTS;
  if (analysis.terminalDG >= terminal3DG.optimal[0] &&
      analysis.terminalDG <= terminal3DG.optimal[1]) {
    priority += 50;
  } else if (analysis.terminalDG >= terminal3DG.acceptable[0] &&
             analysis.terminalDG <= terminal3DG.acceptable[1]) {
    priority += 20;
  } else {
    priority -= 30;
  }

  // Penalize poly-A/T patterns
  if (analysis.patterns.hasPolyAT) {
    priority -= 40;
  }

  // Slight preference for minimal length change
  priority -= Math.abs(lengthDelta) * 5;

  // Prefer optimal length range
  if (analysis.sequence.length >= 18 && analysis.sequence.length <= 24) {
    priority += 15;
  }

  return priority;
}

/**
 * Optimize a single primer for GC clamp and 3' end composition
 */
export function optimizePrimer(
  fullSeq: string,
  startPos: number,
  currentLength: number,
  isFwd: boolean,
  options: {
    targetGcClamp?: boolean;
    target3PrimeDG?: boolean;
    maxLengthChange?: number;
  } = {}
): OptimizedPrimer {
  const {
    targetGcClamp = true,
    target3PrimeDG = true,
    maxLengthChange = 3,
  } = options;

  // Get current primer
  const currentSeq = fullSeq.slice(startPos, startPos + currentLength);
  const currentAnalysis = analyze3PrimeEnd(currentSeq);

  // If already optimal, return as-is
  if (currentAnalysis.quality === 'excellent') {
    return {
      seq: currentSeq,
      length: currentLength,
      optimized: false,
      reason: 'Already optimal',
      analysis: currentAnalysis,
      variants: [],
    };
  }

  // Generate variants
  const variants = generateLengthVariants(fullSeq, startPos, currentLength, isFwd, {
    extendMax: maxLengthChange,
    contractMax: Math.min(2, maxLengthChange),
  });

  // Filter variants that meet minimum requirements
  const viableVariants = variants.filter(v => {
    if (targetGcClamp && v.analysis.gcCounts.last2 === 0) return false;
    if (target3PrimeDG) {
      const { terminal3DG } = SMART_DESIGN_CONSTRAINTS;
      if (v.analysis.terminalDG > terminal3DG.maximum ||
          v.analysis.terminalDG < terminal3DG.minimum) return false;
    }
    return true;
  });

  // If no viable variants, return best of all
  const candidates = viableVariants.length > 0 ? viableVariants : variants;

  if (candidates.length === 0) {
    return {
      seq: currentSeq,
      length: currentLength,
      optimized: false,
      reason: 'No better variants found',
      analysis: currentAnalysis,
      variants: [],
    };
  }

  const best = candidates[0];

  // Check if best is actually better than current
  const currentPriority = calculateVariantPriority(currentAnalysis, 0);
  if (best.priority <= currentPriority) {
    return {
      seq: currentSeq,
      length: currentLength,
      optimized: false,
      reason: 'Current primer is already best option',
      analysis: currentAnalysis,
      variants: variants.slice(0, 5),
    };
  }

  // Return optimized result
  const improvements: string[] = [];
  if (best.analysis.gcCounts.last2 > currentAnalysis.gcCounts.last2) {
    improvements.push(`GC clamp improved (${currentAnalysis.gcCounts.last2} → ${best.analysis.gcCounts.last2})`);
  }
  if (best.analysis.quality !== currentAnalysis.quality) {
    improvements.push(`Quality improved (${currentAnalysis.quality} → ${best.analysis.quality})`);
  }
  if (best.lengthDelta !== 0) {
    const direction = best.lengthDelta > 0 ? 'extended' : 'contracted';
    improvements.push(`Length ${direction} by ${Math.abs(best.lengthDelta)}bp`);
  }

  return {
    seq: best.seq,
    length: best.length,
    originalSeq: currentSeq,
    originalLength: currentLength,
    optimized: true,
    lengthDelta: best.lengthDelta,
    reason: improvements.join('; ') || 'Improved 3\' end composition',
    analysis: best.analysis,
    originalAnalysis: currentAnalysis,
    variants: variants.slice(0, 5),
  };
}

/**
 * Score a primer variant comprehensively
 */
export function scorePrimerVariant(
  seq: string,
  template: string = '',
  isFwd: boolean = true,
  temperature: number = 55
): PrimerScore {
  const primerTm = tm(seq, '', true);
  const gcContent = (seq.match(/[GC]/gi) || []).length / seq.length;
  const terminalDG = calculate3primeTerminalDG(seq).dG;
  const hairpinDG = calculateHairpinDG(seq, temperature);
  const homodimerDG = calculateHomodimerDG(seq, temperature);

  const analysis = analyze3PrimeEnd(seq);

  // Individual scores
  const scores: Record<string, number> = {
    tm: scoreTm(primerTm),
    gc: scoreGc(gcContent),
    length: scoreLength(seq.length),
    gcClamp: scoreGcClamp(seq),
    homopolymer: scoreHomopolymer(seq),
    terminal3DG: scoreTerminal3DG(terminalDG),
    hairpin: Math.exp(-Math.max(0, -3 - hairpinDG) * 0.8),  // Simple hairpin score
    homodimer: Math.exp(-Math.max(0, -6 - homodimerDG) * 0.5),
    // 3' end composition score (new)
    threePrimeComposition: score3PrimeComposition(analysis),
  };

  // Calculate composite
  const composite = calculateCompositeScore({
    [`tm${isFwd ? 'Fwd' : 'Rev'}`]: scores.tm,
    [`gc${isFwd ? 'Fwd' : 'Rev'}`]: scores.gc,
    [`length${isFwd ? 'Fwd' : 'Rev'}`]: scores.length,
    [`gcClamp${isFwd ? 'Fwd' : 'Rev'}`]: scores.gcClamp,
    [`homopolymer${isFwd ? 'Fwd' : 'Rev'}`]: scores.homopolymer,
    [`hairpin${isFwd ? 'Fwd' : 'Rev'}`]: scores.hairpin,
    [`selfDimer${isFwd ? 'Fwd' : 'Rev'}`]: scores.homodimer,
    terminal3DG: scores.terminal3DG,
  });

  const quality = classifyQuality(composite.score);

  return {
    seq,
    length: seq.length,
    tm: Math.round(primerTm * 10) / 10,
    gc: Math.round(gcContent * 100) / 100,
    terminalDG: Math.round(terminalDG * 100) / 100,
    hairpinDG: Math.round(hairpinDG * 100) / 100,
    homodimerDG: Math.round(homodimerDG * 100) / 100,
    scores,
    compositeScore: composite.score,
    qualityTier: quality.tier,
    analysis3Prime: analysis,
  };
}

/**
 * Score 3' end composition quality
 *
 * Combines GC clamp, terminal ΔG, and pattern analysis into a single score.
 */
export function score3PrimeComposition(analysis: Analysis3Prime): number {
  let score = 1.0;

  // GC clamp component (40% weight)
  if (analysis.gcCounts.last2 === 1) {
    // Ideal
  } else if (analysis.gcCounts.last2 === 2) {
    score -= 0.06;  // Slight penalty for strong clamp
  } else if (analysis.gcCounts.last2 === 0) {
    score -= 0.20;  // Significant penalty for no clamp
  }

  // Terminal ΔG component (35% weight)
  const { terminal3DG } = SMART_DESIGN_CONSTRAINTS;
  if (analysis.terminalDG >= terminal3DG.optimal[0] &&
      analysis.terminalDG <= terminal3DG.optimal[1]) {
    // Optimal - no penalty
  } else if (analysis.terminalDG > terminal3DG.optimal[1]) {
    // Too weak
    const excess = analysis.terminalDG - terminal3DG.optimal[1];
    score -= Math.min(0.35, excess * 0.07);
  } else if (analysis.terminalDG < terminal3DG.optimal[0]) {
    // Too tight
    const excess = terminal3DG.optimal[0] - analysis.terminalDG;
    score -= Math.min(0.15, excess * 0.03);  // Less severe penalty
  }

  // Pattern component (25% weight)
  if (analysis.patterns.hasPolyAT) {
    score -= 0.15;
  }
  if (!analysis.endsWithGC) {
    score -= 0.05;
  }
  if (analysis.gcCounts.last5 <= 1) {
    score -= 0.05;  // AT-rich 3' end
  }

  return Math.max(0, Math.min(1, score));
}

/**
 * Generate design suggestions for improving a primer
 */
export function generateDesignSuggestions(primerScore: PrimerScore): DesignSuggestions {
  const suggestions: Array<{
    priority: string;
    issue: string;
    action: string;
    expectedImprovement: string;
  }> = [];
  const { scores, analysis3Prime, compositeScore, qualityTier } = primerScore;

  // Check each component and suggest improvements

  // GC Clamp
  if (scores.gcClamp < 0.8) {
    if (analysis3Prime.gcCounts.last2 === 0) {
      suggestions.push({
        priority: 'high',
        issue: 'No GC clamp',
        action: 'Extend primer by 1-2 bp to capture G or C at 3\' end',
        expectedImprovement: 'GC clamp score +0.15 to +0.50',
      });
    }
  }

  // Terminal ΔG
  if (scores.terminal3DG < 0.7) {
    if (analysis3Prime.terminalDG > -6) {
      suggestions.push({
        priority: 'high',
        issue: '3\' end binding too weak',
        action: 'Extend primer to include more GC-rich sequence at 3\' end',
        expectedImprovement: 'Terminal ΔG score +0.20 to +0.30',
      });
    } else if (analysis3Prime.terminalDG < -14) {
      suggestions.push({
        priority: 'medium',
        issue: '3\' end binding too strong (mispriming risk)',
        action: 'Contract primer or shift position to reduce 3\' stability',
        expectedImprovement: 'Improved specificity',
      });
    }
  }

  // Tm
  if (scores.tm < 0.7) {
    if (primerScore.tm < 55) {
      suggestions.push({
        priority: 'medium',
        issue: 'Tm too low',
        action: 'Extend primer length to increase Tm',
        expectedImprovement: 'Tm score +0.10 to +0.20',
      });
    } else if (primerScore.tm > 65) {
      suggestions.push({
        priority: 'medium',
        issue: 'Tm too high',
        action: 'Reduce primer length or shift to AT-richer region',
        expectedImprovement: 'Tm score +0.10 to +0.20',
      });
    }
  }

  // Homopolymer
  if (scores.homopolymer < 0.8) {
    suggestions.push({
      priority: 'medium',
      issue: 'Homopolymer run detected',
      action: 'Shift primer position to avoid run, or change length',
      expectedImprovement: 'Reduced polymerase slippage risk',
    });
  }

  // Hairpin
  if (scores.hairpin < 0.7) {
    suggestions.push({
      priority: 'medium',
      issue: 'Stable hairpin structure possible',
      action: 'Shift primer position to disrupt hairpin-forming sequence',
      expectedImprovement: 'Hairpin score +0.15 to +0.30',
    });
  }

  // Sort by priority
  const priorityOrder: Record<string, number> = { high: 0, medium: 1, low: 2 };
  suggestions.sort((a, b) => priorityOrder[a.priority] - priorityOrder[b.priority]);

  return {
    currentScore: compositeScore,
    currentTier: qualityTier,
    canImprove: suggestions.length > 0,
    suggestions,
    summary: suggestions.length > 0
      ? `${suggestions.length} improvement(s) identified. Primary: ${suggestions[0].action}`
      : 'Primer design is already optimal',
  };
}

/**
 * Iteratively optimize a primer pair
 *
 * Main entry point for smart design. Takes initial primer candidates and
 * iteratively improves them through length adjustment and position optimization.
 */
export function optimizePrimerPair(
  seq: string,
  initialFwd: { seq: string; startPos: number; length: number },
  initialRev: { seq: string; startPos: number; length: number },
  fullTemplate: string,
  options: {
    maxIterations?: number;
    targetScore?: number;
    temperature?: number;
  } = {}
): OptimizedPair {
  const {
    maxIterations = SMART_DESIGN_CONSTRAINTS.maxIterations,
    targetScore = SMART_DESIGN_CONSTRAINTS.targetScore,
    temperature = 55,
  } = options;

  // Track optimization history
  const history: any[] = [];

  // Optimize forward primer
  const fwdOptimized = optimizePrimer(
    fullTemplate,
    initialFwd.startPos,
    initialFwd.length,
    true,
    options
  );

  // Optimize reverse primer (working on reverse complement)
  const revOptimized = optimizePrimer(
    fullTemplate,
    initialRev.startPos,
    initialRev.length,
    false,
    options
  );

  // Score optimized variants
  const fwdScore = scorePrimerVariant(fwdOptimized.seq, fullTemplate, true, temperature);
  const revScore = scorePrimerVariant(revOptimized.seq, fullTemplate, false, temperature);

  // Generate suggestions for further improvement
  const fwdSuggestions = generateDesignSuggestions(fwdScore);
  const revSuggestions = generateDesignSuggestions(revScore);

  // Calculate pair score (average with Tm diff penalty)
  const tmDiff = Math.abs(fwdScore.tm - revScore.tm);
  const tmDiffPenalty = tmDiff > 3 ? (tmDiff - 3) * 2 : 0;
  const pairScore = Math.round((fwdScore.compositeScore + revScore.compositeScore) / 2 - tmDiffPenalty);

  return {
    forward: {
      ...fwdOptimized,
      scoring: fwdScore,
      suggestions: fwdSuggestions,
    },
    reverse: {
      ...revOptimized,
      scoring: revScore,
      suggestions: revSuggestions,
    },
    pairScore,
    pairQuality: classifyQuality(pairScore),
    tmDiff: Math.round(tmDiff * 10) / 10,
    optimizationSummary: {
      fwdOptimized: fwdOptimized.optimized,
      revOptimized: revOptimized.optimized,
      fwdReason: fwdOptimized.reason,
      revReason: revOptimized.reason,
      improvements: [
        fwdOptimized.optimized ? `FWD: ${fwdOptimized.reason}` : null,
        revOptimized.optimized ? `REV: ${revOptimized.reason}` : null,
      ].filter(Boolean) as string[],
    },
    meetsTargetScore: pairScore >= targetScore,
  };
}

/**
 * Quick analysis of primer quality with improvement potential
 */
export function quickAssess(seq: string): QuickAssessment {
  const analysis = analyze3PrimeEnd(seq);
  const gcContent = (seq.match(/[GC]/gi) || []).length / seq.length;
  const primerTm = tm(seq, '', true);

  const issues = [...analysis.issues];

  // Check additional factors
  if (primerTm < 50) issues.push(`Tm too low: ${primerTm.toFixed(1)}°C`);
  if (primerTm > 68) issues.push(`Tm too high: ${primerTm.toFixed(1)}°C`);
  if (gcContent < 0.35) issues.push(`GC content low: ${(gcContent * 100).toFixed(0)}%`);
  if (gcContent > 0.70) issues.push(`GC content high: ${(gcContent * 100).toFixed(0)}%`);
  if (seq.length < 18) issues.push(`Length short: ${seq.length}bp`);
  if (seq.length > 28) issues.push(`Length long: ${seq.length}bp`);

  // Check for runs
  if (/(.)\1{3,}/.test(seq)) {
    issues.push('Contains homopolymer run (4+ bases)');
  }

  const qualityLevel =
    issues.length === 0 ? 'excellent' :
    issues.length <= 1 ? 'good' :
    issues.length <= 2 ? 'acceptable' :
    issues.length <= 3 ? 'marginal' : 'poor';

  return {
    seq,
    length: seq.length,
    tm: Math.round(primerTm * 10) / 10,
    gc: Math.round(gcContent * 100),
    gcClamp: analysis.gcCounts.last2,
    terminal3DG: analysis.terminalDG,
    quality: qualityLevel,
    issues,
    needsOptimization: issues.length > 0,
    canImproveWith: issues.length > 0 ? suggestOptimizationStrategy(issues) : null,
  };
}

/**
 * Suggest optimization strategy based on identified issues
 */
function suggestOptimizationStrategy(issues: string[]): string {
  const strategies: string[] = [];

  for (const issue of issues) {
    if (issue.includes('GC clamp') || issue.includes('0 G/C')) {
      strategies.push('extend_for_gc_clamp');
    }
    if (issue.includes('ΔG too weak') || issue.includes('Tm too low')) {
      strategies.push('extend_length');
    }
    if (issue.includes('ΔG very strong') || issue.includes('Tm too high')) {
      strategies.push('contract_length');
    }
    if (issue.includes('Poly-A/T') || issue.includes('AT-rich')) {
      strategies.push('shift_position');
    }
    if (issue.includes('homopolymer')) {
      strategies.push('shift_position');
    }
  }

  // Deduplicate and prioritize
  const uniqueStrategies = [...new Set(strategies)];

  if (uniqueStrategies.length === 0) return 'manual_review';
  if (uniqueStrategies.length === 1) return uniqueStrategies[0];

  // If conflicting (extend vs contract), prefer shift
  if (uniqueStrategies.includes('extend_length') && uniqueStrategies.includes('contract_length')) {
    return 'shift_position';
  }

  return uniqueStrategies[0];
}
