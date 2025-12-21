/**
 * Consolidated Overhang Validation Module
 *
 * This module provides unified validation for Golden Gate overhangs,
 * consolidating logic from:
 * - overhang-optimizer.ts (filterOverhangs)
 * - auto-domestication-optimizer.ts (validateDomesticationOverhangs)
 * - goldengate.ts (inline checks in findOptimalOverhangSet*)
 *
 * Single source of truth for all overhang validation.
 */

import {
  calculateFidelity,
  calculateEfficiency,
  reverseComplement,
  isPalindrome,
  isHomopolymer,
  getEnzymeLigationData,
  type FidelityResult,
  type EfficiencyResult,
  type LigationMatrix,
} from './fidelity-core.js';

// ============================================================================
// TYPES
// ============================================================================

export type ValidationSeverity = 'error' | 'warning' | 'info';

export interface ValidationIssue {
  type: string;
  severity: ValidationSeverity;
  message: string;
  penalty?: number;
}

export interface OverhangValidation {
  overhang: string;
  isValid: boolean;
  score: number;  // 0-100 quality score
  issues: ValidationIssue[];
  efficiency: EfficiencyResult;
  details: {
    isPalindrome: boolean;
    isHomopolymer: boolean;
    gcContent: number;
    gcCount: number;
    hasLowCorrectFreq: boolean;
    correctFreq: number;
  };
}

export interface SetValidation {
  overhangs: string[];
  individual: OverhangValidation[];
  isValid: boolean;
  setScore: number;
  fidelity: FidelityResult;
  issues: ValidationIssue[];
  conflicts: OverhangConflict[];
  recommendation: string;
}

export interface OverhangConflict {
  type: 'duplicate' | 'reverse_complement' | 'cross_ligation';
  overhang1: string;
  overhang2: string;
  severity: ValidationSeverity;
  message: string;
  crossLigationRatio?: number;
}

export interface ValidationContext {
  /** Existing overhangs in the set (for conflict detection) */
  existingOverhangs?: string[];

  /** Enzyme for ligation data lookup */
  enzyme?: string;

  /** Ligation matrix for cross-ligation check */
  matrix?: LigationMatrix;

  /** Check for cross-ligation potential */
  checkCrossLigation?: boolean;

  /** Minimum correct frequency threshold */
  minCorrectFreq?: number;

  /** Cross-ligation ratio threshold for warning */
  crossLigationThreshold?: number;

  /** Required patterns (e.g., '..TG' for overhangs ending in TG) */
  requiredPatterns?: string[];

  /** Excluded overhangs */
  excludeOverhangs?: string[];
}

export interface ValidationOptions {
  /** Minimum score to consider overhang valid (default: 50) */
  minValidScore?: number;

  /** Whether to treat warnings as errors (default: false) */
  strictMode?: boolean;

  /** Minimum correct ligation frequency (default: 300) */
  minCorrectFreq?: number;

  /** Cross-ligation warning threshold (default: 0.05 = 5%) */
  crossLigationThreshold?: number;
}

// ============================================================================
// CONSTANTS
// ============================================================================

export const DEFAULT_VALIDATION_OPTIONS: Required<ValidationOptions> = {
  minValidScore: 50,
  strictMode: false,
  minCorrectFreq: 300,
  crossLigationThreshold: 0.05,
};

/** Score penalties for different issues */
const PENALTIES = {
  palindrome: 100,        // Invalid - self-complementary
  homopolymer: 30,        // Reduced efficiency
  extremeGC: 15,          // 0% or 100% GC
  lowGC: 8,               // <25% GC
  highGC: 8,              // >75% GC
  lowCorrectFreq: 20,     // Low ligation frequency
  duplicate: 100,         // Same as existing overhang
  reverseComplement: 100, // RC conflict with existing
  crossLigationHigh: 25,  // >10% cross-ligation
  crossLigationMed: 15,   // 5-10% cross-ligation
  crossLigationLow: 5,    // 1-5% cross-ligation
  patternMismatch: 100,   // Doesn't match required pattern
  excluded: 100,          // Explicitly excluded
};

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Count G and C bases
 */
function countGC(seq: string): number {
  return (seq.match(/[GC]/gi) || []).length;
}

/**
 * Check if overhang matches a pattern
 * Pattern uses '.' as wildcard, e.g., '..TG' matches any overhang ending in TG
 */
function matchesPattern(overhang: string, pattern: string): boolean {
  if (!pattern) return true;
  if (pattern.length !== overhang.length) return false;

  for (let i = 0; i < pattern.length; i++) {
    if (pattern[i] !== '.' && pattern[i].toUpperCase() !== overhang[i].toUpperCase()) {
      return false;
    }
  }
  return true;
}

/**
 * Get correct ligation frequency from matrix
 */
function getCorrectFreq(overhang: string, matrix: LigationMatrix | undefined): number {
  if (!matrix) return 0;
  const oh = overhang.toUpperCase();
  const wc = reverseComplement(oh);
  return matrix[oh]?.[wc] || 0;
}

/**
 * Calculate cross-ligation ratio with another overhang
 */
function getCrossLigationRatio(
  overhang: string,
  otherOverhang: string,
  matrix: LigationMatrix | undefined
): number {
  if (!matrix) return 0;

  const oh = overhang.toUpperCase();
  const otherWc = reverseComplement(otherOverhang.toUpperCase());
  const correctFreq = getCorrectFreq(oh, matrix);

  if (correctFreq === 0) return 0;

  const crossFreq = matrix[oh]?.[otherWc] || 0;
  return crossFreq / correctFreq;
}

// ============================================================================
// SINGLE OVERHANG VALIDATION
// ============================================================================

/**
 * Validate a single overhang
 *
 * @param overhang - 4bp overhang sequence
 * @param context - Validation context
 * @param options - Validation options
 * @returns Validation result with score and issues
 */
export function validateOverhang(
  overhang: string,
  context: ValidationContext = {},
  options: ValidationOptions = {}
): OverhangValidation {
  const opts = { ...DEFAULT_VALIDATION_OPTIONS, ...options };
  const oh = (overhang || '').toUpperCase();

  // Initialize result
  const issues: ValidationIssue[] = [];
  let score = 100;

  // Get efficiency analysis
  const efficiency = calculateEfficiency(oh);

  // Basic validation: length
  if (oh.length !== 4) {
    issues.push({
      type: 'invalid_length',
      severity: 'error',
      message: `Overhang must be 4bp, got ${oh.length}bp`,
      penalty: 100,
    });
    score = 0;
  }

  // Check for invalid characters
  if (!/^[ATGC]+$/i.test(oh)) {
    issues.push({
      type: 'invalid_chars',
      severity: 'error',
      message: 'Overhang contains invalid characters (only ATGC allowed)',
      penalty: 100,
    });
    score = 0;
  }

  // If basic validation failed, return early
  if (score === 0) {
    return {
      overhang: oh,
      isValid: false,
      score: 0,
      issues,
      efficiency,
      details: {
        isPalindrome: false,
        isHomopolymer: false,
        gcContent: 0,
        gcCount: 0,
        hasLowCorrectFreq: false,
        correctFreq: 0,
      },
    };
  }

  // Check palindrome (self-complementary) - CRITICAL
  const palindrome = isPalindrome(oh);
  if (palindrome) {
    issues.push({
      type: 'palindrome',
      severity: 'error',
      message: 'Self-complementary overhang - will cause self-ligation',
      penalty: PENALTIES.palindrome,
    });
    score -= PENALTIES.palindrome;
  }

  // Check homopolymer
  const homopolymer = isHomopolymer(oh);
  if (homopolymer) {
    issues.push({
      type: 'homopolymer',
      severity: 'warning',
      message: 'All same base - poor ligation specificity',
      penalty: PENALTIES.homopolymer,
    });
    score -= PENALTIES.homopolymer;
  }

  // Check GC content
  const gcCount = countGC(oh);
  const gcContent = gcCount / 4;

  if (gcContent === 0) {
    issues.push({
      type: 'no_gc',
      severity: 'warning',
      message: '0% GC - weak Watson-Crick pairing',
      penalty: PENALTIES.extremeGC,
    });
    score -= PENALTIES.extremeGC;
  } else if (gcContent === 1) {
    issues.push({
      type: 'all_gc',
      severity: 'warning',
      message: '100% GC - slow melting kinetics',
      penalty: PENALTIES.extremeGC,
    });
    score -= PENALTIES.extremeGC;
  } else if (gcContent < 0.25) {
    issues.push({
      type: 'low_gc',
      severity: 'info',
      message: 'Low GC content (<25%)',
      penalty: PENALTIES.lowGC,
    });
    score -= PENALTIES.lowGC;
  } else if (gcContent > 0.75) {
    issues.push({
      type: 'high_gc',
      severity: 'info',
      message: 'High GC content (>75%)',
      penalty: PENALTIES.highGC,
    });
    score -= PENALTIES.highGC;
  }

  // Check ligation frequency (if matrix available)
  const correctFreq = getCorrectFreq(oh, context.matrix);
  const hasLowCorrectFreq = context.matrix && correctFreq < opts.minCorrectFreq;

  if (hasLowCorrectFreq) {
    issues.push({
      type: 'low_correct_freq',
      severity: 'warning',
      message: `Low correct ligation frequency: ${correctFreq} (min: ${opts.minCorrectFreq})`,
      penalty: PENALTIES.lowCorrectFreq,
    });
    score -= PENALTIES.lowCorrectFreq;
  }

  // Check against excluded list
  if (context.excludeOverhangs) {
    const excludeSet = new Set(context.excludeOverhangs.map(e => e.toUpperCase()));
    const rc = reverseComplement(oh);

    if (excludeSet.has(oh) || excludeSet.has(rc)) {
      issues.push({
        type: 'excluded',
        severity: 'error',
        message: 'Overhang is in exclusion list',
        penalty: PENALTIES.excluded,
      });
      score -= PENALTIES.excluded;
    }
  }

  // Check required patterns
  if (context.requiredPatterns && context.requiredPatterns.length > 0) {
    const matchesAny = context.requiredPatterns.some(p => matchesPattern(oh, p));
    if (!matchesAny) {
      issues.push({
        type: 'pattern_mismatch',
        severity: 'error',
        message: `Overhang doesn't match required patterns: ${context.requiredPatterns.join(', ')}`,
        penalty: PENALTIES.patternMismatch,
      });
      score -= PENALTIES.patternMismatch;
    }
  }

  // Check against existing overhangs (duplicates and RC conflicts)
  if (context.existingOverhangs && context.existingOverhangs.length > 0) {
    const rc = reverseComplement(oh);

    for (const existing of context.existingOverhangs) {
      const existingUpper = existing.toUpperCase();

      if (existingUpper === oh) {
        issues.push({
          type: 'duplicate',
          severity: 'error',
          message: `Duplicate overhang: ${oh}`,
          penalty: PENALTIES.duplicate,
        });
        score -= PENALTIES.duplicate;
        break;
      }

      if (existingUpper === rc) {
        issues.push({
          type: 'reverse_complement',
          severity: 'error',
          message: `RC conflict with ${existing}: would create palindromic junction`,
          penalty: PENALTIES.reverseComplement,
        });
        score -= PENALTIES.reverseComplement;
        break;
      }
    }
  }

  // Check cross-ligation potential
  if (context.checkCrossLigation && context.existingOverhangs && context.matrix) {
    for (const existing of context.existingOverhangs) {
      const ratio = getCrossLigationRatio(oh, existing, context.matrix);

      if (ratio >= 0.10) {
        issues.push({
          type: 'cross_ligation_high',
          severity: 'warning',
          message: `High cross-ligation risk with ${existing}: ${(ratio * 100).toFixed(1)}%`,
          penalty: PENALTIES.crossLigationHigh,
        });
        score -= PENALTIES.crossLigationHigh;
      } else if (ratio >= 0.05) {
        issues.push({
          type: 'cross_ligation_med',
          severity: 'info',
          message: `Medium cross-ligation risk with ${existing}: ${(ratio * 100).toFixed(1)}%`,
          penalty: PENALTIES.crossLigationMed,
        });
        score -= PENALTIES.crossLigationMed;
      } else if (ratio >= 0.01) {
        issues.push({
          type: 'cross_ligation_low',
          severity: 'info',
          message: `Low cross-ligation with ${existing}: ${(ratio * 100).toFixed(1)}%`,
          penalty: PENALTIES.crossLigationLow,
        });
        score -= PENALTIES.crossLigationLow;
      }
    }
  }

  // Ensure score is in valid range
  score = Math.max(0, Math.min(100, score));

  // Determine validity
  const hasErrors = issues.some(i => i.severity === 'error');
  const isValid = !hasErrors && score >= opts.minValidScore;

  return {
    overhang: oh,
    isValid,
    score,
    issues,
    efficiency,
    details: {
      isPalindrome: palindrome,
      isHomopolymer: homopolymer,
      gcContent,
      gcCount,
      hasLowCorrectFreq: hasLowCorrectFreq || false,
      correctFreq,
    },
  };
}

// ============================================================================
// SET VALIDATION
// ============================================================================

/**
 * Validate a set of overhangs for use together
 *
 * @param overhangs - Array of 4bp overhang sequences
 * @param enzyme - Enzyme name for fidelity calculation
 * @param options - Validation options
 * @returns Set validation result
 */
export function validateOverhangSet(
  overhangs: string[],
  enzyme: string = 'BsaI',
  options: ValidationOptions = {}
): SetValidation {
  const opts = { ...DEFAULT_VALIDATION_OPTIONS, ...options };

  if (!Array.isArray(overhangs) || overhangs.length === 0) {
    return {
      overhangs: [],
      individual: [],
      isValid: false,
      setScore: 0,
      fidelity: calculateFidelity([], enzyme),
      issues: [{ type: 'empty', severity: 'error', message: 'No overhangs provided' }],
      conflicts: [],
      recommendation: 'Provide at least one overhang',
    };
  }

  // Get ligation matrix
  const enzymeData = getEnzymeLigationData(enzyme);
  const matrix = enzymeData?.matrix;

  // Validate each overhang individually, with context of other overhangs
  const individual: OverhangValidation[] = [];
  const conflicts: OverhangConflict[] = [];
  const setIssues: ValidationIssue[] = [];

  for (let i = 0; i < overhangs.length; i++) {
    const oh = overhangs[i];
    const otherOverhangs = overhangs.filter((_, j) => j !== i);

    const validation = validateOverhang(oh, {
      existingOverhangs: otherOverhangs,
      enzyme,
      matrix,
      checkCrossLigation: true,
      minCorrectFreq: opts.minCorrectFreq,
    }, opts);

    individual.push(validation);

    // Extract conflicts for easier analysis
    for (const issue of validation.issues) {
      if (issue.type === 'duplicate') {
        conflicts.push({
          type: 'duplicate',
          overhang1: oh,
          overhang2: otherOverhangs.find(o => o.toUpperCase() === oh.toUpperCase()) || '',
          severity: 'error',
          message: issue.message,
        });
      } else if (issue.type === 'reverse_complement') {
        const rc = reverseComplement(oh);
        conflicts.push({
          type: 'reverse_complement',
          overhang1: oh,
          overhang2: otherOverhangs.find(o => o.toUpperCase() === rc) || '',
          severity: 'error',
          message: issue.message,
        });
      } else if (issue.type.startsWith('cross_ligation')) {
        const matchResult = issue.message.match(/with (\w+):/);
        const other = matchResult ? matchResult[1] : '';
        const ratioMatch = issue.message.match(/([\d.]+)%/);
        const ratio = ratioMatch ? parseFloat(ratioMatch[1]) / 100 : 0;

        conflicts.push({
          type: 'cross_ligation',
          overhang1: oh,
          overhang2: other,
          severity: issue.severity,
          message: issue.message,
          crossLigationRatio: ratio,
        });
      }
    }
  }

  // Calculate fidelity for the set
  const fidelity = calculateFidelity(overhangs, enzyme);

  // Check for set-level issues
  if (fidelity.finalFidelity < 0.90) {
    setIssues.push({
      type: 'low_fidelity',
      severity: 'warning',
      message: `Assembly fidelity ${fidelity.finalFidelityPercent} is below 90%`,
    });
  }

  if (fidelity.hasGTRisks) {
    setIssues.push({
      type: 'gt_mismatch_risk',
      severity: 'warning',
      message: `G:T mismatch risks detected for ${fidelity.gtRisks.length} overhang pair(s)`,
    });
  }

  // Calculate set score (average of individual scores, penalized by conflicts)
  const avgScore = individual.reduce((sum, v) => sum + v.score, 0) / individual.length;
  const conflictPenalty = conflicts.filter(c => c.severity === 'error').length * 20;
  const setScore = Math.max(0, avgScore - conflictPenalty);

  // Determine overall validity
  const hasErrors = individual.some(v => !v.isValid) ||
                    conflicts.some(c => c.severity === 'error');
  const isValid = !hasErrors && setScore >= opts.minValidScore;

  // Generate recommendation
  // Check for palindrome issues in individual validations
  const hasPalindromes = individual.some(v => v.details.isPalindrome);

  let recommendation = '';
  if (isValid && fidelity.finalFidelity >= 0.95) {
    recommendation = 'Excellent overhang set - high fidelity expected';
  } else if (isValid && fidelity.finalFidelity >= 0.90) {
    recommendation = 'Good overhang set - acceptable fidelity';
  } else if (hasPalindromes) {
    recommendation = 'Replace palindromic overhangs - they cause self-ligation';
  } else if (conflicts.some(c => c.type === 'reverse_complement')) {
    recommendation = 'Remove RC conflicts - overhangs and their RCs cannot be in same set';
  } else if (fidelity.finalFidelity < 0.90) {
    recommendation = 'Consider replacing low-fidelity overhangs to improve assembly success';
  } else {
    recommendation = 'Review and address validation issues';
  }

  return {
    overhangs: overhangs.map(o => o.toUpperCase()),
    individual,
    isValid,
    setScore,
    fidelity,
    issues: [...setIssues, ...individual.flatMap(v => v.issues)],
    conflicts: deduplicateConflicts(conflicts),
    recommendation,
  };
}

/**
 * Remove duplicate conflicts (A-B and B-A are the same)
 */
function deduplicateConflicts(conflicts: OverhangConflict[]): OverhangConflict[] {
  const seen = new Set<string>();
  const result: OverhangConflict[] = [];

  for (const conflict of conflicts) {
    const key = [conflict.overhang1, conflict.overhang2].sort().join('-') + conflict.type;
    if (!seen.has(key)) {
      seen.add(key);
      result.push(conflict);
    }
  }

  return result;
}

// ============================================================================
// CANDIDATE FILTERING
// ============================================================================

/**
 * Filter overhang candidates based on validation criteria
 *
 * @param candidates - Array of candidate overhangs
 * @param context - Validation context
 * @param options - Validation options
 * @returns Filtered and scored candidates
 */
export function filterOverhangCandidates(
  candidates: string[],
  context: ValidationContext = {},
  options: ValidationOptions = {}
): Array<{ overhang: string; validation: OverhangValidation }> {
  const opts = { ...DEFAULT_VALIDATION_OPTIONS, ...options };

  // Get enzyme data for matrix
  const enzyme = context.enzyme || 'BsaI';
  const enzymeData = getEnzymeLigationData(enzyme);
  const matrix = enzymeData?.matrix || context.matrix;

  // Validate each candidate
  const validated = candidates.map(oh => ({
    overhang: oh.toUpperCase(),
    validation: validateOverhang(oh, { ...context, matrix }, opts),
  }));

  // Filter to valid candidates only
  const valid = validated.filter(v => v.validation.isValid);

  // Sort by score (descending)
  valid.sort((a, b) => b.validation.score - a.validation.score);

  return valid;
}

/**
 * Check if an overhang has zero cross-ligation with a set
 */
export function hasZeroCrossLigation(
  overhang: string,
  existingSet: string[],
  matrix: LigationMatrix | undefined
): boolean {
  if (!matrix) return true; // Assume OK if no matrix

  const oh = overhang.toUpperCase();
  const ohWc = reverseComplement(oh);

  for (const existing of existingSet) {
    const existingUpper = existing.toUpperCase();
    const existingWc = reverseComplement(existingUpper);

    // Check oh → existingWc (oh would ligate with existing's partner)
    if ((matrix[oh]?.[existingWc] || 0) > 0) return false;

    // Check existing → ohWc (existing would ligate with oh's partner)
    if ((matrix[existingUpper]?.[ohWc] || 0) > 0) return false;
  }

  return true;
}

/**
 * Check if a set has zero cross-ligation between all pairs
 */
export function setHasZeroCrossLigation(
  overhangs: string[],
  matrix: LigationMatrix | undefined
): boolean {
  if (!matrix) return true;

  for (let i = 0; i < overhangs.length; i++) {
    for (let j = i + 1; j < overhangs.length; j++) {
      const oh1 = overhangs[i].toUpperCase();
      const oh2 = overhangs[j].toUpperCase();
      const wc2 = reverseComplement(oh2);

      if ((matrix[oh1]?.[wc2] || 0) > 0) return false;
    }
  }

  return true;
}

// ============================================================================
// SITE RECREATION CHECK
// ============================================================================

/**
 * Recognition sites for common Golden Gate enzymes
 */
const ENZYME_RECOGNITION_SITES: Record<string, { recognition: string; overhangLen: number }> = {
  'BsaI': { recognition: 'GGTCTC', overhangLen: 4 },
  'BsaI-HFv2': { recognition: 'GGTCTC', overhangLen: 4 },
  'BsmBI': { recognition: 'CGTCTC', overhangLen: 4 },
  'BsmBI-v2': { recognition: 'CGTCTC', overhangLen: 4 },
  'BbsI': { recognition: 'GAAGAC', overhangLen: 4 },
  'BbsI-HF': { recognition: 'GAAGAC', overhangLen: 4 },
  'Esp3I': { recognition: 'CGTCTC', overhangLen: 4 },
  'SapI': { recognition: 'GCTCTTC', overhangLen: 3 },
};

export interface SiteRecreationResult {
  /** Whether the junction would recreate a recognition site */
  recreatesSite: boolean;

  /** Direction of the recreated site (if any) */
  direction?: 'forward' | 'reverse';

  /** Position in the scar context where site would appear */
  position?: number;

  /** The junction/scar context that was analyzed */
  context: string;

  /** The specific recognition site found */
  siteFound?: string;

  /** Risk level */
  risk: 'none' | 'low' | 'medium' | 'high';

  /** Recommendation */
  recommendation: string;
}

/**
 * Comprehensive check if assembling with a junction would recreate a recognition site
 *
 * This checks the full junction context (upstream + overhang + downstream)
 * for any occurrence of the enzyme's recognition site.
 *
 * @param upstreamSeq - Sequence upstream of the junction
 * @param downstreamSeq - Sequence downstream of the junction
 * @param overhang - The 4bp overhang at the junction
 * @param enzyme - Enzyme name
 * @returns Comprehensive site recreation analysis
 */
export function checkSiteRecreation(
  upstreamSeq: string,
  downstreamSeq: string,
  overhang: string,
  enzyme: string
): SiteRecreationResult {
  const enzymeInfo = ENZYME_RECOGNITION_SITES[enzyme] || ENZYME_RECOGNITION_SITES['BsaI'];
  const recognition = enzymeInfo.recognition.toUpperCase();
  const recognitionRC = reverseComplement(recognition);

  // Build context: need enough sequence to span a recognition site
  const contextSize = recognition.length + 4;
  const upstream = (upstreamSeq || '').toUpperCase().slice(-contextSize);
  const downstream = (downstreamSeq || '').toUpperCase().slice(0, contextSize);
  const oh = (overhang || '').toUpperCase();

  // The assembled scar region
  // In Golden Gate, the overhang is shared between fragments
  // The scar context is: [end of upstream] + [overhang] + [start of downstream]
  const scarContext = upstream + oh + downstream;

  // Check for forward recognition site
  const fwdMatch = scarContext.indexOf(recognition);

  // Check for reverse recognition site
  const revMatch = scarContext.indexOf(recognitionRC);

  if (fwdMatch >= 0) {
    // Determine risk level based on position
    // If site is entirely in upstream or downstream, it's existing (not created by junction)
    const siteStart = fwdMatch;
    const siteEnd = fwdMatch + recognition.length;
    const ohStart = upstream.length;
    const ohEnd = upstream.length + oh.length;

    // High risk: site spans the junction
    const spansJunction = siteStart < ohEnd && siteEnd > ohStart;

    return {
      recreatesSite: true,
      direction: 'forward',
      position: fwdMatch,
      context: scarContext,
      siteFound: recognition,
      risk: spansJunction ? 'high' : 'medium',
      recommendation: spansJunction
        ? 'Junction would recreate recognition site - choose different junction position'
        : 'Recognition site near junction - verify it does not span the junction',
    };
  }

  if (revMatch >= 0) {
    const siteStart = revMatch;
    const siteEnd = revMatch + recognitionRC.length;
    const ohStart = upstream.length;
    const ohEnd = upstream.length + oh.length;
    const spansJunction = siteStart < ohEnd && siteEnd > ohStart;

    return {
      recreatesSite: true,
      direction: 'reverse',
      position: revMatch,
      context: scarContext,
      siteFound: recognitionRC,
      risk: spansJunction ? 'high' : 'medium',
      recommendation: spansJunction
        ? 'Junction would recreate reverse recognition site - choose different junction position'
        : 'Reverse recognition site near junction - verify it does not span the junction',
    };
  }

  // Check for partial site at junction boundaries (could create site with next fragment)
  let partialRisk: 'none' | 'low' = 'none';
  let partialRecommendation = 'No site recreation detected - junction is safe';

  // Check if end of context could form site with next fragment
  for (let i = 2; i < recognition.length; i++) {
    const ending = scarContext.slice(-i);
    if (recognition.startsWith(ending) || recognitionRC.startsWith(ending)) {
      partialRisk = 'low';
      partialRecommendation = 'Partial recognition site at junction end - verify next fragment does not complete it';
      break;
    }
  }

  return {
    recreatesSite: false,
    context: scarContext,
    risk: partialRisk,
    recommendation: partialRecommendation,
  };
}

/**
 * Batch check multiple junctions for site recreation
 */
export function checkMultipleJunctions(
  junctions: Array<{
    upstreamSeq: string;
    downstreamSeq: string;
    overhang: string;
  }>,
  enzyme: string
): {
  results: SiteRecreationResult[];
  hasHighRisk: boolean;
  hasMediumRisk: boolean;
  summary: string;
} {
  const results = junctions.map(j =>
    checkSiteRecreation(j.upstreamSeq, j.downstreamSeq, j.overhang, enzyme)
  );

  const hasHighRisk = results.some(r => r.risk === 'high');
  const hasMediumRisk = results.some(r => r.risk === 'medium');

  let summary = 'All junctions are safe';
  if (hasHighRisk) {
    const count = results.filter(r => r.risk === 'high').length;
    summary = `${count} junction(s) would recreate recognition site - must fix before assembly`;
  } else if (hasMediumRisk) {
    const count = results.filter(r => r.risk === 'medium').length;
    summary = `${count} junction(s) have recognition sites nearby - review carefully`;
  }

  return {
    results,
    hasHighRisk,
    hasMediumRisk,
    summary,
  };
}

// ============================================================================
// EXPORTS
// ============================================================================

export { reverseComplement, isPalindrome, isHomopolymer };
