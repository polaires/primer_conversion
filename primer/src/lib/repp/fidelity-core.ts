/**
 * Unified Fidelity Calculation Module
 *
 * Consolidates all fidelity calculations into a single source of truth:
 * - Matrix-based experimental fidelity (Pryor et al. 2020)
 * - Static fidelity fallback
 * - G:T mismatch risk analysis
 * - Efficiency penalties (TNNA, palindromes, etc.)
 *
 * This module is the ONLY place where fidelity should be calculated.
 * All other modules should import from here.
 */

import ligationData from './ligation-data.json';
import {
  calculateEfficiency,
  calculateSetEfficiency,
  isPalindrome,
  isHomopolymer,
  type EfficiencyResult,
  type SetEfficiencyResult,
} from './overhang-efficiency.js';

// ============================================================================
// TYPES
// ============================================================================

export interface LigationMatrix {
  [overhang: string]: {
    [partner: string]: number;
  };
}

export interface EnzymeLigationData {
  overhangLength: number;
  overhangs: string[];
  matrix: LigationMatrix;
  optimalSets?: {
    [size: number]: {
      overhangs: string[];
      fidelity: number;
      lowestJunction: { overhang: string; fidelity: number };
    };
  };
}

export interface JunctionFidelity {
  overhang: string;
  wcPartner: string;
  fidelity: number;
  fidelityPercent: string;
  correctFreq: number;
  totalFreq: number;
  efficiency?: EfficiencyResult;
  error?: string;
}

export interface GTMismatchRisk {
  overhang1: string;
  overhang2: string;
  reverseComplement2: string;
  wobblePositions: number[];
  wobbleCount: number;
  matchCount: number;
  positionWeight: number;
  risk: 'critical' | 'high' | 'medium';
  expectedMisLigation: number;
  expectedMisLigationPercent: string;
  description: string;
}

export interface FidelityResult {
  enzyme: string;
  enzymeFullName: string;
  overhangs: string[];
  numJunctions: number;

  // Core fidelity metrics
  assemblyFidelity: number;
  assemblyFidelityPercent: string;

  // Per-junction breakdown
  junctions: JunctionFidelity[];
  sortedJunctions: JunctionFidelity[];
  lowestFidelity: JunctionFidelity;

  // Efficiency analysis
  efficiency: SetEfficiencyResult;
  efficiencyAdjustedFidelity: number;
  efficiencyAdjustedFidelityPercent: string;

  // G:T mismatch analysis
  gtRisks: GTMismatchRisk[];
  hasGTRisks: boolean;
  gtAdjustedFidelity: number;
  gtAdjustedFidelityPercent: string;

  // Final combined score
  finalFidelity: number;
  finalFidelityPercent: string;

  // Metadata
  source: 'matrix' | 'static' | 'fallback';
  warnings: string[];
  calculationMethod: string;
}

export interface FidelityOptions {
  /** Include G:T mismatch risk analysis (default: true) */
  includeGTRisks?: boolean;

  /** Include efficiency penalties (default: true) */
  includeEfficiency?: boolean;

  /** Use experimental matrix data if available (default: true) */
  useMatrixData?: boolean;

  /** Fallback fidelity when no data available (default: 0.85) */
  fallbackFidelity?: number;

  /** G:T mismatch factor (default: 0.20) */
  gtMismatchFactor?: number;

  /** Minimum matches to flag G:T risk (default: 3) */
  gtRiskMatchThreshold?: number;
}

// ============================================================================
// CONSTANTS
// ============================================================================

/** Enzyme full names */
export const ENZYME_FULL_NAMES: Record<string, string> = {
  'BsaI': 'BsaI-HFv2',
  'BsaI-HFv2': 'BsaI-HFv2',
  'BsmBI': 'BsmBI-v2',
  'BsmBI-v2': 'BsmBI-v2',
  'BbsI': 'BbsI-HF',
  'BbsI-HF': 'BbsI-HF',
  'Esp3I': 'Esp3I',
  'SapI': 'SapI',
};

/** Enzymes with experimental ligation data */
export const ENZYMES_WITH_DATA = ['BsaI-HFv2', 'BsmBI-v2', 'BbsI-HF', 'Esp3I', 'SapI'];

/** Default fidelity options */
export const DEFAULT_FIDELITY_OPTIONS: Required<FidelityOptions> = {
  includeGTRisks: true,
  includeEfficiency: true,
  useMatrixData: true,
  fallbackFidelity: 0.85,
  gtMismatchFactor: 0.20,
  gtRiskMatchThreshold: 3,
};

/** Static fidelity data for fallback */
export const STATIC_FIDELITY: Record<string, { fidelity: number; category: string }> = {
  // High fidelity overhangs (>98%)
  'GGAG': { fidelity: 0.99, category: 'excellent' },
  'TACT': { fidelity: 0.99, category: 'excellent' },
  'AATG': { fidelity: 0.98, category: 'excellent' },
  'AGGT': { fidelity: 0.98, category: 'excellent' },
  'GCTT': { fidelity: 0.98, category: 'excellent' },
  'CGCT': { fidelity: 0.98, category: 'excellent' },
  'TGCC': { fidelity: 0.98, category: 'excellent' },
  'ACTA': { fidelity: 0.98, category: 'excellent' },
  'TTCG': { fidelity: 0.98, category: 'excellent' },
  'GCTG': { fidelity: 0.98, category: 'excellent' },
  // Good fidelity overhangs (95-98%)
  'CAGA': { fidelity: 0.97, category: 'good' },
  'TCGA': { fidelity: 0.97, category: 'good' },
  'GTGC': { fidelity: 0.96, category: 'good' },
  'CTAC': { fidelity: 0.96, category: 'good' },
  'GAGT': { fidelity: 0.96, category: 'good' },
  'ATCC': { fidelity: 0.95, category: 'good' },
  'CCGA': { fidelity: 0.95, category: 'good' },
  // Medium fidelity (90-95%)
  'GAAG': { fidelity: 0.93, category: 'medium' },
  'TGGA': { fidelity: 0.92, category: 'medium' },
  'CAGG': { fidelity: 0.91, category: 'medium' },
  'AGCC': { fidelity: 0.90, category: 'medium' },
  // Lower fidelity (<90%)
  'TAAA': { fidelity: 0.85, category: 'low' },
  'TTTA': { fidelity: 0.85, category: 'low' },
  'AAAA': { fidelity: 0.75, category: 'poor' },
  'TTTT': { fidelity: 0.75, category: 'poor' },
};

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Calculate reverse complement of a DNA sequence
 */
export function reverseComplement(seq: string): string {
  const complement: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
  return seq
    .toUpperCase()
    .split('')
    .reverse()
    .map(base => complement[base] || base)
    .join('');
}

/**
 * Normalize enzyme name to match ligation data keys
 */
function normalizeEnzymeName(enzyme: string): string {
  const mapping: Record<string, string> = {
    'BsaI': 'BsaI-HFv2',
    'BsmBI': 'BsmBI-v2',
    'BbsI': 'BbsI-HF',
  };
  return mapping[enzyme] || enzyme;
}

/**
 * Get ligation data for an enzyme
 */
export function getEnzymeLigationData(enzyme: string): EnzymeLigationData | null {
  const normalizedName = normalizeEnzymeName(enzyme);
  const data = (ligationData as any).enzymes?.[normalizedName];

  if (!data) {
    return null;
  }

  return data as EnzymeLigationData;
}

/**
 * Check if two bases can pair (including G:T wobbles)
 */
function canPair(base1: string, base2: string, gtMismatchFactor: number): {
  pairs: boolean;
  isWobble: boolean;
  factor: number;
} {
  const wc: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
  const b1 = (base1 || '').toUpperCase();
  const b2 = (base2 || '').toUpperCase();

  if (wc[b1] === b2) {
    return { pairs: true, isWobble: false, factor: 1.0 };
  }

  // G:T wobble pairs
  if ((b1 === 'G' && b2 === 'T') || (b1 === 'T' && b2 === 'G')) {
    return { pairs: true, isWobble: true, factor: gtMismatchFactor };
  }

  return { pairs: false, isWobble: false, factor: 0 };
}

// ============================================================================
// CORE FIDELITY CALCULATION
// ============================================================================

/**
 * Calculate G:T mismatch risks for an overhang set
 */
export function findGTMismatchRisks(
  overhangs: string[],
  options: Pick<FidelityOptions, 'gtMismatchFactor' | 'gtRiskMatchThreshold'> = {}
): GTMismatchRisk[] {
  const {
    gtMismatchFactor = DEFAULT_FIDELITY_OPTIONS.gtMismatchFactor,
    gtRiskMatchThreshold = DEFAULT_FIDELITY_OPTIONS.gtRiskMatchThreshold,
  } = options;

  if (!Array.isArray(overhangs) || overhangs.length < 2) {
    return [];
  }

  const risks: GTMismatchRisk[] = [];

  for (let i = 0; i < overhangs.length; i++) {
    for (let j = i + 1; j < overhangs.length; j++) {
      const oh1 = (overhangs[i] || '').toUpperCase();
      const oh2 = (overhangs[j] || '').toUpperCase();

      if (!oh1 || !oh2) continue;

      const oh2rc = reverseComplement(oh2);

      let wobbleCount = 0;
      let matchCount = 0;
      const wobblePositions: number[] = [];

      for (let k = 0; k < oh1.length && k < oh2rc.length; k++) {
        const pair = canPair(oh1[k], oh2rc[k], gtMismatchFactor);
        if (pair.isWobble) {
          wobbleCount++;
          wobblePositions.push(k);
        }
        if (pair.pairs) matchCount++;
      }

      // Flag if 3+ matches with at least 1 G:T wobble
      if (matchCount >= gtRiskMatchThreshold && wobbleCount >= 1) {
        const expectedMisLigation = Math.pow(gtMismatchFactor, wobbleCount);
        const positionWeight = wobblePositions.reduce((sum, pos) => sum + (pos + 1), 0) / wobblePositions.length;
        const positionRisk = positionWeight / oh1.length;

        risks.push({
          overhang1: oh1,
          overhang2: oh2,
          reverseComplement2: oh2rc,
          wobblePositions,
          wobbleCount,
          matchCount,
          positionWeight,
          risk: matchCount === oh1.length ? 'critical' :
                matchCount >= gtRiskMatchThreshold && positionRisk > 0.5 ? 'high' : 'medium',
          expectedMisLigation,
          expectedMisLigationPercent: `${(expectedMisLigation * 100).toFixed(1)}%`,
          description: `${oh1} may mis-ligate with RC(${oh2})=${oh2rc} via ${wobbleCount} G:T wobble(s)`,
        });
      }
    }
  }

  return risks;
}

/**
 * Calculate fidelity using experimental ligation matrix
 */
function calculateMatrixFidelity(
  overhangs: string[],
  enzymeData: EnzymeLigationData,
  enzyme: string
): Omit<FidelityResult, 'efficiency' | 'efficiencyAdjustedFidelity' | 'efficiencyAdjustedFidelityPercent' |
                         'gtRisks' | 'hasGTRisks' | 'gtAdjustedFidelity' | 'gtAdjustedFidelityPercent' |
                         'finalFidelity' | 'finalFidelityPercent'> {
  const matrix = enzymeData.matrix;
  const junctions: JunctionFidelity[] = [];
  let assemblyFidelity = 1.0;
  const warnings: string[] = [];

  for (const oh of overhangs) {
    const ohUpper = oh.toUpperCase();
    const wcPartner = reverseComplement(ohUpper);

    // Get correct ligation frequency
    const correctFreq = matrix[ohUpper]?.[wcPartner] || 0;

    if (correctFreq === 0) {
      warnings.push(`No ligation data for overhang ${ohUpper}`);
      junctions.push({
        overhang: ohUpper,
        wcPartner,
        fidelity: 0,
        fidelityPercent: '0.0%',
        correctFreq: 0,
        totalFreq: 0,
        error: 'No experimental data',
      });
      continue;
    }

    // Calculate total competing ligation frequency
    let totalFreq = correctFreq;

    for (const otherOh of overhangs) {
      if (otherOh.toUpperCase() === ohUpper) continue;
      const otherWc = reverseComplement(otherOh.toUpperCase());
      totalFreq += matrix[ohUpper]?.[otherWc] || 0;
    }

    const junctionFidelity = totalFreq > 0 ? correctFreq / totalFreq : 0;
    assemblyFidelity *= junctionFidelity;

    junctions.push({
      overhang: ohUpper,
      wcPartner,
      fidelity: junctionFidelity,
      fidelityPercent: `${(junctionFidelity * 100).toFixed(1)}%`,
      correctFreq,
      totalFreq,
    });

    if (junctionFidelity < 0.95) {
      warnings.push(`Junction ${ohUpper} has ${(junctionFidelity * 100).toFixed(1)}% fidelity`);
    }
  }

  const sortedJunctions = [...junctions].sort((a, b) => a.fidelity - b.fidelity);
  const lowestFidelity = sortedJunctions[0] || junctions[0];

  return {
    enzyme,
    enzymeFullName: ENZYME_FULL_NAMES[enzyme] || enzyme,
    overhangs: overhangs.map(o => o.toUpperCase()),
    numJunctions: overhangs.length,
    assemblyFidelity,
    assemblyFidelityPercent: `${(assemblyFidelity * 100).toFixed(1)}%`,
    junctions,
    sortedJunctions,
    lowestFidelity,
    source: 'matrix',
    warnings,
    calculationMethod: 'experimental-ligation-matrix',
  };
}

/**
 * Calculate fidelity using static data (fallback)
 */
function calculateStaticFidelity(
  overhangs: string[],
  fallbackFidelity: number,
  enzyme: string
): Omit<FidelityResult, 'efficiency' | 'efficiencyAdjustedFidelity' | 'efficiencyAdjustedFidelityPercent' |
                         'gtRisks' | 'hasGTRisks' | 'gtAdjustedFidelity' | 'gtAdjustedFidelityPercent' |
                         'finalFidelity' | 'finalFidelityPercent'> {
  const junctions: JunctionFidelity[] = [];
  let assemblyFidelity = 1.0;
  const warnings: string[] = [];

  for (const oh of overhangs) {
    const ohUpper = oh.toUpperCase();
    const wcPartner = reverseComplement(ohUpper);
    const staticData = STATIC_FIDELITY[ohUpper];
    const junctionFidelity = staticData?.fidelity || fallbackFidelity;

    assemblyFidelity *= junctionFidelity;

    junctions.push({
      overhang: ohUpper,
      wcPartner,
      fidelity: junctionFidelity,
      fidelityPercent: `${(junctionFidelity * 100).toFixed(1)}%`,
      correctFreq: 0,
      totalFreq: 0,
    });

    if (!staticData) {
      warnings.push(`No fidelity data for ${ohUpper}, using fallback ${(fallbackFidelity * 100).toFixed(0)}%`);
    }
  }

  const sortedJunctions = [...junctions].sort((a, b) => a.fidelity - b.fidelity);
  const lowestFidelity = sortedJunctions[0] || junctions[0];

  return {
    enzyme,
    enzymeFullName: ENZYME_FULL_NAMES[enzyme] || enzyme,
    overhangs: overhangs.map(o => o.toUpperCase()),
    numJunctions: overhangs.length,
    assemblyFidelity,
    assemblyFidelityPercent: `${(assemblyFidelity * 100).toFixed(1)}%`,
    junctions,
    sortedJunctions,
    lowestFidelity,
    source: 'static',
    warnings,
    calculationMethod: 'static-fidelity-table',
  };
}

// ============================================================================
// MAIN EXPORT: UNIFIED FIDELITY CALCULATION
// ============================================================================

/**
 * Calculate assembly fidelity with unified interface
 *
 * This is the MAIN function that should be used for all fidelity calculations.
 * It consolidates matrix-based, static, efficiency, and G:T mismatch analysis.
 *
 * @param overhangs - Array of 4bp overhang sequences
 * @param enzyme - Enzyme name (e.g., 'BsaI', 'BsmBI')
 * @param options - Calculation options
 * @returns Comprehensive fidelity result
 */
export function calculateFidelity(
  overhangs: string[],
  enzyme: string = 'BsaI',
  options: FidelityOptions = {}
): FidelityResult {
  const opts = { ...DEFAULT_FIDELITY_OPTIONS, ...options };

  // Handle empty input
  if (!Array.isArray(overhangs) || overhangs.length === 0) {
    const emptyEfficiency = calculateSetEfficiency([]);
    return {
      enzyme,
      enzymeFullName: ENZYME_FULL_NAMES[enzyme] || enzyme,
      overhangs: [],
      numJunctions: 0,
      assemblyFidelity: 0,
      assemblyFidelityPercent: '0.0%',
      junctions: [],
      sortedJunctions: [],
      lowestFidelity: {} as JunctionFidelity,
      efficiency: emptyEfficiency,
      efficiencyAdjustedFidelity: 0,
      efficiencyAdjustedFidelityPercent: '0.0%',
      gtRisks: [],
      hasGTRisks: false,
      gtAdjustedFidelity: 0,
      gtAdjustedFidelityPercent: '0.0%',
      finalFidelity: 0,
      finalFidelityPercent: '0.0%',
      source: 'fallback',
      warnings: ['No overhangs provided'],
      calculationMethod: 'none',
    };
  }

  // Step 1: Calculate base fidelity (matrix or static)
  let baseFidelity: ReturnType<typeof calculateMatrixFidelity>;

  if (opts.useMatrixData) {
    const enzymeData = getEnzymeLigationData(enzyme);
    if (enzymeData) {
      baseFidelity = calculateMatrixFidelity(overhangs, enzymeData, enzyme);
    } else {
      baseFidelity = calculateStaticFidelity(overhangs, opts.fallbackFidelity, enzyme);
    }
  } else {
    baseFidelity = calculateStaticFidelity(overhangs, opts.fallbackFidelity, enzyme);
  }

  // Step 2: Calculate efficiency
  const efficiency = opts.includeEfficiency
    ? calculateSetEfficiency(overhangs)
    : calculateSetEfficiency([]);

  // Apply efficiency adjustment
  const efficiencyAdjustedFidelity = opts.includeEfficiency
    ? baseFidelity.assemblyFidelity * efficiency.averageEfficiency
    : baseFidelity.assemblyFidelity;

  // Step 3: Calculate G:T mismatch risks
  const gtRisks = opts.includeGTRisks
    ? findGTMismatchRisks(overhangs, {
        gtMismatchFactor: opts.gtMismatchFactor,
        gtRiskMatchThreshold: opts.gtRiskMatchThreshold,
      })
    : [];

  // Apply G:T penalty (only for static method - matrix already includes cross-ligation)
  let gtPenalty = 1.0;
  if (baseFidelity.source === 'static' && opts.includeGTRisks) {
    gtRisks.forEach(risk => {
      gtPenalty *= (1 - risk.expectedMisLigation);
    });
  }

  const gtAdjustedFidelity = baseFidelity.assemblyFidelity * gtPenalty;

  // Step 4: Calculate final fidelity (all adjustments combined)
  const finalFidelity = opts.includeEfficiency
    ? gtAdjustedFidelity * efficiency.averageEfficiency
    : gtAdjustedFidelity;

  // Add efficiency info to junctions
  const junctionsWithEfficiency = baseFidelity.junctions.map(j => ({
    ...j,
    efficiency: calculateEfficiency(j.overhang),
  }));

  return {
    ...baseFidelity,
    junctions: junctionsWithEfficiency,
    sortedJunctions: [...junctionsWithEfficiency].sort((a, b) => a.fidelity - b.fidelity),
    efficiency,
    efficiencyAdjustedFidelity,
    efficiencyAdjustedFidelityPercent: `${(efficiencyAdjustedFidelity * 100).toFixed(1)}%`,
    gtRisks,
    hasGTRisks: gtRisks.length > 0,
    gtAdjustedFidelity,
    gtAdjustedFidelityPercent: `${(gtAdjustedFidelity * 100).toFixed(1)}%`,
    finalFidelity,
    finalFidelityPercent: `${(finalFidelity * 100).toFixed(1)}%`,
  };
}

// ============================================================================
// CONVENIENCE EXPORTS
// ============================================================================

/** Re-export efficiency functions for convenience */
export {
  calculateEfficiency,
  calculateSetEfficiency,
  isPalindrome,
  isHomopolymer,
  type EfficiencyResult,
  type SetEfficiencyResult,
};

/**
 * Quick fidelity check (simplified interface)
 */
export function quickFidelityCheck(
  overhangs: string[],
  enzyme: string = 'BsaI'
): { fidelity: number; isPerfect: boolean; warnings: string[] } {
  const result = calculateFidelity(overhangs, enzyme);
  return {
    fidelity: result.finalFidelity,
    isPerfect: result.finalFidelity >= 0.9999,
    warnings: result.warnings,
  };
}

/**
 * Compare fidelity across multiple enzymes
 */
export function compareEnzymeFidelity(
  overhangs: string[]
): { byEnzyme: Record<string, FidelityResult>; recommended: string } {
  const results: Record<string, FidelityResult> = {};

  for (const enzyme of ENZYMES_WITH_DATA) {
    results[enzyme] = calculateFidelity(overhangs, enzyme);
  }

  // Find best enzyme
  let best = ENZYMES_WITH_DATA[0];
  let bestFidelity = 0;

  for (const [enzyme, result] of Object.entries(results)) {
    if (result.finalFidelity > bestFidelity) {
      bestFidelity = result.finalFidelity;
      best = enzyme;
    }
  }

  return {
    byEnzyme: results,
    recommended: best,
  };
}
