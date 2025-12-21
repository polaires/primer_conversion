/**
 * Overhang Efficiency Penalties for Golden Gate Assembly
 *
 * Implements SOFT penalties for efficiency-reducing overhang patterns.
 * Based on NEB research (Potapov et al. 2018, Pryor et al. 2020).
 *
 * Key principle: TNNA and other patterns reduce efficiency but DON'T
 * eliminate ligation. Use soft penalties (multipliers) not hard filters.
 *
 * References:
 * - Potapov et al. 2018 ACS Synth Biol - Overhang fidelity profiling
 * - NEB Technical Guide: Golden Gate Assembly
 */

// ============================================================================
// TYPES
// ============================================================================

interface PatternPenalty {
  regex?: RegExp;
  test?: (overhang: string) => boolean;
  factor: number;
  note: string;
  severity: 'low' | 'medium' | 'high';
}

export interface EfficiencyPenalty {
  type: 'pattern' | 'specific';
  name: string;
  factor: number;
  originalFactor?: number;
  severity?: 'low' | 'medium' | 'high';
  reducedDueToSpecific?: boolean;
}

export interface EfficiencyResult {
  overhang: string;
  efficiency: number;
  efficiencyPercent: string;
  warnings: string[];
  appliedPenalties: EfficiencyPenalty[];
  isOptimal: boolean;
  isAcceptable: boolean;
  isPoor?: boolean;
  gcCount: number;
  gcPercent: number;
  isPalindrome: boolean;
  isHomopolymer: boolean;
  isTNNA: boolean;
}

export interface SetEfficiencyResult {
  overhangs: string[];
  individual: EfficiencyResult[];
  combinedEfficiency: number;
  combinedEfficiencyPercent: string;
  averageEfficiency: number;
  averageEfficiencyPercent: string;
  worstOverhang: EfficiencyResult;
  bestOverhang: EfficiencyResult;
  warnings: string[];
  summary: {
    total: number;
    optimal: number;
    acceptable: number;
    poor: number;
    palindromes: number;
    homopolymers: number;
    tnnaPattern: number;
  };
  recommendation: string;
}

interface EfficiencyCategory {
  category: 'excellent' | 'good' | 'acceptable' | 'poor' | 'avoid';
  color: string;
  description: string;
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Count G and C bases in a sequence
 */
function countGC(seq: string): number {
  return (seq.match(/[GC]/gi) || []).length;
}

/**
 * Count A and T bases in a sequence
 */
function countAT(seq: string): number {
  return (seq.match(/[AT]/gi) || []).length;
}

// ============================================================================
// CONFIGURATION
// ============================================================================

/**
 * Ligation efficiency penalties based on NEB research
 * These are SOFT penalties - reduce score, don't eliminate
 *
 * Pattern-based penalties are multiplicative factors (0-1)
 * where 1.0 = no penalty, lower = reduced efficiency
 */
export const EFFICIENCY_PENALTIES = {
  // Pattern-based penalties (multiplicative)
  patterns: {
    // TNNA overhangs have ~30% reduced ligation efficiency
    // But they CAN work - just slower
    TNNA: {
      regex: /^T..A$/,
      factor: 0.70,
      note: 'TNNA pattern - reduced ligation efficiency (~30%)',
      severity: 'medium' as const,
    },

    // High GC (100%) - slow melting during cycling
    HIGH_GC: {
      test: (oh: string) => countGC(oh) === 4,
      factor: 0.85,
      note: '100% GC - may have slow melting kinetics',
      severity: 'low' as const,
    },

    // Low GC (0%) - weak base pairing
    LOW_GC: {
      test: (oh: string) => countGC(oh) === 0,
      factor: 0.80,
      note: '0% GC - weak Watson-Crick pairing',
      severity: 'medium' as const,
    },

    // Homopolymer runs (3+ of same base)
    HOMOPOLYMER: {
      test: (oh: string) => /(.)\1{2,}/.test(oh),
      factor: 0.65,
      note: 'Homopolymer run - poor ligation specificity',
      severity: 'high' as const,
    },

    // Near-palindromic (1 base from palindrome)
    NEAR_PALINDROME: {
      test: (oh: string) => {
        const rc = oh.split('').reverse().map(b => {
          const comp: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
          return comp[b.toUpperCase()] || b;
        }).join('');
        let mismatches = 0;
        for (let i = 0; i < oh.length; i++) {
          if (oh[i].toUpperCase() !== rc[i]) mismatches++;
        }
        return mismatches === 1; // Near-palindrome: 1 mismatch from perfect
      },
      factor: 0.75,
      note: 'Near-palindromic - may have self-ligation tendency',
      severity: 'medium' as const,
    },
  } as Record<string, PatternPenalty>,

  // Specific overhang penalties (from experimental data)
  // These are known problematic sequences from Pryor et al.
  specific: {
    'TAAA': 0.65,
    'TTTA': 0.65,
    'AAAA': 0.55,  // Homopolymer + TNNA-like
    'TTTT': 0.55,  // Homopolymer
    'CCCC': 0.50,  // Homopolymer - very poor
    'GGGG': 0.45,  // Homopolymer - synthesis issues too
    'ATAT': 0.50,  // Palindrome
    'TATA': 0.50,  // Palindrome
    'GCGC': 0.45,  // Palindrome
    'CGCG': 0.45,  // Palindrome
    'ACGT': 0.40,  // Palindrome
    'CATG': 0.35,  // Strong palindrome (NlaIII site)
    'GATC': 0.30,  // Strong palindrome (Sau3AI/DpnI site)
  } as Record<string, number>,
};

// ============================================================================
// MAIN FUNCTIONS
// ============================================================================

/**
 * Check if an overhang is palindromic (self-complementary)
 */
export function isPalindrome(overhang: string): boolean {
  if (!overhang || overhang.length !== 4) return false;
  const oh = overhang.toUpperCase();
  const comp: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
  const rc = oh.split('').reverse().map(b => comp[b] || b).join('');
  return oh === rc;
}

/**
 * Check if an overhang is a homopolymer (all same base)
 */
export function isHomopolymer(overhang: string): boolean {
  if (!overhang || overhang.length === 0) return false;
  const oh = overhang.toUpperCase();
  return oh.split('').every(b => b === oh[0]);
}

/**
 * Calculate effective ligation efficiency for an overhang
 *
 * This is the core function that applies all efficiency penalties.
 * Returns a factor between 0 and 1 that should be multiplied with
 * the base fidelity score.
 */
export function calculateEfficiency(overhang: string): EfficiencyResult {
  if (!overhang || overhang.length !== 4) {
    return {
      overhang: overhang || '',
      efficiency: 0,
      efficiencyPercent: '0%',
      warnings: ['Invalid overhang length'],
      appliedPenalties: [],
      isOptimal: false,
      isAcceptable: false,
      gcCount: 0,
      gcPercent: 0,
      isPalindrome: false,
      isHomopolymer: false,
      isTNNA: false,
    };
  }

  const oh = overhang.toUpperCase();
  let factor = 1.0;
  const warnings: string[] = [];
  const appliedPenalties: EfficiencyPenalty[] = [];

  // Check specific overhangs first (these have empirically measured penalties)
  const hasSpecificPenalty = EFFICIENCY_PENALTIES.specific[oh] !== undefined;
  if (hasSpecificPenalty) {
    factor *= EFFICIENCY_PENALTIES.specific[oh];
    warnings.push(`Known low-efficiency overhang: ${oh}`);
    appliedPenalties.push({
      type: 'specific',
      name: oh,
      factor: EFFICIENCY_PENALTIES.specific[oh],
    });
  }

  // Check patterns - these are additive characteristics
  // Apply pattern penalties even if specific penalty exists (they compound)
  for (const [name, penalty] of Object.entries(EFFICIENCY_PENALTIES.patterns)) {
    let matches = false;

    if (penalty.regex) {
      matches = penalty.regex.test(oh);
    } else if (penalty.test) {
      matches = penalty.test(oh);
    }

    if (matches) {
      // If we already have a specific penalty, apply a reduced pattern penalty
      // (since specific penalties already account for pattern effects)
      const patternFactor = hasSpecificPenalty
        ? 1.0 - (1.0 - penalty.factor) * 0.3  // Apply 30% of additional pattern penalty
        : penalty.factor;

      factor *= patternFactor;
      warnings.push(penalty.note);
      appliedPenalties.push({
        type: 'pattern',
        name,
        factor: patternFactor,
        originalFactor: penalty.factor,
        severity: penalty.severity,
        reducedDueToSpecific: hasSpecificPenalty,
      });
    }
  }

  // Calculate GC content info
  const gcCount = countGC(oh);
  const gcPercent = (gcCount / 4) * 100;

  return {
    overhang: oh,
    efficiency: factor,
    efficiencyPercent: `${(factor * 100).toFixed(0)}%`,
    warnings,
    appliedPenalties,
    isOptimal: factor >= 0.90,
    isAcceptable: factor >= 0.70,
    isPoor: factor < 0.70,
    gcCount,
    gcPercent,
    isPalindrome: isPalindrome(oh),
    isHomopolymer: isHomopolymer(oh),
    isTNNA: /^T..A$/i.test(oh),
  };
}

/**
 * Calculate efficiency for a set of overhangs
 * Returns individual and combined efficiency metrics
 */
export function calculateSetEfficiency(overhangs: string[]): SetEfficiencyResult {
  if (!Array.isArray(overhangs) || overhangs.length === 0) {
    return {
      overhangs: [],
      individual: [],
      combinedEfficiency: 0,
      combinedEfficiencyPercent: '0%',
      averageEfficiency: 0,
      averageEfficiencyPercent: '0%',
      worstOverhang: {} as EfficiencyResult,
      bestOverhang: {} as EfficiencyResult,
      warnings: ['No overhangs provided'],
      summary: {
        total: 0,
        optimal: 0,
        acceptable: 0,
        poor: 0,
        palindromes: 0,
        homopolymers: 0,
        tnnaPattern: 0,
      },
      recommendation: 'No overhangs to analyze',
    };
  }

  const individual = overhangs.map(oh => calculateEfficiency(oh));

  // Combined efficiency is product of all individual efficiencies
  const combinedEfficiency = individual.reduce((prod, e) => prod * e.efficiency, 1.0);

  // Average efficiency
  const averageEfficiency = individual.reduce((sum, e) => sum + e.efficiency, 0) / individual.length;

  // Find worst overhang
  const sorted = [...individual].sort((a, b) => a.efficiency - b.efficiency);
  const worstOverhang = sorted[0];
  const bestOverhang = sorted[sorted.length - 1];

  // Collect all warnings
  const allWarnings = individual.flatMap(e => e.warnings);

  // Count problematic overhangs
  const problematicCount = individual.filter(e => e.isPoor).length;
  const suboptimalCount = individual.filter(e => !e.isOptimal && !e.isPoor).length;

  return {
    overhangs: overhangs.map(oh => oh.toUpperCase()),
    individual,
    combinedEfficiency,
    combinedEfficiencyPercent: `${(combinedEfficiency * 100).toFixed(1)}%`,
    averageEfficiency,
    averageEfficiencyPercent: `${(averageEfficiency * 100).toFixed(1)}%`,
    worstOverhang,
    bestOverhang,
    warnings: allWarnings,
    summary: {
      total: individual.length,
      optimal: individual.filter(e => e.isOptimal).length,
      acceptable: individual.filter(e => e.isAcceptable && !e.isOptimal).length,
      poor: problematicCount,
      palindromes: individual.filter(e => e.isPalindrome).length,
      homopolymers: individual.filter(e => e.isHomopolymer).length,
      tnnaPattern: individual.filter(e => e.isTNNA).length,
    },
    recommendation: problematicCount > 0
      ? 'Consider replacing poor-efficiency overhangs'
      : suboptimalCount > 2
        ? 'Some overhangs have reduced efficiency - may need more colonies'
        : 'Efficiency looks good',
  };
}

/**
 * Get efficiency category for display
 */
export function getEfficiencyCategory(efficiency: number): EfficiencyCategory {
  if (efficiency >= 0.95) {
    return { category: 'excellent', color: 'green', description: 'Optimal efficiency' };
  }
  if (efficiency >= 0.85) {
    return { category: 'good', color: 'blue', description: 'Good efficiency' };
  }
  if (efficiency >= 0.70) {
    return { category: 'acceptable', color: 'yellow', description: 'Acceptable but reduced efficiency' };
  }
  if (efficiency >= 0.50) {
    return { category: 'poor', color: 'orange', description: 'Poor efficiency - consider alternatives' };
  }
  return { category: 'avoid', color: 'red', description: 'Very poor - strongly recommend alternatives' };
}

export { countGC, countAT };
