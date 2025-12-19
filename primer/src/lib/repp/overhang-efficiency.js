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

/**
 * Count G and C bases in a sequence
 * @param {string} seq - DNA sequence
 * @returns {number} Number of G or C bases
 */
function countGC(seq) {
  return (seq.match(/[GC]/gi) || []).length;
}

/**
 * Count A and T bases in a sequence
 * @param {string} seq - DNA sequence
 * @returns {number} Number of A or T bases
 */
function countAT(seq) {
  return (seq.match(/[AT]/gi) || []).length;
}

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
      severity: 'medium',
    },

    // High GC (100%) - slow melting during cycling
    HIGH_GC: {
      test: (oh) => countGC(oh) === 4,
      factor: 0.85,
      note: '100% GC - may have slow melting kinetics',
      severity: 'low',
    },

    // Low GC (0%) - weak base pairing
    LOW_GC: {
      test: (oh) => countGC(oh) === 0,
      factor: 0.80,
      note: '0% GC - weak Watson-Crick pairing',
      severity: 'medium',
    },

    // Homopolymer runs (3+ of same base)
    HOMOPOLYMER: {
      test: (oh) => /(.)\1{2,}/.test(oh),
      factor: 0.65,
      note: 'Homopolymer run - poor ligation specificity',
      severity: 'high',
    },

    // Near-palindromic (1 base from palindrome)
    NEAR_PALINDROME: {
      test: (oh) => {
        const rc = oh.split('').reverse().map(b => {
          const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
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
      severity: 'medium',
    },
  },

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
  },
};

/**
 * Check if an overhang is palindromic (self-complementary)
 * @param {string} overhang - 4bp overhang sequence
 * @returns {boolean} True if palindromic
 */
export function isPalindrome(overhang) {
  if (!overhang || overhang.length !== 4) return false;
  const oh = overhang.toUpperCase();
  const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
  const rc = oh.split('').reverse().map(b => comp[b] || b).join('');
  return oh === rc;
}

/**
 * Check if an overhang is a homopolymer (all same base)
 * @param {string} overhang - 4bp overhang sequence
 * @returns {boolean} True if homopolymer
 */
export function isHomopolymer(overhang) {
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
 *
 * @param {string} overhang - 4bp overhang sequence
 * @returns {Object} Efficiency factor and warnings
 */
export function calculateEfficiency(overhang) {
  if (!overhang || overhang.length !== 4) {
    return {
      overhang: overhang || '',
      efficiency: 0,
      efficiencyPercent: '0%',
      warnings: ['Invalid overhang length'],
      isOptimal: false,
      isAcceptable: false,
    };
  }

  const oh = overhang.toUpperCase();
  let factor = 1.0;
  const warnings = [];
  const appliedPenalties = [];

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
 *
 * @param {string[]} overhangs - Array of overhang sequences
 * @returns {Object} Set efficiency analysis
 */
export function calculateSetEfficiency(overhangs) {
  if (!Array.isArray(overhangs) || overhangs.length === 0) {
    return {
      overhangs: [],
      individual: [],
      combinedEfficiency: 0,
      averageEfficiency: 0,
      worstOverhang: null,
      warnings: ['No overhangs provided'],
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
    bestOverhang: sorted[sorted.length - 1],
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
 * @param {number} efficiency - Efficiency factor (0-1)
 * @returns {Object} Category info
 */
export function getEfficiencyCategory(efficiency) {
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
