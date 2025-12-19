/**
 * Repp Configuration
 * Default settings for DNA assembly planning
 * Based on: https://github.com/Lattice-Automation/repp
 */

/**
 * Default configuration for assembly planning
 */
export const DEFAULT_CONFIG = {
  // Assembly constraints
  fragmentsMaxCount: 6,           // Maximum fragments in final assembly
  fragmentsMinHomology: 15,       // Minimum junction length between fragments (bp)
  fragmentsMaxHomology: 120,      // Maximum junction length between fragments (bp)
  fragmentsMaxHairpinMelt: 47.0,  // Maximum hairpin Tm at junctions (°C)

  // Gibson Assembly costs
  gibsonAssemblyCost: 12.98,      // Cost per Gibson Assembly reaction ($)
  gibsonAssemblyTimeCost: 0.0,    // Time cost per Gibson Assembly

  // PCR costs
  pcrBpCost: 0.60,                // Cost per bp of primer ($)
  pcrRxnCost: 0.27,               // Cost per PCR reaction ($)
  pcrTimeCost: 0.0,               // Time cost per PCR
  pcrMinLength: 60,               // Minimum PCR fragment length (bp)
  pcrPrimerMaxPairPenalty: 30.0,  // Maximum primer3 pair penalty
  pcrPrimerMaxEmbedLength: 20,    // Maximum embedded bp for junction
  pcrPrimerMaxOfftargetTm: 55.0,  // Maximum off-target Tm (°C)
  pcrBufferLength: 20,            // Buffer for primer selection

  // Synthesis constraints and costs
  syntheticMinLength: 125,        // Minimum synthesis length (bp)
  syntheticMaxLength: 3000,       // Maximum synthesis length (bp)
};

/**
 * Synthetic fragment cost tiers (IDT gBlocks pricing)
 * Key = max length in bp, Value = { fixed: boolean, cost: number }
 */
export const SYNTHETIC_FRAGMENT_COSTS = {
  250: { fixed: true, cost: 89.0 },
  500: { fixed: true, cost: 89.0 },
  750: { fixed: true, cost: 129.0 },
  1000: { fixed: true, cost: 149.0 },
  1250: { fixed: true, cost: 209.0 },
  1500: { fixed: true, cost: 249.0 },
  1750: { fixed: true, cost: 289.0 },
  2000: { fixed: true, cost: 329.0 },
  2250: { fixed: true, cost: 399.0 },
  2500: { fixed: true, cost: 449.0 },
  2750: { fixed: true, cost: 499.0 },
  3000: { fixed: true, cost: 549.0 },
};

/**
 * Synthetic plasmid cost tiers (IDT gene synthesis pricing)
 * For delivery of insert in a plasmid
 */
export const SYNTHETIC_PLASMID_COSTS = {
  500: { fixed: true, cost: 160 },
  3000: { fixed: false, cost: 0.35 },  // per bp
  30000: { fixed: false, cost: 0.60 }, // per bp
};

/**
 * Calculate the cost of synthesizing a fragment
 * @param {number} length - Fragment length in bp
 * @param {Object} config - Configuration object
 * @returns {number} Cost in dollars
 */
export function synthFragmentCost(length, config = DEFAULT_CONFIG) {
  // May need to split into multiple fragments if too long
  const fragCount = Math.ceil(length / config.syntheticMaxLength);
  const fragLength = Math.floor(length / fragCount);

  const costTier = findCostTier(fragLength, SYNTHETIC_FRAGMENT_COSTS);

  if (costTier.fixed) {
    return fragCount * costTier.cost;
  }

  return fragCount * fragLength * costTier.cost;
}

/**
 * Calculate the cost of synthesizing an insert delivered in a plasmid
 * @param {number} insertLength - Insert length in bp
 * @returns {number} Cost in dollars
 */
export function synthPlasmidCost(insertLength) {
  const costTier = findCostTier(insertLength, SYNTHETIC_PLASMID_COSTS);

  if (costTier.fixed) {
    return costTier.cost;
  }

  return insertLength * costTier.cost;
}

/**
 * Find the appropriate cost tier for a given sequence length
 * @param {number} seqLength - Sequence length in bp
 * @param {Object} costs - Cost tier map
 * @returns {Object} Cost tier { fixed, cost }
 */
function findCostTier(seqLength, costs) {
  const lengths = Object.keys(costs).map(Number).sort((a, b) => a - b);

  for (const maxLength of lengths) {
    if (maxLength >= seqLength) {
      return costs[maxLength];
    }
  }

  // Length exceeds all tiers - return very high cost
  return { fixed: true, cost: Number.MAX_SAFE_INTEGER };
}

/**
 * Calculate PCR primer cost
 * @param {string} fwdPrimer - Forward primer sequence
 * @param {string} revPrimer - Reverse primer sequence
 * @param {Object} config - Configuration object
 * @returns {number} Total primer cost
 */
export function primerCost(fwdPrimer, revPrimer, config = DEFAULT_CONFIG) {
  return (fwdPrimer.length + revPrimer.length) * config.pcrBpCost + config.pcrRxnCost;
}

/**
 * Merge user config with defaults
 * @param {Object} userConfig - User-provided configuration
 * @returns {Object} Merged configuration
 */
export function mergeConfig(userConfig = {}) {
  return { ...DEFAULT_CONFIG, ...userConfig };
}
