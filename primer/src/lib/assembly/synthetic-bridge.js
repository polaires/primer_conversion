/**
 * Synthetic Bridge Generation Module
 *
 * Automatically designs synthetic oligonucleotides (bridges) to solve
 * difficult assembly junctions:
 * - Poor overlap Tm
 * - High secondary structure
 * - Repeat sequences
 * - GC-extreme regions
 *
 * Supports ordering as:
 * - Ultramer (up to 200 bp, IDT)
 * - gBlock (up to 3000 bp, IDT)
 * - Eblock (up to 1500 bp, IDT)
 * - Gene fragment (Twist, Genscript, etc.)
 */

import { tm } from '../tm.js';
import { calculateTmQ5, calculateGC } from '../tmQ5.js';
import { dg } from '../dg.js';
import { reverseComplement } from '../sequenceUtils.js';

/**
 * Synthesis product specifications
 */
const SYNTHESIS_PRODUCTS = {
  ultramer: {
    name: 'IDT Ultramer',
    minLength: 60,
    maxLength: 200,
    baseCost: 0.30, // per bp
    turnaround: '1-2 days',
    description: 'High-fidelity synthetic oligo',
  },
  eblock: {
    name: 'IDT Eblock',
    minLength: 200,
    maxLength: 1500,
    baseCost: 0.10, // per bp (volume pricing)
    turnaround: '3-5 days',
    description: 'Double-stranded DNA block',
  },
  gblock: {
    name: 'IDT gBlock',
    minLength: 125,
    maxLength: 3000,
    baseCost: 0.12, // per bp
    turnaround: '3-5 days',
    description: 'Double-stranded gene fragment',
  },
  twist: {
    name: 'Twist Gene Fragment',
    minLength: 300,
    maxLength: 5000,
    baseCost: 0.07, // per bp
    turnaround: '5-7 days',
    description: 'High-throughput gene synthesis',
  },
};

/**
 * Bridge design parameters
 */
const BRIDGE_PARAMS = {
  // Optimal overlap lengths
  overlapLength: {
    min: 20,
    optimal: 30,
    max: 50,
  },

  // Target Tm range
  overlapTm: {
    min: 55,
    optimal: 62,
    max: 68,
  },

  // Acceptable secondary structure
  maxHairpinDG: -3.0, // kcal/mol

  // GC content targets
  gcContent: {
    min: 0.40,
    max: 0.60,
  },
};

/**
 * Find optimal overlap region within a sequence
 * Returns the best position and length for overlap
 */
function findOptimalOverlap(sequence, targetTm = BRIDGE_PARAMS.overlapTm.optimal) {
  const results = [];

  // Scan through possible overlap positions and lengths
  for (let length = BRIDGE_PARAMS.overlapLength.min; length <= BRIDGE_PARAMS.overlapLength.max; length++) {
    for (let start = 0; start <= sequence.length - length; start++) {
      const overlap = sequence.substring(start, start + length);
      const overlapTm = calculateTmQ5(overlap);
      const gc = calculateGC(overlap);
      const hairpinDG = dg(overlap);

      // Score based on how close to targets
      const tmScore = 1 - Math.abs(overlapTm - targetTm) / 20;
      const gcScore = gc >= BRIDGE_PARAMS.gcContent.min && gc <= BRIDGE_PARAMS.gcContent.max ? 1 : 0.5;
      const hairpinScore = hairpinDG > BRIDGE_PARAMS.maxHairpinDG ? 1 : 0.5;

      const totalScore = tmScore * gcScore * hairpinScore;

      results.push({
        start,
        length,
        sequence: overlap,
        tm: overlapTm,
        gc,
        hairpinDG,
        score: totalScore,
      });
    }
  }

  // Sort by score and return best
  results.sort((a, b) => b.score - a.score);
  return results[0];
}

/**
 * Design a synthetic bridge for a problematic junction
 *
 * @param {Object} junction - Junction information
 * @param {string} junction.leftSequence - Sequence ending at junction (3' end)
 * @param {string} junction.rightSequence - Sequence starting at junction (5' end)
 * @param {Object} options - Design options
 * @returns {Object} Bridge design
 */
export function designSyntheticBridge(junction, options = {}) {
  const {
    leftSequence,
    rightSequence,
    leftContext = '',  // Additional context before junction
    rightContext = '', // Additional context after junction
  } = junction;

  const {
    targetTm = BRIDGE_PARAMS.overlapTm.optimal,
    minOverlap = BRIDGE_PARAMS.overlapLength.min,
    maxBridgeLength = 200, // Default to Ultramer-compatible
    includeRestrictionSites = [], // e.g., ['BsaI', 'BsmBI']
  } = options;

  // 1. Find optimal overlap on left side
  const leftRegion = leftContext + leftSequence;
  const leftSearchRegion = leftRegion.slice(-100); // Search last 100 bp
  const leftOverlap = findOptimalOverlap(leftSearchRegion, targetTm);

  // 2. Find optimal overlap on right side
  const rightRegion = rightSequence + rightContext;
  const rightSearchRegion = rightRegion.slice(0, 100); // Search first 100 bp
  const rightOverlap = findOptimalOverlap(rightSearchRegion, targetTm);

  // 3. Calculate bridge sequence
  // Bridge = left overlap + junction region + right overlap
  const leftOverlapSeq = leftOverlap.sequence;
  const rightOverlapSeq = rightOverlap.sequence;

  // The junction region is what connects them
  // For a minimal bridge, we just need the overlaps
  // For a longer bridge, we include more context

  let bridgeSequence;
  let bridgeType;

  // Determine minimum bridge length
  const minBridgeLength = leftOverlap.length + rightOverlap.length;

  if (minBridgeLength <= SYNTHESIS_PRODUCTS.ultramer.maxLength) {
    // Can use Ultramer
    bridgeSequence = leftOverlapSeq + rightOverlapSeq;
    bridgeType = 'ultramer';
  } else if (minBridgeLength <= SYNTHESIS_PRODUCTS.gblock.maxLength) {
    // Need gBlock
    bridgeSequence = leftOverlapSeq + rightOverlapSeq;
    bridgeType = 'gblock';
  } else {
    // Very long bridge - use gene fragment
    bridgeSequence = leftOverlapSeq + rightOverlapSeq;
    bridgeType = 'twist';
  }

  // 4. Validate bridge properties
  const bridgeGC = calculateGC(bridgeSequence);
  const bridgeHairpinDG = dg(bridgeSequence);

  // Check for problematic features
  const warnings = [];

  if (bridgeGC < 0.25 || bridgeGC > 0.75) {
    warnings.push({
      type: 'extreme_gc',
      value: bridgeGC,
      message: `Bridge has extreme GC content (${(bridgeGC * 100).toFixed(1)}%) - may be difficult to synthesize`,
    });
  }

  if (bridgeHairpinDG < -8) {
    warnings.push({
      type: 'strong_secondary_structure',
      value: bridgeHairpinDG,
      message: `Bridge may form strong secondary structure (ΔG = ${bridgeHairpinDG.toFixed(1)} kcal/mol)`,
    });
  }

  // Check for homopolymer runs
  const homopolymerMatch = bridgeSequence.match(/(.)\1{5,}/g);
  if (homopolymerMatch) {
    warnings.push({
      type: 'homopolymer_run',
      value: homopolymerMatch[0],
      message: `Bridge contains long homopolymer run (${homopolymerMatch[0]}) - may cause synthesis issues`,
    });
  }

  // 5. Calculate cost estimate
  const product = SYNTHESIS_PRODUCTS[bridgeType];
  const estimatedCost = bridgeSequence.length * product.baseCost;

  return {
    // Bridge sequence
    bridgeSequence,
    bridgeLength: bridgeSequence.length,
    bridgeGC: (bridgeGC * 100).toFixed(1) + '%',
    bridgeHairpinDG: bridgeHairpinDG.toFixed(1) + ' kcal/mol',

    // Overlap details
    leftOverlap: {
      sequence: leftOverlapSeq,
      length: leftOverlap.length,
      tm: leftOverlap.tm.toFixed(1) + '°C',
      gc: (leftOverlap.gc * 100).toFixed(1) + '%',
      position: 'end of left fragment',
    },
    rightOverlap: {
      sequence: rightOverlapSeq,
      length: rightOverlap.length,
      tm: rightOverlap.tm.toFixed(1) + '°C',
      gc: (rightOverlap.gc * 100).toFixed(1) + '%',
      position: 'start of right fragment',
    },

    // Ordering information
    orderAs: product.name,
    productType: bridgeType,
    estimatedCost: '$' + estimatedCost.toFixed(2),
    turnaround: product.turnaround,

    // Warnings and suggestions
    warnings,
    synthesisComplexity: warnings.length === 0 ? 'low' :
                         warnings.length <= 2 ? 'medium' : 'high',

    // For ordering
    orderingNotes: [
      `Order as: ${product.name}`,
      `Length: ${bridgeSequence.length} bp`,
      `Estimated cost: $${estimatedCost.toFixed(2)}`,
      warnings.length > 0 ? 'Review synthesis warnings before ordering' : '',
    ].filter(Boolean),
  };
}

/**
 * Design bridges for multiple problematic junctions
 *
 * @param {Object[]} junctions - Array of junction objects
 * @param {Object} options - Design options
 * @returns {Object[]} Array of bridge designs
 */
export function designMultipleBridges(junctions, options = {}) {
  return junctions.map((junction, index) => ({
    junctionIndex: index,
    ...designSyntheticBridge(junction, options),
  }));
}

/**
 * Analyze if a junction needs a synthetic bridge
 *
 * @param {Object} junction - Junction information
 * @returns {Object} Analysis result
 */
export function analyzeJunctionForBridge(junction) {
  const {
    overlapTm,
    overlapGC,
    hairpinDG,
    overlapLength,
    hasRepeats = false,
  } = junction;

  const issues = [];
  let needsBridge = false;

  // Check overlap Tm
  if (overlapTm !== undefined && overlapTm < 50) {
    issues.push({
      type: 'low_tm',
      value: overlapTm,
      severity: overlapTm < 45 ? 'high' : 'medium',
      message: `Overlap Tm is low (${overlapTm.toFixed(1)}°C)`,
    });
    needsBridge = true;
  }

  // Check overlap length
  if (overlapLength !== undefined && overlapLength < 18) {
    issues.push({
      type: 'short_overlap',
      value: overlapLength,
      severity: 'high',
      message: `Overlap is too short (${overlapLength} bp)`,
    });
    needsBridge = true;
  }

  // Check secondary structure
  if (hairpinDG !== undefined && hairpinDG < -5) {
    issues.push({
      type: 'secondary_structure',
      value: hairpinDG,
      severity: hairpinDG < -8 ? 'high' : 'medium',
      message: `Strong secondary structure (ΔG = ${hairpinDG.toFixed(1)} kcal/mol)`,
    });
    if (hairpinDG < -8) needsBridge = true;
  }

  // Check GC content
  if (overlapGC !== undefined && (overlapGC < 0.25 || overlapGC > 0.75)) {
    issues.push({
      type: 'extreme_gc',
      value: overlapGC,
      severity: 'medium',
      message: `Extreme GC content (${(overlapGC * 100).toFixed(1)}%)`,
    });
  }

  // Check for repeats
  if (hasRepeats) {
    issues.push({
      type: 'repeat_sequence',
      severity: 'high',
      message: 'Junction contains repeat sequences that may cause mispriming',
    });
    needsBridge = true;
  }

  return {
    needsBridge,
    confidence: issues.length === 0 ? 'high' :
                issues.some(i => i.severity === 'high') ? 'low' : 'medium',
    issues,
    recommendation: needsBridge
      ? 'A synthetic bridge is recommended to ensure reliable assembly'
      : issues.length > 0
      ? 'Junction may work but consider a bridge for higher reliability'
      : 'Junction looks good - no bridge needed',
  };
}

/**
 * Optimize bridge sequence for synthesis
 * Removes problematic features while maintaining function
 */
export function optimizeBridgeForSynthesis(bridgeSequence, options = {}) {
  const {
    targetGC = 0.50,
    maxHomopolymer = 5,
  } = options;

  let optimized = bridgeSequence;
  const changes = [];

  // 1. Break up long homopolymer runs (if possible without affecting overlaps)
  // This is simplified - in practice, we'd be more careful about overlaps
  const homopolymerRegex = /([ATGC])\1{5,}/gi;
  let match;

  while ((match = homopolymerRegex.exec(optimized)) !== null) {
    const run = match[0];
    const base = match[1];
    const position = match.index;

    // Only modify if it's in the middle (not in overlap regions)
    if (position > 30 && position < optimized.length - 30) {
      // Insert a different base to break the run
      const breakBase = base === 'G' ? 'A' : base === 'C' ? 'T' : base === 'A' ? 'G' : 'C';
      const midpoint = Math.floor(run.length / 2);
      const newRun = run.slice(0, midpoint) + breakBase + run.slice(midpoint + 1);

      optimized = optimized.slice(0, position) + newRun + optimized.slice(position + run.length);

      changes.push({
        type: 'break_homopolymer',
        position,
        original: run,
        modified: newRun,
      });
    }
  }

  // 2. Check if optimization changed key properties
  const originalGC = calculateGC(bridgeSequence);
  const optimizedGC = calculateGC(optimized);
  const originalDG = dg(bridgeSequence);
  const optimizedDG = dg(optimized);

  return {
    original: bridgeSequence,
    optimized,
    changes,
    comparison: {
      gcChange: `${(originalGC * 100).toFixed(1)}% → ${(optimizedGC * 100).toFixed(1)}%`,
      dgChange: `${originalDG.toFixed(1)} → ${optimizedDG.toFixed(1)} kcal/mol`,
    },
    improved: changes.length > 0,
  };
}

export { SYNTHESIS_PRODUCTS, BRIDGE_PARAMS };
