/**
 * Assembly Yield Prediction Module
 *
 * Predicts assembly success based on:
 * - Junction fidelity (from ligation data)
 * - Fragment properties (size, GC content)
 * - Assembly method (Golden Gate vs Gibson)
 * - Known risk factors
 *
 * References:
 * - Pryor et al. 2020 PLOS ONE (ligation fidelity data)
 * - Potapov et al. 2018 ACS Synth Biol (optimal overhang sets)
 * - Gibson et al. 2009 Nature Methods (assembly efficiency)
 */

import { reverseComplement } from '../repp/enzymes.js';
import { calculateGC } from '../tmQ5.js';
import ligationData from '../repp/ligation-data.json';

/**
 * Risk factor thresholds based on literature and empirical data
 */
const RISK_THRESHOLDS = {
  // Junction fidelity thresholds
  junctionFidelity: {
    excellent: 0.99,
    good: 0.95,
    acceptable: 0.90,
    marginal: 0.80,
    poor: 0.70,
  },

  // Fragment size thresholds (bp)
  fragmentSize: {
    optimal: { min: 200, max: 3000 },
    acceptable: { min: 100, max: 5000 },
    difficult: { min: 50, max: 8000 },
  },

  // GC content thresholds
  gcContent: {
    optimal: { min: 0.40, max: 0.60 },
    acceptable: { min: 0.30, max: 0.70 },
  },

  // Overlap Tm for Gibson/HiFi (°C)
  overlapTm: {
    optimal: { min: 50, max: 65 },
    acceptable: { min: 45, max: 70 },
  },

  // Total assembly size (bp)
  totalSize: {
    easy: 5000,
    moderate: 10000,
    difficult: 15000,
    veryDifficult: 20000,
  },
};

/**
 * Analyze fragment properties for risk factors
 */
function analyzeFragment(fragment, index) {
  const risks = [];
  const size = fragment.sequence?.length || fragment.size || 0;
  const gc = fragment.gc || calculateGC(fragment.sequence);

  // Size analysis
  if (size < RISK_THRESHOLDS.fragmentSize.acceptable.min) {
    risks.push({
      type: 'fragment_too_short',
      severity: 'high',
      fragment: index,
      value: size,
      threshold: RISK_THRESHOLDS.fragmentSize.acceptable.min,
      message: `Fragment ${index + 1} is very short (${size} bp) - may have low transformation efficiency`,
    });
  } else if (size > RISK_THRESHOLDS.fragmentSize.acceptable.max) {
    risks.push({
      type: 'fragment_too_long',
      severity: size > RISK_THRESHOLDS.fragmentSize.difficult.max ? 'high' : 'medium',
      fragment: index,
      value: size,
      threshold: RISK_THRESHOLDS.fragmentSize.acceptable.max,
      message: `Fragment ${index + 1} is large (${size} bp) - consider splitting`,
    });
  }

  // GC content analysis
  if (gc < RISK_THRESHOLDS.gcContent.acceptable.min) {
    risks.push({
      type: 'low_gc_content',
      severity: gc < 0.25 ? 'high' : 'medium',
      fragment: index,
      value: gc,
      message: `Fragment ${index + 1} has low GC content (${(gc * 100).toFixed(1)}%) - may affect amplification`,
    });
  } else if (gc > RISK_THRESHOLDS.gcContent.acceptable.max) {
    risks.push({
      type: 'high_gc_content',
      severity: gc > 0.75 ? 'high' : 'medium',
      fragment: index,
      value: gc,
      message: `Fragment ${index + 1} has high GC content (${(gc * 100).toFixed(1)}%) - may form secondary structures`,
    });
  }

  return risks;
}

/**
 * Calculate junction fidelity from ligation data
 */
function calculateJunctionFidelity(overhangs, enzyme = 'BsaI') {
  const enzymeKey = enzyme.includes('-') ? enzyme : `${enzyme}-HFv2`;
  const enzymeData = ligationData.enzymes[enzymeKey] || ligationData.enzymes['BsaI-HFv2'];
  const matrix = enzymeData.matrix;

  const junctions = [];
  let overallFidelity = 1.0;

  for (let i = 0; i < overhangs.length; i++) {
    const oh = overhangs[i].toUpperCase();
    const rc = reverseComplement(oh);

    // Correct ligation frequency
    const correct = (matrix[oh]?.[rc] || 0) + (matrix[rc]?.[oh] || 0);

    // Total possible ligations in this set
    let total = 0;
    for (const other of overhangs) {
      const otherUpper = other.toUpperCase();
      const otherRc = reverseComplement(otherUpper);
      total += matrix[oh]?.[otherUpper] || 0;
      total += matrix[oh]?.[otherRc] || 0;
      total += matrix[rc]?.[otherUpper] || 0;
      total += matrix[rc]?.[otherRc] || 0;
    }

    const fidelity = total > 0 ? correct / total : 0;
    overallFidelity *= fidelity;

    // Find worst cross-reactivity
    let worstCross = null;
    let worstCrossFreq = 0;

    for (const other of overhangs) {
      if (other.toUpperCase() === oh) continue;
      const otherRc = reverseComplement(other.toUpperCase());
      const crossFreq = (matrix[oh]?.[otherRc] || 0) + (matrix[rc]?.[other.toUpperCase()] || 0);

      if (crossFreq > worstCrossFreq) {
        worstCrossFreq = crossFreq;
        worstCross = other;
      }
    }

    junctions.push({
      index: i,
      overhang: oh,
      reverseComplement: rc,
      fidelity,
      fidelityPercent: (fidelity * 100).toFixed(1) + '%',
      status: fidelity >= RISK_THRESHOLDS.junctionFidelity.excellent ? 'excellent' :
              fidelity >= RISK_THRESHOLDS.junctionFidelity.good ? 'good' :
              fidelity >= RISK_THRESHOLDS.junctionFidelity.acceptable ? 'acceptable' :
              fidelity >= RISK_THRESHOLDS.junctionFidelity.marginal ? 'marginal' : 'poor',
      worstCrossReactivity: worstCross ? {
        partner: worstCross,
        frequency: worstCrossFreq,
      } : null,
    });
  }

  return { junctions, overallFidelity };
}

/**
 * Predict assembly yield for Golden Gate assembly
 *
 * @param {Object} assembly - Assembly configuration
 * @param {string[]} assembly.overhangs - Array of overhang sequences
 * @param {Object[]} assembly.fragments - Array of fragment objects with sequence/size
 * @param {string} [assembly.enzyme='BsaI'] - Type IIS enzyme
 * @param {string} [assembly.vector] - Vector sequence (optional)
 * @returns {Object} Yield prediction
 */
export function predictGoldenGateYield(assembly) {
  const {
    overhangs = [],
    fragments = [],
    enzyme = 'BsaI',
    vector = null,
  } = assembly;

  const risks = [];

  // 1. Calculate junction fidelities
  const { junctions, overallFidelity } = calculateJunctionFidelity(overhangs, enzyme);

  // Identify low-fidelity junctions
  for (const junction of junctions) {
    if (junction.fidelity < RISK_THRESHOLDS.junctionFidelity.acceptable) {
      risks.push({
        type: 'low_fidelity_junction',
        severity: junction.fidelity < RISK_THRESHOLDS.junctionFidelity.poor ? 'high' : 'medium',
        junction: junction.index,
        overhang: junction.overhang,
        fidelity: junction.fidelity,
        message: `Junction ${junction.index + 1} (${junction.overhang}) has ${junction.status} fidelity (${junction.fidelityPercent})`,
        suggestion: junction.worstCrossReactivity
          ? `Cross-reacts with ${junction.worstCrossReactivity.partner}`
          : 'Consider using a different overhang',
      });
    }
  }

  // 2. Analyze fragments
  for (let i = 0; i < fragments.length; i++) {
    const fragmentRisks = analyzeFragment(fragments[i], i);
    risks.push(...fragmentRisks);
  }

  // 3. Calculate total assembly size
  const totalSize = fragments.reduce((sum, f) => sum + (f.sequence?.length || f.size || 0), 0);

  if (totalSize > RISK_THRESHOLDS.totalSize.difficult) {
    risks.push({
      type: 'large_total_size',
      severity: totalSize > RISK_THRESHOLDS.totalSize.veryDifficult ? 'high' : 'medium',
      value: totalSize,
      message: `Total assembly size (${(totalSize / 1000).toFixed(1)} kb) may reduce transformation efficiency`,
    });
  }

  // 4. Check number of fragments
  const numFragments = fragments.length;
  if (numFragments > 8) {
    risks.push({
      type: 'many_fragments',
      severity: numFragments > 12 ? 'high' : 'medium',
      value: numFragments,
      message: `${numFragments} fragments increases complexity - fidelity drops exponentially`,
    });
  }

  // 5. Calculate expected outcomes
  // Transformation efficiency decreases with size (empirical)
  const sizeEfficiencyFactor = Math.exp(-totalSize / 50000);

  // Colony count estimate (assuming ~10^6 CFU/μg baseline)
  const baseColonies = 1000;
  const expectedColonies = Math.round(baseColonies * overallFidelity * sizeEfficiencyFactor);

  // Colonies to screen for 95% confidence of getting correct clone
  // P(at least one correct in n colonies) = 1 - (1-p)^n >= 0.95
  // n >= log(0.05) / log(1-p)
  const coloniesToScreen = overallFidelity > 0
    ? Math.ceil(Math.log(0.05) / Math.log(1 - overallFidelity))
    : 100;

  // 6. Determine confidence level
  const highSeverityRisks = risks.filter(r => r.severity === 'high').length;
  const mediumSeverityRisks = risks.filter(r => r.severity === 'medium').length;

  let confidence;
  if (highSeverityRisks > 0) {
    confidence = 'low';
  } else if (mediumSeverityRisks > 2) {
    confidence = 'medium';
  } else if (overallFidelity >= 0.80) {
    confidence = 'high';
  } else {
    confidence = 'medium';
  }

  return {
    method: 'Golden Gate',
    enzyme,

    // Core predictions
    overallFidelity,
    overallFidelityPercent: (overallFidelity * 100).toFixed(1) + '%',
    expectedCorrectRate: (overallFidelity * 100).toFixed(0) + '%',
    expectedColonies,
    coloniesToScreen: Math.min(coloniesToScreen, 100),

    // Junction details
    junctions,
    weakestJunction: junctions.reduce((min, j) => j.fidelity < min.fidelity ? j : min, junctions[0]),

    // Assembly stats
    numFragments,
    totalSize,
    totalSizeKb: (totalSize / 1000).toFixed(1) + ' kb',

    // Risk analysis
    risks,
    riskSummary: {
      high: highSeverityRisks,
      medium: mediumSeverityRisks,
      low: risks.filter(r => r.severity === 'low').length,
    },

    // Overall assessment
    confidence,
    recommendation: confidence === 'high'
      ? 'Assembly looks good. Proceed with standard protocol.'
      : confidence === 'medium'
      ? 'Assembly may need optimization. Consider screening more colonies.'
      : 'Assembly has significant risks. Review and address issues before proceeding.',
  };
}

/**
 * Predict assembly yield for Gibson/HiFi assembly
 *
 * @param {Object} assembly - Assembly configuration
 * @param {Object[]} assembly.junctions - Array of junction objects with overlap info
 * @param {Object[]} assembly.fragments - Array of fragment objects
 * @returns {Object} Yield prediction
 */
export function predictGibsonYield(assembly) {
  const {
    junctions = [],
    fragments = [],
  } = assembly;

  const risks = [];

  // 1. Analyze junction overlaps
  const junctionAnalysis = junctions.map((junction, index) => {
    const { overlapTm, overlapLength, overlapGC, hairpinDG } = junction;

    let status = 'good';
    const junctionRisks = [];

    // Check overlap Tm
    if (overlapTm !== undefined) {
      if (overlapTm < RISK_THRESHOLDS.overlapTm.acceptable.min) {
        status = 'poor';
        junctionRisks.push({
          type: 'low_overlap_tm',
          severity: overlapTm < 40 ? 'high' : 'medium',
          junction: index,
          value: overlapTm,
          message: `Junction ${index + 1} has low overlap Tm (${overlapTm.toFixed(1)}°C)`,
        });
      } else if (overlapTm > RISK_THRESHOLDS.overlapTm.acceptable.max) {
        status = 'marginal';
        junctionRisks.push({
          type: 'high_overlap_tm',
          severity: 'low',
          junction: index,
          value: overlapTm,
          message: `Junction ${index + 1} has high overlap Tm (${overlapTm.toFixed(1)}°C) - may be slow to anneal`,
        });
      }
    }

    // Check overlap length
    if (overlapLength !== undefined) {
      if (overlapLength < 15) {
        status = 'poor';
        junctionRisks.push({
          type: 'short_overlap',
          severity: 'high',
          junction: index,
          value: overlapLength,
          message: `Junction ${index + 1} has short overlap (${overlapLength} bp) - minimum 15 bp recommended`,
        });
      } else if (overlapLength > 40) {
        junctionRisks.push({
          type: 'long_overlap',
          severity: 'low',
          junction: index,
          value: overlapLength,
          message: `Junction ${index + 1} has long overlap (${overlapLength} bp) - may increase primer cost`,
        });
      }
    }

    // Check secondary structure
    if (hairpinDG !== undefined && hairpinDG < -5) {
      status = status === 'good' ? 'marginal' : status;
      junctionRisks.push({
        type: 'overlap_secondary_structure',
        severity: hairpinDG < -8 ? 'high' : 'medium',
        junction: index,
        value: hairpinDG,
        message: `Junction ${index + 1} overlap may form secondary structure (ΔG = ${hairpinDG.toFixed(1)} kcal/mol)`,
      });
    }

    risks.push(...junctionRisks);

    return {
      index,
      ...junction,
      status,
      risks: junctionRisks,
    };
  });

  // 2. Analyze fragments
  for (let i = 0; i < fragments.length; i++) {
    const fragmentRisks = analyzeFragment(fragments[i], i);
    risks.push(...fragmentRisks);
  }

  // 3. Calculate total size
  const totalSize = fragments.reduce((sum, f) => sum + (f.sequence?.length || f.size || 0), 0);

  // 4. Estimate efficiency
  // Gibson assembly efficiency drops with more fragments and larger size
  const numJunctions = junctions.length;
  const junctionEfficiency = Math.pow(0.95, numJunctions); // ~5% loss per junction
  const sizeEfficiency = Math.exp(-totalSize / 30000); // Drops faster than GG
  const goodJunctions = junctionAnalysis.filter(j => j.status === 'good').length;
  const junctionQuality = goodJunctions / Math.max(numJunctions, 1);

  const overallEfficiency = junctionEfficiency * sizeEfficiency * junctionQuality;

  // 5. Expected outcomes
  const baseColonies = 500; // Gibson typically gives fewer colonies
  const expectedColonies = Math.round(baseColonies * overallEfficiency);
  const coloniesToScreen = overallEfficiency > 0
    ? Math.ceil(Math.log(0.05) / Math.log(1 - overallEfficiency))
    : 100;

  // 6. Determine confidence
  const highSeverityRisks = risks.filter(r => r.severity === 'high').length;
  const mediumSeverityRisks = risks.filter(r => r.severity === 'medium').length;

  let confidence;
  if (highSeverityRisks > 0) {
    confidence = 'low';
  } else if (mediumSeverityRisks > 2 || overallEfficiency < 0.50) {
    confidence = 'medium';
  } else {
    confidence = 'high';
  }

  return {
    method: 'Gibson/HiFi Assembly',

    // Core predictions
    overallEfficiency,
    overallEfficiencyPercent: (overallEfficiency * 100).toFixed(1) + '%',
    expectedCorrectRate: (overallEfficiency * 100).toFixed(0) + '%',
    expectedColonies,
    coloniesToScreen: Math.min(coloniesToScreen, 100),

    // Junction details
    junctions: junctionAnalysis,
    problematicJunctions: junctionAnalysis.filter(j => j.status !== 'good'),

    // Assembly stats
    numFragments: fragments.length,
    numJunctions,
    totalSize,
    totalSizeKb: (totalSize / 1000).toFixed(1) + ' kb',

    // Risk analysis
    risks,
    riskSummary: {
      high: highSeverityRisks,
      medium: mediumSeverityRisks,
      low: risks.filter(r => r.severity === 'low').length,
    },

    // Overall assessment
    confidence,
    recommendation: confidence === 'high'
      ? 'Assembly looks good. Use standard NEBuilder HiFi protocol.'
      : confidence === 'medium'
      ? 'Consider optimizing junction overlaps or using synthetic bridges.'
      : 'Assembly has significant risks. Consider redesigning problematic junctions.',
  };
}

/**
 * Generate suggestions to improve assembly yield
 */
export function generateImprovementSuggestions(prediction) {
  const suggestions = [];

  for (const risk of prediction.risks) {
    switch (risk.type) {
      case 'low_fidelity_junction':
        suggestions.push({
          priority: 'high',
          type: 'replace_overhang',
          junction: risk.junction,
          action: `Replace overhang ${risk.overhang} with a higher-fidelity alternative`,
          details: 'Use the cross-reactivity module to find better alternatives',
        });
        break;

      case 'fragment_too_long':
        suggestions.push({
          priority: 'medium',
          type: 'split_fragment',
          fragment: risk.fragment,
          action: `Split fragment ${risk.fragment + 1} into smaller pieces`,
          details: 'Optimal fragment size is 200-3000 bp for best efficiency',
        });
        break;

      case 'low_overlap_tm':
      case 'short_overlap':
        suggestions.push({
          priority: 'high',
          type: 'synthetic_bridge',
          junction: risk.junction,
          action: `Use a synthetic bridge oligo for junction ${risk.junction + 1}`,
          details: 'A gBlock or Ultramer can provide optimal overlap sequences',
        });
        break;

      case 'overlap_secondary_structure':
        suggestions.push({
          priority: 'medium',
          type: 'redesign_overlap',
          junction: risk.junction,
          action: `Shift overlap region to avoid secondary structure at junction ${risk.junction + 1}`,
          details: 'Move junction site 10-20 bp in either direction',
        });
        break;

      case 'many_fragments':
        suggestions.push({
          priority: 'medium',
          type: 'reduce_fragments',
          action: 'Consider combining adjacent fragments or using synthetic fragments',
          details: 'Each additional fragment reduces success rate by ~10-20%',
        });
        break;

      case 'large_total_size':
        suggestions.push({
          priority: 'low',
          type: 'optimize_transformation',
          action: 'Use high-efficiency competent cells and extended recovery',
          details: 'NEB 10-beta or similar high-efficiency cells recommended for large plasmids',
        });
        break;
    }
  }

  // Sort by priority
  const priorityOrder = { high: 0, medium: 1, low: 2 };
  suggestions.sort((a, b) => priorityOrder[a.priority] - priorityOrder[b.priority]);

  return suggestions;
}

export { RISK_THRESHOLDS };
