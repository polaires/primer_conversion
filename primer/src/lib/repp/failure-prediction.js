/**
 * Failure Prediction for Golden Gate Assembly
 *
 * Predicts potential failure modes and suggests mitigations.
 * Based on empirical data from assembly troubleshooting literature.
 *
 * Key failure modes:
 * 1. Cross-ligation (wrong fragments ligating)
 * 2. Self-ligation (same fragment circularizing)
 * 3. Missing fragments (incomplete assembly)
 * 4. Multiple bands (mixture of products)
 * 5. No product (complete failure)
 */

import { reverseComplement } from './enzymes.js';
import { getEnzymeLigationData, calculateExperimentalFidelity } from './goldengate.js';
import { findGTMismatchRisks, calculateEnhancedFidelity } from './goldengate-primer-optimizer.js';
import { isPalindrome, isHomopolymer, calculateSetEfficiency } from './overhang-efficiency.js';

/**
 * Failure mode definitions
 */
export const FAILURE_MODES = {
  CROSS_LIGATION: {
    name: 'Cross-ligation',
    description: 'Incorrect fragments ligating due to overhang similarity',
    severity: 'high',
    symptoms: ['Multiple bands', 'Wrong size product', 'Incorrect sequence'],
    mitigation: 'Choose overhangs with higher specificity, increase temperature',
  },
  SELF_LIGATION: {
    name: 'Self-ligation',
    description: 'Fragment circularizing on itself',
    severity: 'medium',
    symptoms: ['Reduced yield', 'Extra bands', 'Shorter products'],
    mitigation: 'Avoid palindromic overhangs, use higher fragment concentrations',
  },
  INCOMPLETE_ASSEMBLY: {
    name: 'Incomplete assembly',
    description: 'Missing one or more fragments',
    severity: 'medium',
    symptoms: ['Ladder pattern', 'Shorter products', 'Low yield'],
    mitigation: 'Balance fragment concentrations, extend ligation time',
  },
  LOW_EFFICIENCY: {
    name: 'Low efficiency ligation',
    description: 'Overhangs with poor ligation kinetics',
    severity: 'medium',
    symptoms: ['Low yield', 'Longer reaction needed', 'Background'],
    mitigation: 'Replace TNNA or low-GC overhangs, increase enzyme',
  },
  G_T_MISMATCH: {
    name: 'G:T wobble mis-ligation',
    description: 'G:T wobble pairing causing mismatch ligation',
    severity: 'high',
    symptoms: ['Incorrect sequences', 'Mixture of correct/incorrect'],
    mitigation: 'Redesign overhangs to avoid G:T pairing potential',
  },
  INTERNAL_SITE: {
    name: 'Internal restriction site',
    description: 'Fragment contains recognition site that gets cut',
    severity: 'critical',
    symptoms: ['Extra bands', 'Small fragments', 'No full-length'],
    mitigation: 'Domesticate sequences or use different enzyme',
  },
  HAIRPIN_FORMATION: {
    name: 'Primer hairpin',
    description: 'Strong secondary structure in primer',
    severity: 'medium',
    symptoms: ['PCR failure', 'Low amplification', 'Artifacts'],
    mitigation: 'Redesign primer with different homology region',
  },
  MISPRIMING: {
    name: 'Mispriming',
    description: 'Primer annealing to wrong location',
    severity: 'medium',
    symptoms: ['Multiple bands', 'Wrong product', 'Smear'],
    mitigation: 'Increase annealing temperature, redesign primer',
  },
};

/**
 * Get cross-ligation matrix for an overhang set
 *
 * @param {string[]} overhangs - Array of overhangs
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Cross-ligation analysis
 */
function getCrossLigationMatrix(overhangs, enzyme) {
  const enzymeData = getEnzymeLigationData(enzyme);
  if (!enzymeData) {
    return { available: false, matrix: null };
  }

  const matrix = [];
  const ohList = overhangs.map(oh => oh.toUpperCase());

  for (let i = 0; i < ohList.length; i++) {
    const row = [];
    const oh1 = ohList[i];

    for (let j = 0; j < ohList.length; j++) {
      const oh2 = ohList[j];
      const rc2 = reverseComplement(oh2);

      // Get ligation frequency
      const freq = enzymeData.matrix?.[oh1]?.[rc2] || 0;
      row.push(freq);
    }
    matrix.push(row);
  }

  return { available: true, matrix, overhangs: ohList };
}

/**
 * Predict cross-ligation risks for an overhang set
 *
 * @param {string[]} overhangs - Array of overhangs
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Cross-ligation risk analysis
 */
export function predictCrossLigation(overhangs, enzyme = 'BsaI') {
  const crossData = getCrossLigationMatrix(overhangs, enzyme);
  const risks = [];

  if (!crossData.available) {
    return {
      enzyme,
      available: false,
      risks: [],
      message: 'No experimental data available for this enzyme',
    };
  }

  const { matrix, overhangs: ohList } = crossData;

  // Find significant off-diagonal values
  for (let i = 0; i < matrix.length; i++) {
    const correctFreq = matrix[i][i]; // Diagonal = correct ligation

    if (correctFreq === 0) continue;

    for (let j = 0; j < matrix[i].length; j++) {
      if (i === j) continue; // Skip diagonal

      const crossFreq = matrix[i][j];
      if (crossFreq > 0) {
        const ratio = crossFreq / correctFreq;

        if (ratio >= 0.01) { // ≥1% cross-ligation
          let severity;
          if (ratio >= 0.10) severity = 'high';
          else if (ratio >= 0.05) severity = 'medium';
          else severity = 'low';

          risks.push({
            source: ohList[i],
            target: ohList[j],
            targetRC: reverseComplement(ohList[j]),
            correctFreq,
            crossFreq,
            ratio,
            ratioPercent: `${(ratio * 100).toFixed(1)}%`,
            severity,
            failureMode: 'CROSS_LIGATION',
          });
        }
      }
    }
  }

  // Sort by ratio (worst first)
  risks.sort((a, b) => b.ratio - a.ratio);

  // Calculate overall cross-ligation risk
  const highRisks = risks.filter(r => r.severity === 'high');
  const mediumRisks = risks.filter(r => r.severity === 'medium');

  let overallRisk;
  if (highRisks.length > 0) overallRisk = 'high';
  else if (mediumRisks.length > 0) overallRisk = 'medium';
  else if (risks.length > 0) overallRisk = 'low';
  else overallRisk = 'minimal';

  return {
    enzyme,
    available: true,
    overhangs: ohList,
    risks,
    riskCount: risks.length,
    highRiskCount: highRisks.length,
    mediumRiskCount: mediumRisks.length,
    overallRisk,
    worstCase: risks[0] || null,
    matrix,
  };
}

/**
 * Predict self-ligation risks
 *
 * @param {string[]} overhangs - Array of overhangs
 * @returns {Object} Self-ligation analysis
 */
export function predictSelfLigation(overhangs) {
  const risks = [];

  for (const oh of overhangs) {
    const ohUpper = oh.toUpperCase();

    // Check if palindromic
    if (isPalindrome(ohUpper)) {
      risks.push({
        overhang: ohUpper,
        reason: 'Palindromic overhang can self-ligate',
        severity: 'high',
        failureMode: 'SELF_LIGATION',
      });
    }

    // Check if near-palindrome (1 mismatch away)
    const rc = reverseComplement(ohUpper);
    let mismatches = 0;
    for (let i = 0; i < ohUpper.length; i++) {
      if (ohUpper[i] !== rc[i]) mismatches++;
    }
    if (mismatches === 1) {
      risks.push({
        overhang: ohUpper,
        reason: 'Near-palindromic overhang may have self-ligation tendency',
        severity: 'medium',
        failureMode: 'SELF_LIGATION',
      });
    }
  }

  return {
    overhangs: overhangs.map(oh => oh.toUpperCase()),
    risks,
    hasRisks: risks.length > 0,
    palindromeCount: risks.filter(r => r.reason.includes('Palindromic')).length,
    nearPalindromeCount: risks.filter(r => r.reason.includes('Near-palindromic')).length,
  };
}

/**
 * Predict efficiency issues
 *
 * @param {string[]} overhangs - Array of overhangs
 * @returns {Object} Efficiency analysis
 */
export function predictEfficiencyIssues(overhangs) {
  const setEfficiency = calculateSetEfficiency(overhangs);
  const risks = [];

  for (const individual of setEfficiency.individual) {
    if (individual.efficiency < 0.70) {
      risks.push({
        overhang: individual.overhang,
        efficiency: individual.efficiency,
        efficiencyPercent: individual.efficiencyPercent,
        reason: individual.warnings.join('; ') || 'Low efficiency overhang',
        severity: individual.efficiency < 0.50 ? 'high' : 'medium',
        failureMode: 'LOW_EFFICIENCY',
      });
    }
  }

  return {
    overhangs: overhangs.map(oh => oh.toUpperCase()),
    risks,
    hasRisks: risks.length > 0,
    combinedEfficiency: setEfficiency.combinedEfficiency,
    averageEfficiency: setEfficiency.averageEfficiency,
    worstOverhang: setEfficiency.worstOverhang,
  };
}

/**
 * Comprehensive failure prediction for an assembly
 *
 * @param {string[]} overhangs - Array of overhangs
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Additional options
 * @returns {Object} Comprehensive failure prediction
 */
export function predictFailureModes(overhangs, enzyme = 'BsaI', options = {}) {
  const {
    primerData = [],          // Array of primer quality data
    internalSites = [],       // Internal sites found
  } = options;

  const predictions = [];
  let overallRisk = 'minimal';

  // 1. Cross-ligation analysis
  const crossLigation = predictCrossLigation(overhangs, enzyme);
  if (crossLigation.risks.length > 0) {
    predictions.push({
      mode: FAILURE_MODES.CROSS_LIGATION,
      details: crossLigation,
      severity: crossLigation.overallRisk,
    });
    if (crossLigation.overallRisk === 'high') overallRisk = 'high';
    else if (crossLigation.overallRisk === 'medium' && overallRisk !== 'high') {
      overallRisk = 'medium';
    }
  }

  // 2. Self-ligation analysis
  const selfLigation = predictSelfLigation(overhangs);
  if (selfLigation.hasRisks) {
    predictions.push({
      mode: FAILURE_MODES.SELF_LIGATION,
      details: selfLigation,
      severity: selfLigation.risks.some(r => r.severity === 'high') ? 'high' : 'medium',
    });
    if (selfLigation.risks.some(r => r.severity === 'high') && overallRisk !== 'high') {
      overallRisk = overallRisk === 'minimal' ? 'medium' : overallRisk;
    }
  }

  // 3. G:T wobble analysis
  const gtRisks = findGTMismatchRisks(overhangs);
  if (gtRisks.length > 0) {
    predictions.push({
      mode: FAILURE_MODES.G_T_MISMATCH,
      details: {
        risks: gtRisks,
        riskCount: gtRisks.length,
        criticalCount: gtRisks.filter(r => r.risk === 'critical').length,
      },
      severity: gtRisks.some(r => r.risk === 'critical') ? 'high' : 'medium',
    });
    if (gtRisks.some(r => r.risk === 'critical')) overallRisk = 'high';
  }

  // 4. Efficiency analysis
  const efficiencyIssues = predictEfficiencyIssues(overhangs);
  if (efficiencyIssues.hasRisks) {
    predictions.push({
      mode: FAILURE_MODES.LOW_EFFICIENCY,
      details: efficiencyIssues,
      severity: efficiencyIssues.risks.some(r => r.severity === 'high') ? 'high' : 'medium',
    });
  }

  // 5. Internal sites
  if (internalSites.length > 0) {
    predictions.push({
      mode: FAILURE_MODES.INTERNAL_SITE,
      details: {
        sites: internalSites,
        count: internalSites.length,
      },
      severity: 'critical',
    });
    overallRisk = 'critical';
  }

  // Calculate overall assembly fidelity
  const fidelityData = calculateEnhancedFidelity(overhangs, enzyme);

  // Sort predictions by severity
  const severityOrder = { critical: 0, high: 1, medium: 2, low: 3 };
  predictions.sort((a, b) => severityOrder[a.severity] - severityOrder[b.severity]);

  // Generate recommendations
  const recommendations = generateRecommendations(predictions);

  return {
    enzyme,
    overhangs: overhangs.map(oh => oh.toUpperCase()),
    numOverhangs: overhangs.length,
    overallRisk,
    predictions,
    predictedIssuesCount: predictions.length,
    fidelity: fidelityData,
    recommendations,
    riskSummary: {
      critical: predictions.filter(p => p.severity === 'critical').length,
      high: predictions.filter(p => p.severity === 'high').length,
      medium: predictions.filter(p => p.severity === 'medium').length,
      low: predictions.filter(p => p.severity === 'low').length,
    },
    expectedSuccessRate: calculateExpectedSuccessRate(predictions, fidelityData),
  };
}

/**
 * Generate recommendations based on predictions
 *
 * @param {Array} predictions - Array of failure predictions
 * @returns {Array} Prioritized recommendations
 */
function generateRecommendations(predictions) {
  const recommendations = [];

  for (const pred of predictions) {
    if (pred.severity === 'critical' || pred.severity === 'high') {
      recommendations.push({
        priority: pred.severity === 'critical' ? 1 : 2,
        issue: pred.mode.name,
        recommendation: pred.mode.mitigation,
        action: getSpecificAction(pred),
      });
    }
  }

  // Sort by priority
  recommendations.sort((a, b) => a.priority - b.priority);

  return recommendations;
}

/**
 * Get specific action for a prediction
 *
 * @param {Object} prediction - Failure prediction
 * @returns {string} Specific action to take
 */
function getSpecificAction(prediction) {
  const details = prediction.details;

  switch (prediction.mode.name) {
    case 'Cross-ligation':
      if (details.worstCase) {
        return `Replace overhang ${details.worstCase.source} or ${details.worstCase.target}`;
      }
      return 'Run findOptimalOverhangSet() to find alternative';

    case 'Self-ligation':
      const palindromes = details.risks.filter(r => r.reason.includes('Palindromic'));
      if (palindromes.length > 0) {
        return `Replace palindromic overhang(s): ${palindromes.map(p => p.overhang).join(', ')}`;
      }
      return 'Review near-palindromic overhangs';

    case 'G:T wobble mis-ligation':
      return 'Redesign overhangs to eliminate G:T wobble pairing potential';

    case 'Internal restriction site':
      return 'Run suggestDomestication() or use different enzyme';

    case 'Low efficiency ligation':
      if (details.worstOverhang) {
        return `Replace low-efficiency overhang: ${details.worstOverhang.overhang}`;
      }
      return 'Replace TNNA or low-GC overhangs';

    default:
      return prediction.mode.mitigation;
  }
}

/**
 * Calculate expected success rate based on predictions
 *
 * Uses an additive penalty model with diminishing returns to avoid
 * overly pessimistic estimates when multiple issues are present.
 *
 * @param {Array} predictions - Failure predictions
 * @param {Object} fidelityData - Fidelity analysis
 * @returns {Object} Expected success rate
 */
function calculateExpectedSuccessRate(predictions, fidelityData) {
  let baseFidelity = fidelityData.gtAdjustedFidelity || fidelityData.baseFidelity || 0.8;

  // Calculate total penalty using additive model with diminishing returns
  // Each severity level contributes a penalty, but penalties don't compound multiplicatively
  let totalPenalty = 0;
  let criticalCount = 0;
  let highCount = 0;
  let mediumCount = 0;
  let lowCount = 0;

  for (const pred of predictions) {
    switch (pred.severity) {
      case 'critical':
        criticalCount++;
        break;
      case 'high':
        highCount++;
        break;
      case 'medium':
        mediumCount++;
        break;
      case 'low':
        lowCount++;
        break;
    }
  }

  // Apply penalties with diminishing returns for multiple issues of same type
  // Critical: First one is severe, additional ones add less
  totalPenalty += criticalCount > 0 ? 0.70 + Math.min(criticalCount - 1, 3) * 0.10 : 0;
  // High: Significant but not fatal
  totalPenalty += highCount > 0 ? 0.25 + Math.min(highCount - 1, 3) * 0.10 : 0;
  // Medium: Moderate concern
  totalPenalty += mediumCount > 0 ? 0.10 + Math.min(mediumCount - 1, 4) * 0.03 : 0;
  // Low: Minor concern
  totalPenalty += lowCount > 0 ? 0.03 + Math.min(lowCount - 1, 4) * 0.01 : 0;

  // Cap total penalty at 95% (always some chance of success if fidelity is non-zero)
  totalPenalty = Math.min(0.95, totalPenalty);

  // Apply penalty to base fidelity
  const adjustedFidelity = baseFidelity * (1 - totalPenalty);

  // Clamp to reasonable range
  const finalFidelity = Math.max(0.01, Math.min(1.0, adjustedFidelity));

  let tier;
  if (finalFidelity >= 0.90) tier = 'excellent';
  else if (finalFidelity >= 0.70) tier = 'good';
  else if (finalFidelity >= 0.50) tier = 'moderate';
  else if (finalFidelity >= 0.20) tier = 'low';
  else tier = 'very_low';

  return {
    rate: finalFidelity,
    ratePercent: `${(finalFidelity * 100).toFixed(0)}%`,
    tier,
    coloniesNeeded: Math.ceil(3 / finalFidelity), // To get 3 correct colonies on average
    penaltyBreakdown: {
      critical: criticalCount,
      high: highCount,
      medium: mediumCount,
      low: lowCount,
      totalPenalty: `${(totalPenalty * 100).toFixed(0)}%`,
    },
  };
}

/**
 * Quick risk assessment for an overhang set
 *
 * @param {string[]} overhangs - Array of overhangs
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Quick risk summary
 */
export function quickRiskAssessment(overhangs, enzyme = 'BsaI') {
  const fidelity = calculateEnhancedFidelity(overhangs, enzyme);
  const efficiency = calculateSetEfficiency(overhangs);
  const selfLigation = predictSelfLigation(overhangs);

  const issues = [];

  if (fidelity.baseFidelity < 0.90) {
    issues.push(`Fidelity ${fidelity.baseFidelityPercent} < 90%`);
  }
  if (efficiency.combinedEfficiency < 0.70) {
    issues.push(`Efficiency ${efficiency.combinedEfficiencyPercent} < 70%`);
  }
  if (selfLigation.palindromeCount > 0) {
    issues.push(`${selfLigation.palindromeCount} palindromic overhang(s)`);
  }
  if (fidelity.hasGTRisks) {
    issues.push(`${fidelity.gtRisks.length} G:T mismatch risk(s)`);
  }

  return {
    overhangs: overhangs.map(oh => oh.toUpperCase()),
    fidelity: fidelity.baseFidelity,
    efficiency: efficiency.combinedEfficiency,
    issues,
    isGood: issues.length === 0,
    riskLevel: issues.length === 0 ? 'low' :
               issues.length <= 2 ? 'medium' : 'high',
  };
}
