/**
 * Cross-Reactivity Analysis for Golden Gate Assembly
 *
 * Provides tools for visualizing and analyzing overhang cross-reactivity:
 * - Heatmap data generation
 * - Fidelity matrix calculation
 * - Problem junction identification
 * - Compatibility scoring
 */

import ligationData from './ligation-data.json';
import { reverseComplement } from './enzymes.js';

/**
 * Get ligation matrix for an enzyme
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Ligation matrix
 */
function getMatrix(enzyme) {
  const keyMap = {
    'BsaI': 'BsaI-HFv2',
    'BsaI-HFv2': 'BsaI-HFv2',
    'BbsI': 'BbsI-HF',
    'BbsI-HF': 'BbsI-HF',
    'BsmBI': 'BsmBI-v2',
    'BsmBI-v2': 'BsmBI-v2',
    'Esp3I': 'Esp3I',
    'SapI': 'SapI',
  };

  const dataKey = keyMap[enzyme];
  if (!dataKey || !ligationData.enzymes[dataKey]) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  return ligationData.enzymes[dataKey].matrix;
}

/**
 * Generate cross-reactivity heatmap data for a set of overhangs
 *
 * Returns a matrix where each cell [i][j] represents the ligation
 * frequency between overhang i and the reverse complement of overhang j.
 *
 * Diagonal cells (i=j) show correct ligation frequency.
 * Off-diagonal cells show cross-ligation (bad) frequencies.
 *
 * @param {string[]} overhangs - Array of overhangs to analyze
 * @param {string} enzyme - Enzyme name (default: BsaI)
 * @returns {Object} Heatmap data
 *
 * @example
 * const data = generateHeatmapData(['GGAG', 'AATG', 'GCTT'], 'BsaI');
 * // Use data.matrix for visualization
 */
export function generateHeatmapData(overhangs, enzyme = 'BsaI') {
  const matrix = getMatrix(enzyme);
  const n = overhangs.length;

  // Normalize overhangs
  const normalizedOH = overhangs.map(oh => oh.toUpperCase());

  // Build heatmap matrix
  // Row i: overhang i trying to ligate
  // Col j: reverse complement of overhang j (potential partner)
  const heatmap = [];
  let maxValue = 0;
  let minValue = Infinity;

  for (let i = 0; i < n; i++) {
    const row = [];
    const ohI = normalizedOH[i];

    for (let j = 0; j < n; j++) {
      const ohJ = normalizedOH[j];
      const rcJ = reverseComplement(ohJ);

      // Ligation frequency of overhang i with reverse complement of j
      const freq = matrix[ohI]?.[rcJ] || 0;

      row.push(freq);

      if (freq > 0) {
        maxValue = Math.max(maxValue, freq);
        minValue = Math.min(minValue, freq);
      }
    }

    heatmap.push(row);
  }

  // Calculate normalized values (0-1) for visualization
  const normalizedHeatmap = heatmap.map(row =>
    row.map(val => val > 0 ? (val - minValue) / (maxValue - minValue || 1) : 0)
  );

  // Identify cross-ligation issues (off-diagonal with significant frequency)
  const crossLigationIssues = [];

  for (let i = 0; i < n; i++) {
    const ohI = normalizedOH[i];
    const correctFreq = heatmap[i][i];

    for (let j = 0; j < n; j++) {
      if (i === j) continue;

      const crossFreq = heatmap[i][j];
      if (crossFreq > 0 && correctFreq > 0) {
        const ratio = crossFreq / correctFreq;

        if (ratio > 0.01) { // >1% cross-reactivity
          crossLigationIssues.push({
            from: ohI,
            to: reverseComplement(normalizedOH[j]),
            crossFreq,
            correctFreq,
            ratio,
            severity: ratio > 0.1 ? 'high' : ratio > 0.05 ? 'medium' : 'low',
          });
        }
      }
    }
  }

  // Sort by severity
  crossLigationIssues.sort((a, b) => b.ratio - a.ratio);

  return {
    enzyme,
    overhangs: normalizedOH,
    reverseComplements: normalizedOH.map(oh => reverseComplement(oh)),
    matrix: heatmap,
    normalizedMatrix: normalizedHeatmap,
    statistics: {
      maxFrequency: maxValue,
      minFrequency: minValue === Infinity ? 0 : minValue,
    },
    crossLigationIssues,
    hasCriticalIssues: crossLigationIssues.some(issue => issue.severity === 'high'),
  };
}

/**
 * Calculate per-junction fidelity for a set of overhangs
 *
 * @param {string[]} overhangs - Array of overhangs
 * @param {string} enzyme - Enzyme name
 * @returns {Object[]} Array of junction fidelity data
 */
export function calculateJunctionFidelities(overhangs, enzyme = 'BsaI') {
  const matrix = getMatrix(enzyme);
  const normalizedOH = overhangs.map(oh => oh.toUpperCase());

  const junctions = normalizedOH.map((oh, index) => {
    const rc = reverseComplement(oh);

    // Correct ligation (both directions)
    const correct = (matrix[oh]?.[rc] || 0) + (matrix[rc]?.[oh] || 0);

    // Total possible in this set
    let total = 0;
    for (const other of normalizedOH) {
      const otherRc = reverseComplement(other);
      total += matrix[oh]?.[other] || 0;
      total += matrix[oh]?.[otherRc] || 0;
      total += matrix[rc]?.[other] || 0;
      total += matrix[rc]?.[otherRc] || 0;
    }

    const fidelity = total > 0 ? correct / total : 0;

    // Find worst cross-reaction partner
    let worstCross = { partner: null, frequency: 0 };
    for (const other of normalizedOH) {
      if (other === oh) continue;
      const otherRc = reverseComplement(other);

      const cross1 = matrix[oh]?.[otherRc] || 0;
      const cross2 = matrix[rc]?.[other] || 0;
      const crossTotal = cross1 + cross2;

      if (crossTotal > worstCross.frequency) {
        worstCross = { partner: other, frequency: crossTotal };
      }
    }

    return {
      index,
      overhang: oh,
      reverseComplement: rc,
      fidelity,
      fidelityPercent: (fidelity * 100).toFixed(1) + '%',
      correctFrequency: correct,
      totalFrequency: total,
      worstCrossReaction: worstCross,
      status: fidelity >= 0.99 ? 'excellent' :
              fidelity >= 0.95 ? 'good' :
              fidelity >= 0.90 ? 'acceptable' :
              fidelity >= 0.80 ? 'marginal' : 'poor',
    };
  });

  // Sort by fidelity to identify weakest links
  const sortedByFidelity = [...junctions].sort((a, b) => a.fidelity - b.fidelity);

  // Calculate overall assembly fidelity
  const overallFidelity = junctions.reduce((prod, j) => prod * j.fidelity, 1);

  return {
    enzyme,
    junctions,
    weakestJunction: sortedByFidelity[0],
    strongestJunction: sortedByFidelity[sortedByFidelity.length - 1],
    overallFidelity,
    overallFidelityPercent: (overallFidelity * 100).toFixed(2) + '%',
    expectedCorrectAssemblies: overallFidelity,
    coloniesToScreen: Math.ceil(10 / overallFidelity),
  };
}

/**
 * Find compatible replacement overhangs for problematic ones
 *
 * @param {string[]} currentOverhangs - Current overhang set
 * @param {number} problemIndex - Index of overhang to replace
 * @param {string} enzyme - Enzyme name
 * @param {number} topN - Number of alternatives to return
 * @returns {Object[]} Sorted alternatives
 */
export function findBetterAlternatives(currentOverhangs, problemIndex, enzyme = 'BsaI', topN = 5) {
  const matrix = getMatrix(enzyme);
  const enzymeData = ligationData.enzymes[enzyme.includes('-') ? enzyme : `${enzyme}-HFv2`] || ligationData.enzymes[enzyme];

  if (!enzymeData) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  const allOverhangs = enzymeData.overhangs;
  const normalizedCurrent = currentOverhangs.map(oh => oh.toUpperCase());
  const problemOH = normalizedCurrent[problemIndex];

  // Get overhangs not in current set
  const used = new Set(normalizedCurrent);
  normalizedCurrent.forEach(oh => used.add(reverseComplement(oh)));

  const candidates = allOverhangs.filter(oh => {
    const rc = reverseComplement(oh);
    // Not already used, not palindromic
    return !used.has(oh) && !used.has(rc) && oh !== rc;
  });

  // Score each candidate
  const scored = candidates.map(candidate => {
    // Create test set with replacement
    const testSet = [...normalizedCurrent];
    testSet[problemIndex] = candidate;

    // Calculate fidelity for this candidate specifically
    const rc = reverseComplement(candidate);
    const correct = (matrix[candidate]?.[rc] || 0) + (matrix[rc]?.[candidate] || 0);

    let total = 0;
    for (const other of testSet) {
      const otherRc = reverseComplement(other);
      total += matrix[candidate]?.[other] || 0;
      total += matrix[candidate]?.[otherRc] || 0;
      total += matrix[rc]?.[other] || 0;
      total += matrix[rc]?.[otherRc] || 0;
    }

    const fidelity = total > 0 ? correct / total : 0;

    // Also calculate overall assembly fidelity with this replacement
    let overallFidelity = 1.0;
    for (let i = 0; i < testSet.length; i++) {
      const oh = testSet[i];
      const ohRc = reverseComplement(oh);
      const ohCorrect = (matrix[oh]?.[ohRc] || 0) + (matrix[ohRc]?.[oh] || 0);

      let ohTotal = 0;
      for (const other of testSet) {
        const otherRc = reverseComplement(other);
        ohTotal += matrix[oh]?.[other] || 0;
        ohTotal += matrix[oh]?.[otherRc] || 0;
        ohTotal += matrix[ohRc]?.[other] || 0;
        ohTotal += matrix[ohRc]?.[otherRc] || 0;
      }

      if (ohTotal > 0) {
        overallFidelity *= ohCorrect / ohTotal;
      }
    }

    return {
      overhang: candidate,
      reverseComplement: rc,
      junctionFidelity: fidelity,
      junctionFidelityPercent: (fidelity * 100).toFixed(1) + '%',
      overallAssemblyFidelity: overallFidelity,
      overallFidelityPercent: (overallFidelity * 100).toFixed(2) + '%',
      improvement: fidelity - (enzymeData.overhangFidelity?.[problemOH] || 0),
    };
  });

  // Sort by overall assembly fidelity
  scored.sort((a, b) => b.overallAssemblyFidelity - a.overallAssemblyFidelity);

  return scored.slice(0, topN);
}

/**
 * Generate SVG heatmap visualization data
 *
 * @param {string[]} overhangs - Array of overhangs
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Visualization options
 * @returns {Object} SVG-ready data
 */
export function generateHeatmapSVGData(overhangs, enzyme = 'BsaI', options = {}) {
  const {
    cellSize = 40,
    padding = 60,
    colorScale = 'viridis', // 'viridis', 'plasma', 'inferno', 'magma'
  } = options;

  const data = generateHeatmapData(overhangs, enzyme);
  const n = overhangs.length;

  // Color scales (simplified - in real app use d3-scale)
  const colorScales = {
    viridis: (t) => {
      // Simplified viridis approximation
      const r = Math.round(255 * Math.min(1, 0.267 + t * 0.733));
      const g = Math.round(255 * Math.min(1, t));
      const b = Math.round(255 * (1 - t * 0.5));
      return `rgb(${r},${g},${b})`;
    },
    redGreen: (t) => {
      // Red (cross-ligation) to Green (correct)
      if (t < 0.5) {
        return `rgb(${Math.round(255 * (1 - t * 2))}, 0, 0)`;
      } else {
        return `rgb(0, ${Math.round(255 * ((t - 0.5) * 2))}, 0)`;
      }
    },
  };

  const getColor = colorScales[colorScale] || colorScales.viridis;

  // Generate cell data
  const cells = [];

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const value = data.normalizedMatrix[i][j];
      const rawValue = data.matrix[i][j];
      const isDiagonal = i === j;

      cells.push({
        row: i,
        col: j,
        x: padding + j * cellSize,
        y: padding + i * cellSize,
        width: cellSize - 2,
        height: cellSize - 2,
        value: rawValue,
        normalizedValue: value,
        color: getColor(value),
        isDiagonal,
        tooltip: isDiagonal
          ? `${data.overhangs[i]} ↔ ${data.reverseComplements[i]}: ${rawValue} (correct)`
          : `${data.overhangs[i]} → ${data.reverseComplements[j]}: ${rawValue} (cross-ligation)`,
      });
    }
  }

  // Generate axis labels
  const rowLabels = data.overhangs.map((oh, i) => ({
    text: oh,
    x: padding - 5,
    y: padding + i * cellSize + cellSize / 2,
    anchor: 'end',
  }));

  const colLabels = data.reverseComplements.map((rc, j) => ({
    text: rc,
    x: padding + j * cellSize + cellSize / 2,
    y: padding - 5,
    anchor: 'middle',
    transform: `rotate(-45, ${padding + j * cellSize + cellSize / 2}, ${padding - 5})`,
  }));

  return {
    width: padding * 2 + n * cellSize,
    height: padding * 2 + n * cellSize,
    cells,
    rowLabels,
    colLabels,
    title: `Cross-Reactivity Matrix (${enzyme})`,
    legend: {
      min: data.statistics.minFrequency,
      max: data.statistics.maxFrequency,
      label: 'Ligation Frequency',
    },
    issues: data.crossLigationIssues,
  };
}

/**
 * Export utilities
 */
export {
  getMatrix,
};
