/**
 * Visualization Utilities for Primer Design
 *
 * Provides analysis and data transformation for visualizing:
 * - GC content hotspots
 * - Secondary structure problems
 * - Primer design conflicts
 */

import { fold, dg, type FoldStructure } from './fold.js';
import { calculateTmQ5, calculateGC } from './tmQ5.js';

// =============================================================================
// Type Definitions
// =============================================================================

/**
 * GC threshold configuration
 */
export interface GCThresholds {
  highGC?: number;
  veryHighGC?: number;
  lowGC?: number;
  veryLowGC?: number;
}

/**
 * Hotspot region
 */
export interface Hotspot {
  start: number;
  end: number;
  type: 'high' | 'low';
  severity: 'warning' | 'severe';
  maxGC: number;
  minGC: number;
}

/**
 * GC hotspot analysis result
 */
export interface GCHotspotAnalysis {
  perPosition: number[];
  hotspots: Hotspot[];
  overallGC: number;
  windowSize: number;
}

/**
 * Structure energy thresholds
 */
export interface StructureThresholds {
  warning?: number;
  severe?: number;
  critical?: number;
}

/**
 * Structure hotspot
 */
export interface StructureHotspot {
  start: number;
  end: number;
  severity: 'warning' | 'severe' | 'critical';
  minDG: number;
}

/**
 * Structure type info
 */
export interface StructureType {
  type: string;
  energy: number;
  pairs: number;
}

/**
 * Structure hotspot analysis result
 */
export interface StructureHotspotAnalysis {
  perPosition: number[];
  hotspots: StructureHotspot[];
  overallDG: number;
  structures: FoldStructure[];
  windowSize: number;
}

/**
 * Base pair info
 */
export interface BasePair {
  i: number;
  j: number;
  type: string;
  energy: number;
  span: number;
}

/**
 * Arc for visualization
 */
export interface Arc {
  start: number;
  end: number;
  height: number;
  type: string;
}

/**
 * Structure visualization data
 */
export interface StructureVisualization {
  sequence: string;
  length: number;
  dotBracket: string[];
  dotBracketString: string;
  pairs: BasePair[];
  pairedPositions: Set<number>;
  unpairedPositions: Set<number>;
  arcs: Arc[];
  structureTypes: StructureType[];
  totalEnergy: number;
}

/**
 * Conflict type
 */
export type ConflictType =
  | 'gc_very_high' | 'gc_high' | 'gc_very_low' | 'gc_low'
  | 'structure_critical' | 'structure_severe' | 'structure_warning'
  | 'repeats' | 'near_start' | 'near_end';

/**
 * Conflict severity
 */
export type ConflictSeverity = 'warning' | 'severe' | 'critical';

/**
 * Design conflict
 */
export interface DesignConflict {
  type: ConflictType;
  severity: ConflictSeverity;
  position: number;
  value: number | RepeatAnalysis;
  message: string;
}

/**
 * Suggestion type
 */
export interface Suggestion {
  type: string;
  message: string;
  positions?: number[];
}

/**
 * Repeat analysis result
 */
export interface RepeatAnalysis {
  hasProblematicRepeats: boolean;
  type: 'homopolymer' | 'dinucleotide' | 'trinucleotide' | null;
  length: number;
  position: number;
}

/**
 * GC data point for visualization
 */
export interface GCDataPoint {
  position: number;
  gc: number;
  severity: 'severe' | 'warning' | 'ok';
}

/**
 * Structure data point for visualization
 */
export interface StructureDataPoint {
  position: number;
  dg: number;
  severity: 'critical' | 'severe' | 'warning' | 'ok';
}

/**
 * Visualization data for conflicts
 */
export interface ConflictVisualization {
  gcData: GCDataPoint[];
  structureData: StructureDataPoint[];
  hotspots: {
    gc: Hotspot[];
    structure: StructureHotspot[];
  };
}

/**
 * Design conflict analysis result
 */
export interface DesignConflictAnalysis {
  position: number;
  regionStart: number;
  regionEnd: number;
  region: string;
  conflicts: DesignConflict[];
  canDesign: boolean;
  suggestions: Suggestion[];
  visualization: ConflictVisualization;
}

/**
 * Design options
 */
export interface DesignOptions {
  windowSize?: number;
  minPrimerLength?: number;
  maxPrimerLength?: number;
}

// =============================================================================
// GC Content Analysis
// =============================================================================

/**
 * Analyze sequence for GC content hotspots
 * Returns regions where GC% is problematically high or low
 *
 * @param sequence - DNA sequence
 * @param windowSize - Sliding window size (default: 20bp)
 * @param thresholds - GC% thresholds
 * @returns Analysis results with hotspot regions
 */
export function analyzeGCHotspots(
  sequence: string,
  windowSize: number = 20,
  thresholds: GCThresholds = {}
): GCHotspotAnalysis {
  const {
    highGC = 0.70,      // >= 70% GC is problematic
    veryHighGC = 0.80,  // >= 80% GC is severe
    lowGC = 0.30,       // <= 30% GC is problematic
    veryLowGC = 0.20,   // <= 20% GC is severe
  } = thresholds;

  const seq = sequence.toUpperCase();
  const n = seq.length;
  const results: GCHotspotAnalysis = {
    perPosition: [],      // GC% at each position (windowed)
    hotspots: [],         // Problematic regions
    overallGC: calculateGC(seq),
    windowSize,
  };

  if (n < windowSize) {
    results.perPosition = new Array(n).fill(results.overallGC);
    return results;
  }

  // Calculate windowed GC% for each position
  let gcCount = 0;

  // Initialize first window
  for (let i = 0; i < windowSize; i++) {
    if (seq[i] === 'G' || seq[i] === 'C') gcCount++;
  }

  // Slide window across sequence
  for (let i = 0; i <= n - windowSize; i++) {
    const gc = gcCount / windowSize;

    // Store GC% centered on window
    const centerPos = i + Math.floor(windowSize / 2);
    results.perPosition[centerPos] = gc;

    // Update window
    if (i + windowSize < n) {
      if (seq[i] === 'G' || seq[i] === 'C') gcCount--;
      if (seq[i + windowSize] === 'G' || seq[i + windowSize] === 'C') gcCount++;
    }
  }

  // Fill edges with nearest value
  const halfWindow = Math.floor(windowSize / 2);
  for (let i = 0; i < halfWindow; i++) {
    results.perPosition[i] = results.perPosition[halfWindow];
  }
  for (let i = n - halfWindow; i < n; i++) {
    results.perPosition[i] = results.perPosition[n - halfWindow - 1];
  }

  // Identify hotspot regions
  let currentHotspot: Hotspot | null = null;

  for (let i = 0; i < n; i++) {
    const gc = results.perPosition[i];
    let severity: 'warning' | 'severe' | null = null;
    let type: 'high' | 'low' | null = null;

    if (gc >= veryHighGC) {
      severity = 'severe';
      type = 'high';
    } else if (gc >= highGC) {
      severity = 'warning';
      type = 'high';
    } else if (gc <= veryLowGC) {
      severity = 'severe';
      type = 'low';
    } else if (gc <= lowGC) {
      severity = 'warning';
      type = 'low';
    }

    if (severity && type) {
      if (currentHotspot && currentHotspot.type === type) {
        // Extend current hotspot
        currentHotspot.end = i;
        if (severity === 'severe') currentHotspot.severity = 'severe';
        currentHotspot.maxGC = Math.max(currentHotspot.maxGC, gc);
        currentHotspot.minGC = Math.min(currentHotspot.minGC, gc);
      } else {
        // Start new hotspot
        if (currentHotspot) results.hotspots.push(currentHotspot);
        currentHotspot = {
          start: i,
          end: i,
          type,
          severity,
          maxGC: gc,
          minGC: gc,
        };
      }
    } else {
      if (currentHotspot) {
        results.hotspots.push(currentHotspot);
        currentHotspot = null;
      }
    }
  }

  if (currentHotspot) results.hotspots.push(currentHotspot);

  return results;
}

// =============================================================================
// Secondary Structure Analysis
// =============================================================================

function severityLevel(severity: 'warning' | 'severe' | 'critical'): number {
  return { warning: 1, severe: 2, critical: 3 }[severity] || 0;
}

/**
 * Analyze sequence for secondary structure hotspots
 * Identifies regions prone to forming hairpins, self-dimers, etc.
 *
 * @param sequence - DNA sequence
 * @param windowSize - Window size for analysis (default: 30bp)
 * @param thresholds - Energy thresholds (kcal/mol)
 * @returns Analysis results with structure hotspots
 */
export function analyzeStructureHotspots(
  sequence: string,
  windowSize: number = 30,
  thresholds: StructureThresholds = {}
): StructureHotspotAnalysis {
  const {
    warning = -3.0,       // ΔG < -3.0 is warning
    severe = -6.0,        // ΔG < -6.0 is severe
    critical = -8.0,      // ΔG < -8.0 is critical
  } = thresholds;

  const seq = sequence.toUpperCase();
  const n = seq.length;
  const results: StructureHotspotAnalysis = {
    perPosition: [],      // ΔG at each position (windowed)
    hotspots: [],         // Problematic regions
    overallDG: n >= 10 ? dg(seq, 37) : 0,
    structures: [],       // Full structure analysis
    windowSize,
  };

  if (n < 10) {
    results.perPosition = new Array(n).fill(0);
    return results;
  }

  // Calculate windowed ΔG for each position
  const halfWindow = Math.floor(windowSize / 2);

  for (let i = 0; i < n; i++) {
    const start = Math.max(0, i - halfWindow);
    const end = Math.min(n, i + halfWindow);
    const subseq = seq.slice(start, end);

    if (subseq.length >= 10) {
      results.perPosition[i] = dg(subseq, 37);
    } else {
      results.perPosition[i] = 0;
    }
  }

  // Get full structure for the sequence
  if (n >= 10) {
    results.structures = fold(seq, 37);
  }

  // Identify hotspot regions
  let currentHotspot: StructureHotspot | null = null;

  for (let i = 0; i < n; i++) {
    const energy = results.perPosition[i];
    let severity: 'warning' | 'severe' | 'critical' | null = null;

    if (energy <= critical) {
      severity = 'critical';
    } else if (energy <= severe) {
      severity = 'severe';
    } else if (energy <= warning) {
      severity = 'warning';
    }

    if (severity) {
      if (currentHotspot) {
        // Extend current hotspot
        currentHotspot.end = i;
        if (severityLevel(severity) > severityLevel(currentHotspot.severity)) {
          currentHotspot.severity = severity;
        }
        currentHotspot.minDG = Math.min(currentHotspot.minDG, energy);
      } else {
        // Start new hotspot
        currentHotspot = {
          start: i,
          end: i,
          severity,
          minDG: energy,
        };
      }
    } else {
      if (currentHotspot) {
        results.hotspots.push(currentHotspot);
        currentHotspot = null;
      }
    }
  }

  if (currentHotspot) results.hotspots.push(currentHotspot);

  return results;
}

// =============================================================================
// Structure Visualization
// =============================================================================

/**
 * Convert fold.js structure output to visualization-ready data
 *
 * @param sequence - DNA sequence
 * @param structures - Output from fold()
 * @returns Visualization data
 */
export function structuresToVisualization(sequence: string, structures: FoldStructure[]): StructureVisualization {
  const n = sequence.length;
  const result: StructureVisualization = {
    sequence,
    length: n,
    dotBracket: new Array(n).fill('.'),  // Dot-bracket notation
    dotBracketString: '',
    pairs: [],                             // [{i, j, type, energy}]
    pairedPositions: new Set(),            // Positions involved in pairing
    unpairedPositions: new Set(),          // Positions not paired
    arcs: [],                              // For arc diagram rendering
    structureTypes: [],                    // Structural motifs
    totalEnergy: 0,
  };

  // Process each structure
  for (const struct of structures) {
    result.totalEnergy += struct.e;

    if (struct.desc) {
      result.structureTypes.push({
        type: struct.desc,
        energy: struct.e,
        pairs: struct.ij.length,
      });
    }

    // Process base pairs
    for (const pair of struct.ij) {
      const i = pair[0];
      const j = pair[1];
      if (i < n && j < n) {
        // Update dot-bracket notation
        result.dotBracket[i] = '(';
        result.dotBracket[j] = ')';

        // Store pair info
        result.pairs.push({
          i,
          j,
          type: struct.desc || 'unknown',
          energy: struct.e / (struct.ij.length || 1),  // Distribute energy
          span: j - i,
        });

        result.pairedPositions.add(i);
        result.pairedPositions.add(j);

        // Create arc for visualization
        result.arcs.push({
          start: i,
          end: j,
          height: (j - i) / 2,  // Arc height proportional to span
          type: struct.desc || 'unknown',
        });
      }
    }
  }

  // Identify unpaired positions
  for (let i = 0; i < n; i++) {
    if (!result.pairedPositions.has(i)) {
      result.unpairedPositions.add(i);
    }
  }

  // Convert dot-bracket to string
  result.dotBracketString = result.dotBracket.join('');

  // Sort arcs by span for proper layering
  result.arcs.sort((a, b) => (b.end - b.start) - (a.end - a.start));

  return result;
}

// =============================================================================
// Design Conflict Analysis
// =============================================================================

/**
 * Analyze sequence for problematic repeats
 */
function analyzeRepeats(sequence: string): RepeatAnalysis {
  const seq = sequence.toUpperCase();
  const result: RepeatAnalysis = {
    hasProblematicRepeats: false,
    type: null,
    length: 0,
    position: 0,
  };

  // Check for homopolymer runs (>= 5 same base)
  const homopolymerMatch = seq.match(/([ATGC])\1{4,}/);
  if (homopolymerMatch) {
    result.hasProblematicRepeats = true;
    result.type = 'homopolymer';
    result.length = homopolymerMatch[0].length;
    result.position = homopolymerMatch.index || 0;
    return result;
  }

  // Check for dinucleotide repeats (>= 4 repeats)
  const dinucMatch = seq.match(/([ATGC]{2})\1{3,}/);
  if (dinucMatch) {
    result.hasProblematicRepeats = true;
    result.type = 'dinucleotide';
    result.length = dinucMatch[0].length;
    result.position = dinucMatch.index || 0;
    return result;
  }

  // Check for trinucleotide repeats (>= 3 repeats)
  const trinucMatch = seq.match(/([ATGC]{3})\1{2,}/);
  if (trinucMatch) {
    result.hasProblematicRepeats = true;
    result.type = 'trinucleotide';
    result.length = trinucMatch[0].length;
    result.position = trinucMatch.index || 0;
    return result;
  }

  return result;
}

/**
 * Find better positions for mutation design
 */
function findBetterPositions(
  sequence: string,
  currentPosition: number,
  gcAnalysis: GCHotspotAnalysis,
  structAnalysis: StructureHotspotAnalysis
): number[] {
  const positions: { pos: number; score: number; gc: number; dg: number }[] = [];
  const searchRadius = 20;

  for (let offset = -searchRadius; offset <= searchRadius; offset++) {
    if (offset === 0) continue;

    const pos = currentPosition + offset;
    if (pos < 15 || pos >= sequence.length - 15) continue;

    // Get metrics at this position (adjusted for analysis window)
    const gcIdx = Math.min(gcAnalysis.perPosition.length - 1, Math.max(0, pos));
    const structIdx = Math.min(structAnalysis.perPosition.length - 1, Math.max(0, pos));

    const gc = gcAnalysis.perPosition[gcIdx] || 0.5;
    const dg = structAnalysis.perPosition[structIdx] || 0;

    // Score this position (lower is better)
    let score = 0;
    if (gc >= 0.70) score += (gc - 0.50) * 100;
    if (gc <= 0.30) score += (0.50 - gc) * 100;
    if (dg <= -5.0) score += Math.abs(dg) * 10;

    if (score < 20) {  // Only include reasonably good positions
      positions.push({ pos, score, gc, dg });
    }
  }

  // Sort by score and return positions
  positions.sort((a, b) => a.score - b.score);
  return positions.map(p => p.pos);
}

/**
 * Analyze why primer design failed at a specific position
 *
 * @param sequence - Full template sequence
 * @param position - Mutation/target position
 * @param options - Design parameters
 * @returns Conflict analysis
 */
export function analyzeDesignConflicts(
  sequence: string,
  position: number,
  options: DesignOptions = {}
): DesignConflictAnalysis {
  const {
    windowSize = 50,
    minPrimerLength = 15,
    maxPrimerLength = 35,
  } = options;

  const seq = sequence.toUpperCase();
  const n = seq.length;

  // Define analysis region around the target position
  const regionStart = Math.max(0, position - windowSize);
  const regionEnd = Math.min(n, position + windowSize);
  const region = seq.slice(regionStart, regionEnd);

  const results: DesignConflictAnalysis = {
    position,
    regionStart,
    regionEnd,
    region,
    conflicts: [],
    canDesign: true,
    suggestions: [],
    visualization: {
      gcData: [],
      structureData: [],
      hotspots: { gc: [], structure: [] },
    },
  };

  // Analyze GC content
  const gcAnalysis = analyzeGCHotspots(region, 15);
  const positionInRegion = position - regionStart;

  // Check GC at mutation position
  const localGC = gcAnalysis.perPosition[positionInRegion];
  if (localGC >= 0.80) {
    results.conflicts.push({
      type: 'gc_very_high',
      severity: 'severe',
      position,
      value: localGC,
      message: `Very high GC content (${(localGC * 100).toFixed(1)}%) - primers will have extremely high Tm`,
    });
    results.canDesign = false;
  } else if (localGC >= 0.70) {
    results.conflicts.push({
      type: 'gc_high',
      severity: 'warning',
      position,
      value: localGC,
      message: `High GC content (${(localGC * 100).toFixed(1)}%) - primers may have high Tm`,
    });
  } else if (localGC <= 0.20) {
    results.conflicts.push({
      type: 'gc_very_low',
      severity: 'severe',
      position,
      value: localGC,
      message: `Very low GC content (${(localGC * 100).toFixed(1)}%) - primers may not bind efficiently`,
    });
  } else if (localGC <= 0.30) {
    results.conflicts.push({
      type: 'gc_low',
      severity: 'warning',
      position,
      value: localGC,
      message: `Low GC content (${(localGC * 100).toFixed(1)}%) - primers may have low Tm`,
    });
  }

  // Analyze secondary structure
  const structAnalysis = analyzeStructureHotspots(region, 25);
  const localDG = structAnalysis.perPosition[positionInRegion];

  if (localDG <= -8.0) {
    results.conflicts.push({
      type: 'structure_critical',
      severity: 'critical',
      position,
      value: localDG,
      message: `Strong secondary structure (ΔG = ${localDG.toFixed(1)} kcal/mol) - primers will form stable hairpins`,
    });
    results.canDesign = false;
  } else if (localDG <= -5.0) {
    results.conflicts.push({
      type: 'structure_severe',
      severity: 'severe',
      position,
      value: localDG,
      message: `Significant secondary structure (ΔG = ${localDG.toFixed(1)} kcal/mol) - primers may form hairpins`,
    });
  } else if (localDG <= -3.0) {
    results.conflicts.push({
      type: 'structure_warning',
      severity: 'warning',
      position,
      value: localDG,
      message: `Moderate secondary structure (ΔG = ${localDG.toFixed(1)} kcal/mol)`,
    });
  }

  // Check for repeats
  const repeatAnalysis = analyzeRepeats(region);
  if (repeatAnalysis.hasProblematicRepeats) {
    results.conflicts.push({
      type: 'repeats',
      severity: 'warning',
      position,
      value: repeatAnalysis,
      message: `Contains ${repeatAnalysis.type} repeat (${repeatAnalysis.length}bp) - may cause mispriming`,
    });
  }

  // Check proximity to edges
  if (position < minPrimerLength) {
    results.conflicts.push({
      type: 'near_start',
      severity: 'severe',
      position,
      value: position,
      message: `Too close to sequence start (${position}bp) - insufficient flanking for primer design`,
    });
    results.canDesign = false;
  }
  if (n - position < minPrimerLength) {
    results.conflicts.push({
      type: 'near_end',
      severity: 'severe',
      position,
      value: n - position,
      message: `Too close to sequence end (${n - position}bp from end) - insufficient flanking for primer design`,
    });
    results.canDesign = false;
  }

  // Generate suggestions
  if (!results.canDesign) {
    // Find better positions nearby
    const betterPositions = findBetterPositions(seq, position, gcAnalysis, structAnalysis);
    if (betterPositions.length > 0) {
      results.suggestions.push({
        type: 'relocate',
        message: `Consider moving mutation to nearby position: ${betterPositions.slice(0, 3).join(', ')}`,
        positions: betterPositions.slice(0, 5),
      });
    }
  }

  // Add visualization data
  results.visualization = {
    gcData: gcAnalysis.perPosition.map((gc, i) => ({
      position: regionStart + i,
      gc,
      severity: gc >= 0.80 ? 'severe' : gc >= 0.70 ? 'warning' : gc <= 0.20 ? 'severe' : gc <= 0.30 ? 'warning' : 'ok',
    })),
    structureData: structAnalysis.perPosition.map((dg, i) => ({
      position: regionStart + i,
      dg,
      severity: dg <= -8.0 ? 'critical' : dg <= -5.0 ? 'severe' : dg <= -3.0 ? 'warning' : 'ok',
    })),
    hotspots: {
      gc: gcAnalysis.hotspots.map(h => ({
        ...h,
        start: h.start + regionStart,
        end: h.end + regionStart,
      })),
      structure: structAnalysis.hotspots.map(h => ({
        ...h,
        start: h.start + regionStart,
        end: h.end + regionStart,
      })),
    },
  };

  return results;
}

// =============================================================================
// Color Utilities
// =============================================================================

/**
 * Generate color for GC content visualization
 */
export function gcToColor(gc: number): string {
  if (gc >= 0.80) return '#ef4444';      // red-500 (severe high)
  if (gc >= 0.70) return '#f97316';      // orange-500 (warning high)
  if (gc >= 0.60) return '#facc15';      // yellow-400 (slightly high)
  if (gc <= 0.20) return '#3b82f6';      // blue-500 (severe low)
  if (gc <= 0.30) return '#60a5fa';      // blue-400 (warning low)
  if (gc <= 0.40) return '#93c5fd';      // blue-300 (slightly low)
  return '#22c55e';                       // green-500 (optimal)
}

/**
 * Generate color for structure energy visualization
 */
export function dgToColor(dg: number): string {
  if (dg <= -8.0) return '#dc2626';      // red-600 (critical)
  if (dg <= -6.0) return '#ef4444';      // red-500 (severe)
  if (dg <= -4.0) return '#f97316';      // orange-500 (warning)
  if (dg <= -2.0) return '#facc15';      // yellow-400 (mild)
  return '#22c55e';                       // green-500 (ok)
}

/**
 * Export all functions
 */
export default {
  analyzeGCHotspots,
  analyzeStructureHotspots,
  structuresToVisualization,
  analyzeDesignConflicts,
  gcToColor,
  dgToColor,
};
