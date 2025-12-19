/**
 * Unified Primer Design Module
 *
 * This module provides a single, unified API for all primer design operations:
 * - Amplification (PCR)
 * - Substitution
 * - Insertion
 * - Deletion
 * - Amino acid mutations (with codon optimization)
 * - Batch operations
 *
 * Core insight: All operations can be expressed as:
 *   REPLACE(template, start, end, replacement)
 *
 * Where:
 * - replacement = null → Amplify (extract region start-end)
 * - replacement = "" → Delete (remove bases start-end)
 * - replacement = "XYZ" → Substitute/Insert
 * - replacement = optimizedCodon → AA mutation (via helper)
 */

import { primers, generateAlternatives } from './primers.js';
import {
  // designSubstitutionPrimers,  // Not implemented yet
  designCodonChangePrimers,
  designInsertionPrimers,
  designDeletionPrimers,
  designRegionSubstitutionPrimers,
  selectOptimalCodon,
  MUTAGENESIS_DEFAULTS,
  MUTATION_TYPES,
} from './mutagenesis.js';

// Stub implementations for functions not yet fully implemented
const designSubstitutionPrimers = (_seq: string, _start: number, _end: number, _replacement: string, _options: any) => ({} as any);
const CODON_TABLE: any = {};
const CODON_TO_AA: any = {};
const AA_NAMES: any = {};
const analyzePrimerPair = (_fwd: any, _rev: any, _options?: any) => ({} as any);
import { reverseComplement, expandAmbiguousBases, hasAmbiguousBases, countCombinations } from './sequenceUtils.js';
import { calculateTmQ5, calculateGC } from './tmQ5.js';
import { calculateEquilibriumEfficiency } from './equilibrium.js';
import { analyzeGQuadruplex } from './scoring.js';

// ============================================================================
// Types and Interfaces
// ============================================================================

export interface AAHelper {
  newAA: string;
  organism?: string;
  orfStart?: number;
  codonPosition?: number;
  originalAA?: string;
}

export interface DesignSpec {
  start: number;
  end: number;
  replacement?: string | null;
  aaHelper?: AAHelper;
  notation?: string;
  notationType?: string;
}

export interface DesignOptions {
  circular?: boolean;
  optimalTm?: number;
  useSmartDesign?: boolean;
  exhaustiveSearch?: boolean;
  maxTmDiff?: number;
  strategy?: string;
  organism?: string;
  orfStart?: number;
  includeEquilibrium?: boolean;
  annealingTemperature?: number;
}

export interface PrimerInfo {
  sequence: string;
  length: number;
  tm: number;
  gc: number;
  gcPercent: string;
  dg: number;
  hasGCClamp: boolean;
  start?: number;
  end?: number;
  gQuadruplex?: any;
}

export interface AmplifiedRegion {
  start: number;
  end: number;
  length: number;
  isWrapped: boolean;
}

export interface DesignResult {
  type: string;
  operation: string;
  design?: string;
  strategy?: string;
  originalSequence: string;
  amplifiedRegion?: AmplifiedRegion;
  description: string;
  forward: PrimerInfo;
  reverse: PrimerInfo;
  annealingTemp?: number;
  tmDifference?: number;
  quality?: string;
  qualityTier?: string;
  compositeScore?: number;
  effectiveScore?: number;
  criticalWarnings?: number;
  piecewiseScores?: any;
  forwardPiecewiseScores?: any;
  reversePiecewiseScores?: any;
  penalty?: number;
  circularWrapped?: boolean;
  alternativePrimers?: any[];
  alternativeCategories?: any;
  hasG4Risk?: boolean;
  warnings?: string[];
  equilibrium?: any;
  isLibrary?: boolean;
  librarySize?: number;
  librarySequence?: string;
}

export interface BatchResult extends DesignResult {
  success: boolean;
  error?: string;
  originalSpec: DesignSpec;
}

export interface ParsedBatch {
  specs: (DesignSpec | { error: string; originalLine: string })[];
  orfStart: number;
}

// ============================================================================
// Circular Sequence Helpers
// ============================================================================

/**
 * Extract a region from a sequence, handling circular wrap-around
 * For circular sequences where start > end, extracts from start to end of sequence
 * plus from beginning of sequence to end position
 *
 * @param seq - The template sequence
 * @param start - Start position (0-based)
 * @param end - End position (0-based, exclusive)
 * @param circular - Whether sequence is circular
 * @returns The extracted region
 */
function extractCircularRegion(seq: string, start: number, end: number, circular: boolean): string {
  if (start <= end) {
    // Normal linear extraction
    return seq.slice(start, end);
  }

  if (!circular) {
    throw new Error(`Invalid region: start=${start} > end=${end} not allowed for linear sequences`);
  }

  // Circular wrap-around: extract from start to end of seq, then from 0 to end
  return seq.slice(start) + seq.slice(0, end);
}

/**
 * Calculate the length of a region, handling circular wrap-around
 *
 * @param start - Start position (0-based)
 * @param end - End position (0-based, exclusive)
 * @param seqLength - Total sequence length
 * @param circular - Whether sequence is circular
 * @returns Length of the region
 */
function getCircularRegionLength(start: number, end: number, seqLength: number, circular: boolean): number {
  if (start <= end) {
    return end - start;
  }

  if (!circular) {
    throw new Error(`Invalid region: start=${start} > end=${end} not allowed for linear sequences`);
  }

  // Circular wrap-around length
  return (seqLength - start) + end;
}

/**
 * Create a working sequence for circular designs by doubling the template
 * This allows standard algorithms to work across the origin
 *
 * @param seq - The template sequence
 * @param start - Start position (0-based)
 * @param end - End position (0-based, exclusive)
 * @param circular - Whether sequence is circular
 * @returns Object with working sequence and adjusted coordinates
 */
function prepareCircularWorkingSequence(
  seq: string,
  start: number,
  end: number,
  circular: boolean
): { workingSeq: string; workingStart: number; workingEnd: number; isWrapped: boolean } {
  if (start <= end) {
    // Normal linear region - no wrapping needed
    return {
      workingSeq: seq,
      workingStart: start,
      workingEnd: end,
      isWrapped: false,
    };
  }

  if (!circular) {
    throw new Error(`Invalid region: start=${start} > end=${end} not allowed for linear sequences`);
  }

  // Circular wrap-around: double the sequence and adjust coordinates
  // Region spans from start to seqLen, then 0 to end
  // In doubled sequence: region is at start to (seqLen + end)
  const seqLen = seq.length;
  return {
    workingSeq: seq + seq,
    workingStart: start,
    workingEnd: seqLen + end,
    isWrapped: true,
  };
}

/**
 * Adjust primer coordinates back to circular range after designing on doubled sequence
 *
 * @param primer - Primer object with start/end coordinates
 * @param seqLength - Original sequence length
 * @returns Primer with adjusted coordinates
 */
function adjustCircularCoordinates(primer: any, seqLength: number): any {
  const result = { ...primer };
  if (result.start !== undefined && result.start >= seqLength) {
    result.start = result.start % seqLength;
  }
  if (result.end !== undefined && result.end > seqLength) {
    result.end = result.end % seqLength;
    if (result.end === 0) result.end = seqLength; // End at seqLength, not 0
  }
  return result;
}

// ============================================================================
// Unified Design API
// ============================================================================

/**
 * Unified primer design function
 *
 * @param template - The template DNA sequence
 * @param spec - Design specification
 * @param options - Design options (Tm, GC, etc.)
 * @returns Design result with primers and analysis
 */
export function designUnified(
  template: string,
  spec: DesignSpec,
  options: DesignOptions = {}
): DesignResult {
  const seq = template.toUpperCase();
  const { start, end, replacement, aaHelper } = spec;
  const opts = { ...MUTAGENESIS_DEFAULTS, ...options };
  const isCircular = opts.circular !== false; // Default to circular (most templates are plasmids)

  // Validate inputs
  if (!seq || seq.length < 50) {
    throw new Error('Template must be at least 50 bp');
  }

  // Validate region bounds
  if (start < 0 || start >= seq.length) {
    throw new Error(`Invalid region: start=${start} out of bounds, template length=${seq.length}`);
  }
  if (end < 0 || end > seq.length) {
    throw new Error(`Invalid region: end=${end} out of bounds, template length=${seq.length}`);
  }

  // For linear sequences, start must be <= end
  // For circular sequences, start > end is allowed (wraps around origin)
  if (start > end && !isCircular) {
    throw new Error(`Invalid region: start=${start} > end=${end} not allowed for linear sequences. Enable circular option for wrap-around regions.`);
  }

  // Determine operation type and dispatch to appropriate handler
  const operationType = determineOperationType(spec);

  switch (operationType) {
    case 'amplify':
      return designAmplification(seq, start, end, opts as any, isCircular);

    case 'delete':
      return designDeletion(seq, start, end, opts as any, isCircular);

    case 'aa_mutation':
      return designAAMutation(seq, start, aaHelper!, opts as any);

    case 'substitute':
    case 'insert':
      return designSubstitution(seq, start, end, replacement!, opts as any, isCircular);

    default:
      throw new Error(`Unknown operation type: ${operationType}`);
  }
}

/**
 * Determine the operation type from the specification
 */
function determineOperationType(spec: DesignSpec): string {
  const { start, end, replacement, aaHelper } = spec;

  // AA mutation helper takes precedence
  if (aaHelper) {
    if (!aaHelper.newAA || aaHelper.newAA.trim() === '') {
      throw new Error('New amino acid must be specified for AA mutation');
    }
    return 'aa_mutation';
  }

  // No replacement = amplify
  if (replacement === null || replacement === undefined) {
    return 'amplify';
  }

  // Empty replacement = delete
  if (replacement === '') {
    return 'delete';
  }

  // Start === end with replacement = insert
  if (start === end) {
    return 'insert';
  }

  // Otherwise = substitute
  return 'substitute';
}

// ============================================================================
// Amplification (PCR)
// ============================================================================

/**
 * Design primers to amplify a region
 * Uses inward-facing primers (→ region ←)
 *
 * Algorithm matches MutagenesisDesigner exactly:
 * 1. Call primers() to get default pair
 * 2. Generate alternatives with generateAlternatives()
 * 3. Filter by MAX_TM_DIFF = 5°C
 * 4. If better alternative exists (higher compositeScore), use it
 * 5. Add original design to alternatives if upgraded
 *
 * For circular plasmids with wrap-around regions (start > end),
 * the region crosses the origin and is extracted as seq[start:] + seq[:end]
 */
function designAmplification(
  template: string,
  start: number,
  end: number,
  options: DesignOptions,
  isCircular: boolean = true
): DesignResult {
  const seqLen = template.length;
  const isWrapped = start > end;

  // Extract the region, handling circular wrap-around
  const regionSeq = extractCircularRegion(template, start, end, isCircular);
  const regionLength = getCircularRegionLength(start, end, seqLen, isCircular);

  if (regionSeq.length < 50) {
    throw new Error('Amplification region must be at least 50 bp');
  }

  // Get primer pair using primers() - this already generates alternatives internally
  // via generateAlternativesInternal(), matching the mutagenesis approach
  const [fwd, rev] = primers(regionSeq, {
    optimalTm: options.optimalTm || 62,
    useCompositeScore: true,
    useSmartDesign: options.useSmartDesign === true,
    exhaustiveSearch: options.exhaustiveSearch === true,
  });

  const currentScore = fwd.scoring.compositeScore ?? 0;

  // Use alternatives from primers() - same approach as MutagenesisDesigner
  // This is the unified approach matching how mutagenesis generates alternatives
  let alternativePrimers: any[] = [];
  const MAX_TM_DIFF = options.maxTmDiff || 5;

  if (fwd.alternatives && fwd.alternatives.length > 0) {
    alternativePrimers = fwd.alternatives
      // Filter to reasonable Tm diff
      .filter((alt: any) => Math.abs(alt.forward.tm - alt.reverse.tm) <= MAX_TM_DIFF)
      // Add comparison info
      .map((alt: any) => ({
        ...alt,
        isBetterThanCurrent: alt.compositeScore > currentScore,
        scoreDelta: alt.compositeScore - currentScore,
      }))
      .slice(0, 10);
  }

  // Also include alternativeOptimized from Smart Design if present
  if (fwd.alternativeOptimized && alternativePrimers.length < 10) {
    const smartAlt = fwd.alternativeOptimized;
    // Check if not already in list
    const isDuplicate = alternativePrimers.some((alt: any) =>
      alt.forward.sequence === smartAlt.forward.seq &&
      alt.reverse.sequence === smartAlt.reverse.seq
    );
    if (!isDuplicate) {
      alternativePrimers.push({
        forward: {
          sequence: smartAlt.forward.seq,
          length: smartAlt.forward.seq.length,
          tm: smartAlt.forward.tm,
          gc: smartAlt.forward.gc,
          dg: smartAlt.forward.dg,
          gcPercent: `${(smartAlt.forward.gc * 100).toFixed(1)}%`,
          hasGCClamp: /[GC]$/.test(smartAlt.forward.seq),
        },
        reverse: {
          sequence: smartAlt.reverse.seq,
          length: smartAlt.reverse.seq.length,
          tm: smartAlt.reverse.tm,
          gc: smartAlt.reverse.gc,
          dg: smartAlt.reverse.dg,
          gcPercent: `${(smartAlt.reverse.gc * 100).toFixed(1)}%`,
          hasGCClamp: /[GC]$/.test(smartAlt.reverse.seq),
        },
        tmDiff: (smartAlt as any).tmDiff || Math.abs(smartAlt.forward.tm - smartAlt.reverse.tm),
        compositeScore: (smartAlt as any).pairScore || 0,
        qualityTier: ((smartAlt as any).tmDiff || 0) <= 2 ? 'excellent' : ((smartAlt as any).tmDiff || 0) <= 5 ? 'good' : 'acceptable',
        label: '⚡ Smart Design',
        explanation: (smartAlt as any).tradeoff || '',
        isBetterThanCurrent: ((smartAlt as any).pairScore || 0) > currentScore,
        scoreDelta: ((smartAlt as any).pairScore || 0) - currentScore,
      });
    }
  }

  const tmDifference = Math.round(Math.abs(fwd.tm - rev.tm) * 10) / 10;

  // For wrapped regions, describe the position as "start-end (wrapping origin)"
  const regionDescription = isWrapped
    ? `${regionLength} bp region (${start + 1}-origin-${end})`
    : `${regionLength} bp region (${start + 1}-${end})`;

  // Calculate forward/reverse primer positions in original template coordinates
  // Forward primer starts at 'start' position in template
  // Reverse primer ends at 'end' position (or wrapped position)
  let fwdStart = start;
  let fwdEnd = (start + fwd.seq.length) % seqLen;
  if (fwdEnd === 0) fwdEnd = seqLen;

  // Reverse primer anneals at the end of the region
  // In the extracted regionSeq, reverse primer is at the 3' end
  let revEndInRegion = regionLength;
  let revStartInRegion = regionLength - rev.seq.length;

  // Map back to original template coordinates
  let revStart: number, revEnd: number;
  if (isWrapped) {
    // Region is template[start:seqLen] + template[0:end]
    // Reverse primer position in original coordinates
    const posInFirst = seqLen - start; // length of first part
    if (revStartInRegion >= posInFirst) {
      // Reverse primer is entirely in the second part (0 to end)
      revStart = revStartInRegion - posInFirst;
      revEnd = end;
    } else if (revEndInRegion <= posInFirst) {
      // Reverse primer is entirely in the first part
      revStart = start + revStartInRegion;
      revEnd = start + revEndInRegion;
    } else {
      // Reverse primer spans the origin
      revStart = start + revStartInRegion;
      revEnd = end;
    }
  } else {
    revStart = start + revStartInRegion;
    revEnd = end;
  }

  // Analyze G-Quadruplex risk for both primers
  const fwdG4 = analyzeGQuadruplex(fwd.seq);
  const revG4 = analyzeGQuadruplex(rev.seq);
  const hasG4Risk = fwdG4.hasG4Motif || revG4.hasG4Motif || fwdG4.hasGGGG || revG4.hasGGGG;

  // Build warnings list
  const warnings: string[] = [];
  if (fwdG4.severity === 'critical') {
    warnings.push(`Forward primer: ${fwdG4.message}`);
  } else if (fwdG4.severity === 'warning') {
    warnings.push(`Forward primer: ${fwdG4.message}`);
  }
  if (revG4.severity === 'critical') {
    warnings.push(`Reverse primer: ${revG4.message}`);
  } else if (revG4.severity === 'warning') {
    warnings.push(`Reverse primer: ${revG4.message}`);
  }

  // Calculate equilibrium efficiency if requested
  let equilibriumAnalysis: any = null;
  if (options.includeEquilibrium) {
    try {
      const fwdInput = { seq: fwd.seq, bindingSite: regionSeq.slice(0, fwd.seq.length) };
      const revInput = { seq: rev.seq, bindingSite: reverseComplement(regionSeq.slice(-rev.seq.length)) };
      equilibriumAnalysis = calculateEquilibriumEfficiency(fwdInput, revInput, regionSeq, {
        temperature: options.annealingTemperature || 55,
        includeOffTarget: true,
      });
    } catch (e: any) {
      // Equilibrium calculation failed, continue without it
      equilibriumAnalysis = { error: e.message };
    }
  }

  return {
    type: 'amplification',
    operation: 'amplify',
    design: 'pcr',
    originalSequence: template,
    amplifiedRegion: { start, end, length: regionLength, isWrapped },
    description: `PCR amplification of ${regionDescription}`,

    // Primers
    forward: {
      sequence: fwd.seq,
      length: fwd.seq.length,
      tm: fwd.tm,
      gc: fwd.gc,
      gcPercent: `${(fwd.gc * 100).toFixed(1)}%`,
      dg: fwd.dg,
      hasGCClamp: fwd.seq[fwd.seq.length - 1] === 'G' || fwd.seq[fwd.seq.length - 1] === 'C',
      start: fwdStart,
      end: fwdEnd,
      gQuadruplex: fwdG4,  // G-Quadruplex analysis
    },
    reverse: {
      sequence: rev.seq,
      length: rev.seq.length,
      tm: rev.tm,
      gc: rev.gc,
      gcPercent: `${(rev.gc * 100).toFixed(1)}%`,
      dg: rev.dg,
      hasGCClamp: rev.seq[rev.seq.length - 1] === 'G' || rev.seq[rev.seq.length - 1] === 'C',
      start: revStart,
      end: revEnd,
      gQuadruplex: revG4,  // G-Quadruplex analysis
    },

    // Q5 annealing temp formula: min(lower_Tm + 1, 72)
    annealingTemp: Math.round(Math.min(Math.min(fwd.tm, rev.tm) + 1, 72)),
    tmDifference,

    // Scoring - use effectiveScore for quality when critical warnings exist
    quality: fwd.scoring?.qualityTier || (currentScore >= 80 ? 'excellent' :
             currentScore >= 70 ? 'good' : 'acceptable'),
    qualityTier: fwd.scoring?.qualityTier || (currentScore >= 80 ? 'excellent' :
             currentScore >= 70 ? 'good' : 'acceptable'),
    compositeScore: currentScore as number,
    effectiveScore: fwd.scoring?.effectiveScore ?? currentScore as number,
    criticalWarnings: fwd.scoring?.criticalWarnings ?? 0,
    piecewiseScores: fwd.scoring?.piecewiseScores,
    forwardPiecewiseScores: fwd.scoring?.piecewiseScores,
    reversePiecewiseScores: rev.scoring?.piecewiseScores,
    penalty: Math.round(((fwd.scoring?.penalty || 0) + (rev.scoring?.penalty || 0)) * 10) / 10,

    // Circular wrap info
    circularWrapped: isWrapped,

    // Alternative primers for comparison (unified with mutagenesis approach)
    alternativePrimers,
    // Alternative categories for grouped view (from exhaustive search)
    alternativeCategories: fwd.alternativeCategories || null,

    // G-Quadruplex risk assessment
    hasG4Risk,
    warnings: warnings.length > 0 ? warnings : undefined,

    // Equilibrium thermodynamic analysis (if requested)
    equilibrium: equilibriumAnalysis,
  };
}

// ============================================================================
// Deletion
// ============================================================================

/**
 * Design primers to delete a region
 * Uses outward-facing primers at the deletion boundaries
 *
 * For circular plasmids with wrap-around regions (start > end),
 * the deletion crosses the origin.
 */
function designDeletion(
  template: string,
  start: number,
  end: number,
  options: DesignOptions,
  isCircular: boolean = true
): DesignResult {
  const seqLen = template.length;
  const isWrapped = start > end;
  const deleteLength = getCircularRegionLength(start, end, seqLen, isCircular);

  if (deleteLength <= 0) {
    throw new Error('Deletion length must be greater than 0');
  }

  let result: any;

  if (isWrapped) {
    // For wrap-around deletions, use the doubled sequence approach
    const { workingSeq, workingStart, workingEnd } = prepareCircularWorkingSequence(
      template, start, end, isCircular
    );

    // Design primers on the doubled sequence
    result = designDeletionPrimers(workingSeq, workingStart, deleteLength, options as any);

    // Adjust coordinates back to circular range
    if (result.forward) {
      result.forward = adjustCircularCoordinates(result.forward, seqLen);
    }
    if (result.reverse) {
      result.reverse = adjustCircularCoordinates(result.reverse, seqLen);
    }

    result.circularWrapped = true;
  } else {
    // Standard linear deletion
    result = designDeletionPrimers(template, start, deleteLength, options as any);
  }

  return {
    ...result,
    operation: 'delete',
    strategy: options.strategy || 'back-to-back',
  };
}

// ============================================================================
// Substitution / Insertion
// ============================================================================

/**
 * Design primers for substitution or insertion
 *
 * For circular plasmids with wrap-around regions (start > end),
 * the substitution crosses the origin.
 */
function designSubstitution(
  template: string,
  start: number,
  end: number,
  replacement: string,
  options: DesignOptions,
  isCircular: boolean = true
): DesignResult {
  const seqLen = template.length;
  const isInsertion = start === end;
  const isWrapped = start > end;
  const cleanReplacement = replacement.toUpperCase();

  // Check for ambiguous bases (library design)
  const isLibrary = hasAmbiguousBases(cleanReplacement);
  const librarySize = isLibrary ? countCombinations(cleanReplacement) : 1;

  let result: any;

  if (isInsertion) {
    // Insertion at a point (no wrap-around needed for insertions)
    result = designInsertionPrimers(template, start, cleanReplacement, options as any);
    result.operation = 'insert';
  } else if (isWrapped) {
    // Substitution with wrap-around: use doubled sequence approach
    const deleteLength = getCircularRegionLength(start, end, seqLen, isCircular);
    const { workingSeq, workingStart, workingEnd } = prepareCircularWorkingSequence(
      template, start, end, isCircular
    );

    // Design primers on the doubled sequence
    result = designRegionSubstitutionPrimers(workingSeq, workingStart, deleteLength, cleanReplacement, options as any);
    result.operation = 'substitute';

    // Adjust coordinates back to circular range
    if (result.forward) {
      result.forward = adjustCircularCoordinates(result.forward, seqLen);
    }
    if (result.reverse) {
      result.reverse = adjustCircularCoordinates(result.reverse, seqLen);
    }

    result.circularWrapped = true;
  } else {
    // Standard linear substitution
    const deleteLength = end - start;
    result = designRegionSubstitutionPrimers(template, start, deleteLength, cleanReplacement, options as any);
    result.operation = 'substitute';
  }

  // Add library info if applicable
  if (isLibrary) {
    result.isLibrary = true;
    result.librarySize = librarySize;
    result.librarySequence = cleanReplacement;
  }

  result.strategy = options.strategy || 'back-to-back';

  return result;
}

// ============================================================================
// Amino Acid Mutation
// ============================================================================

/**
 * Design primers for amino acid mutation with codon optimization
 */
function designAAMutation(
  template: string,
  start: number,
  aaHelper: AAHelper,
  options: DesignOptions
): DesignResult {
  const { newAA, organism = 'ecoli', orfStart = 1 } = aaHelper;

  // Calculate codon position from nucleotide position
  // If start is provided as a codon position (1-based), use it directly
  // Otherwise, calculate from nucleotide position
  let codonPosition: number;

  if (aaHelper.codonPosition) {
    codonPosition = aaHelper.codonPosition;
  } else {
    // Assume start is the nucleotide position of the codon start
    const orfOffset = (orfStart - 1);
    codonPosition = Math.floor((start - orfOffset) / 3) + 1;
  }

  // Use existing codon change primer design
  const result = designCodonChangePrimers(template, codonPosition, newAA, {
    ...options,
    organism,
    orfStart,
  } as any);

  return {
    ...result,
    operation: 'aa_mutation',
    strategy: options.strategy || 'back-to-back',
  };
}

// ============================================================================
// Batch Design
// ============================================================================

/**
 * Design primers for multiple modifications
 *
 * @param template - Template sequence
 * @param specs - Array of design specifications
 * @param options - Shared design options
 * @returns Array of design results
 */
export function designBatch(
  template: string,
  specs: DesignSpec[],
  options: DesignOptions = {}
): BatchResult[] {
  const results: BatchResult[] = [];

  for (const spec of specs) {
    try {
      const result = designUnified(template, spec, options);
      results.push({
        ...result,
        success: true,
        originalSpec: spec,
      } as BatchResult);
    } catch (error: any) {
      results.push({
        success: false,
        error: error.message,
        originalSpec: spec,
      } as any);
    }
  }

  return results;
}

// ============================================================================
// Helper: Parse mutation notation to unified spec
// ============================================================================

/**
 * Parse mutation notation (Y66W, del100-105, ins50_ACGT) to unified spec
 *
 * @param notation - Mutation notation
 * @param context - Context information
 * @returns Unified specification
 */
export function parseNotationToSpec(
  notation: string,
  context: { orfStart?: number } = {}
): DesignSpec {
  const { orfStart = 1 } = context;
  const upper = notation.toUpperCase().trim();

  // AA mutation: Y66W
  const aaMatch = upper.match(/^([ARNDCEQGHILKMFPSTWYV*])(\d+)([ARNDCEQGHILKMFPSTWYV*])$/);
  if (aaMatch) {
    const codonPosition = parseInt(aaMatch[2]);
    const nucPosition = (orfStart - 1) + (codonPosition - 1) * 3;

    return {
      start: nucPosition,
      end: nucPosition + 3,
      replacement: null, // Will be resolved by aaHelper
      aaHelper: {
        newAA: aaMatch[3],
        codonPosition,
        originalAA: aaMatch[1],
      },
      notation: `${aaMatch[1]}${aaMatch[2]}${aaMatch[3]}`,
      notationType: 'aa_mutation',
    };
  }

  // Deletion: del100-105 or del100
  const delMatch = upper.match(/^DEL(\d+)(?:-(\d+))?$/);
  if (delMatch) {
    const start = parseInt(delMatch[1]) - 1; // Convert to 0-based
    const end = delMatch[2] ? parseInt(delMatch[2]) : start + 1;

    return {
      start,
      end,
      replacement: '',
      notation: delMatch[2] ? `del${delMatch[1]}-${delMatch[2]}` : `del${delMatch[1]}`,
      notationType: 'deletion',
    };
  }

  // Insertion: ins50_ACGT
  const insMatch = upper.match(/^INS(\d+)_([ATGC]+)$/);
  if (insMatch) {
    const position = parseInt(insMatch[1]) - 1; // Convert to 0-based

    return {
      start: position,
      end: position,
      replacement: insMatch[2],
      notation: `ins${insMatch[1]}_${insMatch[2]}`,
      notationType: 'insertion',
    };
  }

  // Substitution: sub100-105_ACGT
  const subMatch = upper.match(/^SUB(\d+)-(\d+)_([ATGC]+)$/);
  if (subMatch) {
    const start = parseInt(subMatch[1]) - 1;
    const end = parseInt(subMatch[2]);

    return {
      start,
      end,
      replacement: subMatch[3],
      notation: `sub${subMatch[1]}-${subMatch[2]}_${subMatch[3]}`,
      notationType: 'substitution',
    };
  }

  throw new Error(`Could not parse mutation notation: ${notation}`);
}

/**
 * Parse batch text input to array of specs
 *
 * Format:
 * ORF: 1
 * Y66W
 * S65T
 * del100-105
 * ins50_ACGT
 * # This is a comment
 */
export function parseBatchText(text: string, defaultOptions: { orfStart?: number } = {}): ParsedBatch {
  const lines = text.split('\n');
  const specs: (DesignSpec | { error: string; originalLine: string })[] = [];
  let orfStart = defaultOptions.orfStart || 1;

  for (const line of lines) {
    const trimmed = line.trim();

    // Skip empty lines and comments
    if (!trimmed || trimmed.startsWith('#')) continue;

    // ORF directive
    const orfMatch = trimmed.match(/^ORF:\s*(\d+)/i);
    if (orfMatch) {
      orfStart = parseInt(orfMatch[1]);
      continue;
    }

    try {
      const spec = parseNotationToSpec(trimmed, { orfStart });
      specs.push(spec);
    } catch (error: any) {
      specs.push({
        error: error.message,
        originalLine: trimmed,
      });
    }
  }

  return { specs, orfStart };
}

// ============================================================================
// Exports
// ============================================================================

export {
  MUTAGENESIS_DEFAULTS,
  MUTATION_TYPES,
  CODON_TABLE,
  CODON_TO_AA,
  AA_NAMES,
  selectOptimalCodon,
  analyzePrimerPair,
};

export default designUnified;
