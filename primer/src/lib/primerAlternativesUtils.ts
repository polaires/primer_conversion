/**
 * Primer Alternatives Utilities
 *
 * Shared utilities for comparing, filtering, and displaying primer alternatives.
 * Used by both UnifiedPrimerDesigner and MutagenesisDesigner.
 */

// =============================================================================
// Type Definitions
// =============================================================================

/**
 * Trade-off type
 */
export type TradeOffType = 'better' | 'worse' | 'warn' | 'info';

/**
 * Trade-off comparison object
 */
export interface TradeOff {
  type: TradeOffType;
  label: string;
  delta: string;
  detail?: string;
  priority: number;
}

/**
 * Primer data structure (simplified)
 */
export interface PrimerData {
  sequence?: string;
  seq?: string;
  length?: number;
  tm?: number;
  gc?: number;
  gcPercent?: string;
  dg?: number;
  terminal3DG?: number;
  terminalDG?: number;
  hasGCClamp?: boolean;
}

/**
 * Primer pair alternative
 */
export interface PrimerAlternative {
  forward?: PrimerData;
  reverse?: PrimerData;
  compositeScore?: number;
  score?: number;
  qualityTier?: string;
  penalty?: number;
  tmDiff?: number;
  label?: string | null;
  explanation?: string | null;
  design?: string;
  isBetterThanCurrent?: boolean;
  scoreDelta?: number;
}

/**
 * Normalized alternative with all fields populated
 */
export interface NormalizedAlternative {
  originalIdx: number;
  forward: {
    sequence: string;
    length: number;
    tm: number;
    gc: number;
    gcPercent: string;
    dg: number | null;
    hasGCClamp: boolean;
  };
  reverse: {
    sequence: string;
    length: number;
    tm: number;
    gc: number;
    gcPercent: string;
    dg: number | null;
    hasGCClamp: boolean;
  };
  compositeScore: number;
  qualityTier: string;
  penalty: number;
  tmDiff: number;
  label: string | null;
  explanation: string | null;
  design: string;
  isBetterThanCurrent: boolean;
  scoreDelta: number;
}

/**
 * Alternative with trade-off analysis
 */
export interface AlternativeWithTradeOffs extends NormalizedAlternative {
  tradeOffs: TradeOff[];
}

/**
 * Filter configuration
 */
export interface FilterConfig {
  minScore?: number;
  maxTmDiff?: number;
  requireGcClamp?: boolean;
}

/**
 * Sort configuration
 */
export interface SortConfig {
  field?: 'score' | 'tmDiff' | 'fwdTm' | 'revTm' | 'fwdLen' | 'revLen' | 'totalLen';
  direction?: 'asc' | 'desc';
}

/**
 * Display preparation options
 */
export interface DisplayOptions {
  filters?: FilterConfig;
  sortConfig?: SortConfig;
  maxDisplay?: number;
  showAll?: boolean;
}

/**
 * Display result
 */
export interface DisplayResult {
  displayed: AlternativeWithTradeOffs[];
  totalCount: number;
  filteredCount: number;
  hiddenCount: number;
}

/**
 * CSV export options
 */
export interface CSVExportOptions {
  includeHeader?: boolean;
}

/**
 * Quality tier type
 */
export type QualityTier = 'excellent' | 'good' | 'acceptable' | 'poor';

// =============================================================================
// Trade-off Generation
// =============================================================================

/**
 * Generate trade-off comparisons with delta notation for alternative primers
 * @param alt - Alternative primer pair
 * @param orig - Original/current primer pair
 * @returns Array of trade-off objects with type, label, delta, and priority
 */
export function generateTradeOffs(alt: PrimerAlternative, orig: PrimerAlternative): TradeOff[] {
  const tradeOffs: TradeOff[] = [];

  // Score delta (most important - show first if significant)
  const altScore = alt.compositeScore || alt.score || 0;
  const origScore = orig.compositeScore || orig.score || 0;
  const scoreDelta = altScore - origScore;

  if (Math.abs(scoreDelta) >= 1) {
    tradeOffs.push({
      type: scoreDelta > 0 ? 'better' : 'worse',
      label: 'Score',
      delta: scoreDelta > 0 ? `+${Math.round(scoreDelta)}` : `${Math.round(scoreDelta)}`,
      priority: scoreDelta > 2 ? 0 : 1,
    });
  }

  // Tm difference comparison
  const altTmDiff = Math.abs((alt.forward?.tm || 0) - (alt.reverse?.tm || 0));
  const origTmDiff = Math.abs((orig.forward?.tm || 0) - (orig.reverse?.tm || 0));
  const tmDiffDelta = origTmDiff - altTmDiff; // Positive = alt is better (lower diff)

  if (Math.abs(tmDiffDelta) >= 0.5) {
    tradeOffs.push({
      type: tmDiffDelta > 0 ? 'better' : 'worse',
      label: 'Tm match',
      delta: `${altTmDiff.toFixed(1)}°`,
      detail: tmDiffDelta > 0 ? 'tighter' : 'wider',
      priority: 2,
    });
  }

  // GC Clamp comparison
  const altGcClamps = (alt.forward?.hasGCClamp ? 1 : 0) + (alt.reverse?.hasGCClamp ? 1 : 0);
  const origGcClamps = (orig.forward?.hasGCClamp ? 1 : 0) + (orig.reverse?.hasGCClamp ? 1 : 0);
  const gcDelta = altGcClamps - origGcClamps;

  if (gcDelta !== 0) {
    tradeOffs.push({
      type: gcDelta > 0 ? 'better' : 'worse',
      label: 'GC clamp',
      delta: `${altGcClamps}/2`,
      detail: gcDelta > 0 ? `+${gcDelta}` : `${gcDelta}`,
      priority: 3,
    });
  } else if (altGcClamps < 2) {
    // Show as warning if both missing GC clamps
    tradeOffs.push({
      type: 'warn',
      label: 'GC clamp',
      delta: `${altGcClamps}/2`,
      priority: 6,
    });
  }

  // 3' Terminal ΔG comparison (if available)
  const altFwdDg = alt.forward?.terminal3DG || alt.forward?.terminalDG || alt.forward?.dg;
  const altRevDg = alt.reverse?.terminal3DG || alt.reverse?.terminalDG || alt.reverse?.dg;
  const origFwdDg = orig.forward?.terminal3DG || orig.forward?.terminalDG || orig.forward?.dg;
  const origRevDg = orig.reverse?.terminal3DG || orig.reverse?.terminalDG || orig.reverse?.dg;

  if (altFwdDg && altRevDg && origFwdDg && origRevDg) {
    const altWorstDg = Math.max(altFwdDg, altRevDg);
    const origWorstDg = Math.max(origFwdDg, origRevDg);
    const dgDelta = altWorstDg - origWorstDg;

    if (Math.abs(dgDelta) >= 0.5 && origWorstDg > -11 && origWorstDg < -4) {
      const isBetter = dgDelta < 0 && altWorstDg >= -11;
      const isWorse = dgDelta > 0 && altWorstDg > -6;
      if (isBetter || isWorse) {
        tradeOffs.push({
          type: isBetter ? 'better' : 'worse',
          label: "3' ΔG",
          delta: `${altWorstDg.toFixed(1)}`,
          detail: isBetter ? '↑stable' : '↓stable',
          priority: 4,
        });
      }
    }
  }

  // Length comparison
  const altAvgLen = ((alt.forward?.length || 0) + (alt.reverse?.length || 0)) / 2;
  const origAvgLen = ((orig.forward?.length || 0) + (orig.reverse?.length || 0)) / 2;
  const lenDelta = altAvgLen - origAvgLen;
  const inOptimalRange = altAvgLen >= 18 && altAvgLen <= 24;
  const origInOptimal = origAvgLen >= 18 && origAvgLen <= 24;

  if (Math.abs(lenDelta) >= 2) {
    if (inOptimalRange && !origInOptimal) {
      tradeOffs.push({
        type: 'better',
        label: 'Length',
        delta: `${Math.round(altAvgLen)}bp`,
        detail: '→ optimal',
        priority: 5,
      });
    } else if (!inOptimalRange && origInOptimal) {
      tradeOffs.push({
        type: 'worse',
        label: 'Length',
        delta: `${Math.round(altAvgLen)}bp`,
        detail: '← optimal',
        priority: 5,
      });
    } else if (!inOptimalRange) {
      tradeOffs.push({
        type: 'info',
        label: 'Length',
        delta: `${Math.round(altAvgLen)}bp`,
        priority: 7,
      });
    }
  }

  // Sort by priority (lower = more important)
  tradeOffs.sort((a, b) => a.priority - b.priority);

  return tradeOffs;
}

// =============================================================================
// Filtering & Sorting
// =============================================================================

/**
 * Filter alternatives based on criteria
 * @param alternatives - Array of alternative primer pairs
 * @param filters - Filter criteria
 * @returns Filtered alternatives
 */
export function filterAlternatives(
  alternatives: NormalizedAlternative[],
  filters: FilterConfig = {}
): NormalizedAlternative[] {
  const { minScore = 0, maxTmDiff = 10, requireGcClamp = false } = filters;

  return alternatives.filter(alt => {
    // Score filter
    const score = alt.compositeScore || 0;
    if (minScore > 0 && score < minScore) return false;

    // Tm difference filter
    const tmDiff = Math.abs((alt.forward?.tm || 0) - (alt.reverse?.tm || 0));
    if (tmDiff > maxTmDiff) return false;

    // GC clamp filter
    if (requireGcClamp) {
      const gcClamps = (alt.forward?.hasGCClamp ? 1 : 0) + (alt.reverse?.hasGCClamp ? 1 : 0);
      if (gcClamps < 2) return false;
    }

    return true;
  });
}

/**
 * Sort alternatives by specified field
 * @param alternatives - Array of alternative primer pairs
 * @param sortConfig - Sort configuration
 * @returns Sorted alternatives (new array)
 */
export function sortAlternatives(
  alternatives: NormalizedAlternative[],
  sortConfig: SortConfig = {}
): NormalizedAlternative[] {
  const { field = 'score', direction = 'desc' } = sortConfig;

  return [...alternatives].sort((a, b) => {
    let aVal: number, bVal: number;

    switch (field) {
      case 'score':
        aVal = a.compositeScore || 0;
        bVal = b.compositeScore || 0;
        break;
      case 'tmDiff':
        aVal = Math.abs((a.forward?.tm || 0) - (a.reverse?.tm || 0));
        bVal = Math.abs((b.forward?.tm || 0) - (b.reverse?.tm || 0));
        break;
      case 'fwdTm':
        aVal = a.forward?.tm || 0;
        bVal = b.forward?.tm || 0;
        break;
      case 'revTm':
        aVal = a.reverse?.tm || 0;
        bVal = b.reverse?.tm || 0;
        break;
      case 'fwdLen':
        aVal = a.forward?.length || 0;
        bVal = b.forward?.length || 0;
        break;
      case 'revLen':
        aVal = a.reverse?.length || 0;
        bVal = b.reverse?.length || 0;
        break;
      case 'totalLen':
        aVal = (a.forward?.length || 0) + (a.reverse?.length || 0);
        bVal = (b.forward?.length || 0) + (b.reverse?.length || 0);
        break;
      default:
        aVal = a.compositeScore || 0;
        bVal = b.compositeScore || 0;
    }

    const diff = aVal - bVal;
    return direction === 'desc' ? -diff : diff;
  });
}

// =============================================================================
// Display Helpers
// =============================================================================

/**
 * Thresholds for Tm difference color coding
 */
export const TM_DIFF_THRESHOLDS = {
  EXCELLENT: 1.0,  // Green - ideal
  GOOD: 1.5,       // Light green
  WARNING: 3.0,    // Yellow - acceptable
  // Above WARNING is red/poor
};

/**
 * Get color class for ΔTm value based on thresholds
 * @param tmDiff - Tm difference in °C
 * @returns CSS class name
 */
export function getTmDiffColorClass(tmDiff: number): string {
  if (tmDiff <= TM_DIFF_THRESHOLDS.EXCELLENT) return 'tm-excellent';
  if (tmDiff <= TM_DIFF_THRESHOLDS.GOOD) return 'tm-good';
  if (tmDiff <= TM_DIFF_THRESHOLDS.WARNING) return 'tm-warning';
  return 'tm-poor';
}

/**
 * Get quality tier from composite score
 * @param score - Composite score (0-100)
 * @returns Quality tier name
 */
export function getQualityTier(score: number): QualityTier {
  if (score >= 80) return 'excellent';
  if (score >= 70) return 'good';
  if (score >= 60) return 'acceptable';
  return 'poor';
}

/**
 * Get quality color for a given tier
 * @param tier - Quality tier (excellent, good, acceptable, poor)
 * @returns Hex color code
 */
export function getQualityColor(tier: QualityTier): string {
  switch (tier) {
    case 'excellent': return '#22c55e';
    case 'good': return '#84cc16';
    case 'acceptable': return '#eab308';
    case 'poor': return '#ef4444';
    default: return '#6b7280';
  }
}

/**
 * Format primer sequences for display
 * @param sequence - Full primer sequence
 * @param maxLen - Maximum length before truncation
 * @returns Formatted sequence (possibly truncated)
 */
export function formatSequence(sequence: string, maxLen: number = 30): string {
  if (!sequence) return '';
  if (sequence.length <= maxLen) return sequence;
  return `${sequence.slice(0, maxLen)}...`;
}

/**
 * Calculate Q5 annealing temperature
 * @param fwdTm - Forward primer Tm
 * @param revTm - Reverse primer Tm
 * @returns Recommended annealing temperature
 */
export function calculateAnnealingTemp(fwdTm: number, revTm: number): number {
  const lowerTm = Math.min(fwdTm, revTm);
  return Math.round(Math.min(lowerTm + 1, 72));
}

// =============================================================================
// CSV Export
// =============================================================================

/**
 * Generate CSV content for alternatives export
 * @param alternatives - Array of alternatives to export
 * @param options - Export options
 * @returns CSV content
 */
export function generateAlternativesCSV(
  alternatives: PrimerAlternative[],
  options: CSVExportOptions = {}
): string {
  const { includeHeader = true } = options;

  const rows: string[] = [];

  if (includeHeader) {
    rows.push([
      'Rank',
      'Score',
      'Fwd Sequence',
      'Fwd Tm',
      'Fwd Length',
      'Fwd GC Clamp',
      'Rev Sequence',
      'Rev Tm',
      'Rev Length',
      'Rev GC Clamp',
      'Tm Diff',
      'Label',
    ].join(','));
  }

  alternatives.forEach((alt, idx) => {
    const tmDiff = Math.abs((alt.forward?.tm || 0) - (alt.reverse?.tm || 0));
    rows.push([
      idx + 1,
      alt.compositeScore || alt.score || '',
      alt.forward?.sequence || '',
      alt.forward?.tm?.toFixed(1) || '',
      alt.forward?.length || '',
      alt.forward?.hasGCClamp ? 'Yes' : 'No',
      alt.reverse?.sequence || '',
      alt.reverse?.tm?.toFixed(1) || '',
      alt.reverse?.length || '',
      alt.reverse?.hasGCClamp ? 'Yes' : 'No',
      tmDiff.toFixed(1),
      alt.label || '',
    ].join(','));
  });

  return rows.join('\n');
}

// =============================================================================
// Alternative Normalization
// =============================================================================

/**
 * Normalize alternative data structure for consistent handling
 * Different sources (primers(), mutagenesis) may have slightly different formats
 *
 * @param alt - Raw alternative from any source
 * @param originalIdx - Original index in the source array
 * @returns Normalized alternative object
 */
export function normalizeAlternative(alt: PrimerAlternative, originalIdx: number = 0): NormalizedAlternative {
  // Extract forward/reverse with fallbacks
  const forward = alt.forward || {};
  const reverse = alt.reverse || {};

  // Normalize GC clamp detection
  const fwdHasGCClamp = forward.hasGCClamp ?? /[GC]$/i.test(forward.sequence || '');
  const revHasGCClamp = reverse.hasGCClamp ?? /[GC]$/i.test(reverse.sequence || '');

  return {
    originalIdx,
    forward: {
      sequence: forward.sequence || forward.seq || '',
      length: forward.length || (forward.sequence || forward.seq || '').length,
      tm: forward.tm || 0,
      gc: forward.gc || 0,
      gcPercent: forward.gcPercent || (forward.gc ? `${(forward.gc * 100).toFixed(1)}%` : ''),
      dg: forward.dg || forward.terminal3DG || forward.terminalDG || null,
      hasGCClamp: fwdHasGCClamp,
    },
    reverse: {
      sequence: reverse.sequence || reverse.seq || '',
      length: reverse.length || (reverse.sequence || reverse.seq || '').length,
      tm: reverse.tm || 0,
      gc: reverse.gc || 0,
      gcPercent: reverse.gcPercent || (reverse.gc ? `${(reverse.gc * 100).toFixed(1)}%` : ''),
      dg: reverse.dg || reverse.terminal3DG || reverse.terminalDG || null,
      hasGCClamp: revHasGCClamp,
    },
    compositeScore: alt.compositeScore || alt.score || 0,
    qualityTier: alt.qualityTier || getQualityTier(alt.compositeScore || alt.score || 0),
    penalty: alt.penalty || 0,
    tmDiff: alt.tmDiff || Math.abs((forward.tm || 0) - (reverse.tm || 0)),
    label: alt.label || null,
    explanation: alt.explanation || null,
    design: alt.design || 'overlapping',
    // Comparison flags
    isBetterThanCurrent: alt.isBetterThanCurrent || false,
    scoreDelta: alt.scoreDelta || 0,
  };
}

/**
 * Prepare alternatives list for display
 * Normalizes, filters, sorts, and adds display metadata
 *
 * @param alternatives - Raw alternatives array
 * @param currentDesign - Current/original design for comparison
 * @param options - Processing options
 * @returns Processed alternatives ready for display
 */
export function prepareAlternativesForDisplay(
  alternatives: PrimerAlternative[],
  currentDesign: PrimerAlternative,
  options: DisplayOptions = {}
): DisplayResult {
  const {
    filters = {},
    sortConfig = { field: 'score', direction: 'desc' },
    maxDisplay = 10,
    showAll = false,
  } = options;

  // Normalize all alternatives
  const normalized = alternatives.map((alt, idx) => normalizeAlternative(alt, idx));

  // Filter
  const filtered = filterAlternatives(normalized, filters);

  // Sort
  const sorted = sortAlternatives(filtered, sortConfig);

  // Add trade-offs and comparison info
  const withTradeOffs: AlternativeWithTradeOffs[] = sorted.map(alt => ({
    ...alt,
    tradeOffs: generateTradeOffs(alt, currentDesign),
    isBetterThanCurrent: (alt.compositeScore || 0) > (currentDesign.compositeScore || currentDesign.score || 0),
  }));

  // Limit display count
  const displayed = showAll ? withTradeOffs : withTradeOffs.slice(0, maxDisplay);

  return {
    displayed,
    totalCount: alternatives.length,
    filteredCount: filtered.length,
    hiddenCount: withTradeOffs.length - displayed.length,
  };
}
