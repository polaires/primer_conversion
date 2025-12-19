import React, { useMemo, useState, useRef, useEffect } from 'react';
import { SeqViz } from 'seqviz';
import HairpinDiagram, { HairpinBadge } from './HairpinDiagram.jsx';
import PrimerStructureViewer from './PrimerStructureViewer.jsx';
import { foldSequence } from '../lib/fold.js';
import { classify3PrimeStructureSeverity } from '../lib/smartPrimers.js';

/**
 * PrimerOnTemplateViewer - Visualize designed primers on the target gene
 * Light theme version with improved UI
 */

// Theme type
interface Theme {
  bg: string;
  bgSecondary: string;
  bgTertiary: string;
  border: string;
  borderStrong: string;
  text: string;
  textSecondary: string;
  textMuted: string;
  success: string;
  successBg: string;
  warning: string;
  warningBg: string;
  danger: string;
  dangerBg: string;
  info: string;
  infoBg: string;
  forward: string;
  reverse: string;
  amplicon: string;
  mutation: string;
}

// Light theme colors
const theme: Theme = {
  bg: '#ffffff',
  bgSecondary: '#f8fafc',
  bgTertiary: '#f1f5f9',
  border: '#e2e8f0',
  borderStrong: '#cbd5e1',
  text: '#1e293b',
  textSecondary: '#475569',
  textMuted: '#94a3b8',
  success: '#16a34a',
  successBg: '#dcfce7',
  warning: '#ca8a04',
  warningBg: '#fef9c3',
  danger: '#dc2626',
  dangerBg: '#fee2e2',
  info: '#2563eb',
  infoBg: '#dbeafe',
  forward: '#2563eb',
  reverse: '#7c3aed',
  amplicon: '#16a34a',
  mutation: '#ea580c',
};

// Primer type
interface Primer {
  sequence: string;
  tm?: number | string;
  gc?: number;
  gcPercent?: string;
  dg?: number;
  start?: number;
  end?: number;
}

// Fold result type (from fold.js)
interface FoldResult {
  e: number;
  ij?: number[][];
  structure?: string;
}

// Severity classification result type
interface SeverityResult {
  level: 'none' | 'info' | 'low' | 'moderate' | 'warning' | 'critical';
  shouldWarn: boolean;
  tooltip?: string;
}

// Binding result type
interface BindingResult {
  start: number;
  end: number;
  matchLength: number;
  score: number;
  method: string;
  offset?: number;
}

// Binding options type
interface BindingOptions {
  positionHint?: {
    start: number;
    end: number;
  };
  mutationPosition?: number | null;
  isMutagenesis?: boolean;
}

// SeqViz annotation type
interface Annotation {
  name: string;
  start: number;
  end: number;
  direction: 1 | -1;
  color: string;
  type: string;
}

// Reverse complement helper
function reverseComplement(seq: string): string {
  const comp: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G', a: 't', t: 'a', g: 'c', c: 'g' };
  return seq.split('').reverse().map(c => comp[c] || c).join('');
}

/**
 * IMPROVED: Better primer binding with dual-anchor strategy
 *
 * For mutagenesis primers:
 * - Check BOTH 5' and 3' ends (not just 3')
 * - Use position hints when available
 * - Use mutation position as anchor
 * - Score based on combined anchor matches with 3' weighting
 */
function findPrimerBinding(
  template: string,
  primerSeq: string,
  isReverse: boolean = false,
  options: BindingOptions = {}
): BindingResult | null {
  if (!template || !primerSeq) return null;

  const { positionHint, mutationPosition, isMutagenesis = false } = options;

  const tmpl = template.toUpperCase();
  const primer = primerSeq.toUpperCase();
  const searchSeq = isReverse ? reverseComplement(primer) : primer;

  // Priority 1: Use explicit position hints
  if (positionHint?.start !== undefined && positionHint?.end !== undefined) {
    return {
      start: positionHint.start,
      end: positionHint.end,
      matchLength: primerSeq.length,
      score: 1.0,
      method: 'explicit_position'
    };
  }

  // Priority 2: Exact match
  const exactIdx = tmpl.indexOf(searchSeq);
  if (exactIdx !== -1) {
    return {
      start: exactIdx,
      end: exactIdx + searchSeq.length,
      matchLength: searchSeq.length,
      score: 1.0,
      method: 'exact'
    };
  }

  // Priority 3: For mutagenesis - use mutation position as anchor
  if (isMutagenesis && mutationPosition !== undefined && mutationPosition !== null) {
    // For back-to-back (Q5 SDM) design:
    // - Forward primer: contains mutation, starts ~10bp before mutation, extends downstream
    // - Reverse primer: does NOT contain mutation, ends just before forward starts
    //
    // For overlapping (QuikChange) design:
    // - Both primers contain the mutation in the middle

    let estimatedStart: number;
    if (isReverse) {
      // Reverse primer should END before the mutation site
      // Assuming forward primer starts ~10bp before mutation, reverse ends there
      const forwardFlankBefore = 10; // typical 5' flank on forward primer
      const estimatedEnd = mutationPosition - forwardFlankBefore;
      estimatedStart = estimatedEnd - searchSeq.length;
    } else {
      // Forward primer: mutation is typically ~10bp from 5' end (30-40% into primer)
      estimatedStart = mutationPosition - Math.floor(searchSeq.length * 0.4);
    }

    const clampedStart = Math.max(0, Math.min(tmpl.length - searchSeq.length, estimatedStart));

    return {
      start: clampedStart,
      end: clampedStart + searchSeq.length,
      matchLength: searchSeq.length,
      score: 0.85,
      method: 'mutation_anchor'
    };
  }

  // Priority 4: Dual-anchor fuzzy matching
  // Check both 5' and 3' ends simultaneously
  const minAnchorLen = Math.min(8, Math.floor(primer.length * 0.4));

  for (let anchorLen = Math.floor(primer.length * 0.5); anchorLen >= minAnchorLen; anchorLen--) {
    const anchor5 = searchSeq.slice(0, anchorLen);
    const anchor3 = searchSeq.slice(-anchorLen);

    // Find 5' anchor positions
    const positions5: number[] = [];
    let idx = tmpl.indexOf(anchor5);
    while (idx !== -1) {
      positions5.push(idx);
      idx = tmpl.indexOf(anchor5, idx + 1);
    }

    // Find 3' anchor positions
    const positions3: number[] = [];
    idx = tmpl.indexOf(anchor3);
    while (idx !== -1) {
      positions3.push(idx);
      idx = tmpl.indexOf(anchor3, idx + 1);
    }

    // Look for consistent pairs (5' and 3' at expected distance)
    const expectedDistance = searchSeq.length - anchorLen;

    for (const pos5 of positions5) {
      for (const pos3 of positions3) {
        const actualDistance = pos3 - pos5;

        // Allow some tolerance for indels
        if (Math.abs(actualDistance - expectedDistance) <= 2) {
          return {
            start: pos5,
            end: pos3 + anchorLen,
            matchLength: (pos3 + anchorLen) - pos5,
            score: (anchorLen * 2) / searchSeq.length,
            method: 'dual_anchor'
          };
        }
      }
    }
  }

  // Priority 5: Single 3' anchor (most critical for PCR)
  for (let anchorLen = primer.length - 1; anchorLen >= minAnchorLen; anchorLen--) {
    const anchor3 = searchSeq.slice(-anchorLen);
    const idx3 = tmpl.indexOf(anchor3);

    if (idx3 !== -1) {
      const inferredStart = idx3 - (searchSeq.length - anchorLen);
      const start = Math.max(0, inferredStart);

      return {
        start,
        end: idx3 + anchorLen,
        matchLength: (idx3 + anchorLen) - start,
        score: anchorLen / searchSeq.length,
        method: 'anchor_3prime'
      };
    }
  }

  // Priority 6: Scoring-based alignment with 3' weighting
  let bestMatch: BindingResult | null = null;
  let bestScore = 0;

  for (let i = 0; i <= tmpl.length - searchSeq.length; i++) {
    let score = 0;
    let matches3Prime = 0;

    for (let j = 0; j < searchSeq.length; j++) {
      if (tmpl[i + j] === searchSeq[j]) {
        // Weight 3' end higher (last 10 bases critical for extension)
        const is3Prime = j >= searchSeq.length - 10;
        score += is3Prime ? 2 : 1;
        if (is3Prime) matches3Prime++;
      }
    }

    // Require good 3' match for valid binding
    if (score > bestScore && matches3Prime >= 7) {
      bestScore = score;
      bestMatch = {
        start: i,
        end: i + searchSeq.length,
        matchLength: searchSeq.length,
        score: score / (searchSeq.length + 10), // Normalize
        method: 'weighted_alignment'
      };
    }
  }

  return bestMatch;
}

/**
 * Clean label component - simplified design without connector lines
 */
interface PrimerLabelProps {
  x: number;
  y: number;
  label: string;
  sublabel: string;
  color: string;
  anchor?: 'start' | 'middle' | 'end';
  position?: 'above' | 'below';
}

function PrimerLabel({
  x,
  y,
  label,
  sublabel,
  color,
  anchor = 'middle',
  position = 'above'
}: PrimerLabelProps): React.ReactElement {
  const isAbove = position === 'above';
  const barHeight = 24;

  // For "above": labels above the bar (negative offset from bar top)
  // For "below": labels below the bar (positive offset, accounting for bar height)
  const mainLabelY = isAbove ? y - 22 : y + barHeight + 16;
  const sublabelY = isAbove ? y - 8 : y + barHeight + 30;

  return (
    <g>
      {/* Main label */}
      <text
        x={x}
        y={mainLabelY}
        fill={color}
        fontSize="12"
        fontFamily="system-ui, -apple-system, sans-serif"
        fontWeight="600"
        textAnchor={anchor}
      >
        {label}
      </text>

      {/* Sublabel */}
      <text
        x={x}
        y={sublabelY}
        fill={theme.textSecondary}
        fontSize="10"
        fontFamily="system-ui, sans-serif"
        fontWeight="500"
        textAnchor={anchor}
      >
        {sublabel}
      </text>
    </g>
  );
}

/**
 * Primer block with integrated direction arrow
 */
interface PrimerBlockProps {
  x: number;
  width: number;
  y: number;
  color: string;
  direction: 'forward' | 'reverse';
  onHover?: () => void;
  onLeave?: () => void;
}

function PrimerBlock({
  x,
  width,
  y,
  color,
  direction,
  onHover,
  onLeave
}: PrimerBlockProps): React.ReactElement {
  const arrowSize = 12;
  const height = 24;
  const radius = 6;

  // Build path for rounded rect with arrow
  const isForward = direction === 'forward';

  const path = isForward
    ? `M ${x + radius} ${y}
       h ${width - radius - arrowSize}
       l ${arrowSize} ${height / 2}
       l ${-arrowSize} ${height / 2}
       h ${-(width - radius - arrowSize)}
       a ${radius} ${radius} 0 0 1 ${-radius} ${-radius}
       v ${-(height - radius * 2)}
       a ${radius} ${radius} 0 0 1 ${radius} ${-radius}
       z`
    : `M ${x + arrowSize} ${y}
       h ${width - radius - arrowSize}
       a ${radius} ${radius} 0 0 1 ${radius} ${radius}
       v ${height - radius * 2}
       a ${radius} ${radius} 0 0 1 ${-radius} ${radius}
       h ${-(width - radius - arrowSize)}
       l ${-arrowSize} ${-height / 2}
       z`;

  return (
    <path
      d={path}
      fill={color}
      onMouseEnter={onHover}
      onMouseLeave={onLeave}
      style={{ cursor: 'pointer', filter: 'drop-shadow(0 2px 4px rgba(0,0,0,0.1))' }}
    />
  );
}

/**
 * SeqViz-based Template Visualization
 * Uses the battle-tested SeqViz library for proper primer positioning
 */
interface TemplateVisualizationProps {
  template: string;
  forward?: Primer | null;
  reverse?: Primer | null;
  mutationPosition?: number | null;
  mutationLength?: number;
  width?: number;
  showSequence?: boolean;
  compact?: boolean;
  isMutagenesis?: boolean;
}

function TemplateVisualization({
  template,
  forward,
  reverse,
  mutationPosition = null,
  mutationLength = 1,
  width = 800,
  showSequence = true,
  compact = false,
  isMutagenesis = false
}: TemplateVisualizationProps): React.ReactElement {
  const height = compact ? 200 : 250;

  // Build annotations for SeqViz
  const annotations = useMemo((): Annotation[] => {
    const annots: Annotation[] = [];

    // Get primer positions - use explicit positions if available
    let fwdStart: number | null = null;
    let fwdEnd: number | null = null;
    let revStart: number | null = null;
    let revEnd: number | null = null;

    if (forward?.sequence) {
      if (typeof forward.start === 'number' && typeof forward.end === 'number') {
        fwdStart = forward.start;
        fwdEnd = forward.end;
      } else {
        // Fallback: search for primer binding
        const binding = findPrimerBinding(template, forward.sequence, false, {
          mutationPosition,
          isMutagenesis
        });
        if (binding) {
          fwdStart = binding.start;
          fwdEnd = binding.end;
        }
      }
    }

    if (reverse?.sequence) {
      if (typeof reverse.start === 'number' && typeof reverse.end === 'number') {
        revStart = reverse.start;
        revEnd = reverse.end;
      } else {
        // Fallback: search for primer binding
        const binding = findPrimerBinding(template, reverse.sequence, true, {
          mutationPosition,
          isMutagenesis
        });
        if (binding) {
          revStart = binding.start;
          revEnd = binding.end;
        }
      }
    }

    // Add forward primer annotation
    if (fwdStart !== null && fwdEnd !== null) {
      annots.push({
        name: `Forward (${forward!.sequence.length}bp, Tm ${forward!.tm || '?'}°C)`,
        start: fwdStart,
        end: fwdEnd,
        direction: 1,
        color: theme.forward,
        type: 'primer'
      });
    }

    // Add reverse primer annotation
    if (revStart !== null && revEnd !== null) {
      annots.push({
        name: `Reverse (${reverse!.sequence.length}bp, Tm ${reverse!.tm || '?'}°C)`,
        start: revStart,
        end: revEnd,
        direction: -1,
        color: theme.reverse,
        type: 'primer'
      });
    }

    // Add mutation marker as annotation
    if (mutationPosition !== null && mutationPosition !== undefined) {
      annots.push({
        name: `Mutation (pos ${mutationPosition + 1})`,
        start: mutationPosition,
        end: mutationPosition + (mutationLength || 1),
        direction: 1,
        color: theme.mutation,
        type: 'mutation'
      });
    }

    return annots;
  }, [template, forward, reverse, mutationPosition, mutationLength, isMutagenesis]);

  // Calculate primer region for info display
  const primerInfo = useMemo(() => {
    const fwdAnnot = annotations.find(a => a.type === 'primer' && a.direction === 1);
    const revAnnot = annotations.find(a => a.type === 'primer' && a.direction === -1);

    if (fwdAnnot && revAnnot) {
      const minStart = Math.min(fwdAnnot.start, revAnnot.start);
      const maxEnd = Math.max(fwdAnnot.end, revAnnot.end);
      return { size: maxEnd - minStart, fwd: fwdAnnot, rev: revAnnot };
    }
    return null;
  }, [annotations]);

  if (!template) {
    return (
      <div className="p-10 text-center rounded-xl" style={{
        color: theme.textMuted,
        background: theme.bgSecondary
      }}>
        No template sequence provided
      </div>
    );
  }

  return (
    <div className="rounded-xl overflow-hidden border" style={{
      background: theme.bg,
      borderColor: theme.border
    }}>
      {/* Header */}
      <div className="px-4 py-3 flex justify-between items-center border-b" style={{
        borderColor: theme.border,
        background: theme.bgSecondary
      }}>
        <span className="font-semibold text-sm" style={{ color: theme.text }}>
          Template ({template.length.toLocaleString()} bp)
        </span>
        {primerInfo && (
          <span className="text-xs font-medium" style={{ color: theme.amplicon }}>
            {primerInfo.size.toLocaleString()} bp {isMutagenesis ? 'region' : 'amplicon'}
          </span>
        )}
      </div>

      {/* SeqViz Viewer */}
      <div style={{ height: `${height}px` }}>
        <SeqViz
          name="Template"
          seq={template.toUpperCase()}
          annotations={annotations}
          viewer="linear"
          showComplement={false}
          showIndex={true}
          style={{ height: '100%', width: '100%' }}
          zoom={{ linear: 50 }}
          rotateOnScroll={false}
        />
      </div>

      {/* Primer Info Legend */}
      {primerInfo && (
        <div className="px-4 py-2.5 flex gap-5 text-xs border-t" style={{
          borderColor: theme.border,
          background: theme.bgSecondary
        }}>
          <div className="flex items-center gap-1.5">
            <span className="w-3 h-3 rounded" style={{
              background: theme.forward
            }}></span>
            <span style={{ color: theme.textSecondary }}>
              Forward: {primerInfo.fwd.start + 1}–{primerInfo.fwd.end}
            </span>
          </div>
          <div className="flex items-center gap-1.5">
            <span className="w-3 h-3 rounded" style={{
              background: theme.reverse
            }}></span>
            <span style={{ color: theme.textSecondary }}>
              Reverse: {primerInfo.rev.start + 1}–{primerInfo.rev.end}
            </span>
          </div>
          {mutationPosition !== null && (
            <div className="flex items-center gap-1.5">
              <span className="w-3 h-3 rounded" style={{
                background: theme.mutation
              }}></span>
              <span style={{ color: theme.textSecondary }}>
                Mutation: pos {mutationPosition + 1}
              </span>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

/**
 * Sequence detail view - Light Theme with proper double-stranded DNA visualization
 *
 * For forward primer: binds to top strand, extends 5'→3' (left to right)
 * For reverse primer: binds to bottom strand, extends 5'→3' (right to left on diagram)
 */
interface SequenceDetailViewProps {
  template: string;
  forward?: Primer | null;
  reverse?: Primer | null;
  mutationPosition?: number | null;
  isMutagenesis?: boolean;
}

interface RegionData {
  sequence: string;
  primerStart: number;
  primerEnd: number;
  absoluteStart: number;
  bindingOffset: number;
}

function SequenceDetailView({
  template,
  forward,
  reverse,
  mutationPosition = null,
  isMutagenesis = false
}: SequenceDetailViewProps): React.ReactElement {
  // For mutagenesis primers, use position data if available (forward primer contains mutation
  // and won't be found by sequence search on original template)
  const fwdBinding = useMemo((): BindingResult | null => {
    if (!forward?.sequence) return null;
    if (forward.start !== undefined && forward.end !== undefined) {
      return { start: forward.start, end: forward.end, matchLength: forward.sequence.length, score: 1.0, method: 'explicit' };
    }
    return findPrimerBinding(template, forward.sequence, false, {
      mutationPosition,
      isMutagenesis
    });
  }, [forward, template, mutationPosition, isMutagenesis]);

  const revBinding = useMemo((): BindingResult | null => {
    if (!reverse?.sequence) return null;
    if (reverse.start !== undefined && reverse.end !== undefined) {
      return { start: reverse.start, end: reverse.end, matchLength: reverse.sequence.length, score: 1.0, method: 'explicit' };
    }
    return findPrimerBinding(template, reverse.sequence, true, {
      mutationPosition,
      isMutagenesis
    });
  }, [reverse, template, mutationPosition, isMutagenesis]);

  const getRegionSequence = (
    binding: BindingResult | null,
    primer: Primer | null | undefined,
    contextBefore: number = 10,
    contextAfter: number = 10
  ): RegionData | null => {
    if (!binding) return null;

    const start = Math.max(0, binding.start - contextBefore);
    const end = Math.min(template.length, binding.end + contextAfter);
    const seq = template.slice(start, end);

    return {
      sequence: seq,
      primerStart: binding.start - start,
      primerEnd: binding.end - start,
      absoluteStart: start,
      bindingOffset: binding.offset || 0
    };
  };

  const fwdRegion = getRegionSequence(fwdBinding, forward, 15, 5);
  const revRegion = getRegionSequence(revBinding, reverse, 5, 15);

  // Shared styles for sequence display container
  const gridContainerStyle: React.CSSProperties = {
    background: theme.bgTertiary
  };

  const labelStyle: React.CSSProperties = {
    color: theme.textMuted
  };

  const directionLabelStyle: React.CSSProperties = {
    color: theme.textMuted
  };

  // Template strand color - distinct from primer colors
  const templateColor = '#059669'; // Emerald green for template strands

  // Render forward primer binding to single-stranded template
  // Forward primer (5'→3') binds complementarily to the bottom strand (3'→5')
  // The primer sequence is identical to the top strand, so we show bottom strand as template
  const renderForwardAlignment = (region: RegionData, primerSeq: string, primerColor: string): React.ReactElement => {
    const totalLength = region.sequence.length;
    const cellWidth = 10;

    // Bottom strand (3'→5') is the complement of top strand - this is what forward primer binds to
    const comp: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
    const templateStrand = region.sequence.split('').map(c => comp[c] || c).join('');

    const rowStyle: React.CSSProperties = {
      display: 'grid',
      gridTemplateColumns: `repeat(${totalLength}, ${cellWidth}px)`,
      gap: '0'
    };

    const cellStyle: React.CSSProperties = {
      height: '18px',
      lineHeight: '18px'
    };

    return (
      <div className="flex flex-col gap-0.5">
        {/* Template strand (3'→5') - what the forward primer binds to */}
        <div style={rowStyle}>
          {templateStrand.split('').map((base, i) => {
            const inPrimer = i >= region.primerStart && i < region.primerEnd;
            return (
              <div key={`t-${i}`} className="text-center rounded-sm" style={{
                ...cellStyle,
                color: inPrimer ? templateColor : theme.textMuted,
                fontWeight: inPrimer ? '600' : 'normal',
                background: inPrimer ? `${templateColor}15` : 'transparent'
              }}>
                {base}
              </div>
            );
          })}
        </div>

        {/* Base pairing indicators */}
        <div style={rowStyle}>
          {Array(totalLength).fill(null).map((_, i) => {
            const inPrimerRange = i >= region.primerStart && i < region.primerEnd;
            return (
              <div key={`bp-${i}`} className="text-center text-[10px] h-3.5 leading-[14px] opacity-60" style={{
                color: inPrimerRange ? theme.textMuted : 'transparent'
              }}>
                {inPrimerRange ? '│' : ''}
              </div>
            );
          })}
        </div>

        {/* Primer row (5'→3') - binds complementarily */}
        <div style={rowStyle}>
          {Array(totalLength).fill(null).map((_, i) => {
            const primerIndex = i - region.primerStart;
            const inPrimerRange = primerIndex >= 0 && primerIndex < primerSeq.length;
            return (
              <div key={`p-${i}`} className="text-center font-semibold rounded-sm" style={{
                ...cellStyle,
                color: primerColor,
                background: inPrimerRange ? `${primerColor}15` : 'transparent'
              }}>
                {inPrimerRange ? primerSeq[primerIndex] : ''}
              </div>
            );
          })}
        </div>

        {/* Direction arrows - pointing right (extension direction) */}
        <div style={rowStyle}>
          {Array(totalLength).fill(null).map((_, i) => {
            const primerIndex = i - region.primerStart;
            const inPrimerRange = primerIndex >= 0 && primerIndex < primerSeq.length;
            return (
              <div key={`a-${i}`} className="text-center text-[9px] h-3.5 leading-[14px] opacity-70" style={{
                color: primerColor
              }}>
                {inPrimerRange ? '▸' : ''}
              </div>
            );
          })}
        </div>
      </div>
    );
  };

  // Render reverse primer binding to single-stranded template
  // Reverse primer (5'→3') binds complementarily to the top strand (5'→3')
  // The primer is antiparallel, so its 5' end aligns with template 3' end
  const renderReverseAlignment = (region: RegionData, primerSeq: string, primerColor: string): React.ReactElement => {
    const totalLength = region.sequence.length;
    const cellWidth = 10;

    // Top strand (5'→3') is the template - this is what reverse primer binds to
    const templateStrand = region.sequence;

    const rowStyle: React.CSSProperties = {
      display: 'grid',
      gridTemplateColumns: `repeat(${totalLength}, ${cellWidth}px)`,
      gap: '0'
    };

    const cellStyle: React.CSSProperties = {
      height: '18px',
      lineHeight: '18px'
    };

    return (
      <div className="flex flex-col gap-0.5">
        {/* Template strand (5'→3') - what the reverse primer binds to */}
        <div style={rowStyle}>
          {templateStrand.split('').map((base, i) => {
            const inPrimer = i >= region.primerStart && i < region.primerEnd;
            return (
              <div key={`t-${i}`} className="text-center rounded-sm" style={{
                ...cellStyle,
                color: inPrimer ? templateColor : theme.textMuted,
                fontWeight: inPrimer ? '600' : 'normal',
                background: inPrimer ? `${templateColor}15` : 'transparent'
              }}>
                {base}
              </div>
            );
          })}
        </div>

        {/* Base pairing indicators */}
        <div style={rowStyle}>
          {Array(totalLength).fill(null).map((_, i) => {
            const inPrimerRange = i >= region.primerStart && i < region.primerEnd;
            return (
              <div key={`bp-${i}`} className="text-center text-[10px] h-3.5 leading-[14px] opacity-60" style={{
                color: inPrimerRange ? theme.textMuted : 'transparent'
              }}>
                {inPrimerRange ? '│' : ''}
              </div>
            );
          })}
        </div>

        {/* Primer row (5'→3' right to left) - binds complementarily */}
        {/* Primer is antiparallel: 5' on right, 3' on left */}
        <div style={rowStyle}>
          {Array(totalLength).fill(null).map((_, i) => {
            const primerIndex = i - region.primerStart;
            const inPrimerRange = primerIndex >= 0 && primerIndex < primerSeq.length;
            // Reverse the primer display so 5' end is on the right (antiparallel)
            const displayIndex = inPrimerRange ? (primerSeq.length - 1 - primerIndex) : -1;
            return (
              <div key={`p-${i}`} className="text-center font-semibold rounded-sm" style={{
                ...cellStyle,
                color: primerColor,
                background: inPrimerRange ? `${primerColor}15` : 'transparent'
              }}>
                {inPrimerRange ? primerSeq[displayIndex] : ''}
              </div>
            );
          })}
        </div>

        {/* Direction arrows - pointing left (extension direction) */}
        <div style={rowStyle}>
          {Array(totalLength).fill(null).map((_, i) => {
            const primerIndex = i - region.primerStart;
            const inPrimerRange = primerIndex >= 0 && primerIndex < primerSeq.length;
            return (
              <div key={`a-${i}`} className="text-center text-[9px] h-3.5 leading-[14px] opacity-70" style={{
                color: primerColor
              }}>
                {inPrimerRange ? '◂' : ''}
              </div>
            );
          })}
        </div>
      </div>
    );
  };

  return (
    <div className="grid gap-4 mt-4" style={{
      gridTemplateColumns: 'repeat(auto-fit, minmax(380px, 1fr))'
    }}>
      {/* Forward Primer Detail */}
      {fwdRegion && forward && (
        <div className="rounded-lg p-4 border shadow-sm" style={{
          background: theme.bg,
          borderColor: theme.border
        }}>
          <div className="flex justify-between items-center mb-3.5">
            <span className="font-semibold text-sm" style={{ color: theme.forward }}>
              Forward Primer Binding
            </span>
            <span className="text-xs font-mono px-2 py-0.5 rounded" style={{
              color: theme.textMuted,
              background: theme.bgTertiary
            }}>
              pos {fwdBinding!.start + 1}-{fwdBinding!.end}
            </span>
          </div>

          <div className="font-mono text-xs leading-relaxed overflow-x-auto p-3 rounded-lg shadow-inner" style={gridContainerStyle}>
            <div className="flex items-start gap-2">
              <div className="flex flex-col gap-0.5 min-w-[85px]">
                <span className="pr-2 whitespace-nowrap font-medium h-[18px] leading-[18px]" style={labelStyle}>Template:</span>
                <span className="h-3.5"></span>
                <span className="pr-2 whitespace-nowrap font-medium h-[18px] leading-[18px]" style={labelStyle}>Primer:</span>
                <span className="h-3.5"></span>
              </div>
              <div className="flex flex-col gap-0.5">
                <span className="pr-1 whitespace-nowrap h-[18px] leading-[18px]" style={directionLabelStyle}>3'-</span>
                <span className="h-3.5"></span>
                <span className="pr-1 whitespace-nowrap h-[18px] leading-[18px]" style={directionLabelStyle}>5'-</span>
                <span className="h-3.5"></span>
              </div>
              <div className="flex-1">
                {renderForwardAlignment(fwdRegion, forward.sequence, theme.forward)}
              </div>
              <div className="flex flex-col gap-0.5">
                <span className="pl-1 whitespace-nowrap h-[18px] leading-[18px]" style={directionLabelStyle}>-5'</span>
                <span className="h-3.5"></span>
                <span className="pl-1 whitespace-nowrap h-[18px] leading-[18px]" style={directionLabelStyle}>-3'</span>
                <span className="h-3.5"></span>
              </div>
            </div>
          </div>

          <div className="flex gap-4 mt-3.5 pt-3.5 border-t" style={{
            borderColor: theme.border
          }}>
            <span className="text-xs" style={{ color: theme.textMuted }}>
              Tm: <strong style={{ color: theme.text }}>{forward.tm}°C</strong>
            </span>
            <span className="text-xs" style={{ color: theme.textMuted }}>
              GC: <strong style={{ color: theme.text }}>{forward.gcPercent || `${((forward.gc || 0) * 100).toFixed(0)}%`}</strong>
            </span>
            <span className="text-xs" style={{ color: theme.textMuted }}>
              dG: <strong style={{ color: (forward.dg || 0) < -3 ? theme.mutation : theme.text }}>{forward.dg?.toFixed(1) || 'N/A'} kcal/mol</strong>
            </span>
          </div>
        </div>
      )}

      {/* Reverse Primer Detail - single-stranded view */}
      {revRegion && reverse && (
        <div className="rounded-lg p-4 border shadow-sm" style={{
          background: theme.bg,
          borderColor: theme.border
        }}>
          <div className="flex justify-between items-center mb-3.5">
            <span className="font-semibold text-sm" style={{ color: theme.reverse }}>
              Reverse Primer Binding
            </span>
            <span className="text-xs font-mono px-2 py-0.5 rounded" style={{
              color: theme.textMuted,
              background: theme.bgTertiary
            }}>
              pos {revBinding!.start + 1}-{revBinding!.end}
            </span>
          </div>

          <div className="font-mono text-xs leading-relaxed overflow-x-auto p-3 rounded-lg shadow-inner" style={gridContainerStyle}>
            <div className="flex items-start gap-2">
              <div className="flex flex-col gap-0.5 min-w-[85px]">
                <span className="pr-2 whitespace-nowrap font-medium h-[18px] leading-[18px]" style={labelStyle}>Template:</span>
                <span className="h-3.5"></span>
                <span className="pr-2 whitespace-nowrap font-medium h-[18px] leading-[18px]" style={labelStyle}>Primer:</span>
                <span className="h-3.5"></span>
              </div>
              <div className="flex flex-col gap-0.5">
                <span className="pr-1 whitespace-nowrap h-[18px] leading-[18px]" style={directionLabelStyle}>5'-</span>
                <span className="h-3.5"></span>
                <span className="pr-1 whitespace-nowrap h-[18px] leading-[18px]" style={directionLabelStyle}>3'-</span>
                <span className="h-3.5"></span>
              </div>
              <div className="flex-1">
                {renderReverseAlignment(revRegion, reverse.sequence, theme.reverse)}
              </div>
              <div className="flex flex-col gap-0.5">
                <span className="pl-1 whitespace-nowrap h-[18px] leading-[18px]" style={directionLabelStyle}>-3'</span>
                <span className="h-3.5"></span>
                <span className="pl-1 whitespace-nowrap h-[18px] leading-[18px]" style={directionLabelStyle}>-5'</span>
                <span className="h-3.5"></span>
              </div>
            </div>
          </div>

          <div className="flex gap-4 mt-3.5 pt-3.5 border-t" style={{
            borderColor: theme.border
          }}>
            <span className="text-xs" style={{ color: theme.textMuted }}>
              Tm: <strong style={{ color: theme.text }}>{reverse.tm}°C</strong>
            </span>
            <span className="text-xs" style={{ color: theme.textMuted }}>
              GC: <strong style={{ color: theme.text }}>{reverse.gcPercent || `${((reverse.gc || 0) * 100).toFixed(0)}%`}</strong>
            </span>
            <span className="text-xs" style={{ color: theme.textMuted }}>
              dG: <strong style={{ color: (reverse.dg || 0) < -3 ? theme.mutation : theme.text }}>{reverse.dg?.toFixed(1) || 'N/A'} kcal/mol</strong>
            </span>
          </div>
        </div>
      )}
    </div>
  );
}

/**
 * Main PrimerOnTemplateViewer Component - Light Theme
 */
interface PrimerOnTemplateViewerProps {
  template: string;
  forward?: Primer | null;
  reverse?: Primer | null;
  mutationPosition?: number | null;
  mutationLength?: number;
  showHairpinDiagrams?: boolean;
  showSequenceDetails?: boolean;
  width?: number;
  compact?: boolean;
  isSDMMode?: boolean;
}

export default function PrimerOnTemplateViewer({
  template,
  forward,
  reverse,
  mutationPosition = null,
  mutationLength = 1,
  showHairpinDiagrams = true,
  showSequenceDetails = true,
  width = 800,
  compact = false,
  isSDMMode = false
}: PrimerOnTemplateViewerProps): React.ReactElement {
  const [activeTab, setActiveTab] = useState<'overview' | 'structure'>('overview');
  const [expandedPrimer, setExpandedPrimer] = useState<string | null>(null);

  const fwdFold = useMemo((): FoldResult | null => {
    if (!forward?.sequence) return null;
    return foldSequence(forward.sequence);
  }, [forward?.sequence]);

  const revFold = useMemo((): FoldResult | null => {
    if (!reverse?.sequence) return null;
    return foldSequence(reverse.sequence);
  }, [reverse?.sequence]);

  const hasProblematicStructure = (fwdFold?.e ?? 0) < -3 || (revFold?.e ?? 0) < -3;

  // Use severity classification for each primer
  const fwdSeverity = useMemo((): SeverityResult => {
    if (!fwdFold || !forward?.sequence) return { level: 'none', shouldWarn: false };
    return classify3PrimeStructureSeverity({
      energy: fwdFold.e ?? 0,
      basePairs: fwdFold.ij ?? [],
      seqLength: forward.sequence.length,
    } as any) as any;
  }, [fwdFold, forward?.sequence]);

  const revSeverity = useMemo((): SeverityResult => {
    if (!revFold || !reverse?.sequence) return { level: 'none', shouldWarn: false };
    return classify3PrimeStructureSeverity({
      energy: revFold.e ?? 0,
      basePairs: revFold.ij ?? [],
      seqLength: reverse.sequence.length,
    } as any) as any;
  }, [revFold, reverse?.sequence]);

  // Get the worst severity between the two primers
  const worstSeverity = useMemo((): SeverityResult => {
    const severityOrder = ['none', 'info', 'low', 'moderate', 'warning', 'critical'];
    const fwdIdx = severityOrder.indexOf(fwdSeverity.level);
    const revIdx = severityOrder.indexOf(revSeverity.level);
    return fwdIdx > revIdx ? fwdSeverity : revSeverity;
  }, [fwdSeverity, revSeverity]);

  // has3PrimeIssue now means severity is warning or critical
  const has3PrimeIssue = worstSeverity.level === 'critical' || worstSeverity.level === 'warning';

  if (!template || (!forward && !reverse)) {
    return (
      <div className="p-5 text-center" style={{ color: theme.textMuted }}>
        No primers to display
      </div>
    );
  }

  return (
    <div className="rounded-lg overflow-hidden shadow-sm" style={{
      background: theme.bg,
      border: has3PrimeIssue ? `2px solid ${theme.danger}` : `1px solid ${theme.border}`
    }}>
      {/* Header */}
      <div className="px-4 py-3.5 flex justify-between items-center border-b" style={{
        background: theme.bgSecondary,
        borderColor: theme.border
      }}>
        <div className="flex items-center gap-3">
          <span className="font-semibold text-base" style={{ color: theme.text }}>
            Primer Visualization
          </span>

          {/* Severity-based structure badge */}
          {(hasProblematicStructure || worstSeverity.shouldWarn) && (
            <span
              className="px-3 py-1 rounded-md text-xs font-semibold cursor-help"
              style={{
                background: worstSeverity.level === 'critical' ? theme.dangerBg
                  : worstSeverity.level === 'warning' ? '#fef3c7'
                  : worstSeverity.level === 'moderate' ? theme.warningBg
                  : theme.infoBg,
                color: worstSeverity.level === 'critical' ? theme.danger
                  : worstSeverity.level === 'warning' ? '#d97706'
                  : worstSeverity.level === 'moderate' ? theme.warning
                  : theme.info
              }}
              title={worstSeverity.tooltip}
            >
              {worstSeverity.level === 'critical' ? '🔴 3\' End Blocked!'
                : worstSeverity.level === 'warning' ? '⚠ 3\' Structure Risk'
                : worstSeverity.level === 'moderate' ? '△ 3\' Structure'
                : '⚡ Secondary Structure'}
            </span>
          )}
        </div>

        {/* Tab buttons */}
        <div className="flex gap-1 p-0.5 rounded-lg" style={{ background: theme.bgTertiary }}>
          <button
            onClick={() => setActiveTab('overview')}
            className="px-3.5 py-1.5 border-0 rounded-md text-xs cursor-pointer transition-all duration-150"
            style={{
              background: activeTab === 'overview' ? theme.bg : 'transparent',
              color: activeTab === 'overview' ? theme.text : theme.textMuted,
              fontWeight: activeTab === 'overview' ? '600' : '500',
              boxShadow: activeTab === 'overview' ? '0 1px 2px rgba(0,0,0,0.08)' : 'none'
            }}
          >
            Overview
          </button>
          <button
            onClick={() => setActiveTab('structure')}
            className="px-3.5 py-1.5 border-0 rounded-md text-xs cursor-pointer transition-all duration-150 relative"
            style={{
              background: activeTab === 'structure' ? theme.bg : 'transparent',
              color: activeTab === 'structure' ? theme.text : theme.textMuted,
              fontWeight: activeTab === 'structure' ? '600' : '500',
              boxShadow: activeTab === 'structure' ? '0 1px 2px rgba(0,0,0,0.08)' : 'none'
            }}
          >
            Structure
            {hasProblematicStructure && (
              <span className="absolute top-1 right-1 w-1.5 h-1.5 rounded-full" style={{
                background: theme.danger
              }} />
            )}
          </button>
        </div>
      </div>

      {/* Content */}
      <div className="p-4">
        {activeTab === 'overview' && (
          <>
            <TemplateVisualization
              template={template}
              forward={forward}
              reverse={reverse}
              mutationPosition={mutationPosition}
              mutationLength={mutationLength}
              width={width - 36}
              compact={compact}
              isMutagenesis={isSDMMode || mutationPosition !== null}
            />

            {showSequenceDetails && !compact && (
              <SequenceDetailView
                template={template}
                forward={forward}
                reverse={reverse}
                mutationPosition={mutationPosition}
                isMutagenesis={isSDMMode || mutationPosition !== null}
              />
            )}

            {/* Compact structure summary - detailed values in Analysis section */}
            <div className="flex gap-4 mt-4 px-3.5 py-3 rounded-lg border items-center justify-between" style={{
              background: theme.bgSecondary,
              borderColor: theme.border
            }}>
              <div className="flex gap-5 items-center">
                {forward?.sequence && fwdFold && (
                  <div className="flex items-center gap-1.5">
                    <span className="w-2 h-2 rounded-full" style={{
                      background: fwdFold.e > -3 ? '#22c55e' : fwdFold.e > -6 ? '#f59e0b' : '#ef4444'
                    }}></span>
                    <span className="text-xs font-medium" style={{ color: theme.forward }}>
                      Fwd: {fwdFold.e.toFixed(1)} kcal/mol
                    </span>
                  </div>
                )}
                {reverse?.sequence && revFold && (
                  <div className="flex items-center gap-1.5">
                    <span className="w-2 h-2 rounded-full" style={{
                      background: revFold.e > -3 ? '#22c55e' : revFold.e > -6 ? '#f59e0b' : '#ef4444'
                    }}></span>
                    <span className="text-xs font-medium" style={{ color: theme.reverse }}>
                      Rev: {revFold.e.toFixed(1)} kcal/mol
                    </span>
                  </div>
                )}
              </div>
              <span className="text-xs" style={{ color: theme.textMuted }}>
                See Analysis section for detailed breakdown
              </span>
            </div>
          </>
        )}

        {activeTab === 'structure' && showHairpinDiagrams && (
          <div className="flex flex-col gap-5">
            {/* Introduction */}
            <div className="p-3.5 rounded-lg text-xs leading-relaxed border" style={{
              background: theme.infoBg,
              color: theme.info,
              borderColor: `${theme.info}30`
            }}>
              <strong>"Glass Box" Physics Visualization:</strong> Unlike simple pass/fail scores,
              these diagrams show the <em>actual physics</em> of your primers - how they fold, where they're stable,
              and whether they'll bind to each other.
            </div>

            <PrimerStructureViewer
              forwardSeq={forward?.sequence}
              reverseSeq={reverse?.sequence}
              forwardName="Forward Primer"
              reverseName="Reverse Primer"
              showCrossDimer={true}
              showHeatmap={true}
              showSparkline={false}
              isSDMMode={isSDMMode}
            />

            {/* Legacy diagrams */}
            <details className="mt-2">
              <summary className="cursor-pointer text-xs px-3 py-2.5 rounded-lg border font-medium" style={{
                color: theme.textSecondary,
                background: theme.bgSecondary,
                borderColor: theme.border
              }}>
                View Classic Stem-Loop Diagrams
              </summary>
              <div className="mt-3 flex flex-col gap-4">
                {forward?.sequence && (
                  <HairpinDiagram
                    sequence={forward.sequence}
                    foldResult={fwdFold as any}
                    primerName="Forward Primer"
                    width={Math.min(500, width - 36)}
                    showDetails={true}
                  />
                )}
                {reverse?.sequence && (
                  <HairpinDiagram
                    sequence={reverse.sequence}
                    foldResult={revFold as any}
                    primerName="Reverse Primer"
                    width={Math.min(500, width - 36)}
                    showDetails={true}
                  />
                )}
              </div>
            </details>

            {/* Severity-based Recommendations */}
            {worstSeverity.shouldWarn && (
              <div className="p-4 rounded-lg border" style={{
                background: worstSeverity.level === 'critical' ? theme.dangerBg : '#fef3c7',
                borderColor: worstSeverity.level === 'critical' ? theme.danger : '#f59e0b'
              }}>
                {worstSeverity.level === 'critical' ? (
                  <>
                    <div className="font-semibold text-sm mb-2.5 flex items-center gap-2" style={{
                      color: theme.danger
                    }}>
                      🔴 Redesign Required
                      <span className="px-2 py-0.5 rounded text-[10px] font-semibold uppercase text-white" style={{
                        background: theme.danger
                      }}>
                        Critical
                      </span>
                    </div>
                    <ul className="text-xs m-0 pl-5 leading-relaxed" style={{ color: '#991b1b' }}>
                      <li><strong>The 3' end is trapped in a stable structure</strong></li>
                      <li>DNA polymerase requires a free 3' OH to begin extension</li>
                      <li>PCR will likely <strong>fail</strong> or have very low yield</li>
                    </ul>
                    <div className="mt-3 px-3 py-2.5 rounded-md text-xs font-medium" style={{
                      background: 'rgba(220, 38, 38, 0.1)',
                      color: theme.danger
                    }}>
                      ⚠ <strong>Action Required:</strong> Adjust primer length or position to free the 3' end
                    </div>
                  </>
                ) : (
                  <>
                    <div className="font-semibold text-sm mb-2.5 flex items-center gap-2" style={{
                      color: '#d97706'
                    }}>
                      ⚠ Optimization Recommended
                      <span className="px-2 py-0.5 rounded text-[10px] font-semibold uppercase text-white" style={{
                        background: '#f59e0b'
                      }}>
                        {worstSeverity.level}
                      </span>
                    </div>
                    <ul className="text-xs m-0 pl-5 leading-relaxed" style={{ color: '#92400e' }}>
                      <li>Structure near 3' end may compete with template binding</li>
                      <li>PCR may work but with <strong>reduced efficiency</strong></li>
                      <li>Consider optimization if you experience low yields</li>
                    </ul>
                    <div className="mt-3 px-3 py-2.5 rounded-md text-xs" style={{
                      background: 'rgba(245, 158, 11, 0.1)',
                      color: '#92400e'
                    }}>
                      <strong>Tips:</strong> Adjust primer length by 1-2 bases, use touchdown PCR, or add 2-5% DMSO to reduce secondary structure
                    </div>
                  </>
                )}
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}

/**
 * Compact inline badge - Light Theme
 */
interface PrimerBindingBadgeProps {
  template: string;
  primer?: Primer | null;
  direction?: 'forward' | 'reverse';
  mutationPosition?: number | null;
}

export function PrimerBindingBadge({
  template,
  primer,
  direction = 'forward',
  mutationPosition = null
}: PrimerBindingBadgeProps): React.ReactElement {
  const isReverse = direction === 'reverse';
  const binding = primer?.sequence ? findPrimerBinding(template, primer.sequence, isReverse, {
    mutationPosition,
    isMutagenesis: mutationPosition !== null
  }) : null;
  const color = isReverse ? theme.reverse : theme.forward;

  if (!binding) {
    return (
      <span className="text-xs" style={{ color: theme.textMuted }}>
        Binding not found
      </span>
    );
  }

  return (
    <span className="px-2.5 py-0.5 rounded-md text-xs font-mono font-medium" style={{
      background: `${color}15`,
      color: color
    }}>
      {binding.start + 1}-{binding.end} ({binding.matchLength}bp)
    </span>
  );
}
