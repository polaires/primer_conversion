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

// Light theme colors
const theme = {
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

// Reverse complement helper
function reverseComplement(seq) {
  const comp = { A: 'T', T: 'A', G: 'C', C: 'G', a: 't', t: 'a', g: 'c', c: 'g' };
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
function findPrimerBinding(template, primerSeq, isReverse = false, options = {}) {
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

    let estimatedStart;
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
    const positions5 = [];
    let idx = tmpl.indexOf(anchor5);
    while (idx !== -1) {
      positions5.push(idx);
      idx = tmpl.indexOf(anchor5, idx + 1);
    }

    // Find 3' anchor positions
    const positions3 = [];
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
  let bestMatch = null;
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
function PrimerLabel({
  x,
  y,
  label,
  sublabel,
  color,
  anchor = 'middle',
  position = 'above' // 'above' or 'below'
}) {
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
function PrimerBlock({
  x,
  width,
  y,
  color,
  direction, // 'forward' (→) or 'reverse' (←)
  onHover,
  onLeave
}) {
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
}) {
  const height = compact ? 200 : 250;

  // Build annotations for SeqViz
  const annotations = useMemo(() => {
    const annots = [];

    // Get primer positions - use explicit positions if available
    let fwdStart = null, fwdEnd = null;
    let revStart = null, revEnd = null;

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
        name: `Forward (${forward.sequence.length}bp, Tm ${forward.tm || '?'}°C)`,
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
        name: `Reverse (${reverse.sequence.length}bp, Tm ${reverse.tm || '?'}°C)`,
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
      <div style={{
        padding: '40px',
        textAlign: 'center',
        color: theme.textMuted,
        background: theme.bgSecondary,
        borderRadius: '12px'
      }}>
        No template sequence provided
      </div>
    );
  }

  return (
    <div style={{
      background: theme.bg,
      borderRadius: '12px',
      border: `1px solid ${theme.border}`,
      overflow: 'hidden'
    }}>
      {/* Header */}
      <div style={{
        padding: '12px 16px',
        borderBottom: `1px solid ${theme.border}`,
        background: theme.bgSecondary,
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center'
      }}>
        <span style={{ fontWeight: '600', color: theme.text, fontSize: '13px' }}>
          Template ({template.length.toLocaleString()} bp)
        </span>
        {primerInfo && (
          <span style={{ color: theme.amplicon, fontSize: '12px', fontWeight: '500' }}>
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
        <div style={{
          padding: '10px 16px',
          borderTop: `1px solid ${theme.border}`,
          background: theme.bgSecondary,
          display: 'flex',
          gap: '20px',
          fontSize: '11px'
        }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
            <span style={{
              width: '12px',
              height: '12px',
              borderRadius: '3px',
              background: theme.forward
            }}></span>
            <span style={{ color: theme.textSecondary }}>
              Forward: {primerInfo.fwd.start + 1}–{primerInfo.fwd.end}
            </span>
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
            <span style={{
              width: '12px',
              height: '12px',
              borderRadius: '3px',
              background: theme.reverse
            }}></span>
            <span style={{ color: theme.textSecondary }}>
              Reverse: {primerInfo.rev.start + 1}–{primerInfo.rev.end}
            </span>
          </div>
          {mutationPosition !== null && (
            <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
              <span style={{
                width: '12px',
                height: '12px',
                borderRadius: '3px',
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
function SequenceDetailView({ template, forward, reverse, mutationPosition = null, isMutagenesis = false }) {
  // For mutagenesis primers, use position data if available (forward primer contains mutation
  // and won't be found by sequence search on original template)
  const fwdBinding = useMemo(() => {
    if (!forward?.sequence) return null;
    if (forward.start !== undefined && forward.end !== undefined) {
      return { start: forward.start, end: forward.end, matchLength: forward.sequence.length };
    }
    return findPrimerBinding(template, forward.sequence, false, {
      mutationPosition,
      isMutagenesis
    });
  }, [forward, template, mutationPosition, isMutagenesis]);

  const revBinding = useMemo(() => {
    if (!reverse?.sequence) return null;
    if (reverse.start !== undefined && reverse.end !== undefined) {
      return { start: reverse.start, end: reverse.end, matchLength: reverse.sequence.length };
    }
    return findPrimerBinding(template, reverse.sequence, true, {
      mutationPosition,
      isMutagenesis
    });
  }, [reverse, template, mutationPosition, isMutagenesis]);

  const getRegionSequence = (binding, primer, contextBefore = 10, contextAfter = 10) => {
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
  const gridContainerStyle = {
    fontFamily: 'monospace',
    fontSize: '12px',
    lineHeight: '1.6',
    overflowX: 'auto',
    background: theme.bgTertiary,
    padding: '12px',
    borderRadius: '8px',
    boxShadow: 'inset 0 1px 2px rgba(0,0,0,0.05)'
  };

  const labelStyle = {
    color: theme.textMuted,
    paddingRight: '8px',
    whiteSpace: 'nowrap',
    fontWeight: '500'
  };

  const directionLabelStyle = {
    color: theme.textMuted,
    paddingRight: '4px',
    whiteSpace: 'nowrap'
  };

  // Template strand color - distinct from primer colors
  const templateColor = '#059669'; // Emerald green for template strands

  // Render forward primer binding to single-stranded template
  // Forward primer (5'→3') binds complementarily to the bottom strand (3'→5')
  // The primer sequence is identical to the top strand, so we show bottom strand as template
  const renderForwardAlignment = (region, primerSeq, primerColor) => {
    const totalLength = region.sequence.length;
    const cellWidth = 10;

    // Bottom strand (3'→5') is the complement of top strand - this is what forward primer binds to
    const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
    const templateStrand = region.sequence.split('').map(c => comp[c] || c).join('');

    const rowStyle = {
      display: 'grid',
      gridTemplateColumns: `repeat(${totalLength}, ${cellWidth}px)`,
      gap: '0'
    };

    const cellStyle = {
      textAlign: 'center',
      height: '18px',
      lineHeight: '18px'
    };

    return (
      <div style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
        {/* Template strand (3'→5') - what the forward primer binds to */}
        <div style={rowStyle}>
          {templateStrand.split('').map((base, i) => {
            const inPrimer = i >= region.primerStart && i < region.primerEnd;
            return (
              <div key={`t-${i}`} style={{
                ...cellStyle,
                color: inPrimer ? templateColor : theme.textMuted,
                fontWeight: inPrimer ? '600' : 'normal',
                background: inPrimer ? `${templateColor}15` : 'transparent',
                borderRadius: '2px'
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
              <div key={`bp-${i}`} style={{
                ...cellStyle,
                color: inPrimerRange ? theme.textMuted : 'transparent',
                fontSize: '10px',
                height: '14px',
                lineHeight: '14px',
                opacity: 0.6
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
              <div key={`p-${i}`} style={{
                ...cellStyle,
                color: primerColor,
                fontWeight: '600',
                background: inPrimerRange ? `${primerColor}15` : 'transparent',
                borderRadius: '2px'
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
              <div key={`a-${i}`} style={{
                ...cellStyle,
                color: primerColor,
                fontSize: '9px',
                height: '14px',
                lineHeight: '14px',
                opacity: 0.7
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
  const renderReverseAlignment = (region, primerSeq, primerColor) => {
    const totalLength = region.sequence.length;
    const cellWidth = 10;

    // Top strand (5'→3') is the template - this is what reverse primer binds to
    const templateStrand = region.sequence;

    const rowStyle = {
      display: 'grid',
      gridTemplateColumns: `repeat(${totalLength}, ${cellWidth}px)`,
      gap: '0'
    };

    const cellStyle = {
      textAlign: 'center',
      height: '18px',
      lineHeight: '18px'
    };

    return (
      <div style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
        {/* Template strand (5'→3') - what the reverse primer binds to */}
        <div style={rowStyle}>
          {templateStrand.split('').map((base, i) => {
            const inPrimer = i >= region.primerStart && i < region.primerEnd;
            return (
              <div key={`t-${i}`} style={{
                ...cellStyle,
                color: inPrimer ? templateColor : theme.textMuted,
                fontWeight: inPrimer ? '600' : 'normal',
                background: inPrimer ? `${templateColor}15` : 'transparent',
                borderRadius: '2px'
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
              <div key={`bp-${i}`} style={{
                ...cellStyle,
                color: inPrimerRange ? theme.textMuted : 'transparent',
                fontSize: '10px',
                height: '14px',
                lineHeight: '14px',
                opacity: 0.6
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
              <div key={`p-${i}`} style={{
                ...cellStyle,
                color: primerColor,
                fontWeight: '600',
                background: inPrimerRange ? `${primerColor}15` : 'transparent',
                borderRadius: '2px'
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
              <div key={`a-${i}`} style={{
                ...cellStyle,
                color: primerColor,
                fontSize: '9px',
                height: '14px',
                lineHeight: '14px',
                opacity: 0.7
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
    <div style={{
      display: 'grid',
      gridTemplateColumns: 'repeat(auto-fit, minmax(380px, 1fr))',
      gap: '16px',
      marginTop: '16px'
    }}>
      {/* Forward Primer Detail */}
      {fwdRegion && (
        <div style={{
          background: theme.bg,
          borderRadius: '10px',
          padding: '16px',
          border: `1px solid ${theme.border}`,
          boxShadow: '0 2px 4px rgba(0,0,0,0.04)'
        }}>
          <div style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            marginBottom: '14px'
          }}>
            <span style={{ color: theme.forward, fontWeight: '600', fontSize: '13px' }}>
              Forward Primer Binding
            </span>
            <span style={{
              color: theme.textMuted,
              fontSize: '11px',
              fontFamily: 'monospace',
              background: theme.bgTertiary,
              padding: '2px 8px',
              borderRadius: '4px'
            }}>
              pos {fwdBinding.start + 1}-{fwdBinding.end}
            </span>
          </div>

          <div style={gridContainerStyle}>
            <div style={{ display: 'flex', alignItems: 'flex-start', gap: '8px' }}>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '2px', minWidth: '85px' }}>
                <span style={{ ...labelStyle, height: '18px', lineHeight: '18px' }}>Template:</span>
                <span style={{ height: '14px' }}></span>
                <span style={{ ...labelStyle, height: '18px', lineHeight: '18px' }}>Primer:</span>
                <span style={{ height: '14px' }}></span>
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                <span style={{ ...directionLabelStyle, height: '18px', lineHeight: '18px' }}>3'-</span>
                <span style={{ height: '14px' }}></span>
                <span style={{ ...directionLabelStyle, height: '18px', lineHeight: '18px' }}>5'-</span>
                <span style={{ height: '14px' }}></span>
              </div>
              <div style={{ flex: 1 }}>
                {renderForwardAlignment(fwdRegion, forward.sequence, theme.forward)}
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                <span style={{ ...directionLabelStyle, height: '18px', lineHeight: '18px', paddingLeft: '4px' }}>-5'</span>
                <span style={{ height: '14px' }}></span>
                <span style={{ ...directionLabelStyle, height: '18px', lineHeight: '18px', paddingLeft: '4px' }}>-3'</span>
                <span style={{ height: '14px' }}></span>
              </div>
            </div>
          </div>

          <div style={{
            display: 'flex',
            gap: '16px',
            marginTop: '14px',
            paddingTop: '14px',
            borderTop: `1px solid ${theme.border}`
          }}>
            <span style={{ color: theme.textMuted, fontSize: '11px' }}>
              Tm: <strong style={{ color: theme.text }}>{forward.tm}°C</strong>
            </span>
            <span style={{ color: theme.textMuted, fontSize: '11px' }}>
              GC: <strong style={{ color: theme.text }}>{forward.gcPercent || `${(forward.gc * 100).toFixed(0)}%`}</strong>
            </span>
            <span style={{ color: theme.textMuted, fontSize: '11px' }}>
              dG: <strong style={{ color: forward.dg < -3 ? theme.mutation : theme.text }}>{forward.dg?.toFixed(1) || 'N/A'} kcal/mol</strong>
            </span>
          </div>
        </div>
      )}

      {/* Reverse Primer Detail - single-stranded view */}
      {revRegion && (
        <div style={{
          background: theme.bg,
          borderRadius: '10px',
          padding: '16px',
          border: `1px solid ${theme.border}`,
          boxShadow: '0 2px 4px rgba(0,0,0,0.04)'
        }}>
          <div style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            marginBottom: '14px'
          }}>
            <span style={{ color: theme.reverse, fontWeight: '600', fontSize: '13px' }}>
              Reverse Primer Binding
            </span>
            <span style={{
              color: theme.textMuted,
              fontSize: '11px',
              fontFamily: 'monospace',
              background: theme.bgTertiary,
              padding: '2px 8px',
              borderRadius: '4px'
            }}>
              pos {revBinding.start + 1}-{revBinding.end}
            </span>
          </div>

          <div style={gridContainerStyle}>
            <div style={{ display: 'flex', alignItems: 'flex-start', gap: '8px' }}>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '2px', minWidth: '85px' }}>
                <span style={{ ...labelStyle, height: '18px', lineHeight: '18px' }}>Template:</span>
                <span style={{ height: '14px' }}></span>
                <span style={{ ...labelStyle, height: '18px', lineHeight: '18px' }}>Primer:</span>
                <span style={{ height: '14px' }}></span>
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                <span style={{ ...directionLabelStyle, height: '18px', lineHeight: '18px' }}>5'-</span>
                <span style={{ height: '14px' }}></span>
                <span style={{ ...directionLabelStyle, height: '18px', lineHeight: '18px' }}>3'-</span>
                <span style={{ height: '14px' }}></span>
              </div>
              <div style={{ flex: 1 }}>
                {renderReverseAlignment(revRegion, reverse.sequence, theme.reverse)}
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                <span style={{ ...directionLabelStyle, height: '18px', lineHeight: '18px', paddingLeft: '4px' }}>-3'</span>
                <span style={{ height: '14px' }}></span>
                <span style={{ ...directionLabelStyle, height: '18px', lineHeight: '18px', paddingLeft: '4px' }}>-5'</span>
                <span style={{ height: '14px' }}></span>
              </div>
            </div>
          </div>

          <div style={{
            display: 'flex',
            gap: '16px',
            marginTop: '14px',
            paddingTop: '14px',
            borderTop: `1px solid ${theme.border}`
          }}>
            <span style={{ color: theme.textMuted, fontSize: '11px' }}>
              Tm: <strong style={{ color: theme.text }}>{reverse.tm}°C</strong>
            </span>
            <span style={{ color: theme.textMuted, fontSize: '11px' }}>
              GC: <strong style={{ color: theme.text }}>{reverse.gcPercent || `${(reverse.gc * 100).toFixed(0)}%`}</strong>
            </span>
            <span style={{ color: theme.textMuted, fontSize: '11px' }}>
              dG: <strong style={{ color: reverse.dg < -3 ? theme.mutation : theme.text }}>{reverse.dg?.toFixed(1) || 'N/A'} kcal/mol</strong>
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
}) {
  const [activeTab, setActiveTab] = useState('overview');
  const [expandedPrimer, setExpandedPrimer] = useState(null);

  const fwdFold = useMemo(() => {
    if (!forward?.sequence) return null;
    return foldSequence(forward.sequence);
  }, [forward?.sequence]);

  const revFold = useMemo(() => {
    if (!reverse?.sequence) return null;
    return foldSequence(reverse.sequence);
  }, [reverse?.sequence]);

  const hasProblematicStructure = (fwdFold?.e < -3) || (revFold?.e < -3);

  // Use severity classification for each primer
  const fwdSeverity = useMemo(() => {
    if (!fwdFold || !forward?.sequence) return { level: 'none', shouldWarn: false };
    return classify3PrimeStructureSeverity({
      energy: fwdFold.e ?? 0,
      basePairs: fwdFold.ij ?? [],
      seqLength: forward.sequence.length,
    });
  }, [fwdFold, forward?.sequence]);

  const revSeverity = useMemo(() => {
    if (!revFold || !reverse?.sequence) return { level: 'none', shouldWarn: false };
    return classify3PrimeStructureSeverity({
      energy: revFold.e ?? 0,
      basePairs: revFold.ij ?? [],
      seqLength: reverse.sequence.length,
    });
  }, [revFold, reverse?.sequence]);

  // Get the worst severity between the two primers
  const worstSeverity = useMemo(() => {
    const severityOrder = ['none', 'info', 'low', 'moderate', 'warning', 'critical'];
    const fwdIdx = severityOrder.indexOf(fwdSeverity.level);
    const revIdx = severityOrder.indexOf(revSeverity.level);
    return fwdIdx > revIdx ? fwdSeverity : revSeverity;
  }, [fwdSeverity, revSeverity]);

  // has3PrimeIssue now means severity is warning or critical
  const has3PrimeIssue = worstSeverity.level === 'critical' || worstSeverity.level === 'warning';

  if (!template || (!forward && !reverse)) {
    return (
      <div style={{ color: theme.textMuted, padding: '20px', textAlign: 'center' }}>
        No primers to display
      </div>
    );
  }

  return (
    <div style={{
      background: theme.bg,
      borderRadius: '10px',
      border: has3PrimeIssue ? `2px solid ${theme.danger}` : `1px solid ${theme.border}`,
      overflow: 'hidden',
      boxShadow: '0 1px 3px rgba(0,0,0,0.08)'
    }}>
      {/* Header */}
      <div style={{
        padding: '14px 18px',
        background: theme.bgSecondary,
        borderBottom: `1px solid ${theme.border}`,
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center'
      }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
          <span style={{ color: theme.text, fontWeight: '600', fontSize: '15px' }}>
            Primer Visualization
          </span>

          {/* Severity-based structure badge */}
          {(hasProblematicStructure || worstSeverity.shouldWarn) && (
            <span
              style={{
                padding: '4px 12px',
                background: worstSeverity.level === 'critical' ? theme.dangerBg
                  : worstSeverity.level === 'warning' ? '#fef3c7'
                  : worstSeverity.level === 'moderate' ? theme.warningBg
                  : theme.infoBg,
                color: worstSeverity.level === 'critical' ? theme.danger
                  : worstSeverity.level === 'warning' ? '#d97706'
                  : worstSeverity.level === 'moderate' ? theme.warning
                  : theme.info,
                borderRadius: '6px',
                fontSize: '11px',
                fontWeight: '600',
                cursor: 'help',
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
        <div style={{ display: 'flex', gap: '4px', background: theme.bgTertiary, padding: '3px', borderRadius: '8px' }}>
          <button
            onClick={() => setActiveTab('overview')}
            style={{
              padding: '6px 14px',
              background: activeTab === 'overview' ? theme.bg : 'transparent',
              border: 'none',
              borderRadius: '6px',
              color: activeTab === 'overview' ? theme.text : theme.textMuted,
              fontSize: '12px',
              fontWeight: activeTab === 'overview' ? '600' : '500',
              cursor: 'pointer',
              boxShadow: activeTab === 'overview' ? '0 1px 2px rgba(0,0,0,0.08)' : 'none',
              transition: 'all 0.15s ease'
            }}
          >
            Overview
          </button>
          <button
            onClick={() => setActiveTab('structure')}
            style={{
              padding: '6px 14px',
              background: activeTab === 'structure' ? theme.bg : 'transparent',
              border: 'none',
              borderRadius: '6px',
              color: activeTab === 'structure' ? theme.text : theme.textMuted,
              fontSize: '12px',
              fontWeight: activeTab === 'structure' ? '600' : '500',
              cursor: 'pointer',
              boxShadow: activeTab === 'structure' ? '0 1px 2px rgba(0,0,0,0.08)' : 'none',
              transition: 'all 0.15s ease',
              position: 'relative'
            }}
          >
            Structure
            {hasProblematicStructure && (
              <span style={{
                position: 'absolute',
                top: '4px',
                right: '4px',
                width: '6px',
                height: '6px',
                borderRadius: '50%',
                background: theme.danger
              }} />
            )}
          </button>
        </div>
      </div>

      {/* Content */}
      <div style={{ padding: '18px' }}>
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
            <div style={{
              display: 'flex',
              gap: '16px',
              marginTop: '16px',
              padding: '12px 14px',
              background: theme.bgSecondary,
              borderRadius: '8px',
              border: `1px solid ${theme.border}`,
              alignItems: 'center',
              justifyContent: 'space-between'
            }}>
              <div style={{ display: 'flex', gap: '20px', alignItems: 'center' }}>
                {forward?.sequence && fwdFold && (
                  <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
                    <span style={{
                      width: '8px',
                      height: '8px',
                      borderRadius: '50%',
                      background: fwdFold.e > -3 ? '#22c55e' : fwdFold.e > -6 ? '#f59e0b' : '#ef4444'
                    }}></span>
                    <span style={{ color: theme.forward, fontSize: '12px', fontWeight: '500' }}>
                      Fwd: {fwdFold.e.toFixed(1)} kcal/mol
                    </span>
                  </div>
                )}
                {reverse?.sequence && revFold && (
                  <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
                    <span style={{
                      width: '8px',
                      height: '8px',
                      borderRadius: '50%',
                      background: revFold.e > -3 ? '#22c55e' : revFold.e > -6 ? '#f59e0b' : '#ef4444'
                    }}></span>
                    <span style={{ color: theme.reverse, fontSize: '12px', fontWeight: '500' }}>
                      Rev: {revFold.e.toFixed(1)} kcal/mol
                    </span>
                  </div>
                )}
              </div>
              <span style={{ color: theme.textMuted, fontSize: '11px' }}>
                See Analysis section for detailed breakdown
              </span>
            </div>
          </>
        )}

        {activeTab === 'structure' && showHairpinDiagrams && (
          <div style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
            {/* Introduction */}
            <div style={{
              padding: '14px',
              background: theme.infoBg,
              borderRadius: '8px',
              color: theme.info,
              fontSize: '12px',
              lineHeight: '1.6',
              border: `1px solid ${theme.info}30`
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
            <details style={{ marginTop: '8px' }}>
              <summary style={{
                color: theme.textSecondary,
                cursor: 'pointer',
                fontSize: '12px',
                padding: '10px 12px',
                background: theme.bgSecondary,
                borderRadius: '8px',
                border: `1px solid ${theme.border}`,
                fontWeight: '500'
              }}>
                View Classic Stem-Loop Diagrams
              </summary>
              <div style={{ marginTop: '12px', display: 'flex', flexDirection: 'column', gap: '16px' }}>
                {forward?.sequence && (
                  <HairpinDiagram
                    sequence={forward.sequence}
                    foldResult={fwdFold}
                    primerName="Forward Primer"
                    width={Math.min(500, width - 36)}
                    showDetails={true}
                  />
                )}
                {reverse?.sequence && (
                  <HairpinDiagram
                    sequence={reverse.sequence}
                    foldResult={revFold}
                    primerName="Reverse Primer"
                    width={Math.min(500, width - 36)}
                    showDetails={true}
                  />
                )}
              </div>
            </details>

            {/* Severity-based Recommendations */}
            {worstSeverity.shouldWarn && (
              <div style={{
                padding: '16px',
                background: worstSeverity.level === 'critical' ? theme.dangerBg : '#fef3c7',
                borderRadius: '8px',
                border: `1px solid ${worstSeverity.level === 'critical' ? theme.danger : '#f59e0b'}`
              }}>
                {worstSeverity.level === 'critical' ? (
                  <>
                    <div style={{
                      color: theme.danger,
                      fontWeight: '600',
                      fontSize: '14px',
                      marginBottom: '10px',
                      display: 'flex',
                      alignItems: 'center',
                      gap: '8px'
                    }}>
                      🔴 Redesign Required
                      <span style={{
                        padding: '2px 8px',
                        background: theme.danger,
                        color: '#fff',
                        borderRadius: '4px',
                        fontSize: '10px',
                        fontWeight: '600',
                        textTransform: 'uppercase',
                      }}>
                        Critical
                      </span>
                    </div>
                    <ul style={{ color: '#991b1b', fontSize: '12px', margin: 0, paddingLeft: '20px', lineHeight: '1.7' }}>
                      <li><strong>The 3' end is trapped in a stable structure</strong></li>
                      <li>DNA polymerase requires a free 3' OH to begin extension</li>
                      <li>PCR will likely <strong>fail</strong> or have very low yield</li>
                    </ul>
                    <div style={{
                      marginTop: '12px',
                      padding: '10px 12px',
                      background: 'rgba(220, 38, 38, 0.1)',
                      borderRadius: '6px',
                      fontSize: '12px',
                      color: theme.danger,
                      fontWeight: '500'
                    }}>
                      ⚠ <strong>Action Required:</strong> Adjust primer length or position to free the 3' end
                    </div>
                  </>
                ) : (
                  <>
                    <div style={{
                      color: '#d97706',
                      fontWeight: '600',
                      fontSize: '14px',
                      marginBottom: '10px',
                      display: 'flex',
                      alignItems: 'center',
                      gap: '8px'
                    }}>
                      ⚠ Optimization Recommended
                      <span style={{
                        padding: '2px 8px',
                        background: '#f59e0b',
                        color: '#fff',
                        borderRadius: '4px',
                        fontSize: '10px',
                        fontWeight: '600',
                        textTransform: 'uppercase',
                      }}>
                        {worstSeverity.level}
                      </span>
                    </div>
                    <ul style={{ color: '#92400e', fontSize: '12px', margin: 0, paddingLeft: '20px', lineHeight: '1.7' }}>
                      <li>Structure near 3' end may compete with template binding</li>
                      <li>PCR may work but with <strong>reduced efficiency</strong></li>
                      <li>Consider optimization if you experience low yields</li>
                    </ul>
                    <div style={{
                      marginTop: '12px',
                      padding: '10px 12px',
                      background: 'rgba(245, 158, 11, 0.1)',
                      borderRadius: '6px',
                      fontSize: '12px',
                      color: '#92400e',
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
export function PrimerBindingBadge({ template, primer, direction = 'forward', mutationPosition = null }) {
  const isReverse = direction === 'reverse';
  const binding = primer?.sequence ? findPrimerBinding(template, primer.sequence, isReverse, {
    mutationPosition,
    isMutagenesis: mutationPosition !== null
  }) : null;
  const color = isReverse ? theme.reverse : theme.forward;

  if (!binding) {
    return (
      <span style={{ color: theme.textMuted, fontSize: '11px' }}>
        Binding not found
      </span>
    );
  }

  return (
    <span style={{
      padding: '3px 10px',
      background: `${color}15`,
      color: color,
      borderRadius: '6px',
      fontSize: '11px',
      fontFamily: 'monospace',
      fontWeight: '500'
    }}>
      {binding.start + 1}-{binding.end} ({binding.matchLength}bp)
    </span>
  );
}
