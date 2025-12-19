import React, { useMemo, useState, useEffect, useRef, useCallback } from 'react';
import { foldSequence, fold, dg } from '../lib/fold.js';
import FornaViewer from './FornaViewer.jsx';

/**
 * PrimerStructureViewer - "Glass Box" visualization of primer physics
 *
 * Light theme version with merged charts and interactive features
 */

// ============================================================================
// TYPE DEFINITIONS
// ============================================================================

interface BasePair extends Array<number> {
  0: number;
  1: number;
}

interface FoldResult {
  e: number;
  ij: BasePair[];
}

interface DimerAlignment {
  offset: number;
  score: number;
  pairs: Array<{
    i: number;
    j: number;
    base1: string;
    base2: string;
  }>;
  maxConsecutive: number;
  involves3Prime1: boolean;
  involves3Prime2: boolean;
  involves5Prime1: boolean;
  involves5Prime2: boolean;
}

interface Node {
  id: number;
  base: string;
  x: number;
  y: number;
  vx: number;
  vy: number;
  pairingProb: number;
}

interface ForceDirectedStructureProps {
  sequence: string;
  basePairs?: BasePair[];
  width?: number;
  height?: number;
  label?: string;
  onShift?: ((offset: number) => void) | null;
}

interface CrossDimerZipperProps {
  forwardSeq: string;
  reverseSeq: string;
  width?: number;
  height?: number;
  isSDMMode?: boolean;
}

interface ThermodynamicHeatmapProps {
  sequence: string;
  templateRegion?: any;
  width?: number;
  label?: string;
}

interface PrimerStructureViewerProps {
  forwardSeq?: string;
  reverseSeq?: string;
  forwardName?: string;
  reverseName?: string;
  templateRegion?: any;
  showCrossDimer?: boolean;
  showHeatmap?: boolean;
  showSparkline?: boolean;
  isSDMMode?: boolean;
  onForwardShift?: ((offset: number) => void) | null;
  onReverseShift?: ((offset: number) => void) | null;
  useFornaViewer?: boolean;
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

const complement: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };

function reverseComplement(seq: string): string {
  return seq.split('').reverse().map(c => complement[c] || c).join('');
}

function basePairEnergy(base1: string, base2: string): number {
  const pair = `${base1}${base2}`;
  const energies: Record<string, number> = {
    'GC': -2.4, 'CG': -2.4,
    'AT': -1.5, 'TA': -1.5,
    'GT': -0.5, 'TG': -0.5,
    'GA': 0, 'AG': 0, 'CA': 0, 'AC': 0, 'CT': 0, 'TC': 0,
    'AA': 0.5, 'TT': 0.5, 'GG': 0.5, 'CC': 0.5,
  };
  return energies[pair] || 0;
}

function canPair(base1: string, base2: string): boolean {
  return complement[base1] === base2;
}

function findBestDimerAlignment(seq1: string, seq2: string): DimerAlignment | null {
  const rc2 = reverseComplement(seq2);
  let bestScore = 0;
  let bestAlignment: DimerAlignment | null = null;

  for (let offset = -(seq2.length - 3); offset < seq1.length - 2; offset++) {
    let score = 0;
    const pairs: Array<{ i: number; j: number; base1: string; base2: string }> = [];
    let consecutive = 0;
    let maxConsecutive = 0;

    for (let i = 0; i < seq1.length; i++) {
      const j = i - offset;
      if (j >= 0 && j < rc2.length) {
        if (canPair(seq1[i], seq2[seq2.length - 1 - j])) {
          const energy = Math.abs(basePairEnergy(seq1[i], seq2[seq2.length - 1 - j]));
          score += energy;
          pairs.push({ i, j: seq2.length - 1 - j, base1: seq1[i], base2: seq2[seq2.length - 1 - j] });
          consecutive++;
          maxConsecutive = Math.max(maxConsecutive, consecutive);
        } else {
          consecutive = 0;
        }
      }
    }

    const weightedScore = score + (maxConsecutive * 2);

    if (weightedScore > bestScore) {
      bestScore = weightedScore;
      bestAlignment = {
        offset,
        score: weightedScore,
        pairs,
        maxConsecutive,
        involves3Prime1: pairs.some(p => p.i >= seq1.length - 5),
        involves3Prime2: pairs.some(p => p.j >= seq2.length - 5),
        involves5Prime1: pairs.some(p => p.i < 5),
        involves5Prime2: pairs.some(p => p.j < 5),
      };
    }
  }

  return bestAlignment;
}

function calculatePositionwiseDG(sequence: string, windowSize: number = 10): number[] {
  const energies: number[] = [];
  const halfWindow = Math.floor(windowSize / 2);

  for (let i = 0; i < sequence.length; i++) {
    const start = Math.max(0, i - halfWindow);
    const end = Math.min(sequence.length, i + halfWindow + 1);
    const subseq = sequence.slice(start, end);
    let energy = subseq.length >= 6 ? dg(subseq, 37) : 0;
    // Clamp infinity values to reasonable bounds
    if (!isFinite(energy)) {
      energy = 0;
    }
    energies.push(energy);
  }

  return energies;
}

function calculatePairingProbability(sequence: string, index: number, basePairs: BasePair[]): number {
  const isPaired = basePairs.some(([i, j]) => i === index || j === index);
  if (isPaired) return 1.0;

  const base = sequence[index];
  const windowSize = 5;
  const start = Math.max(0, index - windowSize);
  const end = Math.min(sequence.length, index + windowSize + 1);
  const window = sequence.slice(start, end);

  const gcCount = window.split('').filter(b => b === 'G' || b === 'C').length;
  const gcRatio = gcCount / window.length;

  let pairingPotential = 0;
  for (let i = 0; i < window.length; i++) {
    if (complement[base] === window[i]) pairingPotential += 0.1;
  }

  return Math.min(1.0, gcRatio * 0.5 + pairingPotential);
}

// ============================================================================
// LIGHT THEME COLORS
// ============================================================================

const theme = {
  bg: '#ffffff',
  bgSecondary: '#f8fafc',
  bgTertiary: '#f1f5f9',
  border: '#e2e8f0',
  borderStrong: '#cbd5e1',
  text: '#1e293b',
  textSecondary: '#475569',
  textMuted: '#94a3b8',
  // Semantic colors
  success: '#16a34a',
  successBg: '#dcfce7',
  warning: '#ca8a04',
  warningBg: '#fef9c3',
  danger: '#dc2626',
  dangerBg: '#fee2e2',
  info: '#2563eb',
  infoBg: '#dbeafe',
  // Base colors
  gc: '#16a34a',
  at: '#ea580c',
  prime3: '#2563eb',
};

// ============================================================================
// FORCE-DIRECTED 2D STRUCTURE
// ============================================================================

function ForceDirectedStructure({
  sequence,
  basePairs = [],
  width = 300,
  height = 250,
  label = 'Primer',
  onShift = null
}: ForceDirectedStructureProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [nodes, setNodes] = useState<Node[]>([]);
  const [savedNodes, setSavedNodes] = useState<Node[] | null>(null);
  const [simpleView, setSimpleView] = useState<boolean>(false);
  const [hoveredNode, setHoveredNode] = useState<number | null>(null);

  // Create deterministic layout based on structure
  useEffect(() => {
    if (!sequence || sequence.length < 2) return;

    const n = sequence.length;
    const hasStructure = basePairs && basePairs.length > 0;
    const centerX = width / 2;
    const padding = 30;

    // Build pair map
    const pairMap = new Map<number, number>();
    if (basePairs) {
      basePairs.forEach(([i, j]) => {
        pairMap.set(i, j);
        pairMap.set(j, i);
      });
    }

    const nodeArray: Node[] = [];

    if (hasStructure && basePairs.length > 0) {
      // For hairpin: create proper vertical stem-loop layout
      // Sort pairs by first index to understand structure
      const sortedPairs = [...basePairs].sort((a, b) => a[0] - b[0]);

      // Find the innermost pair (smallest gap = closest to loop)
      const innerPair = sortedPairs.reduce((inner, pair) => {
        const gap = pair[1] - pair[0];
        const innerGap = inner[1] - inner[0];
        return gap < innerGap ? pair : inner;
      }, sortedPairs[0]);

      // Layout constants - generous spacing, viewBox will scale to fit
      const stemWidth = 80; // Horizontal distance between left and right stems
      const loopTopY = 50; // Y position where loop starts
      const stemPairCount = sortedPairs.length;

      // Fixed node spacing - viewBox handles scaling
      const nodeSpacing = 28;

      // Calculate stem bottom Y based on number of pairs
      const stemBottomY = loopTopY + 20 + stemPairCount * nodeSpacing;

      // Calculate positions for all bases
      // Strategy: Walk through sequence, place each base appropriately

      for (let i = 0; i < n; i++) {
        let x: number, y: number;
        const partner = pairMap.get(i);
        const isPaired = partner !== undefined;

        if (isPaired) {
          // This base is part of a pair
          const isLeftStem = i < partner!; // 5' side of pair

          // Find which pair this belongs to
          const pairIdx = sortedPairs.findIndex(([a, b]) => a === i || b === i);

          if (pairIdx !== -1) {
            // Position based on pair index
            // pairIdx 0 = outermost (bottom), higher = closer to loop (top)
            const stemLevel = stemPairCount - 1 - pairIdx;
            y = loopTopY + 20 + stemLevel * nodeSpacing;
            x = isLeftStem ? (centerX - stemWidth / 2) : (centerX + stemWidth / 2);
          } else {
            // Fallback
            x = centerX;
            y = height / 2;
          }
        } else {
          // Unpaired base - could be in loop, 5' tail, or 3' tail
          const loopStart = innerPair[0] + 1;
          const loopEnd = innerPair[1] - 1;
          const firstPairedIdx = sortedPairs[0][0];
          const lastPairedIdx = sortedPairs[sortedPairs.length - 1][1];

          if (i >= loopStart && i <= loopEnd) {
            // Loop base - arrange in arc at top
            const loopLen = loopEnd - loopStart + 1;
            const loopIdx = i - loopStart;
            // Arc from left stem top to right stem top
            // Calculate radius to ensure minimum 22px spacing between bases
            // Arc spans 0.6π radians (108°), so: spacing = radius * 0.6π / loopLen
            // For 22px spacing: radius = 22 * loopLen / (0.6π) ≈ 11.7 * loopLen
            const minRadius = 30;
            const spacingBasedRadius = (loopLen * 22) / (Math.PI * 0.6);
            const loopRadius = Math.max(minRadius, spacingBasedRadius);
            // Angle goes from ~150° to ~30° (top arc)
            const startAngle = Math.PI * 0.8;
            const endAngle = Math.PI * 0.2;
            const angle = startAngle - (startAngle - endAngle) * ((loopIdx + 0.5) / loopLen);

            x = centerX + Math.cos(angle) * loopRadius;
            // Use full radius for Y to prevent vertical compression
            y = loopTopY - 10 + Math.sin(angle) * (loopRadius * 0.6);
          } else if (i < firstPairedIdx) {
            // 5' tail - position along left side going down from bottom of stem
            const tailLen = firstPairedIdx;
            const tailIdx = i;
            const tailSpacing = 22; // Fixed spacing, viewBox scales
            x = centerX - stemWidth / 2 - 25;
            y = stemBottomY + 15 + (tailLen - 1 - tailIdx) * tailSpacing;
          } else if (i > lastPairedIdx) {
            // 3' tail - position along right side going down from bottom of stem
            const tailStart = lastPairedIdx + 1;
            const tailIdx = i - tailStart;
            const tailSpacing = 22; // Fixed spacing, viewBox scales
            x = centerX + stemWidth / 2 + 25;
            y = stemBottomY + 15 + tailIdx * tailSpacing;
          } else {
            // Unpaired base between pairs (bulge) - offset slightly
            // Find nearest paired base to position relative to
            let nearestPairY = height / 2;
            for (const [a, b] of sortedPairs) {
              if (i > a && i < b) {
                const pairIdx = sortedPairs.findIndex(([pa, pb]) => pa === a);
                const stemLevel = stemPairCount - 1 - pairIdx;
                nearestPairY = loopTopY + 20 + stemLevel * nodeSpacing;
                break;
              }
            }
            // Offset bulge bases slightly
            const isLeftSide = i < innerPair[0] + (innerPair[1] - innerPair[0]) / 2;
            x = isLeftSide ? (centerX - stemWidth / 2 - 12) : (centerX + stemWidth / 2 + 12);
            y = nearestPairY;
          }
        }

        nodeArray.push({
          id: i,
          base: sequence[i],
          x,
          y,
          vx: 0,
          vy: 0,
          pairingProb: calculatePairingProbability(sequence, i, basePairs || []),
        });
      }
    } else {
      // No structure: simple arc from 5' to 3'
      for (let i = 0; i < n; i++) {
        const progress = i / Math.max(1, n - 1);
        const x = padding + progress * (width - padding * 2);
        const y = height / 2 + Math.sin(progress * Math.PI) * 40;

        nodeArray.push({
          id: i,
          base: sequence[i],
          x,
          y,
          vx: 0,
          vy: 0,
          pairingProb: calculatePairingProbability(sequence, i, basePairs || []),
        });
      }
    }

    setNodes(nodeArray);
    setSavedNodes(nodeArray);
  }, [sequence, basePairs, width, height]);

  // Calculate viewBox based on actual node positions with padding
  const viewBox = useMemo(() => {
    if (nodes.length === 0) return `0 0 ${width} ${height}`;

    const padding = 25;
    const xs = nodes.map(n => n.x);
    const ys = nodes.map(n => n.y);
    const minX = Math.min(...xs) - padding;
    const maxX = Math.max(...xs) + padding;
    const minY = Math.min(...ys) - padding;
    const maxY = Math.max(...ys) + padding + 30; // Extra space for status bar

    const contentWidth = maxX - minX;
    const contentHeight = maxY - minY;

    return `${minX} ${minY} ${contentWidth} ${contentHeight}`;
  }, [nodes, width, height]);

  const energy = useMemo(() => {
    if (!sequence || sequence.length < 6) return 0;
    return dg(sequence, 37);
  }, [sequence]);

  const getProbabilityColor = (prob: number, index: number): string => {
    const is3Prime = index >= sequence.length - 5;
    if (is3Prime && prob > 0.5) return theme.danger;
    if (prob > 0.7) return '#ea580c';
    if (prob > 0.4) return theme.warning;
    if (is3Prime) return theme.prime3;
    return theme.success;
  };

  const energyColor = energy >= -1 ? theme.success : energy >= -3 ? theme.warning : theme.danger;

  // Toggle to simple view - restore saved positions when switching back
  const handleToggleView = () => {
    if (simpleView && savedNodes) {
      setNodes(savedNodes);
    }
    setSimpleView(!simpleView);
  };

  if (simpleView) {
    return (
      <div style={{ background: theme.bg, borderRadius: '8px', padding: '12px', border: `1px solid ${theme.border}`, overflow: 'hidden' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '8px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
            <span style={{ color: theme.text, fontWeight: 'bold', fontSize: '13px' }}>{label}</span>
            {onShift && (
              <div style={{ display: 'flex', gap: '2px' }}>
                <button onClick={() => onShift(-1)} style={{ background: theme.bgSecondary, border: `1px solid ${theme.border}`, borderRadius: '4px', padding: '2px 6px', cursor: 'pointer', fontSize: '12px' }}>◀</button>
                <button onClick={() => onShift(1)} style={{ background: theme.bgSecondary, border: `1px solid ${theme.border}`, borderRadius: '4px', padding: '2px 6px', cursor: 'pointer', fontSize: '12px' }}>▶</button>
              </div>
            )}
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
            <span style={{ color: energyColor, fontSize: '12px', fontFamily: 'monospace', fontWeight: 'bold' }}>
              ΔG = {energy.toFixed(1)} kcal/mol
            </span>
            <button onClick={handleToggleView} style={{ background: theme.bgSecondary, border: `1px solid ${theme.border}`, borderRadius: '4px', color: theme.text, padding: '4px 8px', fontSize: '10px', cursor: 'pointer' }}>
              ◉ Full View
            </button>
          </div>
        </div>

        <div style={{ background: theme.bgTertiary, borderRadius: '6px', padding: '12px', fontFamily: 'monospace' }}>
          <div style={{ display: 'flex', alignItems: 'center', marginBottom: '8px' }}>
            <span style={{ color: theme.prime3, marginRight: '8px', fontSize: '11px' }}>5'</span>
            {sequence.split('').map((base, i) => {
              const isPaired = basePairs.some(([a, b]) => a === i || b === i);
              const is3Prime = i >= sequence.length - 5;
              const prob = nodes[i]?.pairingProb || 0;

              return (
                <span key={i} style={{
                  display: 'inline-block', width: '16px', height: '24px', lineHeight: '24px',
                  textAlign: 'center', fontSize: '11px', fontWeight: isPaired ? 'bold' : 'normal',
                  color: getProbabilityColor(prob, i),
                  background: is3Prime ? 'rgba(37, 99, 235, 0.1)' : 'transparent',
                  borderRadius: '2px',
                  borderBottom: isPaired ? `2px solid ${getProbabilityColor(prob, i)}` : 'none',
                }} title={`${base} (pos ${i + 1}) - Risk: ${(prob * 100).toFixed(0)}%${is3Prime ? ' [3\']' : ''}`}>
                  {base}
                </span>
              );
            })}
            <span style={{ color: theme.danger, marginLeft: '8px', fontSize: '11px' }}>3'</span>
          </div>
          <div style={{ display: 'flex', gap: '12px', fontSize: '9px', color: theme.textMuted, marginTop: '8px' }}>
            <span><span style={{ color: theme.success }}>●</span> Safe</span>
            <span><span style={{ color: theme.warning }}>●</span> Moderate</span>
            <span><span style={{ color: theme.danger }}>●</span> Critical</span>
            <span style={{ marginLeft: 'auto' }}>{basePairs.length > 0 ? `${basePairs.length} bp` : 'No structure'}</span>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div style={{ background: theme.bg, borderRadius: '8px', padding: '12px', border: `1px solid ${theme.border}`, overflow: 'hidden' }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '8px' }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
          <span style={{ color: theme.text, fontWeight: 'bold', fontSize: '13px' }}>{label}</span>
          {onShift && (
            <div style={{ display: 'flex', gap: '2px' }}>
              <button onClick={() => onShift(-1)} style={{ background: theme.bgSecondary, border: `1px solid ${theme.border}`, borderRadius: '4px', padding: '2px 6px', cursor: 'pointer', fontSize: '12px' }} title="Shift primer 1bp upstream">◀</button>
              <button onClick={() => onShift(1)} style={{ background: theme.bgSecondary, border: `1px solid ${theme.border}`, borderRadius: '4px', padding: '2px 6px', cursor: 'pointer', fontSize: '12px' }} title="Shift primer 1bp downstream">▶</button>
            </div>
          )}
        </div>
        <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
          <span style={{ color: energyColor, fontSize: '12px', fontFamily: 'monospace', fontWeight: 'bold' }}>
            ΔG = {energy.toFixed(1)} kcal/mol
          </span>
          <button onClick={handleToggleView} style={{ background: theme.bgSecondary, border: `1px solid ${theme.border}`, borderRadius: '4px', color: theme.text, padding: '4px 8px', fontSize: '10px', cursor: 'pointer' }}>
            ○ Simple View
          </button>
        </div>
      </div>

      <svg
        ref={svgRef}
        width={width}
        height={height}
        viewBox={viewBox}
        preserveAspectRatio="xMidYMid meet"
        style={{
          background: theme.bgTertiary,
          borderRadius: '6px',
          overflow: 'hidden',
          display: 'block',
          maxWidth: '100%'
        }}
      >
        {/* Gradient definition for backbone direction */}
        <defs>
          <linearGradient id={`backbone-gradient-${label.replace(/\s/g, '')}`} x1="0%" y1="0%" x2="100%" y2="0%">
            <stop offset="0%" stopColor={theme.prime3} />
            <stop offset="100%" stopColor={theme.danger} />
          </linearGradient>
          <marker id={`arrow-${label.replace(/\s/g, '')}`} markerWidth="6" markerHeight="6" refX="3" refY="3" orient="auto">
            <path d="M0,0 L6,3 L0,6 Z" fill={theme.textMuted} />
          </marker>
        </defs>

        {/* Backbone lines with gradient coloring showing 5'→3' direction */}
        {nodes.length > 1 && nodes.map((node, i) => {
          if (i === 0) return null;
          const prev = nodes[i - 1];
          const progress = i / (nodes.length - 1);
          // Interpolate color from blue (5') to red (3')
          const r = Math.round(37 + progress * (220 - 37));
          const g = Math.round(99 - progress * 99);
          const b = Math.round(235 - progress * (235 - 38));
          const color = `rgb(${r}, ${g}, ${b})`;

          return (
            <line
              key={`backbone-${i}`}
              x1={prev.x}
              y1={prev.y}
              x2={node.x}
              y2={node.y}
              stroke={color}
              strokeWidth="3"
              strokeLinecap="round"
              opacity="0.7"
            />
          );
        })}

        {/* Base pair connections - dotted lines between paired bases */}
        {basePairs.map(([i, j], idx) => {
          if (!nodes[i] || !nodes[j]) return null;
          const involves3Prime = i >= sequence.length - 5 || j >= sequence.length - 5;
          const midX = (nodes[i].x + nodes[j].x) / 2;
          const midY = (nodes[i].y + nodes[j].y) / 2;

          return (
            <g key={`pair-${idx}`}>
              <line
                x1={nodes[i].x}
                y1={nodes[i].y}
                x2={nodes[j].x}
                y2={nodes[j].y}
                stroke={involves3Prime ? theme.danger : '#ea580c'}
                strokeWidth={involves3Prime ? 2.5 : 2}
                strokeDasharray={involves3Prime ? '0' : '4,3'}
                opacity={0.8}
              />
              {/* Small indicator for base pair */}
              <circle cx={midX} cy={midY} r="3" fill={involves3Prime ? theme.danger : '#ea580c'} opacity="0.6" />
            </g>
          );
        })}

        {/* Nucleotide nodes */}
        {nodes.map((node, i) => {
          const prob = node.pairingProb || 0;
          const color = getProbabilityColor(prob, i);
          const isHovered = hoveredNode === i;
          const is3Prime = i >= sequence.length - 5;
          const is5Prime = i < 3;
          const nodeRadius = isHovered ? 12 : 10;

          return (
            <g key={i} onMouseEnter={() => setHoveredNode(i)} onMouseLeave={() => setHoveredNode(null)} style={{ cursor: 'pointer' }}>
              {/* Outer ring for 3' region */}
              {is3Prime && (
                <circle cx={node.x} cy={node.y} r={nodeRadius + 3} fill="none" stroke={theme.danger} strokeWidth="1" strokeDasharray="2,2" opacity="0.5" />
              )}
              {/* Main node circle */}
              <circle
                cx={node.x}
                cy={node.y}
                r={nodeRadius}
                fill={theme.bg}
                stroke={color}
                strokeWidth={is3Prime ? 2.5 : 1.5}
                style={{ filter: isHovered ? 'drop-shadow(0 0 4px rgba(0,0,0,0.2))' : 'none' }}
              />
              {/* Base letter */}
              <text x={node.x} y={node.y + 4} fill={color} fontSize="10" fontFamily="monospace" textAnchor="middle" fontWeight="bold">
                {node.base}
              </text>
              {/* Tooltip on hover */}
              {isHovered && (
                <g>
                  <rect x={node.x - 55} y={node.y - 42} width={110} height={26} rx={4} fill={theme.text} />
                  <text x={node.x} y={node.y - 25} fill={theme.bg} fontSize="9" fontFamily="monospace" textAnchor="middle">
                    {`${node.base}${i + 1} | ${(prob * 100).toFixed(0)}% pairing risk`}
                  </text>
                </g>
              )}
            </g>
          );
        })}

        {/* 5' and 3' labels - positioned to the side of first/last nodes */}
        {nodes.length > 0 && (
          <>
            {/* 5' label - to the left of first node */}
            <g>
              <circle cx={nodes[0].x - 18} cy={nodes[0].y} r="12" fill={theme.prime3} opacity="0.15" />
              <text x={nodes[0].x - 18} y={nodes[0].y + 4} fill={theme.prime3} fontSize="10" fontFamily="monospace" textAnchor="middle" fontWeight="bold">5'</text>
            </g>
            {/* 3' label - to the right of last node */}
            <g>
              <circle cx={nodes[nodes.length - 1].x + 18} cy={nodes[nodes.length - 1].y} r="12" fill={theme.danger} opacity="0.15" />
              <text x={nodes[nodes.length - 1].x + 18} y={nodes[nodes.length - 1].y + 4} fill={theme.danger} fontSize="10" fontFamily="monospace" textAnchor="middle" fontWeight="bold">3'</text>
            </g>
            {/* Direction arrow indicator */}
            {nodes.length > 3 && (
              <g opacity="0.4">
                <text x={width / 2} y={18} fill={theme.textMuted} fontSize="9" fontFamily="monospace" textAnchor="middle">
                  5' → 3' direction
                </text>
              </g>
            )}
          </>
        )}

      </svg>

      {/* Status bar - outside SVG so it's not affected by viewBox */}
      <div style={{
        background: basePairs.length > 0 ? theme.warningBg : theme.successBg,
        borderRadius: '4px',
        padding: '6px 12px',
        marginTop: '8px',
        textAlign: 'center'
      }}>
        <span style={{ color: basePairs.length > 0 ? theme.warning : theme.success, fontSize: '10px', fontFamily: 'monospace' }}>
          {basePairs.length > 0 ? `${basePairs.length} base pair${basePairs.length > 1 ? 's' : ''} – Hairpin detected` : 'Linear structure – No secondary structure'}
        </span>
      </div>

      <div style={{ display: 'flex', gap: '12px', marginTop: '8px', fontSize: '9px', color: theme.textMuted, justifyContent: 'center' }}>
        <span><span style={{ color: theme.success }}>●</span> Safe</span>
        <span><span style={{ color: theme.warning }}>●</span> Moderate</span>
        <span><span style={{ color: '#ea580c' }}>●</span> High risk</span>
        <span><span style={{ color: theme.danger }}>●</span> 3' structure</span>
      </div>
    </div>
  );
}

// ============================================================================
// CROSS-DIMER "ZIPPER" VISUALIZATION
// ============================================================================

function CrossDimerZipper({
  forwardSeq,
  reverseSeq,
  width = 600,
  height = 200,
  isSDMMode = false
}: CrossDimerZipperProps) {
  const alignment = useMemo(() => {
    if (!forwardSeq || !reverseSeq) return null;
    return findBestDimerAlignment(forwardSeq, reverseSeq);
  }, [forwardSeq, reverseSeq]);

  if (!alignment) {
    return (
      <div style={{ background: theme.bg, borderRadius: '8px', padding: '16px', textAlign: 'center', color: theme.textMuted, border: `1px solid ${theme.border}` }}>
        Unable to analyze cross-dimer
      </div>
    );
  }

  const { pairs, maxConsecutive, involves3Prime1, involves3Prime2, involves5Prime1, involves5Prime2, offset } = alignment;

  const only5PrimeOverlap = (involves5Prime1 || involves5Prime2) && !involves3Prime1 && !involves3Prime2;
  const hasLigationJunction = isSDMMode && only5PrimeOverlap;

  // Determine severity with better classification
  // Key insight: 3' extensible dimers (where 3' end is hybridized) are CRITICAL
  // Internal/5' dimers are annoying but often PCR-viable
  let severity: 'critical' | 'warning' | 'safe';
  let dimerType = 'none';

  if (involves3Prime1 && involves3Prime2 && maxConsecutive >= 4) {
    severity = 'critical';
    dimerType = '3prime_extensible';
  } else if ((involves3Prime1 || involves3Prime2) && maxConsecutive >= 3) {
    severity = isSDMMode && !involves3Prime1 && !involves3Prime2 ? 'safe' : 'warning';
    dimerType = 'internal_with_3prime';
  } else if (only5PrimeOverlap && isSDMMode) {
    severity = 'safe';
    dimerType = 'ligation_junction';
  } else if (maxConsecutive >= 4) {
    severity = 'warning';
    dimerType = 'internal';
  } else if (pairs.length > 0) {
    severity = 'safe';
    dimerType = 'minor';
  } else {
    severity = 'safe';
    dimerType = 'none';
  }

  const severityStyles = {
    critical: { bg: theme.dangerBg, border: theme.danger, text: '#991b1b', badge: theme.danger, badgeBg: theme.dangerBg },
    warning: { bg: theme.warningBg, border: theme.warning, text: '#854d0e', badge: theme.warning, badgeBg: theme.warningBg },
    safe: { bg: theme.successBg, border: theme.success, text: '#166534', badge: theme.success, badgeBg: theme.successBg }
  };

  const colors = severityStyles[severity];

  const pairedPositions1 = new Set(pairs.map(p => p.i));
  const pairedPositions2 = new Set(pairs.map(p => p.j));

  return (
    <div style={{ background: theme.bg, borderRadius: '8px', padding: '16px', border: `1px solid ${colors.border}` }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '12px' }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
          <span style={{ fontSize: '18px' }}>
            {severity === 'critical' ? '⚠️' : severity === 'warning' ? '⚡' : '✓'}
          </span>
          <span style={{ color: theme.text, fontWeight: 'bold' }}>Cross-Dimer Analysis</span>
        </div>
        <span style={{
          padding: '4px 12px', borderRadius: '4px', background: colors.badgeBg, color: colors.badge,
          fontSize: '11px', fontWeight: 'bold', fontFamily: 'monospace'
        }}>
          {dimerType === '3prime_extensible' ? "3' EXTENSIBLE DIMER!" :
           dimerType === 'ligation_junction' ? '🔗 Ligation Junction (Intended)' :
           dimerType === 'internal_with_3prime' ? 'Internal + 3\' End' :
           dimerType === 'internal' ? 'Internal Binding' :
           dimerType === 'minor' ? 'Minor Binding' : 'No Issues'}
        </span>
      </div>

      {/* Zipper Visualization */}
      {(() => {
        // Calculate alignment padding
        // offset tells us: for forward position i, it pairs with reverse position (reverseSeq.length - 1 - (i - offset))
        // When reverse is displayed reversed, display position d shows reverseSeq[reverseSeq.length - 1 - d]
        // We need forward position i to align with reverse display position where pair.j is shown
        // For pair at (i, j): forward shows at i, reverse shows j at display position (reverseSeq.length - 1 - j)
        // Calculate padding needed so paired positions align
        const fwdPadding = Math.max(0, -offset);
        const revPadding = Math.max(0, offset);

        // Total display length
        const totalLen = Math.max(forwardSeq.length + fwdPadding, reverseSeq.length + revPadding);

        return (
          <div style={{ background: theme.bgTertiary, borderRadius: '6px', padding: '16px', fontFamily: 'monospace', fontSize: '12px', overflowX: 'auto' }}>
            {/* Forward primer */}
            <div style={{ display: 'flex', alignItems: 'center', marginBottom: '4px', whiteSpace: 'nowrap' }}>
              <span style={{ color: theme.info, width: '70px', flexShrink: 0 }}>Forward:</span>
              <span style={{ color: theme.info, width: '50px', flexShrink: 0, textAlign: 'right', marginRight: '4px' }}>5' —</span>
              {/* Padding spacers */}
              {Array(fwdPadding).fill(null).map((_, i) => (
                <span key={`fpad-${i}`} style={{ display: 'inline-block', width: '10px', textAlign: 'center', color: 'transparent' }}>·</span>
              ))}
              {forwardSeq.split('').map((base, i) => {
                const isPaired = pairedPositions1.has(i);
                const is3Prime = i >= forwardSeq.length - 5;
                const is5Prime = i < 5;
                const isFirst3Prime = i === forwardSeq.length - 5;
                const isLast = i === forwardSeq.length - 1;

                let color = theme.textMuted;
                let bg = 'transparent';

                if (isPaired && is3Prime) {
                  color = theme.danger;
                  bg = 'rgba(220, 38, 38, 0.15)';
                } else if (isPaired && is5Prime && hasLigationJunction) {
                  color = theme.success;
                  bg = 'rgba(22, 163, 74, 0.15)';
                } else if (isPaired) {
                  color = theme.warning;
                  bg = 'rgba(202, 138, 4, 0.1)';
                } else if (is3Prime) {
                  color = theme.prime3;
                }

                const boxStyle = is3Prime && involves3Prime1 ? {
                  borderTop: `2px solid ${theme.danger}`,
                  borderBottom: `2px solid ${theme.danger}`,
                  borderLeft: isFirst3Prime ? `2px solid ${theme.danger}` : 'none',
                  borderRight: isLast ? `2px solid ${theme.danger}` : 'none',
                  borderRadius: isFirst3Prime ? '3px 0 0 3px' : isLast ? '0 3px 3px 0' : '0',
                  marginTop: '-2px', marginBottom: '-2px', paddingTop: '2px', paddingBottom: '2px',
                } : is5Prime && hasLigationJunction ? {
                  borderTop: `2px solid ${theme.success}`,
                  borderBottom: `2px solid ${theme.success}`,
                  borderLeft: i === 0 ? `2px solid ${theme.success}` : 'none',
                  borderRight: i === 4 ? `2px solid ${theme.success}` : 'none',
                  borderRadius: i === 0 ? '3px 0 0 3px' : i === 4 ? '0 3px 3px 0' : '0',
                  marginTop: '-2px', marginBottom: '-2px', paddingTop: '2px', paddingBottom: '2px',
                } : {};

                return (
                  <span key={i} style={{ color, background: bg, padding: '0 1px', fontWeight: isPaired ? 'bold' : 'normal', display: 'inline-block', width: '10px', textAlign: 'center', ...boxStyle }}>
                    {base}
                  </span>
                );
              })}
              <span style={{ color: theme.info, marginLeft: '4px' }}>—▶ 3'</span>
            </div>

            {/* Binding lines */}
            <div style={{ display: 'flex', alignItems: 'center', marginBottom: '4px', whiteSpace: 'nowrap' }}>
              <span style={{ width: '70px', flexShrink: 0 }}></span>
              <span style={{ width: '50px', flexShrink: 0, marginRight: '4px' }}></span>
              {/* Padding spacers */}
              {Array(fwdPadding).fill(null).map((_, i) => (
                <span key={`bpad-${i}`} style={{ display: 'inline-block', width: '10px', textAlign: 'center' }}>&nbsp;</span>
              ))}
              {forwardSeq.split('').map((_, i) => {
                const pair = pairs.find(p => p.i === i);
                if (!pair) return <span key={i} style={{ display: 'inline-block', width: '10px', textAlign: 'center' }}>&nbsp;</span>;

                const is3Prime = i >= forwardSeq.length - 5 || pair.j >= reverseSeq.length - 5;
                const is5Prime = i < 5 && pair.j >= reverseSeq.length - 5;
                const isStrong = canPair(forwardSeq[i], reverseSeq[pair.j]) && (forwardSeq[i] === 'G' || forwardSeq[i] === 'C');

                let color = theme.warning;
                if (is3Prime) color = theme.danger;
                else if (is5Prime && hasLigationJunction) color = theme.success;

                return (
                  <span key={i} style={{ color, fontWeight: 'bold', display: 'inline-block', width: '10px', textAlign: 'center' }}>
                    {isStrong ? '║' : '│'}
                  </span>
                );
              })}
            </div>

            {/* Reverse primer - displayed 3' to 5' but aligned by offset */}
            <div style={{ display: 'flex', alignItems: 'center', whiteSpace: 'nowrap' }}>
              <span style={{ color: '#7c3aed', width: '70px', flexShrink: 0 }}>Reverse:</span>
              <span style={{ color: '#7c3aed', width: '50px', flexShrink: 0, textAlign: 'right', marginRight: '4px' }}>3' ◀—</span>
              {/* Padding spacers for reverse */}
              {Array(revPadding).fill(null).map((_, i) => (
                <span key={`rpad-${i}`} style={{ display: 'inline-block', width: '10px', textAlign: 'center', color: 'transparent' }}>·</span>
              ))}
              {reverseSeq.split('').reverse().map((base, i) => {
                const actualIndex = reverseSeq.length - 1 - i;
                const isPaired = pairedPositions2.has(actualIndex);
                const is3Prime = actualIndex >= reverseSeq.length - 5;
                const is5Prime = actualIndex < 5;
                const isFirst3PrimeDisplay = i === 0 && is3Prime;
                const isLast3PrimeDisplay = actualIndex === reverseSeq.length - 5;

                let color = theme.textMuted;
                let bg = 'transparent';

                if (isPaired && is3Prime) {
                  color = theme.danger;
                  bg = 'rgba(220, 38, 38, 0.15)';
                } else if (isPaired && is5Prime && hasLigationJunction) {
                  color = theme.success;
                  bg = 'rgba(22, 163, 74, 0.15)';
                } else if (isPaired) {
                  color = theme.warning;
                  bg = 'rgba(202, 138, 4, 0.1)';
                } else if (is3Prime) {
                  color = '#7c3aed';
                }

                const boxStyle = is3Prime && involves3Prime2 ? {
                  borderTop: `2px solid ${theme.danger}`,
                  borderBottom: `2px solid ${theme.danger}`,
                  borderLeft: isFirst3PrimeDisplay ? `2px solid ${theme.danger}` : 'none',
                  borderRight: isLast3PrimeDisplay ? `2px solid ${theme.danger}` : 'none',
                  borderRadius: isFirst3PrimeDisplay ? '3px 0 0 3px' : isLast3PrimeDisplay ? '0 3px 3px 0' : '0',
                  marginTop: '-2px', marginBottom: '-2px', paddingTop: '2px', paddingBottom: '2px',
                } : {};

                return (
                  <span key={i} style={{ color, background: bg, padding: '0 1px', fontWeight: isPaired ? 'bold' : 'normal', display: 'inline-block', width: '10px', textAlign: 'center', ...boxStyle }}>
                    {base}
                  </span>
                );
              })}
              <span style={{ color: '#7c3aed', marginLeft: '4px' }}>— 5'</span>
            </div>
          </div>
        );
      })()}

      {/* Stats */}
      <div style={{ display: 'flex', gap: '16px', marginTop: '12px', padding: '10px', background: theme.bgSecondary, borderRadius: '6px' }}>
        <div style={{ fontSize: '11px' }}>
          <span style={{ color: theme.textMuted }}>Paired bases: </span>
          <span style={{ color: theme.text, fontWeight: 'bold' }}>{pairs.length}</span>
        </div>
        <div style={{ fontSize: '11px' }}>
          <span style={{ color: theme.textMuted }}>Max consecutive: </span>
          <span style={{ color: theme.text, fontWeight: 'bold' }}>{maxConsecutive}</span>
        </div>
        <div style={{ fontSize: '11px' }}>
          <span style={{ color: theme.textMuted }}>FWD 3' involved: </span>
          <span style={{ color: involves3Prime1 ? theme.danger : theme.success, fontWeight: 'bold' }}>
            {involves3Prime1 ? 'YES ⚠️' : 'No'}
          </span>
        </div>
        <div style={{ fontSize: '11px' }}>
          <span style={{ color: theme.textMuted }}>REV 3' involved: </span>
          <span style={{ color: involves3Prime2 ? theme.danger : theme.success, fontWeight: 'bold' }}>
            {involves3Prime2 ? 'YES ⚠️' : 'No'}
          </span>
        </div>
      </div>

      {/* Dynamic suggestions */}
      <div style={{ marginTop: '12px', padding: '12px', background: colors.bg, borderRadius: '6px', border: `1px solid ${colors.border}`, fontSize: '12px' }}>
        <div style={{ fontWeight: 'bold', color: colors.text, marginBottom: '8px', display: 'flex', alignItems: 'center', gap: '6px' }}>
          {severity === 'critical' ? '⚠️ Critical: 3\' Extensible Dimer!' :
           severity === 'warning' ? '⚡ Attention Recommended' :
           hasLigationJunction ? '🔗 SDM Ligation Junction Detected' : '✓ Good Primer Pair'}
        </div>

        {severity === 'critical' ? (
          <div style={{ color: colors.text }}>
            <p style={{ margin: '0 0 8px 0' }}>
              <strong>3' ends can extend off each other</strong> ({maxConsecutive} consecutive bp).
              This is the most severe type of primer dimer - polymerase can extend both 3' ends, creating artifacts.
            </p>
            <div style={{ marginTop: '10px', paddingTop: '10px', borderTop: `1px solid ${theme.danger}33` }}>
              <strong>💡 Suggestions:</strong>
              <ul style={{ margin: '6px 0 0 16px', padding: 0, color: theme.text, lineHeight: '1.6' }}>
                <li>Try shifting one primer by 2-3 bases</li>
                {involves3Prime1 && <li>Modify forward primer's 3' end</li>}
                {involves3Prime2 && <li>Modify reverse primer's 3' end</li>}
                <li>Use hot-start polymerase</li>
              </ul>
            </div>
          </div>
        ) : severity === 'warning' ? (
          <div style={{ color: colors.text }}>
            <p style={{ margin: '0 0 8px 0' }}>
              {dimerType === 'internal_with_3prime'
                ? `Internal binding with one 3' end involved (${pairs.length} bp). This may cause some dimer formation but 3' end is not fully hybridized.`
                : `Internal/5' binding detected (${pairs.length} bp). The 3' ends are free (dangling), so this is annoying but often PCR-viable.`}
            </p>
            <div style={{ marginTop: '10px', paddingTop: '10px', borderTop: `1px solid ${theme.warning}33` }}>
              <strong>💡 Suggestions:</strong>
              <ul style={{ margin: '6px 0 0 16px', padding: 0, color: theme.text, lineHeight: '1.6' }}>
                <li>Monitor for primer-dimer bands</li>
                <li>Consider increasing annealing temp by 2-3°C</li>
                <li>Reduce primer concentration if dimers observed</li>
              </ul>
            </div>
          </div>
        ) : hasLigationJunction ? (
          <div style={{ color: colors.text }}>
            <p style={{ margin: '0' }}>
              <strong>SDM back-to-back design:</strong> The 5' ends show complementarity which is <em>expected and intended</em>.
              This forms the ligation junction for plasmid circularization after KLD enzyme treatment.
            </p>
            <p style={{ margin: '8px 0 0 0', color: theme.textSecondary, fontSize: '11px' }}>
              ✓ 3' ends are free (not hybridized) - primers will work correctly<br/>
              ✓ 5' junction allows efficient KLD ligation
            </p>
          </div>
        ) : (
          <div style={{ color: colors.text }}>
            <p style={{ margin: '0' }}>
              {pairs.length > 0
                ? `Minor complementarity (${pairs.length} bp) is within acceptable limits.`
                : 'No significant cross-dimer potential detected.'}
            </p>
          </div>
        )}
      </div>
    </div>
  );
}

// ============================================================================
// COMBINED HEATMAP WITH SPARKLINE OVERLAY
// ============================================================================

function ThermodynamicHeatmap({
  sequence,
  templateRegion = null,
  width = 600,
  label = 'Primer'
}: ThermodynamicHeatmapProps) {
  const [hoveredBar, setHoveredBar] = useState<number | null>(null);

  // Calculate position-wise ΔG for sparkline
  const energies = useMemo(() => {
    if (!sequence || sequence.length < 10) return [];
    return calculatePositionwiseDG(sequence, 10);
  }, [sequence]);

  if (!sequence) return null;

  const getBindingStrength = (base: string, index: number): number => {
    const isGC = base === 'G' || base === 'C';
    const position3Prime = index >= sequence.length - 5;
    let score = isGC ? 0.8 : 0.4;
    if (position3Prime) score = Math.min(1, score + 0.2);
    return score;
  };

  const last2 = sequence.slice(-2);
  const hasGCClamp = (last2[0] === 'G' || last2[0] === 'C') || (last2[1] === 'G' || last2[1] === 'C');
  const strongGCClamp = (last2[0] === 'G' || last2[0] === 'C') && (last2[1] === 'G' || last2[1] === 'C');

  const gcCount = sequence.split('').filter(b => b === 'G' || b === 'C').length;
  const gcContent = (gcCount / sequence.length * 100).toFixed(1);

  // Layout constants
  const barHeight = 40;
  const sparklineHeight = 50;
  const chartGap = 24; // Gap between binding strength bars and ΔG sparkline

  // Sparkline calculations - filter out any non-finite values
  const finiteEnergies = energies.filter(e => isFinite(e));
  const minE = finiteEnergies.length > 0 ? Math.min(...finiteEnergies, -5) : -5;
  const maxE = finiteEnergies.length > 0 ? Math.max(...finiteEnergies, 0) : 0;
  const range = maxE - minE || 1;

  return (
    <div style={{ background: theme.bg, borderRadius: '8px', padding: '12px', border: `1px solid ${theme.border}` }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '10px' }}>
        <span style={{ color: theme.text, fontWeight: 'bold', fontSize: '13px' }}>{label} Analysis</span>
        <div style={{ display: 'flex', gap: '12px', fontSize: '11px' }}>
          <span style={{ color: theme.textMuted }}>
            GC: <strong style={{ color: theme.text }}>{gcContent}%</strong>
          </span>
          <span style={{ color: theme.textMuted }}>
            3' Clamp: <strong style={{ color: strongGCClamp ? theme.success : hasGCClamp ? theme.warning : theme.danger }}>
              {strongGCClamp ? 'Strong' : hasGCClamp ? 'Weak' : 'None'}
            </strong>
          </span>
        </div>
      </div>

      <div style={{ background: theme.bgTertiary, borderRadius: '6px', padding: '12px', paddingTop: '40px', overflowX: 'auto' }}>
        {/* Y-axis label for binding strength */}
        <div style={{ display: 'flex' }}>
          <div style={{ display: 'flex', flexDirection: 'column', justifyContent: 'center', marginRight: '8px', height: `${barHeight}px` }}>
            <span style={{ writingMode: 'vertical-rl', transform: 'rotate(180deg)', fontSize: '9px', color: theme.textMuted, fontFamily: 'monospace' }}>
              Strength
            </span>
          </div>

          <div style={{ flex: 1, position: 'relative' }}>
            {/* Binding Strength Bars */}
            <div style={{ display: 'flex', height: `${barHeight}px`, position: 'relative' }}>
              {sequence.split('').map((base, i) => {
                const strength = getBindingStrength(base, i);
                const height = Math.round(strength * barHeight);
                const isGC = base === 'G' || base === 'C';
                const is3Prime = i >= sequence.length - 5;
                const isFirst3Prime = i === sequence.length - 5;
                const isLast = i === sequence.length - 1;
                const isHovered = hoveredBar === i;

                let color = isGC ? theme.gc : theme.at;
                if (is3Prime) color = isGC ? theme.gc : theme.prime3;

                return (
                  <div key={i} style={{
                    display: 'flex', flexDirection: 'column', alignItems: 'center', width: '18px',
                    position: 'relative', cursor: 'pointer',
                    borderTop: is3Prime ? `2px solid ${theme.prime3}50` : 'none',
                    borderBottom: is3Prime ? `2px solid ${theme.prime3}50` : 'none',
                    borderLeft: isFirst3Prime ? `2px solid ${theme.prime3}50` : 'none',
                    borderRight: isLast ? `2px solid ${theme.prime3}50` : 'none',
                    background: is3Prime ? `${theme.prime3}10` : 'transparent',
                    borderRadius: isFirst3Prime ? '4px 0 0 4px' : isLast ? '0 4px 4px 0' : '0',
                  }}
                    onMouseEnter={() => setHoveredBar(i)}
                    onMouseLeave={() => setHoveredBar(null)}
                  >
                    <div style={{
                      width: isHovered ? '16px' : '14px', height: `${height}px`,
                      background: color, opacity: isHovered ? 1 : (0.4 + strength * 0.4),
                      borderRadius: '2px 2px 0 0', marginTop: `${barHeight - height}px`,
                      transition: 'all 0.15s ease',
                      boxShadow: isHovered ? `0 0 8px ${color}` : 'none',
                    }} />

                    {isHovered && (
                      <div style={{
                        position: 'absolute', bottom: '100%', left: '50%', transform: 'translateX(-50%)',
                        background: theme.text, borderRadius: '4px', padding: '4px 8px',
                        fontSize: '9px', color: theme.bg, whiteSpace: 'nowrap', zIndex: 10, marginBottom: '4px', fontFamily: 'monospace'
                      }}>
                        <div>{base} @ pos {i + 1}</div>
                        <div>{(strength * 100).toFixed(0)}% • {energies[i]?.toFixed(1) || 'N/A'} kcal/mol</div>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>

            {/* Sequence labels */}
            <div style={{ display: 'flex', marginTop: '4px' }}>
              {sequence.split('').map((base, i) => {
                const isGC = base === 'G' || base === 'C';
                const is3Prime = i >= sequence.length - 5;
                return (
                  <span key={i} style={{
                    width: '18px', textAlign: 'center', fontFamily: 'monospace', fontSize: '12px',
                    color: is3Prime ? (isGC ? theme.gc : theme.prime3) : (isGC ? theme.gc : theme.textMuted),
                    fontWeight: is3Prime ? 'bold' : 'normal',
                    background: is3Prime ? `${theme.prime3}10` : 'transparent'
                  }}>
                    {base}
                  </span>
                );
              })}
            </div>

            {/* Position markers */}
            <div style={{ display: 'flex' }}>
              {sequence.split('').map((_, i) => (
                <span key={i} style={{ width: '18px', textAlign: 'center', fontSize: '8px', color: theme.textMuted, fontFamily: 'monospace' }}>
                  {(i + 1) % 5 === 0 ? i + 1 : ''}
                </span>
              ))}
            </div>
          </div>
        </div>

        {/* === GAP BETWEEN SECTIONS === */}
        <div style={{ height: `${chartGap}px` }} />

        {/* ΔG Sparkline Section */}
        {energies.length > 0 && (
          <div style={{ display: 'flex' }}>
            <div style={{ display: 'flex', flexDirection: 'column', justifyContent: 'space-between', marginRight: '8px', height: `${sparklineHeight}px` }}>
              <span style={{ fontSize: '8px', color: theme.textMuted, fontFamily: 'monospace' }}>{maxE.toFixed(0)}</span>
              <span style={{ writingMode: 'vertical-rl', transform: 'rotate(180deg)', fontSize: '9px', color: theme.textMuted, fontFamily: 'monospace' }}>
                ΔG
              </span>
              <span style={{ fontSize: '8px', color: theme.textMuted, fontFamily: 'monospace' }}>{minE.toFixed(0)}</span>
            </div>

            <div style={{ flex: 1, position: 'relative' }}>
              <svg width={sequence.length * 18} height={sparklineHeight} style={{ display: 'block' }}>
                {/* Background */}
                <rect x="0" y="0" width={sequence.length * 18} height={sparklineHeight} fill={theme.bgSecondary} rx="4" />

                {/* Zero line if applicable */}
                {minE < 0 && maxE > 0 && (
                  <line
                    x1="0"
                    y1={((maxE - 0) / range) * sparklineHeight}
                    x2={sequence.length * 18}
                    y2={((maxE - 0) / range) * sparklineHeight}
                    stroke={theme.border}
                    strokeWidth="1"
                  />
                )}

                {/* Critical threshold line at -3 kcal/mol */}
                {minE < -3 && (
                  <line
                    x1="0"
                    y1={((maxE - (-3)) / range) * sparklineHeight}
                    x2={sequence.length * 18}
                    y2={((maxE - (-3)) / range) * sparklineHeight}
                    stroke={theme.danger}
                    strokeWidth="1"
                    strokeDasharray="4,4"
                    opacity="0.5"
                  />
                )}

                {/* Sparkline path */}
                <path
                  d={energies.map((e, i) => {
                    const x = i * 18 + 9;
                    const safeE = isFinite(e) ? e : 0;
                    const y = ((maxE - safeE) / range) * sparklineHeight;
                    return `${i === 0 ? 'M' : 'L'} ${x} ${y}`;
                  }).join(' ')}
                  fill="none"
                  stroke={theme.text}
                  strokeWidth="2"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />

                {/* Data points */}
                {energies.map((e, i) => {
                  const x = i * 18 + 9;
                  const safeE = isFinite(e) ? e : 0;
                  const y = ((maxE - safeE) / range) * sparklineHeight;
                  const is3Prime = i >= sequence.length - 5;
                  return (
                    <circle
                      key={i}
                      cx={x}
                      cy={y}
                      r={is3Prime ? 3 : 2}
                      fill={is3Prime ? theme.danger : theme.text}
                      opacity={is3Prime ? 1 : 0.6}
                    />
                  );
                })}
              </svg>

              {/* Sparkline label */}
              <div style={{ marginTop: '4px', fontSize: '9px', color: theme.textMuted, fontFamily: 'monospace' }}>
                Position-wise ΔG (kcal/mol) — <span style={{ color: theme.danger }}>dashed line = -3 threshold</span>
              </div>
            </div>
          </div>
        )}
      </div>

      {/* Legend */}
      <div style={{ display: 'flex', gap: '16px', marginTop: '10px', fontSize: '10px', color: theme.textMuted, alignItems: 'center', flexWrap: 'wrap' }}>
        <span><span style={{ display: 'inline-block', width: '10px', height: '10px', background: theme.gc, borderRadius: '2px', marginRight: '4px' }}></span>G/C</span>
        <span><span style={{ display: 'inline-block', width: '10px', height: '10px', background: theme.at, borderRadius: '2px', marginRight: '4px' }}></span>A/T</span>
        <span><span style={{ display: 'inline-block', width: '10px', height: '10px', background: theme.prime3, opacity: 0.5, borderRadius: '2px', marginRight: '4px' }}></span>3' region</span>
        <span style={{ marginLeft: 'auto' }}><span style={{ display: 'inline-block', width: '16px', height: '2px', background: theme.text, marginRight: '4px', verticalAlign: 'middle' }}></span>ΔG profile</span>
      </div>
    </div>
  );
}

// ============================================================================
// MAIN EXPORT
// ============================================================================

export default function PrimerStructureViewer({
  forwardSeq,
  reverseSeq,
  forwardName = 'Forward Primer',
  reverseName = 'Reverse Primer',
  templateRegion = null,
  showCrossDimer = true,
  showHeatmap = true,
  showSparkline = false, // Deprecated - now merged into heatmap
  isSDMMode = false,
  onForwardShift = null,
  onReverseShift = null,
  useFornaViewer = true, // Use fornac library for smooth structure visualization
}: PrimerStructureViewerProps) {
  const fwdFold = useMemo(() => {
    if (!forwardSeq || forwardSeq.length < 6) return { e: 0, ij: [] };
    return foldSequence(forwardSeq);
  }, [forwardSeq]);

  const revFold = useMemo(() => {
    if (!reverseSeq || reverseSeq.length < 6) return { e: 0, ij: [] };
    return foldSequence(reverseSeq);
  }, [reverseSeq]);

  // Calculate appropriate dimensions based on sequence length
  const calcDimensions = (seq: string | undefined, pairs: BasePair[]) => {
    if (!seq) return { width: 450, height: 420 };
    const len = seq.length;
    const hasPairs = pairs && pairs.length > 0;
    // Scale width based on sequence length - larger for readability
    const width = Math.min(600, Math.max(450, len * 14));
    // Scale height based on sequence and pair count - larger for readability
    const height = hasPairs
      ? Math.min(600, Math.max(420, len * 12 + pairs.length * 25))
      : Math.min(500, Math.max(350, len * 10));
    return { width, height };
  };

  const fwdDims = calcDimensions(forwardSeq, fwdFold.ij);
  const revDims = calcDimensions(reverseSeq, revFold.ij);

  // Choose which structure viewer to use
  const StructureViewer = useFornaViewer ? FornaViewer : ForceDirectedStructure;

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
      {/* Structure visualization panels */}
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(450px, 1fr))', gap: '20px' }}>
        {forwardSeq && (
          <div style={{ overflow: 'auto', maxHeight: '550px' }}>
            <StructureViewer
              sequence={forwardSeq}
              basePairs={fwdFold.ij}
              width={fwdDims.width}
              height={fwdDims.height}
              label={forwardName}
              onShift={onForwardShift}
            />
          </div>
        )}
        {reverseSeq && (
          <div style={{ overflow: 'auto', maxHeight: '550px' }}>
            <StructureViewer
              sequence={reverseSeq}
              basePairs={revFold.ij}
              width={revDims.width}
              height={revDims.height}
              label={reverseName}
              onShift={onReverseShift}
            />
          </div>
        )}
      </div>

      {/* Cross-dimer check */}
      {showCrossDimer && forwardSeq && reverseSeq && (
        <CrossDimerZipper forwardSeq={forwardSeq} reverseSeq={reverseSeq} isSDMMode={isSDMMode} />
      )}

      {/* Combined heatmap with sparkline overlay */}
      {showHeatmap && (
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(350px, 1fr))', gap: '16px' }}>
          {forwardSeq && (
            <div style={{ overflow: 'auto', maxWidth: '100%' }}>
              <ThermodynamicHeatmap sequence={forwardSeq} label={forwardName} />
            </div>
          )}
          {reverseSeq && (
            <div style={{ overflow: 'auto', maxWidth: '100%' }}>
              <ThermodynamicHeatmap sequence={reverseSeq} label={reverseName} />
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export { ForceDirectedStructure, CrossDimerZipper, ThermodynamicHeatmap, FornaViewer };
