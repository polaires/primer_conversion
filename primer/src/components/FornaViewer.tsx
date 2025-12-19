import { useMemo, useState, FC } from 'react';
import { dg } from '../lib/fold.js';

// Theme interface
interface Theme {
  bg: string;
  bgSecondary: string;
  bgTertiary: string;
  border: string;
  text: string;
  textSecondary: string;
  textMuted: string;
  success: string;
  successBg: string;
  warning: string;
  warningBg: string;
  danger: string;
  dangerBg: string;
  prime3: string;
}

// Light theme colors
const theme: Theme = {
  bg: '#ffffff',
  bgSecondary: '#f8fafc',
  bgTertiary: '#f1f5f9',
  border: '#e2e8f0',
  text: '#1e293b',
  textSecondary: '#475569',
  textMuted: '#94a3b8',
  success: '#16a34a',
  successBg: '#dcfce7',
  warning: '#ca8a04',
  warningBg: '#fef9c3',
  danger: '#dc2626',
  dangerBg: '#fee2e2',
  prime3: '#2563eb',
};

// Type definitions
interface Node {
  id: number;
  base: string;
  x: number;
  y: number;
}

type BasePair = [number, number];

interface FornaViewerProps {
  sequence: string;
  basePairs?: BasePair[];
  width?: number;
  height?: number;
  label?: string;
  onShift?: ((direction: number) => void) | null;
}

/**
 * Row-based hairpin layout - proper stem-loop structure with dynamic spacing
 *
 * Key principles:
 * 1. Calculate spacing between pairs dynamically based on bulges
 * 2. Both bases in a pair get EXACTLY the same Y coordinate
 * 3. Bulge bases are distributed evenly in the stretched gap
 */
function calculateClassicHairpinLayout(
  sequence: string,
  basePairs: BasePair[],
  width: number,
  height: number
): Node[] {
  const n = sequence.length;
  const nodes: Node[] = new Array(n);
  const centerX = width / 2;

  console.log('[FornaViewer] === Layout Debug ===');
  console.log('[FornaViewer] sequence:', sequence);
  console.log('[FornaViewer] basePairs:', JSON.stringify(basePairs));
  console.log('[FornaViewer] dimensions:', width, 'x', height);

  // Spacing parameters
  const unitSpacing = 34; // Base unit for vertical spacing
  const stemWidth = 90;   // Horizontal distance between left and right stems

  if (!basePairs || basePairs.length === 0) {
    console.log('[FornaViewer] No base pairs - linear layout');
    // No structure: simple horizontal line with gentle curve
    const startX = 60;
    const endX = width - 60;
    const midY = height / 2;
    for (let i = 0; i < n; i++) {
      const t = n > 1 ? i / (n - 1) : 0.5;
      const x = startX + t * (endX - startX);
      const y = midY + Math.sin(t * Math.PI) * 25;
      nodes[i] = { id: i, base: sequence[i], x, y };
    }
    return nodes;
  }

  // Build pair map (index -> partner index)
  const pairMap = new Map<number, number>();
  basePairs.forEach(([i, j]) => {
    pairMap.set(i, j);
    pairMap.set(j, i);
  });

  // Sort pairs by first index (outermost to innermost on left side)
  const sortedPairs = [...basePairs].sort((a, b) => a[0] - b[0]);
  const numPairs = sortedPairs.length;

  // Find the innermost pair (smallest gap between indices = loop boundary)
  const innerPair = sortedPairs.reduce((inner, pair) =>
    (pair[1] - pair[0]) < (inner[1] - inner[0]) ? pair : inner
  , sortedPairs[0]);

  // Structure boundaries
  // firstPairedIdx = minimum of all left indices (first in sorted pairs)
  const firstPairedIdx = sortedPairs[0][0];
  // lastPairedIdx = maximum of all right indices (NOT just last sorted pair!)
  const lastPairedIdx = Math.max(...sortedPairs.map(p => p[1]));
  const loopStart = innerPair[0] + 1;
  const loopEnd = innerPair[1] - 1;
  const loopLen = Math.max(1, loopEnd - loopStart + 1);

  // Tail length (5' end)
  const tail5Len = firstPairedIdx;

  // X positions
  const leftX = centerX - stemWidth / 2;
  const rightX = centerX + stemWidth / 2;

  // === STEP 1: Calculate dynamic spacing for each segment ===
  // For each pair of adjacent base pairs, count bulges and determine segment height

  // segmentHeights[i] = vertical units between pair i and pair i+1
  // (pair 0 is outermost/bottom, pair numPairs-1 is innermost/top)
  const segmentHeights: number[] = [];

  for (let p = 0; p < numPairs - 1; p++) {
    const [leftA, rightA] = sortedPairs[p];      // Lower pair (bigger indices spread)
    const [leftB, rightB] = sortedPairs[p + 1];  // Upper pair (smaller indices spread)

    // Count unpaired bases on left side between leftA and leftB
    let leftBulgeCount = 0;
    for (let i = leftA + 1; i < leftB; i++) {
      if (!pairMap.has(i)) leftBulgeCount++;
    }

    // Count unpaired bases on right side between rightB and rightA
    let rightBulgeCount = 0;
    for (let i = rightB + 1; i < rightA; i++) {
      if (!pairMap.has(i)) rightBulgeCount++;
    }

    // Segment height = max bulges + 1 (for the pair itself)
    const maxBulges = Math.max(leftBulgeCount, rightBulgeCount);
    segmentHeights.push(maxBulges + 1);
  }

  // === STEP 2: Calculate Y positions for each pair ===
  const loopTopY = 70;
  const stemTopY = loopTopY + 60; // Innermost pair Y position

  // Build Y positions from top (innermost) to bottom (outermost)
  const pairYPositions: number[] = new Array(numPairs);
  pairYPositions[numPairs - 1] = stemTopY; // Innermost pair at top

  for (let p = numPairs - 2; p >= 0; p--) {
    // Distance from pair p+1 down to pair p
    pairYPositions[p] = pairYPositions[p + 1] + segmentHeights[p] * unitSpacing;
  }

  const stemBottomY = pairYPositions[0]; // Outermost pair Y

  // === STEP 3: Create Y map for paired bases ===
  const pairYMap = new Map<number, number>();
  const pairRankMap = new Map<number, number>(); // Maps base index to its pair rank

  console.log('[FornaViewer] sortedPairs:', JSON.stringify(sortedPairs));
  console.log('[FornaViewer] segmentHeights:', JSON.stringify(segmentHeights));
  console.log('[FornaViewer] pairYPositions:', JSON.stringify(pairYPositions));

  for (let p = 0; p < numPairs; p++) {
    const [a, b] = sortedPairs[p];
    const y = pairYPositions[p];
    pairYMap.set(a, y);
    pairYMap.set(b, y); // CRITICAL: Both bases get EXACTLY same Y
    pairRankMap.set(a, p);
    pairRankMap.set(b, p);
    console.log(`[FornaViewer] Pair ${p}: indices [${a}, ${b}] → Y = ${y}`);
  }

  // === STEP 4: Position all nodes ===
  const loopRadius = Math.max(50, loopLen * 18);

  for (let i = 0; i < n; i++) {
    let x: number, y: number;
    const partner = pairMap.get(i);
    const isPaired = partner !== undefined;

    if (i < firstPairedIdx) {
      // === 5' TAIL ===
      const distFromStemBottom = tail5Len - i;
      x = leftX;
      y = stemBottomY + distFromStemBottom * unitSpacing;

    } else if (i > lastPairedIdx) {
      // === 3' TAIL ===
      const tailIdx = i - lastPairedIdx - 1;
      x = rightX;
      y = stemBottomY + (tailIdx + 1) * unitSpacing;

    } else if (i >= loopStart && i <= loopEnd) {
      // === LOOP ===
      const loopIdx = i - loopStart;
      const startAngle = Math.PI * 0.85;
      const endAngle = Math.PI * 0.15;
      const angle = startAngle - (startAngle - endAngle) * ((loopIdx + 0.5) / loopLen);

      x = centerX + Math.cos(angle) * loopRadius;
      y = loopTopY - Math.sin(angle) * loopRadius * 0.65;

    } else if (isPaired) {
      // === STEM (paired base) ===
      // Use pre-calculated Y - never interpolate paired bases!
      const isLeftStem = i < partner!;
      y = pairYMap.get(i)!;
      x = isLeftStem ? leftX : rightX;

    } else {
      // === BULGE (unpaired base within stem region) ===
      const isLeftSide = i < (innerPair[0] + innerPair[1]) / 2;
      x = isLeftSide ? leftX : rightX;

      // Find the segment this bulge belongs to
      // We need the CLOSEST pair above and below the bulge index
      let lowerPairRank = -1;
      let upperPairRank = -1;

      for (let p = 0; p < numPairs; p++) {
        const stemIdx = isLeftSide ? sortedPairs[p][0] : sortedPairs[p][1];

        if (isLeftSide) {
          // Left side: smaller index = lower position (larger Y)
          if (stemIdx < i) {
            // This pair is below the bulge - find closest (largest stemIdx < i)
            if (lowerPairRank === -1 || stemIdx > sortedPairs[lowerPairRank][0]) {
              lowerPairRank = p;
            }
          }
          if (stemIdx > i) {
            // This pair is above the bulge - find closest (smallest stemIdx > i)
            if (upperPairRank === -1 || stemIdx < sortedPairs[upperPairRank][0]) {
              upperPairRank = p;
            }
          }
        } else {
          // Right side: smaller index = UPPER position (smaller Y, near loop)
          if (stemIdx < i) {
            // This pair is ABOVE the bulge in Y - find closest (largest stemIdx < i)
            if (upperPairRank === -1 || stemIdx > sortedPairs[upperPairRank][1]) {
              upperPairRank = p;
            }
          }
          if (stemIdx > i) {
            // This pair is BELOW the bulge in Y - find closest (smallest stemIdx > i)
            if (lowerPairRank === -1 || stemIdx < sortedPairs[lowerPairRank][1]) {
              lowerPairRank = p;
            }
          }
        }
      }

      // Get the boundary indices and Y positions
      const lowerIdx = isLeftSide ? sortedPairs[lowerPairRank][0] : sortedPairs[lowerPairRank][1];
      const upperIdx = isLeftSide ? sortedPairs[upperPairRank][0] : sortedPairs[upperPairRank][1];
      const lowerY = pairYPositions[lowerPairRank];
      const upperY = pairYPositions[upperPairRank];

      // Count total unpaired bases in this segment on this side
      const bulgeIndices: number[] = [];
      const rangeStart = isLeftSide ? lowerIdx + 1 : upperIdx + 1;
      const rangeEnd = isLeftSide ? upperIdx : lowerIdx;

      for (let j = rangeStart; j < rangeEnd; j++) {
        if (!pairMap.has(j)) bulgeIndices.push(j);
      }

      // Find position of current base in the bulge sequence
      const bulgePosition = bulgeIndices.indexOf(i);
      const totalBulges = bulgeIndices.length;

      // Distribute evenly between the two pairs
      if (totalBulges > 0 && bulgePosition >= 0) {
        let t: number;
        if (isLeftSide) {
          // Left side: indices increase going UP (toward loop)
          // Position 0 (smallest index) is closest to lowerY (bottom)
          t = (bulgePosition + 1) / (totalBulges + 1);
        } else {
          // Right side: indices increase going DOWN (toward 3' end)
          // Position 0 (smallest index) is closest to upperY (top)
          // So we need to invert: position 0 → close to top, position N → close to bottom
          t = (totalBulges - bulgePosition) / (totalBulges + 1);
        }
        y = lowerY + (upperY - lowerY) * t;
      } else {
        y = (lowerY + upperY) / 2;
      }
    }

    nodes[i] = { id: i, base: sequence[i], x, y };
  }

  // Log final positions
  console.log('[FornaViewer] Final node positions:');
  nodes.forEach((node, i) => {
    const paired = pairMap.get(i);
    const pairedY = paired !== undefined ? nodes[paired]?.y : null;
    const yMatch = pairedY !== null ? (node.y === pairedY ? '✓' : `✗ partner Y=${pairedY}`) : '';
    console.log(`  [${i}] ${node.base} → (${node.x.toFixed(1)}, ${node.y.toFixed(1)}) ${yMatch}`);
  });

  return nodes;
}

/**
 * Generate simple polyline path - straight lines with rounded corners
 */
function generateSimplePath(nodes: Node[]): string {
  if (nodes.length < 2) return '';

  // Just connect nodes with straight lines - CSS will round the corners
  let path = `M ${nodes[0].x} ${nodes[0].y}`;
  for (let i = 1; i < nodes.length; i++) {
    path += ` L ${nodes[i].x} ${nodes[i].y}`;
  }
  return path;
}

/**
 * FornaViewer - Improved hairpin structure visualization
 */
const FornaViewer: FC<FornaViewerProps> = ({
  sequence,
  basePairs = [],
  width = 420,  // Larger default
  height = 380, // Larger default
  label = 'Primer',
  onShift = null,
}) => {
  const [hoveredNode, setHoveredNode] = useState<number | null>(null);

  // Calculate node positions
  const nodes = useMemo((): Node[] => {
    if (!sequence) return [];
    return calculateClassicHairpinLayout(sequence, basePairs, width, height);
  }, [sequence, basePairs, width, height]);

  // Generate backbone path
  const backbonePath = useMemo(() => generateSimplePath(nodes), [nodes]);

  // Calculate energy
  const energy = useMemo((): number => {
    if (!sequence || sequence.length < 6) return 0;
    return dg(sequence, 37) as number;
  }, [sequence]);

  // Calculate viewBox to fit content
  const viewBox = useMemo((): string => {
    if (nodes.length === 0) return `0 0 ${width} ${height}`;

    const padding = 55; // Extra padding for 5'/3' labels
    const xs = nodes.map(n => n.x);
    const ys = nodes.map(n => n.y);
    const minX = Math.min(...xs) - padding;
    const maxX = Math.max(...xs) + padding;
    const minY = Math.min(...ys) - padding;
    const maxY = Math.max(...ys) + padding;

    return `${minX} ${minY} ${maxX - minX} ${maxY - minY}`;
  }, [nodes, width, height]);

  const energyColor = energy >= -1 ? theme.success : energy >= -3 ? theme.warning : theme.danger;

  const getNodeColor = (base: string, index: number): string => {
    const is3Prime = index >= sequence.length - 5;
    const isGC = base === 'G' || base === 'C';
    if (is3Prime) return isGC ? theme.success : theme.prime3;
    return isGC ? theme.success : '#ea580c';
  };

  if (!sequence) return null;

  const svgHeight = height - 70; // Account for header and footer

  return (
    <div className="rounded-lg p-4 overflow-hidden" style={{
      background: theme.bg,
      border: `1px solid ${theme.border}`,
    }}>
      {/* Header */}
      <div className="flex justify-between items-center mb-3">
        <div className="flex items-center gap-2.5">
          <span className="font-bold text-[15px]" style={{ color: theme.text }}>{label}</span>
          {onShift && (
            <div className="flex gap-0.5">
              <button
                onClick={() => onShift(-1)}
                className="rounded px-2 py-1 cursor-pointer text-sm"
                style={{
                  background: theme.bgSecondary,
                  border: `1px solid ${theme.border}`,
                }}
                title="Shift primer 1bp upstream"
              >
                ◀
              </button>
              <button
                onClick={() => onShift(1)}
                className="rounded px-2 py-1 cursor-pointer text-sm"
                style={{
                  background: theme.bgSecondary,
                  border: `1px solid ${theme.border}`,
                }}
                title="Shift primer 1bp downstream"
              >
                ▶
              </button>
            </div>
          )}
        </div>
        <span className="text-sm font-mono font-bold" style={{ color: energyColor }}>
          ΔG = {energy.toFixed(1)} kcal/mol
        </span>
      </div>

      {/* SVG Visualization */}
      <svg
        width="100%"
        height={svgHeight}
        viewBox={viewBox}
        preserveAspectRatio="xMidYMid meet"
        className="rounded-lg block"
        style={{ background: theme.bgTertiary }}
      >
        <defs>
          <linearGradient
            id={`backbone-grad-${label.replace(/\s/g, '')}`}
            gradientUnits="userSpaceOnUse"
            x1={nodes[0]?.x || 0}
            y1={nodes[0]?.y || 0}
            x2={nodes[nodes.length - 1]?.x || 0}
            y2={nodes[nodes.length - 1]?.y || 0}
          >
            <stop offset="0%" stopColor={theme.prime3} />
            <stop offset="100%" stopColor={theme.danger} />
          </linearGradient>
        </defs>

        {/* Smooth backbone curve */}
        {backbonePath && (
          <path
            d={backbonePath}
            fill="none"
            stroke={`url(#backbone-grad-${label.replace(/\s/g, '')})`}
            strokeWidth="3.5"
            strokeLinecap="round"
            strokeLinejoin="round"
            opacity="0.75"
          />
        )}

        {/* Base pair connections */}
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
                strokeWidth={involves3Prime ? 3 : 2.5}
                strokeDasharray={involves3Prime ? '0' : '5,4'}
                opacity={0.85}
              />
              <circle
                cx={midX}
                cy={midY}
                r="4"
                fill={involves3Prime ? theme.danger : '#ea580c'}
                opacity="0.7"
              />
            </g>
          );
        })}

        {/* Nucleotide nodes */}
        {nodes.map((node, i) => {
          const color = getNodeColor(node.base, i);
          const isHovered = hoveredNode === i;
          const is3Prime = i >= sequence.length - 5;
          const nodeRadius = isHovered ? 15 : 13;

          return (
            <g
              key={i}
              onMouseEnter={() => setHoveredNode(i)}
              onMouseLeave={() => setHoveredNode(null)}
              style={{ cursor: 'pointer' }}
            >
              {/* 3' region indicator ring */}
              {is3Prime && (
                <circle
                  cx={node.x}
                  cy={node.y}
                  r={nodeRadius + 4}
                  fill="none"
                  stroke={theme.danger}
                  strokeWidth="1.5"
                  strokeDasharray="3,2"
                  opacity="0.6"
                />
              )}
              {/* Main node circle */}
              <circle
                cx={node.x}
                cy={node.y}
                r={nodeRadius}
                fill={theme.bg}
                stroke={color}
                strokeWidth={is3Prime ? 3 : 2}
                style={{
                  filter: isHovered ? 'drop-shadow(0 2px 6px rgba(0,0,0,0.25))' : 'none',
                  transition: 'all 0.15s ease'
                }}
              />
              {/* Base letter */}
              <text
                x={node.x}
                y={node.y + 5}
                fill={color}
                fontSize="13"
                fontFamily="ui-monospace, SFMono-Regular, monospace"
                textAnchor="middle"
                fontWeight="bold"
              >
                {node.base}
              </text>
              {/* Tooltip */}
              {isHovered && (
                <g>
                  <rect
                    x={node.x - 50}
                    y={node.y - 45}
                    width={100}
                    height={26}
                    rx={5}
                    fill={theme.text}
                    opacity="0.95"
                  />
                  <text
                    x={node.x}
                    y={node.y - 27}
                    fill={theme.bg}
                    fontSize="11"
                    fontFamily="ui-monospace, monospace"
                    textAnchor="middle"
                  >
                    {`${node.base}${i + 1}${is3Prime ? " (3' region)" : ''}`}
                  </text>
                </g>
              )}
            </g>
          );
        })}

        {/* 5' and 3' labels */}
        {nodes.length > 0 && (
          <>
            <g>
              <circle
                cx={nodes[0].x - 30}
                cy={nodes[0].y}
                r="14"
                fill={theme.prime3}
                opacity="0.2"
              />
              <text
                x={nodes[0].x - 30}
                y={nodes[0].y + 4}
                fill={theme.prime3}
                fontSize="11"
                fontFamily="ui-monospace, monospace"
                textAnchor="middle"
                fontWeight="bold"
              >
                5'
              </text>
            </g>
            <g>
              <circle
                cx={nodes[nodes.length - 1].x + 30}
                cy={nodes[nodes.length - 1].y}
                r="14"
                fill={theme.danger}
                opacity="0.2"
              />
              <text
                x={nodes[nodes.length - 1].x + 30}
                y={nodes[nodes.length - 1].y + 4}
                fill={theme.danger}
                fontSize="11"
                fontFamily="ui-monospace, monospace"
                textAnchor="middle"
                fontWeight="bold"
              >
                3'
              </text>
            </g>
          </>
        )}
      </svg>

      {/* Status bar */}
      <div className="rounded-md p-2 px-3.5 mt-3 text-center" style={{
        background: basePairs.length > 0 ? theme.warningBg : theme.successBg,
      }}>
        <span className="text-xs font-mono font-medium" style={{
          color: basePairs.length > 0 ? theme.warning : theme.success,
        }}>
          {basePairs.length > 0
            ? `${basePairs.length} base pair${basePairs.length > 1 ? 's' : ''} – Hairpin detected`
            : 'Linear structure – No secondary structure'}
        </span>
      </div>

      {/* Legend */}
      <div className="flex gap-4 mt-2.5 text-[11px] justify-center" style={{ color: theme.textMuted }}>
        <span><span style={{ color: theme.success }}>●</span> G/C</span>
        <span><span style={{ color: '#ea580c' }}>●</span> A/T</span>
        <span><span style={{ color: theme.prime3 }}>●</span> 3' region</span>
        <span><span style={{ color: theme.danger }}>●</span> High risk</span>
      </div>
    </div>
  );
};

export default FornaViewer;
export { FornaViewer };
