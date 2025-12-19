import { useMemo, useState, FC } from 'react';
import { foldSequence } from '../lib/fold.js';
import { classify3PrimeStructureSeverity } from '../lib/smartPrimers.js';

/**
 * HairpinDiagram - Draws actual stem-loop/hairpin structures
 *
 * Instead of just showing arc diagrams or ΔG numbers, this component
 * renders the physical structure of the hairpin so users can immediately
 * see problematic secondary structures.
 *
 * A massive stem-loop at the 3' end convinces a user to redesign
 * much faster than just seeing "-5.4 kcal/mol".
 */

// Type definitions
type BasePair = [number, number];

interface FoldResult {
  e?: number;
  ij?: BasePair[];
  [key: string]: unknown;
}

interface StemPair {
  i: number;
  j: number;
  bases: [string, string];
}

interface HairpinStructure {
  stemPairs: StemPair[];
  loopStart: number;
  loopEnd: number;
  loopSequence: string;
  paired: Set<number>;
  pairMap: Map<number, number>;
  has3PrimeStructure: boolean;
  has3PrimeInStem: boolean;
  has3PrimeInLoop: boolean;
  threePrimeRegion: number;
}

interface SeverityResult {
  level: 'none' | 'info' | 'low' | 'moderate' | 'warning' | 'critical';
  label: string;
  shortLabel: string;
  message?: string;
  tooltip: string;
  shouldWarn: boolean;
}

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
  prime5: string;
  prime3: string;
  nodeBg?: string;
  node3PrimeBg?: string;
  backbone?: string;
  prime3Bg?: string;
}

interface StemLoopDiagramProps {
  sequence: string;
  structure: HairpinStructure;
  energy: number;
  width?: number;
  height?: number;
}

interface FlatStructureDiagramProps {
  sequence: string;
  basePairs: BasePair[];
  energy: number;
  width?: number;
  height?: number;
}

interface HairpinDiagramProps {
  sequence: string;
  foldResult?: FoldResult | null;
  primerName?: string;
  width?: number;
  showDetails?: boolean;
  temperature?: number;
  mode?: 'auto' | 'stemloop' | 'flat';
}

interface HairpinBadgeProps {
  sequence: string;
  foldResult?: FoldResult | null;
  lightTheme?: boolean;
  showTooltip?: boolean;
  temperature?: number;
}

// Color scale based on energy (more negative = worse for primers)
function energyToColor(dG: number): string {
  if (dG >= -1) return '#22c55e';      // Green - minimal
  if (dG >= -2) return '#84cc16';      // Light green
  if (dG >= -3) return '#eab308';      // Yellow - moderate
  if (dG >= -4) return '#f97316';      // Orange - problematic
  return '#ef4444';                     // Red - severe
}

function energySeverity(dG: number): { level: string; label: string } {
  if (dG >= -1) return { level: 'minimal', label: 'Minimal Structure' };
  if (dG >= -2) return { level: 'low', label: 'Low Stability' };
  if (dG >= -3) return { level: 'moderate', label: 'Moderate Structure' };
  if (dG >= -4) return { level: 'high', label: 'High Stability - Consider Redesign' };
  return { level: 'severe', label: 'Severe Structure - Redesign Recommended!' };
}

/**
 * Parse fold results to identify hairpin structure elements
 */
function parseHairpinStructure(sequence: string, basePairs: BasePair[]): HairpinStructure | null {
  if (!basePairs || basePairs.length === 0) {
    return null;
  }

  // Sort pairs by position
  const sortedPairs = [...basePairs].sort((a, b) => a[0] - b[0]);

  // Build structure info
  const paired = new Set<number>();
  const pairMap = new Map<number, number>();

  for (const [i, j] of sortedPairs) {
    paired.add(i);
    paired.add(j);
    pairMap.set(i, j);
    pairMap.set(j, i);
  }

  // Find stem and loop regions
  const stemPairs: StemPair[] = [];
  let loopStart = -1;
  let loopEnd = -1;

  // For a simple hairpin, pairs are nested: (0,n), (1,n-1), etc.
  // Find the innermost pair to identify the loop
  for (const [i, j] of sortedPairs) {
    stemPairs.push({ i, j, bases: [sequence[i], sequence[j]] });

    // The loop is between the innermost pair
    const innerUnpaired = j - i - 1;
    if (innerUnpaired > 0 && (loopStart === -1 || i > loopStart)) {
      // Check if there are unpaired bases between this pair
      let hasInnerPairs = false;
      for (let k = i + 1; k < j; k++) {
        const partner = pairMap.get(k);
        if (partner !== undefined && partner > k && partner < j) {
          hasInnerPairs = true;
          break;
        }
      }
      if (!hasInnerPairs) {
        loopStart = i + 1;
        loopEnd = j;
      }
    }
  }

  // Get loop sequence
  const loopSequence = loopStart >= 0 ? sequence.slice(loopStart, loopEnd) : '';

  // Identify 3' region involvement (last 10 bases)
  const seqLen = sequence.length;
  const threePrimeRegion = seqLen - 10;
  const has3PrimeInStem = stemPairs.some(p => p.i >= threePrimeRegion || p.j >= threePrimeRegion);
  const has3PrimeInLoop = loopStart >= threePrimeRegion;

  return {
    stemPairs,
    loopStart,
    loopEnd,
    loopSequence,
    paired,
    pairMap,
    has3PrimeStructure: has3PrimeInStem || has3PrimeInLoop,
    has3PrimeInStem,
    has3PrimeInLoop,
    threePrimeRegion
  };
}

/**
 * SVG Hairpin Structure Renderer
 * Draws the classic stem-loop diagram with the loop at the top
 * White theme version
 */
const StemLoopDiagram: FC<StemLoopDiagramProps> = ({
  sequence,
  structure,
  energy,
  width = 400,
  height = 300
}) => {
  const { stemPairs, loopSequence, threePrimeRegion, has3PrimeStructure } = structure;

  const stemLength = stemPairs.length;
  const loopLength = loopSequence.length;

  // Calculate dimensions
  const centerX = width / 2;
  const stemSpacing = 20;     // Horizontal spacing between stem sides
  const baseHeight = 18;      // Vertical spacing between base pairs
  const loopRadius = Math.max(30, loopLength * 4);
  const fontSize = Math.min(12, baseHeight - 4);

  const stemStartY = height - 40;
  const stemTopY = stemStartY - (stemLength * baseHeight);

  const color = energyToColor(energy);

  // White theme colors
  const theme = {
    bg: '#f8fafc',
    nodeBg: '#ffffff',
    node3PrimeBg: '#fee2e2',
    backbone: '#cbd5e1',
    text: '#1e293b',
    textMuted: '#64748b',
    prime5: '#2563eb',
    prime3: '#dc2626',
  };

  return (
    <svg width={width} height={height} style={{ background: theme.bg, borderRadius: '6px' }}>
      <defs>
        {/* Glow filter for 3' warnings */}
        <filter id="glow3prime">
          <feGaussianBlur stdDeviation="2" result="coloredBlur"/>
          <feMerge>
            <feMergeNode in="coloredBlur"/>
            <feMergeNode in="SourceGraphic"/>
          </feMerge>
        </filter>
      </defs>

      {/* Title */}
      <text x={20} y={25} fill={theme.text} fontSize="13" fontFamily="monospace" fontWeight="bold">
        Hairpin Structure
      </text>
      <text x={width - 20} y={25} fill={color} fontSize="12" fontFamily="monospace" textAnchor="end">
        ΔG = {energy.toFixed(1)} kcal/mol
      </text>

      {/* Draw stem backbone lines */}
      <line
        x1={centerX - stemSpacing}
        y1={stemStartY}
        x2={centerX - stemSpacing}
        y2={stemTopY}
        stroke={theme.backbone}
        strokeWidth="2"
      />
      <line
        x1={centerX + stemSpacing}
        y1={stemStartY}
        x2={centerX + stemSpacing}
        y2={stemTopY}
        stroke={theme.backbone}
        strokeWidth="2"
      />

      {/* Draw stem base pairs */}
      {stemPairs.map((pair, idx) => {
        const y = stemStartY - (idx * baseHeight) - 10;
        const is3Prime = pair.i >= threePrimeRegion || pair.j >= threePrimeRegion;
        const pairColor = is3Prime ? theme.prime3 : color;

        return (
          <g key={`pair-${idx}`} filter={is3Prime ? 'url(#glow3prime)' : undefined}>
            {/* Base pair connecting line */}
            <line
              x1={centerX - stemSpacing + 8}
              y1={y}
              x2={centerX + stemSpacing - 8}
              y2={y}
              stroke={pairColor}
              strokeWidth={is3Prime ? 2 : 1}
              strokeDasharray={is3Prime ? '' : '2,2'}
              opacity={0.6}
            />

            {/* Left base (5' strand) */}
            <circle
              cx={centerX - stemSpacing}
              cy={y}
              r={7}
              fill={is3Prime ? theme.node3PrimeBg : theme.nodeBg}
              stroke={pairColor}
              strokeWidth={is3Prime ? 2 : 1}
            />
            <text
              x={centerX - stemSpacing}
              y={y + 4}
              fill={pairColor}
              fontSize={fontSize}
              fontFamily="monospace"
              textAnchor="middle"
              fontWeight={is3Prime ? 'bold' : 'normal'}
            >
              {pair.bases[0]}
            </text>

            {/* Right base (3' strand) */}
            <circle
              cx={centerX + stemSpacing}
              cy={y}
              r={7}
              fill={is3Prime ? theme.node3PrimeBg : theme.nodeBg}
              stroke={pairColor}
              strokeWidth={is3Prime ? 2 : 1}
            />
            <text
              x={centerX + stemSpacing}
              y={y + 4}
              fill={pairColor}
              fontSize={fontSize}
              fontFamily="monospace"
              textAnchor="middle"
              fontWeight={is3Prime ? 'bold' : 'normal'}
            >
              {pair.bases[1]}
            </text>

            {/* Position labels (every 5th pair) */}
            {idx % 5 === 0 && (
              <>
                <text
                  x={centerX - stemSpacing - 15}
                  y={y + 3}
                  fill={theme.textMuted}
                  fontSize="9"
                  fontFamily="monospace"
                  textAnchor="end"
                >
                  {pair.i + 1}
                </text>
                <text
                  x={centerX + stemSpacing + 15}
                  y={y + 3}
                  fill={theme.textMuted}
                  fontSize="9"
                  fontFamily="monospace"
                  textAnchor="start"
                >
                  {pair.j + 1}
                </text>
              </>
            )}
          </g>
        );
      })}

      {/* Draw loop */}
      {loopSequence && loopLength > 0 && (
        <g>
          {/* Loop arc */}
          <path
            d={`M ${centerX - stemSpacing} ${stemTopY - 10}
                C ${centerX - stemSpacing - loopRadius/2} ${stemTopY - loopRadius}
                  ${centerX + stemSpacing + loopRadius/2} ${stemTopY - loopRadius}
                  ${centerX + stemSpacing} ${stemTopY - 10}`}
            fill="none"
            stroke={color}
            strokeWidth="2"
            opacity="0.5"
          />

          {/* Loop bases */}
          {loopSequence.split('').map((base, idx) => {
            const angle = Math.PI + (Math.PI * (idx + 0.5) / loopLength);
            const x = centerX + Math.cos(angle) * (loopRadius * 0.7);
            const loopY = stemTopY - 10 - loopRadius * 0.5 + Math.sin(angle) * (loopRadius * 0.4);

            return (
              <g key={`loop-${idx}`}>
                <circle cx={x} cy={loopY} r={6} fill={theme.nodeBg} stroke={color} strokeWidth="1" />
                <text
                  x={x}
                  y={loopY + 3}
                  fill={color}
                  fontSize="10"
                  fontFamily="monospace"
                  textAnchor="middle"
                >
                  {base}
                </text>
              </g>
            );
          })}

          {/* Loop label */}
          <text
            x={centerX}
            y={stemTopY - loopRadius - 15}
            fill={theme.textMuted}
            fontSize="10"
            fontFamily="monospace"
            textAnchor="middle"
          >
            Loop ({loopLength}nt)
          </text>
        </g>
      )}

      {/* 5' and 3' labels */}
      <text
        x={centerX - stemSpacing - 15}
        y={stemStartY + 20}
        fill={theme.prime5}
        fontSize="11"
        fontFamily="monospace"
        textAnchor="middle"
        fontWeight="bold"
      >
        5'
      </text>
      <text
        x={centerX + stemSpacing + 15}
        y={stemStartY + 20}
        fill={theme.prime3}
        fontSize="11"
        fontFamily="monospace"
        textAnchor="middle"
        fontWeight="bold"
      >
        3'
      </text>

      {/* 3' warning - severity based */}
      {has3PrimeStructure && (() => {
        // Calculate severity for this specific structure
        const basePairsArray: BasePair[] = stemPairs.map(p => [p.i, p.j]);
        const severity = (classify3PrimeStructureSeverity as any)({
          energy,
          basePairs: basePairsArray,
          seqLength: sequence.length,
          structure,
        }) as SeverityResult;

        // Only show warning in SVG for critical/warning levels
        if (severity.level === 'none' || severity.level === 'info' || severity.level === 'low') {
          return null;
        }

        const warningColors: Record<string, { bg: string; border: string; text: string }> = {
          critical: { bg: '#fee2e2', border: '#dc2626', text: '#dc2626' },
          warning: { bg: '#fef3c7', border: '#d97706', text: '#d97706' },
          moderate: { bg: '#fef9c3', border: '#ca8a04', text: '#ca8a04' },
        };
        const colors = warningColors[severity.level] || warningColors.moderate;

        const warningText = severity.level === 'critical'
          ? "🔴 3' end blocked - redesign needed"
          : severity.level === 'warning'
          ? "⚠ 3' structure risk - may affect efficiency"
          : "△ 3' structure - monitor efficiency";

        return (
          <g>
            <rect
              x={width/2 - 130}
              y={height - 35}
              width={260}
              height={28}
              rx={4}
              fill={colors.bg}
              stroke={colors.border}
              strokeWidth="1"
            />
            <text
              x={width/2}
              y={height - 17}
              fill={colors.text}
              fontSize="10"
              fontFamily="monospace"
              textAnchor="middle"
              fontWeight="bold"
            >
              {warningText}
            </text>
          </g>
        );
      })()}
    </svg>
  );
};

/**
 * Alternative flat view for complex/multi-loop structures
 * White theme version
 */
const FlatStructureDiagram: FC<FlatStructureDiagramProps> = ({ sequence, basePairs, energy, width = 600, height = 180 }) => {
  const seqLength = sequence.length;
  const padding = 40;
  const baseWidth = Math.min(18, (width - 2 * padding) / seqLength);
  const arcMaxHeight = 80;
  const threePrimeRegion = seqLength - 10;

  const color = energyToColor(energy);

  // White theme colors
  const theme = {
    bg: '#f8fafc',
    text: '#1e293b',
    textMuted: '#64748b',
    prime5: '#2563eb',
    prime3: '#dc2626',
    prime3Bg: '#dbeafe',
  };

  // Calculate positions
  const basePositions = sequence.split('').map((_, i) => ({
    x: padding + i * baseWidth + baseWidth / 2,
    y: height - 50
  }));

  // Generate arc paths
  const arcs = basePairs.map(([i, j], idx) => {
    if (i >= basePositions.length || j >= basePositions.length) return null;

    const x1 = basePositions[i].x;
    const x2 = basePositions[j].x;
    const span = Math.abs(j - i);
    const arcHeight = Math.min(arcMaxHeight, span * 3);

    const involves3Prime = i >= threePrimeRegion || j >= threePrimeRegion;

    return {
      path: `M ${x1} ${height - 50} Q ${(x1 + x2) / 2} ${height - 50 - arcHeight} ${x2} ${height - 50}`,
      involves3Prime,
      span,
      key: idx
    };
  }).filter((arc): arc is NonNullable<typeof arc> => arc !== null);

  return (
    <svg width={width} height={height} style={{ background: theme.bg, borderRadius: '6px' }}>
      {/* Energy label */}
      <text x={padding} y={20} fill={theme.text} fontSize="12" fontFamily="monospace">
        Secondary Structure (ΔG = {energy.toFixed(1)} kcal/mol)
      </text>

      {/* Severity badge */}
      <rect x={width - 110} y={8} width={95} height={20} rx={4} fill={color} opacity={0.15} />
      <text x={width - 62} y={22} fill={color} fontSize="10" fontFamily="monospace" textAnchor="middle">
        {energySeverity(energy).level}
      </text>

      {/* Draw arcs */}
      {arcs.map(arc => (
        <path
          key={arc.key}
          d={arc.path}
          fill="none"
          stroke={arc.involves3Prime ? theme.prime3 : color}
          strokeWidth={arc.involves3Prime ? 2.5 : 1.5}
          opacity={arc.involves3Prime ? 0.9 : 0.6}
        />
      ))}

      {/* Draw sequence */}
      {sequence.split('').map((base, i) => {
        const pos = basePositions[i];
        if (!pos) return null;

        const isPaired = basePairs.some(([a, b]) => a === i || b === i);
        const in3Prime = i >= threePrimeRegion;

        let baseColor = theme.textMuted;
        if (isPaired && in3Prime) {
          baseColor = theme.prime3;
        } else if (isPaired) {
          baseColor = color;
        } else if (in3Prime) {
          baseColor = theme.prime5;
        }

        return (
          <g key={i}>
            {in3Prime && (
              <rect
                x={pos.x - baseWidth/2}
                y={pos.y - 10}
                width={baseWidth}
                height={20}
                fill={theme.prime5}
                opacity={0.1}
              />
            )}
            <text
              x={pos.x}
              y={pos.y}
              fill={baseColor}
              fontSize={Math.min(12, baseWidth - 2)}
              fontFamily="monospace"
              textAnchor="middle"
              dominantBaseline="middle"
              fontWeight={isPaired ? 'bold' : 'normal'}
            >
              {base}
            </text>
            {(i + 1) % 10 === 0 && (
              <text x={pos.x} y={pos.y + 18} fill={theme.textMuted} fontSize="8" fontFamily="monospace" textAnchor="middle">
                {i + 1}
              </text>
            )}
          </g>
        );
      })}

      {/* Labels */}
      <text x={basePositions[0]?.x - 12} y={height - 50} fill={theme.prime5} fontSize="10" fontFamily="monospace" textAnchor="end">
        5'
      </text>
      <text x={basePositions[seqLength - 1]?.x + 12} y={height - 50} fill={theme.prime3} fontSize="10" fontFamily="monospace">
        3'
      </text>

      {/* Legend */}
      <g transform={`translate(${padding}, ${height - 20})`}>
        <rect x={0} y={-8} width={8} height={8} fill={theme.prime5} opacity={0.3} />
        <text x={12} y={0} fill={theme.textMuted} fontSize="9" fontFamily="monospace">3' region</text>

        <rect x={80} y={-8} width={8} height={8} fill={color} />
        <text x={92} y={0} fill={theme.textMuted} fontSize="9" fontFamily="monospace">Base paired</text>

        {arcs.some(a => a.involves3Prime) && (
          <>
            <rect x={170} y={-8} width={8} height={8} fill={theme.prime3} />
            <text x={182} y={0} fill={theme.prime3} fontSize="9" fontFamily="monospace" fontWeight="bold">
              3' in structure!
            </text>
          </>
        )}
      </g>
    </svg>
  );
};

/**
 * Main HairpinDiagram Component
 *
 * Automatically chooses the best visualization based on structure complexity:
 * - Simple hairpins: Stem-loop diagram
 * - Complex structures: Flat arc diagram
 */
const HairpinDiagram: FC<HairpinDiagramProps> = ({
  sequence,
  foldResult = null,
  primerName = 'Primer',
  width = 400,
  showDetails = true,
  temperature = 55, // Use PCR annealing temp by default
  mode = 'auto'  // 'auto', 'stemloop', 'flat'
}) => {
  const [viewMode, setViewMode] = useState<'auto' | 'stemloop' | 'flat'>(mode);

  // Compute fold if not provided (use temperature for accurate PCR analysis)
  const computedFold = useMemo((): FoldResult | null => {
    if (foldResult) return foldResult;
    if (!sequence || sequence.length < 6) return null;

    try {
      return foldSequence(sequence, temperature) as FoldResult;
    } catch (e) {
      console.error('Fold failed:', e);
      return null;
    }
  }, [sequence, foldResult, temperature]);

  // White theme colors
  const theme: Theme = {
    bg: '#ffffff',
    bgSecondary: '#f8fafc',
    bgTertiary: '#f1f5f9',
    border: '#e2e8f0',
    text: '#1e293b',
    textSecondary: '#475569',
    textMuted: '#64748b',
    success: '#16a34a',
    successBg: '#dcfce7',
    warning: '#ca8a04',
    warningBg: '#fef9c3',
    danger: '#dc2626',
    dangerBg: '#fee2e2',
    prime5: '#2563eb',
    prime3: '#dc2626',
  };

  if (!sequence) {
    return (
      <div style={{ color: theme.textMuted, padding: '20px', textAlign: 'center', background: theme.bg, borderRadius: '8px' }}>
        No sequence provided
      </div>
    );
  }

  const energy = computedFold?.e ?? 0;
  const basePairs: BasePair[] = (computedFold?.ij ?? []) as BasePair[];

  // No significant structure
  if (energy > -0.5 || basePairs.length === 0) {
    return (
      <div style={{
        background: theme.bg,
        borderRadius: '8px',
        padding: '20px',
        border: `1px solid ${theme.border}`
      }}>
        <div style={{
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          marginBottom: '15px'
        }}>
          <span style={{ color: theme.text, fontWeight: 'bold' }}>{primerName}</span>
          <span style={{ color: theme.success, fontSize: '13px' }}>
            No significant structure (ΔG = {energy.toFixed(1)} kcal/mol)
          </span>
        </div>

        <div style={{
          fontFamily: 'monospace',
          fontSize: '13px',
          color: theme.text,
          background: theme.bgTertiary,
          padding: '12px',
          borderRadius: '4px',
          letterSpacing: '1px',
          wordBreak: 'break-all'
        }}>
          5'-{sequence}-3'
        </div>

        <div style={{
          marginTop: '15px',
          padding: '12px',
          background: theme.successBg,
          borderRadius: '6px',
          border: `1px solid ${theme.success}`,
          display: 'flex',
          alignItems: 'center',
          gap: '10px'
        }}>
          <span style={{ fontSize: '20px' }}>✓</span>
          <div>
            <div style={{ color: theme.success, fontWeight: 'bold' }}>Excellent - No Hairpin Formation</div>
            <div style={{ color: theme.textSecondary, fontSize: '12px' }}>
              This primer has minimal secondary structure and should perform well
            </div>
          </div>
        </div>
      </div>
    );
  }

  // Parse structure for visualization
  const structure = parseHairpinStructure(sequence, basePairs);
  const severity = energySeverity(energy);
  const color = energyToColor(energy);

  // Determine view mode
  const effectiveMode = viewMode === 'auto'
    ? (structure && structure.stemPairs.length <= 15 && structure.stemPairs.length >= 2 ? 'stemloop' : 'flat')
    : viewMode;

  return (
    <div style={{
      background: theme.bg,
      borderRadius: '8px',
      padding: '16px',
      border: `1px solid ${energy < -4 ? theme.danger : theme.border}`
    }}>
      {/* Header */}
      <div style={{
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center',
        marginBottom: '12px'
      }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
          <span style={{ color: theme.text, fontWeight: 'bold', fontSize: '14px' }}>
            {primerName}
          </span>
          <span style={{
            padding: '3px 10px',
            borderRadius: '4px',
            fontSize: '11px',
            fontFamily: 'monospace',
            background: color + '18',
            color: color,
            fontWeight: 'bold'
          }}>
            {severity.label}
          </span>
        </div>

        {/* View mode toggle */}
        <div style={{ display: 'flex', gap: '4px' }}>
          <button
            onClick={() => setViewMode('stemloop')}
            style={{
              padding: '4px 10px',
              background: effectiveMode === 'stemloop' ? theme.prime5 : theme.bgSecondary,
              border: `1px solid ${effectiveMode === 'stemloop' ? theme.prime5 : theme.border}`,
              borderRadius: '4px',
              color: effectiveMode === 'stemloop' ? '#fff' : theme.textMuted,
              fontSize: '11px',
              cursor: 'pointer'
            }}
          >
            Stem-Loop
          </button>
          <button
            onClick={() => setViewMode('flat')}
            style={{
              padding: '4px 10px',
              background: effectiveMode === 'flat' ? theme.prime5 : theme.bgSecondary,
              border: `1px solid ${effectiveMode === 'flat' ? theme.prime5 : theme.border}`,
              borderRadius: '4px',
              color: effectiveMode === 'flat' ? '#fff' : theme.textMuted,
              fontSize: '11px',
              cursor: 'pointer'
            }}
          >
            Arc View
          </button>
        </div>
      </div>

      {/* Structure Visualization */}
      <div style={{ marginBottom: '12px', overflowX: 'auto' }}>
        {effectiveMode === 'stemloop' && structure ? (
          <StemLoopDiagram
            sequence={sequence}
            structure={structure}
            energy={energy}
            width={width}
            height={280}
          />
        ) : (
          <FlatStructureDiagram
            sequence={sequence}
            basePairs={basePairs}
            energy={energy}
            width={Math.max(width, sequence.length * 18 + 80)}
            height={180}
          />
        )}
      </div>

      {/* 3' Structure Warning - Severity-based */}
      {(() => {
        const severityResult = (classify3PrimeStructureSeverity as any)({
          energy,
          basePairs,
          seqLength: sequence.length,
          structure,
        }) as SeverityResult;

        // Only show warning for levels that warrant user attention
        if (!severityResult.shouldWarn && severityResult.level === 'none') return null;

        // Color scheme based on severity
        const severityColors: Record<string, { bg: string; border: string; text: string; icon: string }> = {
          critical: { bg: theme.dangerBg, border: theme.danger, text: theme.danger, icon: '🔴' },
          warning: { bg: '#fef3c7', border: '#f59e0b', text: '#d97706', icon: '🟠' },
          moderate: { bg: '#fef9c3', border: '#eab308', text: '#ca8a04', icon: '🟡' },
          low: { bg: '#e0f2fe', border: '#0ea5e9', text: '#0369a1', icon: 'ℹ️' },
          info: { bg: '#f0f9ff', border: '#0ea5e9', text: '#0369a1', icon: 'ℹ️' },
          none: { bg: '#dcfce7', border: '#22c55e', text: '#16a34a', icon: '✓' },
        };

        const colors = severityColors[severityResult.level] || severityColors.info;

        return (
          <div style={{
            padding: '12px',
            background: colors.bg,
            borderRadius: '6px',
            border: `1px solid ${colors.border}`,
            marginBottom: '12px'
          }}>
            <div style={{
              display: 'flex',
              alignItems: 'flex-start',
              gap: '10px'
            }}>
              <span style={{ fontSize: '20px' }}>{colors.icon}</span>
              <div style={{ flex: 1 }}>
                <div style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: '8px',
                  marginBottom: '4px'
                }}>
                  <div style={{ color: colors.text, fontWeight: 'bold', fontSize: '13px' }}>
                    {severityResult.label}
                  </div>
                  <span style={{
                    padding: '2px 6px',
                    background: colors.border,
                    color: '#fff',
                    borderRadius: '4px',
                    fontSize: '10px',
                    fontWeight: '600',
                    textTransform: 'uppercase',
                  }}>
                    {severityResult.level}
                  </span>
                </div>

                {severityResult.message && (
                  <div style={{ color: theme.text, fontSize: '12px', lineHeight: '1.5', marginBottom: '8px' }}>
                    {severityResult.message}
                  </div>
                )}

                {/* Detailed explanation for critical/warning levels */}
                {severityResult.level === 'critical' && (
                  <>
                    <div style={{ color: theme.text, fontSize: '12px', lineHeight: '1.5', marginBottom: '4px' }}>
                      The 3' end is involved in a {structure?.has3PrimeInStem ? 'stem' : 'loop'} structure:
                    </div>
                    <ul style={{ color: theme.text, fontSize: '12px', margin: '4px 0 8px 0', paddingLeft: '20px', lineHeight: '1.6' }}>
                      <li>DNA polymerase <strong>requires</strong> a free 3' OH to begin extension</li>
                      <li>Hairpin ΔG = {energy.toFixed(1)} kcal/mol (stable at annealing temp)</li>
                      <li>PCR will likely <strong>fail</strong> or have very low yield</li>
                    </ul>
                    <div style={{
                      padding: '8px 10px',
                      background: 'rgba(220, 38, 38, 0.1)',
                      borderRadius: '4px',
                      fontSize: '12px',
                      fontWeight: '600',
                      color: theme.danger
                    }}>
                      ⚠ Action Required: Redesign primer to free the 3' end
                    </div>
                  </>
                )}

                {severityResult.level === 'warning' && (
                  <>
                    <ul style={{ color: theme.text, fontSize: '12px', margin: '4px 0 8px 0', paddingLeft: '20px', lineHeight: '1.6' }}>
                      <li>Structure may compete with template binding</li>
                      <li>ΔG = {energy.toFixed(1)} kcal/mol (moderately stable)</li>
                      <li>PCR may work but with reduced efficiency</li>
                    </ul>
                    <div style={{
                      padding: '8px 10px',
                      background: 'rgba(217, 119, 6, 0.1)',
                      borderRadius: '4px',
                      fontSize: '12px',
                      color: colors.text
                    }}>
                      <strong>Tips:</strong> Adjust length by 1-2 bases, use touchdown PCR, or add 2-5% DMSO
                    </div>
                  </>
                )}

                {(severityResult.level === 'moderate' || severityResult.level === 'low' || severityResult.level === 'info') && (
                  <div style={{
                    padding: '6px 10px',
                    background: 'rgba(3, 105, 161, 0.08)',
                    borderRadius: '4px',
                    fontSize: '11px',
                    color: theme.textMuted,
                    lineHeight: '1.5'
                  }}>
                    ΔG = {energy.toFixed(1)} kcal/mol • {severityResult.level === 'info' || severityResult.level === 'low'
                      ? 'Structure is weak and typically melts at annealing temperature'
                      : 'Monitor PCR efficiency; optimize if yields are low'}
                  </div>
                )}
              </div>
            </div>
          </div>
        );
      })()}

      {/* Details panel */}
      {showDetails && (
        <div style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))',
          gap: '12px',
          padding: '12px',
          background: theme.bgSecondary,
          borderRadius: '6px',
          border: `1px solid ${theme.border}`
        }}>
          <div>
            <div style={{ color: theme.textMuted, fontSize: '11px', marginBottom: '2px' }}>Free Energy (ΔG)</div>
            <div style={{ color: color, fontSize: '16px', fontWeight: 'bold', fontFamily: 'monospace' }}>
              {energy.toFixed(1)} kcal/mol
            </div>
          </div>

          <div>
            <div style={{ color: theme.textMuted, fontSize: '11px', marginBottom: '2px' }}>Base Pairs</div>
            <div style={{ color: theme.text, fontSize: '16px', fontWeight: 'bold' }}>
              {basePairs.length}
            </div>
          </div>

          {structure?.loopSequence && (
            <div>
              <div style={{ color: theme.textMuted, fontSize: '11px', marginBottom: '2px' }}>Loop Size</div>
              <div style={{ color: theme.text, fontSize: '16px', fontWeight: 'bold' }}>
                {structure.loopSequence.length} nt
              </div>
            </div>
          )}

          {(() => {
            const severityResult = (classify3PrimeStructureSeverity as any)({
              energy,
              basePairs,
              seqLength: sequence.length,
              structure,
            }) as SeverityResult;

            const statusColors: Record<string, string> = {
              critical: theme.danger,
              warning: '#d97706',
              moderate: '#ca8a04',
              low: '#0369a1',
              info: '#0369a1',
              none: theme.success,
            };

            return (
              <div style={{ position: 'relative' }}>
                <div style={{ color: theme.textMuted, fontSize: '11px', marginBottom: '2px' }}>3' Status</div>
                <div
                  style={{
                    color: statusColors[severityResult.level] || theme.success,
                    fontSize: '14px',
                    fontWeight: 'bold',
                    cursor: 'help',
                    display: 'inline-flex',
                    alignItems: 'center',
                    gap: '4px',
                  }}
                  title={severityResult.tooltip}
                >
                  {severityResult.shortLabel}
                  <span style={{
                    width: '14px',
                    height: '14px',
                    borderRadius: '50%',
                    background: statusColors[severityResult.level] || theme.success,
                    color: '#fff',
                    fontSize: '9px',
                    display: 'inline-flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    fontWeight: 'normal',
                  }}>
                    ?
                  </span>
                </div>
              </div>
            );
          })()}
        </div>
      )}

      {/* Dot-bracket notation */}
      {showDetails && (
        <div style={{
          marginTop: '12px',
          padding: '10px',
          background: theme.bgTertiary,
          borderRadius: '4px',
          fontFamily: 'monospace',
          fontSize: '11px',
          overflowX: 'auto',
          border: `1px solid ${theme.border}`
        }}>
          <div style={{ color: theme.textMuted, marginBottom: '4px' }}>Dot-bracket notation:</div>
          <div style={{ color: theme.text, letterSpacing: '0.5px', wordBreak: 'break-all' }}>
            {sequence}
          </div>
          <div style={{ color: '#ea580c', letterSpacing: '0.5px' }}>
            {(() => {
              const dots = new Array(sequence.length).fill('.');
              basePairs.forEach(([i, j]) => {
                if (i < dots.length && j < dots.length) {
                  dots[i] = '(';
                  dots[j] = ')';
                }
              });
              return dots.join('');
            })()}
          </div>
        </div>
      )}
    </div>
  );
};

export default HairpinDiagram;

/**
 * Compact badge version for table cells
 * Updated to use severity-based classification for better user guidance
 */
export const HairpinBadge: FC<HairpinBadgeProps> = ({ sequence, foldResult = null, lightTheme = false, showTooltip = true, temperature = 55 }) => {
  const [_isHovered, setIsHovered] = useState<boolean>(false);

  const computedFold = useMemo((): FoldResult | null => {
    if (foldResult) return foldResult;
    if (!sequence || sequence.length < 6) return null;

    try {
      return foldSequence(sequence, temperature) as FoldResult;
    } catch (e) {
      return null;
    }
  }, [sequence, foldResult, temperature]);

  const energy = computedFold?.e ?? 0;
  const basePairs: BasePair[] = (computedFold?.ij ?? []) as BasePair[];
  const seqLength = sequence?.length ?? 0;

  // Use new severity classification
  const severity = useMemo((): SeverityResult => {
    return (classify3PrimeStructureSeverity as any)({
      energy,
      basePairs,
      seqLength,
      structure: null,
    }) as SeverityResult;
  }, [energy, basePairs, seqLength]);

  // Color schemes based on severity and theme
  const colorSchemes: Record<string, Record<string, { bg: string; text: string; border: string }>> = {
    light: {
      none: { bg: '#dcfce7', text: '#16a34a', border: '#bbf7d0' },
      info: { bg: '#f0f9ff', text: '#0369a1', border: '#e0f2fe' },
      low: { bg: '#e0f2fe', text: '#0284c7', border: '#bae6fd' },
      moderate: { bg: '#fef9c3', text: '#ca8a04', border: '#fef08a' },
      warning: { bg: '#fef3c7', text: '#d97706', border: '#fde68a' },
      critical: { bg: '#fee2e2', text: '#dc2626', border: '#fecaca' },
    },
    dark: {
      none: { bg: '#052e16', text: '#4ade80', border: '#166534' },
      info: { bg: '#0c4a6e', text: '#7dd3fc', border: '#075985' },
      low: { bg: '#0c4a6e', text: '#38bdf8', border: '#0369a1' },
      moderate: { bg: '#422006', text: '#fcd34d', border: '#854d0e' },
      warning: { bg: '#451a03', text: '#fbbf24', border: '#92400e' },
      critical: { bg: '#450a0a', text: '#fca5a5', border: '#991b1b' },
    },
  };

  const themeKey = lightTheme ? 'light' : 'dark';
  const colors = colorSchemes[themeKey][severity.level] || colorSchemes[themeKey].none;

  // No structure case
  if (energy > -1 && severity.level === 'none') {
    return (
      <span
        style={{
          padding: '4px 10px',
          background: colors.bg,
          color: colors.text,
          borderRadius: '6px',
          fontSize: '11px',
          fontFamily: 'monospace',
          fontWeight: '500',
          cursor: showTooltip ? 'help' : 'default',
        }}
        title={showTooltip ? severity.tooltip : undefined}
      >
        ✓ No hairpin
      </span>
    );
  }

  // Severity-based badge icons and labels
  const badgeConfig: Record<string, { icon: string; label: string }> = {
    critical: { icon: '🔴', label: "3' blocked" },
    warning: { icon: '⚠', label: "3' risk" },
    moderate: { icon: '△', label: "3' struct" },
    low: { icon: 'ℹ', label: 'Hairpin' },
    info: { icon: '⚡', label: 'Hairpin' },
    none: { icon: '⚡', label: 'Hairpin' },
  };

  const config = badgeConfig[severity.level] || badgeConfig.none;

  return (
    <span
      style={{
        position: 'relative',
        display: 'inline-block',
      }}
      onMouseEnter={() => setIsHovered(true)}
      onMouseLeave={() => setIsHovered(false)}
    >
      <span
        style={{
          padding: '4px 10px',
          background: colors.bg,
          color: colors.text,
          border: `1px solid ${colors.border}`,
          borderRadius: '6px',
          fontSize: '11px',
          fontFamily: 'monospace',
          fontWeight: severity.shouldWarn ? '600' : '500',
          cursor: showTooltip ? 'help' : 'default',
          display: 'inline-flex',
          alignItems: 'center',
          gap: '4px',
        }}
        title={showTooltip ? severity.tooltip : undefined}
      >
        {config.icon} {config.label} {energy.toFixed(1)}
      </span>
    </span>
  );
};
