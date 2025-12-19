import { useMemo } from 'react';
import { foldSequence } from '../lib/fold.js';

/**
 * SecondaryStructureViewer - Draws hairpin/secondary structures visually
 *
 * Shows the primer sequence with arc diagrams connecting base pairs,
 * making it immediately obvious when there's a problematic stem-loop
 * at the 3' end (much more intuitive than just showing "-5.4 kcal/mol")
 */

// Type definitions
type BasePair = [number, number];

interface FoldResult {
  e: number;
  ij: BasePair[];
  desc?: string;
}

interface BasePosition {
  x: number;
  y: number;
}

interface ArcData {
  path: string;
  span: number;
  involves3Prime: boolean;
  i: number;
  j: number;
  key: string;
}

type SeverityLevel = 'none' | 'low' | 'moderate' | 'high' | 'severe';

// Color scale for structure stability (more negative = more stable = worse for primers)
function energyToColor(dG: number): string {
  if (dG >= -1) return '#22c55e';      // Green - minimal structure
  if (dG >= -2) return '#84cc16';      // Light green
  if (dG >= -3) return '#eab308';      // Yellow - moderate concern
  if (dG >= -4) return '#f97316';      // Orange - problematic
  return '#ef4444';                     // Red - severe structure
}

// Get severity level for energy
function energyToSeverity(dG: number): SeverityLevel {
  if (dG >= -1) return 'none';
  if (dG >= -2) return 'low';
  if (dG >= -3) return 'moderate';
  if (dG >= -4) return 'high';
  return 'severe';
}

// Check if a position is in the 3' critical region (last 10bp)
function isIn3PrimeRegion(pos: number, seqLength: number, regionSize: number = 10): boolean {
  return pos >= seqLength - regionSize;
}

interface ArcDiagramProps {
  sequence: string;
  basePairs: BasePair[];
  energy: number;
  width?: number;
  height?: number;
}

/**
 * Arc diagram component - draws arcs connecting base pairs
 */
function ArcDiagram({ sequence, basePairs, energy, width = 600, height = 200 }: ArcDiagramProps) {
  const seqLength = sequence.length;
  const padding = 40;
  const baseWidth = Math.min(20, (width - 2 * padding) / seqLength);

  // Calculate positions for each nucleotide
  const basePositions = useMemo((): BasePosition[] => {
    return sequence.split('').map((_, i) => ({
      x: padding + i * baseWidth + baseWidth / 2,
      y: height - 40
    }));
  }, [sequence, baseWidth, padding, height]);

  // Generate arc paths for base pairs
  const arcs = useMemo((): ArcData[] => {
    if (!basePairs || basePairs.length === 0) return [];

    return basePairs.map(([i, j], idx) => {
      const x1 = basePositions[i]?.x;
      const x2 = basePositions[j]?.x;

      if (x1 === undefined || x2 === undefined) return null;

      const arcRadius = (Math.abs(x2 - x1) / 2);

      // Determine if this base pair involves the 3' region
      const involves3Prime = isIn3PrimeRegion(i, seqLength) || isIn3PrimeRegion(j, seqLength);

      return {
        path: `M ${x1} ${height - 40} A ${arcRadius} ${arcRadius * 0.8} 0 0 1 ${x2} ${height - 40}`,
        span: Math.abs(j - i),
        involves3Prime,
        i,
        j,
        key: `arc-${idx}`
      };
    }).filter((arc): arc is ArcData => arc !== null);
  }, [basePairs, basePositions, height, seqLength]);

  const color = energyToColor(energy);
  const severity = energyToSeverity(energy);

  return (
    <svg width={width} height={height} className="bg-[#1a1a2e]">
      {/* Title and energy label */}
      <text x={padding} y={20} fill="#94a3b8" fontSize="12" fontFamily="monospace">
        Secondary Structure (dG = {energy.toFixed(1)} kcal/mol)
      </text>

      {/* Severity indicator */}
      <rect x={width - 100} y={8} width={80} height={20} rx={4} fill={color} opacity={0.3} />
      <text x={width - 60} y={22} fill={color} fontSize="11" fontFamily="monospace" textAnchor="middle">
        {severity}
      </text>

      {/* Draw arcs for base pairs */}
      {arcs.map(arc => (
        <path
          key={arc.key}
          d={arc.path}
          fill="none"
          stroke={arc.involves3Prime ? '#ef4444' : color}
          strokeWidth={arc.involves3Prime ? 2.5 : 1.5}
          opacity={arc.involves3Prime ? 0.9 : 0.6}
        />
      ))}

      {/* Draw sequence */}
      {sequence.split('').map((base, i) => {
        const pos = basePositions[i];
        const isPaired = basePairs?.some(([a, b]) => a === i || b === i);
        const in3Prime = isIn3PrimeRegion(i, seqLength);

        // Color code: paired bases in structure color, 3' region highlighted
        let baseColor = '#64748b';  // Default gray
        if (isPaired && in3Prime) {
          baseColor = '#ef4444';    // Red for paired 3' bases (dangerous!)
        } else if (isPaired) {
          baseColor = color;        // Structure color for paired bases
        } else if (in3Prime) {
          baseColor = '#60a5fa';    // Blue for unpaired 3' bases
        }

        return (
          <g key={`base-${i}`}>
            {/* Background for 3' region */}
            {in3Prime && (
              <rect
                x={pos.x - baseWidth / 2}
                y={pos.y - 12}
                width={baseWidth}
                height={24}
                fill="#3b82f6"
                opacity={0.1}
              />
            )}

            {/* Nucleotide letter */}
            <text
              x={pos.x}
              y={pos.y}
              fill={baseColor}
              fontSize={Math.min(14, baseWidth - 2)}
              fontFamily="monospace"
              textAnchor="middle"
              dominantBaseline="middle"
              fontWeight={isPaired ? 'bold' : 'normal'}
            >
              {base}
            </text>

            {/* Position number every 5bp */}
            {(i + 1) % 5 === 0 && (
              <text
                x={pos.x}
                y={pos.y + 20}
                fill="#475569"
                fontSize="9"
                fontFamily="monospace"
                textAnchor="middle"
              >
                {i + 1}
              </text>
            )}
          </g>
        );
      })}

      {/* 3' label */}
      <text
        x={basePositions[seqLength - 1]?.x + 15 || width - padding}
        y={height - 40}
        fill="#60a5fa"
        fontSize="11"
        fontFamily="monospace"
        fontWeight="bold"
      >
        3'
      </text>

      {/* 5' label */}
      <text
        x={basePositions[0]?.x - 15 || padding}
        y={height - 40}
        fill="#94a3b8"
        fontSize="11"
        fontFamily="monospace"
        textAnchor="end"
      >
        5'
      </text>

      {/* Legend */}
      <g transform={`translate(${padding}, ${height - 15})`}>
        <rect x={0} y={-8} width={10} height={10} fill="#3b82f6" opacity={0.3} />
        <text x={15} y={0} fill="#64748b" fontSize="9" fontFamily="monospace">3' region (critical)</text>

        <rect x={120} y={-8} width={10} height={10} fill={color} />
        <text x={135} y={0} fill="#64748b" fontSize="9" fontFamily="monospace">Base paired</text>

        {arcs.some(a => a.involves3Prime) && (
          <>
            <rect x={220} y={-8} width={10} height={10} fill="#ef4444" />
            <text x={235} y={0} fill="#ef4444" fontSize="9" fontFamily="monospace" fontWeight="bold">
              3' structure (problematic!)
            </text>
          </>
        )}
      </g>
    </svg>
  );
}

interface BracketNotationProps {
  sequence: string;
  basePairs: BasePair[];
}

/**
 * Bracket notation display (dot-bracket format)
 */
function BracketNotation({ sequence, basePairs }: BracketNotationProps) {
  const notation = useMemo(() => {
    const dots = new Array(sequence.length).fill('.');
    if (basePairs) {
      basePairs.forEach(([i, j]) => {
        if (i < dots.length && j < dots.length) {
          dots[i] = '(';
          dots[j] = ')';
        }
      });
    }
    return dots.join('');
  }, [sequence, basePairs]);

  return (
    <div className="font-mono text-xs bg-slate-950 px-3 py-2 rounded overflow-x-auto">
      <div className="text-slate-400 mb-1">Dot-bracket notation:</div>
      <div className="text-slate-200 tracking-wide">{sequence}</div>
      <div className="text-orange-500 tracking-wide">{notation}</div>
    </div>
  );
}

type StructureType = 'HAIRPIN' | 'INTERNAL_LOOP' | 'BULGE' | 'STACK' | 'MULTI_LOOP';

interface StructureTypeIndicatorProps {
  type: StructureType | string;
  energy: number;
}

/**
 * Structure type indicator with icon
 */
function StructureTypeIndicator({ type, energy }: StructureTypeIndicatorProps) {
  const icons: Record<StructureType, string> = {
    'HAIRPIN': '🔄',
    'INTERNAL_LOOP': '🔃',
    'BULGE': '📎',
    'STACK': '📚',
    'MULTI_LOOP': '🌀'
  };

  const descriptions: Record<StructureType, string> = {
    'HAIRPIN': 'Hairpin loop - stem with a loop at the end',
    'INTERNAL_LOOP': 'Internal loop - asymmetric bulge in stem',
    'BULGE': 'Bulge loop - extra bases on one side of stem',
    'STACK': 'Stacked base pairs',
    'MULTI_LOOP': 'Multi-branch loop - complex junction'
  };

  return (
    <div className="flex items-center gap-2 px-3 py-2 bg-slate-800 rounded-md border-l-[3px]"
         style={{ borderLeftColor: energyToColor(energy) }}>
      <span className="text-xl">{icons[type as StructureType] || '🧬'}</span>
      <div>
        <div className="text-slate-200 font-bold text-[13px]">
          {type?.replace(/_/g, ' ') || 'Unknown Structure'}
        </div>
        <div className="text-slate-400 text-[11px]">
          {descriptions[type as StructureType] || 'Secondary structure element'}
        </div>
      </div>
    </div>
  );
}

interface StructureWarningProps {
  sequence: string;
  basePairs: BasePair[];
  energy: number;
}

/**
 * Warning panel for 3' structure issues
 */
function StructureWarning({ sequence, basePairs }: StructureWarningProps) {
  const seqLength = sequence.length;

  // Check for base pairs involving the 3' region
  const threePrimePairs = useMemo(() => {
    if (!basePairs) return [];
    return basePairs.filter(([i, j]) =>
      isIn3PrimeRegion(i, seqLength) || isIn3PrimeRegion(j, seqLength)
    );
  }, [basePairs, seqLength]);

  if (threePrimePairs.length === 0) {
    return (
      <div className="p-3 bg-green-950 rounded-md border border-green-700 flex items-center gap-2.5">
        <span className="text-xl">✅</span>
        <div>
          <div className="text-green-400 font-bold">3' End Clear</div>
          <div className="text-green-300 text-xs">
            No secondary structure in the critical 3' annealing region
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="p-3 bg-red-950 rounded-md border border-red-600 flex items-center gap-2.5">
      <span className="text-2xl">⚠️</span>
      <div>
        <div className="text-red-300 font-bold">
          3' End Structure Detected!
        </div>
        <div className="text-red-200 text-xs">
          {threePrimePairs.length} base pair(s) involve the 3' annealing region.
          This can severely reduce primer efficiency by preventing proper template binding.
        </div>
        <div className="text-red-400 text-[11px] mt-1">
          Consider repositioning the mutation site or redesigning with different flanking sequences.
        </div>
      </div>
    </div>
  );
}

interface SecondaryStructureViewerProps {
  sequence: string;
  foldResult?: FoldResult | null;
  primerName?: string;
  width?: number;
  showBracketNotation?: boolean;
  showStructureType?: boolean;
  temperature?: number;
}

/**
 * Main SecondaryStructureViewer component
 */
export default function SecondaryStructureViewer({
  sequence,
  foldResult = null,
  primerName = 'Primer',
  width = 600,
  showBracketNotation = true,
  showStructureType = true,
  temperature = 55, // Use PCR annealing temp by default
}: SecondaryStructureViewerProps) {
  // Compute fold if not provided (use temperature for accurate PCR analysis)
  const computedFold = useMemo((): FoldResult | null => {
    if (foldResult) return foldResult;
    if (!sequence || sequence.length < 6) return null;

    try {
      return foldSequence(sequence, temperature) as FoldResult;
    } catch (e) {
      console.error('Fold computation failed:', e);
      return null;
    }
  }, [sequence, foldResult, temperature]);

  if (!sequence) {
    return (
      <div className="text-slate-400 p-5 text-center">
        No sequence provided
      </div>
    );
  }

  const energy = computedFold?.e ?? 0;
  const basePairs = computedFold?.ij ?? [];
  const structureType = computedFold?.desc;

  // No significant structure
  if (energy > -0.5 || basePairs.length === 0) {
    return (
      <div className="bg-slate-950 rounded-lg p-4 border border-slate-800">
        <div className="text-slate-400 text-[13px] mb-3 flex justify-between items-center">
          <span className="font-bold">{primerName}</span>
          <span className="text-green-500">No significant structure (dG = {energy.toFixed(1)} kcal/mol)</span>
        </div>

        <div className="font-mono text-sm text-slate-200 bg-[#1a1a2e] p-3 rounded tracking-wide">
          5'-{sequence}-3'
        </div>

        <div className="mt-3 p-2.5 bg-green-950 rounded-md border border-green-700 text-green-400 text-xs flex items-center gap-2">
          <span>✅</span>
          <span>Minimal secondary structure - good primer design!</span>
        </div>
      </div>
    );
  }

  return (
    <div className="bg-slate-950 rounded-lg p-4 border border-slate-800">
      {/* Header */}
      <div className="flex justify-between items-center mb-3">
        <span className="text-slate-200 font-bold text-sm">
          {primerName}
        </span>
        <span className="text-xs font-mono" style={{ color: energyToColor(energy) }}>
          {sequence.length}bp
        </span>
      </div>

      {/* Structure type indicator */}
      {showStructureType && structureType && (
        <div className="mb-3">
          <StructureTypeIndicator type={structureType} energy={energy} />
        </div>
      )}

      {/* Arc diagram visualization */}
      <div className="mb-3 overflow-x-auto">
        <ArcDiagram
          sequence={sequence}
          basePairs={basePairs}
          energy={energy}
          width={Math.max(width, sequence.length * 18 + 80)}
          height={160}
        />
      </div>

      {/* Warning for 3' structure */}
      <div className="mb-3">
        <StructureWarning
          sequence={sequence}
          basePairs={basePairs}
          energy={energy}
        />
      </div>

      {/* Bracket notation */}
      {showBracketNotation && (
        <div className="mb-3">
          <BracketNotation sequence={sequence} basePairs={basePairs} />
        </div>
      )}

      {/* Base pair details */}
      {basePairs.length > 0 && (
        <div className="text-[11px] text-slate-500 font-mono p-2 bg-slate-800 rounded">
          <div className="mb-1 text-slate-400">
            Base pairs ({basePairs.length}):
          </div>
          <div className="flex flex-wrap gap-1.5">
            {basePairs.map(([i, j], idx) => {
              const in3Prime = isIn3PrimeRegion(i, sequence.length) ||
                               isIn3PrimeRegion(j, sequence.length);
              return (
                <span
                  key={idx}
                  className={`px-1.5 py-0.5 rounded-sm ${
                    in3Prime
                      ? 'bg-red-900 text-red-300'
                      : 'bg-slate-700 text-slate-300'
                  }`}
                >
                  {sequence[i]}{i + 1}:{sequence[j]}{j + 1}
                </span>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
}

interface SecondaryStructureBadgeProps {
  sequence: string;
  foldResult?: FoldResult | null;
  temperature?: number;
}

/**
 * Compact version for inline display
 */
export function SecondaryStructureBadge({ sequence, foldResult = null, temperature = 55 }: SecondaryStructureBadgeProps) {
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
  const basePairs = computedFold?.ij ?? [];
  const seqLength = sequence?.length ?? 0;

  // Check for 3' involvement
  const has3PrimeStructure = basePairs.some(([i, j]) =>
    isIn3PrimeRegion(i, seqLength) || isIn3PrimeRegion(j, seqLength)
  );

  if (energy > -1) {
    return (
      <span className="px-2 py-0.5 bg-green-950 text-green-400 rounded text-[11px] font-mono">
        ✓ No structure
      </span>
    );
  }

  return (
    <span className={`px-2 py-0.5 rounded text-[11px] font-mono ${
      has3PrimeStructure
        ? 'bg-red-950 text-red-300'
        : 'bg-yellow-950 text-yellow-300'
    }`}>
      {has3PrimeStructure ? '⚠️ 3\' structure!' : '⚡'} {energy.toFixed(1)} kcal/mol
    </span>
  );
}

interface MiniStructurePreviewProps {
  sequence: string;
  foldResult?: FoldResult | null;
  width?: number;
  height?: number;
  temperature?: number;
}

/**
 * Mini arc preview for table cells
 */
export function MiniStructurePreview({ sequence, foldResult = null, width = 100, height = 30, temperature = 55 }: MiniStructurePreviewProps) {
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
  const basePairs = computedFold?.ij ?? [];

  if (energy > -1 || basePairs.length === 0) {
    return (
      <svg width={width} height={height}>
        <line x1={5} y1={height - 5} x2={width - 5} y2={height - 5}
              stroke="#334155" strokeWidth={2} />
        <text x={width / 2} y={height / 2} fill="#4ade80" fontSize="8"
              textAnchor="middle" fontFamily="monospace">OK</text>
      </svg>
    );
  }

  const seqLen = sequence.length;
  const scale = (width - 10) / seqLen;
  const color = energyToColor(energy);

  return (
    <svg width={width} height={height}>
      {/* Base line */}
      <line x1={5} y1={height - 5} x2={width - 5} y2={height - 5}
            stroke="#334155" strokeWidth={1} />

      {/* Arcs */}
      {basePairs.map(([i, j], idx) => {
        const x1 = 5 + i * scale;
        const x2 = 5 + j * scale;
        const arcR = Math.abs(x2 - x1) / 2;
        const in3Prime = isIn3PrimeRegion(i, seqLen) || isIn3PrimeRegion(j, seqLen);

        return (
          <path
            key={idx}
            d={`M ${x1} ${height - 5} A ${arcR} ${arcR * 0.7} 0 0 1 ${x2} ${height - 5}`}
            fill="none"
            stroke={in3Prime ? '#ef4444' : color}
            strokeWidth={in3Prime ? 1.5 : 1}
            opacity={0.8}
          />
        );
      })}
    </svg>
  );
}
