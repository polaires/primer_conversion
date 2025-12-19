/**
 * IsothermalAssemblyPanel - Professional Primer Designer UI
 *
 * Professional interface for Gibson/NEBuilder HiFi assembly primer design
 * featuring range-based optimization with position flexibility.
 *
 * Key Features:
 * - Range-based scoring visualization
 * - Position flexibility display (junction sliding)
 * - Cross-junction secondary structure warnings
 * - Self-complementarity detection
 * - Quality tier breakdown with detailed metrics
 * - Alternative designs with selection reasons
 * - Export to IDT/CSV/FASTA formats
 * - Unified construct visualization (circular/linear)
 * - Assembly quality gauge
 * - Junction compatibility matrix
 */

import React, { useState, useMemo, useCallback } from 'react';
import {
  findOptimalOverlapRangeBased,
  designAssemblyPrimersOptimized,
  designAssembly,
  ASSEMBLY_METHODS,
  DEFAULT_ASSEMBLY_CONFIG,
} from '../lib/assemblyCore.js';

// ============================================================================
// TYPE DEFINITIONS
// ============================================================================

type QualityTier = 'excellent' | 'good' | 'acceptable' | 'marginal';
type Direction = 'Forward' | 'Reverse';
type ExportFormat = 'idt' | 'csv' | 'fasta';
type ViewMode = 'circular' | 'linear';
type PartType = 'promoter' | 'rbs' | 'cds' | 'terminator' | 'other';
type WarningSeverity = 'high' | 'medium' | 'low' | 'info';

interface PartTypeDefinition {
  name: string;
  color: string;
  icon: string;
}

interface QualityColorScheme {
  bg: string;
  text: string;
  border: string;
  badge: string;
}

interface PrimerInfo {
  name: string;
  sequence: string;
  direction: Direction;
  fragmentId: string;
  tm?: number;
  gc?: number | string;
  length?: number;
}

interface ScoreBreakdown {
  tmInRange?: number;
  gcInRange?: number;
  lengthInRange?: number;
  gcClamp?: number;
  terminal3DG?: number;
  hairpin?: number;
  patternAvoidance?: number;
  selfComplementarity?: number;
  crossJunctionHairpin?: number;
  balancedGC?: number;
  homopolymer?: number;
  threePrimeComp?: number;
  gQuadruplex?: number;
}

interface ScoreWeights {
  tmInRange?: number;
  gcInRange?: number;
  lengthInRange?: number;
  gcClamp?: number;
  terminal3DG?: number;
  hairpin?: number;
  patternAvoidance?: number;
  selfComplementarity?: number;
  crossJunctionHairpin?: number;
  balancedGC?: number;
  homopolymer?: number;
  threePrimeComp?: number;
  gQuadruplex?: number;
}

interface Primer {
  sequence?: string;
  tm?: number;
  tmFull?: number;
  gcPercent?: number | string;
  length?: number;
  hairpinDG?: number;
  qualityTier?: QualityTier;
  qualityScore?: number;
  homologyTail?: string;
  annealingRegion?: string;
  scoreBreakdown?: ScoreBreakdown;
}

interface PrimerPair {
  forward?: Primer;
  reverse?: Primer;
  pair?: {
    qualityTier?: QualityTier;
    pairScore?: number;
    annealingTemp?: number;
    tmDifference?: number;
    heterodimerDG?: number;
  };
  alternatives?: Array<{
    selectionReason?: string;
    pairScore?: number;
    compositeScore?: number;
  }>;
  warnings?: string[];
}

interface Overlap {
  sequence?: string;
  length?: number;
  tm?: number;
  gcPercent?: number | string;
  terminal3DG?: number | string;
  qualityTier?: QualityTier;
  compositeScore?: number;
  offset?: number;
  warnings?: string[];
  scores?: ScoreBreakdown;
}

interface Junction {
  overlap?: Overlap;
  alternatives?: Overlap[];
  optimizationSummary?: OptimizationSummary;
}

interface Fragment {
  id?: string;
  sequence?: string;
  length?: number;
  originalLength?: number;
  type?: PartType;
  primers?: PrimerPair;
  startAngle?: number;
  sweepAngle?: number;
  midAngle?: number;
  index?: number;
  displayLength?: number;
}

interface AssemblyQuality {
  tier?: QualityTier;
  averageScore?: number;
  minimumScore?: number;
  weakestJunction?: number;
}

interface OptimizationSummary {
  totalEvaluated?: number;
  validCount?: number;
  bestScore?: number;
  scoreDistribution?: {
    excellent: number;
    good: number;
    acceptable: number;
    marginal: number;
  };
}

interface AssemblyResult {
  fragments?: Fragment[];
  junctions?: Junction[];
  quality?: AssemblyQuality;
  warnings?: string[];
  optimization?: {
    description?: string;
  };
  error?: string;
  assembly?: {
    totalLength?: number;
  };
}

interface AssemblyConfig {
  method?: string;
  autoOptimize?: boolean;
  [key: string]: unknown;
}

// ============================================================================
// COMPONENT PROP INTERFACES
// ============================================================================

interface ExportPanelProps {
  fragments?: Fragment[];
  junctions?: Junction[];
  variant?: string;
  onExport?: (format: string, content: string) => void;
}

interface AssemblyQualityGaugeProps {
  quality?: AssemblyQuality;
  junctionCount?: number;
  fragmentCount?: number;
}

interface ConstructVisualizationProps {
  fragments?: Fragment[];
  junctions?: Junction[];
  variant?: string;
  totalLength?: number;
}

interface JunctionCompatibilityMatrixProps {
  junctions?: Junction[];
}

interface EnhancedPrimerCardProps {
  primer?: Primer;
  direction: Direction;
  fragmentId: string;
  onCopy?: (sequence: string) => void;
}

interface ScoreBarProps {
  score: number;
  label: string;
  max?: number;
  showValue?: boolean;
}

interface QualityBadgeProps {
  tier?: QualityTier;
  score?: number;
}

interface ScoreBreakdownPanelProps {
  scores: ScoreBreakdown;
  weights: ScoreWeights;
  title?: string;
}

interface WarningItemProps {
  warning: string;
  severity?: WarningSeverity;
}

interface OverlapCardProps {
  overlap: Overlap;
  index: number;
  isSelected?: boolean;
  onSelect?: (overlap: Overlap) => void;
  showDetails?: boolean;
  isClickable?: boolean;
}

interface PairQualityPanelProps {
  pair?: {
    qualityTier?: QualityTier;
    pairScore?: number;
    annealingTemp?: number;
    tmDifference?: number;
    heterodimerDG?: number;
  };
}

interface AlternativesPanelProps {
  alternatives?: Array<{
    selectionReason?: string;
    pairScore?: number;
    compositeScore?: number;
  }>;
  onSelect?: (alt: unknown) => void;
}

interface OptimizationSummaryProps {
  summary?: OptimizationSummary;
}

interface IsothermalAssemblyPanelProps {
  fragments?: Fragment[];
  variant?: string;
  config?: Partial<AssemblyConfig>;
  onPrimersCopied?: (text: string) => void;
}

// Part type definitions with colors (matching Golden Gate)
const PART_TYPES: Record<PartType, PartTypeDefinition> = {
  promoter: { name: 'Promoter', color: '#ef4444', icon: 'P' },
  rbs: { name: 'RBS', color: '#f97316', icon: 'R' },
  cds: { name: 'CDS', color: '#22c55e', icon: 'C' },
  terminator: { name: 'Terminator', color: '#3b82f6', icon: 'T' },
  other: { name: 'Other', color: '#8b5cf6', icon: 'O' },
};

// Quality tier color mapping
const QUALITY_COLORS: Record<QualityTier, QualityColorScheme> = {
  excellent: { bg: '#dcfce7', text: '#166534', border: '#86efac', badge: '#22c55e' },
  good: { bg: '#dbeafe', text: '#1e40af', border: '#93c5fd', badge: '#3b82f6' },
  acceptable: { bg: '#fef3c7', text: '#92400e', border: '#fcd34d', badge: '#f59e0b' },
  marginal: { bg: '#fee2e2', text: '#991b1b', border: '#fca5a5', badge: '#ef4444' },
};

// Score color gradient
function getScoreColor(score: number): string {
  if (score >= 80) return '#22c55e';
  if (score >= 65) return '#84cc16';
  if (score >= 50) return '#f59e0b';
  if (score >= 35) return '#f97316';
  return '#ef4444';
}

// Get quality status text
function getQualityStatus(score: number): string {
  if (score >= 80) return 'Excellent';
  if (score >= 65) return 'Good';
  if (score >= 50) return 'Acceptable';
  return 'Marginal';
}

// ============================================================================
// EXPORT PANEL - IDT, CSV, FASTA formats
// ============================================================================
function ExportPanel({ fragments, junctions, variant, onExport }: ExportPanelProps): React.ReactElement {
  const [copiedFormat, setCopiedFormat] = useState<ExportFormat | null>(null);
  const [showExportMenu, setShowExportMenu] = useState<boolean>(false);

  // Generate all primers list
  const allPrimers = useMemo<PrimerInfo[]>(() => {
    if (!fragments) return [];
    return fragments.flatMap((frag, i) => [
      {
        name: `${frag.id || `Fragment_${i + 1}`}_F`,
        sequence: frag.primers?.forward?.sequence || '',
        direction: 'Forward' as Direction,
        fragmentId: frag.id || `Fragment ${i + 1}`,
        tm: frag.primers?.forward?.tm,
        gc: frag.primers?.forward?.gcPercent,
        length: frag.primers?.forward?.length,
      },
      {
        name: `${frag.id || `Fragment_${i + 1}`}_R`,
        sequence: frag.primers?.reverse?.sequence || '',
        direction: 'Reverse' as Direction,
        fragmentId: frag.id || `Fragment ${i + 1}`,
        tm: frag.primers?.reverse?.tm,
        gc: frag.primers?.reverse?.gcPercent,
        length: frag.primers?.reverse?.length,
      },
    ]);
  }, [fragments]);

  // IDT format (tab-separated for bulk order)
  const generateIDTFormat = useCallback((): string => {
    const lines = allPrimers.map(p =>
      `${p.name}\t${p.sequence}\t25nm\tSTD`
    );
    return lines.join('\n');
  }, [allPrimers]);

  // CSV format with full details
  const generateCSVFormat = useCallback((): string => {
    const headers = ['Name', 'Sequence', 'Direction', 'Fragment', 'Length (bp)', 'Tm (°C)', 'GC (%)'];
    const rows = allPrimers.map(p => [
      p.name,
      p.sequence,
      p.direction,
      p.fragmentId,
      p.length || '',
      p.tm ? p.tm.toFixed(1) : '',
      p.gc || '',
    ]);
    return [headers.join(','), ...rows.map(r => r.join(','))].join('\n');
  }, [allPrimers]);

  // FASTA format
  const generateFASTAFormat = useCallback((): string => {
    return allPrimers.map(p =>
      `>${p.name} ${p.fragmentId} ${p.direction} Tm=${p.tm?.toFixed(1) || '-'}C\n${p.sequence}`
    ).join('\n');
  }, [allPrimers]);

  const handleExport = async (format: ExportFormat): Promise<void> => {
    let content = '';
    let filename = '';
    let mimeType = 'text/plain';

    switch (format) {
      case 'idt':
        content = generateIDTFormat();
        filename = `isothermal_primers_idt.txt`;
        break;
      case 'csv':
        content = generateCSVFormat();
        filename = `isothermal_primers.csv`;
        mimeType = 'text/csv';
        break;
      case 'fasta':
        content = generateFASTAFormat();
        filename = `isothermal_primers.fasta`;
        break;
      default:
        return;
    }

    // Copy to clipboard
    try {
      await navigator.clipboard.writeText(content);
      setCopiedFormat(format);
      setTimeout(() => setCopiedFormat(null), 2000);
    } catch (err) {
      console.error('Failed to copy:', err);
    }

    onExport?.(format, content);
  };

  const handleDownload = (format: ExportFormat): void => {
    let content = '';
    let filename = '';
    let mimeType = 'text/plain';

    switch (format) {
      case 'idt':
        content = generateIDTFormat();
        filename = `isothermal_primers_idt.txt`;
        break;
      case 'csv':
        content = generateCSVFormat();
        filename = `isothermal_primers.csv`;
        mimeType = 'text/csv';
        break;
      case 'fasta':
        content = generateFASTAFormat();
        filename = `isothermal_primers.fasta`;
        break;
      default:
        return;
    }

    // Download file
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  return (
    <div className="iso-export-panel">
      <div className="iso-export-header">
        <h5>
          <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
            <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
          </svg>
          Export Primers ({allPrimers.length})
        </h5>
      </div>
      <div className="iso-export-buttons">
        <button
          className={`iso-export-btn idt ${copiedFormat === 'idt' ? 'copied' : ''}`}
          onClick={() => handleExport('idt')}
          title="Copy in IDT bulk order format (TSV)"
        >
          <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
            <path d="M16 1H4c-1.1 0-2 .9-2 2v14h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z"/>
          </svg>
          {copiedFormat === 'idt' ? 'Copied!' : 'IDT Format'}
        </button>
        <button
          className={`iso-export-btn csv ${copiedFormat === 'csv' ? 'copied' : ''}`}
          onClick={() => handleExport('csv')}
          title="Copy as CSV spreadsheet"
        >
          <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
            <path d="M14 2H6c-1.1 0-1.99.9-1.99 2L4 20c0 1.1.89 2 1.99 2H18c1.1 0 2-.9 2-2V8l-6-6zm2 16H8v-2h8v2zm0-4H8v-2h8v2zm-3-5V3.5L18.5 9H13z"/>
          </svg>
          {copiedFormat === 'csv' ? 'Copied!' : 'CSV'}
        </button>
        <button
          className={`iso-export-btn fasta ${copiedFormat === 'fasta' ? 'copied' : ''}`}
          onClick={() => handleExport('fasta')}
          title="Copy as FASTA format"
        >
          <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
            <path d="M14 2H6c-1.1 0-1.99.9-1.99 2L4 20c0 1.1.89 2 1.99 2H18c1.1 0 2-.9 2-2V8l-6-6zM6 20V4h7v5h5v11H6z"/>
          </svg>
          {copiedFormat === 'fasta' ? 'Copied!' : 'FASTA'}
        </button>
        <div className="iso-export-dropdown">
          <button
            className="iso-export-btn download"
            onClick={() => setShowExportMenu(!showExportMenu)}
            title="Download file"
          >
            <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
              <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
            </svg>
            Download
            <svg viewBox="0 0 24 24" width="12" height="12" fill="currentColor" className="ml-1">
              <path d="M7.41 8.59L12 13.17l4.59-4.58L18 10l-6 6-6-6 1.41-1.41z"/>
            </svg>
          </button>
          {showExportMenu && (
            <div className="iso-export-menu">
              <button onClick={() => { handleDownload('idt'); setShowExportMenu(false); }}>
                IDT Format (.txt)
              </button>
              <button onClick={() => { handleDownload('csv'); setShowExportMenu(false); }}>
                CSV Spreadsheet (.csv)
              </button>
              <button onClick={() => { handleDownload('fasta'); setShowExportMenu(false); }}>
                FASTA (.fasta)
              </button>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

// ============================================================================
// ASSEMBLY QUALITY GAUGE - Visual circular progress indicator
// ============================================================================
function AssemblyQualityGauge({ quality, junctionCount, fragmentCount }: AssemblyQualityGaugeProps): React.ReactElement {
  const score = quality?.averageScore ?? 0;
  const minScore = quality?.minimumScore ?? 0;
  const tier = quality?.tier || 'marginal';

  const circumference = 2 * Math.PI * 40;
  const progress = (score / 100) * circumference;
  const color = getScoreColor(score);
  const status = getQualityStatus(score);

  return (
    <div className="iso-quality-gauge">
      <div className="iso-gauge-main">
        <div className="iso-gauge-circle">
          <svg viewBox="0 0 100 100">
            {/* Background circle */}
            <circle
              cx="50" cy="50" r="40"
              fill="none"
              stroke="#e5e7eb"
              strokeWidth="8"
            />
            {/* Progress circle */}
            <circle
              cx="50" cy="50" r="40"
              fill="none"
              stroke={color}
              strokeWidth="8"
              strokeLinecap="round"
              strokeDasharray={`${progress} ${circumference}`}
              transform="rotate(-90 50 50)"
              style={{ transition: 'stroke-dasharray 0.5s ease' }}
            />
          </svg>
          <div className="iso-gauge-center">
            <span className="iso-gauge-score">{score.toFixed(0)}</span>
            <span className="iso-gauge-unit">%</span>
          </div>
        </div>
        <div className="iso-gauge-status" style={{ color }}>
          {status}
        </div>
      </div>
      <div className="iso-gauge-details">
        <div className="iso-gauge-stat">
          <span className="iso-gauge-stat-value">{fragmentCount}</span>
          <span className="iso-gauge-stat-label">Fragments</span>
        </div>
        <div className="iso-gauge-stat">
          <span className="iso-gauge-stat-value">{junctionCount}</span>
          <span className="iso-gauge-stat-label">Junctions</span>
        </div>
        <div className="iso-gauge-stat">
          <span className="iso-gauge-stat-value" style={{ color: getScoreColor(minScore) }}>
            {minScore.toFixed(0)}
          </span>
          <span className="iso-gauge-stat-label">Min Score</span>
        </div>
      </div>
    </div>
  );
}

// ============================================================================
// CONSTRUCT VISUALIZATION - Unified circular/linear view
// ============================================================================
function ConstructVisualization({ fragments, junctions, variant, totalLength }: ConstructVisualizationProps): React.ReactElement {
  const [viewMode, setViewMode] = useState<ViewMode>('circular');
  const [hoveredFragment, setHoveredFragment] = useState<number | null>(null);
  const [hoveredJunction, setHoveredJunction] = useState<number | null>(null);

  // Calculate angles for circular view
  const fragmentAngles = useMemo<Fragment[]>(() => {
    if (!fragments?.length) return [];
    const total = fragments.reduce((sum, f) => sum + (f.originalLength || f.length || 100), 0);
    let currentAngle = -90;

    return fragments.map((frag, i) => {
      const length = frag.originalLength || frag.length || 100;
      const sweepAngle = (length / total) * 360;
      const startAngle = currentAngle;
      currentAngle += sweepAngle;
      return {
        ...frag,
        startAngle,
        sweepAngle,
        midAngle: startAngle + sweepAngle / 2,
        index: i,
        displayLength: length,
      };
    });
  }, [fragments]);

  // Generate arc path for SVG
  const describeArc = (cx: number, cy: number, r: number, startAngle: number, endAngle: number): string => {
    const start = {
      x: cx + r * Math.cos(Math.PI * startAngle / 180),
      y: cy + r * Math.sin(Math.PI * startAngle / 180),
    };
    const end = {
      x: cx + r * Math.cos(Math.PI * endAngle / 180),
      y: cy + r * Math.sin(Math.PI * endAngle / 180),
    };
    const largeArc = endAngle - startAngle > 180 ? 1 : 0;
    return `M ${start.x} ${start.y} A ${r} ${r} 0 ${largeArc} 1 ${end.x} ${end.y}`;
  };

  const cx = 150, cy = 150, r = 100;

  // Get color for fragment
  const getFragmentColor = (frag: Fragment, index: number): string => {
    if (frag.type && PART_TYPES[frag.type]) {
      return PART_TYPES[frag.type].color;
    }
    // Default color cycle
    const colors = ['#ef4444', '#f97316', '#22c55e', '#3b82f6', '#8b5cf6', '#ec4899'];
    return colors[index % colors.length];
  };

  return (
    <div className="iso-construct-viz">
      {/* View mode toggle */}
      <div className="iso-viz-toggle">
        <button
          className={viewMode === 'circular' ? 'active' : ''}
          onClick={() => setViewMode('circular')}
        >
          <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
            <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm0 18c-4.42 0-8-3.58-8-8s3.58-8 8-8 8 3.58 8 8-3.58 8-8 8z"/>
          </svg>
          Circular
        </button>
        <button
          className={viewMode === 'linear' ? 'active' : ''}
          onClick={() => setViewMode('linear')}
        >
          <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
            <path d="M3 13h2v-2H3v2zm0 4h2v-2H3v2zm0-8h2V7H3v2zm4 4h14v-2H7v2zm0 4h14v-2H7v2zM7 7v2h14V7H7z"/>
          </svg>
          Linear
        </button>
      </div>

      {/* Circular View */}
      {viewMode === 'circular' && (
        <div className="iso-circular-view">
          <svg viewBox="0 0 300 300" className="iso-plasmid-svg">
            {/* Background circle */}
            <circle cx={cx} cy={cy} r={r} fill="none" stroke="#e5e7eb" strokeWidth="24" />

            {/* Fragment arcs */}
            {fragmentAngles.map((frag, i) => {
              const color = getFragmentColor(frag, i);
              const isHovered = hoveredFragment === i;
              return (
                <g key={i}>
                  <path
                    d={describeArc(cx, cy, r, frag.startAngle || 0, (frag.startAngle || 0) + (frag.sweepAngle || 0) - 2)}
                    fill="none"
                    stroke={color}
                    strokeWidth={isHovered ? 28 : 24}
                    strokeLinecap="round"
                    className="cursor-pointer"
                    style={{ transition: 'stroke-width 0.2s' }}
                    onMouseEnter={() => setHoveredFragment(i)}
                    onMouseLeave={() => setHoveredFragment(null)}
                  />
                </g>
              );
            })}

            {/* Junction markers */}
            {junctions?.map((junction, i) => {
              const frag = fragmentAngles[i];
              if (!frag) return null;
              const angle = ((frag.startAngle || 0) * Math.PI) / 180;
              const jx = cx + r * Math.cos(angle);
              const jy = cy + r * Math.sin(angle);
              const isHovered = hoveredJunction === i;
              const colors = QUALITY_COLORS[junction.overlap?.qualityTier || 'good'] || QUALITY_COLORS.good;

              return (
                <g
                  key={`j-${i}`}
                  onMouseEnter={() => setHoveredJunction(i)}
                  onMouseLeave={() => setHoveredJunction(null)}
                  className="cursor-pointer"
                >
                  <circle
                    cx={jx}
                    cy={jy}
                    r={isHovered ? 8 : 6}
                    fill={colors.badge}
                    stroke="#fff"
                    strokeWidth="2"
                    style={{ transition: 'r 0.2s' }}
                  />
                  {isHovered && (
                    <g>
                      <rect
                        x={jx - 35}
                        y={jy - 35}
                        width="70"
                        height="28"
                        rx="4"
                        fill="rgba(0,0,0,0.85)"
                      />
                      <text x={jx} y={jy - 24} textAnchor="middle" fill="#fff" fontSize="9" fontWeight="600">
                        Junction {i + 1}
                      </text>
                      <text x={jx} y={jy - 12} textAnchor="middle" fill="#a5f3fc" fontSize="8" fontFamily="Monaco, monospace">
                        {junction.overlap?.sequence?.slice(0, 12) || '...'}
                      </text>
                    </g>
                  )}
                </g>
              );
            })}

            {/* Part labels */}
            {fragmentAngles.map((frag, i) => {
              const labelR = r - 40;
              const angle = ((frag.midAngle || 0) * Math.PI) / 180;
              const lx = cx + labelR * Math.cos(angle);
              const ly = cy + labelR * Math.sin(angle);

              return (
                <text
                  key={`label-${i}`}
                  x={lx}
                  y={ly}
                  textAnchor="middle"
                  dominantBaseline="middle"
                  fontSize="10"
                  fontWeight="600"
                  fill={hoveredFragment === i ? getFragmentColor(frag, i) : '#374151'}
                >
                  {(frag.id || `F${i + 1}`).slice(0, 6)}
                </text>
              );
            })}

            {/* Center info */}
            <text x={cx} y={cy - 8} textAnchor="middle" fontSize="20" fontWeight="700" fill="#1f2937">
              {totalLength ? totalLength.toLocaleString() : '—'}
            </text>
            <text x={cx} y={cy + 10} textAnchor="middle" fontSize="10" fill="#6b7280">
              bp total
            </text>
          </svg>

          {/* Legend */}
          <div className="iso-plasmid-legend">
            {fragmentAngles?.map((frag, i) => (
              <div
                key={i}
                className={`iso-legend-item ${hoveredFragment === i ? 'active' : ''}`}
                onMouseEnter={() => setHoveredFragment(i)}
                onMouseLeave={() => setHoveredFragment(null)}
              >
                <span className="iso-legend-dot" style={{ backgroundColor: getFragmentColor(frag, i) }} />
                <span className="iso-legend-name">{frag.id || `Fragment ${i + 1}`}</span>
                <span className="iso-legend-bp">{frag.displayLength?.toLocaleString()} bp</span>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Linear View */}
      {viewMode === 'linear' && (
        <div className="iso-linear-view">
          <div className="iso-linear-track">
            {fragments?.map((frag, i) => {
              const color = getFragmentColor(frag, i);
              const junction = junctions?.[i];
              const junctionColors = QUALITY_COLORS[junction?.overlap?.qualityTier || 'good'] || QUALITY_COLORS.good;

              return (
                <React.Fragment key={i}>
                  {/* Junction node */}
                  <div
                    className={`iso-junction-node ${hoveredJunction === i ? 'highlighted' : ''}`}
                    onMouseEnter={() => setHoveredJunction(i)}
                    onMouseLeave={() => setHoveredJunction(null)}
                    style={{ borderColor: junctionColors.badge }}
                  >
                    <span className="iso-junction-label">J{i + 1}</span>
                    <span className="iso-junction-seq" style={{ color: junctionColors.badge }}>
                      {junction?.overlap?.sequence?.slice(0, 8) || '...'}
                    </span>
                  </div>

                  {/* Fragment block */}
                  <div
                    className={`iso-fragment-block ${hoveredFragment === i ? 'highlighted' : ''}`}
                    onMouseEnter={() => setHoveredFragment(i)}
                    onMouseLeave={() => setHoveredFragment(null)}
                    style={{ borderColor: color }}
                  >
                    <div className="iso-fragment-stripe" style={{ backgroundColor: color }} />
                    <div className="iso-fragment-info">
                      <span className="iso-fragment-name">{frag.id || `Fragment ${i + 1}`}</span>
                      <span className="iso-fragment-details">{(frag.originalLength || frag.length || 0).toLocaleString()} bp</span>
                    </div>
                  </div>
                </React.Fragment>
              );
            })}
            {/* Final junction */}
            {junctions && junctions.length > 0 && (
              <div
                className={`iso-junction-node ${hoveredJunction === junctions.length ? 'highlighted' : ''}`}
                style={{ borderColor: QUALITY_COLORS[junctions[0]?.overlap?.qualityTier || 'good']?.badge || '#3b82f6' }}
              >
                <span className="iso-junction-label">J1</span>
                <span className="iso-junction-seq">circular</span>
              </div>
            )}
          </div>
          <div className="iso-circular-badge">
            <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
              <path d="M12 4V1L8 5l4 4V6c3.31 0 6 2.69 6 6 0 1.01-.25 1.97-.7 2.8l1.46 1.46A7.93 7.93 0 0020 12c0-4.42-3.58-8-8-8zm0 14c-3.31 0-6-2.69-6-6 0-1.01.25-1.97.7-2.8L5.24 7.74A7.93 7.93 0 004 12c0 4.42 3.58 8 8 8v3l4-4-4-4v3z"/>
            </svg>
            Circular Assembly
          </div>
        </div>
      )}

      {/* Hovered fragment tooltip */}
      {hoveredFragment !== null && fragments?.[hoveredFragment] && (
        <div className="iso-viz-tooltip">
          <div className="iso-tooltip-header" style={{ borderColor: getFragmentColor(fragments[hoveredFragment], hoveredFragment) }}>
            {fragments[hoveredFragment].id || `Fragment ${hoveredFragment + 1}`}
          </div>
          <div className="iso-tooltip-content">
            <span>{(fragments[hoveredFragment].originalLength || fragments[hoveredFragment].length || 0).toLocaleString()} bp</span>
            {fragments[hoveredFragment].type && (
              <span className="iso-tooltip-type">{PART_TYPES[fragments[hoveredFragment].type!]?.name || 'Other'}</span>
            )}
          </div>
        </div>
      )}
    </div>
  );
}

// ============================================================================
// OVERLAP TM COMPARISON - Shows Tm differences between junctions
// ============================================================================
function JunctionCompatibilityMatrix({ junctions }: JunctionCompatibilityMatrixProps): React.ReactElement | null {
  if (!junctions || junctions.length < 2) return null;

  // Get junction Tms
  const junctionTms = useMemo<number[]>(() => {
    return junctions.map(j => j.overlap?.tm || 50);
  }, [junctions]);

  // Calculate Tm differences between junctions
  const matrix = useMemo(() => {
    return junctions.map((j1, i) =>
      junctions.map((j2, k) => {
        if (i === k) return { type: 'self', tm: junctionTms[i] };
        const tmDiff = Math.abs(junctionTms[i] - junctionTms[k]);
        return {
          type: tmDiff <= 2 ? 'caution' : tmDiff <= 5 ? 'ok' : 'good',
          tmDiff,
          tm1: junctionTms[i],
          tm2: junctionTms[k],
        };
      })
    );
  }, [junctions, junctionTms]);

  const getCellColor = (cell: { type: string; tm?: number; tmDiff?: number; tm1?: number; tm2?: number }): string => {
    if (cell.type === 'self') return '#f3f4f6';
    if (cell.type === 'caution') return '#fef3c7'; // Yellow - similar Tm, mispriming risk
    if (cell.type === 'ok') return '#dbeafe'; // Blue - moderate
    return '#dcfce7'; // Green - different Tm, good
  };

  return (
    <div className="iso-compatibility-matrix">
      <div className="iso-matrix-header">
        <h5>
          <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
            <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zM9 17H7v-7h2v7zm4 0h-2V7h2v10zm4 0h-2v-4h2v4z"/>
          </svg>
          Overlap Tm Comparison
        </h5>
        <span className="iso-matrix-hint">Shows ΔTm between junction overlaps</span>
      </div>

      {/* Junction Tm summary */}
      <div className="iso-junction-tm-summary">
        {junctions.map((j, i) => (
          <div key={i} className="iso-junction-tm-chip">
            <span className="iso-tm-label">J{i + 1}</span>
            <span className="iso-tm-value">{junctionTms[i]?.toFixed(1)}°C</span>
          </div>
        ))}
      </div>

      <div className="iso-matrix-grid" style={{ gridTemplateColumns: `auto repeat(${junctions.length}, 1fr)` }}>
        {/* Header row */}
        <div className="iso-matrix-cell header corner"></div>
        {junctions.map((_, i) => (
          <div key={`h-${i}`} className="iso-matrix-cell header">J{i + 1}</div>
        ))}

        {/* Data rows */}
        {matrix.map((row, i) => (
          <React.Fragment key={`row-${i}`}>
            <div className="iso-matrix-cell header">J{i + 1}</div>
            {row.map((cell, k) => (
              <div
                key={`cell-${i}-${k}`}
                className={`iso-matrix-cell ${cell.type}`}
                style={{ backgroundColor: getCellColor(cell) }}
                title={cell.type === 'self'
                  ? `J${i + 1} Tm: ${cell.tm?.toFixed(1)}°C`
                  : `J${i + 1} (${cell.tm1?.toFixed(1)}°C) vs J${k + 1} (${cell.tm2?.toFixed(1)}°C)`}
              >
                {cell.type === 'self' ? (
                  <span className="iso-matrix-tm">{cell.tm?.toFixed(0)}°</span>
                ) : (
                  <span className="iso-matrix-value">Δ{cell.tmDiff?.toFixed(0)}°</span>
                )}
              </div>
            ))}
          </React.Fragment>
        ))}
      </div>
      <div className="iso-matrix-legend">
        <span className="iso-legend-item">
          <span className="iso-legend-box bg-green-100" />
          <span>ΔTm &gt;5°C</span>
          <span className="iso-legend-status good">Good</span>
        </span>
        <span className="iso-legend-item">
          <span className="iso-legend-box bg-blue-100" />
          <span>ΔTm 2-5°C</span>
          <span className="iso-legend-status ok">OK</span>
        </span>
        <span className="iso-legend-item">
          <span className="iso-legend-box bg-yellow-100" />
          <span>ΔTm &lt;2°C</span>
          <span className="iso-legend-status caution">Caution</span>
        </span>
      </div>
      <p className="iso-matrix-explainer">
        Similar Tm values between overlaps may increase mispriming risk during assembly.
      </p>
    </div>
  );
}

// ============================================================================
// ENHANCED PRIMER CARD - Better visual structure display
// ============================================================================
function EnhancedPrimerCard({ primer, direction, fragmentId, onCopy }: EnhancedPrimerCardProps): React.ReactElement {
  const [copied, setCopied] = useState<boolean>(false);
  const [showBreakdown, setShowBreakdown] = useState<boolean>(false);
  const safePrimer = primer || {};
  const colors = QUALITY_COLORS[safePrimer.qualityTier || 'marginal'] || QUALITY_COLORS.marginal;

  const handleCopy = (): void => {
    if (safePrimer.sequence) {
      navigator.clipboard.writeText(safePrimer.sequence);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
      onCopy?.(safePrimer.sequence);
    }
  };

  const formatTm = (tm?: number): string => (typeof tm === 'number' && !isNaN(tm)) ? tm.toFixed(1) : '-';
  const formatDG = (dg?: number): string => (typeof dg === 'number' && !isNaN(dg)) ? dg.toFixed(1) : '-';

  // Calculate tail and annealing lengths
  const tailLength = safePrimer.homologyTail?.length || 0;
  const annealingLength = safePrimer.annealingRegion?.length || (safePrimer.sequence?.length || 0) - tailLength;

  return (
    <div className="iso-primer-card-enhanced" style={{ borderColor: colors.border }}>
      {/* Header */}
      <div className="iso-primer-header">
        <div className="iso-primer-direction" style={{ backgroundColor: direction === 'Forward' ? '#3b82f6' : '#8b5cf6' }}>
          {direction === 'Forward' ? '5′→3′' : '3′←5′'}
        </div>
        <span className="iso-primer-name">{fragmentId}_{direction[0]}</span>
        <QualityBadge tier={safePrimer.qualityTier} score={safePrimer.qualityScore} />
        <button className="iso-copy-btn" onClick={handleCopy} title="Copy sequence">
          {copied ? (
            <svg viewBox="0 0 24 24" width="16" height="16" fill="#22c55e">
              <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
            </svg>
          ) : (
            <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
              <path d="M16 1H4c-1.1 0-2 .9-2 2v14h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z" />
            </svg>
          )}
        </button>
      </div>

      {/* Enhanced Sequence Visualization */}
      <div className="iso-primer-structure">
        <div className="iso-structure-bar">
          {safePrimer.homologyTail && (
            <div
              className="iso-structure-segment tail"
              style={{ flex: tailLength }}
              title={`Homology tail: ${tailLength}bp`}
            >
              <span className="iso-segment-label">Overlap</span>
            </div>
          )}
          <div
            className="iso-structure-segment annealing"
            style={{ flex: annealingLength }}
            title={`Annealing region: ${annealingLength}bp`}
          >
            <span className="iso-segment-label">Annealing</span>
          </div>
        </div>
        <div className="iso-structure-lengths">
          {tailLength > 0 && <span className="iso-length-label tail">{tailLength}bp</span>}
          <span className="iso-length-label annealing">{annealingLength}bp</span>
        </div>
      </div>

      {/* Sequence Display */}
      <div className="iso-primer-sequence-display">
        {safePrimer.homologyTail && (
          <span className="iso-seq-tail">{safePrimer.homologyTail}</span>
        )}
        <span className="iso-seq-annealing">{safePrimer.annealingRegion || safePrimer.sequence?.slice(tailLength) || ''}</span>
      </div>

      {/* Metrics Grid */}
      <div className="iso-primer-metrics-grid">
        <div className="iso-metric-item">
          <span className="iso-metric-label">Total</span>
          <span className="iso-metric-value">{safePrimer.length || '-'} bp</span>
        </div>
        <div className="iso-metric-item">
          <span className="iso-metric-label">Anneal Tm</span>
          <span className="iso-metric-value">{formatTm(safePrimer.tm)}°C</span>
        </div>
        <div className="iso-metric-item">
          <span className="iso-metric-label">Full Tm</span>
          <span className="iso-metric-value">{formatTm(safePrimer.tmFull)}°C</span>
        </div>
        <div className="iso-metric-item">
          <span className="iso-metric-label">GC</span>
          <span className="iso-metric-value">{safePrimer.gcPercent || '-'}</span>
        </div>
        <div className="iso-metric-item">
          <span className="iso-metric-label">Hairpin</span>
          <span className="iso-metric-value" style={{ color: (safePrimer.hairpinDG || 0) >= -2 ? '#22c55e' : (safePrimer.hairpinDG || 0) >= -5 ? '#f59e0b' : '#ef4444' }}>
            {formatDG(safePrimer.hairpinDG)}
          </span>
        </div>
        <div className="iso-metric-item">
          <span className="iso-metric-label">3′ Clamp</span>
          <span className="iso-metric-value">
            {safePrimer.sequence ? (() => {
              const last2 = safePrimer.sequence.slice(-2).toUpperCase();
              const gcCount = (last2.match(/[GC]/g) || []).length;
              // 1 G/C in last 2 = ideal, 2 = good, 0 = weak
              if (gcCount >= 1) {
                return <span className="text-green-500">{last2}</span>;
              }
              return <span className="text-amber-500">{last2}</span>;
            })() : '—'}
          </span>
        </div>
      </div>

      {/* Score Breakdown Toggle */}
      {safePrimer.scoreBreakdown && (
        <div className="iso-primer-breakdown">
          <button
            className={`iso-breakdown-toggle ${showBreakdown ? 'expanded' : ''}`}
            onClick={() => setShowBreakdown(!showBreakdown)}
          >
            <span>Score Breakdown</span>
            <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
              <path d={showBreakdown ? "M7.41 15.41L12 10.83l4.59 4.58L18 14l-6-6-6 6z" : "M7.41 8.59L12 13.17l4.59-4.58L18 10l-6 6-6-6 1.41-1.41z"}/>
            </svg>
          </button>
          {showBreakdown && (
            <ScoreBreakdownPanel
              scores={safePrimer.scoreBreakdown}
              weights={{
                tmInRange: 0.15, gcInRange: 0.08, lengthInRange: 0.07,
                gcClamp: 0.16, terminal3DG: 0.18, hairpin: 0.12,
                homopolymer: 0.06, threePrimeComp: 0.06, gQuadruplex: 0.08,
                selfComplementarity: 0.04,
              }}
              title=""
            />
          )}
        </div>
      )}
    </div>
  );
}

// Score bar component
function ScoreBar({ score, label, max = 100, showValue = true }: ScoreBarProps): React.ReactElement {
  const percentage = Math.min(100, (score / max) * 100);
  const color = getScoreColor(score);

  return (
    <div className="iso-score-bar">
      <div className="iso-score-bar-header">
        <span className="iso-score-bar-label">{label}</span>
        {showValue && <span className="iso-score-bar-value" style={{ color }}>{score.toFixed(1)}</span>}
      </div>
      <div className="iso-score-bar-track">
        <div
          className="iso-score-bar-fill"
          style={{ width: `${percentage}%`, backgroundColor: color }}
        />
      </div>
    </div>
  );
}

// Quality badge component
function QualityBadge({ tier = 'marginal', score }: QualityBadgeProps): React.ReactElement {
  const safeTier = tier || 'marginal';
  const colors = QUALITY_COLORS[safeTier] || QUALITY_COLORS.marginal;
  return (
    <span
      className="iso-quality-badge"
      style={{
        backgroundColor: colors.bg,
        color: colors.text,
        borderColor: colors.border,
      }}
    >
      <span className="iso-quality-dot" style={{ backgroundColor: colors.badge }} />
      {safeTier.charAt(0).toUpperCase() + safeTier.slice(1)}
      {score !== undefined && <span className="iso-quality-score">({typeof score === 'number' ? score.toFixed(0) : score})</span>}
    </span>
  );
}

// Detailed score breakdown panel
function ScoreBreakdownPanel({ scores, weights, title }: ScoreBreakdownPanelProps): React.ReactElement {
  const weightedScores = useMemo(() => {
    const entries: Array<{
      key: string;
      label: string;
      score: number;
      weight: number;
      contribution: number;
    }> = [];
    for (const [key, score] of Object.entries(scores)) {
      const weight = weights[key as keyof ScoreWeights] || 0;
      if (weight > 0) {
        entries.push({
          key,
          label: formatScoreLabel(key),
          score: score * 100,
          weight,
          contribution: score * weight * 100,
        });
      }
    }
    return entries.sort((a, b) => b.contribution - a.contribution);
  }, [scores, weights]);

  return (
    <div className="iso-score-breakdown">
      {title && <h5>{title}</h5>}
      <div className="iso-breakdown-list">
        {weightedScores.map(({ key, label, score, weight, contribution }) => (
          <div key={key} className="iso-breakdown-item">
            <div className="iso-breakdown-label">
              <span className="iso-breakdown-name">{label}</span>
              <span className="iso-breakdown-weight">({(weight * 100).toFixed(0)}%)</span>
            </div>
            <div className="iso-breakdown-bar-container">
              <div
                className="iso-breakdown-bar"
                style={{
                  width: `${score}%`,
                  backgroundColor: getScoreColor(score),
                }}
              />
            </div>
            <span className="iso-breakdown-score">{score.toFixed(0)}</span>
          </div>
        ))}
      </div>
    </div>
  );
}

// Format score key to readable label
function formatScoreLabel(key: string): string {
  const labels: Record<string, string> = {
    tmInRange: 'Tm in Range',
    gcInRange: 'GC in Range',
    lengthInRange: 'Length in Range',
    gcClamp: 'GC Clamp',
    terminal3DG: "3' Terminal ΔG",
    hairpin: 'Hairpin',
    patternAvoidance: 'Pattern Avoidance',
    selfComplementarity: 'Self-Complementarity',
    crossJunctionHairpin: 'Cross-Junction',
    balancedGC: 'Balanced GC',
    homopolymer: 'Homopolymer',
    threePrimeComp: "3' Composition",
    gQuadruplex: 'G-Quadruplex',
  };
  return labels[key] || key.replace(/([A-Z])/g, ' $1').trim();
}

// Helper to classify warnings - shift messages are informational, not real warnings
function classifyWarning(warning: string): { isOptimizationNote: boolean; severity: WarningSeverity } {
  const shiftPatterns = [
    /shifted/i,
    /position shifted/i,
    /junction shifted/i,
    /start position/i,
  ];
  const isShiftMessage = shiftPatterns.some(p => p.test(warning));
  return {
    isOptimizationNote: isShiftMessage,
    severity: isShiftMessage ? 'info' :
      (warning.includes('G-quadruplex') || warning.includes('TTTT') || warning.includes('hairpin')) ? 'high' : 'medium'
  };
}

// Filter real warnings from optimization notes
function filterWarnings(warnings?: string[]): { realWarnings: string[]; optimizationNotes: string[] } {
  if (!warnings) return { realWarnings: [], optimizationNotes: [] };
  const realWarnings: string[] = [];
  const optimizationNotes: string[] = [];
  warnings.forEach(w => {
    const { isOptimizationNote } = classifyWarning(w);
    if (isOptimizationNote) {
      optimizationNotes.push(w);
    } else {
      realWarnings.push(w);
    }
  });
  return { realWarnings, optimizationNotes };
}

// Warning item component
function WarningItem({ warning, severity = 'medium' }: WarningItemProps): React.ReactElement {
  const severityColors: Record<WarningSeverity, { bg: string; border: string; icon: string }> = {
    high: { bg: '#fef2f2', border: '#ef4444', icon: '#dc2626' },
    medium: { bg: '#fffbeb', border: '#f59e0b', icon: '#d97706' },
    low: { bg: '#eff6ff', border: '#3b82f6', icon: '#2563eb' },
    info: { bg: '#f0fdf4', border: '#22c55e', icon: '#16a34a' },
  };
  const colors = severityColors[severity] || severityColors.medium;

  const icons: Record<WarningSeverity, React.ReactElement> = {
    high: <path d="M12 5.99L19.53 19H4.47L12 5.99M12 2L1 21h22L12 2zm1 14h-2v2h2v-2zm0-6h-2v4h2v-4z" />,
    medium: <path d="M12 5.99L19.53 19H4.47L12 5.99M12 2L1 21h22L12 2zm1 14h-2v2h2v-2zm0-6h-2v4h2v-4z" />,
    low: <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z" />,
    info: <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />,
  };

  return (
    <div
      className={`iso-warning-item ${severity}`}
      style={{
        backgroundColor: colors.bg,
        borderLeftColor: colors.border,
      }}
    >
      <svg viewBox="0 0 24 24" width="14" height="14" fill={colors.icon}>
        {icons[severity] || icons.medium}
      </svg>
      <span>{warning}</span>
    </div>
  );
}

// Junction info tooltip component
function JunctionInfoTooltip(): React.ReactElement {
  const [show, setShow] = useState<boolean>(false);

  return (
    <div className="iso-info-tooltip-wrapper">
      <button
        className="iso-info-btn"
        onMouseEnter={() => setShow(true)}
        onMouseLeave={() => setShow(false)}
        onClick={() => setShow(!show)}
      >
        <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
          <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 17h-2v-2h2v2zm2.07-7.75l-.9.92C13.45 12.9 13 13.5 13 15h-2v-.5c0-1.1.45-2.1 1.17-2.83l1.24-1.26c.37-.36.59-.86.59-1.41 0-1.1-.9-2-2-2s-2 .9-2 2H8c0-2.21 1.79-4 4-4s4 1.79 4 4c0 .88-.36 1.68-.93 2.25z"/>
        </svg>
      </button>
      {show && (
        <div className="iso-info-tooltip">
          <h6>What is a Junction?</h6>
          <p>
            A <strong>junction</strong> is the overlap region where two DNA fragments
            join during assembly. The enzyme chews back fragment ends, exposing
            complementary sequences that anneal together.
          </p>
          <div className="iso-info-metrics">
            <span><strong>Tm:</strong> Melting temperature - higher = stronger annealing</span>
            <span><strong>Length:</strong> 15-40bp typical for Gibson/NEBuilder</span>
            <span><strong>GC%:</strong> 40-60% is ideal for stability</span>
          </div>
        </div>
      )}
    </div>
  );
}

// Overlap card component
function OverlapCard({ overlap, index, isSelected, onSelect, showDetails = false, isClickable = false }: OverlapCardProps): React.ReactElement {
  const [expanded, setExpanded] = useState<boolean>(false);
  const colors = QUALITY_COLORS[overlap.qualityTier || 'marginal'] || QUALITY_COLORS.marginal;

  return (
    <div
      className={`iso-overlap-card ${isSelected ? 'selected' : ''} ${isClickable ? 'clickable' : ''}`}
      style={{ borderColor: isSelected ? colors.badge : colors.border }}
      onClick={() => onSelect?.(overlap)}
    >
      <div className="iso-overlap-header">
        <div className="iso-overlap-title">
          <span className="iso-overlap-index">Junction {index + 1}</span>
          <QualityBadge tier={overlap.qualityTier} score={overlap.compositeScore} />
        </div>
        {overlap.offset && overlap.offset > 0 && (
          <span className="iso-offset-badge">+{overlap.offset}bp shifted</span>
        )}
      </div>

      <div className="iso-overlap-metrics">
        <div className="iso-metric">
          <span className="iso-metric-label">Length</span>
          <span className="iso-metric-value">{overlap.length} bp</span>
        </div>
        <div className="iso-metric">
          <span className="iso-metric-label">Tm</span>
          <span className="iso-metric-value">{overlap.tm}°C</span>
        </div>
        <div className="iso-metric">
          <span className="iso-metric-label">GC</span>
          <span className="iso-metric-value">{overlap.gcPercent}</span>
        </div>
        <div className="iso-metric">
          <span className="iso-metric-label">3' ΔG</span>
          <span className="iso-metric-value">{overlap.terminal3DG} kcal/mol</span>
        </div>
      </div>

      <div className="iso-overlap-sequence">
        <code>{overlap.sequence}</code>
      </div>

      {overlap.warnings && overlap.warnings.length > 0 && (
        <div className="iso-overlap-warnings">
          {overlap.warnings.slice(0, 2).map((w, i) => (
            <WarningItem
              key={i}
              warning={w}
              severity={w.includes('G-quadruplex') || w.includes('TTTT') ? 'high' : 'medium'}
            />
          ))}
          {overlap.warnings.length > 2 && (
            <button className="iso-more-warnings" onClick={(e) => { e.stopPropagation(); setExpanded(!expanded); }}>
              {expanded ? 'Show less' : `+${overlap.warnings.length - 2} more`}
            </button>
          )}
        </div>
      )}

      {(showDetails || expanded) && overlap.scores && (
        <div className="iso-overlap-details">
          <ScoreBreakdownPanel
            scores={overlap.scores}
            weights={{
              tmInRange: 0.16, gcInRange: 0.10, lengthInRange: 0.07,
              gcClamp: 0.14, terminal3DG: 0.14, hairpin: 0.12,
              patternAvoidance: 0.12, selfComplementarity: 0.08,
              balancedGC: 0.05, crossJunctionHairpin: 0.02,
            }}
            title="Quality Factors"
          />
        </div>
      )}
    </div>
  );
}

// Pair quality panel
function PairQualityPanel({ pair }: PairQualityPanelProps): React.ReactElement | null {
  if (!pair) return null;

  return (
    <div className="iso-pair-panel">
      <div className="iso-pair-header">
        <h5>Primer Pair Quality</h5>
        <QualityBadge tier={pair.qualityTier} score={pair.pairScore} />
      </div>
      <div className="iso-pair-metrics">
        <div className="iso-pair-metric">
          <span className="iso-pair-label">Annealing Temp</span>
          <span className="iso-pair-value">{pair.annealingTemp}°C</span>
        </div>
        <div className="iso-pair-metric">
          <span className="iso-pair-label">Tm Difference</span>
          <span className="iso-pair-value" style={{ color: (pair.tmDifference || 0) <= 3 ? '#22c55e' : (pair.tmDifference || 0) <= 5 ? '#f59e0b' : '#ef4444' }}>
            {pair.tmDifference}°C
          </span>
        </div>
        <div className="iso-pair-metric">
          <span className="iso-pair-label">Heterodimer ΔG</span>
          <span className="iso-pair-value" style={{ color: (pair.heterodimerDG || 0) >= -6 ? '#22c55e' : (pair.heterodimerDG || 0) >= -9 ? '#f59e0b' : '#ef4444' }}>
            {pair.heterodimerDG} kcal/mol
          </span>
        </div>
      </div>
    </div>
  );
}

// Alternatives panel
function AlternativesPanel({ alternatives, onSelect }: AlternativesPanelProps): React.ReactElement | null {
  if (!alternatives || alternatives.length === 0) return null;

  return (
    <div className="iso-alternatives-panel">
      <h5>Alternative Designs</h5>
      <div className="iso-alternatives-list">
        {alternatives.map((alt, i) => (
          <button
            key={i}
            className="iso-alternative-item"
            onClick={() => onSelect?.(alt)}
          >
            <span className="iso-alt-reason">{alt.selectionReason}</span>
            <span className="iso-alt-score">Score: {alt.pairScore || alt.compositeScore?.toFixed(0)}</span>
          </button>
        ))}
      </div>
    </div>
  );
}

// Optimization summary
function OptimizationSummary({ summary }: OptimizationSummaryProps): React.ReactElement | null {
  if (!summary) return null;

  return (
    <div className="iso-optimization-summary">
      <div className="iso-summary-header">
        <svg viewBox="0 0 24 24" width="16" height="16" fill="#7c3aed">
          <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm-5 14H7v-2h7v2zm3-4H7v-2h10v2zm0-4H7V7h10v2z" />
        </svg>
        <span>Optimization Summary</span>
      </div>
      <div className="iso-summary-stats">
        <div className="iso-stat">
          <span className="iso-stat-value">{summary.totalEvaluated}</span>
          <span className="iso-stat-label">Candidates Evaluated</span>
        </div>
        <div className="iso-stat">
          <span className="iso-stat-value">{summary.validCount}</span>
          <span className="iso-stat-label">Valid Candidates</span>
        </div>
        <div className="iso-stat">
          <span className="iso-stat-value" style={{ color: getScoreColor(summary.bestScore || 0) }}>
            {summary.bestScore?.toFixed(1)}
          </span>
          <span className="iso-stat-label">Best Score</span>
        </div>
      </div>
      {summary.scoreDistribution && (
        <div className="iso-distribution">
          <div className="iso-dist-bar">
            <div
              className="iso-dist-segment excellent"
              style={{ flex: summary.scoreDistribution.excellent }}
              title={`${summary.scoreDistribution.excellent} excellent`}
            />
            <div
              className="iso-dist-segment good"
              style={{ flex: summary.scoreDistribution.good }}
              title={`${summary.scoreDistribution.good} good`}
            />
            <div
              className="iso-dist-segment acceptable"
              style={{ flex: summary.scoreDistribution.acceptable }}
              title={`${summary.scoreDistribution.acceptable} acceptable`}
            />
            <div
              className="iso-dist-segment marginal"
              style={{ flex: summary.scoreDistribution.marginal }}
              title={`${summary.scoreDistribution.marginal} marginal`}
            />
          </div>
          <div className="iso-dist-legend">
            <span className="iso-legend-item"><span className="iso-legend-dot excellent" />Excellent</span>
            <span className="iso-legend-item"><span className="iso-legend-dot good" />Good</span>
            <span className="iso-legend-item"><span className="iso-legend-dot acceptable" />Acceptable</span>
            <span className="iso-legend-item"><span className="iso-legend-dot marginal" />Marginal</span>
          </div>
        </div>
      )}
    </div>
  );
}

// ============================================================================
// MAIN COMPONENT
// ============================================================================
export default function IsothermalAssemblyPanel({
  fragments,
  variant = 'nebuilder_hifi',
  config = {},
  onPrimersCopied,
}: IsothermalAssemblyPanelProps): React.ReactElement {
  const [activeTab, setActiveTab] = useState<'overview' | 'construct' | 'primers' | 'junctions' | 'protocol'>('overview');
  const [selectedJunction, setSelectedJunction] = useState<number>(0);
  const [showAdvanced, setShowAdvanced] = useState<boolean>(false);

  // Merge config with defaults and enable autoOptimize
  const assemblyConfig = useMemo<AssemblyConfig>(() => ({
    ...DEFAULT_ASSEMBLY_CONFIG,
    ...config,
    method: variant === 'nebuilder_hifi' ? 'NEBUILDER_HIFI' : 'GIBSON',
    autoOptimize: true,
  }), [variant, config]);

  // Run assembly design
  const assemblyResult = useMemo<AssemblyResult | null>(() => {
    if (!fragments || fragments.length < 2) return null;

    try {
      return designAssembly(fragments, assemblyConfig) as AssemblyResult;
    } catch (error) {
      console.error('Assembly design error:', error);
      return { error: (error as Error).message };
    }
  }, [fragments, assemblyConfig]);

  // Get total length from assembly result (already calculated by designAssembly)
  const totalLength = useMemo<number>(() => {
    // Use pre-calculated totalLength from assembly result, or calculate from fragments
    if (assemblyResult?.assembly?.totalLength) {
      return assemblyResult.assembly.totalLength;
    }
    if (!assemblyResult?.fragments) return 0;
    return assemblyResult.fragments.reduce((sum, f) => sum + (f.originalLength || f.length || 0), 0);
  }, [assemblyResult]);

  // Handle copy
  const handleCopy = useCallback((text: string) => {
    onPrimersCopied?.(text);
  }, [onPrimersCopied]);

  if (!fragments || fragments.length < 2) {
    return (
      <div className="iso-empty-state">
        <svg viewBox="0 0 24 24" width="48" height="48" fill="#9ca3af">
          <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm-5 14H7v-2h7v2zm3-4H7v-2h10v2zm0-4H7V7h10v2z" />
        </svg>
        <h4>Add at least 2 fragments to design assembly primers</h4>
        <p>Enter your DNA sequences to get optimized primers with homology regions.</p>
      </div>
    );
  }

  if (assemblyResult?.error) {
    return (
      <div className="iso-error-state">
        <svg viewBox="0 0 24 24" width="48" height="48" fill="#ef4444">
          <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z" />
        </svg>
        <h4>Design Error</h4>
        <p>{assemblyResult.error}</p>
      </div>
    );
  }

  if (!assemblyResult) {
    return (
      <div className="iso-loading-state">
        <div className="iso-spinner" />
        <p>Optimizing assembly design...</p>
      </div>
    );
  }

  const { fragments: fragmentDesigns, junctions, quality = {}, warnings, optimization } = assemblyResult;
  const safeQuality: AssemblyQuality = {
    tier: quality?.tier || 'good',
    averageScore: quality?.averageScore ?? 0,
    minimumScore: quality?.minimumScore ?? 0,
    weakestJunction: quality?.weakestJunction,
  };

  return (
    <div className="iso-assembly-panel">
      {/* Header */}
      <div className="iso-panel-header">
        <div className="iso-header-title">
          <svg viewBox="0 0 24 24" width="24" height="24" fill="#7c3aed">
            <path d="M12 3C8 3 4 6 4 12s4 9 8 9" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" fill="none" />
            <path d="M12 3c4 0 8 3 8 9s-4 9-8 9" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeDasharray="3 3" fill="none" />
            <circle cx="12" cy="12" r="3" fill="currentColor" />
          </svg>
          <div>
            <h3>{variant === 'nebuilder_hifi' ? 'NEBuilder HiFi' : 'Gibson'} Assembly</h3>
            <span className="iso-subtitle">Optimized Primer Design</span>
          </div>
        </div>
        <div className="iso-header-quality">
          <QualityBadge tier={safeQuality.tier} score={safeQuality.averageScore} />
        </div>
      </div>

      {/* Optimization badge */}
      {optimization && (
        <div className="iso-optimization-badge">
          <svg viewBox="0 0 24 24" width="14" height="14" fill="#7c3aed">
            <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
          </svg>
          <span>{optimization.description}</span>
        </div>
      )}

      {/* Tab navigation */}
      <div className="iso-tabs">
        <button
          className={`iso-tab ${activeTab === 'overview' ? 'active' : ''}`}
          onClick={() => setActiveTab('overview')}
        >
          Overview
        </button>
        <button
          className={`iso-tab ${activeTab === 'construct' ? 'active' : ''}`}
          onClick={() => setActiveTab('construct')}
        >
          Construct
        </button>
        <button
          className={`iso-tab ${activeTab === 'primers' ? 'active' : ''}`}
          onClick={() => setActiveTab('primers')}
        >
          Primers ({(fragmentDesigns?.length || 0) * 2})
        </button>
        <button
          className={`iso-tab ${activeTab === 'junctions' ? 'active' : ''}`}
          onClick={() => setActiveTab('junctions')}
        >
          Junctions ({junctions?.length || 0})
        </button>
        <button
          className={`iso-tab ${activeTab === 'protocol' ? 'active' : ''}`}
          onClick={() => setActiveTab('protocol')}
        >
          Protocol
        </button>
      </div>

      {/* Tab content */}
      <div className="iso-tab-content">
        {/* Overview Tab - Enhanced with Gauge and Export */}
        {activeTab === 'overview' && (
          <div className="iso-overview">
            {/* Quality Gauge + Export Panel */}
            <div className="iso-overview-header">
              <AssemblyQualityGauge
                quality={safeQuality}
                junctionCount={junctions?.length || 0}
                fragmentCount={fragmentDesigns?.length || 0}
              />
              <ExportPanel
                fragments={fragmentDesigns}
                junctions={junctions}
                variant={variant}
                onExport={handleCopy}
              />
            </div>

            {/* Warnings - filtered to exclude optimization notes */}
            {(() => {
              const { realWarnings, optimizationNotes } = filterWarnings(warnings);
              return (
                <>
                  {realWarnings.length > 0 && (
                    <div className="iso-warnings-section">
                      <h5>
                        <svg viewBox="0 0 24 24" width="16" height="16" fill="#f59e0b">
                          <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z" />
                        </svg>
                        Warnings ({realWarnings.length})
                      </h5>
                      <div className="iso-warnings-list">
                        {realWarnings.slice(0, showAdvanced ? realWarnings.length : 3).map((w, i) => {
                          const { severity } = classifyWarning(w);
                          return <WarningItem key={i} warning={w} severity={severity} />;
                        })}
                        {realWarnings.length > 3 && !showAdvanced && (
                          <button className="iso-show-more" onClick={() => setShowAdvanced(true)}>
                            Show {realWarnings.length - 3} more warnings
                          </button>
                        )}
                      </div>
                    </div>
                  )}

                  {/* Optimization notes - collapsible, less prominent */}
                  {optimizationNotes.length > 0 && (
                    <details className="iso-optimization-notes">
                      <summary>
                        <svg viewBox="0 0 24 24" width="14" height="14" fill="#16a34a">
                          <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
                        </svg>
                        <span>Optimization applied ({optimizationNotes.length} adjustments)</span>
                      </summary>
                      <div className="iso-opt-notes-list">
                        {optimizationNotes.map((w, i) => (
                          <WarningItem key={i} warning={w} severity="info" />
                        ))}
                      </div>
                    </details>
                  )}
                </>
              );
            })()}

            {/* Quick primer list */}
            <div className="iso-primer-list-quick">
              <h5>All Primers</h5>
              <div className="iso-primer-table">
                <div className="iso-primer-table-header">
                  <span>Fragment</span>
                  <span>Forward</span>
                  <span>Reverse</span>
                  <span>Quality</span>
                </div>
                {fragmentDesigns?.map((frag, i) => (
                  <div key={i} className="iso-primer-table-row">
                    <span className="iso-frag-name">{frag.id || `Fragment ${i + 1}`}</span>
                    <span className="iso-primer-seq" title={frag.primers?.forward?.sequence || ''}>
                      {(frag.primers?.forward?.sequence || '').slice(0, 20)}...
                    </span>
                    <span className="iso-primer-seq" title={frag.primers?.reverse?.sequence || ''}>
                      {(frag.primers?.reverse?.sequence || '').slice(0, 20)}...
                    </span>
                    <QualityBadge tier={frag.primers?.pair?.qualityTier || 'good'} />
                  </div>
                ))}
              </div>
            </div>
          </div>
        )}

        {/* Construct Tab - Unified Visualization */}
        {activeTab === 'construct' && (
          <div className="iso-construct-view">
            <ConstructVisualization
              fragments={fragmentDesigns}
              junctions={junctions}
              variant={variant}
              totalLength={totalLength}
            />

            {/* Assembly Summary */}
            <div className="iso-assembly-summary">
              <div className="iso-summary-card">
                <span className="iso-summary-label">Total Length</span>
                <span className="iso-summary-value">{totalLength.toLocaleString()} bp</span>
              </div>
              <div className="iso-summary-card">
                <span className="iso-summary-label">Fragments</span>
                <span className="iso-summary-value">{fragmentDesigns?.length || 0}</span>
              </div>
              <div className="iso-summary-card">
                <span className="iso-summary-label">Junctions</span>
                <span className="iso-summary-value">{junctions?.length || 0}</span>
              </div>
              <div className="iso-summary-card">
                <span className="iso-summary-label">Topology</span>
                <span className="iso-summary-value">Circular</span>
              </div>
            </div>

            {/* Junction Compatibility Matrix */}
            {junctions && junctions.length >= 2 && (
              <JunctionCompatibilityMatrix junctions={junctions} />
            )}
          </div>
        )}

        {/* Primers Tab - Enhanced */}
        {activeTab === 'primers' && (
          <div className="iso-primers-view">
            {fragmentDesigns?.map((frag, i) => (
              <div key={i} className="iso-fragment-primers">
                <h5 className="iso-fragment-title">
                  <span className="iso-fragment-badge">{i + 1}</span>
                  {frag.id || `Fragment ${i + 1}`}
                  <span className="iso-fragment-length">{(frag.originalLength || frag.length || 0).toLocaleString()} bp</span>
                </h5>

                <div className="iso-primers-grid">
                  <EnhancedPrimerCard
                    primer={frag.primers?.forward}
                    direction="Forward"
                    fragmentId={frag.id || `Fragment_${i + 1}`}
                    onCopy={handleCopy}
                  />
                  <EnhancedPrimerCard
                    primer={frag.primers?.reverse}
                    direction="Reverse"
                    fragmentId={frag.id || `Fragment_${i + 1}`}
                    onCopy={handleCopy}
                  />
                </div>

                <PairQualityPanel pair={frag.primers?.pair} />

                {frag.primers?.alternatives && frag.primers.alternatives.length > 0 && (
                  <AlternativesPanel alternatives={frag.primers.alternatives} />
                )}

                {frag.primers?.warnings && frag.primers.warnings.length > 0 && (
                  <div className="iso-primer-warnings">
                    {frag.primers.warnings.map((w, j) => (
                      <WarningItem key={j} warning={w} severity="low" />
                    ))}
                  </div>
                )}
              </div>
            ))}
          </div>
        )}

        {/* Junctions Tab */}
        {activeTab === 'junctions' && (
          <div className="iso-junctions-view">
            {/* Header with info tooltip */}
            <div className="iso-junctions-header">
              <h4>Junction Overlaps</h4>
              <JunctionInfoTooltip />
            </div>

            <div className="iso-junction-selector">
              {junctions?.map((j, i) => (
                <button
                  key={i}
                  className={`iso-junction-btn ${selectedJunction === i ? 'active' : ''}`}
                  onClick={() => setSelectedJunction(i)}
                >
                  <span>J{i + 1}</span>
                  <QualityBadge tier={j.overlap?.qualityTier} />
                </button>
              ))}
            </div>

            {junctions?.[selectedJunction] && (
              <div className="iso-junction-detail">
                {/* Current selection label */}
                <div className="iso-current-selection-label">
                  <span>Current Selection</span>
                </div>
                <OverlapCard
                  overlap={junctions[selectedJunction].overlap!}
                  index={selectedJunction}
                  isSelected={true}
                  showDetails={true}
                />

                {junctions[selectedJunction].alternatives && junctions[selectedJunction].alternatives!.length > 0 && (
                  <div className="iso-junction-alternatives">
                    <h5>
                      Alternative Overlaps
                      <span className="iso-alt-hint">Click to use instead</span>
                    </h5>
                    <div className="iso-alternatives-grid">
                      {junctions[selectedJunction].alternatives!.map((alt, i) => (
                        <OverlapCard
                          key={i}
                          overlap={alt}
                          index={i}
                          isSelected={false}
                          onSelect={() => {
                            // Swap the current overlap with the alternative
                            // Note: This would need backend support to actually regenerate primers
                            // For now, show a message that this is the selected alternative
                            alert(`Alternative overlap selected:\n\nSequence: ${alt.sequence}\nTm: ${alt.tm}°C\n\nNote: Swapping overlaps would require regenerating primers. This feature is coming soon.`);
                          }}
                          isClickable={true}
                        />
                      ))}
                    </div>
                  </div>
                )}

                {junctions[selectedJunction].optimizationSummary && (
                  <OptimizationSummary summary={junctions[selectedJunction].optimizationSummary} />
                )}
              </div>
            )}
          </div>
        )}

        {/* Protocol Tab */}
        {activeTab === 'protocol' && (
          <div className="iso-protocol-view">
            <div className="iso-protocol-section">
              <h5>
                <svg viewBox="0 0 24 24" width="18" height="18" fill="#7c3aed">
                  <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm-5 14H7v-2h7v2zm3-4H7v-2h10v2zm0-4H7V7h10v2z"/>
                </svg>
                Assembly Method
              </h5>
              <div className="iso-method-card">
                <div className="iso-method-header">
                  <span className="iso-method-name">{variant === 'nebuilder_hifi' ? 'NEBuilder HiFi DNA Assembly' : 'Gibson Assembly'}</span>
                  <span className="iso-method-badge">{variant === 'nebuilder_hifi' ? 'NEB E5520' : 'NEB E5510'}</span>
                </div>
                <p className="iso-method-desc">
                  {variant === 'nebuilder_hifi'
                    ? 'High-fidelity isothermal assembly with optimized overlap regions. Recommended for complex assemblies and high accuracy requirements.'
                    : 'Standard isothermal assembly using T5 exonuclease, Phusion polymerase, and Taq ligase. Reliable for routine assemblies.'}
                </p>
              </div>
            </div>

            <div className="iso-protocol-section">
              <h5>
                <svg viewBox="0 0 24 24" width="18" height="18" fill="#7c3aed">
                  <path d="M11.99 2C6.47 2 2 6.48 2 12s4.47 10 9.99 10C17.52 22 22 17.52 22 12S17.52 2 11.99 2zM12 20c-4.42 0-8-3.58-8-8s3.58-8 8-8 8 3.58 8 8-3.58 8-8 8zm.5-13H11v6l5.25 3.15.75-1.23-4.5-2.67z"/>
                </svg>
                Reaction Setup
              </h5>
              <div className="iso-protocol-steps">
                <div className="iso-step">
                  <span className="iso-step-num">1</span>
                  <div className="iso-step-content">
                    <span className="iso-step-title">Mix Fragments</span>
                    <span className="iso-step-detail">Combine {fragmentDesigns?.length || 0} PCR-amplified fragments at equimolar ratio (0.03-0.2 pmol each)</span>
                  </div>
                </div>
                <div className="iso-step">
                  <span className="iso-step-num">2</span>
                  <div className="iso-step-content">
                    <span className="iso-step-title">Add Master Mix</span>
                    <span className="iso-step-detail">Add 10 µL of 2X {variant === 'nebuilder_hifi' ? 'NEBuilder HiFi' : 'Gibson Assembly'} Master Mix</span>
                  </div>
                </div>
                <div className="iso-step">
                  <span className="iso-step-num">3</span>
                  <div className="iso-step-content">
                    <span className="iso-step-title">Add Water</span>
                    <span className="iso-step-detail">Add nuclease-free water to 20 µL total volume</span>
                  </div>
                </div>
                <div className="iso-step">
                  <span className="iso-step-num">4</span>
                  <div className="iso-step-content">
                    <span className="iso-step-title">Incubate</span>
                    <span className="iso-step-detail">
                      50°C for {(fragmentDesigns?.length || 0) <= 3 ? '15 minutes' : '60 minutes'}
                      {(fragmentDesigns?.length || 0) > 3 && ' (extended time for >3 fragments)'}
                    </span>
                  </div>
                </div>
                <div className="iso-step">
                  <span className="iso-step-num">5</span>
                  <div className="iso-step-content">
                    <span className="iso-step-title">Transform</span>
                    <span className="iso-step-detail">Transform 2 µL into competent E. coli cells</span>
                  </div>
                </div>
              </div>
            </div>

            <div className="iso-protocol-section">
              <h5>
                <svg viewBox="0 0 24 24" width="18" height="18" fill="#7c3aed">
                  <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
                </svg>
                Optimization Features Used
              </h5>
              <div className="iso-features-grid">
                <div className="iso-feature-item">
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="#22c55e">
                    <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
                  </svg>
                  <span>Range-based scoring</span>
                </div>
                <div className="iso-feature-item">
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="#22c55e">
                    <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
                  </svg>
                  <span>Junction position flexibility (±10bp)</span>
                </div>
                <div className="iso-feature-item">
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="#22c55e">
                    <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
                  </svg>
                  <span>Cross-junction secondary structure</span>
                </div>
                <div className="iso-feature-item">
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="#22c55e">
                    <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
                  </svg>
                  <span>Self-complementarity detection</span>
                </div>
                <div className="iso-feature-item">
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="#22c55e">
                    <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
                  </svg>
                  <span>G-quadruplex avoidance</span>
                </div>
                <div className="iso-feature-item">
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="#22c55e">
                    <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
                  </svg>
                  <span>Heterodimer analysis</span>
                </div>
              </div>
            </div>

            <div className="iso-protocol-section">
              <h5>
                <svg viewBox="0 0 24 24" width="18" height="18" fill="#7c3aed">
                  <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z"/>
                </svg>
                Tips for Success
              </h5>
              <div className="iso-tips-list">
                <div className="iso-tip">
                  <span className="iso-tip-bullet">•</span>
                  <span>Use high-fidelity polymerase (Q5, Phusion) for fragment amplification</span>
                </div>
                <div className="iso-tip">
                  <span className="iso-tip-bullet">•</span>
                  <span>Gel purify or clean up PCR products before assembly</span>
                </div>
                <div className="iso-tip">
                  <span className="iso-tip-bullet">•</span>
                  <span>Quantify fragments accurately using Nanodrop or Qubit</span>
                </div>
                <div className="iso-tip">
                  <span className="iso-tip-bullet">•</span>
                  <span>Store assembled product on ice before transformation</span>
                </div>
                <div className="iso-tip">
                  <span className="iso-tip-bullet">•</span>
                  <span>Use high-efficiency competent cells for complex assemblies</span>
                </div>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
