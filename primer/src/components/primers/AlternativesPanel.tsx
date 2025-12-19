/**
 * Unified Alternatives Panel Component
 *
 * A reusable component for displaying and comparing primer alternatives.
 * Used by both UnifiedPrimerDesigner and MutagenesisDesigner.
 *
 * Features:
 * - Card/Table view toggle
 * - Filtering (score, Tm diff, GC clamp)
 * - Sorting
 * - CSV export
 * - Multi-select comparison
 * - Trade-off analysis
 * - Diversity badges
 */

import React, { useState, useCallback, useMemo, useEffect, useRef } from 'react';
import {
  generateTradeOffs,
  filterAlternatives,
  sortAlternatives,
  getTmDiffColorClass,
  getQualityTier,
  formatSequence,
  calculateAnnealingTemp,
  generateAlternativesCSV,
  normalizeAlternative,
  TM_DIFF_THRESHOLDS,
} from '../../lib/primerAlternativesUtils.js';
import { identifyBadges, identifyStrengths, generateLabel, generateExplanation } from '../../lib/diversitySelection.js';
import { downloadFile } from '../../lib/sequenceUtils.js';
import { get3PrimeStructureBadge, classify3PrimeStructureSeverity } from '../../lib/smartPrimers.js';
import { foldSequence } from '../../lib/fold.js';

// =============================================================================
// Type Definitions
// =============================================================================

interface Primer {
  sequence?: string;
  tm?: number;
  length?: number;
  hasGCClamp?: boolean;
}

interface Design {
  forward?: Primer;
  reverse?: Primer;
  compositeScore?: number;
  score?: number;
}

interface Alternative extends Design {
  originalIdx: number;
  label?: string | { text: string; svgPath: string };
  explanation?: string;
  badges?: string[];
  tradeOffs?: TradeOff[];
  isBetterThanCurrent?: boolean;
}

interface AlternativeCategories {
  [key: string]: Alternative[];
}

interface TradeOff {
  type: string;
  label: string;
  delta: string;
  detail?: string;
}

interface Features {
  viewToggle?: boolean;
  filters?: boolean;
  sorting?: boolean;
  exportCSV?: boolean;
  compare?: boolean;
  expandSequences?: boolean;
  badges?: boolean;
  showCurrentDesign?: boolean;
}

interface CategoryConfig {
  svgPath: string;
  title: string;
  description: string;
  color: string;
}

interface FilterState {
  minScore: number;
  maxTmDiff: number;
  requireGcClamp: boolean;
  maxLength?: number;
}

interface FilterPreset {
  label: string;
  minScore: number;
  maxTmDiff: number;
  requireGcClamp: boolean;
  maxLength?: number;
}

interface SortConfig {
  field: string;
  direction: 'asc' | 'desc';
}

interface StructureSeverity {
  level: string;
  label: string;
  tooltip: string;
}

interface StructureInfo {
  fwd: { severity: StructureSeverity; energy: number } | null;
  rev: { severity: StructureSeverity; energy: number } | null;
  worst: { severity: StructureSeverity; energy: number } | null;
  worstPrimer?: string;
}

interface DesignInsight {
  label: string;
  detail: string;
}

interface DesignInsights {
  pros: DesignInsight[];
  cons: DesignInsight[];
}

// =============================================================================
// Category Configuration
// =============================================================================

/**
 * Category metadata for grouped view display
 * Each category explains WHY these alternatives were selected
 * Professional design with SVG icons for state-of-art primer designer
 */
const CATEGORY_CONFIG: Record<string, CategoryConfig> = {
  highestScore: {
    // Star icon SVG path
    svgPath: 'M12 2l3.09 6.26L22 9.27l-5 4.87 1.18 6.88L12 17.77l-6.18 3.25L7 14.14 2 9.27l6.91-1.01L12 2z',
    title: 'Highest Scoring',
    description: 'Best overall quality across all metrics',
    color: '#10b981', // green
  },
  shortestPrimers: {
    // Minimize/compress icon SVG path
    svgPath: 'M4 14h4v4H6v-2H4v-2zm0-4h2V8h2V6H4v4zm12 6h-2v2h4v-4h-2v2zm-2-6V8h2v2h2V6h-4v4z',
    title: 'Shortest Primers',
    description: 'Lower synthesis cost',
    color: '#3b82f6', // blue
  },
  longestPrimers: {
    // Maximize/expand icon SVG path
    svgPath: 'M21 11V3h-8l3.29 3.29-10 10L3 13v8h8l-3.29-3.29 10-10L21 11z',
    title: 'Longest Primers',
    description: 'Maximum specificity',
    color: '#8b5cf6', // purple
  },
  bestTmMatch: {
    // Target/crosshair icon SVG path
    svgPath: 'M12 2a10 10 0 100 20 10 10 0 000-20zm0 18a8 8 0 110-16 8 8 0 010 16zm0-14a6 6 0 100 12 6 6 0 000-12zm0 10a4 4 0 110-8 4 4 0 010 8zm0-6a2 2 0 100 4 2 2 0 000-4z',
    title: 'Best Tm Match',
    description: 'Most consistent annealing',
    color: '#f59e0b', // amber
  },
  fewestIssues: {
    // Shield check icon SVG path
    svgPath: 'M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10zm-1.5-5.5l-3-3 1.5-1.5 1.5 1.5 4-4 1.5 1.5-5.5 5.5z',
    title: 'Fewest Issues',
    description: 'Most reliable (no warnings)',
    color: '#06b6d4', // cyan
  },
};

/**
 * SVG Icon component for professional rendering
 */
interface SvgIconProps {
  path: string;
  size?: number;
  color?: string;
  className?: string;
}

function SvgIcon({ path, size = 16, color = 'currentColor', className = '' }: SvgIconProps) {
  return (
    <svg
      viewBox="0 0 24 24"
      width={size}
      height={size}
      fill={color}
      className={className}
      style={{ flexShrink: 0 }}
    >
      <path d={path} />
    </svg>
  );
}

/**
 * Label badge icons for alternative strengths
 */
const LABEL_ICONS: Record<string, string> = {
  bestOverall: 'M12 2l3.09 6.26L22 9.27l-5 4.87 1.18 6.88L12 17.77l-6.18 3.25L7 14.14 2 9.27l6.91-1.01L12 2z', // star
  bestTmMatch: 'M12 2a10 10 0 100 20 10 10 0 000-20zm0 18a8 8 0 110-16 8 8 0 010 16zm0-14a6 6 0 100 12 6 6 0 000-12zm0 10a4 4 0 110-8 4 4 0 010 8zm0-6a2 2 0 100 4 2 2 0 000-4z', // target
  safestDimer: 'M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10zm-1.5-5.5l-3-3 1.5-1.5 1.5 1.5 4-4 1.5 1.5-5.5 5.5z', // shield-check
  shortestAmplicon: 'M4 14h4v4H6v-2H4v-2zm0-4h2V8h2V6H4v4zm12 6h-2v2h4v-4h-2v2zm-2-6V8h2v2h2V6h-4v4z', // minimize
  cleanAmplicon: 'M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41L9 16.17z', // check
  optimalSize: 'M3 5v14a2 2 0 002 2h14a2 2 0 002-2V5a2 2 0 00-2-2H5a2 2 0 00-2 2zm16 0v2h-2V5h2zm-4 0v4h-2V5h2zm-4 0v2h-2V5h2zm-4 0v4H5V5h2z', // ruler
  compact: 'M4 14h4v4H6v-2H4v-2zm0-4h2V8h2V6H4v4zm12 6h-2v2h4v-4h-2v2zm-2-6V8h2v2h2V6h-4v4z', // minimize
  specific: 'M21 11V3h-8l3.29 3.29-10 10L3 13v8h8l-3.29-3.29 10-10L21 11z', // maximize/expand
};

// Order categories for display
const CATEGORY_ORDER = ['highestScore', 'fewestIssues', 'bestTmMatch', 'shortestPrimers', 'longestPrimers'];

// Filter presets for quick filtering
const FILTER_PRESETS: Record<string, FilterPreset> = {
  all: { label: 'All', minScore: 0, maxTmDiff: 10, requireGcClamp: false },
  highQuality: { label: 'High Quality', minScore: 70, maxTmDiff: 10, requireGcClamp: false },
  tmMatched: { label: 'Tm-Matched', minScore: 0, maxTmDiff: 2, requireGcClamp: false },
  shortPrimers: { label: 'Short', minScore: 0, maxTmDiff: 10, requireGcClamp: false, maxLength: 22 },
  reliable: { label: 'Reliable', minScore: 60, maxTmDiff: 3, requireGcClamp: true },
};

// =============================================================================
// Sub-components
// =============================================================================

/**
 * Score bar visualization
 */
interface ScoreBarProps {
  score: number;
  size?: 'normal' | 'mini';
}

function ScoreBar({ score, size = 'normal' }: ScoreBarProps) {
  const tier = getQualityTier(score);
  const height = size === 'mini' ? '4px' : '8px';

  return (
    <div className={`score-bar-container ${size}`} style={{ height }}>
      <div
        className={`score-bar ${tier}`}
        style={{ width: `${score}%` }}
      />
    </div>
  );
}

/**
 * Trade-off chip display
 */
interface TradeOffChipProps {
  tradeOff: TradeOff;
}

function TradeOffChip({ tradeOff }: TradeOffChipProps) {
  const { type, label, delta, detail } = tradeOff;

  return (
    <span
      className={`trade-off-chip trade-off-${type}`}
      title={detail ? `${label}: ${delta} (${detail})` : `${label}: ${delta}`}
    >
      <span className="chip-label">{label}</span>
      <span className="chip-value">{delta}</span>
    </span>
  );
}

/**
 * 3' Structure Badge for alternative designs
 * Shows severity-based badge with tooltip explaining the impact
 */
interface Structure3PrimeBadgeProps {
  design: Design;
}

function Structure3PrimeBadge({ design }: Structure3PrimeBadgeProps) {
  const [showTooltip, setShowTooltip] = useState(false);

  // Compute structure severity for both primers
  const structureInfo = useMemo((): StructureInfo => {
    const results: StructureInfo = { fwd: null, rev: null, worst: null };

    // Forward primer
    if (design.forward?.sequence) {
      try {
        const fold = foldSequence(design.forward.sequence);
        const severity = classify3PrimeStructureSeverity({
          energy: fold?.e ?? 0,
          basePairs: fold?.ij ?? [],
          seqLength: design.forward.sequence.length,
        } as any);
        results.fwd = { severity, energy: fold?.e ?? 0 };
      } catch (e) {
        results.fwd = null;
      }
    }

    // Reverse primer
    if (design.reverse?.sequence) {
      try {
        const fold = foldSequence(design.reverse.sequence);
        const severity = classify3PrimeStructureSeverity({
          energy: fold?.e ?? 0,
          basePairs: fold?.ij ?? [],
          seqLength: design.reverse.sequence.length,
        } as any);
        results.rev = { severity, energy: fold?.e ?? 0 };
      } catch (e) {
        results.rev = null;
      }
    }

    // Determine worst severity
    const severityOrder = ['none', 'info', 'low', 'moderate', 'warning', 'critical'];
    const fwdIdx = results.fwd ? severityOrder.indexOf(results.fwd.severity.level) : 0;
    const revIdx = results.rev ? severityOrder.indexOf(results.rev.severity.level) : 0;

    if (fwdIdx > revIdx && results.fwd) {
      results.worst = results.fwd;
      results.worstPrimer = 'forward';
    } else if (results.rev) {
      results.worst = results.rev;
      results.worstPrimer = 'reverse';
    }

    return results;
  }, [design.forward?.sequence, design.reverse?.sequence]);

  // Don't show badge if no significant structure
  if (!structureInfo.worst || structureInfo.worst.severity.level === 'none') {
    return null;
  }

  const { severity, energy } = structureInfo.worst;

  // Badge styling based on severity
  const badgeStyles: Record<string, { bg: string; text: string; border: string; icon: string }> = {
    critical: { bg: '#fee2e2', text: '#dc2626', border: '#fecaca', icon: '🔴' },
    warning: { bg: '#fef3c7', text: '#d97706', border: '#fde68a', icon: '⚠' },
    moderate: { bg: '#fef9c3', text: '#ca8a04', border: '#fef08a', icon: '△' },
    low: { bg: '#e0f2fe', text: '#0284c7', border: '#bae6fd', icon: 'ℹ' },
    info: { bg: '#f0f9ff', text: '#0369a1', border: '#e0f2fe', icon: '⚡' },
  };

  const style = badgeStyles[severity.level] || badgeStyles.info;

  return (
    <span
      style={{
        position: 'relative',
        display: 'inline-block',
      }}
      onMouseEnter={() => setShowTooltip(true)}
      onMouseLeave={() => setShowTooltip(false)}
    >
      <span
        style={{
          padding: '2px 8px',
          background: style.bg,
          color: style.text,
          border: `1px solid ${style.border}`,
          borderRadius: '4px',
          fontSize: '10px',
          fontWeight: '600',
          cursor: 'help',
          display: 'inline-flex',
          alignItems: 'center',
          gap: '3px',
          whiteSpace: 'nowrap',
        }}
        title={severity.tooltip}
      >
        {style.icon} 3' {severity.level === 'critical' ? 'blocked' : severity.level === 'warning' ? 'risk' : 'struct'}
      </span>

      {/* Tooltip */}
      {showTooltip && (
        <div
          style={{
            position: 'absolute',
            bottom: '100%',
            left: '50%',
            transform: 'translateX(-50%)',
            marginBottom: '8px',
            padding: '10px 12px',
            background: '#1f2937',
            color: '#f9fafb',
            borderRadius: '6px',
            fontSize: '11px',
            lineHeight: '1.5',
            whiteSpace: 'pre-wrap',
            width: '280px',
            boxShadow: '0 10px 25px rgba(0,0,0,0.3)',
            zIndex: 1000,
          }}
        >
          <div style={{ fontWeight: '600', marginBottom: '6px', color: style.text.replace('#', 'rgb(').replace(/(..)(..)(..)/, (m, r, g, b) => `${parseInt(r, 16)}, ${parseInt(g, 16)}, ${parseInt(b, 16)}`) + ')' }}>
            {severity.label}
          </div>
          <div style={{ fontSize: '10px', color: '#9ca3af' }}>
            {structureInfo.worstPrimer} primer: ΔG = {energy.toFixed(1)} kcal/mol
          </div>
          <div style={{ marginTop: '6px', fontSize: '10px' }}>
            {severity.level === 'critical'
              ? 'PCR will likely fail. Redesign required.'
              : severity.level === 'warning'
              ? 'May reduce efficiency. Consider optimization.'
              : 'Usually not problematic.'}
          </div>
          <div
            style={{
              position: 'absolute',
              bottom: '-6px',
              left: '50%',
              transform: 'translateX(-50%) rotate(45deg)',
              width: '12px',
              height: '12px',
              background: '#1f2937',
            }}
          />
        </div>
      )}
    </span>
  );
}

/**
 * ΔTm indicator with consistent styling
 */
interface TmDiffIndicatorProps {
  tmDiff: number;
  showLabel?: boolean;
  showWarning?: boolean;
}

function TmDiffIndicator({ tmDiff, showLabel = true, showWarning = true }: TmDiffIndicatorProps) {
  const colorClass = getTmDiffColorClass(tmDiff);

  return (
    <span className={`tm-diff ${colorClass}`}>
      {showLabel && <>ΔTm: {tmDiff.toFixed(1)}°C</>}
      {!showLabel && <>{tmDiff.toFixed(1)}°C</>}
      {showWarning && tmDiff > TM_DIFF_THRESHOLDS.GOOD && (
        <span className="warning-icon" title={`ΔTm > ${TM_DIFF_THRESHOLDS.GOOD}°C may affect annealing`}>!</span>
      )}
    </span>
  );
}

/**
 * Generate pros and cons for a design
 */
function generateDesignInsights(design: Design): DesignInsights {
  const pros: DesignInsight[] = [];
  const cons: DesignInsight[] = [];

  const score = design.compositeScore || design.score || 0;
  const tmDiff = Math.abs((design.forward?.tm || 0) - (design.reverse?.tm || 0));
  const fwdLen = design.forward?.length || 0;
  const revLen = design.reverse?.length || 0;
  const avgLen = (fwdLen + revLen) / 2;
  const gcClamps = (design.forward?.hasGCClamp ? 1 : 0) + (design.reverse?.hasGCClamp ? 1 : 0);

  // Score-based insights
  if (score >= 80) pros.push({ label: 'Excellent score', detail: `${score}/100` });
  else if (score >= 70) pros.push({ label: 'Good score', detail: `${score}/100` });
  else if (score < 60) cons.push({ label: 'Low score', detail: `${score}/100` });

  // Tm matching
  if (tmDiff <= 1.0) pros.push({ label: 'Excellent Tm match', detail: `ΔTm ${tmDiff.toFixed(1)}°C` });
  else if (tmDiff <= 1.5) pros.push({ label: 'Good Tm match', detail: `ΔTm ${tmDiff.toFixed(1)}°C` });
  else if (tmDiff > 3.0) cons.push({ label: 'High Tm difference', detail: `ΔTm ${tmDiff.toFixed(1)}°C` });

  // Length
  if (avgLen >= 18 && avgLen <= 22) pros.push({ label: 'Optimal length', detail: `${Math.round(avgLen)}bp avg` });
  else if (avgLen < 18) cons.push({ label: 'Short primers', detail: `${Math.round(avgLen)}bp avg` });
  else if (avgLen > 25) cons.push({ label: 'Long primers', detail: `${Math.round(avgLen)}bp avg` });

  // GC clamps
  if (gcClamps === 2) pros.push({ label: 'Both GC clamps', detail: '2/2' });
  else if (gcClamps === 0) cons.push({ label: 'No GC clamps', detail: '0/2' });

  // 3' Structure analysis - check for structure issues
  try {
    let worstSeverity = 'none';
    let worstEnergy = 0;

    if (design.forward?.sequence) {
      const fold = foldSequence(design.forward.sequence);
      const severity = classify3PrimeStructureSeverity({
        energy: fold?.e ?? 0,
        basePairs: fold?.ij ?? [],
        seqLength: design.forward.sequence.length,
      } as any);
      if (['critical', 'warning'].includes(severity.level) && severity.level !== 'none') {
        worstSeverity = severity.level;
        worstEnergy = fold?.e ?? 0;
      }
    }

    if (design.reverse?.sequence) {
      const fold = foldSequence(design.reverse.sequence);
      const severity = classify3PrimeStructureSeverity({
        energy: fold?.e ?? 0,
        basePairs: fold?.ij ?? [],
        seqLength: design.reverse.sequence.length,
      } as any);
      if (severity.level === 'critical' || (severity.level === 'warning' && worstSeverity !== 'critical')) {
        worstSeverity = severity.level;
        worstEnergy = fold?.e ?? 0;
      }
    }

    // Add to pros/cons based on 3' structure severity
    if (worstSeverity === 'critical') {
      cons.push({ label: "3' blocked", detail: `ΔG ${worstEnergy.toFixed(1)} - redesign needed` });
    } else if (worstSeverity === 'warning') {
      cons.push({ label: "3' structure risk", detail: `ΔG ${worstEnergy.toFixed(1)} - may affect efficiency` });
    } else if (worstSeverity === 'none') {
      // Only add as pro if no structure issues at all
      pros.push({ label: 'Free 3\' ends', detail: 'No blocking structures' });
    }
  } catch (e) {
    // Silently handle fold errors
  }

  return { pros, cons };
}

/**
 * Current/Original Design Card with Pros/Cons
 */
interface CurrentDesignCardProps {
  design: Design;
  isActive: boolean;
  isAlternativeSelected: boolean;
  isBestOverall: boolean;
  onSelect: () => void;
  onCopy: (fwd?: string, rev?: string) => void;
}

function CurrentDesignCard({
  design,
  isActive,
  isAlternativeSelected,
  isBestOverall,
  onSelect,
  onCopy,
}: CurrentDesignCardProps) {
  const score = design.compositeScore || design.score || 0;
  const tmDiff = Math.abs((design.forward?.tm || 0) - (design.reverse?.tm || 0));
  const { pros, cons } = generateDesignInsights(design);

  // Determine label based on context
  let labelText = 'CURRENT DESIGN';
  if (isAlternativeSelected) {
    labelText = 'ORIGINAL DESIGN';
  } else if (isBestOverall) {
    labelText = 'BEST DESIGN';
  }

  return (
    <div
      className={`current-design-card ${isActive ? 'active' : ''}`}
      onClick={isAlternativeSelected ? onSelect : undefined}
      style={{ cursor: isAlternativeSelected ? 'pointer' : 'default' }}
    >
      <div className="current-design-header">
        <div className="current-design-badge">
          <span className="star-icon">★</span>
          {labelText}
        </div>
        <div className="current-design-score">
          <ScoreBar score={score} />
          <span className="score-value">{score}/100</span>
        </div>
      </div>

      {/* Pros & Cons Badges */}
      <div className="design-insights">
        {pros.length > 0 && (
          <div className="insights-pros">
            {pros.map((pro, i) => (
              <span key={i} className="insight-badge pro" data-tooltip={pro.detail}>
                ✓ {pro.label}
              </span>
            ))}
          </div>
        )}
        {cons.length > 0 && (
          <div className="insights-cons">
            {cons.map((con, i) => (
              <span key={i} className="insight-badge con" data-tooltip={con.detail}>
                ✗ {con.label}
              </span>
            ))}
          </div>
        )}
      </div>

      <div className="current-design-meta">
        <span>{design.forward?.length || 0}bp + {design.reverse?.length || 0}bp</span>
        <span>|</span>
        <TmDiffIndicator tmDiff={tmDiff} showWarning={false} />
        <span>|</span>
        <span>
          Tm: {design.forward?.tm?.toFixed(1) || 0}°C / {design.reverse?.tm?.toFixed(1) || 0}°C
        </span>
        {/* 3' Structure Badge */}
        <Structure3PrimeBadge design={design} />
      </div>
      <div className="current-design-actions">
        <button
          type="button"
          className="action-btn"
          onClick={(e: React.MouseEvent<HTMLButtonElement>) => {
            e.stopPropagation();
            onCopy(design.forward?.sequence, design.reverse?.sequence);
          }}
          title="Copy both primers"
        >
          Copy
        </button>
        {isAlternativeSelected && (
          <span className="click-hint">Click to revert</span>
        )}
      </div>
    </div>
  );
}

/**
 * Why This Alternative? Tooltip
 */
interface WhyThisTooltipProps {
  alternative: Alternative;
  isVisible: boolean;
}

function WhyThisTooltip({ alternative, isVisible }: WhyThisTooltipProps) {
  if (!isVisible || !alternative.explanation) return null;

  return (
    <div className="why-this-tooltip">
      <div className="tooltip-header">Why this alternative?</div>
      <div className="tooltip-content">{alternative.explanation}</div>
      {alternative.badges && alternative.badges.length > 0 && (
        <div className="tooltip-badges">
          {alternative.badges.map((badge, i) => (
            <span key={i} className="tooltip-badge">{badge}</span>
          ))}
        </div>
      )}
    </div>
  );
}

/**
 * Alternative Card Component with Quick Actions on Hover
 */
interface AlternativeCardProps {
  alternative: Alternative;
  displayIndex: number;
  isSelected: boolean;
  isRecommended: boolean;
  isExpanded: boolean;
  isCompareSelected: boolean;
  isFocused: boolean;
  onSelect: () => void;
  onToggleExpand: () => void;
  onToggleCompare: () => void;
  onCopy: (fwd?: string, rev?: string) => void;
  onShowComparison: () => void;
  tradeOffs?: TradeOff[];
  copyFeedback: boolean;
}

function AlternativeCard({
  alternative,
  displayIndex,
  isSelected,
  isRecommended,
  isExpanded,
  isCompareSelected,
  isFocused,
  onSelect,
  onToggleExpand,
  onToggleCompare,
  onCopy,
  onShowComparison,
  tradeOffs,
  copyFeedback,
}: AlternativeCardProps) {
  const [showTooltip, setShowTooltip] = useState(false);
  const score = alternative.compositeScore || 0;
  const tmDiff = Math.abs((alternative.forward?.tm || 0) - (alternative.reverse?.tm || 0));

  return (
    <div
      className={`alt-card ${isSelected ? 'selected' : ''} ${isRecommended ? 'recommended' : ''} ${isFocused ? 'focused' : ''}`}
      tabIndex={0}
      data-index={alternative.originalIdx}
    >
      {/* Main content area with hover overlay */}
      <div className="alt-card-main">
        {/* Quick Actions Overlay - visible on hover */}
        <div className="alt-card-quick-actions">
          <button
            type="button"
            className="quick-action-btn primary"
            onClick={(e: React.MouseEvent<HTMLButtonElement>) => {
              e.stopPropagation();
              onSelect();
            }}
            title="Use this design (Enter)"
          >
            Use This
          </button>
          <button
            type="button"
            className="quick-action-btn"
            onClick={(e: React.MouseEvent<HTMLButtonElement>) => {
              e.stopPropagation();
              onCopy(alternative.forward?.sequence, alternative.reverse?.sequence);
            }}
            title="Copy primers (C)"
          >
            {copyFeedback ? '✓' : 'Copy'}
          </button>
          <button
            type="button"
            className={`quick-action-btn ${isCompareSelected ? 'active' : ''}`}
            onClick={(e: React.MouseEvent<HTMLButtonElement>) => {
              e.stopPropagation();
              onToggleCompare();
            }}
            title="Add to comparison (X)"
          >
            {isCompareSelected ? '✓ Compare' : 'Compare'}
          </button>
        </div>

        <div className="alt-card-header" onClick={onSelect}>
          <div className="alt-card-rank">
            <span className="rank-badge">#{displayIndex + 1}</span>
            {alternative.label && (
              <span
                className="diversity-label"
                onMouseEnter={() => setShowTooltip(true)}
                onMouseLeave={() => setShowTooltip(false)}
              >
                {/* Handle both old string labels and new object labels with SVG */}
                {typeof alternative.label === 'object' && alternative.label.svgPath ? (
                  <>
                    <SvgIcon path={alternative.label.svgPath} size={12} className="label-icon" />
                    <span className="label-text">{alternative.label.text}</span>
                  </>
                ) : (
                  alternative.label as React.ReactNode
                )}
                <WhyThisTooltip alternative={alternative} isVisible={showTooltip} />
              </span>
            )}
            {isSelected && <span className="status-badge selected">Selected</span>}
            {isRecommended && !isSelected && !alternative.label && (
              <span className="status-badge recommended">Better</span>
            )}
          </div>
          <div className="alt-card-score">
            <ScoreBar score={score} size="mini" />
            <span className="score-num">{score}</span>
          </div>
        </div>

        <div className="alt-card-meta" onClick={onSelect}>
          <span>{alternative.forward?.length}bp + {alternative.reverse?.length}bp</span>
          <span>|</span>
          <TmDiffIndicator tmDiff={tmDiff} />
          {/* 3' Structure Badge */}
          <Structure3PrimeBadge design={alternative} />
        </div>

        {tradeOffs && tradeOffs.length > 0 && (
          <div className="alt-card-tags">
            {tradeOffs.slice(0, 3).map((t, i) => (
              <TradeOffChip key={i} tradeOff={t} />
            ))}
            {tradeOffs.length > 3 && (
              <span
                className="trade-off-more"
                title={tradeOffs.slice(3).map(t => `${t.label}: ${t.delta}`).join(', ')}
              >
                +{tradeOffs.length - 3} more
              </span>
            )}
          </div>
        )}
      </div>

      {/* Collapsible Sequences - outside hover overlay */}
      <div className="alt-card-sequences">
        <button
          type="button"
          className="expand-seq-btn"
          onClick={(e: React.MouseEvent<HTMLButtonElement>) => {
            e.stopPropagation();
            onToggleExpand();
          }}
          aria-expanded={isExpanded}
          aria-controls={`seq-details-${alternative.originalIdx}`}
          aria-label={isExpanded ? 'Hide primer sequences' : 'Show primer sequences'}
        >
          {isExpanded ? '- Hide sequences' : '+ View sequences'}
        </button>
        {isExpanded && (
          <div className="seq-details" id={`seq-details-${alternative.originalIdx}`}>
            <div className="seq-row">
              <span className="seq-label">Fwd:</span>
              <code className="seq-code">{alternative.forward?.sequence}</code>
              <span className="seq-tm">{alternative.forward?.tm?.toFixed(1)}°C</span>
            </div>
            <div className="seq-row">
              <span className="seq-label">Rev:</span>
              <code className="seq-code">{alternative.reverse?.sequence}</code>
              <span className="seq-tm">{alternative.reverse?.tm?.toFixed(1)}°C</span>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

/**
 * Table Row for Alternative
 */
interface AlternativeTableRowProps {
  alternative: Alternative;
  displayIndex: number;
  isSelected: boolean;
  isCompareSelected: boolean;
  onSelect: () => void;
  onToggleCompare: () => void;
  sortConfig: SortConfig;
}

function AlternativeTableRow({
  alternative,
  displayIndex,
  isSelected,
  isCompareSelected,
  onSelect,
  onToggleCompare,
  sortConfig,
}: AlternativeTableRowProps) {
  const score = alternative.compositeScore || 0;
  const tmDiff = Math.abs((alternative.forward?.tm || 0) - (alternative.reverse?.tm || 0));

  return (
    <tr
      className={`${isSelected ? 'selected' : ''} ${alternative.isBetterThanCurrent ? 'recommended' : ''}`}
      onClick={onSelect}
    >
      <td className="col-rank">
        <input
          type="checkbox"
          checked={isCompareSelected}
          onChange={(e: React.ChangeEvent<HTMLInputElement>) => {
            e.stopPropagation();
            onToggleCompare();
          }}
          onClick={(e: React.MouseEvent<HTMLInputElement>) => e.stopPropagation()}
        />
        #{displayIndex + 1}
        {alternative.label && (
          <span className="label-badge">
            {typeof alternative.label === 'object' && alternative.label.svgPath ? (
              <>
                <SvgIcon path={alternative.label.svgPath} size={10} className="label-icon" />
                {alternative.label.text}
              </>
            ) : (
              alternative.label as React.ReactNode
            )}
          </span>
        )}
      </td>
      <td className="col-score">
        <ScoreBar score={score} size="mini" />
        <span>{score}</span>
      </td>
      <td className="col-fwd-seq" title={alternative.forward?.sequence}>
        <code>{formatSequence(alternative.forward?.sequence, 20)}</code>
      </td>
      <td className="col-fwd-tm">{alternative.forward?.tm?.toFixed(1)}°C</td>
      <td className="col-fwd-len">{alternative.forward?.length}bp</td>
      <td className="col-rev-seq" title={alternative.reverse?.sequence}>
        <code>{formatSequence(alternative.reverse?.sequence, 20)}</code>
      </td>
      <td className="col-rev-tm">{alternative.reverse?.tm?.toFixed(1)}°C</td>
      <td className="col-rev-len">{alternative.reverse?.length}bp</td>
      <td className="tm-cell">
        <TmDiffIndicator tmDiff={tmDiff} showLabel={false} showWarning={false} />
      </td>
      <td className="col-gc-clamp">
        {(alternative.forward?.hasGCClamp ? 1 : 0) + (alternative.reverse?.hasGCClamp ? 1 : 0)}/2
      </td>
    </tr>
  );
}

// =============================================================================
// Comparison Modal
// =============================================================================

/**
 * Side-by-side comparison modal for selected alternatives
 */
interface ComparisonModalProps {
  alternatives: Alternative[];
  currentDesign: Design;
  onClose: () => void;
  onSelect: (alternative: Alternative) => void;
}

function ComparisonModal({ alternatives, currentDesign, onClose, onSelect }: ComparisonModalProps) {
  if (!alternatives || alternatives.length === 0) return null;

  const getComparisonValue = (alt: Alternative, metric: string): number => {
    switch (metric) {
      case 'score': return alt.compositeScore || 0;
      case 'tmDiff': return Math.abs((alt.forward?.tm || 0) - (alt.reverse?.tm || 0));
      case 'fwdTm': return alt.forward?.tm || 0;
      case 'revTm': return alt.reverse?.tm || 0;
      case 'fwdLen': return alt.forward?.length || 0;
      case 'revLen': return alt.reverse?.length || 0;
      case 'gcClamps': return (alt.forward?.hasGCClamp ? 1 : 0) + (alt.reverse?.hasGCClamp ? 1 : 0);
      default: return 0;
    }
  };

  const getBestValue = (metric: string, lower = false): number => {
    const values = alternatives.map(alt => getComparisonValue(alt, metric));
    return lower ? Math.min(...values) : Math.max(...values);
  };

  const metrics: Array<{ key: string; label: string; unit: string; best: 'high' | 'low' | null }> = [
    { key: 'score', label: 'Score', unit: '', best: 'high' },
    { key: 'tmDiff', label: 'ΔTm', unit: '°C', best: 'low' },
    { key: 'fwdTm', label: 'Fwd Tm', unit: '°C', best: null },
    { key: 'revTm', label: 'Rev Tm', unit: '°C', best: null },
    { key: 'fwdLen', label: 'Fwd Length', unit: 'bp', best: null },
    { key: 'revLen', label: 'Rev Length', unit: 'bp', best: null },
    { key: 'gcClamps', label: 'GC Clamps', unit: '/2', best: 'high' },
  ];

  return (
    <div className="comparison-modal-overlay" onClick={onClose}>
      <div className="comparison-modal" onClick={(e: React.MouseEvent<HTMLDivElement>) => e.stopPropagation()}>
        <div className="comparison-header">
          <h3>Compare Alternatives</h3>
          <button type="button" className="close-btn" onClick={onClose}>×</button>
        </div>
        <div className="comparison-content">
          <table className="comparison-table">
            <thead>
              <tr>
                <th className="metric-col">Metric</th>
                {alternatives.map((alt, idx) => (
                  <th key={idx} className="alt-col">
                    <div className="alt-header">
                      <span className="alt-rank">#{alt.originalIdx + 1}</span>
                      {alt.label && <span className="alt-label">{typeof alt.label === 'object' && 'text' in alt.label ? alt.label.text : alt.label}</span>}
                    </div>
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {metrics.map(({ key, label, unit, best }) => {
                const bestVal = best ? getBestValue(key, best === 'low') : null;
                return (
                  <tr key={key}>
                    <td className="metric-label">{label}</td>
                    {alternatives.map((alt, idx) => {
                      const val = getComparisonValue(alt, key);
                      const isBest = best && val === bestVal;
                      return (
                        <td key={idx} className={`metric-value ${isBest ? 'best' : ''}`}>
                          {key === 'tmDiff' ? val.toFixed(1) : key.includes('Tm') ? val.toFixed(1) : val}{unit}
                        </td>
                      );
                    })}
                  </tr>
                );
              })}
              <tr className="sequence-row">
                <td className="metric-label">Forward</td>
                {alternatives.map((alt, idx) => (
                  <td key={idx} className="seq-cell">
                    <code>{alt.forward?.sequence}</code>
                  </td>
                ))}
              </tr>
              <tr className="sequence-row">
                <td className="metric-label">Reverse</td>
                {alternatives.map((alt, idx) => (
                  <td key={idx} className="seq-cell">
                    <code>{alt.reverse?.sequence}</code>
                  </td>
                ))}
              </tr>
            </tbody>
          </table>
        </div>
        <div className="comparison-footer">
          {alternatives.map((alt, idx) => (
            <button
              key={idx}
              type="button"
              className="select-btn"
              onClick={() => {
                onSelect(alt);
                onClose();
              }}
            >
              Use #{alt.originalIdx + 1}
            </button>
          ))}
        </div>
      </div>
    </div>
  );
}

// =============================================================================
// Main Component
// =============================================================================

/**
 * AlternativesPanel - Unified component for displaying primer alternatives
 */
interface AlternativesPanelProps {
  alternatives?: Alternative[];
  alternativeCategories?: AlternativeCategories | null;
  currentDesign?: Design;
  originalDesign?: Design | null;
  isAlternativeSelected?: boolean;
  selectedIndex?: number;
  onSelectAlternative: (alternative: Alternative) => void;
  onRevert: () => void;
  wasUpgraded?: boolean;
  originalScore?: number;
  features?: Features;
}

export default function AlternativesPanel({
  alternatives = [],
  alternativeCategories = null,
  currentDesign = {},
  originalDesign = null,
  isAlternativeSelected = false,
  selectedIndex = -1,
  onSelectAlternative,
  onRevert,
  wasUpgraded = false,
  originalScore = 0,
  features = {},
}: AlternativesPanelProps) {
  // Default features
  const {
    viewToggle = true,
    filters = true,
    sorting = true,
    exportCSV = true,
    compare = true,
    expandSequences = true,
    badges = true,
    showCurrentDesign = true,
  } = features;

  // Check if we have categorized data for grouped view
  const hasCategories = alternativeCategories && Object.keys(alternativeCategories).length > 0;

  // Refs
  const panelRef = useRef<HTMLDivElement>(null);

  // State - default to 'grouped' if categories available, otherwise 'card'
  const [viewMode, setViewMode] = useState<'card' | 'table' | 'grouped'>(hasCategories ? 'grouped' : 'card');
  const [showAll, setShowAll] = useState(false);
  const [expandedSequences, setExpandedSequences] = useState<Set<number>>(new Set());
  const [compareSelection, setCompareSelection] = useState<Set<number>>(new Set());
  const [dismissedBanner, setDismissedBanner] = useState(false);
  const [copiedFeedback, setCopiedFeedback] = useState<string | null>(null);
  const [sortConfig, setSortConfig] = useState<SortConfig>({ field: 'score', direction: 'desc' });
  const [filterState, setFilterState] = useState<FilterState>({
    minScore: 0,
    maxTmDiff: 10,
    requireGcClamp: false,
  });
  const [activePreset, setActivePreset] = useState('all');
  const [focusedIndex, setFocusedIndex] = useState(-1);
  const [showComparisonModal, setShowComparisonModal] = useState(false);

  // Normalize and process alternatives
  const processedAlternatives = useMemo(() => {
    const normalized = alternatives.map((alt, idx) => {
      const norm = normalizeAlternative(alt, idx);
      // Add badges and labels if enabled
      if (badges) {
        const strengths = identifyStrengths(norm, alternatives, idx);
        norm.label = norm.label || generateLabel(strengths);
        norm.explanation = norm.explanation || generateExplanation(norm, strengths);
        norm.badges = identifyBadges(norm);
      }
      return norm;
    });

    // Apply filters
    const filtered = filterAlternatives(normalized, filterState);

    // Apply sorting
    const sorted = sortAlternatives(filtered, sortConfig);

    // Add trade-offs - compare against the displayed design (originalDesign || currentDesign)
    const displayedDesign = originalDesign || currentDesign;
    const displayedScore = displayedDesign?.compositeScore || displayedDesign?.score || 0;
    const withTradeOffs = sorted.map(alt => ({
      ...alt,
      tradeOffs: generateTradeOffs(alt, displayedDesign),
      // Only mark as "better" if alternative score is strictly higher than displayed design
      isBetterThanCurrent: (alt.compositeScore || 0) > displayedScore,
    }));

    return withTradeOffs;
  }, [alternatives, filterState, sortConfig, currentDesign, originalDesign, badges]);

  const displayedAlternatives = showAll
    ? processedAlternatives
    : processedAlternatives.slice(0, 6);
  const hiddenCount = processedAlternatives.length - displayedAlternatives.length;
  const hasActiveFilters = filterState.minScore > 0 || filterState.maxTmDiff < 10 || filterState.requireGcClamp;

  // Handlers
  const handleCopy = useCallback(async (fwd?: string, rev?: string) => {
    const text = `Forward: ${fwd}\nReverse: ${rev}`;
    try {
      await navigator.clipboard.writeText(text);
      setCopiedFeedback('pair');
      setTimeout(() => setCopiedFeedback(null), 1500);
    } catch (err) {
      console.warn('Copy failed:', err);
    }
  }, []);

  const toggleExpand = useCallback((idx: number) => {
    setExpandedSequences(prev => {
      const next = new Set(prev);
      if (next.has(idx)) {
        next.delete(idx);
      } else {
        next.add(idx);
      }
      return next;
    });
  }, []);

  const toggleCompare = useCallback((idx: number) => {
    setCompareSelection(prev => {
      const next = new Set(prev);
      if (next.has(idx)) {
        next.delete(idx);
      } else {
        next.add(idx);
      }
      return next;
    });
  }, []);

  const handleExport = useCallback(() => {
    const toExport = compareSelection.size > 0
      ? processedAlternatives.filter(alt => compareSelection.has(alt.originalIdx))
      : processedAlternatives;

    const csv = generateAlternativesCSV(toExport);
    downloadFile(csv, 'alternative_primers.csv', 'text/csv');
  }, [processedAlternatives, compareSelection]);

  const handleSort = useCallback((field: string) => {
    setSortConfig(prev => ({
      field,
      direction: prev.field === field && prev.direction === 'desc' ? 'asc' : 'desc',
    }));
  }, []);

  const clearFilters = useCallback(() => {
    setFilterState({ minScore: 0, maxTmDiff: 10, requireGcClamp: false });
    setActivePreset('all');
  }, []);

  // Apply filter preset
  const applyPreset = useCallback((presetKey: string) => {
    const preset = FILTER_PRESETS[presetKey];
    if (preset) {
      setFilterState({
        minScore: preset.minScore,
        maxTmDiff: preset.maxTmDiff,
        requireGcClamp: preset.requireGcClamp,
      });
      setActivePreset(presetKey);
    }
  }, []);

  // Open comparison modal
  const openComparison = useCallback(() => {
    if (compareSelection.size >= 2) {
      setShowComparisonModal(true);
    }
  }, [compareSelection]);

  // Keyboard shortcuts
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      // Only handle when panel is focused or has focus within
      if (!panelRef.current?.contains(document.activeElement) && document.activeElement !== panelRef.current) {
        return;
      }

      // Don't handle if user is typing in an input
      if ((e.target as HTMLElement).tagName === 'INPUT' || (e.target as HTMLElement).tagName === 'SELECT' || (e.target as HTMLElement).tagName === 'TEXTAREA') {
        return;
      }

      const maxIndex = displayedAlternatives.length - 1;

      switch (e.key.toLowerCase()) {
        case 'j': // Next alternative
          e.preventDefault();
          setFocusedIndex(prev => Math.min(prev + 1, maxIndex));
          break;
        case 'k': // Previous alternative
          e.preventDefault();
          setFocusedIndex(prev => Math.max(prev - 1, 0));
          break;
        case 'enter': // Select focused alternative
          e.preventDefault();
          if (focusedIndex >= 0 && focusedIndex <= maxIndex) {
            onSelectAlternative(displayedAlternatives[focusedIndex]);
          }
          break;
        case 'c': // Copy focused alternative
          e.preventDefault();
          if (focusedIndex >= 0 && focusedIndex <= maxIndex) {
            const alt = displayedAlternatives[focusedIndex];
            handleCopy(alt.forward?.sequence, alt.reverse?.sequence);
          }
          break;
        case 'x': // Toggle compare for focused alternative
          e.preventDefault();
          if (focusedIndex >= 0 && focusedIndex <= maxIndex) {
            toggleCompare(displayedAlternatives[focusedIndex].originalIdx);
          }
          break;
        case 'escape': // Clear focus or close modal
          if (showComparisonModal) {
            setShowComparisonModal(false);
          } else {
            setFocusedIndex(-1);
          }
          break;
        default:
          break;
      }
    };

    document.addEventListener('keydown', handleKeyDown);
    return () => document.removeEventListener('keydown', handleKeyDown);
  }, [displayedAlternatives, focusedIndex, handleCopy, toggleCompare, onSelectAlternative, showComparisonModal]);

  // Get alternatives selected for comparison
  const comparisonAlternatives = useMemo(() => {
    return processedAlternatives.filter(alt => compareSelection.has(alt.originalIdx));
  }, [processedAlternatives, compareSelection]);

  // Compute whether current design is the best overall
  const isBestOverall = useMemo(() => {
    const displayedDesign = originalDesign || currentDesign;
    const displayedScore = displayedDesign?.compositeScore || displayedDesign?.score || 0;
    // Current design is best if its score is >= all alternative scores
    return alternatives.every(alt => {
      const altScore = alt.compositeScore || alt.score || 0;
      return displayedScore >= altScore;
    });
  }, [alternatives, currentDesign, originalDesign]);

  // No alternatives to show
  if (alternatives.length === 0) {
    return null;
  }

  return (
    <div className="alternatives-panel" ref={panelRef} tabIndex={-1}>
      {/* Comparison Modal */}
      {showComparisonModal && (
        <ComparisonModal
          alternatives={comparisonAlternatives}
          currentDesign={currentDesign}
          onClose={() => setShowComparisonModal(false)}
          onSelect={onSelectAlternative}
        />
      )}

      {/* Auto-optimization Banner */}
      {wasUpgraded && !dismissedBanner && (
        <div className="upgrade-notice">
          <span className="upgrade-icon">*</span>
          <span className="upgrade-text">
            Auto-optimized: {originalScore} &rarr; {currentDesign.compositeScore || currentDesign.score}
          </span>
          <button
            type="button"
            className="upgrade-info-btn"
            title="The design was automatically improved by optimizing primer 3' ends"
          >
            ?
          </button>
          <button
            type="button"
            className="upgrade-dismiss-btn"
            onClick={() => setDismissedBanner(true)}
            title="Dismiss"
          >
            x
          </button>
        </div>
      )}

      {/* Toolbar */}
      <div className="alternates-toolbar">
        {viewToggle && (
          <div className="view-toggle">
            {hasCategories && (
              <button
                type="button"
                className={`view-toggle-btn ${viewMode === 'grouped' ? 'active' : ''}`}
                onClick={() => setViewMode('grouped')}
                title="Grouped by Purpose"
              >
                By Purpose
              </button>
            )}
            <button
              type="button"
              className={`view-toggle-btn ${viewMode === 'card' ? 'active' : ''}`}
              onClick={() => setViewMode('card')}
              title="Card View"
            >
              Cards
            </button>
            <button
              type="button"
              className={`view-toggle-btn ${viewMode === 'table' ? 'active' : ''}`}
              onClick={() => setViewMode('table')}
              title="Table View"
            >
              Table
            </button>
          </div>
        )}

        <div className="toolbar-actions">
          {compare && compareSelection.size > 0 && (
            <>
              <span className="selection-count">{compareSelection.size} selected</span>
              {compareSelection.size >= 2 && (
                <button
                  type="button"
                  className="toolbar-btn compare"
                  onClick={openComparison}
                  title="Compare selected alternatives side-by-side"
                >
                  Compare
                </button>
              )}
            </>
          )}
          {exportCSV && (
            <button
              type="button"
              className="toolbar-btn"
              onClick={handleExport}
              title={compareSelection.size > 0 ? 'Export selected' : 'Export all'}
            >
              Export CSV
            </button>
          )}
          {isAlternativeSelected && (
            <button
              type="button"
              className="toolbar-btn revert"
              onClick={onRevert}
              title="Revert to original design"
            >
              Revert
            </button>
          )}
        </div>
      </div>

      {/* Filter Presets */}
      {filters && (
        <div className="filter-presets">
          <span className="presets-label">Quick filters:</span>
          {Object.entries(FILTER_PRESETS).map(([key, preset]) => (
            <button
              key={key}
              type="button"
              className={`preset-btn ${activePreset === key ? 'active' : ''}`}
              onClick={() => applyPreset(key)}
            >
              {preset.label}
            </button>
          ))}
        </div>
      )}

      {/* Filters */}
      {filters && (
        <div className="alternates-filters" role="group" aria-label="Filter alternatives">
          <div className="filter-group">
            <label htmlFor="filter-min-score">Min Score:</label>
            <select
              id="filter-min-score"
              value={filterState.minScore}
              onChange={(e: React.ChangeEvent<HTMLSelectElement>) => setFilterState(f => ({ ...f, minScore: parseInt(e.target.value) }))}
              className={filterState.minScore > 0 ? 'filter-active' : ''}
              aria-label="Filter by minimum score"
            >
              <option value="0">All</option>
              <option value="60">&gt;=60</option>
              <option value="70">&gt;=70</option>
              <option value="80">&gt;=80</option>
            </select>
          </div>
          <div className="filter-group">
            <label htmlFor="filter-max-tm">Max ΔTm:</label>
            <select
              id="filter-max-tm"
              value={filterState.maxTmDiff}
              onChange={(e: React.ChangeEvent<HTMLSelectElement>) => setFilterState(f => ({ ...f, maxTmDiff: parseFloat(e.target.value) }))}
              className={filterState.maxTmDiff < 10 ? 'filter-active' : ''}
              aria-label="Filter by maximum Tm difference"
            >
              <option value="10">All</option>
              <option value="5">&lt;=5°C</option>
              <option value="3">&lt;=3°C</option>
              <option value="2">&lt;=2°C</option>
            </select>
          </div>
          <div className="filter-group">
            <label className={`checkbox-label ${filterState.requireGcClamp ? 'filter-active' : ''}`}>
              <input
                type="checkbox"
                checked={filterState.requireGcClamp}
                onChange={(e: React.ChangeEvent<HTMLInputElement>) => setFilterState(f => ({ ...f, requireGcClamp: e.target.checked }))}
                aria-label="Require both primers to have GC clamp"
              />
              Both GC clamps
            </label>
          </div>
          {hasActiveFilters && (
            <div className="filter-feedback" role="status" aria-live="polite">
              <span className="filter-count">
                Showing {processedAlternatives.length} of {alternatives.length}
              </span>
              <button
                type="button"
                className="clear-filters-btn"
                onClick={clearFilters}
                aria-label="Clear all filters"
              >
                Clear
              </button>
            </div>
          )}
        </div>
      )}

      {/* Recommendation Callout */}
      {alternatives.length > 0 && !isAlternativeSelected && (
        <div className="recommendation-callout">
          <span className="recommendation-icon">💡</span>
          <span className="recommendation-text">
            {isBestOverall ? (
              <>
                <strong>Recommended:</strong> The current design has the highest overall score.
                Browse alternatives below if you need specific trade-offs like shorter primers or better Tm matching.
              </>
            ) : (
              <>
                <strong>Alternatives available:</strong> Some designs below have higher scores.
                Review the pros/cons to find the best fit for your needs.
              </>
            )}
          </span>
        </div>
      )}

      {/* Current Design Card */}
      {showCurrentDesign && (
        <CurrentDesignCard
          design={originalDesign || currentDesign}
          isActive={!isAlternativeSelected}
          isAlternativeSelected={isAlternativeSelected}
          isBestOverall={isBestOverall}
          onSelect={onRevert}
          onCopy={handleCopy}
        />
      )}

      {/* Grouped View - Shows alternatives organized by purpose */}
      {viewMode === 'grouped' && hasCategories && (
        <div className="alternates-grouped-view">
          {CATEGORY_ORDER.map(categoryKey => {
            const categoryAlts = alternativeCategories[categoryKey];
            if (!categoryAlts || categoryAlts.length === 0) return null;

            const config = CATEGORY_CONFIG[categoryKey];
            if (!config) return null;

            return (
              <div key={categoryKey} className="category-section">
                <div
                  className="category-header"
                  style={{ borderLeftColor: config.color }}
                >
                  <span className="category-icon">
                    <SvgIcon path={config.svgPath} size={18} color={config.color} />
                  </span>
                  <div className="category-info">
                    <span className="category-title">{config.title}</span>
                    <span className="category-description">{config.description}</span>
                  </div>
                  <span className="category-count">{categoryAlts.length}</span>
                </div>
                <div className="category-alternatives">
                  {categoryAlts.slice(0, 2).map((alt, idx) => {
                    const score = alt.compositeScore || 0;
                    const tmDiff = Math.abs((alt.forward?.tm || 0) - (alt.reverse?.tm || 0));
                    const totalLen = (alt.forward?.length || 0) + (alt.reverse?.length || 0);
                    const altKey = `${alt.forward?.sequence}-${alt.reverse?.sequence}`;

                    return (
                      <div
                        key={altKey}
                        className="category-alt-card"
                        onClick={() => onSelectAlternative(alt)}
                      >
                        <div className="cat-alt-header">
                          <span className="cat-alt-rank">#{idx + 1}</span>
                          <div className="cat-alt-score">
                            <ScoreBar score={score} size="mini" />
                            <span>{score}</span>
                          </div>
                        </div>
                        <div className="cat-alt-meta">
                          <span className="cat-alt-len">{totalLen}bp total</span>
                          <span className="cat-alt-tm">
                            <TmDiffIndicator tmDiff={tmDiff} showWarning={false} />
                          </span>
                          <Structure3PrimeBadge design={alt} />
                        </div>
                        <div className="cat-alt-sequences">
                          <div className="cat-seq-row">
                            <span className="cat-seq-label">F:</span>
                            <code className="cat-seq-code">{formatSequence(alt.forward?.sequence, 18)}</code>
                          </div>
                          <div className="cat-seq-row">
                            <span className="cat-seq-label">R:</span>
                            <code className="cat-seq-code">{formatSequence(alt.reverse?.sequence, 18)}</code>
                          </div>
                        </div>
                      </div>
                    );
                  })}
                </div>
              </div>
            );
          })}

          {/* Show flat list link */}
          <div className="grouped-view-footer">
            <button
              type="button"
              className="show-all-link"
              onClick={() => setViewMode('card')}
            >
              View all {alternatives.length} alternatives as cards
            </button>
          </div>
        </div>
      )}

      {/* Card View */}
      {viewMode === 'card' && (
        <div className="alternates-card-view">
          <div className="keyboard-hint">
            <span>Keyboard: J/K navigate, Enter select, C copy, X compare</span>
          </div>
          {displayedAlternatives.map((alt, displayIdx) => (
            <AlternativeCard
              key={alt.originalIdx}
              alternative={alt}
              displayIndex={displayIdx}
              isSelected={selectedIndex === alt.originalIdx && isAlternativeSelected}
              isRecommended={alt.isBetterThanCurrent || false}
              isExpanded={expandSequences && expandedSequences.has(alt.originalIdx)}
              isCompareSelected={compare && compareSelection.has(alt.originalIdx)}
              isFocused={focusedIndex === displayIdx}
              onSelect={() => onSelectAlternative(alt)}
              onToggleExpand={() => toggleExpand(alt.originalIdx)}
              onToggleCompare={() => toggleCompare(alt.originalIdx)}
              onCopy={handleCopy}
              onShowComparison={openComparison}
              tradeOffs={alt.tradeOffs}
              copyFeedback={copiedFeedback === 'pair'}
            />
          ))}

          {processedAlternatives.length === 0 && (
            <div className="empty-state">
              <span className="empty-icon">?</span>
              <p>No primers match the current filters.</p>
              <button
                type="button"
                className="clear-filters-btn"
                onClick={clearFilters}
              >
                Clear filters
              </button>
            </div>
          )}

          {hiddenCount > 0 && (
            <button
              type="button"
              className="show-more-btn"
              onClick={() => setShowAll(true)}
            >
              Show {hiddenCount} more alternatives
            </button>
          )}

          {showAll && processedAlternatives.length > 6 && (
            <button
              type="button"
              className="show-less-btn"
              onClick={() => setShowAll(false)}
            >
              Show less
            </button>
          )}
        </div>
      )}

      {/* Table View */}
      {viewMode === 'table' && (
        <div className="alternates-table-view">
          <table className="alternates-table">
            <thead>
              <tr>
                <th className="col-rank">#</th>
                <th
                  className={`col-score sortable ${sortConfig.field === 'score' ? 'sorted' : ''}`}
                  onClick={() => handleSort('score')}
                >
                  Score {sortConfig.field === 'score' ? (sortConfig.direction === 'desc' ? '↓' : '↑') : '⇅'}
                </th>
                <th className="col-fwd-seq">Fwd Sequence</th>
                <th
                  className={`col-fwd-tm sortable ${sortConfig.field === 'fwdTm' ? 'sorted' : ''}`}
                  onClick={() => handleSort('fwdTm')}
                >
                  Fwd Tm {sortConfig.field === 'fwdTm' ? (sortConfig.direction === 'desc' ? '↓' : '↑') : '⇅'}
                </th>
                <th
                  className={`col-fwd-len sortable ${sortConfig.field === 'fwdLen' ? 'sorted' : ''}`}
                  onClick={() => handleSort('fwdLen')}
                >
                  Len {sortConfig.field === 'fwdLen' ? (sortConfig.direction === 'desc' ? '↓' : '↑') : '⇅'}
                </th>
                <th className="col-rev-seq">Rev Sequence</th>
                <th
                  className={`col-rev-tm sortable ${sortConfig.field === 'revTm' ? 'sorted' : ''}`}
                  onClick={() => handleSort('revTm')}
                >
                  Rev Tm {sortConfig.field === 'revTm' ? (sortConfig.direction === 'desc' ? '↓' : '↑') : '⇅'}
                </th>
                <th
                  className={`col-rev-len sortable ${sortConfig.field === 'revLen' ? 'sorted' : ''}`}
                  onClick={() => handleSort('revLen')}
                >
                  Len {sortConfig.field === 'revLen' ? (sortConfig.direction === 'desc' ? '↓' : '↑') : '⇅'}
                </th>
                <th
                  className={`col-tm-diff sortable ${sortConfig.field === 'tmDiff' ? 'sorted' : ''}`}
                  onClick={() => handleSort('tmDiff')}
                >
                  ΔTm {sortConfig.field === 'tmDiff' ? (sortConfig.direction === 'desc' ? '↓' : '↑') : '⇅'}
                </th>
                <th className="col-gc-clamp">GC</th>
              </tr>
            </thead>
            <tbody>
              {displayedAlternatives.map((alt, displayIdx) => (
                <AlternativeTableRow
                  key={alt.originalIdx}
                  alternative={alt}
                  displayIndex={displayIdx}
                  isSelected={selectedIndex === alt.originalIdx && isAlternativeSelected}
                  isCompareSelected={compare && compareSelection.has(alt.originalIdx)}
                  onSelect={() => onSelectAlternative(alt)}
                  onToggleCompare={() => toggleCompare(alt.originalIdx)}
                  sortConfig={sortConfig}
                />
              ))}
            </tbody>
          </table>

          {processedAlternatives.length === 0 && (
            <div className="empty-state">
              <p>No primers match the current filters.</p>
              <button type="button" className="clear-filters-btn" onClick={clearFilters}>
                Clear filters
              </button>
            </div>
          )}

          {hiddenCount > 0 && (
            <button
              type="button"
              className="show-more-btn"
              onClick={() => setShowAll(true)}
            >
              Show {hiddenCount} more
            </button>
          )}
        </div>
      )}
    </div>
  );
}
