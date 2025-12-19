import React, { useState, useCallback, useMemo } from 'react';
import {
  designSequencingPrimers,
  designPrimerAtPosition,
  SEQUENCING_DEFAULTS,
} from '../lib/index.js';
import { translateDNA, reverseComplement } from '../lib/sequenceUtils.js';
import { getAvailablePresets, getPreset, autoDetectPreset } from '../lib/presets.js';

/**
 * Arrow component for primer visualization
 */
function PrimerArrow({ direction, color, isSelected }) {
  const arrowHeight = 20;
  const arrowWidth = 40;
  const headWidth = 12;

  if (direction === 'forward') {
    return (
      <svg width={arrowWidth} height={arrowHeight} viewBox={`0 0 ${arrowWidth} ${arrowHeight}`}>
        <polygon
          points={`0,4 ${arrowWidth - headWidth},4 ${arrowWidth - headWidth},0 ${arrowWidth},${arrowHeight / 2} ${arrowWidth - headWidth},${arrowHeight} ${arrowWidth - headWidth},${arrowHeight - 4} 0,${arrowHeight - 4}`}
          fill={color}
          stroke={isSelected ? '#1d4ed8' : 'none'}
          strokeWidth={isSelected ? 2 : 0}
        />
      </svg>
    );
  } else {
    return (
      <svg width={arrowWidth} height={arrowHeight} viewBox={`0 0 ${arrowWidth} ${arrowHeight}`}>
        <polygon
          points={`${arrowWidth},4 ${headWidth},4 ${headWidth},0 0,${arrowHeight / 2} ${headWidth},${arrowHeight} ${headWidth},${arrowHeight - 4} ${arrowWidth},${arrowHeight - 4}`}
          fill={color}
          stroke={isSelected ? '#1d4ed8' : 'none'}
          strokeWidth={isSelected ? 2 : 0}
        />
      </svg>
    );
  }
}

/**
 * Quality badge component with modern styling
 */
function QualityBadge({ quality }) {
  const colors = {
    excellent: { bg: '#dcfce7', text: '#166534', border: '#86efac' },
    good: { bg: '#fef9c3', text: '#854d0e', border: '#fde047' },
    acceptable: { bg: '#fed7aa', text: '#9a3412', border: '#fdba74' },
    poor: { bg: '#fecaca', text: '#991b1b', border: '#fca5a5' },
  };
  const style = colors[quality] || colors.acceptable;

  return (
    <span
      style={{
        display: 'inline-flex',
        alignItems: 'center',
        padding: '3px 10px',
        fontSize: '12px',
        fontWeight: '600',
        textTransform: 'capitalize',
        borderRadius: '12px',
        backgroundColor: style.bg,
        color: style.text,
        border: `1px solid ${style.border}`,
        whiteSpace: 'nowrap',
      }}
    >
      {quality}
    </span>
  );
}


/**
 * Modern Primer Card component - state-of-the-art design
 */
function PrimerCard({ primer, copyToClipboard, isExpanded, onToggleExpand }) {
  const hasAlternatives = primer.alternatives && primer.alternatives.length > 0;

  // Derive quality tier directly from composite score for consistency
  const getQualityFromScore = (score) => {
    if (score >= 90) return 'excellent';
    if (score >= 75) return 'good';
    if (score >= 60) return 'acceptable';
    return 'poor';
  };

  // Use score-based quality for display consistency
  const displayQuality = primer.compositeScore
    ? getQualityFromScore(primer.compositeScore)
    : (primer.qualityTier || primer.quality);

  // Filter warnings to exclude informational/positive messages
  const actualWarnings = (primer.warnings || []).filter(w => {
    const lower = w.toLowerCase();
    // Exclude positive/informational messages
    if (lower.includes('no g-quadruplex') || lower.includes('no risk')) return false;
    if (lower.includes('all primers have compatible')) return false;
    // "Strong 3' end" with 2-3 G/C is actually GOOD for sequencing (GC clamp)
    if (lower.includes('strong 3') && lower.includes('end')) return false;
    return true;
  });

  // Check if G-Quadruplex has an actual warning (not just info)
  const hasGQuadruplexWarning = primer.gQuadruplex?.severity === 'critical' ||
    primer.gQuadruplex?.severity === 'warning';

  // Determine if there are meaningful details to show
  const hasDetails = hasAlternatives ||
    actualWarnings.length > 0 ||
    (primer.issues && primer.issues.length > 0) ||
    hasGQuadruplexWarning ||
    primer.thermodynamics;

  const getScoreGradient = (score) => {
    if (score >= 90) return 'linear-gradient(90deg, #22c55e 0%, #4ade80 100%)';
    if (score >= 75) return 'linear-gradient(90deg, #84cc16 0%, #a3e635 100%)';
    if (score >= 60) return 'linear-gradient(90deg, #eab308 0%, #facc15 100%)';
    if (score >= 40) return 'linear-gradient(90deg, #f97316 0%, #fb923c 100%)';
    return 'linear-gradient(90deg, #ef4444 0%, #f87171 100%)';
  };

  return (
    <div style={{
      backgroundColor: '#fff',
      borderRadius: '12px',
      border: '1px solid #e2e8f0',
      marginBottom: '12px',
      overflow: 'hidden',
      boxShadow: '0 1px 3px rgba(0,0,0,0.05)',
      transition: 'box-shadow 0.2s, border-color 0.2s',
    }}>
      {/* Main card content */}
      <div style={{
        padding: '16px 20px',
        display: 'flex',
        flexDirection: 'column',
        gap: '12px',
      }}>
        {/* Header row: Name + Quality */}
        <div style={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          gap: '12px',
        }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
            <span style={{
              fontWeight: '600',
              fontSize: '15px',
              color: '#1e293b',
            }}>
              {primer.name}
            </span>
            {primer.isRescue && (
              <span style={{
                padding: '2px 8px',
                fontSize: '11px',
                fontWeight: '500',
                backgroundColor: '#fef3c7',
                color: '#92400e',
                borderRadius: '6px',
              }}>
                Rescue
              </span>
            )}
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
            <QualityBadge quality={displayQuality} />
          </div>
        </div>

        {/* Sequence row */}
        <div style={{
          display: 'flex',
          alignItems: 'center',
          gap: '12px',
          padding: '10px 14px',
          backgroundColor: '#f8fafc',
          borderRadius: '8px',
          border: '1px solid #e2e8f0',
        }}>
          <code style={{
            fontFamily: "'SF Mono', 'Monaco', 'Inconsolata', 'Roboto Mono', monospace",
            fontSize: '13px',
            color: '#334155',
            letterSpacing: '0.5px',
            flex: 1,
            wordBreak: 'break-all',
          }}>
            {primer.sequence}
          </code>
          <button
            type="button"
            onClick={() => copyToClipboard(primer.sequence)}
            style={{
              padding: '6px 14px',
              fontSize: '12px',
              fontWeight: '500',
              color: '#3b82f6',
              backgroundColor: '#eff6ff',
              border: '1px solid #bfdbfe',
              borderRadius: '6px',
              cursor: 'pointer',
              whiteSpace: 'nowrap',
              transition: 'all 0.15s',
            }}
          >
            Copy
          </button>
        </div>

        {/* Stats row */}
        <div style={{
          display: 'flex',
          alignItems: 'center',
          gap: '20px',
          flexWrap: 'wrap',
        }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
            <span style={{ color: '#64748b', fontSize: '13px' }}>Length:</span>
            <span style={{ color: '#1e293b', fontWeight: '500', fontSize: '13px' }}>{primer.length} bp</span>
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
            <span style={{ color: '#64748b', fontSize: '13px' }}>Tm:</span>
            <span style={{ color: '#1e293b', fontWeight: '500', fontSize: '13px' }}>{primer.tm}°C</span>
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
            <span style={{ color: '#64748b', fontSize: '13px' }}>GC:</span>
            <span style={{ color: '#1e293b', fontWeight: '500', fontSize: '13px' }}>{primer.gcPercent}</span>
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
            <span style={{ color: '#64748b', fontSize: '13px' }}>Position:</span>
            <span style={{ color: '#1e293b', fontWeight: '500', fontSize: '13px' }}>{primer.position} bp</span>
          </div>

          {/* Score bar */}
          <div style={{
            display: 'flex',
            alignItems: 'center',
            gap: '8px',
            marginLeft: 'auto',
          }}>
            <span style={{ color: '#64748b', fontSize: '13px' }}>Score:</span>
            <div style={{
              width: '80px',
              height: '8px',
              backgroundColor: '#e5e7eb',
              borderRadius: '4px',
              overflow: 'hidden',
            }}>
              <div style={{
                width: `${Math.min(100, primer.compositeScore || 0)}%`,
                height: '100%',
                background: getScoreGradient(primer.compositeScore),
                borderRadius: '4px',
                transition: 'width 0.3s ease',
              }} />
            </div>
            <span style={{
              color: '#1e293b',
              fontWeight: '600',
              fontSize: '13px',
              minWidth: '24px',
            }}>
              {primer.compositeScore?.toFixed(0) || '–'}
            </span>
          </div>
        </div>

        {/* Expandable details toggle */}
        {hasDetails && (
          <button
            type="button"
            onClick={onToggleExpand}
            style={{
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              gap: '6px',
              padding: '8px',
              backgroundColor: 'transparent',
              border: 'none',
              color: '#64748b',
              fontSize: '12px',
              cursor: 'pointer',
              borderTop: '1px solid #e2e8f0',
              marginTop: '4px',
              marginLeft: '-20px',
              marginRight: '-20px',
              marginBottom: '-16px',
            }}
          >
            {isExpanded ? '▲ Hide Details' : '▼ Show Details'}
            {hasAlternatives && ` (${primer.alternatives.length} alternatives)`}
          </button>
        )}
      </div>

      {/* Expanded details panel */}
      {isExpanded && (
        <div style={{
          padding: '16px 20px',
          backgroundColor: '#f8fafc',
          borderTop: '1px solid #e2e8f0',
        }}>
          {/* Critical Issues (shown first, in red) */}
          {primer.issues && primer.issues.length > 0 && (
            <div style={{ marginBottom: '16px' }}>
              <div style={{ fontWeight: '600', fontSize: '13px', color: '#991b1b', marginBottom: '8px' }}>
                Issues
              </div>
              {primer.issues.map((issue, i) => (
                <div key={i} style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: '8px',
                  padding: '8px 12px',
                  backgroundColor: '#fee2e2',
                  borderRadius: '6px',
                  marginBottom: '6px',
                  color: '#991b1b',
                  fontSize: '13px',
                }}>
                  <span>⚠</span>
                  <span>{issue}</span>
                </div>
              ))}
            </div>
          )}

          {/* Warnings (filtered to exclude informational messages) */}
          {actualWarnings.length > 0 && (
            <div style={{ marginBottom: '16px' }}>
              <div style={{ fontWeight: '600', fontSize: '13px', color: '#475569', marginBottom: '8px' }}>
                Warnings
              </div>
              {actualWarnings.map((w, i) => (
                <div key={i} style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: '8px',
                  padding: '8px 12px',
                  backgroundColor: '#fef3c7',
                  borderRadius: '6px',
                  marginBottom: '6px',
                  color: '#92400e',
                  fontSize: '13px',
                }}>
                  <span>⚠</span>
                  <span>{w}</span>
                </div>
              ))}
            </div>
          )}

          {/* G-Quadruplex - only show if there's an actual warning/critical */}
          {hasGQuadruplexWarning && (
            <div style={{
              display: 'flex',
              alignItems: 'center',
              gap: '8px',
              padding: '8px 12px',
              backgroundColor: primer.gQuadruplex.severity === 'critical' ? '#fee2e2' : '#fef3c7',
              borderRadius: '6px',
              marginBottom: '16px',
              color: primer.gQuadruplex.severity === 'critical' ? '#991b1b' : '#92400e',
              fontSize: '13px',
            }}>
              <span>{primer.gQuadruplex.severity === 'critical' ? '⚠' : '⚡'}</span>
              <span>G-Quadruplex: {primer.gQuadruplex.message}</span>
            </div>
          )}

          {/* Thermodynamics */}
          {primer.thermodynamics && (
            <div style={{ marginBottom: '16px' }}>
              <div style={{ fontWeight: '600', fontSize: '13px', color: '#475569', marginBottom: '8px' }}>
                Thermodynamics
              </div>
              <div style={{
                display: 'flex',
                gap: '16px',
                flexWrap: 'wrap',
                fontSize: '13px',
              }}>
                <span style={{ color: '#475569' }}>
                  Hairpin ΔG: <strong>{primer.thermodynamics.hairpinDG?.toFixed(1)} kcal/mol</strong>
                </span>
                <span style={{ color: '#475569' }}>
                  Self-dimer ΔG: <strong>{primer.thermodynamics.homodimerDG?.toFixed(1)} kcal/mol</strong>
                </span>
                <span style={{ color: '#475569' }}>
                  3′ ΔG: <strong>{primer.thermodynamics.terminal3DG?.toFixed(1)} kcal/mol</strong>
                </span>
              </div>
            </div>
          )}

          {/* Alternatives */}
          {hasAlternatives && (
            <div>
              <div style={{ fontWeight: '600', fontSize: '13px', color: '#475569', marginBottom: '10px' }}>
                Alternative Primers ({primer.alternatives.length})
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '8px' }}>
                {primer.alternatives.map((alt, i) => (
                  <div key={i} style={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    gap: '12px',
                    padding: '10px 14px',
                    backgroundColor: '#fff',
                    borderRadius: '8px',
                    border: '1px solid #e2e8f0',
                    fontSize: '13px',
                  }}>
                    <code style={{
                      fontFamily: "'SF Mono', 'Monaco', 'Inconsolata', monospace",
                      fontSize: '12px',
                      color: '#475569',
                      flex: 1,
                      wordBreak: 'break-all',
                    }}>
                      {alt.sequence}
                    </code>
                    <div style={{ display: 'flex', alignItems: 'center', gap: '12px', flexShrink: 0 }}>
                      <span style={{ color: '#64748b' }}>{alt.length}bp</span>
                      <span style={{ color: '#64748b' }}>{alt.tm}°C</span>
                      <span style={{ color: '#64748b' }}>{alt.gcPercent}</span>
                      <div style={{
                        width: '50px',
                        height: '6px',
                        backgroundColor: '#e5e7eb',
                        borderRadius: '3px',
                        overflow: 'hidden',
                      }}>
                        <div style={{
                          width: `${Math.min(100, alt.compositeScore || 0)}%`,
                          height: '100%',
                          background: getScoreGradient(alt.compositeScore),
                        }} />
                      </div>
                      <span style={{ fontWeight: '500', minWidth: '24px' }}>{alt.compositeScore?.toFixed(0)}</span>
                      <button
                        type="button"
                        onClick={() => copyToClipboard(alt.sequence)}
                        style={{
                          padding: '4px 10px',
                          fontSize: '11px',
                          color: '#3b82f6',
                          backgroundColor: '#eff6ff',
                          border: '1px solid #bfdbfe',
                          borderRadius: '4px',
                          cursor: 'pointer',
                        }}
                      >
                        Copy
                      </button>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

/**
 * Interactive Coverage Map Component
 * Shows primer positions and coverage regions on a linear map with arrow-style primers
 */
function CoverageMap({ results, onPrimerSelect, selectedPrimer }) {
  const { templateLength, primers, coverage, template } = results;
  const [hoveredRegion, setHoveredRegion] = useState(null);

  const mapWidth = 800;

  // Comprehensive 6-frame ORF detection
  // Forward frames: +1 (offset 0), +2 (offset 1), +3 (offset 2)
  // Reverse frames: -1, -2, -3 (reverse complement then offset 0, 1, 2)
  const { orfs } = useMemo(() => {
    if (!template) return { orfs: [] };

    const detectedOrfs = [];
    const minOrfLength = 90; // 30 codons minimum

    // Helper to detect ORFs in a translated protein sequence
    const findOrfsInFrame = (protein, frame, isReverse, seqLength) => {
      const starts = [];
      const stops = [];

      protein.split('').forEach((aa, idx) => {
        if (aa === 'M') starts.push(idx);
        if (aa === '*') stops.push(idx);
      });

      // Find ORFs from each start codon to nearest downstream stop
      const usedStops = new Set();
      for (const startIdx of starts) {
        const nextStop = stops.find(s => s > startIdx && !usedStops.has(s));
        if (nextStop) {
          const lengthInAA = nextStop - startIdx;
          const lengthInBP = lengthInAA * 3 + 3; // +3 for stop codon

          if (lengthInBP >= minOrfLength) {
            usedStops.add(nextStop);

            // Convert amino acid position to DNA position
            let dnaStart, dnaEnd;
            if (isReverse) {
              // For reverse strand, positions are from the 3' end
              const offset = Math.abs(frame) - 1;
              dnaEnd = seqLength - (startIdx * 3 + offset);
              dnaStart = seqLength - (nextStop * 3 + 3 + offset);
            } else {
              const offset = frame - 1;
              dnaStart = startIdx * 3 + offset;
              dnaEnd = nextStop * 3 + 3 + offset;
            }

            detectedOrfs.push({
              start: Math.max(0, dnaStart),
              end: Math.min(seqLength, dnaEnd),
              length: lengthInBP,
              proteinLength: lengthInAA,
              frame: isReverse ? -frame : frame,
              strand: isReverse ? '-' : '+',
            });
          }
        }
      }
    };

    // Forward frames (+1, +2, +3)
    for (let frame = 1; frame <= 3; frame++) {
      const protein = translateDNA(template, frame - 1);
      findOrfsInFrame(protein, frame, false, template.length);
    }

    // Reverse frames (-1, -2, -3)
    const revComp = reverseComplement(template);
    for (let frame = 1; frame <= 3; frame++) {
      const protein = translateDNA(revComp, frame - 1);
      findOrfsInFrame(protein, frame, true, template.length);
    }

    // Sort by length (longest first) and then by start position
    detectedOrfs.sort((a, b) => b.length - a.length || a.start - b.start);

    return { orfs: detectedOrfs };
  }, [template]);

  // Calculate coverage depth at each position for heatmap
  // Also track which primers cover each bin for interactivity
  const { coverageData, primersByBin } = useMemo(() => {
    const data = new Array(Math.ceil(templateLength / 10)).fill(0);
    const binPrimers = new Array(data.length).fill(null).map(() => []);

    for (const primer of primers) {
      const start = primer.readStart ?? 0;
      const end = primer.readEnd ?? templateLength;
      for (let i = Math.floor(start / 10); i < Math.ceil(end / 10) && i < data.length; i++) {
        data[i]++;
        binPrimers[i].push(primer);
      }
    }
    return { coverageData: data, primersByBin: binPrimers };
  }, [primers, templateLength]);

  // Generate ruler marks
  const rulerMarks = useMemo(() => {
    const marks = [];
    const step = templateLength <= 500 ? 100 : templateLength <= 2000 ? 200 : 500;
    for (let i = 0; i <= templateLength; i += step) {
      marks.push(i);
    }
    if (marks[marks.length - 1] !== templateLength) {
      marks.push(templateLength);
    }
    return marks;
  }, [templateLength]);

  const getCoverageColor = (depth) => {
    if (depth === 0) return '#fee2e2';
    if (depth === 1) return '#fef3c7';
    if (depth === 2) return '#d1fae5';
    return '#a7f3d0';
  };

  const getPrimerColor = (primer, isSelected) => {
    if (isSelected) return '#3b82f6';
    if (primer.isRescue) return '#f59e0b';
    if (primer.direction === 'forward') return '#22c55e';
    return '#8b5cf6';
  };

  const forwardPrimers = primers.filter(p => p.direction === 'forward');
  const reversePrimers = primers.filter(p => p.direction === 'reverse');

  return (
    <div className="coverage-map-container">
      <div className="coverage-map-header">
        <h4>Interactive Coverage Map</h4>
        <div className="map-legend">
          <span className="legend-item">
            <svg width="24" height="14" viewBox="0 0 24 14">
              <polygon points="0,3 16,3 16,0 24,7 16,14 16,11 0,11" fill="#22c55e" />
            </svg>
            Forward
          </span>
          <span className="legend-item">
            <svg width="24" height="14" viewBox="0 0 24 14">
              <polygon points="24,3 8,3 8,0 0,7 8,14 8,11 24,11" fill="#8b5cf6" />
            </svg>
            Reverse
          </span>
          <span className="legend-item">
            <svg width="24" height="14" viewBox="0 0 24 14">
              <polygon points="0,3 16,3 16,0 24,7 16,14 16,11 0,11" fill="#f59e0b" />
            </svg>
            Rescue
          </span>
          <span style={{ borderLeft: '1px solid #e2e8f0', height: '16px', margin: '0 8px' }} />
          <span className="legend-item">
            <span className="legend-color" style={{ backgroundColor: '#fee2e2' }}></span>
            0x
          </span>
          <span className="legend-item">
            <span className="legend-color" style={{ backgroundColor: '#fef3c7' }}></span>
            1x
          </span>
          <span className="legend-item">
            <span className="legend-color" style={{ backgroundColor: '#d1fae5' }}></span>
            2x+
          </span>
          <span style={{ borderLeft: '1px solid #e2e8f0', height: '16px', margin: '0 8px' }} />
          <span className="legend-item" style={{ display: 'flex', alignItems: 'center', gap: '4px' }}>
            <span style={{
              width: '16px',
              height: '5px',
              background: 'linear-gradient(90deg, #3b82f6 0%, #60a5fa 100%)',
              borderRadius: '2px',
            }} />
            ORF (+)
          </span>
          <span className="legend-item" style={{ display: 'flex', alignItems: 'center', gap: '4px' }}>
            <span style={{
              width: '16px',
              height: '5px',
              background: 'linear-gradient(90deg, #f59e0b 0%, #fbbf24 100%)',
              borderRadius: '2px',
            }} />
            ORF (−)
          </span>
        </div>
      </div>

      <div className="coverage-map" style={{ width: mapWidth + 60 }}>
        {/* Ruler */}
        <div className="map-ruler" style={{ width: mapWidth, marginLeft: 45 }}>
          {rulerMarks.map((mark) => (
            <div
              key={mark}
              className="ruler-mark"
              style={{ left: `${(mark / templateLength) * 100}%` }}
            >
              <div className="ruler-tick"></div>
              <span className="ruler-label">{mark}</span>
            </div>
          ))}
        </div>

        {/* Interactive Coverage Track */}
        <div style={{
          display: 'flex',
          alignItems: 'center',
          marginBottom: '8px',
        }}>
          <span style={{
            width: '45px',
            fontSize: '10px',
            color: '#64748b',
            textAlign: 'right',
            paddingRight: '6px',
            flexShrink: 0,
          }}>
            Cov
          </span>
          <div style={{
            width: mapWidth,
            position: 'relative',
          }}>
            {/* Coverage heatmap bar - interactive cells */}
            <div style={{
              height: '20px',
              position: 'relative',
              borderRadius: '4px',
              overflow: 'hidden',
              border: '1px solid #e2e8f0',
            }}>
              {coverageData.map((depth, idx) => {
                const binStart = idx * 10;
                const binEnd = Math.min((idx + 1) * 10, templateLength);
                const binPrimers = primersByBin[idx] || [];
                const isHovered = hoveredRegion === idx;
                const hasSelection = binPrimers.includes(selectedPrimer);

                return (
                  <div
                    key={idx}
                    style={{
                      position: 'absolute',
                      left: `${(binStart / templateLength) * 100}%`,
                      width: `${((binEnd - binStart) / templateLength) * 100}%`,
                      height: '100%',
                      backgroundColor: getCoverageColor(depth),
                      cursor: depth > 0 ? 'pointer' : 'default',
                      transition: 'filter 0.15s, transform 0.15s',
                      filter: isHovered ? 'brightness(0.85)' : hasSelection ? 'brightness(0.9)' : 'none',
                      boxShadow: isHovered ? 'inset 0 0 0 2px rgba(59, 130, 246, 0.5)' : 'none',
                    }}
                    title={`${binStart}-${binEnd} bp: ${depth}x coverage${depth > 0 ? `\nCovered by: ${binPrimers.map(p => p.name).join(', ')}` : ''}\nClick to highlight primers`}
                    onMouseEnter={() => setHoveredRegion(idx)}
                    onMouseLeave={() => setHoveredRegion(null)}
                    onClick={() => {
                      if (binPrimers.length > 0) {
                        // Select the first primer covering this region
                        onPrimerSelect(binPrimers[0]);
                      }
                    }}
                  />
                );
              })}
            </div>

            {/* Hovered region info tooltip */}
            {hoveredRegion !== null && primersByBin[hoveredRegion]?.length > 0 && (
              <div style={{
                position: 'absolute',
                left: `${((hoveredRegion * 10 + 5) / templateLength) * 100}%`,
                top: '-28px',
                transform: 'translateX(-50%)',
                backgroundColor: 'rgba(30, 41, 59, 0.95)',
                color: 'white',
                padding: '4px 8px',
                borderRadius: '4px',
                fontSize: '10px',
                whiteSpace: 'nowrap',
                pointerEvents: 'none',
                zIndex: 20,
              }}>
                {primersByBin[hoveredRegion].map(p => p.name).join(', ')}
              </div>
            )}
          </div>
        </div>

        {/* ORF Track - Separate bar for 6-frame ORFs */}
        {orfs.length > 0 && (
          <div style={{
            display: 'flex',
            alignItems: 'center',
            marginBottom: '8px',
          }}>
            <span style={{
              width: '45px',
              fontSize: '10px',
              color: '#64748b',
              textAlign: 'right',
              paddingRight: '6px',
              flexShrink: 0,
            }}>
              ORF
            </span>
            <div style={{
              width: mapWidth,
              position: 'relative',
              height: '28px',
              backgroundColor: '#f8fafc',
              borderRadius: '4px',
              border: '1px solid #e2e8f0',
            }}>
              {/* Forward strand ORFs (top half) */}
              {orfs.filter(orf => orf.strand === '+').map((orf, idx) => {
                const frameOffset = (orf.frame - 1) * 3; // Vertical offset based on frame
                return (
                  <div
                    key={`fwd-orf-${idx}`}
                    style={{
                      position: 'absolute',
                      left: `${(orf.start / templateLength) * 100}%`,
                      width: `${((orf.end - orf.start) / templateLength) * 100}%`,
                      top: `${2 + frameOffset}px`,
                      height: '6px',
                      background: 'linear-gradient(90deg, #3b82f6 0%, #60a5fa 100%)',
                      borderRadius: '3px',
                      cursor: 'help',
                      opacity: 0.9,
                    }}
                    title={`ORF +${orf.frame}: ${orf.start + 1}-${orf.end} (${orf.proteinLength} aa)`}
                  />
                );
              })}
              {/* Reverse strand ORFs (bottom half) */}
              {orfs.filter(orf => orf.strand === '-').map((orf, idx) => {
                const frameOffset = (Math.abs(orf.frame) - 1) * 3;
                return (
                  <div
                    key={`rev-orf-${idx}`}
                    style={{
                      position: 'absolute',
                      left: `${(orf.start / templateLength) * 100}%`,
                      width: `${((orf.end - orf.start) / templateLength) * 100}%`,
                      bottom: `${2 + frameOffset}px`,
                      height: '6px',
                      background: 'linear-gradient(90deg, #f59e0b 0%, #fbbf24 100%)',
                      borderRadius: '3px',
                      cursor: 'help',
                      opacity: 0.9,
                    }}
                    title={`ORF ${orf.frame}: ${orf.start + 1}-${orf.end} (${orf.proteinLength} aa)`}
                  />
                );
              })}
              {/* Center line separating forward/reverse */}
              <div style={{
                position: 'absolute',
                left: 0,
                right: 0,
                top: '50%',
                height: '1px',
                backgroundColor: '#e2e8f0',
                pointerEvents: 'none',
              }} />
            </div>
          </div>
        )}

        {/* Gap indicators */}
        {coverage.gaps.map((gap, idx) => (
          <div
            key={`gap-${idx}`}
            className="coverage-gap"
            style={{
              left: `${(gap.start / templateLength) * 100}%`,
              width: `${((gap.end - gap.start + 1) / templateLength) * 100}%`,
              marginLeft: 45,
            }}
            title={`Gap: ${gap.start}-${gap.end} (${gap.length} bp)`}
          >
            <span className="gap-label">Gap</span>
          </div>
        ))}

        {/* Forward primers track with arrows */}
        <div className="primer-track forward-track">
          <span className="track-label">Fwd</span>
          <div className="track-content" style={{ width: mapWidth }}>
            {forwardPrimers.map((primer, idx) => {
              const isSelected = selectedPrimer === primer;
              const leftPos = (primer.position / templateLength) * 100;
              return (
                <div
                  key={`fwd-${idx}`}
                  className={`primer-arrow-container ${isSelected ? 'selected' : ''}`}
                  style={{ left: `${leftPos}%` }}
                  onClick={() => onPrimerSelect(primer)}
                  title={`${primer.name}: ${primer.sequence}\nPos: ${primer.position}, Tm: ${primer.tm}°C`}
                >
                  <span className="primer-arrow-label">{primer.name}</span>
                  <PrimerArrow
                    direction="forward"
                    color={getPrimerColor(primer, isSelected)}
                    isSelected={isSelected}
                  />
                </div>
              );
            })}
          </div>
        </div>

        {/* Read coverage track (forward) */}
        <div className="read-track forward-read-track">
          <span className="track-label"></span>
          <div className="track-content" style={{ width: mapWidth }}>
            {forwardPrimers.map((primer, idx) => {
              const readStart = primer.readStart ?? primer.position;
              const readEnd = primer.readEnd ?? templateLength;
              const isSelected = selectedPrimer === primer;
              return (
                <div
                  key={`fwd-read-${idx}`}
                  className={`read-region ${isSelected ? 'selected' : ''}`}
                  style={{
                    left: `${(readStart / templateLength) * 100}%`,
                    width: `${((readEnd - readStart) / templateLength) * 100}%`,
                    backgroundColor: isSelected ? 'rgba(59, 130, 246, 0.25)' : 'rgba(34, 197, 94, 0.2)',
                    borderLeft: `2px solid ${isSelected ? '#3b82f6' : '#22c55e'}`,
                  }}
                  title={`${primer.name} read: ${readStart}-${readEnd} bp`}
                />
              );
            })}
          </div>
        </div>

        {/* Reverse primers track with arrows */}
        <div className="primer-track reverse-track">
          <span className="track-label">Rev</span>
          <div className="track-content" style={{ width: mapWidth }}>
            {reversePrimers.map((primer, idx) => {
              const isSelected = selectedPrimer === primer;
              const leftPos = (primer.position / templateLength) * 100;
              return (
                <div
                  key={`rev-${idx}`}
                  className={`primer-arrow-container reverse ${isSelected ? 'selected' : ''}`}
                  style={{ left: `${leftPos}%` }}
                  onClick={() => onPrimerSelect(primer)}
                  title={`${primer.name}: ${primer.sequence}\nPos: ${primer.position}, Tm: ${primer.tm}°C`}
                >
                  <span className="primer-arrow-label">{primer.name}</span>
                  <PrimerArrow
                    direction="reverse"
                    color={getPrimerColor(primer, isSelected)}
                    isSelected={isSelected}
                  />
                </div>
              );
            })}
          </div>
        </div>

        {/* Read coverage track (reverse) */}
        <div className="read-track reverse-read-track">
          <span className="track-label"></span>
          <div className="track-content" style={{ width: mapWidth }}>
            {reversePrimers.map((primer, idx) => {
              const readStart = primer.readStart ?? 0;
              const readEnd = primer.readEnd ?? primer.position;
              const isSelected = selectedPrimer === primer;
              return (
                <div
                  key={`rev-read-${idx}`}
                  className={`read-region ${isSelected ? 'selected' : ''}`}
                  style={{
                    left: `${(readStart / templateLength) * 100}%`,
                    width: `${((readEnd - readStart) / templateLength) * 100}%`,
                    backgroundColor: isSelected ? 'rgba(59, 130, 246, 0.25)' : 'rgba(139, 92, 246, 0.2)',
                    borderRight: `2px solid ${isSelected ? '#3b82f6' : '#8b5cf6'}`,
                  }}
                  title={`${primer.name} read: ${readStart}-${readEnd} bp`}
                />
              );
            })}
          </div>
        </div>

        {/* Template bar */}
        <div className="template-bar" style={{ width: mapWidth, marginLeft: 45 }}>
          <span className="bar-label">5'</span>
          <div className="bar-content"></div>
          <span className="bar-label">3'</span>
        </div>
      </div>

      {/* Selected primer details - improved card layout */}
      {selectedPrimer && (() => {
        // Derive quality tier from composite score for consistency
        const getQualityFromScore = (score) => {
          if (score >= 90) return 'excellent';
          if (score >= 75) return 'good';
          if (score >= 60) return 'acceptable';
          return 'poor';
        };
        const displayQuality = selectedPrimer.compositeScore
          ? getQualityFromScore(selectedPrimer.compositeScore)
          : (selectedPrimer.qualityTier || selectedPrimer.quality);

        // Filter warnings to exclude informational/positive messages
        const actualWarnings = (selectedPrimer.warnings || []).filter(w => {
          const lower = w.toLowerCase();
          if (lower.includes('no g-quadruplex') || lower.includes('no risk')) return false;
          if (lower.includes('all primers have compatible')) return false;
          // "Strong 3' end" with G/C is actually GOOD for sequencing (GC clamp)
          if (lower.includes('strong 3') && lower.includes('end')) return false;
          return true;
        });

        return (
        <div className="selected-primer-card">
          <div className="primer-card-header">
            <div className="primer-card-title">
              <span className={`direction-indicator ${selectedPrimer.direction}`}>
                {selectedPrimer.direction === 'forward' ? '→' : '←'}
              </span>
              <h5>{selectedPrimer.name}</h5>
              <QualityBadge quality={displayQuality} />
              {selectedPrimer.isRescue && <span className="rescue-tag">Rescue</span>}
            </div>
            <button
              className="close-btn"
              onClick={() => onPrimerSelect(null)}
              title="Close"
            >
              ×
            </button>
          </div>

          <div className="primer-card-body">
            <div className="sequence-display">
              <label>Sequence (5' → 3')</label>
              <div className="sequence-box">
                <code>{selectedPrimer.sequence}</code>
                <button
                  className="copy-btn"
                  onClick={() => navigator.clipboard.writeText(selectedPrimer.sequence)}
                  title="Copy sequence"
                >
                  Copy
                </button>
              </div>
            </div>

            <div className="primer-stats-grid">
              <div className="stat-item">
                <span className="stat-label">Position</span>
                <span className="stat-value">{selectedPrimer.position} bp</span>
              </div>
              <div className="stat-item">
                <span className="stat-label">Length</span>
                <span className="stat-value">{selectedPrimer.length} bp</span>
              </div>
              <div className="stat-item">
                <span className="stat-label">Tm</span>
                <span className="stat-value">{selectedPrimer.tm}°C</span>
              </div>
              <div className="stat-item">
                <span className="stat-label">GC Content</span>
                <span className="stat-value">{selectedPrimer.gcPercent}</span>
              </div>
              <div className="stat-item span-2">
                <span className="stat-label">Read Coverage</span>
                <span className="stat-value">{selectedPrimer.readStart} → {selectedPrimer.readEnd} bp</span>
              </div>
            </div>

            {actualWarnings.length > 0 && (
              <div className="primer-warnings-box">
                <span className="warnings-label">Warnings</span>
                <ul>
                  {actualWarnings.map((w, i) => (
                    <li key={i}>{w}</li>
                  ))}
                </ul>
              </div>
            )}
          </div>
        </div>
        );
      })()}
    </div>
  );
}

/**
 * Sequencing Primer Designer Component
 *
 * Designs primers optimized for Sanger sequencing with full coverage planning
 */
export default function SequencingDesigner() {
  const [template, setTemplate] = useState('');
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(false);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [selectedPrimer, setSelectedPrimer] = useState(null);
  const [options, setOptions] = useState({
    optimalReadLength: SEQUENCING_DEFAULTS.optimalReadLength,
    minPrimerLength: SEQUENCING_DEFAULTS.minPrimerLength,
    maxPrimerLength: SEQUENCING_DEFAULTS.maxPrimerLength,
    optimalTm: SEQUENCING_DEFAULTS.optimalTm,
    circular: false,
    generateAlternatives: true,
  });
  const [expandedAlternatives, setExpandedAlternatives] = useState({});
  const [selectedPreset, setSelectedPreset] = useState('default');

  // Position-specific primer design
  const [customPosition, setCustomPosition] = useState('');
  const [customDirection, setCustomDirection] = useState('forward');
  const [customPrimer, setCustomPrimer] = useState(null);

  const cleanSequence = useCallback((seq) => {
    return seq.toUpperCase().replace(/[^ATGC]/g, '');
  }, []);

  const handleDesign = useCallback(() => {
    const seq = cleanSequence(template);

    if (seq.length < 100) {
      setError('Template must be at least 100 bp for sequencing primer design');
      return;
    }

    setLoading(true);
    setError(null);
    setResults(null);

    // Use requestAnimationFrame + setTimeout to ensure UI updates before heavy computation
    // This prevents "page unresponsive" warnings from the browser
    requestAnimationFrame(() => {
      setTimeout(() => {
        try {
          const result = designSequencingPrimers(seq, options);
          setResults(result);
        } catch (err) {
          setError(err.message);
        } finally {
          setLoading(false);
        }
      }, 50); // Small delay to let the loading state render
    });
  }, [template, options, cleanSequence]);

  const handleCustomPrimer = useCallback(() => {
    const seq = cleanSequence(template);
    const pos = parseInt(customPosition);

    if (isNaN(pos) || pos < 0 || pos >= seq.length) {
      setError('Invalid position');
      return;
    }

    try {
      const primer = designPrimerAtPosition(seq, pos, customDirection, options);
      setCustomPrimer(primer);
      setError(null);
    } catch (err) {
      setError(err.message);
    }
  }, [template, customPosition, customDirection, options, cleanSequence]);

  const loadExample = useCallback(() => {
    // GFP gene sequence (~720 bp)
    setTemplate(`ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA`);
  }, []);

  const copyToClipboard = useCallback((text) => {
    navigator.clipboard.writeText(text);
  }, []);

  return (
    <div className="sequencing-designer">
      <div className="designer-header">
        <h2>Sequencing Primer Designer</h2>
        <p className="subtitle">
          Design primers optimized for Sanger sequencing with automatic coverage planning
        </p>
      </div>

      <div className="designer-content">
        <div className="input-section">
          <div className="form-group">
            <label htmlFor="template">
              Template Sequence
              <span className="sequence-length">
                {cleanSequence(template).length > 0 && ` (${cleanSequence(template).length} bp)`}
              </span>
            </label>
            <textarea
              id="template"
              value={template}
              onChange={(e) => setTemplate(e.target.value)}
              placeholder="Paste your DNA sequence here (ATGC only)..."
              rows={6}
              spellCheck={false}
            />
            <div className="input-actions">
              <button type="button" className="btn-secondary" onClick={loadExample}>
                Load GFP Example
              </button>
              <button type="button" className="btn-secondary" onClick={() => setTemplate('')}>
                Clear
              </button>
            </div>
          </div>

          <button
            type="button"
            className="toggle-advanced"
            onClick={() => setShowAdvanced(!showAdvanced)}
          >
            {showAdvanced ? '▼ Hide' : '▶ Show'} Advanced Options
          </button>

          {showAdvanced && (
            <div className="advanced-options">
              {/* Application Preset Selector */}
              <div className="preset-section" style={{ marginBottom: '16px', paddingBottom: '16px', borderBottom: '1px solid #e5e7eb' }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
                  <div className="form-group" style={{ flex: 1, marginBottom: 0 }}>
                    <label>Application Preset</label>
                    <select
                      value={selectedPreset}
                      onChange={(e) => {
                        const presetId = e.target.value;
                        setSelectedPreset(presetId);
                        const preset = getPreset(presetId);
                        if (preset && preset.optimalTmRange) {
                          setOptions(prev => ({
                            ...prev,
                            optimalTm: Math.round((preset.optimalTmRange[0] + preset.optimalTmRange[1]) / 2),
                          }));
                        }
                      }}
                      style={{ width: '100%' }}
                    >
                      {getAvailablePresets().map(preset => (
                        <option key={preset.id} value={preset.id}>
                          {preset.name}
                        </option>
                      ))}
                    </select>
                  </div>
                  <button
                    type="button"
                    onClick={() => {
                      const detected = autoDetectPreset(cleanSequence(template));
                      if (detected && detected !== 'default') {
                        setSelectedPreset(detected);
                      }
                    }}
                    disabled={cleanSequence(template).length < 50}
                    style={{
                      padding: '6px 12px',
                      fontSize: '13px',
                      background: '#f3f4f6',
                      border: '1px solid #d1d5db',
                      borderRadius: '4px',
                      cursor: cleanSequence(template).length >= 50 ? 'pointer' : 'not-allowed',
                      opacity: cleanSequence(template).length >= 50 ? 1 : 0.5,
                      marginTop: '20px',
                    }}
                  >
                    Auto-detect
                  </button>
                </div>
                <div style={{ fontSize: '12px', color: '#6b7280', marginTop: '4px' }}>
                  {getAvailablePresets().find(p => p.id === selectedPreset)?.description || ''}
                </div>
              </div>

              <div className="options-grid">
                <div className="form-group">
                  <label>Read Length (bp)</label>
                  <input
                    type="number"
                    value={options.optimalReadLength}
                    onChange={(e) => setOptions({ ...options, optimalReadLength: parseInt(e.target.value) })}
                    min={300}
                    max={1000}
                  />
                </div>
                <div className="form-group">
                  <label>Min Primer Length</label>
                  <input
                    type="number"
                    value={options.minPrimerLength}
                    onChange={(e) => setOptions({ ...options, minPrimerLength: parseInt(e.target.value) })}
                    min={15}
                    max={25}
                  />
                </div>
                <div className="form-group">
                  <label>Max Primer Length</label>
                  <input
                    type="number"
                    value={options.maxPrimerLength}
                    onChange={(e) => setOptions({ ...options, maxPrimerLength: parseInt(e.target.value) })}
                    min={20}
                    max={35}
                  />
                </div>
                <div className="form-group">
                  <label>Target Tm (°C)</label>
                  <input
                    type="number"
                    value={options.optimalTm}
                    onChange={(e) => setOptions({ ...options, optimalTm: parseInt(e.target.value) })}
                    min={45}
                    max={65}
                  />
                </div>
              </div>
              <div style={{ marginTop: '16px', display: 'flex', flexDirection: 'column', gap: '12px' }}>
                <label style={{ display: 'flex', alignItems: 'flex-start', gap: '8px', cursor: 'pointer' }}>
                  <input
                    type="checkbox"
                    checked={options.circular}
                    onChange={(e) => setOptions({ ...options, circular: e.target.checked })}
                    style={{ marginTop: '2px' }}
                  />
                  <div>
                    <div>Circular template (plasmid)</div>
                    <div style={{ fontSize: '12px', color: '#6b7280' }}>Enables primers that wrap around the origin</div>
                  </div>
                </label>
                <label style={{ display: 'flex', alignItems: 'flex-start', gap: '8px', cursor: 'pointer' }}>
                  <input
                    type="checkbox"
                    checked={options.generateAlternatives}
                    onChange={(e) => setOptions({ ...options, generateAlternatives: e.target.checked })}
                    style={{ marginTop: '2px' }}
                  />
                  <div>
                    <div>Generate alternative primers</div>
                    <div style={{ fontSize: '12px', color: '#6b7280' }}>Show up to 5 alternatives per position</div>
                  </div>
                </label>
              </div>
            </div>
          )}

          <div style={{ marginTop: '20px' }}>
            <button
              type="button"
              onClick={handleDesign}
              disabled={loading || cleanSequence(template).length < 100}
              style={{
                display: 'inline-flex',
                alignItems: 'center',
                justifyContent: 'center',
                gap: '10px',
                width: '100%',
                padding: '14px 28px',
                fontSize: '15px',
                fontWeight: '600',
                color: '#fff',
                background: loading
                  ? 'linear-gradient(135deg, #6366f1 0%, #8b5cf6 100%)'
                  : 'linear-gradient(135deg, #3b82f6 0%, #2563eb 100%)',
                border: 'none',
                borderRadius: '10px',
                cursor: loading || cleanSequence(template).length < 100 ? 'not-allowed' : 'pointer',
                opacity: cleanSequence(template).length < 100 ? 0.5 : 1,
                boxShadow: loading ? 'none' : '0 4px 14px rgba(59, 130, 246, 0.35)',
                transition: 'all 0.2s ease',
              }}
            >
              {loading && (
                <>
                  <div style={{
                    width: '18px',
                    height: '18px',
                    border: '2px solid rgba(255,255,255,0.3)',
                    borderTopColor: '#fff',
                    borderRadius: '50%',
                    animation: 'spin 0.8s linear infinite',
                  }} />
                  <span>Calculating Primers...</span>
                  <style>{`@keyframes spin { to { transform: rotate(360deg); } }`}</style>
                </>
              )}
              {!loading && (
                <>
                  <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                    <path d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2" />
                    <path d="M9 14l2 2 4-4" />
                  </svg>
                  <span>Design Sequencing Primers</span>
                </>
              )}
            </button>
          </div>
        </div>

        {error && (
          <div className="error-message">
            <span className="error-icon">!</span>
            {error}
          </div>
        )}

        {results && (
          <div className="results-section">
            <div className="results-summary">
              <h3>Coverage Summary</h3>
              {results.isCircular && (
                <div style={{
                  display: 'inline-flex',
                  alignItems: 'center',
                  gap: '6px',
                  padding: '4px 12px',
                  backgroundColor: '#dbeafe',
                  color: '#1e40af',
                  borderRadius: '16px',
                  fontSize: '13px',
                  marginBottom: '12px',
                }}>
                  <span>○</span>
                  Circular Template
                  {results.circularWrapped && ' (primers wrap origin)'}
                </div>
              )}
              <div className="summary-grid">
                <div className="summary-item">
                  <span className="label">Template Length</span>
                  <span className="value">{results.templateLength} bp</span>
                </div>
                <div className="summary-item">
                  <span className="label">Total Primers</span>
                  <span className="value">{results.primerCount}</span>
                </div>
                <div className="summary-item">
                  <span className="label">Coverage</span>
                  <span className="value highlight">{results.coverage.coveragePercent}</span>
                </div>
                <div className="summary-item">
                  <span className="label">Double Coverage</span>
                  <span className="value">{results.coverage.doubleCoveragePercent}</span>
                </div>
                {results.tmCompatibility && (
                  <div className="summary-item">
                    <span className="label">Tm Range</span>
                    <span className={`value ${results.tmCompatibility.compatible ? '' : 'warning'}`}>
                      {results.tmCompatibility.minTm}–{results.tmCompatibility.maxTm}°C
                    </span>
                  </div>
                )}
              </div>

              {results.coverage.gaps.length > 0 && (
                <div className="coverage-warning">
                  <span className="warning-icon">⚠</span>
                  {results.coverage.gaps.length} gap(s) in coverage detected
                </div>
              )}

              {results.tmCompatibility && !results.tmCompatibility.compatible && (
                <div style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: '8px',
                  padding: '8px 12px',
                  backgroundColor: '#fef3c7',
                  color: '#92400e',
                  borderRadius: '6px',
                  marginTop: '8px',
                  fontSize: '13px',
                }}>
                  <span>⚠</span>
                  {results.tmCompatibility.recommendation}
                </div>
              )}
            </div>

            {/* Interactive Coverage Map */}
            <CoverageMap
              results={results}
              onPrimerSelect={setSelectedPrimer}
              selectedPrimer={selectedPrimer}
            />

            {/* Designed Primers - Modern Card Layout */}
            <div style={{ marginTop: '32px' }}>
              <div style={{
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'space-between',
                marginBottom: '20px',
              }}>
                <h3 style={{
                  fontSize: '18px',
                  fontWeight: '600',
                  color: '#1e293b',
                  margin: 0,
                }}>
                  Designed Primers
                </h3>
                <span style={{
                  fontSize: '13px',
                  color: '#64748b',
                }}>
                  {results.primers.length} total primers
                </span>
              </div>

              {results.forwardPrimers.length > 0 && (
                <div style={{ marginBottom: '28px' }}>
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: '10px',
                    marginBottom: '14px',
                  }}>
                    <div style={{
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                      width: '28px',
                      height: '28px',
                      backgroundColor: '#dcfce7',
                      borderRadius: '6px',
                    }}>
                      <span style={{ fontSize: '14px' }}>→</span>
                    </div>
                    <h4 style={{
                      fontSize: '15px',
                      fontWeight: '600',
                      color: '#166534',
                      margin: 0,
                    }}>
                      Forward Primers ({results.forwardPrimers.length})
                    </h4>
                  </div>
                  {results.forwardPrimers.map((primer, idx) => (
                    <PrimerCard
                      key={idx}
                      primer={primer}
                      copyToClipboard={copyToClipboard}
                      isExpanded={expandedAlternatives[primer.name]}
                      onToggleExpand={() => setExpandedAlternatives(prev => ({
                        ...prev,
                        [primer.name]: !prev[primer.name]
                      }))}
                    />
                  ))}
                </div>
              )}

              {results.reversePrimers.length > 0 && (
                <div>
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: '10px',
                    marginBottom: '14px',
                  }}>
                    <div style={{
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                      width: '28px',
                      height: '28px',
                      backgroundColor: '#ede9fe',
                      borderRadius: '6px',
                    }}>
                      <span style={{ fontSize: '14px' }}>←</span>
                    </div>
                    <h4 style={{
                      fontSize: '15px',
                      fontWeight: '600',
                      color: '#5b21b6',
                      margin: 0,
                    }}>
                      Reverse Primers ({results.reversePrimers.length})
                    </h4>
                  </div>
                  {results.reversePrimers.map((primer, idx) => (
                    <PrimerCard
                      key={idx}
                      primer={primer}
                      copyToClipboard={copyToClipboard}
                      isExpanded={expandedAlternatives[primer.name]}
                      onToggleExpand={() => setExpandedAlternatives(prev => ({
                        ...prev,
                        [primer.name]: !prev[primer.name]
                      }))}
                    />
                  ))}
                </div>
              )}
            </div>

            {/* Export Section */}
            <div style={{
              marginTop: '28px',
              padding: '20px',
              backgroundColor: '#f8fafc',
              borderRadius: '12px',
              border: '1px solid #e2e8f0',
            }}>
              <h4 style={{
                fontSize: '14px',
                fontWeight: '600',
                color: '#475569',
                marginBottom: '14px',
              }}>
                Export All Primers
              </h4>
              <div style={{ display: 'flex', gap: '10px', flexWrap: 'wrap' }}>
                <button
                  type="button"
                  onClick={() => {
                    const text = results.primers
                      .map((p) => `${p.name}\t${p.sequence}`)
                      .join('\n');
                    copyToClipboard(text);
                  }}
                  style={{
                    display: 'inline-flex',
                    alignItems: 'center',
                    gap: '8px',
                    padding: '10px 18px',
                    fontSize: '13px',
                    fontWeight: '500',
                    color: '#475569',
                    backgroundColor: '#fff',
                    border: '1px solid #cbd5e1',
                    borderRadius: '8px',
                    cursor: 'pointer',
                    transition: 'all 0.15s',
                  }}
                >
                  <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                    <path d="M16 4h2a2 2 0 012 2v14a2 2 0 01-2 2H6a2 2 0 01-2-2V6a2 2 0 012-2h2" />
                    <rect x="8" y="2" width="8" height="4" rx="1" />
                  </svg>
                  Copy as TSV
                </button>
                <button
                  type="button"
                  onClick={() => {
                    const fasta = results.primers
                      .map((p) => `>${p.name}\n${p.sequence}`)
                      .join('\n');
                    copyToClipboard(fasta);
                  }}
                  style={{
                    display: 'inline-flex',
                    alignItems: 'center',
                    gap: '8px',
                    padding: '10px 18px',
                    fontSize: '13px',
                    fontWeight: '500',
                    color: '#475569',
                    backgroundColor: '#fff',
                    border: '1px solid #cbd5e1',
                    borderRadius: '8px',
                    cursor: 'pointer',
                    transition: 'all 0.15s',
                  }}
                >
                  <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                    <path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z" />
                    <path d="M14 2v6h6M16 13H8M16 17H8M10 9H8" />
                  </svg>
                  Copy as FASTA
                </button>
              </div>
            </div>
          </div>
        )}

        {/* Custom position primer design */}
        <div className="custom-primer-section">
          <h3>Design Primer at Specific Position</h3>
          <p className="section-description">
            Need a primer at a specific location? Enter the position below.
          </p>
          <div className="custom-primer-form">
            <div className="form-group">
              <label>Position (0-based)</label>
              <input
                type="number"
                value={customPosition}
                onChange={(e) => setCustomPosition(e.target.value)}
                min={0}
                max={cleanSequence(template).length - 1}
                placeholder="e.g., 500"
              />
            </div>
            <div className="form-group">
              <label>Direction</label>
              <select
                value={customDirection}
                onChange={(e) => setCustomDirection(e.target.value)}
              >
                <option value="forward">Forward</option>
                <option value="reverse">Reverse</option>
              </select>
            </div>
            <button
              type="button"
              className="btn-secondary"
              onClick={handleCustomPrimer}
              disabled={!customPosition || cleanSequence(template).length < 100}
            >
              Design
            </button>
          </div>

          {customPrimer && (
            <div style={{ marginTop: '16px' }}>
              <PrimerCard
                primer={{
                  ...customPrimer,
                  name: `Custom ${customPrimer.direction === 'forward' ? 'Fwd' : 'Rev'}`,
                }}
                copyToClipboard={copyToClipboard}
                isExpanded={false}
                onToggleExpand={() => {}}
              />
            </div>
          )}
        </div>
      </div>

      <div className="info-section">
        <h4>About Sequencing Primers</h4>
        <p>
          Sanger sequencing primers have different requirements than PCR primers:
        </p>
        <ul>
          <li>Lower Tm (50-55°C) for Big Dye chemistry</li>
          <li>Readable region starts ~50 bp after primer</li>
          <li>Typical read length is 700-800 bp</li>
          <li>Primers should have G/C at 3' end (GC clamp)</li>
        </ul>
      </div>
    </div>
  );
}
