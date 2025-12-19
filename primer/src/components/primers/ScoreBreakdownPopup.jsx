/**
 * ScoreBreakdownPopup - Interactive score breakdown modal
 *
 * Displays detailed piecewise scores with visual bars and
 * improvement suggestions to help users understand and optimize
 * their primer designs.
 */

import React, { useMemo } from 'react';

// Quality tier colors (matching existing design system)
const QUALITY_COLORS = {
  excellent: '#22c55e',
  good: '#3b82f6',
  acceptable: '#eab308',
  marginal: '#f97316',
  poor: '#ef4444',
};

// Score thresholds for color coding
const getScoreColor = (score) => {
  if (score >= 0.9) return '#22c55e'; // Excellent
  if (score >= 0.7) return '#3b82f6'; // Good
  if (score >= 0.5) return '#eab308'; // Acceptable
  if (score >= 0.3) return '#f97316'; // Marginal
  return '#ef4444'; // Poor
};

const getScoreLabel = (score) => {
  if (score >= 0.9) return 'Excellent';
  if (score >= 0.7) return 'Good';
  if (score >= 0.5) return 'Acceptable';
  if (score >= 0.3) return 'Marginal';
  return 'Poor';
};

// Professional SVG Icons
const Icons = {
  thermometer: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="M14 4v10.54a4 4 0 1 1-4 0V4a2 2 0 0 1 4 0Z" />
    </svg>
  ),
  barChart: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <line x1="12" y1="20" x2="12" y2="10" />
      <line x1="18" y1="20" x2="18" y2="4" />
      <line x1="6" y1="20" x2="6" y2="16" />
    </svg>
  ),
  ruler: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="M21.3 8.7 8.7 21.3c-1 1-2.5 1-3.4 0l-2.6-2.6c-1-1-1-2.5 0-3.4L15.3 2.7c1-1 2.5-1 3.4 0l2.6 2.6c1 1 1 2.5 0 3.4Z" />
      <path d="m7.5 10.5 2 2" />
      <path d="m10.5 7.5 2 2" />
      <path d="m13.5 4.5 2 2" />
      <path d="m4.5 13.5 2 2" />
    </svg>
  ),
  lock: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <rect width="18" height="11" x="3" y="11" rx="2" ry="2" />
      <path d="M7 11V7a5 5 0 0 1 10 0v4" />
    </svg>
  ),
  repeat: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="m17 2 4 4-4 4" />
      <path d="M3 11v-1a4 4 0 0 1 4-4h14" />
      <path d="m7 22-4-4 4-4" />
      <path d="M21 13v1a4 4 0 0 1-4 4H3" />
    </svg>
  ),
  paperclip: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="m21.44 11.05-9.19 9.19a6 6 0 0 1-8.49-8.49l8.57-8.57A4 4 0 1 1 18 8.84l-8.59 8.57a2 2 0 0 1-2.83-2.83l8.49-8.48" />
    </svg>
  ),
  link: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71" />
      <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71" />
    </svg>
  ),
  zap: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <polygon points="13 2 3 14 12 14 11 22 21 10 12 10 13 2" />
    </svg>
  ),
  target: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <circle cx="12" cy="12" r="10" />
      <circle cx="12" cy="12" r="6" />
      <circle cx="12" cy="12" r="2" />
    </svg>
  ),
  crosshair: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <circle cx="12" cy="12" r="10" />
      <line x1="22" y1="12" x2="18" y2="12" />
      <line x1="6" y1="12" x2="2" y2="12" />
      <line x1="12" y1="6" x2="12" y2="2" />
      <line x1="12" y1="22" x2="12" y2="18" />
    </svg>
  ),
  alertTriangle: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="m21.73 18-8-14a2 2 0 0 0-3.48 0l-8 14A2 2 0 0 0 4 21h16a2 2 0 0 0 1.73-3Z" />
      <line x1="12" y1="9" x2="12" y2="13" />
      <line x1="12" y1="17" x2="12.01" y2="17" />
    </svg>
  ),
  scale: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="m16 16 3-8 3 8c-.87.65-1.92 1-3 1s-2.13-.35-3-1Z" />
      <path d="m2 16 3-8 3 8c-.87.65-1.92 1-3 1s-2.13-.35-3-1Z" />
      <path d="M7 21h10" />
      <path d="M12 3v18" />
      <path d="M3 7h2c2 0 5-1 7-2 2 1 5 2 7 2h2" />
    </svg>
  ),
  dna: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="M2 15c6.667-6 13.333 0 20-6" />
      <path d="M9 22c1.798-1.998 2.518-3.995 2.807-5.993" />
      <path d="M15 2c-1.798 1.998-2.518 3.995-2.807 5.993" />
      <path d="m17 6-2.5-2.5" />
      <path d="m14 8-1-1" />
      <path d="m7 18 2.5 2.5" />
      <path d="m3.5 14.5.5.5" />
      <path d="m20 9 .5.5" />
      <path d="m6.5 12.5 1 1" />
      <path d="m16.5 10.5 1 1" />
      <path d="m10 16 1.5 1.5" />
    </svg>
  ),
  clipboard: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <rect width="8" height="4" x="8" y="2" rx="1" ry="1" />
      <path d="M16 4h2a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2H6a2 2 0 0 1-2-2V6a2 2 0 0 1 2-2h2" />
    </svg>
  ),
  lightbulb: (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="M15 14c.2-1 .7-1.7 1.5-2.5 1-.9 1.5-2.2 1.5-3.5A6 6 0 0 0 6 8c0 1 .2 2.2 1.5 3.5.7.7 1.3 1.5 1.5 2.5" />
      <path d="M9 18h6" />
      <path d="M10 22h4" />
    </svg>
  ),
};

// Feature metadata for display and suggestions
const FEATURE_INFO = {
  tm: {
    label: 'Melting Temperature',
    shortLabel: 'Tm',
    description: 'Optimal Tm is 55-65°C for efficient annealing',
    suggestion: 'Adjust primer length to shift Tm into optimal range (55-65°C)',
    icon: Icons.thermometer,
  },
  gc: {
    label: 'GC Content',
    shortLabel: 'GC%',
    description: 'Ideal GC content is 40-60%',
    suggestion: 'Shift primer position or adjust length to achieve 40-60% GC content',
    icon: Icons.barChart,
  },
  length: {
    label: 'Primer Length',
    shortLabel: 'Length',
    description: 'Optimal length is 18-25 bp',
    suggestion: 'Extend or shorten primer to 18-25 bp for optimal specificity',
    icon: Icons.ruler,
  },
  gcClamp: {
    label: 'GC Clamp',
    shortLabel: 'GC Clamp',
    description: '1-2 G/C bases in last 5 bp aids 3\' end stability',
    suggestion: 'Adjust 3\' end to include 1-2 G/C bases in the last 5 positions',
    icon: Icons.lock,
  },
  homopolymer: {
    label: 'Homopolymer',
    shortLabel: 'Runs',
    description: 'Avoid runs of 4+ identical bases',
    suggestion: 'Shift primer position to avoid long runs of identical bases',
    icon: Icons.repeat,
  },
  hairpin: {
    label: 'Hairpin Structure',
    shortLabel: 'Hairpin',
    description: 'Strong hairpins (ΔG < -2 kcal/mol) reduce efficiency',
    suggestion: 'Redesign to reduce self-complementary regions that form hairpins',
    icon: Icons.paperclip,
  },
  homodimer: {
    label: 'Self-Dimer',
    shortLabel: 'Self-Dimer',
    description: 'Self-dimerization reduces available primer',
    suggestion: 'Reduce self-complementary sequences, especially at 3\' end',
    icon: Icons.link,
  },
  heterodimer: {
    label: 'Heterodimer',
    shortLabel: 'Cross-Dimer',
    description: 'Primer-primer binding reduces efficiency',
    suggestion: 'Redesign primers to reduce complementarity between the pair',
    icon: Icons.zap,
  },
  terminal3DG: {
    label: '3\' Terminal Stability',
    shortLabel: '3\' ΔG',
    description: 'Optimal 3\' ΔG is -6 to -11 kcal/mol for efficient extension',
    suggestion: 'Adjust 3\' end sequence: too loose (> -6) needs more G/C, too tight (< -11) needs less',
    icon: Icons.target,
  },
  offTarget: {
    label: 'Off-Target Specificity',
    shortLabel: 'Specificity',
    description: 'Fewer off-target sites = higher specificity',
    suggestion: 'Extend primer or shift position to reduce off-target binding sites',
    icon: Icons.crosshair,
  },
  gQuadruplex: {
    label: 'G-Quadruplex Risk',
    shortLabel: 'G4 Risk',
    description: 'GGGG motifs can form problematic structures',
    suggestion: 'Avoid sequences with GGGG runs or G-rich regions',
    icon: Icons.alertTriangle,
  },
  tmDiff: {
    label: 'Tm Difference',
    shortLabel: 'ΔTm',
    description: 'Ideal Tm difference is < 2°C between primers',
    suggestion: 'Adjust one primer\'s length to match Tm values more closely',
    icon: Icons.scale,
  },
  threePrimeComp: {
    label: '3\' End Composition',
    shortLabel: '3\' Comp',
    description: 'Avoid 3\' end self-complementarity',
    suggestion: 'Redesign 3\' end to reduce self-complementary bases',
    icon: Icons.dna,
  },
};

// Internal weights for sorting (not displayed to user)
const FEATURE_WEIGHTS = {
  offTarget: 0.25,
  terminal3DG: 0.20,
  gQuadruplex: 0.10,
  heterodimer: 0.06,
  tm: 0.05,
  hairpin: 0.05,
  homodimer: 0.04,
  gc: 0.04,
  gcClamp: 0.03,
  tmDiff: 0.03,
  homopolymer: 0.02,
  length: 0.01,
  threePrimeComp: 0.04,
};

/**
 * Score bar component with gradient fill
 */
function ScoreBar({ score, color }) {
  return (
    <div style={{
      width: '100%',
      height: '6px',
      backgroundColor: '#e2e8f0',
      borderRadius: '3px',
      overflow: 'hidden',
    }}>
      <div style={{
        width: `${Math.max(0, Math.min(100, score * 100))}%`,
        height: '100%',
        backgroundColor: color,
        borderRadius: '3px',
        transition: 'width 0.3s ease',
      }} />
    </div>
  );
}

/**
 * Unchecked score row component - for features that weren't analyzed
 */
function UncheckedScoreRow({ feature, label: customLabel, recommendation }) {
  const info = FEATURE_INFO[feature] || {
    label: feature,
    shortLabel: feature,
    description: '',
    icon: Icons.clipboard,
  };

  return (
    <div style={{
      padding: '10px 12px',
      backgroundColor: '#f8fafc',
      borderRadius: '8px',
      marginBottom: '6px',
      border: '1px dashed #cbd5e1',
    }}>
      <div style={{
        display: 'flex',
        alignItems: 'center',
        gap: '10px',
        marginBottom: recommendation ? '8px' : '0',
      }}>
        <div style={{
          width: '28px',
          height: '28px',
          borderRadius: '6px',
          backgroundColor: '#f1f5f9',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          color: '#94a3b8',
          flexShrink: 0,
        }}>
          {info.icon}
        </div>
        <div style={{ flex: 1, minWidth: 0 }}>
          <div style={{
            fontWeight: '600',
            fontSize: '13px',
            color: '#64748b',
            whiteSpace: 'nowrap',
            overflow: 'hidden',
            textOverflow: 'ellipsis',
          }}>
            {customLabel || info.label}
          </div>
          <div style={{
            fontSize: '11px',
            color: '#94a3b8',
          }}>
            {info.description}
          </div>
        </div>
        <div style={{
          textAlign: 'right',
          flexShrink: 0,
        }}>
          <div style={{
            fontWeight: '600',
            fontSize: '12px',
            color: '#94a3b8',
            backgroundColor: '#f1f5f9',
            padding: '4px 8px',
            borderRadius: '4px',
          }}>
            Not checked
          </div>
        </div>
      </div>
      {recommendation && (
        <div style={{
          fontSize: '11px',
          color: '#64748b',
          backgroundColor: '#f1f5f9',
          padding: '6px 8px',
          borderRadius: '4px',
          display: 'flex',
          alignItems: 'center',
          gap: '6px',
        }}>
          <span style={{ color: '#3b82f6' }}>{Icons.lightbulb}</span>
          {recommendation}
        </div>
      )}
    </div>
  );
}

/**
 * Individual score row component
 */
function ScoreRow({ feature, score, label: customLabel }) {
  const info = FEATURE_INFO[feature] || {
    label: feature,
    shortLabel: feature,
    description: '',
    suggestion: '',
    icon: Icons.clipboard,
  };

  const color = getScoreColor(score);
  const statusLabel = getScoreLabel(score);
  const percentage = Math.round(score * 100);

  return (
    <div style={{
      padding: '10px 12px',
      backgroundColor: '#f8fafc',
      borderRadius: '8px',
      marginBottom: '6px',
      border: '1px solid #e2e8f0',
    }}>
      <div style={{
        display: 'flex',
        alignItems: 'center',
        gap: '10px',
        marginBottom: '6px',
      }}>
        <div style={{
          width: '28px',
          height: '28px',
          borderRadius: '6px',
          backgroundColor: `${color}15`,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          color: color,
          flexShrink: 0,
        }}>
          {info.icon}
        </div>
        <div style={{ flex: 1, minWidth: 0 }}>
          <div style={{
            fontWeight: '600',
            fontSize: '13px',
            color: '#1e293b',
            whiteSpace: 'nowrap',
            overflow: 'hidden',
            textOverflow: 'ellipsis',
          }}>
            {customLabel || info.label}
          </div>
          <div style={{
            fontSize: '11px',
            color: '#64748b',
            whiteSpace: 'nowrap',
            overflow: 'hidden',
            textOverflow: 'ellipsis',
          }}>
            {info.description}
          </div>
        </div>
        <div style={{
          textAlign: 'right',
          flexShrink: 0,
        }}>
          <div style={{
            fontWeight: '700',
            fontSize: '15px',
            color: color,
          }}>
            {percentage}%
          </div>
          <div style={{
            fontSize: '10px',
            color: color,
            fontWeight: '500',
          }}>
            {statusLabel}
          </div>
        </div>
      </div>
      <ScoreBar score={score} color={color} />
    </div>
  );
}

/**
 * Main ScoreBreakdownPopup component
 */
export default function ScoreBreakdownPopup({
  compositeScore,
  rawScore,
  criticalWarnings = 0,
  quality,
  forwardScores = {},
  reverseScores = {},
  hasTemplate = false,
  onClose,
}) {
  // Calculate if there's a penalty applied
  const hasPenalty = criticalWarnings > 0 && rawScore !== undefined && rawScore !== compositeScore;
  // Combine and organize scores for display
  const organizedScores = useMemo(() => {
    const scores = [];

    // Helper to add score if available
    // Skip offTarget scores if no template was provided (score would be meaningless default)
    const addScore = (key, value, label, weight) => {
      if (value !== undefined && value !== null && !isNaN(value)) {
        // Skip offTarget if no template - will show as "unchecked" separately
        if (key === 'offTarget' && !hasTemplate) {
          return;
        }
        scores.push({ key, value, label, weight });
      }
    };

    // Forward primer scores
    if (forwardScores) {
      addScore('terminal3DG', forwardScores.terminal3DG, 'Fwd 3\' Stability', 0.20);
      addScore('offTarget', forwardScores.offTarget, 'Fwd Specificity', 0.25);
      addScore('gQuadruplex', forwardScores.gQuadruplex, 'Fwd G-Quadruplex', 0.05);
      addScore('hairpin', forwardScores.hairpin, 'Fwd Hairpin', 0.02);
      addScore('homodimer', forwardScores.homodimer, 'Fwd Self-Dimer', 0.04);
      addScore('tm', forwardScores.tm, 'Fwd Tm', 0.02);
      addScore('gc', forwardScores.gc, 'Fwd GC%', 0.02);
      addScore('gcClamp', forwardScores.gcClamp, 'Fwd GC Clamp', 0.03);
      addScore('homopolymer', forwardScores.homopolymer, 'Fwd Homopolymer', 0.02);
      addScore('length', forwardScores.length, 'Fwd Length', 0.01);
    }

    // Reverse primer scores (if available)
    if (reverseScores && Object.keys(reverseScores).length > 0) {
      addScore('terminal3DG', reverseScores.terminal3DG, 'Rev 3\' Stability', 0.20);
      addScore('offTarget', reverseScores.offTarget, 'Rev Specificity', 0.25);
      addScore('gQuadruplex', reverseScores.gQuadruplex, 'Rev G-Quadruplex', 0.15);
      addScore('hairpin', reverseScores.hairpin, 'Rev Hairpin', 0.05);
      addScore('homodimer', reverseScores.homodimer, 'Rev Self-Dimer', 0.04);
      addScore('tm', reverseScores.tm, 'Rev Tm', 0.05);
      addScore('gc', reverseScores.gc, 'Rev GC%', 0.04);
      addScore('gcClamp', reverseScores.gcClamp, 'Rev GC Clamp', 0.03);
      addScore('homopolymer', reverseScores.homopolymer, 'Rev Homopolymer', 0.02);
      addScore('length', reverseScores.length, 'Rev Length', 0.01);

      // Pair-level scores
      addScore('heterodimer', reverseScores.heterodimer, 'Heterodimer', 0.06);
      addScore('tmDiff', reverseScores.tmDiff, 'Tm Difference', 0.03);
    }

    // Sort by weight (importance) descending
    return scores.sort((a, b) => b.weight - a.weight);
  }, [forwardScores, reverseScores, hasTemplate]);

  // Identify improvement opportunities (scores < 0.7)
  const improvements = useMemo(() => {
    return organizedScores
      .filter(s => s.value < 0.7)
      .sort((a, b) => (b.weight * (1 - b.value)) - (a.weight * (1 - a.value)))
      .slice(0, 3);
  }, [organizedScores]);

  const qualityTier = quality?.tier || quality || 'good';
  const qualityColor = QUALITY_COLORS[qualityTier] || QUALITY_COLORS.good;

  return (
    <div
      className="score-breakdown-overlay"
      onClick={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
      style={{
        position: 'fixed',
        top: 0,
        left: 0,
        right: 0,
        bottom: 0,
        backgroundColor: 'rgba(0, 0, 0, 0.5)',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        zIndex: 1000,
        padding: '20px',
      }}
    >
      <div
        className="score-breakdown-modal"
        style={{
          backgroundColor: 'white',
          borderRadius: '16px',
          boxShadow: '0 20px 40px rgba(0, 0, 0, 0.2)',
          maxWidth: '520px',
          width: '100%',
          maxHeight: '85vh',
          display: 'flex',
          flexDirection: 'column',
          overflow: 'hidden',
        }}
      >
        {/* Header */}
        <div style={{
          padding: '20px 24px',
          borderBottom: '1px solid #e2e8f0',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
        }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '16px' }}>
            {/* Quality Ring */}
            <div style={{
              position: 'relative',
              width: '56px',
              height: '56px',
            }}>
              <svg
                viewBox="0 0 56 56"
                style={{
                  transform: 'rotate(-90deg)',
                  width: '100%',
                  height: '100%',
                }}
              >
                <circle
                  cx="28"
                  cy="28"
                  r="24"
                  fill="none"
                  stroke="#e2e8f0"
                  strokeWidth="5"
                />
                <circle
                  cx="28"
                  cy="28"
                  r="24"
                  fill="none"
                  stroke={qualityColor}
                  strokeWidth="5"
                  strokeLinecap="round"
                  strokeDasharray={`${(compositeScore / 100) * 151} 151`}
                  style={{ transition: 'stroke-dasharray 0.5s ease' }}
                />
              </svg>
              <div style={{
                position: 'absolute',
                top: '50%',
                left: '50%',
                transform: 'translate(-50%, -50%)',
                textAlign: 'center',
              }}>
                <div style={{
                  fontSize: '16px',
                  fontWeight: '700',
                  color: qualityColor,
                }}>
                  {compositeScore}
                </div>
              </div>
            </div>

            <div>
              <h2 style={{ margin: 0, fontSize: '16px', fontWeight: '600', color: '#1e293b' }}>
                Score Breakdown
              </h2>
              <div style={{
                fontSize: '13px',
                color: qualityColor,
                fontWeight: '600',
                textTransform: 'capitalize',
                marginTop: '2px',
              }}>
                {quality?.label || qualityTier} Quality
              </div>
              {hasPenalty && (
                <div style={{
                  fontSize: '11px',
                  color: '#ef4444',
                  marginTop: '4px',
                }}>
                  Raw: {rawScore} − {criticalWarnings * 20} penalty = {compositeScore}
                </div>
              )}
            </div>
          </div>

          <button
            onClick={onClose}
            style={{
              width: '32px',
              height: '32px',
              borderRadius: '8px',
              border: 'none',
              backgroundColor: '#f1f5f9',
              cursor: 'pointer',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              fontSize: '18px',
              color: '#64748b',
              transition: 'background-color 0.2s',
            }}
            onMouseEnter={(e) => e.currentTarget.style.backgroundColor = '#e2e8f0'}
            onMouseLeave={(e) => e.currentTarget.style.backgroundColor = '#f1f5f9'}
          >
            ×
          </button>
        </div>

        {/* Content */}
        <div style={{
          flex: 1,
          overflowY: 'auto',
          padding: '16px 20px',
        }}>
          {/* Top Improvement Opportunities */}
          {improvements.length > 0 && (
            <div style={{
              marginBottom: '16px',
              padding: '14px',
              backgroundColor: '#fffbeb',
              borderRadius: '10px',
              border: '1px solid #fde68a',
            }}>
              <h3 style={{
                margin: '0 0 10px 0',
                fontSize: '13px',
                fontWeight: '600',
                color: '#92400e',
                display: 'flex',
                alignItems: 'center',
                gap: '8px',
              }}>
                <span style={{ color: '#f59e0b' }}>{Icons.lightbulb}</span>
                Improvement Opportunities
              </h3>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '6px' }}>
                {improvements.map((item, index) => {
                  const info = FEATURE_INFO[item.key] || {};
                  return (
                    <div key={index} style={{
                      padding: '10px 12px',
                      backgroundColor: 'white',
                      borderRadius: '6px',
                      fontSize: '12px',
                      border: '1px solid #fde68a',
                    }}>
                      <div style={{
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'space-between',
                        marginBottom: '4px',
                      }}>
                        <span style={{
                          fontWeight: '600',
                          color: '#1e293b',
                        }}>
                          {item.label}
                        </span>
                        <span style={{
                          fontWeight: '600',
                          color: getScoreColor(item.value),
                        }}>
                          {Math.round(item.value * 100)}%
                        </span>
                      </div>
                      <div style={{ color: '#64748b', lineHeight: '1.4' }}>
                        {info.suggestion || 'Consider optimizing this parameter'}
                      </div>
                    </div>
                  );
                })}
              </div>
            </div>
          )}

          {/* All Scores */}
          <div>
            <div style={{
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'space-between',
              marginBottom: '10px',
            }}>
              <h3 style={{
                margin: 0,
                fontSize: '13px',
                fontWeight: '600',
                color: '#334155',
              }}>
                All Scores
              </h3>
              <span style={{
                fontSize: '11px',
                color: '#94a3b8',
                fontStyle: 'italic',
              }}>
                Ranked by importance
              </span>
            </div>

            {/* Show unchecked specificity at top when no template provided */}
            {!hasTemplate && (
              <UncheckedScoreRow
                feature="offTarget"
                label="Specificity"
                recommendation="Add a template sequence to check for off-target binding sites"
              />
            )}

            {organizedScores.map((item, index) => (
              <ScoreRow
                key={`${item.key}-${item.label}-${index}`}
                feature={item.key}
                score={item.value}
                label={item.label}
              />
            ))}

            {organizedScores.length === 0 && hasTemplate && (
              <div style={{
                padding: '40px',
                textAlign: 'center',
                color: '#64748b',
              }}>
                No detailed scores available
              </div>
            )}
          </div>
        </div>

        {/* Footer */}
        <div style={{
          padding: '14px 20px',
          borderTop: '1px solid #e2e8f0',
          backgroundColor: '#f8fafc',
          display: 'flex',
          justifyContent: 'flex-end',
        }}>
          <button
            onClick={onClose}
            style={{
              padding: '8px 20px',
              backgroundColor: '#3b82f6',
              color: 'white',
              border: 'none',
              borderRadius: '8px',
              cursor: 'pointer',
              fontWeight: '500',
              fontSize: '13px',
              transition: 'background-color 0.2s',
            }}
            onMouseEnter={(e) => e.currentTarget.style.backgroundColor = '#2563eb'}
            onMouseLeave={(e) => e.currentTarget.style.backgroundColor = '#3b82f6'}
          >
            Close
          </button>
        </div>
      </div>
    </div>
  );
}
