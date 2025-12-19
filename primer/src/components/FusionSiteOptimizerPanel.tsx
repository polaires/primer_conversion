import React, { useState, useMemo, useCallback } from 'react';
import {
  getDomesticationSummary,
} from '../lib/repp/index.js';

/**
 * Fusion Site Optimizer Panel
 *
 * This component extends the Golden Gate Assembly workflow to support
 * automatic junction site optimization within sequences (SplitSet-like functionality).
 *
 * Features:
 * - Sequence input for fusion site scanning
 * - Candidate site visualization
 * - Multi-factor scoring display
 * - Optimization algorithm selection
 * - Failure mode prediction
 */

// Type Definitions

interface AlgorithmInfo {
  id: string;
  name: string;
  description: string;
  icon: string;
  maxJunctions?: number;
}

interface ScoringWeights {
  overhangQuality: number;
  forwardPrimer: number;
  reversePrimer: number;
  riskFactors: number;
  biologicalContext: number;
}

interface ScoringPreset {
  id: string;
  name: string;
  description: string;
  weights: ScoringWeights;
}

interface ScoreDetail {
  score?: number;
  [key: string]: any;
}

interface CandidateScores {
  overhangQuality?: ScoreDetail;
  forwardPrimer?: ScoreDetail;
  reversePrimer?: ScoreDetail;
  riskFactors?: ScoreDetail;
  biologicalContext?: ScoreDetail;
}

interface FusionCandidate {
  overhang: string;
  position: number;
  composite: number;
  scores?: CandidateScores;
  warnings?: string[];
  isTNNA?: boolean;
  isHighGC?: boolean;
}

interface Constraints {
  minFragmentSize: number;
  maxFragmentSize: number;
  minDistanceFromEnds: number;
  minSetFidelity: number;
  circular: boolean;
}

interface BioContext {
  isCodingSequence: boolean;
  codingFrame: number;
  proteinDomains: any[];
  scarPreferences: string;
}

interface OptimizationSolution {
  junctions?: any[];
  setFidelity?: number;
  overhangs?: string[];
}

interface OptimizationResult {
  solution?: OptimizationSolution;
  algorithm?: string;
  nodesExplored?: number;
  optimal?: boolean;
  failurePrediction?: FailurePrediction;
}

interface FailurePredictionSummary {
  predictedSuccessRate?: number;
  predictedSuccessPercent?: string;
  highRiskCount?: number;
  mediumRiskCount?: number;
  lowRiskCount?: number;
  recommendation?: string;
}

interface Prediction {
  severity: 'high' | 'medium' | 'low';
  type: string;
  fragment?: string;
  probability: number;
  message: string;
  mitigation?: string;
}

interface FailurePrediction {
  summary?: FailurePredictionSummary;
  predictions: Prediction[];
}

interface DomesticationSite {
  sequence: string;
  position: number;
  orientation: string;
  recommendedJunction?: {
    overhang: string;
    quality?: {
      score?: number;
    };
  };
  hasValidOption?: boolean;
}

interface DomesticationSummary {
  status: 'compatible' | 'auto-fixable' | 'needs-attention';
  icon: string;
  title: string;
  description: string;
  sites?: DomesticationSite[];
  additionalFragments?: number;
}

interface AutoDomesticationInfo {
  enabled: boolean;
  sites?: DomesticationSite[];
  additionalFragments?: number;
}

interface OptimizationParams {
  algorithm: string;
  weights: ScoringWeights;
  constraints: Constraints;
  bioContext: BioContext;
  manualCandidates: FusionCandidate[] | null;
  autoDomestication: AutoDomesticationInfo;
  effectiveFragments: number;
}

// Component Props Interfaces

interface FusionCandidateCardProps {
  candidate: FusionCandidate;
  isSelected: boolean;
  onSelect: (candidate: FusionCandidate) => void;
  showDetails: boolean;
}

interface ScoringWeightsEditorProps {
  weights: ScoringWeights;
  onChange: (weights: ScoringWeights) => void;
  preset: string;
  onPresetChange: (preset: string) => void;
}

interface FailurePredictionDisplayProps {
  predictions?: FailurePrediction;
}

interface DomesticationStatusBannerProps {
  summary: DomesticationSummary | null;
  onToggleDetails: () => void;
  showDetails: boolean;
}

interface OptimizerResultSummaryProps {
  result: OptimizationResult | null;
}

interface FusionSiteOptimizerPanelProps {
  sequence: string;
  enzyme?: string;
  numFragments?: number;
  onOptimize?: (params: OptimizationParams) => void;
  result?: OptimizationResult | null;
  candidates?: FusionCandidate[];
  isOptimizing?: boolean;
  globalAutoDomestication?: boolean;
  onAutoDomesticationChange?: (enabled: boolean) => void;
}

// Algorithm selection options
const OPTIMIZER_ALGORITHMS: Record<string, AlgorithmInfo> = {
  AUTO: {
    id: 'auto',
    name: 'Auto Select',
    description: 'Automatically choose best algorithm based on junction count',
    icon: '⚡',
  },
  BRANCH_BOUND: {
    id: 'branch_bound',
    name: 'Branch & Bound',
    description: 'Guarantees global optimum (≤6 junctions)',
    icon: '🎯',
    maxJunctions: 6,
  },
  DP_VALIDATED: {
    id: 'dp_validated',
    name: 'DP + Validation',
    description: 'Dynamic programming with full set validation (6-10 junctions)',
    icon: '📊',
    maxJunctions: 10,
  },
  MONTE_CARLO: {
    id: 'monte_carlo',
    name: 'Monte Carlo',
    description: 'Stochastic optimization for large assemblies (>10 junctions)',
    icon: '🎲',
  },
};

// Scoring weight presets
const SCORING_PRESETS: Record<string, ScoringPreset> = {
  BALANCED: {
    id: 'balanced',
    name: 'Balanced',
    description: 'Equal emphasis on all factors',
    weights: {
      overhangQuality: 0.20,
      forwardPrimer: 0.20,
      reversePrimer: 0.20,
      riskFactors: 0.25,
      biologicalContext: 0.15,
    },
  },
  FIDELITY_FIRST: {
    id: 'fidelity_first',
    name: 'Fidelity First',
    description: 'Maximize assembly fidelity',
    weights: {
      overhangQuality: 0.35,
      forwardPrimer: 0.15,
      reversePrimer: 0.15,
      riskFactors: 0.25,
      biologicalContext: 0.10,
    },
  },
  PRIMER_QUALITY: {
    id: 'primer_quality',
    name: 'Primer Quality',
    description: 'Prioritize PCR success',
    weights: {
      overhangQuality: 0.15,
      forwardPrimer: 0.25,
      reversePrimer: 0.25,
      riskFactors: 0.25,
      biologicalContext: 0.10,
    },
  },
  BIOLOGICAL: {
    id: 'biological',
    name: 'Biology-Aware',
    description: 'Consider coding frames and domains',
    weights: {
      overhangQuality: 0.15,
      forwardPrimer: 0.15,
      reversePrimer: 0.15,
      riskFactors: 0.25,
      biologicalContext: 0.30,
    },
  },
};

// Fusion site candidate card component
function FusionCandidateCard({ candidate, isSelected, onSelect, showDetails }: FusionCandidateCardProps) {
  const getQualityColor = (score: number): string => {
    if (score >= 85) return '#22c55e';
    if (score >= 70) return '#3b82f6';
    if (score >= 55) return '#f59e0b';
    return '#ef4444';
  };

  const getQualityLabel = (score: number): string => {
    if (score >= 85) return 'Excellent';
    if (score >= 70) return 'Good';
    if (score >= 55) return 'Acceptable';
    return 'Poor';
  };

  return (
    <div
      className={`fusion-candidate-card ${isSelected ? 'selected' : ''}`}
      onClick={() => onSelect(candidate)}
    >
      <div className="candidate-header">
        <code className="candidate-overhang">{candidate.overhang}</code>
        <span className="candidate-position">pos {candidate.position}</span>
        {isSelected && <span className="selected-badge">✓</span>}
      </div>

      <div className="candidate-scores">
        <div className="score-bar-group">
          <div className="score-bar-label">Composite</div>
          <div className="score-bar-track">
            <div
              className="score-bar-fill"
              style={{
                width: `${candidate.composite || 0}%`,
                backgroundColor: getQualityColor(candidate.composite || 0),
              }}
            />
          </div>
          <span
            className="score-value"
            style={{ color: getQualityColor(candidate.composite || 0) }}
          >
            {(candidate.composite || 0).toFixed(0)}
          </span>
        </div>
      </div>

      {showDetails && candidate.scores && (
        <div className="candidate-details">
          <div className="detail-row">
            <span className="detail-label">Overhang Quality</span>
            <span className="detail-value">{candidate.scores.overhangQuality?.score || '—'}</span>
          </div>
          <div className="detail-row">
            <span className="detail-label">Forward Primer</span>
            <span className="detail-value">{candidate.scores.forwardPrimer?.score || '—'}</span>
          </div>
          <div className="detail-row">
            <span className="detail-label">Reverse Primer</span>
            <span className="detail-value">{candidate.scores.reversePrimer?.score || '—'}</span>
          </div>
          <div className="detail-row">
            <span className="detail-label">Risk Score</span>
            <span className="detail-value">{candidate.scores.riskFactors?.score || '—'}</span>
          </div>
          {candidate.warnings?.length && candidate.warnings.length > 0 && (
            <div className="candidate-warnings">
              {candidate.warnings.map((w, i) => (
                <span key={i} className="warning-tag">{w}</span>
              ))}
            </div>
          )}
        </div>
      )}

      <div className="candidate-meta">
        <span className={`quality-badge quality-${getQualityLabel(candidate.composite || 0).toLowerCase()}`}>
          {getQualityLabel(candidate.composite || 0)}
        </span>
        {candidate.isTNNA && <span className="efficiency-badge">TNNA</span>}
        {candidate.isHighGC && <span className="efficiency-badge gc">High GC</span>}
      </div>
    </div>
  );
}

// Scoring weights editor component
function ScoringWeightsEditor({ weights, onChange, preset, onPresetChange }: ScoringWeightsEditorProps) {
  const [isCustom, setIsCustom] = useState<boolean>(false);

  const handleWeightChange = (key: keyof ScoringWeights, value: string) => {
    const newWeights = { ...weights, [key]: parseFloat(value) };
    // Normalize weights to sum to 1
    const total = Object.values(newWeights).reduce((sum, w) => sum + w, 0);
    const normalized: ScoringWeights = {} as ScoringWeights;
    for (const k in newWeights) {
      normalized[k as keyof ScoringWeights] = newWeights[k as keyof ScoringWeights] / total;
    }
    onChange(normalized);
    setIsCustom(true);
  };

  return (
    <div className="scoring-weights-editor">
      <div className="preset-selector">
        <label>Scoring Preset</label>
        <div className="preset-chips">
          {Object.entries(SCORING_PRESETS).map(([, p]) => (
            <button
              key={p.id}
              className={`preset-chip ${preset === p.id ? 'active' : ''} ${isCustom && preset !== p.id ? 'inactive' : ''}`}
              onClick={() => {
                onPresetChange(p.id);
                onChange(p.weights);
                setIsCustom(false);
              }}
              title={p.description}
            >
              {p.name}
            </button>
          ))}
          {isCustom && (
            <span className="custom-indicator">Custom</span>
          )}
        </div>
      </div>

      <div className="weights-sliders">
        {Object.entries(weights).map(([key, value]) => (
          <div key={key} className="weight-slider-row">
            <label className="weight-label">
              {key.replace(/([A-Z])/g, ' $1').replace(/^./, s => s.toUpperCase())}
            </label>
            <input
              type="range"
              min="0"
              max="0.5"
              step="0.05"
              value={value}
              onChange={(e: React.ChangeEvent<HTMLInputElement>) => handleWeightChange(key as keyof ScoringWeights, e.target.value)}
              className="weight-slider"
            />
            <span className="weight-value">{(value * 100).toFixed(0)}%</span>
          </div>
        ))}
      </div>
    </div>
  );
}

// Failure prediction display component
function FailurePredictionDisplay({ predictions }: FailurePredictionDisplayProps) {
  if (!predictions || predictions.predictions?.length === 0) {
    return (
      <div className="failure-prediction-empty">
        <span className="success-icon">✓</span>
        <span>No significant failure risks predicted</span>
      </div>
    );
  }

  const getSeverityColor = (severity: string): string => {
    switch (severity) {
      case 'high': return '#ef4444';
      case 'medium': return '#f59e0b';
      case 'low': return '#3b82f6';
      default: return '#6b7280';
    }
  };

  const getSeverityIcon = (severity: string): string => {
    switch (severity) {
      case 'high': return '🔴';
      case 'medium': return '🟠';
      case 'low': return '🔵';
      default: return '⚪';
    }
  };

  return (
    <div className="failure-prediction-display">
      {/* Summary Header */}
      <div className="prediction-summary">
        <div className="summary-gauge">
          <svg viewBox="0 0 100 100" className="success-gauge-svg">
            <circle cx="50" cy="50" r="40" fill="none" stroke="#e5e7eb" strokeWidth="10" />
            <circle
              cx="50"
              cy="50"
              r="40"
              fill="none"
              stroke={predictions.summary?.predictedSuccessRate && predictions.summary.predictedSuccessRate >= 0.8 ? '#22c55e' :
                      predictions.summary?.predictedSuccessRate && predictions.summary.predictedSuccessRate >= 0.6 ? '#f59e0b' : '#ef4444'}
              strokeWidth="10"
              strokeLinecap="round"
              strokeDasharray={`${(predictions.summary?.predictedSuccessRate || 0) * 251.2} 251.2`}
              transform="rotate(-90 50 50)"
            />
          </svg>
          <div className="gauge-center-text">
            <span className="gauge-percent">{predictions.summary?.predictedSuccessPercent || '—'}</span>
            <span className="gauge-label">Success</span>
          </div>
        </div>

        <div className="summary-stats">
          <div className="stat-item high">
            <span className="stat-count">{predictions.summary?.highRiskCount || 0}</span>
            <span className="stat-label">High Risk</span>
          </div>
          <div className="stat-item medium">
            <span className="stat-count">{predictions.summary?.mediumRiskCount || 0}</span>
            <span className="stat-label">Medium Risk</span>
          </div>
          <div className="stat-item low">
            <span className="stat-count">{predictions.summary?.lowRiskCount || 0}</span>
            <span className="stat-label">Low Risk</span>
          </div>
        </div>
      </div>

      {/* Recommendation Banner */}
      {predictions.summary?.recommendation && (
        <div className={`recommendation-banner ${
          predictions.summary.predictedSuccessRate && predictions.summary.predictedSuccessRate >= 0.8 ? 'success' :
          predictions.summary.predictedSuccessRate && predictions.summary.predictedSuccessRate >= 0.6 ? 'warning' : 'danger'
        }`}>
          <span className="recommendation-icon">
            {predictions.summary.predictedSuccessRate && predictions.summary.predictedSuccessRate >= 0.8 ? '✓' :
             predictions.summary.predictedSuccessRate && predictions.summary.predictedSuccessRate >= 0.6 ? '⚠' : '✕'}
          </span>
          <span className="recommendation-text">{predictions.summary.recommendation}</span>
        </div>
      )}

      {/* Prediction List */}
      <div className="predictions-list">
        {predictions.predictions.map((pred, i) => (
          <div key={i} className={`prediction-item severity-${pred.severity}`}>
            <div className="prediction-header">
              <span className="prediction-icon">{getSeverityIcon(pred.severity)}</span>
              <span className="prediction-type">{pred.type.replace(/_/g, ' ')}</span>
              {pred.fragment && <span className="prediction-fragment">{pred.fragment}</span>}
              <span
                className="prediction-probability"
                style={{ color: getSeverityColor(pred.severity) }}
              >
                {(pred.probability * 100).toFixed(0)}% risk
              </span>
            </div>
            <p className="prediction-message">{pred.message}</p>
            {pred.mitigation && (
              <div className="prediction-mitigation">
                <span className="mitigation-label">Mitigation:</span>
                <span className="mitigation-text">{pred.mitigation}</span>
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
}

// Auto-Domestication Status Banner
function DomesticationStatusBanner({ summary, onToggleDetails, showDetails }: DomesticationStatusBannerProps) {
  if (!summary) return null;

  const getStatusStyles = (): { bg: string; border: string; icon: string } => {
    switch (summary.status) {
      case 'compatible':
        return {
          bg: 'linear-gradient(135deg, #dcfce7 0%, #bbf7d0 100%)',
          border: '#22c55e',
          icon: '#16a34a',
        };
      case 'auto-fixable':
        return {
          bg: 'linear-gradient(135deg, #fef3c7 0%, #fde68a 100%)',
          border: '#f59e0b',
          icon: '#d97706',
        };
      case 'needs-attention':
        return {
          bg: 'linear-gradient(135deg, #fee2e2 0%, #fecaca 100%)',
          border: '#ef4444',
          icon: '#dc2626',
        };
      default:
        return {
          bg: '#f3f4f6',
          border: '#9ca3af',
          icon: '#6b7280',
        };
    }
  };

  const styles = getStatusStyles();

  return (
    <div
      className="domestication-status-banner"
      style={{
        background: styles.bg,
        borderLeft: `4px solid ${styles.border}`,
        borderRadius: '8px',
        padding: '12px 16px',
        marginBottom: '16px',
      }}
    >
      <div className="flex items-start gap-3">
        <span className="text-2xl leading-none">{summary.icon}</span>
        <div className="flex-1">
          <div className="font-semibold text-gray-900 mb-1">
            {summary.title}
          </div>
          <div className="text-sm text-gray-700">
            {summary.description}
          </div>

          {summary.sites && summary.sites.length > 0 && (
            <button
              onClick={onToggleDetails}
              style={{
                background: 'none',
                border: 'none',
                color: styles.border,
                fontSize: '12px',
                cursor: 'pointer',
                padding: '4px 0',
                marginTop: '8px',
                display: 'flex',
                alignItems: 'center',
                gap: '4px',
              }}
            >
              {showDetails ? '▼ Hide Details' : '▶ Show Sites'}
            </button>
          )}

          {showDetails && summary.sites && summary.sites.length > 0 && (
            <div className="mt-3">
              {summary.sites.map((site, i) => (
                <div
                  key={i}
                  className="rounded-md p-3 mb-2 text-sm"
                  style={{ background: 'rgba(255,255,255,0.7)' }}
                >
                  <div className="flex items-center gap-3 flex-wrap">
                    <span className="font-medium">Site {i + 1}</span>
                    <code style={{
                      background: '#fef3c7',
                      padding: '2px 6px',
                      borderRadius: '4px',
                      fontFamily: 'monospace',
                      fontSize: '12px',
                    }}>
                      {site.sequence}
                    </code>
                    <span className="text-gray-600 text-xs">
                      pos {site.position} ({site.orientation})
                    </span>
                    {site.recommendedJunction && (
                      <span style={{
                        background: site.hasValidOption ? '#dcfce7' : '#fee2e2',
                        color: site.hasValidOption ? '#16a34a' : '#dc2626',
                        padding: '2px 8px',
                        borderRadius: '4px',
                        fontSize: '11px',
                        fontWeight: 500,
                      }}>
                        Junction: <code>{site.recommendedJunction.overhang}</code>
                        {' '}(score: {site.recommendedJunction.quality?.score || '—'})
                      </span>
                    )}
                  </div>
                </div>
              ))}

              {summary.additionalFragments && summary.additionalFragments > 0 && (
                <div className="mt-3 p-2 rounded-md text-xs text-yellow-900" style={{ background: 'rgba(245, 158, 11, 0.1)' }}>
                  <strong>Auto-Domestication:</strong> Each internal site will become an additional
                  fragment boundary. Your {summary.additionalFragments === 1 ? 'assembly' : 'assemblies'} will
                  have <strong>+{summary.additionalFragments}</strong> extra fragment{summary.additionalFragments > 1 ? 's' : ''}.
                </div>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

// Optimizer result summary component
function OptimizerResultSummary({ result }: OptimizerResultSummaryProps) {
  if (!result) return null;

  return (
    <div className="optimizer-result-summary">
      <div className="result-header">
        <div className="result-badge success">
          <svg viewBox="0 0 24 24" width="24" height="24" fill="currentColor">
            <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
          </svg>
        </div>
        <div className="result-title">
          <h3>Optimization Complete</h3>
          <p>
            {result.solution?.junctions?.length || 0} junctions •
            {result.solution?.setFidelity ? ` ${(result.solution.setFidelity * 100).toFixed(1)}% fidelity` : ''} •
            Algorithm: {result.algorithm || 'auto'}
          </p>
        </div>
      </div>

      {/* Algorithm Info */}
      <div className="algorithm-info-panel">
        <div className="info-item">
          <span className="info-label">Algorithm Used</span>
          <span className="info-value">
            {OPTIMIZER_ALGORITHMS[result.algorithm?.toUpperCase() || '']?.name || result.algorithm || 'Auto'}
          </span>
        </div>
        {result.nodesExplored !== undefined && (
          <div className="info-item">
            <span className="info-label">Nodes Explored</span>
            <span className="info-value">{result.nodesExplored?.toLocaleString()}</span>
          </div>
        )}
        {result.optimal !== undefined && (
          <div className="info-item">
            <span className="info-label">Optimality</span>
            <span className={`info-value ${result.optimal ? 'optimal' : 'heuristic'}`}>
              {result.optimal ? 'Global Optimum' : 'Heuristic Solution'}
            </span>
          </div>
        )}
      </div>

      {/* Junction Summary */}
      {result.solution?.overhangs && (
        <div className="junctions-summary">
          <h4>Selected Overhangs</h4>
          <div className="overhangs-list">
            {result.solution.overhangs.map((oh, i) => (
              <div key={i} className="overhang-chip">
                <span className="junction-num">{i + 1}</span>
                <code className="overhang-seq">{oh}</code>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Set Fidelity Gauge */}
      {result.solution?.setFidelity !== undefined && (
        <div className="set-fidelity-display">
          <div className="fidelity-label">Assembly Fidelity</div>
          <div className="fidelity-bar-large">
            <div
              className="fidelity-fill"
              style={{
                width: `${result.solution.setFidelity * 100}%`,
                backgroundColor: result.solution.setFidelity >= 0.95 ? '#22c55e' :
                                 result.solution.setFidelity >= 0.90 ? '#84cc16' :
                                 result.solution.setFidelity >= 0.80 ? '#f59e0b' : '#ef4444',
              }}
            />
          </div>
          <span className="fidelity-percent">
            {(result.solution.setFidelity * 100).toFixed(1)}%
          </span>
        </div>
      )}
    </div>
  );
}

// Main Fusion Site Optimizer Panel
export default function FusionSiteOptimizerPanel({
  sequence,
  enzyme = 'BsaI',
  numFragments = 4,
  onOptimize,
  result,
  candidates,
  isOptimizing,
  globalAutoDomestication = true,
  onAutoDomesticationChange,
}: FusionSiteOptimizerPanelProps) {
  // State for optimizer settings
  const [algorithm, setAlgorithm] = useState<string>('auto');
  const [scoringPreset, setScoringPreset] = useState<string>('balanced');
  const [weights, setWeights] = useState<ScoringWeights>(SCORING_PRESETS.BALANCED.weights);
  const [selectedCandidates, setSelectedCandidates] = useState<FusionCandidate[]>([]);
  const [showAllCandidates, setShowAllCandidates] = useState<boolean>(false);
  const [activeTab, setActiveTab] = useState<string>('candidates');

  // Auto-domestication state - use global prop if provided
  const [showDomesticationDetails, setShowDomesticationDetails] = useState<boolean>(false);
  const autoDomesticationEnabled = globalAutoDomestication;
  const setAutoDomesticationEnabled = onAutoDomesticationChange || (() => {});

  // Analyze sequence for domestication needs
  const domesticationSummary = useMemo<DomesticationSummary | null>(() => {
    if (!sequence || sequence.length < 50) return null;
    try {
      return getDomesticationSummary(sequence, enzyme) as DomesticationSummary;
    } catch (e) {
      console.warn('Domestication analysis error:', e);
      return null;
    }
  }, [sequence, enzyme]);

  // Calculate effective fragment count (user request + domestication junctions)
  const effectiveFragments = useMemo<number>(() => {
    if (!domesticationSummary || !autoDomesticationEnabled) return numFragments;
    return numFragments + (domesticationSummary.additionalFragments || 0);
  }, [numFragments, domesticationSummary, autoDomesticationEnabled]);

  // Constraint settings
  const [constraints, setConstraints] = useState<Constraints>({
    minFragmentSize: 200,
    maxFragmentSize: 3000,
    minDistanceFromEnds: 50,
    minSetFidelity: 0.90,
    circular: true,
  });

  // Biological context settings
  const [bioContext, setBioContext] = useState<BioContext>({
    isCodingSequence: false,
    codingFrame: 0,
    proteinDomains: [],
    scarPreferences: 'nonCoding',
  });

  // Handle candidate selection
  const handleCandidateSelect = useCallback((candidate: FusionCandidate) => {
    setSelectedCandidates(prev => {
      const exists = prev.find(c => c.position === candidate.position);
      if (exists) {
        return prev.filter(c => c.position !== candidate.position);
      }
      return [...prev, candidate].sort((a, b) => a.position - b.position);
    });
  }, []);

  // Run optimization
  const handleOptimize = useCallback(() => {
    if (onOptimize) {
      // Include auto-domestication info in optimization request
      const domesticationInfo: AutoDomesticationInfo = autoDomesticationEnabled && domesticationSummary?.status !== 'compatible'
        ? {
            enabled: true,
            sites: domesticationSummary?.sites || [],
            additionalFragments: domesticationSummary?.additionalFragments || 0,
          }
        : { enabled: false };

      onOptimize({
        algorithm,
        weights,
        constraints,
        bioContext,
        manualCandidates: selectedCandidates.length > 0 ? selectedCandidates : null,
        autoDomestication: domesticationInfo,
        effectiveFragments,
      });
    }
  }, [algorithm, weights, constraints, bioContext, selectedCandidates, onOptimize, autoDomesticationEnabled, domesticationSummary, effectiveFragments]);

  // Filter and sort candidates for display
  const displayCandidates = useMemo<FusionCandidate[]>(() => {
    if (!candidates || candidates.length === 0) return [];
    const sorted = [...candidates].sort((a, b) => (b.composite || 0) - (a.composite || 0));
    return showAllCandidates ? sorted : sorted.slice(0, 20);
  }, [candidates, showAllCandidates]);

  return (
    <div className="fusion-optimizer-panel">
      {/* Header */}
      <div className="optimizer-header">
        <div className="header-content">
          <h3>
            <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
              <path d="M19.43 12.98c.04-.32.07-.64.07-.98s-.03-.66-.07-.98l2.11-1.65c.19-.15.24-.42.12-.64l-2-3.46c-.12-.22-.39-.3-.61-.22l-2.49 1c-.52-.4-1.08-.73-1.69-.98l-.38-2.65C14.46 2.18 14.25 2 14 2h-4c-.25 0-.46.18-.49.42l-.38 2.65c-.61.25-1.17.59-1.69.98l-2.49-1c-.23-.09-.49 0-.61.22l-2 3.46c-.13.22-.07.49.12.64l2.11 1.65c-.04.32-.07.65-.07.98s.03.66.07.98l-2.11 1.65c-.19.15-.24.42-.12.64l2 3.46c.12.22.39.3.61.22l2.49-1c.52.4 1.08.73 1.69.98l.38 2.65c.03.24.24.42.49.42h4c.25 0 .46-.18.49-.42l.38-2.65c.61-.25 1.17-.59 1.69-.98l2.49 1c.23.09.49 0 .61-.22l2-3.46c.12-.22.07-.49-.12-.64l-2.11-1.65zM12 15.5c-1.93 0-3.5-1.57-3.5-3.5s1.57-3.5 3.5-3.5 3.5 1.57 3.5 3.5-1.57 3.5-3.5 3.5z"/>
            </svg>
            Fusion Site Optimizer
          </h3>
          <p>Automatically find optimal junction positions within your sequence</p>
        </div>
        <div className="header-stats">
          {sequence && <span className="stat">{sequence.length.toLocaleString()} bp</span>}
          <span className="stat">{enzyme}</span>
          <span className="stat">
            {effectiveFragments !== numFragments ? (
              <span title={`${numFragments} requested + ${effectiveFragments - numFragments} for domestication`}>
                {effectiveFragments} fragments
                <span className="text-[10px] text-amber-500 ml-1">
                  (+{effectiveFragments - numFragments})
                </span>
              </span>
            ) : (
              `${numFragments} fragments`
            )}
          </span>
        </div>
      </div>

      {/* Auto-Domestication Status Banner */}
      {domesticationSummary && (
        <div className="px-4">
          <DomesticationStatusBanner
            summary={domesticationSummary}
            showDetails={showDomesticationDetails}
            onToggleDetails={() => setShowDomesticationDetails(!showDomesticationDetails)}
          />
          {domesticationSummary.status !== 'compatible' && (
            <div className="flex items-center gap-2 mb-3 text-sm">
              <label className="flex items-center gap-1.5 cursor-pointer">
                <input
                  type="checkbox"
                  checked={autoDomesticationEnabled}
                  onChange={(e: React.ChangeEvent<HTMLInputElement>) => setAutoDomesticationEnabled(e.target.checked)}
                  style={{ accentColor: '#f59e0b' }}
                />
                <span className="text-gray-700">
                  Enable auto-domestication (add junctions to break internal sites)
                </span>
              </label>
            </div>
          )}
        </div>
      )}

      {/* Main Content */}
      <div className="optimizer-content">
        {/* Settings Panel */}
        <div className="optimizer-settings">
          {/* Algorithm Selection */}
          <div className="settings-section">
            <h4>
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                <path d="M9 4v1.38c-.83-.33-1.72-.5-2.61-.5-1.79 0-3.58.68-4.95 2.05l3.33 3.33h1.11v1.11c.86.86 1.98 1.31 3.11 1.36V15H6v3c0 1.1.9 2 2 2h10c1.1 0 2-.9 2-2V4H9zm-1.11 6.41V8.26H5.61L4.57 7.22a5.07 5.07 0 0 1 1.82-.34c1.34 0 2.59.52 3.54 1.46l1.41 1.41-.2.2c-.51.51-1.19.8-1.92.8-.47 0-.93-.12-1.33-.34zM18 18H8v-1h2v-2h4v2h2v1zm0-4h-4v-1h-2v-2h-1v-1H9v-2h3V5h6v9z"/>
              </svg>
              Optimization Algorithm
            </h4>
            <div className="algorithm-selector">
              {Object.entries(OPTIMIZER_ALGORITHMS).map(([, alg]) => (
                <button
                  key={alg.id}
                  className={`algorithm-btn ${algorithm === alg.id ? 'active' : ''}`}
                  onClick={() => setAlgorithm(alg.id)}
                  title={alg.description}
                  disabled={!!(alg.maxJunctions && numFragments - 1 > alg.maxJunctions)}
                >
                  <span className="alg-icon">{alg.icon}</span>
                  <span className="alg-name">{alg.name}</span>
                  {alg.maxJunctions && (
                    <span className="alg-limit">≤{alg.maxJunctions} junctions</span>
                  )}
                </button>
              ))}
            </div>
          </div>

          {/* Scoring Weights */}
          <div className="settings-section">
            <h4>
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                <path d="M12 3c-4.97 0-9 4.03-9 9s4.03 9 9 9c.83 0 1.5-.67 1.5-1.5 0-.39-.15-.74-.39-1.01-.23-.26-.38-.61-.38-.99 0-.83.67-1.5 1.5-1.5H16c2.76 0 5-2.24 5-5 0-4.42-4.03-8-9-8zm-5.5 9c-.83 0-1.5-.67-1.5-1.5S5.67 9 6.5 9 8 9.67 8 10.5 7.33 12 6.5 12zm3-4C8.67 8 8 7.33 8 6.5S8.67 5 9.5 5s1.5.67 1.5 1.5S10.33 8 9.5 8zm5 0c-.83 0-1.5-.67-1.5-1.5S13.67 5 14.5 5s1.5.67 1.5 1.5S15.33 8 14.5 8zm3 4c-.83 0-1.5-.67-1.5-1.5S16.67 9 17.5 9s1.5.67 1.5 1.5-.67 1.5-1.5 1.5z"/>
              </svg>
              Scoring Weights
            </h4>
            <ScoringWeightsEditor
              weights={weights}
              onChange={setWeights}
              preset={scoringPreset}
              onPresetChange={setScoringPreset}
            />
          </div>

          {/* Constraints */}
          <div className="settings-section">
            <h4>
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                <path d="M12 1L3 5v6c0 5.55 3.84 10.74 9 12 5.16-1.26 9-6.45 9-12V5l-9-4zm0 10.99h7c-.53 4.12-3.28 7.79-7 8.94V12H5V6.3l7-3.11v8.8z"/>
              </svg>
              Constraints
            </h4>
            <div className="constraints-grid">
              <div className="constraint-row">
                <label>Min Fragment Size</label>
                <input
                  type="number"
                  value={constraints.minFragmentSize}
                  onChange={(e: React.ChangeEvent<HTMLInputElement>) => setConstraints(c => ({ ...c, minFragmentSize: parseInt(e.target.value) || 200 }))}
                  min={100}
                  max={1000}
                />
                <span className="constraint-unit">bp</span>
              </div>
              <div className="constraint-row">
                <label>Max Fragment Size</label>
                <input
                  type="number"
                  value={constraints.maxFragmentSize}
                  onChange={(e: React.ChangeEvent<HTMLInputElement>) => setConstraints(c => ({ ...c, maxFragmentSize: parseInt(e.target.value) || 3000 }))}
                  min={500}
                  max={10000}
                />
                <span className="constraint-unit">bp</span>
              </div>
              <div className="constraint-row">
                <label>Min Set Fidelity</label>
                <input
                  type="range"
                  value={constraints.minSetFidelity}
                  onChange={(e: React.ChangeEvent<HTMLInputElement>) => setConstraints(c => ({ ...c, minSetFidelity: parseFloat(e.target.value) }))}
                  min={0.7}
                  max={0.99}
                  step={0.01}
                />
                <span className="constraint-value">{(constraints.minSetFidelity * 100).toFixed(0)}%</span>
              </div>
              <div className="constraint-row checkbox">
                <label>
                  <input
                    type="checkbox"
                    checked={constraints.circular}
                    onChange={(e: React.ChangeEvent<HTMLInputElement>) => setConstraints(c => ({ ...c, circular: e.target.checked }))}
                  />
                  Circular Assembly
                </label>
              </div>
            </div>
          </div>

          {/* Biological Context */}
          <div className="settings-section">
            <h4>
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                <path d="M18 4l-4 4h3v7c0 1.1-.9 2-2 2s-2-.9-2-2V8c0-2.21-1.79-4-4-4S5 5.79 5 8v7H2l4 4 4-4H7V8c0-1.1.9-2 2-2s2 .9 2 2v7c0 2.21 1.79 4 4 4s4-1.79 4-4V8h3l-4-4z"/>
              </svg>
              Biological Context
            </h4>
            <div className="bio-context-settings">
              <div className="constraint-row checkbox">
                <label>
                  <input
                    type="checkbox"
                    checked={bioContext.isCodingSequence}
                    onChange={(e: React.ChangeEvent<HTMLInputElement>) => setBioContext(c => ({ ...c, isCodingSequence: e.target.checked }))}
                  />
                  Coding Sequence
                </label>
              </div>
              {bioContext.isCodingSequence && (
                <>
                  <div className="constraint-row">
                    <label>Reading Frame</label>
                    <select
                      value={bioContext.codingFrame}
                      onChange={(e: React.ChangeEvent<HTMLSelectElement>) => setBioContext(c => ({ ...c, codingFrame: parseInt(e.target.value) }))}
                    >
                      <option value={0}>Frame 1 (0)</option>
                      <option value={1}>Frame 2 (+1)</option>
                      <option value={2}>Frame 3 (+2)</option>
                    </select>
                  </div>
                  <div className="constraint-row">
                    <label>Scar Preferences</label>
                    <select
                      value={bioContext.scarPreferences}
                      onChange={(e: React.ChangeEvent<HTMLSelectElement>) => setBioContext(c => ({ ...c, scarPreferences: e.target.value }))}
                    >
                      <option value="coding">Coding (avoid stops)</option>
                      <option value="linker">Linker regions</option>
                      <option value="nonCoding">Non-coding</option>
                    </select>
                  </div>
                </>
              )}
            </div>
          </div>

          {/* Optimize Button */}
          <button
            className="optimize-btn"
            onClick={handleOptimize}
            disabled={!sequence || isOptimizing}
          >
            {isOptimizing ? (
              <>
                <span className="spinner"></span>
                Optimizing...
              </>
            ) : (
              <>
                <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
                  <path d="M17.65 6.35C16.2 4.9 14.21 4 12 4c-4.42 0-7.99 3.58-7.99 8s3.57 8 7.99 8c3.73 0 6.84-2.55 7.73-6h-2.08c-.82 2.33-3.04 4-5.65 4-3.31 0-6-2.69-6-6s2.69-6 6-6c1.66 0 3.14.69 4.22 1.78L13 11h7V4l-2.35 2.35z"/>
                </svg>
                Find Optimal Junctions
              </>
            )}
          </button>
        </div>

        {/* Results Area */}
        <div className="optimizer-results">
          {/* Tabs */}
          <div className="results-tabs">
            <button
              className={`tab-btn ${activeTab === 'candidates' ? 'active' : ''}`}
              onClick={() => setActiveTab('candidates')}
            >
              Candidates
              {candidates && <span className="tab-count">{candidates.length}</span>}
            </button>
            <button
              className={`tab-btn ${activeTab === 'result' ? 'active' : ''}`}
              onClick={() => setActiveTab('result')}
              disabled={!result}
            >
              Optimization Result
            </button>
            <button
              className={`tab-btn ${activeTab === 'failures' ? 'active' : ''}`}
              onClick={() => setActiveTab('failures')}
              disabled={!result?.failurePrediction}
            >
              Failure Analysis
              {result?.failurePrediction?.summary?.highRiskCount && result.failurePrediction.summary.highRiskCount > 0 && (
                <span className="tab-badge danger">{result.failurePrediction.summary.highRiskCount}</span>
              )}
            </button>
          </div>

          {/* Tab Content */}
          <div className="tab-content">
            {activeTab === 'candidates' && (
              <div className="candidates-panel">
                {!candidates || candidates.length === 0 ? (
                  <div className="empty-state">
                    <svg viewBox="0 0 24 24" width="48" height="48" fill="currentColor" opacity="0.3">
                      <path d="M15.5 14h-.79l-.28-.27C15.41 12.59 16 11.11 16 9.5 16 5.91 13.09 3 9.5 3S3 5.91 3 9.5 5.91 16 9.5 16c1.61 0 3.09-.59 4.23-1.57l.27.28v.79l5 4.99L20.49 19l-4.99-5zm-6 0C7.01 14 5 11.99 5 9.5S7.01 5 9.5 5 14 7.01 14 9.5 11.99 14 9.5 14z"/>
                    </svg>
                    <p>No fusion site candidates scanned yet</p>
                    <small>Enter a sequence and run optimization to scan for candidates</small>
                  </div>
                ) : (
                  <>
                    <div className="candidates-header">
                      <span className="candidates-count">{candidates.length} candidates found</span>
                      <button
                        className="show-all-btn"
                        onClick={() => setShowAllCandidates(!showAllCandidates)}
                      >
                        {showAllCandidates ? 'Show Top 20' : 'Show All'}
                      </button>
                      {selectedCandidates.length > 0 && (
                        <span className="selected-count">{selectedCandidates.length} selected</span>
                      )}
                    </div>
                    <div className="candidates-grid">
                      {displayCandidates.map((candidate) => (
                        <FusionCandidateCard
                          key={candidate.position}
                          candidate={candidate}
                          isSelected={selectedCandidates.some(c => c.position === candidate.position)}
                          onSelect={handleCandidateSelect}
                          showDetails={selectedCandidates.some(c => c.position === candidate.position)}
                        />
                      ))}
                    </div>
                  </>
                )}
              </div>
            )}

            {activeTab === 'result' && (
              <div className="result-panel">
                {result ? (
                  <OptimizerResultSummary result={result} />
                ) : (
                  <div className="empty-state">
                    <svg viewBox="0 0 24 24" width="48" height="48" fill="currentColor" opacity="0.3">
                      <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm-5 14H7v-2h7v2zm3-4H7v-2h10v2zm0-4H7V7h10v2z"/>
                    </svg>
                    <p>No optimization result yet</p>
                    <small>Run the optimizer to find the best junction positions</small>
                  </div>
                )}
              </div>
            )}

            {activeTab === 'failures' && (
              <div className="failures-panel">
                {result?.failurePrediction ? (
                  <FailurePredictionDisplay predictions={result.failurePrediction} />
                ) : (
                  <div className="empty-state">
                    <svg viewBox="0 0 24 24" width="48" height="48" fill="currentColor" opacity="0.3">
                      <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
                    </svg>
                    <p>No failure analysis available</p>
                    <small>Run optimization to see predicted failure modes</small>
                  </div>
                )}
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}
