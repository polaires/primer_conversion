import { useState } from 'react';
import SequenceViewer from './SequenceViewer';
import ScoreBreakdownPopup from './primers/ScoreBreakdownPopup.jsx';

// Types for piecewise scoring
interface PiecewiseScores {
  tm: number;
  gc: number;
  terminal3DG: number;
  offTarget: number;
  hairpin: number;
  homodimer: number;
  heterodimer?: number;
  tmDiff?: number;
  [key: string]: number | undefined;
}

// Types for scoring object
interface PrimerScoring {
  compositeScore?: number;
  qualityTier?: QualityTier;
  penalty?: number;
  piecewiseScores?: PiecewiseScores;
}

// Primer result from library
interface PrimerResult {
  seq: string;
  tm: number;
  tmTotal: number;
  gc: number;
  dg: number;
  len: number;
  offTargetCount: number;
  scoring: PrimerScoring;
  dict: () => Record<string, unknown>;
}

type QualityTier = 'excellent' | 'good' | 'acceptable' | 'marginal' | 'poor';

interface QualityInfoItem {
  color: string;
  label: string;
  description: string;
}

// Quality tier colors and descriptions
const QUALITY_INFO: Record<QualityTier, QualityInfoItem> = {
  excellent: {
    color: '#22c55e',
    label: 'Excellent',
    description: 'High confidence primer pair - expected >90% PCR success rate',
  },
  good: {
    color: '#3b82f6',
    label: 'Good',
    description: 'Reliable primer pair - expected ~80% PCR success rate',
  },
  acceptable: {
    color: '#eab308',
    label: 'Acceptable',
    description: 'May work but consider alternatives - expected ~65% success rate',
  },
  marginal: {
    color: '#f97316',
    label: 'Marginal',
    description: 'Likely to have issues - expected ~50% success rate',
  },
  poor: {
    color: '#ef4444',
    label: 'Poor',
    description: 'High failure risk - redesign recommended',
  },
};

type ViewerMode = 'linear' | 'circular' | 'both';

interface PrimerResultsProps {
  results: PrimerResult[] | null;
  error: string | null;
  sequence: string;
  options?: {
    addFwd?: string;
    addRev?: string;
    [key: string]: unknown;
  };
}

export default function PrimerResults({ results, error, sequence, options = {} }: PrimerResultsProps) {
  const [showJson, setShowJson] = useState<boolean>(false);
  const [viewerMode, setViewerMode] = useState<ViewerMode>('linear');
  const [showScoreBreakdown, setShowScoreBreakdown] = useState<boolean>(false);

  if (error) {
    return (
      <div className="results-container error">
        <h3>Error</h3>
        <p>{error}</p>
      </div>
    );
  }

  if (!results || results.length === 0) {
    return null;
  }

  const fwd = results[0];
  const rev = results[1];

  return (
    <div className="results-container">
      <div className="results-header">
        <h3>Designed Primers</h3>
        <button
          type="button"
          className="toggle-json"
          onClick={() => setShowJson(!showJson)}
        >
          {showJson ? 'Table View' : 'JSON View'}
        </button>
      </div>

      {showJson ? (
        <pre className="json-output">
          {JSON.stringify(
            results.filter(Boolean).map((p) => p.dict()),
            null,
            2
          )}
        </pre>
      ) : (
        <>
          {/* Quality Score Badge - Clickable for breakdown */}
          {fwd?.scoring?.compositeScore !== undefined && (
            <>
              <button
                type="button"
                className="quality-badge"
                onClick={() => setShowScoreBreakdown(true)}
                style={{
                  backgroundColor: fwd.scoring.qualityTier ? QUALITY_INFO[fwd.scoring.qualityTier]?.color : '#6b7280',
                  color: 'white',
                  padding: '8px 16px',
                  borderRadius: '8px',
                  marginBottom: '16px',
                  display: 'inline-flex',
                  alignItems: 'center',
                  gap: '8px',
                  cursor: 'pointer',
                  border: 'none',
                  fontSize: '14px',
                  transition: 'all 0.2s ease',
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.transform = 'scale(1.02)';
                  e.currentTarget.style.boxShadow = '0 4px 12px rgba(0,0,0,0.2)';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.transform = 'scale(1)';
                  e.currentTarget.style.boxShadow = 'none';
                }}
                title="Click to see detailed score breakdown"
              >
                <strong>Quality Score: {fwd.scoring.compositeScore}/100</strong>
                {' '}
                ({fwd.scoring.qualityTier ? QUALITY_INFO[fwd.scoring.qualityTier]?.label : fwd.scoring.qualityTier})
                <span style={{ fontSize: '12px', opacity: 0.8 }}>• Click for details</span>
              </button>

              {/* Score Breakdown Popup */}
              {showScoreBreakdown && (
                <ScoreBreakdownPopup
                  compositeScore={fwd.scoring.compositeScore}
                  quality={{ tier: fwd.scoring.qualityTier, label: fwd.scoring.qualityTier ? QUALITY_INFO[fwd.scoring.qualityTier]?.label : undefined }}
                  forwardScores={fwd.scoring?.piecewiseScores}
                  reverseScores={rev?.scoring?.piecewiseScores}
                  hasTemplate={!!sequence}
                  onClose={() => setShowScoreBreakdown(false)}
                />
              )}
            </>
          )}

          <table className="results-table">
            <thead>
              <tr>
                <th>Dir</th>
                <th>Tm</th>
                <th>Total Tm</th>
                <th>GC</th>
                <th>dG</th>
                <th
                  title="Lower penalty = better primer (legacy scoring)"
                  style={{ cursor: 'help' }}
                >
                  Penalty
                </th>
                <th
                  title="Calibrated quality score (0-100). Higher = better. Based on empirical validation with 81.9% F1 score."
                  style={{ cursor: 'help' }}
                >
                  Score
                </th>
                <th>Sequence (5' to 3')</th>
              </tr>
            </thead>
            <tbody>
              {fwd && (
                <tr>
                  <td className="dir fwd">FWD</td>
                  <td>{fwd.tm}</td>
                  <td>{fwd.tmTotal}</td>
                  <td>{fwd.gc}</td>
                  <td>{fwd.dg}</td>
                  <td>{fwd.scoring.penalty?.toFixed(1)}</td>
                  <td
                    style={{
                      color: fwd.scoring?.qualityTier ? QUALITY_INFO[fwd.scoring.qualityTier]?.color : undefined,
                      fontWeight: 'bold',
                    }}
                    title={fwd.scoring?.qualityTier ? QUALITY_INFO[fwd.scoring.qualityTier]?.description : undefined}
                  >
                    {fwd.scoring?.compositeScore ?? 'N/A'}
                  </td>
                  <td className="sequence">{fwd.seq}</td>
                </tr>
              )}
              {rev && (
                <tr>
                  <td className="dir rev">REV</td>
                  <td>{rev.tm}</td>
                  <td>{rev.tmTotal}</td>
                  <td>{rev.gc}</td>
                  <td>{rev.dg}</td>
                  <td>{rev.scoring.penalty?.toFixed(1)}</td>
                  <td
                    style={{
                      color: rev.scoring?.qualityTier ? QUALITY_INFO[rev.scoring.qualityTier]?.color : undefined,
                      fontWeight: 'bold',
                    }}
                    title={rev.scoring?.qualityTier ? QUALITY_INFO[rev.scoring.qualityTier]?.description : undefined}
                  >
                    {rev.scoring?.compositeScore ?? 'N/A'}
                  </td>
                  <td className="sequence">{rev.seq}</td>
                </tr>
              )}
            </tbody>
          </table>

          <div className="scoring-details">
            <h4
              title="Feature scores calibrated on Döring immunoglobulin PCR dataset (829 pairs). Higher = better for each metric."
              style={{ cursor: 'help' }}
            >
              Scoring Breakdown
            </h4>
            <div className="scoring-grid">
              {fwd && (
                <div className="scoring-card">
                  <h5>Forward Primer</h5>
                  <dl>
                    <dt>Length</dt>
                    <dd>{fwd.len} bp</dd>
                    {/* Piecewise scores (calibrated, 0-1 scale) */}
                    {fwd.scoring?.piecewiseScores && (
                      <>
                        <dt
                          title="Melting temperature score. Optimal: 55-60°C"
                          style={{ cursor: 'help' }}
                        >
                          Tm Score
                        </dt>
                        <dd>{(fwd.scoring.piecewiseScores.tm * 100).toFixed(0)}%</dd>
                        <dt
                          title="GC content score. Optimal: 40-60%"
                          style={{ cursor: 'help' }}
                        >
                          GC Score
                        </dt>
                        <dd>{(fwd.scoring.piecewiseScores.gc * 100).toFixed(0)}%</dd>
                        <dt
                          title="3' terminal binding stability. Critical for priming efficiency."
                          style={{ cursor: 'help' }}
                        >
                          3' ΔG Score
                        </dt>
                        <dd>{(fwd.scoring.piecewiseScores.terminal3DG * 100).toFixed(0)}%</dd>
                        <dt
                          title="Off-target specificity. Most predictive of PCR success (+0.515 correlation)."
                          style={{ cursor: 'help' }}
                        >
                          Specificity
                        </dt>
                        <dd>{(fwd.scoring.piecewiseScores.offTarget * 100).toFixed(0)}%</dd>
                        <dt
                          title="Hairpin formation risk. Lower stability = better."
                          style={{ cursor: 'help' }}
                        >
                          Hairpin Score
                        </dt>
                        <dd>{(fwd.scoring.piecewiseScores.hairpin * 100).toFixed(0)}%</dd>
                        <dt
                          title="Self-dimer formation risk."
                          style={{ cursor: 'help' }}
                        >
                          Homodimer Score
                        </dt>
                        <dd>{(fwd.scoring.piecewiseScores.homodimer * 100).toFixed(0)}%</dd>
                      </>
                    )}
                    <dt>Off-target Count</dt>
                    <dd>{fwd.offTargetCount}</dd>
                    <dt className="total">Legacy Penalty</dt>
                    <dd className="total">{fwd.scoring.penalty?.toFixed(1)}</dd>
                  </dl>
                </div>
              )}
              {rev && (
                <div className="scoring-card">
                  <h5>Reverse Primer</h5>
                  <dl>
                    <dt>Length</dt>
                    <dd>{rev.len} bp</dd>
                    {/* Piecewise scores (calibrated, 0-1 scale) */}
                    {rev.scoring?.piecewiseScores && (
                      <>
                        <dt
                          title="Melting temperature score. Optimal: 55-60°C"
                          style={{ cursor: 'help' }}
                        >
                          Tm Score
                        </dt>
                        <dd>{(rev.scoring.piecewiseScores.tm * 100).toFixed(0)}%</dd>
                        <dt
                          title="GC content score. Optimal: 40-60%"
                          style={{ cursor: 'help' }}
                        >
                          GC Score
                        </dt>
                        <dd>{(rev.scoring.piecewiseScores.gc * 100).toFixed(0)}%</dd>
                        <dt
                          title="3' terminal binding stability. Critical for priming efficiency."
                          style={{ cursor: 'help' }}
                        >
                          3' ΔG Score
                        </dt>
                        <dd>{(rev.scoring.piecewiseScores.terminal3DG * 100).toFixed(0)}%</dd>
                        <dt
                          title="Off-target specificity. Most predictive of PCR success (+0.515 correlation)."
                          style={{ cursor: 'help' }}
                        >
                          Specificity
                        </dt>
                        <dd>{(rev.scoring.piecewiseScores.offTarget * 100).toFixed(0)}%</dd>
                        <dt
                          title="Hairpin formation risk. Lower stability = better."
                          style={{ cursor: 'help' }}
                        >
                          Hairpin Score
                        </dt>
                        <dd>{(rev.scoring.piecewiseScores.hairpin * 100).toFixed(0)}%</dd>
                        <dt
                          title="Self-dimer formation risk."
                          style={{ cursor: 'help' }}
                        >
                          Homodimer Score
                        </dt>
                        <dd>{(rev.scoring.piecewiseScores.homodimer * 100).toFixed(0)}%</dd>
                        {rev.scoring.piecewiseScores.heterodimer !== undefined && (
                          <>
                            <dt
                              title="Primer-primer dimer risk between forward and reverse."
                              style={{ cursor: 'help' }}
                            >
                              Heterodimer Score
                            </dt>
                            <dd>{(rev.scoring.piecewiseScores.heterodimer * 100).toFixed(0)}%</dd>
                          </>
                        )}
                        {rev.scoring.piecewiseScores.tmDiff !== undefined && (
                          <>
                            <dt
                              title="Tm difference between primers. Has 'little effect' on PCR success."
                              style={{ cursor: 'help' }}
                            >
                              Tm Match Score
                            </dt>
                            <dd>{(rev.scoring.piecewiseScores.tmDiff * 100).toFixed(0)}%</dd>
                          </>
                        )}
                      </>
                    )}
                    <dt>Off-target Count</dt>
                    <dd>{rev.offTargetCount}</dd>
                    <dt className="total">Legacy Penalty</dt>
                    <dd className="total">{rev.scoring.penalty?.toFixed(1)}</dd>
                  </dl>
                </div>
              )}
            </div>
          </div>
        </>
      )}

      <div className="copy-section">
        <h4>Quick Copy</h4>
        <div className="copy-buttons">
          {fwd && (
            <button
              type="button"
              onClick={() => navigator.clipboard.writeText(fwd.seq)}
              title="Copy FWD primer"
            >
              Copy FWD: {fwd.seq.slice(0, 15)}...
            </button>
          )}
          {rev && (
            <button
              type="button"
              onClick={() => navigator.clipboard.writeText(rev.seq)}
              title="Copy REV primer"
            >
              Copy REV: {rev.seq.slice(0, 15)}...
            </button>
          )}
        </div>
      </div>

      {sequence && (
        <div className="viewer-section">
          <div className="viewer-header">
            <h4>Sequence Visualization</h4>
            <div className="viewer-controls">
              <button
                type="button"
                className={viewerMode === 'linear' ? 'active' : ''}
                onClick={() => setViewerMode('linear')}
              >
                Linear
              </button>
              <button
                type="button"
                className={viewerMode === 'circular' ? 'active' : ''}
                onClick={() => setViewerMode('circular')}
              >
                Circular
              </button>
              <button
                type="button"
                className={viewerMode === 'both' ? 'active' : ''}
                onClick={() => setViewerMode('both')}
              >
                Both
              </button>
            </div>
          </div>
          <SequenceViewer
            sequence={sequence}
            name="Template"
            fwdPrimer={fwd.seq as any}
            revPrimer={rev?.seq as any}
            addFwd={options.addFwd as any}
            addRev={options.addRev as any}
            viewer={viewerMode}
            height={viewerMode === 'both' ? 400 : 300}
          />
        </div>
      )}
    </div>
  );
}
