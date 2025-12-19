/**
 * EnhancedAnalysisSection - Shared component for detailed primer analysis display
 *
 * Extracted from UnifiedPrimerDesigner for reuse across:
 * - UnifiedPrimerDesigner (amplification/mutagenesis)
 * - EnhancedScorer (standalone scoring mode)
 *
 * Displays:
 * - Heterodimer analysis with 3' involvement detection
 * - Secondary structure (hairpin/self-dimer) with Zuker algorithm
 * - Off-target binding summary
 * - Tm methods comparison
 * - Warnings with severity levels
 */

import { useState } from 'react';

// Type definitions for analysis data structures
interface ThermodynamicsData {
  hairpinDG?: number;
  homodimerDG?: number;
}

interface AnnealingRegion {
  tm?: number;
  length?: number;
}

interface PrimerAnalysis {
  hairpinDG?: number;
  selfDimerDG?: number;
  thermodynamics?: ThermodynamicsData;
  isAssemblyPrimer?: boolean;
  annealingRegion?: AnnealingRegion;
  fullPrimerTm?: number;
  tm?: number;
  tmComparison?: TmComparison;
}

interface HeterodimerWarning {
  message: string;
  severity: string;
  tooltip?: string;
}

interface HeterodimerAnalysis {
  isExpectedOverlap?: boolean;
  dimerType?: string;
  severity?: 'critical' | 'warning' | 'ok';
  overlapLength?: number;
  heterodimerDG?: number;
  involves3Prime1?: boolean;
  involves3Prime2?: boolean;
  maxConsecutive?: number;
  warnings?: HeterodimerWarning[];
}

interface OffTargetData {
  offTargetCount?: number;
}

interface OffTargets {
  forward?: OffTargetData;
  reverse?: OffTargetData;
}

interface TmMethods {
  q5?: number;
  general?: number;
  simple?: number;
}

interface TmComparison {
  methods?: TmMethods;
}

interface Warning {
  primer: string;
  type: string;
  message: string;
  severity: string;
}

interface PairAnalysis {
  heterodimer?: HeterodimerAnalysis;
}

interface AnalysisData {
  forward?: PrimerAnalysis;
  reverse?: PrimerAnalysis;
  heterodimer?: HeterodimerAnalysis;
  pair?: PairAnalysis;
  offTargets?: OffTargets;
  fwdTmComparison?: TmComparison;
  warnings?: Warning[];
}

interface EnhancedAnalysisSectionProps {
  analysis: AnalysisData | null;
  collapsible?: boolean;
  defaultCollapsed?: boolean;
  sectionId?: string;
}

export default function EnhancedAnalysisSection({
  analysis,
  collapsible = true,
  defaultCollapsed = false,
  sectionId = 'section-analysis',
}: EnhancedAnalysisSectionProps) {
  const [isCollapsed, setIsCollapsed] = useState<boolean>(defaultCollapsed);
  const [showThermodynamics, setShowThermodynamics] = useState<boolean>(false);

  if (!analysis) {
    return null;
  }

  // Support both enhancedAnalysis format and direct analyzePrimers() output
  const forward = analysis.forward || {};
  const reverse = analysis.reverse || {};
  const heterodimer = analysis.heterodimer || analysis.pair?.heterodimer || {};
  const offTargets = analysis.offTargets || {};
  const fwdTmComparison = analysis.fwdTmComparison || analysis.forward?.tmComparison || null;
  const warnings = analysis.warnings || [];

  const toggleCollapse = (): void => {
    if (collapsible) {
      setIsCollapsed(!isCollapsed);
    }
  };

  return (
    <div id={sectionId} className="enhanced-analysis-section collapsible-section">
      {collapsible ? (
        <button
          type="button"
          className="section-header collapsible"
          onClick={toggleCollapse}
        >
          <h4>Analysis Details</h4>
          <span className="collapse-icon">{isCollapsed ? '▶' : '▼'}</span>
        </button>
      ) : (
        <div className="section-header">
          <h4>Analysis Details</h4>
        </div>
      )}

      {!isCollapsed && (
        <div className="section-content">
          {/* Heterodimer Analysis */}
          <div className="analysis-card">
            <div className="analysis-header">
              <span
                className="analysis-title"
                title="Analysis of potential binding between forward and reverse primers"
              >
                {heterodimer.isExpectedOverlap
                  ? 'Primer Overlap (QuikChange)'
                  : heterodimer.dimerType === '3prime_extensible'
                  ? "3' Extensible Dimer Analysis"
                  : 'Heterodimer Analysis'}
              </span>
              <span
                className={`status-badge ${
                  heterodimer.isExpectedOverlap
                    ? 'good'
                    : heterodimer.severity === 'critical'
                    ? 'critical'
                    : heterodimer.severity === 'warning'
                    ? 'warning'
                    : 'good'
                }`}
              >
                {heterodimer.isExpectedOverlap
                  ? 'Expected'
                  : heterodimer.severity === 'critical'
                  ? 'CRITICAL'
                  : heterodimer.severity === 'warning'
                  ? 'Warning'
                  : 'OK'}
              </span>
            </div>
            <div className="analysis-content">
              {heterodimer.isExpectedOverlap ? (
                <>
                  <div
                    className="stat-row"
                    title="For QuikChange mutagenesis, primers are designed as reverse complements. Full overlap is expected and correct."
                  >
                    <span className="label">Overlap:</span>
                    <span className="value good">
                      {heterodimer.overlapLength} bp (100%)
                    </span>
                  </div>
                  <div className="info-note">
                    Full primer complementarity is expected for QuikChange-style
                    overlapping mutagenesis.
                  </div>
                </>
              ) : (
                <>
                  <div
                    className="stat-row"
                    title="Gibbs free energy of heterodimer formation. 3' end involvement is more critical than ΔG alone."
                  >
                    <span className="label">Heterodimer ΔG:</span>
                    <span
                      className={`value ${
                        heterodimer.severity === 'critical'
                          ? 'critical'
                          : (heterodimer.heterodimerDG ?? 0) < -6
                          ? 'warning'
                          : ''
                      }`}
                    >
                      {heterodimer.heterodimerDG?.toFixed(1) || 'N/A'} kcal/mol
                    </span>
                  </div>
                  {/* Show 3' involvement details */}
                  {(heterodimer.involves3Prime1 || heterodimer.involves3Prime2) && (
                    <div
                      className="stat-row"
                      title="3' end involvement is the most critical factor for primer-dimer artifacts"
                    >
                      <span className="label">3' Involvement:</span>
                      <span
                        className={`value ${
                          heterodimer.severity === 'critical' ? 'critical' : 'warning'
                        }`}
                      >
                        {heterodimer.involves3Prime1 && heterodimer.involves3Prime2
                          ? 'Both primers'
                          : heterodimer.involves3Prime1
                          ? 'Forward primer'
                          : 'Reverse primer'}
                        {heterodimer.maxConsecutive &&
                          ` (${heterodimer.maxConsecutive} consecutive bp)`}
                      </span>
                    </div>
                  )}
                  {heterodimer.warnings && heterodimer.warnings.length > 0 && (
                    <div className="mini-warnings">
                      {heterodimer.warnings.map((w, i) => (
                        <div
                          key={i}
                          className={`mini-warning ${w.severity}`}
                          title={w.tooltip || ''}
                        >
                          {w.message}
                        </div>
                      ))}
                    </div>
                  )}
                </>
              )}
            </div>
          </div>

          {/* Secondary Structure */}
          <div className="analysis-card">
            <div className="analysis-header">
              <span
                className="analysis-title"
                title="Self-folding analysis using Zuker algorithm"
              >
                Secondary Structure (Zuker Algorithm)
              </span>
            </div>
            <div className="analysis-content two-column">
              <div className="column">
                <span className="column-title">Forward</span>
                <div
                  className="stat-row"
                  title="Hairpin/folding energy using accurate Zuker algorithm. Values > -3 kcal/mol are ideal."
                >
                  <span className="label">Hairpin ΔG:</span>
                  <span
                    className={`value ${
                      (forward.hairpinDG ?? forward.thermodynamics?.hairpinDG ?? 0) < -3
                        ? (forward.hairpinDG ?? forward.thermodynamics?.hairpinDG ?? 0) < -6
                          ? 'critical'
                          : 'warning'
                        : ''
                    }`}
                  >
                    {(forward.hairpinDG ?? forward.thermodynamics?.hairpinDG)?.toFixed(1) ?? 'N/A'} kcal/mol
                  </span>
                </div>
                <div
                  className="stat-row"
                  title="Energy of self-dimer (primer binding to itself). Values > -6 kcal/mol are acceptable."
                >
                  <span className="label">Self-dimer ΔG:</span>
                  <span
                    className={`value ${
                      (forward.selfDimerDG ?? forward.thermodynamics?.homodimerDG ?? 0) < -6
                        ? (forward.selfDimerDG ?? forward.thermodynamics?.homodimerDG ?? 0) < -8
                          ? 'critical'
                          : 'warning'
                        : ''
                    }`}
                  >
                    {(forward.selfDimerDG ?? forward.thermodynamics?.homodimerDG)?.toFixed(1) ?? 'N/A'} kcal/mol
                  </span>
                </div>
              </div>
              {reverse && Object.keys(reverse).length > 0 && (
                <div className="column">
                  <span className="column-title">Reverse</span>
                  <div
                    className="stat-row"
                    title="Hairpin/folding energy using accurate Zuker algorithm. Values > -3 kcal/mol are ideal."
                  >
                    <span className="label">Hairpin ΔG:</span>
                    <span
                      className={`value ${
                        (reverse.hairpinDG ?? reverse.thermodynamics?.hairpinDG ?? 0) < -3
                          ? (reverse.hairpinDG ?? reverse.thermodynamics?.hairpinDG ?? 0) < -6
                            ? 'critical'
                            : 'warning'
                          : ''
                      }`}
                    >
                      {(reverse.hairpinDG ?? reverse.thermodynamics?.hairpinDG)?.toFixed(1) ?? 'N/A'} kcal/mol
                    </span>
                  </div>
                  <div
                    className="stat-row"
                    title="Energy of self-dimer (primer binding to itself). Values > -6 kcal/mol are acceptable."
                  >
                    <span className="label">Self-dimer ΔG:</span>
                    <span
                      className={`value ${
                        (reverse.selfDimerDG ?? reverse.thermodynamics?.homodimerDG ?? 0) < -6
                          ? (reverse.selfDimerDG ?? reverse.thermodynamics?.homodimerDG ?? 0) < -8
                            ? 'critical'
                            : 'warning'
                          : ''
                      }`}
                    >
                      {(reverse.selfDimerDG ?? reverse.thermodynamics?.homodimerDG)?.toFixed(1) ?? 'N/A'} kcal/mol
                    </span>
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Assembly Primer Tm Breakdown */}
          {(forward.isAssemblyPrimer || reverse?.isAssemblyPrimer) && (
            <div className="analysis-card">
              <div className="analysis-header">
                <span
                  className="analysis-title"
                  title="For assembly primers, we show both the annealing region Tm (what matters for PCR) and full primer Tm"
                >
                  Tm Analysis (Assembly Primer)
                </span>
                <span className="status-badge good">Key Info</span>
              </div>
              <div className="analysis-content two-column">
                {forward.isAssemblyPrimer && forward.annealingRegion && (
                  <div className="column">
                    <span className="column-title">Forward</span>
                    <div
                      className="stat-row"
                      title="This is the Tm of the annealing region - use this for calculating your annealing temperature"
                    >
                      <span className="label">Annealing Tm:</span>
                      <span className="value highlight">
                        {forward.annealingRegion.tm?.toFixed(1) || forward.tm?.toFixed(1)}°C
                      </span>
                    </div>
                    <div
                      className="stat-row"
                      title="Full primer Tm is higher but NOT relevant for annealing temperature selection"
                    >
                      <span className="label">Full Primer Tm:</span>
                      <span className="value" style={{ color: '#6b7280' }}>
                        {forward.fullPrimerTm?.toFixed(1) || forward.tm?.toFixed(1)}°C
                      </span>
                    </div>
                    <div className="stat-row">
                      <span className="label">Annealing Length:</span>
                      <span className="value">
                        {forward.annealingRegion.length} bp
                      </span>
                    </div>
                  </div>
                )}
                {reverse?.isAssemblyPrimer && reverse?.annealingRegion && (
                  <div className="column">
                    <span className="column-title">Reverse</span>
                    <div
                      className="stat-row"
                      title="This is the Tm of the annealing region - use this for calculating your annealing temperature"
                    >
                      <span className="label">Annealing Tm:</span>
                      <span className="value highlight">
                        {reverse.annealingRegion.tm?.toFixed(1) || reverse.tm?.toFixed(1)}°C
                      </span>
                    </div>
                    <div
                      className="stat-row"
                      title="Full primer Tm is higher but NOT relevant for annealing temperature selection"
                    >
                      <span className="label">Full Primer Tm:</span>
                      <span className="value" style={{ color: '#6b7280' }}>
                        {reverse.fullPrimerTm?.toFixed(1) || reverse.tm?.toFixed(1)}°C
                      </span>
                    </div>
                    <div className="stat-row">
                      <span className="label">Annealing Length:</span>
                      <span className="value">
                        {reverse.annealingRegion.length} bp
                      </span>
                    </div>
                  </div>
                )}
              </div>
              <div className="info-note" style={{ marginTop: '8px', fontSize: '12px', color: '#6b7280' }}>
                Use the <strong>Annealing Tm</strong> (not full primer Tm) for calculating your PCR annealing temperature.
                Recommended: Annealing Tm minus 5°C.
              </div>
            </div>
          )}

          {/* Off-Target Summary */}
          {(offTargets.forward || offTargets.reverse) && (
            <div className="analysis-card">
              <div className="analysis-header">
                <span className="analysis-title">Off-Target Binding</span>
              </div>
              <div className="analysis-content two-column">
                <div className="column">
                  <span className="column-title">Forward</span>
                  <div className="stat-row">
                    <span className="label">Sites found:</span>
                    <span
                      className={`value ${
                        (offTargets.forward?.offTargetCount ?? 0) > 2 ? 'warning' : ''
                      }`}
                    >
                      {offTargets.forward?.offTargetCount || 0}
                    </span>
                  </div>
                </div>
                <div className="column">
                  <span className="column-title">Reverse</span>
                  <div className="stat-row">
                    <span className="label">Sites found:</span>
                    <span
                      className={`value ${
                        (offTargets.reverse?.offTargetCount ?? 0) > 2 ? 'warning' : ''
                      }`}
                    >
                      {offTargets.reverse?.offTargetCount || 0}
                    </span>
                  </div>
                </div>
              </div>
            </div>
          )}

          {/* Tm Comparison */}
          {fwdTmComparison && (
            <div className="analysis-card">
              <button
                type="button"
                className="analysis-header clickable"
                onClick={() => setShowThermodynamics(!showThermodynamics)}
              >
                <span className="analysis-title">Tm Methods Comparison</span>
                <span className="toggle-icon">{showThermodynamics ? '▼' : '▶'}</span>
              </button>
              {showThermodynamics && (
                <div className="analysis-content">
                  <div className="tm-comparison">
                    <div className="stat-row">
                      <span className="label">Q5 Tm:</span>
                      <span className="value highlight">
                        {fwdTmComparison.methods?.q5}°C
                      </span>
                    </div>
                    <div className="stat-row">
                      <span className="label">General NN Tm:</span>
                      <span className="value">
                        {fwdTmComparison.methods?.general?.toFixed(1)}°C
                      </span>
                    </div>
                    <div className="stat-row">
                      <span className="label">Simple Formula Tm:</span>
                      <span className="value">{fwdTmComparison.methods?.simple}°C</span>
                    </div>
                  </div>
                </div>
              )}
            </div>
          )}

          {/* Warnings */}
          {warnings && warnings.length > 0 && (
            <div className="all-warnings">
              <h5>Analysis Warnings ({warnings.length})</h5>
              {warnings.map((warning, idx) => (
                <div key={idx} className={`warning-item ${warning.severity}`}>
                  <span className="warning-badge">{warning.primer}</span>
                  <span className="warning-type">{warning.type}</span>
                  <span className="warning-message">{warning.message}</span>
                </div>
              ))}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
