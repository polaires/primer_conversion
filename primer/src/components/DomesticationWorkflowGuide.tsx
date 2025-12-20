/**
 * Domestication Workflow Guide Component
 *
 * Provides a visual step-by-step workflow guide for Golden Gate domestication.
 * Integrates with the state-of-the-art domestication primer workflow system.
 *
 * Features:
 * - Strategy comparison (Silent Mutation vs Mutagenic Junction)
 * - NEB fidelity-based overhang scoring
 * - Ready-to-order primer specifications
 * - Interactive workflow steps with progress tracking
 */

import { useState, useMemo, useCallback, FC } from 'react';
import {
  runDomesticationWorkflow,
  analyzeSequenceForDomestication,
} from '../lib/repp/domestication-primer-workflow.js';
import {
  GOLDEN_GATE_ENZYMES,
  findInternalSites,
} from '../lib/repp/goldengate.js';

// ============================================================================
// TYPE DEFINITIONS
// ============================================================================

type WorkflowStepType = 'overview' | 'analysis' | 'strategy' | 'primers' | 'summary';
type StrategyType = 'none' | 'silent_mutation' | 'mutagenic_junction' | 'hybrid' | 'alternative_enzyme';

interface Site {
  position: number;
  [key: string]: unknown;
}

interface InternalSites {
  hasSites: boolean;
  count: number;
  sites: Site[];
}

interface Analysis {
  error?: string;
  siteCount?: number;
  silentMutationPossible?: boolean;
  recommendedStrategy?: string;
  [key: string]: unknown;
}

interface MutationData {
  mutation?: {
    sequencePosition?: number;
    position?: number;
    codonStart?: number;
    originalCodon?: string;
    newCodon?: string;
    aminoAcid?: string;
    score?: number;
  };
  sequencePosition?: number;
  position?: number;
  codonStart?: number;
  originalCodon?: string;
  newCodon?: string;
  aminoAcid?: string;
  score?: number;
}

interface JunctionData {
  junctionPosition?: number;
  position?: number;
  overhang?: string;
  fidelity?: {
    singleOverhang?: number;
  } | number;
}

interface FidelityData {
  assembly?: number;
  [key: string]: unknown;
}

interface DomesticationData {
  mutations?: MutationData[];
  junctions?: JunctionData[];
  fidelity?: FidelityData | number;
  additionalFragments?: number;
}

interface ValidationData {
  isValid: boolean;
  warnings?: Array<{ message?: string } | string>;
}

interface ThermodynamicsData {
  hairpinDG?: number;
  homodimerDG?: number;
  terminal3DG?: number;
}

interface PrimerData {
  sequence?: string;
  length?: number;
  homologyTm?: number;
  gcPercent?: string;
  gc?: number;
  thermodynamics?: ThermodynamicsData;
  score?: number;
  qualityTier?: string;
  issues?: string[];
  warnings?: string[];
}

interface HeterodimerData {
  dG?: number;
  risk?: string;
}

interface PrimerPairData {
  type?: string;
  name?: string;
  forward?: PrimerData;
  reverse?: PrimerData;
  junctionPosition?: number;
  overhang?: string;
  pairQuality?: number;
  tmDifference?: number;
  heterodimer?: HeterodimerData;
}

interface OrderListItem {
  name: string;
  sequence: string;
  length: number;
  tm?: number;
}

interface PrimerSummaryData {
  totalPrimers: number;
  totalLength: number;
  estimatedCost?: number;
  totalCost?: number;
  orderList?: OrderListItem[];
}

interface WorkflowStep {
  title: string;
  description: string;
  materials?: string[];
}

interface WorkflowGuideData {
  steps?: WorkflowStep[];
}

interface WorkflowResult {
  success: boolean;
  message?: string;
  strategy: StrategyType;
  domestication?: DomesticationData;
  analysis?: Analysis;
  validation?: ValidationData;
  primers?: PrimerPairData[];
  primerSummary?: PrimerSummaryData;
  workflowGuide?: WorkflowGuideData;
}

interface DomesticationWorkflowGuideProps {
  sequence: string;
  enzyme?: string;
  existingOverhangs?: string[];
  onWorkflowComplete?: (result: WorkflowResult) => void;
  onCancel?: () => void;
}

interface WorkflowProgressProps {
  steps: string[];
  currentStep: WorkflowStepType;
}

interface OverviewStepProps {
  sequence: string;
  enzyme: string;
  internalSites: InternalSites;
  onContinue: () => void;
}

interface AnalysisStepProps {
  sequence: string;
  analysis: Analysis | null;
  selectedFrame: number;
  onFrameChange: (frame: number) => void;
  selectedOrganism: string;
  onOrganismChange: (organism: string) => void;
  isCoding: boolean;
  onCodingChange: (isCoding: boolean) => void;
  onRunWorkflow: () => void;
  isProcessing: boolean;
  onBack: () => void;
}

interface StrategyStepProps {
  result: WorkflowResult;
  onContinue: () => void;
  onBack: () => void;
}

interface PrimersStepProps {
  result: WorkflowResult;
  onContinue: () => void;
  onBack: () => void;
}

interface PrimerPairCardProps {
  primerPair: PrimerPairData;
  index: number;
}

interface SummaryStepProps {
  result: WorkflowResult;
  onAccept: () => void;
  onBack: () => void;
}

interface EnzymeInfo {
  recognition?: string;
  [key: string]: unknown;
}

// ============================================================================
// MAIN COMPONENT
// ============================================================================

export const DomesticationWorkflowGuide: FC<DomesticationWorkflowGuideProps> = ({
  sequence,
  enzyme = 'BsaI',
  existingOverhangs = [],
  onWorkflowComplete,
  onCancel,
}) => {
  // State
  const [step, setStep] = useState<WorkflowStepType>('overview');
  const [workflowResult, setWorkflowResult] = useState<WorkflowResult | null>(null);
  const [selectedFrame, setSelectedFrame] = useState<number>(0);
  const [selectedOrganism, setSelectedOrganism] = useState<string>('ecoli');
  const [isCoding, setIsCoding] = useState<boolean>(true);
  const [isProcessing, setIsProcessing] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);

  // Initial analysis
  const analysis = useMemo((): Analysis | null => {
    if (!sequence) return null;
    try {
      return (analyzeSequenceForDomestication as any)(sequence, enzyme, {
        frame: selectedFrame,
        isCodingSequence: isCoding,
        organism: selectedOrganism,
        existingOverhangs,
      }) as Analysis;
    } catch (err) {
      return { error: (err as Error).message };
    }
  }, [sequence, enzyme, selectedFrame, isCoding, selectedOrganism, existingOverhangs]);

  // Run workflow
  const runWorkflow = useCallback(async () => {
    if (!sequence) return;

    setIsProcessing(true);
    setError(null);

    try {
      const result = (runDomesticationWorkflow as any)(sequence, enzyme, {
        frame: selectedFrame,
        isCodingSequence: isCoding,
        organism: selectedOrganism,
        existingOverhangs,
        includeWorkflowGuide: true,
      }) as WorkflowResult;

      setWorkflowResult(result);
      setStep('strategy');
    } catch (err) {
      setError((err as Error).message);
    } finally {
      setIsProcessing(false);
    }
  }, [sequence, enzyme, selectedFrame, isCoding, selectedOrganism, existingOverhangs]);

  // Accept results and proceed
  const acceptAndContinue = useCallback(() => {
    if (workflowResult?.success && onWorkflowComplete) {
      onWorkflowComplete(workflowResult);
    }
  }, [workflowResult, onWorkflowComplete]);

  // Render
  // Detect internal enzyme sites in the sequence
  const internalSites = useMemo((): InternalSites => {
    if (!sequence) return { hasSites: false, count: 0, sites: [] };
    const result = findInternalSites(sequence, enzyme);
    return {
      hasSites: result.hasSites,
      count: result.count,
      sites: result.sites.map(s => ({
        position: s.position,
        sequence: s.sequence,
        orientation: s.orientation,
      })),
    };
  }, [sequence, enzyme]);

  if (!sequence) {
    return (
      <div className="workflow-guide empty">
        <p>No sequence provided</p>
      </div>
    );
  }

  if (!internalSites?.hasSites) {
    return (
      <div className="workflow-guide no-sites">
        <div className="success-icon-large">✓</div>
        <h3>Sequence Compatible</h3>
        <p>No internal {enzyme} sites found. Sequence is ready for Golden Gate assembly.</p>
        <button className="btn-primary" onClick={onCancel}>Continue</button>
      </div>
    );
  }

  return (
    <div className="domestication-workflow-guide">
      {/* Progress Steps */}
      <WorkflowProgress
        steps={['Overview', 'Analysis', 'Strategy', 'Primers', 'Summary']}
        currentStep={step}
      />

      {/* Error Display */}
      {error && (
        <div className="workflow-error">
          <span className="error-icon">⚠</span>
          <span>{error}</span>
          <button onClick={() => setError(null)}>×</button>
        </div>
      )}

      {/* Step Content */}
      {step === 'overview' && (
        <OverviewStep
          sequence={sequence}
          enzyme={enzyme}
          internalSites={internalSites}
          onContinue={() => setStep('analysis')}
        />
      )}

      {step === 'analysis' && (
        <AnalysisStep
          sequence={sequence}
          analysis={analysis}
          selectedFrame={selectedFrame}
          onFrameChange={setSelectedFrame}
          selectedOrganism={selectedOrganism}
          onOrganismChange={setSelectedOrganism}
          isCoding={isCoding}
          onCodingChange={setIsCoding}
          onRunWorkflow={runWorkflow}
          isProcessing={isProcessing}
          onBack={() => setStep('overview')}
        />
      )}

      {step === 'strategy' && workflowResult && (
        <StrategyStep
          result={workflowResult}
          onContinue={() => setStep('primers')}
          onBack={() => setStep('analysis')}
        />
      )}

      {step === 'primers' && workflowResult && (
        <PrimersStep
          result={workflowResult}
          onContinue={() => setStep('summary')}
          onBack={() => setStep('strategy')}
        />
      )}

      {step === 'summary' && workflowResult && (
        <SummaryStep
          result={workflowResult}
          onAccept={acceptAndContinue}
          onBack={() => setStep('primers')}
        />
      )}
    </div>
  );
};

// ============================================================================
// STEP COMPONENTS
// ============================================================================

const WorkflowProgress: FC<WorkflowProgressProps> = ({ steps, currentStep }) => {
  const stepMap: Record<WorkflowStepType, number> = {
    overview: 0,
    analysis: 1,
    strategy: 2,
    primers: 3,
    summary: 4,
  };
  const current = stepMap[currentStep] ?? 0;

  return (
    <div className="workflow-progress">
      {steps.map((label, i) => (
        <div
          key={label}
          className={`progress-item ${i < current ? 'done' : ''} ${i === current ? 'active' : ''}`}
        >
          <div className="progress-dot">{i < current ? '✓' : i + 1}</div>
          <div className="progress-label">{label}</div>
        </div>
      ))}
    </div>
  );
};

const OverviewStep: FC<OverviewStepProps> = ({ sequence, enzyme, internalSites, onContinue }) => {
  const enzymeInfo = (GOLDEN_GATE_ENZYMES as any)[enzyme];

  return (
    <div className="workflow-step overview-step">
      <h2>Domestication Required</h2>
      <p className="step-intro">
        Your sequence contains internal {enzyme} recognition sites that must be removed
        for successful Golden Gate assembly.
      </p>

      <div className="info-grid">
        <div className="info-card warning">
          <div className="info-icon">⚠</div>
          <div className="info-content">
            <div className="info-title">{internalSites.count} Internal Site{internalSites.count !== 1 ? 's' : ''}</div>
            <div className="info-detail">
              Positions: {internalSites.sites.map(s => s.position + 1).join(', ')}
            </div>
          </div>
        </div>

        <div className="info-card">
          <div className="info-icon">🧬</div>
          <div className="info-content">
            <div className="info-title">Sequence Length</div>
            <div className="info-detail">{sequence.length.toLocaleString()} bp</div>
          </div>
        </div>

        <div className="info-card">
          <div className="info-icon">✂️</div>
          <div className="info-content">
            <div className="info-title">{enzyme} Enzyme</div>
            <div className="info-detail">
              Recognition: {enzymeInfo?.recognition || 'N/A'}
            </div>
          </div>
        </div>
      </div>

      <div className="workflow-preview">
        <h3>Available Strategies</h3>
        <div className="strategy-cards">
          <div className="strategy-card preferred">
            <div className="strategy-badge">Preferred</div>
            <h4>Silent Mutation</h4>
            <p>Introduce synonymous codon changes to eliminate sites while preserving protein sequence.</p>
            <ul>
              <li>One-pot compatible</li>
              <li>No additional fragments</li>
              <li>Highest fidelity</li>
            </ul>
          </div>

          <div className="strategy-card">
            <h4>Mutagenic Junction</h4>
            <p>Create assembly junctions at site locations, splitting into additional fragments.</p>
            <ul>
              <li>Works for non-coding regions</li>
              <li>NEB fidelity-optimized overhangs</li>
              <li>Wider search radius (±50bp)</li>
            </ul>
          </div>
        </div>
      </div>

      <div className="step-actions">
        <button className="btn-primary" onClick={onContinue}>
          Configure Options →
        </button>
      </div>
    </div>
  );
};

const AnalysisStep: FC<AnalysisStepProps> = ({
  sequence: _sequence,
  analysis,
  selectedFrame,
  onFrameChange,
  selectedOrganism,
  onOrganismChange,
  isCoding,
  onCodingChange,
  onRunWorkflow,
  isProcessing,
  onBack,
}) => {
  return (
    <div className="workflow-step analysis-step">
      <h2>Configuration</h2>
      <p className="step-intro">
        Configure the domestication parameters for optimal results.
      </p>

      <div className="config-section">
        <label className="config-checkbox">
          <input
            type="checkbox"
            checked={isCoding}
            onChange={(e: React.ChangeEvent<HTMLInputElement>) => onCodingChange(e.target.checked)}
          />
          <span>Coding Sequence (CDS)</span>
        </label>
        <p className="config-help">
          Enable for protein-coding sequences to use synonymous mutations
        </p>
      </div>

      {isCoding && (
        <>
          <div className="config-section">
            <label>Reading Frame</label>
            <div className="frame-options">
              {[0, 1, 2].map(frame => (
                <button
                  key={frame}
                  className={`frame-btn ${selectedFrame === frame ? 'active' : ''}`}
                  onClick={() => onFrameChange(frame)}
                >
                  Frame {frame}
                </button>
              ))}
            </div>
            <p className="config-help">
              Select the correct reading frame for accurate codon analysis
            </p>
          </div>

          <div className="config-section">
            <label>Target Organism</label>
            <select
              value={selectedOrganism}
              onChange={(e: React.ChangeEvent<HTMLSelectElement>) => onOrganismChange(e.target.value)}
            >
              <option value="ecoli">E. coli</option>
              <option value="yeast">S. cerevisiae (Yeast)</option>
            </select>
            <p className="config-help">
              Codon usage tables for synonymous mutation selection
            </p>
          </div>
        </>
      )}

      {/* Analysis Preview */}
      {analysis && !analysis.error && (
        <div className="analysis-preview">
          <h3>Pre-Analysis</h3>
          <div className="analysis-stats">
            <div className="stat">
              <span className="stat-value">{analysis.siteCount}</span>
              <span className="stat-label">Sites to Remove</span>
            </div>
            {analysis.silentMutationPossible !== undefined && (
              <div className="stat">
                <span className="stat-value">
                  {analysis.silentMutationPossible ? '✓' : '✗'}
                </span>
                <span className="stat-label">Silent Mutation Available</span>
              </div>
            )}
            <div className="stat">
              <span className="stat-value">{analysis.recommendedStrategy || 'TBD'}</span>
              <span className="stat-label">Recommended Strategy</span>
            </div>
          </div>
        </div>
      )}

      {analysis?.error && (
        <div className="analysis-error">
          <span className="error-icon">⚠</span>
          <span>{analysis.error}</span>
        </div>
      )}

      <div className="step-actions">
        <button className="btn-secondary" onClick={onBack}>← Back</button>
        <button
          className="btn-primary"
          onClick={onRunWorkflow}
          disabled={isProcessing}
        >
          {isProcessing ? 'Analyzing...' : 'Run Workflow →'}
        </button>
      </div>
    </div>
  );
};

const StrategyStep: FC<StrategyStepProps> = ({ result, onContinue, onBack }) => {
  const { strategy, domestication, analysis: _analysis } = result;

  const strategyInfo: Record<string, { name: string; desc: string; icon: string }> = {
    none: { name: 'No Domestication', desc: 'Sequence is already compatible', icon: '✓' },
    silent_mutation: { name: 'Silent Mutations', desc: 'Synonymous codon changes', icon: '🧬' },
    mutagenic_junction: { name: 'Mutagenic Junctions', desc: 'Assembly junction splitting', icon: '✂️' },
    hybrid: { name: 'Hybrid Approach', desc: 'Combination of strategies', icon: '⚗️' },
    alternative_enzyme: { name: 'Alternative Enzyme', desc: 'Use different restriction enzyme', icon: '🔄' },
  };

  const info = strategyInfo[strategy] || strategyInfo.none;

  return (
    <div className="workflow-step strategy-step">
      <h2>Selected Strategy</h2>

      <div className="selected-strategy">
        <div className="strategy-icon">{info.icon}</div>
        <div className="strategy-info">
          <h3>{info.name}</h3>
          <p>{info.desc}</p>
        </div>
      </div>

      {/* Strategy Details */}
      {strategy === 'silent_mutation' && domestication?.mutations && (
        <div className="strategy-details">
          <h4>Mutations Applied</h4>
          <table className="mutations-table">
            <thead>
              <tr>
                <th>Position</th>
                <th>Original</th>
                <th>New</th>
                <th>Amino Acid</th>
                <th>Score</th>
              </tr>
            </thead>
            <tbody>
              {domestication.mutations.map((mut, i) => {
                // Handle both nested (mut.mutation.x) and flat (mut.x) structures
                const m = mut.mutation || mut;
                const position = m.sequencePosition ?? m.position ?? m.codonStart;
                const originalCodon = m.originalCodon || '';
                const newCodon = m.newCodon || '';
                const aminoAcid = m.aminoAcid || '';
                const score = m.score ?? mut.score;

                return (
                  <tr key={i}>
                    <td>{position !== undefined ? position + 1 : 'N/A'}</td>
                    <td><code>{originalCodon}</code></td>
                    <td><code>{newCodon}</code></td>
                    <td>{aminoAcid}</td>
                    <td>
                      <span className={`score-badge ${(score ?? 0) >= 80 ? 'good' : (score ?? 0) >= 60 ? 'ok' : 'poor'}`}>
                        {score?.toFixed?.(0) || 'N/A'}
                      </span>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      )}

      {(strategy === 'mutagenic_junction' || strategy === 'hybrid') && domestication?.junctions && domestication.junctions.length > 0 && (
        <div className="strategy-details">
          <h4>Junction Positions</h4>
          <table className="junctions-table">
            <thead>
              <tr>
                <th>Junction</th>
                <th>Position</th>
                <th>Overhang</th>
                <th>Fidelity</th>
              </tr>
            </thead>
            <tbody>
              {domestication.junctions.map((jnc, i) => {
                // Handle both position and junctionPosition field names
                const pos = jnc.junctionPosition ?? jnc.position ?? 0;
                // Handle fidelity from junction or from nested fidelity object
                const fidelity = typeof jnc.fidelity === 'object'
                  ? jnc.fidelity?.singleOverhang ?? 0.9
                  : jnc.fidelity ?? 0.9;
                return (
                  <tr key={i}>
                    <td>J{i + 1}</td>
                    <td>{pos + 1}</td>
                    <td><code>{jnc.overhang}</code></td>
                    <td>
                      <span className={`fidelity-badge ${fidelity >= 0.95 ? 'high' : fidelity >= 0.85 ? 'medium' : 'low'}`}>
                        {(fidelity * 100).toFixed(0)}%
                      </span>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>

          {domestication.fidelity && (
            <div className="assembly-fidelity">
              <strong>Overall Assembly Fidelity:</strong>{' '}
              <span className={`fidelity-value ${getFidelityValue(domestication.fidelity) >= 0.95 ? 'excellent' : getFidelityValue(domestication.fidelity) >= 0.85 ? 'good' : 'acceptable'}`}>
                {(getFidelityValue(domestication.fidelity) * 100).toFixed(1)}%
              </span>
              <span className="fidelity-source">
                (NEB experimental data)
              </span>
            </div>
          )}
        </div>
      )}

      {/* Hybrid strategy: show mutations too */}
      {strategy === 'hybrid' && domestication?.mutations && domestication.mutations.length > 0 && (
        <div className="strategy-details">
          <h4>Silent Mutations Applied</h4>
          <table className="mutations-table">
            <thead>
              <tr>
                <th>Position</th>
                <th>Original</th>
                <th>New</th>
                <th>Amino Acid</th>
                <th>Score</th>
              </tr>
            </thead>
            <tbody>
              {domestication.mutations.map((mut, i) => {
                const m = mut.mutation || mut;
                const position = m.sequencePosition ?? m.position ?? m.codonStart;
                const originalCodon = m.originalCodon || '';
                const newCodon = m.newCodon || '';
                const aminoAcid = m.aminoAcid || '';
                const score = m.score ?? mut.score;

                return (
                  <tr key={i}>
                    <td>{position !== undefined ? position + 1 : 'N/A'}</td>
                    <td><code>{originalCodon}</code></td>
                    <td><code>{newCodon}</code></td>
                    <td>{aminoAcid}</td>
                    <td>
                      <span className={`score-badge ${(score ?? 0) >= 80 ? 'good' : (score ?? 0) >= 60 ? 'ok' : 'poor'}`}>
                        {score?.toFixed?.(0) || 'N/A'}
                      </span>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      )}

      {/* Validation Status */}
      {result.validation && (
        <div className="validation-summary">
          <h4>Validation</h4>
          <div className={`validation-status ${result.validation?.isValid ? 'valid' : 'invalid'}`}>
            <span className="status-icon">{result.validation?.isValid ? '✓' : '⚠'}</span>
            <span>{result.validation?.isValid ? 'All checks passed' : 'Issues detected'}</span>
          </div>
          {result.validation?.warnings && result.validation.warnings.length > 0 && (
            <ul className="validation-warnings">
              {result.validation.warnings.map((w, i) => (
                <li key={i}>{typeof w === 'object' ? w.message : w}</li>
              ))}
            </ul>
          )}
        </div>
      )}

      <div className="step-actions">
        <button className="btn-secondary" onClick={onBack}>← Back</button>
        <button className="btn-primary" onClick={onContinue}>
          View Primers →
        </button>
      </div>
    </div>
  );
};

// Helper function to extract fidelity value
function getFidelityValue(fidelity: FidelityData | number): number {
  if (typeof fidelity === 'number') return fidelity;
  return fidelity.assembly ?? 0.9;
}

const PrimersStep: FC<PrimersStepProps> = ({ result, onContinue, onBack }) => {
  const { primers, primerSummary } = result;

  if (!primers || primers.length === 0) {
    return (
      <div className="workflow-step primers-step">
        <h2>Primer Design</h2>
        <div className="no-primers">
          <p>No primers needed for this strategy.</p>
        </div>
        <div className="step-actions">
          <button className="btn-secondary" onClick={onBack}>← Back</button>
          <button className="btn-primary" onClick={onContinue}>Continue →</button>
        </div>
      </div>
    );
  }

  return (
    <div className="workflow-step primers-step">
      <h2>Designed Primers</h2>

      {/* Primer Pairs */}
      <div className="primer-list">
        {primers.map((primerPair, i) => (
          <PrimerPairCard key={i} primerPair={primerPair} index={i} />
        ))}
      </div>

      {/* Ordering Summary */}
      {primerSummary && (
        <div className="ordering-summary">
          <h3>Ordering Summary</h3>
          <div className="summary-stats">
            <div className="stat">
              <span className="stat-value">{primerSummary.totalPrimers}</span>
              <span className="stat-label">Total Primers</span>
            </div>
            <div className="stat">
              <span className="stat-value">{primerSummary.totalLength} nt</span>
              <span className="stat-label">Total Length</span>
            </div>
            {(primerSummary.estimatedCost || primerSummary.totalCost) && (
              <div className="stat">
                <span className="stat-value">${(primerSummary.estimatedCost || primerSummary.totalCost || 0).toFixed(2)}</span>
                <span className="stat-label">Est. Cost</span>
              </div>
            )}
          </div>

          {primerSummary.orderList && primerSummary.orderList.length > 0 && (
            <div className="order-table-container">
              <h4>Ready-to-Order List</h4>
              <table className="order-table">
                <thead>
                  <tr>
                    <th>Name</th>
                    <th>Sequence (5' → 3')</th>
                    <th>Length</th>
                    <th>Tm</th>
                  </tr>
                </thead>
                <tbody>
                  {primerSummary.orderList.map((item, i) => (
                    <tr key={i}>
                      <td>{item.name}</td>
                      <td><code className="primer-seq">{item.sequence}</code></td>
                      <td>{item.length} nt</td>
                      <td>{item.tm?.toFixed(1) || 'N/A'}°C</td>
                    </tr>
                  ))}
                </tbody>
              </table>
              <button
                className="btn-copy"
                onClick={() => {
                  const text = primerSummary.orderList!
                    .map(p => `${p.name}\t${p.sequence}`)
                    .join('\n');
                  navigator.clipboard.writeText(text);
                }}
              >
                Copy for Ordering
              </button>
            </div>
          )}
        </div>
      )}

      <div className="step-actions">
        <button className="btn-secondary" onClick={onBack}>← Back</button>
        <button className="btn-primary" onClick={onContinue}>
          View Summary →
        </button>
      </div>
    </div>
  );
};

const PrimerPairCard: FC<PrimerPairCardProps> = ({ primerPair, index }) => {
  const [expanded, setExpanded] = useState<boolean>(false);

  // Handle silent_mutation_info type - no primers to show
  if (primerPair.type === 'silent_mutation_info') {
    return (
      <div className="primer-card info-card">
        <div className="primer-header">
          <span className="primer-name">Silent Mutation Applied</span>
          <span className="primer-badge info">No additional primers</span>
        </div>
        <div className="primer-details always-visible">
          <p className="info-text">
            Mutation applied directly to template sequence. Use standard primers for this region.
          </p>
        </div>
      </div>
    );
  }

  const { forward, reverse, type, junctionPosition, overhang, pairQuality, tmDifference } = primerPair;
  const pairName = primerPair.name || (type === 'junction' ? `Junction ${index + 1}` : `Fragment ${index + 1}`);

  return (
    <div className={`primer-card primer-pair ${expanded ? 'expanded' : ''}`}>
      <div className="primer-header" onClick={() => setExpanded(!expanded)}>
        <span className="primer-name">{pairName}</span>
        {type === 'junction' && junctionPosition !== undefined && (
          <span className="primer-badge junction">pos {junctionPosition + 1}</span>
        )}
        {overhang && (
          <span className="primer-badge overhang">{overhang}</span>
        )}
        {pairQuality !== undefined && (
          <span className={`quality-indicator ${pairQuality >= 80 ? 'excellent' : pairQuality >= 65 ? 'good' : 'acceptable'}`}>
            {pairQuality.toFixed(0)}
          </span>
        )}
        <button className="expand-btn">{expanded ? '−' : '+'}</button>
      </div>

      {expanded && (
        <div className="primer-details">
          {/* Forward Primer */}
          {forward && (
            <div className="individual-primer">
              <div className="primer-direction">
                <span className="direction-label">Forward (5' → 3')</span>
                <span className="primer-metrics">
                  {forward.length} nt | Tm: {forward.homologyTm?.toFixed(1) || 'N/A'}°C | GC: {forward.gcPercent || `${forward.gc?.toFixed(1)}%` || 'N/A'}
                </span>
              </div>
              <div className="sequence-display">
                <code>{forward.sequence}</code>
              </div>
              {forward.thermodynamics && (
                <div className="thermo-grid">
                  <div className="thermo-item">
                    <span className="label">Hairpin ΔG:</span>
                    <span className="value">{forward.thermodynamics.hairpinDG?.toFixed(1) || 'N/A'}</span>
                  </div>
                  <div className="thermo-item">
                    <span className="label">Homodimer ΔG:</span>
                    <span className="value">{forward.thermodynamics.homodimerDG?.toFixed(1) || 'N/A'}</span>
                  </div>
                  <div className="thermo-item">
                    <span className="label">3' Terminal ΔG:</span>
                    <span className="value">{forward.thermodynamics.terminal3DG?.toFixed(1) || 'N/A'}</span>
                  </div>
                </div>
              )}
              {forward.score !== undefined && (
                <div className="quality-score">
                  <span className={`score ${forward.qualityTier || (forward.score >= 80 ? 'excellent' : forward.score >= 65 ? 'good' : 'acceptable')}`}>
                    {forward.score.toFixed(0)} / 100
                  </span>
                  <span className="classification">{forward.qualityTier || 'N/A'}</span>
                </div>
              )}
              {forward.issues && forward.issues.length > 0 && (
                <div className="primer-issues">
                  {forward.issues.map((issue, i) => (
                    <span key={i} className="issue">{issue}</span>
                  ))}
                </div>
              )}
              {forward.warnings && forward.warnings.length > 0 && (
                <div className="primer-warnings">
                  {forward.warnings.map((warning, i) => (
                    <span key={i} className="warning">{warning}</span>
                  ))}
                </div>
              )}
            </div>
          )}

          {/* Reverse Primer */}
          {reverse && (
            <div className="individual-primer">
              <div className="primer-direction">
                <span className="direction-label">Reverse (5' → 3')</span>
                <span className="primer-metrics">
                  {reverse.length} nt | Tm: {reverse.homologyTm?.toFixed(1) || 'N/A'}°C | GC: {reverse.gcPercent || `${reverse.gc?.toFixed(1)}%` || 'N/A'}
                </span>
              </div>
              <div className="sequence-display">
                <code>{reverse.sequence}</code>
              </div>
              {reverse.thermodynamics && (
                <div className="thermo-grid">
                  <div className="thermo-item">
                    <span className="label">Hairpin ΔG:</span>
                    <span className="value">{reverse.thermodynamics.hairpinDG?.toFixed(1) || 'N/A'}</span>
                  </div>
                  <div className="thermo-item">
                    <span className="label">Homodimer ΔG:</span>
                    <span className="value">{reverse.thermodynamics.homodimerDG?.toFixed(1) || 'N/A'}</span>
                  </div>
                  <div className="thermo-item">
                    <span className="label">3' Terminal ΔG:</span>
                    <span className="value">{reverse.thermodynamics.terminal3DG?.toFixed(1) || 'N/A'}</span>
                  </div>
                </div>
              )}
              {reverse.score !== undefined && (
                <div className="quality-score">
                  <span className={`score ${reverse.qualityTier || (reverse.score >= 80 ? 'excellent' : reverse.score >= 65 ? 'good' : 'acceptable')}`}>
                    {reverse.score.toFixed(0)} / 100
                  </span>
                  <span className="classification">{reverse.qualityTier || 'N/A'}</span>
                </div>
              )}
              {reverse.issues && reverse.issues.length > 0 && (
                <div className="primer-issues">
                  {reverse.issues.map((issue, i) => (
                    <span key={i} className="issue">{issue}</span>
                  ))}
                </div>
              )}
              {reverse.warnings && reverse.warnings.length > 0 && (
                <div className="primer-warnings">
                  {reverse.warnings.map((warning, i) => (
                    <span key={i} className="warning">{warning}</span>
                  ))}
                </div>
              )}
            </div>
          )}

          {/* Pair metrics */}
          {(tmDifference !== undefined || primerPair.heterodimer) && (
            <div className="pair-metrics">
              {tmDifference !== undefined && (
                <div className="metric-item">
                  <span className="label">Tm Difference:</span>
                  <span className={`value ${tmDifference <= 3 ? 'good' : tmDifference <= 5 ? 'ok' : 'poor'}`}>
                    {tmDifference.toFixed(1)}°C
                  </span>
                </div>
              )}
              {primerPair.heterodimer && (
                <div className="metric-item">
                  <span className="label">Heterodimer ΔG:</span>
                  <span className={`value ${primerPair.heterodimer.risk === 'low' ? 'good' : primerPair.heterodimer.risk === 'moderate' ? 'ok' : 'poor'}`}>
                    {primerPair.heterodimer.dG?.toFixed(1) || 'N/A'} ({primerPair.heterodimer.risk})
                  </span>
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
};

const SummaryStep: FC<SummaryStepProps> = ({ result, onAccept, onBack }) => {
  const { workflowGuide, validation: _validation, domestication, strategy: _strategy } = result;

  return (
    <div className="workflow-step summary-step">
      <h2>Workflow Summary</h2>

      {/* Success Status */}
      {result.success ? (
        <div className="success-banner">
          <span className="success-icon">✓</span>
          <div>
            <strong>Domestication Complete</strong>
            <p>All internal sites have been addressed successfully.</p>
          </div>
        </div>
      ) : (
        <div className="warning-banner">
          <span className="warning-icon">⚠</span>
          <div>
            <strong>Workflow Incomplete</strong>
            <p>{result.message || 'Some issues need attention.'}</p>
          </div>
        </div>
      )}

      {/* Workflow Guide Steps */}
      {workflowGuide && workflowGuide.steps && (
        <div className="workflow-steps">
          <h3>Laboratory Workflow</h3>
          <ol className="lab-steps">
            {workflowGuide.steps.map((step, i) => (
              <li key={i} className="lab-step">
                <div className="step-number">{i + 1}</div>
                <div className="step-content">
                  <strong>{step.title}</strong>
                  <p>{step.description}</p>
                  {step.materials && (
                    <div className="step-materials">
                      <span className="materials-label">Materials:</span>
                      <span>{step.materials.join(', ')}</span>
                    </div>
                  )}
                </div>
              </li>
            ))}
          </ol>
        </div>
      )}

      {/* Key Metrics */}
      <div className="metrics-summary">
        <h3>Key Metrics</h3>
        <div className="metrics-grid">
          {domestication?.fidelity && (
            <div className="metric">
              <span className="metric-value">
                {(getFidelityValue(domestication.fidelity) * 100).toFixed(1)}%
              </span>
              <span className="metric-label">Assembly Fidelity</span>
            </div>
          )}
          {result.primers && (
            <div className="metric">
              <span className="metric-value">{result.primers.length}</span>
              <span className="metric-label">Primers Designed</span>
            </div>
          )}
          {domestication?.additionalFragments !== undefined && (
            <div className="metric">
              <span className="metric-value">{domestication.additionalFragments}</span>
              <span className="metric-label">Additional Fragments</span>
            </div>
          )}
        </div>
      </div>

      {/* Final Actions */}
      <div className="step-actions">
        <button className="btn-secondary" onClick={onBack}>← Back</button>
        <button
          className="btn-primary btn-accept"
          onClick={onAccept}
          disabled={!result.success}
        >
          ✓ Accept & Apply
        </button>
      </div>
    </div>
  );
};

// ============================================================================
// STYLES
// ============================================================================

export const DomesticationWorkflowStyles = `
.domestication-workflow-guide {
  max-width: 800px;
  margin: 0 auto;
  padding: 20px;
  font-family: system-ui, -apple-system, sans-serif;
}

/* Progress Bar */
.workflow-progress {
  display: flex;
  justify-content: space-between;
  margin-bottom: 24px;
  padding: 0 16px;
}

.progress-item {
  display: flex;
  flex-direction: column;
  align-items: center;
  flex: 1;
  position: relative;
}

.progress-item::after {
  content: '';
  position: absolute;
  top: 14px;
  left: 50%;
  width: 100%;
  height: 2px;
  background: #e5e7eb;
}

.progress-item:last-child::after {
  display: none;
}

.progress-item.done::after {
  background: #22c55e;
}

.progress-dot {
  width: 28px;
  height: 28px;
  border-radius: 50%;
  background: #e5e7eb;
  color: #6b7280;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 12px;
  font-weight: 600;
  z-index: 1;
}

.progress-item.done .progress-dot {
  background: #22c55e;
  color: white;
}

.progress-item.active .progress-dot {
  background: #3b82f6;
  color: white;
}

.progress-label {
  margin-top: 6px;
  font-size: 11px;
  color: #6b7280;
}

.progress-item.active .progress-label {
  color: #3b82f6;
  font-weight: 500;
}

/* Step Content */
.workflow-step {
  background: white;
  border-radius: 12px;
  padding: 24px;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.workflow-step h2 {
  margin: 0 0 8px;
  font-size: 1.4rem;
  color: #111827;
}

.step-intro {
  color: #6b7280;
  margin-bottom: 20px;
  line-height: 1.5;
}

/* Info Grid */
.info-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 12px;
  margin-bottom: 20px;
}

.info-card {
  display: flex;
  gap: 10px;
  padding: 12px;
  background: #f9fafb;
  border-radius: 8px;
}

.info-card.warning {
  background: #fef3c7;
}

.info-icon {
  font-size: 20px;
}

.info-title {
  font-weight: 600;
  font-size: 14px;
}

.info-detail {
  font-size: 12px;
  color: #6b7280;
}

/* Strategy Cards */
.strategy-cards {
  display: grid;
  grid-template-columns: repeat(2, 1fr);
  gap: 16px;
  margin-top: 16px;
}

.strategy-card {
  padding: 16px;
  border: 2px solid #e5e7eb;
  border-radius: 10px;
  position: relative;
}

.strategy-card.preferred {
  border-color: #22c55e;
  background: #f0fdf4;
}

.strategy-badge {
  position: absolute;
  top: -10px;
  right: 12px;
  background: #22c55e;
  color: white;
  font-size: 10px;
  padding: 2px 8px;
  border-radius: 10px;
  font-weight: 600;
}

.strategy-card h4 {
  margin: 0 0 8px;
  font-size: 1rem;
}

.strategy-card p {
  font-size: 13px;
  color: #6b7280;
  margin-bottom: 10px;
}

.strategy-card ul {
  margin: 0;
  padding-left: 18px;
  font-size: 12px;
}

.strategy-card li {
  margin-bottom: 4px;
}

/* Configuration */
.config-section {
  margin-bottom: 20px;
}

.config-section label {
  display: block;
  font-weight: 500;
  margin-bottom: 8px;
}

.config-checkbox {
  display: flex;
  align-items: center;
  gap: 8px;
  cursor: pointer;
}

.config-checkbox input {
  width: 18px;
  height: 18px;
}

.config-help {
  font-size: 12px;
  color: #6b7280;
  margin-top: 4px;
}

.frame-options {
  display: flex;
  gap: 8px;
}

.frame-btn {
  padding: 8px 16px;
  border: 2px solid #e5e7eb;
  border-radius: 6px;
  background: white;
  cursor: pointer;
  transition: all 0.2s;
}

.frame-btn:hover {
  border-color: #93c5fd;
}

.frame-btn.active {
  border-color: #3b82f6;
  background: #eff6ff;
}

/* Tables */
.mutations-table,
.junctions-table,
.order-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
  margin-top: 12px;
}

.mutations-table th,
.mutations-table td,
.junctions-table th,
.junctions-table td,
.order-table th,
.order-table td {
  padding: 8px 10px;
  text-align: left;
  border-bottom: 1px solid #e5e7eb;
}

.mutations-table th,
.junctions-table th,
.order-table th {
  background: #f9fafb;
  font-weight: 500;
}

code {
  background: #f3f4f6;
  padding: 2px 6px;
  border-radius: 4px;
  font-family: 'SF Mono', Monaco, monospace;
  font-size: 12px;
}

.primer-seq {
  word-break: break-all;
  font-size: 11px;
}

/* Badges */
.score-badge, .fidelity-badge {
  padding: 2px 8px;
  border-radius: 10px;
  font-size: 11px;
  font-weight: 500;
}

.score-badge.good, .fidelity-badge.high {
  background: #dcfce7;
  color: #166534;
}

.score-badge.ok, .fidelity-badge.medium {
  background: #fef3c7;
  color: #92400e;
}

.score-badge.poor, .fidelity-badge.low {
  background: #fef2f2;
  color: #dc2626;
}

/* Primer Cards */
.primer-list {
  display: flex;
  flex-direction: column;
  gap: 8px;
  margin-bottom: 20px;
}

.primer-card {
  border: 1px solid #e5e7eb;
  border-radius: 8px;
  overflow: hidden;
}

.primer-header {
  display: flex;
  align-items: center;
  gap: 12px;
  padding: 10px 12px;
  background: #f9fafb;
  cursor: pointer;
}

.primer-name {
  font-weight: 500;
  flex: 1;
}

.primer-length, .primer-tm {
  font-size: 12px;
  color: #6b7280;
}

.expand-btn {
  width: 24px;
  height: 24px;
  border: none;
  background: #e5e7eb;
  border-radius: 4px;
  cursor: pointer;
}

.primer-details {
  padding: 12px;
  background: white;
  border-top: 1px solid #e5e7eb;
}

.sequence-display {
  margin-bottom: 12px;
}

.sequence-display label {
  display: block;
  font-size: 11px;
  color: #6b7280;
  margin-bottom: 4px;
}

.sequence-display code {
  display: block;
  padding: 8px;
  word-break: break-all;
}

.thermo-grid {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 8px;
}

.thermo-item {
  font-size: 12px;
}

.thermo-item .label {
  color: #6b7280;
}

/* Primer Pair Card Enhancements */
.primer-card.info-card {
  background: #f0f9ff;
  border-color: #93c5fd;
}

.primer-card.info-card .primer-header {
  background: transparent;
  cursor: default;
}

.primer-details.always-visible {
  border-top: none;
  padding-top: 0;
}

.info-text {
  color: #6b7280;
  font-size: 13px;
  margin: 0;
}

.primer-badge {
  font-size: 10px;
  padding: 2px 8px;
  border-radius: 10px;
  font-weight: 500;
}

.primer-badge.info {
  background: #dbeafe;
  color: #1e40af;
}

.primer-badge.junction {
  background: #fef3c7;
  color: #92400e;
}

.primer-badge.overhang {
  background: #dcfce7;
  color: #166534;
  font-family: 'SF Mono', Monaco, monospace;
}

.quality-indicator {
  font-size: 11px;
  padding: 2px 8px;
  border-radius: 10px;
  font-weight: 600;
}

.quality-indicator.excellent {
  background: #dcfce7;
  color: #166534;
}

.quality-indicator.good {
  background: #dbeafe;
  color: #1e40af;
}

.quality-indicator.acceptable {
  background: #fef3c7;
  color: #92400e;
}

.individual-primer {
  padding: 12px;
  background: #f9fafb;
  border-radius: 8px;
  margin-bottom: 12px;
}

.individual-primer:last-child {
  margin-bottom: 0;
}

.primer-direction {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 8px;
}

.direction-label {
  font-weight: 600;
  font-size: 13px;
  color: #374151;
}

.primer-metrics {
  font-size: 11px;
  color: #6b7280;
}

.quality-score {
  display: flex;
  align-items: center;
  gap: 8px;
  margin-top: 8px;
}

.quality-score .score {
  padding: 4px 10px;
  border-radius: 6px;
  font-size: 12px;
  font-weight: 600;
}

.quality-score .score.excellent {
  background: #dcfce7;
  color: #166534;
}

.quality-score .score.good {
  background: #dbeafe;
  color: #1e40af;
}

.quality-score .score.acceptable {
  background: #fef3c7;
  color: #92400e;
}

.quality-score .score.poor {
  background: #fef2f2;
  color: #dc2626;
}

.quality-score .classification {
  font-size: 11px;
  color: #6b7280;
  text-transform: capitalize;
}

.pair-metrics {
  display: flex;
  gap: 16px;
  padding: 10px 12px;
  background: #f3f4f6;
  border-radius: 6px;
  margin-top: 12px;
}

.pair-metrics .metric-item {
  display: flex;
  align-items: center;
  gap: 6px;
  font-size: 12px;
}

.pair-metrics .metric-item .label {
  color: #6b7280;
}

.pair-metrics .metric-item .value {
  font-weight: 500;
}

.pair-metrics .metric-item .value.good {
  color: #166534;
}

.pair-metrics .metric-item .value.ok {
  color: #92400e;
}

.pair-metrics .metric-item .value.poor {
  color: #dc2626;
}

.primer-issues, .primer-warnings {
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  margin-top: 8px;
}

.primer-issues .issue {
  font-size: 11px;
  padding: 2px 8px;
  background: #fef2f2;
  color: #dc2626;
  border-radius: 4px;
}

.primer-warnings .warning {
  font-size: 11px;
  padding: 2px 8px;
  background: #fef3c7;
  color: #92400e;
  border-radius: 4px;
}

/* Summary */
.success-banner, .warning-banner {
  display: flex;
  align-items: center;
  gap: 12px;
  padding: 16px;
  border-radius: 8px;
  margin-bottom: 20px;
}

.success-banner {
  background: #dcfce7;
  color: #166534;
}

.warning-banner {
  background: #fef3c7;
  color: #92400e;
}

.success-icon, .warning-icon {
  font-size: 24px;
}

.lab-steps {
  list-style: none;
  padding: 0;
  margin: 0;
}

.lab-step {
  display: flex;
  gap: 12px;
  padding: 12px 0;
  border-bottom: 1px solid #e5e7eb;
}

.lab-step:last-child {
  border-bottom: none;
}

.lab-step .step-number {
  width: 28px;
  height: 28px;
  background: #3b82f6;
  color: white;
  border-radius: 50%;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 12px;
  font-weight: 600;
  flex-shrink: 0;
}

.lab-step .step-content strong {
  display: block;
  margin-bottom: 4px;
}

.lab-step .step-content p {
  font-size: 13px;
  color: #6b7280;
  margin: 0;
}

.step-materials {
  margin-top: 8px;
  font-size: 12px;
  color: #6b7280;
}

.metrics-grid {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 16px;
  margin-top: 12px;
}

.metric {
  text-align: center;
  padding: 16px;
  background: #f9fafb;
  border-radius: 8px;
}

.metric-value {
  display: block;
  font-size: 1.5rem;
  font-weight: 600;
  color: #111827;
}

.metric-label {
  font-size: 12px;
  color: #6b7280;
}

/* Actions */
.step-actions {
  display: flex;
  justify-content: flex-end;
  gap: 10px;
  margin-top: 24px;
  padding-top: 20px;
  border-top: 1px solid #e5e7eb;
}

.btn-primary, .btn-secondary {
  padding: 10px 20px;
  border-radius: 6px;
  font-weight: 500;
  cursor: pointer;
  transition: all 0.2s;
}

.btn-primary {
  background: #3b82f6;
  color: white;
  border: none;
}

.btn-primary:hover {
  background: #2563eb;
}

.btn-primary:disabled {
  background: #93c5fd;
  cursor: not-allowed;
}

.btn-primary.btn-accept {
  background: #22c55e;
}

.btn-primary.btn-accept:hover {
  background: #16a34a;
}

.btn-secondary {
  background: white;
  color: #374151;
  border: 1px solid #d1d5db;
}

.btn-secondary:hover {
  background: #f9fafb;
}

.btn-copy {
  margin-top: 12px;
  padding: 8px 16px;
  background: #f3f4f6;
  border: 1px solid #d1d5db;
  border-radius: 6px;
  font-size: 13px;
  cursor: pointer;
}

.btn-copy:hover {
  background: #e5e7eb;
}

/* Error */
.workflow-error, .analysis-error {
  background: #fef2f2;
  border: 1px solid #fecaca;
  color: #dc2626;
  padding: 10px 14px;
  border-radius: 8px;
  margin-bottom: 16px;
  display: flex;
  align-items: center;
  gap: 8px;
}

.workflow-error button {
  margin-left: auto;
  background: none;
  border: none;
  color: #dc2626;
  cursor: pointer;
  font-size: 16px;
}

/* No sites */
.workflow-guide.no-sites {
  text-align: center;
  padding: 40px;
}

.success-icon-large {
  width: 60px;
  height: 60px;
  background: #dcfce7;
  color: #22c55e;
  font-size: 30px;
  border-radius: 50%;
  display: flex;
  align-items: center;
  justify-content: center;
  margin: 0 auto 16px;
}

/* Validation */
.validation-summary {
  margin-top: 20px;
  padding: 12px;
  background: #f9fafb;
  border-radius: 8px;
}

.validation-status {
  display: flex;
  align-items: center;
  gap: 8px;
  padding: 8px 12px;
  border-radius: 6px;
  font-weight: 500;
}

.validation-status.valid {
  background: #dcfce7;
  color: #166534;
}

.validation-status.invalid {
  background: #fef2f2;
  color: #dc2626;
}

.validation-warnings {
  margin: 8px 0 0;
  padding-left: 20px;
  font-size: 13px;
  color: #92400e;
}

/* Analysis Preview */
.analysis-preview {
  background: #f9fafb;
  border-radius: 8px;
  padding: 16px;
  margin-top: 20px;
}

.analysis-preview h3 {
  margin: 0 0 12px;
  font-size: 14px;
}

.analysis-stats {
  display: flex;
  gap: 24px;
}

.analysis-stats .stat {
  text-align: center;
}

.analysis-stats .stat-value {
  display: block;
  font-size: 1.2rem;
  font-weight: 600;
}

.analysis-stats .stat-label {
  font-size: 11px;
  color: #6b7280;
}

/* Selected Strategy */
.selected-strategy {
  display: flex;
  align-items: center;
  gap: 16px;
  padding: 16px;
  background: #f0f9ff;
  border: 2px solid #3b82f6;
  border-radius: 10px;
  margin-bottom: 20px;
}

.strategy-icon {
  font-size: 32px;
}

.selected-strategy h3 {
  margin: 0;
  font-size: 1.1rem;
}

.selected-strategy p {
  margin: 4px 0 0;
  font-size: 13px;
  color: #6b7280;
}

/* Strategy Details */
.strategy-details {
  background: white;
  border: 1px solid #e5e7eb;
  border-radius: 8px;
  padding: 16px;
  margin-bottom: 16px;
}

.strategy-details h4 {
  margin: 0 0 12px;
  font-size: 14px;
}

.assembly-fidelity {
  margin-top: 16px;
  padding: 12px;
  background: #f0fdf4;
  border-radius: 6px;
  font-size: 13px;
}

.fidelity-value {
  font-weight: 600;
  margin: 0 4px;
}

.fidelity-value.excellent {
  color: #166534;
}

.fidelity-value.good {
  color: #3b82f6;
}

.fidelity-value.acceptable {
  color: #92400e;
}

.fidelity-source {
  font-size: 11px;
  color: #6b7280;
}
`;

export default DomesticationWorkflowGuide;
