/**
 * Enhanced Domestication Panel
 *
 * A comprehensive UI component for state-of-the-art Golden Gate domestication.
 * Provides:
 * - Strategy selection (Mutagenic Junction PRIMARY, Silent Mutation, Alt Enzyme)
 * - Reading frame detection with user confirmation
 * - Visual protein alignment preview
 * - Codon optimization mode selection
 * - Site-by-site mutation selection
 * - Pre-flight validation before committing changes
 *
 * STRATEGY HIERARCHY (user can override):
 * 1. Mutagenic Junction (PRIMARY) - One-pot Golden Gate compatible
 * 2. Silent Mutation - No extra fragments, but needs pre-modified template
 * 3. Alternative Enzyme - No mutations needed
 *
 * DESIGN PRINCIPLE: Never auto-commit changes without explicit user confirmation.
 */

import React, { useState, useMemo, useCallback } from 'react';

// Import the enhanced domestication functions
import {
  createDomesticationPlan,
  executeDomesticationPlan,
  ENHANCED_CONFIG,
  detectORFs,
  validateReadingFrame,
  translateSequence,
} from '../lib/repp/enhanced-domestication.js';
import { GOLDEN_GATE_ENZYMES, findInternalSites } from '../lib/repp/goldengate.js';
import { designAllMutagenicJunctions } from '../lib/repp/mutagenic-junction-domesticator.js';
import { recommendAlternativeEnzymes } from '../lib/repp/auto-domestication-optimizer.js';
import { designIntegratedPrimers, generatePrimerSummary } from '../lib/repp/domestication-primer-workflow.js';

// ============================================================================
// Type Definitions
// ============================================================================

type StepType = 'analyze' | 'strategy' | 'frame' | 'mutations' | 'preview' | 'primers' | 'complete';

interface InternalSite {
  position: number;
  sequence: string;
  orientation: string;
}

interface InternalSitesResult {
  hasSites: boolean;
  count: number;
  sites: InternalSite[];
}

interface ORF {
  frame: number;
  start: number;
  end: number;
  proteinLength: number;
  strand: 'forward' | 'reverse';
  avgCodonUsage: number;
}

interface ORFDetectionResult {
  hasOrfs: boolean;
  totalFound: number;
  orfs: ORF[];
  bestOrf: ORF | null;
}

interface ReadingFrameValidation {
  valid: boolean;
  internalStops: number[];
  warnings: string[];
}

interface MutationOption {
  type: 'silent_mutation' | 'mutagenic_junction';
  score: number;
  codonChange?: string;
  aminoAcid?: string;
  frequencyChange?: string;
  warnings?: string[];
  junction?: {
    junctionPosition: number;
  };
  overhang?: string;
}

interface SiteOption {
  site: InternalSite;
  options: MutationOption[];
  recommended: MutationOption | null;
}

interface MutationStep {
  step: string;
  siteOptions: SiteOption[];
}

interface DomesticationPlan {
  steps: MutationStep[];
  preview: {
    validation: {
      proteinPreserved: boolean;
      noRemainingSites: boolean;
    };
    original: {
      protein: string;
      proteinLength: number;
    };
    domesticated: {
      protein: string;
      proteinLength: number;
    };
    comparison: {
      proteinIdentical: boolean;
      differences: any[];
    };
    totalMutations: number;
    mutations: Array<{
      position: number;
      change: string;
      codon: string;
      aminoAcid: string;
    }>;
    totalJunctions: number;
    junctions: Array<{
      sitePosition?: number;
      position?: number;
      overhang?: string;
      description?: string;
    }>;
  };
}

interface ExecutionResult {
  success: boolean;
  strategy: string;
  domesticatedSequence?: string;
  message?: string;
  mutations: any[];
  junctions: any[];
  switchToEnzyme?: string;
}

interface PrimerComponents {
  fivePrime: string;
  mutation: string;
  threePrime: string;
}

interface PrimerDetail {
  sequence: string;
  length: number;
  homologyTm?: number;
  components?: PrimerComponents;
}

interface PrimerPairDesign {
  name?: string;
  type: string;
  forward?: PrimerDetail;
  reverse?: PrimerDetail;
  codonChange?: string;
  overhang?: string;
  pairQuality?: number;
  mutationIndex?: number;
  instructions?: string;
  pcrProtocol?: {
    annealingTemp: number;
    extensionTime: string;
  };
}

interface PrimerSummary {
  totalPrimers: number;
  totalLength: number;
  estimatedCost?: number;
  orderList?: Array<{
    name: string;
    sequence: string;
    length: number;
    tm?: number;
  }>;
}

interface PrimerDesignResult {
  primers: PrimerPairDesign[];
  primerSummary: PrimerSummary | null;
  strategy: string;
}

interface MutagenicAnalysisResult {
  success: boolean;
  junctions: any[];
  additionalFragments: number;
}

interface AlternativeEnzyme {
  enzyme: string;
  fullName?: string;
  recognition: string;
  overhangLength: number;
  internalSites: number;
  isCompatible: boolean;
  isCurrent: boolean;
}

interface FrameValidationOption {
  frame: number;
  validation: ReadingFrameValidation;
  orfs: ORF[];
  bestOrf: ORF | undefined;
  isRecommended: boolean;
  proteinPreview: string;
}

interface EnhancedDomesticationPanelProps {
  sequence: string;
  enzyme?: string;
  onDomesticationComplete?: (result: ExecutionResult & { primers: PrimerPairDesign[]; primerSummary: PrimerSummary | null }) => void;
  onCancel?: () => void;
}

// ============================================================================
// SVG ICONS
// ============================================================================

const Icons = {
  dna: (
    <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
      <path d="M4 2h2v2c0 1.44.68 2.61 1.88 3.78.86.83 2.01 1.63 3.12 2.38 1.11-.75 2.26-1.55 3.12-2.38C15.32 6.61 16 5.44 16 4V2h2v2c0 2.08-.96 3.71-2.22 5 1.26 1.29 2.22 2.92 2.22 5v2h-2v-2c0-1.44-.68-2.61-1.88-3.78-.86-.83-2.01-1.63-3.12-2.38-1.11.75-2.26 1.55-3.12 2.38C6.68 11.39 6 12.56 6 14v2H4v-2c0-2.08.96-3.71 2.22-5C4.96 7.71 4 6.08 4 4V2zm5 6c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5zm4 1c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5zm-4 4c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5zm4 1c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5z"/>
    </svg>
  ),
  scissors: (
    <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
      <path d="M9.64 7.64c.23-.5.36-1.05.36-1.64 0-2.21-1.79-4-4-4S2 3.79 2 6s1.79 4 4 4c.59 0 1.14-.13 1.64-.36L10 12l-2.36 2.36C7.14 14.13 6.59 14 6 14c-2.21 0-4 1.79-4 4s1.79 4 4 4 4-1.79 4-4c0-.59-.13-1.14-.36-1.64L12 14l7 7h3v-1L9.64 7.64zM6 8c-1.1 0-2-.89-2-2s.9-2 2-2 2 .89 2 2-.9 2-2 2zm0 12c-1.1 0-2-.89-2-2s.9-2 2-2 2 .89 2 2-.9 2-2 2zm6-7.5c-.28 0-.5-.22-.5-.5s.22-.5.5-.5.5.22.5.5-.22.5-.5.5zM19 3l-6 6 2 2 7-7V3h-3z"/>
    </svg>
  ),
  refresh: (
    <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
      <path d="M17.65 6.35C16.2 4.9 14.21 4 12 4c-4.42 0-7.99 3.58-7.99 8s3.57 8 7.99 8c3.73 0 6.84-2.55 7.73-6h-2.08c-.82 2.33-3.04 4-5.65 4-3.31 0-6-2.69-6-6s2.69-6 6-6c1.66 0 3.14.69 4.22 1.78L13 11h7V4l-2.35 2.35z"/>
    </svg>
  ),
  flask: (
    <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
      <path d="M19.8 18.4L14 10.67V6.5l1.35-1.69c.26-.33.03-.81-.39-.81H9.04c-.42 0-.65.48-.39.81L10 6.5v4.17L4.2 18.4c-.49.66-.02 1.6.8 1.6h14c.82 0 1.29-.94.8-1.6z"/>
    </svg>
  ),
  check: (
    <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
      <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
    </svg>
  ),
  warning: (
    <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
      <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
    </svg>
  ),
  close: (
    <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
      <path d="M19 6.41L17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12z"/>
    </svg>
  ),
  info: (
    <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
      <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-6h2v6zm0-8h-2V7h2v2z"/>
    </svg>
  ),
  copy: (
    <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
      <path d="M16 1H4c-1.1 0-2 .9-2 2v14h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z"/>
    </svg>
  ),
  checkCircle: (
    <svg viewBox="0 0 24 24" width="24" height="24" fill="currentColor">
      <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-2 15l-5-5 1.41-1.41L10 14.17l7.59-7.59L19 8l-9 9z"/>
    </svg>
  ),
  download: (
    <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
      <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
    </svg>
  ),
};

// ============================================================================
// SEQUENCE RULER COMPONENT - Visual sequence context
// ============================================================================

interface SequenceAnnotation {
  start: number;
  end: number;
  label: string;
  type: 'site' | 'mutation' | 'highlight';
  orientation?: string;
  details?: string;
}

interface SequenceRulerProps {
  sequence: string;
  annotations: SequenceAnnotation[];
  title?: string;
  showFullSequence?: boolean;
  contextSize?: number;
}

function SequenceRuler({
  sequence,
  annotations,
  title = 'Sequence Map',
  showFullSequence = false,
  contextSize = 25
}: SequenceRulerProps) {
  if (!sequence || !annotations.length) return null;

  // For full sequence view, show segments around each annotation
  const segments = useMemo(() => {
    if (showFullSequence) {
      // Show entire sequence with annotations
      return [{
        start: 0,
        end: sequence.length,
        annotations: annotations
      }];
    }

    // Create focused segments around each annotation
    return annotations.map((ann, idx) => {
      const start = Math.max(0, ann.start - contextSize);
      const end = Math.min(sequence.length, ann.end + contextSize);
      return {
        start,
        end,
        annotations: [ann],
        index: idx
      };
    });
  }, [sequence, annotations, showFullSequence, contextSize]);

  // Generate position markers for a segment
  const getPositionMarkers = (start: number, end: number) => {
    const markers: { pos: number; label: string }[] = [];
    const firstMarker = Math.ceil((start + 1) / 10) * 10;
    for (let pos = firstMarker; pos <= end; pos += 10) {
      markers.push({ pos: pos - 1, label: String(pos) });
    }
    return markers;
  };

  // Render a single segment
  const renderSegment = (segment: typeof segments[0], segIndex: number) => {
    const { start, end, annotations: segAnnotations } = segment;
    const segmentSeq = sequence.slice(start, end);
    const markers = getPositionMarkers(start, end);

    return (
      <div key={segIndex} className="sequence-segment">
        {/* Position ruler */}
        <div className="position-ruler">
          <span className="ruler-prefix">{start > 0 ? '...' : ''}</span>
          <span className="ruler-track">
            {segmentSeq.split('').map((_, i) => {
              const absPos = start + i;
              const marker = markers.find(m => m.pos === absPos);
              return (
                <span key={i} className="ruler-tick">
                  {marker ? <span className="ruler-marker">{marker.label}</span> : ''}
                </span>
              );
            })}
          </span>
          <span className="ruler-suffix">{end < sequence.length ? '...' : ''}</span>
        </div>

        {/* Tick marks */}
        <div className="tick-row">
          <span className="ruler-prefix"></span>
          <span className="tick-track">
            {segmentSeq.split('').map((_, i) => {
              const absPos = start + i;
              const isMajor = (absPos + 1) % 10 === 0;
              const isMinor = (absPos + 1) % 5 === 0;
              return (
                <span key={i} className={`tick ${isMajor ? 'major' : isMinor ? 'minor' : ''}`}>
                  {isMajor ? '|' : isMinor ? '·' : ' '}
                </span>
              );
            })}
          </span>
          <span className="ruler-suffix"></span>
        </div>

        {/* Sequence with highlights */}
        <div className="sequence-row">
          <span className="seq-prefix">{start > 0 ? '...' : '   '}</span>
          <code className="sequence-track">
            {segmentSeq.split('').map((nt, i) => {
              const absPos = start + i;
              // Check if this position is within any annotation
              const matchingAnn = segAnnotations.find(
                a => absPos >= a.start && absPos < a.end
              );
              const isHighlight = !!matchingAnn;
              const highlightClass = matchingAnn ? `highlight-${matchingAnn.type}` : '';

              return (
                <span
                  key={i}
                  className={`nt nt-${nt.toUpperCase()} ${isHighlight ? highlightClass : ''}`}
                >
                  {nt}
                </span>
              );
            })}
          </code>
          <span className="seq-suffix">{end < sequence.length ? '...' : '   '}</span>
        </div>

        {/* Annotation indicators */}
        <div className="annotation-row">
          <span className="ann-prefix"></span>
          <span className="annotation-track">
            {segmentSeq.split('').map((_, i) => {
              const absPos = start + i;
              const matchingAnn = segAnnotations.find(
                a => absPos >= a.start && absPos < a.end
              );
              if (!matchingAnn) return <span key={i} className="ann-space"> </span>;

              // Show caret at center of annotation
              const annCenter = Math.floor((matchingAnn.start + matchingAnn.end) / 2);
              if (absPos === annCenter) {
                return <span key={i} className="ann-caret">↑</span>;
              }
              if (absPos >= matchingAnn.start && absPos < matchingAnn.end) {
                return <span key={i} className="ann-underline">─</span>;
              }
              return <span key={i} className="ann-space"> </span>;
            })}
          </span>
          <span className="ann-suffix"></span>
        </div>

        {/* Annotation labels */}
        {segAnnotations.map((ann, annIdx) => (
          <div key={annIdx} className={`annotation-label type-${ann.type}`}>
            <span className="ann-icon">
              {ann.type === 'site' ? '◆' : ann.type === 'mutation' ? '→' : '●'}
            </span>
            <span className="ann-text">
              <strong>{ann.label}</strong>
              {ann.orientation && <span className="ann-orient"> ({ann.orientation})</span>}
              {ann.details && <span className="ann-details"> — {ann.details}</span>}
            </span>
            <span className="ann-pos">pos {ann.start + 1}</span>
          </div>
        ))}
      </div>
    );
  };

  return (
    <div className="sequence-ruler">
      <div className="ruler-header">
        <span className="ruler-icon">{Icons.dna}</span>
        <span className="ruler-title">{title}</span>
        <span className="ruler-legend">
          <span className="legend-item nt-A">A</span>
          <span className="legend-item nt-T">T</span>
          <span className="legend-item nt-G">G</span>
          <span className="legend-item nt-C">C</span>
        </span>
      </div>
      <div className="ruler-content">
        {segments.map((seg, idx) => renderSegment(seg, idx))}
      </div>
    </div>
  );
}

// ============================================================================
// MAIN COMPONENT
// ============================================================================

export function EnhancedDomesticationPanel({
  sequence,
  enzyme = 'BsaI',
  onDomesticationComplete,
  onCancel,
}: EnhancedDomesticationPanelProps) {
  // State
  const [step, setStep] = useState<StepType>('analyze');
  const [plan, setPlan] = useState<DomesticationPlan | null>(null);
  const [selectedStrategy, setSelectedStrategy] = useState<string>(ENHANCED_CONFIG.defaultStrategy);
  const [selectedFrame, setSelectedFrame] = useState<number | null>(null);
  const [codonMode, setCodonMode] = useState<string>(ENHANCED_CONFIG.codonModes.CONSERVATIVE);
  const [mutationSelections, setMutationSelections] = useState<Record<string, MutationOption>>({});
  const [error, setError] = useState<string | null>(null);
  const [isProcessing, setIsProcessing] = useState<boolean>(false);
  const [executionResult, setExecutionResult] = useState<ExecutionResult | null>(null);
  const [primerDesign, setPrimerDesign] = useState<PrimerDesignResult | null>(null);
  const [selectedAlternativeEnzyme, setSelectedAlternativeEnzyme] = useState<string | null>(null);

  // Detect internal enzyme sites in the sequence
  const internalSites = useMemo<InternalSitesResult | null>(() => {
    if (!sequence) return null;
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

  // ORF detection
  const orfDetection = useMemo<ORFDetectionResult | null>(() => {
    if (!sequence) return null;
    return detectORFs(sequence) as any;
  }, [sequence]);

  // Create domestication plan when we have all inputs
  const createPlan = useCallback(() => {
    if (!sequence || selectedFrame === null) return;

    setIsProcessing(true);
    setError(null);

    try {
      const newPlan = createDomesticationPlan(sequence, enzyme, {
        frame: selectedFrame,
        codonMode,
        preferredStrategy: selectedStrategy,
      }) as DomesticationPlan;

      setPlan(newPlan);

      // Initialize mutation selections with recommended options
      const initialSelections: Record<string, MutationOption> = {};
      const mutationStep = newPlan.steps.find((s: MutationStep) => s.step === 'MUTATION_OPTIONS');
      if (mutationStep) {
        for (const siteOption of mutationStep.siteOptions) {
          if (siteOption.recommended) {
            initialSelections[`site_${siteOption.site.position}`] = siteOption.recommended;
          }
        }
      }
      setMutationSelections(initialSelections);

      setStep('mutations');
    } catch (err) {
      setError((err as Error).message);
    } finally {
      setIsProcessing(false);
    }
  }, [sequence, enzyme, selectedFrame, codonMode, selectedStrategy]);

  // Execute the plan and design primers
  const executePlan = useCallback(() => {
    if (!plan) return;

    setIsProcessing(true);
    setError(null);

    try {
      const result = executeDomesticationPlan(plan as any, {
        frame: selectedFrame ?? undefined,
        strategy: selectedStrategy,
        ...mutationSelections,
      }) as ExecutionResult;

      if (result.success) {
        setExecutionResult(result);

        // Design primers for the domestication result
        try {
          const primerResult = designIntegratedPrimers(
            result.domesticatedSequence || sequence,
            result as any,  // strategyResult - contains strategy, mutations, junctions
            enzyme,
            {
              frame: selectedFrame ?? 0,
              organism: 'ecoli',
              targetTm: 60,
              existingOverhangs: [],
            }
          );

          // Extract primers array from the result object
          const primersArray = (primerResult?.primers || []) as any[];
          const summaryResult = (primerResult?.summary || generatePrimerSummary(primersArray, enzyme)) as any;

          // Ensure totalLength is present in summary
          const primerSummary = summaryResult ? {
            ...summaryResult,
            totalLength: summaryResult.totalLength ||
              (summaryResult.orderList?.reduce((sum: number, p: any) => sum + (p.length || 0), 0) || 0)
          } as PrimerSummary : null;

          setPrimerDesign({
            primers: primersArray as PrimerPairDesign[],
            primerSummary,
            strategy: result.strategy,
          });
        } catch (primerErr) {
          console.warn('Primer design warning:', (primerErr as Error).message);
          // Continue even if primer design fails - user can still apply domestication
          setPrimerDesign(null);
        }

        setStep('primers');
      } else {
        setError(result.message || 'Execution failed');
      }
    } catch (err) {
      setError((err as Error).message);
    } finally {
      setIsProcessing(false);
    }
  }, [plan, selectedFrame, selectedStrategy, mutationSelections, sequence, enzyme]);

  // Complete the workflow and return result
  const completeWorkflow = useCallback(() => {
    if (executionResult && onDomesticationComplete) {
      // Include primer design in the result
      const finalResult = {
        ...executionResult,
        primers: primerDesign?.primers || [],
        primerSummary: primerDesign?.primerSummary || null,
      };
      onDomesticationComplete(finalResult);
    }
    setStep('complete');
  }, [executionResult, primerDesign, onDomesticationComplete]);

  // Frame selection handler
  const handleFrameSelect = useCallback((frame: number) => {
    setSelectedFrame(frame);
  }, []);

  // Auto-select recommended frame when entering frame step
  const handleEnterFrameStep = useCallback(() => {
    // Auto-select the recommended frame if one exists and none is selected
    if (selectedFrame === null && orfDetection?.bestOrf?.frame !== undefined) {
      setSelectedFrame(orfDetection.bestOrf.frame);
    }
    setStep('frame');
  }, [selectedFrame, orfDetection]);

  // Handle strategy continue - different behavior based on strategy
  const handleStrategyContinue = useCallback(() => {
    if (selectedStrategy === ENHANCED_CONFIG.strategies.ALTERNATIVE_ENZYME) {
      // For alternative enzyme, complete workflow with enzyme switch
      if (selectedAlternativeEnzyme && onDomesticationComplete) {
        onDomesticationComplete({
          success: true,
          strategy: 'ALTERNATIVE_ENZYME',
          switchToEnzyme: selectedAlternativeEnzyme,
          message: `Switch to ${selectedAlternativeEnzyme} which has no internal recognition sites`,
          domesticatedSequence: sequence, // Sequence unchanged
          mutations: [],
          junctions: [],
          primers: [],
          primerSummary: null,
        });
      }
      setStep('complete');
    } else {
      // For other strategies, proceed to frame selection
      handleEnterFrameStep();
    }
  }, [selectedStrategy, selectedAlternativeEnzyme, sequence, onDomesticationComplete, handleEnterFrameStep]);

  // Mutation selection handler
  const handleMutationSelect = useCallback((siteKey: string, option: MutationOption) => {
    setMutationSelections(prev => ({
      ...prev,
      [siteKey]: option,
    }));
  }, []);

  // Render based on current step
  if (!sequence) {
    return <div className="domestication-panel empty">No sequence provided</div>;
  }

  if (!internalSites?.hasSites) {
    return (
      <div className="domestication-panel no-sites">
        <div className="success-icon">{Icons.checkCircle}</div>
        <h3>No Domestication Needed</h3>
        <p>Sequence has no internal {enzyme} sites and is already compatible.</p>
        <button className="btn-primary" onClick={onCancel}>Continue</button>
      </div>
    );
  }

  return (
    <div className="enhanced-domestication-panel">
      {/* Progress indicator */}
      <ProgressIndicator
        steps={['Analyze', 'Strategy', 'Frame', 'Mutations', 'Preview', 'Primers', 'Complete']}
        currentStep={step}
      />

      {/* Error display */}
      {error && (
        <div className="error-banner">
          <span className="error-icon">{Icons.warning}</span>
          <span>{error}</span>
          <button onClick={() => setError(null)}>{Icons.close}</button>
        </div>
      )}

      {/* Step content */}
      {step === 'analyze' && (
        <AnalyzeStep
          sequence={sequence}
          enzyme={enzyme}
          internalSites={internalSites}
          orfDetection={orfDetection}
          onContinue={() => setStep('strategy')}
        />
      )}

      {step === 'strategy' && (
        <StrategySelectionStep
          sequence={sequence}
          enzyme={enzyme}
          internalSites={internalSites}
          selectedStrategy={selectedStrategy}
          onStrategySelect={setSelectedStrategy}
          selectedAlternativeEnzyme={selectedAlternativeEnzyme}
          onAlternativeEnzymeSelect={setSelectedAlternativeEnzyme}
          onContinue={handleStrategyContinue}
          onBack={() => setStep('analyze')}
        />
      )}

      {step === 'frame' && (
        <FrameSelectionStep
          sequence={sequence}
          orfDetection={orfDetection}
          selectedFrame={selectedFrame}
          onFrameSelect={handleFrameSelect}
          codonMode={codonMode}
          onCodonModeChange={setCodonMode}
          onContinue={createPlan}
          onBack={() => setStep('strategy')}
          isProcessing={isProcessing}
        />
      )}

      {step === 'mutations' && plan && (
        <MutationSelectionStep
          sequence={sequence}
          enzyme={enzyme}
          plan={plan}
          selectedStrategy={selectedStrategy}
          mutationSelections={mutationSelections}
          onMutationSelect={handleMutationSelect}
          onContinue={() => setStep('preview')}
          onBack={() => setStep('frame')}
        />
      )}

      {step === 'preview' && plan && (
        <PreviewStep
          plan={plan}
          sequence={sequence}
          enzyme={enzyme}
          frame={selectedFrame}
          mutationSelections={mutationSelections}
          onApprove={executePlan}
          onBack={() => setStep('mutations')}
          isProcessing={isProcessing}
        />
      )}

      {step === 'primers' && (
        <PrimersStep
          executionResult={executionResult}
          primerDesign={primerDesign}
          enzyme={enzyme}
          onContinue={completeWorkflow}
          onBack={() => setStep('preview')}
        />
      )}

      {step === 'complete' && (
        <CompleteStep
          onClose={onCancel}
        />
      )}
    </div>
  );
}

// ============================================================================
// STEP COMPONENTS
// ============================================================================

interface ProgressIndicatorProps {
  steps: string[];
  currentStep: StepType;
}

function ProgressIndicator({ steps, currentStep }: ProgressIndicatorProps) {
  const stepIndex: Record<string, number> = {
    'analyze': 0,
    'strategy': 1,
    'frame': 2,
    'mutations': 3,
    'preview': 4,
    'primers': 5,
    'complete': 6,
  };

  const current = stepIndex[currentStep] ?? 0;

  return (
    <div className="progress-indicator">
      {steps.map((label, i) => (
        <div
          key={label}
          className={`progress-step ${i < current ? 'complete' : ''} ${i === current ? 'current' : ''}`}
        >
          <div className="step-number">{i < current ? Icons.check : i + 1}</div>
          <div className="step-label">{label}</div>
        </div>
      ))}
    </div>
  );
}

interface AnalyzeStepProps {
  sequence: string;
  enzyme: string;
  internalSites: InternalSitesResult;
  orfDetection: ORFDetectionResult | null;
  onContinue: () => void;
}

function AnalyzeStep({ sequence, enzyme, internalSites, orfDetection, onContinue }: AnalyzeStepProps) {
  return (
    <div className="step-content analyze-step">
      <h2>Sequence Analysis</h2>
      <p className="step-description">
        Your sequence contains internal {enzyme} recognition sites that need to be removed for Golden Gate assembly.
      </p>

      <div className="analysis-summary">
        <div className="summary-card warning">
          <div className="card-icon">{Icons.warning}</div>
          <div className="card-content">
            <div className="card-label">INTERNAL SITES</div>
            <div className="card-value">{internalSites.count}</div>
            <div className="card-detail">
              at position{internalSites.count !== 1 ? 's' : ''} {internalSites.sites.map(s => s.position + 1).join(', ')}
            </div>
          </div>
        </div>

        <div className="summary-card info">
          <div className="card-icon">{Icons.dna}</div>
          <div className="card-content">
            <div className="card-label">SEQUENCE LENGTH</div>
            <div className="card-value">{sequence.length.toLocaleString()}</div>
            <div className="card-detail">base pairs</div>
          </div>
        </div>

        <div className="summary-card success">
          <div className="card-icon">{orfDetection?.hasOrfs ? Icons.check : Icons.info}</div>
          <div className="card-content">
            <div className="card-label">ORF DETECTION</div>
            <div className="card-value">{orfDetection?.totalFound || 0}</div>
            <div className="card-detail">
              {orfDetection?.hasOrfs
                ? `best in frame ${orfDetection.bestOrf?.frame} (${orfDetection.bestOrf?.proteinLength} aa)`
                : 'no significant ORFs found'}
            </div>
          </div>
        </div>
      </div>

      {/* Sequence Map Visualization */}
      <SequenceRuler
        sequence={sequence}
        annotations={internalSites.sites.map((site, i) => ({
          start: site.position,
          end: site.position + (site.sequence?.length || 6),
          label: `${enzyme} site ${i + 1}`,
          type: 'site' as const,
          orientation: site.orientation,
          details: site.sequence
        }))}
        title={`Internal ${enzyme} Sites`}
        contextSize={20}
      />

      <div className="site-details">
        <h3>Site Details</h3>
        <table className="sites-table">
          <thead>
            <tr>
              <th>#</th>
              <th>Position</th>
              <th>Recognition Sequence</th>
              <th>Orientation</th>
            </tr>
          </thead>
          <tbody>
            {internalSites.sites.map((site, i) => {
              const siteSeq = site.sequence || '';
              return (
                <tr key={i}>
                  <td className="site-num">{i + 1}</td>
                  <td>{site.position + 1}</td>
                  <td><code className="site-seq">{siteSeq}</code></td>
                  <td className={`orient-${site.orientation?.toLowerCase()}`}>{site.orientation}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      <div className="step-actions">
        <button className="btn-primary" onClick={onContinue}>
          Continue to Strategy Selection →
        </button>
      </div>
    </div>
  );
}

/**
 * Strategy Selection Step
 *
 * Allows user to choose between domestication approaches:
 * 1. Mutagenic Junction (PRIMARY) - One-pot Golden Gate compatible
 * 2. Silent Mutation - No extra fragments
 * 3. Alternative Enzyme - No mutations needed
 */
interface StrategySelectionStepProps {
  sequence: string;
  enzyme: string;
  internalSites: InternalSitesResult;
  selectedStrategy: string;
  onStrategySelect: (strategy: string) => void;
  selectedAlternativeEnzyme: string | null;
  onAlternativeEnzymeSelect: (enzyme: string) => void;
  onContinue: () => void;
  onBack: () => void;
}

function StrategySelectionStep({
  sequence,
  enzyme,
  internalSites,
  selectedStrategy,
  onStrategySelect,
  selectedAlternativeEnzyme,
  onAlternativeEnzymeSelect,
  onContinue,
  onBack,
}: StrategySelectionStepProps) {
  // Analyze available strategies
  const mutagenicAnalysis = useMemo<MutagenicAnalysisResult>(() => {
    try {
      return designAllMutagenicJunctions(sequence, enzyme, { frame: 0 }) as MutagenicAnalysisResult;
    } catch (e) {
      return { success: false, junctions: [], additionalFragments: internalSites?.count || 0 };
    }
  }, [sequence, enzyme, internalSites]);

  const alternativeEnzymes = useMemo<AlternativeEnzyme[]>(() => {
    try {
      const alternatives = recommendAlternativeEnzymes(sequence, enzyme) as AlternativeEnzyme[];
      return alternatives.filter((e: AlternativeEnzyme) => e.isCompatible && !e.isCurrent);
    } catch (e) {
      return [];
    }
  }, [sequence, enzyme]);

  const strategies = [
    {
      id: ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION,
      name: 'Mutagenic Junction',
      tagline: 'RECOMMENDED for One-Pot Golden Gate',
      icon: Icons.scissors,
      description: 'Primers introduce silent mutations during PCR. The mutation becomes permanent after amplification.',
      benefits: [
        'One-pot Golden Gate compatible (no re-cutting)',
        'No pre-modified template needed',
        'High-fidelity overhang selection',
        'Mutations are permanent after PCR',
      ],
      tradeoffs: [
        `Adds ${mutagenicAnalysis?.additionalFragments || internalSites?.count || 1} additional fragment(s)`,
        'Requires additional primers',
      ],
      feasible: mutagenicAnalysis?.success || (mutagenicAnalysis?.junctions?.length > 0),
      details: {
        fragments: (mutagenicAnalysis?.additionalFragments || internalSites?.count || 1) + 1,
        junctionsDesigned: mutagenicAnalysis?.junctions?.length || 0,
      },
    },
    {
      id: ENHANCED_CONFIG.strategies.SILENT_MUTATION,
      name: 'Silent Mutation',
      tagline: 'Simpler Assembly (Fewer Fragments)',
      icon: Icons.dna,
      description: 'Directly modify the template sequence with synonymous codon changes before assembly.',
      benefits: [
        'No additional fragments needed',
        'Simpler assembly workflow',
        'Can use existing primer designs',
      ],
      tradeoffs: [
        'Requires pre-modified template (synthesis or site-directed mutagenesis)',
        'NOT one-pot compatible if sites re-cut during reaction',
        'Template must be prepared separately',
      ],
      feasible: true,
      details: {
        fragments: 1,
        requiresModification: true,
      },
    },
    {
      id: ENHANCED_CONFIG.strategies.ALTERNATIVE_ENZYME,
      name: 'Alternative Enzyme',
      tagline: 'No Mutations Required',
      icon: Icons.refresh,
      description: 'Switch to a different Type IIS enzyme that has no internal recognition sites in your sequence.',
      benefits: [
        'Preserves original sequence exactly',
        'No mutations needed',
        'One-pot compatible',
      ],
      tradeoffs: alternativeEnzymes.length === 0
        ? ['No compatible alternatives available']
        : [
            'May require different overhang designs',
            'Different enzyme characteristics',
          ],
      feasible: alternativeEnzymes.length > 0,
      details: {
        availableEnzymes: alternativeEnzymes,
      },
    },
  ];

  // Get all enzymes for comparison (including non-compatible ones)
  const allAlternativeEnzymes = useMemo<AlternativeEnzyme[]>(() => {
    try {
      return recommendAlternativeEnzymes(sequence, enzyme) as AlternativeEnzyme[];
    } catch (e) {
      return [];
    }
  }, [sequence, enzyme]);

  // Determine if continue is allowed for alternative enzyme strategy
  const canContinueAlternativeEnzyme = selectedStrategy === ENHANCED_CONFIG.strategies.ALTERNATIVE_ENZYME
    && selectedAlternativeEnzyme
    && alternativeEnzymes.some((e: AlternativeEnzyme) => e.enzyme === selectedAlternativeEnzyme);

  return (
    <div className="step-content strategy-step">
      <h2>Choose Domestication Strategy</h2>
      <p className="step-description">
        Select how to handle the {internalSites?.count || 0} internal {enzyme} site(s).
        <strong> Mutagenic Junction</strong> is recommended for one-pot Golden Gate assemblies.
      </p>

      <div className="strategy-options">
        {strategies.map(strategy => {
          // ONE-POT is only for Mutagenic Junction and Alternative Enzyme (NOT Silent Mutation)
          const isOnePot = strategy.id === ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION ||
                           strategy.id === ENHANCED_CONFIG.strategies.ALTERNATIVE_ENZYME;
          const isSelected = selectedStrategy === strategy.id;

          return (
            <div
              key={strategy.id}
              className={`strategy-card ${isSelected ? 'selected' : ''} ${!strategy.feasible ? 'disabled' : ''} ${strategy.id === ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION ? 'recommended' : ''}`}
              onClick={() => strategy.feasible && onStrategySelect(strategy.id)}
            >
              <div className="strategy-content">
                <div className="strategy-header">
                  <span className="strategy-title">{strategy.name}</span>
                  <div className="strategy-badges">
                    {strategy.id === ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION && (
                      <span className="strategy-badge recommended">Recommended</span>
                    )}
                    {isOnePot && strategy.feasible && (
                      <span className="strategy-badge one-pot">One-pot</span>
                    )}
                  </div>
                </div>

                <p className="strategy-description">{strategy.description}</p>

                <div className="strategy-meta">
                  <div className="meta-benefits">
                    {strategy.benefits.slice(0, 2).map((b, i) => (
                      <span key={i} className="meta-item benefit">✓ {b}</span>
                    ))}
                  </div>
                  <div className="meta-tradeoffs">
                    {strategy.tradeoffs.slice(0, 1).map((t, i) => (
                      <span key={i} className="meta-item tradeoff">⚠ {t}</span>
                    ))}
                  </div>
                </div>
              </div>

              <div className="strategy-icon-container">
                <div className="strategy-icon">{strategy.icon}</div>
              </div>
            </div>
          );
        })}
      </div>

      {/* Enzyme Selection Section - shown when Alternative Enzyme is selected */}
      {selectedStrategy === ENHANCED_CONFIG.strategies.ALTERNATIVE_ENZYME && (
        <div className="enzyme-selection-section">
          <h3>Select Alternative Enzyme</h3>
          <p className="section-description">
            Choose an enzyme with no internal recognition sites in your sequence.
            Enzymes are sorted by compatibility and ligation data availability.
          </p>

          <div className="enzyme-comparison-table">
            <table className="enzyme-table">
              <thead>
                <tr>
                  <th>Enzyme</th>
                  <th>Recognition</th>
                  <th>Overhang</th>
                  <th>Internal Sites</th>
                  <th>Status</th>
                </tr>
              </thead>
              <tbody>
                {/* Current enzyme row */}
                <tr className="current-enzyme">
                  <td>
                    <strong>{enzyme}</strong>
                    <span className="current-badge">Current</span>
                  </td>
                  <td><code>{(GOLDEN_GATE_ENZYMES as any)[enzyme]?.recognition || 'N/A'}</code></td>
                  <td>{(GOLDEN_GATE_ENZYMES as any)[enzyme]?.overhangLength || 4} nt</td>
                  <td className="site-count warning">{internalSites?.count || 0} sites</td>
                  <td><span className="status-badge incompatible">Needs domestication</span></td>
                </tr>

                {/* Alternative enzymes */}
                {allAlternativeEnzymes.filter((e: AlternativeEnzyme) => !e.isCurrent).map((alt) => (
                  <tr
                    key={alt.enzyme}
                    className={`enzyme-option ${alt.isCompatible ? 'compatible' : 'incompatible'} ${selectedAlternativeEnzyme === alt.enzyme ? 'selected' : ''}`}
                    onClick={() => alt.isCompatible && onAlternativeEnzymeSelect(alt.enzyme)}
                  >
                    <td>
                      <strong>{alt.enzyme}</strong>
                      {alt.fullName && <span className="enzyme-fullname">{alt.fullName}</span>}
                    </td>
                    <td><code>{alt.recognition}</code></td>
                    <td>{alt.overhangLength} nt</td>
                    <td className={`site-count ${alt.internalSites === 0 ? 'success' : 'warning'}`}>
                      {alt.internalSites === 0 ? 'None' : `${alt.internalSites} sites`}
                    </td>
                    <td>
                      {alt.isCompatible ? (
                        <span className="status-badge compatible">
                          {selectedAlternativeEnzyme === alt.enzyme ? '✓ Selected' : 'Compatible'}
                        </span>
                      ) : (
                        <span className="status-badge incompatible">Has internal sites</span>
                      )}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          {alternativeEnzymes.length === 0 && (
            <div className="no-alternatives-warning">
              <span className="warning-icon">{Icons.warning}</span>
              <span>No compatible alternative enzymes found. All available Type IIS enzymes have internal recognition sites in your sequence.</span>
            </div>
          )}
        </div>
      )}

      <div className="step-actions">
        <button className="btn-secondary" onClick={onBack}>← Back</button>
        <button
          className="btn-primary"
          onClick={onContinue}
          disabled={
            !selectedStrategy ||
            (selectedStrategy === ENHANCED_CONFIG.strategies.ALTERNATIVE_ENZYME && !canContinueAlternativeEnzyme)
          }
        >
          {selectedStrategy === ENHANCED_CONFIG.strategies.ALTERNATIVE_ENZYME
            ? `Switch to ${selectedAlternativeEnzyme || '...'} →`
            : 'Continue to Frame Selection →'}
        </button>
      </div>
    </div>
  );
}

interface FrameSelectionStepProps {
  sequence: string;
  orfDetection: ORFDetectionResult | null;
  selectedFrame: number | null;
  onFrameSelect: (frame: number) => void;
  codonMode: string;
  onCodonModeChange: (mode: string) => void;
  onContinue: () => void;
  onBack: () => void;
  isProcessing: boolean;
}

function FrameSelectionStep({
  sequence,
  orfDetection,
  selectedFrame,
  onFrameSelect,
  codonMode,
  onCodonModeChange,
  onContinue,
  onBack,
  isProcessing,
}: FrameSelectionStepProps) {
  const frameOptions: FrameValidationOption[] = [0, 1, 2].map(frame => {
    const validationResult = validateReadingFrame(sequence, frame) as any;
    const validation: ReadingFrameValidation = {
      valid: validationResult.isValid || false,
      internalStops: validationResult.internalStops || [],
      warnings: validationResult.validations?.filter((v: any) => !v.passed).map((v: any) => v.message) || [],
    };
    const orfs = orfDetection?.orfs.filter((o: ORF) => o.frame === frame && o.strand === 'forward') || [];
    const bestOrf = orfs[0];

    const translationResult = translateSequence(sequence, frame) as any;
    return {
      frame,
      validation,
      orfs,
      bestOrf,
      isRecommended: orfDetection?.bestOrf?.frame === frame,
      proteinPreview: translationResult.protein.slice(0, 30),
    };
  });

  return (
    <div className="step-content frame-step">
      <h2>Reading Frame Selection</h2>
      <p className="step-description">
        <strong>CRITICAL:</strong> Select the correct reading frame to ensure mutations remain silent (synonymous).
        Wrong frame selection will alter the protein sequence!
      </p>

      <div className="frame-options">
        {frameOptions.map(opt => {
          // Determine severity based on stop codons
          const stopCount = opt.validation?.internalStops?.length || 0;
          const severityClass = stopCount === 0 ? 'severity-ok' : stopCount <= 5 ? 'severity-warning' : 'severity-danger';

          return (
            <div
              key={opt.frame}
              className={`frame-option ${selectedFrame === opt.frame ? 'selected' : ''} ${severityClass}`}
              onClick={() => onFrameSelect(opt.frame)}
            >
              <div className="frame-header">
                <div className="frame-number">Frame {opt.frame}</div>
                <div className="frame-badges">
                  {opt.isRecommended && <span className="frame-badge recommended">Best</span>}
                  {selectedFrame === opt.frame && <span className="frame-badge selected">✓</span>}
                </div>
              </div>

              <div className="frame-stats">
                {opt.bestOrf ? (
                  <>
                    <div className="frame-stat">
                      <span className="label">ORF Length</span>
                      <span className="value">{opt.bestOrf.proteinLength} aa</span>
                    </div>
                    <div className="frame-stat">
                      <span className="label">Codon Usage</span>
                      <span className="value">
                        {typeof opt.bestOrf.avgCodonUsage === 'number'
                          ? `${opt.bestOrf.avgCodonUsage.toFixed(0)}/1000`
                          : 'N/A'}
                      </span>
                    </div>
                  </>
                ) : (
                  <div className="frame-stat">
                    <span className="label">No significant ORF detected</span>
                  </div>
                )}
              </div>

              <div className="protein-preview">
                <code>{opt.proteinPreview}...</code>
              </div>

              {/* Status indicator - always show for consistent card heights */}
              <div className={`frame-status ${stopCount === 0 ? 'success' : stopCount <= 5 ? 'warning' : 'danger'}`}>
                {stopCount === 0 ? (
                  <>{Icons.check} No internal stops — Recommended</>
                ) : (
                  <>{Icons.warning} {stopCount} internal stop{stopCount !== 1 ? 's' : ''} — {stopCount > 5 ? 'Avoid this frame' : 'Use with caution'}</>
                )}
              </div>
            </div>
          );
        })}
      </div>

      <div className="codon-mode-section">
        <h3>Codon Optimization Mode</h3>
        <div className="mode-options">
          <label className={codonMode === ENHANCED_CONFIG.codonModes.CONSERVATIVE ? 'selected' : ''}>
            <input
              type="radio"
              name="codonMode"
              value={ENHANCED_CONFIG.codonModes.CONSERVATIVE}
              checked={codonMode === ENHANCED_CONFIG.codonModes.CONSERVATIVE}
              onChange={(e: React.ChangeEvent<HTMLInputElement>) => onCodonModeChange(e.target.value)}
            />
            <div className="mode-content">
              <div className="mode-title">Conservative</div>
              <div className="mode-desc">Prefer codons with similar usage frequency</div>
            </div>
          </label>

          <label className={codonMode === ENHANCED_CONFIG.codonModes.OPTIMIZED ? 'selected' : ''}>
            <input
              type="radio"
              name="codonMode"
              value={ENHANCED_CONFIG.codonModes.OPTIMIZED}
              checked={codonMode === ENHANCED_CONFIG.codonModes.OPTIMIZED}
              onChange={(e: React.ChangeEvent<HTMLInputElement>) => onCodonModeChange(e.target.value)}
            />
            <div className="mode-content">
              <div className="mode-title">Optimized</div>
              <div className="mode-desc">Prefer codons optimal for target organism</div>
            </div>
          </label>

          <label className={codonMode === ENHANCED_CONFIG.codonModes.CUSTOM ? 'selected' : ''}>
            <input
              type="radio"
              name="codonMode"
              value={ENHANCED_CONFIG.codonModes.CUSTOM}
              checked={codonMode === ENHANCED_CONFIG.codonModes.CUSTOM}
              onChange={(e: React.ChangeEvent<HTMLInputElement>) => onCodonModeChange(e.target.value)}
            />
            <div className="mode-content">
              <div className="mode-title">Custom</div>
              <div className="mode-desc">Manually select mutations for each site</div>
            </div>
          </label>
        </div>
      </div>

      <div className="step-actions">
        <button className="btn-secondary" onClick={onBack}>← Back</button>
        <button
          className="btn-primary"
          onClick={onContinue}
          disabled={selectedFrame === null || isProcessing}
        >
          {isProcessing ? 'Processing...' : 'Generate Mutation Options →'}
        </button>
      </div>
    </div>
  );
}

interface MutationSelectionStepProps {
  sequence: string;
  enzyme: string;
  plan: DomesticationPlan;
  selectedStrategy: string;
  mutationSelections: Record<string, MutationOption>;
  onMutationSelect: (siteKey: string, option: MutationOption) => void;
  onContinue: () => void;
  onBack: () => void;
}

function MutationSelectionStep({
  sequence,
  enzyme,
  plan,
  selectedStrategy,
  mutationSelections,
  onMutationSelect,
  onContinue,
  onBack,
}: MutationSelectionStepProps) {
  const mutationStep = plan.steps.find((s: MutationStep) => s.step === 'MUTATION_OPTIONS');

  // Build annotations for the sequence ruler
  const rulerAnnotations: SequenceAnnotation[] = useMemo(() => {
    if (!mutationStep) return [];
    return mutationStep.siteOptions.map((siteOption: SiteOption, i: number) => {
      const siteKey = `site_${siteOption.site.position}`;
      const selected = mutationSelections[siteKey];
      const siteSeq = siteOption.site.sequence || '';

      return {
        start: siteOption.site.position,
        end: siteOption.site.position + siteSeq.length,
        label: `Site ${i + 1}`,
        type: selected?.type === 'silent_mutation' ? 'mutation' as const : 'site' as const,
        orientation: siteOption.site.orientation,
        details: selected?.type === 'silent_mutation'
          ? `${selected.codonChange} (silent)`
          : selected?.type === 'mutagenic_junction'
          ? 'junction'
          : siteSeq
      };
    });
  }, [mutationStep, mutationSelections]);

  // Determine which option type matches the selected strategy
  const isMutagenicStrategy = selectedStrategy === ENHANCED_CONFIG.strategies.MUTAGENIC_JUNCTION;
  const preferredType = isMutagenicStrategy ? 'mutagenic_junction' : 'silent_mutation';
  const strategyName = isMutagenicStrategy ? 'Mutagenic Junction' : 'Silent Mutation';

  // Auto-select appropriate options based on strategy when component mounts or strategy changes
  React.useEffect(() => {
    if (!mutationStep) return;

    mutationStep.siteOptions.forEach((siteOption: SiteOption) => {
      const siteKey = `site_${siteOption.site.position}`;
      const current = mutationSelections[siteKey];
      const filteredOptions = siteOption.options.filter((o: MutationOption) => o.type === preferredType);
      const displayOptions = filteredOptions.length > 0 ? filteredOptions : siteOption.options;

      // If current selection is not in display options, auto-select first matching option
      if (!current || !displayOptions.includes(current)) {
        onMutationSelect(siteKey, displayOptions[0]);
      }
    });
  }, [mutationStep, preferredType, mutationSelections, onMutationSelect]);

  if (!mutationStep) return null;

  return (
    <div className="step-content mutations-step">
      <h2>Site Options ({strategyName} Strategy)</h2>
      <p className="step-description">
        {isMutagenicStrategy
          ? 'Each site will be handled by a mutagenic junction (one-pot Golden Gate compatible).'
          : 'Each site will be handled by a silent mutation applied to the template sequence.'}
      </p>

      {/* Sequence overview showing all sites */}
      <SequenceRuler
        sequence={sequence}
        annotations={rulerAnnotations}
        title={`All ${enzyme} Sites to Modify`}
        contextSize={15}
      />

      <div className="site-mutations">
        {mutationStep.siteOptions.map((siteOption: SiteOption, siteIndex: number) => {
          const siteKey = `site_${siteOption.site.position}`;
          const selected = mutationSelections[siteKey];

          // Filter options based on strategy - only show relevant options
          const filteredOptions = siteOption.options.filter((o: MutationOption) => o.type === preferredType);
          // Fall back to all options if no matching options
          const displayOptions = filteredOptions.length > 0 ? filteredOptions : siteOption.options;

          // Determine effective selection - use first option if current selection isn't in displayed options
          const effectiveSelected = displayOptions.includes(selected) ? selected : displayOptions[0];
          const isSelected = (option: MutationOption) => option === effectiveSelected;

          // Create color-coded sequence display
          const siteSeq = siteOption.site.sequence || '';
          const coloredSeq = siteSeq.split('').map((nt, i) => (
            <span key={i} className={`nt-${nt}`}>{nt}</span>
          ));

          return (
            <div key={siteIndex} className="site-mutation-card">
              <div className="site-header">
                <span className="site-position">Site at position {siteOption.site.position + 1}</span>
                <code className="site-sequence">{siteOption.site.sequence}</code>
              </div>

              {/* Sequence minimap showing the site in context */}
              <div className="sequence-minimap">
                <div className="minimap-label">Recognition Site</div>
                <div className="minimap-sequence">
                  <span className="site-highlight">{coloredSeq}</span>
                  {effectiveSelected?.type === 'mutagenic_junction' && effectiveSelected.junction && (
                    <span className="junction-marker"> ← junction at {effectiveSelected.junction.junctionPosition}</span>
                  )}
                </div>
              </div>

              <div className="mutation-options">
                {displayOptions.map((option: MutationOption, optIndex: number) => {
                  const isRecommended = option === siteOption.recommended;
                  const isOptionSelected = isSelected(option);

                  // Format warnings nicely
                  const formatWarning = (w: string) => {
                    if (w === 'RARE_CODON' || w.includes('rare') || w.includes('Rare')) {
                      return '⚠️ Low frequency';
                    }
                    return `⚠️ ${w}`;
                  };

                  return (
                    <div
                      key={optIndex}
                      className={`mutation-option ${isOptionSelected ? 'selected' : ''} ${isRecommended ? 'recommended' : ''}`}
                      onClick={() => onMutationSelect(siteKey, option)}
                    >
                      <div className="option-header">
                        <span className="option-type">
                          {option.type === 'silent_mutation' ? 'Silent Mutation' : 'Mutagenic Junction'}
                        </span>
                        <div className="option-badges">
                          {isRecommended && (
                            <span className="option-badge recommended">Recommended</span>
                          )}
                          <span className={`option-badge ${option.type === 'silent_mutation' ? 'silent' : 'junction'}`}>
                            {option.type === 'silent_mutation' ? 'No fragments' : '+1 fragment'}
                          </span>
                        </div>
                      </div>

                      <div className="option-details">
                        {option.type === 'silent_mutation' ? (
                          <>
                            <div className="detail">
                              <span className="label">Codon:</span>
                              <code>{option.codonChange}</code>
                              <span className="aa">({option.aminoAcid})</span>
                            </div>
                            <div className="detail">
                              <span className="label">Frequency:</span>
                              <span>{option.frequencyChange}</span>
                            </div>
                            {option.warnings && option.warnings.length > 0 && (
                              <div className="warnings">
                                {option.warnings.map((w: string, i: number) => (
                                  <span key={i} className="warning-tag">{formatWarning(w)}</span>
                                ))}
                              </div>
                            )}
                          </>
                        ) : (
                          <>
                            <div className="detail">
                              <span className="label">Junction at:</span>
                              <span>position {option.junction?.junctionPosition}</span>
                            </div>
                            <div className="detail">
                              <span className="label">Overhang:</span>
                              <code>{option.overhang}</code>
                            </div>
                          </>
                        )}
                      </div>
                    </div>
                  );
                })}
              </div>
            </div>
          );
        })}
      </div>

      <div className="step-actions">
        <button className="btn-secondary" onClick={onBack}>← Back</button>
        <button
          className="btn-primary"
          onClick={onContinue}
          disabled={Object.keys(mutationSelections).length !== mutationStep.siteOptions.length}
        >
          Preview Changes →
        </button>
      </div>
    </div>
  );
}

interface PreviewStepProps {
  plan: DomesticationPlan;
  sequence: string;
  enzyme: string;
  frame: number | null;
  mutationSelections: Record<string, MutationOption>;
  onApprove: () => void;
  onBack: () => void;
  isProcessing: boolean;
}

function PreviewStep({
  plan,
  sequence,
  enzyme,
  mutationSelections,
  onApprove,
  onBack,
  isProcessing,
}: PreviewStepProps) {
  const preview = plan.preview;
  const mutationStep = plan.steps.find((s: MutationStep) => s.step === 'MUTATION_OPTIONS');

  // Build annotations showing what will be modified
  const previewAnnotations: SequenceAnnotation[] = useMemo(() => {
    if (!mutationStep) return [];
    return mutationStep.siteOptions.map((siteOption: SiteOption, i: number) => {
      const siteKey = `site_${siteOption.site.position}`;
      const selected = mutationSelections[siteKey];
      const siteSeq = siteOption.site.sequence || '';

      return {
        start: siteOption.site.position,
        end: siteOption.site.position + siteSeq.length,
        label: selected?.type === 'silent_mutation'
          ? `${selected.codonChange}`
          : selected?.type === 'mutagenic_junction'
          ? `Junction ${i + 1}`
          : `Site ${i + 1}`,
        type: 'mutation' as const,
        details: selected?.type === 'silent_mutation' ? 'silent mutation' : 'mutagenic junction'
      };
    });
  }, [mutationStep, mutationSelections]);

  if (!preview) {
    return <div>Generating preview...</div>;
  }

  const allPassed = preview.validation.proteinPreserved && preview.validation.noRemainingSites;

  return (
    <div className="step-content preview-step">
      <h2>Pre-Flight Check</h2>

      {/* Success/Warning Banner */}
      <div className={`status-banner ${allPassed ? 'success' : 'warning'}`}>
        <div className="status-icon">
          {allPassed ? Icons.checkCircle : Icons.warning}
        </div>
        <div className="status-content">
          <div className="status-title">
            {allPassed ? 'All Checks Passed!' : 'Review Required'}
          </div>
          <div className="status-subtitle">
            {allPassed
              ? 'Your domestication preserves the protein sequence and removes all internal sites.'
              : 'Some checks did not pass. Please review before proceeding.'}
          </div>
        </div>
        <div className="status-checks">
          <span className={`status-check ${preview.validation.proteinPreserved ? '' : 'failed'}`}>
            {preview.validation.proteinPreserved ? '✓' : '✗'} Protein preserved
          </span>
          <span className={`status-check ${preview.validation.noRemainingSites ? '' : 'failed'}`}>
            {preview.validation.noRemainingSites ? '✓' : '✗'} Sites removed
          </span>
        </div>
      </div>

      {/* Sequence Map showing modifications */}
      <SequenceRuler
        sequence={sequence}
        annotations={previewAnnotations}
        title="Modifications to Apply"
        contextSize={18}
      />

      {/* Protein Verification - Compact Display */}
      <div className="protein-verification">
        <div className="verification-header">
          <h3>Protein Verification</h3>
          <div className={`verification-result ${preview.comparison?.proteinIdentical ? 'success' : 'error'}`}>
            {preview.comparison?.proteinIdentical ? (
              <>{Icons.check} Identical — No amino acid changes</>
            ) : (
              <>{Icons.warning} {preview.comparison?.differences?.length || 0} difference(s) detected</>
            )}
          </div>
        </div>
        <div className="protein-stats">
          <div className="stat">
            <span className="stat-label">Original</span>
            <span className="stat-value">{preview.original.proteinLength} aa</span>
          </div>
          <div className="stat">
            <span className="stat-label">Domesticated</span>
            <span className="stat-value">{preview.domesticated.proteinLength} aa</span>
          </div>
          <div className="stat">
            <span className="stat-label">Identity</span>
            <span className="stat-value">{preview.comparison?.proteinIdentical ? '100%' : '<100%'}</span>
          </div>
        </div>
      </div>

      {/* Mutation Summary - show when using silent mutation strategy */}
      {preview.totalMutations > 0 && preview.mutations && (
        <div className="mutation-summary">
          <h3>Silent Mutations to Apply ({preview.totalMutations})</h3>
          <table className="mutations-table">
            <thead>
              <tr>
                <th>Position</th>
                <th>Context</th>
                <th>Codon Change</th>
                <th>Amino Acid</th>
              </tr>
            </thead>
            <tbody>
              {preview.mutations.map((mut, i) => {
                // Show mutation in sequence context
                const pos = mut.position;
                const before = sequence.slice(Math.max(0, pos - 3), pos);
                const after = sequence.slice(pos + 1, Math.min(sequence.length, pos + 4));
                const oldNt = mut.change?.split('→')[0] || '';
                const newNt = mut.change?.split('→')[1] || '';

                return (
                  <tr key={i}>
                    <td className="pos-col">{pos + 1}</td>
                    <td className="context-col">
                      <code className="mutation-context">
                        <span className="ctx-before">{before}</span>
                        <span className="ctx-old">{oldNt}</span>
                        <span className="ctx-arrow">→</span>
                        <span className="ctx-new">{newNt}</span>
                        <span className="ctx-after">{after}</span>
                      </code>
                    </td>
                    <td><code className="codon-change">{mut.codon}</code></td>
                    <td className="aa-col">{mut.aminoAcid}</td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      )}

      {/* Junction Summary - show when using mutagenic junction strategy */}
      {preview.totalJunctions > 0 && preview.junctions && (
        <div className="junction-summary">
          <h3>Mutagenic Junctions ({preview.totalJunctions})</h3>
          <p className="junction-info">
            Sites will be removed via PCR primers during assembly. No pre-modified template needed.
          </p>
          <table className="junctions-table">
            <thead>
              <tr>
                <th>Site Position</th>
                <th>Junction Position</th>
                <th>Overhang</th>
                <th>Description</th>
              </tr>
            </thead>
            <tbody>
              {preview.junctions.map((jnc, i) => (
                <tr key={i}>
                  <td>{(jnc.sitePosition ?? 0) + 1}</td>
                  <td>{(jnc.position ?? 0) + 1}</td>
                  <td><code>{jnc.overhang || 'N/A'}</code></td>
                  <td>{jnc.description || 'Junction primer'}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      {/* Show message when no changes needed */}
      {preview.totalMutations === 0 && preview.totalJunctions === 0 && (
        <div className="no-changes-message">
          <p>No mutations or junctions selected. Please go back and select options for each site.</p>
        </div>
      )}

      {/* Final Approval */}
      <div className="approval-section">
        {preview.validation?.proteinPreserved && preview.validation?.noRemainingSites ? (
          <div className="approval-ready">
            <span className="approval-icon">{Icons.check}</span>
            <span className="approval-text">All checks passed. Ready to apply domestication.</span>
          </div>
        ) : (
          <div className="approval-warning">
            <span className="approval-icon">{Icons.warning}</span>
            <span className="approval-text">
              Some checks failed. Review carefully before proceeding.
            </span>
          </div>
        )}
      </div>

      <div className="step-actions">
        <button className="btn-secondary" onClick={onBack}>← Back</button>
        <button
          className="btn-primary btn-approve"
          onClick={onApprove}
          disabled={isProcessing || !preview.validation?.proteinPreserved}
        >
          {isProcessing ? 'Applying...' : <>{Icons.check} Approve & Apply Domestication</>}
        </button>
      </div>
    </div>
  );
}

/**
 * Calculate GC% for a sequence
 */
function calculateGC(seq: string): number {
  if (!seq) return 0;
  const gc = (seq.match(/[GCgc]/g) || []).length;
  return (gc / seq.length) * 100;
}

/**
 * Copy text to clipboard with visual feedback
 */
function copyToClipboard(text: string, setCopied: (id: string) => void, id: string) {
  navigator.clipboard.writeText(text);
  setCopied(id);
  setTimeout(() => setCopied(''), 2000);
}

/**
 * Render colored nucleotide sequence
 */
function ColoredSequence({ seq }: { seq: string }) {
  return (
    <span className="colored-sequence">
      {seq.split('').map((nt, i) => (
        <span key={i} className={`nt-${nt.toUpperCase()}`}>{nt}</span>
      ))}
    </span>
  );
}

/**
 * Reusable primer card component for displaying primer pairs
 */
interface PrimerCardProps {
  primerPair: PrimerPairDesign;
  index: number;
  expandedPrimer: number | null;
  setExpandedPrimer: (index: number | null) => void;
  primerType: string;
}

function PrimerCard({ primerPair, index, expandedPrimer, setExpandedPrimer, primerType }: PrimerCardProps) {
  const [copiedId, setCopiedId] = useState<string>('');
  const isExpanded = expandedPrimer === index;
  const isPCR = primerType === 'pcr';

  // Determine display name based on primer type
  let pairName = primerPair.name;
  if (!pairName) {
    if (primerPair.type === 'junction') {
      pairName = `Junction ${index + 1}`;
    } else if (primerPair.type === 'pcr_mutagenesis') {
      pairName = `Mutagenesis Site ${(primerPair.mutationIndex || 0) + 1}`;
    } else if (primerPair.type === 'golden_gate') {
      pairName = 'Assembly Primers';
    } else {
      pairName = `Primer Pair ${index + 1}`;
    }
  }

  // Helper to render a single primer row
  const renderPrimer = (primer: any, direction: 'forward' | 'reverse', label: string) => {
    if (!primer) return null;
    const seq = primer.sequence || '';
    const gcPercent = calculateGC(seq);
    const tm = primer.homologyTm || primer.tm;
    const copyId = `${index}-${direction}`;
    const isCopied = copiedId === copyId;

    return (
      <div className={`primer-row ${direction}`}>
        <div className="primer-row-header">
          <span className={`primer-direction-badge ${direction}`}>
            {direction === 'forward' ? 'FWD' : 'REV'}
          </span>
          <span className="primer-label">{label}</span>
          <div className="primer-metrics-inline">
            <span className="metric"><strong>{seq.length}</strong> nt</span>
            <span className="metric-sep">·</span>
            <span className="metric">Tm <strong>{tm?.toFixed(1) || '—'}°C</strong></span>
            <span className="metric-sep">·</span>
            <span className={`metric ${gcPercent < 40 || gcPercent > 60 ? 'warn' : ''}`}>
              GC <strong>{gcPercent.toFixed(0)}%</strong>
            </span>
          </div>
        </div>
        <div className="primer-sequence-row">
          <code className="primer-sequence-box">
            {isPCR && primer.components ? (
              <>
                <span className="seq-segment fiveprime">{primer.components.fivePrime}</span>
                <span className="seq-segment mutation-highlight">{primer.components.mutation}</span>
                <span className="seq-segment threeprime">{primer.components.threePrime}</span>
              </>
            ) : (
              <ColoredSequence seq={seq} />
            )}
          </code>
          <button
            className={`btn-copy-inline ${isCopied ? 'copied' : ''}`}
            onClick={(e) => {
              e.stopPropagation();
              copyToClipboard(seq, setCopiedId, copyId);
            }}
            title="Copy sequence"
          >
            {isCopied ? Icons.check : Icons.copy}
          </button>
        </div>
      </div>
    );
  };

  return (
    <div className={`primer-card-v2 ${isExpanded ? 'expanded' : ''} ${primerType}`}>
      <div className="primer-card-header-v2" onClick={() => setExpandedPrimer(isExpanded ? null : index)}>
        <div className="primer-card-title-section">
          <span className="primer-pair-icon">{Icons.dna}</span>
          <span className="primer-pair-name">{pairName}</span>
          {primerPair.codonChange && (
            <span className="primer-tag mutation">{primerPair.codonChange}</span>
          )}
          {primerPair.overhang && (
            <span className="primer-tag overhang">{primerPair.overhang}</span>
          )}
        </div>
        <div className="primer-card-actions">
          {primerPair.pairQuality !== undefined && (
            <span className={`quality-badge ${primerPair.pairQuality >= 80 ? 'excellent' : primerPair.pairQuality >= 65 ? 'good' : 'acceptable'}`}>
              {primerPair.pairQuality >= 80 ? '★' : primerPair.pairQuality >= 65 ? '●' : '○'} {primerPair.pairQuality.toFixed(0)}%
            </span>
          )}
          <button className="expand-toggle">{isExpanded ? '▲' : '▼'}</button>
        </div>
      </div>

      {/* Always show primer sequences */}
      <div className="primer-sequences-preview">
        {renderPrimer(primerPair.forward, 'forward', `${pairName}_F`)}
        {renderPrimer(primerPair.reverse, 'reverse', `${pairName}_R`)}
      </div>

      {isExpanded && (
        <div className="primer-expanded-details">
          {/* Instructions for PCR primers */}
          {isPCR && primerPair.instructions && (
            <div className="primer-protocol-note">
              <span className="protocol-icon">{Icons.info}</span>
              <span><strong>Protocol:</strong> {primerPair.instructions}</span>
            </div>
          )}

          {/* PCR protocol info */}
          {isPCR && primerPair.pcrProtocol && (
            <div className="pcr-conditions">
              <span className="condition">
                <span className="condition-label">Annealing:</span>
                <span className="condition-value">{primerPair.pcrProtocol.annealingTemp}°C</span>
              </span>
              <span className="condition">
                <span className="condition-label">Extension:</span>
                <span className="condition-value">{primerPair.pcrProtocol.extensionTime}</span>
              </span>
            </div>
          )}

          {/* Copy Both Button */}
          <button
            className="btn-copy-both"
            onClick={(e) => {
              e.stopPropagation();
              const text = [
                primerPair.forward ? `>${pairName}_F\n${primerPair.forward.sequence}` : '',
                primerPair.reverse ? `>${pairName}_R\n${primerPair.reverse.sequence}` : ''
              ].filter(Boolean).join('\n');
              copyToClipboard(text, setCopiedId, `${index}-both`);
            }}
          >
            {copiedId === `${index}-both` ? <>{Icons.check} Copied!</> : <>{Icons.copy} Copy Both (FASTA)</>}
          </button>
        </div>
      )}
    </div>
  );
}

interface PrimersStepProps {
  executionResult: ExecutionResult | null;
  primerDesign: PrimerDesignResult | null;
  enzyme: string;
  onContinue: () => void;
  onBack: () => void;
}

function PrimersStep({ executionResult, primerDesign, onContinue, onBack }: PrimersStepProps) {
  const [expandedPrimer, setExpandedPrimer] = useState<number | null>(null);
  const [copiedState, setCopiedState] = useState<string>('');

  // Download helper functions
  const downloadFile = (content: string, filename: string, type: string) => {
    const blob = new Blob([content], { type });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const generateFASTA = (orderList: any[]) => {
    return orderList.map(p => `>${p.name}\n${p.sequence}`).join('\n\n');
  };

  const generateCSV = (orderList: any[]) => {
    const header = 'Name,Sequence,Length,Tm,GC%';
    const rows = orderList.map(p => {
      const gc = calculateGC(p.sequence);
      return `${p.name},${p.sequence},${p.length},${p.tm?.toFixed(1) || 'N/A'},${gc.toFixed(1)}`;
    });
    return [header, ...rows].join('\n');
  };

  if (!executionResult) {
    return (
      <div className="step-content primers-step">
        <h2>Primer Design</h2>
        <div className="no-result">
          <p>No domestication result available.</p>
        </div>
        <div className="step-actions">
          <button className="btn-secondary" onClick={onBack}>← Back</button>
        </div>
      </div>
    );
  }

  const { primers, primerSummary } = primerDesign || {};
  const hasJunctions = executionResult.junctions && executionResult.junctions.length > 0;
  const hasMutations = executionResult.mutations && executionResult.mutations.length > 0;
  const totalPrimers = primerSummary?.totalPrimers || 0;

  return (
    <div className="step-content primers-step-v2">
      {/* Success Header */}
      <div className="primers-header">
        <div className="primers-title-section">
          <div className="primers-icon">{Icons.dna}</div>
          <div>
            <h2>Primers Ready</h2>
            <p className="primers-subtitle">
              {hasJunctions
                ? `${executionResult.junctions.length} junction primer pair${executionResult.junctions.length !== 1 ? 's' : ''} designed`
                : hasMutations
                ? `${executionResult.mutations.length} silent mutation${executionResult.mutations.length !== 1 ? 's' : ''} applied`
                : 'Domestication complete'}
            </p>
          </div>
        </div>
        {totalPrimers > 0 && (
          <div className="primers-quick-stats">
            <div className="quick-stat">
              <span className="quick-stat-value">{totalPrimers}</span>
              <span className="quick-stat-label">Primers</span>
            </div>
            <div className="quick-stat">
              <span className="quick-stat-value">{primerSummary?.totalLength || 0}</span>
              <span className="quick-stat-label">Total nt</span>
            </div>
            {primerSummary?.estimatedCost && (
              <div className="quick-stat highlight">
                <span className="quick-stat-value">${primerSummary.estimatedCost.toFixed(0)}</span>
                <span className="quick-stat-label">Est. Cost</span>
              </div>
            )}
          </div>
        )}
      </div>

      {/* Strategy Badge */}
      <div className="strategy-badge-inline">
        <span className="strategy-badge-icon">{hasJunctions ? Icons.scissors : Icons.dna}</span>
        <span className="strategy-badge-text">
          {hasJunctions ? 'Mutagenic Junction Strategy' : 'Silent Mutation Strategy'}
        </span>
        {hasJunctions && (
          <span className="strategy-badge-detail">• {executionResult.junctions.length + 1} fragments</span>
        )}
      </div>

      {/* Primer List - Group by step for silent mutation strategy */}
      {primers && primers.length > 0 ? (
        <div className="primer-list">
          {/* Check if we have multi-step workflow (silent mutation) */}
          {(() => {
            const pcrPrimers = primers.filter((p: PrimerPairDesign) => p.type === 'pcr_mutagenesis');
            const ggPrimers = primers.filter((p: PrimerPairDesign) => p.type === 'golden_gate' || p.type === 'junction');
            const hasMultiStep = pcrPrimers.length > 0 && ggPrimers.length > 0;

            if (hasMultiStep) {
              return (
                <>
                  {/* Step 1: PCR Mutagenesis */}
                  <div className="primer-step-section">
                    <div className="step-header pcr">
                      <span className="step-badge">Step 1</span>
                      <h3>Site-Directed Mutagenesis (PCR)</h3>
                      <p className="step-description">
                        Use these primers with overlap extension PCR or QuikChange to introduce silent mutations
                      </p>
                    </div>
                    {pcrPrimers.map((primerPair: PrimerPairDesign, i: number) => (
                      <PrimerCard
                        key={`pcr-${i}`}
                        primerPair={primerPair}
                        index={i}
                        expandedPrimer={expandedPrimer}
                        setExpandedPrimer={setExpandedPrimer}
                        primerType="pcr"
                      />
                    ))}
                  </div>

                  {/* Step 2: Golden Gate Assembly */}
                  <div className="primer-step-section">
                    <div className="step-header gg">
                      <span className="step-badge">Step 2</span>
                      <h3>Golden Gate Assembly</h3>
                      <p className="step-description">
                        After obtaining mutated template, use these primers for final assembly
                      </p>
                    </div>
                    {ggPrimers.map((primerPair: PrimerPairDesign, i: number) => (
                      <PrimerCard
                        key={`gg-${i}`}
                        primerPair={primerPair}
                        index={pcrPrimers.length + i}
                        expandedPrimer={expandedPrimer}
                        setExpandedPrimer={setExpandedPrimer}
                        primerType="golden_gate"
                      />
                    ))}
                  </div>
                </>
              );
            }

            // Single-step workflow (mutagenic junction - one-pot)
            return (
              <>
                <h3>{hasJunctions ? 'Mutagenic Junction Primers' : 'Assembly Primers'}</h3>
                {hasJunctions && (
                  <p className="one-pot-note">
                    {Icons.check} One-pot compatible: Mutations are introduced during Golden Gate assembly
                  </p>
                )}
                {primers.map((primerPair: PrimerPairDesign, i: number) => {
                  if (primerPair.type === 'silent_mutation_info') {
                    return (
                      <div key={i} className="primer-card info-card">
                        <div className="primer-header">
                          <span className="primer-name">Silent Mutation Applied</span>
                          <span className="primer-badge info">Template modified</span>
                        </div>
                      </div>
                    );
                  }
                  return (
                    <PrimerCard
                      key={i}
                      primerPair={primerPair}
                      index={i}
                      expandedPrimer={expandedPrimer}
                      setExpandedPrimer={setExpandedPrimer}
                      primerType="junction"
                    />
                  );
                })}
              </>
            );
          })()}
        </div>
      ) : (
        <div className="no-primers-needed">
          <div className="info-icon">{Icons.info}</div>
          <p>
            {hasMutations
              ? 'Silent mutations applied directly to template. Use your existing primers or design standard amplification primers.'
              : 'No additional primers needed for this domestication strategy.'}
          </p>
        </div>
      )}

      {/* Order & Export Section */}
      {primerSummary && primerSummary.orderList && primerSummary.orderList.length > 0 && (
        <div className="order-export-section">
          <div className="order-section-header">
            <h3>{Icons.download} Ready to Order</h3>
            <div className="export-buttons">
              <button
                className={`btn-export ${copiedState === 'clipboard' ? 'copied' : ''}`}
                onClick={() => {
                  const text = primerSummary.orderList!
                    .map(p => `${p.name}\t${p.sequence}`)
                    .join('\n');
                  navigator.clipboard.writeText(text);
                  setCopiedState('clipboard');
                  setTimeout(() => setCopiedState(''), 2000);
                }}
              >
                {copiedState === 'clipboard' ? <>{Icons.check} Copied!</> : <>{Icons.copy} Copy All</>}
              </button>
              <button
                className="btn-export"
                onClick={() => {
                  const csv = generateCSV(primerSummary.orderList!);
                  downloadFile(csv, 'domestication_primers.csv', 'text/csv');
                }}
              >
                {Icons.download} CSV
              </button>
              <button
                className="btn-export"
                onClick={() => {
                  const fasta = generateFASTA(primerSummary.orderList!);
                  downloadFile(fasta, 'domestication_primers.fasta', 'text/plain');
                }}
              >
                {Icons.download} FASTA
              </button>
            </div>
          </div>

          <div className="order-table-v2">
            <table>
              <thead>
                <tr>
                  <th className="col-name">Name</th>
                  <th className="col-sequence">Sequence (5' → 3')</th>
                  <th className="col-length">Length</th>
                  <th className="col-tm">Tm</th>
                  <th className="col-gc">GC%</th>
                  <th className="col-copy"></th>
                </tr>
              </thead>
              <tbody>
                {primerSummary.orderList.map((item, i) => {
                  const gc = calculateGC(item.sequence);
                  const rowCopyId = `row-${i}`;
                  return (
                    <tr key={i} className={i % 2 === 0 ? 'even' : 'odd'}>
                      <td className="col-name">
                        <span className="primer-name-cell">{item.name}</span>
                      </td>
                      <td className="col-sequence">
                        <code className="sequence-mono">
                          <ColoredSequence seq={item.sequence} />
                        </code>
                      </td>
                      <td className="col-length">{item.length} nt</td>
                      <td className="col-tm">{item.tm?.toFixed(1) || '—'}°C</td>
                      <td className={`col-gc ${gc < 40 || gc > 60 ? 'warn' : ''}`}>{gc.toFixed(0)}%</td>
                      <td className="col-copy">
                        <button
                          className={`btn-copy-row ${copiedState === rowCopyId ? 'copied' : ''}`}
                          onClick={() => {
                            navigator.clipboard.writeText(item.sequence);
                            setCopiedState(rowCopyId);
                            setTimeout(() => setCopiedState(''), 2000);
                          }}
                          title="Copy sequence"
                        >
                          {copiedState === rowCopyId ? Icons.check : Icons.copy}
                        </button>
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>

          <div className="order-summary-footer">
            <span className="footer-stat">{primerSummary.totalPrimers} primers</span>
            <span className="footer-sep">•</span>
            <span className="footer-stat">{primerSummary.totalLength} nt total</span>
            {primerSummary.estimatedCost && (
              <>
                <span className="footer-sep">•</span>
                <span className="footer-stat highlight">~${primerSummary.estimatedCost.toFixed(2)} estimated</span>
              </>
            )}
          </div>
        </div>
      )}

      <div className="step-actions">
        <button className="btn-secondary" onClick={onBack}>← Back</button>
        <button className="btn-primary" onClick={onContinue}>
          Complete Workflow →
        </button>
      </div>
    </div>
  );
}

interface CompleteStepProps {
  onClose?: () => void;
}

function CompleteStep({ onClose }: CompleteStepProps) {
  return (
    <div className="step-content complete-step">
      <div className="success-animation">
        <div className="success-icon">{Icons.checkCircle}</div>
        <h2>Domestication Complete!</h2>
        <p>All internal sites have been successfully removed while preserving the protein sequence. Primers have been designed and are ready for ordering.</p>
      </div>

      <div className="step-actions">
        <button className="btn-primary" onClick={onClose}>
          Done
        </button>
      </div>
    </div>
  );
}

// ============================================================================
// STYLES (CSS-in-JS for portability)
// ============================================================================

export const EnhancedDomesticationStyles = `
.enhanced-domestication-panel {
  max-width: 900px;
  margin: 0 auto;
  padding: 20px;
  font-family: system-ui, -apple-system, sans-serif;
}

.progress-indicator {
  display: flex;
  justify-content: space-between;
  margin-bottom: 30px;
  padding: 0 20px;
}

.progress-step {
  display: flex;
  flex-direction: column;
  align-items: center;
  flex: 1;
  position: relative;
}

.progress-step::after {
  content: '';
  position: absolute;
  top: 14px;
  left: 50%;
  width: 100%;
  height: 3px;
  background: #e5e7eb;
  border-radius: 2px;
}

.progress-step:last-child::after {
  display: none;
}

.progress-step.complete::after {
  background: linear-gradient(90deg, #22c55e 0%, #4ade80 100%);
}

.step-number {
  width: 28px;
  height: 28px;
  border-radius: 50%;
  background: #f3f4f6;
  color: #9ca3af;
  display: flex;
  align-items: center;
  justify-content: center;
  font-weight: 600;
  font-size: 12px;
  position: relative;
  z-index: 1;
  border: 2px solid #e5e7eb;
  transition: all 0.2s;
}

.progress-step.complete .step-number {
  background: #22c55e;
  color: white;
  border-color: #22c55e;
  box-shadow: 0 2px 4px rgba(34, 197, 94, 0.3);
}

.progress-step.current .step-number {
  background: #3b82f6;
  color: white;
  border-color: #3b82f6;
  box-shadow: 0 2px 8px rgba(59, 130, 246, 0.4);
  transform: scale(1.1);
}

.step-label {
  margin-top: 6px;
  font-size: 11px;
  color: #9ca3af;
  font-weight: 500;
}

.progress-step.complete .step-label {
  color: #22c55e;
}

.progress-step.current .step-label {
  color: #3b82f6;
  font-weight: 600;
}

.step-content {
  background: white;
  border-radius: 12px;
  padding: 24px;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.step-content h2 {
  margin: 0 0 8px;
  font-size: 1.5rem;
  color: #111827;
}

.step-description {
  color: #6b7280;
  margin-bottom: 24px;
}

.error-banner {
  background: #fef2f2;
  border: 1px solid #fecaca;
  color: #dc2626;
  padding: 12px 16px;
  border-radius: 8px;
  margin-bottom: 16px;
  display: flex;
  align-items: center;
  gap: 8px;
}

.error-banner button {
  margin-left: auto;
  background: none;
  border: none;
  cursor: pointer;
  color: #dc2626;
}

.analysis-summary {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
  gap: 16px;
  margin-bottom: 24px;
}

.summary-card {
  display: flex;
  gap: 12px;
  padding: 16px;
  background: #f9fafb;
  border-radius: 8px;
}

.card-icon {
  width: 40px;
  height: 40px;
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 20px;
}

.card-icon.warning { background: #fef3c7; }
.card-icon.info { background: #dbeafe; }
.card-icon.success { background: #dcfce7; }

.card-title {
  font-weight: 600;
  color: #111827;
}

.card-detail {
  font-size: 13px;
  color: #6b7280;
}

.sites-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 14px;
}

.sites-table th,
.sites-table td {
  padding: 10px 12px;
  text-align: left;
  border-bottom: 1px solid #e5e7eb;
}

.sites-table th {
  background: #f9fafb;
  font-weight: 500;
  color: #374151;
}

.sites-table code {
  background: #f3f4f6;
  padding: 2px 6px;
  border-radius: 4px;
  font-family: 'SF Mono', Monaco, monospace;
}

.context .highlight {
  background: #fef3c7;
  color: #92400e;
  font-weight: 600;
}

.frame-options {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 16px;
  margin-bottom: 24px;
}

.frame-option {
  border: 2px solid #e5e7eb;
  border-radius: 12px;
  padding: 16px;
  cursor: pointer;
  transition: all 0.2s;
}

.frame-option:hover {
  border-color: #93c5fd;
}

.frame-option.selected {
  border-color: #3b82f6;
  background: #eff6ff;
}

.frame-option.recommended {
  border-color: #22c55e;
}

.frame-header {
  display: flex;
  align-items: center;
  gap: 8px;
  margin-bottom: 12px;
}

.frame-number {
  font-weight: 600;
  font-size: 1.1rem;
}

.recommended-badge,
.selected-badge {
  font-size: 11px;
  padding: 2px 8px;
  border-radius: 12px;
}

.recommended-badge {
  background: #dcfce7;
  color: #166534;
}

.selected-badge {
  background: #dbeafe;
  color: #1d4ed8;
}

.frame-details .detail-row {
  display: flex;
  justify-content: space-between;
  font-size: 13px;
  margin-bottom: 4px;
}

.frame-details .label {
  color: #6b7280;
}

.protein-preview {
  margin-top: 8px;
  font-size: 12px;
}

.protein-preview code {
  display: block;
  background: #f3f4f6;
  padding: 4px 8px;
  border-radius: 4px;
  margin-top: 4px;
  overflow: hidden;
  text-overflow: ellipsis;
}

.warning-row {
  color: #dc2626;
  font-size: 12px;
  margin-top: 8px;
}

.codon-mode-section {
  margin-top: 24px;
  padding-top: 24px;
  border-top: 1px solid #e5e7eb;
}

.codon-mode-section h3 {
  margin: 0 0 16px;
  font-size: 1rem;
}

.mode-options {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 12px;
}

.mode-options label {
  display: flex;
  align-items: flex-start;
  gap: 12px;
  padding: 12px;
  border: 2px solid #e5e7eb;
  border-radius: 8px;
  cursor: pointer;
  transition: all 0.2s;
}

.mode-options label:hover {
  border-color: #93c5fd;
}

.mode-options label.selected {
  border-color: #3b82f6;
  background: #eff6ff;
}

.mode-options input[type="radio"] {
  margin-top: 2px;
}

.mode-title {
  font-weight: 500;
}

.mode-desc {
  font-size: 12px;
  color: #6b7280;
}

.site-mutation-card {
  border: 1px solid #e5e7eb;
  border-radius: 12px;
  margin-bottom: 16px;
  overflow: hidden;
}

.site-header {
  background: #f9fafb;
  padding: 12px 16px;
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.site-position {
  font-weight: 500;
}

.site-sequence {
  background: #fef3c7;
  padding: 4px 8px;
  border-radius: 4px;
}

.mutation-options {
  padding: 16px;
  display: grid;
  gap: 12px;
}

.mutation-option {
  border: 2px solid #e5e7eb;
  border-radius: 8px;
  padding: 12px;
  cursor: pointer;
  transition: all 0.2s;
}

.mutation-option:hover {
  border-color: #93c5fd;
}

.mutation-option.selected {
  border-color: #3b82f6;
  background: #eff6ff;
}

.mutation-option.recommended {
  position: relative;
}

.option-header {
  display: flex;
  align-items: center;
  gap: 8px;
  margin-bottom: 8px;
}

.option-type {
  font-weight: 500;
}

.option-score {
  font-size: 12px;
  color: #6b7280;
}

.recommended-tag {
  font-size: 11px;
  background: #dcfce7;
  color: #166534;
  padding: 2px 8px;
  border-radius: 12px;
  margin-left: auto;
}

.option-details .detail {
  display: flex;
  gap: 8px;
  font-size: 13px;
  margin-bottom: 4px;
}

.option-details .label {
  color: #6b7280;
}

.option-details .aa {
  color: #6b7280;
}

.warnings {
  margin-top: 8px;
  display: flex;
  gap: 4px;
  flex-wrap: wrap;
}

.warning-badge {
  font-size: 11px;
  background: #fef3c7;
  color: #92400e;
  padding: 2px 6px;
  border-radius: 4px;
}

.validation-status {
  display: flex;
  gap: 24px;
  margin-bottom: 24px;
}

.validation-item {
  display: flex;
  align-items: center;
  gap: 8px;
  padding: 12px 16px;
  border-radius: 8px;
  font-weight: 500;
}

.validation-item.passed {
  background: #dcfce7;
  color: #166534;
}

.validation-item.failed {
  background: #fef2f2;
  color: #dc2626;
}

.protein-alignment {
  background: #f9fafb;
  border-radius: 8px;
  padding: 16px;
  margin-bottom: 24px;
}

.protein-alignment h3 {
  margin: 0 0 12px;
  font-size: 1rem;
}

.alignment-row {
  display: flex;
  align-items: center;
  gap: 12px;
  margin-bottom: 8px;
}

.row-label {
  width: 100px;
  font-weight: 500;
  color: #6b7280;
}

.protein-sequence {
  flex: 1;
  background: white;
  padding: 8px 12px;
  border-radius: 4px;
  font-family: 'SF Mono', Monaco, monospace;
  font-size: 12px;
  overflow-x: auto;
  white-space: nowrap;
}

.alignment-status {
  margin-top: 12px;
  padding-top: 12px;
  border-top: 1px solid #e5e7eb;
}

.alignment-status .identical {
  color: #166534;
  font-weight: 600;
}

.alignment-status .different {
  color: #dc2626;
  font-weight: 600;
}

.mutation-summary h3,
.junction-summary h3 {
  margin: 0 0 12px;
  font-size: 1rem;
}

.junction-summary {
  margin-top: 20px;
}

.junction-info {
  font-size: 13px;
  color: #6b7280;
  margin-bottom: 12px;
  padding: 8px 12px;
  background: #f0f9ff;
  border-radius: 6px;
  border-left: 3px solid #3b82f6;
}

.no-changes-message {
  padding: 16px;
  background: #fef3c7;
  border-radius: 8px;
  color: #92400e;
  text-align: center;
}

.no-changes-message p {
  margin: 0;
}

.mutations-table,
.junctions-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 14px;
  margin-bottom: 24px;
}

.mutations-table th,
.mutations-table td,
.junctions-table th,
.junctions-table td {
  padding: 8px 12px;
  text-align: left;
  border-bottom: 1px solid #e5e7eb;
}

.mutations-table th,
.junctions-table th {
  background: #f9fafb;
}

.approval-section {
  padding: 16px;
  border-radius: 8px;
  margin-bottom: 24px;
}

.approval-ready {
  background: #dcfce7;
  padding: 16px;
  border-radius: 8px;
  display: flex;
  align-items: center;
  gap: 12px;
  color: #166534;
}

.approval-warning {
  background: #fef3c7;
  padding: 16px;
  border-radius: 8px;
  display: flex;
  align-items: center;
  gap: 12px;
  color: #92400e;
}

.approval-icon {
  font-size: 24px;
}

.step-actions {
  display: flex;
  justify-content: flex-end;
  gap: 12px;
  margin-top: 24px;
  padding-top: 24px;
  border-top: 1px solid #e5e7eb;
}

.btn-primary,
.btn-secondary {
  padding: 10px 20px;
  border-radius: 8px;
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

.btn-primary.btn-approve {
  background: #22c55e;
}

.btn-primary.btn-approve:hover {
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

.complete-step {
  text-align: center;
  padding: 48px 24px;
}

.success-animation .success-icon {
  width: 80px;
  height: 80px;
  background: #dcfce7;
  color: #22c55e;
  font-size: 40px;
  border-radius: 50%;
  display: flex;
  align-items: center;
  justify-content: center;
  margin: 0 auto 24px;
}

.complete-step h2 {
  color: #166534;
}

.complete-step .step-actions {
  justify-content: center;
}

.domestication-panel.no-sites {
  text-align: center;
  padding: 48px;
}

.domestication-panel.no-sites .success-icon {
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

/* Strategy Selection Step Styles */
.strategy-step h2 {
  margin-bottom: 8px;
}

.strategy-options {
  display: flex;
  flex-direction: column;
  gap: 12px;
  margin-bottom: 24px;
}

.strategy-card {
  border: 2px solid #e5e7eb;
  border-radius: 12px;
  padding: 16px 20px;
  cursor: pointer;
  transition: all 0.2s;
  position: relative;
  overflow: hidden;
}

.strategy-card::before {
  content: '';
  position: absolute;
  top: 0;
  left: 0;
  width: 4px;
  height: 100%;
  background: transparent;
  transition: background 0.2s;
}

.strategy-card:hover:not(.disabled) {
  border-color: #93c5fd;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
}

.strategy-card.selected {
  border-color: #3b82f6;
  background: linear-gradient(135deg, #eff6ff 0%, #f0f9ff 100%);
}

.strategy-card.selected::before {
  background: #3b82f6;
}

.strategy-card.recommended {
  border-color: #86efac;
  background: linear-gradient(135deg, #f0fdf4 0%, #ecfdf5 100%);
}

.strategy-card.recommended::before {
  background: #22c55e;
}

.strategy-card.recommended.selected {
  border-color: #22c55e;
}

.strategy-card.disabled {
  opacity: 0.5;
  cursor: not-allowed;
}

.strategy-header {
  display: flex;
  align-items: flex-start;
  gap: 16px;
  margin-bottom: 12px;
}

.strategy-icon {
  font-size: 24px;
  width: 40px;
  height: 40px;
  display: flex;
  align-items: center;
  justify-content: center;
  background: #f3f4f6;
  border-radius: 10px;
  flex-shrink: 0;
}

.strategy-card.recommended .strategy-icon {
  background: #dcfce7;
}

.strategy-card.selected .strategy-icon {
  background: #dbeafe;
}

.strategy-title {
  flex: 1;
}

.strategy-title h3 {
  margin: 0 0 2px;
  font-size: 1rem;
  font-weight: 600;
  display: flex;
  align-items: center;
  gap: 8px;
}

.strategy-tagline {
  font-size: 11px;
  color: #6b7280;
  font-weight: 500;
  text-transform: uppercase;
  letter-spacing: 0.02em;
}

.strategy-card.recommended .strategy-tagline {
  color: #15803d;
}

.recommended-badge {
  font-size: 10px;
  background: #22c55e;
  color: white;
  padding: 2px 8px;
  border-radius: 10px;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 0.02em;
}

.selected-badge {
  width: 24px;
  height: 24px;
  background: #3b82f6;
  border-radius: 50%;
  display: flex;
  align-items: center;
  justify-content: center;
  color: white;
  flex-shrink: 0;
}

.selected-badge svg {
  width: 14px;
  height: 14px;
}

.strategy-description {
  color: #4b5563;
  margin-bottom: 12px;
  font-size: 13px;
  line-height: 1.5;
  padding-left: 56px;
}

.strategy-details {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 12px;
  padding-left: 56px;
}

.strategy-details h4 {
  font-size: 11px;
  margin: 0 0 6px;
  color: #6b7280;
  text-transform: uppercase;
  letter-spacing: 0.03em;
  font-weight: 600;
}

.strategy-details ul {
  list-style: none;
  padding: 0;
  margin: 0;
}

.strategy-details li {
  font-size: 12px;
  margin-bottom: 3px;
  display: flex;
  align-items: flex-start;
  gap: 4px;
  line-height: 1.4;
}

.strategy-details .benefit {
  color: #166534;
}

.strategy-details .benefit svg {
  width: 12px;
  height: 12px;
  flex-shrink: 0;
  margin-top: 2px;
}

.strategy-details .tradeoff {
  color: #92400e;
}

.strategy-unavailable {
  margin-top: 12px;
  margin-left: 56px;
  padding: 8px 12px;
  background: #fef2f2;
  border-radius: 6px;
  color: #dc2626;
  font-size: 12px;
  text-align: center;
}

/* Enzyme Selection Section Styles */
.enzyme-selection-section {
  margin-top: 24px;
  padding: 20px;
  background: #f8fafc;
  border: 1px solid #e2e8f0;
  border-radius: 12px;
}

.enzyme-selection-section h3 {
  margin: 0 0 8px;
  font-size: 16px;
  color: #1e293b;
}

.enzyme-selection-section .section-description {
  margin: 0 0 16px;
  font-size: 13px;
  color: #64748b;
}

.enzyme-comparison-table {
  overflow-x: auto;
}

.enzyme-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}

.enzyme-table th,
.enzyme-table td {
  padding: 10px 12px;
  text-align: left;
  border-bottom: 1px solid #e2e8f0;
}

.enzyme-table th {
  background: #f1f5f9;
  font-weight: 600;
  color: #475569;
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 0.5px;
}

.enzyme-table tr.current-enzyme {
  background: #fef3c7;
}

.enzyme-table tr.current-enzyme td {
  border-bottom-color: #fde68a;
}

.enzyme-table tr.enzyme-option {
  cursor: pointer;
  transition: background 0.15s;
}

.enzyme-table tr.enzyme-option.compatible:hover {
  background: #ecfdf5;
}

.enzyme-table tr.enzyme-option.incompatible {
  opacity: 0.6;
  cursor: not-allowed;
}

.enzyme-table tr.enzyme-option.selected {
  background: #dbeafe;
}

.enzyme-table tr.enzyme-option.selected td {
  border-bottom-color: #93c5fd;
}

.current-badge {
  display: inline-block;
  margin-left: 8px;
  padding: 2px 6px;
  background: #fbbf24;
  color: #78350f;
  font-size: 10px;
  font-weight: 600;
  border-radius: 4px;
  text-transform: uppercase;
}

.enzyme-fullname {
  display: block;
  font-size: 11px;
  color: #64748b;
  font-weight: normal;
}

.site-count {
  font-weight: 600;
}

.site-count.success {
  color: #16a34a;
}

.site-count.warning {
  color: #dc2626;
}

.status-badge {
  display: inline-block;
  padding: 3px 8px;
  border-radius: 12px;
  font-size: 11px;
  font-weight: 500;
}

.status-badge.compatible {
  background: #dcfce7;
  color: #166534;
}

.status-badge.incompatible {
  background: #fee2e2;
  color: #dc2626;
}

.no-alternatives-warning {
  display: flex;
  align-items: center;
  gap: 12px;
  padding: 16px;
  background: #fef3c7;
  border: 1px solid #fde68a;
  border-radius: 8px;
  color: #92400e;
  margin-top: 16px;
}

.no-alternatives-warning .warning-icon {
  flex-shrink: 0;
}

/* Primers Step Styles */
.primers-step .strategy-summary-card {
  background: #f9fafb;
  border: 1px solid #e5e7eb;
  border-radius: 12px;
  padding: 16px;
  margin-bottom: 24px;
}

.primers-step .summary-header {
  display: flex;
  align-items: center;
  gap: 12px;
  margin-bottom: 12px;
}

.primers-step .strategy-icon {
  font-size: 24px;
}

.primers-step .strategy-name {
  font-size: 16px;
  font-weight: 600;
  color: #111827;
}

.primers-step .summary-details {
  display: flex;
  gap: 24px;
}

.primers-step .detail-item {
  display: flex;
  align-items: center;
  gap: 6px;
  font-size: 14px;
}

.primers-step .detail-item .label {
  color: #6b7280;
}

.primers-step .detail-item .value {
  font-weight: 500;
  color: #111827;
}

.primers-step .primer-list {
  margin-bottom: 24px;
}

.primers-step .primer-list h3 {
  margin: 0 0 12px;
  font-size: 16px;
}

.primers-step .primer-card {
  border: 1px solid #e5e7eb;
  border-radius: 8px;
  margin-bottom: 8px;
  overflow: hidden;
}

.primers-step .primer-card.info-card {
  background: #f0f9ff;
  border-color: #93c5fd;
}

.primers-step .primer-header {
  display: flex;
  align-items: center;
  gap: 12px;
  padding: 12px 16px;
  background: #f9fafb;
  cursor: pointer;
}

.primers-step .primer-card.info-card .primer-header {
  cursor: default;
  background: transparent;
}

.primers-step .primer-name {
  flex: 1;
  font-weight: 500;
}

.primers-step .primer-badge {
  font-size: 11px;
  padding: 2px 8px;
  border-radius: 10px;
  font-weight: 500;
}

.primers-step .primer-badge.info {
  background: #dbeafe;
  color: #1e40af;
}

.primers-step .primer-badge.overhang {
  background: #dcfce7;
  color: #166534;
  font-family: monospace;
}

.primers-step .quality-indicator {
  font-size: 11px;
  padding: 2px 8px;
  border-radius: 10px;
  font-weight: 600;
}

.primers-step .quality-indicator.excellent {
  background: #dcfce7;
  color: #166534;
}

.primers-step .quality-indicator.good {
  background: #dbeafe;
  color: #1e40af;
}

.primers-step .quality-indicator.acceptable {
  background: #fef3c7;
  color: #92400e;
}

.primers-step .expand-btn {
  width: 24px;
  height: 24px;
  border: none;
  background: #e5e7eb;
  border-radius: 4px;
  cursor: pointer;
  font-size: 14px;
}

.primers-step .primer-details {
  padding: 16px;
  background: white;
  border-top: 1px solid #e5e7eb;
}

.primers-step .individual-primer {
  padding: 12px;
  background: #f9fafb;
  border-radius: 6px;
  margin-bottom: 12px;
}

.primers-step .individual-primer:last-child {
  margin-bottom: 0;
}

.primers-step .primer-direction {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 8px;
}

.primers-step .direction-label {
  font-weight: 600;
  font-size: 13px;
  color: #374151;
}

.primers-step .primer-metrics {
  font-size: 12px;
  color: #6b7280;
}

.primers-step .sequence-display code {
  display: block;
  background: #f3f4f6;
  padding: 8px;
  border-radius: 4px;
  font-size: 12px;
  word-break: break-all;
}

/* Mutation highlighting in primer sequences */
.seq-segment {
  font-family: 'SF Mono', Monaco, 'Consolas', monospace;
}

.seq-segment.fiveprime {
  color: #6b7280;
}

.seq-segment.mutation-highlight {
  background: #fef08a;
  color: #92400e;
  font-weight: 600;
  padding: 1px 2px;
  border-radius: 2px;
}

.seq-segment.threeprime {
  color: #6b7280;
}

/* PCR protocol info */
.pcr-protocol {
  padding: 8px 12px;
  background: #fffbeb;
  border: 1px solid #fde68a;
  border-radius: 6px;
  font-size: 12px;
  color: #92400e;
  margin-top: 12px;
}

.pcr-protocol .protocol-label {
  font-weight: 600;
  margin-right: 4px;
}

/* Primer instructions */
.primer-instructions {
  padding: 8px 12px;
  background: #f0f9ff;
  border-left: 3px solid #3b82f6;
  font-size: 13px;
  color: #1e40af;
  margin-bottom: 12px;
}

/* Step sections for multi-step workflow (Silent Mutation) */
.primer-step-section {
  margin-bottom: 24px;
  border: 1px solid #e5e7eb;
  border-radius: 12px;
  overflow: hidden;
}

.primer-step-section .step-header {
  padding: 16px 20px;
  border-bottom: 1px solid #e5e7eb;
}

.primer-step-section .step-header.pcr {
  background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%);
}

.primer-step-section .step-header.gg {
  background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%);
}

.primer-step-section .step-badge {
  display: inline-block;
  padding: 4px 10px;
  background: rgba(0, 0, 0, 0.15);
  border-radius: 12px;
  font-size: 11px;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 0.5px;
  margin-bottom: 8px;
}

.primer-step-section .step-header.pcr .step-badge {
  background: rgba(217, 119, 6, 0.2);
  color: #92400e;
}

.primer-step-section .step-header.gg .step-badge {
  background: rgba(59, 130, 246, 0.2);
  color: #1e40af;
}

.primer-step-section .step-header h3 {
  margin: 0 0 4px;
  font-size: 16px;
  font-weight: 600;
}

.primer-step-section .step-header .step-description {
  margin: 0;
  font-size: 13px;
  color: #6b7280;
}

.primer-step-section .step-header.pcr .step-description {
  color: #78350f;
}

.primer-step-section .step-header.gg .step-description {
  color: #1e3a8a;
}

/* One-pot note styling */
.one-pot-note {
  display: flex;
  align-items: center;
  gap: 8px;
  padding: 12px 16px;
  background: #ecfdf5;
  border: 1px solid #a7f3d0;
  border-radius: 8px;
  color: #065f46;
  font-size: 14px;
  margin-bottom: 16px;
}

.no-primers-needed {
  padding: 24px;
  background: #f0f9ff;
  border-radius: 8px;
  text-align: center;
  margin-bottom: 24px;
}

.no-primers-needed .info-icon {
  font-size: 32px;
  margin-bottom: 8px;
}

.no-primers-needed p {
  margin: 0;
  color: #6b7280;
}

.ordering-summary {
  background: #f9fafb;
  border: 1px solid #e5e7eb;
  border-radius: 12px;
  padding: 20px;
  margin-bottom: 24px;
}

.ordering-summary h3 {
  margin: 0 0 16px;
  font-size: 16px;
}

.ordering-summary h4 {
  margin: 16px 0 12px;
  font-size: 14px;
}

.ordering-summary .summary-stats {
  display: flex;
  gap: 24px;
  margin-bottom: 16px;
}

.ordering-summary .stat {
  display: flex;
  flex-direction: column;
  align-items: center;
}

.ordering-summary .stat-value {
  font-size: 24px;
  font-weight: 600;
  color: #111827;
}

.ordering-summary .stat-label {
  font-size: 12px;
  color: #6b7280;
}

.order-table-container {
  margin-top: 16px;
}

.order-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
  margin-bottom: 12px;
}

.order-table th,
.order-table td {
  padding: 8px 12px;
  text-align: left;
  border-bottom: 1px solid #e5e7eb;
}

.order-table th {
  background: white;
  font-weight: 500;
}

.order-table .primer-seq {
  font-size: 11px;
  word-break: break-all;
  max-width: 300px;
  display: inline-block;
}

.btn-copy {
  padding: 8px 16px;
  background: #f3f4f6;
  color: #374151;
  border: 1px solid #d1d5db;
  border-radius: 6px;
  cursor: pointer;
  font-size: 13px;
  font-weight: 500;
  display: inline-flex;
  align-items: center;
  gap: 6px;
  transition: all 0.2s;
}

.btn-copy:hover {
  background: #e5e7eb;
  border-color: #9ca3af;
  color: #111827;
}

.btn-copy:active {
  background: #d1d5db;
}
`;

export default EnhancedDomesticationPanel;
