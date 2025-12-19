/**
 * EnhancedScorer - Standalone primer scoring with full analysis & visualization
 *
 * This component provides the same analysis depth as UnifiedPrimerDesigner
 * but for user-provided primers (instead of designing new ones).
 *
 * Features:
 * - Input form for forward/reverse primers and optional template
 * - Unified analysis via analyzePrimers() from primerAnalysis.js
 * - SummaryStatusPanel for quick pass/warn/fail overview
 * - EnhancedAnalysisSection for detailed thermodynamic analysis
 * - HairpinDiagram visualization for secondary structure
 * - PrimerOnTemplateViewer for template binding visualization
 */

import { useState, useCallback, FC, FormEvent } from 'react';
import { analyzePrimers } from '../lib/primerAnalysis.js';
import { checkHeterodimer, compareTmMethods } from '../lib/mutagenesis.js';
import { calculateTmQ5, calculateGC } from '../lib/tmQ5.js';
import { offTargets } from '../lib/offTargets.js';
import EnhancedAnalysisSection from './primers/EnhancedAnalysisSection';
import SummaryStatusPanel from './primers/SummaryStatusPanel';
import ScoreBreakdownPopup from './primers/ScoreBreakdownPopup';
import HairpinDiagram from './HairpinDiagram.jsx';
import PrimerOnTemplateViewer from './PrimerOnTemplateViewer.jsx';

// Type definitions
interface Thermodynamics {
  hairpinDG?: number;
  homodimerDG?: number;
  [key: string]: unknown;
}

interface Scores {
  [key: string]: number | undefined;
}

interface GQuadruplex {
  hasGQuadruplex?: boolean;
  [key: string]: unknown;
}

interface AnnealingRegion {
  sequence: string;
  length: number;
  tm: number;
  gcPercent: string;
}

interface GoldenGateDetection {
  primaryEnzyme?: string;
  [key: string]: unknown;
}

interface PrimerAnalysis {
  sequence: string;
  tm: number;
  gc: number;
  hairpinDG?: number;
  selfDimerDG?: number;
  thermodynamics?: Thermodynamics;
  scores?: Scores;
  gQuadruplex?: GQuadruplex;
  isAssemblyPrimer?: boolean;
  annealingRegion?: AnnealingRegion;
  tailRegion?: string;
  fullPrimerTm?: number;
  templateBinding?: unknown;
}

interface HeterodimerDetails {
  dg?: number;
  [key: string]: unknown;
}

interface TmComparison {
  [key: string]: unknown;
}

interface OffTargetResult {
  [key: string]: unknown;
}

interface OffTargetAnalysis {
  forward: OffTargetResult | null;
  reverse: OffTargetResult | null;
}

interface EnhancedAnalysis {
  forward: PrimerAnalysis;
  reverse: PrimerAnalysis | null;
  heterodimer: HeterodimerDetails | null;
  fwdTmComparison: TmComparison;
  revTmComparison: TmComparison | null;
  offTargets: OffTargetAnalysis | null;
  warnings: string[];
  composite: CompositeScore | null;
  quality: Quality;
  isAssemblyMode: boolean;
}

interface CompositeScore {
  score?: number;
  [key: string]: unknown;
}

interface Quality {
  tier?: string;
  label?: string;
}

interface PrimerResult {
  sequence: string;
  tm: number;
  gc: number;
  length: number;
  isAssemblyPrimer?: boolean;
  annealingRegion?: AnnealingRegion;
  fullPrimerTm?: number;
  goldenGateDetection?: GoldenGateDetection;
  thermodynamics?: Thermodynamics;
}

interface Results {
  forward: PrimerResult;
  reverse: PrimerResult | null;
  template: string | null;
  mode: string;
  analysis: EnhancedAnalysis;
  composite: CompositeScore | null;
  quality: Quality;
  effectiveScore?: number;
  criticalWarnings?: number;
}

interface CollapsedSections {
  summary: boolean;
  primers: boolean;
  analysis: boolean;
  visualization: boolean;
  hairpins: boolean;
}

// Quality tier colors
const QUALITY_COLORS: Record<string, string> = {
  excellent: '#22c55e',
  good: '#3b82f6',
  acceptable: '#eab308',
  marginal: '#f97316',
  poor: '#ef4444',
};

const EnhancedScorer: FC = () => {
  // Form state
  const [fwdPrimer, setFwdPrimer] = useState<string>('');
  const [revPrimer, setRevPrimer] = useState<string>('');
  const [template, setTemplate] = useState<string>('');
  const [showAdvanced, setShowAdvanced] = useState<boolean>(false);
  const [mode, setMode] = useState<string>('amplification');

  // Results state
  const [results, setResults] = useState<Results | null>(null);
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);

  // Score breakdown popup state
  const [showScoreBreakdown, setShowScoreBreakdown] = useState<boolean>(false);

  // Collapsible sections
  const [collapsedSections, setCollapsedSections] = useState<CollapsedSections>({
    summary: false,
    primers: false,
    analysis: false,
    visualization: false,
    hairpins: false,
  });

  const toggleSection = useCallback((sectionId: keyof CollapsedSections): void => {
    setCollapsedSections((prev) => ({
      ...prev,
      [sectionId]: !prev[sectionId],
    }));
  }, []);

  const scrollToSection = useCallback((sectionId: string): void => {
    const element = document.getElementById(`section-${sectionId}`);
    if (element) {
      element.scrollIntoView({ behavior: 'smooth', block: 'start' });
      setCollapsedSections((prev) => ({
        ...prev,
        [sectionId as keyof CollapsedSections]: false,
      }));
    }
  }, []);

  // Clean sequence input
  const cleanSequence = (seq: string): string => seq.toUpperCase().replace(/[^ATGC]/gi, '');

  // Analyze primers
  const handleAnalyze = useCallback((): void => {
    const fwdSeq = cleanSequence(fwdPrimer);
    if (!fwdSeq || fwdSeq.length < 10) {
      setError('Forward primer must be at least 10 bases');
      return;
    }

    setLoading(true);
    setError(null);
    setResults(null);

    // Use setTimeout to allow UI to update
    setTimeout(() => {
      try {
        const revSeq = revPrimer ? cleanSequence(revPrimer) : null;
        const templateSeq = template ? cleanSequence(template) : null;

        // Calculate basic properties
        const fwdTm = calculateTmQ5(fwdSeq) as number;
        const fwdGc = calculateGC(fwdSeq) as number;
        const revTm = revSeq ? (calculateTmQ5(revSeq) as number) : null;
        const revGc = revSeq ? (calculateGC(revSeq) as number) : null;

        // Run unified analysis
        const analysis = analyzePrimers(
          { seq: fwdSeq, tm: fwdTm, gc: fwdGc },
          revSeq ? { seq: revSeq, tm: revTm, gc: revGc } : undefined,
          { mode, template: templateSeq }
        ) as any;

        // Get heterodimer details if we have both primers
        let heterodimerDetails: HeterodimerDetails | null = null;
        if (revSeq) {
          heterodimerDetails = checkHeterodimer(fwdSeq, revSeq) as HeterodimerDetails;
        }

        // Get Tm comparison
        const fwdTmComparison = compareTmMethods(fwdSeq) as TmComparison;
        const revTmComparison = revSeq ? (compareTmMethods(revSeq) as TmComparison) : null;

        // Check off-targets if template provided
        let offTargetAnalysis: OffTargetAnalysis | null = null;
        if (templateSeq) {
          offTargetAnalysis = {
            forward: offTargets(fwdSeq, templateSeq) as unknown as OffTargetResult,
            reverse: revSeq ? (offTargets(revSeq, templateSeq) as unknown as OffTargetResult) : null,
          };
        }

        // Build enhanced analysis object (matching UnifiedPrimerDesigner format)
        const enhancedAnalysis: EnhancedAnalysis = {
          forward: {
            sequence: fwdSeq,
            tm: fwdTm,
            gc: fwdGc,
            hairpinDG: analysis.forward?.thermodynamics?.hairpinDG,
            selfDimerDG: analysis.forward?.thermodynamics?.homodimerDG,
            thermodynamics: analysis.forward?.thermodynamics,
            scores: analysis.forward?.scores,
            gQuadruplex: analysis.forward?.gQuadruplex,
            // Assembly-specific fields
            isAssemblyPrimer: analysis.forward?.isAssemblyPrimer,
            annealingRegion: analysis.forward?.annealingRegion,
            tailRegion: analysis.forward?.tailRegion,
            fullPrimerTm: analysis.forward?.fullPrimerTm,
            templateBinding: analysis.forward?.templateBinding,
          },
          reverse: revSeq
            ? {
                sequence: revSeq,
                tm: revTm!,
                gc: revGc!,
                hairpinDG: analysis.reverse?.thermodynamics?.hairpinDG,
                selfDimerDG: analysis.reverse?.thermodynamics?.homodimerDG,
                thermodynamics: analysis.reverse?.thermodynamics,
                scores: analysis.reverse?.scores,
                gQuadruplex: analysis.reverse?.gQuadruplex,
                // Assembly-specific fields
                isAssemblyPrimer: analysis.reverse?.isAssemblyPrimer,
                annealingRegion: analysis.reverse?.annealingRegion,
                tailRegion: analysis.reverse?.tailRegion,
                fullPrimerTm: analysis.reverse?.fullPrimerTm,
                templateBinding: analysis.reverse?.templateBinding,
              }
            : null,
          heterodimer: heterodimerDetails,
          fwdTmComparison,
          revTmComparison,
          offTargets: offTargetAnalysis,
          warnings: analysis.warnings || [],
          composite: analysis.composite,
          quality: analysis.quality,
          // Assembly mode flag
          isAssemblyMode: mode === 'assembly' && !!templateSeq,
        };

        // Use analysis results for Tm/GC when available (they may be for annealing region)
        const fwdResultTm = analysis.forward?.tm ?? fwdTm;
        const fwdResultGc = analysis.forward?.gc ?? fwdGc;
        const revResultTm = analysis.reverse?.tm ?? revTm;
        const revResultGc = analysis.reverse?.gc ?? revGc;

        setResults({
          forward: {
            sequence: fwdSeq,
            tm: Math.round(fwdResultTm * 10) / 10,
            gc: Math.round((fwdResultGc || 0) * 100),  // Guard against NaN
            length: fwdSeq.length,
            // Assembly-specific: include annealing region data
            isAssemblyPrimer: analysis.forward?.isAssemblyPrimer,
            annealingRegion: analysis.forward?.annealingRegion,
            fullPrimerTm: analysis.forward?.fullPrimerTm,
            // GG detection info
            goldenGateDetection: analysis.forward?.goldenGateDetection,
            // Thermodynamics for display
            thermodynamics: analysis.forward?.thermodynamics,
          },
          reverse: revSeq
            ? {
                sequence: revSeq,
                tm: Math.round((revResultTm || 0) * 10) / 10,
                gc: Math.round((revResultGc || 0) * 100),  // Guard against NaN
                length: revSeq.length,
                // Assembly-specific: include annealing region data
                isAssemblyPrimer: analysis.reverse?.isAssemblyPrimer,
                annealingRegion: analysis.reverse?.annealingRegion,
                fullPrimerTm: analysis.reverse?.fullPrimerTm,
                // GG detection info
                goldenGateDetection: analysis.reverse?.goldenGateDetection,
                // Thermodynamics for display
                thermodynamics: analysis.reverse?.thermodynamics,
              }
            : null,
          template: templateSeq,
          mode,
          analysis: enhancedAnalysis,
          composite: analysis.composite,
          quality: analysis.quality,
          // Include effective score with critical warning penalties
          effectiveScore: analysis.effectiveScore,
          criticalWarnings: analysis.criticalWarnings,
        });
      } catch (err) {
        console.error('Analysis error:', err);
        const errorMessage = err instanceof Error ? err.message : 'Analysis failed';
        setError(errorMessage);
      } finally {
        setLoading(false);
      }
    }, 50);
  }, [fwdPrimer, revPrimer, template, mode]);

  // Handle form submit
  const handleSubmit = (e: FormEvent<HTMLFormElement>): void => {
    e.preventDefault();
    handleAnalyze();
  };

  const getQualityTier = (): string => {
    if (!results?.quality) return 'good';
    return typeof results.quality === 'string' ? results.quality : (results.quality.tier || 'good');
  };

  const getQualityLabel = (): string => {
    if (!results?.quality) return 'Good';
    if (typeof results.quality === 'string') return results.quality;
    return results.quality.label || results.quality.tier || 'Good';
  };

  return (
    <div className="enhanced-scorer">
      {/* Input Form */}
      <div className="scorer-form-section">
        <form onSubmit={handleSubmit} className="scorer-form">
          <div className="form-group">
            <label htmlFor="fwdPrimer">
              Forward Primer <span className="required">*</span>
            </label>
            <input
              type="text"
              id="fwdPrimer"
              value={fwdPrimer}
              onChange={(e) => setFwdPrimer(e.target.value.toUpperCase().replace(/[^ATGC]/gi, ''))}
              placeholder="Enter forward primer sequence (5' to 3')..."
              required
              className="font-mono"
            />
            {fwdPrimer && (
              <small className="char-count">{cleanSequence(fwdPrimer).length} bp</small>
            )}
          </div>

          <div className="form-group">
            <label htmlFor="revPrimer">Reverse Primer (optional)</label>
            <input
              type="text"
              id="revPrimer"
              value={revPrimer}
              onChange={(e) => setRevPrimer(e.target.value.toUpperCase().replace(/[^ATGC]/gi, ''))}
              placeholder="Enter reverse primer sequence (5' to 3')..."
              className="font-mono"
            />
            {revPrimer && (
              <small className="char-count">{cleanSequence(revPrimer).length} bp</small>
            )}
          </div>

          <div className="form-group">
            <label htmlFor="template">Template Sequence (optional)</label>
            <textarea
              id="template"
              value={template}
              onChange={(e) => setTemplate(e.target.value.toUpperCase().replace(/[^ATGC]/gi, ''))}
              placeholder="Enter template sequence for off-target checking and visualization..."
              rows={3}
              className="font-mono"
            />
            {template && (
              <small className="char-count">{cleanSequence(template).length} bp</small>
            )}
          </div>

          <div className="form-group">
            <button
              type="button"
              className="toggle-advanced"
              onClick={() => setShowAdvanced(!showAdvanced)}
            >
              {showAdvanced ? 'Hide' : 'Show'} Advanced Options
            </button>
          </div>

          {showAdvanced && (
            <div className="advanced-options">
              <div className="form-group">
                <label htmlFor="mode">Analysis Mode</label>
                <select
                  id="mode"
                  value={mode}
                  onChange={(e) => setMode(e.target.value)}
                >
                  <option value="amplification">PCR Amplification</option>
                  <option value="mutagenesis">Mutagenesis (QuikChange)</option>
                  <option value="sequencing">Sanger Sequencing</option>
                  <option value="goldengate">Golden Gate Assembly</option>
                  <option value="assembly">Gibson/NEBuilder Assembly</option>
                </select>
                <small>
                  Different modes use different optimal ranges and thresholds.
                  {mode === 'goldengate' || mode === 'assembly' ? (
                    <span className="text-sky-700 block mt-1">
                      Assembly modes score the annealing region (3' template-binding portion) separately from the full primer.
                    </span>
                  ) : null}
                </small>
              </div>
            </div>
          )}

          <button type="submit" className="submit-btn" disabled={loading}>
            {loading ? 'Analyzing...' : 'Analyze Primers'}
          </button>
        </form>
      </div>

      {/* Error Display */}
      {error && (
        <div className="error-message p-3 bg-red-100 text-red-600 rounded-md mt-4">
          {error}
        </div>
      )}

      {/* Loading State */}
      {loading && (
        <div className="loading-state text-center py-10 text-slate-500">
          <div className="spinner w-10 h-10 border-4 border-slate-200 border-t-blue-500 rounded-full mx-auto mb-4" style={{
            animation: 'spin 1s linear infinite',
          }} />
          <p>Analyzing primers...</p>
        </div>
      )}

      {/* Results */}
      {results && !loading && (
        <div className="scorer-results mt-6">
          {/* Quality Badge - Clickable for score breakdown */}
          <div className="quality-header mb-5">
            <button
              type="button"
              className="quality-badge inline-flex items-center gap-3 px-5 py-3 rounded-lg cursor-pointer transition-all duration-200"
              onClick={() => setShowScoreBreakdown(true)}
              style={{
                backgroundColor: QUALITY_COLORS[getQualityTier()] + '15',
                border: `2px solid ${QUALITY_COLORS[getQualityTier()]}`,
              }}
              onMouseEnter={(e) => {
                e.currentTarget.style.transform = 'scale(1.02)';
                e.currentTarget.style.boxShadow = '0 4px 12px rgba(0,0,0,0.15)';
              }}
              onMouseLeave={(e) => {
                e.currentTarget.style.transform = 'scale(1)';
                e.currentTarget.style.boxShadow = 'none';
              }}
              title="Click to see score breakdown"
            >
              <span
                className="text-3xl font-bold"
                style={{
                  color: QUALITY_COLORS[getQualityTier()],
                }}
              >
                {results.effectiveScore ?? results.composite?.score ?? 'N/A'}
              </span>
              <div>
                <div
                  className="text-base font-semibold capitalize"
                  style={{
                    color: QUALITY_COLORS[getQualityTier()],
                  }}
                >
                  {getQualityLabel()}
                </div>
                <div className="text-xs text-slate-500 flex items-center gap-1">
                  <span>Composite Score</span>
                  {(results.criticalWarnings ?? 0) > 0 && (
                    <span className="text-[10px] text-red-500">
                      ({results.criticalWarnings} critical)
                    </span>
                  )}
                  <span className="text-[10px] opacity-70">• Click for details</span>
                </div>
              </div>
            </button>
          </div>

          {/* Score Breakdown Popup */}
          {showScoreBreakdown && (
            <ScoreBreakdownPopup
              compositeScore={results.effectiveScore ?? results.composite?.score ?? 0}
              rawScore={results.composite?.score}
              criticalWarnings={results.criticalWarnings}
              quality={results.quality}
              forwardScores={results.analysis?.forward?.scores}
              reverseScores={results.analysis?.reverse?.scores}
              hasTemplate={!!results.template}
              onClose={() => setShowScoreBreakdown(false)}
            />
          )}

          {/* Golden Gate Detection Warning - Mode Mismatch */}
          {(results.forward.goldenGateDetection || results.reverse?.goldenGateDetection) &&
           results.mode !== 'goldengate' && (
            <div className="mb-4 p-3 px-4 bg-red-50 border-2 border-red-400 rounded-lg flex items-start gap-3">
              <span className="text-xl">&#9888;</span>
              <div>
                <strong className="text-red-700 text-sm">
                  Golden Gate Primer Detected - Wrong Analysis Mode!
                </strong>
                <p className="mt-1.5 text-[13px] text-red-800">
                  <strong>{results.forward.goldenGateDetection?.primaryEnzyme || results.reverse?.goldenGateDetection?.primaryEnzyme}</strong> recognition site detected in your primer.
                  You are currently using <strong>{results.mode}</strong> mode, which scores the full primer sequence.
                </p>
                <p className="mt-2 text-[13px] text-red-800">
                  <strong>This will give inaccurate results!</strong> The Tm shown ({results.forward.tm}°C) is for the full primer,
                  not the annealing region that actually binds to your template during PCR.
                </p>
                <div className="mt-2.5 flex gap-2 flex-wrap">
                  <button
                    type="button"
                    onClick={() => {
                      setMode('goldengate');
                      setShowAdvanced(true);
                    }}
                    className="px-3 py-1.5 bg-red-600 text-white border-none rounded cursor-pointer text-xs font-semibold"
                  >
                    Switch to Golden Gate Mode
                  </button>
                  {!results.template && (
                    <span className="text-xs text-red-900 self-center">
                      + Add template sequence for annealing region analysis
                    </span>
                  )}
                </div>
              </div>
            </div>
          )}

          {/* Golden Gate Detection Reminder - In correct mode but no template */}
          {(results.forward.goldenGateDetection || results.reverse?.goldenGateDetection) &&
           results.mode === 'goldengate' && !results.template && (
            <div className="mb-4 p-3 px-4 bg-amber-100 border border-amber-300 rounded-lg flex items-start gap-3">
              <span className="text-lg">&#128161;</span>
              <div>
                <strong className="text-amber-900 text-sm">
                  Add Template for Better Analysis
                </strong>
                <p className="mt-1 text-[13px] text-amber-950">
                  {results.forward.goldenGateDetection?.primaryEnzyme || results.reverse?.goldenGateDetection?.primaryEnzyme} site detected.
                  Add a template sequence to analyze only the annealing region (the 3' portion that binds to your template).
                  Without template, the full primer Tm ({results.forward.tm}°C) is shown, which is higher than your actual annealing Tm.
                </p>
              </div>
            </div>
          )}

          {/* Summary Status Panel */}
          <div id="section-summary" className="collapsible-section mb-5">
            <button
              type="button"
              className="section-header collapsible flex justify-between items-center w-full p-3 px-4 bg-slate-50 border border-slate-200 rounded-lg cursor-pointer"
              onClick={() => toggleSection('summary')}
            >
              <h3 className="m-0 text-base">Quick Status</h3>
              <span>{collapsedSections.summary ? '▶' : '▼'}</span>
            </button>
            {!collapsedSections.summary && (
              <div className="p-4 border border-slate-200 border-t-0 rounded-b-lg">
                <SummaryStatusPanel
                  forward={results.forward}
                  reverse={results.reverse}
                  analysis={results.analysis as any}
                  quality={getQualityTier()}
                  onSectionClick={scrollToSection}
                />
              </div>
            )}
          </div>

          {/* Primer Details */}
          <div id="section-primers" className="collapsible-section mb-5">
            <button
              type="button"
              className="section-header collapsible flex justify-between items-center w-full p-3 px-4 bg-slate-50 border border-slate-200 rounded-lg cursor-pointer"
              onClick={() => toggleSection('primers')}
            >
              <h3 className="m-0 text-base">Primer Details</h3>
              <span>{collapsedSections.primers ? '▶' : '▼'}</span>
            </button>
            {!collapsedSections.primers && (
              <div className="p-4 border border-slate-200 border-t-0 rounded-b-lg">
                {/* Assembly Mode Info Banner */}
                {(results.mode === 'assembly' || results.mode === 'goldengate') && results.template && (
                  <div className="mb-4 p-3 bg-green-50 border border-green-300 rounded-md text-[13px]">
                    <strong className="text-green-600">
                      {results.mode === 'goldengate' ? 'Golden Gate' : 'Assembly'} Mode Active
                    </strong>
                    <p className="mt-1.5 mb-0 text-green-800">
                      <strong>Annealing Tm</strong> shown below is for the template-binding region only (what matters for PCR).
                      The full primer Tm is higher but not relevant for annealing temperature selection.
                      Secondary structure is checked on the full primer sequence.
                    </p>
                  </div>
                )}

                {/* Warning: Assembly mode without template */}
                {(results.mode === 'assembly' || results.mode === 'goldengate') && !results.template && (
                  <div className="mb-4 p-3 bg-amber-50 border border-amber-300 rounded-md text-[13px]">
                    <strong className="text-amber-700">No Template Provided</strong>
                    <p className="mt-1.5 mb-0 text-amber-900">
                      Without a template, the <strong>full primer Tm</strong> is shown ({results.forward.tm}°C).
                      This is <strong>higher than your actual annealing Tm</strong>!
                      Add a template sequence to see the correct annealing region Tm.
                    </p>
                  </div>
                )}

                {/* Warning: Template provided but annealing region not found */}
                {(results.mode === 'assembly' || results.mode === 'goldengate') && results.template &&
                 !results.forward.isAssemblyPrimer && (
                  <div className="mb-4 p-3 bg-red-50 border-2 border-red-400 rounded-md text-[13px]">
                    <strong className="text-red-700">Template Mismatch - Annealing Region Not Found!</strong>
                    <p className="mt-1.5 mb-0 text-red-800">
                      The primer's 3' end does not match the template sequence. This means:
                    </p>
                    <ul className="mt-2 mb-0 pl-5 text-red-800">
                      <li><strong>Tm shown ({results.forward.tm}°C)</strong> is for the full primer, NOT the annealing region</li>
                      <li>The actual annealing Tm would be significantly lower</li>
                      <li>Check that you provided the correct template sequence</li>
                      <li>Verify the primer's 3' region matches the template</li>
                    </ul>
                  </div>
                )}

                <table className="primer-table w-full border-collapse">
                  <thead>
                    <tr className="bg-slate-50">
                      <th className="p-2.5 text-left border-b-2 border-slate-200">Direction</th>
                      <th className="p-2.5 text-left border-b-2 border-slate-200">Sequence (5' → 3')</th>
                      <th className="p-2.5 text-center border-b-2 border-slate-200">Length</th>
                      <th className="p-2.5 text-center border-b-2 border-slate-200">
                        {results.forward.isAssemblyPrimer && results.forward.annealingRegion ? (
                          <span title="Annealing region Tm - the temperature at which the primer-binding portion anneals to template">
                            Annealing Tm (°C)
                          </span>
                        ) : (
                          <span title="Melting temperature of the full primer sequence">
                            Tm (°C)
                          </span>
                        )}
                      </th>
                      <th className="p-2.5 text-center border-b-2 border-slate-200">GC%</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      <td className="p-2.5 font-semibold text-blue-600">Forward</td>
                      <td className="p-2.5 font-mono text-[13px] break-all">
                        {results.forward.isAssemblyPrimer && results.forward.annealingRegion ? (
                          <span>
                            <span className="text-slate-400">
                              {results.forward.sequence.slice(0, results.forward.sequence.length - results.forward.annealingRegion.length)}
                            </span>
                            <span className="text-blue-600 font-semibold">
                              {results.forward.annealingRegion.sequence}
                            </span>
                          </span>
                        ) : (
                          results.forward.sequence
                        )}
                      </td>
                      <td className="p-2.5 text-center">{results.forward.length} bp</td>
                      <td className="p-2.5 text-center font-semibold">
                        {results.forward.isAssemblyPrimer && results.forward.annealingRegion ? (
                          <span title={`Full primer Tm: ${results.forward.fullPrimerTm}°C`}>
                            {results.forward.tm}
                            <span className="text-[10px] text-slate-500 block">
                              (annealing)
                            </span>
                          </span>
                        ) : (
                          results.forward.tm
                        )}
                      </td>
                      <td className="p-2.5 text-center">{results.forward.gc}%</td>
                    </tr>
                    {results.reverse && (
                      <tr>
                        <td className="p-2.5 font-semibold text-red-600">Reverse</td>
                        <td className="p-2.5 font-mono text-[13px] break-all">
                          {results.reverse.isAssemblyPrimer && results.reverse.annealingRegion ? (
                            <span>
                              <span className="text-slate-400">
                                {results.reverse.sequence.slice(0, results.reverse.sequence.length - results.reverse.annealingRegion.length)}
                              </span>
                              <span className="text-red-600 font-semibold">
                                {results.reverse.annealingRegion.sequence}
                              </span>
                            </span>
                          ) : (
                            results.reverse.sequence
                          )}
                        </td>
                        <td className="p-2.5 text-center">{results.reverse.length} bp</td>
                        <td className="p-2.5 text-center font-semibold">
                          {results.reverse.isAssemblyPrimer && results.reverse.annealingRegion ? (
                            <span title={`Full primer Tm: ${results.reverse.fullPrimerTm}°C`}>
                              {results.reverse.tm}
                              <span className="text-[10px] text-slate-500 block">
                                (annealing)
                              </span>
                            </span>
                          ) : (
                            results.reverse.tm
                          )}
                        </td>
                        <td className="p-2.5 text-center">{results.reverse.gc}%</td>
                      </tr>
                    )}
                  </tbody>
                </table>

                {/* Assembly Primer Annealing Details */}
                {results.forward.isAssemblyPrimer && results.forward.annealingRegion && (
                  <div className="mt-4 p-3 bg-blue-50 rounded-md border border-blue-200">
                    <h4 className="m-0 mb-2.5 text-sm text-blue-800">
                      Forward Primer - Annealing Region Analysis
                    </h4>
                    <div className="grid grid-cols-[repeat(auto-fit,minmax(140px,1fr))] gap-2.5">
                      <div>
                        <span className="text-[11px] text-slate-500">Annealing Sequence</span>
                        <div className="font-mono text-xs text-blue-800 font-semibold">
                          {results.forward.annealingRegion.sequence}
                        </div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Annealing Length</span>
                        <div className="font-semibold">{results.forward.annealingRegion.length} bp</div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Annealing Tm</span>
                        <div className="font-semibold">{results.forward.annealingRegion.tm}°C</div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Annealing GC</span>
                        <div className="font-semibold">{results.forward.annealingRegion.gcPercent}</div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Full Primer Tm</span>
                        <div className="font-semibold text-slate-500">{results.forward.fullPrimerTm}°C</div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Tail Length</span>
                        <div className="font-semibold text-slate-500">
                          {results.forward.sequence.length - results.forward.annealingRegion.length} bp
                        </div>
                      </div>
                    </div>
                  </div>
                )}

                {results.reverse?.isAssemblyPrimer && results.reverse?.annealingRegion && (
                  <div className="mt-3 p-3 bg-red-50 rounded-md border border-red-200">
                    <h4 className="m-0 mb-2.5 text-sm text-red-800">
                      Reverse Primer - Annealing Region Analysis
                    </h4>
                    <div className="grid grid-cols-[repeat(auto-fit,minmax(140px,1fr))] gap-2.5">
                      <div>
                        <span className="text-[11px] text-slate-500">Annealing Sequence</span>
                        <div className="font-mono text-xs text-red-800 font-semibold">
                          {results.reverse.annealingRegion.sequence}
                        </div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Annealing Length</span>
                        <div className="font-semibold">{results.reverse.annealingRegion.length} bp</div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Annealing Tm</span>
                        <div className="font-semibold">{results.reverse.annealingRegion.tm}°C</div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Annealing GC</span>
                        <div className="font-semibold">{results.reverse.annealingRegion.gcPercent}</div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Full Primer Tm</span>
                        <div className="font-semibold text-slate-500">{results.reverse.fullPrimerTm}°C</div>
                      </div>
                      <div>
                        <span className="text-[11px] text-slate-500">Tail Length</span>
                        <div className="font-semibold text-slate-500">
                          {results.reverse.sequence.length - results.reverse.annealingRegion.length} bp
                        </div>
                      </div>
                    </div>
                  </div>
                )}

                {results.reverse && (
                  <div className="mt-3 p-2.5 bg-slate-50 rounded-md">
                    <strong>Tm Difference:</strong>{' '}
                    <span className="font-semibold" style={{
                      color: Math.abs(results.forward.tm - results.reverse.tm) <= 2 ? '#22c55e' :
                             Math.abs(results.forward.tm - results.reverse.tm) <= 5 ? '#eab308' : '#ef4444'
                    }}>
                      {Math.abs(results.forward.tm - results.reverse.tm).toFixed(1)}°C
                    </span>
                    {Math.abs(results.forward.tm - results.reverse.tm) <= 2 && ' (ideal)'}
                    {Math.abs(results.forward.tm - results.reverse.tm) > 2 && Math.abs(results.forward.tm - results.reverse.tm) <= 5 && ' (acceptable)'}
                    {Math.abs(results.forward.tm - results.reverse.tm) > 5 && ' (may cause issues)'}
                    {results.forward.isAssemblyPrimer && (
                      <span className="text-[11px] text-slate-500 ml-2">
                        (annealing region Tm comparison)
                      </span>
                    )}
                  </div>
                )}
              </div>
            )}
          </div>

          {/* Enhanced Analysis Section */}
          <div className="mb-5">
            <EnhancedAnalysisSection
              analysis={results.analysis as any}
              collapsible={true}
              defaultCollapsed={collapsedSections.analysis}
            />
          </div>

          {/* Hairpin Visualization */}
          <div id="section-hairpins" className="collapsible-section mb-5">
            <button
              type="button"
              className="section-header collapsible flex justify-between items-center w-full p-3 px-4 bg-slate-50 border border-slate-200 rounded-lg cursor-pointer"
              onClick={() => toggleSection('hairpins')}
            >
              <h3 className="m-0 text-base">Secondary Structure Visualization</h3>
              <span>{collapsedSections.hairpins ? '▶' : '▼'}</span>
            </button>
            {!collapsedSections.hairpins && (
              <div className="p-4 border border-slate-200 border-t-0 rounded-b-lg grid gap-4" style={{
                gridTemplateColumns: results.reverse ? '1fr 1fr' : '1fr'
              }}>
                <HairpinDiagram
                  sequence={results.forward.sequence}
                  primerName="Forward Primer"
                  width={380}
                  showDetails={true}
                />
                {results.reverse && (
                  <HairpinDiagram
                    sequence={results.reverse.sequence}
                    primerName="Reverse Primer"
                    width={380}
                    showDetails={true}
                  />
                )}
              </div>
            )}
          </div>

          {/* Template Visualization */}
          {results.template && (
            <div id="section-visualization" className="collapsible-section mb-5">
              <button
                type="button"
                className="section-header collapsible flex justify-between items-center w-full p-3 px-4 bg-slate-50 border border-slate-200 rounded-lg cursor-pointer"
                onClick={() => toggleSection('visualization')}
              >
                <h3 className="m-0 text-base">Template Binding Visualization</h3>
                <span>{collapsedSections.visualization ? '▶' : '▼'}</span>
              </button>
              {!collapsedSections.visualization && (
                <div className="p-4 border border-slate-200 border-t-0 rounded-b-lg">
                  <PrimerOnTemplateViewer
                    template={results.template}
                    forward={{
                      sequence: results.forward.sequence,
                      tm: results.forward.tm,
                      gc: results.forward.gc / 100,  // Convert from percentage to fraction
                      gcPercent: `${results.forward.gc}%`,
                      dg: results.forward.thermodynamics?.hairpinDG,
                    }}
                    reverse={results.reverse ? {
                      sequence: results.reverse.sequence,
                      tm: results.reverse.tm,
                      gc: results.reverse.gc / 100,  // Convert from percentage to fraction
                      gcPercent: `${results.reverse.gc}%`,
                      dg: results.reverse.thermodynamics?.hairpinDG,
                    } : null}
                    showHairpinDiagrams={true}
                    showSequenceDetails={true}
                    width={800}
                  />
                </div>
              )}
            </div>
          )}

          {/* Copy Buttons */}
          <div className="copy-section mt-5 p-4 bg-slate-50 rounded-lg flex gap-3 flex-wrap">
            <button
              type="button"
              onClick={() => navigator.clipboard.writeText(results.forward.sequence)}
              className="px-4 py-2 bg-blue-600 text-white border-none rounded-md cursor-pointer text-[13px]"
            >
              Copy Forward
            </button>
            {results.reverse && (
              <button
                type="button"
                onClick={() => navigator.clipboard.writeText(results.reverse!.sequence)}
                className="px-4 py-2 bg-red-600 text-white border-none rounded-md cursor-pointer text-[13px]"
              >
                Copy Reverse
              </button>
            )}
            {results.reverse && (
              <button
                type="button"
                onClick={() => navigator.clipboard.writeText(
                  `Forward: ${results.forward.sequence}\nReverse: ${results.reverse!.sequence}`
                )}
                className="px-4 py-2 bg-slate-500 text-white border-none rounded-md cursor-pointer text-[13px]"
              >
                Copy Both
              </button>
            )}
          </div>
        </div>
      )}

      <style>{`
        @keyframes spin {
          to { transform: rotate(360deg); }
        }
        .enhanced-scorer .form-group {
          margin-bottom: 16px;
        }
        .enhanced-scorer .form-group label {
          display: block;
          margin-bottom: 6px;
          font-weight: 500;
          color: #374151;
        }
        .enhanced-scorer .form-group .required {
          color: #ef4444;
        }
        .enhanced-scorer .form-group input,
        .enhanced-scorer .form-group textarea,
        .enhanced-scorer .form-group select {
          width: 100%;
          padding: 10px 12px;
          border: 1px solid #d1d5db;
          border-radius: 6px;
          font-size: 14px;
        }
        .enhanced-scorer .form-group input:focus,
        .enhanced-scorer .form-group textarea:focus,
        .enhanced-scorer .form-group select:focus {
          outline: none;
          border-color: #3b82f6;
          box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }
        .enhanced-scorer .char-count {
          display: block;
          margin-top: 4px;
          color: #6b7280;
          font-size: 12px;
        }
        .enhanced-scorer .toggle-advanced {
          background: none;
          border: none;
          color: #3b82f6;
          cursor: pointer;
          font-size: 13px;
          padding: 0;
        }
        .enhanced-scorer .toggle-advanced:hover {
          text-decoration: underline;
        }
        .enhanced-scorer .advanced-options {
          padding: 16px;
          background: #f8fafc;
          border-radius: 8px;
          margin-bottom: 16px;
        }
        .enhanced-scorer .submit-btn {
          width: 100%;
          padding: 12px;
          background: #3b82f6;
          color: white;
          border: none;
          border-radius: 8px;
          font-size: 16px;
          font-weight: 600;
          cursor: pointer;
          transition: background 0.2s;
        }
        .enhanced-scorer .submit-btn:hover:not(:disabled) {
          background: #2563eb;
        }
        .enhanced-scorer .submit-btn:disabled {
          background: #9ca3af;
          cursor: not-allowed;
        }
      `}</style>
    </div>
  );
};

export default EnhancedScorer;
