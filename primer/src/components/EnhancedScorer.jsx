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

import React, { useState, useMemo, useCallback } from 'react';
import { analyzePrimers } from '../lib/primerAnalysis.js';
import { checkHeterodimer, compareTmMethods, DIMER_THRESHOLDS } from '../lib/mutagenesis.js';
import { calculateTmQ5, calculateGC } from '../lib/tmQ5.js';
import { offTargets } from '../lib/offTargets.js';
import EnhancedAnalysisSection from './primers/EnhancedAnalysisSection';
import SummaryStatusPanel from './primers/SummaryStatusPanel';
import ScoreBreakdownPopup from './primers/ScoreBreakdownPopup.jsx';
import HairpinDiagram from './HairpinDiagram.jsx';
import PrimerOnTemplateViewer from './PrimerOnTemplateViewer.jsx';

// Quality tier colors
const QUALITY_COLORS = {
  excellent: '#22c55e',
  good: '#3b82f6',
  acceptable: '#eab308',
  marginal: '#f97316',
  poor: '#ef4444',
};

export default function EnhancedScorer() {
  // Form state
  const [fwdPrimer, setFwdPrimer] = useState('');
  const [revPrimer, setRevPrimer] = useState('');
  const [template, setTemplate] = useState('');
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [mode, setMode] = useState('amplification');

  // Results state
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  // Score breakdown popup state
  const [showScoreBreakdown, setShowScoreBreakdown] = useState(false);

  // Collapsible sections
  const [collapsedSections, setCollapsedSections] = useState({
    summary: false,
    primers: false,
    analysis: false,
    visualization: false,
    hairpins: false,
  });

  const toggleSection = useCallback((sectionId) => {
    setCollapsedSections((prev) => ({
      ...prev,
      [sectionId]: !prev[sectionId],
    }));
  }, []);

  const scrollToSection = useCallback((sectionId) => {
    const element = document.getElementById(`section-${sectionId}`);
    if (element) {
      element.scrollIntoView({ behavior: 'smooth', block: 'start' });
      setCollapsedSections((prev) => ({
        ...prev,
        [sectionId]: false,
      }));
    }
  }, []);

  // Clean sequence input
  const cleanSequence = (seq) => seq.toUpperCase().replace(/[^ATGC]/gi, '');

  // Analyze primers
  const handleAnalyze = useCallback(() => {
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
        const fwdTm = calculateTmQ5(fwdSeq);
        const fwdGc = calculateGC(fwdSeq);
        const revTm = revSeq ? calculateTmQ5(revSeq) : null;
        const revGc = revSeq ? calculateGC(revSeq) : null;

        // Run unified analysis
        const analysis = analyzePrimers(
          { seq: fwdSeq, tm: fwdTm, gc: fwdGc },
          revSeq ? { seq: revSeq, tm: revTm, gc: revGc } : null,
          { mode, template: templateSeq }
        );

        // Get heterodimer details if we have both primers
        let heterodimerDetails = null;
        if (revSeq) {
          heterodimerDetails = checkHeterodimer(fwdSeq, revSeq);
        }

        // Get Tm comparison
        const fwdTmComparison = compareTmMethods(fwdSeq);
        const revTmComparison = revSeq ? compareTmMethods(revSeq) : null;

        // Check off-targets if template provided
        let offTargetAnalysis = null;
        if (templateSeq) {
          offTargetAnalysis = {
            forward: offTargets(fwdSeq, templateSeq),
            reverse: revSeq ? offTargets(revSeq, templateSeq) : null,
          };
        }

        // Build enhanced analysis object (matching UnifiedPrimerDesigner format)
        const enhancedAnalysis = {
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
                tm: revTm,
                gc: revGc,
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
                tm: Math.round(revResultTm * 10) / 10,
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
        setError(err.message || 'Analysis failed');
      } finally {
        setLoading(false);
      }
    }, 50);
  }, [fwdPrimer, revPrimer, template, mode]);

  // Handle form submit
  const handleSubmit = (e) => {
    e.preventDefault();
    handleAnalyze();
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
              style={{ fontFamily: 'monospace' }}
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
              style={{ fontFamily: 'monospace' }}
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
              style={{ fontFamily: 'monospace' }}
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
                    <span style={{ color: '#0369a1', display: 'block', marginTop: '4px' }}>
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
        <div className="error-message" style={{
          padding: '12px',
          backgroundColor: '#fee2e2',
          color: '#dc2626',
          borderRadius: '6px',
          marginTop: '16px'
        }}>
          {error}
        </div>
      )}

      {/* Loading State */}
      {loading && (
        <div className="loading-state" style={{
          textAlign: 'center',
          padding: '40px',
          color: '#64748b'
        }}>
          <div className="spinner" style={{
            width: '40px',
            height: '40px',
            border: '3px solid #e2e8f0',
            borderTopColor: '#3b82f6',
            borderRadius: '50%',
            animation: 'spin 1s linear infinite',
            margin: '0 auto 16px',
          }} />
          <p>Analyzing primers...</p>
        </div>
      )}

      {/* Results */}
      {results && !loading && (
        <div className="scorer-results" style={{ marginTop: '24px' }}>
          {/* Quality Badge - Clickable for score breakdown */}
          <div className="quality-header" style={{ marginBottom: '20px' }}>
            <button
              type="button"
              className="quality-badge"
              onClick={() => setShowScoreBreakdown(true)}
              style={{
                display: 'inline-flex',
                alignItems: 'center',
                gap: '12px',
                padding: '12px 20px',
                backgroundColor: QUALITY_COLORS[results.quality?.tier || results.quality] + '15',
                border: `2px solid ${QUALITY_COLORS[results.quality?.tier || results.quality]}`,
                borderRadius: '8px',
                cursor: 'pointer',
                transition: 'all 0.2s ease',
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
                style={{
                  fontSize: '28px',
                  fontWeight: 'bold',
                  color: QUALITY_COLORS[results.quality?.tier || results.quality],
                }}
              >
                {results.effectiveScore ?? results.composite?.score ?? 'N/A'}
              </span>
              <div>
                <div
                  style={{
                    fontSize: '16px',
                    fontWeight: '600',
                    color: QUALITY_COLORS[results.quality?.tier || results.quality],
                    textTransform: 'capitalize',
                  }}
                >
                  {results.quality?.label || results.quality?.tier || results.quality}
                </div>
                <div style={{ fontSize: '12px', color: '#64748b', display: 'flex', alignItems: 'center', gap: '4px' }}>
                  <span>Composite Score</span>
                  {results.criticalWarnings > 0 && (
                    <span style={{ fontSize: '10px', color: '#ef4444' }}>
                      ({results.criticalWarnings} critical)
                    </span>
                  )}
                  <span style={{ fontSize: '10px', opacity: 0.7 }}>• Click for details</span>
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
            <div style={{
              marginBottom: '16px',
              padding: '12px 16px',
              backgroundColor: '#fef2f2',
              border: '2px solid #f87171',
              borderRadius: '8px',
              display: 'flex',
              alignItems: 'flex-start',
              gap: '12px',
            }}>
              <span style={{ fontSize: '20px' }}>&#9888;</span>
              <div>
                <strong style={{ color: '#b91c1c', fontSize: '14px' }}>
                  Golden Gate Primer Detected - Wrong Analysis Mode!
                </strong>
                <p style={{ margin: '6px 0 0', fontSize: '13px', color: '#991b1b' }}>
                  <strong>{results.forward.goldenGateDetection?.primaryEnzyme || results.reverse?.goldenGateDetection?.primaryEnzyme}</strong> recognition site detected in your primer.
                  You are currently using <strong>{results.mode}</strong> mode, which scores the full primer sequence.
                </p>
                <p style={{ margin: '8px 0 0', fontSize: '13px', color: '#991b1b' }}>
                  <strong>This will give inaccurate results!</strong> The Tm shown ({results.forward.tm}°C) is for the full primer,
                  not the annealing region that actually binds to your template during PCR.
                </p>
                <div style={{ marginTop: '10px', display: 'flex', gap: '8px', flexWrap: 'wrap' }}>
                  <button
                    type="button"
                    onClick={() => {
                      setMode('goldengate');
                      setShowAdvanced(true);
                    }}
                    style={{
                      padding: '6px 12px',
                      backgroundColor: '#dc2626',
                      color: 'white',
                      border: 'none',
                      borderRadius: '4px',
                      cursor: 'pointer',
                      fontSize: '12px',
                      fontWeight: '600',
                    }}
                  >
                    Switch to Golden Gate Mode
                  </button>
                  {!results.template && (
                    <span style={{ fontSize: '12px', color: '#7f1d1d', alignSelf: 'center' }}>
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
            <div style={{
              marginBottom: '16px',
              padding: '12px 16px',
              backgroundColor: '#fef3c7',
              border: '1px solid #fcd34d',
              borderRadius: '8px',
              display: 'flex',
              alignItems: 'flex-start',
              gap: '12px',
            }}>
              <span style={{ fontSize: '18px' }}>&#128161;</span>
              <div>
                <strong style={{ color: '#92400e', fontSize: '14px' }}>
                  Add Template for Better Analysis
                </strong>
                <p style={{ margin: '4px 0 0', fontSize: '13px', color: '#78350f' }}>
                  {results.forward.goldenGateDetection?.primaryEnzyme || results.reverse?.goldenGateDetection?.primaryEnzyme} site detected.
                  Add a template sequence to analyze only the annealing region (the 3' portion that binds to your template).
                  Without template, the full primer Tm ({results.forward.tm}°C) is shown, which is higher than your actual annealing Tm.
                </p>
              </div>
            </div>
          )}

          {/* Summary Status Panel */}
          <div id="section-summary" className="collapsible-section" style={{ marginBottom: '20px' }}>
            <button
              type="button"
              className="section-header collapsible"
              onClick={() => toggleSection('summary')}
              style={{
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                width: '100%',
                padding: '12px 16px',
                backgroundColor: '#f8fafc',
                border: '1px solid #e2e8f0',
                borderRadius: '8px',
                cursor: 'pointer',
              }}
            >
              <h3 style={{ margin: 0, fontSize: '16px' }}>Quick Status</h3>
              <span>{collapsedSections.summary ? '▶' : '▼'}</span>
            </button>
            {!collapsedSections.summary && (
              <div style={{ padding: '16px', border: '1px solid #e2e8f0', borderTop: 'none', borderRadius: '0 0 8px 8px' }}>
                <SummaryStatusPanel
                  forward={results.forward}
                  reverse={results.reverse}
                  analysis={results.analysis}
                  quality={results.quality?.tier || results.quality}
                  onSectionClick={scrollToSection}
                />
              </div>
            )}
          </div>

          {/* Primer Details */}
          <div id="section-primers" className="collapsible-section" style={{ marginBottom: '20px' }}>
            <button
              type="button"
              className="section-header collapsible"
              onClick={() => toggleSection('primers')}
              style={{
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                width: '100%',
                padding: '12px 16px',
                backgroundColor: '#f8fafc',
                border: '1px solid #e2e8f0',
                borderRadius: '8px',
                cursor: 'pointer',
              }}
            >
              <h3 style={{ margin: 0, fontSize: '16px' }}>Primer Details</h3>
              <span>{collapsedSections.primers ? '▶' : '▼'}</span>
            </button>
            {!collapsedSections.primers && (
              <div style={{ padding: '16px', border: '1px solid #e2e8f0', borderTop: 'none', borderRadius: '0 0 8px 8px' }}>
                {/* Assembly Mode Info Banner */}
                {(results.mode === 'assembly' || results.mode === 'goldengate') && results.template && (
                  <div style={{
                    marginBottom: '16px',
                    padding: '12px',
                    backgroundColor: '#f0fdf4',
                    border: '1px solid #bbf7d0',
                    borderRadius: '6px',
                    fontSize: '13px',
                  }}>
                    <strong style={{ color: '#16a34a' }}>
                      {results.mode === 'goldengate' ? 'Golden Gate' : 'Assembly'} Mode Active
                    </strong>
                    <p style={{ margin: '6px 0 0', color: '#15803d' }}>
                      <strong>Annealing Tm</strong> shown below is for the template-binding region only (what matters for PCR).
                      The full primer Tm is higher but not relevant for annealing temperature selection.
                      Secondary structure is checked on the full primer sequence.
                    </p>
                  </div>
                )}

                {/* Warning: Assembly mode without template */}
                {(results.mode === 'assembly' || results.mode === 'goldengate') && !results.template && (
                  <div style={{
                    marginBottom: '16px',
                    padding: '12px',
                    backgroundColor: '#fffbeb',
                    border: '1px solid #fcd34d',
                    borderRadius: '6px',
                    fontSize: '13px',
                  }}>
                    <strong style={{ color: '#b45309' }}>No Template Provided</strong>
                    <p style={{ margin: '6px 0 0', color: '#92400e' }}>
                      Without a template, the <strong>full primer Tm</strong> is shown ({results.forward.tm}°C).
                      This is <strong>higher than your actual annealing Tm</strong>!
                      Add a template sequence to see the correct annealing region Tm.
                    </p>
                  </div>
                )}

                {/* Warning: Template provided but annealing region not found */}
                {(results.mode === 'assembly' || results.mode === 'goldengate') && results.template &&
                 !results.forward.isAssemblyPrimer && (
                  <div style={{
                    marginBottom: '16px',
                    padding: '12px',
                    backgroundColor: '#fef2f2',
                    border: '2px solid #f87171',
                    borderRadius: '6px',
                    fontSize: '13px',
                  }}>
                    <strong style={{ color: '#b91c1c' }}>Template Mismatch - Annealing Region Not Found!</strong>
                    <p style={{ margin: '6px 0 0', color: '#991b1b' }}>
                      The primer's 3' end does not match the template sequence. This means:
                    </p>
                    <ul style={{ margin: '8px 0 0', paddingLeft: '20px', color: '#991b1b' }}>
                      <li><strong>Tm shown ({results.forward.tm}°C)</strong> is for the full primer, NOT the annealing region</li>
                      <li>The actual annealing Tm would be significantly lower</li>
                      <li>Check that you provided the correct template sequence</li>
                      <li>Verify the primer's 3' region matches the template</li>
                    </ul>
                  </div>
                )}

                <table className="primer-table" style={{ width: '100%', borderCollapse: 'collapse' }}>
                  <thead>
                    <tr style={{ backgroundColor: '#f8fafc' }}>
                      <th style={{ padding: '10px', textAlign: 'left', borderBottom: '2px solid #e2e8f0' }}>Direction</th>
                      <th style={{ padding: '10px', textAlign: 'left', borderBottom: '2px solid #e2e8f0' }}>Sequence (5' → 3')</th>
                      <th style={{ padding: '10px', textAlign: 'center', borderBottom: '2px solid #e2e8f0' }}>Length</th>
                      <th style={{ padding: '10px', textAlign: 'center', borderBottom: '2px solid #e2e8f0' }}>
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
                      <th style={{ padding: '10px', textAlign: 'center', borderBottom: '2px solid #e2e8f0' }}>GC%</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      <td style={{ padding: '10px', fontWeight: '600', color: '#2563eb' }}>Forward</td>
                      <td style={{ padding: '10px', fontFamily: 'monospace', fontSize: '13px', wordBreak: 'break-all' }}>
                        {results.forward.isAssemblyPrimer && results.forward.annealingRegion ? (
                          <span>
                            <span style={{ color: '#9ca3af' }}>
                              {results.forward.sequence.slice(0, results.forward.sequence.length - results.forward.annealingRegion.length)}
                            </span>
                            <span style={{ color: '#2563eb', fontWeight: '600' }}>
                              {results.forward.annealingRegion.sequence}
                            </span>
                          </span>
                        ) : (
                          results.forward.sequence
                        )}
                      </td>
                      <td style={{ padding: '10px', textAlign: 'center' }}>{results.forward.length} bp</td>
                      <td style={{ padding: '10px', textAlign: 'center', fontWeight: '600' }}>
                        {results.forward.isAssemblyPrimer && results.forward.annealingRegion ? (
                          <span title={`Full primer Tm: ${results.forward.fullPrimerTm}°C`}>
                            {results.forward.tm}
                            <span style={{ fontSize: '10px', color: '#6b7280', display: 'block' }}>
                              (annealing)
                            </span>
                          </span>
                        ) : (
                          results.forward.tm
                        )}
                      </td>
                      <td style={{ padding: '10px', textAlign: 'center' }}>{results.forward.gc}%</td>
                    </tr>
                    {results.reverse && (
                      <tr>
                        <td style={{ padding: '10px', fontWeight: '600', color: '#dc2626' }}>Reverse</td>
                        <td style={{ padding: '10px', fontFamily: 'monospace', fontSize: '13px', wordBreak: 'break-all' }}>
                          {results.reverse.isAssemblyPrimer && results.reverse.annealingRegion ? (
                            <span>
                              <span style={{ color: '#9ca3af' }}>
                                {results.reverse.sequence.slice(0, results.reverse.sequence.length - results.reverse.annealingRegion.length)}
                              </span>
                              <span style={{ color: '#dc2626', fontWeight: '600' }}>
                                {results.reverse.annealingRegion.sequence}
                              </span>
                            </span>
                          ) : (
                            results.reverse.sequence
                          )}
                        </td>
                        <td style={{ padding: '10px', textAlign: 'center' }}>{results.reverse.length} bp</td>
                        <td style={{ padding: '10px', textAlign: 'center', fontWeight: '600' }}>
                          {results.reverse.isAssemblyPrimer && results.reverse.annealingRegion ? (
                            <span title={`Full primer Tm: ${results.reverse.fullPrimerTm}°C`}>
                              {results.reverse.tm}
                              <span style={{ fontSize: '10px', color: '#6b7280', display: 'block' }}>
                                (annealing)
                              </span>
                            </span>
                          ) : (
                            results.reverse.tm
                          )}
                        </td>
                        <td style={{ padding: '10px', textAlign: 'center' }}>{results.reverse.gc}%</td>
                      </tr>
                    )}
                  </tbody>
                </table>

                {/* Assembly Primer Annealing Details */}
                {results.forward.isAssemblyPrimer && results.forward.annealingRegion && (
                  <div style={{
                    marginTop: '16px',
                    padding: '12px',
                    backgroundColor: '#eff6ff',
                    borderRadius: '6px',
                    border: '1px solid #bfdbfe',
                  }}>
                    <h4 style={{ margin: '0 0 10px', fontSize: '14px', color: '#1e40af' }}>
                      Forward Primer - Annealing Region Analysis
                    </h4>
                    <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(140px, 1fr))', gap: '10px' }}>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Annealing Sequence</span>
                        <div style={{ fontFamily: 'monospace', fontSize: '12px', color: '#1e40af', fontWeight: '600' }}>
                          {results.forward.annealingRegion.sequence}
                        </div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Annealing Length</span>
                        <div style={{ fontWeight: '600' }}>{results.forward.annealingRegion.length} bp</div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Annealing Tm</span>
                        <div style={{ fontWeight: '600' }}>{results.forward.annealingRegion.tm}°C</div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Annealing GC</span>
                        <div style={{ fontWeight: '600' }}>{results.forward.annealingRegion.gcPercent}</div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Full Primer Tm</span>
                        <div style={{ fontWeight: '600', color: '#6b7280' }}>{results.forward.fullPrimerTm}°C</div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Tail Length</span>
                        <div style={{ fontWeight: '600', color: '#6b7280' }}>
                          {results.forward.sequence.length - results.forward.annealingRegion.length} bp
                        </div>
                      </div>
                    </div>
                  </div>
                )}

                {results.reverse?.isAssemblyPrimer && results.reverse?.annealingRegion && (
                  <div style={{
                    marginTop: '12px',
                    padding: '12px',
                    backgroundColor: '#fef2f2',
                    borderRadius: '6px',
                    border: '1px solid #fecaca',
                  }}>
                    <h4 style={{ margin: '0 0 10px', fontSize: '14px', color: '#991b1b' }}>
                      Reverse Primer - Annealing Region Analysis
                    </h4>
                    <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(140px, 1fr))', gap: '10px' }}>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Annealing Sequence</span>
                        <div style={{ fontFamily: 'monospace', fontSize: '12px', color: '#991b1b', fontWeight: '600' }}>
                          {results.reverse.annealingRegion.sequence}
                        </div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Annealing Length</span>
                        <div style={{ fontWeight: '600' }}>{results.reverse.annealingRegion.length} bp</div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Annealing Tm</span>
                        <div style={{ fontWeight: '600' }}>{results.reverse.annealingRegion.tm}°C</div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Annealing GC</span>
                        <div style={{ fontWeight: '600' }}>{results.reverse.annealingRegion.gcPercent}</div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Full Primer Tm</span>
                        <div style={{ fontWeight: '600', color: '#6b7280' }}>{results.reverse.fullPrimerTm}°C</div>
                      </div>
                      <div>
                        <span style={{ fontSize: '11px', color: '#6b7280' }}>Tail Length</span>
                        <div style={{ fontWeight: '600', color: '#6b7280' }}>
                          {results.reverse.sequence.length - results.reverse.annealingRegion.length} bp
                        </div>
                      </div>
                    </div>
                  </div>
                )}

                {results.reverse && (
                  <div style={{ marginTop: '12px', padding: '10px', backgroundColor: '#f8fafc', borderRadius: '6px' }}>
                    <strong>Tm Difference:</strong>{' '}
                    <span style={{
                      color: Math.abs(results.forward.tm - results.reverse.tm) <= 2 ? '#22c55e' :
                             Math.abs(results.forward.tm - results.reverse.tm) <= 5 ? '#eab308' : '#ef4444',
                      fontWeight: '600'
                    }}>
                      {Math.abs(results.forward.tm - results.reverse.tm).toFixed(1)}°C
                    </span>
                    {Math.abs(results.forward.tm - results.reverse.tm) <= 2 && ' (ideal)'}
                    {Math.abs(results.forward.tm - results.reverse.tm) > 2 && Math.abs(results.forward.tm - results.reverse.tm) <= 5 && ' (acceptable)'}
                    {Math.abs(results.forward.tm - results.reverse.tm) > 5 && ' (may cause issues)'}
                    {results.forward.isAssemblyPrimer && (
                      <span style={{ fontSize: '11px', color: '#6b7280', marginLeft: '8px' }}>
                        (annealing region Tm comparison)
                      </span>
                    )}
                  </div>
                )}
              </div>
            )}
          </div>

          {/* Enhanced Analysis Section */}
          <div style={{ marginBottom: '20px' }}>
            <EnhancedAnalysisSection
              analysis={results.analysis}
              collapsible={true}
              defaultCollapsed={collapsedSections.analysis}
            />
          </div>

          {/* Hairpin Visualization */}
          <div id="section-hairpins" className="collapsible-section" style={{ marginBottom: '20px' }}>
            <button
              type="button"
              className="section-header collapsible"
              onClick={() => toggleSection('hairpins')}
              style={{
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                width: '100%',
                padding: '12px 16px',
                backgroundColor: '#f8fafc',
                border: '1px solid #e2e8f0',
                borderRadius: '8px',
                cursor: 'pointer',
              }}
            >
              <h3 style={{ margin: 0, fontSize: '16px' }}>Secondary Structure Visualization</h3>
              <span>{collapsedSections.hairpins ? '▶' : '▼'}</span>
            </button>
            {!collapsedSections.hairpins && (
              <div style={{
                padding: '16px',
                border: '1px solid #e2e8f0',
                borderTop: 'none',
                borderRadius: '0 0 8px 8px',
                display: 'grid',
                gridTemplateColumns: results.reverse ? '1fr 1fr' : '1fr',
                gap: '16px'
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
            <div id="section-visualization" className="collapsible-section" style={{ marginBottom: '20px' }}>
              <button
                type="button"
                className="section-header collapsible"
                onClick={() => toggleSection('visualization')}
                style={{
                  display: 'flex',
                  justifyContent: 'space-between',
                  alignItems: 'center',
                  width: '100%',
                  padding: '12px 16px',
                  backgroundColor: '#f8fafc',
                  border: '1px solid #e2e8f0',
                  borderRadius: '8px',
                  cursor: 'pointer',
                }}
              >
                <h3 style={{ margin: 0, fontSize: '16px' }}>Template Binding Visualization</h3>
                <span>{collapsedSections.visualization ? '▶' : '▼'}</span>
              </button>
              {!collapsedSections.visualization && (
                <div style={{ padding: '16px', border: '1px solid #e2e8f0', borderTop: 'none', borderRadius: '0 0 8px 8px' }}>
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
          <div className="copy-section" style={{
            marginTop: '20px',
            padding: '16px',
            backgroundColor: '#f8fafc',
            borderRadius: '8px',
            display: 'flex',
            gap: '12px',
            flexWrap: 'wrap'
          }}>
            <button
              type="button"
              onClick={() => navigator.clipboard.writeText(results.forward.sequence)}
              style={{
                padding: '8px 16px',
                backgroundColor: '#2563eb',
                color: 'white',
                border: 'none',
                borderRadius: '6px',
                cursor: 'pointer',
                fontSize: '13px',
              }}
            >
              Copy Forward
            </button>
            {results.reverse && (
              <button
                type="button"
                onClick={() => navigator.clipboard.writeText(results.reverse.sequence)}
                style={{
                  padding: '8px 16px',
                  backgroundColor: '#dc2626',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '13px',
                }}
              >
                Copy Reverse
              </button>
            )}
            {results.reverse && (
              <button
                type="button"
                onClick={() => navigator.clipboard.writeText(
                  `Forward: ${results.forward.sequence}\nReverse: ${results.reverse.sequence}`
                )}
                style={{
                  padding: '8px 16px',
                  backgroundColor: '#6b7280',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '13px',
                }}
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
}
