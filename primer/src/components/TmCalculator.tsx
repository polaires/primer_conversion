import { useState, useCallback, useMemo, CSSProperties } from 'react';
import {
  calculateTm,
  calculateAnnealing,
  calculateGC,
  getPolymeraseList,
  getPolymeraseInfo,
} from '../lib/polymerases.js';

// Type definitions
interface Polymerase {
  id: string;
  name: string;
  type: string;
  color: string;
}

interface ProtocolStep {
  temp: number;
  time: string;
}

interface PolymeraseProtocol {
  initialDenature: ProtocolStep;
  denature: ProtocolStep;
  extension: ProtocolStep;
  finalExtension: ProtocolStep;
  cycles: number;
}

interface PolymeraseInfo {
  id: string;
  name: string;
  type: string;
  color: string;
  fidelity: string;
  extension: string;
  hotStart: boolean;
  proofReading: boolean;
  method: string;
  notes: string;
  defaults: {
    primerConc: number;
    maxAnneal: number;
  };
  protocol: PolymeraseProtocol;
}

interface PrimerStats {
  length: number;
  gc: string;
}

interface TmQuality {
  status: 'optimal' | 'acceptable' | 'warning';
  message: string;
}

interface PrimerResult {
  sequence: string;
  length: number;
  gc: number;
  tm: number;
}

interface CalculationResult {
  mode: 'single' | 'pair';
  primer1: PrimerResult;
  primer2?: PrimerResult;
  annealingTemp?: number;
  tmDifference?: number;
  isCapped?: boolean;
  method?: string;
}

// Nucleotide colors for sequence display
const NT_COLORS: Record<string, string> = {
  A: '#22c55e',
  T: '#ef4444',
  G: '#f59e0b',
  C: '#3b82f6',
};

interface ColorizedSequenceInputProps {
  sequence: string;
  stats: PrimerStats;
}

/**
 * Colorized sequence display component - integrated into input
 */
function ColorizedSequenceInput({ sequence, stats }: ColorizedSequenceInputProps) {
  if (!sequence) return null;

  return (
    <div className="tm-seq-input-preview">
      <div className="tm-seq-colored">
        {sequence.split('').map((nt, i) => (
          <span key={i} style={{ color: NT_COLORS[nt] || '#6b7280' }}>{nt}</span>
        ))}
      </div>
      {stats && (
        <div className="tm-seq-stats-inline">
          <span className="tm-stat-chip">{stats.length} bp</span>
          <span className="tm-stat-chip gc">GC {stats.gc}%</span>
        </div>
      )}
    </div>
  );
}

interface TemperatureGaugeProps {
  value: number;
  min?: number;
  max?: number;
  label: string;
  quality?: TmQuality;
}

/**
 * Visual temperature gauge component
 */
function TemperatureGauge({ value, min = 40, max = 85, label, quality }: TemperatureGaugeProps) {
  const percentage = Math.max(0, Math.min(100, ((value - min) / (max - min)) * 100));
  const optimalStart = ((55 - min) / (max - min)) * 100;
  const optimalEnd = ((72 - min) / (max - min)) * 100;

  return (
    <div className={`tm-gauge ${quality?.status || ''}`}>
      <div className="tm-gauge-label">{label}</div>
      <div className="tm-gauge-value">{value}°C</div>
      <div className="tm-gauge-bar">
        <div className="tm-gauge-optimal" style={{
          left: `${optimalStart}%`,
          width: `${optimalEnd - optimalStart}%`
        }} />
        <div className="tm-gauge-indicator" style={{ left: `${percentage}%` }}>
          <div className="tm-gauge-needle" />
        </div>
      </div>
      <div className="tm-gauge-scale">
        <span>{min}°C</span>
        <span className="optimal-label">Optimal</span>
        <span>{max}°C</span>
      </div>
      {quality && (
        <div className={`tm-gauge-quality ${quality.status}`}>
          {quality.message}
        </div>
      )}
    </div>
  );
}

interface PolymeraseSelectorProps {
  selected: string;
  onChange: (id: string) => void;
  showMethod: boolean;
  onToggleMethod: () => void;
}

/**
 * Collapsible Polymerase selector with integrated method info
 */
function PolymeraseSelector({ selected, onChange, showMethod, onToggleMethod }: PolymeraseSelectorProps) {
  const polymerases = getPolymeraseList() as Polymerase[];
  const [expanded, setExpanded] = useState<boolean>(false);
  const selectedInfo = getPolymeraseInfo(selected) as PolymeraseInfo | null;

  return (
    <div className={`tm-polymerase-selector ${expanded ? 'expanded' : 'collapsed'}`}>
      {/* Compact header showing current selection */}
      <button
        className="tm-poly-compact-header"
        onClick={() => setExpanded(!expanded)}
      >
        <div className="tm-poly-current">
          <span className="tm-poly-indicator" style={{ backgroundColor: selectedInfo?.color }} />
          <span className="tm-poly-current-name">{selectedInfo?.name}</span>
          <span className="tm-poly-current-type">{selectedInfo?.type}</span>
        </div>
        <div className="tm-poly-expand-icon">
          <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
            <path d={expanded ? "M7.41 15.41L12 10.83l4.59 4.58L18 14l-6-6-6 6z" : "M7.41 8.59L12 13.17l4.59-4.58L18 10l-6 6-6-6z"} />
          </svg>
        </div>
      </button>

      {/* Expanded polymerase selection */}
      {expanded && (
        <div className="tm-poly-expanded-content">
          <div className="tm-polymerase-grid">
            {polymerases.map(p => (
              <button
                key={p.id}
                className={`tm-polymerase-btn ${selected === p.id ? 'active' : ''}`}
                onClick={() => {
                  onChange(p.id);
                  setExpanded(false);
                }}
                style={{ '--poly-color': p.color } as CSSProperties}
              >
                <span className="tm-poly-indicator" style={{ backgroundColor: p.color }} />
                <span className="tm-poly-name">{p.name}</span>
                <span className="tm-poly-type">{p.type}</span>
              </button>
            ))}
          </div>

          {/* Polymerase details */}
          {selectedInfo && (
            <div className="tm-polymerase-details">
              <div className="tm-poly-specs-row">
                <div className="tm-poly-spec-item">
                  <span className="tm-spec-label">Fidelity</span>
                  <span className="tm-spec-value">{selectedInfo.fidelity}</span>
                </div>
                <div className="tm-poly-spec-item">
                  <span className="tm-spec-label">Extension</span>
                  <span className="tm-spec-value">{selectedInfo.extension}</span>
                </div>
                <div className="tm-poly-spec-item">
                  <span className="tm-spec-label">Hot Start</span>
                  <span className={`tm-spec-value ${selectedInfo.hotStart ? 'yes' : 'no'}`}>
                    {selectedInfo.hotStart ? 'Yes' : 'No'}
                  </span>
                </div>
                <div className="tm-poly-spec-item">
                  <span className="tm-spec-label">Proofreading</span>
                  <span className={`tm-spec-value ${selectedInfo.proofReading ? 'yes' : 'no'}`}>
                    {selectedInfo.proofReading ? 'Yes' : 'No'}
                  </span>
                </div>
              </div>

              {/* Method toggle moved here */}
              <div className="tm-method-toggle-row">
                <label className="tm-method-toggle-label">
                  <div className={`tm-toggle ${showMethod ? 'active' : ''}`} onClick={(e) => { e.stopPropagation(); onToggleMethod(); }}>
                    <div className="tm-toggle-slider" />
                  </div>
                  <span>Show calculation method in results</span>
                </label>
              </div>

              {showMethod && (
                <div className="tm-poly-method-preview">
                  <span className="tm-method-label">Tm Method:</span>
                  <span className="tm-method-value">{selectedInfo.method}</span>
                </div>
              )}

              <p className="tm-poly-notes">{selectedInfo.notes}</p>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

/**
 * Multi-Polymerase Tm Calculator Component - Modern UI
 */
export default function TmCalculator() {
  const [primer1, setPrimer1] = useState<string>('');
  const [primer2, setPrimer2] = useState<string>('');
  const [polymerase, setPolymerase] = useState<string>('q5'); // Default to Q5
  const [primerConc, setPrimerConc] = useState<number>(500);
  const [showMethod, setShowMethod] = useState<boolean>(true); // Show by default
  const [result, setResult] = useState<CalculationResult | null>(null);
  const [error, setError] = useState<string | null>(null);

  const polymeraseInfo = useMemo(() => getPolymeraseInfo(polymerase) as PolymeraseInfo | null, [polymerase]);

  // Update primer concentration when polymerase changes
  const handlePolymeraseChange = useCallback((newPolymerase: string): void => {
    setPolymerase(newPolymerase);
    const info = getPolymeraseInfo(newPolymerase) as PolymeraseInfo | null;
    if (info) {
      setPrimerConc(info.defaults.primerConc);
    }
    setResult(null);
  }, []);

  // Clean sequence input
  const cleanSequence = useCallback((seq: string): string => {
    return seq.toUpperCase().replace(/[^ATGCN]/g, '').replace(/N/g, '');
  }, []);

  // Live cleaned sequences
  const cleanedPrimer1 = useMemo(() => cleanSequence(primer1), [primer1, cleanSequence]);
  const cleanedPrimer2 = useMemo(() => cleanSequence(primer2), [primer2, cleanSequence]);

  // Calculate result when primers change
  const calculate = useCallback((): void => {
    setError(null);
    setResult(null);

    const seq1 = cleanedPrimer1;
    const seq2 = cleanedPrimer2;

    if (!seq1) {
      setError('Please enter at least one primer sequence');
      return;
    }

    try {
      if (seq1 && seq2) {
        // Two primers - calculate annealing temperature
        const annealResult = calculateAnnealing(seq1, seq2, polymerase, { primerConc }) as {
          tm1: number;
          tm2: number;
          annealingTemp: number;
          tmDifference: number;
          isCapped: boolean;
          method: string;
        };

        setResult({
          mode: 'pair',
          primer1: {
            sequence: seq1,
            length: seq1.length,
            gc: calculateGC(seq1) as number,
            tm: annealResult.tm1,
          },
          primer2: {
            sequence: seq2,
            length: seq2.length,
            gc: calculateGC(seq2) as number,
            tm: annealResult.tm2,
          },
          annealingTemp: annealResult.annealingTemp,
          tmDifference: annealResult.tmDifference,
          isCapped: annealResult.isCapped,
          method: annealResult.method,
        });
      } else {
        // Single primer - calculate Tm only
        const tm = calculateTm(seq1, polymerase, { primerConc }) as number;

        setResult({
          mode: 'single',
          primer1: {
            sequence: seq1,
            length: seq1.length,
            gc: calculateGC(seq1) as number,
            tm,
          },
        });
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'An error occurred');
    }
  }, [cleanedPrimer1, cleanedPrimer2, polymerase, primerConc]);

  // Load example primers
  const loadExample = useCallback((): void => {
    setPrimer1('ACGTACGTACGTACGTACGT');
    setPrimer2('TGCATGCATGCATGCATGCA');
  }, []);

  // Clear all
  const clearAll = useCallback((): void => {
    setPrimer1('');
    setPrimer2('');
    setResult(null);
    setError(null);
  }, []);

  // Get Tm quality assessment
  const getTmQuality = useCallback((tm: number): TmQuality => {
    if (tm >= 55 && tm <= 72) return { status: 'optimal', message: 'Optimal range' };
    if (tm >= 50 && tm < 55) return { status: 'acceptable', message: 'Acceptable (low)' };
    if (tm > 72 && tm <= 80) return { status: 'acceptable', message: 'Acceptable (high)' };
    if (tm < 50) return { status: 'warning', message: 'Too low' };
    return { status: 'warning', message: 'Too high' };
  }, []);

  // Get concentration options based on polymerase
  const concentrationOptions = useMemo((): number[] => {
    if (polymeraseInfo) {
      const defaultConc = polymeraseInfo.defaults.primerConc;
      const options = [100, 200, 300, 500, 1000].filter(c =>
        c === defaultConc || [200, 500, 1000].includes(c)
      );
      if (!options.includes(defaultConc)) {
        options.push(defaultConc);
        options.sort((a, b) => a - b);
      }
      return options;
    }
    return [200, 500, 1000];
  }, [polymeraseInfo]);

  // Theme color from polymerase
  const themeColor = polymeraseInfo?.color || '#0d9488';

  return (
    <div className="tm-calculator-modern" style={{ '--tm-theme-color': themeColor } as CSSProperties}>
      {/* Header */}
      <div className="tm-header-modern">
        <div className="tm-header-icon">
          <svg viewBox="0 0 24 24" width="32" height="32" fill="none" stroke="currentColor" strokeWidth="2">
            <path d="M14 4v10.54a4 4 0 1 1-4 0V4a2 2 0 0 1 4 0Z" />
            <path d="M12 14a1 1 0 1 0 0 2 1 1 0 0 0 0-2z" fill="currentColor" />
          </svg>
        </div>
        <div className="tm-header-text">
          <h2>Tm Calculator</h2>
          <p>Calculate melting and annealing temperatures for PCR primers</p>
        </div>
      </div>

      {/* Collapsible Polymerase Selector */}
      <PolymeraseSelector
        selected={polymerase}
        onChange={handlePolymeraseChange}
        showMethod={showMethod}
        onToggleMethod={() => setShowMethod(!showMethod)}
      />

      {/* Input Section */}
      <div className="tm-input-section-modern">
        {/* Forward Primer */}
        <div className="tm-input-card">
          <div className="tm-input-header">
            <div className="tm-input-badge fwd">FWD</div>
            <label>Forward Primer (5' → 3')</label>
          </div>
          <div className="tm-input-with-preview">
            <textarea
              value={primer1}
              onChange={(e) => setPrimer1(e.target.value)}
              placeholder="Enter primer sequence..."
              rows={1}
              spellCheck={false}
              className={cleanedPrimer1 ? 'has-content' : ''}
            />
            {cleanedPrimer1 && (
              <ColorizedSequenceInput
                sequence={cleanedPrimer1}
                stats={{
                  length: cleanedPrimer1.length,
                  gc: (calculateGC(cleanedPrimer1) as number * 100).toFixed(0)
                }}
              />
            )}
          </div>
        </div>

        {/* Reverse Primer */}
        <div className="tm-input-card">
          <div className="tm-input-header">
            <div className="tm-input-badge rev">REV</div>
            <label>Reverse Primer (5' → 3')</label>
            <span className="tm-optional-tag">Optional</span>
          </div>
          <div className="tm-input-with-preview">
            <textarea
              value={primer2}
              onChange={(e) => setPrimer2(e.target.value)}
              placeholder="Enter for annealing temp..."
              rows={1}
              spellCheck={false}
              className={cleanedPrimer2 ? 'has-content' : ''}
            />
            {cleanedPrimer2 && (
              <ColorizedSequenceInput
                sequence={cleanedPrimer2}
                stats={{
                  length: cleanedPrimer2.length,
                  gc: (calculateGC(cleanedPrimer2) as number * 100).toFixed(0)
                }}
              />
            )}
          </div>
        </div>

        {/* Concentration Options */}
        <div className="tm-options-modern">
          <div className="tm-option-group">
            <label>Primer Concentration</label>
            <div className="tm-conc-pills">
              {concentrationOptions.map(conc => (
                <button
                  key={conc}
                  className={`tm-conc-pill ${primerConc === conc ? 'active' : ''}`}
                  onClick={() => setPrimerConc(conc)}
                >
                  {conc} nM
                </button>
              ))}
            </div>
          </div>
        </div>

        {/* Actions */}
        <div className="tm-actions-modern">
          <button className="tm-btn-calculate" onClick={calculate}>
            <svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" strokeWidth="2">
              <path d="M14 4v10.54a4 4 0 1 1-4 0V4a2 2 0 0 1 4 0Z" />
            </svg>
            Calculate Tm
          </button>
          <button className="tm-btn-secondary" onClick={loadExample}>
            Example
          </button>
          <button className="tm-btn-secondary" onClick={clearAll}>
            Clear
          </button>
        </div>
      </div>

      {/* Error Display */}
      {error && (
        <div className="tm-error-modern">
          <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
            <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z"/>
          </svg>
          <span>{error}</span>
        </div>
      )}

      {/* Results */}
      {result && (
        <div className="tm-results-modern">
          {/* Method Badge - shown by default */}
          {showMethod && (
            <div className="tm-method-badge">
              <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                <path d="M9.4 16.6L4.8 12l4.6-4.6L8 6l-6 6 6 6 1.4-1.4zm5.2 0l4.6-4.6-4.6-4.6L16 6l6 6-6 6-1.4-1.4z"/>
              </svg>
              <span>{polymeraseInfo?.method || 'Nearest-neighbor thermodynamics'}</span>
            </div>
          )}

          {/* Annealing Temperature Hero */}
          {result.mode === 'pair' && result.annealingTemp && (
            <div className="tm-annealing-hero">
              <div className="tm-annealing-content">
                <div className="tm-annealing-label">Recommended Annealing Temperature</div>
                <div className="tm-annealing-value">
                  {result.annealingTemp}<span className="tm-degree">°C</span>
                </div>
                {result.isCapped && (
                  <div className="tm-capped-badge">Capped at {polymeraseInfo?.defaults.maxAnneal || 72}°C max</div>
                )}
              </div>
              <div className="tm-annealing-visual">
                <div className="tm-thermometer">
                  <div className="tm-thermometer-fill" style={{
                    height: `${Math.min(100, ((result.annealingTemp - 40) / 45) * 100)}%`
                  }} />
                  <div className="tm-thermometer-mark" style={{ bottom: '33%' }}>55°C</div>
                  <div className="tm-thermometer-mark optimal" style={{ bottom: '71%' }}>72°C</div>
                </div>
              </div>
            </div>
          )}

          {/* Primer Results - Compact without redundant sequence display */}
          <div className="tm-primers-grid">
            {/* Primer 1 */}
            <div className="tm-primer-card">
              <div className="tm-primer-header">
                <div className={`tm-primer-badge ${result.mode === 'pair' ? 'fwd' : ''}`}>
                  {result.mode === 'pair' ? 'Forward' : 'Primer'}
                </div>
                <div className="tm-primer-quick-stats">
                  <span>{result.primer1.length} bp</span>
                  <span>GC {(result.primer1.gc * 100).toFixed(0)}%</span>
                </div>
              </div>

              <TemperatureGauge
                value={result.primer1.tm}
                label="Melting Temperature"
                quality={getTmQuality(result.primer1.tm)}
              />
            </div>

            {/* Primer 2 */}
            {result.mode === 'pair' && result.primer2 && (
              <div className="tm-primer-card">
                <div className="tm-primer-header">
                  <div className="tm-primer-badge rev">Reverse</div>
                  <div className="tm-primer-quick-stats">
                    <span>{result.primer2.length} bp</span>
                    <span>GC {(result.primer2.gc * 100).toFixed(0)}%</span>
                  </div>
                </div>

                <TemperatureGauge
                  value={result.primer2.tm}
                  label="Melting Temperature"
                  quality={getTmQuality(result.primer2.tm)}
                />
              </div>
            )}
          </div>

          {/* Tm Difference Warning */}
          {result.mode === 'pair' && result.tmDifference && result.tmDifference > 5 && (
            <div className="tm-warning-modern">
              <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
                <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
              </svg>
              <div>
                <strong>Tm Difference: {result.tmDifference}°C</strong>
                <p>Consider redesigning primers for Tm within 5°C of each other.</p>
              </div>
            </div>
          )}

          {/* Protocol Timeline */}
          {polymeraseInfo && (
            <div className="tm-protocol-modern">
              <h4>
                <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
                  <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm0 16H5V5h14v14z"/>
                  <path d="M7 12h2v5H7zm4-3h2v8h-2zm4-3h2v11h-2z"/>
                </svg>
                {polymeraseInfo.name} PCR Protocol
              </h4>
              <div className="tm-protocol-timeline">
                <div className="tm-protocol-step">
                  <div className="tm-step-marker">1</div>
                  <div className="tm-step-content">
                    <div className="tm-step-label">Initial Denaturation</div>
                    <div className="tm-step-value">
                      {polymeraseInfo.protocol.initialDenature.temp}°C
                      <span className="tm-step-time">{polymeraseInfo.protocol.initialDenature.time}</span>
                    </div>
                  </div>
                </div>
                <div className="tm-protocol-step cycling">
                  <div className="tm-step-marker cycle">
                    <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                      <path d="M12 4V1L8 5l4 4V6c3.31 0 6 2.69 6 6 0 1.01-.25 1.97-.7 2.8l1.46 1.46C19.54 15.03 20 13.57 20 12c0-4.42-3.58-8-8-8zm0 14c-3.31 0-6-2.69-6-6 0-1.01.25-1.97.7-2.8L5.24 7.74C4.46 8.97 4 10.43 4 12c0 4.42 3.58 8 8 8v3l4-4-4-4v3z"/>
                    </svg>
                  </div>
                  <div className="tm-step-content">
                    <div className="tm-step-label">Cycling <span className="tm-cycle-count">{polymeraseInfo.protocol.cycles}×</span></div>
                    <div className="tm-cycle-steps">
                      <div className="tm-cycle-step">
                        <span className="tm-cycle-temp denature">{polymeraseInfo.protocol.denature.temp}°C</span>
                        <span className="tm-cycle-time">{polymeraseInfo.protocol.denature.time}</span>
                        <span className="tm-cycle-desc">Denature</span>
                      </div>
                      <div className="tm-cycle-step highlight">
                        <span className="tm-cycle-temp anneal">{result.annealingTemp || result.primer1.tm}°C</span>
                        <span className="tm-cycle-time">20-30 sec</span>
                        <span className="tm-cycle-desc">Anneal</span>
                      </div>
                      <div className="tm-cycle-step">
                        <span className="tm-cycle-temp extend">{polymeraseInfo.protocol.extension.temp}°C</span>
                        <span className="tm-cycle-time">{polymeraseInfo.protocol.extension.time}</span>
                        <span className="tm-cycle-desc">Extend</span>
                      </div>
                    </div>
                  </div>
                </div>
                <div className="tm-protocol-step">
                  <div className="tm-step-marker">3</div>
                  <div className="tm-step-content">
                    <div className="tm-step-label">Final Extension</div>
                    <div className="tm-step-value">
                      {polymeraseInfo.protocol.finalExtension.temp}°C
                      <span className="tm-step-time">{polymeraseInfo.protocol.finalExtension.time}</span>
                    </div>
                  </div>
                </div>
              </div>
            </div>
          )}
        </div>
      )}

      {/* Info Footer - Minimal */}
      <div className="tm-info-minimal">
        <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
          <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-6h2v6zm0-8h-2V7h2v2z"/>
        </svg>
        <span>Supports {(getPolymeraseList() as Polymerase[]).length} polymerases with manufacturer-specific Tm calculations</span>
      </div>
    </div>
  );
}
