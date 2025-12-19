import { useState, useCallback, useEffect, FC, ChangeEvent } from 'react';
import {
  planAssembly,
  parseFasta,
  digest,
  DEFAULT_CONFIG,
  PART_SOURCES,
  getAvailablePartSources,
  fetchPartDatabase,
  loadCachedDatabase,
  saveCachedDatabase,
  mergeDatabases,
} from '../lib/repp/index.js';

// Type definitions
interface Primer {
  seq: string;
  tm?: number;
  gc?: number;
}

interface Fragment {
  id: string;
  seq: string;
  pcrSeq?: string;
  type: string;
  cost?: number;
  url?: string;
  source?: string;
  primers?: Primer[];
  start?: number;
  end?: number;
}

interface PartSource {
  id: string;
  name: string;
  color: string;
  icon: string;
  cost: number;
}

interface Solution {
  count: number;
  cost: number;
  fragments: Fragment[];
}

interface AssemblyResult {
  target: string;
  seq: string;
  solutions: Solution[];
}

interface Config {
  fragmentsMaxCount: number;
  fragmentsMinHomology: number;
  fragmentsMaxHomology: number;
  pcrMinLength: number;
  pcrBpCost: number;
  pcrRxnCost: number;
  gibsonAssemblyCost: number;
  [key: string]: number;
}

interface FragmentDetailProps {
  frag: Fragment;
  index: number;
  onCopy?: (message: string) => void;
}

interface AssemblyVisualizationProps {
  fragments: Fragment[];
  targetLength: number;
}

interface CostBreakdownProps {
  solution: Solution;
  config: Config;
}

interface PartSourceSelectorProps {
  sources: PartSource[];
  enabledSources: string[];
  onToggle: (sourceId: string) => void;
  loadingStates: Record<string, boolean>;
  loadedCounts: Record<string, number>;
  onLoad: (sourceId: string) => void;
}

const COMMON_ENZYMES: string[] = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'XbaI', 'SpeI', 'PstI', 'SalI', 'NotI', 'NcoI', 'NdeI', 'BglII', 'KpnI', 'SacI', 'BsaI', 'BsmBI'];

// Fragment detail component with expandable sequence view
const FragmentDetail: FC<FragmentDetailProps> = ({ frag, index, onCopy }) => {
  const [expanded, setExpanded] = useState<boolean>(false);
  const [seqView, setSeqView] = useState<'full' | 'start' | 'end'>('full');

  const formatSequence = (seq: string, chunkSize: number = 60): string[] => {
    if (!seq) return [];
    const chunks: string[] = [];
    for (let i = 0; i < seq.length; i += chunkSize) {
      chunks.push(seq.slice(i, i + chunkSize));
    }
    return chunks;
  };

  const getDisplaySeq = (): string => {
    const seq = frag.pcrSeq || frag.seq;
    if (!seq) return '';
    if (seqView === 'start') return seq.slice(0, 120);
    if (seqView === 'end') return seq.slice(-120);
    return seq;
  };

  const copyToClipboard = (text: string, label: string): void => {
    navigator.clipboard.writeText(text).then(() => {
      onCopy && onCopy(`Copied ${label}`);
    });
  };

  const fullSeq = frag.pcrSeq || frag.seq;
  const partSourceInfo = frag.source ? (PART_SOURCES as Record<string, PartSource>)[frag.source] : undefined;

  return (
    <div className={`fragment-detail ${expanded ? 'expanded' : ''}`}>
      <div className="fragment-header" onClick={() => setExpanded(!expanded)}>
        <div className="fragment-info">
          <span className="fragment-number">{index + 1}</span>
          <span className="fragment-id">
            {frag.url ? (
              <a
                href={frag.url}
                target="_blank"
                rel="noopener noreferrer"
                onClick={(e) => e.stopPropagation()}
                className="fragment-link"
              >
                {frag.id}
              </a>
            ) : frag.id}
          </span>
          <span className={`type-badge ${frag.type}`}>{frag.type}</span>
          {partSourceInfo && (
            <span
              className="source-badge"
              style={{ backgroundColor: partSourceInfo.color }}
              title={partSourceInfo.name}
            >
              {partSourceInfo.icon}
            </span>
          )}
        </div>
        <div className="fragment-meta">
          <span className="fragment-length">{fullSeq.length} bp</span>
          <span className="fragment-cost">${(frag.cost || 0).toFixed(2)}</span>
          <span className="expand-icon">{expanded ? '▼' : '▶'}</span>
        </div>
      </div>

      {expanded && (
        <div className="fragment-body">
          {/* Sequence Section */}
          <div className="sequence-section">
            <div className="section-header">
              <h5>Sequence</h5>
              <div className="seq-controls">
                <button
                  className={seqView === 'start' ? 'active' : ''}
                  onClick={() => setSeqView('start')}
                >
                  First 120bp
                </button>
                <button
                  className={seqView === 'end' ? 'active' : ''}
                  onClick={() => setSeqView('end')}
                >
                  Last 120bp
                </button>
                <button
                  className={seqView === 'full' ? 'active' : ''}
                  onClick={() => setSeqView('full')}
                >
                  Full
                </button>
                <button
                  className="copy-btn"
                  onClick={() => copyToClipboard(fullSeq, 'sequence')}
                >
                  Copy
                </button>
              </div>
            </div>
            <div className="sequence-display">
              {formatSequence(getDisplaySeq()).map((chunk, i) => (
                <div key={i} className="seq-line">
                  <span className="seq-pos">{(seqView === 'end' ? fullSeq.length - 120 : 0) + i * 60 + 1}</span>
                  <span className="seq-text">{chunk}</span>
                </div>
              ))}
              {seqView !== 'full' && fullSeq.length > 120 && (
                <div className="seq-truncated">... {fullSeq.length - 120} bp not shown</div>
              )}
            </div>
          </div>

          {/* Primers Section */}
          {frag.primers && frag.primers.length >= 2 && (
            <div className="primers-detail-section">
              <h5>PCR Primers</h5>
              <div className="primer-detail">
                <div className="primer-row">
                  <span className="primer-dir fwd">FWD</span>
                  <div className="primer-info">
                    <code className="primer-seq">{frag.primers[0].seq}</code>
                    <div className="primer-stats">
                      <span>Length: {frag.primers[0].seq.length} bp</span>
                      <span>Tm: {frag.primers[0].tm?.toFixed(1) || 'N/A'}°C</span>
                      <span>GC: {frag.primers[0].gc?.toFixed(1) || 'N/A'}%</span>
                    </div>
                  </div>
                  <button
                    className="copy-btn"
                    onClick={() => copyToClipboard(frag.primers![0].seq, 'FWD primer')}
                  >
                    Copy
                  </button>
                </div>
                <div className="primer-row">
                  <span className="primer-dir rev">REV</span>
                  <div className="primer-info">
                    <code className="primer-seq">{frag.primers[1].seq}</code>
                    <div className="primer-stats">
                      <span>Length: {frag.primers[1].seq.length} bp</span>
                      <span>Tm: {frag.primers[1].tm?.toFixed(1) || 'N/A'}°C</span>
                      <span>GC: {frag.primers[1].gc?.toFixed(1) || 'N/A'}%</span>
                    </div>
                  </div>
                  <button
                    className="copy-btn"
                    onClick={() => copyToClipboard(frag.primers![1].seq, 'REV primer')}
                  >
                    Copy
                  </button>
                </div>
              </div>
            </div>
          )}

          {/* Fragment position info */}
          {(frag.start !== undefined || frag.end !== undefined) && (
            <div className="position-info">
              <span>Position: {frag.start} - {frag.end}</span>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

// Assembly visualization component
const AssemblyVisualization: FC<AssemblyVisualizationProps> = ({ fragments, targetLength }) => {
  return (
    <div className="assembly-viz">
      <div className="viz-track target-track">
        <span className="track-label">Target</span>
        <div className="track-bar">
          <div className="target-segment" style={{ width: '100%' }}></div>
        </div>
        <span className="track-length">{targetLength} bp</span>
      </div>
      <div className="viz-track fragments-track">
        <span className="track-label">Fragments</span>
        <div className="track-bar">
          {fragments.map((frag, i) => {
            const width = ((frag.pcrSeq || frag.seq).length / targetLength) * 100;
            const colors = ['#3b82f6', '#8b5cf6', '#10b981', '#f59e0b', '#ef4444', '#06b6d4'];
            return (
              <div
                key={i}
                className={`fragment-segment ${frag.type}`}
                style={{
                  width: `${Math.max(width, 2)}%`,
                  backgroundColor: colors[i % colors.length],
                }}
                title={`${frag.id}: ${(frag.pcrSeq || frag.seq).length} bp`}
              >
                <span className="segment-label">{i + 1}</span>
              </div>
            );
          })}
        </div>
      </div>
    </div>
  );
};

// Cost breakdown component
const CostBreakdown: FC<CostBreakdownProps> = ({ solution, config }) => {
  const pcrFragments = solution.fragments.filter(f => f.type === 'pcr');
  const synthFragments = solution.fragments.filter(f => f.type === 'synthetic');

  const primerCost = pcrFragments.reduce((sum, f) => {
    if (!f.primers) return sum;
    return sum + (f.primers[0].seq.length + f.primers[1].seq.length) * config.pcrBpCost;
  }, 0);

  const pcrRxnCost = pcrFragments.length * config.pcrRxnCost;
  const synthCost = synthFragments.reduce((sum, f) => sum + (f.cost || 0), 0);
  const gibsonCost = solution.fragments.length > 1 ? config.gibsonAssemblyCost : 0;

  return (
    <div className="cost-breakdown">
      <h5>Cost Breakdown</h5>
      <div className="cost-items">
        {primerCost > 0 && (
          <div className="cost-item">
            <span>Primers ({pcrFragments.length} pairs)</span>
            <span>${primerCost.toFixed(2)}</span>
          </div>
        )}
        {pcrRxnCost > 0 && (
          <div className="cost-item">
            <span>PCR Reactions</span>
            <span>${pcrRxnCost.toFixed(2)}</span>
          </div>
        )}
        {synthCost > 0 && (
          <div className="cost-item">
            <span>Synthesis ({synthFragments.length} fragments)</span>
            <span>${synthCost.toFixed(2)}</span>
          </div>
        )}
        {gibsonCost > 0 && (
          <div className="cost-item">
            <span>Gibson Assembly</span>
            <span>${gibsonCost.toFixed(2)}</span>
          </div>
        )}
        <div className="cost-item total">
          <span>Total</span>
          <span>${solution.cost.toFixed(2)}</span>
        </div>
      </div>
    </div>
  );
};

// Part Source Selector Component
const PartSourceSelector: FC<PartSourceSelectorProps> = ({ sources, enabledSources, onToggle, loadingStates, loadedCounts, onLoad }) => {
  return (
    <div className="part-sources">
      {sources.map(source => (
        <div key={source.id} className="part-source-item">
          <label className="part-source-toggle">
            <input
              type="checkbox"
              checked={enabledSources.includes(source.id)}
              onChange={() => onToggle(source.id)}
              disabled={loadingStates[source.id]}
            />
            <span
              className="source-icon"
              style={{ backgroundColor: source.color }}
            >
              {source.icon}
            </span>
            <span className="source-name">{source.name}</span>
          </label>
          <div className="source-info">
            <span className="source-cost">${source.cost.toFixed(0)}/part</span>
            {loadedCounts[source.id] > 0 && (
              <span className="source-count">{loadedCounts[source.id].toLocaleString()} parts</span>
            )}
            {loadingStates[source.id] && (
              <span className="source-loading">Loading...</span>
            )}
          </div>
          {enabledSources.includes(source.id) && loadedCounts[source.id] === 0 && !loadingStates[source.id] && (
            <button
              className="load-btn"
              onClick={() => onLoad(source.id)}
            >
              Load Database
            </button>
          )}
        </div>
      ))}
      <p className="source-help">
        <a href="https://lattice-automation.github.io/repp/" target="_blank" rel="noopener noreferrer">
          Part databases from REPP
        </a>
      </p>
    </div>
  );
};

const AssemblyDesigner: FC = () => {
  // Input state
  const [targetSeq, setTargetSeq] = useState<string>('');
  const [targetName, setTargetName] = useState<string>('');
  const [database, setDatabase] = useState<Fragment[]>([]);
  const [userDatabase, setUserDatabase] = useState<Fragment[]>([]); // User-uploaded fragments
  const [backboneName, setBackboneName] = useState<string>('');
  const [selectedEnzymes, setSelectedEnzymes] = useState<string[]>([]);

  // Part sources state
  const [enabledSources, setEnabledSources] = useState<string[]>([]);
  const [partDatabases, setPartDatabases] = useState<Record<string, Fragment[]>>({}); // { addgene: [...], igem: [...] }
  const [sourceLoadingStates, setSourceLoadingStates] = useState<Record<string, boolean>>({});
  const [sourceErrors, setSourceErrors] = useState<Record<string, string | null>>({});

  // Config state
  const [showConfig, setShowConfig] = useState<boolean>(false);
  const [config, setConfig] = useState<Config>({ ...DEFAULT_CONFIG } as Config);

  // Results state
  const [results, setResults] = useState<AssemblyResult | null>(null);
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);
  const [copyMessage, setCopyMessage] = useState<string>('');

  // Available part sources
  const availableSources = getAvailablePartSources() as PartSource[];

  // Merge user database with enabled part sources
  useEffect(() => {
    const enabledDatabases: Record<string, Fragment[]> = {};
    enabledSources.forEach(sourceId => {
      if (partDatabases[sourceId]) {
        enabledDatabases[sourceId] = partDatabases[sourceId];
      }
    });

    // Merge all enabled sources with user database
    const merged: Fragment[] = [
      ...userDatabase,
      ...(mergeDatabases(enabledDatabases) as Fragment[]),
    ];
    setDatabase(merged);
  }, [userDatabase, enabledSources, partDatabases]);

  // Toggle part source
  const togglePartSource = useCallback((sourceId: string): void => {
    setEnabledSources(prev =>
      prev.includes(sourceId)
        ? prev.filter(s => s !== sourceId)
        : [...prev, sourceId]
    );
  }, []);

  // Load part database from source
  const loadPartDatabase = useCallback(async (sourceId: string): Promise<void> => {
    // Check local cache first
    const cached = loadCachedDatabase(sourceId) as Fragment[] | null;
    if (cached && cached.length > 0) {
      setPartDatabases(prev => ({ ...prev, [sourceId]: cached }));
      return;
    }

    setSourceLoadingStates(prev => ({ ...prev, [sourceId]: true }));
    setSourceErrors(prev => ({ ...prev, [sourceId]: null }));

    try {
      const fragments = await fetchPartDatabase(sourceId) as Fragment[];
      setPartDatabases(prev => ({ ...prev, [sourceId]: fragments }));
      saveCachedDatabase(sourceId, fragments);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Unknown error';
      setSourceErrors(prev => ({ ...prev, [sourceId]: errorMessage }));
      console.error(`Failed to load ${sourceId} database:`, err);
    } finally {
      setSourceLoadingStates(prev => ({ ...prev, [sourceId]: false }));
    }
  }, []);

  // Auto-load databases when enabled
  useEffect(() => {
    for (const sourceId of enabledSources) {
      // If source is enabled but not loaded and not currently loading, load it
      if (!partDatabases[sourceId] && !sourceLoadingStates[sourceId]) {
        loadPartDatabase(sourceId);
      }
    }
  }, [enabledSources, partDatabases, sourceLoadingStates, loadPartDatabase]);

  // Get loaded counts for each source
  const loadedCounts: Record<string, number> = Object.fromEntries(
    Object.entries(partDatabases).map(([id, frags]) => [id, frags?.length || 0])
  );

  // Handle target file upload
  const handleTargetUpload = useCallback((e: ChangeEvent<HTMLInputElement>): void => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      try {
        const content = event.target?.result as string;
        const frags = parseFasta(content) as Fragment[];
        if (frags.length > 0) {
          setTargetSeq(frags[0].seq);
          setTargetName(frags[0].id);
          setError(null);
        }
      } catch (err) {
        const errorMessage = err instanceof Error ? err.message : 'Unknown error';
        setError('Failed to parse target file: ' + errorMessage);
      }
    };
    reader.readAsText(file);
  }, []);

  // Handle database file upload (user's own fragments)
  const handleDatabaseUpload = useCallback((e: ChangeEvent<HTMLInputElement>): void => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      try {
        const content = event.target?.result as string;
        const frags = parseFasta(content) as Fragment[];
        const dbFrags: Fragment[] = frags.map((f) => ({
          ...f,
          cost: 65,
          url: '',
          source: 'custom',
        }));
        setUserDatabase(dbFrags);
        setError(null);
      } catch (err) {
        const errorMessage = err instanceof Error ? err.message : 'Unknown error';
        setError('Failed to parse database file: ' + errorMessage);
      }
    };
    reader.readAsText(file);
  }, []);

  // Handle enzyme selection
  const toggleEnzyme = useCallback((enzyme: string): void => {
    setSelectedEnzymes(prev =>
      prev.includes(enzyme)
        ? prev.filter(e => e !== enzyme)
        : [...prev, enzyme]
    );
  }, []);

  // Run assembly planning
  const runAssembly = useCallback((): void => {
    if (!targetSeq) {
      setError('Please provide a target sequence');
      return;
    }

    setLoading(true);
    setError(null);
    setResults(null);

    setTimeout(() => {
      try {
        let finalTarget = targetSeq;

        if (backboneName && selectedEnzymes.length > 0) {
          const backbone = database.find(f => f.id === backboneName);
          if (backbone) {
            const digested = digest(backbone.seq, selectedEnzymes) as { linearized: string } | null;
            if (digested) {
              finalTarget = targetSeq + digested.linearized;
            } else {
              throw new Error(`No cut sites found for ${selectedEnzymes.join(', ')} in ${backboneName}`);
            }
          }
        }

        const result = planAssembly(finalTarget, database, {
          targetName: targetName || 'target',
          config,
        }) as AssemblyResult;

        setResults(result);
      } catch (err) {
        const errorMessage = err instanceof Error ? err.message : 'Unknown error';
        setError('Assembly planning failed: ' + errorMessage);
      } finally {
        setLoading(false);
      }
    }, 50);
  }, [targetSeq, targetName, database, backboneName, selectedEnzymes, config]);

  // Update config
  const updateConfig = useCallback((key: string, value: string): void => {
    const numValue = parseFloat(value);
    setConfig(prev => ({
      ...prev,
      [key]: isNaN(numValue) ? 0 : numValue,
    }));
  }, []);

  // Handle copy feedback
  const handleCopy = useCallback((message: string): void => {
    setCopyMessage(message);
    setTimeout(() => setCopyMessage(''), 2000);
  }, []);

  // Export all primers
  const exportAllPrimers = useCallback((solution: Solution): void => {
    const lines: string[] = [];
    solution.fragments.forEach((frag) => {
      if (frag.primers) {
        lines.push(`>${frag.id}_FWD`);
        lines.push(frag.primers[0].seq);
        lines.push(`>${frag.id}_REV`);
        lines.push(frag.primers[1].seq);
      }
    });
    const text = lines.join('\n');
    navigator.clipboard.writeText(text).then(() => {
      handleCopy('Copied all primers to clipboard');
    });
  }, [handleCopy]);

  return (
    <div className="assembly-designer">
      <div className="assembly-panels">
        {/* Input Panel */}
        <div className="assembly-input-panel">
          <div className="control-section">
            <h3>Target Sequence</h3>
            <textarea
              placeholder="Paste target plasmid sequence (ATGC) or upload FASTA..."
              value={targetSeq}
              onChange={(e) => setTargetSeq(e.target.value.toUpperCase().replace(/[^ATGC]/g, ''))}
              rows={4}
            />
            <div className="file-import">
              <input
                type="file"
                accept=".fa,.fasta,.txt"
                onChange={handleTargetUpload}
              />
            </div>
            {targetName && <small>Target: {targetName} ({targetSeq.length} bp)</small>}
          </div>

          <div className="control-section">
            <h3>Part Sources</h3>
            <p className="help-text">Enable part repositories to search for assembly fragments</p>
            <PartSourceSelector
              sources={availableSources}
              enabledSources={enabledSources}
              onToggle={togglePartSource}
              loadingStates={sourceLoadingStates}
              loadedCounts={loadedCounts}
              onLoad={loadPartDatabase}
            />
            {Object.entries(sourceErrors).map(([sourceId, err]) =>
              err ? (
                <div key={sourceId} className="source-error">
                  {(PART_SOURCES as Record<string, PartSource>)[sourceId]?.name || sourceId}: {err}
                </div>
              ) : null
            )}
          </div>

          <div className="control-section">
            <h3>Custom Fragments (Optional)</h3>
            <p className="help-text">Upload a FASTA file with your own DNA fragments</p>
            <div className="file-import">
              <input
                type="file"
                accept=".fa,.fasta,.txt"
                onChange={handleDatabaseUpload}
              />
            </div>
            {userDatabase.length > 0 && (
              <div className="database-summary">
                <small>{userDatabase.length} custom fragments loaded</small>
                <div className="fragment-list">
                  {userDatabase.slice(0, 5).map(f => (
                    <div key={f.id} className="fragment-item">
                      <span>{f.id}</span>
                      <span>{f.seq.length} bp</span>
                    </div>
                  ))}
                  {userDatabase.length > 5 && <small>...and {userDatabase.length - 5} more</small>}
                </div>
              </div>
            )}
            {database.length > 0 && (
              <div className="database-total">
                <strong>Total: {database.length.toLocaleString()} fragments available</strong>
              </div>
            )}
          </div>

          <div className="control-section">
            <h3>Backbone (Optional)</h3>
            <select
              value={backboneName}
              onChange={(e) => setBackboneName(e.target.value)}
            >
              <option value="">No backbone</option>
              {database.map(f => (
                <option key={f.id} value={f.id}>{f.id}</option>
              ))}
            </select>

            {backboneName && (
              <>
                <h4>Linearization Enzymes</h4>
                <div className="enzyme-select">
                  {COMMON_ENZYMES.map(enzyme => (
                    <label key={enzyme} className="enzyme-option">
                      <input
                        type="checkbox"
                        checked={selectedEnzymes.includes(enzyme)}
                        onChange={() => toggleEnzyme(enzyme)}
                      />
                      {enzyme}
                    </label>
                  ))}
                </div>
              </>
            )}
          </div>

          <div className="control-section">
            <button
              type="button"
              className="toggle-config"
              onClick={() => setShowConfig(!showConfig)}
            >
              {showConfig ? 'Hide' : 'Show'} Advanced Settings
            </button>

            {showConfig && (
              <div className="config-options">
                <div className="config-row">
                  <label>Max Fragments:</label>
                  <input
                    type="number"
                    value={config.fragmentsMaxCount}
                    onChange={(e) => updateConfig('fragmentsMaxCount', e.target.value)}
                    min={2}
                    max={10}
                  />
                </div>
                <div className="config-row">
                  <label>Min Junction (bp):</label>
                  <input
                    type="number"
                    value={config.fragmentsMinHomology}
                    onChange={(e) => updateConfig('fragmentsMinHomology', e.target.value)}
                    min={10}
                    max={50}
                  />
                </div>
                <div className="config-row">
                  <label>Max Junction (bp):</label>
                  <input
                    type="number"
                    value={config.fragmentsMaxHomology}
                    onChange={(e) => updateConfig('fragmentsMaxHomology', e.target.value)}
                    min={50}
                    max={200}
                  />
                </div>
                <div className="config-row">
                  <label>Min PCR Length (bp):</label>
                  <input
                    type="number"
                    value={config.pcrMinLength}
                    onChange={(e) => updateConfig('pcrMinLength', e.target.value)}
                    min={30}
                    max={200}
                  />
                </div>
                <div className="config-row">
                  <label>Primer Cost ($/bp):</label>
                  <input
                    type="number"
                    value={config.pcrBpCost}
                    onChange={(e) => updateConfig('pcrBpCost', e.target.value)}
                    step={0.1}
                    min={0}
                  />
                </div>
              </div>
            )}
          </div>

          <button
            type="button"
            className="run-btn"
            onClick={runAssembly}
            disabled={loading || !targetSeq}
          >
            {loading ? 'Planning...' : 'Plan Assembly'}
          </button>

          {error && <div className="error-message">{error}</div>}
        </div>

        {/* Results Panel */}
        <div className="assembly-results-panel">
          {copyMessage && <div className="copy-toast">{copyMessage}</div>}

          {loading && (
            <div className="loading">
              <div className="spinner"></div>
              <p>Finding optimal assembly...</p>
            </div>
          )}

          {!loading && !results && (
            <div className="empty-results">
              <p>Upload a target sequence and fragment database to plan your assembly.</p>
              <small>Results will appear here</small>
            </div>
          )}

          {results && results.solutions && (
            <div className="solutions">
              <div className="solutions-header">
                <h3>Assembly Solutions</h3>
                <p className="target-info">
                  Target: {results.target} ({results.seq.length} bp)
                </p>
              </div>

              {results.solutions.length === 0 ? (
                <div className="no-solutions">
                  <p>No assembly solutions found.</p>
                  <small>Try uploading more fragments or adjusting settings.</small>
                </div>
              ) : (
                results.solutions.map((solution, idx) => (
                  <div key={idx} className="solution-card">
                    <div className="solution-header">
                      <h4>Solution {idx + 1}</h4>
                      <div className="solution-stats">
                        <span className="stat-badge">{solution.count} fragments</span>
                        <span className="stat-badge cost">${solution.cost.toFixed(2)}</span>
                      </div>
                    </div>

                    {/* Assembly Visualization */}
                    <AssemblyVisualization
                      fragments={solution.fragments}
                      targetLength={results.seq.length}
                    />

                    {/* Cost Breakdown */}
                    <CostBreakdown solution={solution} config={config} />

                    {/* Fragments List */}
                    <div className="fragments-list">
                      <div className="fragments-list-header">
                        <h5>Fragments</h5>
                        {solution.fragments.some(f => f.primers) && (
                          <button
                            className="export-primers-btn"
                            onClick={() => exportAllPrimers(solution)}
                          >
                            Export All Primers
                          </button>
                        )}
                      </div>
                      {solution.fragments.map((frag, i) => (
                        <FragmentDetail
                          key={i}
                          frag={frag}
                          index={i}
                          onCopy={handleCopy}
                        />
                      ))}
                    </div>
                  </div>
                ))
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default AssemblyDesigner;
