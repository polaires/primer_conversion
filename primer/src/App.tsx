import { useState } from 'react';
import PrimerForm from './components/PrimerForm';
import PrimerResults from './components/PrimerResults';
import StandaloneViewer from './components/StandaloneViewer';
import AssemblyDesigner from './components/AssemblyDesigner';
import GoldenGateDesigner from './components/GoldenGateDesigner';
import TmCalculator from './components/TmCalculator';
import SequencingDesigner from './components/SequencingDesigner';
import UnifiedPrimerDesigner from './components/UnifiedPrimerDesigner';
import EnhancedScorer from './components/EnhancedScorer';
import ErrorBoundary from './components/ErrorBoundary';
import { primers, score } from './lib/index.js';
import type { Tool, ToolMode } from './types';

// Library returns generic objects, we'll type them loosely during migration
// eslint-disable-next-line @typescript-eslint/no-explicit-any
type PrimerResult = any;
// eslint-disable-next-line @typescript-eslint/no-explicit-any
type LibOptions = any;

// Tool definitions with metadata
const TOOLS: Tool[] = [
  { id: 'primer-designer', label: 'Primer Designer', category: 'design', description: 'Unified PCR & mutagenesis primer design' },
  { id: 'sequencing', label: 'Sequencing', category: 'design', description: 'Sanger sequencing primer design' },
  { id: 'score', label: 'Score Primers', category: 'analysis', description: 'Evaluate existing primers' },
  { id: 'tm', label: 'Tm Calculator', category: 'analysis', description: 'NEB Q5 melting temperature' },
  { id: 'viewer', label: 'Sequence Viewer', category: 'analysis', description: 'Visualize DNA sequences' },
  { id: 'assembly-studio', label: 'Assembly Studio', category: 'assembly', description: 'Golden Gate, Gibson & NEBuilder HiFi assembly' },
  { id: 'fragment-planner', label: 'Fragment Planner', category: 'assembly', description: 'Cost-optimized assembly from fragment library' },
];

interface CreateSubmitData {
  sequence: string;
  options: Record<string, unknown>;
}

interface ScoreSubmitData {
  fwdPrimer: string;
  revPrimer: string;
  sequence: string;
  options: Record<string, unknown> & { offtargetCheck?: boolean };
}

export default function App() {
  const [mode, setMode] = useState<ToolMode>('primer-designer');
  const [results, setResults] = useState<[PrimerResult, PrimerResult] | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState<boolean>(false);
  const [currentSequence, setCurrentSequence] = useState<string>('');
  const [currentOptions, setCurrentOptions] = useState<Record<string, unknown>>({});

  const handleCreate = ({ sequence, options }: CreateSubmitData): void => {
    setLoading(true);
    setError(null);
    setResults(null);
    setCurrentSequence(sequence);
    setCurrentOptions(options);

    setTimeout(() => {
      try {
        const [fwd, rev] = primers(sequence, options as LibOptions) as [PrimerResult, PrimerResult];
        setResults([fwd, rev]);
      } catch (err) {
        setError(err instanceof Error ? err.message : 'An error occurred');
      } finally {
        setLoading(false);
      }
    }, 10);
  };

  const handleScore = ({ fwdPrimer, revPrimer, sequence, options }: ScoreSubmitData): void => {
    setLoading(true);
    setError(null);
    setResults(null);
    setCurrentSequence(sequence);
    setCurrentOptions(options);

    setTimeout(() => {
      try {
        const [fwd, rev] = score(fwdPrimer, revPrimer, sequence, options.offtargetCheck as LibOptions, options as LibOptions) as [PrimerResult, PrimerResult];
        setResults([fwd, rev]);
      } catch (err) {
        setError(err instanceof Error ? err.message : 'An error occurred');
      } finally {
        setLoading(false);
      }
    }, 10);
  };

  const switchMode = (newMode: ToolMode): void => {
    setMode(newMode);
    setResults(null);
    setError(null);
  };

  // Group tools by category
  const designTools = TOOLS.filter((t) => t.category === 'design');
  const analysisTools = TOOLS.filter((t) => t.category === 'analysis');
  const assemblyTools = TOOLS.filter((t) => t.category === 'assembly');

  return (
    <div className="app">
      <a href="#main-content" className="skip-link">
        Skip to main content
      </a>
      <header className="app-header">
        <div className="header-content">
          <h1>Primers</h1>
          <p className="tagline">Complete PCR Primer Design Suite</p>
        </div>
      </header>

      <nav className="app-nav" aria-label="Main navigation">
        <div className="nav-section" role="group" aria-labelledby="design-tools-label">
          <span id="design-tools-label" className="nav-label">Design</span>
          <div className="nav-buttons" role="tablist">
            {designTools.map((tool) => (
              <button
                key={tool.id}
                type="button"
                role="tab"
                aria-selected={mode === tool.id}
                aria-controls="main-content"
                className={mode === tool.id ? 'active' : ''}
                onClick={() => switchMode(tool.id)}
                title={tool.description}
              >
                {tool.label}
              </button>
            ))}
          </div>
        </div>
        <div className="nav-section" role="group" aria-labelledby="analysis-tools-label">
          <span id="analysis-tools-label" className="nav-label">Analysis</span>
          <div className="nav-buttons" role="tablist">
            {analysisTools.map((tool) => (
              <button
                key={tool.id}
                type="button"
                role="tab"
                aria-selected={mode === tool.id}
                aria-controls="main-content"
                className={mode === tool.id ? 'active' : ''}
                onClick={() => switchMode(tool.id)}
                title={tool.description}
              >
                {tool.label}
              </button>
            ))}
          </div>
        </div>
        <div className="nav-section" role="group" aria-labelledby="assembly-tools-label">
          <span id="assembly-tools-label" className="nav-label">Assembly</span>
          <div className="nav-buttons" role="tablist">
            {assemblyTools.map((tool) => (
              <button
                key={tool.id}
                type="button"
                role="tab"
                aria-selected={mode === tool.id}
                aria-controls="main-content"
                className={mode === tool.id ? 'active' : ''}
                onClick={() => switchMode(tool.id)}
                title={tool.description}
              >
                {tool.label}
              </button>
            ))}
          </div>
        </div>
      </nav>

      <main className="app-main" id="main-content" role="tabpanel" aria-label="Tool content">
        <ErrorBoundary>
          {mode === 'viewer' ? (
            <StandaloneViewer />
          ) : mode === 'assembly-studio' ? (
            <GoldenGateDesigner />
          ) : mode === 'fragment-planner' ? (
            <AssemblyDesigner />
          ) : mode === 'tm' ? (
            <TmCalculator />
          ) : mode === 'sequencing' ? (
            <SequencingDesigner />
          ) : mode === 'primer-designer' ? (
            <UnifiedPrimerDesigner />
          ) : mode === 'score' ? (
            <EnhancedScorer />
          ) : (
          <div className="content">
            <section className="form-section">
              <PrimerForm
                mode={mode}
                onSubmit={mode === 'create' ? handleCreate : handleScore}
              />
            </section>

            <section className="results-section">
              {loading && (
                <div className="loading">
                  <div className="spinner"></div>
                  <p>Designing optimal primers...</p>
                </div>
              )}
              <PrimerResults
                results={results}
                error={error}
                sequence={currentSequence}
                options={currentOptions}
              />
            </section>
          </div>
          )}
        </ErrorBoundary>
      </main>

      <footer className="app-footer">
        <p>
          Based on{' '}
          <a
            href="https://github.com/Lattice-Automation/primers"
            target="_blank"
            rel="noopener noreferrer"
          >
            primers
          </a>,{' '}
          <a
            href="https://github.com/Lattice-Automation/repp"
            target="_blank"
            rel="noopener noreferrer"
          >
            repp
          </a>, and{' '}
          <a
            href="https://github.com/Lattice-Automation/synbio"
            target="_blank"
            rel="noopener noreferrer"
          >
            synbio
          </a>{' '}
          by Lattice Automation
        </p>
        <p className="algorithm-info">
          Algorithms: NEB Q5 Tm | SantaLucia Nearest-Neighbor | Monte Carlo Optimization
        </p>
      </footer>
    </div>
  );
}
