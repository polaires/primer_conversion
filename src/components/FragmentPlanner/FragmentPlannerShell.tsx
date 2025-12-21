/**
 * Fragment Planner Shell
 * Main container component for the enhanced plasmid assembly tool
 */

import { useEffect, useCallback } from 'react';
import { useDesignStore, useHistoryKeyboardShortcuts } from '../../stores/designStore';
import { DesignTimeline } from './DesignTimeline';
import { FragmentList } from './FragmentList';
import { ConstraintPanel } from './ConstraintPanel';
import { FidelityDashboard } from './FidelityDashboard';
import { CircularPlasmidView } from './AssemblyViz/CircularPlasmidView';
import { LinearAssemblyView } from './AssemblyViz/LinearAssemblyView';
import { ActionBar } from './ActionBar';

export function FragmentPlannerShell() {
  // Initialize keyboard shortcuts
  useEffect(() => {
    const cleanup = useHistoryKeyboardShortcuts();
    return cleanup;
  }, []);

  const {
    fragments,
    overhangs,
    fidelity,
    viewMode,
    isOptimizing,
    warnings,
    errors,
    undo,
    redo,
    canUndo,
    canRedo,
    setViewMode,
    optimizeOverhangs,
  } = useDesignStore();

  const totalLength = fragments.reduce((sum, f) => sum + f.length, 0);

  const handleOptimize = useCallback(async () => {
    await optimizeOverhangs();
  }, [optimizeOverhangs]);

  return (
    <div className="min-h-screen bg-slate-50 dark:bg-slate-900">
      {/* Header */}
      <header className="bg-white dark:bg-slate-800 border-b border-slate-200 dark:border-slate-700 px-6 py-4">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4">
            <div className="flex items-center gap-2">
              <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-amber-400 to-orange-500 flex items-center justify-center">
                <svg
                  className="w-5 h-5 text-white"
                  fill="none"
                  viewBox="0 0 24 24"
                  stroke="currentColor"
                >
                  <path
                    strokeLinecap="round"
                    strokeLinejoin="round"
                    strokeWidth={2}
                    d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z"
                  />
                </svg>
              </div>
              <h1 className="text-xl font-semibold text-slate-900 dark:text-white">
                Fragment Planner
              </h1>
            </div>
            <span className="px-2 py-1 text-xs font-medium bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400 rounded-full">
              Golden Gate Assembly
            </span>
          </div>

          {/* Undo/Redo Controls */}
          <div className="flex items-center gap-2">
            <button
              onClick={undo}
              disabled={!canUndo()}
              className="p-2 rounded-lg text-slate-500 hover:text-slate-700 hover:bg-slate-100 dark:text-slate-400 dark:hover:text-slate-200 dark:hover:bg-slate-700 disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
              title="Undo (Ctrl+Z)"
            >
              <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M3 10h10a8 8 0 018 8v2M3 10l6 6m-6-6l6-6"
                />
              </svg>
            </button>
            <button
              onClick={redo}
              disabled={!canRedo()}
              className="p-2 rounded-lg text-slate-500 hover:text-slate-700 hover:bg-slate-100 dark:text-slate-400 dark:hover:text-slate-200 dark:hover:bg-slate-700 disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
              title="Redo (Ctrl+Shift+Z)"
            >
              <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M21 10h-10a8 8 0 00-8 8v2M21 10l-6 6m6-6l-6-6"
                />
              </svg>
            </button>

            <div className="w-px h-6 bg-slate-200 dark:bg-slate-700 mx-2" />

            {/* View Mode Toggle */}
            <div className="flex items-center bg-slate-100 dark:bg-slate-700 rounded-lg p-1">
              <button
                onClick={() => setViewMode('circular')}
                className={`px-3 py-1.5 text-sm font-medium rounded-md transition-colors ${
                  viewMode === 'circular'
                    ? 'bg-white dark:bg-slate-600 text-slate-900 dark:text-white shadow-sm'
                    : 'text-slate-600 dark:text-slate-300 hover:text-slate-900 dark:hover:text-white'
                }`}
              >
                Circular
              </button>
              <button
                onClick={() => setViewMode('linear')}
                className={`px-3 py-1.5 text-sm font-medium rounded-md transition-colors ${
                  viewMode === 'linear'
                    ? 'bg-white dark:bg-slate-600 text-slate-900 dark:text-white shadow-sm'
                    : 'text-slate-600 dark:text-slate-300 hover:text-slate-900 dark:hover:text-white'
                }`}
              >
                Linear
              </button>
            </div>
          </div>
        </div>
      </header>

      {/* Main Content */}
      <div className="flex h-[calc(100vh-73px)]">
        {/* Left Sidebar - Fragments & Constraints */}
        <aside className="w-80 bg-white dark:bg-slate-800 border-r border-slate-200 dark:border-slate-700 flex flex-col">
          {/* Fragment List */}
          <div className="flex-1 overflow-hidden flex flex-col">
            <div className="p-4 border-b border-slate-200 dark:border-slate-700">
              <h2 className="text-sm font-semibold text-slate-900 dark:text-white uppercase tracking-wide">
                Fragments
              </h2>
              <p className="text-xs text-slate-500 dark:text-slate-400 mt-1">
                {fragments.length} fragments | {totalLength.toLocaleString()} bp total
              </p>
            </div>
            <div className="flex-1 overflow-y-auto">
              <FragmentList />
            </div>
          </div>

          {/* Constraints Panel */}
          <div className="border-t border-slate-200 dark:border-slate-700">
            <ConstraintPanel />
          </div>
        </aside>

        {/* Main View */}
        <main className="flex-1 flex flex-col overflow-hidden">
          {/* Visualization Area */}
          <div className="flex-1 p-6 overflow-auto">
            {fragments.length === 0 ? (
              <EmptyState />
            ) : (
              <div className="h-full flex flex-col">
                {/* Assembly Visualization */}
                <div className="flex-1 bg-white dark:bg-slate-800 rounded-xl border border-slate-200 dark:border-slate-700 p-6 min-h-[400px]">
                  {viewMode === 'circular' ? (
                    <CircularPlasmidView
                      fragments={fragments}
                      overhangs={overhangs}
                      totalLength={totalLength}
                    />
                  ) : (
                    <LinearAssemblyView
                      fragments={fragments}
                      overhangs={overhangs}
                    />
                  )}
                </div>

                {/* Fidelity Dashboard */}
                <div className="mt-4">
                  <FidelityDashboard
                    fidelity={fidelity}
                    overhangs={overhangs}
                    isOptimizing={isOptimizing}
                    warnings={warnings}
                    errors={errors}
                    onOptimize={handleOptimize}
                  />
                </div>
              </div>
            )}
          </div>

          {/* Action Bar */}
          <ActionBar />
        </main>

        {/* Right Sidebar - History Timeline */}
        <aside className="w-72 bg-white dark:bg-slate-800 border-l border-slate-200 dark:border-slate-700">
          <DesignTimeline />
        </aside>
      </div>
    </div>
  );
}

function EmptyState() {
  const { addFragment } = useDesignStore();

  const handleAddExample = () => {
    // Add example fragments
    addFragment({
      name: 'pTrc Promoter',
      sequence: 'TTGACAATTAATCATCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCT',
      type: 'promoter',
    });
    addFragment({
      name: 'RBS B0034',
      sequence: 'AAAGAGGAGAAA',
      type: 'rbs',
    });
    addFragment({
      name: 'GFP',
      sequence: 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA',
      type: 'cds',
    });
    addFragment({
      name: 'T7 Terminator',
      sequence: 'CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG',
      type: 'terminator',
    });
  };

  return (
    <div className="h-full flex items-center justify-center">
      <div className="text-center max-w-md">
        <div className="w-16 h-16 rounded-2xl bg-slate-100 dark:bg-slate-700 flex items-center justify-center mx-auto mb-4">
          <svg
            className="w-8 h-8 text-slate-400"
            fill="none"
            viewBox="0 0 24 24"
            stroke="currentColor"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={1.5}
              d="M19 11H5m14 0a2 2 0 012 2v6a2 2 0 01-2 2H5a2 2 0 01-2-2v-6a2 2 0 012-2m14 0V9a2 2 0 00-2-2M5 11V9a2 2 0 012-2m0 0V5a2 2 0 012-2h6a2 2 0 012 2v2M7 7h10"
            />
          </svg>
        </div>
        <h3 className="text-lg font-semibold text-slate-900 dark:text-white mb-2">
          No fragments yet
        </h3>
        <p className="text-slate-500 dark:text-slate-400 mb-6">
          Add DNA fragments to start designing your Golden Gate assembly. You can paste
          sequences, import from files, or load from part registries.
        </p>
        <div className="flex flex-col gap-3">
          <button
            onClick={handleAddExample}
            className="px-4 py-2 bg-amber-500 hover:bg-amber-600 text-white font-medium rounded-lg transition-colors"
          >
            Load Example Design
          </button>
          <button className="px-4 py-2 border border-slate-300 dark:border-slate-600 text-slate-700 dark:text-slate-300 font-medium rounded-lg hover:bg-slate-50 dark:hover:bg-slate-700 transition-colors">
            Import from File
          </button>
        </div>
      </div>
    </div>
  );
}

export default FragmentPlannerShell;
