/**
 * Constraint Panel
 * Configure assembly constraints for the optimizer
 */

import { useState } from 'react';
import { useDesignStore } from '../../stores/designStore';
import { GoldenGateEnzyme, GOLDEN_GATE_ENZYMES } from '../../types/fragmentPlanner';

export function ConstraintPanel() {
  const { constraints, updateConstraints, fragments, optimizeOverhangs, isOptimizing } =
    useDesignStore();
  const [isExpanded, setIsExpanded] = useState(true);

  const handleEnzymeChange = (enzyme: GoldenGateEnzyme) => {
    updateConstraints({ enzyme });
  };

  const handleFidelityChange = (value: number) => {
    updateConstraints({ minFidelity: value / 100 });
  };

  const handleMaxFragmentsChange = (value: number) => {
    updateConstraints({ maxFragments: value });
  };

  const enzymes = Object.entries(GOLDEN_GATE_ENZYMES) as [GoldenGateEnzyme, typeof GOLDEN_GATE_ENZYMES['BsaI']][];

  return (
    <div className="border-t border-slate-200 dark:border-slate-700">
      {/* Header */}
      <button
        onClick={() => setIsExpanded(!isExpanded)}
        className="w-full flex items-center justify-between p-4 text-left hover:bg-slate-50 dark:hover:bg-slate-700/50 transition-colors"
      >
        <div className="flex items-center gap-2">
          <svg
            className="w-5 h-5 text-slate-500 dark:text-slate-400"
            fill="none"
            viewBox="0 0 24 24"
            stroke="currentColor"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M12 6V4m0 2a2 2 0 100 4m0-4a2 2 0 110 4m-6 8a2 2 0 100-4m0 4a2 2 0 110-4m0 4v2m0-6V4m6 6v10m6-2a2 2 0 100-4m0 4a2 2 0 110-4m0 4v2m0-6V4"
            />
          </svg>
          <span className="text-sm font-semibold text-slate-900 dark:text-white uppercase tracking-wide">
            Constraints
          </span>
        </div>
        <svg
          className={`w-4 h-4 text-slate-400 transition-transform ${isExpanded ? 'rotate-180' : ''}`}
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Content */}
      {isExpanded && (
        <div className="px-4 pb-4 space-y-4">
          {/* Enzyme Selection */}
          <div>
            <label className="block text-xs font-medium text-slate-600 dark:text-slate-400 mb-2">
              Enzyme
            </label>
            <div className="grid grid-cols-2 gap-2">
              {enzymes.map(([enzyme, info]) => (
                <button
                  key={enzyme}
                  onClick={() => handleEnzymeChange(enzyme)}
                  className={`flex flex-col items-start p-2 rounded-lg border transition-all ${
                    constraints.enzyme === enzyme
                      ? 'border-amber-500 bg-amber-50 dark:bg-amber-900/20 dark:border-amber-400'
                      : 'border-slate-200 dark:border-slate-600 hover:border-slate-300 dark:hover:border-slate-500'
                  }`}
                >
                  <span
                    className={`text-sm font-medium ${
                      constraints.enzyme === enzyme
                        ? 'text-amber-700 dark:text-amber-400'
                        : 'text-slate-700 dark:text-slate-300'
                    }`}
                  >
                    {enzyme}
                  </span>
                  <span className="text-[10px] font-mono text-slate-500 dark:text-slate-400">
                    {info.recognition}
                  </span>
                </button>
              ))}
            </div>
          </div>

          {/* Minimum Fidelity */}
          <div>
            <div className="flex items-center justify-between mb-2">
              <label className="text-xs font-medium text-slate-600 dark:text-slate-400">
                Minimum Fidelity
              </label>
              <span className="text-xs font-semibold text-slate-900 dark:text-white">
                {(constraints.minFidelity * 100).toFixed(0)}%
              </span>
            </div>
            <input
              type="range"
              min={50}
              max={99}
              value={constraints.minFidelity * 100}
              onChange={(e) => handleFidelityChange(parseInt(e.target.value))}
              className="w-full h-2 bg-slate-200 dark:bg-slate-600 rounded-lg appearance-none cursor-pointer accent-amber-500"
            />
            <div className="flex justify-between text-[10px] text-slate-400 mt-1">
              <span>50%</span>
              <span>75%</span>
              <span>99%</span>
            </div>
          </div>

          {/* Max Fragments */}
          <div>
            <div className="flex items-center justify-between mb-2">
              <label className="text-xs font-medium text-slate-600 dark:text-slate-400">
                Max Fragments
              </label>
              <span className="text-xs font-semibold text-slate-900 dark:text-white">
                {constraints.maxFragments}
              </span>
            </div>
            <input
              type="range"
              min={2}
              max={8}
              value={constraints.maxFragments}
              onChange={(e) => handleMaxFragmentsChange(parseInt(e.target.value))}
              className="w-full h-2 bg-slate-200 dark:bg-slate-600 rounded-lg appearance-none cursor-pointer accent-amber-500"
            />
            <div className="flex justify-between text-[10px] text-slate-400 mt-1">
              <span>2</span>
              <span>5</span>
              <span>8</span>
            </div>
          </div>

          {/* Target Organism */}
          <div>
            <label className="block text-xs font-medium text-slate-600 dark:text-slate-400 mb-2">
              Target Organism
            </label>
            <select
              value={constraints.targetOrganism}
              onChange={(e) =>
                updateConstraints({
                  targetOrganism: e.target.value as 'ecoli' | 'yeast' | 'mammalian',
                })
              }
              className="w-full px-3 py-2 text-sm border border-slate-200 dark:border-slate-600 rounded-lg bg-white dark:bg-slate-700 text-slate-900 dark:text-white focus:ring-2 focus:ring-amber-500 focus:border-transparent"
            >
              <option value="ecoli">E. coli</option>
              <option value="yeast">Yeast</option>
              <option value="mammalian">Mammalian</option>
            </select>
          </div>

          {/* Advanced Options */}
          <div className="space-y-2">
            <label className="flex items-center gap-2 cursor-pointer">
              <input
                type="checkbox"
                checked={constraints.preferHighFidelityOverhangs}
                onChange={(e) =>
                  updateConstraints({ preferHighFidelityOverhangs: e.target.checked })
                }
                className="w-4 h-4 rounded border-slate-300 dark:border-slate-600 text-amber-500 focus:ring-amber-500"
              />
              <span className="text-sm text-slate-700 dark:text-slate-300">
                Prefer NEB-validated overhangs
              </span>
            </label>

            <label className="flex items-center gap-2 cursor-pointer">
              <input
                type="checkbox"
                checked={constraints.allowDomestication}
                onChange={(e) =>
                  updateConstraints({ allowDomestication: e.target.checked })
                }
                className="w-4 h-4 rounded border-slate-300 dark:border-slate-600 text-amber-500 focus:ring-amber-500"
              />
              <span className="text-sm text-slate-700 dark:text-slate-300">
                Allow sequence domestication
              </span>
            </label>
          </div>

          {/* Optimize Button */}
          <button
            onClick={() => optimizeOverhangs()}
            disabled={isOptimizing || fragments.length < 2}
            className="w-full py-2.5 bg-gradient-to-r from-amber-500 to-orange-500 hover:from-amber-600 hover:to-orange-600 text-white font-medium rounded-lg shadow-sm disabled:opacity-50 disabled:cursor-not-allowed transition-all flex items-center justify-center gap-2"
          >
            {isOptimizing ? (
              <>
                <svg className="w-4 h-4 animate-spin" fill="none" viewBox="0 0 24 24">
                  <circle
                    className="opacity-25"
                    cx="12"
                    cy="12"
                    r="10"
                    stroke="currentColor"
                    strokeWidth="4"
                  />
                  <path
                    className="opacity-75"
                    fill="currentColor"
                    d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                  />
                </svg>
                Optimizing...
              </>
            ) : (
              <>
                <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path
                    strokeLinecap="round"
                    strokeLinejoin="round"
                    strokeWidth={2}
                    d="M13 10V3L4 14h7v7l9-11h-7z"
                  />
                </svg>
                Optimize Overhangs
              </>
            )}
          </button>

          {fragments.length < 2 && (
            <p className="text-xs text-center text-slate-500 dark:text-slate-400">
              Add at least 2 fragments to optimize
            </p>
          )}
        </div>
      )}
    </div>
  );
}

export default ConstraintPanel;
