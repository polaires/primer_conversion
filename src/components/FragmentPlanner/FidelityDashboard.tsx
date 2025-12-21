/**
 * Fidelity Dashboard
 * Live metrics display for assembly fidelity and overhang status
 */

import { getFidelityLevel } from '../../types/fragmentPlanner';

interface FidelityDashboardProps {
  fidelity: number;
  overhangs: string[];
  isOptimizing: boolean;
  warnings: string[];
  errors: string[];
  onOptimize: () => void;
}

export function FidelityDashboard({
  fidelity,
  overhangs,
  isOptimizing,
  warnings,
  errors,
  onOptimize,
}: FidelityDashboardProps) {
  const fidelityLevel = getFidelityLevel(fidelity);
  const fidelityPercent = (fidelity * 100).toFixed(1);

  return (
    <div className="bg-white dark:bg-slate-800 rounded-xl border border-slate-200 dark:border-slate-700 p-4">
      <div className="flex items-start gap-6">
        {/* Fidelity Gauge */}
        <div className="flex-shrink-0">
          <div className="relative w-24 h-24">
            {/* Background circle */}
            <svg className="w-24 h-24 transform -rotate-90">
              <circle
                cx="48"
                cy="48"
                r="40"
                stroke="currentColor"
                strokeWidth="8"
                fill="none"
                className="text-slate-200 dark:text-slate-700"
              />
              {/* Fidelity arc */}
              <circle
                cx="48"
                cy="48"
                r="40"
                stroke={fidelityLevel.color}
                strokeWidth="8"
                fill="none"
                strokeLinecap="round"
                strokeDasharray={`${fidelity * 251.2} 251.2`}
                className="transition-all duration-500"
              />
            </svg>
            {/* Center text */}
            <div className="absolute inset-0 flex flex-col items-center justify-center">
              <span className="text-xl font-bold text-slate-900 dark:text-white">
                {fidelity > 0 ? fidelityPercent : '--'}
              </span>
              <span className="text-xs text-slate-500 dark:text-slate-400">%</span>
            </div>
          </div>
          <div className="text-center mt-2">
            <span
              className="text-xs font-medium px-2 py-0.5 rounded-full"
              style={{
                backgroundColor: `${fidelityLevel.color}20`,
                color: fidelityLevel.color,
              }}
            >
              {fidelityLevel.label}
            </span>
          </div>
        </div>

        {/* Overhangs Display */}
        <div className="flex-1 min-w-0">
          <h3 className="text-sm font-semibold text-slate-900 dark:text-white mb-2">
            Assigned Overhangs
          </h3>
          {overhangs.length > 0 ? (
            <div className="flex flex-wrap gap-2">
              {overhangs.map((oh, idx) => (
                <div key={idx} className="flex items-center gap-1">
                  <span className="px-2 py-1 text-xs font-mono font-medium bg-slate-100 dark:bg-slate-700 text-slate-700 dark:text-slate-300 rounded">
                    {oh}
                  </span>
                  {idx < overhangs.length - 1 && (
                    <svg
                      className="w-3 h-3 text-slate-400"
                      fill="none"
                      viewBox="0 0 24 24"
                      stroke="currentColor"
                    >
                      <path
                        strokeLinecap="round"
                        strokeLinejoin="round"
                        strokeWidth={2}
                        d="M9 5l7 7-7 7"
                      />
                    </svg>
                  )}
                </div>
              ))}
            </div>
          ) : (
            <div className="flex items-center gap-2 text-slate-500 dark:text-slate-400">
              <svg
                className="w-4 h-4"
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"
                />
              </svg>
              <span className="text-sm">
                Click "Optimize Overhangs" to assign junction sequences
              </span>
            </div>
          )}

          {/* Cross-ligation info */}
          {overhangs.length > 0 && (
            <div className="mt-3 flex items-center gap-4 text-xs text-slate-500 dark:text-slate-400">
              <span className="flex items-center gap-1">
                <span className="w-2 h-2 rounded-full bg-green-500" />
                No self-complementary
              </span>
              <span className="flex items-center gap-1">
                <span className="w-2 h-2 rounded-full bg-green-500" />
                No palindromes
              </span>
            </div>
          )}
        </div>

        {/* Status/Actions */}
        <div className="flex-shrink-0">
          {isOptimizing ? (
            <div className="flex items-center gap-2 text-amber-600 dark:text-amber-400">
              <svg className="w-5 h-5 animate-spin" fill="none" viewBox="0 0 24 24">
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
              <span className="text-sm font-medium">Optimizing...</span>
            </div>
          ) : overhangs.length > 0 ? (
            <button
              onClick={onOptimize}
              className="flex items-center gap-2 px-3 py-1.5 text-sm text-slate-600 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 rounded-lg transition-colors"
            >
              <svg
                className="w-4 h-4"
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15"
                />
              </svg>
              Re-optimize
            </button>
          ) : null}
        </div>
      </div>

      {/* Warnings & Errors */}
      {(warnings.length > 0 || errors.length > 0) && (
        <div className="mt-4 space-y-2">
          {errors.map((error, idx) => (
            <div
              key={`error-${idx}`}
              className="flex items-start gap-2 p-2 bg-red-50 dark:bg-red-900/20 rounded-lg"
            >
              <svg
                className="w-4 h-4 text-red-500 flex-shrink-0 mt-0.5"
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M12 8v4m0 4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"
                />
              </svg>
              <span className="text-sm text-red-700 dark:text-red-400">{error}</span>
            </div>
          ))}
          {warnings.map((warning, idx) => (
            <div
              key={`warning-${idx}`}
              className="flex items-start gap-2 p-2 bg-amber-50 dark:bg-amber-900/20 rounded-lg"
            >
              <svg
                className="w-4 h-4 text-amber-500 flex-shrink-0 mt-0.5"
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z"
                />
              </svg>
              <span className="text-sm text-amber-700 dark:text-amber-400">{warning}</span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

export default FidelityDashboard;
