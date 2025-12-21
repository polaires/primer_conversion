/**
 * Linear Assembly View
 * Linear visualization of fragment assembly with junctions
 */

import React from 'react';
import { Fragment, FRAGMENT_TYPES } from '../../../types/fragmentPlanner';

interface LinearAssemblyViewProps {
  fragments: Fragment[];
  overhangs: string[];
}

export function LinearAssemblyView({
  fragments,
  overhangs,
}: LinearAssemblyViewProps) {
  if (fragments.length === 0) {
    return (
      <div className="h-full flex items-center justify-center">
        <p className="text-slate-400 dark:text-slate-500">No fragments to display</p>
      </div>
    );
  }

  const totalLength = fragments.reduce((sum, f) => sum + f.length, 0);

  return (
    <div className="h-full flex flex-col">
      {/* Title */}
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-sm font-medium text-slate-700 dark:text-slate-300">
          Linear Assembly View
        </h3>
        <span className="text-xs text-slate-500 dark:text-slate-400">
          {totalLength.toLocaleString()} bp total
        </span>
      </div>

      {/* Assembly Track */}
      <div className="flex-1 flex flex-col justify-center">
        {/* Main assembly bar */}
        <div className="relative">
          {/* Connection line */}
          <div className="absolute left-0 right-0 top-1/2 h-0.5 bg-slate-300 dark:bg-slate-600 -translate-y-1/2" />

          {/* Circular indicator at start */}
          <div className="absolute -left-2 top-1/2 -translate-y-1/2 w-4 h-4 rounded-full border-2 border-slate-400 dark:border-slate-500 bg-white dark:bg-slate-800 flex items-center justify-center">
            <div className="w-1.5 h-1.5 rounded-full bg-slate-400 dark:bg-slate-500" />
          </div>

          {/* Fragment segments */}
          <div className="flex gap-0 relative z-10 ml-4">
            {fragments.map((fragment, idx) => {
              const typeInfo = FRAGMENT_TYPES[fragment.type];
              const widthPercent = (fragment.length / totalLength) * 100;
              const minWidth = 80;

              return (
                <div key={fragment.id} className="flex items-center">
                  {/* Junction (overhang) before fragment */}
                  {overhangs[idx] && (
                    <div className="flex flex-col items-center mx-1 z-20">
                      <div className="w-px h-3 bg-slate-400 dark:bg-slate-500" />
                      <div className="px-1.5 py-0.5 bg-white dark:bg-slate-800 border border-slate-300 dark:border-slate-600 rounded text-[10px] font-mono font-medium text-slate-600 dark:text-slate-400">
                        {overhangs[idx]}
                      </div>
                      <div className="w-px h-3 bg-slate-400 dark:bg-slate-500" />
                    </div>
                  )}

                  {/* Fragment block */}
                  <div
                    className="group relative flex-shrink-0 cursor-pointer"
                    style={{ minWidth: `${minWidth}px` }}
                  >
                    {/* Main block */}
                    <div
                      className="h-16 rounded-lg border-2 border-white dark:border-slate-700 shadow-sm flex flex-col items-center justify-center px-2 transition-transform hover:scale-105"
                      style={{ backgroundColor: typeInfo.color }}
                    >
                      {/* Type badge */}
                      <div className="absolute -top-2.5 left-1/2 -translate-x-1/2 w-5 h-5 rounded-full bg-white dark:bg-slate-800 border-2 flex items-center justify-center text-[10px] font-bold"
                        style={{ borderColor: typeInfo.color, color: typeInfo.color }}
                      >
                        {typeInfo.icon}
                      </div>

                      {/* Name */}
                      <span className="text-xs font-medium text-white truncate max-w-full">
                        {fragment.name}
                      </span>

                      {/* Length */}
                      <span className="text-[10px] text-white/80">
                        {fragment.length.toLocaleString()} bp
                      </span>
                    </div>

                    {/* Hover tooltip */}
                    <div className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 px-2 py-1 bg-slate-900 dark:bg-slate-700 text-white text-xs rounded opacity-0 group-hover:opacity-100 transition-opacity whitespace-nowrap pointer-events-none z-30">
                      <div className="font-medium">{fragment.name}</div>
                      <div className="text-slate-300">
                        {fragment.length.toLocaleString()} bp | GC: {fragment.gcContent.toFixed(1)}%
                      </div>
                    </div>
                  </div>
                </div>
              );
            })}

            {/* Final junction (closes the circle) */}
            {overhangs[fragments.length] && (
              <div className="flex flex-col items-center mx-1 z-20">
                <div className="w-px h-3 bg-slate-400 dark:bg-slate-500" />
                <div className="px-1.5 py-0.5 bg-white dark:bg-slate-800 border border-slate-300 dark:border-slate-600 rounded text-[10px] font-mono font-medium text-slate-600 dark:text-slate-400">
                  {overhangs[fragments.length]}
                </div>
                <div className="w-px h-3 bg-slate-400 dark:bg-slate-500" />
              </div>
            )}
          </div>

          {/* Circular indicator at end */}
          <div className="absolute -right-2 top-1/2 -translate-y-1/2 w-4 h-4 rounded-full border-2 border-slate-400 dark:border-slate-500 bg-white dark:bg-slate-800 flex items-center justify-center">
            <div className="w-1.5 h-1.5 rounded-full bg-slate-400 dark:bg-slate-500" />
          </div>
        </div>

        {/* Scale bar */}
        <div className="mt-8 flex items-center justify-center">
          <div className="flex items-center gap-2 text-xs text-slate-500 dark:text-slate-400">
            <div className="w-20 h-0.5 bg-slate-300 dark:bg-slate-600" />
            <span>{Math.round(totalLength / fragments.length).toLocaleString()} bp avg</span>
          </div>
        </div>
      </div>

      {/* Fragment summary table */}
      <div className="mt-6 overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="text-left text-xs text-slate-500 dark:text-slate-400 border-b border-slate-200 dark:border-slate-700">
              <th className="pb-2 font-medium">#</th>
              <th className="pb-2 font-medium">Fragment</th>
              <th className="pb-2 font-medium">Type</th>
              <th className="pb-2 font-medium text-right">Length</th>
              <th className="pb-2 font-medium text-right">GC%</th>
              <th className="pb-2 font-medium text-center">Left OH</th>
              <th className="pb-2 font-medium text-center">Right OH</th>
            </tr>
          </thead>
          <tbody>
            {fragments.map((fragment, idx) => {
              const typeInfo = FRAGMENT_TYPES[fragment.type];
              return (
                <tr
                  key={fragment.id}
                  className="border-b border-slate-100 dark:border-slate-800 hover:bg-slate-50 dark:hover:bg-slate-800/50"
                >
                  <td className="py-2 text-slate-400 dark:text-slate-500">
                    {idx + 1}
                  </td>
                  <td className="py-2 font-medium text-slate-900 dark:text-white">
                    {fragment.name}
                  </td>
                  <td className="py-2">
                    <span
                      className="inline-flex items-center gap-1 px-1.5 py-0.5 rounded text-xs font-medium text-white"
                      style={{ backgroundColor: typeInfo.color }}
                    >
                      {typeInfo.icon} {typeInfo.name}
                    </span>
                  </td>
                  <td className="py-2 text-right font-mono text-slate-600 dark:text-slate-300">
                    {fragment.length.toLocaleString()}
                  </td>
                  <td className="py-2 text-right text-slate-600 dark:text-slate-300">
                    {fragment.gcContent.toFixed(1)}%
                  </td>
                  <td className="py-2 text-center">
                    {overhangs[idx] ? (
                      <span className="font-mono text-xs bg-blue-100 dark:bg-blue-900/30 text-blue-700 dark:text-blue-400 px-1.5 py-0.5 rounded">
                        {overhangs[idx]}
                      </span>
                    ) : (
                      <span className="text-slate-300 dark:text-slate-600">-</span>
                    )}
                  </td>
                  <td className="py-2 text-center">
                    {overhangs[idx + 1] ? (
                      <span className="font-mono text-xs bg-violet-100 dark:bg-violet-900/30 text-violet-700 dark:text-violet-400 px-1.5 py-0.5 rounded">
                        {overhangs[idx + 1]}
                      </span>
                    ) : (
                      <span className="text-slate-300 dark:text-slate-600">-</span>
                    )}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}

export default LinearAssemblyView;
