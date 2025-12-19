/**
 * SequenceConflictMap Component
 *
 * Visualizes WHY primer design failed by highlighting:
 * - High/Low GC content regions (Red = High GC, Blue = Low GC)
 * - Strong secondary structure regions (Yellow/Orange)
 * - Repeat sequences
 *
 * Unlike NEB BaseChanger which just says "No primers found",
 * this shows the user exactly what's wrong and where.
 */

import React, { useMemo, useState } from 'react';
import {
  analyzeGCHotspots,
  analyzeStructureHotspots,
  analyzeDesignConflicts,
  gcToColor,
  dgToColor,
} from '../lib/visualization.js';

const SEVERITY_COLORS = {
  critical: '#dc2626',  // red-600
  severe: '#f97316',    // orange-500
  warning: '#facc15',   // yellow-400
  ok: '#22c55e',        // green-500
};

const SEVERITY_BG = {
  critical: 'bg-red-100 border-red-300',
  severe: 'bg-orange-100 border-orange-300',
  warning: 'bg-yellow-100 border-yellow-300',
  ok: 'bg-green-100 border-green-300',
};

/**
 * Main conflict visualization component
 */
export function SequenceConflictMap({
  sequence,
  mutationPosition,
  error,
  showDetails = true,
  height = 120,
}) {
  const [activeTab, setActiveTab] = useState('gc');
  const [hoveredPosition, setHoveredPosition] = useState(null);

  // Analyze conflicts
  const analysis = useMemo(() => {
    if (!sequence || sequence.length < 20) return null;
    return analyzeDesignConflicts(sequence, mutationPosition);
  }, [sequence, mutationPosition]);

  // Separate GC and structure analyses for full sequence
  const gcAnalysis = useMemo(() => {
    if (!sequence || sequence.length < 20) return null;
    return analyzeGCHotspots(sequence, 15);
  }, [sequence]);

  const structureAnalysis = useMemo(() => {
    if (!sequence || sequence.length < 20) return null;
    return analyzeStructureHotspots(sequence, 25);
  }, [sequence]);

  if (!sequence || sequence.length < 20) {
    return (
      <div className="p-4 bg-gray-100 rounded text-gray-500 text-sm">
        Sequence too short for analysis (minimum 20bp)
      </div>
    );
  }

  const conflicts = analysis?.conflicts || [];
  const hasCriticalConflicts = conflicts.some(c => c.severity === 'critical' || c.severity === 'severe');

  return (
    <div className="space-y-4">
      {/* Conflict Summary */}
      {conflicts.length > 0 && (
        <ConflictSummary conflicts={conflicts} canDesign={analysis?.canDesign} />
      )}

      {/* Tab Navigation */}
      <div className="flex border-b border-gray-200">
        <TabButton
          active={activeTab === 'gc'}
          onClick={() => setActiveTab('gc')}
          label="GC Content"
          icon={
            <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
            </svg>
          }
        />
        <TabButton
          active={activeTab === 'structure'}
          onClick={() => setActiveTab('structure')}
          label="Secondary Structure"
          icon={
            <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 5a1 1 0 011-1h14a1 1 0 011 1v2a1 1 0 01-1 1H5a1 1 0 01-1-1V5zM4 13a1 1 0 011-1h6a1 1 0 011 1v6a1 1 0 01-1 1H5a1 1 0 01-1-1v-6zM16 13a1 1 0 011-1h2a1 1 0 011 1v6a1 1 0 01-1 1h-2a1 1 0 01-1-1v-6z" />
            </svg>
          }
        />
        <TabButton
          active={activeTab === 'combined'}
          onClick={() => setActiveTab('combined')}
          label="Combined View"
          icon={
            <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2V6zM14 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2V6zM4 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2v-2zM14 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2v-2z" />
            </svg>
          }
        />
      </div>

      {/* Visualization */}
      <div className="relative">
        {activeTab === 'gc' && (
          <GCVisualization
            sequence={sequence}
            analysis={gcAnalysis}
            mutationPosition={mutationPosition}
            height={height}
            onHover={setHoveredPosition}
          />
        )}
        {activeTab === 'structure' && (
          <StructureVisualization
            sequence={sequence}
            analysis={structureAnalysis}
            mutationPosition={mutationPosition}
            height={height}
            onHover={setHoveredPosition}
          />
        )}
        {activeTab === 'combined' && (
          <CombinedVisualization
            sequence={sequence}
            gcAnalysis={gcAnalysis}
            structureAnalysis={structureAnalysis}
            mutationPosition={mutationPosition}
            height={height * 2}
            onHover={setHoveredPosition}
          />
        )}

        {/* Position Tooltip */}
        {hoveredPosition !== null && (
          <PositionTooltip
            position={hoveredPosition}
            sequence={sequence}
            gcAnalysis={gcAnalysis}
            structureAnalysis={structureAnalysis}
          />
        )}
      </div>

      {/* Suggestions */}
      {analysis?.suggestions?.length > 0 && (
        <Suggestions suggestions={analysis.suggestions} />
      )}

      {/* Legend */}
      {showDetails && <Legend activeTab={activeTab} />}
    </div>
  );
}

/**
 * Tab button component
 */
function TabButton({ active, onClick, label, icon }) {
  return (
    <button
      onClick={onClick}
      className={`flex items-center gap-2 px-4 py-2 text-sm font-medium border-b-2 transition-colors ${
        active
          ? 'border-blue-500 text-blue-600'
          : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
      }`}
    >
      {icon}
      {label}
    </button>
  );
}

/**
 * Conflict summary cards
 */
function ConflictSummary({ conflicts, canDesign }) {
  const grouped = conflicts.reduce((acc, c) => {
    acc[c.severity] = acc[c.severity] || [];
    acc[c.severity].push(c);
    return acc;
  }, {});

  return (
    <div className={`p-4 rounded-lg border ${canDesign ? 'bg-yellow-50 border-yellow-200' : 'bg-red-50 border-red-200'}`}>
      <div className="flex items-center gap-2 mb-2">
        {!canDesign ? (
          <svg className="w-5 h-5 text-red-500" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
          </svg>
        ) : (
          <svg className="w-5 h-5 text-yellow-500" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
          </svg>
        )}
        <span className={`font-semibold ${canDesign ? 'text-yellow-800' : 'text-red-800'}`}>
          {canDesign ? 'Design Warnings' : 'Cannot Design Primers'}
        </span>
      </div>
      <ul className="space-y-1 text-sm">
        {conflicts.map((conflict, i) => (
          <li key={i} className="flex items-start gap-2">
            <span
              className="w-2 h-2 rounded-full mt-1.5 flex-shrink-0"
              style={{ backgroundColor: SEVERITY_COLORS[conflict.severity] }}
            />
            <span className={conflict.severity === 'critical' ? 'text-red-700' : 'text-gray-700'}>
              {conflict.message}
            </span>
          </li>
        ))}
      </ul>
    </div>
  );
}

/**
 * GC content visualization - heatmap style
 */
function GCVisualization({ sequence, analysis, mutationPosition, height, onHover }) {
  if (!analysis) return null;

  const n = sequence.length;
  const barWidth = Math.max(1, 800 / n);

  return (
    <div className="border rounded-lg p-4 bg-white">
      <div className="text-sm font-medium text-gray-700 mb-2">
        GC Content Map ({(analysis.overallGC * 100).toFixed(1)}% overall)
      </div>
      <svg
        width="100%"
        height={height}
        viewBox={`0 0 800 ${height}`}
        preserveAspectRatio="none"
        className="border rounded"
        onMouseLeave={() => onHover(null)}
      >
        {/* Background */}
        <rect x="0" y="0" width="800" height={height} fill="#f9fafb" />

        {/* GC bars */}
        {analysis.perPosition.map((gc, i) => {
          const x = (i / n) * 800;
          const barHeight = gc * (height - 20);
          return (
            <rect
              key={i}
              x={x}
              y={height - 20 - barHeight}
              width={Math.max(barWidth, 1)}
              height={barHeight}
              fill={gcToColor(gc)}
              opacity={0.8}
              onMouseEnter={() => onHover(i)}
            />
          );
        })}

        {/* Optimal range indicator */}
        <rect x="0" y={height - 20 - 0.5 * (height - 20)} width="800" height="1" fill="#22c55e" opacity="0.5" />
        <text x="805" y={height - 20 - 0.5 * (height - 20) + 4} fontSize="8" fill="#22c55e">50%</text>

        {/* Mutation position marker */}
        {mutationPosition !== undefined && (
          <g>
            <line
              x1={(mutationPosition / n) * 800}
              y1="0"
              x2={(mutationPosition / n) * 800}
              y2={height - 20}
              stroke="#000"
              strokeWidth="2"
              strokeDasharray="4,2"
            />
            <polygon
              points={`${(mutationPosition / n) * 800 - 5},${height - 15} ${(mutationPosition / n) * 800 + 5},${height - 15} ${(mutationPosition / n) * 800},${height - 5}`}
              fill="#000"
            />
          </g>
        )}

        {/* X-axis */}
        <line x1="0" y1={height - 20} x2="800" y2={height - 20} stroke="#e5e7eb" />
        <text x="0" y={height - 5} fontSize="10" fill="#6b7280">1</text>
        <text x="790" y={height - 5} fontSize="10" fill="#6b7280" textAnchor="end">{n}</text>
      </svg>

      {/* Hotspots list */}
      {analysis.hotspots.length > 0 && (
        <div className="mt-2 text-xs text-gray-500">
          <span className="font-medium">Hotspots: </span>
          {analysis.hotspots.map((h, i) => (
            <span key={i} className="inline-block mr-2">
              <span className={`px-1 rounded ${h.severity === 'severe' ? 'bg-red-100 text-red-700' : 'bg-yellow-100 text-yellow-700'}`}>
                {h.start + 1}-{h.end + 1} ({h.type === 'high' ? 'High' : 'Low'} GC)
              </span>
            </span>
          ))}
        </div>
      )}
    </div>
  );
}

/**
 * Secondary structure visualization
 */
function StructureVisualization({ sequence, analysis, mutationPosition, height, onHover }) {
  if (!analysis) return null;

  const n = sequence.length;
  const barWidth = Math.max(1, 800 / n);

  return (
    <div className="border rounded-lg p-4 bg-white">
      <div className="text-sm font-medium text-gray-700 mb-2">
        Secondary Structure Stability (ΔG = {analysis.overallDG.toFixed(1)} kcal/mol)
      </div>
      <svg
        width="100%"
        height={height}
        viewBox={`0 0 800 ${height}`}
        preserveAspectRatio="none"
        className="border rounded"
        onMouseLeave={() => onHover(null)}
      >
        {/* Background */}
        <rect x="0" y="0" width="800" height={height} fill="#f9fafb" />

        {/* Zero line */}
        <line x1="0" y1="20" x2="800" y2="20" stroke="#d1d5db" strokeDasharray="4,4" />
        <text x="805" y="24" fontSize="8" fill="#9ca3af">0</text>

        {/* Structure bars (inverted - more negative = taller bar going down) */}
        {analysis.perPosition.map((dg, i) => {
          const x = (i / n) * 800;
          const barHeight = Math.min(Math.abs(dg) * 8, height - 40);  // Scale factor
          return (
            <rect
              key={i}
              x={x}
              y={20}
              width={Math.max(barWidth, 1)}
              height={barHeight}
              fill={dgToColor(dg)}
              opacity={0.8}
              onMouseEnter={() => onHover(i)}
            />
          );
        })}

        {/* Threshold lines */}
        <line x1="0" y1={20 + 3 * 8} x2="800" y2={20 + 3 * 8} stroke="#facc15" strokeWidth="1" opacity="0.5" />
        <line x1="0" y1={20 + 6 * 8} x2="800" y2={20 + 6 * 8} stroke="#f97316" strokeWidth="1" opacity="0.5" />
        <line x1="0" y1={20 + 8 * 8} x2="800" y2={20 + 8 * 8} stroke="#dc2626" strokeWidth="1" opacity="0.5" />

        {/* Mutation position marker */}
        {mutationPosition !== undefined && (
          <g>
            <line
              x1={(mutationPosition / n) * 800}
              y1="0"
              x2={(mutationPosition / n) * 800}
              y2={height - 20}
              stroke="#000"
              strokeWidth="2"
              strokeDasharray="4,2"
            />
            <polygon
              points={`${(mutationPosition / n) * 800 - 5},${height - 15} ${(mutationPosition / n) * 800 + 5},${height - 15} ${(mutationPosition / n) * 800},${height - 5}`}
              fill="#000"
            />
          </g>
        )}

        {/* X-axis */}
        <line x1="0" y1={height - 20} x2="800" y2={height - 20} stroke="#e5e7eb" />
        <text x="0" y={height - 5} fontSize="10" fill="#6b7280">1</text>
        <text x="790" y={height - 5} fontSize="10" fill="#6b7280" textAnchor="end">{n}</text>
      </svg>

      {/* Structure types */}
      {analysis.structures.length > 0 && (
        <div className="mt-2 text-xs text-gray-500">
          <span className="font-medium">Structures found: </span>
          {analysis.structures.filter(s => s.desc).map((s, i) => (
            <span key={i} className="inline-block mr-2">
              <span className="px-1 rounded bg-purple-100 text-purple-700">
                {s.desc} ({s.e.toFixed(1)} kcal/mol)
              </span>
            </span>
          )).slice(0, 5)}
        </div>
      )}
    </div>
  );
}

/**
 * Combined GC and Structure visualization
 */
function CombinedVisualization({ sequence, gcAnalysis, structureAnalysis, mutationPosition, height, onHover }) {
  if (!gcAnalysis || !structureAnalysis) return null;

  return (
    <div className="space-y-2">
      <GCVisualization
        sequence={sequence}
        analysis={gcAnalysis}
        mutationPosition={mutationPosition}
        height={height / 2 - 10}
        onHover={onHover}
      />
      <StructureVisualization
        sequence={sequence}
        analysis={structureAnalysis}
        mutationPosition={mutationPosition}
        height={height / 2 - 10}
        onHover={onHover}
      />
    </div>
  );
}

/**
 * Position tooltip showing details
 */
function PositionTooltip({ position, sequence, gcAnalysis, structureAnalysis }) {
  const gc = gcAnalysis?.perPosition[position];
  const dg = structureAnalysis?.perPosition[position];
  const base = sequence[position];

  return (
    <div className="absolute top-0 right-0 bg-white border rounded shadow-lg p-3 text-sm z-10">
      <div className="font-medium mb-1">Position {position + 1}</div>
      <div className="text-gray-600">Base: <span className="font-mono">{base}</span></div>
      {gc !== undefined && (
        <div className="text-gray-600">
          GC: <span style={{ color: gcToColor(gc) }}>{(gc * 100).toFixed(1)}%</span>
        </div>
      )}
      {dg !== undefined && (
        <div className="text-gray-600">
          ΔG: <span style={{ color: dgToColor(dg) }}>{dg.toFixed(1)} kcal/mol</span>
        </div>
      )}
    </div>
  );
}

/**
 * Suggestions panel
 */
function Suggestions({ suggestions }) {
  return (
    <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
      <div className="flex items-center gap-2 mb-2">
        <svg className="w-5 h-5 text-blue-500" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
        </svg>
        <span className="font-semibold text-blue-800">Suggestions</span>
      </div>
      <ul className="space-y-1 text-sm text-blue-700">
        {suggestions.map((s, i) => (
          <li key={i}>
            {s.message}
            {s.positions && (
              <span className="ml-1 text-xs">
                (positions: {s.positions.map(p => p + 1).join(', ')})
              </span>
            )}
          </li>
        ))}
      </ul>
    </div>
  );
}

/**
 * Legend component
 */
function Legend({ activeTab }) {
  return (
    <div className="flex flex-wrap gap-4 text-xs text-gray-600 pt-2 border-t">
      {(activeTab === 'gc' || activeTab === 'combined') && (
        <>
          <div className="flex items-center gap-1">
            <span className="w-3 h-3 rounded" style={{ backgroundColor: '#ef4444' }} />
            <span>High GC (≥70%)</span>
          </div>
          <div className="flex items-center gap-1">
            <span className="w-3 h-3 rounded" style={{ backgroundColor: '#22c55e' }} />
            <span>Optimal GC (40-60%)</span>
          </div>
          <div className="flex items-center gap-1">
            <span className="w-3 h-3 rounded" style={{ backgroundColor: '#3b82f6' }} />
            <span>Low GC (≤30%)</span>
          </div>
        </>
      )}
      {(activeTab === 'structure' || activeTab === 'combined') && (
        <>
          <div className="flex items-center gap-1">
            <span className="w-3 h-3 rounded" style={{ backgroundColor: '#dc2626' }} />
            <span>Strong hairpin (≤-8 kcal)</span>
          </div>
          <div className="flex items-center gap-1">
            <span className="w-3 h-3 rounded" style={{ backgroundColor: '#f97316' }} />
            <span>Moderate (≤-5 kcal)</span>
          </div>
          <div className="flex items-center gap-1">
            <span className="w-3 h-3 rounded" style={{ backgroundColor: '#22c55e' }} />
            <span>Stable (&gt;-3 kcal)</span>
          </div>
        </>
      )}
      <div className="flex items-center gap-1">
        <span className="w-3 h-0.5 bg-black" style={{ borderStyle: 'dashed' }} />
        <span>Mutation position</span>
      </div>
    </div>
  );
}

export default SequenceConflictMap;
