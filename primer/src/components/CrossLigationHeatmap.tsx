import { useState } from 'react';

// Types for heatmap data structure
interface RowStat {
  overhang: string;
  fidelity: number;
  fidelityPercent: string;
}

interface Hotspot {
  row: number;
  col: number;
  source: string;
  target: string;
  ratioPercent: string;
  severity: 'high' | 'medium' | 'low';
}

interface Stats {
  overallFidelity: number;
  overallFidelityPercent: string;
}

interface Visualization {
  title?: string;
  xAxisLabel?: string;
  yAxisLabel?: string;
}

interface HeatmapData {
  error?: string;
  overhangs: string[];
  matrix: number[][];
  normalizedMatrix: number[][];
  rowStats: RowStat[];
  hotspots: Hotspot[];
  stats: Stats;
  visualization: Visualization;
}

interface HoveredCell {
  row: number;
  col: number;
  value: number;
  normalized: number;
}

interface CrossLigationHeatmapProps {
  heatmapData: HeatmapData | null;
  enzyme?: string;
  compact?: boolean;
}

/**
 * CrossLigationHeatmap - Visualizes cross-ligation frequencies between overhangs
 *
 * Shows a color-coded matrix where:
 * - Diagonal (green) = correct ligations
 * - Off-diagonal (red) = cross-ligations (problematic)
 * - White = no ligation
 */
export function CrossLigationHeatmap({ heatmapData, enzyme = 'BsaI', compact = false }: CrossLigationHeatmapProps) {
  const [hoveredCell, setHoveredCell] = useState<HoveredCell | null>(null);
  const [showLegend] = useState<boolean>(!compact);

  if (!heatmapData || heatmapData.error) {
    return (
      <div className="heatmap-error">
        <span className="error-icon">!</span>
        <span>{heatmapData?.error || 'No heatmap data available'}</span>
      </div>
    );
  }

  const { overhangs, matrix, normalizedMatrix, rowStats, hotspots, stats, visualization } = heatmapData;

  // Ensure overhangs array exists
  if (!overhangs || !Array.isArray(overhangs) || overhangs.length === 0) {
    return (
      <div className="heatmap-error">
        <span className="error-icon">!</span>
        <span>No overhang data available for heatmap</span>
      </div>
    );
  }

  const n = overhangs.length;

  // Get color for a cell based on normalized value and position
  const getCellColor = (value: number, row: number, col: number, normalized: number): string => {
    if (value === 0) return '#f8f9fa'; // No ligation - light gray

    // Diagonal = correct ligation (green scale)
    if (row === col) {
      const intensity = Math.min(1, normalized);
      return `rgba(34, 197, 94, ${0.3 + intensity * 0.7})`; // Green
    }

    // Off-diagonal = cross-ligation (red scale)
    const intensity = Math.min(1, normalized * 5); // Amplify for visibility
    return `rgba(239, 68, 68, ${0.2 + intensity * 0.8})`; // Red
  };

  // Cell size based on compact mode and number of overhangs
  const cellSize = compact ? Math.min(24, 120 / n) : Math.min(40, 200 / n);
  const labelWidth = compact ? 40 : 50;

  return (
    <div className={`cross-ligation-heatmap ${compact ? 'compact' : ''}`}>
      {!compact && (
        <div className="heatmap-header">
          <h4>{visualization?.title || `Cross-Ligation Heatmap (${enzyme})`}</h4>
          <div className="heatmap-stats">
            <span className={`fidelity-badge ${(stats?.overallFidelity ?? 0) >= 0.95 ? 'excellent' : (stats?.overallFidelity ?? 0) >= 0.90 ? 'good' : 'warning'}`}>
              Fidelity: {stats?.overallFidelityPercent || 'N/A'}
            </span>
            {hotspots && hotspots.length > 0 && (
              <span className="hotspot-badge">
                {hotspots.length} hotspot{hotspots.length > 1 ? 's' : ''}
              </span>
            )}
          </div>
        </div>
      )}

      <div className="heatmap-container">
        {/* Y-axis label */}
        {!compact && (
          <div className="y-axis-label" style={{ width: labelWidth }}>
            <span>{visualization?.yAxisLabel || 'Source'}</span>
          </div>
        )}

        <div className="heatmap-grid-wrapper">
          {/* X-axis labels (top) */}
          <div className="x-axis-labels" style={{ marginLeft: labelWidth }}>
            {overhangs.map((oh, i) => (
              <div
                key={`x-${i}`}
                className="axis-label x-label"
                style={{ width: cellSize }}
                title={oh}
              >
                {compact ? oh.slice(0, 2) : oh}
              </div>
            ))}
          </div>

          <div className="heatmap-with-y-labels">
            {/* Y-axis labels + Grid */}
            <div className="y-axis-labels">
              {overhangs.map((oh, i) => (
                <div
                  key={`y-${i}`}
                  className="axis-label y-label"
                  style={{ height: cellSize, width: labelWidth }}
                  title={oh}
                >
                  {compact ? oh.slice(0, 2) : oh}
                </div>
              ))}
            </div>

            {/* Heatmap grid */}
            <div
              className="heatmap-grid"
              style={{
                gridTemplateColumns: `repeat(${n}, ${cellSize}px)`,
                gridTemplateRows: `repeat(${n}, ${cellSize}px)`,
              }}
            >
              {matrix.map((row, i) =>
                row.map((value, j) => {
                  const normalized = normalizedMatrix[i][j];
                  const isHovered = hoveredCell?.row === i && hoveredCell?.col === j;
                  const isDiagonal = i === j;
                  const isHotspot = hotspots?.some(h => h.row === i && h.col === j) ?? false;

                  return (
                    <div
                      key={`${i}-${j}`}
                      className={`heatmap-cell ${isDiagonal ? 'diagonal' : ''} ${isHotspot ? 'hotspot' : ''} ${isHovered ? 'hovered' : ''}`}
                      style={{
                        backgroundColor: getCellColor(value, i, j, normalized),
                        width: cellSize,
                        height: cellSize,
                      }}
                      onMouseEnter={() => setHoveredCell({ row: i, col: j, value, normalized })}
                      onMouseLeave={() => setHoveredCell(null)}
                    >
                      {!compact && value > 0 && cellSize >= 30 && (
                        <span className="cell-value">
                          {value >= 100 ? Math.round(value) : value.toFixed(0)}
                        </span>
                      )}
                    </div>
                  );
                })
              )}
            </div>
          </div>
        </div>

        {/* X-axis label (bottom) */}
        {!compact && (
          <div className="x-axis-label-bottom">
            <span>{visualization?.xAxisLabel || 'Target (RC)'}</span>
          </div>
        )}
      </div>

      {/* Tooltip */}
      {hoveredCell && (
        <div className="heatmap-tooltip">
          <div className="tooltip-row">
            <span className="tooltip-label">
              {overhangs[hoveredCell.row]} → {overhangs[hoveredCell.col]}
            </span>
          </div>
          <div className="tooltip-row">
            <span className="tooltip-type">
              {hoveredCell.row === hoveredCell.col ? 'Correct ligation' : 'Cross-ligation'}
            </span>
          </div>
          <div className="tooltip-row">
            <span className="tooltip-freq">Frequency: {hoveredCell.value}</span>
          </div>
        </div>
      )}

      {/* Legend */}
      {showLegend && (
        <div className="heatmap-legend">
          <div className="legend-item">
            <div className="legend-color" style={{ backgroundColor: 'rgba(34, 197, 94, 0.8)' }}></div>
            <span>Correct (diagonal)</span>
          </div>
          <div className="legend-item">
            <div className="legend-color" style={{ backgroundColor: 'rgba(239, 68, 68, 0.6)' }}></div>
            <span>Cross-ligation</span>
          </div>
          <div className="legend-item">
            <div className="legend-color" style={{ backgroundColor: '#f8f9fa', border: '1px solid #e5e7eb' }}></div>
            <span>No ligation</span>
          </div>
        </div>
      )}

      {/* Junction Fidelity Bars */}
      {!compact && rowStats && (
        <div className="junction-fidelity-section">
          <h5>Per-Junction Fidelity</h5>
          <div className="fidelity-bars">
            {rowStats.map((stat, i) => (
              <div key={i} className="fidelity-bar-row">
                <span className="fidelity-bar-label">{stat.overhang}</span>
                <div className="fidelity-bar-track">
                  <div
                    className={`fidelity-bar-fill ${stat.fidelity >= 0.99 ? 'excellent' : stat.fidelity >= 0.95 ? 'good' : stat.fidelity >= 0.90 ? 'medium' : 'low'}`}
                    style={{ width: `${stat.fidelity * 100}%` }}
                  ></div>
                </div>
                <span className="fidelity-bar-value">{stat.fidelityPercent}</span>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Hotspots Warning */}
      {!compact && hotspots && hotspots.length > 0 && (
        <div className="hotspots-section">
          <h5>Cross-Ligation Hotspots</h5>
          <div className="hotspots-list">
            {hotspots.slice(0, 5).map((hotspot, i) => (
              <div key={i} className={`hotspot-item ${hotspot.severity}`}>
                <span className="hotspot-pair">
                  {hotspot.source} → {hotspot.target}
                </span>
                <span className="hotspot-ratio">{hotspot.ratioPercent} of correct</span>
                <span className={`hotspot-severity ${hotspot.severity}`}>
                  {hotspot.severity}
                </span>
              </div>
            ))}
            {hotspots.length > 5 && (
              <div className="hotspots-more">
                +{hotspots.length - 5} more
              </div>
            )}
          </div>
        </div>
      )}

      <style>{`
        .cross-ligation-heatmap {
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
          padding: 16px;
          background: white;
          border-radius: 8px;
          border: 1px solid #e5e7eb;
        }

        .cross-ligation-heatmap.compact {
          padding: 8px;
        }

        .heatmap-error {
          display: flex;
          align-items: center;
          gap: 8px;
          padding: 12px;
          background: #fef2f2;
          border: 1px solid #fecaca;
          border-radius: 6px;
          color: #dc2626;
        }

        .heatmap-header {
          display: flex;
          justify-content: space-between;
          align-items: center;
          margin-bottom: 16px;
        }

        .heatmap-header h4 {
          margin: 0;
          font-size: 14px;
          font-weight: 600;
          color: #1f2937;
        }

        .heatmap-stats {
          display: flex;
          gap: 8px;
        }

        .fidelity-badge, .hotspot-badge {
          padding: 4px 8px;
          border-radius: 4px;
          font-size: 12px;
          font-weight: 500;
        }

        .fidelity-badge.excellent { background: #dcfce7; color: #166534; }
        .fidelity-badge.good { background: #fef3c7; color: #92400e; }
        .fidelity-badge.warning { background: #fecaca; color: #dc2626; }

        .hotspot-badge {
          background: #fef3c7;
          color: #92400e;
        }

        .heatmap-container {
          display: flex;
          flex-direction: column;
          align-items: center;
        }

        .y-axis-label {
          writing-mode: vertical-rl;
          text-orientation: mixed;
          transform: rotate(180deg);
          font-size: 11px;
          color: #6b7280;
          display: flex;
          align-items: center;
          justify-content: center;
        }

        .heatmap-grid-wrapper {
          display: flex;
          flex-direction: column;
        }

        .x-axis-labels {
          display: flex;
        }

        .axis-label {
          font-size: 10px;
          font-family: 'Monaco', 'Menlo', monospace;
          color: #374151;
          display: flex;
          align-items: center;
          justify-content: center;
          overflow: hidden;
          text-overflow: ellipsis;
        }

        .x-label {
          height: 20px;
        }

        .heatmap-with-y-labels {
          display: flex;
        }

        .y-axis-labels {
          display: flex;
          flex-direction: column;
        }

        .y-label {
          justify-content: flex-end;
          padding-right: 6px;
        }

        .heatmap-grid {
          display: grid;
          gap: 1px;
          background: #e5e7eb;
          border: 1px solid #d1d5db;
          border-radius: 4px;
          overflow: hidden;
        }

        .heatmap-cell {
          display: flex;
          align-items: center;
          justify-content: center;
          cursor: pointer;
          transition: all 0.15s ease;
          position: relative;
        }

        .heatmap-cell.hovered {
          outline: 2px solid #3b82f6;
          outline-offset: -1px;
          z-index: 1;
        }

        .heatmap-cell.diagonal {
          border: 1px solid rgba(34, 197, 94, 0.5);
        }

        .heatmap-cell.hotspot:not(.diagonal) {
          border: 1px solid rgba(239, 68, 68, 0.8);
        }

        .cell-value {
          font-size: 8px;
          color: #374151;
          font-weight: 500;
        }

        .x-axis-label-bottom {
          margin-top: 8px;
          font-size: 11px;
          color: #6b7280;
        }

        .heatmap-tooltip {
          position: fixed;
          background: rgba(0, 0, 0, 0.85);
          color: white;
          padding: 8px 12px;
          border-radius: 6px;
          font-size: 12px;
          pointer-events: none;
          z-index: 1000;
          transform: translate(-50%, -100%);
          margin-top: -10px;
        }

        .tooltip-row {
          margin: 2px 0;
        }

        .tooltip-label {
          font-family: 'Monaco', monospace;
          font-weight: 600;
        }

        .tooltip-type {
          color: #9ca3af;
        }

        .heatmap-legend {
          display: flex;
          justify-content: center;
          gap: 16px;
          margin-top: 12px;
          padding-top: 12px;
          border-top: 1px solid #e5e7eb;
        }

        .legend-item {
          display: flex;
          align-items: center;
          gap: 6px;
          font-size: 11px;
          color: #6b7280;
        }

        .legend-color {
          width: 16px;
          height: 16px;
          border-radius: 3px;
        }

        .junction-fidelity-section {
          margin-top: 16px;
          padding-top: 16px;
          border-top: 1px solid #e5e7eb;
        }

        .junction-fidelity-section h5 {
          margin: 0 0 12px 0;
          font-size: 12px;
          font-weight: 600;
          color: #374151;
        }

        .fidelity-bars {
          display: flex;
          flex-direction: column;
          gap: 6px;
        }

        .fidelity-bar-row {
          display: flex;
          align-items: center;
          gap: 8px;
        }

        .fidelity-bar-label {
          width: 50px;
          font-size: 10px;
          font-family: 'Monaco', monospace;
          color: #374151;
          text-align: right;
        }

        .fidelity-bar-track {
          flex: 1;
          height: 8px;
          background: #f3f4f6;
          border-radius: 4px;
          overflow: hidden;
        }

        .fidelity-bar-fill {
          height: 100%;
          border-radius: 4px;
          transition: width 0.3s ease;
        }

        .fidelity-bar-fill.excellent { background: #22c55e; }
        .fidelity-bar-fill.good { background: #84cc16; }
        .fidelity-bar-fill.medium { background: #f59e0b; }
        .fidelity-bar-fill.low { background: #ef4444; }

        .fidelity-bar-value {
          width: 45px;
          font-size: 10px;
          color: #6b7280;
          text-align: right;
        }

        .hotspots-section {
          margin-top: 16px;
          padding-top: 16px;
          border-top: 1px solid #e5e7eb;
        }

        .hotspots-section h5 {
          margin: 0 0 12px 0;
          font-size: 12px;
          font-weight: 600;
          color: #374151;
        }

        .hotspots-list {
          display: flex;
          flex-direction: column;
          gap: 6px;
        }

        .hotspot-item {
          display: flex;
          align-items: center;
          gap: 12px;
          padding: 6px 10px;
          background: #fef2f2;
          border-radius: 4px;
          font-size: 11px;
        }

        .hotspot-pair {
          font-family: 'Monaco', monospace;
          font-weight: 500;
          color: #1f2937;
        }

        .hotspot-ratio {
          color: #6b7280;
        }

        .hotspot-severity {
          padding: 2px 6px;
          border-radius: 3px;
          font-size: 10px;
          font-weight: 500;
        }

        .hotspot-severity.high { background: #dc2626; color: white; }
        .hotspot-severity.medium { background: #f59e0b; color: white; }
        .hotspot-severity.low { background: #6b7280; color: white; }

        .hotspots-more {
          font-size: 11px;
          color: #6b7280;
          text-align: center;
          padding: 4px;
        }
      `}</style>
    </div>
  );
}

interface CrossLigationHeatmapCompactProps {
  heatmapData: HeatmapData | null;
  enzyme?: string;
}

/**
 * Compact version for embedding in other components
 */
export function CrossLigationHeatmapCompact({ heatmapData, enzyme }: CrossLigationHeatmapCompactProps) {
  return <CrossLigationHeatmap heatmapData={heatmapData} enzyme={enzyme} compact={true} />;
}

export default CrossLigationHeatmap;
