/**
 * Circular Plasmid View
 * Circular visualization of assembled plasmid with fragments
 */

import React, { useMemo } from 'react';
import { Fragment, FRAGMENT_TYPES } from '../../../types/fragmentPlanner';

interface CircularPlasmidViewProps {
  fragments: Fragment[];
  overhangs: string[];
  totalLength: number;
}

export function CircularPlasmidView({
  fragments,
  overhangs,
  totalLength,
}: CircularPlasmidViewProps) {
  const size = 320;
  const center = size / 2;
  const radius = 120;
  const innerRadius = 80;
  const overhangRadius = radius + 30;

  // Calculate arc positions for each fragment
  const fragmentArcs = useMemo(() => {
    if (fragments.length === 0 || totalLength === 0) return [];

    let currentAngle = -90; // Start at 12 o'clock

    return fragments.map((fragment, idx) => {
      const sweepAngle = (fragment.length / totalLength) * 360;
      const startAngle = currentAngle;
      const endAngle = currentAngle + sweepAngle;
      const midAngle = startAngle + sweepAngle / 2;

      currentAngle = endAngle;

      return {
        fragment,
        startAngle,
        endAngle,
        midAngle,
        sweepAngle,
      };
    });
  }, [fragments, totalLength]);

  // Convert angle to SVG arc coordinates
  const polarToCartesian = (cx: number, cy: number, r: number, angle: number) => {
    const angleRad = (angle * Math.PI) / 180;
    return {
      x: cx + r * Math.cos(angleRad),
      y: cy + r * Math.sin(angleRad),
    };
  };

  // Create arc path
  const describeArc = (
    cx: number,
    cy: number,
    r: number,
    startAngle: number,
    endAngle: number
  ) => {
    const start = polarToCartesian(cx, cy, r, endAngle);
    const end = polarToCartesian(cx, cy, r, startAngle);
    const largeArcFlag = endAngle - startAngle <= 180 ? 0 : 1;

    return `M ${start.x} ${start.y} A ${r} ${r} 0 ${largeArcFlag} 0 ${end.x} ${end.y}`;
  };

  // Create donut segment path
  const describeDonutSegment = (
    cx: number,
    cy: number,
    innerR: number,
    outerR: number,
    startAngle: number,
    endAngle: number
  ) => {
    const outerStart = polarToCartesian(cx, cy, outerR, startAngle);
    const outerEnd = polarToCartesian(cx, cy, outerR, endAngle);
    const innerStart = polarToCartesian(cx, cy, innerR, startAngle);
    const innerEnd = polarToCartesian(cx, cy, innerR, endAngle);
    const largeArcFlag = endAngle - startAngle <= 180 ? 0 : 1;

    return `
      M ${outerStart.x} ${outerStart.y}
      A ${outerR} ${outerR} 0 ${largeArcFlag} 1 ${outerEnd.x} ${outerEnd.y}
      L ${innerEnd.x} ${innerEnd.y}
      A ${innerR} ${innerR} 0 ${largeArcFlag} 0 ${innerStart.x} ${innerStart.y}
      Z
    `;
  };

  if (fragments.length === 0) {
    return (
      <div className="h-full flex items-center justify-center">
        <p className="text-slate-400 dark:text-slate-500">No fragments to display</p>
      </div>
    );
  }

  return (
    <div className="h-full flex items-center justify-center">
      <svg
        viewBox={`0 0 ${size} ${size}`}
        className="w-full max-w-md h-auto"
      >
        {/* Background circle */}
        <circle
          cx={center}
          cy={center}
          r={radius}
          fill="none"
          stroke="#e2e8f0"
          strokeWidth="1"
          className="dark:stroke-slate-700"
        />

        {/* Fragment segments */}
        {fragmentArcs.map(({ fragment, startAngle, endAngle, midAngle }, idx) => {
          const typeInfo = FRAGMENT_TYPES[fragment.type];
          const labelPos = polarToCartesian(center, center, radius + 45, midAngle);
          const overhangPos = polarToCartesian(center, center, overhangRadius, startAngle);

          return (
            <g key={fragment.id}>
              {/* Segment */}
              <path
                d={describeDonutSegment(
                  center,
                  center,
                  innerRadius,
                  radius,
                  startAngle,
                  endAngle
                )}
                fill={typeInfo.color}
                fillOpacity={0.8}
                stroke="white"
                strokeWidth="2"
                className="cursor-pointer hover:fill-opacity-100 transition-all"
              />

              {/* Junction marker (overhang) */}
              {overhangs[idx] && (
                <g>
                  <circle
                    cx={overhangPos.x}
                    cy={overhangPos.y}
                    r="14"
                    fill="white"
                    stroke="#94a3b8"
                    strokeWidth="1"
                    className="dark:fill-slate-800 dark:stroke-slate-600"
                  />
                  <text
                    x={overhangPos.x}
                    y={overhangPos.y}
                    textAnchor="middle"
                    dominantBaseline="central"
                    className="text-[8px] font-mono font-medium fill-slate-700 dark:fill-slate-300"
                  >
                    {overhangs[idx]}
                  </text>
                </g>
              )}

              {/* Fragment label */}
              <text
                x={labelPos.x}
                y={labelPos.y}
                textAnchor="middle"
                dominantBaseline="central"
                className="text-xs font-medium fill-slate-700 dark:fill-slate-300 pointer-events-none"
              >
                {fragment.name.length > 12
                  ? fragment.name.slice(0, 10) + '...'
                  : fragment.name}
              </text>
            </g>
          );
        })}

        {/* Center info */}
        <g>
          <circle
            cx={center}
            cy={center}
            r={innerRadius - 10}
            fill="white"
            className="dark:fill-slate-800"
          />
          <text
            x={center}
            y={center - 10}
            textAnchor="middle"
            className="text-lg font-bold fill-slate-900 dark:fill-white"
          >
            {(totalLength / 1000).toFixed(1)} kb
          </text>
          <text
            x={center}
            y={center + 10}
            textAnchor="middle"
            className="text-xs fill-slate-500 dark:fill-slate-400"
          >
            {fragments.length} fragments
          </text>
        </g>

        {/* Legend */}
        <g transform={`translate(10, ${size - 70})`}>
          {fragments.slice(0, 4).map((fragment, idx) => {
            const typeInfo = FRAGMENT_TYPES[fragment.type];
            return (
              <g key={fragment.id} transform={`translate(0, ${idx * 16})`}>
                <rect
                  width="10"
                  height="10"
                  rx="2"
                  fill={typeInfo.color}
                />
                <text
                  x="14"
                  y="8"
                  className="text-[10px] fill-slate-600 dark:fill-slate-400"
                >
                  {fragment.name.length > 15
                    ? fragment.name.slice(0, 13) + '...'
                    : fragment.name}
                </text>
              </g>
            );
          })}
          {fragments.length > 4 && (
            <text
              y={4 * 16 + 8}
              className="text-[10px] fill-slate-500 dark:fill-slate-500"
            >
              +{fragments.length - 4} more
            </text>
          )}
        </g>
      </svg>
    </div>
  );
}

export default CircularPlasmidView;
