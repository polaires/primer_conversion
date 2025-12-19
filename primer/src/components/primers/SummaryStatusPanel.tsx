/**
 * SummaryStatusPanel - Quick status overview for primer pair quality
 *
 * Displays pass/warning/fail status badges for key metrics:
 * - Primer Quality (composite score)
 * - Tm Difference
 * - GC Clamp
 * - Heterodimer
 * - Secondary Structure
 * - Self-Dimer
 * - Off-Targets
 *
 * Extracted from UnifiedPrimerDesigner for reuse in:
 * - UnifiedPrimerDesigner
 * - EnhancedScorer (standalone scoring mode)
 */

import { useMemo } from 'react';

// Evidence-based thresholds (from IDT, Premier Biosoft)
const DIMER_THRESHOLDS = {
  hairpin: { internal: { ideal: -3.0, critical: -6.0 } },
  selfDimer: { internal: { ideal: -6.0, critical: -8.0 } },
};

type StatusType = 'pass' | 'warning' | 'fail';
type QualityTier = 'excellent' | 'good' | 'acceptable' | 'marginal' | 'poor';
type SectionType = 'primers' | 'pcr' | 'analysis';

interface StatusCheck {
  label: string;
  status: StatusType;
  value: string;
  section: SectionType;
  threshold: string;
}

interface PrimerData {
  tm?: number;
  gc?: number;
  sequence?: string;
}

interface ThermodynamicsData {
  hairpinDG?: number;
  homodimerDG?: number;
}

interface PrimerAnalysis {
  hairpinDG?: number;
  selfDimerDG?: number;
  thermodynamics?: ThermodynamicsData;
}

interface HeterodimerData {
  isExpectedOverlap?: boolean;
  severity?: string;
  dimerType?: string;
  overlapLength?: number;
  heterodimerDG?: number;
}

interface PairData {
  heterodimer?: HeterodimerData;
}

interface OffTargetData {
  offTargetCount?: number;
}

interface OffTargetsData {
  forward?: OffTargetData;
  reverse?: OffTargetData;
}

interface AnalysisData {
  heterodimer?: HeterodimerData;
  pair?: PairData;
  forward?: PrimerAnalysis;
  reverse?: PrimerAnalysis;
  offTargets?: OffTargetsData;
}

interface SummaryStatusPanelProps {
  forward: PrimerData | null;
  reverse?: PrimerData | null;
  analysis?: AnalysisData | null;
  quality?: QualityTier | string;
  onSectionClick?: (section: SectionType) => void;
}

interface SummaryStatusResult {
  checks: Record<string, StatusCheck>;
  overall: StatusType;
}

export default function SummaryStatusPanel({
  forward,
  reverse,
  analysis,
  quality,
  onSectionClick,
}: SummaryStatusPanelProps) {
  const summaryStatus = useMemo<SummaryStatusResult | null>(() => {
    if (!forward) return null;

    const checks: Record<string, StatusCheck> = {};

    // Quality check
    checks.quality = {
      label: 'Primer Quality',
      status:
        quality === 'excellent' || quality === 'good'
          ? 'pass'
          : quality === 'acceptable'
          ? 'warning'
          : 'fail',
      value: quality || 'N/A',
      section: 'primers',
      threshold:
        'Excellent/Good = optimal, Acceptable = may need optimization, Poor = redesign recommended',
    };

    // Tm Difference check (if we have both primers)
    if (reverse) {
      const tmDiff = Math.abs((forward.tm || 0) - (reverse.tm || 0));
      checks.tmDiff = {
        label: 'Tm Difference',
        status: tmDiff <= 2 ? 'pass' : tmDiff <= 5 ? 'warning' : 'fail',
        value: `${tmDiff.toFixed(1)}°C`,
        section: 'pcr',
        threshold:
          '≤2°C = ideal, 2-5°C = acceptable, >5°C = may cause uneven amplification',
      };
    }

    // GC Clamp check
    const fwdHasGCClamp =
      forward.sequence &&
      (forward.sequence[forward.sequence.length - 1] === 'G' ||
        forward.sequence[forward.sequence.length - 1] === 'C');
    const revHasGCClamp =
      reverse?.sequence &&
      (reverse.sequence[reverse.sequence.length - 1] === 'G' ||
        reverse.sequence[reverse.sequence.length - 1] === 'C');

    checks.gcClamp = {
      label: 'GC Clamp',
      status: fwdHasGCClamp && (!reverse || revHasGCClamp) ? 'pass' : 'warning',
      value:
        fwdHasGCClamp && (!reverse || revHasGCClamp)
          ? reverse
            ? 'Both'
            : 'Yes'
          : fwdHasGCClamp || revHasGCClamp
          ? 'Partial'
          : 'None',
      section: 'primers',
      threshold:
        "G or C at 3' end improves specificity. Both primers having GC clamp is ideal.",
    };

    // Enhanced analysis checks
    if (analysis) {
      const heterodimer = analysis.heterodimer || analysis.pair?.heterodimer || {};
      const fwdAnalysis = analysis.forward || {};
      const revAnalysis = analysis.reverse || {};
      const offTargets = analysis.offTargets || {};

      // Heterodimer check
      const isExpectedOverlap = heterodimer.isExpectedOverlap;
      const heteroSeverity = heterodimer.severity || 'safe';
      const heterodimerStatus: StatusType = isExpectedOverlap
        ? 'pass'
        : heteroSeverity === 'critical'
        ? 'fail'
        : heteroSeverity === 'warning'
        ? 'warning'
        : 'pass';
      const heterodimerType = heterodimer.dimerType || 'none';
      const heteroLabel = isExpectedOverlap
        ? 'Primer Overlap (QuikChange)'
        : heterodimerType === '3prime_extensible'
        ? "3' Extensible Dimer"
        : heterodimerType === 'internal_with_3prime'
        ? "Heterodimer (3' involved)"
        : 'Heterodimer';
      const heteroValue = isExpectedOverlap
        ? `${heterodimer.overlapLength || 'Full'} bp overlap`
        : (heterodimer.heterodimerDG?.toFixed(1) || 'N/A') + ' kcal/mol';
      const heteroThreshold = isExpectedOverlap
        ? 'Full primer complementarity is expected for QuikChange-style mutagenesis'
        : heterodimerType === '3prime_extensible'
        ? "CRITICAL: Both primers' 3' ends hybridize - causes primer-dimer artifacts"
        : ">-6 kcal/mol = OK. 3' end involvement is more critical than ΔG alone.";

      checks.heterodimer = {
        label: heteroLabel,
        status: heterodimerStatus,
        value: heteroValue,
        section: 'analysis',
        threshold: heteroThreshold,
      };

      // Hairpin/Secondary Structure check
      const fwdHairpin = fwdAnalysis.hairpinDG || fwdAnalysis.thermodynamics?.hairpinDG;
      const revHairpin = revAnalysis.hairpinDG || revAnalysis.thermodynamics?.hairpinDG;
      const worstHairpin = Math.min(fwdHairpin || 0, revHairpin || 0);
      const hairpinThresh = DIMER_THRESHOLDS.hairpin.internal;

      checks.secondaryStructure = {
        label: 'Hairpin',
        status:
          worstHairpin >= hairpinThresh.ideal
            ? 'pass'
            : worstHairpin >= hairpinThresh.critical
            ? 'warning'
            : 'fail',
        value: `${worstHairpin?.toFixed(1) || 'N/A'} kcal/mol`,
        section: 'analysis',
        threshold: `>${hairpinThresh.ideal} = ideal (Premier Biosoft), ${hairpinThresh.ideal} to ${hairpinThresh.critical} = monitor, <${hairpinThresh.critical} = problematic`,
      };

      // Self-dimer check
      const fwdSelfDimer = fwdAnalysis.selfDimerDG || fwdAnalysis.thermodynamics?.homodimerDG;
      const revSelfDimer = revAnalysis.selfDimerDG || revAnalysis.thermodynamics?.homodimerDG;
      const worstSelfDimer = Math.min(fwdSelfDimer || 0, revSelfDimer || 0);
      const selfDimerThresh = DIMER_THRESHOLDS.selfDimer.internal;

      checks.selfDimer = {
        label: 'Self-Dimer',
        status:
          worstSelfDimer >= selfDimerThresh.ideal
            ? 'pass'
            : worstSelfDimer >= selfDimerThresh.critical
            ? 'warning'
            : 'fail',
        value: `${worstSelfDimer?.toFixed(1) || 'N/A'} kcal/mol`,
        section: 'analysis',
        threshold: `>${selfDimerThresh.ideal} = OK (Premier Biosoft), ${selfDimerThresh.ideal} to ${selfDimerThresh.critical} = monitor, <${selfDimerThresh.critical} = IDT critical threshold`,
      };

      // Off-targets check
      if (offTargets.forward || offTargets.reverse) {
        const maxOffTargets = Math.max(
          offTargets.forward?.offTargetCount || 0,
          offTargets.reverse?.offTargetCount || 0
        );
        checks.offTargets = {
          label: 'Off-Targets',
          status: maxOffTargets === 0 ? 'pass' : maxOffTargets <= 2 ? 'warning' : 'fail',
          value: maxOffTargets.toString(),
          section: 'analysis',
          threshold:
            '0 = specific binding only, 1-2 = minor risk of non-specific products, >2 = may need redesign',
        };
      }
    }

    // Calculate overall status
    const statuses = Object.values(checks).map((c) => c.status);
    const overall: StatusType = statuses.includes('fail')
      ? 'fail'
      : statuses.includes('warning')
      ? 'warning'
      : 'pass';

    return { checks, overall };
  }, [forward, reverse, analysis, quality]);

  if (!summaryStatus) {
    return null;
  }

  const statusColors: Record<StatusType, string> = {
    pass: '#22c55e',
    warning: '#eab308',
    fail: '#ef4444',
  };

  const statusIcons: Record<StatusType, string> = {
    pass: '✓',
    warning: '⚠',
    fail: '✕',
  };

  return (
    <div className="summary-status-panel">
      <div className="summary-status-header">
        <span
          className={`overall-status ${summaryStatus.overall}`}
          style={{
            backgroundColor: statusColors[summaryStatus.overall] + '20',
            color: statusColors[summaryStatus.overall],
            padding: '4px 12px',
            borderRadius: '4px',
            fontWeight: 'bold',
            fontSize: '13px',
          }}
        >
          {statusIcons[summaryStatus.overall]} Overall:{' '}
          {summaryStatus.overall === 'pass'
            ? 'Good'
            : summaryStatus.overall === 'warning'
            ? 'Needs Review'
            : 'Issues Found'}
        </span>
      </div>

      <div
        className="summary-status-grid"
        style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fit, minmax(140px, 1fr))',
          gap: '8px',
          marginTop: '12px',
        }}
      >
        {Object.entries(summaryStatus.checks).map(([key, check]) => (
          <div
            key={key}
            className={`status-item ${check.status}`}
            style={{
              padding: '8px 12px',
              borderRadius: '6px',
              backgroundColor: statusColors[check.status] + '10',
              border: `1px solid ${statusColors[check.status]}40`,
              cursor: onSectionClick ? 'pointer' : 'default',
            }}
            onClick={() => onSectionClick?.(check.section)}
            title={check.threshold}
          >
            <div
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: '6px',
                marginBottom: '4px',
              }}
            >
              <span
                style={{
                  color: statusColors[check.status],
                  fontWeight: 'bold',
                  fontSize: '12px',
                }}
              >
                {statusIcons[check.status]}
              </span>
              <span
                style={{
                  fontSize: '11px',
                  color: '#64748b',
                  fontWeight: '500',
                }}
              >
                {check.label}
              </span>
            </div>
            <div
              style={{
                fontSize: '13px',
                fontWeight: '600',
                color: '#1e293b',
                textTransform: 'capitalize',
              }}
            >
              {check.value}
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}
