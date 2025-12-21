import React, { useState, useCallback, useMemo, useEffect } from 'react';
import { SeqViz } from 'seqviz';
import {
  GOLDEN_GATE_ENZYMES,
  ENZYMES_WITH_DATA,
  parseFasta,
  // Experimental fidelity functions
  calculateExperimentalFidelity,
  compareEnzymeFidelity,
  findProblematicPairs,
  // Optimized Golden Gate primer design (UI-compatible adapter)
  designOptimizedGoldenGateAssemblyForUI,
  // Fusion Site Optimizer
  optimizeFusionSites,
  scanForFusionSites,
  scoreFusionSiteComposite,
  OPTIMIZER_DEFAULTS,
  // Auto-Domestication Optimizer (legacy - kept for compatibility)
  optimizeWithDomestication,
  analyzeForDomestication,
  DOMESTICATION_DEFAULTS,
  // Unified Domestication Optimizer (PREFERRED - one-pot compatible strategies)
  analyzeDomesticationOptions,
  DOMESTICATION_STRATEGY,
  designEnhancedMutagenicJunction,
} from '../lib/repp/index.js';

import {
  optimizeAssemblyBoundaries,
  assessBoundaryOptimizationPotential,
} from '../lib/repp/primer-boundary-optimizer.js';

import { CrossLigationHeatmapCompact } from './CrossLigationHeatmap';

import FusionSiteOptimizerPanel from './FusionSiteOptimizerPanel';
import { EnhancedDomesticationPanel, EnhancedDomesticationStyles } from './EnhancedDomesticationPanel';
import { DomesticationWorkflowGuide, DomesticationWorkflowStyles } from './DomesticationWorkflowGuide';
import {
  designNEBuilderAssembly,
  NEBUILDER_PARAMS,
} from '../lib/nebuilder.js';
import {
  simulateAssembly,
  exportToGenBank,
  exportToFasta,
  exportProject,
  importProject,
} from '../lib/assemblyCore.js';
import IsothermalAssemblyPanel from './IsothermalAssemblyPanel';
import {
  findOptimalOverhangSet,
  getLigationFrequency,
  getEnzymeLigationData,
  getRecommendedOverhangs,
  generateCrossLigationHeatmap,
} from '../lib/repp/goldengate.js';
import { reverseComplement } from '../lib/repp/enzymes.js';

import { designAllMutagenicJunctions } from '../lib/repp/mutagenic-junction-domesticator.js';

// TypeScript Interfaces
interface Part {
  id: string;
  seq: string;
  type: 'promoter' | 'rbs' | 'cds' | 'terminator' | 'other' | 'fragment';
  rangeStart?: number;
  rangeEnd?: number;
  originalSeq?: string;
  index?: number;
  length?: number;
  _originalSequence?: string;
  leftOverhang?: string;
  rightOverhang?: string;
  primers?: {
    forward: {
      sequence: string;
      length: number;
      tm: number;
      gc: number;
      qualityScore?: number;
      qualityTier?: string;
      breakdown?: any;
      structure?: {
        extra?: string;
        bsaISite?: string;
        recognitionSite?: string;
        spacer?: string;
        overhang?: string;
        homology?: string;
        annealing?: string;
      };
    };
    reverse: {
      sequence: string;
      length: number;
      tm: number;
      gc: number;
      qualityScore?: number;
      qualityTier?: string;
      breakdown?: any;
      structure?: {
        extra?: string;
        bsaISite?: string;
        recognitionSite?: string;
        spacer?: string;
        overhang?: string;
        homology?: string;
        annealing?: string;
      };
    };
    pairQuality?: {
      tmDifference?: { value: number };
      heterodimer?: any;
      score?: number;
      quality?: any;
    };
    pcr?: {
      annealingTemp: number;
      extensionTime: number;
      tmDifference?: number;
    };
  };
  _domesticationApproved?: boolean;
  _domesticationStrategy?: string;
  _domesticationJunctions?: Array<{
    sitePosition?: number;
    junctionPosition?: number;
    position?: number;
    overhang: string;
    site?: { position: number };
  }>;
  _domesticationMutations?: Array<{
    position: number;
    original: string;
    mutated: string;
  }>;
  _domesticationPrimers?: any[];
  _workflowApplied?: boolean;
  _workflowResult?: any;
  _optimized?: {
    leftOverhang: string | null;
    rightOverhang: string | null;
    junctionQuality: number | null;
  };
  _domesticated?: boolean;
  _userConfigured?: boolean;
  _parentPart?: string;
  _fragmentIndex?: number;
  _totalFragments?: number;
  _overhang?: string | null;
  _strategy?: string;
  _onePotCompatible?: boolean;
  _configuredPrimers?: any[];
}


interface PartWithAngles extends Part {
  startAngle: number;
  sweepAngle: number;
  midAngle: number;
  index: number;
}

interface AssemblyMethod {
  id: string;
  name: string;
  description: string;
  icon: React.ReactElement;
  color: string;
}

interface IsothermalVariant {
  id: string;
  name: string;
  description: string;
  badge: string;
  optimalOverlap: { min: number; max: number; ideal: number };
  optimalTm: { min: number; max: number; ideal: number };
  recommended?: boolean;
}

interface PartTypeInfo {
  name: string;
  color: string;
  icon: string;
}

interface OverlapSettings {
  overlapLength: number;
  targetTm: number;
}

interface DomesticationAnalysis {
  error?: boolean | { type?: string; details?: any[]; violations?: any[] };
  message?: string;
  needsDomestication?: boolean;
  totalFragments?: number;
  additionalFragmentsNeeded?: number;
  domesticationJunctions?: any[];
  alternativeEnzymes?: any[];
  recommendedEnzyme?: string;
}

interface FusionCandidate {
  position: number;
  overhang: string;
  score: any;
}

interface AssemblyResult {
  parts: any[];
  overhangs?: string[];
  assembledSequence?: string;
  assembledLength?: number;
  fidelity?: {
    percentage: string;
    individual?: any[];
    baseFidelityPercent?: string;
    gtAdjusted?: number;
  };
  warnings?: string[];
  primers?: any[];
  simulation?: any;
  overlapAnalysis?: any;
  protocol?: string;
  hasInternalSites?: boolean;
  internalSiteIssues?: any[];
  _autoDomestication?: {
    applied: boolean;
    details: any[];
  };
  _optimizedData?: any;
}

interface PartCardProps {
  part: Part;
  index: number;
  onChange: (index: number, newPart: Part) => void;
  onRemove: (index: number) => void;
  canRemove: boolean;
  enzyme: string;
  onDragStart: (e: React.DragEvent, index: number) => void;
  onDragOver: (e: React.DragEvent, index: number) => void;
  onDrop: (e: React.DragEvent, index: number) => void;
  isDragging: boolean;
  onSplitSequence?: (index: number, sequence: string, numFragments: number, enzyme: string) => Promise<any>;
  useEnhancedWorkflow?: boolean;
  onOpenEnhancedDomestication?: (part: Part & { index: number }) => void;
  onOpenWorkflowGuide?: (part: Part) => void;
}

interface CircularPlasmidViewProps {
  parts: Part[];
  overhangs: string[];
  totalLength: number;
}

interface LinearAssemblyDiagramProps {
  parts: Part[];
  overhangs: string[];
}

interface FidelityGaugeProps {
  overhangs: string[];
  enzyme: string;
  staticFidelity: number;
}

interface PrimerResultsProps {
  result: AssemblyResult;
  onCopy?: (message: string) => void;
  method: string;
  enzyme: string;
  isGoldenGate?: boolean;
}

interface ExperimentalFidelity {
  assemblyFidelity: number;
  source: string;
}

interface EnzymeComparison {
  ranked: Array<{
    enzyme: string;
    assemblyFidelity: number;
    [key: string]: any;
  }>;
}

interface HeatmapData {
  error?: boolean;
  [key: string]: any;
}

interface CrossLigationData {
  error?: boolean;
  [key: string]: any;
}

interface DomesticationData {
  domesticatedSequence?: string;
  [key: string]: any;
}


// Professional SVG icons for assembly methods
const MethodIcons = {
  goldenGate: (
    <svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M12 2L2 7l10 5 10-5-10-5z" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M2 17l10 5 10-5" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M2 12l10 5 10-5" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
    </svg>
  ),
  isothermal: (
    <svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M12 3C8 3 4 6 4 12s4 9 8 9" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M12 3c4 0 8 3 8 9s-4 9-8 9" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeDasharray="3 3"/>
      <circle cx="7" cy="8" r="1.5" fill="currentColor"/>
      <circle cx="17" cy="8" r="1.5" fill="currentColor"/>
      <circle cx="7" cy="16" r="1.5" fill="currentColor"/>
      <circle cx="17" cy="16" r="1.5" fill="currentColor"/>
      <path d="M9 12h6" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
    </svg>
  ),
};

// Assembly method definitions - consolidated to 2 main methods
const ASSEMBLY_METHODS = {
  GOLDEN_GATE: {
    id: 'golden_gate',
    name: 'Golden Gate',
    description: 'Type IIS restriction enzyme-based assembly',
    icon: MethodIcons.goldenGate,
    color: '#f59e0b',
  },
  ISOTHERMAL: {
    id: 'isothermal',
    name: 'Isothermal Assembly',
    description: 'Overlap-based assembly (Gibson / NEBuilder HiFi)',
    icon: MethodIcons.isothermal,
    color: '#7c3aed',
  },
};

// Protocol variants for isothermal assembly
const ISOTHERMAL_VARIANTS = {
  GIBSON: {
    id: 'gibson',
    name: 'Gibson Assembly',
    description: 'Standard isothermal assembly',
    badge: 'DIY / NEB E5510',
    optimalOverlap: { min: 20, max: 40, ideal: 25 },
    optimalTm: { min: 50, max: 60, ideal: 55 },
  },
  NEBUILDER: {
    id: 'nebuilder_hifi',
    name: 'NEBuilder HiFi',
    description: 'High-fidelity assembly',
    badge: 'NEB E5520',
    optimalOverlap: { min: 15, max: 35, ideal: 20 },
    optimalTm: { min: 48, max: 65, ideal: 55 },
    recommended: true,
  },
};

const MAX_FRAGMENTS = 12; // Increased for Gibson/HiFi methods
const MAX_GOLDEN_GATE_FRAGMENTS = 5;

// Part type definitions with colors
const PART_TYPES = {
  promoter: { name: 'Promoter', color: '#ef4444', icon: 'P' },
  rbs: { name: 'RBS/5\'UTR', color: '#f97316', icon: 'R' },
  cds: { name: 'CDS', color: '#22c55e', icon: 'C' },
  terminator: { name: 'Terminator', color: '#3b82f6', icon: 'T' },
  other: { name: 'Other', color: '#8b5cf6', icon: 'O' },
  fragment: { name: 'Fragment', color: '#6366f1', icon: 'F' },
};

// Compact Part Card - collapsed by default, shows summary
function PartCard({ part, index, onChange, onRemove, canRemove, enzyme, onDragStart, onDragOver, onDrop, isDragging, onSplitSequence, useEnhancedWorkflow, onOpenEnhancedDomestication, onOpenWorkflowGuide }: PartCardProps) {
  const [isExpanded, setIsExpanded] = useState<boolean>(false);
  const [isEditing, setIsEditing] = useState<boolean>(!part.seq);
  const [showSplitPanel, setShowSplitPanel] = useState<boolean>(false);
  const [splitNumFragments, setSplitNumFragments] = useState<number>(3);
  const [isSplitting, setIsSplitting] = useState<boolean>(false);
  const [splitError, setSplitError] = useState<string | null>(null);

  // Check if sequence is large enough to split (>500bp)
  const canSplit = part.seq && part.seq.length >= 500;

  const handleSeqChange = (e: React.ChangeEvent<HTMLTextAreaElement>) => {
    const seq = e.target.value.toUpperCase().replace(/[^ATGCNRYSWKMBDHV]/gi, '');
    // Reset range when sequence changes
    onChange(index, { ...part, seq, rangeStart: 1, rangeEnd: seq.length || 0 });
  };

  const handleNameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    onChange(index, { ...part, id: e.target.value });
  };

  const handleTypeChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    onChange(index, { ...part, type: e.target.value as Part['type'] });
  };

  const handleRangeChange = (field: string, value: string | number) => {
    // Allow free typing - store raw value, validate on blur
    const numValue = typeof value === 'string' ? parseInt(value) || 1 : value;
    if (field === 'rangeStart') {
      onChange(index, { ...part, rangeStart: numValue });
    } else {
      onChange(index, { ...part, rangeEnd: numValue });
    }
  };

  const handleRangeBlur = (field: string) => {
    const seqLen = part.seq?.length || 0;
    if (field === 'rangeStart') {
      const numVal = parseInt(String(part.rangeStart)) || 1;
      const newStart = Math.max(1, Math.min(numVal, seqLen));
      onChange(index, { ...part, rangeStart: newStart });
    } else {
      // Allow end < start for circular plasmids (wrap-around selection)
      const numVal = parseInt(String(part.rangeEnd)) || seqLen;
      const newEnd = Math.max(1, Math.min(numVal, seqLen));
      onChange(index, { ...part, rangeEnd: newEnd });
    }
  };

  const setFullRange = () => {
    onChange(index, { ...part, rangeStart: 1, rangeEnd: part.seq?.length || 0 });
  };

  // Handle splitting sequence into optimized fragments
  const handleSplitSequence = async () => {
    if (!canSplit || !onSplitSequence) return;

    setIsSplitting(true);
    setSplitError(null);

    try {
      // Call the parent's split handler which runs the optimizer
      await onSplitSequence(index, part.seq, splitNumFragments, enzyme);
      setShowSplitPanel(false);
    } catch (err: unknown) {
      setSplitError(err instanceof Error ? err.message : String(err));
    } finally {
      setIsSplitting(false);
    }
  };

  const typeInfo = PART_TYPES[part.type as keyof typeof PART_TYPES] || PART_TYPES.other;

  // Check for internal enzyme sites (only for Golden Gate when enzyme is specified)
  const hasSites = enzyme && part.seq && (() => {
    const enzSite = GOLDEN_GATE_ENZYMES[enzyme as keyof typeof GOLDEN_GATE_ENZYMES]?.recognition;
    if (!enzSite) return false;
    const complement: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
    const enzSiteRc = enzSite.split('').reverse().map(c =>
      complement[c] || c
    ).join('');
    return part.seq.toUpperCase().includes(enzSite) ||
           part.seq.toUpperCase().includes(enzSiteRc);
  })();

  // Get selected range
  const rangeStart = typeof part.rangeStart === 'number' ? part.rangeStart : (part.rangeStart ? parseInt(String(part.rangeStart)) : 1);
  const rangeEnd = typeof part.rangeEnd === 'number' ? part.rangeEnd : (part.rangeEnd ? parseInt(String(part.rangeEnd)) : part.seq?.length || 0);
  // Support circular wrap-around: if end < start, selection spans origin
  const selectedSeq = part.seq ? (
    rangeEnd >= rangeStart
      ? part.seq.slice(rangeStart - 1, rangeEnd)
      : part.seq.slice(rangeStart - 1) + part.seq.slice(0, rangeEnd)
  ) : '';
  const isSubrange = part.seq && (rangeStart > 1 || rangeEnd < part.seq.length);
  const isCircularWrap = rangeEnd < rangeStart;

  const gcContent = selectedSeq ? (
    ((selectedSeq.toUpperCase().match(/[GC]/g) || []).length / selectedSeq.length) * 100
  ).toFixed(1) : 0;

  return (
    <div
      className={`part-card-compact ${isDragging ? 'dragging' : ''} ${isExpanded ? 'expanded' : ''} ${hasSites ? 'has-warning' : ''}`}
      draggable
      onDragStart={(e: React.DragEvent) => onDragStart(e, index)}
      onDragOver={(e: React.DragEvent) => onDragOver(e, index)}
      onDrop={(e: React.DragEvent) => onDrop(e, index)}
    >
      {/* Main Row */}
      <div className="part-card-main">
        {/* Drag Handle */}
        <div className="drag-handle" title="Drag to reorder">
          <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
            <path d="M11 18c0 1.1-.9 2-2 2s-2-.9-2-2 .9-2 2-2 2 .9 2 2zm-2-8c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm0-6c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm6 4c1.1 0 2-.9 2-2s-.9-2-2-2-2 .9-2 2 .9 2 2 2zm0 2c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm0 6c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2z"/>
          </svg>
        </div>

        {/* Type Color Tag with Warning Indicator */}
        <div
          className={`part-type-tag ${hasSites ? 'has-sites' : ''}`}
          style={{ backgroundColor: typeInfo.color }}
          title={hasSites ? `${typeInfo.name} - has internal ${enzyme} sites` : typeInfo.name}
        >
          {typeInfo.icon}
        </div>

        {/* Part Info */}
        <div className="part-main-info">
          <input
            type="text"
            className="part-name-compact"
            placeholder={`Part ${index + 1}`}
            value={part.id}
            onChange={handleNameChange}
          />
          {part.seq && (
            <span className="part-summary">
              <span className={`bp-count ${isSubrange ? 'has-range' : ''}`}>
                {isSubrange ? `${selectedSeq.length}/${part.seq.length} bp` : `${part.seq.length} bp`}
              </span>
              {isSubrange && <span className="range-badge">{rangeStart}-{rangeEnd}</span>}
              <span className="gc-count">{gcContent}% GC</span>
            </span>
          )}
        </div>

        {/* Type Selector */}
        <select
          className="part-type-compact"
          value={part.type || 'other'}
          onChange={handleTypeChange}
        >
          {Object.entries(PART_TYPES).map(([key, val]) => (
            <option key={key} value={key}>{val.name}</option>
          ))}
        </select>

        {/* Actions */}
        <div className="part-actions">
          {part.seq ? (
            <button
              className="action-icon-btn"
              onClick={() => setIsExpanded(!isExpanded)}
              title={isExpanded ? 'Collapse' : 'View/Edit Sequence'}
            >
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                {isExpanded ? (
                  <path d="M7.41 15.41L12 10.83l4.59 4.58L18 14l-6-6-6 6z"/>
                ) : (
                  <path d="M7.41 8.59L12 13.17l4.59-4.58L18 10l-6 6-6-6 1.41-1.41z"/>
                )}
              </svg>
            </button>
          ) : (
            <button
              className="action-icon-btn primary"
              onClick={() => setIsExpanded(true)}
              title="Add Sequence"
            >
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                <path d="M19 13h-6v6h-2v-6H5v-2h6V5h2v6h6v2z"/>
              </svg>
            </button>
          )}
          {canRemove && (
            <button
              className="action-icon-btn danger"
              onClick={() => onRemove(index)}
              title="Remove part"
            >
              <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                <path d="M19 6.41L17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12z"/>
              </svg>
            </button>
          )}
        </div>
      </div>

      {/* Status Bar - shows when there are internal sites OR when domestication is configured */}
      {part.seq && (hasSites || part._domesticationApproved) && (
        <div className={`part-status-bar domestication ${part._domesticationApproved ? 'configured' : ''}`}>
          <div className="status-content">
            {part._domesticationApproved ? (
              <>
                <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                  <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
                </svg>
                <span>
                  Domestication configured: {part._domesticationStrategy === 'mutagenic_junction'
                    ? `${part._domesticationJunctions?.length || 0} mutagenic junction(s)`
                    : `${part._domesticationMutations?.length || 0} silent mutation(s)`}
                </span>
              </>
            ) : (
              <>
                <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                  <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z"/>
                </svg>
                <span>Internal {enzyme} sites detected</span>
              </>
            )}
          </div>
          <div className="status-actions">
            {useEnhancedWorkflow && onOpenEnhancedDomestication && (
              <button
                className="domestication-workflow-btn"
                onClick={() => onOpenEnhancedDomestication({ index, ...part })}
                title={part._domesticationApproved
                  ? "Reconfigure domestication strategy and primers"
                  : "Review mutations, select strategy, and design primers for domestication"}
              >
                <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                  <path d="M19.14 12.94c.04-.31.06-.63.06-.94 0-.31-.02-.63-.06-.94l2.03-1.58c.18-.14.23-.41.12-.61l-1.92-3.32c-.12-.22-.37-.29-.59-.22l-2.39.96c-.5-.38-1.03-.7-1.62-.94l-.36-2.54c-.04-.24-.24-.41-.48-.41h-3.84c-.24 0-.43.17-.47.41l-.36 2.54c-.59.24-1.13.57-1.62.94l-2.39-.96c-.22-.08-.47 0-.59.22L2.74 8.87c-.12.21-.08.47.12.61l2.03 1.58c-.04.31-.06.63-.06.94s.02.63.06.94l-2.03 1.58c-.18.14-.23.41-.12.61l1.92 3.32c.12.22.37.29.59.22l2.39-.96c.5.38 1.03.7 1.62.94l.36 2.54c.05.24.24.41.48.41h3.84c.24 0 .44-.17.47-.41l.36-2.54c.59-.24 1.13-.56 1.62-.94l2.39.96c.22.08.47 0 .59-.22l1.92-3.32c.12-.22.07-.47-.12-.61l-2.01-1.58zM12 15.6c-1.98 0-3.6-1.62-3.6-3.6s1.62-3.6 3.6-3.6 3.6 1.62 3.6 3.6-1.62 3.6-3.6 3.6z"/>
                </svg>
                {part._domesticationApproved ? 'Reconfigure' : 'Configure'}
              </button>
            )}
          </div>
        </div>
      )}

      {/* Expanded Sequence Editor */}
      {isExpanded && (
        <div className="part-expanded-content">
          <textarea
            className="sequence-textarea"
            placeholder={`Paste ${typeInfo.name.toLowerCase()} sequence (without ${enzyme} sites)...`}
            value={part.seq}
            onChange={handleSeqChange}
            rows={4}
            autoFocus={!part.seq}
          />

          {/* Range Selection */}
          {part.seq && part.seq.length > 0 && (
            <div className="range-selection-panel">
              <div className="range-header">
                <span className="range-title">
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                    <path d="M9.5 3A6.5 6.5 0 0 1 16 9.5c0 1.61-.59 3.09-1.56 4.23l.27.27h.79l5 5-1.5 1.5-5-5v-.79l-.27-.27A6.516 6.516 0 0 1 9.5 16 6.5 6.5 0 0 1 3 9.5 6.5 6.5 0 0 1 9.5 3m0 2C7 5 5 7 5 9.5S7 14 9.5 14 14 12 14 9.5 12 5 9.5 5z"/>
                  </svg>
                  Select Range
                </span>
                <button className="btn-link" onClick={setFullRange} title="Use full sequence">
                  Full Sequence
                </button>
              </div>
              <div className="range-controls">
                <div className="range-input-group">
                  <label>Start</label>
                  <input
                    type="text"
                    value={rangeStart}
                    onChange={(e: React.ChangeEvent<HTMLInputElement>) => handleRangeChange('rangeStart', e.target.value)}
                    onBlur={() => handleRangeBlur('rangeStart')}
                    className="range-input"
                    placeholder="1"
                  />
                </div>
                <div className="range-separator">
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                    <path d="M8.59 16.59L13.17 12 8.59 7.41 10 6l6 6-6 6-1.41-1.41z"/>
                  </svg>
                </div>
                <div className="range-input-group">
                  <label>End</label>
                  <input
                    type="text"
                    value={rangeEnd}
                    onChange={(e: React.ChangeEvent<HTMLInputElement>) => handleRangeChange('rangeEnd', e.target.value)}
                    onBlur={() => handleRangeBlur('rangeEnd')}
                    className="range-input"
                    placeholder={part.seq.length.toString()}
                  />
                </div>
                <div className="range-result">
                  <span className="range-length">{selectedSeq.length} bp selected</span>
                  {rangeEnd < rangeStart && <span className="circular-badge" title="Circular wrap-around selection">circular</span>}
                </div>
              </div>

              {/* Minimized Sequence Viewers */}
              <div className="range-preview-panels">
                <div className="range-preview start">
                  <div className="preview-label">
                    <span className="label-text">5' Start Region</span>
                    <span className="label-pos">pos {rangeStart}</span>
                  </div>
                  <div className="mini-sequence-viewer">
                    <span className="seq-context">{rangeStart > 10 ? '...' : ''}</span>
                    <span className="seq-before">
                      {part.seq.slice(Math.max(0, rangeStart - 11), rangeStart - 1)}
                    </span>
                    <span className="seq-selected">
                      {part.seq.slice(rangeStart - 1, Math.min(rangeStart + 19, rangeEnd))}
                    </span>
                    <span className="seq-after">
                      {rangeEnd > rangeStart + 20 ? '...' : ''}
                    </span>
                  </div>
                </div>
                <div className="range-preview end">
                  <div className="preview-label">
                    <span className="label-text">3' End Region</span>
                    <span className="label-pos">pos {rangeEnd}</span>
                  </div>
                  <div className="mini-sequence-viewer">
                    <span className="seq-context">
                      {rangeEnd - rangeStart > 20 ? '...' : ''}
                    </span>
                    <span className="seq-selected">
                      {part.seq.slice(Math.max(rangeStart - 1, rangeEnd - 20), rangeEnd)}
                    </span>
                    <span className="seq-after">
                      {part.seq.slice(rangeEnd, Math.min(rangeEnd + 10, part.seq.length))}
                    </span>
                    <span className="seq-context">{rangeEnd < part.seq.length - 10 ? '...' : ''}</span>
                  </div>
                </div>
              </div>
            </div>
          )}

          <div className="sequence-actions">
            <button className="btn-small" onClick={() => setIsExpanded(false)}>
              Done
            </button>
          </div>
        </div>
      )}

      {/* Split into Fragments Panel - inline optimizer */}
      {showSplitPanel && canSplit && (
        <div className="split-sequence-panel p-4 mt-2 border-t-2 border-amber-500" style={{
          background: 'linear-gradient(135deg, #fffbeb 0%, #fef3c7 100%)',
        }}>
          <div className="flex items-center justify-between mb-3">
            <h4 className="m-0 text-sm text-amber-800 flex items-center gap-2">
              <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
                <path d="M14 4l2.29 2.29-2.88 2.88 1.42 1.42 2.88-2.88L20 10V4h-6zm-4 0H4v6l2.29-2.29 4.71 4.7V20h2v-8.41l-5.29-5.3L10 4z"/>
              </svg>
              Split into Optimized Fragments
            </h4>
            <button
              onClick={() => setShowSplitPanel(false)}
              className="bg-none border-none cursor-pointer p-1"
            >
              <svg viewBox="0 0 24 24" width="18" height="18" fill="#92400e">
                <path d="M19 6.41L17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12z"/>
              </svg>
            </button>
          </div>

          <p className="text-xs text-amber-900 m-0 mb-4">
            This will find optimal junction positions to split your {part.seq.length.toLocaleString()} bp sequence
            into fragments with high-fidelity overhangs for {enzyme} Golden Gate assembly.
          </p>

          <div className="flex items-center gap-4 mb-4">
            <label className="text-sm text-amber-800 font-medium">
              Number of fragments:
            </label>
            <div className="flex gap-1">
              {[2, 3, 4, 5, 6].map(n => (
                <button
                  key={n}
                  onClick={() => setSplitNumFragments(n)}
                  className={`py-1.5 px-3 rounded-md cursor-pointer text-sm ${
                    splitNumFragments === n
                      ? 'border-2 border-amber-500 bg-amber-500 text-white font-semibold'
                      : 'border border-amber-600 bg-white text-amber-800 font-normal'
                  }`}
                >
                  {n}
                </button>
              ))}
            </div>
            <span className="text-xs text-amber-700">
              ({splitNumFragments - 1} junctions, ~{Math.round(part.seq.length / splitNumFragments).toLocaleString()} bp each)
            </span>
          </div>

          {splitError && (
            <div className="py-2 px-3 bg-red-50 border border-red-300 rounded-md text-red-600 text-xs mb-3">
              {splitError}
            </div>
          )}

          <div className="flex justify-end gap-2">
            <button
              onClick={() => setShowSplitPanel(false)}
              className="py-2 px-4 border border-amber-600 rounded-md bg-white text-amber-800 cursor-pointer text-sm"
            >
              Cancel
            </button>
            <button
              onClick={handleSplitSequence}
              disabled={isSplitting}
              className={`py-2 px-4 border-none rounded-md text-white text-sm font-semibold flex items-center gap-2 ${
                isSplitting ? 'bg-amber-600 cursor-not-allowed' : 'bg-amber-500 cursor-pointer'
              }`}
            >
              {isSplitting ? (
                <>
                  <span className="animate-spin">⟳</span>
                  Optimizing...
                </>
              ) : (
                <>
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                    <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
                  </svg>
                  Split & Optimize
                </>
              )}
            </button>
          </div>
        </div>
      )}
    </div>
  );
}

// Enhanced Circular Plasmid Visualization
function CircularPlasmidView({ parts, overhangs, totalLength }: CircularPlasmidViewProps) {
  const [hoveredJunction, setHoveredJunction] = useState<number | null>(null);
  const [hoveredPart, setHoveredPart] = useState<number | null>(null);

  // Calculate angles for each part
  const partAngles = useMemo((): PartWithAngles[] => {
    if (!parts.length) return [];

    const total = parts.reduce((sum: number, p: Part) => sum + (p.seq?.length || 100), 0);
    let currentAngle = -90; // Start at top

    return parts.map((part: Part, i: number): PartWithAngles => {
      const partLength = part.seq?.length || 100;
      const sweepAngle = (partLength / total) * 360;
      const startAngle = currentAngle;
      currentAngle += sweepAngle;
      return {
        ...part,
        startAngle,
        sweepAngle,
        midAngle: startAngle + sweepAngle / 2,
        index: i,
      };
    });
  }, [parts]);

  // Generate arc path
  const describeArc = (cx: number, cy: number, r: number, startAngle: number, endAngle: number) => {
    const start = {
      x: cx + r * Math.cos(Math.PI * startAngle / 180),
      y: cy + r * Math.sin(Math.PI * startAngle / 180)
    };
    const end = {
      x: cx + r * Math.cos(Math.PI * endAngle / 180),
      y: cy + r * Math.sin(Math.PI * endAngle / 180)
    };
    const largeArc = endAngle - startAngle > 180 ? 1 : 0;
    return `M ${start.x} ${start.y} A ${r} ${r} 0 ${largeArc} 1 ${end.x} ${end.y}`;
  };

  const cx = 150, cy = 150, r = 100;

  return (
    <div className="circular-plasmid-view">
      <svg viewBox="0 0 300 300" className="plasmid-svg">
        {/* Background circle */}
        <circle cx={cx} cy={cy} r={r} fill="none" stroke="#e5e7eb" strokeWidth="20" />

        {/* Part arcs */}
        {partAngles.map((part: PartWithAngles, i: number) => {
          const typeInfo = PART_TYPES[part.type as keyof typeof PART_TYPES] || PART_TYPES.other;
          const isHovered = hoveredPart === i;
          return (
            <g key={i}>
              <path
                d={describeArc(cx, cy, r, part.startAngle, part.startAngle + part.sweepAngle - 1)}
                fill="none"
                stroke={typeInfo.color}
                strokeWidth={isHovered ? 24 : 20}
                strokeLinecap="round"
                className="cursor-pointer transition-all duration-200"
                onMouseEnter={() => setHoveredPart(i)}
                onMouseLeave={() => setHoveredPart(null)}
              />
            </g>
          );
        })}

        {/* Junction markers */}
        {partAngles.map((part: PartWithAngles, i: number) => {
          const angle = (part.startAngle * Math.PI) / 180;
          const jx = cx + (r + 20) * Math.cos(angle);
          const jy = cy + (r + 20) * Math.sin(angle);
          const isHovered = hoveredJunction === i;

          return (
            <g
              key={`junction-${i}`}
              onMouseEnter={() => setHoveredJunction(i)}
              onMouseLeave={() => setHoveredJunction(null)}
              className="cursor-pointer"
            >
              <circle
                cx={cx + r * Math.cos(angle)}
                cy={cy + r * Math.sin(angle)}
                r={isHovered ? 6 : 4}
                fill="#10b981"
                className="transition-all duration-200"
              />
              {isHovered && (
                <g>
                  <rect
                    x={jx - 20}
                    y={jy - 10}
                    width="40"
                    height="20"
                    rx="4"
                    fill="rgba(0,0,0,0.8)"
                  />
                  <text
                    x={jx}
                    y={jy + 4}
                    textAnchor="middle"
                    fill="white"
                    fontSize="10"
                    fontFamily="Monaco, monospace"
                  >
                    {overhangs[i]}
                  </text>
                </g>
              )}
            </g>
          );
        })}

        {/* Part labels */}
        {partAngles.map((part: PartWithAngles, i: number) => {
          const labelR = r - 35;
          const angle = (part.midAngle * Math.PI) / 180;
          const lx = cx + labelR * Math.cos(angle);
          const ly = cy + labelR * Math.sin(angle);
          const typeInfo = PART_TYPES[part.type as keyof typeof PART_TYPES] || PART_TYPES.other;

          return (
            <g key={`label-${i}`}>
              <text
                x={lx}
                y={ly}
                textAnchor="middle"
                dominantBaseline="middle"
                fontSize="10"
                fontWeight="600"
                fill={hoveredPart === i ? typeInfo.color : '#374151'}
              >
                {(part.id || `Part ${i + 1}`).slice(0, 8)}
              </text>
            </g>
          );
        })}

        {/* Center info */}
        <text x={cx} y={cy - 10} textAnchor="middle" fontSize="24" fontWeight="700" fill="#1f2937">
          {totalLength || '—'}
        </text>
        <text x={cx} y={cy + 12} textAnchor="middle" fontSize="11" fill="#6b7280">
          bp total
        </text>
      </svg>

      {/* Hovered part tooltip */}
      {hoveredPart !== null && (
        <div className="plasmid-tooltip">
          <div className="tooltip-header" style={{ borderColor: PART_TYPES[partAngles[hoveredPart]?.type as keyof typeof PART_TYPES]?.color }}>
            <span className="tooltip-type">{PART_TYPES[partAngles[hoveredPart]?.type as keyof typeof PART_TYPES]?.name}</span>
            <span className="tooltip-name">{partAngles[hoveredPart]?.id || `Part ${hoveredPart + 1}`}</span>
          </div>
          <div className="tooltip-stats">
            <span>{partAngles[hoveredPart]?.seq?.length || 0} bp</span>
            <span>{overhangs[hoveredPart]} → {overhangs[hoveredPart + 1]}</span>
          </div>
        </div>
      )}

      {/* Legend */}
      <div className="plasmid-legend">
        {parts.map((part: Part, i: number) => {
          const typeInfo = PART_TYPES[part.type as keyof typeof PART_TYPES] || PART_TYPES.other;
          return (
            <div
              key={i}
              className={`legend-item ${hoveredPart === i ? 'active' : ''}`}
              onMouseEnter={() => setHoveredPart(i)}
              onMouseLeave={() => setHoveredPart(null)}
            >
              <span className="legend-dot" style={{ backgroundColor: typeInfo.color }}></span>
              <span className="legend-name">{part.id || `Part ${i + 1}`}</span>
            </div>
          );
        })}
      </div>
    </div>
  );
}

// Linear Assembly Diagram (enhanced)
function LinearAssemblyDiagram({ parts, overhangs }: LinearAssemblyDiagramProps) {
  const [hoveredOH, setHoveredOH] = useState<string | null>(null);

  return (
    <div className="linear-assembly-diagram">
      <div className="assembly-track">
        {parts.map((part: Part, i: number) => {
          const typeInfo = PART_TYPES[part.type] || PART_TYPES.other;
          const leftOH = overhangs[i];
          const rightOH = overhangs[i + 1];

          return (
            <React.Fragment key={i}>
              {/* Left junction */}
              <div
                className={`junction-node ${hoveredOH === leftOH ? 'highlighted' : ''}`}
                onMouseEnter={() => setHoveredOH(leftOH)}
                onMouseLeave={() => setHoveredOH(null)}
              >
                <span className="junction-seq">{leftOH}</span>
              </div>

              {/* Part block */}
              <div className="part-block" style={{ borderColor: typeInfo.color }}>
                <div className="part-color-stripe" style={{ backgroundColor: typeInfo.color }}></div>
                <div className="part-block-content">
                  <span className="part-block-name">{part.id || `Part ${i + 1}`}</span>
                  <span className="part-block-info">{typeInfo.name} • {part.seq?.length || 0} bp</span>
                </div>
              </div>

              {/* Right junction (only for last part) */}
              {i === parts.length - 1 && (
                <div
                  className={`junction-node ${hoveredOH === rightOH ? 'highlighted' : ''}`}
                  onMouseEnter={() => setHoveredOH(rightOH)}
                  onMouseLeave={() => setHoveredOH(null)}
                >
                  <span className="junction-seq">{rightOH}</span>
                </div>
              )}
            </React.Fragment>
          );
        })}
      </div>
      <div className="circular-indicator-badge">
        <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
          <path d="M12 4V1L8 5l4 4V6c3.31 0 6 2.69 6 6 0 1.01-.25 1.97-.7 2.8l1.46 1.46A7.93 7.93 0 0020 12c0-4.42-3.58-8-8-8zm0 14c-3.31 0-6-2.69-6-6 0-1.01.25-1.97.7-2.8L5.24 7.74A7.93 7.93 0 004 12c0 4.42 3.58 8 8 8v3l4-4-4-4v3z"/>
        </svg>
        Circular Assembly
      </div>
    </div>
  );
}

// Fidelity Gauge (enhanced with enzyme comparison)
function FidelityGauge({ overhangs, enzyme, staticFidelity }: FidelityGaugeProps) {
  const [showDetails, setShowDetails] = useState<boolean>(false);

  const experimentalFidelity = useMemo((): ExperimentalFidelity | null => {
    if (!overhangs || overhangs.length === 0) return null;
    try {
      return calculateExperimentalFidelity(overhangs, enzyme) as ExperimentalFidelity;
    } catch (e: unknown) {
      return null;
    }
  }, [overhangs, enzyme]);

  const enzymeComparison = useMemo((): EnzymeComparison | null => {
    if (!overhangs || overhangs.length === 0) return null;
    try {
      return compareEnzymeFidelity(overhangs) as EnzymeComparison;
    } catch (e: unknown) {
      return null;
    }
  }, [overhangs]);

  const problematicPairs = useMemo((): any[] => {
    if (!overhangs || overhangs.length === 0) return [];
    try {
      return findProblematicPairs(overhangs, enzyme, 0.05) as any[];
    } catch (e: unknown) {
      return [];
    }
  }, [overhangs, enzyme]);

  const fidelity = experimentalFidelity?.assemblyFidelity || staticFidelity;
  const fidelityPercent = (fidelity * 100).toFixed(1);
  const isExperimental = experimentalFidelity?.source === 'experimental';

  const getFidelityColor = (f: number) => {
    if (f >= 0.95) return '#22c55e';
    if (f >= 0.90) return '#84cc16';
    if (f >= 0.80) return '#f59e0b';
    return '#ef4444';
  };

  const getFidelityStatus = (f: number) => {
    if (f >= 0.95) return 'Excellent';
    if (f >= 0.90) return 'Good';
    if (f >= 0.80) return 'Moderate';
    return 'Low';
  };

  return (
    <div className="fidelity-gauge-compact">
      <div className="gauge-main">
        <div className="gauge-circle">
          <svg viewBox="0 0 100 100">
            <circle cx="50" cy="50" r="40" fill="none" stroke="#e5e7eb" strokeWidth="8" />
            <circle
              cx="50"
              cy="50"
              r="40"
              fill="none"
              stroke={getFidelityColor(fidelity)}
              strokeWidth="8"
              strokeLinecap="round"
              strokeDasharray={`${fidelity * 251.2} 251.2`}
              transform="rotate(-90 50 50)"
            />
          </svg>
          <div className="gauge-center">
            <span className="gauge-percent">{fidelityPercent}</span>
            <span className="gauge-unit">%</span>
          </div>
        </div>
        <div className="gauge-label">
          <span className="status-text" style={{ color: getFidelityColor(fidelity) }}>
            {getFidelityStatus(fidelity)}
          </span>
          {isExperimental && <span className="exp-badge">Exp. Data</span>}
        </div>
      </div>

      {problematicPairs.length > 0 && (
        <div className="warning-banner">
          <span className="warning-icon">⚠️</span>
          <span>{problematicPairs.length} potential cross-ligation pair(s)</span>
        </div>
      )}

      <div className="fidelity-buttons">
        <button className="details-toggle" onClick={() => setShowDetails(!showDetails)}>
          {showDetails ? '▲ Hide Enzymes' : '▼ Compare Enzymes'}
        </button>
      </div>

      {showDetails && enzymeComparison && (
        <div className="enzyme-comparison-list">
          {enzymeComparison.ranked.map((e: any, i: number) => (
            <div key={e.enzyme} className={`enzyme-row ${e.enzyme === enzyme ? 'current' : ''}`}>
              <span className="enzyme-rank">#{i + 1}</span>
              <span className="enzyme-label">{e.enzyme}</span>
              <div className="enzyme-bar-track">
                <div
                  className="enzyme-bar-fill"
                  style={{ width: `${e.assemblyFidelity * 100}%`, backgroundColor: getFidelityColor(e.assemblyFidelity) }}
                />
              </div>
              <span className="enzyme-pct" style={{ color: getFidelityColor(e.assemblyFidelity) }}>
                {e.fidelityPercent}
              </span>
              {e.enzyme === enzyme && <span className="current-badge">✓</span>}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

// ============================================================================
// GG PRIMER CARD - Clean card-based primer display (adapted from Isothermal)
// ============================================================================
interface GGPrimerCardProps {
  primer: any;
  direction: 'Forward' | 'Reverse';
  partId: string;
  overhang?: string;
  onCopy: (text: string, label: string) => void;
  copiedId: string | null;
}

function GGPrimerCard({ primer, direction, partId, overhang, onCopy, copiedId }: GGPrimerCardProps) {
  const [showDetails, setShowDetails] = useState(false);
  const safePrimer = primer || {};
  const structure = safePrimer.structure || {};

  const isForward = direction === 'Forward';
  const dirColor = isForward ? '#3b82f6' : '#8b5cf6';
  const copyId = `${partId}-${isForward ? 'fwd' : 'rev'}`;

  // Calculate segment lengths for structure bar
  const extraLen = structure.extra?.length || 0;
  const enzymeLen = (structure.bsaISite || structure.recognitionSite)?.length || 0;
  const spacerLen = structure.spacer?.length || 0;
  const overhangLen = structure.overhang?.length || 0;
  // For Golden Gate, homology IS the annealing region (template-binding)
  const annealingLen = structure.annealing?.length || structure.homology?.length || 0;
  const totalLen = safePrimer.length || (extraLen + enzymeLen + spacerLen + overhangLen + annealingLen);

  // Get quality color
  const getQualityColor = (score: number) => {
    if (score >= 85) return '#22c55e';
    if (score >= 70) return '#3b82f6';
    if (score >= 55) return '#f59e0b';
    return '#ef4444';
  };

  // Check for issues - use homology or annealing for 3' end check
  const annealingSeq = structure.annealing || structure.homology || '';
  const last2 = annealingSeq.slice(-2).toUpperCase();
  const gcClamp = (last2.match(/[GC]/g) || []).length;
  const hasGGGG = /GGGG/i.test(annealingSeq);

  const breakdown = safePrimer.breakdown;

  return (
    <div className="gg-primer-card">
      {/* Header */}
      <div className="gg-primer-card-header">
        <span className="gg-primer-direction" style={{ backgroundColor: dirColor }}>
          {isForward ? '5′→3′' : '3′←5′'}
        </span>
        <span className="gg-primer-name">{partId}_{isForward ? 'F' : 'R'}</span>
        {safePrimer.qualityScore != null && (
          <span className="gg-primer-quality" style={{
            backgroundColor: getQualityColor(safePrimer.qualityScore) + '15',
            color: getQualityColor(safePrimer.qualityScore),
            borderColor: getQualityColor(safePrimer.qualityScore) + '40'
          }}>
            {safePrimer.qualityScore}/100
          </span>
        )}
        <button
          className={`gg-primer-copy ${copiedId === copyId ? 'copied' : ''}`}
          onClick={() => onCopy(safePrimer.sequence || '', copyId)}
          title="Copy sequence"
        >
          {copiedId === copyId ? (
            <svg viewBox="0 0 24 24" width="14" height="14" fill="#22c55e">
              <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
            </svg>
          ) : (
            <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
              <path d="M16 1H4c-1.1 0-2 .9-2 2v14h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z" />
            </svg>
          )}
        </button>
      </div>

      {/* Structure Visualization Bar */}
      <div className="gg-primer-structure">
        <div className="gg-structure-bar">
          {extraLen > 0 && (
            <div className="gg-segment extra" style={{ flex: extraLen }} title={`Extra: ${extraLen}bp`} />
          )}
          {enzymeLen > 0 && (
            <div className="gg-segment enzyme" style={{ flex: enzymeLen }} title={`BsaI site: ${enzymeLen}bp`} />
          )}
          {spacerLen > 0 && (
            <div className="gg-segment spacer" style={{ flex: spacerLen }} title={`Spacer: ${spacerLen}bp`} />
          )}
          {overhangLen > 0 && (
            <div className="gg-segment overhang" style={{ flex: overhangLen }} title={`Overhang: ${overhangLen}bp`} />
          )}
          {annealingLen > 0 && (
            <div className="gg-segment annealing" style={{ flex: annealingLen }} title={`Annealing: ${annealingLen}bp`} />
          )}
        </div>
        <div className="gg-structure-legend">
          {enzymeLen > 0 && <span className="legend-item enzyme">BsaI</span>}
          {overhangLen > 0 && <span className="legend-item overhang">{overhang || structure.overhang || 'OH'}</span>}
          {annealingLen > 0 && <span className="legend-item annealing">{annealingLen}bp anneal</span>}
        </div>
      </div>

      {/* Sequence Display - with underline for annealing region */}
      <div className="gg-primer-sequence">
        <code>
          {structure.extra && <span className="seq-extra">{structure.extra}</span>}
          {(structure.bsaISite || structure.recognitionSite) && (
            <span className="seq-enzyme">{structure.bsaISite || structure.recognitionSite}</span>
          )}
          {structure.spacer && <span className="seq-spacer">{structure.spacer}</span>}
          {structure.overhang && <span className="seq-overhang">{structure.overhang}</span>}
          {(structure.homology || structure.annealing) && (
            <span className="seq-annealing">{structure.homology || structure.annealing}</span>
          )}
        </code>
      </div>

      {/* Metrics Grid */}
      <div className="gg-primer-metrics">
        <div className="gg-metric">
          <span className="gg-metric-value">{totalLen}</span>
          <span className="gg-metric-label">bp</span>
        </div>
        <div className="gg-metric">
          <span className="gg-metric-value">{safePrimer.tm?.toFixed(1) || '—'}</span>
          <span className="gg-metric-label">Tm °C</span>
        </div>
        <div className="gg-metric">
          <span className="gg-metric-value">{safePrimer.gc?.toFixed(0) || '—'}</span>
          <span className="gg-metric-label">% GC</span>
        </div>
        <div className="gg-metric">
          <span className={`gg-metric-value ${gcClamp >= 1 ? 'good' : 'warning'}`}>{last2 || '—'}</span>
          <span className="gg-metric-label">3′ end</span>
        </div>
        {breakdown?.hairpin && (
          <div className="gg-metric">
            <span className={`gg-metric-value ${breakdown.hairpin.score >= 70 ? 'good' : breakdown.hairpin.score >= 50 ? 'warning' : 'bad'}`}>
              {breakdown.hairpin.value?.toFixed(1) || '—'}
            </span>
            <span className="gg-metric-label">Hairpin</span>
          </div>
        )}
        {hasGGGG && (
          <div className="gg-metric warning">
            <span className="gg-metric-value bad">G4</span>
            <span className="gg-metric-label">Risk</span>
          </div>
        )}
      </div>

      {/* Expandable Details */}
      {breakdown && (
        <button className="gg-details-toggle" onClick={() => setShowDetails(!showDetails)}>
          {showDetails ? '▼ Hide details' : '▶ Thermodynamics'}
        </button>
      )}

      {showDetails && breakdown && (
        <div className="gg-primer-details">
          <div className="gg-detail-row">
            <span>Hairpin ΔG</span>
            <span className={breakdown.hairpin?.score >= 70 ? 'good' : breakdown.hairpin?.score >= 50 ? 'warning' : 'bad'}>
              {breakdown.hairpin?.value?.toFixed(1) || '0.0'} kcal/mol
            </span>
          </div>
          <div className="gg-detail-row">
            <span>Self-dimer ΔG</span>
            <span className={breakdown.homodimer?.score >= 70 ? 'good' : breakdown.homodimer?.score >= 50 ? 'warning' : 'bad'}>
              {breakdown.homodimer?.value?.toFixed(1) || '0.0'} kcal/mol
            </span>
          </div>
          {breakdown.terminal3DG && (
            <div className="gg-detail-row">
              <span>3′ Duplex ΔG</span>
              <span className={breakdown.terminal3DG?.score >= 70 ? 'good' : 'warning'}>
                {breakdown.terminal3DG?.value?.toFixed(1) || '—'} kcal/mol
              </span>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

// Pair Quality Panel for Golden Gate
interface GGPairQualityProps {
  pairQuality: any;
  fwdTm?: number;
  revTm?: number;
}

function GGPairQualityPanel({ pairQuality, fwdTm, revTm }: GGPairQualityProps) {
  if (!pairQuality) return null;

  const tmDiff = pairQuality.tmDifference?.value || Math.abs((fwdTm || 0) - (revTm || 0));
  const heterodimer = pairQuality.heterodimer;

  const getQualityColor = (quality: string) => {
    if (quality === 'excellent') return '#22c55e';
    if (quality === 'good') return '#3b82f6';
    if (quality === 'acceptable') return '#f59e0b';
    return '#ef4444';
  };

  return (
    <div className="gg-pair-panel">
      <div className="gg-pair-header">
        <span>Primer Pair Quality</span>
        {pairQuality.score != null && (
          <span className="gg-pair-score" style={{
            backgroundColor: getQualityColor(pairQuality.quality) + '15',
            color: getQualityColor(pairQuality.quality)
          }}>
            {pairQuality.score}/100 • {pairQuality.quality}
          </span>
        )}
      </div>
      <div className="gg-pair-metrics">
        <div className="gg-pair-metric">
          <span className="gg-pair-label">ΔTm</span>
          <span className={`gg-pair-value ${tmDiff <= 2 ? 'good' : tmDiff <= 5 ? 'warning' : 'bad'}`}>
            {tmDiff.toFixed(1)}°C
          </span>
        </div>
        {heterodimer && (
          <div className="gg-pair-metric">
            <span className="gg-pair-label">Heterodimer</span>
            <span className={`gg-pair-value ${heterodimer.score >= 70 ? 'good' : heterodimer.score >= 50 ? 'warning' : 'bad'}`}>
              {heterodimer.fullPrimers?.toFixed(1) || '—'} kcal/mol
            </span>
          </div>
        )}
      </div>
      {heterodimer && heterodimer.fullPrimers < -8 && (
        <div className="gg-pair-warning">
          <svg viewBox="0 0 24 24" width="14" height="14" fill="#f59e0b">
            <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
          </svg>
          Strong heterodimer may cause primer-dimer artifacts
        </div>
      )}
    </div>
  );
}

// Enhanced Primer Results with better contrast and CSV export
function PrimerResults({ result, onCopy, method, enzyme, isGoldenGate: isGoldenGateProp }: PrimerResultsProps) {
  // Use prop if provided, otherwise derive from method
  const isGoldenGate = isGoldenGateProp !== undefined ? isGoldenGateProp : method === 'golden_gate';
  const [activeTab, setActiveTab] = useState<string>('overview');
  const [copiedId, setCopiedId] = useState<string | null>(null);
  const [expandedThermo, setExpandedThermo] = useState<Record<string, boolean>>({}); // Track expanded thermo panels per part

  // Generate cross-ligation heatmap data from result's overhangs
  const crossLigationData = useMemo(() => {
    if (!result?.parts || method !== 'golden_gate') return null;
    try {
      // Extract unique overhangs from parts (each part has left and right overhang)
      // Also include _overhang for domestication fragments and check result.overhangs directly
      const overhangs: string[] = [];

      // First, try to use the result.overhangs array directly if available
      // This is the authoritative source from the optimizer
      if (result.overhangs && Array.isArray(result.overhangs)) {
        result.overhangs.forEach((oh: string) => {
          if (oh && !overhangs.includes(oh)) {
            overhangs.push(oh);
          }
        });
      }

      // If no overhangs from result.overhangs, extract from parts
      if (overhangs.length === 0) {
        result.parts.forEach((part: any) => {
          if (part.leftOverhang && !overhangs.includes(part.leftOverhang)) {
            overhangs.push(part.leftOverhang);
          }
          if (part.rightOverhang && !overhangs.includes(part.rightOverhang)) {
            overhangs.push(part.rightOverhang);
          }
          // Include domestication junction overhangs (_overhang property)
          if (part._overhang && !overhangs.includes(part._overhang)) {
            overhangs.push(part._overhang);
          }
        });
      }

      if (overhangs.length < 2) return null;
      return generateCrossLigationHeatmap(overhangs, enzyme || 'BsaI') as CrossLigationData;
    } catch (e: unknown) {
      return null;
    }
  }, [result, enzyme, method]);

  // Toggle thermo panel for a specific part
  const toggleThermo = (partId: string) => {
    setExpandedThermo(prev => ({
      ...prev,
      [partId]: !prev[partId]
    }));
  };

  if (!result) return null;

  const copyToClipboard = (text: string, id: string) => {
    navigator.clipboard.writeText(text);
    setCopiedId(id);
    setTimeout(() => setCopiedId(null), 1500);
    onCopy && onCopy(`Copied ${id}`);
  };

  const exportAllPrimersFasta = () => {
    const lines: string[] = [];
    result.parts.forEach(part => {
      if (part.primers) {
        lines.push(`>${part.id}_FWD`);
        lines.push(part.primers.forward.sequence);
        lines.push(`>${part.id}_REV`);
        lines.push(part.primers.reverse.sequence);
      }
    });
    copyToClipboard(lines.join('\n'), 'all primers (FASTA)');
  };

  const downloadCSV = () => {
    const headers = ['Part', 'Direction', 'Sequence', 'Length', 'Tm (°C)', 'GC (%)', 'Left Overhang', 'Right Overhang'];
    const rows: string[] = [];

    result.parts.forEach(part => {
      if (part.primers) {
        rows.push([
          part.id,
          'Forward',
          part.primers.forward.sequence,
          part.primers.forward.length,
          part.primers.forward.tm.toFixed(1),
          part.primers.forward.gc.toFixed(1),
          part.leftOverhang,
          part.rightOverhang
        ].join(','));
        rows.push([
          part.id,
          'Reverse',
          part.primers.reverse.sequence,
          part.primers.reverse.length,
          part.primers.reverse.tm.toFixed(1),
          part.primers.reverse.gc.toFixed(1),
          part.leftOverhang,
          part.rightOverhang
        ].join(','));
      }
    });

    const csv = [headers.join(','), ...rows].join('\n');
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'golden_gate_primers.csv';
    a.click();
    URL.revokeObjectURL(url);
  };

  const downloadFasta = () => {
    const fasta = `>assembled_${result.parts.map((p: any) => p.id).join('+')}\n${(result.assembledSequence || '').match(/.{1,60}/g)?.join('\n') || ''}`;
    const blob = new Blob([fasta], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'assembled_sequence.fasta';
    a.click();
    URL.revokeObjectURL(url);
  };

  const copyIDTFormat = () => {
    const lines: string[] = [];
    result.parts.forEach(part => {
      if (part.primers) {
        // Forward primer: Name, Sequence, Notes (with Ta), Scale, Purification
        const fwdTa = part.primers.forward.tm.toFixed(1);
        const revTa = part.primers.reverse.tm.toFixed(1);
        lines.push([
          `${part.id}_FWD`,
          part.primers.forward.sequence,
          `Ta=${fwdTa}C`,
          '25nm',
          'STD'
        ].join('\t'));
        lines.push([
          `${part.id}_REV`,
          part.primers.reverse.sequence,
          `Ta=${revTa}C`,
          '25nm',
          'STD'
        ].join('\t'));
      }
    });
    copyToClipboard(lines.join('\n'), 'primers (IDT format)');
  };

  return (
    <div className="results-panel-enhanced">
      {/* Results Header */}
      <div className="results-header-card">
        <div className="success-badge">
          <svg viewBox="0 0 24 24" width="24" height="24" fill="currentColor">
            <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
          </svg>
        </div>
        <div className="results-title">
          <h3>Assembly Design Complete</h3>
          <p>
            {result.parts.length} {isGoldenGate ? 'parts' : 'fragments'} • {result.assembledLength?.toLocaleString()} bp
            {isGoldenGate ? (
              <> • {result.fidelity?.percentage} fidelity</>
            ) : (
              <> • {result.parts.length} junctions</>
            )}
          </p>
        </div>
        <div className="export-buttons">
          <button className="export-btn idt" onClick={copyIDTFormat} title="Copy all primers in IDT bulk order format">
            <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
              <path d="M16 1H4c-1.1 0-2 .9-2 2v14h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z"/>
            </svg>
            IDT
          </button>
          <button className="export-btn csv" onClick={downloadCSV} title="Download all primers as CSV">
            <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
              <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
            </svg>
            CSV
          </button>
          <button className="export-btn fasta" onClick={exportAllPrimersFasta} title="Copy all primers as FASTA">
            <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
              <path d="M16 1H4c-1.1 0-2 .9-2 2v14h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z"/>
            </svg>
            FASTA
          </button>
        </div>
      </div>

      {/* Warnings */}
      {result.warnings && result.warnings.length > 0 && (
        <div className="warnings-card">
          {result.warnings.map((w: string, i: number) => (
            <div key={i} className="warning-item">
              <svg className="warning-icon" viewBox="0 0 24 24" width="16" height="16" fill="#f59e0b">
                <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
              </svg>
              {w}
            </div>
          ))}
        </div>
      )}

      {/* Workflow Steps - shown when user-configured domestication with multi-step workflow */}
      {(() => { console.log('[PrimerResults] _autoDomestication:', result._autoDomestication); return null; })()}
      {result._autoDomestication?.details?.some(d => d.workflowSteps?.length > 0) && (
        <div className="workflow-steps-card">
          <div className="workflow-header">
            <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
              <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm-5 14H7v-2h7v2zm3-4H7v-2h10v2zm0-4H7V7h10v2z"/>
            </svg>
            <span className="workflow-title">Multi-Step Workflow Required</span>
            <span className="workflow-badge">Silent Mutation Strategy</span>
          </div>
          <div className="workflow-steps-list">
            {result._autoDomestication.details
              .filter(d => d.workflowSteps?.length > 0)
              .flatMap((d: any, di: number) => d.workflowSteps.map((step: any, si: number) => (
                <div key={`${di}-${si}`} className={`workflow-step-item ${step.type}`}>
                  <div className="step-indicator">
                    <span className="step-number">{step.step}</span>
                  </div>
                  <div className="step-content">
                    <div className="step-title">{step.title}</div>
                    <div className="step-description">{step.description}</div>
                    {step.primers?.length > 0 && (
                      <div className="step-primer-count">
                        {step.primers.length} primer pair{step.primers.length !== 1 ? 's' : ''}
                      </div>
                    )}
                  </div>
                </div>
              )))}
          </div>
          <div className="workflow-note">
            <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
              <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-6h2v6zm0-8h-2V7h2v2z"/>
            </svg>
            <span>
              Complete Step 1 (mutagenesis) before proceeding to Step 2 (assembly).
              See Primers tab for detailed sequences.
            </span>
          </div>
        </div>
      )}

      {/* Tabs - Streamlined with Overview first */}
      <div className="results-tabs-bar">
        <button className={`tab ${activeTab === 'overview' ? 'active' : ''}`} onClick={() => setActiveTab('overview')}>
          <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor" style={{ marginRight: 6 }}>
            <path d="M3 13h8V3H3v10zm0 8h8v-6H3v6zm10 0h8V11h-8v10zm0-18v6h8V3h-8z"/>
          </svg>
          Overview
        </button>
        <button className={`tab ${activeTab === 'viewer' ? 'active' : ''}`} onClick={() => setActiveTab('viewer')}>
          <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor" style={{ marginRight: 6 }}>
            <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-2 15l-5-5 1.41-1.41L10 14.17l7.59-7.59L19 8l-9 9z"/>
          </svg>
          Construct
        </button>
        <button className={`tab ${activeTab === 'primers' ? 'active' : ''}`} onClick={() => setActiveTab('primers')}>
          <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor" style={{ marginRight: 6 }}>
            <path d="M9.4 16.6L4.8 12l4.6-4.6L8 6l-6 6 6 6 1.4-1.4zm5.2 0l4.6-4.6-4.6-4.6L16 6l6 6-6 6-1.4-1.4z"/>
          </svg>
          Primers ({result.parts.length * 2})
        </button>
        <button className={`tab ${activeTab === 'fidelity' ? 'active' : ''}`} onClick={() => setActiveTab('fidelity')}>
          <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor" style={{ marginRight: 6 }}>
            <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zM9 17H7v-7h2v7zm4 0h-2V7h2v10zm4 0h-2v-4h2v4z"/>
          </svg>
          Junctions
        </button>
        <button className={`tab ${activeTab === 'protocol' ? 'active' : ''}`} onClick={() => setActiveTab('protocol')}>
          <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor" style={{ marginRight: 6 }}>
            <path d="M14 2H6c-1.1 0-1.99.9-1.99 2L4 20c0 1.1.89 2 1.99 2H18c1.1 0 2-.9 2-2V8l-6-6zm2 16H8v-2h8v2zm0-4H8v-2h8v2zm-3-5V3.5L18.5 9H13z"/>
          </svg>
          Protocol
        </button>
      </div>

      {/* Tab Content */}
      <div className="tab-panel">
        {/* Overview Tab - Dashboard with all key information */}
        {activeTab === 'overview' && (
          <div className="gg-overview-panel">
            {/* Top Row: Quality Gauge + Export + Summary */}
            <div className="gg-overview-header">
              {/* Quality Gauge */}
              <div className="gg-quality-gauge">
                {(() => {
                  // Parse fidelity percentage from string like "95.2%" to decimal 0.952
                  const fidelityValue = result.fidelity?.percentage
                    ? parseFloat(result.fidelity.percentage.replace('%', '')) / 100
                    : 0;
                  const fidelityColor = fidelityValue >= 0.95 ? '#22c55e' :
                                       fidelityValue >= 0.90 ? '#84cc16' :
                                       fidelityValue >= 0.80 ? '#f59e0b' : '#ef4444';
                  const fidelityStatus = fidelityValue >= 0.95 ? 'Excellent' :
                                        fidelityValue >= 0.90 ? 'Good' :
                                        fidelityValue >= 0.80 ? 'Moderate' : 'Low';
                  return (
                    <>
                      <svg viewBox="0 0 120 120" className="gg-gauge-svg">
                        {/* Background circle */}
                        <circle cx="60" cy="60" r="50" fill="none" stroke="#e5e7eb" strokeWidth="10" />
                        {/* Progress circle */}
                        <circle
                          cx="60" cy="60" r="50"
                          fill="none"
                          stroke={fidelityColor}
                          strokeWidth="10"
                          strokeLinecap="round"
                          strokeDasharray={`${fidelityValue * 314.16} 314.16`}
                          transform="rotate(-90 60 60)"
                          style={{ transition: 'stroke-dasharray 0.5s ease' }}
                        />
                      </svg>
                      <div className="gg-gauge-center">
                        <span className="gg-gauge-value">{(fidelityValue * 100).toFixed(1)}</span>
                        <span className="gg-gauge-unit">%</span>
                      </div>
                      <div className="gg-gauge-label">
                        <span className="gg-gauge-status" style={{ color: fidelityColor }}>
                          {fidelityStatus}
                        </span>
                        <span className="gg-gauge-sublabel">Assembly Fidelity</span>
                      </div>
                    </>
                  );
                })()}
              </div>

              {/* Export Panel */}
              <div className="gg-export-panel">
                <div className="gg-export-title">
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                    <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                  </svg>
                  Export Primers
                </div>
                <div className="gg-export-buttons">
                  <button className="gg-export-btn" onClick={copyIDTFormat} title="Copy all primers in IDT bulk order format">
                    <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                      <path d="M16 1H4c-1.1 0-2 .9-2 2v14h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z"/>
                    </svg>
                    IDT
                  </button>
                  <button className="gg-export-btn" onClick={downloadCSV} title="Download all primers as CSV">
                    <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                      <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                    </svg>
                    CSV
                  </button>
                  <button className="gg-export-btn" onClick={exportAllPrimersFasta} title="Copy all primers as FASTA">
                    <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                      <path d="M14 2H6c-1.1 0-1.99.9-1.99 2L4 20c0 1.1.89 2 1.99 2H18c1.1 0 2-.9 2-2V8l-6-6zm2 16H8v-2h8v2zm0-4H8v-2h8v2zm-3-5V3.5L18.5 9H13z"/>
                    </svg>
                    FASTA
                  </button>
                </div>
              </div>

              {/* Assembly Summary Cards */}
              <div className="gg-summary-cards">
                <div className="gg-summary-card">
                  <svg viewBox="0 0 24 24" width="20" height="20" fill="#6366f1">
                    <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm0 18c-4.41 0-8-3.59-8-8s3.59-8 8-8 8 3.59 8 8-3.59 8-8 8z"/>
                  </svg>
                  <span className="gg-summary-value">{result.assembledLength?.toLocaleString() || '—'}</span>
                  <span className="gg-summary-label">bp Total</span>
                </div>
                <div className="gg-summary-card">
                  <svg viewBox="0 0 24 24" width="20" height="20" fill="#8b5cf6">
                    <path d="M4 6H2v14c0 1.1.9 2 2 2h14v-2H4V6zm16-4H8c-1.1 0-2 .9-2 2v12c0 1.1.9 2 2 2h12c1.1 0 2-.9 2-2V4c0-1.1-.9-2-2-2zm0 14H8V4h12v12z"/>
                  </svg>
                  <span className="gg-summary-value">{result.parts?.length || 0}</span>
                  <span className="gg-summary-label">Parts</span>
                </div>
                <div className="gg-summary-card">
                  <svg viewBox="0 0 24 24" width="20" height="20" fill="#06b6d4">
                    <path d="M17 3H7c-1.1 0-2 .9-2 2v16l7-3 7 3V5c0-1.1-.9-2-2-2z"/>
                  </svg>
                  <span className="gg-summary-value">{result.fidelity?.individual?.length || 0}</span>
                  <span className="gg-summary-label">Junctions</span>
                </div>
                <div className="gg-summary-card">
                  <svg viewBox="0 0 24 24" width="20" height="20" fill="#10b981">
                    <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-2 15l-5-5 1.41-1.41L10 14.17l7.59-7.59L19 8l-9 9z"/>
                  </svg>
                  <span className="gg-summary-value">Circular</span>
                  <span className="gg-summary-label">Topology</span>
                </div>
              </div>
            </div>

            {/* Warnings Section */}
            {result.warnings && result.warnings.length > 0 && (
              <div className="gg-warnings-section">
                <div className="gg-section-header">
                  <svg viewBox="0 0 24 24" width="18" height="18" fill="#f59e0b">
                    <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
                  </svg>
                  <span>Notifications ({result.warnings.length})</span>
                </div>
                <div className="gg-warnings-list">
                  {result.warnings.slice(0, 5).map((w: string, i: number) => {
                    const isSuccess = w.startsWith('✅') || w.includes('auto-domesticated');
                    const isWarning = w.startsWith('⚠️') || w.includes('WARNING');
                    return (
                      <div key={i} className={`gg-warning-item ${isSuccess ? 'success' : isWarning ? 'warning' : 'info'}`}>
                        <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                          {isSuccess ? (
                            <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-2 15l-5-5 1.41-1.41L10 14.17l7.59-7.59L19 8l-9 9z"/>
                          ) : isWarning ? (
                            <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
                          ) : (
                            <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-6h2v6zm0-8h-2V7h2v2z"/>
                          )}
                        </svg>
                        <span>{w.replace(/^[✅⚠️🔴ℹ️]\s*/, '')}</span>
                      </div>
                    );
                  })}
                  {result.warnings.length > 5 && (
                    <button className="gg-show-more" onClick={() => setActiveTab('fidelity')}>
                      View all {result.warnings.length} notifications
                    </button>
                  )}
                </div>
              </div>
            )}

            {/* Quick Primer Reference */}
            <div className="gg-primer-quick">
              <div className="gg-section-header">
                <svg viewBox="0 0 24 24" width="18" height="18" fill="#6366f1">
                  <path d="M9.4 16.6L4.8 12l4.6-4.6L8 6l-6 6 6 6 1.4-1.4zm5.2 0l4.6-4.6-4.6-4.6L16 6l6 6-6 6-1.4-1.4z"/>
                </svg>
                <span>Quick Primer Reference</span>
                <button className="gg-view-details" onClick={() => setActiveTab('primers')}>
                  View Details
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                    <path d="M8.59 16.59L13.17 12 8.59 7.41 10 6l6 6-6 6-1.41-1.41z"/>
                  </svg>
                </button>
              </div>
              <div className="gg-primer-table">
                <div className="gg-primer-table-header">
                  <span className="col-part">Part</span>
                  <span className="col-fwd">Forward Primer</span>
                  <span className="col-rev">Reverse Primer</span>
                  <span className="col-oh">Overhangs</span>
                </div>
                {result.parts?.map((part: any, i: number) => (
                  <div key={i} className="gg-primer-table-row">
                    <span className="col-part">
                      <span className="part-color" style={{ backgroundColor: (PART_TYPES as Record<string, { color: string }>)[part.type]?.color || '#8b5cf6' }} />
                      {part.id || `Part ${i + 1}`}
                    </span>
                    <span className="col-fwd" title={part.primers?.forward?.sequence || ''}>
                      <code>{(part.primers?.forward?.sequence || '').slice(0, 25)}...</code>
                    </span>
                    <span className="col-rev" title={part.primers?.reverse?.sequence || ''}>
                      <code>{(part.primers?.reverse?.sequence || '').slice(0, 25)}...</code>
                    </span>
                    <span className="col-oh">
                      <code className="overhang-badge">{part.leftOverhang || '—'}</code>
                      <svg viewBox="0 0 24 24" width="12" height="12" fill="#9ca3af">
                        <path d="M8.59 16.59L13.17 12 8.59 7.41 10 6l6 6-6 6-1.41-1.41z"/>
                      </svg>
                      <code className="overhang-badge">{part.rightOverhang || '—'}</code>
                    </span>
                  </div>
                ))}
              </div>
            </div>

            {/* Junction Fidelity Summary */}
            <div className="gg-junction-summary">
              <div className="gg-section-header">
                <svg viewBox="0 0 24 24" width="18" height="18" fill="#06b6d4">
                  <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zM9 17H7v-7h2v7zm4 0h-2V7h2v10zm4 0h-2v-4h2v4z"/>
                </svg>
                <span>Junction Fidelity</span>
                <button className="gg-view-details" onClick={() => setActiveTab('fidelity')}>
                  View Matrix
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                    <path d="M8.59 16.59L13.17 12 8.59 7.41 10 6l6 6-6 6-1.41-1.41z"/>
                  </svg>
                </button>
              </div>
              <div className="gg-junction-bars">
                {result.fidelity?.individual?.map((item: any, i: number) => (
                  <div key={i} className="gg-junction-bar">
                    <span className="junction-label">J{i + 1}</span>
                    <code className="junction-overhang">{item.overhang}</code>
                    <div className="junction-bar-track">
                      <div
                        className="junction-bar-fill"
                        style={{
                          width: `${item.fidelity * 100}%`,
                          backgroundColor: item.fidelity >= 0.95 ? '#22c55e' :
                                          item.fidelity >= 0.90 ? '#84cc16' :
                                          item.fidelity >= 0.80 ? '#f59e0b' : '#ef4444'
                        }}
                      />
                    </div>
                    <span className="junction-pct" style={{
                      color: item.fidelity >= 0.95 ? '#22c55e' :
                             item.fidelity >= 0.90 ? '#84cc16' :
                             item.fidelity >= 0.80 ? '#f59e0b' : '#ef4444'
                    }}>
                      {(item.fidelity * 100).toFixed(0)}%
                    </span>
                  </div>
                ))}
              </div>
            </div>
          </div>
        )}

        {activeTab === 'viewer' && (
          <AssemblyViewer result={result} />
        )}

        {activeTab === 'primers' && (
          <div className="primers-panel">
            {/* Show mutagenesis primers first if silent mutation workflow */}
            {result._autoDomestication?.details?.some(d => d.workflowSteps?.length > 0) && (
              <div className="gg-mutagenesis-section">
                {/* Step 1 Header */}
                <div className="gg-step-header step1">
                  <div className="gg-step-badge">
                    <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                      <path d="M19.14 12.94c.04-.31.06-.63.06-.94 0-.31-.02-.63-.06-.94l2.03-1.58c.18-.14.23-.41.12-.61l-1.92-3.32c-.12-.22-.37-.29-.59-.22l-2.39.96c-.5-.38-1.03-.7-1.62-.94l-.36-2.54c-.04-.24-.24-.41-.48-.41h-3.84c-.24 0-.43.17-.47.41l-.36 2.54c-.59.24-1.13.57-1.62.94l-2.39-.96c-.22-.08-.47 0-.59.22L2.74 8.87c-.12.21-.08.47.12.61l2.03 1.58c-.04.31-.06.63-.06.94s.02.63.06.94l-2.03 1.58c-.18.14-.23.41-.12.61l1.92 3.32c.12.22.37.29.59.22l2.39-.96c.5.38 1.03.7 1.62.94l.36 2.54c.05.24.24.41.48.41h3.84c.24 0 .44-.17.47-.41l.36-2.54c.59-.24 1.13-.56 1.62-.94l2.39.96c.22.08.47 0 .59-.22l1.92-3.32c.12-.22.07-.47-.12-.61l-2.01-1.58zM12 15.6c-1.98 0-3.6-1.62-3.6-3.6s1.62-3.6 3.6-3.6 3.6 1.62 3.6 3.6-1.62 3.6-3.6 3.6z"/>
                    </svg>
                    <span>Step 1</span>
                  </div>
                  <h4>Site-Directed Mutagenesis</h4>
                  <span className="gg-step-subtitle">Remove internal restriction sites with silent mutations</span>
                </div>

                {/* Mutagenesis Primer Cards */}
                <div className="gg-mutagenesis-list">
                  {result._autoDomestication.details
                    .filter((d: any) => d.workflowSteps?.length > 0)
                    .flatMap((detail: any, di: number) =>
                      detail.workflowSteps
                        .filter((step: any) => step.type === 'pcr_mutagenesis' && step.primers?.length > 0)
                        .flatMap((step: any, si: number) => step.primers.map((primer: any, pi: number) => (
                          <div key={`mut-${di}-${si}-${pi}`} className="gg-mutagenesis-card">
                            {/* Card Header */}
                            <div className="gg-mut-card-header">
                              <span className="gg-mut-badge">
                                <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                                  <path d="M17.66 7.93L12 2.27 6.34 7.93c-3.12 3.12-3.12 8.19 0 11.31C7.9 20.8 9.95 21.58 12 21.58s4.1-.78 5.66-2.34c3.12-3.12 3.12-8.19 0-11.31zM12 19.59c-1.6 0-3.11-.62-4.24-1.76C6.62 16.69 6 15.19 6 13.59s.62-3.11 1.76-4.24L12 5.1v14.49z"/>
                                </svg>
                                M{di + 1}
                              </span>
                              <span className="gg-mut-name">{primer.name || `Mutagenesis ${pi + 1}`}</span>
                              {primer.codonChange && (
                                <span className="gg-mut-change">
                                  <code>{primer.codonChange}</code>
                                </span>
                              )}
                              {primer.pairQuality && (
                                <span className="gg-mut-quality" style={{
                                  backgroundColor: primer.pairQuality >= 85 ? '#22c55e15' : primer.pairQuality >= 70 ? '#3b82f615' : '#f59e0b15',
                                  color: primer.pairQuality >= 85 ? '#22c55e' : primer.pairQuality >= 70 ? '#3b82f6' : '#f59e0b',
                                }}>
                                  {primer.pairQuality}/100
                                </span>
                              )}
                            </div>

                            {/* Primer Sequences */}
                            <div className="gg-mut-primers">
                              <div className="gg-mut-primer-row">
                                <span className="gg-mut-dir fwd">FWD</span>
                                <div className="gg-mut-seq-wrapper">
                                  <code className="gg-mut-sequence" onClick={() => copyToClipboard(primer.forward?.sequence, 'Forward primer')}>
                                    {primer.forward?.sequence}
                                  </code>
                                </div>
                                <div className="gg-mut-stats">
                                  <span>{primer.forward?.length}bp</span>
                                  <span>Tm {primer.forward?.homologyTm || primer.overlapTm || '—'}°C</span>
                                </div>
                              </div>
                              <div className="gg-mut-primer-row">
                                <span className="gg-mut-dir rev">REV</span>
                                <div className="gg-mut-seq-wrapper">
                                  <code className="gg-mut-sequence" onClick={() => copyToClipboard(primer.reverse?.sequence, 'Reverse primer')}>
                                    {primer.reverse?.sequence}
                                  </code>
                                </div>
                                <div className="gg-mut-stats">
                                  <span>{primer.reverse?.length}bp</span>
                                  <span>Tm {primer.reverse?.homologyTm || primer.overlapTm || '—'}°C</span>
                                </div>
                              </div>
                            </div>

                            {/* Instructions */}
                            {primer.instructions && (
                              <div className="gg-mut-instructions">
                                <svg viewBox="0 0 24 24" width="14" height="14" fill="#f59e0b">
                                  <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-6h2v6zm0-8h-2V7h2v2z"/>
                                </svg>
                                <span>{primer.instructions}</span>
                              </div>
                            )}
                          </div>
                        )))
                    )}
                </div>

                {/* Step 2 Divider */}
                <div className="gg-step-divider">
                  <div className="gg-divider-line"></div>
                  <span className="gg-divider-text">After completing mutagenesis PCR, proceed to assembly</span>
                  <div className="gg-divider-line"></div>
                </div>

                {/* Step 2 Header */}
                <div className="gg-step-header step2">
                  <div className="gg-step-badge">
                    <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                      <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm-5 14H7v-2h7v2zm3-4H7v-2h10v2zm0-4H7V7h10v2z"/>
                    </svg>
                    <span>Step 2</span>
                  </div>
                  <h4>Golden Gate Assembly</h4>
                  <span className="gg-step-subtitle">Assemble domesticated fragments</span>
                </div>
              </div>
            )}

            {/* Part Primers - Card Layout */}
            <div className="gg-primers-list">
              {result.parts.filter((part: Part) => part.primers).map((part: Part, i: number) => (
                <div key={i} className="gg-part-primers">
                  {/* Part Header */}
                  <div className="gg-part-header">
                    <span className="gg-part-badge" style={{ backgroundColor: (PART_TYPES as Record<string, {color: string}>)[part.type]?.color || '#6b7280' }}>
                      {i + 1}
                    </span>
                    <span className="gg-part-name">{part.id}</span>
                    <span className="gg-part-length">{(part.length || part.seq?.length || 0).toLocaleString()} bp</span>
                    <div className="gg-part-overhangs">
                      <code>{part.leftOverhang}</code>
                      <svg viewBox="0 0 24 24" width="14" height="14" fill="#9ca3af">
                        <path d="M8.59 16.59L13.17 12 8.59 7.41 10 6l6 6-6 6-1.41-1.41z"/>
                      </svg>
                      <code>{part.rightOverhang}</code>
                    </div>
                  </div>

                  {/* Primer Cards Grid */}
                  <div className="gg-primers-grid">
                    <GGPrimerCard
                      primer={part.primers?.forward}
                      direction="Forward"
                      partId={part.id}
                      overhang={part.leftOverhang}
                      onCopy={copyToClipboard}
                      copiedId={copiedId}
                    />
                    <GGPrimerCard
                      primer={part.primers?.reverse}
                      direction="Reverse"
                      partId={part.id}
                      overhang={part.rightOverhang}
                      onCopy={copyToClipboard}
                      copiedId={copiedId}
                    />
                  </div>

                  {/* Pair Quality Panel */}
                  <GGPairQualityPanel
                    pairQuality={part.primers?.pairQuality}
                    fwdTm={part.primers?.forward?.tm}
                    revTm={part.primers?.reverse?.tm}
                  />
                </div>
              ))}
            </div>
          </div>
        )}

        {activeTab === 'fidelity' && (
          <div className="fidelity-panel">
            {/* Golden Gate: Fidelity & Cross-ligation Analysis */}
            {isGoldenGate ? (
              <>
                <div className="fidelity-list">
                  {result.fidelity?.individual?.map((item: any, i: number) => {
                    const getColor = (f: number) => f >= 0.95 ? '#22c55e' : f >= 0.90 ? '#84cc16' : f >= 0.80 ? '#f59e0b' : '#ef4444';
                    // Find cross-ligation hotspots for this overhang
                    const overhangUpper = item.overhang?.toUpperCase();
                    const junctionHotspots = crossLigationData?.hotspots?.filter((h: any) =>
                      h.source?.toUpperCase() === overhangUpper || h.target?.toUpperCase() === overhangUpper
                    ) || [];
                    const hasHotspot = junctionHotspots.length > 0;

                    return (
                      <div key={i} className={`fidelity-row-item ${hasHotspot ? 'has-hotspot' : ''}`}>
                        <div className="junction-main-row">
                          <span className="junction-label">Junction {i + 1}</span>
                          <code className="junction-code">{item.overhang}</code>
                          <div className="fidelity-bar-track">
                            <div
                              className="fidelity-bar-fill"
                              style={{ width: `${item.fidelity * 100}%`, backgroundColor: getColor(item.fidelity) }}
                            />
                          </div>
                          <span className="fidelity-pct" style={{ color: getColor(item.fidelity) }}>
                            {(item.fidelity * 100).toFixed(0)}%
                          </span>
                          {hasHotspot && (
                            <span className="crosslig-warning" title="Cross-ligation risk">
                              <svg viewBox="0 0 24 24" width="14" height="14" fill="#f59e0b">
                                <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
                              </svg>
                            </span>
                          )}
                        </div>
                        {hasHotspot && (
                          <div className="crosslig-details">
                            <span className="crosslig-label">Cross-ligation risk with:</span>
                            {junctionHotspots.slice(0, 3).map((h: any, j: number) => {
                              const pairOH = h.source?.toUpperCase() === overhangUpper ? h.target : h.source;
                              return (
                                <span key={j} className={`crosslig-pair severity-${h.severity}`}>
                                  {pairOH} ({h.ratioPercent})
                                </span>
                              );
                            })}
                            {junctionHotspots.length > 3 && (
                              <span className="crosslig-more">+{junctionHotspots.length - 3} more</span>
                            )}
                          </div>
                        )}
                      </div>
                    );
                  })}
                </div>

                {/* Cross-ligation Heatmap Section */}
                {crossLigationData && !crossLigationData.error && (
                  <div className="crosslig-heatmap-section">
                    <div className="crosslig-heatmap-header">
                      <h5>Cross-Ligation Matrix</h5>
                      <span className="crosslig-heatmap-info">
                        Diagonal = correct, Off-diagonal = cross-ligation risk
                      </span>
                    </div>
                    <CrossLigationHeatmapCompact heatmapData={crossLigationData as any} enzyme={enzyme} />
                  </div>
                )}
              </>
            ) : (
              /* Gibson/HiFi: Overlap Junction Analysis */
              <div className="overlap-junction-analysis">
                <div className="junction-analysis-header">
                  <h4>
                    <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
                      <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm-7 14h-2v-4H8l4-4 4 4h-2v4z"/>
                    </svg>
                    Overlap Junction Analysis
                  </h4>
                  <span className="junction-count-badge">{result.fidelity?.individual?.length || result.parts?.length || 0} junctions</span>
                </div>

                <div className="overlap-junction-list">
                  {(result.overlapAnalysis || result.fidelity?.individual || []).map((junction: any, i: number) => {
                    const overlap = junction;
                    const score = overlap.score || overlap.fidelity || 0.8;
                    const tm = overlap.tm || 55;
                    const gc = overlap.gc || 50;
                    const length = overlap.length || overlap.overhang?.length || 20;
                    const sequence = overlap.sequence || overlap.overhang || '';
                    const warnings = overlap.warnings || [];

                    const getScoreColor = (s: number) => s >= 0.8 ? '#22c55e' : s >= 0.6 ? '#84cc16' : s >= 0.4 ? '#f59e0b' : '#ef4444';
                    const getScoreTier = (s: number) => s >= 0.8 ? 'Excellent' : s >= 0.6 ? 'Good' : s >= 0.4 ? 'Acceptable' : 'Poor';
                    const getTmStatus = (t: number) => Math.abs(t - 55) <= 3 ? 'optimal' : Math.abs(t - 55) <= 7 ? 'acceptable' : 'suboptimal';
                    const getGcStatus = (g: number) => g >= 40 && g <= 60 ? 'optimal' : g >= 30 && g <= 70 ? 'acceptable' : 'suboptimal';

                    return (
                      <div key={i} className="overlap-junction-card">
                        <div className="junction-card-header">
                          <div className="junction-title">
                            <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                              <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-2 15l-5-5 1.41-1.41L10 14.17l7.59-7.59L19 8l-9 9z"/>
                            </svg>
                            <span>Junction {i + 1}</span>
                            <span className="fragment-connection">
                              Fragment {i + 1} → Fragment {i + 2 > result.parts?.length ? 1 : i + 2}
                            </span>
                          </div>
                          <div className="junction-score-badge" style={{ backgroundColor: getScoreColor(score) + '20', color: getScoreColor(score) }}>
                            {Math.round(score * 100)}% {getScoreTier(score)}
                          </div>
                        </div>

                        <div className="junction-card-body">
                          {/* Overlap Sequence */}
                          {sequence && (
                            <div className="overlap-sequence-row">
                              <span className="seq-label">Overlap</span>
                              <code className="overlap-seq" title="Click to copy" onClick={() => navigator.clipboard.writeText(sequence)}>
                                {sequence}
                              </code>
                              <span className="seq-length">{length} bp</span>
                            </div>
                          )}

                          {/* Metrics Grid */}
                          <div className="junction-metrics-grid">
                            {/* Tm Metric */}
                            <div className="junction-metric">
                              <div className="metric-header">
                                <svg viewBox="0 0 24 24" width="12" height="12" fill="currentColor">
                                  <path d="M15 13V5c0-1.66-1.34-3-3-3S9 3.34 9 5v8c-1.21.91-2 2.37-2 4 0 2.76 2.24 5 5 5s5-2.24 5-5c0-1.63-.79-3.09-2-4zm-4-8c0-.55.45-1 1-1s1 .45 1 1h-1v1h1v2h-1v1h1v2h-2V5z"/>
                                </svg>
                                <span>Melting Temp</span>
                              </div>
                              <div className="metric-value-row">
                                <span className={`metric-value tm-${getTmStatus(tm)}`}>{tm?.toFixed(1) || '—'}°C</span>
                                <span className="metric-target">target: 55°C</span>
                              </div>
                              <div className="metric-bar-container">
                                <div className="metric-bar-track">
                                  <div className="metric-bar-fill" style={{
                                    width: `${Math.min(100, Math.max(0, (tm / 70) * 100))}%`,
                                    backgroundColor: getTmStatus(tm) === 'optimal' ? '#22c55e' : getTmStatus(tm) === 'acceptable' ? '#f59e0b' : '#ef4444'
                                  }} />
                                  <div className="metric-bar-target" style={{ left: `${(55 / 70) * 100}%` }} />
                                </div>
                              </div>
                            </div>

                            {/* GC Content Metric */}
                            <div className="junction-metric">
                              <div className="metric-header">
                                <svg viewBox="0 0 24 24" width="12" height="12" fill="currentColor">
                                  <path d="M19 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h14c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zM9 17H7v-7h2v7zm4 0h-2V7h2v10zm4 0h-2v-4h2v4z"/>
                                </svg>
                                <span>GC Content</span>
                              </div>
                              <div className="metric-value-row">
                                <span className={`metric-value gc-${getGcStatus(gc)}`}>{gc?.toFixed(0) || '—'}%</span>
                                <span className="metric-target">optimal: 40-60%</span>
                              </div>
                              <div className="metric-bar-container">
                                <div className="metric-bar-track gc-track">
                                  <div className="metric-bar-fill" style={{
                                    width: `${gc || 50}%`,
                                    backgroundColor: getGcStatus(gc) === 'optimal' ? '#22c55e' : getGcStatus(gc) === 'acceptable' ? '#f59e0b' : '#ef4444'
                                  }} />
                                  <div className="gc-optimal-zone" />
                                </div>
                              </div>
                            </div>
                          </div>

                          {/* Warnings */}
                          {warnings.length > 0 && (
                            <div className="junction-warnings">
                              {warnings.map((w: any, wi: number) => (
                                <div key={wi} className="junction-warning-item">
                                  <svg viewBox="0 0 24 24" width="12" height="12" fill="#f59e0b">
                                    <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
                                  </svg>
                                  <span>{w}</span>
                                </div>
                              ))}
                            </div>
                          )}
                        </div>
                      </div>
                    );
                  })}
                </div>

                {/* Summary Stats */}
                <div className="overlap-summary-stats">
                  <div className="summary-stat">
                    <span className="stat-label">Average Tm</span>
                    <span className="stat-value">
                      {((result.overlapAnalysis || result.fidelity?.individual || []).reduce((sum: any, j: any) => sum + (j.tm || 55), 0) /
                        Math.max(1, (result.overlapAnalysis || result.fidelity?.individual || []).length)).toFixed(1)}°C
                    </span>
                  </div>
                  <div className="summary-stat">
                    <span className="stat-label">Average GC</span>
                    <span className="stat-value">
                      {((result.overlapAnalysis || result.fidelity?.individual || []).reduce((sum: any, j: any) => sum + (j.gc || 50), 0) /
                        Math.max(1, (result.overlapAnalysis || result.fidelity?.individual || []).length)).toFixed(0)}%
                    </span>
                  </div>
                  <div className="summary-stat">
                    <span className="stat-label">Overall Quality</span>
                    <span className="stat-value">
                      {Math.round(((result.overlapAnalysis || result.fidelity?.individual || []).reduce((sum: any, j: any) => sum + (j.score || j.fidelity || 0.8), 0) /
                        Math.max(1, (result.overlapAnalysis || result.fidelity?.individual || []).length)) * 100)}%
                    </span>
                  </div>
                </div>
              </div>
            )}
          </div>
        )}

        {activeTab === 'protocol' && result.protocol && (
          <div className="protocol-panel">
            <h4>{(result.protocol as any).title}</h4>
            <ol className="protocol-steps">
              {(result.protocol as any).steps.map((step: any, i: number) => (
                <li key={i}>
                  <strong>{step.title}</strong>
                  <ul>
                    {step.details.map((d: any, j: number) => <li key={j}>{d}</li>)}
                  </ul>
                </li>
              ))}
            </ol>
            <div className="protocol-notes">
              <h5>Notes</h5>
              <ul>
                {(result.protocol as any).notes.map((n: any, i: number) => <li key={i}>{n}</li>)}
              </ul>
            </div>
          </div>
        )}
      </div>

      {/* Assembled Sequence */}
      <div className="sequence-export-card">
        <div className="export-card-header">
          <h4>Assembled Sequence</h4>
          <div className="export-card-actions">
            <button onClick={() => copyToClipboard(result.assembledSequence || '', 'sequence')}>Copy</button>
            <button onClick={downloadFasta}>FASTA</button>
          </div>
        </div>
        <div className="sequence-preview-text">
          <code>{result.assembledSequence?.slice(0, 100)}...</code>
          <span className="seq-length-badge">{result.assembledSequence?.length} bp</span>
        </div>
      </div>
    </div>
  );
}

// Assembly Viewer with SeqViz
function AssemblyViewer({ result }: { result: any }) {
  const [viewMode, setViewMode] = useState<string>('circular');

  const annotations = useMemo(() => {
    if (!result?.parts) return [];
    const annots: any[] = [];
    let pos = 0;
    result.parts.forEach((part: any, i: number) => {
      const len = part.length || part.seq?.length || 0;
      if (len > 0) {
        annots.push({
          name: part.id || `Part ${i + 1}`,
          start: pos,
          end: pos + len,
          direction: 1,
          color: (PART_TYPES as any)[part.type]?.color || '#8b5cf6',
        });
        pos += len;
      }
    });
    return annots;
  }, [result]);

  if (!result?.assembledSequence) return null;

  return (
    <div className="assembly-viewer-panel">
      <div className="viewer-controls-bar">
        <button className={`view-btn ${viewMode === 'linear' ? 'active' : ''}`} onClick={() => setViewMode('linear')}>
          Linear
        </button>
        <button className={`view-btn ${viewMode === 'circular' ? 'active' : ''}`} onClick={() => setViewMode('circular')}>
          Circular
        </button>
      </div>
      <div className="seqviz-container">
        <SeqViz
          name="Assembled Construct"
          seq={result.assembledSequence.toUpperCase()}
          annotations={annotations}
          viewer={viewMode as 'linear' | 'circular' | 'both' | 'both_flip'}
          showComplement={false}
          showIndex={true}
          zoom={{ linear: 25 }}
          rotateOnScroll={false}
          style={{ height: '100%', width: '100%' }}
        />
      </div>
      <div className="parts-legend-bar">
        {result.parts.map((part: Part, i: number) => (
          <span key={i} className="legend-chip">
            <span className="chip-dot" style={{ backgroundColor: PART_TYPES[part.type]?.color || '#8b5cf6' }}></span>
            {part.id || `Part ${i + 1}`}
          </span>
        ))}
      </div>
    </div>
  );
}

// Overlap Settings Component for Gibson/HiFi methods
function OverlapSettings({ settings, onChange, variant }: { settings: any; onChange: any; variant: any }) {
  const isHiFi = variant === 'nebuilder_hifi';
  const defaults = isHiFi ? NEBUILDER_PARAMS : { MIN_OVERLAP: 15, MAX_OVERLAP: 40, TARGET_TM: 55 };

  // Extract min/max overlap values from NEBUILDER_PARAMS or fallback defaults
  const minOverlap = (defaults as any).overlap?.minLength || (defaults as any).MIN_OVERLAP || 15;
  const maxOverlap = (defaults as any).overlap?.maxLength || (defaults as any).MAX_OVERLAP || 40;

  return (
    <div className="overlap-settings">
      <div className="setting-row">
        <label>Target Overlap</label>
        <div className="setting-control">
          <input
            type="range"
            min={minOverlap}
            max={maxOverlap}
            value={settings.overlapLength}
            onChange={(e: React.ChangeEvent<HTMLInputElement>) => onChange({ ...settings, overlapLength: parseInt(e.target.value) })}
          />
          <span className="setting-value">{settings.overlapLength} bp</span>
        </div>
      </div>
      <div className="setting-row">
        <label>Target Tm</label>
        <div className="setting-control">
          <input
            type="range"
            min="48"
            max="65"
            value={settings.targetTm}
            onChange={(e: React.ChangeEvent<HTMLInputElement>) => onChange({ ...settings, targetTm: parseInt(e.target.value) })}
          />
          <span className="setting-value">{settings.targetTm}°C</span>
        </div>
      </div>
      {isHiFi && (
        <div className="setting-info">
          <span className="info-badge">NEB E5520</span>
          <span className="info-text">Recommended: 15-35bp overlaps, Tm 48-65°C</span>
        </div>
      )}
      {!isHiFi && (
        <div className="setting-info">
          <span className="info-badge">Gibson</span>
          <span className="info-text">Recommended: 20-40bp overlaps, Tm 50-60°C</span>
        </div>
      )}
    </div>
  );
}

// Overlap Quality Gauge for Gibson/HiFi Assembly Plan preview
function OverlapQualityGauge({ overlapLength, targetTm, fragmentCount, variant }: { overlapLength: any; targetTm: any; fragmentCount: any; variant: any }) {
  const isHiFi = variant === 'nebuilder_hifi';

  // Calculate estimated quality based on parameters
  const getOverlapQuality = () => {
    // Optimal ranges for Gibson/NEBuilder
    const optimalLength = isHiFi ? { min: 18, max: 30, ideal: 22 } : { min: 20, max: 35, ideal: 25 };
    const optimalTm = isHiFi ? { min: 50, max: 62, ideal: 55 } : { min: 52, max: 60, ideal: 55 };

    // Score overlap length (0-100)
    let lengthScore = 100;
    if (overlapLength < optimalLength.min) {
      lengthScore = Math.max(0, 100 - (optimalLength.min - overlapLength) * 10);
    } else if (overlapLength > optimalLength.max) {
      lengthScore = Math.max(0, 100 - (overlapLength - optimalLength.max) * 5);
    }

    // Score Tm (0-100)
    let tmScore = 100;
    if (targetTm < optimalTm.min) {
      tmScore = Math.max(0, 100 - (optimalTm.min - targetTm) * 8);
    } else if (targetTm > optimalTm.max) {
      tmScore = Math.max(0, 100 - (targetTm - optimalTm.max) * 8);
    }

    // Fragment count penalty (more fragments = harder assembly)
    const fragmentPenalty = Math.max(0, (fragmentCount - 3) * 3);

    // Combined score
    return Math.max(0, Math.min(100, (lengthScore * 0.4 + tmScore * 0.5 + (100 - fragmentPenalty) * 0.1)));
  };

  const quality = getOverlapQuality();
  const qualityPercent = quality.toFixed(0);

  const getQualityColor = (q: number) => {
    if (q >= 85) return '#22c55e';
    if (q >= 70) return '#84cc16';
    if (q >= 55) return '#f59e0b';
    return '#ef4444';
  };

  const getQualityStatus = (q: number) => {
    if (q >= 85) return 'Excellent';
    if (q >= 70) return 'Good';
    if (q >= 55) return 'Acceptable';
    return 'Suboptimal';
  };

  const getTmStatus = () => {
    const optimal = isHiFi ? 55 : 55;
    const diff = Math.abs(targetTm - optimal);
    if (diff <= 3) return { status: 'optimal', color: '#22c55e' };
    if (diff <= 7) return { status: 'acceptable', color: '#f59e0b' };
    return { status: 'suboptimal', color: '#ef4444' };
  };

  const tmStatus = getTmStatus();

  return (
    <div className="overlap-quality-gauge">
      <div className="gauge-main">
        <div className="gauge-circle">
          <svg viewBox="0 0 100 100">
            <circle cx="50" cy="50" r="40" fill="none" stroke="#e5e7eb" strokeWidth="8" />
            <circle
              cx="50"
              cy="50"
              r="40"
              fill="none"
              stroke={getQualityColor(quality)}
              strokeWidth="8"
              strokeLinecap="round"
              strokeDasharray={`${(quality / 100) * 251.2} 251.2`}
              transform="rotate(-90 50 50)"
            />
          </svg>
          <div className="gauge-center">
            <span className="gauge-percent">{qualityPercent}</span>
            <span className="gauge-unit">%</span>
          </div>
        </div>
        <div className="gauge-info">
          <span className="gauge-label" style={{ color: getQualityColor(quality) }}>
            {getQualityStatus(quality)}
          </span>
          <span className="gauge-sublabel">Assembly Readiness</span>
        </div>
      </div>

      <div className="gauge-metrics">
        <div className="gauge-metric">
          <div className="metric-icon">
            <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
              <path d="M4 6h16v2H4zm0 5h16v2H4zm0 5h16v2H4z"/>
            </svg>
          </div>
          <div className="metric-content">
            <span className="metric-value">{overlapLength} bp</span>
            <span className="metric-label">Overlap</span>
          </div>
        </div>
        <div className="gauge-metric">
          <div className="metric-icon" style={{ color: tmStatus.color }}>
            <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
              <path d="M15 13V5c0-1.66-1.34-3-3-3S9 3.34 9 5v8c-1.21.91-2 2.37-2 4 0 2.76 2.24 5 5 5s5-2.24 5-5c0-1.63-.79-3.09-2-4zm-4-8c0-.55.45-1 1-1s1 .45 1 1h-1v1h1v2h-1v1h1v2h-2V5z"/>
            </svg>
          </div>
          <div className="metric-content">
            <span className="metric-value" style={{ color: tmStatus.color }}>{targetTm}°C</span>
            <span className="metric-label">Target Tm</span>
          </div>
        </div>
        <div className="gauge-metric">
          <div className="metric-icon">
            <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
              <path d="M4 8h4V4H4v4zm6 12h4v-4h-4v4zm-6 0h4v-4H4v4zm0-6h4v-4H4v4zm6 0h4v-4h-4v4zm6-10v4h4V4h-4zm-6 4h4V4h-4v4zm6 6h4v-4h-4v4zm0 6h4v-4h-4v4z"/>
            </svg>
          </div>
          <div className="metric-content">
            <span className="metric-value">{fragmentCount}</span>
            <span className="metric-label">Fragments</span>
          </div>
        </div>
      </div>

      <div className="gauge-footer">
        <span className="method-badge">{isHiFi ? 'NEBuilder HiFi' : 'Gibson'}</span>
        <span className="assembly-type">Isothermal Assembly</span>
      </div>
    </div>
  );
}

// Overlap Analysis Display for Gibson/HiFi results
function OverlapAnalysisDisplay({ overlaps }: { overlaps: any }) {
  if (!overlaps || overlaps.length === 0) return null;

  const getScoreColor = (score: number) => {
    if (score >= 0.8) return '#22c55e';
    if (score >= 0.6) return '#84cc16';
    if (score >= 0.4) return '#f59e0b';
    return '#ef4444';
  };

  return (
    <div className="overlap-analysis">
      <h4>Junction Analysis</h4>
      <div className="overlap-list">
        {overlaps.map((overlap: any, i: number) => (
          <div key={i} className="overlap-item">
            <span className="overlap-label">Junction {i + 1}</span>
            <div className="overlap-details">
              <span className="overlap-length">{overlap.length} bp</span>
              <span className="overlap-tm">Tm: {overlap.tm?.toFixed(1) || '—'}°C</span>
              <span className="overlap-gc">GC: {overlap.gc?.toFixed(0) || '—'}%</span>
            </div>
            <div className="overlap-score-bar">
              <div
                className="score-fill"
                style={{
                  width: `${(overlap.score || 0.5) * 100}%`,
                  backgroundColor: getScoreColor(overlap.score || 0.5)
                }}
              />
            </div>
            {overlap.warnings && overlap.warnings.length > 0 && (
              <div className="overlap-warnings">
                {overlap.warnings.map((w: any, j: number) => (
                  <span key={j} className="warning-tag">{w}</span>
                ))}
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
}

// Main Golden Gate Designer Component (now unified as Assembly Studio)
export default function GoldenGateDesigner() {
  // Method selection - now uses 2 main methods
  const [method, setMethod] = useState<string>('golden_gate');
  // Protocol variant for isothermal assembly (gibson or nebuilder_hifi)
  const [isothermalVariant, setIsothermalVariant] = useState<string>('nebuilder_hifi');

  // Computed helpers
  const isGoldenGate = method === 'golden_gate';
  const isIsothermal = method === 'isothermal';
  // For backwards compatibility with existing code that checks specific methods
  const effectiveMethod = isIsothermal ? isothermalVariant : method;

  // Fusion Site Optimizer mode - enables automatic junction position finding
  // (Advanced feature - not commonly used)
  const [fusionMode, setFusionMode] = useState<boolean>(false);
  const [fusionSequence, setFusionSequence] = useState<string>('');
  const [fusionNumFragments, setFusionNumFragments] = useState<number>(4);
  const [fusionCandidates, setFusionCandidates] = useState<FusionCandidate[]>([]);
  const [fusionResult, setFusionResult] = useState<any>(null);
  const [isOptimizing, setIsOptimizing] = useState<boolean>(false);

  // Auto-Domestication - automatically breaks internal restriction sites
  const [autoDomesticationEnabled, setAutoDomesticationEnabled] = useState<boolean>(true);
  const [minFragmentSize, setMinFragmentSize] = useState<number>(DOMESTICATION_DEFAULTS.minFragmentSize);
  const [domesticationError, setDomesticationError] = useState<DomesticationAnalysis | null>(null);

  // Enhanced Domestication - state-of-the-art workflow with user confirmation
  const [showEnhancedDomestication, setShowEnhancedDomestication] = useState<boolean>(false);
  const [enhancedDomesticationPart, setEnhancedDomesticationPart] = useState<Part | null>(null);
  const [useEnhancedWorkflow, setUseEnhancedWorkflow] = useState<boolean>(true); // New workflow enabled by default
  const [needsRedesignAfterDomestication, setNeedsRedesignAfterDomestication] = useState<boolean>(false);

  // Domestication Workflow Guide - complete workflow with primer design
  const [showWorkflowGuide, setShowWorkflowGuide] = useState<boolean>(false);
  const [workflowGuidePart, setWorkflowGuidePart] = useState<Part | null>(null);

  // Primer Boundary Optimizer - adjusts boundaries for better primers
  const [boundaryAssessment, setBoundaryAssessment] = useState<any>(null);
  const [isOptimizingBoundaries, setIsOptimizingBoundaries] = useState<boolean>(false);
  const [showAdvancedOptions, setShowAdvancedOptions] = useState<boolean>(false);

  // Minimum sequence length for PCR-based Golden Gate assembly:
  // - Need ~18-25bp homology region for each primer (for proper Tm)
  // - Both primers need non-overlapping binding sites
  // - NEB recommends fragments >100bp for reliable PCR amplification
  // - Fragments <50bp should use synthetic oligo approach instead
  const MIN_SEQUENCE_LENGTH = 50;
  const RECOMMENDED_MIN_LENGTH = 100;

  // Start with empty parts - user can click "Load Example" to load sample data
  const [parts, setParts] = useState<Part[]>([
    { id: '', seq: '', type: 'promoter' },
    { id: '', seq: '', type: 'cds' },
  ]);
  const [enzyme, setEnzyme] = useState<string>('BsaI');
  const [result, setResult] = useState<AssemblyResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [copyMessage, setCopyMessage] = useState<string>('');
  const [draggedIndex, setDraggedIndex] = useState<number | null>(null);
  const [viewMode, setViewMode] = useState<string>('linear'); // 'linear' or 'circular'

  // Gibson/HiFi specific settings - defaults optimized for NEBuilder HiFi
  const [overlapSettings, setOverlapSettings] = useState<OverlapSettings>({
    overlapLength: 20,  // NEB recommends 15-20bp for 2-3 fragments, 20-30bp for 4-6
    targetTm: 52,       // Must be ≥48°C (reaction at 50°C)
  });
  const [showIsothermalAdvanced, setShowIsothermalAdvanced] = useState<boolean>(false);
  // Use advanced range-based optimization for isothermal assembly
  const [useAdvancedOptimization, setUseAdvancedOptimization] = useState<boolean>(true);

  // Computed values - maxParts based on method
  const maxParts = isGoldenGate ? MAX_GOLDEN_GATE_FRAGMENTS : MAX_FRAGMENTS;

  const recommendedOverhangs = useMemo(() => {
    if (!isGoldenGate) return { overhangs: [], fidelity: 1.0 };
    return getRecommendedOverhangs(parts.length);
  }, [parts.length, isGoldenGate]);

  const totalLength = useMemo(() => {
    return parts.reduce((sum: number, p: Part) => sum + (p.seq?.length || 0), 0);
  }, [parts]);

  const handlePartChange = useCallback((index: number, newPart: Part) => {
    setParts(prev => {
      const updated = [...prev];
      updated[index] = newPart;
      return updated;
    });
    setResult(null);
  }, []);

  const addPart = useCallback(() => {
    if (parts.length < maxParts) {
      const types: Array<'promoter' | 'rbs' | 'cds' | 'terminator' | 'other'> = ['promoter', 'rbs', 'cds', 'terminator', 'other'];
      setParts(prev => [...prev, { id: '', seq: '', type: types[prev.length % types.length] || 'other' }] as Part[]);
      setResult(null);
    }
  }, [parts.length, maxParts]);

  const removePart = useCallback((index: any) => {
    // Allow deletion of any part - user can add parts back or load example
    setParts(prev => prev.filter((_: any, i: number) => i !== index));
    setResult(null);
  }, []);

  /**
   * Handle splitting a single large sequence into optimized fragments
   * This runs the Fusion Site Optimizer and replaces the part with fragment parts
   * Includes auto-domestication support to break internal enzyme sites
   */
  const handleSplitSequenceFromPart = useCallback((partIndex: any, sequence: any, numFragments: any, enzymeToUse: any) => {
    return new Promise((resolve, reject) => {
      try {
        // First, analyze for domestication and check for errors
        const domesticationAnalysis = analyzeForDomestication(sequence, enzymeToUse, {
          minFragmentSize,
          minSiteDistance: minFragmentSize, // Use same value for site distance
        });

        // Check for domestication errors (adjacent sites, fragment size violations, etc.)
        if (domesticationAnalysis.error) {
          setDomesticationError(domesticationAnalysis as any);
          reject(new Error(domesticationAnalysis.message));
          return;
        }

        // Clear any previous domestication errors
        setDomesticationError(null);

        // Check for auto-domestication needs
        const domesticationInfo = optimizeWithDomestication(
          sequence,
          numFragments,
          enzymeToUse,
          { constraints: { minFragmentSize, maxFragmentSize: Math.max(3000, Math.ceil(sequence.length / 2)) } }
        );

        // Use effective fragment count (includes domestication junctions)
        const effectiveFragments = domesticationInfo.needsDomestication
          ? domesticationInfo.totalFragments
          : numFragments;

        // Run the optimizer with effective fragment count
        const result = optimizeFusionSites(sequence, effectiveFragments, {
          enzyme: enzymeToUse,
          algorithm: 'auto',
          minFragmentSize,
          maxFragmentSize: Math.max(3000, Math.ceil(sequence.length / 2)),
          minDistanceFromEnds: 50,
          domesticationJunctions: domesticationInfo.domesticationJunctions || [],
        } as any);

        if (!result.success) {
          reject(new Error(result.error || 'Failed to find valid junction positions'));
          return;
        }

        // Log domestication info if applied
        if (domesticationInfo.needsDomestication) {
          console.log(`Auto-domestication: Split includes ${domesticationInfo.additionalFragmentsNeeded} ` +
            `additional junction(s) to break internal ${enzymeToUse} sites`);
        }

        // Convert junctions to fragment parts
        const junctions = [...(result.junctions || [])].sort((a: any, b: any) => a.position - b.position);
        const fragmentParts: any[] = [];

        // Get enzyme-specific overhang length (4 for BsaI/BsmBI/BbsI, 3 for SapI)
        const overhangLength = GOLDEN_GATE_ENZYMES[enzymeToUse]?.overhangLength || 4;

        for (let i = 0; i <= junctions.length; i++) {
          const isFirst = i === 0;
          const isLast = i === junctions.length;

          // Calculate fragment boundaries
          const fragStart = isFirst ? 0 : junctions[i - 1].position;
          const fragEnd = isLast ? sequence.length : junctions[i].position + overhangLength;

          // Extract fragment sequence
          const fragSeq = sequence.slice(fragStart, fragEnd);

          fragmentParts.push({
            id: `Fragment_${i + 1}`,
            seq: fragSeq,
            type: 'fragment',
            // Store optimization metadata
            _optimized: {
              leftOverhang: isFirst ? null : junctions[i - 1].overhang,
              rightOverhang: isLast ? null : junctions[i].overhang,
              junctionQuality: isFirst ? null : junctions[i - 1].score?.composite,
            },
          });
        }

        // Replace the single part with the fragment parts
        setParts(prev => {
          const updated = [...prev];
          // Remove the original part and insert fragments at that position
          updated.splice(partIndex, 1, ...fragmentParts);
          return updated;
        });

        setResult(null);
        const domesticationNote = domesticationInfo.needsDomestication
          ? ` (includes ${domesticationInfo.additionalFragmentsNeeded} domestication junction${domesticationInfo.additionalFragmentsNeeded > 1 ? 's' : ''})`
          : '';
        setCopyMessage(`Split into ${effectiveFragments} optimized fragments${domesticationNote} (${((result.score?.fidelity || 0.9) * 100).toFixed(0)}% predicted fidelity)`);
        resolve({ ...result, domesticationInfo });
      } catch (err: unknown) {
        reject(err);
      }
    });
  }, [minFragmentSize]);

  // Drag and drop handlers
  const handleDragStart = (e: any, index: any) => {
    setDraggedIndex(index);
    e.dataTransfer.effectAllowed = 'move';
  };

  const handleDragOver = (e: any, index: any) => {
    e.preventDefault();
    e.dataTransfer.dropEffect = 'move';
  };

  const handleDrop = (e: any, dropIndex: any) => {
    e.preventDefault();
    if (draggedIndex === null || draggedIndex === dropIndex) return;

    setParts(prev => {
      const newParts = [...prev];
      const [removed] = newParts.splice(draggedIndex, 1);
      newParts.splice(dropIndex, 0, removed);
      return newParts;
    });
    setDraggedIndex(null);
    setResult(null);
  };

  const handleFileUpload = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const file = (e.target as HTMLInputElement).files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event: ProgressEvent<FileReader>) => {
      try {
        const content = event.target?.result as string;
        const parsed = parseFasta(content);

        if (parsed.length === 0) {
          setError('No sequences found in file');
          return;
        }

        const types: Array<'promoter' | 'rbs' | 'cds' | 'terminator' | 'other'> = ['promoter', 'rbs', 'cds', 'terminator', 'other'];
        const newParts = parsed.slice(0, MAX_FRAGMENTS).map((p: any, i: number) => ({
          id: p.id,
          seq: p.seq.toUpperCase(),
          type: (types[i] || 'other') as 'promoter' | 'rbs' | 'cds' | 'terminator' | 'other',
        }));

        setParts((newParts.length >= 2 ? newParts : [...newParts, { id: '', seq: '', type: 'other' }]) as Part[]);
        setResult(null);
        setError(null);
      } catch (err: unknown) {
        setError('Failed to parse file: ' + (err as Error).message);
      }
    };
    reader.readAsText(file);
  }, []);

  const designAssembly = useCallback(() => {
    setError(null);
    try {
      // Apply range selection to each part
      const partsWithRange = parts.map(p => {
        if (!p.seq) return p;
        const start = p.rangeStart || 1;
        const end = p.rangeEnd || p.seq.length;
        const selectedSeq = p.seq.slice(start - 1, end);
        return {
          ...p,
          seq: selectedSeq,
          originalSeq: p.seq,
          rangeStart: start,
          rangeEnd: end,
        };
      });

      // Validate sequence lengths for PCR-based assembly
      // Golden Gate requires primers to bind to both ends, so need sufficient length
      const tooShortParts = partsWithRange.filter(p => p.seq && p.seq.length > 0 && p.seq.length < MIN_SEQUENCE_LENGTH);
      if (tooShortParts.length > 0) {
        const names = tooShortParts.map(p => p.id || 'unnamed').join(', ');
        setError(`Parts too short for PCR amplification (minimum ${MIN_SEQUENCE_LENGTH}bp): ${names}. For fragments <50bp, consider ordering synthetic oligos with overhangs instead.`);
        return;
      }

      // Warn about short sequences that may have issues
      const shortParts = partsWithRange.filter(p => p.seq && p.seq.length >= MIN_SEQUENCE_LENGTH && p.seq.length < RECOMMENDED_MIN_LENGTH);

      const validParts = partsWithRange.filter(p => p.seq && p.seq.length >= MIN_SEQUENCE_LENGTH);

      if (validParts.length < 2) {
        setError(`At least 2 parts with sequences (minimum ${MIN_SEQUENCE_LENGTH}bp each) required`);
        return;
      }

      let assemblyResult;
      let domesticationApplied = false;
      let expandedParts = validParts;
      let domesticationDetails = []; // Track mutation info for UI display

      if (isGoldenGate && autoDomesticationEnabled) {
        // Auto-domestication using unified optimizer (one-pot compatible strategies)
        // Priority: Direct primer mutation > Mutagenic junction > Alternative enzyme
        const partsToExpand = [];

        for (const part of validParts) {
          // Check if user has already configured domestication for this part
          console.log('[DesignAssembly] Checking part:', part.id, {
            _domesticationApproved: part._domesticationApproved,
            _domesticationMutations: part._domesticationMutations,
            _domesticationJunctions: part._domesticationJunctions,
          });
          if (part._domesticationApproved) {
            // User has configured domestication through EnhancedDomesticationPanel
            // Use their configured strategy and primers instead of auto-detecting
            // Check the actual data (junctions/mutations) to determine the approach - more robust than string matching
            const configuredJunctions = part._domesticationJunctions || [];
            const configuredMutations = part._domesticationMutations || [];
            const configuredPrimers = part._domesticationPrimers || [];
            const hasJunctions = configuredJunctions.length > 0;
            const hasMutations = configuredMutations.length > 0;

            // Debug log to help trace the issue
            console.log('[Domestication] Part:', part.id, {
              strategy: part._domesticationStrategy,
              hasJunctions,
              hasMutations,
              junctionsCount: configuredJunctions.length,
              mutationsCount: configuredMutations.length,
            });

            if (hasJunctions) {
              // User chose mutagenic junction strategy - create fragments
              const junctions = [...configuredJunctions].sort((a, b) =>
                (a.junctionPosition ?? a.position ?? 0) - (b.junctionPosition ?? b.position ?? 0)
              );
              const domesticationOverhangLen = GOLDEN_GATE_ENZYMES[enzyme]?.overhangLength || 4;

              // Track mutations for display
              const partMutations = junctions.map(junction => ({
                sitePosition: junction.sitePosition ?? junction.site?.position,
                junctionPosition: junction.junctionPosition ?? junction.position,
                overhang: junction.overhang,
                strategy: 'mutagenic_junction',
              }));

              // Create fragment parts from configured junctions
              for (let i = 0; i <= junctions.length; i++) {
                const isFirst = i === 0;
                const isLast = i === junctions.length;
                const juncPos = (j: any) => j.junctionPosition ?? j.position ?? 0;
                const fragStart = isFirst ? 0 : juncPos(junctions[i - 1]);
                const fragEnd = isLast ? part.seq.length : juncPos(junctions[i]) + domesticationOverhangLen;
                const fragSeq = part.seq.slice(fragStart, fragEnd);

                partsToExpand.push({
                  id: `${part.id || 'Part'}_frag${i + 1}`,
                  seq: fragSeq,
                  type: part.type,
                  _domesticated: true,
                  _userConfigured: true,
                  _parentPart: part.id,
                  _fragmentIndex: i,
                  _totalFragments: junctions.length + 1,
                  _overhang: isLast ? null : junctions[i].overhang,
                  _strategy: 'mutagenic_junction',
                  _onePotCompatible: true,
                  _configuredPrimers: configuredPrimers,
                });
              }

              domesticationDetails.push({
                partId: part.id,
                siteCount: junctions.length,
                strategy: 'mutagenic_junction',
                mutations: partMutations,
                fragmentCount: junctions.length + 1,
                onePotCompatible: true,
                userConfigured: true,
              });

              domesticationApplied = true;
              continue;
            } else if (hasMutations) {
              // User chose silent mutation strategy - sequence is already modified
              // The part.seq should already contain the domesticated sequence
              console.log('[Domestication] Entering silent mutation branch for part:', part.id);
              partsToExpand.push({
                ...part,
                _domesticated: true,
                _userConfigured: true,
                _strategy: 'silent_mutation',
                _mutations: configuredMutations,
                _onePotCompatible: true,
                _configuredPrimers: configuredPrimers,
              });

              domesticationDetails.push({
                partId: part.id,
                siteCount: configuredMutations.length,
                strategy: 'silent_mutation',
                mutations: configuredMutations.map((m: any) => ({
                  // Handle nested structure: mutation data can be in m.mutation or directly in m
                  position: m.site?.position ?? m.mutation?.position ?? m.position ?? m.sequencePosition,
                  codonChange: m.codonChange || m.mutation?.codonChange || `${m.mutation?.originalCodon || m.originalCodon || ''}→${m.mutation?.newCodon || m.newCodon || ''}`,
                  aminoAcid: m.aminoAcid || m.mutation?.aminoAcid,
                  strategy: 'silent_mutation',
                })),
                fragmentCount: 1,
                onePotCompatible: true,
                userConfigured: true,
                workflowSteps: [
                  {
                    step: 1,
                    title: 'Site-Directed Mutagenesis',
                    description: `Introduce ${configuredMutations.length} silent mutation(s) to remove internal ${enzyme} sites`,
                    type: 'pcr_mutagenesis',
                    primers: configuredPrimers.filter(p => p.type === 'pcr_mutagenesis'),
                  },
                  {
                    step: 2,
                    title: 'Golden Gate Assembly',
                    description: 'Assemble the domesticated template with other parts',
                    type: 'golden_gate',
                    primers: configuredPrimers.filter(p => p.type === 'golden_gate' || p.type === 'junction'),
                  },
                ],
              });

              domesticationApplied = true;
              continue;
            } else {
              // Part was approved but no junctions or mutations - sequence doesn't need changes
              console.log('[Domestication] Part approved but no changes needed:', part.id);
              partsToExpand.push(part);
              continue;
            }
          }

          // Use new unified domestication analysis (auto-detection)
          const analysis = analyzeDomesticationOptions(part.seq, enzyme, {
            frame: 0, // Assume frame 0, could be made configurable
            organism: 'ecoli',
          });

          if (analysis.needsDomestication) {
            try {
              // Use mutagenic junction approach - one-pot compatible!
              const mutagenicResult = designAllMutagenicJunctions(part.seq, enzyme, {
                frame: 0,
                organism: 'ecoli',
              });

              if (mutagenicResult.success && mutagenicResult.junctions && mutagenicResult.junctions.length > 0) {
                // Create fragments with mutagenic primers
                // Note: junctions have junctionPosition property, not position
                const junctions = mutagenicResult.junctions.sort((a: any, b: any) => a.junctionPosition - b.junctionPosition);
                const domesticationOverhangLen = GOLDEN_GATE_ENZYMES[enzyme]?.overhangLength || 4;

                // Track mutations for this part
                const partMutations = [];
                for (const junction of junctions) {
                  if (junction.mutation) {
                    partMutations.push({
                      sitePosition: junction.site?.position,
                      junctionPosition: junction.junctionPosition,
                      overhang: junction.overhang,
                      mutations: [{
                        position: junction.mutation.sequencePosition,
                        originalBase: junction.mutation.originalBase,
                        newBase: junction.mutation.newBase,
                        codon: junction.mutation.originalCodon,
                        inPrimer: junction.mutation.inFragment1Primer ? 'reverse' : 'forward',
                      }],
                      strategy: DOMESTICATION_STRATEGY.MUTAGENIC_JUNCTION,
                    });
                  }
                }

                // Create fragment parts from junctions
                for (let i = 0; i <= junctions.length; i++) {
                  const isFirst = i === 0;
                  const isLast = i === junctions.length;
                  const fragStart = isFirst ? 0 : junctions[i - 1].junctionPosition;
                  const fragEnd = isLast ? part.seq.length : (junctions[i]?.junctionPosition ?? 0) + domesticationOverhangLen;
                  const fragSeq = part.seq.slice(fragStart, fragEnd);

                  partsToExpand.push({
                    id: `${part.id || 'Part'}_frag${i + 1}`,
                    seq: fragSeq,
                    type: part.type,
                    _domesticated: true,
                    _parentPart: part.id,
                    _fragmentIndex: i,
                    _totalFragments: junctions.length + 1,
                    _overhang: isLast ? null : junctions[i].overhang,
                    _strategy: DOMESTICATION_STRATEGY.MUTAGENIC_JUNCTION,
                    _mutation: !isLast ? junctions[i].mutation : null,
                    _primers: !isLast ? junctions[i].primers : null,
                    _onePotCompatible: true, // Mutagenic junctions are one-pot compatible!
                  });
                }

                domesticationDetails.push({
                  partId: part.id,
                  siteCount: analysis.siteCount,
                  strategy: DOMESTICATION_STRATEGY.MUTAGENIC_JUNCTION,
                  mutations: partMutations,
                  fragmentCount: junctions.length + 1,
                  onePotCompatible: true,
                });

                domesticationApplied = true;
                continue; // Skip adding original part
              }

              // Fallback: Check if direct primer mutations work (sites near existing junctions)
              const allSilentAvailable = analysis.siteAnalyses?.every(a => a.silentMutation?.available);
              if (allSilentAvailable && analysis.siteAnalyses) {
                // All sites can be handled with direct primer mutations - no extra fragments needed!
                const mutations = analysis.siteAnalyses
                  .filter(a => a.silentMutation?.bestCandidate)
                  .map(a => ({
                    sitePosition: a.site.position,
                    mutation: a.silentMutation.bestCandidate,
                    strategy: DOMESTICATION_STRATEGY.DIRECT_PRIMER_MUTATION,
                  }));

                // Part stays as-is, mutations will be applied during primer design
                partsToExpand.push({
                  ...part,
                  _domesticated: true,
                  _strategy: DOMESTICATION_STRATEGY.DIRECT_PRIMER_MUTATION,
                  _mutations: mutations,
                  _onePotCompatible: true,
                });

                domesticationDetails.push({
                  partId: part.id,
                  siteCount: analysis.siteCount,
                  strategy: DOMESTICATION_STRATEGY.DIRECT_PRIMER_MUTATION,
                  mutations,
                  fragmentCount: 1, // No additional fragments
                  onePotCompatible: true,
                });

                domesticationApplied = true;
                continue;
              }

              // If mutagenic junction failed, fall back to legacy approach but warn user
              const legacyAnalysis = analyzeForDomestication(part.seq, enzyme);
              if (legacyAnalysis.needsDomestication && (legacyAnalysis.domesticationOptions?.length || 0) > 0) {
                const domesticationJunctions = (legacyAnalysis.domesticationOptions || [])
                  .filter((opt: any) => opt.hasValidOption && opt.recommended)
                  .map((opt: any) => ({
                    position: opt.recommended?.position,
                    overhang: opt.recommended?.overhang,
                  }))
                  .sort((a, b) => a.position - b.position);

                if (domesticationJunctions.length > 0) {
                  const domesticationOverhangLen = GOLDEN_GATE_ENZYMES[enzyme]?.overhangLength || 4;

                  for (let i = 0; i <= domesticationJunctions.length; i++) {
                    const isFirst = i === 0;
                    const isLast = i === domesticationJunctions.length;
                    const fragStart = isFirst ? 0 : domesticationJunctions[i - 1].position;
                    const fragEnd = isLast ? part.seq.length : domesticationJunctions[i].position + domesticationOverhangLen;
                    const fragSeq = part.seq.slice(fragStart, fragEnd);

                    partsToExpand.push({
                      id: `${part.id || 'Part'}_frag${i + 1}`,
                      seq: fragSeq,
                      type: part.type,
                      _domesticated: true,
                      _parentPart: part.id,
                      _fragmentIndex: i,
                      _totalFragments: domesticationJunctions.length + 1,
                      _overhang: isLast ? null : domesticationJunctions[i].overhang,
                      _strategy: DOMESTICATION_STRATEGY.LEGACY_JUNCTION,
                      _onePotCompatible: false, // Legacy is NOT one-pot compatible
                    });
                  }

                  domesticationDetails.push({
                    partId: part.id,
                    siteCount: (legacyAnalysis as any).internalSites?.length || 1,
                    strategy: DOMESTICATION_STRATEGY.LEGACY_JUNCTION,
                    fragmentCount: domesticationJunctions.length + 1,
                    onePotCompatible: false,
                    warning: 'Using legacy junction-based domestication. NOT compatible with one-pot reactions.',
                  });

                  domesticationApplied = true;
                  continue;
                }
              }
            } catch (err: unknown) {
              console.warn(`Auto-domestication failed for ${part.id}:`, err);
              // Fall through to add original part
            }
          }

          // No internal sites or domestication failed - keep original part
          partsToExpand.push(part);
        }

        expandedParts = partsToExpand;

        if (domesticationApplied) {
          const onePotCount = domesticationDetails.filter((d: any) => d.onePotCompatible).length;
          const legacyCount = domesticationDetails.filter((d: any) => !d.onePotCompatible).length;
          console.log(`Auto-domestication: Expanded ${validParts.length} parts to ${expandedParts.length} fragments`);
          console.log(`  One-pot compatible: ${onePotCount}, Legacy (requires sequential): ${legacyCount}`);
        }
      }

      if (isGoldenGate) {
        // Golden Gate Assembly - using optimized primer design with:
        // - 6bp flanking sequences (NEB recommendation)
        // - G:T mismatch detection
        // - Automatic overhang optimization for maximum fidelity

        // Collect domestication overhangs from expanded fragments
        // These must be preserved - they are specifically designed to break internal restriction sites
        const requiredOverhangPatterns: string[] = [];
        const requiredOverhangIndices: number[] = [];
        expandedParts.forEach((part: Part, i: number) => {
          // _overhang is the right-side overhang of this fragment (junction to next fragment)
          // This corresponds to overhang index i+1 in the assembly
          if ((part as any)._overhang) {
            requiredOverhangPatterns.push((part as any)._overhang);
            requiredOverhangIndices.push(i + 1);
          }
        });

        assemblyResult = designOptimizedGoldenGateAssemblyForUI(expandedParts, {
          enzyme,
          circular: true,
          useStandardOverhangs: requiredOverhangPatterns.length === 0, // Don't use standard if we have domestication overhangs
          autoOptimize: true,
          autoDomestication: domesticationApplied, // Tell primer optimizer about domestication
          requiredOverhangPatterns: requiredOverhangPatterns.length > 0 ? requiredOverhangPatterns : undefined,
          requiredOverhangIndices: requiredOverhangIndices.length > 0 ? requiredOverhangIndices : undefined,
        });

        // Add domestication info to result
        if (domesticationApplied) {
          console.log('[Domestication] domesticationDetails:', JSON.stringify(domesticationDetails, null, 2));
          console.log('[Domestication] Has workflowSteps:', domesticationDetails.some((d: any) => d.workflowSteps?.length > 0));

          const onePotCompatible = domesticationDetails.every((d: any) => d.onePotCompatible);
          const strategies = [...new Set(domesticationDetails.map((d: any) => d.strategy))];

          (assemblyResult as any)._autoDomestication = {
            applied: true,
            originalParts: validParts.length,
            expandedParts: expandedParts.length,
            additionalFragments: expandedParts.length - validParts.length,
            onePotCompatible,
            strategies,
            details: domesticationDetails,
            // Summarize mutations for display
            mutationSummary: domesticationDetails
              .filter((d: any) => d.mutations && d.mutations.length > 0)
              .flatMap((d: any) => d.mutations)
              .map((m: any) => ({
                partId: (m as any).partId,
                sitePosition: (m as any).sitePosition,
                strategy: (m as any).strategy,
                baseChanges: (m as any).mutations ? (m as any).mutations.map((mut: any) =>
                  `${mut.originalBase}${mut.position}${mut.newBase}`
                ).join(', ') : null,
              })),
          };

          // Add warning if legacy fallback was used
          if (!onePotCompatible) {
            console.warn('Auto-domestication: Some parts use legacy junction-based approach. ' +
                        'Consider using alternative enzyme or sequential protocol.');
          }
        }
      } else {
        // Gibson or NEBuilder HiFi Assembly
        const fragments = validParts.map((p: Part, i: number) => ({
          id: p.id || `Fragment_${i + 1}`,
          seq: p.seq,
          type: p.type,
        }));

        const designResult = designNEBuilderAssembly(fragments, {
          circular: true,
          targetOverlapLength: overlapSettings.overlapLength,
          targetOverlapTm: overlapSettings.targetTm,
          method: effectiveMethod === 'nebuilder_hifi' ? 'nebuilder' : 'gibson',
        } as any);

        // Simulate the assembly - pass fragments (with seq) and junctions
        const simulation = simulateAssembly(
          fragments,
          designResult.junctions as any,
          { circular: true }
        );

        // Transform to unified result format for display
        // Use designResult.fragments for primer info, validParts for sequences
        assemblyResult = {
          method: effectiveMethod,
          parts: designResult.fragments.map((frag: any, i: number) => {
            const part = validParts[i];
            return {
              id: frag.id || part?.id || `Part ${i + 1}`,
              type: part?.type || 'other',
              seq: part?.seq || '',
              length: part?.seq?.length || frag.length || 0,
              leftOverhang: frag.forward?.homologyTail || '',
              rightOverhang: frag.reverse?.homologyTail || '',
              primers: {
                forward: {
                  sequence: frag.forward?.sequence || '',
                  length: frag.forward?.length || 0,
                  tm: frag.forward?.tm || 60,
                  gc: 50,
                  structure: {
                    extra: '',
                    bsaISite: '',
                    spacer: '',
                    overhang: '',
                    homology: frag.forward?.homologyTail || '',
                    annealing: frag.forward?.annealingRegion || '',
                  },
                },
                reverse: {
                  sequence: frag.reverse?.sequence || '',
                  length: frag.reverse?.length || 0,
                  tm: frag.reverse?.tm || 60,
                  gc: 50,
                  structure: {
                    extra: '',
                    bsaISite: '',
                    spacer: '',
                    overhang: '',
                    homology: frag.reverse?.homologyTail || '',
                    annealing: frag.reverse?.annealingRegion || '',
                  },
                },
                pcr: {
                  annealingTemp: frag.pair?.annealingTemp || 58,
                  extensionTime: Math.ceil((part?.seq?.length || 1000) / 1000) * 30,
                  tmDifference: frag.pair?.tmDiff || 0,
                },
              },
            };
          }),
          assembledSequence: simulation?.assembledSequence || validParts.map(p => p.seq).join(''),
          assembledLength: simulation?.assembledLength || validParts.reduce((sum: number, p: Part) => sum + (p.seq?.length || 0), 0),
          fidelity: {
            overall: 0.95,
            percentage: '95.0%',
            individual: (designResult.junctions || []).map((j, i) => ({
              overhang: j.optimal?.sequence || '',
              fidelity: j.optimal?.compositeScore || 0.9,
            })),
          },
          overlapAnalysis: (designResult.junctions || []).map(j => ({
            length: j.optimal?.length || overlapSettings.overlapLength,
            tm: j.optimal?.tm || overlapSettings.targetTm,
            gc: j.optimal?.gc || 50,
            score: j.optimal?.compositeScore || 0.9,
            warnings: j.optimal?.warnings || [],
          })),
          warnings: designResult.warnings || [],
          protocol: generateProtocol(effectiveMethod, validParts.length),
          simulation: simulation,
        };
      }

      setResult(assemblyResult as any);

      // After design, assess if boundary optimization could help (Golden Gate only)
      if (isGoldenGate && validParts.length >= 2) {
        try {
          const assessment = assessBoundaryOptimizationPotential(validParts);
          setBoundaryAssessment(assessment);
        } catch (e: unknown) {
          console.warn('Boundary assessment failed:', e);
          setBoundaryAssessment(null);
        }
      } else {
        setBoundaryAssessment(null);
      }
    } catch (err: unknown) {
      console.error('Assembly design error:', err);
      setError((err as Error).message);
      setResult(null);
      setBoundaryAssessment(null);
    }
  }, [parts, enzyme, method, isothermalVariant, overlapSettings, autoDomesticationEnabled, isGoldenGate, effectiveMethod]);

  // Auto-redesign after domestication configuration is saved
  useEffect(() => {
    if (needsRedesignAfterDomestication) {
      console.log('[useEffect] Triggering redesign after domestication config saved');
      setNeedsRedesignAfterDomestication(false);
      designAssembly();
    }
  }, [needsRedesignAfterDomestication, designAssembly]);

  // Handler for boundary optimization
  const handleOptimizeBoundaries = useCallback(async () => {
    if (!result || !isGoldenGate || parts.length < 2) return;

    setIsOptimizingBoundaries(true);
    setError(null);

    try {
      const validParts = parts.filter(p => p.seq && p.seq.length >= MIN_SEQUENCE_LENGTH);
      const optimized = optimizeAssemblyBoundaries(validParts, {
        maxShift: 50,
        minFragmentSize: MIN_SEQUENCE_LENGTH,
      });

      if (optimized.success && (optimized.summary?.boundariesOptimized || 0) > 0) {
        // Update parts with optimized sequences
        const newParts = optimized.optimizedFragments.map((frag: any, i: number) => ({
          ...parts[i],
          seq: frag.seq,
          id: frag.id || parts[i].id,
          _boundaryOptimized: frag.lengthChange !== 0,
          _originalSeq: frag.originalSeq,
        }));

        setParts(newParts);

        // Re-run primer design with optimized boundaries
        setTimeout(() => {
          designAssembly();
        }, 100);
      } else {
        setError('No boundary improvements found. Current boundaries are already optimal.');
      }
    } catch (err: unknown) {
      console.error('Boundary optimization error:', err);
      setError('Failed to optimize boundaries: ' + (err as Error).message);
    } finally {
      setIsOptimizingBoundaries(false);
    }
  }, [parts, result, isGoldenGate, designAssembly]);

  // Generate protocol based on method
  const generateProtocol = (methodOrVariant: any, partCount: any) => {
    if (methodOrVariant === 'golden_gate') {
      return null; // Golden Gate uses its own protocol from designGoldenGateAssembly
    }

    const isHiFi = methodOrVariant === 'nebuilder_hifi';
    return {
      title: isHiFi ? 'NEBuilder HiFi DNA Assembly Protocol' : 'Gibson Assembly Protocol',
      steps: [
        {
          title: 'PCR Amplification',
          details: [
            `Amplify ${partCount} fragments with the designed primers`,
            'Use high-fidelity polymerase (Q5 or equivalent)',
            'Verify PCR products by gel electrophoresis',
            'Gel purify or PCR cleanup the products',
          ],
        },
        {
          title: 'Assembly Reaction',
          details: isHiFi ? [
            'Add 0.03-0.2 pmol of each fragment (equimolar amounts)',
            'Add NEBuilder HiFi DNA Assembly Master Mix (2X)',
            'Bring total volume to 20 uL with water',
            'Incubate at 50°C for 15 minutes (2-3 fragments)',
            'For 4-6 fragments, incubate for 60 minutes',
          ] : [
            'Add 0.02-0.5 pmol of each fragment (equimolar amounts)',
            'Add Gibson Assembly Master Mix (2X)',
            'Bring total volume to 20 uL with water',
            'Incubate at 50°C for 15-60 minutes',
          ],
        },
        {
          title: 'Transformation',
          details: [
            'Transform 2 uL of assembly reaction into competent cells',
            'Use NEB 5-alpha or equivalent competent cells',
            'Plate on selective media',
            'Incubate overnight at 37°C',
          ],
        },
      ],
      notes: isHiFi ? [
        'NEBuilder HiFi is optimized for fragments 15-35 bp overlaps',
        'Can assemble up to 5+ fragments in a single reaction',
        'Total DNA should not exceed 0.2 pmol for optimal efficiency',
        'Overlapping sequences should have Tm between 48-65°C',
      ] : [
        'Gibson Assembly works best with 15-40 bp overlaps',
        'Tm of overlapping regions should be 48-65°C',
        'Use equimolar amounts of each fragment',
        'For large assemblies, consider increasing incubation time',
      ],
    };
  };

  const handleCopy = useCallback((message: any) => {
    setCopyMessage(message);
    setTimeout(() => setCopyMessage(''), 2000);
  }, []);

  const loadExample = useCallback(() => {
    setParts([
      {
        id: 'J23100_promoter',
        seq: 'TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC',
        type: 'promoter',
      },
      {
        id: 'B0034_rbs',
        seq: 'AAAGAGGAGAAA',
        type: 'rbs',
      },
      {
        id: 'GFP_cds',
        seq: 'ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA',
        type: 'cds',
      },
      {
        id: 'B0015_terminator',
        seq: 'CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA',
        type: 'terminator',
      },
    ]);
    setResult(null);
    setError(null);
  }, []);

  // Fusion Site Optimizer handler - uses the real library functions
  const handleFusionOptimize = useCallback((options: any) => {
    if (!fusionSequence || fusionSequence.length < 200) {
      setError('Please enter a sequence of at least 200bp for fusion site optimization');
      return;
    }

    setIsOptimizing(true);
    setError(null);

    // Run optimization asynchronously to avoid blocking UI
    setTimeout(() => {
      try {
        // Map UI algorithm names to library algorithm names
        const algorithmMap = {
          'auto': 'auto',
          'branch_bound': 'branchBound',
          'dp_validated': 'hybrid',  // Use hybrid for DP-validated
          'monte_carlo': 'monteCarlo',
        };

        // Check for auto-domestication
        const autoDomestication = options.autoDomestication || { enabled: false };
        let domesticationInfo = null;
        let effectiveNumFragments = fusionNumFragments;

        if (autoDomestication.enabled) {
          // Get domestication analysis to determine required junctions
          domesticationInfo = optimizeWithDomestication(
            fusionSequence,
            fusionNumFragments,
            enzyme,
            { constraints: options.constraints }
          );

          if (domesticationInfo.needsDomestication) {
            effectiveNumFragments = domesticationInfo.totalFragments;
            console.log(`Auto-domestication: Adding ${domesticationInfo.additionalFragmentsNeeded} ` +
              `junction(s) to break internal ${enzyme} sites. Total fragments: ${effectiveNumFragments}`);
          }
        }

        // Build optimizer options from UI settings
        const optimizerOptions = {
          enzyme,
          algorithm: (algorithmMap as any)[options.algorithm] || 'auto',
          minFragmentSize: options.constraints?.minFragmentSize || OPTIMIZER_DEFAULTS.minFragmentSize,
          maxFragmentSize: options.constraints?.maxFragmentSize || OPTIMIZER_DEFAULTS.maxFragmentSize,
          minDistanceFromEnds: OPTIMIZER_DEFAULTS.minDistanceFromEnds,
          searchRadius: OPTIMIZER_DEFAULTS.searchRadius,
          // Biological context from UI
          codingFrame: options.bioContext?.isCodingSequence ? options.bioContext.codingFrame : null,
          scarContext: options.bioContext?.scarPreferences || 'nonCoding',
          // Scoring weights from UI
          fidelityWeight: options.weights?.overhangQuality || 0.20,
          efficiencyWeight: 0.15, // Fixed
          primerQualityWeight: (options.weights?.forwardPrimer || 0.20) + (options.weights?.reversePrimer || 0.20),
          positionWeight: 0.15, // Fixed
          verbose: false,
          // Pass domestication info if available
          domesticationJunctions: domesticationInfo?.domesticationJunctions || [],
        };

        // First, scan for all candidates to show in the UI
        const scanOptions = {
          enzyme,
          minDistanceFromEnds: 50,
          includeEfficiency: true,
          maxCandidates: 500,
        };

        const scannedCandidates = scanForFusionSites(fusionSequence, scanOptions);

        // Score each candidate with full composite scoring
        const scoredCandidates = scannedCandidates.map(candidate => {
          try {
            const scored = scoreFusionSiteComposite(fusionSequence, candidate.position, enzyme, {
              codingFrame: options.bioContext?.isCodingSequence ? options.bioContext.codingFrame : null,
              scarContext: options.bioContext?.scarPreferences || 'nonCoding',
            });
            return {
              ...candidate,
              composite: scored.composite * 100, // Convert to 0-100 scale
              scores: {
                overhangQuality: { score: Math.round(((scored as any).breakdown?.overhangQuality || 0.7) * 100) },
                forwardPrimer: { score: Math.round(((scored as any).breakdown?.forwardPrimer || 0.7) * 100) },
                reversePrimer: { score: Math.round(((scored as any).breakdown?.reversePrimer || 0.7) * 100) },
                riskFactors: { score: Math.round(((scored as any).breakdown?.riskFactors || 0.7) * 100) },
                biologicalContext: { score: Math.round(((scored as any).breakdown?.biologicalContext || 0.7) * 100) },
              },
              warnings: scored.warnings || [],
              isTNNA: candidate.overhang.startsWith('T') && candidate.overhang.endsWith('A'),
              isHighGC: (candidate.overhang.match(/[GC]/g) || []).length >= 3,
            };
          } catch (e: unknown) {
            // Fallback if scoring fails
            return {
              ...candidate,
              composite: candidate.efficiency ? candidate.efficiency * 100 : 70,
              scores: {
                overhangQuality: { score: 70 },
                forwardPrimer: { score: 70 },
                reversePrimer: { score: 70 },
                riskFactors: { score: 70 },
              },
              warnings: [],
            };
          }
        });

        // Sort by composite score
        scoredCandidates.sort((a, b) => b.composite - a.composite);
        setFusionCandidates(scoredCandidates as any);

        // Run the main optimization with effective fragment count (includes domestication junctions)
        const optimizationResult = optimizeFusionSites(fusionSequence, effectiveNumFragments, optimizerOptions);

        if (!optimizationResult.success) {
          setError(optimizationResult.error || 'Optimization failed - no valid junction set found');
          setFusionResult(null);
          setIsOptimizing(false);
          return;
        }

        // Format the junctions for UI display
        const formattedJunctions = (optimizationResult.junctions || []).map((j: any, i: number) => {
          const detailed = optimizationResult.detailedJunctions?.[i];
          return {
            position: j.position,
            overhang: j.overhang,
            composite: (detailed?.composite || j.score?.composite || 0.75) * 100,
            scores: (detailed as any)?.breakdown ? {
              overhangQuality: { score: Math.round(((detailed as any).breakdown.overhangQuality || 0.7) * 100) },
              forwardPrimer: { score: Math.round(((detailed as any).breakdown.forwardPrimer || 0.7) * 100) },
              reversePrimer: { score: Math.round(((detailed as any).breakdown.reversePrimer || 0.7) * 100) },
              riskFactors: { score: Math.round(((detailed as any).breakdown.riskFactors || 0.7) * 100) },
            } : (j as any).scores,
            warnings: detailed?.warnings || (j as any).warnings || [],
          };
        });

        // Format failure prediction for UI
        const failurePred = optimizationResult.failurePrediction || {};
        const successRate = failurePred.expectedSuccessRate || optimizationResult.score?.fidelity || 0.85;

        const formattedFailurePrediction = {
          summary: {
            predictedSuccessRate: successRate,
            predictedSuccessPercent: `${(successRate * 100).toFixed(0)}%`,
            highRiskCount: (failurePred.predictions || []).filter((p: any) => p.severity === 'high').length,
            mediumRiskCount: (failurePred.predictions || []).filter((p: any) => p.severity === 'medium').length,
            lowRiskCount: (failurePred.predictions || []).filter((p: any) => p.severity === 'low').length,
            recommendation: failurePred.recommendations?.[0] ||
              (successRate >= 0.9
                ? 'Assembly design meets quality thresholds. Proceed with confidence.'
                : successRate >= 0.7
                ? 'Consider optimizing low-scoring junctions for better reliability.'
                : 'Multiple high-risk junctions detected. Consider revising the design.'),
          },
          predictions: (failurePred.predictions || []).map((p: any) => ({
            type: p.type || 'risk',
            severity: p.severity || 'medium',
            fragment: p.fragment || '',
            probability: p.probability || 0.2,
            message: p.message || p.description || 'Potential issue detected',
            mitigation: p.mitigation || p.recommendation || 'Review junction design',
          })),
        };

        // Set the result
        setFusionResult({
          algorithm: optimizationResult.algorithm,
          optimal: optimizationResult.algorithm === 'branchBound',
          nodesExplored: (optimizationResult as any).nodesExplored || 0,
          solution: {
            junctions: formattedJunctions,
            overhangs: optimizationResult.overhangs,
            setFidelity: optimizationResult.score?.fidelity || successRate,
          },
          failurePrediction: formattedFailurePrediction,
          // Include additional info for display
          fragmentSizes: optimizationResult.fragmentSizes,
          internalSites: optimizationResult.internalSites,
          summary: optimizationResult.summary,
          // Auto-domestication info
          autoDomestication: domesticationInfo ? {
            enabled: true,
            needsDomestication: domesticationInfo.needsDomestication,
            internalSites: domesticationInfo.internalSites || [],
            domesticationJunctions: domesticationInfo.domesticationJunctions || [],
            additionalFragments: domesticationInfo.additionalFragmentsNeeded || 0,
            message: domesticationInfo.message,
          } : { enabled: false },
        });

        setCopyMessage('Optimization complete');
      } catch (err: unknown) {
        console.error('Fusion optimization error:', err);
        setError('Optimization failed: ' + (err as Error).message);
        setFusionResult(null);
      } finally {
        setIsOptimizing(false);
      }
    }, 100); // Small delay to allow UI to update
  }, [fusionSequence, fusionNumFragments, enzyme]);

  // Handler to accept fusion optimization results and design primers
  const handleAcceptAndDesignPrimers = useCallback(() => {
    if (fusionResult?.success && fusionResult?.fragments) {
      // The fusion results are already in the state, just trigger design
      designAssembly();
    }
  }, [fusionResult, designAssembly]);

  // Export handlers
  const handleExportGenBank = useCallback(() => {
    if (!result || !result.assembledSequence) {
      setError('No assembly result to export. Design primers first.');
      return;
    }

    try {
      const genbank = exportToGenBank(
        result.simulation || { sequence: result.assembledSequence, features: result.parts },
        {
          name: `Assembly_${method}`,
          description: `DNA assembly designed using ${(ASSEMBLY_METHODS as any)[method.toUpperCase()]?.name || method}`,
          molType: 'DNA',
          topology: 'circular',
        } as any
      );

      const blob = new Blob([genbank], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `assembly_${method}_${Date.now()}.gb`;
      a.click();
      URL.revokeObjectURL(url);
      setCopyMessage('GenBank file downloaded');
    } catch (err: unknown) {
      setError('Failed to export GenBank: ' + (err as Error).message);
    }
  }, [result, method]);

  const handleExportFasta = useCallback(() => {
    if (!result || !result.assembledSequence) {
      setError('No assembly result to export. Design primers first.');
      return;
    }

    try {
      const fasta = exportToFasta(
        result.simulation || { sequence: result.assembledSequence, features: result.parts },
        { name: `Assembly_${method}` }
      );

      const blob = new Blob([fasta], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `assembly_${method}_${Date.now()}.fasta`;
      a.click();
      URL.revokeObjectURL(url);
      setCopyMessage('FASTA file downloaded');
    } catch (err: unknown) {
      setError('Failed to export FASTA: ' + (err as Error).message);
    }
  }, [result, method]);

  const handleSaveProject = useCallback(() => {
    const projectData = {
      method,
      parts: parts.map(p => ({
        id: p.id,
        seq: p.seq,
        type: p.type,
        rangeStart: p.rangeStart,
        rangeEnd: p.rangeEnd,
      })),
      enzyme: isGoldenGate ? enzyme : null,
      overlapSettings: !isGoldenGate ? overlapSettings : null,
      result,
    };

    try {
      const content = exportProject(projectData);
      const blob = new Blob([content], { type: 'application/json' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `assembly_project_${Date.now()}.json`;
      a.click();
      URL.revokeObjectURL(url);
      setCopyMessage('Project saved');
    } catch (err: unknown) {
      setError('Failed to save project: ' + (err as Error).message);
    }
  }, [method, parts, enzyme, overlapSettings, result, isGoldenGate]);

  const handleLoadProject = useCallback((e: any) => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event: ProgressEvent<FileReader>) => {
      try {
        const project = importProject(event.target?.result as string);

        // Restore project state
        if (project.method) setMethod(project.method);
        if (project.parts) setParts(project.parts);
        if (project.enzyme) setEnzyme(project.enzyme);
        if (project.overlapSettings) setOverlapSettings(project.overlapSettings);
        if (project.result) setResult(project.result);

        setCopyMessage('Project loaded');
        setError(null);
      } catch (err: unknown) {
        setError('Failed to load project: ' + (err as Error).message);
      }
    };
    reader.readAsText(file);
    // Reset file input
    e.target.value = '';
  }, []);

  const validPartsCount = parts.filter(p => p.seq && p.seq.length > 10).length;

  return (
    <div className="gg-designer-v2">
      {copyMessage && <div className="toast-notification">{copyMessage}</div>}

      {/* Method Selector Cards */}
      <div className="method-selector-pro">
        <div className="method-cards">
          {Object.entries(ASSEMBLY_METHODS).map(([key, m]) => (
            <button
              key={m.id}
              className={`method-card ${method === m.id ? 'active' : ''}`}
              onClick={() => { setMethod(m.id); setResult(null); }}
              style={{ '--method-color': m.color } as React.CSSProperties}
            >
              <div className="method-card-icon">
                {m.icon}
              </div>
              <div className="method-card-content">
                <span className="method-card-name">{m.name}</span>
                <span className="method-card-desc">{m.description}</span>
              </div>
              {method === m.id && (
                <div className="method-card-check">
                  <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
                    <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
                  </svg>
                </div>
              )}
            </button>
          ))}
        </div>
      </div>

      {/* Advanced Options (Golden Gate) - Collapsible */}
      {isGoldenGate && (
        <div className="advanced-options-section">
          <button
            className="advanced-options-toggle"
            onClick={() => setShowAdvancedOptions(!showAdvancedOptions)}
          >
            <span className="advanced-options-toggle-label">
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                <path d="M19.14 12.94c.04-.31.06-.63.06-.94 0-.31-.02-.63-.06-.94l2.03-1.58c.18-.14.23-.41.12-.61l-1.92-3.32c-.12-.22-.37-.29-.59-.22l-2.39.96c-.5-.38-1.03-.7-1.62-.94l-.36-2.54c-.04-.24-.24-.41-.48-.41h-3.84c-.24 0-.43.17-.47.41l-.36 2.54c-.59.24-1.13.57-1.62.94l-2.39-.96c-.22-.08-.47 0-.59.22L2.74 8.87c-.12.21-.08.47.12.61l2.03 1.58c-.04.31-.06.63-.06.94s.02.63.06.94l-2.03 1.58c-.18.14-.23.41-.12.61l1.92 3.32c.12.22.37.29.59.22l2.39-.96c.5.38 1.03.7 1.62.94l.36 2.54c.05.24.24.41.48.41h3.84c.24 0 .44-.17.47-.41l.36-2.54c.59-.24 1.13-.56 1.62-.94l2.39.96c.22.08.47 0 .59-.22l1.92-3.32c.12-.22.07-.47-.12-.61l-2.01-1.58zM12 15.6c-1.98 0-3.6-1.62-3.6-3.6s1.62-3.6 3.6-3.6 3.6 1.62 3.6 3.6-1.62 3.6-3.6 3.6z"/>
              </svg>
              Advanced Options
            </span>
            <svg
              className={`advanced-options-chevron ${showAdvancedOptions ? 'expanded' : ''}`}
              viewBox="0 0 24 24"
              width="16"
              height="16"
              fill="currentColor"
            >
              <path d="M7.41 8.59L12 13.17l4.59-4.58L18 10l-6 6-6-6 1.41-1.41z"/>
            </svg>
          </button>

          {showAdvancedOptions && (
            <div className="advanced-options-content">
              {/* Auto-Domestication Toggle */}
              <div className="fusion-mode-toggle">
                <div className="toggle-label">
                  <span>Auto-Domestication</span>
                  <span>
                    Automatically detect and break internal {enzyme} sites by adding junction positions
                  </span>
                </div>
                <div
                  className={`fusion-mode-switch ${autoDomesticationEnabled ? 'active' : ''}`}
                  onClick={() => {
                    setAutoDomesticationEnabled(!autoDomesticationEnabled);
                    setDomesticationError(null);
                  }}
                />
              </div>

              {/* Min Fragment Size (shown when auto-domestication is enabled) */}
              {autoDomesticationEnabled && (
                <div className="advanced-option-row">
                  <div className="option-label">
                    <span>Minimum Fragment Size</span>
                    <span className="option-hint">
                      Minimum bp for fragments created by domestication junctions
                    </span>
                  </div>
                  <div className="option-input">
                    <input
                      type="number"
                      min="20"
                      max="500"
                      step="10"
                      value={minFragmentSize}
                      onChange={(e: React.ChangeEvent<HTMLInputElement>) => {
                        const val = Math.max(20, Math.min(500, parseInt(e.target.value) || 50));
                        setMinFragmentSize(val);
                        setDomesticationError(null);
                      }}
                      className="fragment-size-input"
                    />
                    <span className="input-unit">bp</span>
                  </div>
                </div>
              )}

              {/* Enhanced Domestication Workflow Toggle */}
              {autoDomesticationEnabled && (
                <div className="fusion-mode-toggle enhanced-workflow-toggle">
                  <div className="toggle-label">
                    <span>🧬 Enhanced Workflow (Recommended)</span>
                    <span>
                      Frame validation, protein preview & user confirmation before mutations
                    </span>
                  </div>
                  <div
                    className={`fusion-mode-switch ${useEnhancedWorkflow ? 'active' : ''}`}
                    onClick={() => setUseEnhancedWorkflow(!useEnhancedWorkflow)}
                  />
                </div>
              )}

              {/* Fusion Site Optimizer Toggle */}
              <div className="fusion-mode-toggle">
                <div className="toggle-label">
                  <span>Fusion Site Optimizer</span>
                  <span>
                    Find optimal junction positions within a single sequence (splits a gene into fragments)
                  </span>
                </div>
                <div
                  className={`fusion-mode-switch ${fusionMode ? 'active' : ''}`}
                  onClick={() => {
                    setFusionMode(!fusionMode);
                    if (!fusionMode) {
                      setResult(null);
                      setError(null);
                    } else {
                      setFusionResult(null);
                      setFusionCandidates([]);
                    }
                  }}
                />
              </div>
            </div>
          )}
        </div>
      )}

      {/* Advanced Options (Isothermal Assembly) - Collapsible */}
      {isIsothermal && (
        <div className="advanced-options-section isothermal-advanced">
          <button
            className="advanced-options-toggle"
            onClick={() => setShowIsothermalAdvanced(!showIsothermalAdvanced)}
          >
            <span className="advanced-options-toggle-label">
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                <path d="M19.14 12.94c.04-.31.06-.63.06-.94 0-.31-.02-.63-.06-.94l2.03-1.58c.18-.14.23-.41.12-.61l-1.92-3.32c-.12-.22-.37-.29-.59-.22l-2.39.96c-.5-.38-1.03-.7-1.62-.94l-.36-2.54c-.04-.24-.24-.41-.48-.41h-3.84c-.24 0-.43.17-.47.41l-.36 2.54c-.59.24-1.13.57-1.62.94l-2.39-.96c-.22-.08-.47 0-.59.22L2.74 8.87c-.12.21-.08.47.12.61l2.03 1.58c-.04.31-.06.63-.06.94s.02.63.06.94l-2.03 1.58c-.18.14-.23.41-.12.61l1.92 3.32c.12.22.37.29.59.22l2.39-.96c.5.38 1.03.7 1.62.94l.36 2.54c.05.24.24.41.48.41h3.84c.24 0 .44-.17.47-.41l.36-2.54c.59-.24 1.13-.56 1.62-.94l2.39.96c.22.08.47 0 .59-.22l1.92-3.32c.12-.22.07-.47-.12-.61l-2.01-1.58zM12 15.6c-1.98 0-3.6-1.62-3.6-3.6s1.62-3.6 3.6-3.6 3.6 1.62 3.6 3.6-1.62 3.6-3.6 3.6z"/>
              </svg>
              Advanced Options
              <span className="current-settings-badge">
                {isothermalVariant === 'nebuilder_hifi' ? 'NEBuilder HiFi' : 'Gibson'} • {overlapSettings.overlapLength}bp • {overlapSettings.targetTm}°C
              </span>
            </span>
            <svg
              className={`advanced-options-chevron ${showIsothermalAdvanced ? 'expanded' : ''}`}
              viewBox="0 0 24 24"
              width="16"
              height="16"
              fill="currentColor"
            >
              <path d="M7.41 8.59L12 13.17l4.59-4.58L18 10l-6 6-6-6 1.41-1.41z"/>
            </svg>
          </button>

          {showIsothermalAdvanced && (
            <div className="advanced-options-content isothermal-options-content">
              {/* Protocol Variant */}
              <div className="isothermal-option-group">
                <div className="option-group-header">
                  <h4>Protocol Variant</h4>
                  <span className="option-group-desc">Choose the assembly kit/protocol</span>
                </div>
                <div className="variant-options-compact">
                  {Object.entries(ISOTHERMAL_VARIANTS).map(([key, v]) => (
                    <button
                      key={v.id}
                      className={`variant-option-compact ${isothermalVariant === v.id ? 'active' : ''}`}
                      onClick={() => {
                        setIsothermalVariant(v.id);
                        // Update overlap settings to optimal for this variant
                        setOverlapSettings({
                          overlapLength: v.optimalOverlap.ideal,
                          targetTm: v.optimalTm.ideal,
                        });
                      }}
                    >
                      <div className="variant-radio-compact"></div>
                      <div className="variant-info-compact">
                        <strong>{v.name}</strong>
                        <span className="variant-badge-compact">{v.badge}</span>
                        {(v as any).recommended && <span className="recommended-tag-compact">Recommended</span>}
                      </div>
                    </button>
                  ))}
                </div>
              </div>

              {/* Overlap Settings */}
              <div className="isothermal-option-group">
                <div className="option-group-header">
                  <h4>Overlap Parameters</h4>
                  <span className="option-group-desc">
                    {isothermalVariant === 'nebuilder_hifi'
                      ? 'NEB recommends 15-20bp for 2-3 fragments, 20-30bp for 4-6 fragments'
                      : 'Recommended 20-40bp overlaps with Tm 50-60°C'}
                  </span>
                </div>
                <OverlapSettings
                  settings={overlapSettings}
                  onChange={setOverlapSettings}
                  variant={isothermalVariant}
                />
              </div>

              {/* Advanced Optimization Mode Toggle */}
              <div className="isothermal-option-group">
                <div className="option-group-header">
                  <h4>Optimization Mode</h4>
                  <span className="option-group-desc">
                    Choose between standard and advanced primer optimization
                  </span>
                </div>
                <div className="optimization-mode-toggle">
                  <label className={`optimization-option ${!useAdvancedOptimization ? 'active' : ''}`}>
                    <input
                      type="radio"
                      name="optimization-mode"
                      checked={!useAdvancedOptimization}
                      onChange={() => setUseAdvancedOptimization(false)}
                    />
                    <div className="optimization-info">
                      <strong>Standard</strong>
                      <span>Target-based optimization with fixed Tm/length</span>
                    </div>
                  </label>
                  <label className={`optimization-option ${useAdvancedOptimization ? 'active' : ''}`}>
                    <input
                      type="radio"
                      name="optimization-mode"
                      checked={useAdvancedOptimization}
                      onChange={() => setUseAdvancedOptimization(true)}
                    />
                    <div className="optimization-info">
                      <strong>Advanced (Recommended)</strong>
                      <span>Range-based scoring with position flexibility, cross-junction analysis</span>
                    </div>
                    <span className="new-badge">NEW</span>
                  </label>
                </div>
              </div>
            </div>
          )}
        </div>
      )}

      {/* Domestication Error Banner with Enzyme Recommendations */}
      {domesticationError && isGoldenGate && autoDomesticationEnabled && (
        <div className="domestication-error-banner">
          <div className="error-header">
            <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor" className="error-icon">
              <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z"/>
            </svg>
            <span className="error-title">Auto-Domestication Error</span>
          </div>
          <p className="error-message">{domesticationError.message}</p>

          {(domesticationError.error as any)?.type === 'ADJACENT_SITES_TOO_CLOSE' && (
            <div className="error-details">
              <p className="detail-label">Sites too close together:</p>
              <ul className="site-list">
                {(domesticationError.error as any).details?.map((detail: any, i: number) => (
                  <li key={i}>
                    Positions {detail.site1Position + 1} and {detail.site2Position + 1} are only {detail.distance}bp apart
                    (minimum: {detail.requiredDistance}bp)
                  </li>
                ))}
              </ul>
            </div>
          )}

          {(domesticationError.error as any)?.type === 'FRAGMENT_SIZE_VIOLATION' && (
            <div className="error-details">
              <p className="detail-label">Fragment size violations:</p>
              <ul className="site-list">
                {(domesticationError.error as any).violations?.map((v: any, i: number) => (
                  <li key={i}>{v.message}</li>
                ))}
              </ul>
            </div>
          )}

          {domesticationError.alternativeEnzymes && domesticationError.alternativeEnzymes.length > 0 && (
            <div className="enzyme-recommendations">
              <p className="recommendation-label">Alternative Enzymes (click to switch):</p>
              <div className="enzyme-buttons">
                {domesticationError.alternativeEnzymes
                  .filter(alt => !alt.isCurrent) // Don't show current enzyme as alternative
                  .slice(0, 4)
                  .map((alt: any) => (
                  <button
                    key={alt.enzyme}
                    className={`enzyme-recommendation-btn ${alt.isCompatible ? 'compatible' : 'partial'}`}
                    onClick={() => {
                      setEnzyme(alt.enzyme);
                      setDomesticationError(null);
                      setResult(null);
                    }}
                    title={`${alt.fullName}: ${alt.internalSites} internal site${alt.internalSites !== 1 ? 's' : ''}`}
                  >
                    <span className="enzyme-name">{alt.enzyme}</span>
                    <span className="enzyme-status">
                      {alt.isCompatible ? (
                        <span className="compatible-badge">No sites</span>
                      ) : (
                        <span className="site-count">{alt.internalSites} site{alt.internalSites !== 1 ? 's' : ''}</span>
                      )}
                    </span>
                    {alt.enzyme === domesticationError.recommendedEnzyme && (
                      <span className="recommended-badge">Recommended</span>
                    )}
                  </button>
                ))}
              </div>
            </div>
          )}

          <button
            className="dismiss-error-btn"
            onClick={() => setDomesticationError(null)}
          >
            Dismiss
          </button>
        </div>
      )}

      {/* Fusion Mode Content */}
      {fusionMode && isGoldenGate ? (
        <div className="fusion-mode-content">
          {/* Sequence Input for Fusion Mode */}
          <div className="setup-card fusion-sequence-card">
            <div className="card-header">
              <h3>
                <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
                  <path d="M14 2H6c-1.1 0-1.99.9-1.99 2L4 20c0 1.1.89 2 1.99 2H18c1.1 0 2-.9 2-2V8l-6-6zm4 18H6V4h7v5h5v11z"/>
                </svg>
                Input Sequence
              </h3>
              <div className="card-actions">
                <span className="seq-length-indicator">
                  {fusionSequence.length.toLocaleString()} bp
                </span>
              </div>
            </div>
            <div className="card-body">
              <textarea
                className="fusion-sequence-input"
                placeholder="Paste your DNA sequence here (minimum 500bp). The optimizer will find optimal junction positions to split this into fragments."
                value={fusionSequence}
                onChange={(e: React.ChangeEvent<HTMLTextAreaElement>) => {
                  const seq = e.target.value.replace(/[^ACGTacgt\n\r]/g, '').toUpperCase();
                  setFusionSequence(seq);
                  setFusionResult(null);
                  setFusionCandidates([]);
                }}
                rows={6}
              />
              <div className="fusion-input-options">
                <div className="option-row">
                  <label>Number of Fragments</label>
                  <div className="fragment-selector">
                    {[2, 3, 4, 5, 6, 8, 10].map(n => (
                      <button
                        key={n}
                        className={`fragment-btn ${fusionNumFragments === n ? 'active' : ''}`}
                        onClick={() => {
                          setFusionNumFragments(n);
                          setFusionResult(null);
                          setFusionCandidates([]);
                        }}
                      >
                        {n}
                      </button>
                    ))}
                  </div>
                  <span className="option-hint">{fusionNumFragments - 1} junctions needed</span>
                </div>
              </div>
            </div>
          </div>

          {/* Fusion Optimizer Panel */}
          <FusionSiteOptimizerPanel
            sequence={fusionSequence}
            enzyme={enzyme}
            numFragments={fusionNumFragments}
            onOptimize={handleFusionOptimize}
            result={fusionResult}
            candidates={fusionCandidates as any}
            isOptimizing={isOptimizing}
            globalAutoDomestication={autoDomesticationEnabled}
            onAutoDomesticationChange={setAutoDomesticationEnabled}
          />

          {/* Accept & Design Primers Action Bar - shown when optimization succeeds */}
          {fusionResult?.success && (
            <div className="fusion-action-bar mt-4 p-4 rounded-lg border border-green-300 flex items-center justify-between gap-4" style={{
              background: 'linear-gradient(135deg, #f0fdf4 0%, #dcfce7 100%)',
            }}>
              <div className="fusion-summary flex gap-4 items-center">
                <span className="py-1 px-3 bg-green-500 text-white rounded-full text-sm font-semibold">
                  {((fusionResult.solution?.setFidelity || fusionResult._fullResult?.score?.fidelity || 0.9) * 100).toFixed(1)}% Fidelity
                </span>
                <span className="text-green-800 text-sm">
                  {fusionResult.numFragments || fusionNumFragments} Fragments
                </span>
                <span className={`py-1 px-2 rounded text-xs ${
                  fusionResult.failurePrediction?.summary?.highRiskCount > 0
                    ? 'bg-red-50 text-red-600'
                    : 'bg-green-50 text-green-800'
                }`}>
                  {fusionResult.failurePrediction?.summary?.highRiskCount > 0
                    ? `${fusionResult.failurePrediction.summary.highRiskCount} high risk`
                    : 'Low risk'}
                </span>
              </div>

              <button
                onClick={handleAcceptAndDesignPrimers}
                disabled={isOptimizing}
                className={`flex items-center gap-2 py-3 px-6 bg-green-500 text-white border-none rounded-lg text-base font-semibold transition-all duration-200 hover:bg-green-600 ${
                  isOptimizing ? 'cursor-not-allowed opacity-70' : 'cursor-pointer'
                }`}
              >
                <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
                  <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
                </svg>
                {isOptimizing ? 'Designing...' : 'Accept & Design Primers'}
              </button>
            </div>
          )}

          {error && (
            <div className="error-message-box">
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
              </svg>
              <span>{error}</span>
            </div>
          )}
        </div>
      ) : (
        <>
          {/* Export/Import Toolbar */}
          <div className="assembly-toolbar">
            <div className="toolbar-group">
              <span className="toolbar-label">Project</span>
              <button className="toolbar-btn" onClick={handleSaveProject} title="Save project">
                <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                  <path d="M17 3H5c-1.11 0-2 .9-2 2v14c0 1.1.89 2 2 2h14c1.1 0 2-.9 2-2V7l-4-4zm-5 16c-1.66 0-3-1.34-3-3s1.34-3 3-3 3 1.34 3 3-1.34 3-3 3zm3-10H5V5h10v4z"/>
                </svg>
                Save
              </button>
              <label className="toolbar-btn" title="Load project">
                <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                  <path d="M9 16h6v-6h4l-7-7-7 7h4zm-4 2h14v2H5z"/>
                </svg>
                Load
                <input type="file" accept=".json" onChange={handleLoadProject} className="hidden" />
              </label>
            </div>
            <div className="toolbar-divider" />
            <div className="toolbar-group">
              <span className="toolbar-label">Export</span>
              <button
                className="toolbar-btn"
                onClick={handleExportGenBank}
                disabled={!result}
                title="Export as GenBank"
              >
                <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                  <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                </svg>
                GenBank
              </button>
              <button
                className="toolbar-btn"
                onClick={handleExportFasta}
                disabled={!result}
                title="Export as FASTA"
              >
                <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                  <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                </svg>
                FASTA
              </button>
            </div>
          </div>

          <div className="gg-main-layout">
        {/* Setup Column */}
        <div className="setup-column">
          {/* Card 1: Enzyme Selection (Golden Gate only) */}
          {isGoldenGate && (
            <div className="setup-card">
              <div className="card-header">
                <h3>
                  <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
                    <path d="M7 14l5-5 5 5z"/>
                  </svg>
                  Type IIS Enzyme
                </h3>
              </div>
              <div className="card-body">
                <div className="enzyme-chips">
                  {Object.entries(GOLDEN_GATE_ENZYMES).map(([name, enz]) => {
                    const hasData = ENZYMES_WITH_DATA.includes(name);
                    return (
                      <button
                        key={name}
                        className={`enzyme-chip ${enzyme === name ? 'selected' : ''}`}
                        onClick={() => { setEnzyme(name); setResult(null); setDomesticationError(null); }}
                      >
                        <span className="chip-name">{name}</span>
                        <span className="chip-site">{enz.recognition}</span>
                        {hasData && <span className="chip-badge">✓</span>}
                      </button>
                    );
                  })}
                </div>
              </div>
            </div>
          )}

          {/* Card 2: Parts */}
          <div className="setup-card">
            <div className="card-header">
              <h3>
                <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
                  <path d="M19 13h-6v6h-2v-6H5v-2h6V5h2v6h6v2z"/>
                </svg>
                {isGoldenGate ? 'Insert Sequences' : 'DNA Fragments'}
                <span className="count-badge">{parts.length}/{maxParts}</span>
              </h3>
              <div className="card-actions">
                <button className="text-btn" onClick={loadExample}>Load Example</button>
                <label className="text-btn file-label">
                  Upload
                  <input type="file" accept=".fa,.fasta,.txt" onChange={handleFileUpload} />
                </label>
              </div>
            </div>
            <div className="card-body">
              <div className="parts-list-compact">
                {parts.map((part: Part, i: number) => (
                  <PartCard
                    key={i}
                    part={part}
                    index={i}
                    onChange={handlePartChange}
                    onRemove={removePart}
                    canRemove={true}
                    enzyme={isGoldenGate ? enzyme : null as any}
                    onDragStart={handleDragStart}
                    onDragOver={handleDragOver}
                    onDrop={handleDrop}
                    isDragging={draggedIndex === i}
                    onSplitSequence={isGoldenGate ? handleSplitSequenceFromPart : null as any}
                    useEnhancedWorkflow={useEnhancedWorkflow}
                    onOpenEnhancedDomestication={(partData) => {
                      setEnhancedDomesticationPart(partData);
                      setShowEnhancedDomestication(true);
                    }}
                    onOpenWorkflowGuide={(partData) => {
                      setWorkflowGuidePart(partData);
                      setShowWorkflowGuide(true);
                    }}
                  />
                ))}
              </div>

              {parts.length < maxParts && (
                <button className="add-part-btn-compact" onClick={addPart}>
                  <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                    <path d="M19 13h-6v6h-2v-6H5v-2h6V5h2v6h6v2z"/>
                  </svg>
                  Add {isGoldenGate ? 'Part' : 'Fragment'}
                </button>
              )}
            </div>
          </div>

          {/* Card 3: Assembly Plan */}
          {validPartsCount >= 2 && (
            <div className="setup-card">
              <div className="card-header">
                <h3>
                  <svg viewBox="0 0 24 24" width="18" height="18" fill="currentColor">
                    <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-1 17.93c-3.95-.49-7-3.85-7-7.93 0-.62.08-1.21.21-1.79L9 15v1c0 1.1.9 2 2 2v1.93zm6.9-2.54c-.26-.81-1-1.39-1.9-1.39h-1v-3c0-.55-.45-1-1-1H8v-2h2c.55 0 1-.45 1-1V7h2c1.1 0 2-.9 2-2v-.41c2.93 1.19 5 4.06 5 7.41 0 2.08-.8 3.97-2.1 5.39z"/>
                  </svg>
                  Assembly Plan
                </h3>
                <div className="view-mode-toggle">
                  <button
                    className={`toggle-btn ${viewMode === 'linear' ? 'active' : ''}`}
                    onClick={() => setViewMode('linear')}
                  >
                    Linear
                  </button>
                  <button
                    className={`toggle-btn ${viewMode === 'circular' ? 'active' : ''}`}
                    onClick={() => setViewMode('circular')}
                  >
                    Circular
                  </button>
                </div>
              </div>
              <div className="card-body">
                {viewMode === 'circular' ? (
                  <CircularPlasmidView
                    parts={parts.filter(p => p.seq)}
                    overhangs={isGoldenGate ? recommendedOverhangs.overhangs : parts.filter(p => p.seq).map((_, i) => `J${i + 1}`)}
                    totalLength={totalLength}
                  />
                ) : (
                  <LinearAssemblyDiagram
                    parts={parts.filter(p => p.seq)}
                    overhangs={isGoldenGate ? recommendedOverhangs.overhangs : parts.filter(p => p.seq).map((_, i) => `${overlapSettings.overlapLength}bp`)}
                  />
                )}

                {isGoldenGate ? (
                  <FidelityGauge
                    overhangs={recommendedOverhangs.overhangs}
                    enzyme={enzyme}
                    staticFidelity={recommendedOverhangs.fidelity}
                  />
                ) : (
                  <OverlapQualityGauge
                    overlapLength={overlapSettings.overlapLength}
                    targetTm={overlapSettings.targetTm}
                    fragmentCount={validPartsCount}
                    variant={isothermalVariant}
                  />
                )}
              </div>
            </div>
          )}

          {/* Design Button */}
          <button
            className="design-btn-primary"
            onClick={designAssembly}
            disabled={validPartsCount < 2}
          >
            <svg viewBox="0 0 24 24" width="20" height="20" fill="currentColor">
              <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
            </svg>
            Design Primers
          </button>

          {error && (
            <div className="error-message-box">
              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z"/>
              </svg>
              <span>{error}</span>
            </div>
          )}
        </div>

        {/* Results Column */}
        <div className="results-column">
          {!result ? (
            <div className="empty-results-state">
              <div className="empty-icon">
                {isGoldenGate ? (
                  <svg viewBox="0 0 24 24" width="64" height="64" fill="currentColor" opacity="0.2">
                    <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-1 17.93c-3.95-.49-7-3.85-7-7.93 0-.62.08-1.21.21-1.79L9 15v1c0 1.1.9 2 2 2v1.93zm6.9-2.54c-.26-.81-1-1.39-1.9-1.39h-1v-3c0-.55-.45-1-1-1H8v-2h2c.55 0 1-.45 1-1V7h2c1.1 0 2-.9 2-2v-.41c2.93 1.19 5 4.06 5 7.41 0 2.08-.8 3.97-2.1 5.39z"/>
                  </svg>
                ) : (
                  <svg viewBox="0 0 24 24" width="64" height="64" fill="currentColor" opacity="0.2">
                    <path d="M4 6h16v2H4zm4 5h12v2H8zm-4 5h16v2H4z"/>
                    <circle cx="4" cy="7" r="2" />
                    <circle cx="8" cy="12" r="2" />
                    <circle cx="4" cy="17" r="2" />
                  </svg>
                )}
              </div>
              <h3>
                {isGoldenGate ? 'Golden Gate Assembly' :
                 isothermalVariant === 'nebuilder_hifi' ? 'NEBuilder HiFi Assembly' : 'Gibson Assembly'}
              </h3>
              <p>
                {isGoldenGate
                  ? 'Enter your insert sequences to design PCR primers with high-fidelity overhangs.'
                  : isothermalVariant === 'nebuilder_hifi'
                    ? 'Design primers with optimized homology regions for high-fidelity isothermal assembly.'
                    : 'Enter your DNA fragments to design primers with overlapping homology tails.'}
              </p>
              <div className="quick-steps">
                {isGoldenGate ? (
                  <>
                    <div className="step"><span className="step-num">1</span> Add DNA sequences</div>
                    <div className="step"><span className="step-num">2</span> Select enzyme</div>
                    <div className="step"><span className="step-num">3</span> Get optimized primers</div>
                  </>
                ) : (
                  <>
                    <div className="step"><span className="step-num">1</span> Add 2-6 DNA fragments</div>
                    <div className="step"><span className="step-num">2</span> Set overlap length (20-25bp)</div>
                    <div className="step"><span className="step-num">3</span> Target Tm ~55°C for overlaps</div>
                    <div className="step"><span className="step-num">4</span> Get primers with homology tails</div>
                  </>
                )}
              </div>
              {!isGoldenGate && (
                <div className="assembly-tips">
                  <div className="tip-header">
                    <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                      <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-6h2v6zm0-8h-2V7h2v2z"/>
                    </svg>
                    <span>Tips for optimal results</span>
                  </div>
                  <ul className="tip-list">
                    <li>Use fragments &gt;200bp for reliable PCR amplification</li>
                    <li>Keep overlap GC content between 40-60%</li>
                    <li>Avoid secondary structures in overlap regions</li>
                  </ul>
                </div>
              )}
            </div>
          ) : (
            <>
              {/* Advanced Isothermal Assembly Panel - state-of-the-art optimization */}
              {isIsothermal && useAdvancedOptimization ? (
                <IsothermalAssemblyPanel
                  fragments={parts.filter(p => p.seq && p.seq.length > 0).map((p: Part, i: number) => ({
                    id: p.id || `Fragment_${i + 1}`,
                    seq: p.seq,
                    type: p.type,
                  })) as any}
                  variant={isothermalVariant}
                  config={{
                    targetOverlap: overlapSettings.overlapLength,
                    targetTm: overlapSettings.targetTm,
                    circular: true,
                  }}
                  onPrimersCopied={handleCopy}
                />
              ) : (
                <>
                  <PrimerResults result={result} onCopy={handleCopy} method={method} enzyme={enzyme} />
                  {!isGoldenGate && result.overlapAnalysis && (
                    <OverlapAnalysisDisplay overlaps={result.overlapAnalysis} />
                  )}
                </>
              )}

              {/* Boundary Optimization Panel - shown when assessment indicates potential improvements */}
              {isGoldenGate && boundaryAssessment && (
                <div className="boundary-optimization-panel mt-4 p-4 rounded-lg" style={{
                  background: boundaryAssessment.needsOptimization
                    ? 'linear-gradient(135deg, #fffbeb 0%, #fef3c7 100%)'
                    : 'linear-gradient(135deg, #f0fdf4 0%, #dcfce7 100%)',
                  border: `1px solid ${boundaryAssessment.needsOptimization ? '#fcd34d' : '#86efac'}`,
                }}>
                  <div className="flex items-start gap-3">
                    <div className={`w-8 h-8 rounded-full flex items-center justify-center flex-shrink-0 ${
                      boundaryAssessment.needsOptimization ? 'bg-amber-400' : 'bg-green-500'
                    }`}>
                      <svg viewBox="0 0 24 24" width="18" height="18" fill="white">
                        {boundaryAssessment.needsOptimization ? (
                          <path d="M12 5.99L19.53 19H4.47L12 5.99M12 2L1 21h22L12 2zm1 14h-2v2h2v-2zm0-6h-2v4h2v-4z"/>
                        ) : (
                          <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
                        )}
                      </svg>
                    </div>
                    <div className="flex-1">
                      <h4 className={`m-0 mb-1 text-[0.938rem] font-semibold ${
                        boundaryAssessment.needsOptimization ? 'text-amber-800' : 'text-green-800'
                      }`}>
                        {boundaryAssessment.needsOptimization
                          ? 'Primer Quality Could Be Improved'
                          : 'Primer Quality is Good'}
                      </h4>
                      <p className={`m-0 mb-3 text-[0.813rem] ${
                        boundaryAssessment.needsOptimization ? 'text-amber-700' : 'text-green-700'
                      }`}>
                        {boundaryAssessment.recommendation}
                      </p>

                      {boundaryAssessment.needsOptimization && boundaryAssessment.issues?.length > 0 && (
                        <div className="mb-3 p-2 bg-white/70 rounded text-xs">
                          <strong className="block mb-1">Issues found:</strong>
                          {boundaryAssessment.issues.slice(0, 3).map((issue: any, i: number) => (
                            <div key={i} className="text-amber-900 ml-2">
                              • Junction {issue.junction + 1} ({issue.side}): {issue.quality} quality
                              {issue.issues?.[0] && ` - ${issue.issues[0]}`}
                            </div>
                          ))}
                          {boundaryAssessment.issues.length > 3 && (
                            <div className="text-amber-900 ml-2 italic">
                              + {boundaryAssessment.issues.length - 3} more issues
                            </div>
                          )}
                        </div>
                      )}

                      {boundaryAssessment.needsOptimization && (
                        <button
                          onClick={handleOptimizeBoundaries}
                          disabled={isOptimizingBoundaries}
                          className={`inline-flex items-center gap-2 py-2 px-4 bg-amber-500 text-white border-none rounded-md text-sm font-medium ${
                            isOptimizingBoundaries ? 'cursor-not-allowed opacity-70' : 'cursor-pointer'
                          }`}
                        >
                          {isOptimizingBoundaries ? (
                            <>
                              <svg className="animate-spin" viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" strokeWidth="2">
                                <circle cx="12" cy="12" r="10" opacity="0.25"/>
                                <path d="M12 2a10 10 0 0 1 10 10" opacity="0.75"/>
                              </svg>
                              Optimizing...
                            </>
                          ) : (
                            <>
                              <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                                <path d="M7.5 5.6L10 7 8.6 4.5 10 2 7.5 3.4 5 2l1.4 2.5L5 7zm12 9.8L17 14l1.4 2.5L17 19l2.5-1.4L22 19l-1.4-2.5L22 14zM22 2l-2.5 1.4L17 2l1.4 2.5L17 7l2.5-1.4L22 7l-1.4-2.5zm-7.63 5.29a.996.996 0 0 0-1.41 0L1.29 18.96a.996.996 0 0 0 0 1.41l2.34 2.34c.39.39 1.02.39 1.41 0L16.7 11.05a.996.996 0 0 0 0-1.41l-2.33-2.35zm-1.03 5.49l-2.12-2.12 2.44-2.44 2.12 2.12-2.44 2.44z"/>
                              </svg>
                              Optimize Boundaries
                            </>
                          )}
                        </button>
                      )}
                    </div>
                  </div>
                </div>
              )}
            </>
          )}
        </div>
      </div>
        </>
      )}

      {/* Enhanced Domestication Panel Modal */}
      {showEnhancedDomestication && enhancedDomesticationPart && (
        <div className="enhanced-domestication-modal-overlay">
          <div className="enhanced-domestication-modal">
            <div className="modal-header">
              <h2>Enhanced Domestication Review</h2>
              <button
                className="modal-close-btn"
                onClick={() => {
                  setShowEnhancedDomestication(false);
                  setEnhancedDomesticationPart(null);
                }}
              >
                <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                  <path d="M19 6.41L17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12z"/>
                </svg>
              </button>
            </div>
            <EnhancedDomesticationPanel
              sequence={enhancedDomesticationPart._originalSequence || enhancedDomesticationPart.seq}
              enzyme={enzyme}
              onDomesticationComplete={(result) => {
                // Update the part with the domesticated sequence
                console.log('[EnhancedDomestication] onComplete result:', {
                  success: result.success,
                  strategy: result.strategy,
                  mutationsCount: result.mutations?.length,
                  junctionsCount: result.junctions?.length,
                  primersCount: result.primers?.length,
                  mutations: result.mutations,
                });
                if (result.success) {
                  const updatedParts = [...parts];
                  const mutationCount = result.mutations?.length || 0;
                  const junctionCount = result.junctions?.length || 0;

                  // Store original sequence if not already stored (for reconfiguration)
                  const originalSeq = enhancedDomesticationPart._originalSequence || enhancedDomesticationPart.seq;

                  const updatedPart = {
                    ...updatedParts[enhancedDomesticationPart.index || 0],
                    seq: result.domesticatedSequence || enhancedDomesticationPart.seq,
                    _originalSequence: originalSeq,
                    _domesticationApproved: true,
                    _domesticationMutations: result.mutations || [],
                    _domesticationJunctions: result.junctions || [],
                    _domesticationPrimers: result.primers || [],
                    _domesticationStrategy: result.strategy,
                  };
                  console.log('[EnhancedDomestication] Storing updated part:', {
                    index: enhancedDomesticationPart.index,
                    _domesticationApproved: updatedPart._domesticationApproved,
                    _domesticationMutations: updatedPart._domesticationMutations,
                    _domesticationJunctions: updatedPart._domesticationJunctions,
                    _domesticationStrategy: updatedPart._domesticationStrategy,
                  });
                  updatedParts[enhancedDomesticationPart.index || 0] = updatedPart;
                  setParts(updatedParts);

                  // Set flag to trigger redesign after state update
                  setNeedsRedesignAfterDomestication(true);

                  // Build appropriate message based on strategy
                  let message;
                  if (junctionCount > 0) {
                    message = `Domestication approved: ${junctionCount} mutagenic junction(s) designed`;
                    if (mutationCount > 0) {
                      message += ` + ${mutationCount} silent mutation(s)`;
                    }
                  } else if (mutationCount > 0) {
                    message = `Domestication approved: ${mutationCount} silent mutation(s) applied`;
                  } else {
                    message = 'Domestication approved: No changes needed';
                  }
                  setCopyMessage(message);
                }
                setShowEnhancedDomestication(false);
                setEnhancedDomesticationPart(null);
              }}
              onCancel={() => {
                setShowEnhancedDomestication(false);
                setEnhancedDomesticationPart(null);
              }}
            />
          </div>
          <style>{EnhancedDomesticationStyles}</style>
        </div>
      )}

      {/* Domestication Workflow Guide Modal */}
      {showWorkflowGuide && workflowGuidePart && (
        <div className="domestication-workflow-modal-overlay">
          <div className="domestication-workflow-modal">
            <DomesticationWorkflowGuide
              sequence={workflowGuidePart.seq}
              enzyme={enzyme}
              existingOverhangs={result?.overhangs || []}
              onWorkflowComplete={(workflowResult: any) => {
                // Apply the workflow results
                if (workflowResult.success && workflowResult.domestication?.domesticatedSequence) {
                  const updatedParts = [...parts];
                  const partIndex = parts.findIndex(p => p === workflowGuidePart);
                  if (partIndex >= 0) {
                    updatedParts[partIndex] = {
                      ...workflowGuidePart,
                      seq: workflowResult.domestication.domesticatedSequence,
                      _workflowApplied: true,
                      _workflowResult: workflowResult,
                    };
                    setParts(updatedParts);
                  }
                }
                setShowWorkflowGuide(false);
                setWorkflowGuidePart(null);
              }}
              onCancel={() => {
                setShowWorkflowGuide(false);
                setWorkflowGuidePart(null);
              }}
            />
          </div>
          <style>{DomesticationWorkflowStyles}</style>
        </div>
      )}
    </div>
  );
}
