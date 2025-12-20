import React, { useState, useCallback, useMemo, useRef, useEffect } from 'react';
import {
  designUnified,
  designBatch,
  parseNotationToSpec,
  parseBatchText,
  MUTAGENESIS_DEFAULTS,
  CODON_TABLE,
  AA_NAMES,
  analyzePrimerPair,
} from '../lib/unifiedPrimerDesign.js';
import {
  checkHeterodimer,
  MUTATION_TYPES,
  DIMER_THRESHOLDS,  // Evidence-based thresholds (IDT, Premier Biosoft)
} from '../lib/index.js';
import { identifyBadges } from '../lib/diversitySelection.js';
import { getTmDiffColorClass, TM_DIFF_THRESHOLDS } from '../lib/primerAlternativesUtils.js';
import {
  COMMON_VECTORS,
  parseSequenceFile,
  reverseComplement,
  get3FrameTranslation,
  hasAmbiguousBases,
  countCombinations,
  downloadFile,
  generateCSVSummary,
  generatePrimerOrderList,
  generateFasta,
} from '../lib/sequenceUtils.js';
import PrimerOnTemplateViewer from './PrimerOnTemplateViewer.jsx';
import HairpinDiagram from './HairpinDiagram.jsx';
import { AlternativesPanel, EnhancedAnalysisSection } from './primers/index.js';
import ScoreBreakdownPopup from './primers/ScoreBreakdownPopup.jsx';

// Stub for missing export
const compareTmMethods = (...args: any[]): any => ({});

// Type definitions
interface SequenceInfo {
  format: string;
  id: string;
  description: string;
}

interface PrimerData {
  sequence: string;
  tm: number;
  gc: number;
  hasGCClamp?: boolean;
  gcPercent?: string;
  length?: number;
  start?: number;
  end?: number;
  terminal3DG?: number;
  hairpinDG?: number;
  selfDimerDG?: number;
  gcContent?: number;
}

interface DesignResult {
  forward: PrimerData;
  reverse: PrimerData;
  quality?: string;
  qualityTier?: string;
  compositeScore?: number;
  effectiveScore?: number;
  criticalWarnings?: number;
  tmDifference?: number;
  annealingTemp?: number;
  description?: string;
  type?: string;
  design?: string;
  strategy?: string;
  isLibrary?: boolean;
  protocol?: Protocol;
  alternateDesigns?: AlternateDesign[];
  alternativePrimers?: AlternateDesign[];
  alternativeCategories?: any;
  isAlternativeSelected?: boolean;
  selectedAlternateIdx?: number;
  wasUpgraded?: boolean;
  originalScore?: number;
  originalSequence?: string;
  mutatedSequence?: string;
  position?: number;
  nucleotidePosition?: number;
  deleteLength?: number;
  deletionLength?: number;
  deletedSequence?: string;
  insertSequence?: string;
  replacement?: string;
  oldCodon?: string;
  newCodon?: string;
  codonUsage?: number;
  forwardPiecewiseScores?: any;
  reversePiecewiseScores?: any;
}

interface AlternateDesign {
  forward: PrimerData;
  reverse: PrimerData;
  compositeScore?: number;
  score?: number;
  qualityTier?: string;
  design?: string;
  heterodimerDG?: number;
  originalIdx?: number;
}

interface Protocol {
  name: string;
  steps?: ProtocolStep[];
  notes?: string[];
}

interface ProtocolStep {
  name: string;
  temp?: string;
  time?: string;
  substeps?: SubStep[];
}

interface SubStep {
  name: string;
  temp: string;
  time: string;
}

interface EnhancedAnalysis {
  forward?: {
    hairpinDG?: number;
    selfDimerDG?: number;
  };
  reverse?: {
    hairpinDG?: number;
    selfDimerDG?: number;
  };
  heterodimer?: {
    heterodimerDG?: number;
    severity?: string;
    dimerType?: string;
    isExpectedOverlap?: boolean;
    overlapLength?: number;
  };
  offTargets?: {
    forward?: {
      offTargetCount?: number;
    };
    reverse?: {
      offTargetCount?: number;
    };
  } | null;
  fwdTmComparison?: any;
}

interface CollapsedSections {
  summary: boolean;
  primers: boolean;
  analysis: boolean;
  pcr: boolean;
  visualization: boolean;
  protocol: boolean;
  alternates: boolean;
}

interface DesignOptions {
  strategy: string;
  optimalTm: number;
  minTm: number;
  maxTm: number;
  minPrimerLength: number;
  maxPrimerLength: number;
  minGC: number;
  maxGC: number;
  circular: boolean;
  useSmartDesign: boolean;
  confineTo5Tails: boolean;
  exhaustiveSearch: boolean;
}

interface LibraryInfo {
  isLibrary: boolean;
  size: number;
}

interface SelectionInfo {
  start: number;
  end: number;
  length: number;
  sequence: string;
  isWrapped: boolean;
}

interface BatchItem {
  id: number;
  label: string;
  start: number;
  end: number;
  replacement?: string | null;
  aaHelper?: {
    newAA: string;
    codonPosition: number;
    organism: string;
  };
}

interface BatchResult {
  success: boolean;
  forward?: PrimerData;
  reverse?: PrimerData;
  originalSpec?: BatchItem;
  description?: string;
  error?: string;
}

interface StatusCheck {
  label: string;
  status: string;
  value: string;
  section: string;
  threshold: string;
}

interface SummaryStatus {
  checks: Record<string, StatusCheck>;
  overall: string;
}

interface AlternatesSort {
  field: string;
  direction: string;
}

interface AlternatesFilters {
  minScore: number;
  maxTmDiff: number;
  requireGcClamp: boolean;
}

interface Section {
  id: string;
  label: string;
  icon: string;
}

/**
 * Unified Primer Designer Component
 *
 * A streamlined UI based on the unified design model:
 * - All operations are SELECT region → REPLACE with something
 * - AA mutation is a helper for codon optimization
 * - Batch is built-in (add more modifications)
 */
export default function UnifiedPrimerDesigner() {
  // Template state
  const [template, setTemplate] = useState<string>('');
  const [sequenceInfo, setSequenceInfo] = useState<SequenceInfo | null>(null);
  const [selectedVector, setSelectedVector] = useState<string>('');
  const [useReverseComplement, setUseReverseComplement] = useState<boolean>(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Selection state
  const [selectionStart, setSelectionStart] = useState<string>('');
  const [selectionEnd, setSelectionEnd] = useState<string>('');

  // Replacement mode
  const [replacementMode, setReplacementMode] = useState<string>('amplify'); // amplify | delete | aa | direct
  const [directSequence, setDirectSequence] = useState<string>('');
  const [aaOldAA, setAaOldAA] = useState<string>('');
  const [aaPosition, setAaPosition] = useState<string>('');
  const [aaNewAA, setAaNewAA] = useState<string>('');
  const [organism, setOrganism] = useState<string>('ecoli');
  const [orfStart, setOrfStart] = useState<string>('1');

  // Batch state
  const [batchItems, setBatchItems] = useState<BatchItem[]>([]);
  const [showBatch, setShowBatch] = useState<boolean>(false);

  // Results
  const [results, setResults] = useState<DesignResult | null>(null);
  const [originalDesign, setOriginalDesign] = useState<DesignResult | null>(null);
  const [batchResults, setBatchResults] = useState<BatchResult[]>([]);
  const [enhancedAnalysis, setEnhancedAnalysis] = useState<EnhancedAnalysis | null>(null);
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);

  // Collapsible sections (matches legacy MutagenesisDesigner)
  const [collapsedSections, setCollapsedSections] = useState<CollapsedSections>({
    summary: false,
    primers: false,
    analysis: false,
    pcr: false,
    visualization: false,
    protocol: true,  // Collapsed by default
    alternates: false,
  });
  const [activeSection, setActiveSection] = useState<string | null>(null);
  const [showThermodynamics, setShowThermodynamics] = useState<boolean>(false);
  const [showScoreBreakdown, setShowScoreBreakdown] = useState<boolean>(false);
  const resultsRef = useRef<HTMLDivElement>(null);

  // Alternative options view state
  const [alternatesViewMode, setAlternatesViewMode] = useState<string>('card'); // 'card' | 'table'
  const [showAllAlternates, setShowAllAlternates] = useState<boolean>(false);
  const [expandedSequences, setExpandedSequences] = useState<Set<number>>(new Set());
  const [compareSelection, setCompareSelection] = useState<Set<number>>(new Set());
  const [alternatesSort, setAlternatesSort] = useState<AlternatesSort>({ field: 'score', direction: 'desc' });
  const [alternatesFilters, setAlternatesFilters] = useState<AlternatesFilters>({
    minScore: 0,
    maxTmDiff: 10,
    requireGcClamp: false,
  });
  const [showCompareModal, setShowCompareModal] = useState<boolean>(false);
  const [toastMessage, setToastMessage] = useState<string | null>(null);
  const [dismissedBanner, setDismissedBanner] = useState<boolean>(false);
  const [copiedFeedback, setCopiedFeedback] = useState<string | null>(null); // 'fwd' | 'rev' | 'pair' | null

  // Advanced options
  const [showAdvanced, setShowAdvanced] = useState<boolean>(false);
  const [designOptions, setDesignOptions] = useState<DesignOptions>({
    strategy: 'back-to-back',
    optimalTm: 62,
    minTm: 55,
    maxTm: 72,
    minPrimerLength: 15,
    maxPrimerLength: 60,
    minGC: 40,  // Percentage (will convert to decimal for API)
    maxGC: 60,
    circular: true,  // Most templates are plasmids
    useSmartDesign: true,
    confineTo5Tails: false,  // Keep mutations in 5' non-annealing tail
    exhaustiveSearch: true,  // Exhaustive search enabled by default for optimal results
  });
  const [loadingMessage, setLoadingMessage] = useState<string>('');
  const [isPreviewResult, setIsPreviewResult] = useState<boolean>(false);  // Track if showing preview vs final result

  // Computed template sequence
  const templateSeq = useMemo(() => {
    const cleaned = template.toUpperCase().replace(/[^ATGCRYSWKMBDHVN]/g, '');
    return useReverseComplement ? reverseComplement(cleaned) : cleaned;
  }, [template, useReverseComplement]);

  // Computed selection info - handles circular wrap-around when start > end
  const selectionInfo = useMemo((): SelectionInfo | null => {
    const start = parseInt(selectionStart);
    const end = parseInt(selectionEnd);

    if (isNaN(start) || isNaN(end) || !templateSeq) return null;

    const seqLen = templateSeq.length;
    const clampedStart = Math.max(0, Math.min(seqLen - 1, start));
    const clampedEnd = Math.max(0, Math.min(seqLen, end));
    const isCircular = designOptions.circular;
    const isWrapped = clampedStart > clampedEnd && isCircular;

    // Calculate length based on whether selection wraps around origin
    let length: number;
    let sequence: string;
    if (isWrapped) {
      // Circular wrap-around: length is (seqLen - start) + end
      length = (seqLen - clampedStart) + clampedEnd;
      sequence = templateSeq.slice(clampedStart) + templateSeq.slice(0, clampedEnd);
    } else if (clampedStart > clampedEnd && !isCircular) {
      // Linear sequence with invalid selection (start > end) - show 0
      length = 0;
      sequence = '';
    } else {
      // Normal linear selection
      length = clampedEnd - clampedStart;
      sequence = templateSeq.slice(clampedStart, clampedEnd);
    }

    return {
      start: clampedStart,
      end: clampedEnd,
      length,
      sequence,
      isWrapped,
    };
  }, [selectionStart, selectionEnd, templateSeq, designOptions.circular]);

  // Library info for direct mode
  const libraryInfo = useMemo((): LibraryInfo | null => {
    if (!directSequence || !hasAmbiguousBases(directSequence)) return null;
    return {
      isLibrary: true,
      size: countCombinations(directSequence),
    };
  }, [directSequence]);

  // Compute summary status for key metrics
  const summaryStatus = useMemo((): SummaryStatus | null => {
    if (!results || !results.forward) return null;

    const checks: Record<string, StatusCheck> = {
      quality: {
        label: 'Primer Quality',
        status: results.quality === 'excellent' || results.quality === 'good' ? 'pass' :
                results.quality === 'acceptable' ? 'warning' : 'fail',
        value: results.quality || results.qualityTier || 'N/A',
        section: 'primers',
        threshold: 'Excellent/Good = optimal, Acceptable = may need optimization, Poor = redesign recommended',
      },
      tmDiff: {
        label: 'Tm Difference',
        status: (results.tmDifference || Math.abs(results.forward.tm - results.reverse.tm)) <= 2 ? 'pass' :
                (results.tmDifference || Math.abs(results.forward.tm - results.reverse.tm)) <= 5 ? 'warning' : 'fail',
        value: `${(results.tmDifference || Math.abs(results.forward.tm - results.reverse.tm)).toFixed(1)}°C`,
        section: 'pcr',
        threshold: '≤2°C = ideal, 2-5°C = acceptable, >5°C = may cause uneven amplification',
      },
      gcClamp: {
        label: 'GC Clamp',
        status: results.forward.hasGCClamp && results.reverse.hasGCClamp ? 'pass' : 'warning',
        value: results.forward.hasGCClamp && results.reverse.hasGCClamp ? 'Both' :
               results.forward.hasGCClamp || results.reverse.hasGCClamp ? 'Partial' : 'None',
        section: 'primers',
        threshold: 'G or C at 3\' end improves specificity. Both primers having GC clamp is ideal.',
      },
    };

    // Check for GC clamp on forward/reverse
    if (results.forward.sequence) {
      const fwdLast = results.forward.sequence[results.forward.sequence.length - 1];
      results.forward.hasGCClamp = fwdLast === 'G' || fwdLast === 'C';
    }
    if (results.reverse.sequence) {
      const revLast = results.reverse.sequence[results.reverse.sequence.length - 1];
      results.reverse.hasGCClamp = revLast === 'G' || revLast === 'C';
    }

    // Add enhanced analysis checks if available
    if (enhancedAnalysis) {
      // Heterodimer check - now uses unified severity classification including 3' involvement
      // Special handling for overlapping (QuikChange) designs where full complementarity is expected
      const isExpectedOverlap = enhancedAnalysis.heterodimer?.isExpectedOverlap;
      const heteroSeverity = enhancedAnalysis.heterodimer?.severity || 'safe';
      const heterodimerStatus = isExpectedOverlap ? 'pass' :
                                heteroSeverity === 'critical' ? 'fail' :
                                heteroSeverity === 'warning' ? 'warning' : 'pass';
      const heterodimerType = enhancedAnalysis.heterodimer?.dimerType || 'none';
      const heteroLabel = isExpectedOverlap ? 'Primer Overlap (QuikChange)' :
                          heterodimerType === '3prime_extensible' ? "3' Extensible Dimer" :
                          heterodimerType === 'internal_with_3prime' ? "Heterodimer (3' involved)" :
                          'Heterodimer';
      // For overlapping designs, show overlap info instead of ΔG (which is expected to be very negative)
      const heteroValue = isExpectedOverlap
        ? `${enhancedAnalysis.heterodimer?.overlapLength || 'Full'} bp overlap`
        : enhancedAnalysis.heterodimer?.heterodimerDG?.toFixed(1) + ' kcal/mol';
      const heteroThreshold = isExpectedOverlap
        ? 'Full primer complementarity is expected for QuikChange-style mutagenesis'
        : heterodimerType === '3prime_extensible'
          ? "CRITICAL: Both primers' 3' ends hybridize - causes primer-dimer artifacts"
          : '>-6 kcal/mol = OK. 3\' end involvement is more critical than ΔG alone.';
      checks.heterodimer = {
        label: heteroLabel,
        status: heterodimerStatus,
        value: heteroValue,
        section: 'analysis',
        threshold: heteroThreshold,
      };

      // Hairpin/Secondary Structure check (now uses accurate Zuker foldDG)
      // Evidence-based thresholds from Premier Biosoft (internal hairpin)
      const fwdHairpin = enhancedAnalysis.forward?.hairpinDG;
      const revHairpin = enhancedAnalysis.reverse?.hairpinDG;
      const worstHairpin = Math.min(fwdHairpin || 0, revHairpin || 0);
      const hairpinThresh = DIMER_THRESHOLDS.hairpin.internal;
      checks.secondaryStructure = {
        label: 'Hairpin',
        status: worstHairpin >= hairpinThresh.ideal ? 'pass' :
                worstHairpin >= hairpinThresh.critical ? 'warning' : 'fail',
        value: `${worstHairpin?.toFixed(1) || 'N/A'} kcal/mol`,
        section: 'analysis',
        threshold: `>${hairpinThresh.ideal} = ideal (Premier Biosoft), ${hairpinThresh.ideal} to ${hairpinThresh.critical} = monitor, <${hairpinThresh.critical} = problematic`,
      };

      // Self-dimer check - Evidence-based thresholds from IDT/Premier Biosoft (internal)
      const fwdSelfDimer = enhancedAnalysis.forward?.selfDimerDG;
      const revSelfDimer = enhancedAnalysis.reverse?.selfDimerDG;
      const worstSelfDimer = Math.min(fwdSelfDimer || 0, revSelfDimer || 0);
      const selfDimerThresh = DIMER_THRESHOLDS.selfDimer.internal;
      checks.selfDimer = {
        label: 'Self-Dimer',
        status: worstSelfDimer >= selfDimerThresh.ideal ? 'pass' :
                worstSelfDimer >= selfDimerThresh.critical ? 'warning' : 'fail',
        value: `${worstSelfDimer?.toFixed(1) || 'N/A'} kcal/mol`,
        section: 'analysis',
        threshold: `>${selfDimerThresh.ideal} = OK (Premier Biosoft), ${selfDimerThresh.ideal} to ${selfDimerThresh.critical} = monitor, <${selfDimerThresh.critical} = IDT critical threshold`,
      };

      if (enhancedAnalysis.offTargets) {
        const maxOffTargets = Math.max(
          enhancedAnalysis.offTargets.forward?.offTargetCount || 0,
          enhancedAnalysis.offTargets.reverse?.offTargetCount || 0
        );
        checks.offTargets = {
          label: 'Off-Targets',
          status: maxOffTargets === 0 ? 'pass' : maxOffTargets <= 2 ? 'warning' : 'fail',
          value: maxOffTargets.toString(),
          section: 'analysis',
          threshold: '0 = specific binding only, 1-2 = minor risk of non-specific products, >2 = may need redesign',
        };
      }
    }

    // Calculate overall status
    const statuses = Object.values(checks).map(c => c.status);
    const overall = statuses.includes('fail') ? 'fail' :
                    statuses.includes('warning') ? 'warning' : 'pass';

    return { checks, overall };
  }, [results, enhancedAnalysis]);

  // Toggle section collapse state
  const toggleSection = useCallback((sectionId: string) => {
    setCollapsedSections(prev => ({
      ...prev,
      [sectionId]: !prev[sectionId as keyof CollapsedSections]
    }));
  }, []);

  // Scroll to section
  const scrollToSection = useCallback((sectionId: string) => {
    const element = document.getElementById(`section-${sectionId}`);
    if (element) {
      element.scrollIntoView({ behavior: 'smooth', block: 'start' });
      setActiveSection(sectionId);
      // Expand section if collapsed
      setCollapsedSections(prev => ({
        ...prev,
        [sectionId]: false
      }));
    }
  }, []);

  // Define available sections for navigation (matches legacy MutagenesisDesigner)
  const getAvailableSections = useCallback((): Section[] => {
    if (!results || results.isLibrary) return [];

    const sections: Section[] = [
      { id: 'summary', label: 'Summary', icon: '📊' },
      { id: 'primers', label: 'Designed Primers', icon: '🧬' },
    ];

    if (enhancedAnalysis) {
      sections.push({ id: 'analysis', label: 'Analysis', icon: '🔬' });
    }

    sections.push({ id: 'pcr', label: 'PCR Parameters', icon: '🌡️' });

    if (templateSeq && results.forward && results.reverse) {
      sections.push({ id: 'visualization', label: 'Visualization', icon: '📐' });
    }

    if (results.protocol) {
      sections.push({ id: 'protocol', label: 'Protocol', icon: '📋' });
    }

    // Sequence comparison is now integrated into the primers section

    if (results.alternateDesigns?.length && results.alternateDesigns.length > 0) {
      sections.push({ id: 'alternates', label: 'Alternate Designs', icon: '🔄' });
    }

    if (results.alternativePrimers?.length && results.alternativePrimers.length > 0) {
      sections.push({ id: 'alternates', label: 'Alternative Primers', icon: '🔄' });
    }

    return sections;
  }, [results, enhancedAnalysis, templateSeq]);

  // Track active section on scroll (matches legacy MutagenesisDesigner)
  useEffect(() => {
    if (!results) return;

    const handleScroll = () => {
      const sections = ['summary', 'primers', 'analysis', 'pcr', 'visualization', 'protocol', 'alternates'];
      let currentSection: string | null = null;

      for (const sectionId of sections) {
        const element = document.getElementById(`section-${sectionId}`);
        if (element) {
          const rect = element.getBoundingClientRect();
          if (rect.top <= 150) {
            currentSection = sectionId;
          }
        }
      }

      if (currentSection !== activeSection) {
        setActiveSection(currentSection);
      }
    };

    window.addEventListener('scroll', handleScroll);
    return () => window.removeEventListener('scroll', handleScroll);
  }, [results, activeSection]);

  // Handle file upload
  const handleFileUpload = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e: ProgressEvent<FileReader>) => {
      try {
        const content = e.target?.result as string;
        const parsed = parseSequenceFile(content);

        if (parsed.entries.length > 0) {
          const entry = parsed.entries[0];
          setTemplate(entry.sequence);
          setSequenceInfo({
            format: parsed.format,
            id: entry.id,
            description: entry.description,
          });
          setError(null);
        }
      } catch (err) {
        setError('Failed to parse file: ' + (err as Error).message);
      }
    };
    reader.readAsText(file);
  }, []);

  // Handle vector selection
  const handleVectorSelect = useCallback((vectorName: string) => {
    setSelectedVector(vectorName);
    if (vectorName && COMMON_VECTORS[vectorName]?.sequence) {
      setTemplate(COMMON_VECTORS[vectorName].sequence);
      setSequenceInfo({
        format: 'vector',
        id: vectorName,
        description: COMMON_VECTORS[vectorName].description,
      });
    }
  }, []);

  // Build specification from current UI state
  const buildSpec = useCallback(() => {
    const start = parseInt(selectionStart) || 0;
    const end = parseInt(selectionEnd) || templateSeq.length;

    switch (replacementMode) {
      case 'amplify':
        return { start, end, replacement: null };

      case 'delete':
        return { start, end, replacement: '' };

      case 'aa': {
        const codonPos = parseInt(aaPosition);
        const orfOffset = (parseInt(orfStart) || 1) - 1;
        const nucStart = orfOffset + (codonPos - 1) * 3;

        return {
          start: nucStart,
          end: nucStart + 3,
          aaHelper: {
            newAA: aaNewAA.toUpperCase(),
            codonPosition: codonPos,
            organism,
          },
        };
      }

      case 'direct':
        return { start, end, replacement: directSequence.toUpperCase() };

      default:
        throw new Error(`Unknown replacement mode: ${replacementMode}`);
    }
  }, [replacementMode, selectionStart, selectionEnd, templateSeq, aaPosition, aaNewAA, orfStart, organism, directSequence]);

  // Helper to process design result and add metadata
  const processDesignResult = useCallback((result: DesignResult): DesignResult => {
    // Add GC clamp info
    if (result.forward?.sequence) {
      const fwdLast = result.forward.sequence[result.forward.sequence.length - 1];
      result.forward.hasGCClamp = fwdLast === 'G' || fwdLast === 'C';
      result.forward.gcPercent = `${((result.forward.gc || 0) * 100).toFixed(1)}%`;
    }
    if (result.reverse?.sequence) {
      const revLast = result.reverse.sequence[result.reverse.sequence.length - 1];
      result.reverse.hasGCClamp = revLast === 'G' || revLast === 'C';
      result.reverse.gcPercent = `${((result.reverse.gc || 0) * 100).toFixed(1)}%`;
    }

    // Calculate Tm difference and annealing temp
    if (result.forward?.tm && result.reverse?.tm) {
      result.tmDifference = Math.abs(result.forward.tm - result.reverse.tm);
      result.annealingTemp = Math.round(Math.min(Math.min(result.forward.tm, result.reverse.tm) + 1, 72));
    }

    return result;
  }, []);

  // Helper to run enhanced analysis
  const runEnhancedAnalysis = useCallback((result: DesignResult) => {
    if (result.forward && result.reverse) {
      try {
        const analysis = analyzePrimerPair(
          result.forward.sequence,
          result.reverse.sequence,
          templateSeq
        );
        const heterodimerCheck = checkHeterodimer(
          result.forward.sequence,
          result.reverse.sequence,
          { designType: result.design } as any
        );
        let fwdTmComparison = null;
        if (result.forward.start !== undefined) {
          try {
            fwdTmComparison = compareTmMethods(
              result.forward.sequence,
              templateSeq.slice(result.forward.start, result.forward.start + result.forward.sequence.length)
            );
          } catch (e) {
            // Tm comparison is optional
          }
        }
        setEnhancedAnalysis({
          ...analysis,
          heterodimer: heterodimerCheck,
          ...(fwdTmComparison && { fwdTmComparison }),
        });
      } catch (analysisErr) {
        console.warn('Enhanced analysis failed:', analysisErr);
      }
    }
  }, [templateSeq]);

  // Handle single design with progressive results
  // When exhaustive search is enabled, first show quick preview, then update with optimal result
  const handleDesign = useCallback(() => {
    if (!templateSeq || templateSeq.length < 50) {
      setError('Template must be at least 50 bp');
      return;
    }

    setLoading(true);
    setError(null);
    setResults(null);
    setEnhancedAnalysis(null);
    setIsPreviewResult(false);

    const spec = buildSpec();
    // Convert GC percentage to decimal for API
    const apiOptions = {
      ...designOptions,
      minGC: designOptions.minGC / 100,
      maxGC: designOptions.maxGC / 100,
    };

    // Progressive search: show quick preview first, then run exhaustive search
    if (designOptions.exhaustiveSearch) {
      setLoadingMessage('Quick preview...');

      // Step 1: Quick search for immediate feedback
      setTimeout(() => {
        try {
          const quickOptions = { ...apiOptions, exhaustiveSearch: false };
          const quickResult = designUnified(templateSeq, spec as any, quickOptions);
          const processedQuick = processDesignResult(quickResult);

          // Show preview results immediately
          setResults(processedQuick);
          setOriginalDesign(processedQuick);
          setIsPreviewResult(true);
          setLoadingMessage('Optimizing... (finding best design)');

          // Run enhanced analysis on preview
          runEnhancedAnalysis(processedQuick);

          // Step 2: Run exhaustive search in background
          setTimeout(() => {
            try {
              const exhaustiveResult = designUnified(templateSeq, spec as any, apiOptions);
              const processedExhaustive = processDesignResult(exhaustiveResult);

              // Update with final optimal result
              setResults(processedExhaustive);
              setOriginalDesign(processedExhaustive);
              setIsPreviewResult(false);

              // Update enhanced analysis with final result
              runEnhancedAnalysis(processedExhaustive);
            } catch (err) {
              // Keep preview results if exhaustive fails
              console.warn('Exhaustive search failed, keeping preview:', err);
              setIsPreviewResult(false);  // Still mark as final
            } finally {
              setLoading(false);
              setLoadingMessage('');
            }
          }, 10);
        } catch (err) {
          setError((err as Error).message);
          setLoading(false);
          setLoadingMessage('');
        }
      }, 10);
    } else {
      // Standard mode: single pass
      setLoadingMessage('Designing primers...');

      setTimeout(() => {
        try {
          const result = designUnified(templateSeq, spec as any, apiOptions);
          const processedResult = processDesignResult(result);

          setResults(processedResult);
          setOriginalDesign(processedResult);
          runEnhancedAnalysis(processedResult);
        } catch (err) {
          setError((err as Error).message);
        } finally {
          setLoading(false);
          setLoadingMessage('');
        }
      }, 10);
    }
  }, [templateSeq, buildSpec, designOptions, processDesignResult, runEnhancedAnalysis]);

  // Add to batch
  const handleAddToBatch = useCallback(() => {
    try {
      const spec = buildSpec();

      // Create a label for this item
      let label: string;
      if ((spec as any).aaHelper) {
        label = `${aaOldAA || '?'}${aaPosition}${aaNewAA}`;
      } else if (spec.replacement === null) {
        label = `Amplify ${spec.start + 1}-${spec.end}`;
      } else if (spec.replacement === '') {
        label = `Delete ${spec.start + 1}-${spec.end}`;
      } else if (spec.start === spec.end) {
        label = `Insert at ${spec.start + 1}`;
      } else {
        label = `Replace ${spec.start + 1}-${spec.end}`;
      }

      setBatchItems(prev => [...prev, { ...spec, label, id: Date.now() } as BatchItem]);
      setShowBatch(true);
    } catch (err) {
      setError((err as Error).message);
    }
  }, [buildSpec, aaOldAA, aaPosition, aaNewAA]);

  // Remove from batch
  const handleRemoveFromBatch = useCallback((id: number) => {
    setBatchItems(prev => prev.filter(item => item.id !== id));
  }, []);

  // Design all batch items
  const handleDesignBatch = useCallback(() => {
    if (!templateSeq || templateSeq.length < 50) {
      setError('Template must be at least 50 bp');
      return;
    }

    if (batchItems.length === 0) {
      setError('No items in batch');
      return;
    }

    setLoading(true);
    setError(null);
    setBatchResults([]);

    setTimeout(() => {
      try {
        const results = designBatch(templateSeq, batchItems, designOptions);
        setBatchResults(results as any);
      } catch (err) {
        setError((err as Error).message);
      } finally {
        setLoading(false);
      }
    }, 10);
  }, [templateSeq, batchItems, designOptions]);

  // Download batch CSV
  const handleDownloadCSV = useCallback(() => {
    const successfulResults = batchResults.filter(r => r.success);
    if (successfulResults.length === 0) return;

    const rows = ['Mutation,Direction,Sequence,Tm,GC,Length'];
    for (const result of successfulResults) {
      const label = result.originalSpec?.label || result.description || '';
      rows.push(`${label},Forward,${result.forward!.sequence},${result.forward!.tm?.toFixed(1) || ''},${((result.forward!.gc || 0) * 100).toFixed(1)}%,${result.forward!.sequence.length}`);
      rows.push(`${label},Reverse,${result.reverse!.sequence},${result.reverse!.tm?.toFixed(1) || ''},${((result.reverse!.gc || 0) * 100).toFixed(1)}%,${result.reverse!.sequence.length}`);
    }

    downloadFile(rows.join('\n'), 'batch_primers.csv', 'text/csv');
  }, [batchResults]);

  // Copy primer to clipboard
  const copyToClipboard = useCallback((text: string) => {
    navigator.clipboard.writeText(text);
  }, []);

  // Download single result handlers
  const handleDownloadPrimers = useCallback(() => {
    if (!results?.forward) return;
    const primers = [
      { name: 'Forward', sequence: results.forward.sequence, note: results.description || '' },
      { name: 'Reverse', sequence: results.reverse.sequence, note: '' },
    ];
    const content = generatePrimerOrderList(primers);
    downloadFile(content, 'primers.txt');
  }, [results]);

  const handleDownloadFasta = useCallback(() => {
    if (!results?.forward) return;
    const sequences = [
      { id: 'Forward_primer', description: results.description || '', sequence: results.forward.sequence },
      { id: 'Reverse_primer', description: '', sequence: results.reverse.sequence },
    ];
    if (results.mutatedSequence) {
      sequences.push({ id: 'Mutated_sequence', description: results.description || '', sequence: results.mutatedSequence });
    }
    const content = generateFasta(sequences);
    downloadFile(content, 'primers.fasta');
  }, [results]);

  const handleDownloadSingleCSV = useCallback(() => {
    if (!results?.forward) return;
    const content = generateCSVSummary([results as any]);
    downloadFile(content, 'primers.csv', 'text/csv');
  }, [results]);

  // Show toast notification
  const showToast = useCallback((message: string, duration: number = 3000) => {
    setToastMessage(message);
    setTimeout(() => setToastMessage(null), duration);
  }, []);

  // Select an alternative design
  const handleSelectAlternate = useCallback((alt: AlternateDesign) => {
    if (!results) return;

    const fwdTm = alt.forward.tm;
    const revTm = alt.reverse.tm;
    const lowerTm = Math.min(fwdTm, revTm);
    const newAnnealingTemp = Math.round(Math.min(lowerTm + 1, 72));
    const newTmDifference = Math.abs(fwdTm - revTm);

    const score = alt.compositeScore || 70;
    const qualityTier = alt.qualityTier || (score >= 80 ? 'excellent' :
                        score >= 70 ? 'good' :
                        score >= 60 ? 'acceptable' : 'poor');

    const updatedResults: DesignResult = {
      ...results,
      forward: alt.forward,
      reverse: alt.reverse,
      quality: qualityTier,
      compositeScore: alt.compositeScore,
      annealingTemp: newAnnealingTemp,
      tmDifference: newTmDifference,
      isAlternativeSelected: true,
    };

    setResults(updatedResults);

    // Show toast notification for feedback
    showToast(`Selected alternative primer pair (Score: ${alt.compositeScore})`);

    // Run enhanced analysis on new primers
    try {
      const analysis = analyzePrimerPair(
        alt.forward.sequence,
        alt.reverse.sequence,
        templateSeq
      );
      const heterodimerCheck = checkHeterodimer(
        alt.forward.sequence,
        alt.reverse.sequence,
        { designType: alt.design || results.design } as any
      );
      setEnhancedAnalysis({
        ...analysis,
        heterodimer: heterodimerCheck,
      });
    } catch (analysisErr) {
      console.warn('Enhanced analysis failed:', analysisErr);
    }
  }, [results, templateSeq, showToast]);

  // Revert to original design
  const handleRevertDesign = useCallback(() => {
    if (!originalDesign) return;

    setResults(originalDesign);

    // Run enhanced analysis on original primers
    try {
      const analysis = analyzePrimerPair(
        originalDesign.forward.sequence,
        originalDesign.reverse.sequence,
        templateSeq
      );
      const heterodimerCheck = checkHeterodimer(
        originalDesign.forward.sequence,
        originalDesign.reverse.sequence,
        { designType: originalDesign.design } as any
      );
      setEnhancedAnalysis({
        ...analysis,
        heterodimer: heterodimerCheck,
      });
    } catch (analysisErr) {
      console.warn('Enhanced analysis failed:', analysisErr);
    }

    scrollToSection('primers');
  }, [originalDesign, templateSeq, scrollToSection]);

  // Copy sequence to clipboard with feedback
  const handleCopySequence = useCallback(async (sequence: string, label: string = 'Sequence') => {
    try {
      await navigator.clipboard.writeText(sequence);
      setCopiedFeedback(label);
      setTimeout(() => setCopiedFeedback(null), 1500);
    } catch (err) {
      console.warn('Copy failed:', err);
    }
  }, []);

  // Copy primer pair to clipboard with feedback
  const handleCopyPrimerPair = useCallback(async (fwd: string, rev: string) => {
    const text = `Forward: ${fwd}\nReverse: ${rev}`;
    try {
      await navigator.clipboard.writeText(text);
      setCopiedFeedback('pair');
      setTimeout(() => setCopiedFeedback(null), 1500);
    } catch (err) {
      console.warn('Copy failed:', err);
    }
  }, []);

  // Toggle sequence expansion for card view
  const toggleSequenceExpansion = useCallback((idx: number) => {
    setExpandedSequences(prev => {
      const next = new Set(prev);
      if (next.has(idx)) {
        next.delete(idx);
      } else {
        next.add(idx);
      }
      return next;
    });
  }, []);

  // Toggle compare selection
  const toggleCompareSelection = useCallback((idx: number) => {
    setCompareSelection(prev => {
      const next = new Set(prev);
      if (next.has(idx)) {
        next.delete(idx);
      } else {
        next.add(idx);
      }
      return next;
    });
  }, []);

  // Export selected alternatives to CSV
  const handleExportAlternatives = useCallback(() => {
    if (!results?.alternativePrimers) return;
    const toExport = compareSelection.size > 0
      ? results.alternativePrimers.filter((_, i) => compareSelection.has(i))
      : results.alternativePrimers;

    const csvRows = [
      ['Rank', 'Score', 'Fwd Sequence', 'Fwd Tm', 'Fwd Length', 'Rev Sequence', 'Rev Tm', 'Rev Length', 'ΔTm', 'GC Clamps'].join(','),
      ...toExport.map((alt, i) => [
        i + 1,
        alt.compositeScore || '',
        alt.forward.sequence,
        alt.forward.tm?.toFixed(1) || '',
        alt.forward.length || '',
        alt.reverse.sequence,
        alt.reverse.tm?.toFixed(1) || '',
        alt.reverse.length || '',
        Math.abs(alt.forward.tm - alt.reverse.tm).toFixed(1),
        `${(alt.forward.hasGCClamp ? 1 : 0) + (alt.reverse.hasGCClamp ? 1 : 0)}/2`,
      ].join(','))
    ];
    const content = csvRows.join('\n');
    downloadFile(content, 'alternative_primers.csv', 'text/csv');
  }, [results, compareSelection]);

  // Sort alternatives by field
  const sortAlternatives = useCallback((alternatives: AlternateDesign[], sortConfig: AlternatesSort): AlternateDesign[] => {
    return [...alternatives].sort((a, b) => {
      let aVal: number, bVal: number;
      switch (sortConfig.field) {
        case 'score':
          aVal = a.compositeScore || 0;
          bVal = b.compositeScore || 0;
          break;
        case 'tmDiff':
          aVal = Math.abs(a.forward.tm - a.reverse.tm);
          bVal = Math.abs(b.forward.tm - b.reverse.tm);
          break;
        case 'fwdTm':
          aVal = a.forward.tm || 0;
          bVal = b.forward.tm || 0;
          break;
        case 'revTm':
          aVal = a.reverse.tm || 0;
          bVal = b.reverse.tm || 0;
          break;
        case 'fwdLen':
          aVal = a.forward.length || 0;
          bVal = b.forward.length || 0;
          break;
        case 'revLen':
          aVal = a.reverse.length || 0;
          bVal = b.reverse.length || 0;
          break;
        default:
          aVal = a.compositeScore || 0;
          bVal = b.compositeScore || 0;
      }
      const diff = aVal - bVal;
      return sortConfig.direction === 'desc' ? -diff : diff;
    });
  }, []);

  // Filter alternatives
  const filterAlternatives = useCallback((alternatives: AlternateDesign[], filters: AlternatesFilters): AlternateDesign[] => {
    return alternatives.filter(alt => {
      if (filters.minScore > 0 && (alt.compositeScore || 0) < filters.minScore) return false;
      const tmDiff = Math.abs(alt.forward.tm - alt.reverse.tm);
      if (tmDiff > filters.maxTmDiff) return false;
      if (filters.requireGcClamp) {
        const gcClamps = (alt.forward.hasGCClamp ? 1 : 0) + (alt.reverse.hasGCClamp ? 1 : 0);
        if (gcClamps < 2) return false;
      }
      return true;
    });
  }, []);

  // Update AA info when position changes
  const handleAAPositionChange = useCallback((pos: string) => {
    setAaPosition(pos);
    const codonPos = parseInt(pos);
    if (!isNaN(codonPos) && templateSeq) {
      const orfOffset = (parseInt(orfStart) || 1) - 1;
      const nucPos = orfOffset + (codonPos - 1) * 3;
      if (nucPos >= 0 && nucPos + 3 <= templateSeq.length) {
        const codon = templateSeq.slice(nucPos, nucPos + 3);
        const aa = Object.entries(CODON_TABLE).find(([key, codons]) =>
          (codons as any).includes(codon)
        )?.[0] || '?';
        setAaOldAA(aa);
        setSelectionStart(String(nucPos));
        setSelectionEnd(String(nucPos + 3));
      }
    }
  }, [templateSeq, orfStart]);

  return (
    <div className="unified-primer-designer">
      <div className="designer-header">
        <h2>Unified Primer Designer</h2>
        <p className="subtitle">
          One tool for all primer design: amplify, delete, insert, substitute, or mutate
        </p>
      </div>

      {/* Step 1: Template */}
      <section className="design-section">
        <h3>1. Template Sequence</h3>

        <div className="file-input-row">
          <input
            type="file"
            ref={fileInputRef}
            onChange={handleFileUpload}
            accept=".fasta,.fa,.fna,.gb,.gbk,.genbank,.txt"
            className="hidden"
          />
          <button
            type="button"
            className="btn-secondary"
            onClick={() => fileInputRef.current?.click()}
          >
            Upload File
          </button>
          <select
            value={selectedVector}
            onChange={(e) => handleVectorSelect(e.target.value)}
            className="vector-select"
          >
            <option value="">-- Common Vectors --</option>
            {Object.entries(COMMON_VECTORS).map(([key, vec]) => (
              <option key={key} value={key} disabled={!vec.sequence}>
                {vec.name} ({vec.length} bp)
              </option>
            ))}
          </select>
        </div>

        {sequenceInfo && (
          <div className="sequence-info-badge">
            <span className="format-badge">{sequenceInfo.format.toUpperCase()}</span>
            <span className="seq-id">{sequenceInfo.id}</span>
          </div>
        )}

        <textarea
          value={template}
          onChange={(e) => { setTemplate(e.target.value); setSequenceInfo(null); }}
          placeholder="Paste DNA sequence (FASTA, GenBank, or raw)..."
          rows={4}
          spellCheck={false}
        />

        {templateSeq.length > 0 && (
          <div className="template-info">
            {templateSeq.length} bp
            {useReverseComplement && ' (reverse complement)'}
          </div>
        )}

        <label className="checkbox-label">
          <input
            type="checkbox"
            checked={useReverseComplement}
            onChange={(e) => setUseReverseComplement(e.target.checked)}
          />
          Use Reverse Complement
        </label>
      </section>

      {/* Step 2: Select Region */}
      <section className="design-section">
        <h3>2. Select Target Region</h3>

        {templateSeq.length > 0 && (
          <>
            <div className="sequence-bar" onClick={(e: React.MouseEvent<HTMLDivElement>) => {
              const rect = e.currentTarget.getBoundingClientRect();
              const pos = Math.round((e.clientX - rect.left) / rect.width * templateSeq.length);
              if (!selectionStart || selectionEnd) {
                setSelectionStart(String(pos));
                setSelectionEnd('');
              } else {
                // For circular plasmids, allow setting end < start (wrap-around)
                // Just set the end position as clicked, don't auto-swap
                if (designOptions.circular) {
                  setSelectionEnd(String(pos));
                } else {
                  // For linear, swap if needed to ensure start <= end
                  setSelectionEnd(String(Math.max(pos, parseInt(selectionStart))));
                  if (pos < parseInt(selectionStart)) {
                    setSelectionStart(String(pos));
                    setSelectionEnd(selectionStart);
                  }
                }
              }
            }}>
              <div className="sequence-bar-bg" />
              {selectionInfo && !selectionInfo.isWrapped && (
                <div
                  className="sequence-bar-selection"
                  style={{
                    left: `${(selectionInfo.start / templateSeq.length) * 100}%`,
                    width: `${(selectionInfo.length / templateSeq.length) * 100}%`,
                  }}
                />
              )}
              {/* For wrap-around selections, show two bars */}
              {selectionInfo && selectionInfo.isWrapped && (
                <>
                  {/* First part: from start to end of sequence */}
                  <div
                    className="sequence-bar-selection sequence-bar-selection-wrapped"
                    style={{
                      left: `${(selectionInfo.start / templateSeq.length) * 100}%`,
                      width: `${((templateSeq.length - selectionInfo.start) / templateSeq.length) * 100}%`,
                    }}
                  />
                  {/* Second part: from beginning to end position */}
                  <div
                    className="sequence-bar-selection sequence-bar-selection-wrapped"
                    style={{
                      left: '0%',
                      width: `${(selectionInfo.end / templateSeq.length) * 100}%`,
                    }}
                  />
                </>
              )}
              <div className="sequence-bar-scale">
                {[0, 0.25, 0.5, 0.75, 1].map(frac => (
                  <span key={frac} style={{ left: `${frac * 100}%` }}>
                    {Math.round(frac * templateSeq.length)}
                  </span>
                ))}
              </div>
            </div>

            <div className="position-inputs">
              <div className="form-group">
                <label>Start</label>
                <input
                  type="number"
                  value={selectionStart}
                  onChange={(e) => setSelectionStart(e.target.value)}
                  min={0}
                  max={templateSeq.length}
                  placeholder="0"
                />
              </div>
              <div className="form-group">
                <label>End</label>
                <input
                  type="number"
                  value={selectionEnd}
                  onChange={(e) => setSelectionEnd(e.target.value)}
                  min={0}
                  max={templateSeq.length}
                  placeholder={String(templateSeq.length)}
                />
              </div>
              {selectionInfo && (
                <div className="selection-info">
                  {selectionInfo.length} bp selected
                  {selectionInfo.isWrapped && (
                    <span className="wrap-indicator" title="Selection wraps around origin (circular plasmid)">
                      {' '}(wrapping origin)
                    </span>
                  )}
                </div>
              )}
            </div>
          </>
        )}
      </section>

      {/* Step 3: Choose Operation */}
      <section className="design-section">
        <h3>3. What to Do</h3>

        <div className="replacement-modes">
          <label className={`mode-card ${replacementMode === 'amplify' ? 'active' : ''}`}>
            <input
              type="radio"
              name="replacementMode"
              value="amplify"
              checked={replacementMode === 'amplify'}
              onChange={(e) => setReplacementMode(e.target.value)}
            />
            <div className="mode-card-icon amplify">
              <svg width="32" height="32" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round">
                <path d="M4 12h6M14 12h6" />
                <path d="M7 8v8M17 8v8" />
                <path d="M10 12l2-2 2 2" />
                <circle cx="7" cy="8" r="1.5" fill="currentColor" stroke="none" />
                <circle cx="7" cy="16" r="1.5" fill="currentColor" stroke="none" />
                <circle cx="17" cy="8" r="1.5" fill="currentColor" stroke="none" />
                <circle cx="17" cy="16" r="1.5" fill="currentColor" stroke="none" />
              </svg>
            </div>
            <div className="mode-card-content">
              <span className="mode-card-label">Amplify Region</span>
              <span className="mode-card-desc">PCR primers to extract this region</span>
            </div>
          </label>

          <label className={`mode-card ${replacementMode === 'delete' ? 'active' : ''}`}>
            <input
              type="radio"
              name="replacementMode"
              value="delete"
              checked={replacementMode === 'delete'}
              onChange={(e) => setReplacementMode(e.target.value)}
            />
            <div className="mode-card-icon delete">
              <svg width="32" height="32" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round">
                <path d="M4 7h16M10 11v6M14 11v6" />
                <path d="M5 7l1 12a2 2 0 002 2h8a2 2 0 002-2l1-12" />
                <path d="M9 7V4a1 1 0 011-1h4a1 1 0 011 1v3" />
              </svg>
            </div>
            <div className="mode-card-content">
              <span className="mode-card-label">Delete Region</span>
              <span className="mode-card-desc">Remove selected bases</span>
            </div>
          </label>

          <label className={`mode-card ${replacementMode === 'aa' ? 'active' : ''}`}>
            <input
              type="radio"
              name="replacementMode"
              value="aa"
              checked={replacementMode === 'aa'}
              onChange={(e) => setReplacementMode(e.target.value)}
            />
            <div className="mode-card-icon aa">
              <svg width="32" height="32" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round">
                <circle cx="8" cy="8" r="3" />
                <circle cx="16" cy="16" r="3" />
                <path d="M10.5 10.5l3 3" />
                <path d="M16 8l2 2-2 2" />
                <path d="M8 14l-2 2 2 2" />
              </svg>
            </div>
            <div className="mode-card-content">
              <span className="mode-card-label">Change Amino Acid</span>
              <span className="mode-card-desc">Codon-optimized mutation</span>
            </div>
          </label>

          <label className={`mode-card ${replacementMode === 'direct' ? 'active' : ''}`}>
            <input
              type="radio"
              name="replacementMode"
              value="direct"
              checked={replacementMode === 'direct'}
              onChange={(e) => setReplacementMode(e.target.value)}
            />
            <div className="mode-card-icon direct">
              <svg width="32" height="32" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round">
                <path d="M4 6h16M4 12h10M4 18h14" />
                <path d="M18 14l3 3-3 3" />
              </svg>
            </div>
            <div className="mode-card-content">
              <span className="mode-card-label">Direct Sequence</span>
              <span className="mode-card-desc">Insert or replace with nucleotides</span>
            </div>
          </label>
        </div>

        {/* AA Mutation Options */}
        {replacementMode === 'aa' && (
          <div className="aa-options">
            <div className="aa-input-row">
              <div className="form-group">
                <label>ORF Start</label>
                <input
                  type="number"
                  value={orfStart}
                  onChange={(e) => setOrfStart(e.target.value)}
                  min={1}
                  placeholder="1"
                />
              </div>
              <div className="form-group">
                <label>Codon Position</label>
                <input
                  type="number"
                  value={aaPosition}
                  onChange={(e) => handleAAPositionChange(e.target.value)}
                  min={1}
                  placeholder="e.g., 66"
                />
              </div>
              <div className="form-group">
                <label>Current AA</label>
                <input
                  type="text"
                  value={aaOldAA}
                  readOnly
                  className="readonly"
                  placeholder="?"
                />
              </div>
              <span className="arrow">→</span>
              <div className="form-group">
                <label>New AA</label>
                <input
                  type="text"
                  value={aaNewAA}
                  onChange={(e) => setAaNewAA(e.target.value.toUpperCase())}
                  maxLength={1}
                  placeholder="W"
                />
              </div>
              <div className="form-group">
                <label>Organism</label>
                <select value={organism} onChange={(e) => setOrganism(e.target.value)}>
                  <option value="ecoli">E. coli</option>
                  <option value="human">Human</option>
                  <option value="yeast">Yeast</option>
                </select>
              </div>
            </div>
            {aaOldAA && aaNewAA && aaPosition && (
              <div className="aa-preview">
                {AA_NAMES[aaOldAA] || aaOldAA}{aaPosition} → {AA_NAMES[aaNewAA] || aaNewAA}
              </div>
            )}
          </div>
        )}

        {/* Direct Sequence Options */}
        {replacementMode === 'direct' && (
          <div className="direct-options">
            <div className="form-group">
              <label>
                Replacement Sequence
                <span className="label-hint"> (leave empty for deletion, IUPAC codes supported)</span>
              </label>
              <input
                type="text"
                value={directSequence}
                onChange={(e) => setDirectSequence(e.target.value.toUpperCase())}
                placeholder="e.g., ACGT or NNK for library"
              />
            </div>
            {libraryInfo && (
              <div className="library-info">
                Library design: {libraryInfo.size} variants
              </div>
            )}
          </div>
        )}
      </section>

      {/* Actions */}
      <section className="design-actions">
        <button
          type="button"
          className="btn-primary"
          onClick={handleDesign}
          disabled={loading || templateSeq.length < 50}
        >
          {loading ? (
            <span className="btn-loading">
              <span className="btn-spinner"></span>
              {loadingMessage || 'Designing...'}
            </span>
          ) : 'Design Primers'}
        </button>

        <button
          type="button"
          className="btn-secondary"
          onClick={handleAddToBatch}
          disabled={templateSeq.length < 50}
        >
          + Add to Batch
        </button>

        <button
          type="button"
          className="btn-link"
          onClick={() => setShowAdvanced(!showAdvanced)}
        >
          {showAdvanced ? 'Hide' : 'Show'} Advanced Options
        </button>
      </section>

      {/* Advanced Options */}
      {showAdvanced && (
        <section className="advanced-options">
          {/* Row 1: Strategy and Smart Design */}
          <div className="option-row">
            <div className="form-group">
              <label>Strategy</label>
              <select
                value={designOptions.strategy}
                onChange={(e) => setDesignOptions(prev => ({ ...prev, strategy: e.target.value }))}
              >
                <option value="back-to-back">Q5 SDM (Back-to-back)</option>
                <option value="overlapping">QuikChange (Overlapping)</option>
              </select>
            </div>
            <label className="checkbox-label">
              <input
                type="checkbox"
                checked={designOptions.useSmartDesign}
                onChange={(e) => setDesignOptions(prev => ({ ...prev, useSmartDesign: e.target.checked }))}
              />
              <span className="checkbox-text">
                Smart Design
                <span className="option-hint-inline"> - Optimize for GC clamp and 3' stability</span>
              </span>
            </label>
            <label className="checkbox-label">
              <input
                type="checkbox"
                checked={designOptions.circular}
                onChange={(e) => setDesignOptions(prev => ({ ...prev, circular: e.target.checked }))}
              />
              <span className="checkbox-text">
                Circular Plasmid
                <span className="option-hint-inline"> - Allow primers to wrap around origin</span>
              </span>
            </label>
            <label className="checkbox-label">
              <input
                type="checkbox"
                checked={designOptions.exhaustiveSearch}
                onChange={(e) => setDesignOptions(prev => ({ ...prev, exhaustiveSearch: e.target.checked }))}
              />
              <span className="checkbox-text">
                Exhaustive Search
                <span className="option-hint-inline"> - Find best design (shows preview first)</span>
              </span>
            </label>
          </div>

          {/* Row 1.5: 5' Tail Confinement (mutagenesis-specific) */}
          <div className="option-row">
            <label className="checkbox-label" title="When enabled, mutations are placed only in the 5' non-annealing overhang. The 3' portion binds perfectly to the original template, improving PCR efficiency.">
              <input
                type="checkbox"
                checked={designOptions.confineTo5Tails}
                onChange={(e) => setDesignOptions(prev => ({ ...prev, confineTo5Tails: e.target.checked }))}
              />
              <span className="checkbox-text">
                Confine Mutations to 5' Overhang
              </span>
            </label>
            <div className="option-explanation">
              <details>
                <summary>When to use this option?</summary>
                <div className="explanation-content">
                  <p><strong>5' Overhang Design:</strong> The mutation is placed in the non-annealing 5' tail, while the 3' region binds perfectly to the original template.</p>
                  <pre className="primer-diagram">
{`Standard:     5'─[flank]─[MUTATION]─[annealing region]─3'
                         ↑ mutation may be in annealing region

5' Confined:  5'─[MUTATION in tail]─[perfect 3' annealing]─3'
                                    ↑ 100% template match`}
                  </pre>
                  <p><strong>Benefits:</strong></p>
                  <ul>
                    <li>Higher PCR efficiency - 3' end binds perfectly to template</li>
                    <li>Better specificity - no mismatches in critical annealing region</li>
                    <li>Recommended for difficult templates (high GC, secondary structures)</li>
                  </ul>
                  <p><strong>Trade-offs:</strong></p>
                  <ul>
                    <li>Longer primers (mutation + full annealing length)</li>
                    <li>Higher synthesis cost</li>
                    <li>May not work for large insertions (&gt;30bp)</li>
                  </ul>
                </div>
              </details>
            </div>
          </div>

          {/* Row 2: Tm Parameters */}
          <div className="option-row primer-params">
            <div className="param-group">
              <label className="param-group-label">Tm Parameters (°C)</label>
              <div className="param-inputs">
                <div className="form-group compact">
                  <label>Target</label>
                  <input
                    type="number"
                    value={designOptions.optimalTm}
                    onChange={(e) => setDesignOptions(prev => ({ ...prev, optimalTm: parseFloat(e.target.value) || 62 }))}
                    min={45}
                    max={75}
                  />
                </div>
                <div className="form-group compact">
                  <label>Min</label>
                  <input
                    type="number"
                    value={designOptions.minTm}
                    onChange={(e) => setDesignOptions(prev => ({ ...prev, minTm: parseFloat(e.target.value) || 55 }))}
                    min={40}
                    max={65}
                  />
                </div>
                <div className="form-group compact">
                  <label>Max</label>
                  <input
                    type="number"
                    value={designOptions.maxTm}
                    onChange={(e) => setDesignOptions(prev => ({ ...prev, maxTm: parseFloat(e.target.value) || 72 }))}
                    min={55}
                    max={80}
                  />
                </div>
              </div>
            </div>

            {/* Primer Length */}
            <div className="param-group">
              <label className="param-group-label">Primer Length (bp)</label>
              <div className="param-inputs">
                <div className="form-group compact">
                  <label>Min</label>
                  <input
                    type="number"
                    value={designOptions.minPrimerLength}
                    onChange={(e) => setDesignOptions(prev => ({ ...prev, minPrimerLength: parseInt(e.target.value) || 15 }))}
                    min={15}
                    max={40}
                  />
                </div>
                <div className="form-group compact">
                  <label>Max</label>
                  <input
                    type="number"
                    value={designOptions.maxPrimerLength}
                    onChange={(e) => setDesignOptions(prev => ({ ...prev, maxPrimerLength: parseInt(e.target.value) || 60 }))}
                    min={25}
                    max={60}
                  />
                </div>
              </div>
            </div>

            {/* GC Content */}
            <div className="param-group">
              <label className="param-group-label">GC Content (%)</label>
              <div className="param-inputs">
                <div className="form-group compact">
                  <label>Min</label>
                  <input
                    type="number"
                    value={designOptions.minGC}
                    onChange={(e) => setDesignOptions(prev => ({ ...prev, minGC: parseInt(e.target.value) || 40 }))}
                    min={20}
                    max={50}
                  />
                </div>
                <div className="form-group compact">
                  <label>Max</label>
                  <input
                    type="number"
                    value={designOptions.maxGC}
                    onChange={(e) => setDesignOptions(prev => ({ ...prev, maxGC: parseInt(e.target.value) || 60 }))}
                    min={50}
                    max={80}
                  />
                </div>
              </div>
            </div>
          </div>
        </section>
      )}

      {/* Batch Panel */}
      {showBatch && batchItems.length > 0 && (
        <section className="batch-panel">
          <h3>Batch Queue ({batchItems.length} items)</h3>
          <div className="batch-list">
            {batchItems.map((item, idx) => (
              <div key={item.id} className="batch-item">
                <span className="batch-num">{idx + 1}</span>
                <span className="batch-label">{item.label}</span>
                <button
                  type="button"
                  className="btn-icon"
                  onClick={() => handleRemoveFromBatch(item.id)}
                >
                  ✕
                </button>
              </div>
            ))}
          </div>
          <div className="batch-actions">
            <button
              type="button"
              className="btn-primary"
              onClick={handleDesignBatch}
              disabled={loading}
            >
              {loading ? 'Processing...' : 'Design All Primers'}
            </button>
            <button
              type="button"
              className="btn-secondary"
              onClick={() => setBatchItems([])}
            >
              Clear Batch
            </button>
          </div>
        </section>
      )}

      {/* Error Display */}
      {error && (
        <div className="error-message">
          <span className="error-icon">⚠</span>
          {error}
        </div>
      )}

      {/* Single Result Display */}
      {results && !results.isLibrary && (
        <section className="results-section" ref={resultsRef}>
              {/* Preview indicator when still optimizing */}
              {isPreviewResult && (
                <div className="preview-banner">
                  <span className="preview-icon">⏳</span>
                  <span className="preview-text">Preview result shown while finding optimal design...</span>
                </div>
              )}
              {/* Quick Status Summary Card (matches legacy MutagenesisDesigner) */}
              {summaryStatus && (
                <div id="section-summary" className={`status-summary-card status-${summaryStatus.overall}`}>
                  <button
                    type="button"
                    className="section-header collapsible"
                    onClick={() => toggleSection('summary')}
                  >
                    <div className="section-header-content">
                      <span className={`overall-status-badge status-${summaryStatus.overall}`}>
                        {summaryStatus.overall === 'pass' ? 'All Checks Passed' :
                         summaryStatus.overall === 'warning' ? 'Some Warnings' : 'Issues Found'}
                      </span>
                      <span className="section-title">Quick Status</span>
                    </div>
                    <span className="collapse-icon">{collapsedSections.summary ? '▶' : '▼'}</span>
                  </button>
                  {!collapsedSections.summary && (
                    <div className="status-checks-grid">
                      {Object.entries(summaryStatus.checks).map(([key, check]) => (
                        <button
                          key={key}
                          type="button"
                          className={`status-check status-${check.status} clickable`}
                          onClick={() => check.section && scrollToSection(check.section)}
                          title={check.threshold}
                        >
                          <span className="check-icon">
                            {check.status === 'pass' ? '✓' : check.status === 'warning' ? '⚠' : '✗'}
                          </span>
                          <div className="check-info">
                            <span className="check-label">{check.label}</span>
                            <span className="check-value">{check.value}</span>
                          </div>
                        </button>
                      ))}
                    </div>
                  )}
                </div>
              )}

              {/* Mutation Summary Section - description only, badges moved to sequence comparison */}
              {results.description && (
                <div className="mutation-summary">
                  <div className="summary-card">
                    <div className="mutation-description">
                      {results.description}
                    </div>
                  </div>
                </div>
              )}

              {/* Designed Primers Section */}
              <div id="section-primers" className="primers-section collapsible-section">
                <button
                  type="button"
                  className="section-header collapsible"
                  onClick={() => toggleSection('primers')}
                >
                  <h3>Designed Primers</h3>
                  <div className="section-header-actions">
                    <div className="download-buttons" onClick={(e) => e.stopPropagation()}>
                      <button type="button" className="btn-download" onClick={handleDownloadPrimers}>TXT</button>
                      <button type="button" className="btn-download" onClick={handleDownloadFasta}>FASTA</button>
                      <button type="button" className="btn-download" onClick={handleDownloadSingleCSV}>CSV</button>
                    </div>
                    <span className="collapse-icon">{collapsedSections.primers ? '▶' : '▼'}</span>
                  </div>
                </button>
                {!collapsedSections.primers && (
                  <>
                    {/* Quality Score Ring - Clickable for Score Breakdown */}
                    <button
                      type="button"
                      className="quality-score-container bg-transparent border-0 p-0 cursor-pointer w-full text-left transition-transform duration-150 hover:scale-[1.02]"
                      onClick={() => setShowScoreBreakdown(true)}
                      title="Click to see detailed score breakdown"
                    >
                      <div className={`quality-ring quality-${results.quality || results.qualityTier || 'good'}`}>
                        <svg viewBox="0 0 36 36" className="quality-svg">
                          <path
                            className="quality-bg"
                            d="M18 2.0845 a 15.9155 15.9155 0 0 1 0 31.831 a 15.9155 15.9155 0 0 1 0 -31.831"
                          />
                          <path
                            className="quality-fill"
                            strokeDasharray={`${results.effectiveScore ?? results.compositeScore ?? 75}, 100`}
                            d="M18 2.0845 a 15.9155 15.9155 0 0 1 0 31.831 a 15.9155 15.9155 0 0 1 0 -31.831"
                          />
                        </svg>
                        <div className="quality-text">
                          <span className="quality-score">{results.effectiveScore ?? results.compositeScore ?? ''}</span>
                        </div>
                      </div>
                      <div className="quality-details">
                        <span className="quality-title">
                          Quality Score: <span className={`quality-tier tier-${results.quality || results.qualityTier || 'good'}`}>{results.quality || results.qualityTier || 'good'}</span>
                          {(results.criticalWarnings ?? 0) > 0 && (
                            <span className="text-[11px] text-red-500 ml-1.5">
                              ({results.criticalWarnings} critical)
                            </span>
                          )}
                        </span>
                        <span className="quality-desc" title={(results.quality || results.qualityTier) === 'excellent' ? 'High confidence - expected >90% PCR success' :
                           (results.quality || results.qualityTier) === 'good' ? 'Reliable - expected ~80% PCR success' :
                           (results.quality || results.qualityTier) === 'acceptable' ? 'May work - consider alternatives' :
                           'High failure risk - redesign recommended'}>
                          {(results.quality || results.qualityTier) === 'excellent' ? 'High confidence - expected >90% PCR success' :
                           (results.quality || results.qualityTier) === 'good' ? 'Reliable - expected ~80% PCR success' :
                           (results.quality || results.qualityTier) === 'acceptable' ? 'May work - consider alternatives' :
                           'High failure risk - redesign recommended'}
                        </span>
                        <span className="text-xs opacity-80">• Click for details</span>
                      </div>
                    </button>

                    {/* Score Breakdown Popup */}
                    {showScoreBreakdown && (
                      <ScoreBreakdownPopup
                        compositeScore={results.effectiveScore ?? results.compositeScore ?? 0}
                        rawScore={results.compositeScore}
                        criticalWarnings={results.criticalWarnings || 0}
                        quality={{ tier: results.quality || results.qualityTier, label: results.quality || results.qualityTier }}
                        forwardScores={results.forwardPiecewiseScores}
                        reverseScores={results.reversePiecewiseScores}
                        hasTemplate={!!template}
                        onClose={() => setShowScoreBreakdown(false)}
                      />
                    )}

                    <div className="primer-pair">
                      {/* Forward Primer */}
                      <div className="primer-card">
                        <div className="primer-header">
                          <span className="direction-badge forward">
                            {results.type === 'amplification' ? 'Forward' : 'Forward (Mutagenic)'}
                          </span>
                        </div>
                        <div className="primer-sequence-display">
                          <code>{results.forward.sequence}</code>
                          <button
                            type="button"
                            className="btn-copy-inline"
                            onClick={() => copyToClipboard(results.forward.sequence)}
                            title="Copy sequence"
                          >
                            Copy
                          </button>
                        </div>
                        <div className="primer-stats">
                          <div className="stat">
                            <span className="label">Length</span>
                            <span className="value">{results.forward.sequence?.length} bp</span>
                          </div>
                          <div className="stat">
                            <span className="label">Tm</span>
                            <span className="value">{results.forward.tm?.toFixed(1)}°C</span>
                          </div>
                          <div className="stat">
                            <span className="label">GC</span>
                            <span className="value">{results.forward.gcPercent || `${((results.forward.gc || 0) * 100).toFixed(1)}%`}</span>
                          </div>
                          <div className="stat">
                            <span className="label">GC Clamp</span>
                            <span className="value">{results.forward.hasGCClamp ? 'Yes' : 'No'}</span>
                          </div>
                        </div>
                      </div>

                      {/* Reverse Primer */}
                      <div className="primer-card">
                        <div className="primer-header">
                          <span className="direction-badge reverse">Reverse</span>
                        </div>
                        <div className="primer-sequence-display">
                          <code>{results.reverse.sequence}</code>
                          <button
                            type="button"
                            className="btn-copy-inline"
                            onClick={() => copyToClipboard(results.reverse.sequence)}
                            title="Copy sequence"
                          >
                            Copy
                          </button>
                        </div>
                        <div className="primer-stats">
                          <div className="stat">
                            <span className="label">Length</span>
                            <span className="value">{results.reverse.sequence?.length} bp</span>
                          </div>
                          <div className="stat">
                            <span className="label">Tm</span>
                            <span className="value">{results.reverse.tm?.toFixed(1)}°C</span>
                          </div>
                          <div className="stat">
                            <span className="label">GC</span>
                            <span className="value">{results.reverse.gcPercent || `${((results.reverse.gc || 0) * 100).toFixed(1)}%`}</span>
                          </div>
                          <div className="stat">
                            <span className="label">GC Clamp</span>
                            <span className="value">{results.reverse.hasGCClamp ? 'Yes' : 'No'}</span>
                          </div>
                        </div>
                      </div>
                    </div>

                    {/* PCR Parameters Summary - BELOW primer cards */}
                    {results.forward && results.reverse && (
                      <div className="pcr-params-inline">
                        {(() => {
                          const templateLen = templateSeq?.length || 3000;
                          const templateKb = templateLen / 1000;
                          const extensionSec = Math.max(30, Math.ceil(templateKb * 25));
                          const lowerTm = Math.min(results.forward.tm, results.reverse.tm);
                          const annealingTemp = Math.round(Math.min(lowerTm + 1, 72));
                          const extensionTime = extensionSec >= 60
                            ? `${Math.floor(extensionSec / 60)}:${(extensionSec % 60).toString().padStart(2, '0')}`
                            : `${extensionSec}s`;
                          return (
                            <>
                              <div className="pcr-param">
                                <span className="param-label">Annealing Temp</span>
                                <span className="param-value">{annealingTemp}°C</span>
                                <span className="param-note">Q5: lower Tm + 1°C</span>
                              </div>
                              <div className="pcr-param">
                                <span className="param-label">Extension Time</span>
                                <span className="param-value">{extensionTime}</span>
                                <span className="param-note">~25s/kb for {templateKb.toFixed(1)} kb template</span>
                              </div>
                            </>
                          );
                        })()}
                      </div>
                    )}

                    {/* Sequence Comparison - for all modes */}
                    {/* For mutation modes: show original vs mutated */}
                    {/* For amplify mode: show the amplified region */}
                    {(results.originalSequence && results.mutatedSequence) || results.type === 'amplification' ? (
                      <div id="section-comparison" className="sequence-comparison-inline">
                        {/* Amplification mode - show the amplified region in summary format */}
                        {results.type === 'amplification' && templateSeq && (
                          <div className="amplification-summary">
                            {(() => {
                              const start = results.forward?.start ?? (parseInt(selectionStart) || 0);
                              const end = results.reverse?.end ?? (parseInt(selectionEnd) || templateSeq.length);
                              const ampLen = end - start;
                              const ampSeq = templateSeq.slice(start, end);
                              const displaySeq = ampSeq.length > 50
                                ? `${ampSeq.slice(0, 25)}...${ampSeq.slice(-25)}`
                                : ampSeq;
                              return (
                                <>
                                  <div className="summary-label-row">
                                    <span className="amplified-sequence-label">Amplified region ({ampLen} bp, position {start + 1}-{end}):</span>
                                    <div className="summary-badges">
                                      <div className="mutation-type-badge">
                                        {(results.type || 'amplification').replace('_', ' ').toUpperCase()}
                                      </div>
                                      <div className="design-badge" data-design={results.design || 'pcr'}>
                                        {results.design === 'pcr' ? 'PCR' : results.design === 'nebuilder' ? 'NEBuilder' : 'PCR'}
                                      </div>
                                    </div>
                                  </div>
                                  <code className="amplified-sequence">{displaySeq}</code>
                                </>
                              );
                            })()}
                          </div>
                        )}

                        {/* Mutation modes - show original vs mutated */}
                        {results.originalSequence && results.mutatedSequence && (
                          <div className="comparison-display">
                            <div className="seq-row">
                              <span className="seq-label">Original</span>
                              <code className="seq-snippet">
                                {(() => {
                                  const pos = results.position ?? 0;
                                  const delLen = results.deleteLength || results.deletionLength || 1;
                                  const start = Math.max(0, pos - 15);
                                  const end = Math.min(results.originalSequence?.length || 0, pos + delLen + 15);
                                  const prefix = results.originalSequence?.slice(start, pos) || '';
                                  const changed = results.originalSequence?.slice(pos, pos + delLen) || '';
                                  const suffix = results.originalSequence?.slice(pos + delLen, end) || '';
                                  const displayChanged = changed.length > 20
                                    ? `${changed.slice(0, 10)}...${changed.slice(-10)}`
                                    : changed;
                                  return (
                                    <>
                                      ...{prefix}<span className="highlight-deleted">{displayChanged}</span>{suffix}...
                                    </>
                                  );
                                })()}
                              </code>
                            </div>
                            <div className="seq-row">
                              <span className="seq-label">Mutated</span>
                              <code className="seq-snippet mutated">
                                {(() => {
                                  const pos = results.position ?? 0;
                                  const insLen = (results.insertSequence || results.replacement || '').length || 1;
                                  const start = Math.max(0, pos - 15);
                                  const end = Math.min(results.mutatedSequence?.length || 0, pos + insLen + 15);
                                  const prefix = results.mutatedSequence?.slice(start, pos) || '';
                                  const suffix = results.mutatedSequence?.slice(pos + insLen, end) || '';
                                  if (results.type === 'deletion') {
                                    const delSuffix = results.mutatedSequence?.slice(pos, pos + 15) || '';
                                    return (
                                      <>
                                        ...{prefix}<span className="highlight-junction">|</span>{delSuffix}...
                                      </>
                                    );
                                  } else {
                                    const changed = results.insertSequence || results.replacement || results.newCodon || results.mutatedSequence?.slice(pos, pos + insLen) || '';
                                    const displayChanged = changed.length > 20
                                      ? `${changed.slice(0, 10)}...${changed.slice(-10)}`
                                      : changed;
                                    return (
                                      <>
                                        ...{prefix}<span className="highlight-inserted">{displayChanged}</span>{suffix}...
                                      </>
                                    );
                                  }
                                })()}
                              </code>
                            </div>
                          </div>
                        )}

                        {/* Deletion summary with badges */}
                        {results.type === 'deletion' && results.deletedSequence && (
                          <div className="deletion-summary">
                            <div className="summary-label-row">
                              <span className="deleted-sequence-label">Deleted sequence ({results.deleteLength || results.deletionLength} bp):</span>
                              <div className="summary-badges">
                                <div className="mutation-type-badge">DELETION</div>
                                <div className="design-badge" data-design={results.design || results.strategy}>
                                  {results.design === 'back-to-back' || results.strategy === 'back-to-back' ? 'Q5 SDM' : 'QuikChange'}
                                </div>
                                {results.position !== undefined && (
                                  <span className="position-info">Position: {results.position + 1}</span>
                                )}
                              </div>
                            </div>
                            <code className="deleted-sequence">
                              {results.deletedSequence.length > 50
                                ? `${results.deletedSequence.slice(0, 25)}...${results.deletedSequence.slice(-25)}`
                                : results.deletedSequence}
                            </code>
                          </div>
                        )}

                        {/* Substitution/Insertion summary with badges */}
                        {results.type === 'substitution' && results.insertSequence && (
                          <div className="substitution-summary">
                            <div className="summary-label-row">
                              <span className="substitution-sequence-label">Inserted sequence ({results.insertSequence.length} bp):</span>
                              <div className="summary-badges">
                                <div className="mutation-type-badge">SUBSTITUTION</div>
                                <div className="design-badge" data-design={results.design || results.strategy}>
                                  {results.design === 'back-to-back' || results.strategy === 'back-to-back' ? 'Q5 SDM' : 'QuikChange'}
                                </div>
                                {results.position !== undefined && (
                                  <span className="position-info">Position: {results.position + 1}</span>
                                )}
                              </div>
                            </div>
                            <code className="substitution-sequence">
                              {results.insertSequence.length > 50
                                ? `${results.insertSequence.slice(0, 25)}...${results.insertSequence.slice(-25)}`
                                : results.insertSequence}
                            </code>
                          </div>
                        )}

                        {/* Codon change details with badges */}
                        {results.type === 'codon_change' && (
                          <div className="codon-change-summary">
                            <div className="summary-label-row">
                              <span className="codon-info-item">
                                <span className="info-label">Codon:</span>
                                <code className="codon-value">{results.oldCodon}</code>
                                <span className="arrow">→</span>
                                <code className="codon-value">{results.newCodon}</code>
                                {results.codonUsage && (
                                  <span className="codon-usage">({(results.codonUsage * 100).toFixed(0)}% usage)</span>
                                )}
                              </span>
                              <div className="summary-badges">
                                <div className="mutation-type-badge">CODON CHANGE</div>
                                <div className="design-badge" data-design={results.design || results.strategy}>
                                  {results.design === 'back-to-back' || results.strategy === 'back-to-back' ? 'Q5 SDM' : 'QuikChange'}
                                </div>
                                {results.position !== undefined && (
                                  <span className="position-info">Position: {results.position + 1}</span>
                                )}
                              </div>
                            </div>
                          </div>
                        )}
                      </div>
                    ) : null}

                  </>
                )}
              </div>

              {/* Enhanced Analysis Section - Using shared component */}
              {enhancedAnalysis && (
                <EnhancedAnalysisSection
                  analysis={enhancedAnalysis as any}
                  collapsible={true}
                  defaultCollapsed={collapsedSections.analysis}
                />
              )}

              {/* Primer Visualization Section */}
              {templateSeq && results.forward && results.reverse && (
                <div id="section-visualization" className="primer-visualization-section collapsible-section">
                  <button
                    type="button"
                    className="section-header collapsible"
                    onClick={() => toggleSection('visualization')}
                  >
                    <h3>Primer Visualization</h3>
                    <span className="collapse-icon">{collapsedSections.visualization ? '▶' : '▼'}</span>
                  </button>
                  {!collapsedSections.visualization && (
                    <div className="section-content">
                      <PrimerOnTemplateViewer
                        template={templateSeq}
                        forward={results.forward}
                        reverse={results.reverse}
                        mutationPosition={results.nucleotidePosition ?? results.position}
                        mutationLength={1}
                        showHairpinDiagrams={true}
                        showSequenceDetails={true}
                        width={800}
                        isSDMMode={results.design === 'back-to-back' || results.strategy === 'back-to-back'}
                      />
                    </div>
                  )}
                </div>
              )}

              {/* Protocol Section - shows for all modes */}
              {results.forward && results.reverse && (
                <div id="section-protocol" className="protocol-section collapsible-section">
                  <button
                    type="button"
                    className="section-header collapsible"
                    onClick={() => toggleSection('protocol')}
                  >
                    <h3>PCR Protocol</h3>
                    <span className="collapse-icon">{collapsedSections.protocol ? '▶' : '▼'}</span>
                  </button>
                  {!collapsedSections.protocol && (
                    <div className="section-content">
                      {/* Key Parameters Highlight (matches legacy) */}
                      {(() => {
                        const templateLen = templateSeq?.length || 3000;
                        const templateKb = templateLen / 1000;
                        // Use 20-30s per kb for Q5 polymerase
                        const extensionSec = Math.max(30, Math.ceil(templateKb * 25));
                        // Q5 annealing temp: min(lower_Tm + 1, 72) - NEB recommendation
                        const lowerTm = Math.min(results.forward.tm, results.reverse.tm);
                        const annealingTemp = Math.round(Math.min(lowerTm + 1, 72));
                        const extensionTime = extensionSec >= 60
                          ? `${Math.floor(extensionSec / 60)}:${(extensionSec % 60).toString().padStart(2, '0')}`
                          : `${extensionSec}s`;

                        return (
                          <>
                            {results.protocol ? (
                              <div className="protocol-card">
                                <h4>{results.protocol.name}</h4>
                                <details open className="protocol-details">
                                  <summary className="protocol-summary">Detailed Steps</summary>
                                  <div className="protocol-steps">
                                    {results.protocol.steps?.map((step, idx) => (
                                      <div key={idx} className="protocol-step">
                                        <span className="step-name">{step.name}</span>
                                        {step.substeps ? (
                                          <div className="substeps">
                                            {step.substeps.map((sub, sidx) => {
                                              const isAnneal = sub.name.toLowerCase().includes('anneal');
                                              const isExtend = sub.name.toLowerCase().includes('extend');
                                              // Replace generic times with calculated values
                                              let displayTemp = sub.temp;
                                              let displayTime = sub.time;
                                              if (isAnneal) {
                                                displayTemp = `${annealingTemp}°C`;
                                              }
                                              if (isExtend && sub.time.includes('/kb')) {
                                                displayTime = extensionTime;
                                              }
                                              return (
                                                <span key={sidx} className={`substep ${isAnneal ? 'highlight-anneal' : ''} ${isExtend ? 'highlight-extend' : ''}`}>
                                                  {sub.name}: {displayTemp} × {displayTime}
                                                </span>
                                              );
                                            })}
                                          </div>
                                        ) : (
                                          <span className="step-params">{step.temp} × {step.time}</span>
                                        )}
                                      </div>
                                    ))}
                                  </div>
                                </details>
                                {results.protocol.notes && (
                                  <div className="protocol-notes">
                                    <h5>Notes</h5>
                                    <ul>
                                      {results.protocol.notes.map((note, idx) => (
                                        <li key={idx}>{note}</li>
                                      ))}
                                    </ul>
                                  </div>
                                )}
                              </div>
                            ) : (
                              <div className="protocol-card">
                                <h4>Standard PCR Protocol</h4>
                                <details open className="protocol-details">
                                  <summary className="protocol-summary">Detailed Steps</summary>
                                  <div className="protocol-steps">
                                    <div className="protocol-step">
                                      <span className="step-name">Initial Denaturation</span>
                                      <span className="step-params">98°C × 30s</span>
                                    </div>
                                    <div className="protocol-step">
                                      <span className="step-name">Cycling (25-35×)</span>
                                      <div className="substeps">
                                        <span className="substep">Denature: 98°C × 10s</span>
                                        <span className="substep highlight-anneal">Anneal: {annealingTemp}°C × 20s</span>
                                        <span className="substep highlight-extend">Extend: 72°C × {extensionTime}</span>
                                      </div>
                                    </div>
                                    <div className="protocol-step">
                                      <span className="step-name">Final Extension</span>
                                      <span className="step-params">72°C × 2:00</span>
                                    </div>
                                    <div className="protocol-step">
                                      <span className="step-name">Hold</span>
                                      <span className="step-params">4°C × ∞</span>
                                    </div>
                                  </div>
                                </details>
                              </div>
                            )}
                          </>
                        );
                      })()}
                    </div>
                  )}
                </div>
              )}

              {/* Alternate Designs for Mutations (using unified AlternativesPanel) */}
              {results.design === 'back-to-back' && results.alternateDesigns && results.alternateDesigns.length > 0 && (
                <div id="section-alternates" className="alternate-designs collapsible-section">
                  <button
                    type="button"
                    className="section-header collapsible"
                    onClick={() => toggleSection('alternates')}
                  >
                    <h3>Alternative Designs (QuikChange-Style)</h3>
                    <span className="collapse-icon">{collapsedSections.alternates ? '▶' : '▼'}</span>
                  </button>
                  {!collapsedSections.alternates && (
                    <AlternativesPanel
                      alternatives={results.alternateDesigns as any}
                      currentDesign={{
                        forward: results.forward,
                        reverse: results.reverse,
                        compositeScore: results.compositeScore,
                        score: results.compositeScore,
                      }}
                      originalDesign={originalDesign as any}
                      isAlternativeSelected={results.isAlternativeSelected || false}
                      selectedIndex={results.selectedAlternateIdx || -1}
                      onSelectAlternative={handleSelectAlternate as any}
                      onRevert={handleRevertDesign}
                      features={{
                        viewToggle: true,
                        filters: false,  // Fewer alternatives, filters not needed
                        sorting: true,
                        exportCSV: true,
                        compare: true,
                        expandSequences: true,
                        badges: true,
                        showCurrentDesign: true,
                      }}
                    />
                  )}
                </div>
              )}

              {/* Alternative Primers for PCR Amplification (using unified AlternativesPanel) */}
              {results.type === 'amplification' && results.alternativePrimers && results.alternativePrimers.length > 0 && (
                <div id="section-alternates" className="alternate-designs collapsible-section">
                  <button
                    type="button"
                    className="section-header collapsible"
                    onClick={() => toggleSection('alternates')}
                  >
                    <h3>Alternative Primer Pairs</h3>
                    <span className="collapse-icon">{collapsedSections.alternates ? '▶' : '▼'}</span>
                  </button>
                  {!collapsedSections.alternates && (
                    <AlternativesPanel
                      alternatives={results.alternativePrimers as any}
                      alternativeCategories={results.alternativeCategories as any}
                      currentDesign={{
                        forward: results.forward,
                        reverse: results.reverse,
                        compositeScore: results.compositeScore,
                        score: results.compositeScore,
                      }}
                      originalDesign={originalDesign as any}
                      isAlternativeSelected={results.isAlternativeSelected || false}
                      selectedIndex={results.selectedAlternateIdx || -1}
                      onSelectAlternative={(alt) => handleSelectAlternate(alt as any)}
                      onRevert={handleRevertDesign}
                      wasUpgraded={results.wasUpgraded}
                      originalScore={results.originalScore}
                      features={{
                        viewToggle: true,
                        filters: true,
                        sorting: true,
                        exportCSV: true,
                        compare: true,
                        expandSequences: true,
                        badges: true,
                        showCurrentDesign: true,
                      }}
                    />
                  )}
                </div>
              )}
        </section>
      )}


      {/* Batch Results Display */}
      {batchResults.length > 0 && (
        <section className="batch-results-section">
          <div className="batch-results-header">
            <h3>
              Batch Results ({batchResults.filter(r => r.success).length}/{batchResults.length} successful)
            </h3>
            <button
              type="button"
              className="btn-secondary"
              onClick={handleDownloadCSV}
            >
              Download CSV
            </button>
          </div>

          <div className="batch-results-list">
            {batchResults.map((result, idx) => (
              <div key={idx} className={`batch-result-item ${result.success ? '' : 'error'}`}>
                {result.success ? (
                  <>
                    <span className="batch-label">{result.originalSpec?.label || result.description}</span>
                    <div className="batch-primers">
                      <code title="Forward">{result.forward!.sequence}</code>
                      <code title="Reverse">{result.reverse!.sequence}</code>
                    </div>
                  </>
                ) : (
                  <span className="batch-error">{result.originalSpec?.label}: {result.error}</span>
                )}
              </div>
            ))}
          </div>
        </section>
      )}

      {/* Comparison Modal */}
      {showCompareModal && results?.alternativePrimers && (
        <div className="compare-modal-overlay" onClick={() => setShowCompareModal(false)}>
          <div className="compare-modal" onClick={(e) => e.stopPropagation()}>
            <div className="compare-modal-header">
              <h3>Primer Comparison</h3>
              <button
                type="button"
                className="compare-modal-close"
                onClick={() => setShowCompareModal(false)}
              >
                ✕
              </button>
            </div>
            <div className="compare-modal-body">
              <div className="compare-grid">
                {(() => {
                  const selectedPrimers = results.alternativePrimers
                    .map((alt, idx) => ({ ...alt, originalIdx: idx }))
                    .filter((_, idx) => compareSelection.has(idx));
                  const orig = originalDesign || results;

                  return selectedPrimers.map((alt, colIdx) => {
                    const tmDiff = Math.abs(alt.forward.tm - alt.reverse.tm);
                    const gcClamps = (alt.forward.hasGCClamp ? 1 : 0) + (alt.reverse.hasGCClamp ? 1 : 0);
                    const scoreDelta = (alt.compositeScore || 0) - (orig.compositeScore || 0);
                    const fwdGC = alt.forward.gcContent || ((alt.forward.sequence.match(/[GC]/g) || []).length / alt.forward.sequence.length * 100);
                    const revGC = alt.reverse.gcContent || ((alt.reverse.sequence.match(/[GC]/g) || []).length / alt.reverse.sequence.length * 100);

                    return (
                      <div key={alt.originalIdx} className="compare-column">
                        <div className="compare-column-header">
                          <span className="compare-rank">#{colIdx + 1}</span>
                          <span className={`compare-score ${alt.compositeScore! >= 80 ? 'excellent' : alt.compositeScore! >= 70 ? 'good' : alt.compositeScore! >= 60 ? 'acceptable' : 'poor'}`}>
                            Score: {alt.compositeScore}
                            {scoreDelta !== 0 && (
                              <span className={`score-delta ${scoreDelta > 0 ? 'better' : 'worse'}`}>
                                ({scoreDelta > 0 ? '+' : ''}{scoreDelta})
                              </span>
                            )}
                          </span>
                        </div>

                        <div className="compare-section">
                          <h4>Forward Primer</h4>
                          <div className="compare-sequence">
                            <code>{alt.forward.sequence}</code>
                            <button
                              type="button"
                              className="copy-btn-mini"
                              onClick={() => handleCopySequence(alt.forward.sequence)}
                              title="Copy sequence"
                            >
                              📋
                            </button>
                          </div>
                          <div className="compare-metrics">
                            <div className="metric-row">
                              <span className="metric-label">Length:</span>
                              <span className="metric-value">{alt.forward.length} bp</span>
                            </div>
                            <div className="metric-row">
                              <span className="metric-label">Tm:</span>
                              <span className="metric-value">{alt.forward.tm?.toFixed(1)}°C</span>
                            </div>
                            <div className="metric-row">
                              <span className="metric-label">GC%:</span>
                              <span className="metric-value">{fwdGC.toFixed(1)}%</span>
                            </div>
                            <div className="metric-row">
                              <span className="metric-label">GC Clamp:</span>
                              <span className={`metric-value ${alt.forward.hasGCClamp ? 'good' : 'warn'}`}>
                                {alt.forward.hasGCClamp ? 'Yes' : 'No'}
                              </span>
                            </div>
                            {alt.forward.terminal3DG !== undefined && (
                              <div className="metric-row">
                                <span className="metric-label">3' ΔG:</span>
                                <span className="metric-value">{alt.forward.terminal3DG?.toFixed(1)} kcal/mol</span>
                              </div>
                            )}
                            {alt.forward.hairpinDG !== undefined && (
                              <div className="metric-row">
                                <span className="metric-label">Hairpin ΔG:</span>
                                <span className={`metric-value ${alt.forward.hairpinDG > -3 ? 'good' : alt.forward.hairpinDG > -6 ? 'warn' : 'poor'}`}>
                                  {alt.forward.hairpinDG?.toFixed(1)} kcal/mol
                                </span>
                              </div>
                            )}
                          </div>
                        </div>

                        <div className="compare-section">
                          <h4>Reverse Primer</h4>
                          <div className="compare-sequence">
                            <code>{alt.reverse.sequence}</code>
                            <button
                              type="button"
                              className="copy-btn-mini"
                              onClick={() => handleCopySequence(alt.reverse.sequence)}
                              title="Copy sequence"
                            >
                              📋
                            </button>
                          </div>
                          <div className="compare-metrics">
                            <div className="metric-row">
                              <span className="metric-label">Length:</span>
                              <span className="metric-value">{alt.reverse.length} bp</span>
                            </div>
                            <div className="metric-row">
                              <span className="metric-label">Tm:</span>
                              <span className="metric-value">{alt.reverse.tm?.toFixed(1)}°C</span>
                            </div>
                            <div className="metric-row">
                              <span className="metric-label">GC%:</span>
                              <span className="metric-value">{revGC.toFixed(1)}%</span>
                            </div>
                            <div className="metric-row">
                              <span className="metric-label">GC Clamp:</span>
                              <span className={`metric-value ${alt.reverse.hasGCClamp ? 'good' : 'warn'}`}>
                                {alt.reverse.hasGCClamp ? 'Yes' : 'No'}
                              </span>
                            </div>
                            {alt.reverse.terminal3DG !== undefined && (
                              <div className="metric-row">
                                <span className="metric-label">3' ΔG:</span>
                                <span className="metric-value">{alt.reverse.terminal3DG?.toFixed(1)} kcal/mol</span>
                              </div>
                            )}
                            {alt.reverse.hairpinDG !== undefined && (
                              <div className="metric-row">
                                <span className="metric-label">Hairpin ΔG:</span>
                                <span className={`metric-value ${alt.reverse.hairpinDG > -3 ? 'good' : alt.reverse.hairpinDG > -6 ? 'warn' : 'poor'}`}>
                                  {alt.reverse.hairpinDG?.toFixed(1)} kcal/mol
                                </span>
                              </div>
                            )}
                          </div>
                        </div>

                        <div className="compare-section">
                          <h4>Pair Analysis</h4>
                          <div className="compare-metrics">
                            <div className="metric-row">
                              <span className="metric-label">ΔTm:</span>
                              <span className={`metric-value ${getTmDiffColorClass(tmDiff)}`}>
                                {tmDiff.toFixed(1)}°C
                                {tmDiff > TM_DIFF_THRESHOLDS.GOOD && ' ⚠️'}
                              </span>
                            </div>
                            <div className="metric-row">
                              <span className="metric-label">GC Clamps:</span>
                              <span className={`metric-value ${gcClamps === 2 ? 'good' : gcClamps === 1 ? 'warn' : 'poor'}`}>
                                {gcClamps}/2
                              </span>
                            </div>
                            {alt.heterodimerDG !== undefined && (
                              <div className="metric-row">
                                <span className="metric-label">Heterodimer ΔG:</span>
                                <span className={`metric-value ${alt.heterodimerDG > -6 ? 'good' : alt.heterodimerDG > -9 ? 'warn' : 'poor'}`}>
                                  {alt.heterodimerDG?.toFixed(1)} kcal/mol
                                </span>
                              </div>
                            )}
                          </div>
                        </div>

                        <div className="compare-actions">
                          <button
                            type="button"
                            className="compare-use-btn"
                            onClick={() => {
                              handleSelectAlternate(alt);
                              setShowCompareModal(false);
                            }}
                          >
                            Use This Pair
                          </button>
                          <button
                            type="button"
                            className="compare-copy-btn"
                            onClick={() => handleCopyPrimerPair(alt.forward.sequence, alt.reverse.sequence)}
                          >
                            📋 Copy Pair
                          </button>
                        </div>
                      </div>
                    );
                  });
                })()}
              </div>
            </div>
            <div className="compare-modal-footer">
              <button
                type="button"
                className="btn-secondary"
                onClick={handleExportAlternatives}
              >
                📥 Export Selected to CSV
              </button>
              <button
                type="button"
                className="btn-primary"
                onClick={() => setShowCompareModal(false)}
              >
                Close
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Toast Notification */}
      {toastMessage && (
        <div className="toast-notification">
          <span className="toast-icon">✓</span>
          <span className="toast-message">{toastMessage}</span>
        </div>
      )}
    </div>
  );
}
