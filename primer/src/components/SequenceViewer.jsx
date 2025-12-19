import React, { useMemo, useState, useCallback, useRef, useEffect } from 'react';
import { SeqViz } from 'seqviz';

// Codon table for amino acid translation
const CODON_TABLE = {
  'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
  'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
  'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
  'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
  'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
  'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
  'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
  'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
  'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
  'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
  'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
  'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
  'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
  'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
};

// Pastel color palette for annotations
const PASTEL_COLORS = {
  promoter: '#fca5a5',
  rbs: '#fdba74',
  cds: '#86efac',
  terminator: '#93c5fd',
  ori: '#c4b5fd',
  gene: '#a5f3fc',
  misc: '#fcd34d',
  primer_fwd: '#60a5fa',
  primer_rev: '#a78bfa',
  amplicon: '#34d399',
  selection: '#fbbf24',
};

// Amino acid color scheme by biochemical properties
const AA_COLORS = {
  // Hydrophobic (warm earth tones)
  A: { bg: '#fef3c7', fg: '#92400e', group: 'hydrophobic' }, // Alanine
  V: { bg: '#fef3c7', fg: '#92400e', group: 'hydrophobic' }, // Valine
  I: { bg: '#fef3c7', fg: '#92400e', group: 'hydrophobic' }, // Isoleucine
  L: { bg: '#fef3c7', fg: '#92400e', group: 'hydrophobic' }, // Leucine
  M: { bg: '#fef3c7', fg: '#92400e', group: 'hydrophobic' }, // Methionine
  F: { bg: '#fed7aa', fg: '#9a3412', group: 'aromatic' },    // Phenylalanine
  W: { bg: '#fed7aa', fg: '#9a3412', group: 'aromatic' },    // Tryptophan
  Y: { bg: '#fed7aa', fg: '#9a3412', group: 'aromatic' },    // Tyrosine
  // Polar uncharged (cool blues)
  S: { bg: '#dbeafe', fg: '#1e40af', group: 'polar' },       // Serine
  T: { bg: '#dbeafe', fg: '#1e40af', group: 'polar' },       // Threonine
  N: { bg: '#dbeafe', fg: '#1e40af', group: 'polar' },       // Asparagine
  Q: { bg: '#dbeafe', fg: '#1e40af', group: 'polar' },       // Glutamine
  C: { bg: '#fce7f3', fg: '#9d174d', group: 'special' },     // Cysteine
  // Positively charged (greens)
  K: { bg: '#dcfce7', fg: '#166534', group: 'positive' },    // Lysine
  R: { bg: '#dcfce7', fg: '#166534', group: 'positive' },    // Arginine
  H: { bg: '#d1fae5', fg: '#065f46', group: 'positive' },    // Histidine
  // Negatively charged (purples)
  D: { bg: '#ede9fe', fg: '#5b21b6', group: 'negative' },    // Aspartic acid
  E: { bg: '#ede9fe', fg: '#5b21b6', group: 'negative' },    // Glutamic acid
  // Special
  G: { bg: '#f3f4f6', fg: '#374151', group: 'special' },     // Glycine
  P: { bg: '#fce7f3', fg: '#9d174d', group: 'special' },     // Proline
  // Stop codon
  '*': { bg: '#fecaca', fg: '#dc2626', group: 'stop' },      // Stop
  '?': { bg: '#e5e7eb', fg: '#6b7280', group: 'unknown' },   // Unknown
};

// Amino acid 3-letter codes
const AA_THREE_LETTER = {
  A: 'Ala', V: 'Val', I: 'Ile', L: 'Leu', M: 'Met',
  F: 'Phe', W: 'Trp', Y: 'Tyr', S: 'Ser', T: 'Thr',
  N: 'Asn', Q: 'Gln', C: 'Cys', K: 'Lys', R: 'Arg',
  H: 'His', D: 'Asp', E: 'Glu', G: 'Gly', P: 'Pro',
  '*': 'Stp', '?': '???',
};

// Nucleotide color schemes for SeqViz bpColors prop
const NUCLEOTIDE_COLORS = {
  classic: {
    A: '#22c55e', // Green - Adenine
    T: '#ef4444', // Red - Thymine
    G: '#f59e0b', // Amber - Guanine
    C: '#3b82f6', // Blue - Cytosine
    U: '#ef4444', // Red - Uracil (RNA)
    N: '#9ca3af', // Gray - Unknown
  },
  muted: {
    A: '#86efac', // Light green
    T: '#fca5a5', // Light red
    G: '#fcd34d', // Light amber
    C: '#93c5fd', // Light blue
    U: '#fca5a5', // Light red
    N: '#d1d5db', // Light gray
  },
  monochrome: null, // No coloring
};

// Calculate Tm using nearest-neighbor method (simplified)
const calculateTm = (seq) => {
  if (!seq || seq.length < 10) return 0;
  const gc = (seq.match(/[GC]/gi) || []).length;
  const at = seq.length - gc;
  if (seq.length < 14) {
    return 2 * at + 4 * gc;
  }
  return 64.9 + 41 * (gc - 16.4) / seq.length;
};

// Calculate GC content
const calculateGC = (seq) => {
  if (!seq) return 0;
  const gc = (seq.match(/[GC]/gi) || []).length;
  return ((gc / seq.length) * 100).toFixed(1);
};

// Reverse complement
const reverseComplement = (seq) => {
  const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
  return seq.split('').reverse().map(c => comp[c.toUpperCase()] || c).join('');
};

// Translate DNA to amino acids
const translateDNA = (seq, frame = 0) => {
  const result = [];
  const s = seq.toUpperCase().slice(frame);
  for (let i = 0; i < s.length - 2; i += 3) {
    const codon = s.slice(i, i + 3);
    result.push(CODON_TABLE[codon] || '?');
  }
  return result.join('');
};

// Translate reverse complement strand
const translateReverseStrand = (seq, frame = 0) => {
  const rc = reverseComplement(seq);
  return translateDNA(rc, frame);
};

// Find ORFs (Open Reading Frames) in translation
const findORFs = (aaSequence, minLength = 20) => {
  const orfs = [];
  let start = -1;

  for (let i = 0; i < aaSequence.length; i++) {
    const aa = aaSequence[i];
    if (aa === 'M' && start === -1) {
      start = i;
    } else if (aa === '*' && start !== -1) {
      const length = i - start;
      if (length >= minLength) {
        orfs.push({ start, end: i, length });
      }
      start = -1;
    }
  }

  return orfs;
};

/**
 * Enhanced SequenceViewer with improved UX
 */
export default function SequenceViewer({
  sequence,
  name = 'Sequence',
  fwdPrimer,
  revPrimer,
  addFwd = '',
  addRev = '',
  viewer: initialViewer = 'linear',
  height = 700,
  annotations: externalAnnotations = [],
  showControls = true,
  showInfo = true,
  showLegend = true,
  compact = false,
  enzymes = [],
  highlightRegions = [],
  circular = false,
  onSequenceChange,
  editable = false,
}) {
  // State
  const [viewer, setViewer] = useState(initialViewer);
  const [zoom, setZoom] = useState(50);
  const [sidebarCollapsed, setSidebarCollapsed] = useState(true); // Start collapsed for more viewer space
  const [splitPosition, setSplitPosition] = useState(50);
  const [isDraggingSplit, setIsDraggingSplit] = useState(false);
  const [showSequenceModal, setShowSequenceModal] = useState(false);
  const [showTranslation, setShowTranslation] = useState(false);
  const [translationFrames, setTranslationFrames] = useState([0]); // 0, 1, 2 for forward; 3, 4, 5 for reverse
  const [showThreeLetter, setShowThreeLetter] = useState(false);
  const [showColoredAA, setShowColoredAA] = useState(true);
  const [nucleotideColorScheme, setNucleotideColorScheme] = useState('classic'); // 'classic', 'muted', 'monochrome'
  const [selection, setSelection] = useState(null);
  const [hoveredAnnotation, setHoveredAnnotation] = useState(null);
  const [viewportStart, setViewportStart] = useState(0);
  const [searchQuery, setSearchQuery] = useState('');
  const [searchResults, setSearchResults] = useState([]);
  const [currentSearchIndex, setCurrentSearchIndex] = useState(0);

  // Refs
  const containerRef = useRef(null);
  const linearViewRef = useRef(null);

  // Calculate GC content
  const gcContent = useMemo(() => calculateGC(sequence), [sequence]);

  // Search functionality
  useEffect(() => {
    if (!searchQuery || !sequence) {
      setSearchResults([]);
      return;
    }
    const results = [];
    const query = searchQuery.toUpperCase();
    const seq = sequence.toUpperCase();
    let idx = seq.indexOf(query);
    while (idx !== -1) {
      results.push({ start: idx, end: idx + query.length });
      idx = seq.indexOf(query, idx + 1);
    }
    setSearchResults(results);
    setCurrentSearchIndex(0);
  }, [searchQuery, sequence]);

  // Handle scroll-to-zoom
  const handleWheel = useCallback((e) => {
    if (e.ctrlKey || e.metaKey) {
      e.preventDefault();
      const delta = e.deltaY > 0 ? -5 : 5;
      setZoom(prev => Math.max(10, Math.min(100, prev + delta)));
    }
  }, []);

  // Handle split pane dragging
  const handleSplitDrag = useCallback((e) => {
    if (!isDraggingSplit || !containerRef.current) return;
    const rect = containerRef.current.getBoundingClientRect();
    const newPos = ((e.clientY - rect.top) / rect.height) * 100;
    setSplitPosition(Math.max(20, Math.min(80, newPos)));
  }, [isDraggingSplit]);

  const handleSplitDragEnd = useCallback(() => {
    setIsDraggingSplit(false);
  }, []);

  useEffect(() => {
    if (isDraggingSplit) {
      document.addEventListener('mousemove', handleSplitDrag);
      document.addEventListener('mouseup', handleSplitDragEnd);
      return () => {
        document.removeEventListener('mousemove', handleSplitDrag);
        document.removeEventListener('mouseup', handleSplitDragEnd);
      };
    }
  }, [isDraggingSplit, handleSplitDrag, handleSplitDragEnd]);

  // Calculate primer positions and create annotations
  const { annotations, primers } = useMemo(() => {
    if (!sequence) {
      return { annotations: [...externalAnnotations, ...highlightRegions], primers: [] };
    }

    const seqUpper = sequence.toUpperCase();
    const annots = externalAnnotations.map(a => ({
      ...a,
      color: a.color || PASTEL_COLORS[a.type] || PASTEL_COLORS.misc,
    }));

    // Add highlight regions
    highlightRegions.forEach(r => {
      annots.push({
        ...r,
        color: r.color || PASTEL_COLORS.selection,
      });
    });

    // Add search results as highlights
    searchResults.forEach((r, i) => {
      annots.push({
        name: `Match ${i + 1}`,
        start: r.start,
        end: r.end,
        color: i === currentSearchIndex ? '#fbbf24' : '#fef3c7',
      });
    });

    const primersList = [];

    // Find FWD primer binding site
    if (fwdPrimer) {
      let fwdSeq = fwdPrimer.seq.toUpperCase();
      let bindingStart = 0;
      let fwdMatch = -1;

      for (let i = 0; i <= fwdSeq.length - 10; i++) {
        const testSeq = fwdSeq.slice(i);
        const idx = seqUpper.indexOf(testSeq);
        if (idx !== -1) {
          fwdMatch = idx;
          bindingStart = i;
          break;
        }
      }

      if (fwdMatch === -1) {
        fwdMatch = seqUpper.indexOf(fwdSeq);
        if (fwdMatch === -1) fwdMatch = 0;
      }

      const fwdEnd = fwdMatch + (fwdSeq.length - bindingStart);

      primersList.push({
        name: `FWD (${fwdSeq.length}bp)`,
        start: fwdMatch,
        end: Math.min(fwdEnd, seqUpper.length),
        direction: 1,
        color: PASTEL_COLORS.primer_fwd,
        tm: fwdPrimer.tm || calculateTm(fwdSeq),
        gc: calculateGC(fwdSeq),
      });

      annots.push({
        name: 'FWD Primer',
        start: fwdMatch,
        end: Math.min(fwdEnd, seqUpper.length),
        direction: 1,
        color: PASTEL_COLORS.primer_fwd,
        type: 'primer',
      });
    }

    // Find REV primer binding site
    if (revPrimer) {
      let revSeq = revPrimer.seq.toUpperCase();
      const revRc = reverseComplement(revSeq);

      let revMatch = -1;
      let bindingLen = revSeq.length;

      for (let i = 0; i <= revRc.length - 10; i++) {
        const testSeq = revRc.slice(0, revRc.length - i);
        const idx = seqUpper.indexOf(testSeq);
        if (idx !== -1) {
          revMatch = idx;
          bindingLen = testSeq.length;
          break;
        }
      }

      if (revMatch === -1) {
        revMatch = seqUpper.indexOf(revRc);
        if (revMatch === -1) revMatch = Math.max(0, seqUpper.length - revSeq.length);
      }

      const revEnd = revMatch + bindingLen;

      primersList.push({
        name: `REV (${revSeq.length}bp)`,
        start: revMatch,
        end: Math.min(revEnd, seqUpper.length),
        direction: -1,
        color: PASTEL_COLORS.primer_rev,
        tm: revPrimer.tm || calculateTm(revSeq),
        gc: calculateGC(revSeq),
      });

      annots.push({
        name: 'REV Primer',
        start: revMatch,
        end: Math.min(revEnd, seqUpper.length),
        direction: -1,
        color: PASTEL_COLORS.primer_rev,
        type: 'primer',
      });
    }

    // Add amplicon annotation if both primers exist
    if (primersList.length === 2) {
      const ampStart = primersList[0].start;
      const ampEnd = primersList[1].end;
      annots.unshift({
        name: `Amplicon (${ampEnd - ampStart}bp)`,
        start: ampStart,
        end: ampEnd,
        direction: 1,
        color: PASTEL_COLORS.amplicon,
        type: 'amplicon',
      });
    }

    // Add selection as annotation
    if (selection) {
      annots.push({
        name: 'Selection',
        start: selection.start,
        end: selection.end,
        color: PASTEL_COLORS.selection,
        type: 'selection',
      });
    }

    return { annotations: annots, primers: primersList };
  }, [sequence, fwdPrimer, revPrimer, externalAnnotations, highlightRegions, searchResults, currentSearchIndex, selection]);

  // Selection stats
  const selectionStats = useMemo(() => {
    if (!selection || !sequence) return null;
    const selectedSeq = sequence.slice(selection.start, selection.end);
    return {
      length: selectedSeq.length,
      gc: calculateGC(selectedSeq),
      tm: calculateTm(selectedSeq).toFixed(1),
      sequence: selectedSeq,
    };
  }, [selection, sequence]);

  // Minimap click handler
  const handleMinimapClick = useCallback((e) => {
    if (!sequence) return;
    const rect = e.currentTarget.getBoundingClientRect();
    const clickPos = (e.clientX - rect.left) / rect.width;
    const newStart = Math.floor(clickPos * sequence.length);
    // Clamp to valid bounds to prevent negative values
    setViewportStart(Math.max(0, Math.min(newStart, sequence.length - 1)));
    // Could integrate with SeqViz scrolling if API supports it
  }, [sequence]);

  // Navigate to annotation
  const navigateToAnnotation = useCallback((annot) => {
    // Clamp to valid bounds
    setViewportStart(Math.max(0, annot.start));
    setHoveredAnnotation(annot);
  }, []);

  if (!sequence) {
    return (
      <div className="sv-empty-state">
        <div className="sv-empty-icon">
          <svg viewBox="0 0 24 24" width="48" height="48" fill="none" stroke="currentColor" strokeWidth="1.5">
            <path d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
          </svg>
        </div>
        <p>No sequence loaded</p>
        <button className="sv-import-btn" onClick={() => setShowSequenceModal(true)}>
          Import Sequence
        </button>
      </div>
    );
  }

  const viewerHeight = viewer === 'both' ? height : height;

  return (
    <div
      className={`sv-enhanced ${compact ? 'compact' : ''} ${sidebarCollapsed ? 'sidebar-collapsed' : ''}`}
      ref={containerRef}
      onWheel={handleWheel}
    >
      {/* Collapsible Sidebar */}
      {!sidebarCollapsed && (
        <aside className="sv-sidebar">
          {/* Sequence Metadata */}
          <div className="sv-metadata-card">
            <div className="sv-metadata-header">
              <h3 className="sv-seq-name">{name}</h3>
              <span className={`sv-seq-type ${circular ? 'circular' : 'linear'}`}>
                {circular ? 'Circular' : 'Linear'}
              </span>
            </div>
            <div className="sv-metadata-stats">
              <div className="sv-stat">
                <span className="sv-stat-value">{sequence.length.toLocaleString()}</span>
                <span className="sv-stat-label">bp</span>
              </div>
              <div className="sv-stat">
                <span className="sv-stat-value">{gcContent}%</span>
                <span className="sv-stat-label">GC</span>
              </div>
              <div className="sv-stat">
                <span className="sv-stat-value">{annotations.length}</span>
                <span className="sv-stat-label">Features</span>
              </div>
            </div>
            <div className="sv-metadata-actions">
              {editable && (
                <button className="sv-action-btn" onClick={() => setShowSequenceModal(true)}>
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                    <path d="M3 17.25V21h3.75L17.81 9.94l-3.75-3.75L3 17.25zM20.71 7.04c.39-.39.39-1.02 0-1.41l-2.34-2.34c-.39-.39-1.02-.39-1.41 0l-1.83 1.83 3.75 3.75 1.83-1.83z"/>
                  </svg>
                  Edit
                </button>
              )}
              <button className="sv-action-btn">
                <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                  <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                </svg>
                Export
              </button>
            </div>
          </div>

          {/* Search */}
          <div className="sv-search-section">
            <div className="sv-search-input-wrap">
              <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor" className="sv-search-icon">
                <path d="M15.5 14h-.79l-.28-.27C15.41 12.59 16 11.11 16 9.5 16 5.91 13.09 3 9.5 3S3 5.91 3 9.5 5.91 16 9.5 16c1.61 0 3.09-.59 4.23-1.57l.27.28v.79l5 4.99L20.49 19l-4.99-5zm-6 0C7.01 14 5 11.99 5 9.5S7.01 5 9.5 5 14 7.01 14 9.5 11.99 14 9.5 14z"/>
              </svg>
              <input
                type="text"
                className="sv-search-input"
                placeholder="Search sequence..."
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value.toUpperCase())}
              />
              {searchResults.length > 0 && (
                <span className="sv-search-count">
                  {currentSearchIndex + 1}/{searchResults.length}
                </span>
              )}
            </div>
            {searchResults.length > 1 && (
              <div className="sv-search-nav">
                <button onClick={() => setCurrentSearchIndex(prev => prev > 0 ? prev - 1 : searchResults.length - 1)}>
                  ↑
                </button>
                <button onClick={() => setCurrentSearchIndex(prev => prev < searchResults.length - 1 ? prev + 1 : 0)}>
                  ↓
                </button>
              </div>
            )}
          </div>

          {/* Annotations List */}
          <div className="sv-annotations-section">
            <div className="sv-section-header sticky">
              <h4>Features</h4>
              <span className="sv-count-badge">{annotations.filter(a => a.type !== 'selection').length}</span>
            </div>
            <div className="sv-annotations-list">
              {annotations.filter(a => a.type !== 'selection').map((annot, i) => (
                <div
                  key={i}
                  className={`sv-annotation-item ${hoveredAnnotation === annot ? 'active' : ''}`}
                  onClick={() => navigateToAnnotation(annot)}
                  onMouseEnter={() => setHoveredAnnotation(annot)}
                  onMouseLeave={() => setHoveredAnnotation(null)}
                >
                  <span className="sv-annot-color" style={{ backgroundColor: annot.color }}></span>
                  <div className="sv-annot-info">
                    <span className="sv-annot-name">{annot.name}</span>
                    <span className="sv-annot-range">{annot.start + 1}..{annot.end}</span>
                  </div>
                  <span className="sv-annot-length">{annot.end - annot.start} bp</span>
                </div>
              ))}
            </div>
          </div>

          {/* Selection Info */}
          {selectionStats && (
            <div className="sv-selection-card">
              <div className="sv-section-header">
                <h4>Selection</h4>
                <button className="sv-clear-btn" onClick={() => setSelection(null)}>×</button>
              </div>
              <div className="sv-selection-stats">
                <div className="sv-sel-stat">
                  <span className="sv-sel-label">Length</span>
                  <span className="sv-sel-value">{selectionStats.length} bp</span>
                </div>
                <div className="sv-sel-stat">
                  <span className="sv-sel-label">GC</span>
                  <span className="sv-sel-value">{selectionStats.gc}%</span>
                </div>
                <div className="sv-sel-stat">
                  <span className="sv-sel-label">Tm</span>
                  <span className="sv-sel-value">{selectionStats.tm}°C</span>
                </div>
              </div>
              <div className="sv-selection-seq">
                <code>{selectionStats.sequence.slice(0, 50)}{selectionStats.sequence.length > 50 ? '...' : ''}</code>
              </div>
            </div>
          )}
        </aside>
      )}

      {/* Main Viewer Area */}
      <main className="sv-main">
        {/* Top Controls */}
        {showControls && (
          <div className="sv-toolbar">
            <div className="sv-toolbar-left">
              {/* Sidebar Toggle Button */}
              <button
                className="sv-sidebar-toggle"
                onClick={() => setSidebarCollapsed(!sidebarCollapsed)}
                title={sidebarCollapsed ? 'Expand sidebar' : 'Collapse sidebar'}
              >
                <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                  {sidebarCollapsed ? (
                    <path d="M3 18h18v-2H3v2zm0-5h18v-2H3v2zm0-7v2h18V6H3z"/>
                  ) : (
                    <path d="M3 18h13v-2H3v2zm0-5h10v-2H3v2zm0-7v2h13V6H3zm18 9.59L17.42 12 21 8.41 19.59 7l-5 5 5 5L21 15.59z"/>
                  )}
                </svg>
              </button>

              <div className="sv-toolbar-divider"></div>

              <div className="sv-view-toggle">
                <button
                  className={`sv-toggle-btn ${viewer === 'linear' ? 'active' : ''}`}
                  onClick={() => setViewer('linear')}
                >
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                    <rect x="2" y="10" width="20" height="4" rx="1" />
                  </svg>
                  Linear
                </button>
                <button
                  className={`sv-toggle-btn ${viewer === 'circular' ? 'active' : ''}`}
                  onClick={() => setViewer('circular')}
                >
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="none" stroke="currentColor" strokeWidth="2">
                    <circle cx="12" cy="12" r="7" />
                  </svg>
                  Circular
                </button>
                <button
                  className={`sv-toggle-btn ${viewer === 'both' ? 'active' : ''}`}
                  onClick={() => setViewer('both')}
                >
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                    <rect x="2" y="3" width="20" height="3" rx="1" />
                    <circle cx="12" cy="14" r="5" fill="none" stroke="currentColor" strokeWidth="1.5" />
                  </svg>
                  Both
                </button>
              </div>

              <div className="sv-toolbar-divider"></div>

              <button
                className={`sv-tool-btn ${showTranslation ? 'active' : ''}`}
                onClick={() => setShowTranslation(!showTranslation)}
                title="Toggle amino acid translation"
              >
                <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                  <path d="M12.87 15.07l-2.54-2.51.03-.03c1.74-1.94 2.98-4.17 3.71-6.53H17V4h-7V2H8v2H1v2h11.17C11.5 7.92 10.44 9.75 9 11.35 8.07 10.32 7.3 9.19 6.69 8h-2c.73 1.63 1.73 3.17 2.98 4.56l-5.09 5.02L4 19l5-5 3.11 3.11.76-2.04zM18.5 10h-2L12 22h2l1.12-3h4.75L21 22h2l-4.5-12zm-2.62 7l1.62-4.33L19.12 17h-3.24z"/>
                </svg>
                AA
              </button>

              <div className="sv-toolbar-divider"></div>

              {/* Nucleotide color scheme selector */}
              <div className="sv-nucleotide-colors">
                <button
                  className={`sv-tool-btn ${nucleotideColorScheme === 'classic' ? 'active' : ''}`}
                  onClick={() => setNucleotideColorScheme(nucleotideColorScheme === 'classic' ? 'monochrome' : 'classic')}
                  title="Toggle nucleotide coloring"
                >
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                    <path d="M12 2C6.49 2 2 6.49 2 12s4.49 10 10 10 10-4.49 10-10S17.51 2 12 2zm0 18c-4.41 0-8-3.59-8-8s3.59-8 8-8 8 3.59 8 8-3.59 8-8 8zm3-8c0 1.66-1.34 3-3 3s-3-1.34-3-3 1.34-3 3-3 3 1.34 3 3z"/>
                  </svg>
                </button>
                {nucleotideColorScheme !== 'monochrome' && (
                  <div className="sv-color-scheme-pills">
                    <button
                      className={`sv-scheme-pill ${nucleotideColorScheme === 'classic' ? 'active' : ''}`}
                      onClick={() => setNucleotideColorScheme('classic')}
                      title="Classic nucleotide colors"
                    >
                      <span className="sv-nt-preview">
                        <span style={{ color: '#22c55e' }}>A</span>
                        <span style={{ color: '#ef4444' }}>T</span>
                        <span style={{ color: '#f59e0b' }}>G</span>
                        <span style={{ color: '#3b82f6' }}>C</span>
                      </span>
                    </button>
                    <button
                      className={`sv-scheme-pill ${nucleotideColorScheme === 'muted' ? 'active' : ''}`}
                      onClick={() => setNucleotideColorScheme('muted')}
                      title="Muted nucleotide colors"
                    >
                      <span className="sv-nt-preview muted">
                        <span style={{ color: '#86efac' }}>A</span>
                        <span style={{ color: '#fca5a5' }}>T</span>
                        <span style={{ color: '#fcd34d' }}>G</span>
                        <span style={{ color: '#93c5fd' }}>C</span>
                      </span>
                    </button>
                  </div>
                )}
              </div>
            </div>

            <div className="sv-toolbar-right">
              <div className="sv-zoom-control">
                <button
                  className="sv-zoom-btn"
                  onClick={() => setZoom(Math.max(10, zoom - 10))}
                >
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                    <path d="M19 13H5v-2h14v2z"/>
                  </svg>
                </button>
                <div className="sv-zoom-track">
                  <input
                    type="range"
                    min="10"
                    max="100"
                    value={zoom}
                    onChange={(e) => setZoom(parseInt(e.target.value))}
                    className="sv-zoom-slider"
                  />
                </div>
                <button
                  className="sv-zoom-btn"
                  onClick={() => setZoom(Math.min(100, zoom + 10))}
                >
                  <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor">
                    <path d="M19 13h-6v6h-2v-6H5v-2h6V5h2v6h6v2z"/>
                  </svg>
                </button>
                <span className="sv-zoom-value">{zoom}%</span>
              </div>
            </div>
          </div>
        )}

        {/* Viewer Container */}
        <div className={`sv-viewer-wrap ${viewer}`}>
          {viewer === 'both' ? (
            <>
              <div className="sv-split-top" style={{ height: `${splitPosition}%` }}>
                <SeqViz
                  name={name}
                  seq={sequence.toUpperCase()}
                  annotations={annotations}
                  primers={primers}
                  enzymes={enzymes}
                  viewer="circular"
                  showComplement={false}
                  showIndex={true}
                  style={{ height: '100%', width: '100%' }}
                  rotateOnScroll={false}
                  bpColors={NUCLEOTIDE_COLORS[nucleotideColorScheme]}
                  onSelection={(sel) => sel && setSelection({ start: sel.start, end: sel.end })}
                />
              </div>
              <div
                className="sv-split-handle"
                onMouseDown={() => setIsDraggingSplit(true)}
              >
                <div className="sv-split-handle-line"></div>
              </div>
              <div className="sv-split-bottom" style={{ height: `${100 - splitPosition}%` }}>
                <SeqViz
                  name={name}
                  seq={sequence.toUpperCase()}
                  annotations={annotations}
                  primers={primers}
                  enzymes={enzymes}
                  viewer="linear"
                  showComplement={true}
                  showIndex={true}
                  style={{ height: '100%', width: '100%' }}
                  zoom={{ linear: zoom }}
                  rotateOnScroll={false}
                  bpColors={NUCLEOTIDE_COLORS[nucleotideColorScheme]}
                  onSelection={(sel) => sel && setSelection({ start: sel.start, end: sel.end })}
                />
              </div>
            </>
          ) : (
            <div className="sv-single-view" ref={linearViewRef}>
              <SeqViz
                name={name}
                seq={sequence.toUpperCase()}
                annotations={annotations}
                primers={primers}
                enzymes={enzymes}
                viewer={viewer}
                showComplement={viewer === 'linear'}
                showIndex={true}
                style={{ height: `${viewerHeight}px`, width: '100%' }}
                zoom={{ linear: zoom }}
                rotateOnScroll={false}
                bpColors={NUCLEOTIDE_COLORS[nucleotideColorScheme]}
                onSelection={(sel) => sel && setSelection({ start: sel.start, end: sel.end })}
              />
            </div>
          )}
        </div>

        {/* Amino Acid Translation Panel - Modern Design */}
        {showTranslation && viewer !== 'circular' && (
          <div className="sv-translation-panel-modern">
            {/* Header with controls */}
            <div className="sv-translation-header-modern">
              <div className="sv-translation-title">
                <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
                  <path d="M12.87 15.07l-2.54-2.51.03-.03c1.74-1.94 2.98-4.17 3.71-6.53H17V4h-7V2H8v2H1v2h11.17C11.5 7.92 10.44 9.75 9 11.35 8.07 10.32 7.3 9.19 6.69 8h-2c.73 1.63 1.73 3.17 2.98 4.56l-5.09 5.02L4 19l5-5 3.11 3.11.76-2.04zM18.5 10h-2L12 22h2l1.12-3h4.75L21 22h2l-4.5-12zm-2.62 7l1.62-4.33L19.12 17h-3.24z"/>
                </svg>
                <span>Translation</span>
              </div>

              <div className="sv-translation-controls">
                {/* Display options */}
                <div className="sv-translation-options">
                  <button
                    className={`sv-option-btn ${showColoredAA ? 'active' : ''}`}
                    onClick={() => setShowColoredAA(!showColoredAA)}
                    title="Color by amino acid properties"
                  >
                    <svg viewBox="0 0 24 24" width="12" height="12" fill="currentColor">
                      <path d="M12 3c-4.97 0-9 4.03-9 9s4.03 9 9 9c.83 0 1.5-.67 1.5-1.5 0-.39-.15-.74-.39-1.01-.23-.26-.38-.61-.38-.99 0-.83.67-1.5 1.5-1.5H16c2.76 0 5-2.24 5-5 0-4.42-4.03-8-9-8zm-5.5 9c-.83 0-1.5-.67-1.5-1.5S5.67 9 6.5 9 8 9.67 8 10.5 7.33 12 6.5 12zm3-4C8.67 8 8 7.33 8 6.5S8.67 5 9.5 5s1.5.67 1.5 1.5S10.33 8 9.5 8zm5 0c-.83 0-1.5-.67-1.5-1.5S13.67 5 14.5 5s1.5.67 1.5 1.5S15.33 8 14.5 8zm3 4c-.83 0-1.5-.67-1.5-1.5S16.67 9 17.5 9s1.5.67 1.5 1.5-.67 1.5-1.5 1.5z"/>
                    </svg>
                  </button>
                  <button
                    className={`sv-option-btn ${showThreeLetter ? 'active' : ''}`}
                    onClick={() => setShowThreeLetter(!showThreeLetter)}
                    title="Show 3-letter amino acid codes"
                  >
                    3L
                  </button>
                </div>

                {/* Frame selector - Segmented control */}
                <div className="sv-frame-selector-modern">
                  <div className="sv-frame-group">
                    <span className="sv-frame-group-label">5'→3'</span>
                    <div className="sv-frame-pills">
                      {[0, 1, 2].map(f => (
                        <button
                          key={f}
                          className={`sv-frame-pill ${translationFrames.includes(f) ? 'active' : ''}`}
                          onClick={() => {
                            setTranslationFrames(prev =>
                              prev.includes(f) ? prev.filter(x => x !== f) : [...prev, f]
                            );
                          }}
                        >
                          +{f + 1}
                        </button>
                      ))}
                    </div>
                  </div>
                  <div className="sv-frame-divider"></div>
                  <div className="sv-frame-group">
                    <span className="sv-frame-group-label">3'→5'</span>
                    <div className="sv-frame-pills reverse">
                      {[3, 4, 5].map(f => (
                        <button
                          key={f}
                          className={`sv-frame-pill reverse ${translationFrames.includes(f) ? 'active' : ''}`}
                          onClick={() => {
                            setTranslationFrames(prev =>
                              prev.includes(f) ? prev.filter(x => x !== f) : [...prev, f]
                            );
                          }}
                        >
                          −{f - 2}
                        </button>
                      ))}
                    </div>
                  </div>
                </div>
              </div>
            </div>

            {/* Translation tracks */}
            <div className="sv-translation-tracks-modern">
              {translationFrames.length === 0 ? (
                <div className="sv-translation-empty">
                  Select a reading frame above to view translation
                </div>
              ) : (
                translationFrames.sort((a, b) => a - b).map(frame => {
                  const isReverse = frame >= 3;
                  const actualFrame = isReverse ? frame - 3 : frame;
                  const translation = isReverse
                    ? translateReverseStrand(sequence, actualFrame)
                    : translateDNA(sequence, actualFrame);
                  const orfs = findORFs(translation, 10);
                  const displayLength = Math.min(80, translation.length);

                  return (
                    <div key={frame} className={`sv-translation-track-modern ${isReverse ? 'reverse' : 'forward'}`}>
                      <div className="sv-frame-indicator">
                        <span className={`sv-frame-badge ${isReverse ? 'reverse' : 'forward'}`}>
                          {isReverse ? `−${actualFrame + 1}` : `+${actualFrame + 1}`}
                        </span>
                        {orfs.length > 0 && (
                          <span className="sv-orf-badge" title={`${orfs.length} ORF${orfs.length > 1 ? 's' : ''} found`}>
                            {orfs.length} ORF
                          </span>
                        )}
                      </div>
                      <div className="sv-aa-track-container">
                        <div className={`sv-aa-track ${showThreeLetter ? 'three-letter' : 'one-letter'}`}>
                          {translation.slice(0, displayLength).split('').map((aa, i) => {
                            const colors = AA_COLORS[aa] || AA_COLORS['?'];
                            const inOrf = orfs.some(orf => i >= orf.start && i <= orf.end);
                            return (
                              <span
                                key={i}
                                className={`sv-aa-residue ${colors.group} ${inOrf ? 'in-orf' : ''} ${aa === 'M' ? 'start' : ''} ${aa === '*' ? 'stop' : ''}`}
                                style={showColoredAA ? { backgroundColor: colors.bg, color: colors.fg } : {}}
                                title={`${AA_THREE_LETTER[aa] || aa} (${aa}) - Position ${isReverse ? 'rev ' : ''}${i + 1}`}
                              >
                                {showThreeLetter ? (AA_THREE_LETTER[aa] || aa) : aa}
                              </span>
                            );
                          })}
                          {translation.length > displayLength && (
                            <span className="sv-aa-more">+{translation.length - displayLength}</span>
                          )}
                        </div>
                      </div>
                    </div>
                  );
                })
              )}
            </div>

            {/* Legend */}
            <div className="sv-aa-legend">
              <div className="sv-aa-legend-item">
                <span className="sv-aa-legend-dot hydrophobic"></span>
                <span>Hydrophobic</span>
              </div>
              <div className="sv-aa-legend-item">
                <span className="sv-aa-legend-dot polar"></span>
                <span>Polar</span>
              </div>
              <div className="sv-aa-legend-item">
                <span className="sv-aa-legend-dot positive"></span>
                <span>Basic (+)</span>
              </div>
              <div className="sv-aa-legend-item">
                <span className="sv-aa-legend-dot negative"></span>
                <span>Acidic (−)</span>
              </div>
              <div className="sv-aa-legend-item">
                <span className="sv-aa-legend-dot stop"></span>
                <span>Stop</span>
              </div>
            </div>
          </div>
        )}

        {/* Minimap Navigator */}
        {viewer === 'linear' && sequence.length > 500 && (
          <div className="sv-minimap" onClick={handleMinimapClick}>
            <div className="sv-minimap-track">
              {annotations.slice(0, 20).map((annot, i) => (
                <div
                  key={i}
                  className="sv-minimap-feature"
                  style={{
                    left: `${(annot.start / sequence.length) * 100}%`,
                    width: `${Math.max(0.5, ((annot.end - annot.start) / sequence.length) * 100)}%`,
                    backgroundColor: annot.color,
                  }}
                  title={annot.name}
                />
              ))}
              <div
                className="sv-minimap-viewport"
                style={{
                  left: `${Math.max(0, Math.min(100, (viewportStart / sequence.length) * 100))}%`,
                  width: `${Math.min(100, (zoom / 100) * 50)}%`,
                }}
              />
            </div>
            <div className="sv-minimap-labels">
              <span>1</span>
              <span>{Math.floor(sequence.length / 2).toLocaleString()}</span>
              <span>{sequence.length.toLocaleString()}</span>
            </div>
          </div>
        )}

        {/* Legend */}
        {showLegend && annotations.length > 0 && (
          <div className="sv-legend-bar">
            {primers.map((p, i) => (
              <span key={`p-${i}`} className="sv-legend-item">
                <span className="sv-legend-dot" style={{ backgroundColor: p.color }}></span>
                {p.name}
              </span>
            ))}
            {annotations.filter(a => a.type !== 'selection' && a.type !== 'primer').slice(0, 5).map((a, i) => (
              <span key={`a-${i}`} className="sv-legend-item">
                <span className="sv-legend-dot" style={{ backgroundColor: a.color }}></span>
                {a.name}
              </span>
            ))}
            {annotations.length > 7 && (
              <span className="sv-legend-more">+{annotations.length - 7} more</span>
            )}
          </div>
        )}
      </main>

      {/* Annotation Tooltip */}
      {hoveredAnnotation && (
        <div className="sv-tooltip" style={{ position: 'fixed', top: '50%', right: sidebarCollapsed ? '20px' : '280px' }}>
          <div className="sv-tooltip-header">
            <span className="sv-tooltip-color" style={{ backgroundColor: hoveredAnnotation.color }}></span>
            <strong>{hoveredAnnotation.name}</strong>
          </div>
          <div className="sv-tooltip-body">
            <div className="sv-tooltip-row">
              <span>Position</span>
              <span>{hoveredAnnotation.start + 1}..{hoveredAnnotation.end}</span>
            </div>
            <div className="sv-tooltip-row">
              <span>Length</span>
              <span>{hoveredAnnotation.end - hoveredAnnotation.start} bp</span>
            </div>
            {hoveredAnnotation.direction && (
              <div className="sv-tooltip-row">
                <span>Strand</span>
                <span>{hoveredAnnotation.direction > 0 ? 'Forward (+)' : 'Reverse (−)'}</span>
              </div>
            )}
            {hoveredAnnotation.tm && (
              <div className="sv-tooltip-row">
                <span>Tm</span>
                <span>{hoveredAnnotation.tm.toFixed(1)}°C</span>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Sequence Edit Modal */}
      {showSequenceModal && (
        <div className="sv-modal-overlay" onClick={() => setShowSequenceModal(false)}>
          <div className="sv-modal" onClick={e => e.stopPropagation()}>
            <div className="sv-modal-header">
              <h3>Edit Sequence</h3>
              <button className="sv-modal-close" onClick={() => setShowSequenceModal(false)}>×</button>
            </div>
            <div className="sv-modal-body">
              <textarea
                className="sv-seq-textarea"
                value={sequence}
                onChange={(e) => onSequenceChange?.(e.target.value.toUpperCase().replace(/[^ATGCN]/g, ''))}
                placeholder="Paste sequence here (FASTA or raw)..."
                rows={12}
              />
            </div>
            <div className="sv-modal-footer">
              <button className="sv-modal-btn secondary" onClick={() => setShowSequenceModal(false)}>
                Cancel
              </button>
              <button className="sv-modal-btn primary" onClick={() => setShowSequenceModal(false)}>
                Apply
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

/**
 * Compact Assembly Viewer for Golden Gate results
 */
export function AssemblyViewer({
  parts = [],
  assembledSequence = '',
  circular = true,
  height = 250,
}) {
  const [viewer, setViewer] = useState(circular ? 'circular' : 'linear');

  const annotations = useMemo(() => {
    const annots = [];
    let position = 0;

    parts.forEach((part, i) => {
      const partLength = part.seq?.length || part.length || 0;
      if (partLength > 0) {
        annots.push({
          name: part.id || `Part ${i + 1}`,
          start: position,
          end: position + partLength,
          direction: 1,
          color: PASTEL_COLORS[part.type] || PASTEL_COLORS.misc,
        });

        if (i < parts.length - 1) {
          annots.push({
            name: `Junction ${i + 1}`,
            start: position + partLength - 4,
            end: position + partLength,
            direction: 1,
            color: '#fcd34d',
          });
        }

        position += partLength;
      }
    });

    return annots;
  }, [parts]);

  if (!assembledSequence || parts.length === 0) {
    return (
      <div className="av-empty">
        <p>No assembled sequence to display</p>
      </div>
    );
  }

  return (
    <div className="av-container">
      <div className="av-header">
        <h4>Assembled Construct</h4>
        <div className="av-toggle">
          <button className={viewer === 'linear' ? 'active' : ''} onClick={() => setViewer('linear')}>
            Linear
          </button>
          <button className={viewer === 'circular' ? 'active' : ''} onClick={() => setViewer('circular')}>
            Circular
          </button>
        </div>
      </div>

      <div className="av-viewer">
        <SeqViz
          name="Assembled Construct"
          seq={assembledSequence.toUpperCase()}
          annotations={annotations}
          viewer={viewer}
          showComplement={false}
          showIndex={true}
          style={{ height: `${height}px`, width: '100%' }}
          zoom={{ linear: 30 }}
          rotateOnScroll={false}
        />
      </div>

      <div className="av-legend">
        {parts.map((part, i) => (
          <span key={i} className="av-legend-item">
            <span className="av-legend-dot" style={{ backgroundColor: PASTEL_COLORS[part.type] || PASTEL_COLORS.misc }}></span>
            {part.id || `Part ${i + 1}`}
          </span>
        ))}
      </div>
    </div>
  );
}
