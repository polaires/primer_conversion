import { useState } from 'react';
import SequenceViewer from './SequenceViewer';

// Common restriction enzymes for quick selection
const COMMON_ENZYMES: string[] = [
  'EcoRI', 'BamHI', 'HindIII', 'XhoI', 'NotI', 'XbaI', 'SpeI', 'PstI',
  'SalI', 'KpnI', 'SacI', 'NcoI', 'NdeI', 'BglII', 'ApaI', 'SmaI',
  'BsaI', 'BpiI', 'AarI', 'SapI', 'BbsI'
];

interface Annotation {
  name: string;
  start: number;
  end: number;
  color: string;
  direction: number;
  type?: string;
}

interface NewAnnotationInput {
  name: string;
  start: string;
  end: string;
  color: string;
  direction: number | string;
}

export default function StandaloneViewer() {
  const [sequence, setSequence] = useState<string>('');
  const [seqName, setSeqName] = useState<string>('My Sequence');
  const [isCircular, setIsCircular] = useState<boolean>(false);

  // Annotations
  const [annotations, setAnnotations] = useState<Annotation[]>([]);
  const [newAnnotation, setNewAnnotation] = useState<NewAnnotationInput>({
    name: '',
    start: '',
    end: '',
    color: '#3b82f6',
    direction: 1
  });

  // Enzymes
  const [selectedEnzymes, setSelectedEnzymes] = useState<string[]>([]);

  // File import
  const handleFileImport = (e: React.ChangeEvent<HTMLInputElement>): void => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      const content = event.target?.result as string;
      // Try to parse as FASTA
      if (content.startsWith('>')) {
        const lines = content.split('\n');
        const name = lines[0].substring(1).trim();
        const seq = lines.slice(1).join('').replace(/\s/g, '').toUpperCase();
        setSeqName(name || 'Imported Sequence');
        setSequence(seq);
      } else {
        // Plain sequence
        setSequence(content.replace(/\s/g, '').toUpperCase());
      }
    };
    reader.readAsText(file);
  };

  // Handle sequence change from the viewer
  const handleSequenceChange = (newSeq: string): void => {
    setSequence(newSeq);
  };

  // Add annotation
  const handleAddAnnotation = (): void => {
    if (!newAnnotation.name || !newAnnotation.start || !newAnnotation.end) return;
    setAnnotations([...annotations, {
      ...newAnnotation,
      start: parseInt(newAnnotation.start, 10),
      end: parseInt(newAnnotation.end, 10),
      direction: typeof newAnnotation.direction === 'string'
        ? parseInt(newAnnotation.direction, 10)
        : newAnnotation.direction
    }]);
    setNewAnnotation({ name: '', start: '', end: '', color: '#3b82f6', direction: 1 });
  };

  // Remove annotation
  const handleRemoveAnnotation = (index: number): void => {
    setAnnotations(annotations.filter((_, i) => i !== index));
  };

  // Toggle enzyme
  const handleToggleEnzyme = (enzyme: string): void => {
    if (selectedEnzymes.includes(enzyme)) {
      setSelectedEnzymes(selectedEnzymes.filter(e => e !== enzyme));
    } else {
      setSelectedEnzymes([...selectedEnzymes, enzyme]);
    }
  };

  // Export as FASTA
  const handleExportFasta = (): void => {
    const fasta = `>${seqName}\n${sequence.match(/.{1,80}/g)?.join('\n') || ''}`;
    const blob = new Blob([fasta], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${seqName.replace(/\s+/g, '_')}.fasta`;
    a.click();
    URL.revokeObjectURL(url);
  };

  // Sample sequences
  const loadSampleSequence = (type: 'plasmid' | 'gene'): void => {
    if (type === 'plasmid') {
      setSeqName('pUC19');
      setIsCircular(true);
      setSequence('TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTC');
      setAnnotations([
        { name: 'lacZ', start: 146, end: 467, color: '#93c5fd', direction: 1, type: 'gene' },
        { name: 'AmpR', start: 1629, end: 2489, color: '#fca5a5', direction: -1, type: 'gene' },
        { name: 'ori', start: 927, end: 1515, color: '#86efac', direction: 1, type: 'ori' }
      ]);
    } else if (type === 'gene') {
      setSeqName('GFP');
      setIsCircular(false);
      setSequence('ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA');
      setAnnotations([
        { name: 'Start Codon', start: 0, end: 3, color: '#86efac', direction: 1, type: 'cds' },
        { name: 'Chromophore', start: 195, end: 207, color: '#fcd34d', direction: 1, type: 'misc' },
        { name: 'Stop Codon', start: 714, end: 717, color: '#fca5a5', direction: 1, type: 'cds' }
      ]);
    }
  };

  return (
    <div className="standalone-viewer-v2">
      {/* Quick Actions Bar */}
      <div className="sv-quick-actions">
        <div className="sv-quick-left">
          <label className="sv-file-btn">
            <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
              <path d="M9 16h6v-6h4l-7-7-7 7h4v6zm-4 2h14v2H5v-2z"/>
            </svg>
            Import
            <input type="file" accept=".fasta,.fa,.txt,.gb,.gbk" onChange={handleFileImport} hidden />
          </label>
          <button type="button" className="sv-action-btn" onClick={handleExportFasta} disabled={!sequence}>
            <svg viewBox="0 0 24 24" width="16" height="16" fill="currentColor">
              <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
            </svg>
            Export
          </button>
          <div className="sv-samples">
            <span>Load:</span>
            <button type="button" onClick={() => loadSampleSequence('plasmid')}>Plasmid</button>
            <button type="button" onClick={() => loadSampleSequence('gene')}>Gene</button>
          </div>
        </div>
        <div className="sv-quick-right">
          <label className="sv-checkbox">
            <input
              type="checkbox"
              checked={isCircular}
              onChange={(e) => setIsCircular(e.target.checked)}
            />
            Circular
          </label>
        </div>
      </div>

      {/* Main Content */}
      {sequence ? (
        <SequenceViewer
          sequence={sequence}
          name={seqName}
          annotations={annotations}
          enzymes={selectedEnzymes}
          circular={isCircular}
          height={700}
          editable={true}
          onSequenceChange={handleSequenceChange}
        />
      ) : (
        <div className="sv-empty-placeholder">
          <div className="sv-empty-content">
            <svg viewBox="0 0 24 24" width="64" height="64" fill="none" stroke="currentColor" strokeWidth="1">
              <path d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
            </svg>
            <h3>No Sequence Loaded</h3>
            <p>Import a file or paste a sequence to get started</p>
            <div className="sv-empty-actions">
              <label className="sv-primary-btn">
                Import File
                <input type="file" accept=".fasta,.fa,.txt,.gb,.gbk" onChange={handleFileImport} hidden />
              </label>
              <span>or try a sample:</span>
              <button type="button" onClick={() => loadSampleSequence('plasmid')}>pUC19 Plasmid</button>
              <button type="button" onClick={() => loadSampleSequence('gene')}>GFP Gene</button>
            </div>
          </div>
        </div>
      )}

      {/* Bottom Panel - Annotations & Enzymes */}
      {sequence && (
        <div className="sv-bottom-panel">
          {/* Add Annotation */}
          <div className="sv-panel-section">
            <h4>Add Annotation</h4>
            <div className="sv-add-annotation">
              <input
                type="text"
                value={newAnnotation.name}
                onChange={(e) => setNewAnnotation({...newAnnotation, name: e.target.value})}
                placeholder="Name"
              />
              <input
                type="number"
                value={newAnnotation.start}
                onChange={(e) => setNewAnnotation({...newAnnotation, start: e.target.value})}
                placeholder="Start"
                min="0"
              />
              <input
                type="number"
                value={newAnnotation.end}
                onChange={(e) => setNewAnnotation({...newAnnotation, end: e.target.value})}
                placeholder="End"
                min="0"
              />
              <input
                type="color"
                value={newAnnotation.color}
                onChange={(e) => setNewAnnotation({...newAnnotation, color: e.target.value})}
              />
              <select
                value={newAnnotation.direction}
                onChange={(e) => setNewAnnotation({...newAnnotation, direction: e.target.value})}
              >
                <option value={1}>→</option>
                <option value={-1}>←</option>
              </select>
              <button type="button" onClick={handleAddAnnotation}>Add</button>
            </div>
            {annotations.length > 0 && (
              <div className="sv-annotations-chips">
                {annotations.map((ann, i) => (
                  <span key={i} className="sv-annotation-chip" style={{ borderColor: ann.color }}>
                    <span className="chip-color" style={{ background: ann.color }}></span>
                    {ann.name}
                    <button type="button" onClick={() => handleRemoveAnnotation(i)}>×</button>
                  </span>
                ))}
              </div>
            )}
          </div>

          {/* Restriction Enzymes */}
          <div className="sv-panel-section">
            <h4>
              Restriction Enzymes
              {selectedEnzymes.length > 0 && (
                <button type="button" className="sv-clear-link" onClick={() => setSelectedEnzymes([])}>
                  Clear
                </button>
              )}
            </h4>
            <div className="sv-enzyme-grid">
              {COMMON_ENZYMES.map(enzyme => (
                <label key={enzyme} className={`sv-enzyme-chip ${selectedEnzymes.includes(enzyme) ? 'selected' : ''}`}>
                  <input
                    type="checkbox"
                    checked={selectedEnzymes.includes(enzyme)}
                    onChange={() => handleToggleEnzyme(enzyme)}
                    hidden
                  />
                  {enzyme}
                </label>
              ))}
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
