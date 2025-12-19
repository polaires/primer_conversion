/**
 * Sequence Utility Functions
 *
 * Handles file parsing (FASTA, GenBank), sequence manipulation,
 * ambiguous base expansion, and common vectors.
 */

// =============================================================================
// Types and Interfaces
// =============================================================================

export interface Vector {
  name: string;
  length: number;
  description: string;
  sequence: string | null;
}

export interface VectorDatabase {
  [key: string]: Vector;
}

export interface IupacCodes {
  [key: string]: string[];
}

export interface FastaEntry {
  id: string;
  description: string;
  sequence: string;
}

export interface GenBankFeature {
  type: string;
  start: number;
  end: number;
  complement: boolean;
}

export interface GenBankEntry {
  id: string;
  description: string;
  sequence: string;
  features: GenBankFeature[];
  length: number;
}

export type SequenceEntry = FastaEntry | GenBankEntry;

export type SequenceFormat = 'fasta' | 'genbank' | 'raw' | 'unknown';

export interface ParsedSequenceFile {
  format: SequenceFormat;
  entries: SequenceEntry[];
}

export interface FrameTranslation {
  frame: number;
  protein: string;
}

export interface AAMutation {
  type: 'aa_mutation';
  oldAA: string;
  position: number;
  newAA: string;
  orfStart: number;
  original: string;
}

export interface IndelMutation {
  type: 'indel';
  id: string;
  start: number;
  end: number;
  replacement: string;
  original: string;
}

export interface ErrorMutation {
  type: 'error';
  message: string;
  original: string;
}

export type BatchMutation = AAMutation | IndelMutation | ErrorMutation;

export interface Primer {
  name?: string;
  sequence: string;
  note?: string;
  length?: number;
  tm?: number;
  gcPercent?: number;
  dg?: number | string;
}

export interface PrimerPair {
  forward?: Primer;
  reverse?: Primer;
  id?: string;
  description?: string;
}

export interface ComplementMap {
  [key: string]: string;
}

export interface CodonTable {
  [key: string]: string;
}

// =============================================================================
// Common Vectors Database
// =============================================================================

export const COMMON_VECTORS: VectorDatabase = {
  'pUC19': {
    name: 'pUC19',
    length: 2686,
    description: 'High-copy cloning vector',
    sequence: 'TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTC',
  },
  'pET-28a': {
    name: 'pET-28a(+)',
    length: 5369,
    description: 'T7 expression vector with N-terminal His-tag',
    sequence: null, // Would include full sequence
  },
  'pGEX-4T-1': {
    name: 'pGEX-4T-1',
    length: 4969,
    description: 'GST fusion expression vector',
    sequence: null,
  },
  'pcDNA3.1': {
    name: 'pcDNA3.1(+)',
    length: 5428,
    description: 'Mammalian expression vector',
    sequence: null,
  },
};

// =============================================================================
// Ambiguous Base Codes (IUPAC)
// =============================================================================

export const IUPAC_CODES: IupacCodes = {
  'A': ['A'],
  'T': ['T'],
  'G': ['G'],
  'C': ['C'],
  'R': ['A', 'G'],        // puRine
  'Y': ['C', 'T'],        // pYrimidine
  'S': ['G', 'C'],        // Strong
  'W': ['A', 'T'],        // Weak
  'K': ['G', 'T'],        // Keto
  'M': ['A', 'C'],        // aMino
  'B': ['C', 'G', 'T'],   // not A
  'D': ['A', 'G', 'T'],   // not C
  'H': ['A', 'C', 'T'],   // not G
  'V': ['A', 'C', 'G'],   // not T
  'N': ['A', 'C', 'G', 'T'],
};

/**
 * Expand ambiguous sequence into all possible combinations
 * @param sequence - Sequence with IUPAC codes
 * @returns Array of expanded sequences
 */
export function expandAmbiguousBases(sequence: string): string[] {
  const seq = sequence.toUpperCase();
  let combinations: string[] = [''];

  for (const base of seq) {
    const expansions = IUPAC_CODES[base] || [base];
    const newCombinations: string[] = [];

    for (const combo of combinations) {
      for (const exp of expansions) {
        newCombinations.push(combo + exp);
      }
    }
    combinations = newCombinations;

    // Limit to prevent memory issues (e.g., NNNNNN = 4096 combinations)
    if (combinations.length > 10000) {
      throw new Error('Too many combinations (>10000). Reduce ambiguous bases.');
    }
  }

  return combinations;
}

/**
 * Check if sequence contains ambiguous bases
 */
export function hasAmbiguousBases(sequence: string): boolean {
  return /[RYSWKMBDHVN]/i.test(sequence);
}

/**
 * Count expected combinations from ambiguous sequence
 */
export function countCombinations(sequence: string): number {
  let count = 1;
  for (const base of sequence.toUpperCase()) {
    const expansions = IUPAC_CODES[base];
    if (expansions) {
      count *= expansions.length;
    }
  }
  return count;
}

// =============================================================================
// File Format Parsing
// =============================================================================

/**
 * Detect file format from content
 */
export function detectFormat(content: string): SequenceFormat {
  const trimmed = content.trim();

  if (trimmed.startsWith('>')) {
    return 'fasta';
  }
  if (trimmed.startsWith('LOCUS') || /^LOCUS\s+/m.test(trimmed)) {
    return 'genbank';
  }
  if (/^[ATGCatgcNnRYSWKMBDHV\s\d]+$/.test(trimmed.replace(/[\r\n]/g, ''))) {
    return 'raw';
  }
  return 'unknown';
}

/**
 * Parse FASTA format
 * Returns { id, description, sequence }
 */
export function parseFasta(content: string): FastaEntry[] {
  const lines = content.trim().split(/[\r\n]+/);
  const results: FastaEntry[] = [];
  let currentEntry: FastaEntry | null = null;

  for (const line of lines) {
    if (line.startsWith('>')) {
      if (currentEntry) {
        results.push(currentEntry);
      }
      const header = line.slice(1).trim();
      const [id, ...descParts] = header.split(/\s+/);
      currentEntry = {
        id: id || 'Unknown',
        description: descParts.join(' '),
        sequence: '',
      };
    } else if (currentEntry) {
      currentEntry.sequence += line.replace(/\s/g, '').toUpperCase();
    }
  }

  if (currentEntry) {
    results.push(currentEntry);
  }

  return results;
}

/**
 * Parse GenBank format (simplified)
 * Returns { id, description, sequence, features }
 */
export function parseGenBank(content: string): GenBankEntry[] {
  const results: GenBankEntry[] = [];
  const entries = content.split(/\/\/\s*[\r\n]+/).filter(e => e.trim());

  for (const entry of entries) {
    const result: GenBankEntry = {
      id: '',
      description: '',
      sequence: '',
      features: [],
      length: 0,
    };

    // Parse LOCUS line
    const locusMatch = entry.match(/LOCUS\s+(\S+)\s+(\d+)\s+bp/);
    if (locusMatch) {
      result.id = locusMatch[1];
      result.length = parseInt(locusMatch[2]);
    }

    // Parse DEFINITION
    const defMatch = entry.match(/DEFINITION\s+(.+?)(?=\n[A-Z])/s);
    if (defMatch) {
      result.description = defMatch[1].replace(/\s+/g, ' ').trim();
    }

    // Parse FEATURES section for CDS/gene locations
    const featuresMatch = entry.match(/FEATURES\s+Location\/Qualifiers\s*([\s\S]*?)(?=ORIGIN|$)/);
    if (featuresMatch) {
      const featuresText = featuresMatch[1];

      // Find CDS features
      const cdsRegex = /CDS\s+(?:complement\()?(\d+)\.\.(\d+)\)?/g;
      let match: RegExpExecArray | null;
      while ((match = cdsRegex.exec(featuresText)) !== null) {
        result.features.push({
          type: 'CDS',
          start: parseInt(match[1]),
          end: parseInt(match[2]),
          complement: featuresText.slice(match.index - 20, match.index).includes('complement'),
        });
      }

      // Find gene features
      const geneRegex = /gene\s+(?:complement\()?(\d+)\.\.(\d+)\)?/g;
      while ((match = geneRegex.exec(featuresText)) !== null) {
        result.features.push({
          type: 'gene',
          start: parseInt(match[1]),
          end: parseInt(match[2]),
          complement: featuresText.slice(match.index - 20, match.index).includes('complement'),
        });
      }
    }

    // Parse ORIGIN sequence
    const originMatch = entry.match(/ORIGIN\s*([\s\S]*?)$/);
    if (originMatch) {
      result.sequence = originMatch[1]
        .replace(/[\d\s\/]+/g, '')
        .toUpperCase();
    }

    if (result.sequence) {
      results.push(result);
    }
  }

  return results;
}

/**
 * Parse raw sequence (just nucleotides)
 */
export function parseRaw(content: string): FastaEntry[] {
  const sequence = content
    .replace(/[\d\s\r\n]+/g, '')
    .toUpperCase();

  return [{
    id: 'Sequence',
    description: '',
    sequence,
  }];
}

/**
 * Auto-detect and parse sequence file
 */
export function parseSequenceFile(content: string): ParsedSequenceFile {
  const format = detectFormat(content);

  switch (format) {
    case 'fasta':
      return { format, entries: parseFasta(content) };
    case 'genbank':
      return { format, entries: parseGenBank(content) };
    case 'raw':
      return { format, entries: parseRaw(content) };
    default:
      throw new Error('Unknown file format. Supported: FASTA, GenBank, or raw sequence.');
  }
}

// =============================================================================
// Sequence Manipulation
// =============================================================================

/**
 * Get reverse complement of sequence
 */
export function reverseComplement(sequence: string): string {
  const complement: ComplementMap = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D', 'N': 'N',
  };

  return sequence
    .toUpperCase()
    .split('')
    .reverse()
    .map(base => complement[base] || base)
    .join('');
}

/**
 * Translate DNA to protein (single reading frame)
 */
export function translateDNA(sequence: string, frame: number = 0): string {
  const codonTable: CodonTable = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
  };

  const seq = sequence.toUpperCase().slice(frame);
  let protein = '';

  for (let i = 0; i < seq.length - 2; i += 3) {
    const codon = seq.slice(i, i + 3);
    protein += codonTable[codon] || 'X';
  }

  return protein;
}

/**
 * Get 3-frame translation
 */
export function get3FrameTranslation(sequence: string): FrameTranslation[] {
  return [
    { frame: 1, protein: translateDNA(sequence, 0) },
    { frame: 2, protein: translateDNA(sequence, 1) },
    { frame: 3, protein: translateDNA(sequence, 2) },
  ];
}

// =============================================================================
// Batch Mutation Parsing
// =============================================================================

/**
 * Parse batch mutation commands
 * Supports:
 * - ORF: [position] - Set ORF start
 * - [oldAA][codon#][newAA] - AA mutation (e.g., M1Y, A42S)
 * - [ID] [start] [end] [replacement] - Indel operation
 */
export function parseBatchCommands(text: string, defaultOrfStart: number = 1): BatchMutation[] {
  const lines = text.trim().split(/[\r\n]+/).filter(l => l.trim() && !l.trim().startsWith('#'));
  const mutations: BatchMutation[] = [];
  let currentOrfStart = defaultOrfStart;

  for (const line of lines) {
    const trimmed = line.trim();

    // ORF directive
    if (trimmed.toUpperCase().startsWith('ORF:')) {
      const pos = parseInt(trimmed.slice(4).trim());
      if (!isNaN(pos)) {
        currentOrfStart = pos;
      }
      continue;
    }

    // AA mutation format: [old][codon#][new] e.g., M1Y, A42S, F32H
    const aaMutationMatch = trimmed.match(/^([A-Z])(\d+)([A-Z*])$/i);
    if (aaMutationMatch) {
      mutations.push({
        type: 'aa_mutation',
        oldAA: aaMutationMatch[1].toUpperCase(),
        position: parseInt(aaMutationMatch[2]),
        newAA: aaMutationMatch[3].toUpperCase(),
        orfStart: currentOrfStart,
        original: trimmed,
      });
      continue;
    }

    // Indel format: [ID] [start] [end] [replacement?]
    const indelMatch = trimmed.match(/^(\S+)\s+(\d+)\s+(\d+)(?:\s+([A-Z]+))?$/i);
    if (indelMatch) {
      mutations.push({
        type: 'indel',
        id: indelMatch[1],
        start: parseInt(indelMatch[2]),
        end: parseInt(indelMatch[3]),
        replacement: indelMatch[4]?.toUpperCase() || '',
        original: trimmed,
      });
      continue;
    }

    // Unknown format
    mutations.push({
      type: 'error',
      message: `Unknown format: ${trimmed}`,
      original: trimmed,
    });
  }

  return mutations;
}

// =============================================================================
// Download Utilities
// =============================================================================

/**
 * Generate tab-delimited primer list for ordering
 */
export function generatePrimerOrderList(primers: Primer[], projectName: string = 'Primers'): string {
  const lines = ['Name\tSequence\tNote'];

  primers.forEach((primer, index) => {
    const name = primer.name || `${projectName}_${index + 1}`;
    lines.push(`${name}\t${primer.sequence}\t${primer.note || ''}`);
  });

  return lines.join('\n');
}

/**
 * Generate FASTA format
 */
export function generateFasta(sequences: FastaEntry[]): string {
  return sequences
    .map(seq => `>${seq.id}${seq.description ? ' ' + seq.description : ''}\n${seq.sequence}`)
    .join('\n\n');
}

/**
 * Generate CSV summary
 */
export function generateCSVSummary(results: PrimerPair[]): string {
  const headers = ['Primer Name', 'Sequence', 'Length', 'Tm (°C)', 'GC%', 'ΔG', 'Note'];
  const lines = [headers.join(',')];

  results.forEach(r => {
    if (r.forward) {
      lines.push([
        `"${r.id || 'Forward'}_F"`,
        `"${r.forward.sequence}"`,
        r.forward.length,
        r.forward.tm,
        r.forward.gcPercent,
        r.forward.dg || '',
        `"${r.description || ''}"`,
      ].join(','));
    }
    if (r.reverse) {
      lines.push([
        `"${r.id || 'Reverse'}_R"`,
        `"${r.reverse.sequence}"`,
        r.reverse.length,
        r.reverse.tm,
        r.reverse.gcPercent,
        r.reverse.dg || '',
        `""`,
      ].join(','));
    }
  });

  return lines.join('\n');
}

/**
 * Trigger file download in browser
 */
export function downloadFile(content: string, filename: string, mimeType: string = 'text/plain'): void {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}
