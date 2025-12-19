/**
 * Sequence Utilities Module
 * Provides sequence parsing, alignment, and matching for DNA assembly
 * Based on: https://github.com/Lattice-Automation/repp
 */

import { reverseComplement } from './enzymes.js';

/**
 * Fragment representation from FASTA file
 */
export interface Fragment {
  id: string;
  seq: string;
  circular: boolean;
}

/**
 * GenBank feature
 */
export interface GenbankFeature {
  start: number;
  end: number;
  label: string;
}

/**
 * GenBank parse result
 */
export interface GenbankResult {
  id: string;
  seq: string;
  features: GenbankFeature[];
}

/**
 * Sequence match result
 */
export interface SequenceMatch {
  queryStart: number;
  queryEnd: number;
  subjectStart: number;
  subjectEnd: number;
  length: number;
  seq: string;
  mismatches: number;
  forward: boolean;
}

/**
 * Extended match result with additional properties
 */
interface ExtendedMatch {
  length: number;
  seq: string;
  mismatches: number;
}

/**
 * Find matches options
 */
export interface FindMatchesOptions {
  minLength?: number;
  maxMismatches?: number;
  circular?: boolean;
}

/**
 * FASTA format options
 */
export interface FastaOptions {
  lineWidth?: number;
  circular?: boolean;
}

/**
 * Parse a FASTA format string into fragments
 * @param content - FASTA content
 * @returns Array of { id, seq, circular } objects
 */
export function parseFasta(content: string): Fragment[] {
  const fragments: Fragment[] = [];
  const lines = content.split('\n');

  let currentId: string | null = null;
  let currentSeq: string[] = [];
  let isCircular = false;

  for (const line of lines) {
    const trimmed = line.trim();

    if (trimmed.startsWith('>')) {
      // Save previous fragment
      if (currentId !== null) {
        fragments.push({
          id: currentId,
          seq: currentSeq.join('').toUpperCase().replace(/[^ATGC]/gi, ''),
          circular: isCircular,
        });
      }

      // Parse new header
      currentId = trimmed.slice(1).trim();
      currentSeq = [];
      isCircular = trimmed.toLowerCase().includes('circular');
    } else if (trimmed && currentId !== null) {
      currentSeq.push(trimmed);
    }
  }

  // Save last fragment
  if (currentId !== null && currentSeq.length > 0) {
    fragments.push({
      id: currentId,
      seq: currentSeq.join('').toUpperCase().replace(/[^ATGC]/gi, ''),
      circular: isCircular,
    });
  }

  return fragments;
}

/**
 * Parse a GenBank format string
 * @param content - GenBank content
 * @returns { id, seq, features }
 */
export function parseGenbank(content: string): GenbankResult {
  const parts = content.split('ORIGIN');

  if (parts.length !== 2) {
    throw new Error('Invalid GenBank format: missing ORIGIN section');
  }

  // Parse sequence
  const seqPart = parts[1];
  const seq = seqPart.toUpperCase().replace(/[^ATGC]/g, '');

  // Parse ID from LOCUS
  const locusMatch = parts[0].match(/LOCUS\s+(\S+)/);
  const id = locusMatch ? locusMatch[1] : 'unknown';

  // Parse features (simplified)
  const features: GenbankFeature[] = [];
  const featureMatch = parts[0].match(/FEATURES[\s\S]*?(?=ORIGIN|$)/);

  if (featureMatch) {
    const featureRegex = /(\d+)\.\.(\d+)[\s\S]*?\/label="?([^"\n]+)"?/g;
    let match;
    while ((match = featureRegex.exec(featureMatch[0])) !== null) {
      features.push({
        start: parseInt(match[1]) - 1,  // Convert to 0-indexed
        end: parseInt(match[2]) - 1,
        label: match[3].trim(),
      });
    }
  }

  return { id, seq, features };
}

/**
 * Find matching regions between a query and subject sequence
 * Uses a simplified local alignment approach suitable for browser
 * @param query - Query sequence
 * @param subject - Subject sequence (from fragment database)
 * @param options - Search options
 * @returns Array of match objects
 */
export function findMatches(query: string, subject: string, options: FindMatchesOptions = {}): SequenceMatch[] {
  const {
    minLength = 60,
    maxMismatches = 2,
    circular = true,
  } = options;

  const matches: SequenceMatch[] = [];
  const queryUpper = query.toUpperCase();
  const subjectUpper = subject.toUpperCase();

  // For circular targets, double the query to find matches across zero-index
  const searchQuery = circular ? queryUpper + queryUpper : queryUpper;

  // Use k-mer based approach for efficiency
  const kmerSize = 15;
  const kmerIndex = buildKmerIndex(subjectUpper, kmerSize);

  // Find seed matches and extend them
  for (let i = 0; i < searchQuery.length - kmerSize; i++) {
    const kmer = searchQuery.slice(i, i + kmerSize);
    const subjectPositions = kmerIndex.get(kmer);

    if (!subjectPositions) continue;

    for (const subjectPos of subjectPositions) {
      // Try to extend the match
      const extended = extendMatch(searchQuery, subjectUpper, i, subjectPos, maxMismatches);

      if (extended && extended.length >= minLength) {
        // Normalize query positions for circular sequences
        const queryStart = i % queryUpper.length;
        const queryEnd = (i + extended.length - 1) % queryUpper.length;

        matches.push({
          queryStart,
          queryEnd: queryEnd < queryStart ? queryEnd + queryUpper.length : queryEnd,
          subjectStart: subjectPos,
          subjectEnd: subjectPos + extended.length - 1,
          length: extended.length,
          seq: extended.seq,
          mismatches: extended.mismatches,
          forward: true,
        });
      }
    }
  }

  // Also check reverse complement of subject
  const rcSubject = reverseComplement(subjectUpper);
  const rcKmerIndex = buildKmerIndex(rcSubject, kmerSize);

  for (let i = 0; i < searchQuery.length - kmerSize; i++) {
    const kmer = searchQuery.slice(i, i + kmerSize);
    const rcPositions = rcKmerIndex.get(kmer);

    if (!rcPositions) continue;

    for (const rcPos of rcPositions) {
      const extended = extendMatch(searchQuery, rcSubject, i, rcPos, maxMismatches);

      if (extended && extended.length >= minLength) {
        const queryStart = i % queryUpper.length;
        const queryEnd = (i + extended.length - 1) % queryUpper.length;

        // Convert RC position back to forward strand
        const subjectStart = subjectUpper.length - rcPos - extended.length;
        const subjectEnd = subjectUpper.length - rcPos - 1;

        matches.push({
          queryStart,
          queryEnd: queryEnd < queryStart ? queryEnd + queryUpper.length : queryEnd,
          subjectStart,
          subjectEnd,
          length: extended.length,
          seq: extended.seq,
          mismatches: extended.mismatches,
          forward: false,
        });
      }
    }
  }

  // Remove duplicate and overlapping matches
  return cullMatches(matches, minLength);
}

/**
 * Build k-mer index for fast lookups
 * @param seq - Sequence to index
 * @param k - K-mer size
 * @returns Map from kmer to array of positions
 */
function buildKmerIndex(seq: string, k: number): Map<string, number[]> {
  const index = new Map<string, number[]>();

  for (let i = 0; i <= seq.length - k; i++) {
    const kmer = seq.slice(i, i + k);
    if (!index.has(kmer)) {
      index.set(kmer, []);
    }
    index.get(kmer)!.push(i);
  }

  return index;
}

/**
 * Extend a seed match in both directions
 * @param query - Query sequence
 * @param subject - Subject sequence
 * @param queryPos - Query start position
 * @param subjectPos - Subject start position
 * @param maxMismatches - Maximum allowed mismatches
 * @returns Extended match or null
 */
function extendMatch(
  query: string,
  subject: string,
  queryPos: number,
  subjectPos: number,
  maxMismatches: number
): ExtendedMatch | null {
  let qStart = queryPos;
  let sStart = subjectPos;
  let qEnd = queryPos;
  let sEnd = subjectPos;
  let mismatches = 0;

  // Extend right
  while (qEnd < query.length && sEnd < subject.length) {
    if (query[qEnd] !== subject[sEnd]) {
      mismatches++;
      if (mismatches > maxMismatches) break;
    }
    qEnd++;
    sEnd++;
  }

  // Extend left
  while (qStart > 0 && sStart > 0) {
    if (query[qStart - 1] !== subject[sStart - 1]) {
      mismatches++;
      if (mismatches > maxMismatches) {
        mismatches--;  // Undo the increment
        break;
      }
    }
    qStart--;
    sStart--;
  }

  const length = qEnd - qStart;
  if (length < 1) return null;

  return {
    length,
    seq: query.slice(qStart, qEnd),
    mismatches,
  };
}

/**
 * Remove overlapping and duplicate matches
 * @param matches - Array of matches
 * @param minLength - Minimum match length
 * @returns Culled matches
 */
function cullMatches(matches: SequenceMatch[], minLength: number): SequenceMatch[] {
  if (matches.length === 0) return [];

  // Sort by query start, then by length (longest first)
  matches.sort((a, b) => {
    if (a.queryStart !== b.queryStart) {
      return a.queryStart - b.queryStart;
    }
    return b.length - a.length;
  });

  const culled: SequenceMatch[] = [];

  for (const match of matches) {
    if (match.length < minLength) continue;

    // Check if this match significantly overlaps with a previous one
    // We consider matches overlapping if they share >80% of their length
    let overlapping = false;
    for (const prev of culled) {
      // Calculate overlap between prev and match on the query
      const overlapStart = Math.max(prev.queryStart, match.queryStart);
      const overlapEnd = Math.min(prev.queryEnd, match.queryEnd);
      const overlapLength = Math.max(0, overlapEnd - overlapStart + 1);

      // Check if overlap is significant (>80% of either match)
      const prevOverlapRatio = overlapLength / (prev.queryEnd - prev.queryStart + 1);
      const matchOverlapRatio = overlapLength / (match.queryEnd - match.queryStart + 1);

      if (prevOverlapRatio > 0.8 || matchOverlapRatio > 0.8) {
        overlapping = true;
        break;
      }
    }

    if (!overlapping) {
      culled.push(match);
    }
  }

  return culled;
}

/**
 * Find the junction (overlapping region) between two sequences
 * @param seq1 - First sequence (end checked)
 * @param seq2 - Second sequence (start checked)
 * @param minHomology - Minimum overlap length
 * @param maxHomology - Maximum overlap length
 * @returns Junction sequence or empty string
 */
export function findJunction(seq1: string, seq2: string, minHomology: number, maxHomology: number): string {
  const s1 = seq1.toUpperCase();
  const s2 = seq2.toUpperCase();

  const start = Math.max(0, s1.length - maxHomology);
  const end = Math.max(0, s1.length - minHomology);

  for (let i = start; i <= end; i++) {
    const suffix = s1.slice(i);
    if (s2.startsWith(suffix)) {
      return suffix;
    }
  }

  return '';
}

/**
 * Calculate GC content of a sequence
 * @param seq - DNA sequence
 * @returns GC content as percentage (0-100)
 */
export function gcContent(seq: string): number {
  const upper = seq.toUpperCase();
  const gc = (upper.match(/[GC]/g) || []).length;
  return (gc / upper.length) * 100;
}

/**
 * Generate FASTA format string
 * @param id - Sequence ID
 * @param seq - Sequence
 * @param options - Options
 * @returns FASTA formatted string
 */
export function toFasta(id: string, seq: string, options: FastaOptions = {}): string {
  const { lineWidth = 60, circular = false } = options;

  const header = `>${id}${circular ? ' circular' : ''}`;
  const lines: string[] = [];

  for (let i = 0; i < seq.length; i += lineWidth) {
    lines.push(seq.slice(i, i + lineWidth));
  }

  return `${header}\n${lines.join('\n')}`;
}
