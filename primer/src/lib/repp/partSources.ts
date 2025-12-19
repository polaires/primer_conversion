/**
 * Part Sources Configuration
 * Support for Addgene and iGEM part repositories
 * Based on: https://github.com/Lattice-Automation/repp
 */

import { parseFasta } from './sequence.js';

// ============================================================================
// TYPES
// ============================================================================

interface PartSource {
  id: string;
  name: string;
  description: string;
  url: string;
  cost: number;
  databaseUrl: string;
  color: string;
  icon: string;
}

interface PartFragment {
  id: string;
  seq: string;
  cost: number;
  source: string;
  url: string;
  description?: string;
}

interface ParsedPartId {
  source: string;
  id: string;
  url: string;
}

interface SearchOptions {
  caseSensitive?: boolean;
  maxResults?: number;
}

interface CachedData {
  data: PartFragment[];
  timestamp: number;
  version: string;
}

// ============================================================================
// CONFIGURATION
// ============================================================================

/**
 * Part source definitions
 * Costs based on REPP paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0223935
 */
export const PART_SOURCES: Record<string, PartSource> = {
  addgene: {
    id: 'addgene',
    name: 'Addgene',
    description: 'Nonprofit plasmid repository with over 100,000 plasmids',
    url: 'https://www.addgene.org',
    cost: 65.0, // $65 per plasmid
    databaseUrl: 'https://repp.s3.amazonaws.com/addgene.fa.gz',
    color: '#0073b1',
    icon: 'A',
  },
  igem: {
    id: 'igem',
    name: 'iGEM Registry',
    description: 'Registry of Standard Biological Parts with over 70,000 parts',
    url: 'https://parts.igem.org',
    cost: 0.0, // Free for iGEM participants
    databaseUrl: 'https://repp.s3.amazonaws.com/igem.fa.gz',
    color: '#009245',
    icon: 'i',
  },
  dnasu: {
    id: 'dnasu',
    name: 'DNASU',
    description: 'DNA repository from Arizona State University',
    url: 'https://dnasu.org',
    cost: 55.0, // $55 per plasmid
    databaseUrl: 'https://repp.s3.amazonaws.com/dnasu.fa.gz',
    color: '#8C1D40',
    icon: 'D',
  },
};

/**
 * Sample parts for demo/testing purposes when databases cannot be loaded
 * These are well-known parts from each registry
 */
export const SAMPLE_PARTS: Record<string, PartFragment[]> = {
  addgene: [
    { id: 'pUC19 (Addgene #50005)', seq: '', cost: 65, source: 'addgene', url: 'https://www.addgene.org/50005/' },
    { id: 'pET-28a(+) (Addgene #69864)', seq: '', cost: 65, source: 'addgene', url: 'https://www.addgene.org/69864/' },
    { id: 'pCMV-dR8.2 dvpr (Addgene #8455)', seq: '', cost: 65, source: 'addgene', url: 'https://www.addgene.org/8455/' },
  ],
  igem: [
    { id: 'BBa_J23100', seq: 'TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC', cost: 0, source: 'igem', url: 'https://parts.igem.org/Part:BBa_J23100', description: 'Constitutive promoter' },
    { id: 'BBa_B0034', seq: 'AAAGAGGAGAAA', cost: 0, source: 'igem', url: 'https://parts.igem.org/Part:BBa_B0034', description: 'RBS' },
    { id: 'BBa_E0040', seq: '', cost: 0, source: 'igem', url: 'https://parts.igem.org/Part:BBa_E0040', description: 'GFP' },
    { id: 'BBa_B0010', seq: 'CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTC', cost: 0, source: 'igem', url: 'https://parts.igem.org/Part:BBa_B0010', description: 'Terminator T1' },
    { id: 'BBa_B0012', seq: 'TCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA', cost: 0, source: 'igem', url: 'https://parts.igem.org/Part:BBa_B0012', description: 'Terminator TE' },
  ],
};

// ============================================================================
// MAIN FUNCTIONS
// ============================================================================

/**
 * Get all available part sources
 */
export function getAvailablePartSources(): PartSource[] {
  return Object.values(PART_SOURCES);
}

/**
 * Get part source by ID
 */
export function getPartSource(sourceId: string): PartSource | null {
  return PART_SOURCES[sourceId] || null;
}

/**
 * Parse part ID to extract source information
 * Addgene parts: "Addgene #12345" or "pXXX (Addgene #12345)"
 * iGEM parts: "BBa_XXXXX" prefix
 */
export function parsePartId(partId: string): ParsedPartId {
  // Check for Addgene format
  const addgeneMatch = partId.match(/Addgene[#\s]+(\d+)/i);
  if (addgeneMatch) {
    return {
      source: 'addgene',
      id: addgeneMatch[1],
      url: `https://www.addgene.org/${addgeneMatch[1]}/`,
    };
  }

  // Check for iGEM BBa_ format
  const igemMatch = partId.match(/^(BBa_[A-Za-z0-9_]+)/i);
  if (igemMatch) {
    return {
      source: 'igem',
      id: igemMatch[1],
      url: `https://parts.igem.org/Part:${igemMatch[1]}`,
    };
  }

  // Check for iGEM pSB format (common backbones)
  const psbMatch = partId.match(/^(pSB\d+[A-Za-z]\d*)/i);
  if (psbMatch) {
    return {
      source: 'igem',
      id: psbMatch[1],
      url: `https://parts.igem.org/Part:${psbMatch[1]}`,
    };
  }

  return { source: 'custom', id: partId, url: '' };
}

/**
 * Fetch part database from REPP S3 bucket
 */
export async function fetchPartDatabase(
  sourceId: string,
  onProgress: ((loaded: number, total: number) => void) | null = null
): Promise<PartFragment[]> {
  const source = PART_SOURCES[sourceId];
  if (!source) {
    throw new Error(`Unknown part source: ${sourceId}`);
  }

  try {
    const response = await fetch(source.databaseUrl);
    if (!response.ok) {
      // Provide helpful error message for common HTTP errors
      if (response.status === 403) {
        throw new Error(
          `${source.name} database is temporarily unavailable (access denied). ` +
          `Please upload a custom FASTA database instead, or try again later.`
        );
      } else if (response.status === 404) {
        throw new Error(
          `${source.name} database not found. The database URL may have changed. ` +
          `Please upload a custom FASTA database instead.`
        );
      }
      throw new Error(`Failed to fetch ${source.name} database: HTTP ${response.status}`);
    }

    // Handle gzipped content
    const blob = await response.blob();
    const ds = new DecompressionStream('gzip');
    const decompressedStream = blob.stream().pipeThrough(ds);
    const decompressedBlob = await new Response(decompressedStream).blob();
    const text = await decompressedBlob.text();

    // Parse FASTA and add source metadata
    const fragments = parseFasta(text);
    return fragments.map(frag => ({
      ...frag,
      cost: source.cost,
      source: sourceId,
      url: generatePartUrl(frag.id, sourceId),
    }));
  } catch (error: any) {
    // Handle network errors with helpful message
    if (error.name === 'TypeError' && error.message === 'Failed to fetch') {
      throw new Error(
        `Unable to connect to ${source.name} database. Check your internet connection ` +
        `or upload a custom FASTA database instead.`
      );
    }
    console.error(`Error fetching ${source.name} database:`, error);
    throw error;
  }
}

/**
 * Generate URL for a part based on its ID and source
 */
export function generatePartUrl(partId: string, sourceId: string): string {
  const parsed = parsePartId(partId);
  if (parsed.url) {
    return parsed.url;
  }

  switch (sourceId) {
    case 'addgene':
      // Try to extract Addgene ID from part name
      const addgeneNum = partId.match(/\d{4,}/);
      return addgeneNum ? `https://www.addgene.org/${addgeneNum[0]}/` : 'https://www.addgene.org';
    case 'igem':
      return `https://parts.igem.org/Part:${partId}`;
    case 'dnasu':
      return `https://dnasu.org/DNASU/AdvancedSearchOptions.do?search=${partId}`;
    default:
      return '';
  }
}

/**
 * Search parts by name/ID in a database
 */
export function searchParts(
  database: PartFragment[],
  query: string,
  options: SearchOptions = {}
): PartFragment[] {
  const { caseSensitive = false, maxResults = 100 } = options;

  if (!query || query.length < 2) {
    return [];
  }

  const normalizedQuery = caseSensitive ? query : query.toLowerCase();

  const results = database.filter(frag => {
    const id = caseSensitive ? frag.id : frag.id.toLowerCase();
    const description = caseSensitive ? (frag.description || '') : (frag.description || '').toLowerCase();
    return id.includes(normalizedQuery) || description.includes(normalizedQuery);
  });

  return results.slice(0, maxResults);
}

/**
 * Load part database from local storage cache
 */
export function loadCachedDatabase(sourceId: string): PartFragment[] | null {
  try {
    const cacheKey = `repp_db_${sourceId}`;
    const cached = localStorage.getItem(cacheKey);
    if (!cached) return null;

    const { data, timestamp, version }: CachedData = JSON.parse(cached);

    // Cache expires after 7 days
    const CACHE_TTL = 7 * 24 * 60 * 60 * 1000;
    if (Date.now() - timestamp > CACHE_TTL) {
      localStorage.removeItem(cacheKey);
      return null;
    }

    return data;
  } catch (error) {
    console.error('Error loading cached database:', error);
    return null;
  }
}

/**
 * Save part database to local storage cache
 */
export function saveCachedDatabase(sourceId: string, fragments: PartFragment[]): void {
  try {
    const cacheKey = `repp_db_${sourceId}`;
    const cacheData: CachedData = {
      data: fragments,
      timestamp: Date.now(),
      version: '1.0',
    };
    localStorage.setItem(cacheKey, JSON.stringify(cacheData));
  } catch (error) {
    // Storage might be full or unavailable
    console.warn('Could not cache database:', error);
  }
}

/**
 * Clear cached database
 */
export function clearCachedDatabase(sourceId: string = 'all'): void {
  if (sourceId === 'all') {
    Object.keys(PART_SOURCES).forEach(id => {
      localStorage.removeItem(`repp_db_${id}`);
    });
  } else {
    localStorage.removeItem(`repp_db_${sourceId}`);
  }
}

/**
 * Merge multiple databases into one, tagging each fragment with its source
 */
export function mergeDatabases(databases: Record<string, PartFragment[]>): PartFragment[] {
  const merged: PartFragment[] = [];

  for (const [sourceId, fragments] of Object.entries(databases)) {
    const source = PART_SOURCES[sourceId] || { cost: 0 };

    fragments.forEach(frag => {
      merged.push({
        ...frag,
        source: sourceId,
        cost: frag.cost ?? source.cost,
        url: frag.url || generatePartUrl(frag.id, sourceId),
      });
    });
  }

  return merged;
}

export default {
  PART_SOURCES,
  SAMPLE_PARTS,
  getAvailablePartSources,
  getPartSource,
  parsePartId,
  fetchPartDatabase,
  generatePartUrl,
  searchParts,
  loadCachedDatabase,
  saveCachedDatabase,
  clearCachedDatabase,
  mergeDatabases,
};
