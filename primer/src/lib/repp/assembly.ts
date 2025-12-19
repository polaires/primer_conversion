/**
 * Assembly Module
 * Core algorithm for DNA assembly planning - finds optimal fragment combinations
 * Based on: https://github.com/Lattice-Automation/repp
 */

import { DEFAULT_CONFIG, synthFragmentCost, primerCost, Config } from './config.js';
import { findMatches, findJunction, parseFasta, Match } from './sequence.js';
import { primers as designPrimersLib, PrimerDesignResult } from '../primers.js';
import { tmCache, gcCache } from '../tm.js';

/**
 * Fragment types for assembly
 */
export const FragType = {
  LINEAR: 'linear',
  CIRCULAR: 'plasmid',
  PCR: 'pcr',
  SYNTHETIC: 'synthetic',
} as const;

export type FragmentType = typeof FragType[keyof typeof FragType];

/**
 * Primer information
 */
export interface Primer {
  seq: string;
  tm: number;
  gc: number;
  strand: boolean;
}

/**
 * Fragment object
 */
export interface Fragment {
  id: string;
  seq: string;
  type: FragmentType;
  start: number;
  end: number;
  cost: number;
  url: string;
  primers: [Primer, Primer] | null;
  pcrSeq: string | null;
  assemblies?: Assembly[];
}

/**
 * Assembly object
 */
export interface Assembly {
  frags: Fragment[];
  cost: number;
  synths: number;
}

/**
 * Solution object
 */
export interface Solution {
  count: number;
  cost: number;
  fragments: Fragment[];
}

/**
 * Assembly planning result
 */
export interface AssemblyPlanResult {
  target: string;
  seq: string;
  solutions: Solution[];
}

/**
 * Database fragment
 */
export interface DatabaseFragment {
  id: string;
  seq: string;
  cost?: number;
  url?: string;
  circular?: boolean;
}

/**
 * Fragment creation parameters
 */
interface CreateFragmentParams {
  id: string;
  seq: string;
  type?: FragmentType;
  start?: number;
  end?: number | null;
  cost?: number;
  url?: string;
  primers?: [Primer, Primer] | null;
}

/**
 * Create a fragment object
 */
export function createFragment(params: CreateFragmentParams): Fragment {
  const {
    id,
    seq,
    type = FragType.LINEAR,
    start = 0,
    end = null,
    cost = 0,
    url = '',
    primers = null,
  } = params;

  return {
    id,
    seq: seq.toUpperCase(),
    type,
    start,
    end: end !== null ? end : start + seq.length - 1,
    cost,
    url,
    primers,
    pcrSeq: null,
  };
}

/**
 * Calculate fragment cost based on type and preparation
 */
export function fragmentCost(frag: Fragment, config: Config = DEFAULT_CONFIG, includeProcurement: boolean = true): number {
  let cost = 0;

  if (includeProcurement && frag.cost) {
    cost += frag.cost;
  }

  if (frag.type === FragType.PCR && frag.primers) {
    cost += primerCost(frag.primers[0].seq, frag.primers[1].seq, config);
  } else if (frag.type === FragType.SYNTHETIC) {
    cost += synthFragmentCost(frag.seq.length, config);
  }

  return cost;
}

/**
 * Check if fragment can overlap with another via existing homology
 */
function overlapsViaHomology(frag1: Fragment, frag2: Fragment, config: Config, targetLength: number | null = null): boolean {
  // Calculate gap between fragments
  let gap = frag2.start - frag1.end;

  // Handle circular wrapping: if frag1 is near the end and frag2 is at the start
  // This happens when frag1.end is near targetLength and frag2.start is near 0
  if (targetLength && gap < -(targetLength / 2)) {
    // Circular case: gap should be calculated as wrapping around
    // Actual gap = distance from frag1.end to end of target + frag2.start
    gap = (targetLength - frag1.end - 1) + frag2.start;
  }

  return gap <= -config.fragmentsMinHomology;
}

/**
 * Check if fragment can overlap with another via PCR
 */
function overlapsViaPCR(frag1: Fragment, frag2: Fragment, config: Config): boolean {
  return (frag2.start - frag1.end) <= config.pcrPrimerMaxEmbedLength;
}

/**
 * Calculate number of synthetic fragments needed between two fragments
 */
function synthDist(frag1: Fragment, frag2: Fragment, config: Config): number {
  if (overlapsViaPCR(frag1, frag2, config)) {
    return 0;
  }

  const dist = Math.max(1, frag2.start - frag1.end);
  return Math.ceil(dist / config.syntheticMaxLength);
}

/**
 * Estimate cost to connect two fragments
 */
function costTo(frag1: Fragment, frag2: Fragment, config: Config): number {
  const pcrNoHomology = 50 * config.pcrBpCost;
  const pcrHomology = (50 + config.fragmentsMinHomology) * config.pcrBpCost;

  if (frag1 === frag2 || frag1.id === frag2.id) {
    const needsPCR = frag1.type === FragType.PCR || frag1.type === FragType.CIRCULAR;
    return needsPCR ? pcrNoHomology : 0;
  }

  if (overlapsViaPCR(frag1, frag2, config)) {
    return overlapsViaHomology(frag1, frag2, config) ? pcrNoHomology : pcrHomology;
  }

  // Need synthetic fragment
  const dist = frag2.start - frag1.end + config.fragmentsMinHomology * 2;
  const synthCost = synthFragmentCost(dist, config);

  const needsPCR = frag1.type === FragType.PCR || frag1.type === FragType.CIRCULAR;
  return needsPCR ? synthCost + pcrNoHomology : synthCost;
}

/**
 * Create assemblies from fragments using DAG traversal
 */
export function createAssemblies(fragments: Fragment[], target: string, config: Config = DEFAULT_CONFIG): Assembly[] {
  const targetLength = target.length;
  const maxFragments = config.fragmentsMaxCount;

  // Sort fragments by start position
  const sortedFrags = [...fragments].sort((a, b) => a.start - b.start);

  // Initialize assemblies on each fragment
  for (const frag of sortedFrags) {
    frag.assemblies = [{
      frags: [frag],
      cost: costTo(frag, frag, config),
      synths: 0,
    }];
  }

  const completedAssemblies: Assembly[] = [];

  // Build assemblies using dynamic programming
  for (let i = 0; i < sortedFrags.length; i++) {
    const frag = sortedFrags[i];

    // Find fragments this one can reach
    const reachable = findReachable(sortedFrags, i, config);

    for (const j of reachable) {
      const nextFrag = sortedFrags[j];

      for (const assembly of frag.assemblies!) {
        // Try to add nextFrag to this assembly
        const result = tryAddFragment(assembly, nextFrag, maxFragments, targetLength, config);

        if (!result.created) continue;

        if (result.circularized) {
          completedAssemblies.push(result.assembly);
        } else {
          if (!nextFrag.assemblies) nextFrag.assemblies = [];
          nextFrag.assemblies.push(result.assembly);
        }
      }
    }
  }

  // Add fully synthetic solution as fallback
  const synthSolution = createSyntheticAssembly(target, config);
  completedAssemblies.push(synthSolution);

  return completedAssemblies;
}

/**
 * Find reachable fragments from a given position
 */
function findReachable(frags: Fragment[], i: number, config: Config): number[] {
  const reachable: number[] = [];
  const current = frags[i];

  for (let j = i + 1; j < frags.length; j++) {
    const next = frags[j];

    // Skip if fully contained
    if (next.end < current.end) continue;

    reachable.push(j);

    // Stop if we've found the same fragment wrapping around
    if (next.id === current.id) break;
  }

  return reachable;
}

/**
 * Result of trying to add a fragment
 */
interface TryAddResult {
  created: boolean;
  circularized?: boolean;
  assembly?: Assembly;
}

/**
 * Try to add a fragment to an assembly
 */
function tryAddFragment(assembly: Assembly, frag: Fragment, maxFragments: number, targetLength: number, config: Config): TryAddResult {
  const firstFrag = assembly.frags[0];
  const lastFrag = assembly.frags[assembly.frags.length - 1];

  // Check if this would complete the circle
  const circularized = frag.end >= firstFrag.start + targetLength - 1;
  const selfAnnealing = frag.id === firstFrag.id;

  // Calculate synthetic fragments needed
  const synths = synthDist(lastFrag, frag, config);
  const newCount = assembly.frags.length + synths + (selfAnnealing ? 0 : 1);

  // Check if within limits
  if (newCount > maxFragments) {
    return { created: false };
  }

  // Calculate cost
  let annealCost = costTo(lastFrag, frag, config);
  if (selfAnnealing && synths === 0) {
    annealCost = 0;
  }

  // Check if fragment already in assembly (don't count procurement twice)
  const alreadyIncluded = assembly.frags.some(f => f.id === frag.id && f.type === frag.type);
  annealCost += fragmentCost(frag, config, !alreadyIncluded);

  // Create new assembly
  const newFrags = [...assembly.frags.map(f => ({ ...f }))];
  if (!selfAnnealing) {
    newFrags.push({ ...frag });
  }

  return {
    created: true,
    circularized,
    assembly: {
      frags: newFrags,
      cost: assembly.cost + annealCost,
      synths: assembly.synths + synths,
    },
  };
}

/**
 * Create a fully synthetic assembly as fallback
 */
function createSyntheticAssembly(target: string, config: Config): Assembly {
  // For circular assemblies, we need at least 2 fragments to create proper junctions
  // Even if the target is small, we split it to ensure overlapping homology regions
  const minFrags = 2;
  const numFrags = Math.max(minFrags, Math.ceil(target.length / config.syntheticMaxLength));

  // Calculate fragment length with overlap regions
  // Each fragment needs homology on both ends for Gibson assembly
  const homology = config.fragmentsMinHomology;
  const totalCoverage = target.length + numFrags * homology; // Account for overlaps
  const fragLength = Math.ceil(totalCoverage / numFrags);

  const frags: Fragment[] = [];
  let pos = 0;

  for (let i = 0; i < numFrags; i++) {
    // Start with overlap from previous fragment (except first)
    const start = i === 0 ? 0 : pos - homology;

    // End position - ensure we don't exceed target length
    let end = Math.min(target.length, start + fragLength);

    // For the last fragment, extend to include homology with the first fragment
    // by wrapping around to the beginning of the target sequence
    let seq;
    if (i === numFrags - 1) {
      // Last fragment: include homology to first fragment
      seq = target.slice(start, target.length) + target.slice(0, homology);
      end = target.length + homology; // Virtual end position for tracking
    } else {
      seq = target.slice(start, end);
    }

    frags.push({
      id: `synthetic-${i + 1}`,
      seq,
      type: FragType.SYNTHETIC,
      start,
      end: i === numFrags - 1 ? target.length - 1 : end - 1,
      cost: synthFragmentCost(seq.length, config),
      url: '',
      primers: null,
      pcrSeq: null,
    });

    pos = end;
  }

  return {
    frags,
    cost: frags.reduce((sum, f) => sum + f.cost, 0),
    synths: numFrags,
  };
}

/**
 * Group assemblies result
 */
interface GroupAssembliesResult {
  counts: number[];
  byCount: Map<number, Assembly[]>;
}

/**
 * Group assemblies by fragment count and sort by cost
 */
export function groupAssemblies(assemblies: Assembly[]): GroupAssembliesResult {
  const byCount = new Map<number, Assembly[]>();

  for (const assembly of assemblies) {
    const count = assembly.frags.length + assembly.synths;
    if (!byCount.has(count)) {
      byCount.set(count, []);
    }
    byCount.get(count)!.push(assembly);
  }

  // Sort each group by cost
  for (const [count, group] of byCount) {
    group.sort((a, b) => a.cost - b.cost);
  }

  // Get sorted counts
  const counts = [...byCount.keys()].sort((a, b) => a - b);

  return { counts, byCount };
}

/**
 * Fill in assemblies with primers and synthetic fragments
 */
export function fillAssemblies(assemblies: Assembly[], target: string, config: Config = DEFAULT_CONFIG): Solution[] {
  const { counts, byCount } = groupAssemblies(assemblies);

  const solutions: Solution[] = [];
  let minCost = Infinity;

  for (const count of counts) {
    for (const assembly of byCount.get(count)!) {
      // Skip if already more expensive than best
      if (assembly.cost > minCost) break;

      try {
        const filled = fillAssembly(assembly, target, config);
        if (!filled) continue;

        const totalCost = calculateTotalCost(filled, config);
        if (totalCost >= minCost) continue;

        minCost = totalCost;
        solutions.push({
          count: filled.length,
          cost: Math.round(totalCost * 100) / 100,
          fragments: filled,
        });
      } catch (e) {
        // Assembly failed (e.g., can't design primers)
        continue;
      }
    }
  }

  // Sort by fragment count
  solutions.sort((a, b) => a.count - b.count);

  return solutions;
}

/**
 * Fill a single assembly with primers and synthetic fragments
 */
function fillAssembly(assembly: Assembly, target: string, config: Config, debug: boolean = false): Fragment[] | null {
  const filled: Fragment[] = [];
  const frags = assembly.frags;
  const targetLength = target.length;

  for (let i = 0; i < frags.length; i++) {
    const frag = frags[i];
    const prev = frags[(i - 1 + frags.length) % frags.length];
    const next = frags[(i + 1) % frags.length];

    // Check if primers are needed (pass targetLength for circular assembly handling)
    const needsPrimers = frag.type === FragType.PCR ||
                        frag.type === FragType.CIRCULAR ||
                        !overlapsViaHomology(prev, frag, config, targetLength) ||
                        !overlapsViaHomology(frag, next, config, targetLength);

    if (needsPrimers && frag.type !== FragType.SYNTHETIC) {
      // Design primers for this fragment
      const withPrimers = designFragmentPrimers(frag, prev, next, target, config);
      if (!withPrimers) {
        if (debug) console.log(`  fillAssembly: primer design failed for ${frag.id}`);
        return null;
      }
      filled.push(withPrimers);
    } else {
      filled.push({ ...frag });
    }

    // Check if synthetic fragments needed to reach next
    const synthFrags = createSyntheticBridge(frag, next, target, config);
    if (synthFrags) {
      filled.push(...synthFrags);
    }
  }

  // Validate junctions
  if (!validateJunctions(filled, config)) {
    if (debug) {
      console.log(`  fillAssembly: junction validation failed`);
      for (let i = 0; i < filled.length; i++) {
        const f = filled[i];
        const n = filled[(i + 1) % filled.length];
        const fSeq = f.pcrSeq || f.seq;
        const nSeq = n.pcrSeq || n.seq;
        const junction = findJunction(fSeq, nSeq, config.fragmentsMinHomology, config.fragmentsMaxHomology);
        console.log(`    Junction ${f.id} -> ${n.id}: ${junction ? junction.length + 'bp' : 'NONE'}`);
      }
    }
    return null;
  }

  return filled;
}

/**
 * Design primers for a fragment
 */
function designFragmentPrimers(frag: Fragment, prev: Fragment, next: Fragment, target: string, config: Config): Fragment | null {
  try {
    // Get the fragment sequence
    const fragSeq = frag.seq || target.slice(frag.start, frag.end + 1);

    // Use the primers library to design primers
    const [fwdResult, revResult] = designPrimersLib(fragSeq, {
      tmTarget: 60,
      tmRange: 5,
    });

    if (!fwdResult || !revResult || !fwdResult.seq || !revResult.seq) {
      return null;
    }

    // Start with the designed primers
    let fwdPrimer = fwdResult.seq;
    let revPrimer = revResult.seq;
    let fwdTm = fwdResult.tm;
    let revTm = revResult.tm;

    // Track the homology tails added to compute the full PCR product
    let fwdTailLen = 0;
    let revTailLen = 0;

    const targetLength = target.length;

    // Add homology tails if needed
    if (!overlapsViaHomology(prev, frag, config, targetLength)) {
      // Add homology to previous fragment
      const homologyLen = config.fragmentsMinHomology;

      // Handle circular case: if prev.end is near the end and we need to get homology
      let homology: string;
      if (prev.end >= targetLength - homologyLen) {
        // Previous fragment ends near the circular boundary
        // Get homology from the end of the sequence, wrapping around if needed
        const fromEnd = target.slice(Math.max(0, prev.end - homologyLen + 1), prev.end + 1);
        homology = fromEnd;
      } else {
        const prevEnd = Math.max(0, prev.end - homologyLen + 1);
        homology = target.slice(prevEnd, prev.end + 1);
      }
      fwdPrimer = homology + fwdPrimer;
      fwdTailLen = homology.length;
    }

    if (!overlapsViaHomology(frag, next, config, targetLength)) {
      // Add homology to next fragment
      const homologyLen = config.fragmentsMinHomology;

      // Handle circular case: if frag ends near the end and next starts near 0
      let homology: string;
      if (frag.end >= targetLength - config.fragmentsMinHomology && next.start < homologyLen) {
        // This is the circular junction - next fragment starts at the beginning
        // Get homology from the start of the sequence
        homology = target.slice(0, homologyLen);
      } else {
        const nextStart = next.start % targetLength;
        homology = target.slice(nextStart, nextStart + homologyLen);
      }
      revPrimer = reverseComplementDNA(homology) + revPrimer;
      revTailLen = homology.length;
    }

    // Compute the full PCR product sequence including the homology tails
    // The PCR product = fwdTail + fragSeq + revTail (reverse complement of the rev primer tail)
    let pcrProductSeq = fragSeq;

    // Add 5' homology tail (from fwd primer)
    if (fwdTailLen > 0) {
      // The tail sequence is at the beginning of the fwd primer
      const fwdTail = fwdPrimer.slice(0, fwdTailLen);
      pcrProductSeq = fwdTail + pcrProductSeq;
    }

    // Add 3' homology tail (from rev primer - need to get the complement)
    if (revTailLen > 0) {
      // The tail on the rev primer corresponds to the 3' end of the product
      // Handle circular junction case
      let revTail: string;
      if (frag.end >= targetLength - config.fragmentsMinHomology && next.start < revTailLen) {
        // Circular junction - get tail from the beginning of the sequence
        revTail = target.slice(0, revTailLen);
      } else {
        const nextStart = next.start % targetLength;
        revTail = target.slice(nextStart, nextStart + revTailLen);
      }
      pcrProductSeq = pcrProductSeq + revTail;
    }

    return {
      ...frag,
      type: FragType.PCR,
      primers: [
        { seq: fwdPrimer, tm: fwdTm, gc: fwdResult.gc, strand: true },
        { seq: revPrimer, tm: revTm, gc: revResult.gc, strand: false },
      ],
      pcrSeq: pcrProductSeq,
    };
  } catch (e) {
    return null;
  }
}

/**
 * Create synthetic fragments to bridge a gap
 */
function createSyntheticBridge(frag1: Fragment, frag2: Fragment, target: string, config: Config): Fragment[] | null {
  const numSynths = synthDist(frag1, frag2, config);
  if (numSynths === 0) return null;

  const synths: Fragment[] = [];
  const gap = frag2.start - frag1.end;
  const fragLen = Math.ceil(gap / numSynths) + config.fragmentsMinHomology * 2;

  let pos = frag1.end - config.fragmentsMinHomology + 1;

  for (let i = 0; i < numSynths; i++) {
    const start = pos;
    const end = Math.min(start + fragLen, frag2.start + config.fragmentsMinHomology);
    const seq = target.slice(start, end);

    synths.push({
      id: `${frag1.id}-${frag2.id}-synth-${i + 1}`,
      seq,
      type: FragType.SYNTHETIC,
      start,
      end: end - 1,
      cost: synthFragmentCost(seq.length, config),
      url: '',
      primers: null,
      pcrSeq: null,
    });

    pos = end - config.fragmentsMinHomology;
  }

  return synths;
}

/**
 * Validate junctions between fragments
 */
function validateJunctions(frags: Fragment[], config: Config): boolean {
  for (let i = 0; i < frags.length; i++) {
    const frag = frags[i];
    const next = frags[(i + 1) % frags.length];

    const fragSeq = frag.pcrSeq || frag.seq;
    const nextSeq = next.pcrSeq || next.seq;

    const junction = findJunction(fragSeq, nextSeq, config.fragmentsMinHomology, config.fragmentsMaxHomology);

    if (!junction || junction.length < config.fragmentsMinHomology) {
      return false;
    }
  }

  return true;
}

/**
 * Calculate total cost of filled fragments
 */
function calculateTotalCost(frags: Fragment[], config: Config): number {
  const seen = new Set<string>();
  let total = 0;
  let hasGibson = false;
  let hasPCR = false;

  for (const frag of frags) {
    const includeProcurement = !seen.has(frag.id);
    total += fragmentCost(frag, config, includeProcurement);
    seen.add(frag.id);

    if (frag.type !== FragType.LINEAR && frag.type !== FragType.CIRCULAR) {
      hasGibson = true;
    }
    if (frag.type === FragType.PCR) {
      hasPCR = true;
    }
  }

  if (hasGibson) {
    total += config.gibsonAssemblyCost + config.gibsonAssemblyTimeCost;
  }
  if (hasPCR) {
    total += config.pcrTimeCost;
  }

  return total;
}

/**
 * Reverse complement DNA sequence
 */
function reverseComplementDNA(seq: string): string {
  const comp: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
  return seq.split('').reverse().map(c => comp[c] || c).join('');
}

/**
 * Planning options
 */
interface PlanAssemblyOptions {
  config?: Partial<Config>;
  targetName?: string;
}

/**
 * Main assembly planning function
 */
export function planAssembly(targetSeq: string, database: DatabaseFragment[], options: PlanAssemblyOptions = {}): AssemblyPlanResult {
  const config = { ...DEFAULT_CONFIG, ...options.config };
  const target = targetSeq.toUpperCase();

  // Find matches from database
  const allMatches: Fragment[] = [];

  for (const dbFrag of database) {
    const matches = findMatches(target, dbFrag.seq, {
      minLength: config.pcrMinLength,
      circular: true,
    });

    for (const match of matches) {
      allMatches.push(createFragment({
        id: dbFrag.id,
        seq: match.seq,
        type: dbFrag.circular ? FragType.CIRCULAR : FragType.PCR,
        start: match.queryStart,
        end: match.queryEnd,
        cost: dbFrag.cost || 0,
        url: dbFrag.url || '',
      }));
    }
  }

  if (allMatches.length === 0) {
    // No matches - return synthetic-only solution
    const synthAssembly = createSyntheticAssembly(target, config);
    return {
      target: options.targetName || 'target',
      seq: target,
      solutions: [{
        count: synthAssembly.frags.length,
        cost: synthAssembly.cost,
        fragments: synthAssembly.frags,
      }],
    };
  }

  // Create and fill assemblies
  const assemblies = createAssemblies(allMatches, target, config);
  const solutions = fillAssemblies(assemblies, target, config);

  return {
    target: options.targetName || 'target',
    seq: target,
    solutions,
  };
}
