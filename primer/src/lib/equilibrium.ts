/**
 * Equilibrium Efficiency Calculation for Primer Scoring
 *
 * Based on the Pythia algorithm (Mann et al., 2009, Nucleic Acids Research)
 * "A thermodynamic approach to PCR primer design"
 *
 * Core principle: Calculate equilibrium concentrations of all chemical species
 * and determine what fraction of primer is bound to the intended target.
 *
 * Equilibrium Efficiency = min(η_fwd, η_rev)
 * where η = [Primer bound to target] / [Total Primer]
 *
 * Species modeled:
 * - Primer (unfolded) - free primer in solution
 * - Primer (hairpin) - self-folded primer
 * - Primer-Primer (homodimer) - self-dimer
 * - Primer-Template (target) - desired binding
 * - Primer-Template (off-target) - mispriming elsewhere
 * - Fwd-Rev (heterodimer) - cross-dimer between primers
 */

import { dg as foldDG } from './fold.js';
import { tm } from './tm.js';
import { DNA24_ENERGIES, DNA24_COMPLEMENT } from './dna24.js';
import { reverseComplement } from './sequenceUtils.js';

// Gas constant in kcal/(mol·K)
const R = 1.987e-3;

// Default concentrations (M)
const DEFAULT_PRIMER_CONC = 0.5e-6;  // 500 nM (typical PCR)
const DEFAULT_TEMPLATE_CONC = 1e-9;  // 1 nM (late-stage PCR assumption)

// ============================================================================
// PERFORMANCE OPTIMIZATION: Global DG Cache
// These calculations are expensive (O(n²) to O(n⁴)) and the same sequences
// are often scored multiple times during exhaustive search.
// Cache key format: "seq:temp" for single-seq, "seq1:seq2:temp" for pairs
// ============================================================================
const hairpinDGCache = new Map<string, number>();
const homodimerDGCache = new Map<string, number>();
const heterodimerDGCache = new Map<string, number>();

// TypeScript Interfaces

export interface Primer {
  seq?: string;
  sequence?: string;
  bindingSite?: string;
}

export interface EquilibriumOptions {
  temperature?: number;
  primerConc?: number;
  templateConc?: number;
  includeOffTarget?: boolean;
}

export interface CacheStats {
  hairpinHits: number;
  hairpinMisses: number;
  homodimerHits: number;
  homodimerMisses: number;
  heterodimerHits: number;
  heterodimerMisses: number;
}

export interface DGCacheStats extends CacheStats {
  hairpinCacheSize: number;
  homodimerCacheSize: number;
  heterodimerCacheSize: number;
  totalCacheSize: number;
}

export interface SpeciesDG {
  fwdHairpin: number;
  revHairpin: number;
  fwdHomodimer: number;
  revHomodimer: number;
  fwdTarget: number;
  revTarget: number;
  fwdOffTarget: number;
  revOffTarget: number;
  heterodimer: number;
}

export interface EquilibriumConstants {
  [key: string]: number;
}

export interface Concentrations {
  fwdFree: number;
  fwdHairpin: number;
  fwdHomodimer: number;
  fwdBoundTarget: number;
  fwdOffTarget: number;
  revFree: number;
  revHairpin: number;
  revHomodimer: number;
  revBoundTarget: number;
  revOffTarget: number;
  heterodimer: number;
  templateFree: number;
}

export interface LossBreakdown {
  hairpin: number;
  homodimer: number;
  heterodimer: number;
  offTarget: number;
  free: number;
}

export interface EquilibriumResult {
  efficiency: number;
  efficiencyFwd: number;
  efficiencyRev: number;
  bottleneck: 'forward' | 'reverse';
  species: SpeciesDG;
  equilibriumConstants: EquilibriumConstants;
  concentrations: Concentrations;
  losses: {
    fwd: LossBreakdown;
    rev: LossBreakdown;
  };
  quality: 'excellent' | 'good' | 'acceptable' | 'marginal' | 'poor';
}

export interface EfficiencyScoreOptions {
  optimal?: number;
  acceptable?: number;
  steepness?: number;
}

// Cache statistics for debugging
let cacheStats: CacheStats = {
  hairpinHits: 0,
  hairpinMisses: 0,
  homodimerHits: 0,
  homodimerMisses: 0,
  heterodimerHits: 0,
  heterodimerMisses: 0
};

/**
 * Clear all DG caches (call between design sessions)
 */
export function clearDGCache(): void {
  hairpinDGCache.clear();
  homodimerDGCache.clear();
  heterodimerDGCache.clear();
  cacheStats = {
    hairpinHits: 0,
    hairpinMisses: 0,
    homodimerHits: 0,
    homodimerMisses: 0,
    heterodimerHits: 0,
    heterodimerMisses: 0
  };
}

/**
 * Get cache statistics for debugging
 */
export function getDGCacheStats(): DGCacheStats {
  return {
    ...cacheStats,
    hairpinCacheSize: hairpinDGCache.size,
    homodimerCacheSize: homodimerDGCache.size,
    heterodimerCacheSize: heterodimerDGCache.size,
    totalCacheSize: hairpinDGCache.size + homodimerDGCache.size + heterodimerDGCache.size,
  };
}

/**
 * Calculate the equilibrium efficiency for a primer pair
 *
 * @param fwd - Forward primer object with sequence and binding info
 * @param rev - Reverse primer object with sequence and binding info
 * @param template - Full template sequence
 * @param options - Configuration options
 * @returns Equilibrium efficiency results
 */
export function calculateEquilibriumEfficiency(
  fwd: Primer | string,
  rev: Primer | string,
  template: string,
  options: EquilibriumOptions = {}
): EquilibriumResult {
  const {
    temperature = 55,           // Annealing temperature (°C)
    primerConc = DEFAULT_PRIMER_CONC,
    templateConc = DEFAULT_TEMPLATE_CONC,
    includeOffTarget = true,
  } = options;

  const tempK = temperature + 273.15;

  // Calculate ΔG for all species
  const species = calculateAllSpeciesDG(fwd, rev, template, temperature, includeOffTarget);

  // Convert ΔG to equilibrium constants
  const K: EquilibriumConstants = {};
  for (const [key, value] of Object.entries(species)) {
    if (typeof value === 'number' && !isNaN(value)) {
      // K = exp(-ΔG / RT)
      K[key] = Math.exp(-value / (R * tempK));
    }
  }

  // Solve equilibrium concentrations
  const concentrations = solveEquilibrium(K, primerConc, templateConc);

  // Calculate efficiency for each primer
  const η_fwd = concentrations.fwdBoundTarget / primerConc;
  const η_rev = concentrations.revBoundTarget / primerConc;
  const equilibriumEfficiency = Math.min(η_fwd, η_rev);

  // Identify the bottleneck
  const bottleneck: 'forward' | 'reverse' = η_fwd < η_rev ? 'forward' : 'reverse';

  // Calculate loss breakdown
  const fwdLosses: LossBreakdown = {
    hairpin: concentrations.fwdHairpin / primerConc,
    homodimer: concentrations.fwdHomodimer / primerConc,
    heterodimer: concentrations.heterodimer / primerConc / 2,
    offTarget: concentrations.fwdOffTarget / primerConc,
    free: concentrations.fwdFree / primerConc,
  };

  const revLosses: LossBreakdown = {
    hairpin: concentrations.revHairpin / primerConc,
    homodimer: concentrations.revHomodimer / primerConc,
    heterodimer: concentrations.heterodimer / primerConc / 2,
    offTarget: concentrations.revOffTarget / primerConc,
    free: concentrations.revFree / primerConc,
  };

  return {
    efficiency: equilibriumEfficiency,
    efficiencyFwd: η_fwd,
    efficiencyRev: η_rev,
    bottleneck,
    species,
    equilibriumConstants: K,
    concentrations,
    losses: {
      fwd: fwdLosses,
      rev: revLosses,
    },
    // Quality classification based on efficiency
    quality: classifyEfficiency(equilibriumEfficiency),
  };
}

/**
 * Calculate ΔG for all chemical species
 *
 * @param fwd - Forward primer
 * @param rev - Reverse primer
 * @param template - Template sequence
 * @param temperature - Temperature in Celsius
 * @param includeOffTarget - Whether to include off-target analysis
 * @returns ΔG values for all species
 */
export function calculateAllSpeciesDG(
  fwd: Primer | string,
  rev: Primer | string,
  template: string,
  temperature: number = 55,
  includeOffTarget: boolean = true
): SpeciesDG {
  const fwdSeq = typeof fwd === 'string' ? fwd : fwd.seq || fwd.sequence || '';
  const revSeq = typeof rev === 'string' ? rev : rev.seq || rev.sequence || '';

  // Get binding site sequences from template
  const fwdObj = typeof fwd === 'string' ? null : fwd;
  const revObj = typeof rev === 'string' ? null : rev;

  const fwdBindingSite = fwdObj?.bindingSite || template.slice(0, fwdSeq.length);
  const revBindingSite = revObj?.bindingSite || reverseComplement(template.slice(-revSeq.length));

  return {
    // Hairpin formation (self-folding)
    fwdHairpin: calculateHairpinDG(fwdSeq, temperature),
    revHairpin: calculateHairpinDG(revSeq, temperature),

    // Homodimer formation (self-dimer)
    fwdHomodimer: calculateHomodimerDG(fwdSeq, temperature),
    revHomodimer: calculateHomodimerDG(revSeq, temperature),

    // Target binding (desired)
    fwdTarget: calculateDuplexDG(fwdSeq, fwdBindingSite, temperature),
    revTarget: calculateDuplexDG(revSeq, revBindingSite, temperature),

    // Off-target binding (mispriming)
    ...(includeOffTarget ? {
      fwdOffTarget: calculateOffTargetDG(fwdSeq, template, fwdSeq.length, temperature),
      revOffTarget: calculateOffTargetDG(revSeq, template, template.length - revSeq.length, temperature),
    } : {
      fwdOffTarget: 0,
      revOffTarget: 0,
    }),

    // Heterodimer formation (cross-dimer)
    heterodimer: calculateHeterodimerDG(fwdSeq, revSeq, temperature),
  };
}

/**
 * Calculate hairpin ΔG using Zuker fold algorithm (with caching)
 */
export function calculateHairpinDG(seq: string, temperature: number = 55): number {
  if (!seq || seq.length < 6) return 0;

  // Check cache first
  const cacheKey = `${seq}:${temperature}`;
  if (hairpinDGCache.has(cacheKey)) {
    cacheStats.hairpinHits++;
    return hairpinDGCache.get(cacheKey)!;
  }
  cacheStats.hairpinMisses++;

  try {
    // Use the full Zuker algorithm from fold.js
    const dgVal = foldDG(seq, temperature);
    // Return 0 if positive (no stable structure forms)
    const result = Math.min(0, dgVal);
    hairpinDGCache.set(cacheKey, result);
    return result;
  } catch (e) {
    // Fall back to simple estimation
    const result = estimateHairpinDG(seq, temperature);
    hairpinDGCache.set(cacheKey, result);
    return result;
  }
}

/**
 * Simple hairpin estimation when full fold is not available
 */
function estimateHairpinDG(seq: string, temperature: number = 55): number {
  const tempK = temperature + 273.15;
  let minDG = 0;

  // Check for potential hairpin loops (min 4nt loop)
  for (let loopStart = 4; loopStart < seq.length - 4; loopStart++) {
    for (let loopLen = 4; loopLen <= 8 && loopStart + loopLen < seq.length; loopLen++) {
      let dH = 0;
      let dS = 0;
      let stemLen = 0;

      // Loop penalty (Jacobson-Stockmayer)
      dS -= 1.75 * R * 1000 * Math.log(loopLen);
      dH += 0;  // Loop enthalpy is typically small

      // Check stem formation
      for (let i = 0; i < Math.min(loopStart, seq.length - loopStart - loopLen); i++) {
        const base5 = seq[loopStart - 1 - i];
        const base3 = seq[loopStart + loopLen + i];

        if (isComplement(base5, base3)) {
          stemLen++;
          // Add nearest-neighbor contribution
          if (i > 0) {
            const prevBase5 = seq[loopStart - i];
            const prevBase3 = seq[loopStart + loopLen + i - 1];
            const pair = `${prevBase5}${base5}/${complement(prevBase3)}${complement(base3)}`;
            const nn = DNA24_ENERGIES.NN[pair];
            if (nn) {
              dH += nn[0];
              dS += nn[1];
            }
          }
        } else {
          break;
        }
      }

      if (stemLen >= 3) {
        const dG = dH - tempK * (dS / 1000);
        if (dG < minDG) {
          minDG = dG;
        }
      }
    }
  }

  return minDG;
}

/**
 * Calculate homodimer (self-dimer) ΔG
 *
 * A homodimer forms when a sequence binds to another copy of itself.
 * This requires the sequence to be partially self-complementary.
 *
 * Two copies of the same sequence align antiparallel:
 *   Strand 1: 5'-ABCDEF-3'
 *   Strand 2: 3'-ABCDEF-5' (same sequence, flipped)
 *
 * For binding: seq[i] must be complementary to seq[len-1-j] at aligned positions.
 */
export function calculateHomodimerDG(seq: string, temperature: number = 55): number {
  if (!seq || seq.length < 6) return 0;

  seq = seq.toUpperCase();

  // Check cache first
  const cacheKey = `${seq}:${temperature}`;
  if (homodimerDGCache.has(cacheKey)) {
    cacheStats.homodimerHits++;
    return homodimerDGCache.get(cacheKey)!;
  }
  cacheStats.homodimerMisses++;

  const tempK = temperature + 273.15;
  const len = seq.length;

  // seqRev is the sequence reversed (NOT reverse-complemented)
  // This represents the second strand aligned antiparallel
  const seqRev = seq.split('').reverse().join('');
  let minDG = 0;

  // Check all possible alignments of the sequence with its reversed copy
  for (let offset = -(len - 4); offset <= len - 4; offset++) {
    let dH = 0;
    let dS = 0;
    let consecutiveMatches = 0;
    let maxConsecutive = 0;

    for (let i = 0; i < len; i++) {
      const j = i + offset;
      if (j >= 0 && j < len) {
        // Check if seq[i] pairs with seqRev[j]
        // seqRev[j] = seq[len-1-j], so we need seq[i] complementary to seq[len-1-j]
        if (isComplement(seq[i], seqRev[j])) {
          consecutiveMatches++;
          maxConsecutive = Math.max(maxConsecutive, consecutiveMatches);
          // Add nearest-neighbor contribution if we have >= 2 consecutive matches
          if (consecutiveMatches >= 2) {
            const prevI = i - 1;
            const prevJ = j - 1;
            if (prevI >= 0 && prevJ >= 0 && isComplement(seq[prevI], seqRev[prevJ])) {
              // NN pair: strand1 5'->3' paired with strand2 3'->5'
              // seqRev contains the actual bases on the bottom strand
              // Format: "XY/WZ" where XY is top strand, WZ is bottom strand 3'->5'
              const pair = `${seq[prevI]}${seq[i]}/${seqRev[prevJ]}${seqRev[j]}`;
              const nn = DNA24_ENERGIES.NN[pair];
              if (nn) {
                dH += nn[0];
                dS += nn[1];
              }
            }
          }
        } else {
          consecutiveMatches = 0;
        }
      }
    }

    // Only count as homodimer if we have at least 3 consecutive matches
    if (maxConsecutive >= 3 && dH < 0) {
      // Add initiation penalty
      dH += 0.2;  // DNA24 init
      dS += -5.7;
      // Symmetry correction for homodimer
      dS += -1.4;

      const dG = dH - tempK * (dS / 1000);
      if (dG < minDG) {
        minDG = dG;
      }
    }
  }

  // Cache result before returning
  homodimerDGCache.set(cacheKey, minDG);
  return minDG;
}

/**
 * Calculate heterodimer ΔG between two primers
 *
 * A heterodimer forms when two different sequences bind to each other.
 * This requires complementary regions between the sequences.
 *
 * Two different sequences align antiparallel:
 *   Strand 1: 5'-seq1-3'
 *   Strand 2: 3'-seq2-5' (seq2 reversed for antiparallel alignment)
 *
 * For binding: seq1[i] must be complementary to seq2[len2-1-j] at aligned positions.
 */
export function calculateHeterodimerDG(seq1: string, seq2: string, temperature: number = 55): number {
  if (!seq1 || !seq2 || seq1.length < 4 || seq2.length < 4) return 0;

  seq1 = seq1.toUpperCase();
  seq2 = seq2.toUpperCase();

  // Check cache first (use sorted key for symmetry - order shouldn't matter)
  const sortedSeqs = [seq1, seq2].sort();
  const cacheKey = `${sortedSeqs[0]}:${sortedSeqs[1]}:${temperature}`;
  if (heterodimerDGCache.has(cacheKey)) {
    cacheStats.heterodimerHits++;
    return heterodimerDGCache.get(cacheKey)!;
  }
  cacheStats.heterodimerMisses++;

  const tempK = temperature + 273.15;

  // seq2Rev is seq2 reversed (NOT reverse-complemented)
  // This represents seq2 aligned antiparallel
  const seq2Rev = seq2.split('').reverse().join('');
  let minDG = 0;

  // Check all possible alignments
  for (let offset = -(seq1.length - 4); offset <= seq2.length - 4; offset++) {
    let dH = 0;
    let dS = 0;
    let consecutiveMatches = 0;
    let maxConsecutive = 0;

    for (let i = 0; i < seq1.length; i++) {
      const j = i + offset;
      if (j >= 0 && j < seq2Rev.length) {
        // Check if seq1[i] pairs with seq2Rev[j]
        // seq2Rev[j] = seq2[len2-1-j]
        if (isComplement(seq1[i], seq2Rev[j])) {
          consecutiveMatches++;
          maxConsecutive = Math.max(maxConsecutive, consecutiveMatches);
          if (consecutiveMatches >= 2) {
            const prevI = i - 1;
            const prevJ = j - 1;
            if (prevI >= 0 && prevJ >= 0 && isComplement(seq1[prevI], seq2Rev[prevJ])) {
              // NN pair: strand1 5'->3' paired with strand2 3'->5'
              // seq2Rev contains the actual bases on the bottom strand
              // Format: "XY/WZ" where XY is top strand, WZ is bottom strand 3'->5'
              const pair = `${seq1[prevI]}${seq1[i]}/${seq2Rev[prevJ]}${seq2Rev[j]}`;
              const nn = DNA24_ENERGIES.NN[pair];
              if (nn) {
                dH += nn[0];
                dS += nn[1];
              }
            }
          }
        } else {
          consecutiveMatches = 0;
        }
      }
    }

    // Only count as heterodimer if we have at least 3 consecutive matches
    if (maxConsecutive >= 3 && dH < 0) {
      // Add initiation
      dH += 0.2;
      dS += -5.7;

      const dG = dH - tempK * (dS / 1000);
      if (dG < minDG) {
        minDG = dG;
      }
    }
  }

  // Cache result before returning
  heterodimerDGCache.set(cacheKey, minDG);
  return minDG;
}

/**
 * Calculate duplex ΔG between primer and its antiparallel target strand
 *
 * The target should be the strand that the primer hybridizes with, given in 5'->3' direction.
 * For a perfect match, target = reverseComplement(primer).
 *
 * Example:
 *   primer = 5'-ATGC-3'
 *   target = 5'-GCAT-3' (= reverseComplement(primer))
 *
 * When aligned antiparallel:
 *   Primer:  5'-A-T-G-C-3' (positions 0,1,2,3)
 *   Target:  3'-T-A-C-G-5' (target reversed = positions 3,2,1,0)
 *
 * NN pair format: "XY/WZ" where:
 * - XY is the primer dinucleotide 5'->3' (positions i, i+1)
 * - WZ is the target strand 3'->5' (target reversed, positions i, i+1)
 */
export function calculateDuplexDG(primer: string, target: string, temperature: number = 55): number {
  if (!primer || !target) return 0;

  primer = primer.toUpperCase();
  target = target.toUpperCase();

  if (primer.length !== target.length) {
    // Truncate to shorter length
    const len = Math.min(primer.length, target.length);
    primer = primer.slice(0, len);
    target = target.slice(0, len);
  }

  const tempK = temperature + 273.15;

  // Reverse target for antiparallel alignment
  // This gives us the target strand read 3'->5' (left to right in NN notation)
  const targetRev = target.split('').reverse().join('');

  let dH = 0.2;   // Initiation (DNA24)
  let dS = -5.7;

  // Add terminal A/T penalty based on the primer ends
  const firstBase = primer[0];
  const lastBase = primer[primer.length - 1];
  if (firstBase === 'A' || firstBase === 'T') {
    dH += 2.2;
    dS += 6.9;
  }
  if (lastBase === 'A' || lastBase === 'T') {
    dH += 2.2;
    dS += 6.9;
  }

  // Sum nearest-neighbor contributions
  // For aligned duplex:
  //   Primer:      5'-X-Y-3' (positions i, i+1)
  //   Target rev:  3'-W-Z-5' (positions i, i+1, 3'->5' direction)
  // NN pair format = "XY/WZ"
  for (let i = 0; i < primer.length - 1; i++) {
    const pair = `${primer[i]}${primer[i + 1]}/${targetRev[i]}${targetRev[i + 1]}`;
    const nn = DNA24_ENERGIES.NN[pair];
    if (nn) {
      dH += nn[0];
      dS += nn[1];
    }
  }

  return dH - tempK * (dS / 1000);
}

/**
 * Calculate the most stable off-target binding ΔG
 *
 * @param primer - Primer sequence
 * @param template - Full template sequence
 * @param excludeStart - Start position to exclude (intended binding site)
 * @param temperature - Temperature in Celsius
 * @returns Most stable off-target ΔG (most negative)
 */
export function calculateOffTargetDG(
  primer: string,
  template: string,
  excludeStart: number,
  temperature: number = 55
): number {
  if (!primer || !template || primer.length < 8) return 0;

  primer = primer.toUpperCase();
  template = template.toUpperCase();

  const minMatchLen = 8;  // Minimum 8bp match for potential off-target
  const excludeEnd = excludeStart + primer.length;
  let minDG = 0;

  // Get the 3' anchor of the primer (most important for mispriming)
  const anchor = primer.slice(-minMatchLen);

  // Check forward strand for 3' anchor matches
  for (let pos = 0; pos <= template.length - minMatchLen; pos++) {
    // Skip positions that overlap with the intended binding site
    if (pos >= excludeStart && pos < excludeEnd) continue;
    if (pos + minMatchLen > excludeStart && pos + minMatchLen <= excludeEnd) continue;

    const templateRegion = template.slice(pos, pos + minMatchLen);

    // Check if primer's 3' end matches this region
    // The primer binds antiparallel to the template, so we check if the anchor
    // (3' end of primer) equals reverseComplement of template region
    if (anchor === reverseComplement(templateRegion)) {
      // Calculate full binding ΔG for this site
      // The bindingSite IS the antiparallel partner strand (what primer hybridizes with)
      const bindingLen = Math.min(primer.length, template.length - pos);
      const bindingSite = template.slice(pos, pos + bindingLen);
      // For calculateDuplexDG: target = antiparallel partner = bindingSite (already correct orientation)
      const dG = calculateDuplexDG(primer.slice(-bindingLen), bindingSite, temperature);
      if (dG < minDG) {
        minDG = dG;
      }
    }
  }

  // Check reverse strand (primer binding directly to template as if it were antisense)
  for (let pos = 0; pos <= template.length - minMatchLen; pos++) {
    // Skip positions that overlap with the intended binding site
    if (pos >= excludeStart && pos < excludeEnd) continue;

    const templateRegion = template.slice(pos, pos + minMatchLen);

    // Direct match (primer binds to this strand directly)
    if (anchor === templateRegion) {
      const bindingLen = Math.min(primer.length, template.length - pos);
      const bindingSite = template.slice(pos, pos + bindingLen);
      // This would be primer binding to template directly (unusual but possible)
      const dG = calculateDuplexDG(primer.slice(-bindingLen), bindingSite, temperature);
      if (dG < minDG) {
        minDG = dG;
      }
    }
  }

  return minDG;
}

/**
 * Solve the equilibrium concentrations for all species
 *
 * This uses a simplified model based on partition functions.
 * For each primer, we calculate how it partitions between:
 * - Free (unbound)
 * - Hairpin (self-folded)
 * - Homodimer (self-dimer)
 * - Target-bound
 * - Off-target bound
 * - Heterodimer
 *
 * The partition is based on the relative Boltzmann weights (exp(-ΔG/RT))
 * of each state.
 */
export function solveEquilibrium(
  K: EquilibriumConstants,
  primerConc: number,
  templateConc: number
): Concentrations {
  // Compute effective concentrations for bimolecular reactions
  // For bimolecular: effective_K = K * [partner]
  // For unimolecular: effective_K = K

  // Fwd primer partition function
  // Q_fwd = 1 (free) + K_hairpin (unimolecular) + K_homodimer*[P] + K_target*[T] + K_offtarget*[T] + K_hetero*[Rev]
  const fwdHairpinK = K.fwdHairpin || 1e-10;
  const fwdHomodimerK = K.fwdHomodimer || 1e-10;
  const fwdTargetK = K.fwdTarget || 1e-10;
  const fwdOffTargetK = K.fwdOffTarget || 1e-10;
  const heteroK = K.heterodimer || 1e-10;

  const revHairpinK = K.revHairpin || 1e-10;
  const revHomodimerK = K.revHomodimer || 1e-10;
  const revTargetK = K.revTarget || 1e-10;
  const revOffTargetK = K.revOffTarget || 1e-10;

  // Compute partition functions
  // For simplicity, we assume template and partner primer concentrations are approximately constant
  // This is valid when template >> product and when heterodimer is not dominant

  const Q_fwd = 1 +
    fwdHairpinK +
    fwdHomodimerK * primerConc +
    fwdTargetK * templateConc +
    fwdOffTargetK * templateConc +
    heteroK * primerConc;

  const Q_rev = 1 +
    revHairpinK +
    revHomodimerK * primerConc +
    revTargetK * templateConc +
    revOffTargetK * templateConc +
    heteroK * primerConc;

  // Compute concentrations for each species
  // [X] = [P_total] * K_X / Q
  const fwdFree = primerConc * 1 / Q_fwd;
  const fwdHairpin = primerConc * fwdHairpinK / Q_fwd;
  const fwdHomodimer = primerConc * fwdHomodimerK * primerConc / Q_fwd / 2;  // Divide by 2 for stoichiometry
  const fwdTarget = primerConc * fwdTargetK * templateConc / Q_fwd;
  const fwdOffTarget = primerConc * fwdOffTargetK * templateConc / Q_fwd;
  const fwdHeterodimer = primerConc * heteroK * primerConc / Q_fwd / 2;

  const revFree = primerConc * 1 / Q_rev;
  const revHairpin = primerConc * revHairpinK / Q_rev;
  const revHomodimer = primerConc * revHomodimerK * primerConc / Q_rev / 2;
  const revTarget = primerConc * revTargetK * templateConc / Q_rev;
  const revOffTarget = primerConc * revOffTargetK * templateConc / Q_rev;
  const revHeterodimer = primerConc * heteroK * primerConc / Q_rev / 2;

  // Average heterodimer concentration (should be same for both primers in symmetric case)
  const heterodimer = (fwdHeterodimer + revHeterodimer) / 2;

  // Remaining template
  const templateFree = Math.max(0, templateConc - fwdTarget - revTarget);

  return {
    fwdFree,
    fwdHairpin,
    fwdHomodimer,
    fwdBoundTarget: fwdTarget,
    fwdOffTarget,
    revFree,
    revHairpin,
    revHomodimer,
    revBoundTarget: revTarget,
    revOffTarget,
    heterodimer,
    templateFree,
  };
}

/**
 * Classify efficiency into quality tiers
 */
function classifyEfficiency(efficiency: number): 'excellent' | 'good' | 'acceptable' | 'marginal' | 'poor' {
  if (efficiency >= 0.95) return 'excellent';
  if (efficiency >= 0.85) return 'good';
  if (efficiency >= 0.70) return 'acceptable';
  if (efficiency >= 0.50) return 'marginal';
  return 'poor';
}

/**
 * Convert efficiency to a 0-100 score using piecewise logistic
 *
 * @param efficiency - Equilibrium efficiency (0-1)
 * @param options - Scoring options
 * @returns Score from 0-100
 */
export function efficiencyToScore(efficiency: number, options: EfficiencyScoreOptions = {}): number {
  const {
    optimal = 0.95,
    acceptable = 0.70,
    steepness = 10,
  } = options;

  if (efficiency >= optimal) {
    return 100;
  }

  if (efficiency >= acceptable) {
    // Linear interpolation in acceptable range
    return 70 + 30 * (efficiency - acceptable) / (optimal - acceptable);
  }

  // Logistic decay below acceptable
  const excess = acceptable - efficiency;
  return 70 / (1 + Math.exp(steepness * excess));
}

// Helper functions

function isComplement(base1: string, base2: string): boolean {
  return DNA24_COMPLEMENT[base1] === base2;
}

function complement(base: string): string {
  return DNA24_COMPLEMENT[base] || base;
}

// Re-export reverseComplement from sequenceUtils for backward compatibility
export { reverseComplement, foldDG };
