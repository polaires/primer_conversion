/**
 * Off-Target Type A-F Classification
 *
 * Implements comprehensive off-target detection for plasmid PCR context.
 * Based on Section 2.4 of SCORING_MANUSCRIPT_STRATEGY.md
 *
 * Types:
 * - Type A: Alternative binding site (exact/near-exact match on template)
 * - Type B: Antisense binding (primer binds to reverse complement)
 * - Type C: Partial 3' homology (3' end matches elsewhere)
 * - Type D: Internal homology (non-3' homology, low risk)
 * - Type E: Self-complementarity (hairpin/homodimer) - integrated from equilibrium.js
 * - Type F: Primer-primer (heterodimer) - integrated from equilibrium.js
 *
 * Key insight from GM1: off-target is the "dominant failure factor"
 */

import {
  calculateDuplexDG,
  calculateHairpinDG,
  calculateHomodimerDG,
  calculateHeterodimerDG,
} from "./equilibrium.js";

// ============================================================================
// Types and Interfaces
// ============================================================================

export type RiskLevel = 'high' | 'medium' | 'low';

export interface OffTargetSite {
  type: string;
  subtype: string;
  position?: number;
  sequence?: string;
  mismatches?: number;
  dG?: number;
  risk: RiskLevel;
  description: string;
  anchorLength?: number;
  anchorSequence?: string;
  primerStart?: number;
  matchLength?: number;
  matchSequence?: string;
  threshold?: number;
}

export interface ThresholdConfig {
  exactMatchLength?: number;
  nearMatchLength?: number;
  anchorLength?: number;
  antisenseAnchorLength?: number;
  maxMismatches?: number;
  terminal3Mismatches?: number;
  stabilityThreshold?: number;
  temperature?: number;
  minInternalLength?: number;
  hairpinThreshold?: number;
  homodimerThreshold?: number;
  heterodimerThreshold?: number;
}

export interface OffTargetCounts {
  total: number;
  typeA: number;
  typeB: number;
  typeC: number;
  typeD: number;
  highRisk: number;
  mediumRisk: number;
  lowRisk: number;
}

export interface OffTargetClassification {
  status: 'pass' | 'warning' | 'critical';
  counts: OffTargetCounts;
  sites: {
    typeA: OffTargetSite[];
    typeB: OffTargetSite[];
    typeC: OffTargetSite[];
    typeD: OffTargetSite[];
  };
  allSites: OffTargetSite[];
  hasCritical: boolean;
  summary: string;
}

export interface TypeEAnalysis {
  hairpinDG: number;
  homodimerDG: number;
  sites: OffTargetSite[];
  hasIssues: boolean;
}

export interface TypeFAnalysis {
  heterodimerDG: number;
  sites: OffTargetSite[];
  hasIssues: boolean;
}

export interface PairOffTargetAnalysis {
  status: 'pass' | 'warning' | 'critical';
  fwd: OffTargetClassification & { typeE: TypeEAnalysis };
  rev: OffTargetClassification & { typeE: TypeEAnalysis };
  typeF: TypeFAnalysis;
  combined: {
    totalSites: number;
    highRisk: number;
    mediumRisk: number;
    lowRisk: number;
    typeE: number;
    typeF: number;
  };
  hasCritical: boolean;
  score: number;
  scoreBreakdown: {
    fwdBase: number;
    revBase: number;
    fwdEnhanced: number;
    revEnhanced: number;
  };
  thermodynamics: {
    fwdHairpinDG: number;
    fwdHomodimerDG: number;
    revHairpinDG: number;
    revHomodimerDG: number;
    heterodimerDG: number;
  };
}

// ============================================================================
// Constants
// ============================================================================

/**
 * Risk levels for off-target types
 */
export const RISK_LEVELS = {
  HIGH: 'high' as RiskLevel,
  MEDIUM: 'medium' as RiskLevel,
  LOW: 'low' as RiskLevel,
};

/**
 * Default thresholds for off-target detection
 */
export const DEFAULT_THRESHOLDS: Required<ThresholdConfig> = {
  // Minimum match lengths
  exactMatchLength: 15,        // Minimum length for exact match concern
  nearMatchLength: 12,         // Minimum length for near-exact match
  anchorLength: 8,             // Minimum 3' anchor for partial homology
  antisenseAnchorLength: 10,   // Minimum 3' anchor for antisense binding

  // Mismatch tolerance
  maxMismatches: 2,            // Maximum mismatches for near-exact
  terminal3Mismatches: 0,      // No mismatches in last 3bp for near-exact

  // ΔG threshold for viable priming
  stabilityThreshold: -11,     // kcal/mol - binding must be this stable to prime

  // Temperature for ΔG calculations
  temperature: 55,

  // Internal homology
  minInternalLength: 10,

  // Type E/F thresholds
  hairpinThreshold: -3,
  homodimerThreshold: -6,
  heterodimerThreshold: -6,
};

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Reverse complement a DNA sequence
 */
function reverseComplement(seq: string): string {
  const rc: Record<string, string> = { A: 'T', T: 'A', G: 'C', C: 'G' };
  return seq.split('').reverse().map(c => rc[c] || c).join('');
}

/**
 * Count mismatches between two sequences
 */
function countMismatches(seq1: string, seq2: string): number {
  if (seq1.length !== seq2.length) return Infinity;
  let count = 0;
  for (let i = 0; i < seq1.length; i++) {
    if (seq1[i] !== seq2[i]) count++;
  }
  return count;
}

/**
 * Check if the last n bases match (3' end alignment)
 */
function terminal3Match(seq1: string, seq2: string, n: number = 3): boolean {
  const end1 = seq1.slice(-n);
  const end2 = seq2.slice(-n);
  return end1 === end2;
}

/**
 * Find all positions where a substring occurs in a string
 */
function findAllOccurrences(text: string, pattern: string): number[] {
  const positions: number[] = [];
  let pos = 0;
  while ((pos = text.indexOf(pattern, pos)) !== -1) {
    positions.push(pos);
    pos++;
  }
  return positions;
}

// ============================================================================
// Type A-D Detection Functions
// ============================================================================

/**
 * Type A: Find exact and near-exact matches on template
 *
 * Detects alternative binding sites where the primer could bind
 * elsewhere on the plasmid, causing wrong products.
 *
 * @param primer - Primer sequence
 * @param template - Template sequence
 * @param intendedPosition - Position where primer is intended to bind
 * @param options - Detection options
 * @returns Array of Type A off-target sites
 */
export function findTypeAOffTargets(
  primer: string,
  template: string,
  intendedPosition: number = 0,
  options: ThresholdConfig = {}
): OffTargetSite[] {
  const {
    maxMismatches = DEFAULT_THRESHOLDS.maxMismatches,
    terminal3Mismatches = DEFAULT_THRESHOLDS.terminal3Mismatches,
    stabilityThreshold = DEFAULT_THRESHOLDS.stabilityThreshold,
    temperature = DEFAULT_THRESHOLDS.temperature,
  } = options;

  const sites: OffTargetSite[] = [];
  primer = primer.toUpperCase();
  template = template.toUpperCase();

  // The primer binds to the reverse complement of a template region
  // So we search for rc(primer's 3' anchor) in the template
  const anchor = primer.slice(-10);
  const rcAnchor = reverseComplement(anchor);
  const anchorPositions = findAllOccurrences(template, rcAnchor);

  for (const pos of anchorPositions) {
    // The primer's 3' end binds at the END of the rc region
    // So the full binding site starts at: pos - (primer.length - 10)
    // But since rcAnchor matches at the 3' end of binding, start is: pos
    // The full rc(primer) would be at pos - (primer.length - 10) to pos + 10
    const endPos = pos + 10;  // End of rc(anchor) in template
    const startPos = endPos - primer.length;

    // Skip the intended binding position (where rc(primer) is expected)
    if (Math.abs(startPos - intendedPosition) < 5) continue;

    if (startPos < 0 || endPos > template.length) continue;

    // Get the template region that would pair with the primer
    const templateRegion = template.slice(startPos, endPos);
    if (templateRegion.length < primer.length) continue;

    // The primer binds to rc(templateRegion), so compare primer to rc(templateRegion)
    const bindingPartner = reverseComplement(templateRegion);
    const mismatches = countMismatches(primer, bindingPartner);

    // Exact match
    if (mismatches === 0) {
      sites.push({
        type: 'A',
        subtype: 'exact',
        position: startPos,
        sequence: templateRegion,
        mismatches: 0,
        risk: RISK_LEVELS.HIGH,
        description: 'Exact match - alternative amplification site',
      });
    }
    // Near-exact match (≤2 mismatches, none in last 3bp)
    else if (mismatches <= maxMismatches) {
      // Check terminal 3bp of binding partner
      const terminal3Mm = countMismatches(primer.slice(-3), bindingPartner.slice(-3));
      if (terminal3Mm <= terminal3Mismatches) {
        // Calculate binding stability
        const dG = calculateDuplexDG(primer, bindingPartner, temperature);

        if (dG < stabilityThreshold) {
          sites.push({
            type: 'A',
            subtype: 'near_exact',
            position: startPos,
            sequence: templateRegion,
            mismatches,
            dG: Math.round(dG * 100) / 100,
            risk: RISK_LEVELS.HIGH,
            description: `Near-exact match (${mismatches} mismatches) - potential mispriming`,
          });
        }
      }
    }
  }

  return sites;
}

/**
 * Type B: Find antisense binding sites
 *
 * Detects where the primer could bind directly to the template strand
 * (same polarity as template), which would create primer-template hybrids
 * that prevent proper amplification.
 *
 * Normal binding: primer binds to rc(template_region)
 * Antisense binding: primer sequence appears directly in template (same strand)
 *
 * @param primer - Primer sequence
 * @param template - Template sequence
 * @param options - Detection options
 * @returns Array of Type B off-target sites
 */
export function findTypeBOffTargets(
  primer: string,
  template: string,
  options: ThresholdConfig = {}
): OffTargetSite[] {
  const {
    antisenseAnchorLength = DEFAULT_THRESHOLDS.antisenseAnchorLength,
    stabilityThreshold = DEFAULT_THRESHOLDS.stabilityThreshold,
    temperature = DEFAULT_THRESHOLDS.temperature,
  } = options;

  const sites: OffTargetSite[] = [];
  primer = primer.toUpperCase();
  template = template.toUpperCase();

  // For antisense binding, look for primer sequence appearing directly in template
  // This means the primer could bind to the complementary strand at this location
  // (which is effectively "antisense" to normal primer function)
  const anchor = primer.slice(-antisenseAnchorLength);

  // Find where primer's 3' anchor appears directly in template
  const positions = findAllOccurrences(template, anchor);

  for (const pos of positions) {
    // If the primer appears directly, it means there's a complementary site on opposite strand
    const endPos = pos + antisenseAnchorLength;
    const startPos = Math.max(0, endPos - primer.length);
    const templateRegion = template.slice(startPos, endPos);

    // For antisense, the primer would bind to rc(templateRegion) on the opposite strand
    // This creates a conflict with normal template replication
    const bindingPartner = reverseComplement(templateRegion);
    const dG = calculateDuplexDG(primer.slice(-(endPos - startPos)), bindingPartner, temperature);

    if (dG < stabilityThreshold) {
      sites.push({
        type: 'B',
        subtype: 'antisense',
        position: startPos,
        sequence: templateRegion,
        dG: Math.round(dG * 100) / 100,
        risk: RISK_LEVELS.HIGH,
        description: 'Antisense binding - primer appears in template strand',
      });
    }
  }

  return sites;
}

/**
 * Type C: Find partial 3' homology sites
 *
 * Detects where just the 3' end of the primer matches elsewhere,
 * which can cause mispriming even without full-length match.
 *
 * @param primer - Primer sequence
 * @param template - Template sequence
 * @param intendedPosition - Position where primer is intended to bind
 * @param options - Detection options
 * @returns Array of Type C off-target sites
 */
export function findTypeCOffTargets(
  primer: string,
  template: string,
  intendedPosition: number = 0,
  options: ThresholdConfig = {}
): OffTargetSite[] {
  const {
    anchorLength = DEFAULT_THRESHOLDS.anchorLength,
    stabilityThreshold = DEFAULT_THRESHOLDS.stabilityThreshold,
    temperature = DEFAULT_THRESHOLDS.temperature,
  } = options;

  const sites: OffTargetSite[] = [];
  primer = primer.toUpperCase();
  template = template.toUpperCase();

  // Check multiple anchor lengths (8, 10, 12bp)
  for (const len of [8, 10, 12]) {
    if (len > primer.length) continue;

    const anchor = primer.slice(-len);
    const rcAnchor = reverseComplement(anchor);

    // Find where the anchor's reverse complement appears in template
    // (this is where the 3' end could bind)
    const positions = findAllOccurrences(template, rcAnchor);

    for (const pos of positions) {
      // Skip if this is part of the intended binding region
      if (Math.abs(pos - intendedPosition) < primer.length) continue;

      // Calculate binding stability for just the 3' anchor
      const dG = calculateDuplexDG(anchor, rcAnchor, temperature);

      // Only flag if binding is stable enough to potentially prime
      if (dG < stabilityThreshold + 3) {  // Slightly relaxed threshold for partial
        // Check if we already have this site from a longer anchor
        const existingSite = sites.find(s => Math.abs(s.position! - pos) < 5);
        if (existingSite && existingSite.anchorLength! >= len) continue;

        if (existingSite) {
          // Update with longer anchor
          existingSite.anchorLength = len;
          existingSite.dG = Math.round(dG * 100) / 100;
        } else {
          sites.push({
            type: 'C',
            subtype: 'partial_3prime',
            position: pos,
            anchorLength: len,
            anchorSequence: anchor,
            dG: Math.round(dG * 100) / 100,
            risk: RISK_LEVELS.MEDIUM,
            description: `Partial 3' homology (${len}bp) - potential mispriming`,
          });
        }
      }
    }
  }

  return sites;
}

/**
 * Type D: Find internal homology sites
 *
 * Detects homology that doesn't involve the 3' end.
 * Lower risk since extension can't initiate from non-3' regions.
 *
 * @param primer - Primer sequence
 * @param template - Template sequence
 * @param intendedPosition - Position where primer is intended to bind
 * @param options - Detection options
 * @returns Array of Type D off-target sites
 */
export function findTypeDOffTargets(
  primer: string,
  template: string,
  intendedPosition: number = 0,
  options: ThresholdConfig = {}
): OffTargetSite[] {
  const {
    minInternalLength = 10,
    stabilityThreshold = DEFAULT_THRESHOLDS.stabilityThreshold,
    temperature = DEFAULT_THRESHOLDS.temperature,
  } = options;

  const sites: OffTargetSite[] = [];
  primer = primer.toUpperCase();
  template = template.toUpperCase();

  // Check for internal regions (not including 3' end)
  const internalPrimer = primer.slice(0, -5);  // Exclude last 5bp

  for (let len = Math.min(12, internalPrimer.length); len >= minInternalLength; len--) {
    for (let start = 0; start <= internalPrimer.length - len; start++) {
      const region = internalPrimer.slice(start, start + len);
      const rcRegion = reverseComplement(region);

      const positions = findAllOccurrences(template, rcRegion);

      for (const pos of positions) {
        // Skip intended binding region
        if (Math.abs(pos - intendedPosition) < primer.length) continue;

        // Check if already covered by a longer match
        const existingSite = sites.find(s =>
          pos >= s.position! && pos < s.position! + s.matchLength!
        );
        if (existingSite) continue;

        sites.push({
          type: 'D',
          subtype: 'internal',
          position: pos,
          primerStart: start,
          matchLength: len,
          matchSequence: region,
          risk: RISK_LEVELS.LOW,
          description: 'Internal homology - usually tolerable',
        });
      }
    }
  }

  return sites;
}

// ============================================================================
// Comprehensive Classification
// ============================================================================

/**
 * Comprehensive off-target analysis
 *
 * Runs all Type A-D checks and returns classified results.
 * Types E (self-complementarity) and F (primer-primer) are handled
 * by the equilibrium module's hairpin/dimer calculations.
 *
 * @param primer - Primer sequence
 * @param template - Template sequence
 * @param intendedPosition - Position where primer is intended to bind
 * @param options - Detection options
 * @returns Classified off-target analysis results
 */
export function classifyOffTargets(
  primer: string,
  template: string,
  intendedPosition: number = 0,
  options: ThresholdConfig = {}
): OffTargetClassification {
  const typeA = findTypeAOffTargets(primer, template, intendedPosition, options);
  const typeB = findTypeBOffTargets(primer, template, options);
  const typeC = findTypeCOffTargets(primer, template, intendedPosition, options);
  const typeD = findTypeDOffTargets(primer, template, intendedPosition, options);

  const allSites = [...typeA, ...typeB, ...typeC, ...typeD];

  // Count by risk level
  const highRisk = allSites.filter(s => s.risk === RISK_LEVELS.HIGH);
  const mediumRisk = allSites.filter(s => s.risk === RISK_LEVELS.MEDIUM);
  const lowRisk = allSites.filter(s => s.risk === RISK_LEVELS.LOW);

  // Determine overall status
  let status: 'pass' | 'warning' | 'critical' = 'pass';
  if (highRisk.length > 0) {
    status = 'critical';  // High risk of failure
  } else if (mediumRisk.length > 0) {
    status = 'warning';   // May have issues
  }

  return {
    status,
    counts: {
      total: allSites.length,
      typeA: typeA.length,
      typeB: typeB.length,
      typeC: typeC.length,
      typeD: typeD.length,
      highRisk: highRisk.length,
      mediumRisk: mediumRisk.length,
      lowRisk: lowRisk.length,
    },
    sites: {
      typeA,
      typeB,
      typeC,
      typeD,
    },
    allSites,
    hasCritical: highRisk.length > 0,
    summary: generateSummary(typeA, typeB, typeC, typeD),
  };
}

/**
 * Generate human-readable summary of off-target findings
 */
function generateSummary(
  typeA: OffTargetSite[],
  typeB: OffTargetSite[],
  typeC: OffTargetSite[],
  typeD: OffTargetSite[]
): string {
  const parts: string[] = [];

  if (typeA.length > 0) {
    const exact = typeA.filter(s => s.subtype === 'exact').length;
    const near = typeA.filter(s => s.subtype === 'near_exact').length;
    if (exact > 0) parts.push(`${exact} exact match${exact > 1 ? 'es' : ''} (HIGH RISK)`);
    if (near > 0) parts.push(`${near} near-exact match${near > 1 ? 'es' : ''} (HIGH RISK)`);
  }

  if (typeB.length > 0) {
    parts.push(`${typeB.length} antisense binding site${typeB.length > 1 ? 's' : ''} (HIGH RISK)`);
  }

  if (typeC.length > 0) {
    parts.push(`${typeC.length} partial 3' homology site${typeC.length > 1 ? 's' : ''} (MEDIUM RISK)`);
  }

  if (typeD.length > 0) {
    parts.push(`${typeD.length} internal homology site${typeD.length > 1 ? 's' : ''} (LOW RISK)`);
  }

  if (parts.length === 0) {
    return 'No significant off-target sites detected';
  }

  return parts.join('; ');
}

// ============================================================================
// Scoring Functions
// ============================================================================

/**
 * Calculate off-target score (0-1, higher is better)
 *
 * Based on GM1 finding that off-target is dominant failure factor.
 * Uses exponential penalty scaling.
 *
 * @param classification - Result from classifyOffTargets
 * @returns Score from 0 to 1
 */
export function scoreOffTargetClassification(classification: OffTargetClassification): number {
  const { counts } = classification;

  // High risk sites are critical - each one severely impacts score
  if (counts.highRisk >= 3) return 0.0;  // Disqualify
  if (counts.highRisk === 2) return 0.1;
  if (counts.highRisk === 1) return 0.3;

  // Medium risk sites - moderate impact
  let score = 1.0;
  score -= counts.mediumRisk * 0.15;

  // Low risk sites - minor impact
  score -= counts.lowRisk * 0.02;

  return Math.max(0, Math.min(1, score));
}

// ============================================================================
// Type E and F Detection (Self-complementarity and Heterodimers)
// ============================================================================

/**
 * Type E: Detect self-complementarity issues (hairpin and homodimer)
 *
 * Based on manuscript: "Type E: Self-complementarity (hairpin) - MEDIUM RISK"
 * Only penalize if ΔG is stable enough to compete with target binding.
 *
 * @param primer - Primer sequence
 * @param options - Detection options
 * @returns Type E analysis results
 */
export function findTypeEOffTargets(primer: string, options: ThresholdConfig = {}): TypeEAnalysis {
  const {
    hairpinThreshold = -3,      // ΔG threshold for significant hairpin
    homodimerThreshold = -6,    // ΔG threshold for significant homodimer
    temperature = DEFAULT_THRESHOLDS.temperature,
  } = options;

  primer = primer.toUpperCase();
  const sites: OffTargetSite[] = [];

  // Calculate hairpin ΔG
  const hairpinDG = calculateHairpinDG(primer, temperature);
  if (hairpinDG < hairpinThreshold) {
    sites.push({
      type: 'E',
      subtype: 'hairpin',
      dG: Math.round(hairpinDG * 100) / 100,
      threshold: hairpinThreshold,
      risk: hairpinDG < hairpinThreshold - 3 ? RISK_LEVELS.HIGH : RISK_LEVELS.MEDIUM,
      description: `Stable hairpin (ΔG=${hairpinDG.toFixed(1)} kcal/mol) - competes with target binding`,
    });
  }

  // Calculate homodimer ΔG
  const homodimerDG = calculateHomodimerDG(primer, temperature);
  if (homodimerDG < homodimerThreshold) {
    sites.push({
      type: 'E',
      subtype: 'homodimer',
      dG: Math.round(homodimerDG * 100) / 100,
      threshold: homodimerThreshold,
      risk: homodimerDG < homodimerThreshold - 4 ? RISK_LEVELS.HIGH : RISK_LEVELS.MEDIUM,
      description: `Stable homodimer (ΔG=${homodimerDG.toFixed(1)} kcal/mol) - primer self-dimer`,
    });
  }

  return {
    hairpinDG: Math.round(hairpinDG * 100) / 100,
    homodimerDG: Math.round(homodimerDG * 100) / 100,
    sites,
    hasIssues: sites.length > 0,
  };
}

/**
 * Type F: Detect primer-primer heterodimer issues
 *
 * Based on manuscript: "Type F: Primer-primer (heterodimer) - MEDIUM RISK"
 * Cross-dimer between forward and reverse primers.
 *
 * @param fwd - Forward primer sequence
 * @param rev - Reverse primer sequence
 * @param options - Detection options
 * @returns Type F analysis results
 */
export function findTypeFOffTargets(fwd: string, rev: string, options: ThresholdConfig = {}): TypeFAnalysis {
  const {
    heterodimerThreshold = -6,  // ΔG threshold for significant heterodimer
    temperature = DEFAULT_THRESHOLDS.temperature,
  } = options;

  fwd = fwd.toUpperCase();
  rev = rev.toUpperCase();
  const sites: OffTargetSite[] = [];

  // Calculate heterodimer ΔG
  const heterodimerDG = calculateHeterodimerDG(fwd, rev, temperature);
  if (heterodimerDG < heterodimerThreshold) {
    const risk = heterodimerDG < heterodimerThreshold - 4 ? RISK_LEVELS.HIGH : RISK_LEVELS.MEDIUM;
    sites.push({
      type: 'F',
      subtype: 'heterodimer',
      dG: Math.round(heterodimerDG * 100) / 100,
      threshold: heterodimerThreshold,
      risk,
      description: `Stable heterodimer (ΔG=${heterodimerDG.toFixed(1)} kcal/mol) - fwd-rev dimer`,
    });
  }

  return {
    heterodimerDG: Math.round(heterodimerDG * 100) / 100,
    sites,
    hasIssues: sites.length > 0,
  };
}

/**
 * Enhanced off-target scoring including Type E and F
 *
 * Integrates hairpin, homodimer, and heterodimer penalties into the score.
 * Only penalizes structures stable enough to compete with target binding
 * (ΔG < threshold).
 *
 * @param classification - Result from classifyOffTargets
 * @param typeE - Type E analysis (hairpin/homodimer)
 * @param typeF - Type F analysis (heterodimer)
 * @returns Score from 0 to 1
 */
export function scoreOffTargetClassificationEnhanced(
  classification: OffTargetClassification,
  typeE: TypeEAnalysis | null = null,
  typeF: TypeFAnalysis | null = null
): number {
  let score = scoreOffTargetClassification(classification);

  // Type E: Self-complementarity penalties
  if (typeE && typeE.sites) {
    for (const site of typeE.sites) {
      if (site.risk === RISK_LEVELS.HIGH) {
        score -= 0.25;  // Severe penalty for very stable structures
      } else if (site.risk === RISK_LEVELS.MEDIUM) {
        score -= 0.10;  // Moderate penalty
      }
    }
  }

  // Type F: Heterodimer penalties
  if (typeF && typeF.sites) {
    for (const site of typeF.sites) {
      if (site.risk === RISK_LEVELS.HIGH) {
        score -= 0.30;  // Severe penalty for very stable heterodimer
      } else if (site.risk === RISK_LEVELS.MEDIUM) {
        score -= 0.15;  // Moderate penalty
      }
    }
  }

  return Math.max(0, Math.min(1, score));
}

/**
 * Quick off-target check for a primer pair
 *
 * Now includes Type E (hairpin/homodimer) and Type F (heterodimer) analysis.
 *
 * @param fwd - Forward primer sequence
 * @param rev - Reverse primer sequence
 * @param template - Template sequence
 * @param options - Detection options
 * @returns Combined off-target analysis for the pair
 */
export function analyzeOffTargetsPair(
  fwd: string,
  rev: string,
  template: string,
  options: ThresholdConfig = {}
): PairOffTargetAnalysis {
  // Find intended binding positions
  const fwdPos = template.indexOf(fwd.slice(-10));
  const revPos = template.indexOf(reverseComplement(rev.slice(-10)));

  // Type A-D: Template-based off-targets
  const fwdAnalysis = classifyOffTargets(fwd, template, fwdPos, options);
  const revAnalysis = classifyOffTargets(rev, template, revPos, options);

  // Type E: Self-complementarity (hairpin/homodimer) for each primer
  const fwdTypeE = findTypeEOffTargets(fwd, options);
  const revTypeE = findTypeEOffTargets(rev, options);

  // Type F: Heterodimer between forward and reverse
  const typeF = findTypeFOffTargets(fwd, rev, options);

  // Combined status
  let status: 'pass' | 'warning' | 'critical' = 'pass';
  if (fwdAnalysis.status === 'critical' || revAnalysis.status === 'critical') {
    status = 'critical';
  } else if (fwdAnalysis.status === 'warning' || revAnalysis.status === 'warning') {
    status = 'warning';
  }
  // Escalate status for severe Type E/F issues
  if (typeF.sites.some(s => s.risk === RISK_LEVELS.HIGH) ||
      fwdTypeE.sites.some(s => s.risk === RISK_LEVELS.HIGH) ||
      revTypeE.sites.some(s => s.risk === RISK_LEVELS.HIGH)) {
    status = status === 'pass' ? 'warning' : status;
  }

  // Calculate enhanced score including Type E/F
  const fwdBaseScore = scoreOffTargetClassification(fwdAnalysis);
  const revBaseScore = scoreOffTargetClassification(revAnalysis);
  const fwdEnhancedScore = scoreOffTargetClassificationEnhanced(fwdAnalysis, fwdTypeE, null);
  const revEnhancedScore = scoreOffTargetClassificationEnhanced(revAnalysis, revTypeE, typeF);

  // Count all issues
  const allTypeESites = [...fwdTypeE.sites, ...revTypeE.sites];
  const typeEHighRisk = allTypeESites.filter(s => s.risk === RISK_LEVELS.HIGH).length;
  const typeEMediumRisk = allTypeESites.filter(s => s.risk === RISK_LEVELS.MEDIUM).length;
  const typeFHighRisk = typeF.sites.filter(s => s.risk === RISK_LEVELS.HIGH).length;
  const typeFMediumRisk = typeF.sites.filter(s => s.risk === RISK_LEVELS.MEDIUM).length;

  return {
    status,
    fwd: {
      ...fwdAnalysis,
      typeE: fwdTypeE,
    },
    rev: {
      ...revAnalysis,
      typeE: revTypeE,
    },
    typeF,
    combined: {
      totalSites: fwdAnalysis.counts.total + revAnalysis.counts.total + allTypeESites.length + typeF.sites.length,
      highRisk: fwdAnalysis.counts.highRisk + revAnalysis.counts.highRisk + typeEHighRisk + typeFHighRisk,
      mediumRisk: fwdAnalysis.counts.mediumRisk + revAnalysis.counts.mediumRisk + typeEMediumRisk + typeFMediumRisk,
      lowRisk: fwdAnalysis.counts.lowRisk + revAnalysis.counts.lowRisk,
      typeE: allTypeESites.length,
      typeF: typeF.sites.length,
    },
    hasCritical: fwdAnalysis.hasCritical || revAnalysis.hasCritical ||
                 typeEHighRisk > 0 || typeFHighRisk > 0,
    // Use enhanced score that includes Type E/F
    score: Math.min(fwdEnhancedScore, revEnhancedScore),
    // Also provide breakdown
    scoreBreakdown: {
      fwdBase: Math.round(fwdBaseScore * 1000) / 1000,
      revBase: Math.round(revBaseScore * 1000) / 1000,
      fwdEnhanced: Math.round(fwdEnhancedScore * 1000) / 1000,
      revEnhanced: Math.round(revEnhancedScore * 1000) / 1000,
    },
    // Thermodynamic values for reference
    thermodynamics: {
      fwdHairpinDG: fwdTypeE.hairpinDG,
      fwdHomodimerDG: fwdTypeE.homodimerDG,
      revHairpinDG: revTypeE.hairpinDG,
      revHomodimerDG: revTypeE.homodimerDG,
      heterodimerDG: typeF.heterodimerDG,
    },
  };
}
