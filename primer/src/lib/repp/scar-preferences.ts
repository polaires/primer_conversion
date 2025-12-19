/**
 * Scar Sequence Preferences for Golden Gate Assembly
 *
 * The 4bp overhang becomes a "scar" in the final assembled sequence.
 * This module provides scoring and recommendations based on biological context.
 *
 * Key considerations:
 * 1. For coding sequences - avoid stop codons, prefer flexible amino acids
 * 2. For linker regions - prefer Gly/Ser-rich sequences
 * 3. For non-coding - focus on assembly efficiency
 */

// ============================================================================
// TYPES
// ============================================================================

interface AminoAcidProperties {
  flexible?: boolean;
  small?: boolean;
  hydrophobic?: boolean;
  charged?: boolean;
  name: string;
  kinkInducing?: boolean;
  reactive?: boolean;
  stop?: boolean;
  charge?: number;
}

interface ScarPreference {
  preferred?: string[];
  avoid?: string[];
  preferredBonus?: number;
  avoidPenalty?: number;
}

interface StopCodonInfo {
  frame: number;
  codon: string;
  position: number;
  note?: string;
}

interface StopCodonResult {
  overhang: string;
  hasStopCodon: boolean;
  stopCodons: StopCodonInfo[];
  isSafe: boolean;
  checkedFrame2: boolean;
}

interface TranslationFrame {
  codons: string[];
  aminoAcids: string[];
  description: string;
  note?: string;
}

interface TranslationResult {
  overhang: string;
  translations: Record<string, TranslationFrame>;
  hasStopInAnyFrame: boolean;
}

interface ScarNote {
  type: 'warning' | 'error' | 'info';
  message: string;
}

interface ScarScoreResult {
  overhang: string;
  context: string;
  score: number;
  scorePercent: string;
  isPreferred: boolean;
  isAvoided: boolean;
  notes: ScarNote[];
  quality: 'excellent' | 'good' | 'acceptable' | 'poor';
}

interface ScarSetResult {
  context: string;
  overhangs: string[];
  individual: ScarScoreResult[];
  averageScore: number;
  worstScore: number;
  problematicCount: number;
  problematic: string[];
  allPreferred: boolean;
  hasAvoided: boolean;
  recommendation: string;
}

interface ScarCandidate extends ScarScoreResult {
  overhang: string;
}

interface BestScarResult {
  context: string;
  candidates: ScarCandidate[];
  best: ScarCandidate | null;
  preferred: ScarCandidate[];
  acceptable: ScarCandidate[];
  poor: ScarCandidate[];
}

interface AssemblyInfo {
  isCodingSequence?: boolean;
  hasStartCodon?: boolean;
  hasStopCodon?: boolean;
  isLinkerRegion?: boolean;
  junctionPositions?: number[];
}

interface JunctionRecommendation {
  junctionIndex: number;
  position: number;
  recommendedContext: string;
  preferredOverhangs: string[];
  avoidOverhangs: string[];
}

interface ContextRecommendation {
  assemblyInfo: AssemblyInfo;
  recommendations: JunctionRecommendation[];
  summary: string;
}

type CodonContext = 'coding' | 'linker' | 'nonCoding' | 'startCodon' | 'stopCodon';

// ============================================================================
// CODON TABLES
// ============================================================================

/**
 * Codon translation table (DNA -> amino acid)
 */
const CODON_TO_AA: Record<string, string> = {
  // Phenylalanine
  'TTT': 'F', 'TTC': 'F',
  // Leucine
  'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
  // Isoleucine
  'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
  // Methionine (start)
  'ATG': 'M',
  // Valine
  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
  // Serine
  'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
  // Proline
  'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
  // Threonine
  'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
  // Alanine
  'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
  // Tyrosine
  'TAT': 'Y', 'TAC': 'Y',
  // Stop codons
  'TAA': '*', 'TAG': '*', 'TGA': '*',
  // Histidine
  'CAT': 'H', 'CAC': 'H',
  // Glutamine
  'CAA': 'Q', 'CAG': 'Q',
  // Asparagine
  'AAT': 'N', 'AAC': 'N',
  // Lysine
  'AAA': 'K', 'AAG': 'K',
  // Aspartic acid
  'GAT': 'D', 'GAC': 'D',
  // Glutamic acid
  'GAA': 'E', 'GAG': 'E',
  // Cysteine
  'TGT': 'C', 'TGC': 'C',
  // Tryptophan
  'TGG': 'W',
  // Arginine
  'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
  // Glycine
  'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
};

/**
 * Amino acid properties for linker design
 */
const AA_PROPERTIES: Record<string, AminoAcidProperties> = {
  'G': { flexible: true, small: true, hydrophobic: false, charged: false, name: 'Glycine' },
  'A': { flexible: true, small: true, hydrophobic: true, charged: false, name: 'Alanine' },
  'S': { flexible: true, small: true, hydrophobic: false, charged: false, name: 'Serine' },
  'T': { flexible: false, small: false, hydrophobic: false, charged: false, name: 'Threonine' },
  'P': { flexible: false, small: false, hydrophobic: false, charged: false, name: 'Proline', kinkInducing: true },
  'V': { flexible: false, small: false, hydrophobic: true, charged: false, name: 'Valine' },
  'L': { flexible: false, small: false, hydrophobic: true, charged: false, name: 'Leucine' },
  'I': { flexible: false, small: false, hydrophobic: true, charged: false, name: 'Isoleucine' },
  'M': { flexible: false, small: false, hydrophobic: true, charged: false, name: 'Methionine' },
  'F': { flexible: false, small: false, hydrophobic: true, charged: false, name: 'Phenylalanine' },
  'Y': { flexible: false, small: false, hydrophobic: true, charged: false, name: 'Tyrosine' },
  'W': { flexible: false, small: false, hydrophobic: true, charged: false, name: 'Tryptophan' },
  'C': { flexible: false, small: true, hydrophobic: true, charged: false, name: 'Cysteine', reactive: true },
  'N': { flexible: true, small: false, hydrophobic: false, charged: false, name: 'Asparagine' },
  'Q': { flexible: true, small: false, hydrophobic: false, charged: false, name: 'Glutamine' },
  'D': { flexible: false, small: false, hydrophobic: false, charged: true, name: 'Aspartic acid', charge: -1 },
  'E': { flexible: false, small: false, hydrophobic: false, charged: true, name: 'Glutamic acid', charge: -1 },
  'K': { flexible: true, small: false, hydrophobic: false, charged: true, name: 'Lysine', charge: +1 },
  'R': { flexible: false, small: false, hydrophobic: false, charged: true, name: 'Arginine', charge: +1 },
  'H': { flexible: false, small: false, hydrophobic: false, charged: true, name: 'Histidine', charge: +1 },
  '*': { stop: true, name: 'Stop' },
};

/**
 * Scar sequence preferences for different contexts
 */
export const SCAR_PREFERENCES: Record<CodonContext, ScarPreference> = {
  // For coding sequences - prefer overhangs encoding useful amino acids
  coding: {
    preferred: [
      // Overhangs that encode flexible/common amino acids
      'GGAG',  // Gly-X (frame 0)
      'GGTG',  // Gly-X (frame 0)
      'GCAG',  // Ala-X (frame 0)
      'AATG',  // Asn-Met (good for start context) - MoClo standard
      'GCTT',  // Ala-X - MoClo standard
      'TACT',  // Tyr-X - MoClo standard
      'AGGT',  // Frame-dependent
      'GCTG',  // Ala-X
      'TCTG',  // Ser-X
    ],
    avoid: [
      // Stop codons in any frame
      'TGAT',  // Contains TGA (stop) in frame 0
      'TAAT',  // Contains TAA (stop) in frame 0
      'TAGT',  // Contains TAG (stop) in frame 0
      'ATGA',  // Contains TGA in frame 1
      'ATAA',  // Contains TAA in frame 1
      'CTAG',  // Contains TAG in frame 1
      'GTAA',  // Contains TAA in frame 2
      'GTAG',  // Contains TAG in frame 2
      'GTGA',  // Contains TGA in frame 2
    ],
    // Penalty for using avoided overhangs
    avoidPenalty: 50,
    preferredBonus: 10,
  },

  // For linker regions between domains
  linker: {
    preferred: [
      'GGAG',  // Gly-Gly linker context
      'GGTG',  // Gly-X
      'GGCG',  // Gly-X
      'GCAG',  // Ala-X
      'TCTG',  // Ser-X
      'TCAG',  // Ser-X
    ],
    preferredBonus: 15,
  },

  // For non-coding regions (promoters, terminators)
  nonCoding: {
    // No specific preferences - just avoid creating problems
    avoid: [],
    preferredBonus: 0,
    avoidPenalty: 0,
  },

  // For start codon context (ATG should be in overhang)
  startCodon: {
    preferred: [
      'AATG',  // Standard MoClo - contains ATG
      'CATG',  // Contains ATG
      'GATG',  // Contains ATG
      'TATG',  // Contains ATG
    ],
    preferredBonus: 20,
  },

  // For stop codon context
  stopCodon: {
    preferred: [
      'AGGT',  // MoClo standard for CDS end
      'GCTT',  // After stop
    ],
    preferredBonus: 10,
  },
};

// ============================================================================
// MAIN FUNCTIONS
// ============================================================================

/**
 * Check if an overhang contains a stop codon in any reading frame
 */
export function checkStopCodons(overhang: string, contextBefore: string = ''): StopCodonResult {
  const oh = overhang.toUpperCase();
  const before = contextBefore.toUpperCase();
  const stopCodons = ['TAA', 'TAG', 'TGA'];
  const found: StopCodonInfo[] = [];

  // Frame 0: positions 0-2 of overhang
  const codon0 = oh.slice(0, 3);
  if (stopCodons.includes(codon0)) {
    found.push({ frame: 0, codon: codon0, position: 0 });
  }

  // Frame 1: positions 1-3 of overhang
  const codon1 = oh.slice(1, 4);
  if (stopCodons.includes(codon1)) {
    found.push({ frame: 1, codon: codon1, position: 1 });
  }

  // Frame 2: requires context - last 2 bases of upstream + first base of overhang
  if (before.length >= 2) {
    const codon2 = before.slice(-2) + oh.slice(0, 1);
    if (stopCodons.includes(codon2)) {
      found.push({ frame: 2, codon: codon2, position: -2, note: 'Requires upstream context' });
    }
  }

  return {
    overhang: oh,
    hasStopCodon: found.length > 0,
    stopCodons: found,
    isSafe: found.length === 0,
    checkedFrame2: before.length >= 2,
  };
}

/**
 * Translate an overhang to amino acids in all frames
 */
export function translateOverhang(
  overhang: string,
  contextBefore: string = '',
  contextAfter: string = ''
): TranslationResult {
  const oh = overhang.toUpperCase();
  const before = contextBefore.toUpperCase();
  const after = contextAfter.toUpperCase();

  const fullContext = before + oh + after;
  const translations: Record<string, TranslationFrame> = {};

  // Frame 0: starts at position 0 of overhang
  if (oh.length >= 3) {
    const codon1 = oh.slice(0, 3);
    const aa1 = CODON_TO_AA[codon1] || '?';
    translations.frame0 = {
      codons: [codon1],
      aminoAcids: [aa1],
      description: `${codon1} -> ${aa1}`,
    };
  }

  // Frame 1: starts at position 1 of overhang
  if (oh.length >= 4) {
    const codon = oh.slice(1, 4);
    const aa = CODON_TO_AA[codon] || '?';
    translations.frame1 = {
      codons: [codon],
      aminoAcids: [aa],
      description: `${codon} -> ${aa}`,
    };
  }

  // Frame 2: requires context (last 2 of before + first 1 of oh)
  if (before.length >= 2 && oh.length >= 1) {
    const codon = before.slice(-2) + oh.slice(0, 1);
    if (codon.length === 3) {
      const aa = CODON_TO_AA[codon] || '?';
      translations.frame2 = {
        codons: [codon],
        aminoAcids: [aa],
        description: `${codon} -> ${aa}`,
        note: 'Requires upstream context',
      };
    }
  }

  return {
    overhang: oh,
    translations,
    hasStopInAnyFrame: Object.values(translations).some(t =>
      t.aminoAcids.includes('*')
    ),
  };
}

/**
 * Score an overhang based on scar preferences
 */
export function scoreScarSequence(overhang: string, context: CodonContext = 'nonCoding'): ScarScoreResult {
  const prefs = SCAR_PREFERENCES[context] || SCAR_PREFERENCES.nonCoding;
  const oh = overhang.toUpperCase();

  let score = 80; // Base score
  const notes: ScarNote[] = [];

  // Check if in avoid list
  if (prefs.avoid?.includes(oh)) {
    score -= (prefs.avoidPenalty || 50);
    notes.push({
      type: 'warning',
      message: `Overhang ${oh} may create problematic scar sequence in ${context} context`,
    });
  }

  // Check if in preferred list
  if (prefs.preferred?.includes(oh)) {
    score += (prefs.preferredBonus || 10);
    notes.push({
      type: 'info',
      message: `Overhang ${oh} creates favorable scar sequence for ${context}`,
    });
  }

  // Additional checks for coding context
  if (context === 'coding' || context === 'linker') {
    const stopCheck = checkStopCodons(oh);
    if (stopCheck.hasStopCodon) {
      score -= 40;
      notes.push({
        type: 'error',
        message: `Contains stop codon(s) in frame ${stopCheck.stopCodons.map(s => s.frame).join(', ')}`,
      });
    }
  }

  // Normalize score to 0-100
  score = Math.max(0, Math.min(100, score));

  return {
    overhang: oh,
    context,
    score,
    scorePercent: `${score}%`,
    isPreferred: prefs.preferred?.includes(oh) || false,
    isAvoided: prefs.avoid?.includes(oh) || false,
    notes,
    quality: score >= 90 ? 'excellent' :
             score >= 70 ? 'good' :
             score >= 50 ? 'acceptable' : 'poor',
  };
}

/**
 * Score multiple overhangs for a specific context
 */
export function scoreScarSet(overhangs: string[], context: CodonContext = 'nonCoding'): ScarSetResult {
  const individual = overhangs.map(oh => scoreScarSequence(oh, context));

  const avgScore = individual.reduce((s, i) => s + i.score, 0) / individual.length;
  const worstScore = Math.min(...individual.map(i => i.score));
  const problematic = individual.filter(i => i.score < 70);

  return {
    context,
    overhangs: overhangs.map(oh => oh.toUpperCase()),
    individual,
    averageScore: Math.round(avgScore),
    worstScore,
    problematicCount: problematic.length,
    problematic: problematic.map(p => p.overhang),
    allPreferred: individual.every(i => i.isPreferred),
    hasAvoided: individual.some(i => i.isAvoided),
    recommendation: problematic.length > 0
      ? `${problematic.length} overhang(s) have poor scar scores - consider alternatives`
      : 'Scar sequences are acceptable for this context',
  };
}

/**
 * Find best overhang for a specific context from candidates
 */
export function findBestScarCandidate(candidates: string[], context: CodonContext = 'nonCoding'): BestScarResult {
  const scored: ScarCandidate[] = candidates.map(oh => ({
    overhang: oh.toUpperCase(),
    ...scoreScarSequence(oh, context),
  }));

  scored.sort((a, b) => b.score - a.score);

  return {
    context,
    candidates: scored,
    best: scored[0] || null,
    preferred: scored.filter(s => s.isPreferred),
    acceptable: scored.filter(s => s.quality === 'acceptable' || s.quality === 'good'),
    poor: scored.filter(s => s.quality === 'poor'),
  };
}

/**
 * Get context-specific recommendations for an assembly
 */
export function getContextRecommendations(assemblyInfo: AssemblyInfo): ContextRecommendation {
  const {
    isCodingSequence = false,
    hasStartCodon = false,
    hasStopCodon = false,
    isLinkerRegion = false,
    junctionPositions = [],
  } = assemblyInfo;

  const recommendations: JunctionRecommendation[] = junctionPositions.map((pos, idx) => {
    let context: CodonContext = 'nonCoding';

    if (isCodingSequence) {
      if (idx === 0 && hasStartCodon) {
        context = 'startCodon';
      } else if (idx === junctionPositions.length - 1 && hasStopCodon) {
        context = 'stopCodon';
      } else if (isLinkerRegion) {
        context = 'linker';
      } else {
        context = 'coding';
      }
    }

    return {
      junctionIndex: idx,
      position: pos,
      recommendedContext: context,
      preferredOverhangs: SCAR_PREFERENCES[context]?.preferred || [],
      avoidOverhangs: SCAR_PREFERENCES[context]?.avoid || [],
    };
  });

  return {
    assemblyInfo,
    recommendations,
    summary: `${recommendations.filter(r => r.recommendedContext !== 'nonCoding').length} junctions have specific scar preferences`,
  };
}

export { CODON_TO_AA, AA_PROPERTIES };
