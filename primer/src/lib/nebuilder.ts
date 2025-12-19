/**
 * NEBuilder HiFi DNA Assembly Module
 *
 * Specialized module for NEB's NEBuilder HiFi DNA Assembly Kit (E5520).
 * Implements advanced overlap optimization based on NEB's recommendations.
 *
 * References:
 * - NEB NEBuilder HiFi DNA Assembly Master Mix (E5520) manual
 * - NEB Online Assembly Tool algorithms
 *
 * This module provides:
 * - NEBuilder-specific Tm calculations
 * - Advanced overlap optimization with multiple scoring criteria
 * - Assembly simulation and yield prediction
 * - Protocol generation matching NEB format
 */

import { calculateTmQ5, calculateGC } from './tmQ5.js';
import { dg as foldDG } from './fold.js';
import { calculateHairpinDG, calculateHomodimerDG, calculateHeterodimerDG } from './equilibrium.js';
import { primers } from './primers.js';
import { reverseComplement } from './sequenceUtils.js';
import { scoreHairpin as unifiedScoreHairpin } from './scoring.js';
import { ANALYSIS_PRESETS } from './presets.js';

// =============================================================================
// Types and Interfaces
// =============================================================================

export interface NEBuilderParams {
  overlap: {
    minLength: number;
    maxLength: number;
    optimalLength: number;
    minTm: number;
    maxTm: number;
    optimalTm: number;
  };
  fragmentAmounts: {
    '2-3': { min: number; max: number; recommended: number };
    '4-6': { min: number; max: number; recommended: number };
  };
  incubation: {
    '2-3': number;
    '4+': number;
    max: number;
  };
  molarRatios: Record<string, string>;
  reactionTemp: number;
  masterMixVolume: number;
}

export interface OverlapAnalysis {
  sequence: string;
  length: number;
  tm: number;
  gc: number;
  gcPercent: string;
  hairpinDG: number;
  hasGcClamp: boolean;
  patterns: {
    polyT: boolean;
    polyA: boolean;
    polyG: boolean;
    polyC: boolean;
    dinucRepeat: boolean;
    lowComplexity: boolean;
    palindrome: boolean;
  };
  scores: {
    tm: number;
    gcContent: number;
    length: number;
    hairpin: number;
    gcClamp: number;
    patterns: number;
  };
  compositeScore: number;
  quality: string;
  warnings: string[];
  recommendation: string;
}

export interface OverlapOptimization {
  optimal: OverlapAnalysis;
  alternatives: OverlapAnalysis[];
  totalCandidatesEvaluated: number;
  searchSpace: {
    lengthRange: string;
    positionRange: string;
  };
}

export interface Fragment {
  id: string;
  seq: string;
}

export interface NEBuilderDesign {
  method: string;
  version: string;
  fragments: Array<{
    id: string;
    index: number;
    length: number;
    forward: {
      sequence: string;
      length: number;
      annealingRegion: string;
      annealingLength: number;
      homologyTail: string;
      tm: number;
      hairpinDG: number;
    };
    reverse: {
      sequence: string;
      length: number;
      annealingRegion: string;
      annealingLength: number;
      homologyTail: string;
      tm: number;
      hairpinDG: number;
    };
    pair: {
      tmDiff: number;
      heterodimerDG: number;
      annealingTemp: number;
    };
    warnings: string[];
  }>;
  junctions: Array<{
    index: number;
    from: { id: string; position: string };
    to: { id: string; position: string };
    optimal: OverlapAnalysis;
    alternatives: OverlapAnalysis[];
  }>;
  assembly: {
    type: string;
    totalLength: number;
    fragmentCount: number;
  };
  quality: {
    avgJunctionScore: number;
    minJunctionScore: number;
    tier: string;
    recommendation: string;
  };
  protocol: Protocol;
  warnings: string[];
  exportFormats: string[];
}

export interface Protocol {
  title: string;
  kit: string;
  version: string;
  materials: string[];
  steps: Array<{
    step: number;
    title: string;
    substeps: string[];
    tips?: string[];
    table?: {
      headers: string[];
      rows: string[][];
    };
    critical?: string;
  }>;
  troubleshooting: Array<{
    problem: string;
    solutions: string[];
  }>;
  references: string[];
}

// =============================================================================
// NEBuilder-Specific Constants
// =============================================================================

/**
 * NEBuilder HiFi recommended parameters (from NEB manual)
 */
export const NEBUILDER_PARAMS: NEBuilderParams = {
  // Overlap requirements
  overlap: {
    minLength: 15,
    maxLength: 35,
    optimalLength: 20,
    minTm: 48,
    maxTm: 65,
    optimalTm: 55,
  },

  // Fragment input amounts (pmol)
  fragmentAmounts: {
    '2-3': { min: 0.03, max: 0.2, recommended: 0.1 },
    '4-6': { min: 0.2, max: 0.5, recommended: 0.3 },
  },

  // Incubation times (minutes)
  incubation: {
    '2-3': 15,
    '4+': 60,
    max: 60,
  },

  // Vector:Insert molar ratios
  molarRatios: {
    '1:1': 'Standard',
    '1:2': 'Recommended for difficult inserts',
    '1:3': 'For very small inserts',
  },

  // Temperature
  reactionTemp: 50, // °C

  // Master mix volume
  masterMixVolume: 10, // µL per 20 µL reaction
};

/**
 * Scoring weights for overlap optimization
 * Based on experimental data from NEB and published literature
 */
export const OVERLAP_SCORING_WEIGHTS: Record<string, number> = {
  tm: 0.30,           // Tm near target is critical
  gcContent: 0.15,    // Balanced GC important
  length: 0.10,       // Optimal length preference
  hairpin: 0.15,      // Secondary structure avoidance
  gcClamp: 0.10,      // 3' end stability
  patterns: 0.20,     // Avoid problematic sequences
};

// =============================================================================
// Advanced Overlap Analysis
// =============================================================================

/**
 * Comprehensive overlap analysis for NEBuilder
 */
export function analyzeOverlap(
  seq: string,
  options: { targetTm?: number; targetLength?: number } = {}
): OverlapAnalysis {
  const { targetTm = 55, targetLength = 20 } = options;

  seq = seq.toUpperCase();
  const len = seq.length;

  // Calculate basic properties
  const tm = calculateTmQ5(seq);
  const gc = calculateGC(seq);

  // Calculate secondary structure
  let hairpinDG = 0;
  try {
    hairpinDG = foldDG(seq, 50);
  } catch (e) {
    hairpinDG = 0;
  }

  // Analyze 3' end (GC clamp)
  const last2 = seq.slice(-2);
  const gcLast2 = (last2.match(/[GC]/g) || []).length;
  const hasGcClamp = gcLast2 >= 1;

  // Check for problematic patterns
  const patterns = {
    polyT: /TTTT/.test(seq),          // Poor 3' extension
    polyA: /AAAA/.test(seq),          // AT-rich run
    polyG: /GGGG/.test(seq),          // G-quadruplex risk
    polyC: /CCCC/.test(seq),          // GC-rich run
    dinucRepeat: /(.)\1{3,}/.test(seq), // Any 4+ repeat
    lowComplexity: /(.)\1\1/.test(seq.slice(-4)), // Repeat at 3' end
    palindrome: isPalindromic(seq),   // Self-complementary
  };

  // Calculate individual scores (0-1 scale)
  const scores = {
    // Tm score: optimal at target, penalty for deviation
    tm: scoreTmDeviation(tm, targetTm, {
      optimal: 0, acceptable: 5, steep: 10
    }),

    // GC content score: optimal 40-60%
    gcContent: scoreGcContent(gc),

    // Length score: optimal at target
    length: scoreLengthDeviation(len, targetLength, {
      optimalRange: 5, acceptable: 10
    }),

    // Hairpin score: penalty for stable structures
    hairpin: scoreHairpin(hairpinDG),

    // GC clamp score: prefer 1 G/C in last 2 bp
    gcClamp: gcLast2 === 1 ? 1.0 : gcLast2 === 2 ? 0.85 : 0.5,

    // Pattern score: penalty for problematic sequences
    patterns: scorePatterns(patterns),
  };

  // Calculate weighted composite score
  let compositeScore = 0;
  for (const [key, score] of Object.entries(scores)) {
    compositeScore += score * (OVERLAP_SCORING_WEIGHTS[key] || 0);
  }
  compositeScore = Math.round(compositeScore * 100);

  // Generate warnings
  const warnings: string[] = [];
  if (tm < NEBUILDER_PARAMS.overlap.minTm) {
    warnings.push(`Tm ${tm.toFixed(1)}°C below minimum ${NEBUILDER_PARAMS.overlap.minTm}°C`);
  }
  if (tm > NEBUILDER_PARAMS.overlap.maxTm) {
    warnings.push(`Tm ${tm.toFixed(1)}°C above maximum ${NEBUILDER_PARAMS.overlap.maxTm}°C`);
  }
  if (patterns.polyT) warnings.push('Contains poly-T run - may cause poor extension');
  if (patterns.polyG) warnings.push('Contains poly-G run - G-quadruplex risk');
  if (patterns.palindrome) warnings.push('Palindromic - may cause self-annealing');
  if (hairpinDG < -3) warnings.push(`Stable hairpin structure (ΔG = ${hairpinDG.toFixed(1)} kcal/mol)`);
  if (!hasGcClamp) warnings.push('No GC clamp at 3\' end');

  // Classify quality
  const quality = compositeScore >= 80 ? 'excellent' :
                  compositeScore >= 60 ? 'good' :
                  compositeScore >= 40 ? 'acceptable' : 'poor';

  return {
    sequence: seq,
    length: len,
    tm: Math.round(tm * 10) / 10,
    gc: Math.round(gc * 1000) / 10,
    gcPercent: `${(gc * 100).toFixed(1)}%`,
    hairpinDG: Math.round(hairpinDG * 10) / 10,
    hasGcClamp,
    patterns,
    scores,
    compositeScore,
    quality,
    warnings,
    recommendation: quality === 'poor' ?
      'Consider adjusting overlap position or length' :
      quality === 'acceptable' ?
      'May work but consider alternatives' :
      'Good for assembly',
  };
}

/**
 * Find optimal overlap with sliding window
 */
export function optimizeOverlap(
  junction: string,
  options: {
    minLen?: number;
    maxLen?: number;
    targetTm?: number;
    targetLength?: number;
  } = {}
): OverlapOptimization {
  const {
    minLen = NEBUILDER_PARAMS.overlap.minLength,
    maxLen = NEBUILDER_PARAMS.overlap.maxLength,
    targetTm = NEBUILDER_PARAMS.overlap.optimalTm,
    targetLength = NEBUILDER_PARAMS.overlap.optimalLength,
  } = options;

  junction = junction.toUpperCase();
  const junctionCenter = Math.floor(junction.length / 2);

  const candidates: OverlapAnalysis[] = [];

  // Try different lengths
  for (let len = minLen; len <= maxLen; len++) {
    // Try different positions centered around junction
    const halfLen = Math.floor(len / 2);
    const startRange = Math.max(0, junctionCenter - halfLen - 5);
    const endRange = Math.min(junction.length - len, junctionCenter - halfLen + 5);

    for (let start = startRange; start <= endRange; start++) {
      const overlap = junction.slice(start, start + len);
      const analysis = analyzeOverlap(overlap, { targetTm, targetLength });

      candidates.push({
        ...analysis,
        position: start,
        offset: start - (junctionCenter - halfLen),
      } as OverlapAnalysis & { position: number; offset: number });
    }
  }

  // Sort by composite score
  candidates.sort((a, b) => b.compositeScore - a.compositeScore);

  const best = candidates[0];
  const alternatives = candidates.slice(1, 5).filter(c => c.compositeScore > 50);

  return {
    optimal: best,
    alternatives,
    totalCandidatesEvaluated: candidates.length,
    searchSpace: {
      lengthRange: `${minLen}-${maxLen} bp`,
      positionRange: 'centered ±5 bp',
    },
  };
}

// =============================================================================
// Assembly Design with NEBuilder Optimization
// =============================================================================

/**
 * Design NEBuilder HiFi assembly with optimized overlaps
 */
export function designNEBuilderAssembly(
  fragments: Fragment[],
  options: {
    circular?: boolean;
    targetOverlapTm?: number;
    targetOverlapLength?: number;
  } = {}
): NEBuilderDesign {
  const {
    circular = true,
    targetOverlapTm = 55,
    targetOverlapLength = 20,
  } = options;

  if (fragments.length < 2) {
    throw new Error('Assembly requires at least 2 fragments');
  }
  if (fragments.length > 6) {
    throw new Error('NEBuilder HiFi supports maximum 6 fragments efficiently');
  }

  // Analyze and optimize each junction
  const junctions: NEBuilderDesign['junctions'] = [];

  for (let i = 0; i < fragments.length; i++) {
    const frag1 = fragments[i];
    const frag2 = fragments[(i + 1) % fragments.length];

    // Skip last junction for linear assembly
    if (!circular && i === fragments.length - 1) continue;

    // Create junction sequence (end of frag1 + start of frag2)
    const junctionSeq = frag1.seq.slice(-50) + frag2.seq.slice(0, 50);

    const optimization = optimizeOverlap(junctionSeq, {
      targetTm: targetOverlapTm,
      targetLength: targetOverlapLength,
    });

    junctions.push({
      index: i + 1,
      from: { id: frag1.id, position: 'end' },
      to: { id: frag2.id, position: 'start' },
      optimal: optimization.optimal,
      alternatives: optimization.alternatives,
    });
  }

  // Design primers for each fragment
  const fragmentDesigns: NEBuilderDesign['fragments'] = [];

  for (let i = 0; i < fragments.length; i++) {
    const frag = fragments[i];

    // Get overlap sequences
    let leftOverlap = '';
    let rightOverlap = '';

    if (circular || i > 0) {
      const prevJunction = junctions[(i - 1 + junctions.length) % junctions.length];
      leftOverlap = prevJunction.optimal.sequence;
    }

    if (circular || i < fragments.length - 1) {
      const nextJunction = junctions[i % junctions.length];
      rightOverlap = nextJunction.optimal.sequence;
    }

    // Design base primers using existing primers() function
    const [fwd, rev] = primers(frag.seq, {
      optimalTm: 62,
      useCompositeScore: true,
    });

    // Add homology tails
    const fwdWithTail = leftOverlap + fwd.seq;
    const revWithTail = reverseComplement(rightOverlap) + rev.seq;

    // Calculate properties of full primers
    const fwdTm = calculateTmQ5(fwd.seq);
    const revTm = calculateTmQ5(rev.seq);
    const fwdHairpin = calculateHairpinDG(fwdWithTail, 55);
    const revHairpin = calculateHairpinDG(revWithTail, 55);
    const heterodimer = calculateHeterodimerDG(fwdWithTail, revWithTail, 55);

    // Check for issues
    const warnings: string[] = [];
    if (fwdWithTail.length > 60) warnings.push('Forward primer >60bp - consider shorter overlap');
    if (revWithTail.length > 60) warnings.push('Reverse primer >60bp - consider shorter overlap');
    if (fwdHairpin < -4) warnings.push(`Forward primer hairpin risk (ΔG = ${fwdHairpin.toFixed(1)})`);
    if (revHairpin < -4) warnings.push(`Reverse primer hairpin risk (ΔG = ${revHairpin.toFixed(1)})`);
    if (heterodimer < -6) warnings.push(`Primer dimer risk (ΔG = ${heterodimer.toFixed(1)})`);

    fragmentDesigns.push({
      id: frag.id,
      index: i + 1,
      length: frag.seq.length,
      forward: {
        sequence: fwdWithTail,
        length: fwdWithTail.length,
        annealingRegion: fwd.seq,
        annealingLength: fwd.seq.length,
        homologyTail: leftOverlap,
        tm: fwdTm,
        hairpinDG: fwdHairpin,
      },
      reverse: {
        sequence: revWithTail,
        length: revWithTail.length,
        annealingRegion: rev.seq,
        annealingLength: rev.seq.length,
        homologyTail: reverseComplement(rightOverlap),
        tm: revTm,
        hairpinDG: revHairpin,
      },
      pair: {
        tmDiff: Math.abs(fwdTm - revTm),
        heterodimerDG: heterodimer,
        annealingTemp: Math.min(fwdTm, revTm),
      },
      warnings,
    });
  }

  // Calculate assembly quality metrics
  const junctionScores = junctions.map(j => j.optimal.compositeScore);
  const avgJunctionScore = junctionScores.reduce((a, b) => a + b, 0) / junctionScores.length;
  const minJunctionScore = Math.min(...junctionScores);

  // Generate protocol
  const protocol = generateNEBuilderProtocol(fragmentDesigns.length);

  // Collect all warnings
  const allWarnings = [
    ...junctions.flatMap(j => j.optimal.warnings.map(w =>
      `Junction ${j.index} (${j.from.id}→${j.to.id}): ${w}`
    )),
    ...fragmentDesigns.flatMap(f => f.warnings.map(w =>
      `Fragment ${f.index} (${f.id}): ${w}`
    )),
  ];

  return {
    method: 'NEBuilder HiFi DNA Assembly',
    version: '2.0',
    fragments: fragmentDesigns,
    junctions,
    assembly: {
      type: circular ? 'circular' : 'linear',
      totalLength: fragments.reduce((sum, f) => sum + f.seq.length, 0),
      fragmentCount: fragments.length,
    },
    quality: {
      avgJunctionScore: Math.round(avgJunctionScore),
      minJunctionScore,
      tier: avgJunctionScore >= 70 ? 'excellent' :
            avgJunctionScore >= 50 ? 'good' : 'acceptable',
      recommendation: avgJunctionScore >= 70 ?
        'High probability of successful assembly' :
        avgJunctionScore >= 50 ?
        'Assembly should work, check warnings' :
        'Consider redesigning problematic junctions',
    },
    protocol,
    warnings: allWarnings,
    exportFormats: ['tsv', 'csv', 'json', 'genbank'],
  };
}

// =============================================================================
// Protocol Generation
// =============================================================================

/**
 * Generate NEBuilder HiFi protocol matching NEB format
 */
export function generateNEBuilderProtocol(fragmentCount: number): Protocol {
  const incubationTime = fragmentCount <= 3 ? 15 : 60;
  const pmolRange = fragmentCount <= 3 ?
    NEBUILDER_PARAMS.fragmentAmounts['2-3'] :
    NEBUILDER_PARAMS.fragmentAmounts['4-6'];

  return {
    title: 'NEBuilder HiFi DNA Assembly Protocol',
    kit: 'NEB #E5520',
    version: 'Rev. 2024',

    materials: [
      'NEBuilder HiFi DNA Assembly Master Mix (2X)',
      'Purified DNA fragments (PCR products or restriction digests)',
      'Nuclease-free water',
      'PCR tubes or strips',
    ],

    steps: [
      {
        step: 1,
        title: 'Prepare DNA fragments',
        substeps: [
          'Amplify each fragment using high-fidelity PCR (Q5 recommended)',
          'Verify correct size by gel electrophoresis',
          'Purify PCR products using spin columns',
          'Quantify DNA concentration using Qubit or Nanodrop',
          'Calculate molar amounts: pmol = (ng × 1000) / (bp × 650)',
        ],
        tips: [
          'Use ≥15 bp overlap between adjacent fragments',
          'Optimal overlap Tm is 48-65°C',
          'DpnI digest can remove template plasmid',
        ],
      },
      {
        step: 2,
        title: 'Set up assembly reaction',
        substeps: [
          `Mix ${pmolRange.min}-${pmolRange.max} pmol total DNA (recommended: ${pmolRange.recommended} pmol)`,
          'For vector:insert, use 1:2 molar ratio',
          'Add 10 µL NEBuilder HiFi DNA Assembly Master Mix (2X)',
          'Add nuclease-free water to 20 µL total',
          'Mix gently by pipetting',
        ],
        table: {
          headers: ['Component', 'Volume'],
          rows: [
            ['DNA fragments (total)', 'X µL'],
            ['NEBuilder HiFi Master Mix (2X)', '10 µL'],
            ['Nuclease-free water', 'to 20 µL'],
          ],
        },
      },
      {
        step: 3,
        title: 'Incubate',
        substeps: [
          `Incubate at 50°C for ${incubationTime} minutes`,
          fragmentCount <= 3 ?
            'For 2-3 fragments: 15 minutes is sufficient' :
            'For 4+ fragments: incubate full 60 minutes',
          'Reactions can be left at 50°C up to 60 minutes',
        ],
        critical: 'Do not exceed 60 minutes incubation',
      },
      {
        step: 4,
        title: 'Transform',
        substeps: [
          'Place reaction on ice',
          'Transform 2 µL into NEB 10-beta competent cells',
          'Follow standard transformation protocol',
          'Plate on appropriate antibiotic selection',
          'Incubate overnight at 37°C',
        ],
        tips: [
          'Store remaining reaction at -20°C',
          'Expect 10²-10⁴ colonies depending on fragment number',
        ],
      },
    ],

    troubleshooting: [
      {
        problem: 'No colonies',
        solutions: [
          'Verify DNA fragment concentration',
          'Check overlap sequences for errors',
          'Ensure overlaps have Tm ≥48°C',
          'Verify competent cell efficiency',
        ],
      },
      {
        problem: 'Low colony count',
        solutions: [
          'Increase DNA amount',
          'Extend incubation to 60 minutes',
          'Optimize molar ratios',
          'Use longer overlaps (25-30 bp)',
        ],
      },
      {
        problem: 'Wrong assembly',
        solutions: [
          'Screen more colonies',
          'Increase overlap length',
          'Check for internal homologous sequences',
          'Use unique overlaps at each junction',
        ],
      },
    ],

    references: [
      'NEB NEBuilder HiFi DNA Assembly Master Mix Manual (E5520)',
      'Gibson et al. 2009, Nature Methods 6:343-345',
    ],
  };
}

// =============================================================================
// Helper Functions
// =============================================================================

function isPalindromic(seq: string): boolean {
  if (seq.length < 4) return false;
  const rc = reverseComplement(seq);
  return seq === rc;
}

function scoreTmDeviation(tm: number, target: number, { optimal = 0, acceptable = 5, steep = 10 }): number {
  const diff = Math.abs(tm - target);
  if (diff <= optimal) return 1.0;
  if (diff <= acceptable) return 1.0 - (diff - optimal) / (acceptable - optimal) * 0.3;
  return Math.max(0, 0.7 - (diff - acceptable) / (steep - acceptable) * 0.7);
}

function scoreGcContent(gc: number): number {
  if (gc >= 0.40 && gc <= 0.60) return 1.0;
  if (gc >= 0.30 && gc <= 0.70) return 0.7;
  return Math.max(0, 0.4 - Math.abs(gc - 0.5) * 0.8);
}

function scoreLengthDeviation(len: number, target: number, { optimalRange = 5, acceptable = 10 }): number {
  const diff = Math.abs(len - target);
  if (diff <= optimalRange) return 1.0;
  if (diff <= acceptable) return 0.8;
  return Math.max(0.5, 1.0 - diff * 0.05);
}

function scoreHairpin(dg: number): number {
  // Use unified scoring with assembly preset threshold
  return unifiedScoreHairpin(dg, { threshold: ANALYSIS_PRESETS.assembly.hairpinThreshold });
}

function scorePatterns(patterns: Record<string, boolean>): number {
  let score = 1.0;
  if (patterns.polyT) score -= 0.3;
  if (patterns.polyA) score -= 0.2;
  if (patterns.polyG) score -= 0.4;
  if (patterns.polyC) score -= 0.2;
  if (patterns.palindrome) score -= 0.3;
  if (patterns.lowComplexity) score -= 0.2;
  return Math.max(0, score);
}
