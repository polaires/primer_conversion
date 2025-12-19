/**
 * Fusion Site Composite Scorer
 *
 * Comprehensive scoring function that evaluates fusion sites considering ALL factors:
 * 1. Overhang quality (fidelity + efficiency)
 * 2. Primer quality at junction position
 * 3. Context quality (flanking, secondary structure)
 * 4. Risk factors (site creation, mispriming)
 * 5. Biological context (coding frames, domains, scars)
 *
 * This is the core scoring engine for the Fusion Site Optimizer.
 */

import { reverseComplement } from './enzymes.js';
import { getEnzymeLigationData, getOverhangFidelityExperimental } from './goldengate.js';
import { calculateEfficiency, countGC } from './overhang-efficiency.js';
import { checkSiteCreation } from './site-creation-check.js';
import { scoreScarSequence } from './scar-preferences.js';
import { selectOptimalFlankingSequence } from './goldengate-primer-optimizer.js';
import { calculateHairpinDG, calculateHomodimerDG } from '../equilibrium.js';
import { calculateTmQ5, calculate3primeTerminalDG } from '../tmQ5.js';
import {
  scoreTm,
  scoreGc,
  scoreHairpin,
  scoreHomodimer,
  scoreLength,
  scoreGcClamp,
  scoreHomopolymer,
  scoreTerminal3DG,
  scoreGQuadruplex,
} from '../scoring.js';

/**
 * Default weights for fusion site scoring
 */
export const DEFAULT_FUSION_WEIGHTS = {
  overhangQuality: 0.20,
  forwardPrimer: 0.20,
  reversePrimer: 0.20,
  riskFactors: 0.25,
  biologicalContext: 0.15,
};

/**
 * Length scoring options
 */
interface LengthOptions {
  optimalLow: number;
  optimalHigh: number;
  acceptableLow: number;
  acceptableHigh: number;
}

/**
 * Score breakdown
 */
interface ScoreBreakdown {
  tm: number;
  gc: number;
  length: number;
  gcClamp: number;
  homopolymer: number;
  hairpin: number;
  homodimer: number;
  terminal3DG: number;
  gQuadruplex: number;
}

/**
 * Homology region analysis result
 */
interface HomologyAnalysis {
  score: number;
  tm: number;
  gc: number;
  length: number;
  gcClamp: boolean;
  hairpinDG: number;
  homodimerDG: number;
  terminal3DG: number;
  hasG4: boolean;
  issues: string[];
  breakdown: ScoreBreakdown;
  region?: string;
  weight?: number;
}

/**
 * Get overhang fidelity from experimental data
 * Uses centralized function from goldengate.js
 * @param overhang - 4bp overhang
 * @param enzyme - Enzyme name
 * @returns Fidelity (0-1)
 */
function getOverhangFidelity(overhang: string, enzyme: string): number {
  const result = getOverhangFidelityExperimental(overhang, enzyme);
  return result.fidelity;
}

/**
 * Analyze a homology region for primer quality
 *
 * @param region - DNA sequence of the homology region
 * @param isForward - True if forward primer (5'->3'), false for reverse
 * @returns Primer quality analysis
 */
function analyzeHomologyRegion(region: string, isForward: boolean): HomologyAnalysis {
  if (!region || region.length < 15) {
    return {
      score: 50,
      tm: 0,
      gc: 0,
      length: 0,
      gcClamp: false,
      hairpinDG: 0,
      homodimerDG: 0,
      terminal3DG: 0,
      hasG4: false,
      issues: ['Region too short'],
      breakdown: {
        tm: 0,
        gc: 0,
        length: 0,
        gcClamp: 0,
        homopolymer: 0,
        hairpin: 0,
        homodimer: 0,
        terminal3DG: 0,
        gQuadruplex: 0,
      },
    };
  }

  const seq = region.toUpperCase();

  // Calculate basic properties
  const gcCount = countGC(seq);
  const gc = (gcCount / seq.length) * 100;

  // Calculate Tm using proper NEB Q5 nearest-neighbor algorithm
  let tm: number;
  try {
    tm = calculateTmQ5(seq);
  } catch (e) {
    // Fallback to Wallace rule for very short sequences
    if (seq.length < 14) {
      const at = seq.length - gcCount;
      tm = 2 * at + 4 * gcCount;
    } else {
      tm = 64.9 + 41 * (gcCount - 16.4) / seq.length;
    }
  }

  // Calculate thermodynamic properties
  let hairpinDG = 0;
  let homodimerDG = 0;
  let terminal3DG = -8;

  try {
    hairpinDG = calculateHairpinDG(seq);
    homodimerDG = calculateHomodimerDG(seq);
    const termResult = calculate3primeTerminalDG(seq);
    terminal3DG = termResult?.dG ?? -8;
  } catch (e) {
    // Use defaults if calculation fails
  }

  // Score components
  const tmScore = scoreTm(tm);
  const gcScore = scoreGc(gc / 100);
  const lengthScore = scoreLength(seq.length, {
    optimalLow: 18,
    optimalHigh: 25,
    acceptableLow: 15,
    acceptableHigh: 30,
  });
  const gcClampScore = scoreGcClamp(seq);
  const homopolymerScore = scoreHomopolymer(seq);
  const hairpinScore = scoreHairpin(hairpinDG);
  const homodimerScore = scoreHomodimer(homodimerDG);
  const terminal3DGScore = scoreTerminal3DG(terminal3DG);
  const gQuadruplexScore = scoreGQuadruplex(seq);

  // Weighted composite
  const composite = (
    tmScore * 0.15 +
    gcScore * 0.10 +
    lengthScore * 0.10 +
    gcClampScore * 0.10 +
    homopolymerScore * 0.05 +
    hairpinScore * 0.15 +
    homodimerScore * 0.10 +
    terminal3DGScore * 0.15 +
    gQuadruplexScore * 0.10
  );

  // Collect issues
  const issues: string[] = [];
  if (tmScore < 0.7) issues.push(`Tm ${tm.toFixed(1)}°C outside optimal range`);
  if (gcScore < 0.7) issues.push(`GC content ${gc.toFixed(0)}% outside optimal range`);
  if (hairpinScore < 0.5) issues.push(`Strong hairpin potential (ΔG=${hairpinDG.toFixed(1)})`);
  if (gQuadruplexScore < 0.5) issues.push('G-quadruplex risk detected');
  if (terminal3DGScore < 0.5) issues.push(`Weak 3' binding (ΔG=${terminal3DG.toFixed(1)})`);

  return {
    score: Math.round(composite * 100),
    tm,
    gc,
    length: seq.length,
    gcClamp: scoreGcClamp(seq) >= 0.9,
    hairpinDG,
    homodimerDG,
    terminal3DG,
    hasG4: gQuadruplexScore < 0.8,
    issues,
    breakdown: {
      tm: Math.round(tmScore * 100),
      gc: Math.round(gcScore * 100),
      length: Math.round(lengthScore * 100),
      gcClamp: Math.round(gcClampScore * 100),
      homopolymer: Math.round(homopolymerScore * 100),
      hairpin: Math.round(hairpinScore * 100),
      homodimer: Math.round(homodimerScore * 100),
      terminal3DG: Math.round(terminal3DGScore * 100),
      gQuadruplex: Math.round(gQuadruplexScore * 100),
    },
  };
}

/**
 * Calculate risk score from site creation and flanking analysis
 *
 * @param siteCreation - Result from checkSiteCreation()
 * @param flanking - Result from selectOptimalFlankingSequence()
 * @returns Risk score (0-100)
 */
function calculateRiskScore(siteCreation: any, flanking: any): number {
  let score = 100;

  // Site creation risks
  if (siteCreation.hasRisk) {
    const highRisks = siteCreation.risks.filter((r: any) => r.severity === 'high').length;
    const mediumRisks = siteCreation.risks.filter((r: any) => r.severity === 'medium').length;
    score -= highRisks * 30;
    score -= mediumRisks * 15;
  }

  // Flanking quality
  if (flanking?.best) {
    const flankingScore = flanking.best.score || 80;
    score -= (100 - flankingScore) * 0.3;

    // Mispriming risk
    if (flanking.best.mispriming) {
      const mispriming = flanking.best.mispriming;
      if (mispriming.risk === 'high') score -= 25;
      else if (mispriming.risk === 'medium') score -= 15;
      else if (mispriming.risk === 'low') score -= 5;
    }
  }

  return Math.max(0, Math.min(100, score));
}

/**
 * Score codon boundary alignment
 *
 * @param position - Junction position
 * @param codingFrame - Reading frame (0, 1, or 2)
 * @returns Score (0-100)
 */
function scoreCodonBoundary(position: number, codingFrame: number | null): number {
  if (codingFrame === null || codingFrame === undefined) return 100;

  // Calculate position within codon
  const posInCodon = (position - codingFrame) % 3;

  // Best: junction at codon boundary (position 0)
  // Acceptable: position 1 (splits codon but no frameshift)
  // Avoid: position 2 (awkward split)
  if (posInCodon === 0) return 100;
  if (posInCodon === 1) return 80;
  return 60;
}

/**
 * Protein domain
 */
interface ProteinDomain {
  start: number;
  end: number;
  name?: string;
}

/**
 * Score domain integrity (junction should not split domains)
 *
 * @param position - Junction position
 * @param proteinDomains - Array of protein domains
 * @returns Score (0-100)
 */
function scoreDomainIntegrity(position: number, proteinDomains: ProteinDomain[]): number {
  if (!proteinDomains || proteinDomains.length === 0) return 100;

  for (const domain of proteinDomains) {
    if (position > domain.start && position < domain.end) {
      // Junction splits a domain - bad
      return 30;
    }
    // Check if close to domain boundary (within 10bp)
    if (Math.abs(position - domain.start) < 10 || Math.abs(position - domain.end) < 10) {
      return 70; // Acceptable but not ideal
    }
  }

  return 100; // Not in any domain
}

/**
 * Component score with weight
 */
interface ComponentScore {
  score: number;
  weight?: number;
  warnings?: string[];
  issues?: string[];
  [key: string]: any;
}

/**
 * Calculate weighted composite score from component scores
 *
 * @param scores - Object with component scores
 * @returns Weighted composite (0-100)
 */
function calculateWeightedComposite(scores: Record<string, ComponentScore>): number {
  let totalWeight = 0;
  let weightedSum = 0;

  // Use weights stored in each score component, fallback to defaults
  for (const [key, component] of Object.entries(scores)) {
    if (component && typeof component.score === 'number') {
      const weight = component.weight || (DEFAULT_FUSION_WEIGHTS as any)[key] || 0;
      weightedSum += component.score * weight;
      totalWeight += weight;
    }
  }

  if (totalWeight === 0) return 0;
  return weightedSum / totalWeight;
}

/**
 * Collect all warnings from score components
 *
 * @param scores - Object with component scores
 * @returns Array of warning strings
 */
function collectAllWarnings(scores: Record<string, ComponentScore>): string[] {
  const warnings: string[] = [];

  for (const [key, component] of Object.entries(scores)) {
    if (component.warnings) {
      warnings.push(...component.warnings);
    }
    if (component.issues) {
      warnings.push(...component.issues);
    }
  }

  return warnings;
}

/**
 * Scoring options
 */
export interface ScoringOptions {
  weights?: typeof DEFAULT_FUSION_WEIGHTS;
  template?: string | null;
  codingFrame?: number | null;
  proteinDomains?: ProteinDomain[];
  scarContext?: string;
  homologyLength?: number;
}

/**
 * Composite score result
 */
export interface CompositeScore {
  position: number;
  overhang: string;
  reverseComplement: string;
  scores: Record<string, ComponentScore>;
  composite: number;
  compositePercent: string;
  quality: 'excellent' | 'good' | 'acceptable' | 'poor' | 'invalid';
  warnings: string[];
  enzyme: string;
  metadata: {
    homologyLength: number;
    codingFrame: number | null;
    scarContext: string;
    hasDomainData: boolean;
  };
  error?: string;
}

/**
 * Score a fusion site considering ALL factors
 *
 * This is the core scoring function that combines:
 * 1. Overhang quality (fidelity + efficiency)
 * 2. Primer quality at this position
 * 3. Context quality (flanking, secondary structure)
 * 4. Risk factors (site creation, mispriming)
 * 5. Biological context (coding frames, domains)
 *
 * @param sequence - Full DNA sequence
 * @param position - Junction position (start of 4bp overhang)
 * @param enzyme - Enzyme name
 * @param options - Scoring options
 * @returns Comprehensive score breakdown
 */
export function scoreFusionSiteComposite(
  sequence: string,
  position: number,
  enzyme: string,
  options: ScoringOptions = {}
): CompositeScore {
  const {
    weights = DEFAULT_FUSION_WEIGHTS,
    template = null,
    codingFrame = null,
    proteinDomains = [],
    scarContext = 'nonCoding',
    homologyLength = 25,
  } = options;

  const seq = sequence.toUpperCase();
  const overhang = seq.slice(position, position + 4);

  if (overhang.length !== 4) {
    return {
      position,
      overhang: '',
      reverseComplement: '',
      scores: {},
      composite: 0,
      compositePercent: '0%',
      quality: 'invalid',
      warnings: [],
      enzyme,
      metadata: {
        homologyLength,
        codingFrame,
        scarContext,
        hasDomainData: proteinDomains.length > 0,
      },
      error: 'Invalid position - cannot extract 4bp overhang',
    };
  }

  const scores: Record<string, ComponentScore> = {};
  const warnings: string[] = [];

  // Component 1: Overhang Intrinsic Quality
  const fidelityData = getOverhangFidelity(overhang, enzyme);
  const efficiencyData = calculateEfficiency(overhang);

  const overhangScore = Math.round(fidelityData * (efficiencyData as any).efficiency * 100);

  scores.overhangQuality = {
    fidelity: fidelityData,
    efficiency: (efficiencyData as any).efficiency,
    combined: fidelityData * (efficiencyData as any).efficiency,
    score: overhangScore,
    warnings: (efficiencyData as any).warnings,
    weight: weights.overhangQuality || 0.20,
  };

  if ((efficiencyData as any).warnings.length > 0) {
    warnings.push(...(efficiencyData as any).warnings);
  }

  // Component 2: Forward Primer Quality
  const fwdRegionEnd = Math.min(seq.length, position + 4 + homologyLength);
  const fwdRegion = seq.slice(position + 4, fwdRegionEnd);

  if (fwdRegion.length >= 15) {
    const fwdAnalysis = analyzeHomologyRegion(fwdRegion, true);
    scores.forwardPrimer = {
      ...fwdAnalysis,
      region: fwdRegion,
      weight: weights.forwardPrimer || 0.20,
    };
    if (fwdAnalysis.issues.length > 0) {
      warnings.push(...fwdAnalysis.issues.map(i => `Fwd primer: ${i}`));
    }
  } else {
    scores.forwardPrimer = {
      score: 50,
      issues: ['Forward homology region too short'],
      weight: weights.forwardPrimer || 0.20,
    };
  }

  // Component 3: Reverse Primer Quality
  const revRegionStart = Math.max(0, position - homologyLength);
  const revRegion = seq.slice(revRegionStart, position);

  if (revRegion.length >= 15) {
    // Reverse primer binds to complement, so analyze RC
    const revAnalysis = analyzeHomologyRegion(reverseComplement(revRegion), false);
    scores.reversePrimer = {
      ...revAnalysis,
      region: revRegion,
      weight: weights.reversePrimer || 0.20,
    };
    if (revAnalysis.issues.length > 0) {
      warnings.push(...revAnalysis.issues.map(i => `Rev primer: ${i}`));
    }
  } else {
    scores.reversePrimer = {
      score: 50,
      issues: ['Reverse homology region too short'],
      weight: weights.reversePrimer || 0.20,
    };
  }

  // Component 4: Risk Factors
  const siteCreation = checkSiteCreation(seq, position, enzyme);

  // Get optimal flanking
  const downstreamContext = seq.slice(position + 4, Math.min(seq.length, position + 24));
  const flanking = selectOptimalFlankingSequence(enzyme, overhang, downstreamContext, {
    template: template || seq,
  } as any);

  const riskScore = calculateRiskScore(siteCreation, flanking);

  scores.riskFactors = {
    siteCreationRisk: (siteCreation as any).hasRisk ? 0 : 100,
    siteCreationDetails: siteCreation,
    flankingQuality: (flanking as any)?.best?.score || 80,
    mispriming: (flanking as any)?.best?.mispriming || null,  // FIXED: Type assertion
    score: riskScore,
    warnings: (siteCreation as any).risks?.map((r: any) => r.message) || [],
    weight: weights.riskFactors || 0.25,
  };

  if ((siteCreation as any).hasRisk) {
    warnings.push(...(siteCreation as any).risks.map((r: any) => r.message));
  }

  // Component 5: Biological Context
  const codonScore = scoreCodonBoundary(position, codingFrame);
  const domainScore = scoreDomainIntegrity(position, proteinDomains);
  const scarResult = scoreScarSequence(overhang, scarContext as any);  // FIXED: Type assertion for CodonContext

  const bioScore = Math.round(
    codonScore * 0.30 +
    domainScore * 0.35 +
    (scarResult as any).score * 0.35
  );

  scores.biologicalContext = {
    codonBoundary: codonScore,
    domainIntegrity: domainScore,
    scarSequence: (scarResult as any).score,
    scarDetails: scarResult,
    score: bioScore,
    weight: weights.biologicalContext || 0.15,
  };

  if (codonScore < 80) {
    warnings.push(`Junction not on codon boundary (frame ${codingFrame})`);
  }
  if (domainScore < 50) {
    warnings.push('Junction splits a protein domain');
  }
  if ((scarResult as any).notes.length > 0) {
    warnings.push(...(scarResult as any).notes.filter((n: any) => n.type === 'warning').map((n: any) => n.message));
  }

  // Calculate weighted composite
  const composite = calculateWeightedComposite(scores);

  // Determine quality tier
  let quality: 'excellent' | 'good' | 'acceptable' | 'poor';
  if (composite >= 85) quality = 'excellent';
  else if (composite >= 70) quality = 'good';
  else if (composite >= 55) quality = 'acceptable';
  else quality = 'poor';

  return {
    position,
    overhang,
    reverseComplement: reverseComplement(overhang),
    scores,
    composite: Math.round(composite),
    compositePercent: `${Math.round(composite)}%`,
    quality,
    warnings: [...new Set(warnings)], // Deduplicate
    enzyme,
    metadata: {
      homologyLength,
      codingFrame,
      scarContext,
      hasDomainData: proteinDomains.length > 0,
    },
  };
}

/**
 * Multiple fusion sites result
 */
export interface MultipleFusionSitesResult {
  enzyme: string;
  totalPositions: number;
  validPositions: number;
  ranked: CompositeScore[];
  best: CompositeScore;
  worst: CompositeScore;
  statistics: {
    averageScore: number;
    excellent: number;
    good: number;
    acceptable: number;
    poor: number;
  };
}

/**
 * Score multiple fusion sites and rank them
 *
 * @param sequence - DNA sequence
 * @param positions - Array of positions to score
 * @param enzyme - Enzyme name
 * @param options - Scoring options
 * @returns Ranked scores and statistics
 */
export function scoreMultipleFusionSites(
  sequence: string,
  positions: number[],
  enzyme: string,
  options: ScoringOptions = {}
): MultipleFusionSitesResult {
  const scored = positions.map(pos =>
    scoreFusionSiteComposite(sequence, pos, enzyme, options)
  );

  // Sort by composite score
  const ranked = [...scored].sort((a, b) => b.composite - a.composite);

  // Statistics
  const valid = scored.filter(s => !s.error);
  const avgScore = valid.reduce((s, v) => s + v.composite, 0) / valid.length || 0;

  return {
    enzyme,
    totalPositions: positions.length,
    validPositions: valid.length,
    ranked,
    best: ranked[0],
    worst: ranked[ranked.length - 1],
    statistics: {
      averageScore: Math.round(avgScore),
      excellent: valid.filter(s => s.quality === 'excellent').length,
      good: valid.filter(s => s.quality === 'good').length,
      acceptable: valid.filter(s => s.quality === 'acceptable').length,
      poor: valid.filter(s => s.quality === 'poor').length,
    },
  };
}

/**
 * Quick score result
 */
export interface QuickScoreResult {
  position: number;
  overhang: string;
  score: number;
  fidelity: number;
  efficiency: number;
  valid: boolean;
  isPalindrome?: boolean;
  isHomopolymer?: boolean;
  isTNNA?: boolean;
}

/**
 * Quick score a position (faster, less detailed)
 *
 * @param sequence - DNA sequence
 * @param position - Position to score
 * @param enzyme - Enzyme name
 * @returns Quick score result
 */
export function quickScoreFusionSite(sequence: string, position: number, enzyme: string): QuickScoreResult {
  const seq = sequence.toUpperCase();
  const overhang = seq.slice(position, position + 4);

  if (overhang.length !== 4) {
    return { position, overhang: '', score: 0, fidelity: 0, efficiency: 0, valid: false };
  }

  const fidelity = getOverhangFidelity(overhang, enzyme);
  const efficiency = calculateEfficiency(overhang);

  const score = Math.round(fidelity * (efficiency as any).efficiency * 100);

  return {
    position,
    overhang,
    score,
    fidelity,
    efficiency: (efficiency as any).efficiency,
    valid: true,
    isPalindrome: (efficiency as any).isPalindrome,
    isHomopolymer: (efficiency as any).isHomopolymer,
    isTNNA: (efficiency as any).isTNNA,
  };
}

export {
  analyzeHomologyRegion,
  calculateRiskScore,
  scoreCodonBoundary,
  scoreDomainIntegrity,
};
