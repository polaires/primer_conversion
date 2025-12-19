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
  overhangQuality: 0.20,     // Fidelity + efficiency
  forwardPrimer: 0.20,       // Primer quality downstream
  reversePrimer: 0.20,       // Primer quality upstream
  riskFactors: 0.25,         // Site creation, mispriming
  biologicalContext: 0.15,   // Codons, domains, scars
};

/**
 * Get overhang fidelity from experimental data
 * Uses centralized function from goldengate.js
 * @param {string} overhang - 4bp overhang
 * @param {string} enzyme - Enzyme name
 * @returns {number} Fidelity (0-1)
 */
function getOverhangFidelity(overhang, enzyme) {
  const result = getOverhangFidelityExperimental(overhang, enzyme);
  return result.fidelity;
}

/**
 * Analyze a homology region for primer quality
 *
 * @param {string} region - DNA sequence of the homology region
 * @param {boolean} isForward - True if forward primer (5'->3'), false for reverse
 * @returns {Object} Primer quality analysis
 */
function analyzeHomologyRegion(region, isForward) {
  if (!region || region.length < 15) {
    return {
      score: 50,
      tm: 0,
      gc: 0,
      issues: ['Region too short'],
    };
  }

  const seq = region.toUpperCase();

  // Calculate basic properties
  const gcCount = countGC(seq);
  const gc = (gcCount / seq.length) * 100;

  // Calculate Tm using proper NEB Q5 nearest-neighbor algorithm
  let tm;
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
  const issues = [];
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
 * @param {Object} siteCreation - Result from checkSiteCreation()
 * @param {Object} flanking - Result from selectOptimalFlankingSequence()
 * @returns {number} Risk score (0-100)
 */
function calculateRiskScore(siteCreation, flanking) {
  let score = 100;

  // Site creation risks
  if (siteCreation.hasRisk) {
    const highRisks = siteCreation.risks.filter(r => r.severity === 'high').length;
    const mediumRisks = siteCreation.risks.filter(r => r.severity === 'medium').length;
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
 * @param {number} position - Junction position
 * @param {number} codingFrame - Reading frame (0, 1, or 2)
 * @returns {number} Score (0-100)
 */
function scoreCodonBoundary(position, codingFrame) {
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
 * Score domain integrity (junction should not split domains)
 *
 * @param {number} position - Junction position
 * @param {Array} proteinDomains - Array of {start, end, name} objects
 * @returns {number} Score (0-100)
 */
function scoreDomainIntegrity(position, proteinDomains) {
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
 * Calculate weighted composite score from component scores
 *
 * @param {Object} scores - Object with component scores
 * @returns {number} Weighted composite (0-100)
 */
function calculateWeightedComposite(scores) {
  let totalWeight = 0;
  let weightedSum = 0;

  // Use weights stored in each score component, fallback to defaults
  for (const [key, component] of Object.entries(scores)) {
    if (component && typeof component.score === 'number') {
      const weight = component.weight || DEFAULT_FUSION_WEIGHTS[key] || 0;
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
 * @param {Object} scores - Object with component scores
 * @returns {Array} Array of warning strings
 */
function collectAllWarnings(scores) {
  const warnings = [];

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
 * Score a fusion site considering ALL factors
 *
 * This is the core scoring function that combines:
 * 1. Overhang quality (fidelity + efficiency)
 * 2. Primer quality at this position
 * 3. Context quality (flanking, secondary structure)
 * 4. Risk factors (site creation, mispriming)
 * 5. Biological context (coding frames, domains)
 *
 * @param {string} sequence - Full DNA sequence
 * @param {number} position - Junction position (start of 4bp overhang)
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Scoring options
 * @returns {Object} Comprehensive score breakdown
 */
export function scoreFusionSiteComposite(sequence, position, enzyme, options = {}) {
  const {
    weights = DEFAULT_FUSION_WEIGHTS,
    template = null,              // Full template for context checks
    codingFrame = null,           // For codon boundary scoring
    proteinDomains = [],          // For domain integrity scoring
    scarContext = 'nonCoding',    // For scar preferences
    homologyLength = 25,          // Length of primer homology region
  } = options;

  const seq = sequence.toUpperCase();
  const overhang = seq.slice(position, position + 4);

  if (overhang.length !== 4) {
    return {
      position,
      overhang: '',
      error: 'Invalid position - cannot extract 4bp overhang',
      composite: 0,
      quality: 'invalid',
    };
  }

  const scores = {};
  const warnings = [];

  // ═══════════════════════════════════════════════════════════════════
  // COMPONENT 1: Overhang Intrinsic Quality (Weight: 20%)
  // ═══════════════════════════════════════════════════════════════════
  const fidelityData = getOverhangFidelity(overhang, enzyme);
  const efficiencyData = calculateEfficiency(overhang);

  const overhangScore = Math.round(fidelityData * efficiencyData.efficiency * 100);

  scores.overhangQuality = {
    fidelity: fidelityData,
    efficiency: efficiencyData.efficiency,
    combined: fidelityData * efficiencyData.efficiency,
    score: overhangScore,
    warnings: efficiencyData.warnings,
    weight: weights.overhangQuality || 0.20,
  };

  if (efficiencyData.warnings.length > 0) {
    warnings.push(...efficiencyData.warnings);
  }

  // ═══════════════════════════════════════════════════════════════════
  // COMPONENT 2: Forward Primer Quality (Weight: 20%)
  // ═══════════════════════════════════════════════════════════════════
  // Forward primer binds downstream of junction
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

  // ═══════════════════════════════════════════════════════════════════
  // COMPONENT 3: Reverse Primer Quality (Weight: 20%)
  // ═══════════════════════════════════════════════════════════════════
  // Reverse primer binds upstream of junction
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

  // ═══════════════════════════════════════════════════════════════════
  // COMPONENT 4: Risk Factors (Weight: 25%)
  // ═══════════════════════════════════════════════════════════════════
  const siteCreation = checkSiteCreation(seq, position, enzyme);

  // Get optimal flanking
  const downstreamContext = seq.slice(position + 4, Math.min(seq.length, position + 24));
  const flanking = selectOptimalFlankingSequence(enzyme, overhang, downstreamContext, {
    template: template || seq,
  });

  const riskScore = calculateRiskScore(siteCreation, flanking);

  scores.riskFactors = {
    siteCreationRisk: siteCreation.hasRisk ? 0 : 100,
    siteCreationDetails: siteCreation,
    flankingQuality: flanking?.best?.score || 80,
    mispriming: flanking?.best?.mispriming || null,
    score: riskScore,
    warnings: siteCreation.risks?.map(r => r.message) || [],
    weight: weights.riskFactors || 0.25,
  };

  if (siteCreation.hasRisk) {
    warnings.push(...siteCreation.risks.map(r => r.message));
  }

  // ═══════════════════════════════════════════════════════════════════
  // COMPONENT 5: Biological Context (Weight: 15%)
  // ═══════════════════════════════════════════════════════════════════
  const codonScore = scoreCodonBoundary(position, codingFrame);
  const domainScore = scoreDomainIntegrity(position, proteinDomains);
  const scarResult = scoreScarSequence(overhang, scarContext);

  const bioScore = Math.round(
    codonScore * 0.30 +
    domainScore * 0.35 +
    scarResult.score * 0.35
  );

  scores.biologicalContext = {
    codonBoundary: codonScore,
    domainIntegrity: domainScore,
    scarSequence: scarResult.score,
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
  if (scarResult.notes.length > 0) {
    warnings.push(...scarResult.notes.filter(n => n.type === 'warning').map(n => n.message));
  }

  // ═══════════════════════════════════════════════════════════════════
  // Calculate weighted composite
  // ═══════════════════════════════════════════════════════════════════
  const composite = calculateWeightedComposite(scores);

  // Determine quality tier
  let quality;
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
 * Score multiple fusion sites and rank them
 *
 * @param {string} sequence - DNA sequence
 * @param {number[]} positions - Array of positions to score
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Scoring options
 * @returns {Object} Ranked scores and statistics
 */
export function scoreMultipleFusionSites(sequence, positions, enzyme, options = {}) {
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
 * Quick score a position (faster, less detailed)
 *
 * @param {string} sequence - DNA sequence
 * @param {number} position - Position to score
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Quick score result
 */
export function quickScoreFusionSite(sequence, position, enzyme) {
  const seq = sequence.toUpperCase();
  const overhang = seq.slice(position, position + 4);

  if (overhang.length !== 4) {
    return { position, score: 0, valid: false };
  }

  const fidelity = getOverhangFidelity(overhang, enzyme);
  const efficiency = calculateEfficiency(overhang);

  const score = Math.round(fidelity * efficiency.efficiency * 100);

  return {
    position,
    overhang,
    score,
    fidelity,
    efficiency: efficiency.efficiency,
    valid: true,
    isPalindrome: efficiency.isPalindrome,
    isHomopolymer: efficiency.isHomopolymer,
    isTNNA: efficiency.isTNNA,
  };
}

export {
  analyzeHomologyRegion,
  calculateRiskScore,
  scoreCodonBoundary,
  scoreDomainIntegrity,
};
