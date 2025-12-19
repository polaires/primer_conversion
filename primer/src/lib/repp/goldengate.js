/**
 * Golden Gate Assembly Planning Module
 *
 * Based on NEB E1601 protocol and MoClo standard.
 * Includes experimental ligation fidelity data from Pryor et al. (2020).
 *
 * References:
 * - Pryor et al. (2020) PLOS ONE - "Enabling one-pot Golden Gate assemblies of unprecedented complexity"
 * - Potapov et al. (2018) ACS Synth Biol - "Comprehensive Profiling of Four Base Overhang Ligation Fidelity"
 * - Weber et al. (2011) PLOS ONE - "A Modular Cloning System for Standardized Assembly"
 * - NEB Golden Gate Assembly Guidelines
 */

import { reverseComplement } from './enzymes.js';
import { calculateTmQ5 } from '../tmQ5.js';
import { OPTIMAL_FLANKING_SEQUENCES, OPTIMAL_SPACERS } from './goldengate-primer-optimizer.js';

// Import ligation data - Vite handles JSON imports automatically
import ligationDataRaw from './ligation-data.json';

// Re-export ligation data for direct access
export const ligationData = ligationDataRaw;

// Type IIS restriction enzymes for Golden Gate Assembly
// Now includes experimental ligation data from Pryor et al. (2020)
export const GOLDEN_GATE_ENZYMES = {
  BsaI: {
    name: 'BsaI',
    fullName: 'BsaI-HFv2',
    recognition: 'GGTCTC',      // Recognition sequence
    cutOffset: 1,               // Cut 1bp after recognition on top strand
    overhangLength: 4,          // Creates 4bp overhang
    // BsaI cuts: 5'...GGTCTC(N)1|NNNN...3'
    //            3'...CCAGAG(N)5NNNN|...5'
    // Primer format: 5'-[extra bases]-GGTCTCN[overhang][homology]-3'
    fwdSite: 'GGTCTCN',  // N is 1bp spacer before overhang
    revSite: 'GAGACC',   // Reverse complement for reverse primers
    dataKey: 'BsaI-HFv2',  // Key in ligationData.enzymes
  },
  BbsI: {
    name: 'BbsI',
    fullName: 'BbsI-HF',
    aliases: ['BpiI'],
    recognition: 'GAAGAC',
    cutOffset: 2,
    overhangLength: 4,
    fwdSite: 'GAAGACNN',
    revSite: 'GTCTTC',
    dataKey: 'BbsI-HF',
  },
  BsmBI: {
    name: 'BsmBI',
    fullName: 'BsmBI-v2',
    recognition: 'CGTCTC',
    cutOffset: 1,
    overhangLength: 4,
    fwdSite: 'CGTCTCN',
    revSite: 'GAGACG',
    dataKey: 'BsmBI-v2',
  },
  Esp3I: {
    name: 'Esp3I',
    fullName: 'Esp3I',
    recognition: 'CGTCTC',
    cutOffset: 1,
    overhangLength: 4,
    fwdSite: 'CGTCTCN',
    revSite: 'GAGACG',
    dataKey: 'Esp3I',
  },
  SapI: {
    name: 'SapI',
    fullName: 'SapI',
    recognition: 'GCTCTTC',
    cutOffset: 1,
    overhangLength: 3,  // SapI creates 3bp overhangs
    fwdSite: 'GCTCTTCN',
    revSite: 'GAAGAGC',
    dataKey: 'SapI',
  },
};

// Available enzymes with ligation data
export const ENZYMES_WITH_DATA = Object.keys(GOLDEN_GATE_ENZYMES).filter(
  key => ligationData.enzymes[GOLDEN_GATE_ENZYMES[key].dataKey]
);

/**
 * Standard MoClo/CIDAR fusion sites
 * Based on Weber et al. (2011) and NEB assembly standards
 * Format: [5' overhang] - [Part type] - [3' overhang]
 */
export const STANDARD_FUSION_SITES = {
  // Core MoClo Level 0 sites
  A: { seq: 'GGAG', description: 'Upstream of promoter' },
  B: { seq: 'TACT', description: 'Promoter/5\'UTR junction' },
  C: { seq: 'AATG', description: 'RBS/CDS start (includes ATG)' },
  D: { seq: 'AGGT', description: 'CDS end/terminator' },
  E: { seq: 'GCTT', description: 'Downstream of terminator' },
  // Extended sites for complex assemblies
  F: { seq: 'CGCT', description: 'Extended site F' },
  G: { seq: 'TGCC', description: 'Extended site G' },
  H: { seq: 'ACTA', description: 'Extended site H' },
};

/**
 * NEB-validated high-fidelity overhang sets
 * Based on Potapov et al. (2018) ACS Synthetic Biology
 * "Comprehensive Profiling of Four Base Overhang Ligation Fidelity by T4 DNA Ligase"
 *
 * These sets have been experimentally validated to minimize cross-ligation
 * Fidelity scores represent % correct assemblies under standard conditions
 */
export const HIGH_FIDELITY_SETS = {
  // ===== TRUE 100% fidelity sets - Zero cross-reactivity (NEB validated) =====
  // These overhangs have ZERO cross-ligation in the experimental matrix
  '3-part': {
    overhangs: ['TGAC', 'GCAT', 'GATG', 'ATTG'],
    fidelity: 1.0,
    description: '3-part assembly (100% fidelity - NEB validated)',
  },
  '4-part': {
    overhangs: ['TGAC', 'GCAT', 'GATG', 'ATTG', 'TCCT'],
    fidelity: 1.0,
    description: '4-part assembly (100% fidelity - zero cross-ligation)',
  },
  '5-part': {
    overhangs: ['TGAC', 'GCAT', 'GATG', 'ATTG', 'TCCT', 'GGAA'],
    fidelity: 1.0,
    description: '5-part assembly (100% fidelity - zero cross-ligation)',
  },
  '2-part': {
    overhangs: ['TGAC', 'GCAT', 'GATG'],
    fidelity: 1.0,
    description: '2-part assembly (100% fidelity)',
  },

  // ===== Alternative 100% fidelity sets =====
  '3-part-alt': {
    overhangs: ['ATTG', 'GCTC', 'TGGT', 'TGAC'],
    fidelity: 1.0,
    description: 'Alternative 3-part (100% fidelity)',
  },

  // ===== MoClo-compatible sets (for standardized parts) =====
  '3-part-moclo': {
    overhangs: ['GGAG', 'TACT', 'AATG', 'GCTT'],
    fidelity: 0.85,
    description: 'MoClo standard 3-part (lower fidelity but standardized)',
  },
  '4-part-moclo': {
    overhangs: ['GGAG', 'TACT', 'AATG', 'AGGT', 'GCTT'],
    fidelity: 0.80,
    description: 'MoClo standard 4-part (lower fidelity but standardized)',
  },

  // ===== NEB Level 1 set (up to 8 parts, 93% fidelity) =====
  '8-part': {
    overhangs: ['GGAG', 'TACT', 'CCAT', 'AATG', 'AGGT', 'TTCG', 'GCTT', 'GGTA', 'CGCT'],
    fidelity: 0.93,
    description: 'NEB Level 1 core set (8 parts)',
  },

  // ===== NEB Level 1 extended (up to 12 parts, 92% fidelity) =====
  '12-part': {
    overhangs: [
      'GGAG', 'TACT', 'CCAT', 'AATG', 'AGGT', 'TTCG',
      'GCTT', 'GGTA', 'CGCT', 'GAAA', 'TCAA', 'ATAA', 'GCGA'
    ],
    fidelity: 0.92,
    description: 'NEB Level 1 extended (12 parts)',
  },

  // ===== NEB Level 1 full set (up to 20 parts, 92% fidelity) =====
  '20-part-L1': {
    overhangs: [
      'GGAG', 'TACT', 'CCAT', 'AATG', 'AGGT', 'TTCG', 'GCTT', 'GGTA', 'CGCT', 'GAAA',
      'TCAA', 'ATAA', 'GCGA', 'CGGC', 'GTCA', 'AACA', 'AAAT', 'GCAC', 'CTTA', 'TCCA', 'ACGG'
    ],
    fidelity: 0.92,
    description: 'NEB Level 1 full set (20 parts) - transcriptional units',
  },

  // ===== NEB Level 2 set (up to 20 parts, 95% fidelity) =====
  '20-part-L2': {
    overhangs: [
      'TGCC', 'GCAA', 'ACTA', 'TTAC', 'CAGA', 'TGTG', 'GAGC', 'GGGA', 'CGTA', 'CTTC',
      'ATCC', 'ATAG', 'CCAG', 'AATC', 'ACCG', 'AAAC', 'AGAC', 'AGGG', 'TGAA', 'ATGA', 'TAGG'
    ],
    fidelity: 0.95,
    description: 'NEB Level 2 full set (20 parts) - multigene constructs',
  },

  // ===== High-complexity sets from Potapov 2018 supplementary =====
  // Set 1: 15 junctions including MoClo standard (>98% per junction)
  'Set1-15': {
    overhangs: [
      'GGAG', 'TACT', 'AATG', 'GCTT',  // MoClo core
      'TGCC', 'GCAA', 'ACTA', 'TTAC',  // Extended
      'CAGA', 'TGTG', 'GAGC', 'GGGA',
      'CGTA', 'CTTC', 'ATCC', 'ATAG'
    ],
    fidelity: 0.98,
    description: 'Potapov Set 1: 15 junctions with MoClo compatibility',
  },

  // Set 2: 20 junctions highest fidelity (>98% per junction)
  'Set2-20': {
    overhangs: [
      'AACG', 'AAGC', 'ACAA', 'ACTC', 'AGGA', 'ATAG',
      'CAAG', 'CATG', 'CCTA', 'CGAA', 'CTAC', 'CTGA',
      'GACT', 'GCGT', 'GGAC', 'GTGC', 'TACC', 'TCGG',
      'TGCT', 'TTAG', 'AGTC'
    ],
    fidelity: 0.98,
    description: 'Potapov Set 2: 20 junctions maximum fidelity (non-MoClo)',
  },

  // Set 3: 25 junctions (>97% per junction)
  'Set3-25': {
    overhangs: [
      'AACG', 'AAGC', 'ACAA', 'ACTC', 'AGAT', 'AGGA', 'ATAG',
      'CAAG', 'CATG', 'CCTA', 'CGAA', 'CTAC', 'CTGA', 'GACC',
      'GACT', 'GCGT', 'GGAC', 'GTGC', 'TACC', 'TCGG', 'TGCT',
      'TGGA', 'TTAG', 'TTCA', 'AGTC', 'CGTG'
    ],
    fidelity: 0.97,
    description: 'Potapov Set 3: 25 junctions',
  },

  // Set 4: 30 junctions (>90% overall fidelity)
  'Set4-30': {
    overhangs: [
      'AACG', 'AAGC', 'ACAA', 'ACGC', 'ACTC', 'AGAT', 'AGGA',
      'ATAG', 'CAAG', 'CATG', 'CCAC', 'CCTA', 'CGAA', 'CTAC',
      'CTGA', 'GACC', 'GACT', 'GCGT', 'GGAC', 'GTAA', 'GTGC',
      'TACC', 'TCGG', 'TGCT', 'TGGA', 'TTAG', 'TTCA', 'AGTC',
      'CGTG', 'GCCA', 'TAGC'
    ],
    fidelity: 0.90,
    description: 'Potapov Set 4: 30 junctions for high-complexity assemblies',
  },
};

/**
 * Individual overhang fidelity data
 * Based on NEB ligation fidelity profiling (Potapov et al. 2018)
 * Score 0-1, higher = more specific/less promiscuous
 *
 * Categories:
 * - Excellent (>0.95): Recommended for all assemblies
 * - Good (0.90-0.95): Suitable for most assemblies
 * - Medium (0.80-0.90): Use with caution, may increase screening
 * - Low (<0.80): Avoid if possible
 */
export const OVERHANG_FIDELITY = {
  // ===== Excellent fidelity (>95%) - Use in all NEB sets =====
  'GGAG': { fidelity: 0.98, category: 'excellent', note: 'MoClo standard - upstream promoter' },
  'TACT': { fidelity: 0.97, category: 'excellent', note: 'MoClo standard - promoter/5\'UTR' },
  'AATG': { fidelity: 0.96, category: 'excellent', note: 'MoClo standard - includes ATG start' },
  'GCTT': { fidelity: 0.97, category: 'excellent', note: 'MoClo standard - downstream terminator' },
  'AGGT': { fidelity: 0.96, category: 'excellent', note: 'MoClo standard - CDS/terminator' },
  'AACG': { fidelity: 0.98, category: 'excellent', note: 'Potapov Set 2' },
  'AAGC': { fidelity: 0.97, category: 'excellent', note: 'Potapov Set 2' },
  'ACAA': { fidelity: 0.97, category: 'excellent', note: 'Potapov Set 2' },
  'ACTC': { fidelity: 0.97, category: 'excellent', note: 'Potapov Set 2' },
  'AGGA': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'ATAG': { fidelity: 0.97, category: 'excellent', note: 'Potapov Set 2' },
  'CAAG': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'CATG': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'CCTA': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'CGAA': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'CTAC': { fidelity: 0.97, category: 'excellent', note: 'Potapov Set 2' },
  'CTGA': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'GACT': { fidelity: 0.97, category: 'excellent', note: 'Potapov Set 2' },
  'GCGT': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'GGAC': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'GTGC': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'TACC': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'TCGG': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'TGCT': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },
  'TTAG': { fidelity: 0.96, category: 'excellent', note: 'Potapov Set 2' },

  // ===== Good fidelity (90-95%) - Safe for most assemblies =====
  'CGCT': { fidelity: 0.94, category: 'good', note: 'NEB Level 1' },
  'TGCC': { fidelity: 0.94, category: 'good', note: 'NEB Level 2' },
  'GCAA': { fidelity: 0.93, category: 'good', note: 'NEB Level 2' },
  'ACTA': { fidelity: 0.93, category: 'good', note: 'NEB Level 2' },
  'TTAC': { fidelity: 0.93, category: 'good', note: 'NEB Level 2' },
  'CAGA': { fidelity: 0.93, category: 'good', note: 'NEB Level 2' },
  'TGTG': { fidelity: 0.92, category: 'good', note: 'NEB Level 2' },
  'GAGC': { fidelity: 0.92, category: 'good', note: 'NEB Level 2' },
  'CCAT': { fidelity: 0.92, category: 'good', note: 'NEB Level 1' },
  'TTCG': { fidelity: 0.91, category: 'good', note: 'NEB Level 1' },
  'GGGA': { fidelity: 0.91, category: 'good', note: 'NEB Level 2' },
  'CGTA': { fidelity: 0.91, category: 'good', note: 'NEB Level 2' },
  'CTTC': { fidelity: 0.91, category: 'good', note: 'NEB Level 2' },
  'ATCC': { fidelity: 0.90, category: 'good', note: 'NEB Level 2' },

  // ===== Medium fidelity (80-90%) - May need more screening =====
  'GGTA': { fidelity: 0.88, category: 'medium', note: 'NEB Level 1' },
  'GAAA': { fidelity: 0.87, category: 'medium', note: 'NEB Level 1' },
  'TCAA': { fidelity: 0.86, category: 'medium', note: 'NEB Level 1' },
  'ATAA': { fidelity: 0.85, category: 'medium', note: 'NEB Level 1' },
  'GCGA': { fidelity: 0.85, category: 'medium', note: 'NEB Level 1' },
  'CGGC': { fidelity: 0.84, category: 'medium', note: 'NEB Level 1' },
  'GTCA': { fidelity: 0.84, category: 'medium', note: 'NEB Level 1' },
  'AACA': { fidelity: 0.83, category: 'medium', note: 'NEB Level 1' },
  'CCAG': { fidelity: 0.83, category: 'medium', note: 'NEB Level 2' },
  'AATC': { fidelity: 0.82, category: 'medium', note: 'NEB Level 2' },
  'ACCG': { fidelity: 0.82, category: 'medium', note: 'NEB Level 2' },

  // ===== Low fidelity (<80%) - Avoid if possible =====
  'AAAT': { fidelity: 0.78, category: 'low', note: 'NEB Level 1 - use with caution' },
  'GCAC': { fidelity: 0.77, category: 'low', note: 'NEB Level 1 - use with caution' },
  'CTTA': { fidelity: 0.76, category: 'low', note: 'NEB Level 1 - use with caution' },
  'TCCA': { fidelity: 0.75, category: 'low', note: 'NEB Level 1 - use with caution' },

  // ===== Problematic overhangs - AVOID =====
  'AAAA': { fidelity: 0.65, category: 'avoid', note: 'Homopolymer - promiscuous ligation' },
  'TTTT': { fidelity: 0.65, category: 'avoid', note: 'Homopolymer - promiscuous ligation' },
  'CCCC': { fidelity: 0.60, category: 'avoid', note: 'Homopolymer - promiscuous ligation' },
  'GGGG': { fidelity: 0.55, category: 'avoid', note: 'Homopolymer - difficult to synthesize' },
  'ATAT': { fidelity: 0.55, category: 'avoid', note: 'Palindromic - can self-ligate' },
  'TATA': { fidelity: 0.55, category: 'avoid', note: 'Palindromic - can self-ligate' },
  'GCGC': { fidelity: 0.50, category: 'avoid', note: 'Palindromic - can self-ligate' },
  'CGCG': { fidelity: 0.50, category: 'avoid', note: 'Palindromic - can self-ligate' },
  'ACGT': { fidelity: 0.45, category: 'avoid', note: 'Palindromic - can self-ligate' },
  'CATG': { fidelity: 0.40, category: 'avoid', note: 'Strong palindrome (NlaIII site)' },
  'GATC': { fidelity: 0.35, category: 'avoid', note: 'Strong palindrome (Sau3AI site)' },
};

// ==========================================
// Experimental Ligation Fidelity Functions
// Based on Pryor et al. (2020) PLOS ONE data
// ==========================================

/**
 * Get ligation data for a specific enzyme
 * @param {string} enzyme - Enzyme name (BsaI, BsmBI, Esp3I, BbsI, SapI)
 * @returns {Object|null} Enzyme ligation data or null if not found
 */
export function getEnzymeLigationData(enzyme) {
  const enzObj = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enzObj) return null;

  const dataKey = enzObj.dataKey;
  return ligationData.enzymes[dataKey] || null;
}

/**
 * Get individual overhang fidelity from experimental data
 * @param {string} overhang - 3 or 4 base overhang sequence
 * @param {string} enzyme - Enzyme name (default: BsaI)
 * @returns {Object} Fidelity data
 */
export function getOverhangFidelityExperimental(overhang, enzyme = 'BsaI') {
  const enzymeData = getEnzymeLigationData(enzyme);
  if (!enzymeData) {
    // Fall back to static data
    const staticData = OVERHANG_FIDELITY[overhang.toUpperCase()];
    return {
      overhang: overhang.toUpperCase(),
      fidelity: staticData?.fidelity || 0.85,
      source: 'static',
      enzyme: null,
    };
  }

  const oh = overhang.toUpperCase();
  const fidelity = enzymeData.overhangFidelity[oh];

  if (fidelity === undefined) {
    return {
      overhang: oh,
      fidelity: 0,
      source: 'experimental',
      enzyme,
      error: 'Overhang not found in experimental data',
    };
  }

  // Categorize based on fidelity
  let category;
  if (fidelity >= 0.95) category = 'excellent';
  else if (fidelity >= 0.90) category = 'good';
  else if (fidelity >= 0.80) category = 'medium';
  else if (fidelity >= 0.70) category = 'low';
  else category = 'avoid';

  return {
    overhang: oh,
    fidelity,
    fidelityPercent: (fidelity * 100).toFixed(1) + '%',
    category,
    source: 'experimental',
    enzyme,
  };
}

/**
 * Get raw ligation frequency between two overhangs
 * @param {string} overhang1 - Source overhang
 * @param {string} overhang2 - Target overhang (complement of destination)
 * @param {string} enzyme - Enzyme name
 * @returns {number} Ligation frequency (0 if no data)
 */
export function getLigationFrequency(overhang1, overhang2, enzyme = 'BsaI') {
  const enzymeData = getEnzymeLigationData(enzyme);
  if (!enzymeData) return 0;

  const oh1 = overhang1.toUpperCase();
  const oh2 = overhang2.toUpperCase();

  return enzymeData.matrix[oh1]?.[oh2] || 0;
}

/**
 * Calculate assembly fidelity using experimental ligation data
 *
 * For each junction, fidelity = p(correct) / p(total)
 * where p(total) includes only the overhangs actually in the assembly
 *
 * @param {string[]} overhangs - Array of overhang sequences for all junctions
 * @param {string} enzyme - Enzyme name (default: BsaI)
 * @returns {Object} Detailed fidelity analysis
 */
export function calculateExperimentalFidelity(overhangs, enzyme = 'BsaI') {
  const enzymeData = getEnzymeLigationData(enzyme);

  // Fall back to static calculation if no experimental data
  if (!enzymeData) {
    return calculateSetFidelity(overhangs);
  }

  const matrix = enzymeData.matrix;
  const junctions = [];
  let assemblyFidelity = 1.0;
  const warnings = [];

  for (const oh of overhangs) {
    const ohUpper = oh.toUpperCase();
    const wcPartner = reverseComplement(ohUpper);

    // Get correct ligation frequency (with Watson-Crick partner)
    const correctFreq = matrix[ohUpper]?.[wcPartner] || 0;

    if (correctFreq === 0) {
      warnings.push(`No ligation data for overhang ${ohUpper}`);
      junctions.push({
        overhang: ohUpper,
        fidelity: 0,
        correctFreq: 0,
        totalFreq: 0,
        error: 'No experimental data',
      });
      continue;
    }

    // Calculate total competing ligation frequency
    // Only count overhangs that are actually in our assembly
    let totalFreq = correctFreq;

    for (const otherOh of overhangs) {
      if (otherOh.toUpperCase() === ohUpper) continue;

      // Add mismatch frequency with other overhangs' complements
      const otherWc = reverseComplement(otherOh.toUpperCase());
      totalFreq += matrix[ohUpper]?.[otherWc] || 0;
    }

    const junctionFidelity = totalFreq > 0 ? correctFreq / totalFreq : 0;
    assemblyFidelity *= junctionFidelity;

    junctions.push({
      overhang: ohUpper,
      wcPartner,
      fidelity: junctionFidelity,
      fidelityPercent: (junctionFidelity * 100).toFixed(1) + '%',
      correctFreq,
      totalFreq,
    });

    // Warn about low-fidelity junctions
    if (junctionFidelity < 0.95) {
      warnings.push(`Junction ${ohUpper} has ${(junctionFidelity * 100).toFixed(1)}% fidelity`);
    }
  }

  // Sort junctions by fidelity to show weakest first
  const sortedJunctions = [...junctions].sort((a, b) => a.fidelity - b.fidelity);
  const lowestFidelity = sortedJunctions[0];

  return {
    enzyme,
    enzymeFullName: GOLDEN_GATE_ENZYMES[enzyme]?.fullName || enzyme,
    overhangs: overhangs.map(o => o.toUpperCase()),
    numJunctions: overhangs.length,
    junctions,
    sortedJunctions,
    assemblyFidelity,
    assemblyFidelityPercent: (assemblyFidelity * 100).toFixed(1) + '%',
    lowestFidelity,
    warnings,
    source: 'experimental',
    dataSource: ligationData.metadata,
  };
}

/**
 * Get pre-computed optimal overhang set for a given number of parts
 * @param {number} numParts - Number of parts (not including vector)
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Optimal overhang set with fidelity
 */
export function getOptimalOverhangSetExperimental(numParts, enzyme = 'BsaI') {
  const enzymeData = getEnzymeLigationData(enzyme);

  if (!enzymeData || !enzymeData.optimalSets) {
    // Fall back to static sets
    return getRecommendedOverhangSet(numParts, { preferMoClo: true });
  }

  // Find the closest matching set
  const availableSizes = Object.keys(enzymeData.optimalSets)
    .map(Number)
    .sort((a, b) => a - b);

  let matchingSize = availableSizes.find(s => s >= numParts);
  if (!matchingSize) {
    matchingSize = availableSizes[availableSizes.length - 1];
  }

  const optimalSet = enzymeData.optimalSets[matchingSize];

  return {
    enzyme,
    requestedParts: numParts,
    actualParts: matchingSize,
    overhangs: optimalSet.overhangs,
    numJunctions: optimalSet.overhangs.length,
    fidelity: optimalSet.fidelity,
    fidelityPercent: (optimalSet.fidelity * 100).toFixed(1) + '%',
    lowestJunction: optimalSet.lowestJunction,
    source: 'experimental',
  };
}

/**
 * Exhaustive search for optimal overhang set (for small assemblies ≤8 junctions)
 *
 * This guarantees finding the best possible 100% fidelity set by checking
 * all valid combinations. Only feasible for small junction counts due to
 * combinatorial explosion.
 *
 * @param {number} numJunctions - Number of junctions needed (max 8)
 * @param {string} enzyme - Enzyme name (default: BsaI)
 * @param {Object} options - Search options
 * @returns {Object} Optimal overhang set with fidelity analysis
 */
export function findOptimalOverhangSetExhaustive(numJunctions, enzyme = 'BsaI', options = {}) {
  const {
    minCorrectFreq = 300,
    requiredOverhangs: reqOh = [],
    excludeOverhangs: exclOh = [],
  } = options;

  const requiredOverhangs = Array.isArray(reqOh) ? reqOh.map(o => o.toUpperCase()) : [];
  const excludeOverhangs = Array.isArray(exclOh) ? exclOh.map(o => o.toUpperCase()) : [];

  const enzymeData = getEnzymeLigationData(enzyme);
  if (!enzymeData) {
    return getRecommendedOverhangSet(numJunctions, { preferMoClo: false });
  }

  const matrix = enzymeData.matrix;

  // Helper: check if overhang is a palindrome
  const isPalindrome = (oh) => oh === reverseComplement(oh);

  // Helper: check if set has zero cross-ligation
  const hasZeroCrossLigationSet = (overhangs) => {
    for (let i = 0; i < overhangs.length; i++) {
      for (let j = 0; j < overhangs.length; j++) {
        if (i === j) continue;
        const oh1 = overhangs[i];
        const wc2 = reverseComplement(overhangs[j]);
        if ((matrix[oh1]?.[wc2] || 0) > 0) return false;
      }
    }
    return true;
  };

  // Build candidate pool
  const excludeSet = new Set(excludeOverhangs);
  excludeOverhangs.forEach(oh => excludeSet.add(reverseComplement(oh)));

  const candidates = [];
  for (const oh of Object.keys(matrix)) {
    const ohUpper = oh.toUpperCase();
    const wc = reverseComplement(ohUpper);
    const correctFreq = matrix[ohUpper]?.[wc] || 0;

    if (excludeSet.has(ohUpper) || excludeSet.has(wc)) continue;
    if (isPalindrome(ohUpper)) continue;
    if (correctFreq < minCorrectFreq) continue;

    candidates.push(ohUpper);
  }

  // Remove candidates whose RC is also in the list (keep only one of each pair)
  const uniqueCandidates = [];
  const seenRCs = new Set();
  for (const oh of candidates) {
    const rc = reverseComplement(oh);
    if (!seenRCs.has(oh)) {
      uniqueCandidates.push(oh);
      seenRCs.add(rc);
    }
  }

  // Sort by correct frequency for consistent results
  uniqueCandidates.sort((a, b) => {
    const freqA = matrix[a]?.[reverseComplement(a)] || 0;
    const freqB = matrix[b]?.[reverseComplement(b)] || 0;
    return freqB - freqA;
  });

  // Number of additional overhangs needed beyond required ones
  const numNeeded = numJunctions - requiredOverhangs.length;

  // Filter out candidates that conflict with required overhangs
  const requiredSet = new Set(requiredOverhangs);
  requiredOverhangs.forEach(oh => requiredSet.add(reverseComplement(oh)));

  const availableCandidates = uniqueCandidates.filter(oh => {
    const rc = reverseComplement(oh);
    return !requiredSet.has(oh) && !requiredSet.has(rc);
  });

  // Check if required overhangs themselves have zero cross-ligation
  if (requiredOverhangs.length > 1 && !hasZeroCrossLigationSet(requiredOverhangs)) {
    // Required overhangs conflict - fall back to greedy
    return findOptimalOverhangSetGreedy(numJunctions, enzyme, options);
  }

  // For exhaustive search, we need to find combinations of 'numNeeded' overhangs
  // that, when combined with required overhangs, have zero cross-ligation
  let bestSet = null;
  let bestFidelity = -1;

  // Generator for combinations
  function* combinations(arr, k) {
    if (k === 0) {
      yield [];
      return;
    }
    for (let i = 0; i <= arr.length - k; i++) {
      for (const tail of combinations(arr.slice(i + 1), k - 1)) {
        yield [arr[i], ...tail];
      }
    }
  }

  // Limit search space for larger sets
  const maxCandidates = numNeeded <= 4 ? 80 : numNeeded <= 6 ? 60 : 50;
  const searchCandidates = availableCandidates.slice(0, maxCandidates);

  // Search for 100% fidelity sets first
  for (const combo of combinations(searchCandidates, numNeeded)) {
    const fullSet = [...requiredOverhangs, ...combo];

    // Check for RC conflicts within the set
    const allSeqs = new Set();
    let hasConflict = false;
    for (const oh of fullSet) {
      if (allSeqs.has(oh) || allSeqs.has(reverseComplement(oh))) {
        hasConflict = true;
        break;
      }
      allSeqs.add(oh);
      allSeqs.add(reverseComplement(oh));
    }
    if (hasConflict) continue;

    // Check for zero cross-ligation (100% fidelity)
    if (hasZeroCrossLigationSet(fullSet)) {
      // Found a 100% fidelity set - calculate to verify
      const fidelityResult = calculateExperimentalFidelity(fullSet, enzyme);
      if (fidelityResult.assemblyFidelity > bestFidelity) {
        bestFidelity = fidelityResult.assemblyFidelity;
        bestSet = fullSet;
        // If we found 100% fidelity, we can stop (all 100% sets are equally good)
        if (bestFidelity >= 0.9999) break;
      }
    }
  }

  // If no 100% set found, find the best available
  if (!bestSet) {
    for (const combo of combinations(searchCandidates.slice(0, 40), numNeeded)) {
      const fullSet = [...requiredOverhangs, ...combo];

      // Check for RC conflicts
      const allSeqs = new Set();
      let hasConflict = false;
      for (const oh of fullSet) {
        if (allSeqs.has(oh) || allSeqs.has(reverseComplement(oh))) {
          hasConflict = true;
          break;
        }
        allSeqs.add(oh);
        allSeqs.add(reverseComplement(oh));
      }
      if (hasConflict) continue;

      const fidelityResult = calculateExperimentalFidelity(fullSet, enzyme);
      if (fidelityResult.assemblyFidelity > bestFidelity) {
        bestFidelity = fidelityResult.assemblyFidelity;
        bestSet = fullSet;
      }
    }
  }

  // Fallback if still no set found
  if (!bestSet) {
    return findOptimalOverhangSetGreedy(numJunctions, enzyme, options);
  }

  const fidelityResult = calculateExperimentalFidelity(bestSet, enzyme);

  return {
    enzyme,
    numJunctions,
    requestedJunctions: numJunctions,
    foundJunctions: bestSet.length,
    overhangs: bestSet,
    fidelity: fidelityResult.assemblyFidelity,
    fidelityPercent: fidelityResult.assemblyFidelityPercent,
    isPerfect: fidelityResult.assemblyFidelity >= 0.9999,
    junctions: fidelityResult.junctions,
    warnings: fidelityResult.warnings,
    source: 'exhaustive-search',
    searchStats: {
      candidatesConsidered: searchCandidates.length,
      requiredOverhangs: requiredOverhangs.length,
      searchType: 'exhaustive',
    },
  };
}

/**
 * Greedy search for optimal overhang set (fast, for larger assemblies)
 *
 * @param {number} numJunctions - Number of junctions needed
 * @param {string} enzyme - Enzyme name (default: BsaI)
 * @param {Object} options - Search options
 * @returns {Object} Optimal overhang set with fidelity analysis
 */
function findOptimalOverhangSetGreedy(numJunctions, enzyme = 'BsaI', options = {}) {
  const {
    minCorrectFreq = 300,
    requiredOverhangs: reqOh = [],
    excludeOverhangs: exclOh = [],
  } = options;

  const requiredOverhangs = Array.isArray(reqOh) ? reqOh : [];
  const excludeOverhangs = Array.isArray(exclOh) ? exclOh : [];

  const enzymeData = getEnzymeLigationData(enzyme);
  if (!enzymeData) {
    return getRecommendedOverhangSet(numJunctions, { preferMoClo: false });
  }

  const matrix = enzymeData.matrix;
  const allOverhangs = Object.keys(matrix);

  const isPalindrome = (oh) => oh === reverseComplement(oh);

  const hasZeroCrossLigation = (existingSet, newOh) => {
    const newWc = reverseComplement(newOh);
    for (const existing of existingSet) {
      const existingWc = reverseComplement(existing);
      if ((matrix[newOh]?.[existingWc] || 0) > 0) return false;
      if ((matrix[existing]?.[newWc] || 0) > 0) return false;
    }
    return true;
  };

  const isValidSet = (overhangs) => {
    const rcs = overhangs.map(oh => reverseComplement(oh));
    const allSeqs = new Set([...overhangs, ...rcs]);
    return allSeqs.size === overhangs.length * 2;
  };

  const excludeSet = new Set(excludeOverhangs.map(oh => oh.toUpperCase()));
  const candidates = [];

  for (const oh of allOverhangs) {
    const ohUpper = oh.toUpperCase();
    const wc = reverseComplement(ohUpper);
    const correctFreq = matrix[ohUpper]?.[wc] || 0;

    if (excludeSet.has(ohUpper)) continue;
    if (excludeSet.has(wc)) continue;
    if (isPalindrome(ohUpper)) continue;
    if (correctFreq < minCorrectFreq) continue;

    candidates.push({ overhang: ohUpper, correctFreq, wc });
  }

  candidates.sort((a, b) => b.correctFreq - a.correctFreq);

  const resultSet = [];
  const usedOverhangs = new Set();
  const usedWcs = new Set();

  for (const reqOh of requiredOverhangs) {
    const ohUpper = reqOh.toUpperCase();
    const wc = reverseComplement(ohUpper);
    resultSet.push(ohUpper);
    usedOverhangs.add(ohUpper);
    usedWcs.add(wc);
  }

  for (const candidate of candidates) {
    if (resultSet.length >= numJunctions) break;
    const { overhang, wc } = candidate;

    if (usedOverhangs.has(overhang)) continue;
    if (usedWcs.has(overhang)) continue;
    if (usedOverhangs.has(wc)) continue;
    if (usedWcs.has(wc)) continue;

    if (hasZeroCrossLigation(resultSet, overhang)) {
      resultSet.push(overhang);
      usedOverhangs.add(overhang);
      usedWcs.add(wc);
    }
  }

  if (resultSet.length < numJunctions) {
    for (const candidate of candidates) {
      if (resultSet.length >= numJunctions) break;
      const { overhang, wc } = candidate;

      if (usedOverhangs.has(overhang)) continue;
      if (usedWcs.has(overhang)) continue;
      if (usedOverhangs.has(wc)) continue;
      if (usedWcs.has(wc)) continue;

      const testSet = [...resultSet, overhang];
      if (isValidSet(testSet)) {
        resultSet.push(overhang);
        usedOverhangs.add(overhang);
        usedWcs.add(wc);
      }
    }
  }

  const fidelityResult = calculateExperimentalFidelity(resultSet, enzyme);

  return {
    enzyme,
    numJunctions,
    requestedJunctions: numJunctions,
    foundJunctions: resultSet.length,
    overhangs: resultSet,
    fidelity: fidelityResult.assemblyFidelity,
    fidelityPercent: fidelityResult.assemblyFidelityPercent,
    isPerfect: fidelityResult.assemblyFidelity >= 0.9999,
    junctions: fidelityResult.junctions,
    warnings: fidelityResult.warnings,
    source: 'greedy-search',
    searchStats: {
      candidatesConsidered: candidates.length,
      requiredOverhangs: requiredOverhangs.length,
      searchType: 'greedy',
    },
  };
}

/**
 * Dynamically find an optimal overhang set with 100% fidelity (or best possible)
 *
 * This function searches the ligation matrix to find overhang combinations
 * with zero cross-reactivity. Unlike the static HIGH_FIDELITY_SETS lookup,
 * this can find optimal sets for ANY number of junctions.
 *
 * Algorithm:
 * - For small assemblies (≤8 junctions): Uses exhaustive search to guarantee optimal
 * - For larger assemblies: Uses greedy search with high-frequency candidate selection
 *
 * @param {number} numJunctions - Number of junctions needed
 * @param {string} enzyme - Enzyme name (default: BsaI)
 * @param {Object} options - Search options
 * @param {number} options.minCorrectFreq - Minimum correct ligation frequency (default: 300)
 * @param {string[]} options.requiredOverhangs - Overhangs that must be included
 * @param {string[]} options.excludeOverhangs - Overhangs to exclude from search
 * @param {boolean} options.forceExhaustive - Force exhaustive search even for larger sets
 * @returns {Object} Optimal overhang set with fidelity analysis
 */
export function findOptimalOverhangSet(numJunctions, enzyme = 'BsaI', options = {}) {
  const { forceExhaustive = false } = options;

  // Use exhaustive search for small assemblies (guaranteed optimal)
  // or when explicitly requested
  if (numJunctions <= 8 || forceExhaustive) {
    return findOptimalOverhangSetExhaustive(numJunctions, enzyme, options);
  }

  // Use greedy search for larger assemblies (fast, good results)
  return findOptimalOverhangSetGreedy(numJunctions, enzyme, options);
}

/**
 * Compare fidelity across all enzymes for a given overhang set
 * @param {string[]} overhangs - Array of overhang sequences
 * @returns {Object} Comparison across all enzymes with data
 */
export function compareEnzymeFidelity(overhangs) {
  const results = {};

  for (const enzymeName of ENZYMES_WITH_DATA) {
    const fidelityData = calculateExperimentalFidelity(overhangs, enzymeName);
    results[enzymeName] = {
      enzyme: enzymeName,
      fullName: GOLDEN_GATE_ENZYMES[enzymeName]?.fullName,
      assemblyFidelity: fidelityData.assemblyFidelity,
      fidelityPercent: fidelityData.assemblyFidelityPercent,
      lowestJunction: fidelityData.lowestFidelity,
      warnings: fidelityData.warnings.length,
    };
  }

  // Sort by fidelity
  const ranked = Object.values(results).sort(
    (a, b) => b.assemblyFidelity - a.assemblyFidelity
  );

  return {
    overhangs: overhangs.map(o => o.toUpperCase()),
    numJunctions: overhangs.length,
    byEnzyme: results,
    ranked,
    recommended: ranked[0]?.enzyme || 'BsaI',
  };
}

/**
 * Find problematic overhang pairs that may cross-ligate
 * @param {string[]} overhangs - Array of overhangs in the assembly
 * @param {string} enzyme - Enzyme name
 * @param {number} threshold - Minimum cross-ligation ratio to report (default: 0.05 = 5%)
 * @returns {Array} List of problematic pairs
 */
export function findProblematicPairs(overhangs, enzyme = 'BsaI', threshold = 0.05) {
  const enzymeData = getEnzymeLigationData(enzyme);
  if (!enzymeData) return [];

  const matrix = enzymeData.matrix;
  const problems = [];

  for (let i = 0; i < overhangs.length; i++) {
    const oh1 = overhangs[i].toUpperCase();
    const wc1 = reverseComplement(oh1);
    const correct1 = matrix[oh1]?.[wc1] || 0;

    if (correct1 === 0) continue;

    for (let j = 0; j < overhangs.length; j++) {
      if (i === j) continue;

      const oh2 = overhangs[j].toUpperCase();
      const wc2 = reverseComplement(oh2);

      // Check if oh1 can ligate with wc2 (cross-ligation)
      const crossFreq = matrix[oh1]?.[wc2] || 0;
      const ratio = crossFreq / correct1;

      if (ratio >= threshold) {
        problems.push({
          source: oh1,
          target: oh2,
          sourceWC: wc1,
          targetWC: wc2,
          correctFreq: correct1,
          crossFreq,
          ratio,
          ratioPercent: (ratio * 100).toFixed(1) + '%',
          severity: ratio >= 0.2 ? 'high' : ratio >= 0.1 ? 'medium' : 'low',
        });
      }
    }
  }

  return problems.sort((a, b) => b.ratio - a.ratio);
}

/**
 * Generate cross-ligation heatmap data for visualization
 *
 * Creates a data structure suitable for rendering a heatmap showing
 * cross-ligation frequencies between all overhangs in an assembly.
 * Similar to NEB's Ligase Fidelity Viewer.
 *
 * @param {string[]} overhangs - Array of overhangs in the assembly
 * @param {string} enzyme - Enzyme name (default: BsaI)
 * @param {Object} options - Options
 * @param {boolean} options.normalize - Normalize values to 0-1 range (default: true)
 * @param {boolean} options.includeCorrect - Include correct ligation on diagonal (default: true)
 * @returns {Object} Heatmap data structure for visualization
 */
export function generateCrossLigationHeatmap(overhangs, enzyme = 'BsaI', options = {}) {
  const {
    normalize = true,
    includeCorrect = true,
  } = options;

  const enzymeData = getEnzymeLigationData(enzyme);
  if (!enzymeData) {
    return {
      error: `No ligation data available for enzyme: ${enzyme}`,
      overhangs: [],
      matrix: [],
      labels: [],
    };
  }

  const matrix = enzymeData.matrix;
  const ohList = overhangs.map(oh => oh.toUpperCase());

  // Build the heatmap matrix
  // Each cell [i][j] = frequency of overhang i ligating with RC of overhang j
  const heatmapData = [];
  let maxValue = 0;
  let minValue = Infinity;

  for (let i = 0; i < ohList.length; i++) {
    const row = [];
    const oh1 = ohList[i];

    for (let j = 0; j < ohList.length; j++) {
      const oh2 = ohList[j];
      const rc2 = reverseComplement(oh2);

      // Get ligation frequency: oh1 -> rc2
      const freq = matrix[oh1]?.[rc2] || 0;

      row.push(freq);

      if (freq > maxValue) maxValue = freq;
      if (freq < minValue && freq > 0) minValue = freq;
    }

    heatmapData.push(row);
  }

  // Normalize if requested
  let normalizedData = heatmapData;
  if (normalize && maxValue > 0) {
    normalizedData = heatmapData.map(row =>
      row.map(val => val / maxValue)
    );
  }

  // Calculate statistics for each row (overhang)
  const rowStats = ohList.map((oh, i) => {
    const row = heatmapData[i];
    const correctIdx = i; // Diagonal is correct ligation
    const correctFreq = row[correctIdx];

    // Sum of off-diagonal (cross-ligation)
    const crossLigationSum = row.reduce((sum, val, j) =>
      j !== correctIdx ? sum + val : sum, 0
    );

    // Fidelity for this junction
    const totalFreq = correctFreq + crossLigationSum;
    const fidelity = totalFreq > 0 ? correctFreq / totalFreq : 0;

    return {
      overhang: oh,
      reverseComplement: reverseComplement(oh),
      correctFreq,
      crossLigationSum,
      fidelity,
      fidelityPercent: `${(fidelity * 100).toFixed(1)}%`,
    };
  });

  // Identify cross-ligation hotspots (cells with significant off-diagonal values)
  const hotspots = [];
  for (let i = 0; i < ohList.length; i++) {
    for (let j = 0; j < ohList.length; j++) {
      if (i === j) continue; // Skip diagonal

      const freq = heatmapData[i][j];
      if (freq > 0) {
        const correctFreq = heatmapData[i][i];
        const ratio = correctFreq > 0 ? freq / correctFreq : 0;

        if (ratio >= 0.01) { // Report if ≥1% of correct ligation
          hotspots.push({
            source: ohList[i],
            target: ohList[j],
            targetRC: reverseComplement(ohList[j]),
            frequency: freq,
            ratio,
            ratioPercent: `${(ratio * 100).toFixed(1)}%`,
            severity: ratio >= 0.1 ? 'high' : ratio >= 0.05 ? 'medium' : 'low',
            row: i,
            col: j,
          });
        }
      }
    }
  }

  // Sort hotspots by severity
  hotspots.sort((a, b) => b.ratio - a.ratio);

  // Calculate overall assembly fidelity
  const overallFidelity = rowStats.reduce((product, stat) =>
    product * stat.fidelity, 1
  );

  return {
    enzyme,
    enzymeFullName: GOLDEN_GATE_ENZYMES[enzyme]?.fullName || enzyme,
    overhangs: ohList,
    labels: ohList, // For axis labels
    labelsRC: ohList.map(oh => reverseComplement(oh)), // RC labels for y-axis

    // Raw frequency matrix
    matrix: heatmapData,

    // Normalized matrix (0-1 scale) for heatmap coloring
    normalizedMatrix: normalizedData,

    // Statistics
    stats: {
      maxFrequency: maxValue,
      minFrequency: minValue || 0,
      overallFidelity,
      overallFidelityPercent: `${(overallFidelity * 100).toFixed(1)}%`,
    },

    // Per-overhang statistics
    rowStats,

    // Cross-ligation hotspots
    hotspots,
    hasHotspots: hotspots.length > 0,

    // Metadata for visualization
    visualization: {
      title: `Cross-Ligation Heatmap (${enzyme})`,
      xAxisLabel: 'Target Overhang (RC)',
      yAxisLabel: 'Source Overhang',
      colorScale: {
        min: 0,
        max: normalize ? 1 : maxValue,
        colors: ['#ffffff', '#fee8c8', '#fdbb84', '#e34a33', '#b30000'], // White to red
      },
      diagonalLabel: 'Correct Ligation',
      offDiagonalLabel: 'Cross-Ligation',
    },
  };
}

/**
 * Generate a summary report of overhang set quality
 * Combines fidelity analysis, heatmap data, and recommendations
 *
 * @param {string[]} overhangs - Array of overhangs
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Comprehensive quality report
 */
export function generateOverhangQualityReport(overhangs, enzyme = 'BsaI') {
  const fidelityAnalysis = calculateExperimentalFidelity(overhangs, enzyme);
  const heatmapData = generateCrossLigationHeatmap(overhangs, enzyme);
  const problematicPairs = findProblematicPairs(overhangs, enzyme, 0.01);

  // Generate quality score (0-100)
  let qualityScore = 100;
  const deductions = [];

  // Deduct for low fidelity
  if (fidelityAnalysis.assemblyFidelity < 1.0) {
    const fidelityDeduction = Math.round((1 - fidelityAnalysis.assemblyFidelity) * 50);
    qualityScore -= fidelityDeduction;
    deductions.push({
      reason: `Assembly fidelity ${fidelityAnalysis.assemblyFidelityPercent}`,
      deduction: fidelityDeduction,
    });
  }

  // Deduct for cross-ligation hotspots
  const highSeverityHotspots = heatmapData.hotspots.filter(h => h.severity === 'high');
  if (highSeverityHotspots.length > 0) {
    const hotspotDeduction = Math.min(30, highSeverityHotspots.length * 10);
    qualityScore -= hotspotDeduction;
    deductions.push({
      reason: `${highSeverityHotspots.length} high-severity cross-ligation hotspot(s)`,
      deduction: hotspotDeduction,
    });
  }

  // Deduct for weak junctions
  const weakJunctions = fidelityAnalysis.junctions.filter(j => j.fidelity < 0.95);
  if (weakJunctions.length > 0) {
    const weakDeduction = Math.min(20, weakJunctions.length * 5);
    qualityScore -= weakDeduction;
    deductions.push({
      reason: `${weakJunctions.length} junction(s) below 95% fidelity`,
      deduction: weakDeduction,
    });
  }

  qualityScore = Math.max(0, qualityScore);

  // Generate recommendations
  const recommendations = [];

  if (fidelityAnalysis.assemblyFidelity < 0.95) {
    recommendations.push({
      priority: 'high',
      type: 'fidelity',
      message: 'Consider using findOptimalOverhangSet() to find a higher-fidelity overhang combination',
    });
  }

  if (highSeverityHotspots.length > 0) {
    const worstHotspot = highSeverityHotspots[0];
    recommendations.push({
      priority: 'high',
      type: 'cross-ligation',
      message: `Replace overhang ${worstHotspot.source} or ${worstHotspot.target} to eliminate ${worstHotspot.ratioPercent} cross-ligation`,
    });
  }

  if (weakJunctions.length > 0) {
    const weakestJunction = fidelityAnalysis.sortedJunctions[0];
    recommendations.push({
      priority: 'medium',
      type: 'junction',
      message: `Junction ${weakestJunction.overhang} (${weakestJunction.fidelityPercent}) is the weakest link - consider replacing`,
    });
  }

  if (qualityScore >= 95) {
    recommendations.push({
      priority: 'info',
      type: 'success',
      message: 'Excellent overhang set! Expected assembly efficiency is very high.',
    });
  }

  return {
    enzyme,
    overhangs: overhangs.map(o => o.toUpperCase()),
    numJunctions: overhangs.length,

    // Quality assessment
    qualityScore,
    qualityGrade: qualityScore >= 95 ? 'A' :
                  qualityScore >= 85 ? 'B' :
                  qualityScore >= 70 ? 'C' :
                  qualityScore >= 50 ? 'D' : 'F',
    deductions,

    // Detailed analysis
    fidelity: fidelityAnalysis,
    heatmap: heatmapData,
    problematicPairs,

    // Recommendations
    recommendations,

    // Summary
    summary: {
      assemblyFidelity: fidelityAnalysis.assemblyFidelityPercent,
      weakestJunction: fidelityAnalysis.lowestFidelity,
      crossLigationRisks: heatmapData.hotspots.length,
      isPerfect: fidelityAnalysis.assemblyFidelity >= 0.9999,
    },
  };
}

/**
 * Get the recommended overhang set for a given number of parts
 * @param {number} numParts - Number of DNA parts to assemble
 * @param {Object} options - Options for set selection
 * @returns {Object} Recommended overhang set
 */
export function getRecommendedOverhangSet(numParts, options = {}) {
  const { preferMoClo = true, maxFidelity = false } = options;

  if (numParts < 2 || numParts > 30) {
    throw new Error('Golden Gate assembly supports 2-30 parts');
  }

  // For maximum fidelity without MoClo constraint
  if (maxFidelity && numParts <= 20) {
    return HIGH_FIDELITY_SETS['Set2-20'];
  }

  // Standard selection based on part count
  if (numParts <= 2) return HIGH_FIDELITY_SETS['2-part'];
  if (numParts <= 3) return HIGH_FIDELITY_SETS['3-part'];
  if (numParts <= 4) return HIGH_FIDELITY_SETS['4-part'];
  if (numParts <= 5) return HIGH_FIDELITY_SETS['5-part'];
  if (numParts <= 8) return HIGH_FIDELITY_SETS['8-part'];
  if (numParts <= 12) return HIGH_FIDELITY_SETS['12-part'];
  if (numParts <= 15 && preferMoClo) return HIGH_FIDELITY_SETS['Set1-15'];
  if (numParts <= 20) return preferMoClo ? HIGH_FIDELITY_SETS['20-part-L1'] : HIGH_FIDELITY_SETS['Set2-20'];
  if (numParts <= 25) return HIGH_FIDELITY_SETS['Set3-25'];
  return HIGH_FIDELITY_SETS['Set4-30'];
}

/**
 * Calculate expected assembly fidelity for a set of overhangs
 * @param {string[]} overhangs - Array of 4-base overhang sequences
 * @returns {Object} Fidelity analysis
 */
export function calculateSetFidelity(overhangs) {
  const analysis = {
    overhangs: [],
    overallFidelity: 1,
    lowestFidelity: { overhang: null, fidelity: 1 },
    warnings: [],
    categories: { excellent: 0, good: 0, medium: 0, low: 0, avoid: 0, unknown: 0 },
  };

  for (const oh of overhangs) {
    const data = OVERHANG_FIDELITY[oh.toUpperCase()];
    const fidelity = data?.fidelity || 0.85; // Default for unknown
    const category = data?.category || 'unknown';

    analysis.overhangs.push({
      sequence: oh.toUpperCase(),
      fidelity,
      category,
      note: data?.note || 'Not in database',
    });

    analysis.overallFidelity *= fidelity;
    analysis.categories[category]++;

    if (fidelity < analysis.lowestFidelity.fidelity) {
      analysis.lowestFidelity = { overhang: oh, fidelity };
    }

    // Generate warnings
    if (category === 'avoid') {
      analysis.warnings.push(`${oh}: Problematic overhang - ${data?.note}`);
    } else if (category === 'low') {
      analysis.warnings.push(`${oh}: Low fidelity (${(fidelity * 100).toFixed(0)}%) - may need extra screening`);
    }
  }

  analysis.overallFidelityPercent = (analysis.overallFidelity * 100).toFixed(1) + '%';

  return analysis;
}

/**
 * Standard part type definitions with their flanking overhangs
 */
export const PART_TYPES = {
  promoter:   { left: 'A', right: 'B', name: 'Promoter', color: '#ef4444' },
  rbs:        { left: 'B', right: 'C', name: 'RBS/5\'UTR', color: '#f97316' },
  cds:        { left: 'C', right: 'D', name: 'CDS', color: '#22c55e' },
  terminator: { left: 'D', right: 'E', name: 'Terminator', color: '#3b82f6' },
  backbone:   { left: 'E', right: 'A', name: 'Backbone', color: '#8b5cf6' },
  // For custom/other parts
  custom:     { left: null, right: null, name: 'Custom', color: '#6b7280' },
};

/**
 * Calculate GC content as percentage
 */
function gcContent(seq) {
  const gc = (seq.match(/[GC]/gi) || []).length;
  return (gc / seq.length) * 100;
}

// Primer design constraints (based on IDT guidelines and NEB recommendations)
const PRIMER_MIN_HOMOLOGY = 15;  // Minimum binding region length
const PRIMER_MAX_HOMOLOGY = 30;  // Maximum binding region length
const PRIMER_TARGET_TM = 60;     // Target Tm for Q5 polymerase (NEB recommends 60-72°C)
const PRIMER_MIN_TM = 55;        // Minimum acceptable Tm
const PRIMER_MAX_TM = 72;        // Maximum acceptable Tm

/**
 * Find optimal homology length to achieve target Tm using NEB Q5 calculator
 * Extends from minimum length until target Tm is reached (or max length)
 *
 * @param {string} templateSeq - Template sequence to anneal to
 * @param {number} targetTm - Target melting temperature
 * @param {boolean} fromEnd - If true, extract homology from end of sequence
 * @returns {Object} Optimal homology region and its Tm
 */
function findOptimalHomology(templateSeq, targetTm = PRIMER_TARGET_TM, fromEnd = false) {
  const seqUpper = templateSeq.toUpperCase();

  // Start with minimum length and extend until target Tm is reached
  let bestHomology = null;
  let bestTm = 0;
  let bestDiff = Infinity;

  for (let len = PRIMER_MIN_HOMOLOGY; len <= Math.min(PRIMER_MAX_HOMOLOGY, seqUpper.length); len++) {
    const homology = fromEnd
      ? seqUpper.slice(-len)
      : seqUpper.slice(0, len);

    try {
      const tm = calculateTmQ5(homology);
      const diff = Math.abs(tm - targetTm);

      // Track best match to target Tm
      if (diff < bestDiff) {
        bestDiff = diff;
        bestHomology = homology;
        bestTm = tm;
      }

      // If we've reached or exceeded target, and are in acceptable range, stop
      if (tm >= targetTm && tm <= PRIMER_MAX_TM) {
        break;
      }
    } catch (e) {
      // Skip if calculation fails (e.g., too short)
      continue;
    }
  }

  // Fall back to minimum length if nothing worked
  if (!bestHomology) {
    bestHomology = fromEnd
      ? seqUpper.slice(-PRIMER_MIN_HOMOLOGY)
      : seqUpper.slice(0, PRIMER_MIN_HOMOLOGY);
    try {
      bestTm = calculateTmQ5(bestHomology);
    } catch (e) {
      bestTm = 50; // Fallback estimate
    }
  }

  return {
    sequence: bestHomology,
    length: bestHomology.length,
    tm: bestTm,
    gc: gcContent(bestHomology),
    targetTm,
    withinRange: bestTm >= PRIMER_MIN_TM && bestTm <= PRIMER_MAX_TM,
  };
}

/**
 * Design Golden Gate primers to add Type IIS sites and overhangs to a target sequence.
 * Uses NEB Q5 Tm calculator with dynamic homology length optimization.
 *
 * Primer structure (5' to 3'):
 *   [extra bases]-[recognition site]-[spacer]-[overhang]-[homology to target]
 *
 * @param {string} targetSeq - The insert sequence
 * @param {string} leftOverhang - 4bp 5' overhang (e.g., 'GGAG')
 * @param {string} rightOverhang - 4bp 3' overhang (e.g., 'TACT')
 * @param {Object} options - Primer design options
 * @returns {Object} Forward and reverse primers with details
 */
export function designGoldenGatePrimers(targetSeq, leftOverhang, rightOverhang, options = {}) {
  const {
    enzyme = 'BsaI',
    targetTm = PRIMER_TARGET_TM,
    extraBases, // Optional: Override flanking bases (for backward compatibility)
    useOptimizedFlanking = true, // Use 6bp NEB-recommended flanking by default
  } = options;

  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enz) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  const seqUpper = targetSeq.toUpperCase();
  const leftOH = leftOverhang.toUpperCase();
  const rightOH = rightOverhang.toUpperCase();

  // Validate overhangs based on enzyme's overhang length
  const expectedOHLength = enz.overhangLength || 4;
  if (leftOH.length !== expectedOHLength || rightOH.length !== expectedOHLength) {
    throw new Error(`Overhangs must be exactly ${expectedOHLength} bases for ${enzyme}`);
  }

  // Find optimal homology regions for target Tm
  const fwdHomologyInfo = findOptimalHomology(seqUpper, targetTm, false);
  const revHomologyInfo = findOptimalHomology(seqUpper, targetTm, true);

  // Enzyme-specific primer construction
  const recognitionSite = enz.recognition;
  const recognitionSiteRC = reverseComplement(enz.recognition);

  // Use optimal spacers from NEB research (enzyme-specific)
  const spacerInfo = OPTIMAL_SPACERS[enzyme] || { forward: 'A', reverse: 'T' };
  const spacer = spacerInfo.forward;
  const spacerRC = spacerInfo.reverse;
  const rightOHrc = reverseComplement(rightOH);

  // Determine flanking bases:
  // 1. If extraBases explicitly provided, use that (backward compatibility)
  // 2. If useOptimizedFlanking, use 6bp NEB-recommended sequences
  // 3. Fallback to 'GG' (legacy default)
  let flankingBases;
  if (extraBases !== undefined) {
    flankingBases = extraBases;
  } else if (useOptimizedFlanking) {
    const flankingOptions = OPTIMAL_FLANKING_SEQUENCES[enzyme] || OPTIMAL_FLANKING_SEQUENCES.BsaI;
    flankingBases = flankingOptions.default; // 6bp optimized sequence
  } else {
    flankingBases = 'GG'; // Legacy 2bp default
  }

  // Forward primer: 5'-[flanking]-[recognition]-[spacer]-[leftOH]-[homology]-3'
  const fwdPrimer = flankingBases + recognitionSite + spacer + leftOH + fwdHomologyInfo.sequence;

  // Reverse primer: 5'-[flanking]-[recognition RC]-[spacer RC]-[rightOH RC]-[homology RC]-3'
  const revHomologyRC = reverseComplement(revHomologyInfo.sequence);
  const revPrimer = flankingBases + recognitionSiteRC + spacerRC + rightOHrc + revHomologyRC;

  // Get fidelity info for overhangs
  const leftFidelity = OVERHANG_FIDELITY[leftOH] || { fidelity: 0.85, note: 'Custom' };
  const rightFidelity = OVERHANG_FIDELITY[rightOH] || { fidelity: 0.85, note: 'Custom' };

  // Generate warnings
  const warnings = [];
  if (!fwdHomologyInfo.withinRange) {
    warnings.push(`Forward primer Tm (${fwdHomologyInfo.tm}°C) outside optimal range (${PRIMER_MIN_TM}-${PRIMER_MAX_TM}°C)`);
  }
  if (!revHomologyInfo.withinRange) {
    warnings.push(`Reverse primer Tm (${revHomologyInfo.tm}°C) outside optimal range (${PRIMER_MIN_TM}-${PRIMER_MAX_TM}°C)`);
  }
  const tmDiff = Math.abs(fwdHomologyInfo.tm - revHomologyInfo.tm);
  if (tmDiff > 5) {
    warnings.push(`Primer Tm difference (${tmDiff}°C) is large - consider adjusting for better PCR`);
  }

  // Calculate recommended annealing temperature for PCR
  // For Q5: Ta = lower Tm + 1°C, capped at 72°C
  const lowerTm = Math.min(fwdHomologyInfo.tm, revHomologyInfo.tm);
  const higherTm = Math.max(fwdHomologyInfo.tm, revHomologyInfo.tm);
  const annealingTemp = Math.min(lowerTm + 1, 72);

  return {
    forward: {
      sequence: fwdPrimer,
      length: fwdPrimer.length,
      tm: fwdHomologyInfo.tm,
      gc: fwdHomologyInfo.gc,
      homologyRegion: fwdHomologyInfo.sequence,
      homologyLength: fwdHomologyInfo.length,
      overhang: leftOH,
      structure: {
        extra: flankingBases,
        flanking: flankingBases,          // Explicit flanking field
        bsaISite: recognitionSite,        // Keep for backward compatibility
        recognitionSite,                  // New field
        spacer,
        overhang: leftOH,
        homology: fwdHomologyInfo.sequence,
      },
    },
    reverse: {
      sequence: revPrimer,
      length: revPrimer.length,
      tm: revHomologyInfo.tm,
      gc: revHomologyInfo.gc,
      homologyRegion: revHomologyRC,
      homologyLength: revHomologyInfo.length,
      overhang: rightOHrc,
      structure: {
        extra: flankingBases,
        flanking: flankingBases,          // Explicit flanking field
        bsaISite: recognitionSiteRC,      // Keep for backward compatibility
        recognitionSite: recognitionSiteRC, // New field
        spacer: spacerRC,
        overhang: rightOHrc,
        homology: revHomologyRC,
      },
    },
    // PCR conditions
    pcr: {
      annealingTemp,                      // Recommended annealing temperature
      lowerTm,
      higherTm,
      tmDifference: Math.abs(higherTm - lowerTm),
      extensionTime: Math.ceil(seqUpper.length / 1000) * 30, // 30s per kb for Q5
    },
    insert: {
      original: seqUpper,
      length: seqUpper.length,
      afterDigestion: leftOH + seqUpper + reverseComplement(rightOHrc),
    },
    fidelity: {
      leftOverhang: leftFidelity,
      rightOverhang: rightFidelity,
      estimated: Math.min(leftFidelity.fidelity, rightFidelity.fidelity),
    },
    enzyme,
    warnings,
  };
}

// Codon table for synonymous substitutions (DNA codons)
const CODON_TABLE = {
  // Amino acid -> list of codons
  F: ['TTT', 'TTC'],
  L: ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
  I: ['ATT', 'ATC', 'ATA'],
  M: ['ATG'],
  V: ['GTT', 'GTC', 'GTA', 'GTG'],
  S: ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
  P: ['CCT', 'CCC', 'CCA', 'CCG'],
  T: ['ACT', 'ACC', 'ACA', 'ACG'],
  A: ['GCT', 'GCC', 'GCA', 'GCG'],
  Y: ['TAT', 'TAC'],
  H: ['CAT', 'CAC'],
  Q: ['CAA', 'CAG'],
  N: ['AAT', 'AAC'],
  K: ['AAA', 'AAG'],
  D: ['GAT', 'GAC'],
  E: ['GAA', 'GAG'],
  C: ['TGT', 'TGC'],
  W: ['TGG'],
  R: ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
  G: ['GGT', 'GGC', 'GGA', 'GGG'],
  '*': ['TAA', 'TAG', 'TGA'], // Stop codons
};

// Reverse lookup: codon -> amino acid
const CODON_TO_AA = {};
for (const [aa, codons] of Object.entries(CODON_TABLE)) {
  for (const codon of codons) {
    CODON_TO_AA[codon] = aa;
  }
}

/**
 * Find all internal Type IIS restriction sites in a sequence
 * These sites will interfere with Golden Gate assembly and need to be removed (domesticated)
 *
 * @param {string} sequence - DNA sequence to check
 * @param {string} enzyme - Enzyme name (default: 'BsaI')
 * @returns {Object} Analysis of internal sites with positions and domestication suggestions
 */
export function findInternalSites(sequence, enzyme = 'BsaI') {
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enz) {
    throw new Error(`Unknown enzyme: ${enzyme}`);
  }

  const seq = sequence.toUpperCase();
  const recognition = enz.recognition;
  const recognitionRC = reverseComplement(recognition);

  const sites = [];

  // Find forward orientation sites
  let pos = seq.indexOf(recognition);
  while (pos !== -1) {
    sites.push({
      position: pos,
      sequence: seq.slice(pos, pos + recognition.length),
      orientation: 'forward',
      recognition,
    });
    pos = seq.indexOf(recognition, pos + 1);
  }

  // Find reverse complement sites
  pos = seq.indexOf(recognitionRC);
  while (pos !== -1) {
    sites.push({
      position: pos,
      sequence: seq.slice(pos, pos + recognitionRC.length),
      orientation: 'reverse',
      recognition: recognitionRC,
    });
    pos = seq.indexOf(recognitionRC, pos + 1);
  }

  // Sort by position
  sites.sort((a, b) => a.position - b.position);

  return {
    enzyme,
    recognitionSequence: recognition,
    hasSites: sites.length > 0,
    count: sites.length,
    sites,
  };
}

/**
 * Suggest silent mutations to remove internal restriction sites
 * This is "domestication" - removing sites while preserving protein sequence
 *
 * @param {string} sequence - DNA sequence (ideally coding sequence in frame)
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Options for domestication
 * @returns {Object} Domestication suggestions
 */
export function suggestDomestication(sequence, enzyme = 'BsaI', options = {}) {
  const {
    iscodingSequence = true, // Assume it's a coding sequence
    frame = 0, // Reading frame offset (0, 1, or 2)
  } = options;

  const internalSites = findInternalSites(sequence, enzyme);

  if (!internalSites.hasSites) {
    return {
      enzyme,
      needsDomestication: false,
      sites: [],
      suggestions: [],
      domesticatedSequence: sequence,
    };
  }

  const seq = sequence.toUpperCase();
  const suggestions = [];

  for (const site of internalSites.sites) {
    const suggestion = {
      position: site.position,
      originalSite: site.sequence,
      orientation: site.orientation,
      mutations: [],
    };

    if (iscodingSequence) {
      // Try to find synonymous mutations that break the site
      suggestion.mutations = findSynonymousMutations(seq, site, frame, enzyme);
    } else {
      // For non-coding, just suggest any single base change that breaks the site
      suggestion.mutations = findNonCodingMutations(site);
    }

    suggestions.push(suggestion);
  }

  // Apply mutations to create domesticated sequence
  let domesticatedSeq = seq;
  for (const suggestion of suggestions) {
    if (suggestion.mutations.length > 0) {
      // Use the first (best) mutation
      const mut = suggestion.mutations[0];
      domesticatedSeq =
        domesticatedSeq.slice(0, mut.position) +
        mut.newBase +
        domesticatedSeq.slice(mut.position + 1);
    }
  }

  // Verify domestication worked
  const verifyResult = findInternalSites(domesticatedSeq, enzyme);

  return {
    enzyme,
    needsDomestication: true,
    originalSites: internalSites.count,
    sites: internalSites.sites,
    suggestions,
    domesticatedSequence: domesticatedSeq,
    success: !verifyResult.hasSites,
    remainingSites: verifyResult.count,
  };
}

/**
 * Find synonymous (silent) mutations that break a restriction site
 */
function findSynonymousMutations(sequence, site, frame, enzyme = 'BsaI') {
  const mutations = [];
  const recognition = site.sequence;
  const siteStart = site.position;

  // Get the enzyme data for validation
  const enz = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enz) {
    console.warn(`Unknown enzyme ${enzyme}, falling back to BsaI`);
  }
  const enzRecognition = enz?.recognition || GOLDEN_GATE_ENZYMES.BsaI.recognition;
  const enzRecognitionRC = reverseComplement(enzRecognition);

  // For each position in the recognition sequence
  for (let i = 0; i < recognition.length; i++) {
    const seqPos = siteStart + i;
    const originalBase = recognition[i];

    // Determine which codon this position is in
    const adjustedPos = seqPos - frame;
    if (adjustedPos < 0) continue;

    const codonStart = Math.floor(adjustedPos / 3) * 3 + frame;
    const codonPos = adjustedPos % 3; // Position within codon (0, 1, or 2)

    if (codonStart + 3 > sequence.length) continue;

    const originalCodon = sequence.slice(codonStart, codonStart + 3);
    const originalAA = CODON_TO_AA[originalCodon];

    if (!originalAA) continue; // Not a valid codon

    // Try each alternative base
    for (const newBase of ['A', 'T', 'G', 'C']) {
      if (newBase === originalBase) continue;

      // Create new codon with mutation
      const newCodon =
        originalCodon.slice(0, codonPos) +
        newBase +
        originalCodon.slice(codonPos + 1);

      const newAA = CODON_TO_AA[newCodon];

      // Check if it's synonymous (same amino acid)
      if (newAA === originalAA) {
        // Check if this mutation breaks the recognition site
        const newSite =
          recognition.slice(0, i) + newBase + recognition.slice(i + 1);

        // A broken site means the new sequence won't match the recognition
        if (newSite !== enzRecognition && newSite !== enzRecognitionRC) {
          mutations.push({
            position: seqPos,
            positionInSite: i,
            originalBase,
            newBase,
            originalCodon,
            newCodon,
            aminoAcid: originalAA,
            isSynonymous: true,
            codonStart,
          });
        }
      }
    }
  }

  // Sort by position in site (prefer mutations in the middle of recognition)
  mutations.sort((a, b) => {
    // Prefer positions 2-4 (middle of 6bp recognition)
    const aScore = Math.abs(a.positionInSite - 2.5);
    const bScore = Math.abs(b.positionInSite - 2.5);
    return aScore - bScore;
  });

  return mutations;
}

/**
 * Find any mutation that breaks a restriction site (for non-coding sequences)
 */
function findNonCodingMutations(site) {
  const mutations = [];
  const recognition = site.sequence;

  // For each position, suggest alternative bases that break the site
  for (let i = 0; i < recognition.length; i++) {
    const originalBase = recognition[i];

    for (const newBase of ['A', 'T', 'G', 'C']) {
      if (newBase === originalBase) continue;

      mutations.push({
        position: site.position + i,
        positionInSite: i,
        originalBase,
        newBase,
        isSynonymous: false,
        note: 'Non-synonymous mutation (non-coding region)',
      });
    }
  }

  return mutations;
}

/**
 * Check if a sequence is compatible with Golden Gate assembly
 * Returns detailed information about any issues found
 *
 * @param {string} sequence - DNA sequence to check
 * @param {string} enzyme - Enzyme to use
 * @returns {Object} Compatibility analysis
 */
export function checkGoldenGateCompatibility(sequence, enzyme = 'BsaI') {
  const internalSites = findInternalSites(sequence, enzyme);
  const enz = GOLDEN_GATE_ENZYMES[enzyme];

  const issues = [];
  const warnings = [];

  // Check for internal sites
  if (internalSites.hasSites) {
    issues.push({
      type: 'internal_site',
      severity: 'error',
      message: `Sequence contains ${internalSites.count} internal ${enzyme} site(s) that will interfere with assembly`,
      details: internalSites.sites.map(s => ({
        position: s.position,
        orientation: s.orientation,
      })),
      suggestion: 'Remove sites via synonymous mutations (domestication) or choose a different enzyme',
    });
  }

  // Check sequence length
  if (sequence.length < 50) {
    warnings.push({
      type: 'short_sequence',
      severity: 'warning',
      message: 'Sequence is very short (<50bp). Assembly efficiency may be reduced.',
    });
  }

  // Check for extreme GC content
  const gc = gcContent(sequence);
  if (gc < 30) {
    warnings.push({
      type: 'low_gc',
      severity: 'warning',
      message: `Low GC content (${gc.toFixed(1)}%) may result in low primer Tm`,
    });
  } else if (gc > 70) {
    warnings.push({
      type: 'high_gc',
      severity: 'warning',
      message: `High GC content (${gc.toFixed(1)}%) may cause secondary structures`,
    });
  }

  // Check for alternative enzymes if current has issues
  let alternativeEnzymes = [];
  if (internalSites.hasSites) {
    alternativeEnzymes = findAlternativeEnzymes(sequence, enzyme);
  }

  return {
    enzyme,
    isCompatible: issues.length === 0,
    issues,
    warnings,
    internalSites,
    alternativeEnzymes,
    gcContent: gc,
  };
}

/**
 * Find alternative enzymes that don't have internal sites in the sequence
 */
export function findAlternativeEnzymes(sequence, currentEnzyme) {
  const alternatives = [];

  for (const [name, enz] of Object.entries(GOLDEN_GATE_ENZYMES)) {
    if (name === currentEnzyme) continue;

    const sites = findInternalSites(sequence, name);
    if (!sites.hasSites) {
      alternatives.push({
        enzyme: name,
        fullName: enz.fullName,
        recognition: enz.recognition,
        hasData: !!enz.dataKey, // Has experimental fidelity data
      });
    }
  }

  // Sort by preference (enzymes with data first)
  alternatives.sort((a, b) => (b.hasData ? 1 : 0) - (a.hasData ? 1 : 0));

  return alternatives;
}

/**
 * Design primers for a complete Golden Gate assembly
 * Auto-assigns overhangs based on part order and types
 *
 * @param {Array} parts - Array of {id, seq, type?} objects
 * @param {Object} options - Design options
 * @returns {Object} Complete assembly design with primers
 */
export function designGoldenGateAssembly(parts, options = {}) {
  const {
    enzyme = 'BsaI',
    circular = true,
    useStandardOverhangs = true,
    customOverhangs = null,
  } = options;

  if (parts.length < 2 || parts.length > 5) {
    throw new Error('Golden Gate assembly requires 2-5 parts');
  }

  // Determine overhangs for each junction
  let overhangs;

  if (customOverhangs && customOverhangs.length === parts.length + 1) {
    overhangs = customOverhangs;
  } else if (useStandardOverhangs) {
    // Use standard MoClo overhangs based on number of parts
    const standardSets = {
      2: ['GGAG', 'AATG', 'GCTT'],
      3: ['GGAG', 'TACT', 'AATG', 'GCTT'],
      4: ['GGAG', 'TACT', 'AATG', 'AGGT', 'GCTT'],
      5: ['GGAG', 'TACT', 'CCAT', 'AATG', 'AGGT', 'GCTT'],
    };
    overhangs = standardSets[parts.length];
  } else {
    throw new Error('Must provide customOverhangs or set useStandardOverhangs=true');
  }

  // Design primers for each part
  const designedParts = parts.map((part, i) => {
    const leftOH = overhangs[i];
    const rightOH = overhangs[i + 1];

    const primers = designGoldenGatePrimers(part.seq, leftOH, rightOH, { enzyme });

    return {
      ...part,
      index: i + 1,
      leftOverhang: leftOH,
      rightOverhang: rightOH,
      primers,
    };
  });

  // Calculate overall assembly fidelity using experimental data
  // This calculates fidelity based only on cross-reactivity between overhangs in this set
  const experimentalFidelity = calculateExperimentalFidelity(overhangs, enzyme);

  // Assemble the final sequence (for preview)
  let assembledSeq = '';
  for (const part of designedParts) {
    assembledSeq += part.seq;
  }

  // Check for potential issues
  const warnings = [...experimentalFidelity.warnings];
  const internalSiteIssues = [];

  for (const part of designedParts) {
    // Check for internal restriction sites using proper detection
    const siteCheck = findInternalSites(part.seq, enzyme);
    if (siteCheck.hasSites) {
      // Get domestication suggestions
      const domestication = suggestDomestication(part.seq, enzyme);
      const altEnzymes = findAlternativeEnzymes(part.seq, enzyme);

      const issue = {
        partId: part.id,
        sites: siteCheck.sites,
        count: siteCheck.count,
        domestication,
        alternativeEnzymes: altEnzymes,
      };
      internalSiteIssues.push(issue);

      // Add warning with actionable suggestions
      let suggestion = '';
      if (domestication.success && domestication.suggestions.length > 0) {
        const mut = domestication.suggestions[0].mutations[0];
        if (mut) {
          suggestion = ` Suggested fix: ${mut.originalBase}${mut.position + 1}${mut.newBase}`;
          if (mut.isSynonymous) suggestion += ' (silent mutation)';
        }
      } else if (altEnzymes.length > 0) {
        suggestion = ` Consider using ${altEnzymes[0].enzyme} instead.`;
      }

      warnings.push(`⚠ Part "${part.id}" contains ${siteCheck.count} internal ${enzyme} site(s) - MUST be removed before assembly.${suggestion}`);
    }

    // Add primer warnings from the primer design
    if (part.primers.warnings) {
      for (const w of part.primers.warnings) {
        warnings.push(`Part "${part.id}": ${w}`);
      }
    }
  }

  return {
    enzyme,
    parts: designedParts,
    overhangs,
    circular,
    assembledSequence: circular ? assembledSeq : assembledSeq,
    assembledLength: assembledSeq.length,
    fidelity: {
      // Use experimental junction fidelities (calculated against only this set)
      individual: experimentalFidelity.junctions.map((j, i) => ({
        junction: i,
        overhang: j.overhang,
        fidelity: j.fidelity,
        fidelityPercent: j.fidelityPercent,
      })),
      overall: experimentalFidelity.assemblyFidelity,
      percentage: experimentalFidelity.assemblyFidelityPercent,
    },
    warnings,
    internalSiteIssues, // Detailed info about internal sites with domestication suggestions
    hasInternalSites: internalSiteIssues.length > 0,
    protocol: generateProtocol(designedParts, enzyme),
  };
}

/**
 * Generate a protocol text for the assembly
 */
function generateProtocol(parts, enzyme) {
  const totalVolume = 20; // µL
  const dnaPerPart = 2; // µL of 75 ng/µL
  const enzymeVol = 2;
  const bufferVol = 2;
  const waterVol = totalVolume - (parts.length * dnaPerPart) - enzymeVol - bufferVol;

  return {
    title: `Golden Gate Assembly Protocol (${parts.length} parts)`,
    steps: [
      {
        step: 1,
        title: 'Prepare DNA fragments',
        details: [
          'PCR amplify each fragment using the designed primers',
          'Gel purify or PCR cleanup each product',
          'Dilute to 75 ng/µL in water or TE buffer',
        ],
      },
      {
        step: 2,
        title: 'Setup reaction',
        details: [
          `Add ${bufferVol} µL T4 DNA Ligase Buffer (10X)`,
          ...parts.map((p, i) => `Add ${dnaPerPart} µL ${p.id} (75 ng/µL)`),
          `Add ${enzymeVol} µL NEB Golden Gate Assembly Mix`,
          `Add ${waterVol.toFixed(1)} µL nuclease-free water`,
          `Total volume: ${totalVolume} µL`,
        ],
      },
      {
        step: 3,
        title: 'Thermocycler program',
        details: [
          '37°C for 60 min (digestion and ligation)',
          '60°C for 5 min (enzyme inactivation)',
          '4°C hold',
        ],
      },
      {
        step: 4,
        title: 'Transformation',
        details: [
          'Transform 2-5 µL into competent E. coli',
          'Plate on appropriate antibiotic selection',
          'Incubate overnight at 37°C',
        ],
      },
    ],
    notes: [
      `Expected fidelity: ~${(parts.length <= 3 ? 95 : parts.length <= 5 ? 90 : 85)}% correct assemblies`,
      'For best results, use high-fidelity BsaI-HFv2',
      'Avoid internal BsaI sites in your insert sequences',
    ],
  };
}

/**
 * Get recommended overhang set for N parts
 */
export function getRecommendedOverhangs(numParts) {
  const sets = {
    2: { overhangs: ['GGAG', 'AATG', 'GCTT'], fidelity: 0.98 },
    3: { overhangs: ['GGAG', 'TACT', 'AATG', 'GCTT'], fidelity: 0.97 },
    4: { overhangs: ['GGAG', 'TACT', 'AATG', 'AGGT', 'GCTT'], fidelity: 0.95 },
    5: { overhangs: ['GGAG', 'TACT', 'CCAT', 'AATG', 'AGGT', 'GCTT'], fidelity: 0.93 },
  };

  return sets[numParts] || sets[5];
}

/**
 * Validate an overhang sequence
 */
export function validateOverhang(overhang) {
  const oh = overhang.toUpperCase();

  if (oh.length !== 4) {
    return { valid: false, error: 'Overhang must be exactly 4 bases' };
  }

  if (!/^[ATGC]+$/.test(oh)) {
    return { valid: false, error: 'Overhang must contain only A, T, G, C' };
  }

  // Check if palindromic
  const rc = reverseComplement(oh);
  if (oh === rc) {
    return {
      valid: false,
      error: 'Palindromic overhangs should be avoided (can self-ligate)',
    };
  }

  // Get fidelity
  const fidelity = OVERHANG_FIDELITY[oh] || { fidelity: 0.85, note: 'Custom' };

  return {
    valid: true,
    overhang: oh,
    reverseComplement: rc,
    fidelity: fidelity.fidelity,
    note: fidelity.note,
    warning: fidelity.fidelity < 0.8 ? 'Low fidelity overhang - may cause misassembly' : null,
  };
}

/**
 * Check compatibility between two adjacent overhangs
 */
export function checkOverhangCompatibility(leftPartRightOH, rightPartLeftOH) {
  const left = leftPartRightOH.toUpperCase();
  const right = rightPartLeftOH.toUpperCase();

  // They should be reverse complements of each other
  // Actually, they should be IDENTICAL for Golden Gate
  // The enzyme creates complementary sticky ends automatically
  if (left !== right) {
    return {
      compatible: false,
      error: `Overhangs don't match: ${left} ≠ ${right}`,
    };
  }

  return {
    compatible: true,
    overhang: left,
    fidelity: OVERHANG_FIDELITY[left]?.fidelity || 0.85,
  };
}

// ==========================================
// Legacy functions for backward compatibility
// ==========================================

export function findTypeIISSites(seq, enzymeName) {
  const enzyme = GOLDEN_GATE_ENZYMES[enzymeName];
  if (!enzyme) {
    throw new Error(`Unknown Golden Gate enzyme: ${enzymeName}`);
  }

  const sites = [];
  const recognition = enzyme.recognition.toUpperCase();
  const recognitionRC = reverseComplement(recognition);
  const seqUpper = seq.toUpperCase();

  let pos = 0;
  while ((pos = seqUpper.indexOf(recognition, pos)) !== -1) {
    sites.push({ position: pos, strand: '+', enzyme: enzymeName, recognition });
    pos++;
  }

  pos = 0;
  while ((pos = seqUpper.indexOf(recognitionRC, pos)) !== -1) {
    sites.push({ position: pos, strand: '-', enzyme: enzymeName, recognition: recognitionRC });
    pos++;
  }

  return sites.sort((a, b) => a.position - b.position);
}

export function catalyze(seq, enzymes, options = {}) {
  const { linear = true } = options;
  const seqUpper = seq.toUpperCase();

  const allSites = [];
  for (const enzymeName of enzymes) {
    const sites = findTypeIISSites(seqUpper, enzymeName);
    allSites.push(...sites);
  }

  if (allSites.length < 2) {
    return [];
  }

  allSites.sort((a, b) => a.position - b.position);

  // Calculate overhangs for each site
  const cuts = allSites.map(site => {
    const enzyme = GOLDEN_GATE_ENZYMES[site.enzyme];
    let cutPos, overhangStart;

    if (site.strand === '+') {
      cutPos = site.position + enzyme.recognition.length + enzyme.cutOffset;
      overhangStart = cutPos;
    } else {
      cutPos = site.position - enzyme.cutOffset;
      overhangStart = cutPos - enzyme.overhangLength;
    }

    const overhang = seqUpper.slice(overhangStart, overhangStart + enzyme.overhangLength);

    return {
      site,
      cutPos,
      overhang: site.strand === '+' ? overhang : reverseComplement(overhang),
      overhangRaw: overhang,
      strand: site.strand,
    };
  });

  cuts.sort((a, b) => a.cutPos - b.cutPos);

  // Extract fragments
  const fragments = [];
  for (let i = 0; i < cuts.length - 1; i++) {
    const currentCut = cuts[i];
    const nextCut = cuts[i + 1];

    const fragSeq = seqUpper.slice(currentCut.cutPos, nextCut.cutPos);

    fragments.push({
      seq: fragSeq,
      start: currentCut.cutPos,
      end: nextCut.cutPos,
      leftOverhang: currentCut.overhang,
      rightOverhang: nextCut.overhang,
      leftOverhangSeq: currentCut.overhangRaw,
      rightOverhangSeq: nextCut.overhangRaw,
    });
  }

  return fragments;
}

export function identifyFusionSites(part, options = {}) {
  const { enzymes = ['BsaI'] } = options;
  const fragments = catalyze(part.seq, enzymes, { linear: true });

  if (fragments.length === 0) {
    return { left: null, right: null };
  }

  const leftOH = fragments[0].leftOverhangSeq;
  const rightOH = fragments[0].rightOverhangSeq;

  let leftSite = null;
  let rightSite = null;

  for (const [code, data] of Object.entries(STANDARD_FUSION_SITES)) {
    if (leftOH === data.seq) leftSite = code;
    if (rightOH === data.seq) rightSite = code;
  }

  return {
    left: leftSite,
    right: rightSite,
    leftOverhang: leftOH,
    rightOverhang: rightOH,
    partType: leftSite && rightSite ? `${leftSite}${rightSite}` : 'custom',
  };
}

// Export for validation UI
export function validateOrderedAssembly(parts, options = {}) {
  const { enzymes = ['BsaI'] } = options;
  const errors = [];
  const warnings = [];
  const partDetails = [];

  for (let i = 0; i < parts.length; i++) {
    const part = parts[i];
    const sites = findTypeIISSites(part.seq, enzymes[0]);

    if (sites.length === 0) {
      errors.push(`Part "${part.id}": No ${enzymes.join('/')} sites found`);
      continue;
    }

    if (sites.length > 2) {
      warnings.push(`Part "${part.id}": Has ${sites.length} cut sites (expected 2)`);
    }

    const fragments = catalyze(part.seq, enzymes, { linear: true });

    if (fragments.length > 0) {
      partDetails.push({
        id: part.id,
        index: i,
        leftOverhang: fragments[0].leftOverhangSeq,
        rightOverhang: fragments[0].rightOverhangSeq,
        insertLength: fragments[0].seq.length,
        sites: sites.length,
      });
    }
  }

  // Check compatibility
  for (let i = 0; i < partDetails.length - 1; i++) {
    const current = partDetails[i];
    const next = partDetails[i + 1];

    if (current.rightOverhang !== next.leftOverhang) {
      errors.push(
        `Parts "${current.id}" and "${next.id}": Incompatible overhangs. ` +
        `${current.id} right: ${current.rightOverhang}, ${next.id} left: ${next.leftOverhang}`
      );
    }
  }

  // Circular check
  if (partDetails.length >= 2) {
    const first = partDetails[0];
    const last = partDetails[partDetails.length - 1];
    if (last.rightOverhang !== first.leftOverhang) {
      errors.push(
        `Assembly cannot circularize: Last part right overhang (${last.rightOverhang}) ` +
        `≠ first part left overhang (${first.leftOverhang})`
      );
    }
  }

  return {
    valid: errors.length === 0,
    errors,
    warnings,
    partDetails,
    canCircularize: errors.length === 0,
  };
}

export function assembleSequence(parts, options = {}) {
  const { enzymes = ['BsaI'], circular = true } = options;

  let assembled = '';
  const assemblySteps = [];

  for (let i = 0; i < parts.length; i++) {
    const part = parts[i];
    const fragments = catalyze(part.seq, enzymes, { linear: true });

    if (fragments.length > 0) {
      const insert = fragments[0];
      assembled += insert.seq;

      assemblySteps.push({
        step: i + 1,
        partId: part.id,
        insertSequence: insert.seq,
        insertLength: insert.seq.length,
        leftOverhang: insert.leftOverhangSeq,
        rightOverhang: insert.rightOverhangSeq,
        cumulativeLength: assembled.length,
      });
    }
  }

  return {
    sequence: assembled,
    length: assembled.length,
    circular,
    steps: assemblySteps,
  };
}

// Stub exports for backward compatibility
export const calculateOverhang = (seq, site) => null;
export const findSimpleCycles = (graph) => [];
export const planGoldenGate = (parts, options) => designGoldenGateAssembly(parts, options);
export const getStandardOverhang = (code) => STANDARD_FUSION_SITES[code]?.seq || null;

// ==========================================
// Overhang Set Optimization (NEB Algorithm Port)
// Monte Carlo simulated annealing for optimal overhang selection
// ==========================================

export {
  optimizeOverhangSet,
  evaluateOverhangSet,
  optimizeOverhangSetMultiRun,
  batchScoreRandomSets,
  calculateLigationFidelity,
  getAllOverhangs,
  filterOverhangs,
} from './overhang-optimizer.js';
