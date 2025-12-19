/**
 * Unified Primer Analysis Presets
 *
 * Single source of truth for domain-specific analysis configurations.
 * Used by scoring.js, primerAnalysis.js, and other analysis modules.
 *
 * Each preset defines:
 * - Thresholds for thermodynamic parameters (hairpin, dimer ΔG)
 * - Optimal ranges for Tm, GC content, length
 * - Tm calculation method
 * - Scoring weights
 */

import { DEFAULT_WEIGHTS } from './weightCalibration.js';

/**
 * Domain-specific analysis presets
 *
 * Based on:
 * - Published literature (IDT, Premier Biosoft, NEB guidelines)
 * - Empirical validation from the Döring dataset
 * - Domain-specific requirements
 */
export const ANALYSIS_PRESETS = {
  amplification: {
    name: 'Amplification (PCR)',
    description: 'Standard PCR amplification primers',
    thresholds: {
      hairpin: -3.0,
      hairpin3Prime: -2.0,
      homodimer: -6.0,
      homodimer3Prime: -5.0,
      heterodimer: -6.0,
      heterodimer3Prime: -5.0,
    },
    optimalTmRange: [55, 60],
    optimalGcRange: [40, 60],
    optimalLengthRange: [18, 24],
    tmMethod: 'santaLucia',
    weights: DEFAULT_WEIGHTS,
    // Scoring function options (for backward compatibility with SCORING_PRESETS format)
    tmOptions: { optimalLow: 55, optimalHigh: 60 },
    gcOptions: { optimalLow: 40, optimalHigh: 60 },
    lengthOptions: { optimalLow: 18, optimalHigh: 24 },
    hairpinThreshold: -3.0,
    homodimerThreshold: -6.0,
    heterodimerThreshold: -6.0,
  },

  mutagenesis: {
    name: 'Site-Directed Mutagenesis',
    description: 'QuikChange-style mutagenesis primers',
    thresholds: {
      hairpin: -3.0,
      hairpin3Prime: -2.0,
      homodimer: -6.0,
      homodimer3Prime: -5.0,
      heterodimer: -6.0,
      heterodimer3Prime: -5.0,
      // Critical thresholds for severe warnings
      hairpinCritical: -6.0,
      homodimerCritical: -9.0,
      heterodimerCritical: -9.0,
    },
    optimalTmRange: [50, 72],
    optimalGcRange: [40, 60],
    optimalLengthRange: [25, 45],
    tmMethod: 'q5',
    weights: {
      ...DEFAULT_WEIGHTS,
      terminal3DG: 0.25,  // 3' end stability is critical
      heterodimer: 0.15,  // Heterodimer more important in QuikChange
    },
    // Scoring function options
    tmOptions: { optimalLow: 50, optimalHigh: 72 },
    gcOptions: { optimalLow: 40, optimalHigh: 60 },
    lengthOptions: { optimalLow: 25, optimalHigh: 45 },
    hairpinThreshold: -3.0,
    homodimerThreshold: -6.0,
    heterodimerThreshold: -6.0,
  },

  sequencing: {
    name: 'Sanger Sequencing',
    description: 'Primers for Sanger sequencing',
    thresholds: {
      hairpin: -2.5,
      hairpin3Prime: -1.5,
      homodimer: -5.0,
      homodimer3Prime: -4.0,
      heterodimer: -5.0,
      heterodimer3Prime: -4.0,
    },
    optimalTmRange: [55, 60],
    optimalGcRange: [40, 60],
    optimalLengthRange: [18, 24],
    tmMethod: 'santaLucia',
    weights: {
      ...DEFAULT_WEIGHTS,
      offTarget: 0.30,  // More important for sequencing specificity
    },
    // Scoring function options
    tmOptions: { optimalLow: 55, optimalHigh: 60 },
    gcOptions: { optimalLow: 40, optimalHigh: 60 },
    lengthOptions: { optimalLow: 18, optimalHigh: 24 },
    hairpinThreshold: -2.5,
    homodimerThreshold: -5.0,
    heterodimerThreshold: -5.0,
  },

  assembly: {
    name: 'Gibson/NEBuilder Assembly',
    description: 'Gibson Assembly or NEBuilder HiFi primers with overlap regions',
    thresholds: {
      hairpin: -3.0,
      hairpin3Prime: -2.0,
      homodimer: -6.0,
      homodimer3Prime: -5.0,
      heterodimer: -6.0,
      heterodimer3Prime: -5.0,
    },
    optimalTmRange: [48, 65],
    optimalGcRange: [40, 60],
    optimalLengthRange: [18, 30],
    tmMethod: 'q5',
    weights: {
      ...DEFAULT_WEIGHTS,
      heterodimer: 0.15,  // Important for multi-fragment assemblies
    },
    // Scoring function options
    tmOptions: { optimalLow: 48, optimalHigh: 65 },
    gcOptions: { optimalLow: 40, optimalHigh: 60 },
    lengthOptions: { optimalLow: 18, optimalHigh: 30 },
    hairpinThreshold: -3.0,
    homodimerThreshold: -6.0,
    heterodimerThreshold: -6.0,
  },

  goldengate: {
    name: 'Golden Gate Assembly',
    description: 'Golden Gate primers with Type IIS enzyme sites (BsaI, BsmBI, etc.)',
    thresholds: {
      hairpin: -3.0,
      hairpin3Prime: -2.0,
      homodimer: -6.0,
      homodimer3Prime: -5.0,
      heterodimer: -6.0,
      heterodimer3Prime: -5.0,
    },
    // Golden Gate uses specific annealing conditions - NEB recommends 50-60°C for annealing region
    optimalTmRange: [50, 60],
    optimalGcRange: [40, 60],
    // Annealing region typically 18-25bp for Golden Gate
    optimalLengthRange: [18, 25],
    tmMethod: 'q5',
    weights: {
      ...DEFAULT_WEIGHTS,
      heterodimer: 0.15,  // Important for multi-fragment assemblies
    },
    // Scoring function options - score based on annealing region only
    tmOptions: { optimalLow: 50, optimalHigh: 60 },
    gcOptions: { optimalLow: 40, optimalHigh: 60 },
    lengthOptions: { optimalLow: 18, optimalHigh: 25 },
    hairpinThreshold: -3.0,
    homodimerThreshold: -6.0,
    heterodimerThreshold: -6.0,
    // Golden Gate specific flag
    isGoldenGate: true,
  },
};

/**
 * Get preset by mode name
 * @param {string} mode - Preset mode name
 * @returns {Object} Preset configuration
 */
export function getPreset(mode) {
  return ANALYSIS_PRESETS[mode] || ANALYSIS_PRESETS.amplification;
}

/**
 * Legacy alias for backward compatibility
 * @deprecated Use ANALYSIS_PRESETS instead
 */
export const SCORING_PRESETS = ANALYSIS_PRESETS;

/**
 * Application-specific presets for different PCR use cases
 * These augment the base ANALYSIS_PRESETS with application-specific weight adjustments
 */
export const APPLICATION_PRESETS = {
  default: {
    id: 'default',
    name: 'Standard PCR',
    description: 'Balanced settings for general PCR applications',
    weightMultipliers: {},
    icon: 'beaker',
  },
  qpcr: {
    id: 'qpcr',
    name: 'qPCR / RT-qPCR',
    description: 'Optimized for quantitative PCR with strict efficiency requirements',
    weightMultipliers: {
      tm: 1.3,           // Precise Tm matching critical
      gcContent: 1.2,    // Avoid GC-rich regions
      heterodimer: 1.4,  // No primer-dimers allowed
      homodimer: 1.3,    // Self-complementarity affects efficiency
    },
    icon: 'chart',
  },
  colonyPcr: {
    id: 'colonyPcr',
    name: 'Colony PCR',
    description: 'Robust amplification from crude colony lysates',
    weightMultipliers: {
      tm: 0.9,           // Slightly relaxed Tm
      hairpin: 0.8,      // Less stringent secondary structure
      homodimer: 0.8,
    },
    icon: 'bacteria',
  },
  longRange: {
    id: 'longRange',
    name: 'Long-Range PCR',
    description: 'Primers for amplicons >5kb with specialized polymerases',
    weightMultipliers: {
      tm: 1.2,           // Consistent Tm important
      terminal3DG: 1.4,  // 3\' stability critical for processivity
      hairpin: 1.1,
    },
    icon: 'ruler',
  },
  gcRich: {
    id: 'gcRich',
    name: 'GC-Rich Templates',
    description: 'Templates with >65% GC content requiring special handling',
    weightMultipliers: {
      hairpin: 1.5,      // Secondary structures more likely
      homodimer: 1.3,
      tm: 1.2,           // Tm accuracy affected by high GC
    },
    icon: 'dna',
  },
  atRich: {
    id: 'atRich',
    name: 'AT-Rich Templates',
    description: 'Templates with >65% AT content',
    weightMultipliers: {
      terminal3DG: 1.3,  // 3\' stability more challenging
      tm: 1.1,
    },
    icon: 'dna',
  },
  multiplex: {
    id: 'multiplex',
    name: 'Multiplex PCR',
    description: 'Multiple primer pairs in single reaction',
    weightMultipliers: {
      heterodimer: 1.6,  // Cross-dimerization critical
      tm: 1.4,           // All primers must match Tm
      homodimer: 1.2,
    },
    icon: 'layers',
  },
  highFidelity: {
    id: 'highFidelity',
    name: 'High-Fidelity',
    description: 'Maximum specificity for cloning and sensitive applications',
    weightMultipliers: {
      offTarget: 1.5,    // Specificity paramount
      terminal3DG: 1.3,
      heterodimer: 1.2,
    },
    icon: 'shield',
  },
};

/**
 * Get list of available application presets for UI selection
 * @returns {Array} Array of preset objects with id, name, description
 */
export function getAvailablePresets() {
  return Object.values(APPLICATION_PRESETS).map(preset => ({
    id: preset.id,
    name: preset.name,
    description: preset.description,
    icon: preset.icon,
  }));
}

/**
 * Auto-detect the best preset based on template characteristics
 * @param {Object} params - Parameters for detection
 * @param {string} params.sequence - Template sequence
 * @param {number} params.ampliconLength - Expected amplicon length
 * @param {string} params.applicationHint - User hint about application
 * @returns {string} Recommended preset ID
 */
export function autoDetectPreset({ sequence = '', ampliconLength = 0, applicationHint = '' } = {}) {
  const hint = applicationHint.toLowerCase();

  // Check for explicit hints
  if (hint.includes('qpcr') || hint.includes('quantitative')) return 'qpcr';
  if (hint.includes('colony')) return 'colonyPcr';
  if (hint.includes('multiplex')) return 'multiplex';
  if (hint.includes('fidelity') || hint.includes('cloning')) return 'highFidelity';

  // Check amplicon length for long-range
  if (ampliconLength > 5000) return 'longRange';

  // Check GC content of sequence
  if (sequence && sequence.length > 50) {
    const gcCount = (sequence.match(/[GC]/gi) || []).length;
    const gcPercent = (gcCount / sequence.length) * 100;
    if (gcPercent > 65) return 'gcRich';
    if (gcPercent < 35) return 'atRich';
  }

  return 'default';
}
