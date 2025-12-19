/**
 * Repp - DNA Assembly Planning Library
 * JavaScript port of https://github.com/Lattice-Automation/repp
 *
 * Finds optimal fragment combinations for plasmid construction via Gibson Assembly
 */

export { DEFAULT_CONFIG, mergeConfig, synthFragmentCost, synthPlasmidCost, primerCost } from './config.js';
export { ENZYMES, digest, getEnzymeNames, parseEnzyme, reverseComplement } from './enzymes.js';
export { parseFasta, parseGenbank, findMatches, findJunction, gcContent, toFasta } from './sequence.js';
export { FragType, createFragment, fragmentCost, planAssembly, createAssemblies, fillAssemblies } from './assembly.js';

// Part Sources (Addgene, iGEM, DNASU)
export {
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
} from './partSources.js';

// Golden Gate Assembly
export {
  GOLDEN_GATE_ENZYMES,
  STANDARD_FUSION_SITES,
  HIGH_FIDELITY_SETS,
  OVERHANG_FIDELITY,
  // PART_TYPES,  // COMMENTED: Not exported from goldengate.js
  // findTypeIISSites,  // COMMENTED: Not exported from goldengate.js
  // catalyze,  // COMMENTED: Not exported from goldengate.js
  // identifyFusionSites,  // COMMENTED: Not exported from goldengate.js
  // validateOrderedAssembly,  // COMMENTED: Not exported from goldengate.js
  // assembleSequence,  // COMMENTED: Not exported from goldengate.js
  // Primer design functions
  // designGoldenGatePrimers,  // COMMENTED: Not exported from goldengate.js
  // designGoldenGateAssembly,  // COMMENTED: Not exported from goldengate.js
  // getRecommendedOverhangs,  // COMMENTED: Not exported from goldengate.js
  // validateOverhang,  // COMMENTED: Not exported from goldengate.js
  // checkOverhangCompatibility,  // COMMENTED: Not exported from goldengate.js
  // NEB-validated overhang set functions
  // getRecommendedOverhangSet,  // COMMENTED: Not exported from goldengate.js
  // calculateSetFidelity,  // COMMENTED: Not exported from goldengate.js
  // Experimental ligation data (Pryor et al. 2020)
  ligationData,
  ENZYMES_WITH_DATA,
  getEnzymeLigationData,
  getOverhangFidelityExperimental,
  getLigationFrequency,
  calculateExperimentalFidelity,
  getOptimalOverhangSetExperimental,
  findOptimalOverhangSet,
  compareEnzymeFidelity,
  findProblematicPairs,
  // Cross-ligation analysis and visualization
  // generateCrossLigationHeatmap,  // COMMENTED: Not exported from goldengate.js
  // generateOverhangQualityReport,  // COMMENTED: Not exported from goldengate.js
  // Internal site detection and domestication
  findInternalSites,
  // suggestDomestication,  // COMMENTED: Not exported from goldengate.js
  // checkGoldenGateCompatibility,  // COMMENTED: Not exported from goldengate.js
  // findAlternativeEnzymes,  // COMMENTED: Not exported from goldengate.js
  // Legacy/compatibility
  // calculateOverhang,  // COMMENTED: Not exported from goldengate.js
  // findSimpleCycles,  // COMMENTED: Not exported from goldengate.js
  // planGoldenGate,  // COMMENTED: Not exported from goldengate.js
  // getStandardOverhang,  // COMMENTED: Not exported from goldengate.js
} from './goldengate.js';

// Golden Gate Primer Optimizer (state-of-the-art optimizations)
export {
  // Optimal flanking sequences (6bp NEB recommendation)
  OPTIMAL_FLANKING_SEQUENCES,
  OPTIMAL_SPACERS,
  selectOptimalFlankingSequence,
  // G:T mismatch detection
  canPair,
  findGTMismatchRisks,
  calculateEnhancedFidelity,
  // Automatic overhang optimization
  autoOptimizeOverhangs,
  // Golden Gate primer quality scoring
  scoreGoldenGatePrimer,
  scoreGoldenGatePrimerPair,
  // Optimized design functions (recommended)
  designOptimizedGoldenGatePrimers,
  designOptimizedGoldenGateAssembly,
  // UI-compatible adapter
  designOptimizedGoldenGateAssemblyForUI,
} from './goldengate-primer-optimizer.js';

// Fusion Site Optimizer (automatic junction position finding)
export {
  // Main optimization entry point
  optimizeFusionSites,
  // Algorithm-specific functions
  optimizeGreedy,
  optimizeBranchBound,
  optimizeMonteCarlo,
  // Quick optimization for simple cases
  quickOptimize,
  // Configuration defaults
  OPTIMIZER_DEFAULTS,
} from './fusion-site-optimizer.js';

// Fusion Site Scanner (candidate generation)
export {
  scanForFusionSites,
  generateTargetPositions,
  filterByDistance,
  assessFeasibility,
} from './fusion-site-scanner.js';

// Fusion Site Scorer (composite scoring)
export {
  scoreFusionSiteComposite,
  scoreMultipleFusionSites,
  quickScoreFusionSite,
  DEFAULT_FUSION_WEIGHTS,
} from './fusion-site-scorer.js';

// Failure Prediction
export {
  predictFailureModes,
  quickRiskAssessment,
} from './failure-prediction.js';

// Auto-Domestication Optimizer
export {
  analyzeForDomestication,
  findDomesticationJunctions,
  optimizeWithDomestication,
  getDomesticationSummary,
  validateDomesticationOverhangs,
  // Adjacent site detection and enzyme recommendations
  detectAdjacentSites,
  recommendAlternativeEnzymes,
  validateFragmentSizes,
  // Global overhang optimization
  optimizeGlobalOverhangSet,
  // Adjacent site merging
  mergeAdjacentSites,
  // Post-domestication validation
  validatePostDomestication,
  DOMESTICATION_DEFAULTS,
} from './auto-domestication-optimizer.js';

// Silent Mutation Domestication (preferred for one-pot Golden Gate)
export {
  // Main entry point
  domesticateWithSilentMutations,
  // Candidate finding and scoring
  findAllSilentMutationCandidates,
  scoreMutationCandidates,
  // Verification
  verifyProteinSequence,
  validateDomestication as validateSilentDomestication,
  // Utilities
  applyMutation,
  applyMutations,
  // Codon tables
  CODON_TABLE,
  CODON_TO_AA,
  ECOLI_CODON_USAGE,
  YEAST_CODON_USAGE,
  // Strategy comparison
  compareDomesticationStrategies,
} from './silent-mutation-domesticator.js';

// Unified Domestication Optimizer (recommended entry point)
export {
  // Main optimization function (auto-selects best strategy)
  optimizeDomestication,
  // Detailed analysis
  analyzeDomesticationOptions,
  generateDomesticationReport,
  // Strategy constants
  DOMESTICATION_STRATEGY,
  UNIFIED_DOMESTICATION_CONFIG,
  // Mutagenic junction functions (PREFERRED for sites in middle of fragments)
  // designMutagenicJunction,  // COMMENTED: Not exported from domestication-optimizer.js
  // designAllMutagenicJunctions,  // COMMENTED: Not exported from domestication-optimizer.js
  // selectDomesticationStrategy,  // COMMENTED: Not exported from domestication-optimizer.js (use DOMESTICATION_STRATEGY instead)
  // MUTAGENIC_JUNCTION_CONFIG,  // COMMENTED: Not exported from domestication-optimizer.js
} from './domestication-optimizer.js';

// Enhanced Domestication System (STATE-OF-THE-ART with user confirmation)
export {
  // Main workflow functions
  createDomesticationPlan,
  executeDomesticationPlan,
  // Configuration
  ENHANCED_CONFIG,
} from './enhanced-domestication.js';

// Enhanced Mutagenic Junction (STATE-OF-THE-ART with NEB fidelity data)
export {
  // Main design functions
  designEnhancedMutagenicJunction,
  designAllEnhancedJunctions,
  // Scoring functions
  scoreOverhangFidelity,
  scorePrimerQuality,
  classifyQuality as classifyJunctionQuality,
  // Configuration
  ENHANCED_JUNCTION_CONFIG,
} from './enhanced-mutagenic-junction.js';

// Domestication Primer Workflow (COMPLETE INTEGRATED WORKFLOW)
export {
  // Main workflow entry point (RECOMMENDED)
  runDomesticationWorkflow,
  // Analysis functions
  analyzeSequenceForDomestication,
  // Strategy execution
  selectAndExecuteStrategy,
  // Primer design
  designIntegratedPrimers,
  // Validation
  validateCompleteDesign,
  // Workflow guide generation
  generateWorkflowGuide,
  // Configuration
  WORKFLOW_CONFIG,
} from './domestication-primer-workflow.js';

// ORF Detection and Frame Validation
export {
  detectORFs,
  validateReadingFrame,
  translateSequence,
  analyzeSiteCodonContext,
  preFlightCheck,
  compareProteins,
  formatProteinForDisplay,
} from './orf-detector.js';
