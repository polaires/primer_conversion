/**
 * Create or score PCR primers
 * Ported from primers Python library
 */

import { gcCache, tmCache } from "./tm.js";
import { dgCache } from "./fold.js";
import { offTargets } from "./offTargets.js";
import { calculate3primeTerminalDG } from "./tmQ5.js";
import {
  calculateEquilibriumEfficiency,
  calculateHairpinDG,
  calculateHomodimerDG,
  calculateHeterodimerDG,
  efficiencyToScore,
  reverseComplement as rcSeq,
} from "./equilibrium.js";
import { analyzeOffTargetsPair } from "./offTargetClassification.js";
import {
  scoreTm,
  scoreGc,
  scoreTerminal3DG,
  scoreTmDiff,
  scoreHairpin,
  scoreHomodimer,
  scoreHeterodimer,
  scoreOffTarget,
  scoreLength,
  scoreGcClamp,
  scoreHomopolymer,
  scoreOffTargetClassification,
  score3PrimeComposition,
  scoreGQuadruplex,
  analyzeGQuadruplex,
  calculateCompositeScore,
  classifyQuality,
} from "./scoring.js";
import { ANALYSIS_PRESETS } from "./presets.js";
import {
  analyze3PrimeEnd,
  optimizePrimer,
  optimizePrimerPair,
  scorePrimerVariant,
  generateDesignSuggestions,
  quickAssess,
  SMART_DESIGN_CONSTRAINTS,
} from "./smartPrimers.js";

export const LEN_MIN = 15; // min length of the annealing portion of primers
export const LEN_MAX = 32; // default max length for optimal PCR efficiency
export const LEN_MAX_EXTENDED = 60; // extended max length for challenging sequences (IDT supports up to 60bp standard)

/**
 * IDT Primer Length Guidelines:
 * - Standard synthesis: up to 60bp (no extra cost)
 * - Extended synthesis: up to 100bp (special order)
 *
 * The default LEN_MAX of 32bp is a conservative choice for optimal PCR efficiency:
 * - Shorter primers have faster annealing kinetics
 * - Less risk of secondary structure formation
 * - Lower synthesis cost per base
 *
 * However, for challenging sequences (high/low GC content, Tm mismatches),
 * extending primers up to 60bp can significantly improve Tm matching.
 * The allowExtendedPrimers option (default: true) enables this adaptive behavior.
 */

// Type definitions
export interface EquilibriumLosses {
  hairpin?: number;
  homodimer?: number;
  heterodimer?: number;
  heterodimerDG?: number;
  offTarget?: number;
  free?: number;
}

export interface PiecewiseScores {
  tm?: number;
  gc?: number;
  length?: number;
  gcClamp?: number;
  homopolymer?: number;
  hairpin?: number;
  homodimer?: number;
  heterodimer?: number;
  offTarget?: number;
  terminal3DG?: number;
  tmDiff?: number;
  gQuadruplex?: number;
  threePrimeComp?: number;
  selfDimerFwd?: number;
  selfDimerRev?: number;
  gQuadruplexFwd?: number;
  gQuadruplexRev?: number;
}

export interface GQuadruplexAnalysis {
  hasG4Motif: boolean;
  score: number;
  details?: any;
}

export interface Scoring {
  penalty: number;
  penaltyTm: number;
  penaltyTmDiff: number;
  penaltyGc: number;
  penaltyLen: number;
  penaltyDg: number;
  penaltyOffTarget: number;
  equilibriumEfficiency: number | null;
  equilibriumScore: number | null;
  equilibriumLosses: EquilibriumLosses | null;
  piecewiseScores: PiecewiseScores | null;
  compositeScore: number | null;
  effectiveScore?: number;
  criticalWarnings?: number;
  qualityTier: string | null;
  gQuadruplex?: GQuadruplexAnalysis;
}

export interface SmartDesignMetadata {
  optimized: boolean;
  reason?: string;
  originalAnalysis?: any;
  newAnalysis?: any;
}

export interface Primer {
  seq: string;
  len: number;
  tm: number;
  tmTotal: number;
  gc: number;
  dg: number;
  fwd: boolean;
  offTargetCount: number;
  scoring: Scoring;
  penalty: number;
  smartDesign?: SmartDesignMetadata;
  alternatives?: AlternativePair[];
  alternativeCategories?: Record<string, AlternativePair[]>;
  alternativeOptimized?: OptimizationResult;
  dict(): PrimerDict;
}

export interface PrimerDict {
  seq: string;
  len: number;
  tm: number;
  tm_total: number;
  gc: number;
  dg: number;
  fwd: boolean;
  off_target_count: number;
  scoring: {
    penalty: number;
    penalty_tm: number;
    penalty_tm_diff: number;
    penalty_gc: number;
    penalty_len: number;
    penalty_dg: number;
    penalty_off_target: number;
    equilibrium_efficiency: number | null;
    equilibrium_score: number | null;
    equilibrium_losses: EquilibriumLosses | null;
    piecewise_scores: PiecewiseScores | null;
    composite_score: number | null;
    quality_tier: string | null;
  };
}

export interface PrimerFactory {
  optimalTm: number;
  optimalGc: number;
  optimalLen: number;
  penaltyTm: number;
  penaltyTmDiff: number;
  penaltyGc: number;
  penaltyLen: number;
  penaltyDg: number;
  penaltyOffTarget: number;
  build(
    seq: string,
    tm: number,
    tmTotal: number,
    gc: number,
    dg: number,
    fwd: boolean,
    offTargetCount: number
  ): Primer;
  buildPair(fwd: Primer, rev: Primer): [Primer, Primer];
  withOptimalLen(newOptimalLen: number): PrimerFactory;
  addEquilibriumScoring(
    fwd: Primer,
    rev: Primer,
    template: string | null,
    temperature?: number
  ): [Primer, Primer];
  addPiecewiseScoring(
    fwd: Primer,
    rev: Primer | null,
    temperature?: number,
    template?: string | null
  ): [Primer, Primer | null];
}

export interface ScoreOptions {
  optimalTm?: number;
  optimalGc?: number;
  optimalLen?: number;
  penaltyTm?: number;
  penaltyGc?: number;
  penaltyLen?: number;
  penaltyTmDiff?: number;
  penaltyDg?: number;
  penaltyOffTarget?: number;
  includeEquilibrium?: boolean;
  annealingTemperature?: number;
}

export interface PrimersOptions extends ScoreOptions {
  addFwd?: string;
  addRev?: string;
  addFwdLen?: [number, number];
  addRevLen?: [number, number];
  offtargetCheck?: string;
  useCompositeScore?: boolean;
  useSmartDesign?: boolean;
  smartDesignTargetScore?: number;
  allowExtendedPrimers?: boolean;
  exhaustiveSearch?: boolean;
}

export interface AlternativePrimerInfo {
  sequence: string;
  length: number;
  tm: number;
  gc: number;
  dg: number;
  gcPercent: string;
  hasGCClamp: boolean;
}

export interface AlternativePair {
  forward: AlternativePrimerInfo;
  reverse: AlternativePrimerInfo;
  tmDiff: number;
  compositeScore: number;
  qualityTier: string;
  penalty?: number;
  category?: string;
}

export interface OptimizationResult {
  forward: Primer;
  reverse: Primer;
  pairScore: number;
  originalScore?: number;
  optimized: boolean;
  scoreImproved?: boolean;
  improvements?: string[];
  tmDiffPreserved?: boolean;
  scorePreserved?: boolean;
  alternativeOptimized?: {
    forward: Primer;
    reverse: Primer;
    pairScore: number;
    tmDiff: number;
    tradeoff: string;
  };
}

export interface PairCandidate {
  fwd: Primer;
  rev: Primer;
  score: number;
  penalty: number;
  tmDiff?: number;
  compositeScore?: number;
  qualityTier?: string;
  heterodimerDG?: number;
}

export interface PrimerTiers {
  tier1: PairCandidate[];
  tier2: PairCandidate[];
  tier3: PairCandidate[];
  tier4: PairCandidate[];
}

interface ChooseBestOptions {
  useCompositeScore?: boolean;
  annealingTemp?: number;
  template?: string | null;
  exhaustiveSearch?: boolean;
  returnCandidates?: boolean;
}

interface ChooseBestResult {
  best: [Primer, Primer];
  allCandidates: PairCandidate[];
}

/**
 * Create a Scoring object
 */
function createScoring({
  penalty = 0,
  penaltyTm = 0,
  penaltyTmDiff = 0,
  penaltyGc = 0,
  penaltyLen = 0,
  penaltyDg = 0,
  penaltyOffTarget = 0,
  equilibriumEfficiency = null,
  equilibriumScore = null,
  equilibriumLosses = null,
  piecewiseScores = null,
  compositeScore = null,
  qualityTier = null,
}: Partial<Scoring> = {}): Scoring {
  return {
    penalty,
    penaltyTm,
    penaltyTmDiff,
    penaltyGc,
    penaltyLen,
    penaltyDg,
    penaltyOffTarget,
    equilibriumEfficiency,
    equilibriumScore,
    equilibriumLosses,
    piecewiseScores,
    compositeScore,
    qualityTier,
  };
}

/**
 * Create a Primer object
 */
function createPrimer({
  seq,
  len,
  tm,
  tmTotal,
  gc,
  dg,
  fwd,
  offTargetCount,
  scoring,
}: {
  seq: string;
  len: number;
  tm: number;
  tmTotal: number;
  gc: number;
  dg: number;
  fwd: boolean;
  offTargetCount: number;
  scoring: Scoring;
}): Primer {
  return {
    seq,
    len,
    tm,
    tmTotal,
    gc,
    dg,
    fwd,
    offTargetCount,
    scoring,
    get penalty(): number {
      return this.scoring.penalty;
    },
    dict(): PrimerDict {
      return {
        seq: this.seq,
        len: this.len,
        tm: this.tm,
        tm_total: this.tmTotal,
        gc: this.gc,
        dg: this.dg,
        fwd: this.fwd,
        off_target_count: this.offTargetCount,
        scoring: {
          penalty: this.scoring.penalty,
          penalty_tm: this.scoring.penaltyTm,
          penalty_tm_diff: this.scoring.penaltyTmDiff,
          penalty_gc: this.scoring.penaltyGc,
          penalty_len: this.scoring.penaltyLen,
          penalty_dg: this.scoring.penaltyDg,
          penalty_off_target: this.scoring.penaltyOffTarget,
          equilibrium_efficiency: this.scoring.equilibriumEfficiency,
          equilibrium_score: this.scoring.equilibriumScore,
          equilibrium_losses: this.scoring.equilibriumLosses,
          piecewise_scores: this.scoring.piecewiseScores,
          composite_score: this.scoring.compositeScore,
          quality_tier: this.scoring.qualityTier,
        },
      };
    },
  };
}

/**
 * A factory for creating Primers with penalties.
 */
function createPrimerFactory({
  optimalTm = 62.0,
  optimalGc = 0.5,
  optimalLen = 22,
  penaltyTm = 1.0,
  penaltyTmDiff = 1.0,
  penaltyGc = 0.2,
  penaltyLen = 0.5,
  penaltyDg = 2.0,
  penaltyOffTarget = 20.0,
}: Partial<PrimerFactory> = {}): PrimerFactory {
  return {
    optimalTm,
    optimalGc,
    optimalLen,
    penaltyTm,
    penaltyTmDiff,
    penaltyGc,
    penaltyLen,
    penaltyDg,
    penaltyOffTarget,

    /**
     * Create a Primer with a scored penalty.
     */
    build(
      seq: string,
      tm: number,
      tmTotal: number,
      gc: number,
      dg: number,
      fwd: boolean,
      offTargetCount: number
    ): Primer {
      dg = Math.min(dg, 0);
      const penaltyTmVal = Math.abs(tm - this.optimalTm) * this.penaltyTm;
      const penaltyGcVal = Math.abs(gc - this.optimalGc) * this.penaltyGc * 100;
      const penaltyLenVal = Math.abs(seq.length - this.optimalLen) * this.penaltyLen;
      const penaltyDgVal = Math.abs(dg) * this.penaltyDg;
      const penaltyOffTargetVal = offTargetCount * this.penaltyOffTarget;
      const penalty =
        penaltyTmVal + penaltyGcVal + penaltyLenVal + penaltyDgVal + penaltyOffTargetVal;

      return createPrimer({
        seq,
        len: seq.length,
        tm,
        tmTotal,
        gc: Math.round(gc * 100) / 100,
        dg: Math.round(dg * 100) / 100,
        fwd,
        offTargetCount,
        scoring: createScoring({
          penaltyTm: Math.round(penaltyTmVal * 10) / 10,
          penaltyTmDiff: 0, // unknown at this point
          penaltyLen: penaltyLenVal,
          penaltyGc: Math.round(penaltyGcVal * 10) / 10,
          penaltyDg: Math.round(penaltyDgVal * 10) / 10,
          penaltyOffTarget: penaltyOffTargetVal,
          penalty: Math.round(penalty * 10) / 10,
        }),
      });
    },

    /**
     * Create a pair of Primers with a tm_diff penalty added.
     */
    buildPair(fwd: Primer, rev: Primer): [Primer, Primer] {
      const penaltyTmDiffVal = Math.abs(fwd.tm - rev.tm) * this.penaltyTmDiff;

      const newFwd = createPrimer({
        ...fwd,
        scoring: createScoring({
          ...fwd.scoring,
          penalty: Math.round((fwd.scoring.penalty + penaltyTmDiffVal) * 10) / 10,
          penaltyTmDiff: Math.round(penaltyTmDiffVal * 10) / 10,
        }),
      });

      const newRev = createPrimer({
        ...rev,
        scoring: createScoring({
          ...rev.scoring,
          penalty: Math.round((rev.scoring.penalty + penaltyTmDiffVal) * 10) / 10,
          penaltyTmDiff: Math.round(penaltyTmDiffVal * 10) / 10,
        }),
      });

      return [newFwd, newRev];
    },

    /**
     * Create a copy with modified optimal length
     */
    withOptimalLen(newOptimalLen: number): PrimerFactory {
      return createPrimerFactory({
        ...this,
        optimalLen: newOptimalLen,
      });
    },

    /**
     * Add equilibrium efficiency scoring to a primer pair.
     */
    addEquilibriumScoring(
      fwd: Primer,
      rev: Primer,
      template: string | null,
      temperature: number = 55
    ): [Primer, Primer] {
      if (!template) {
        const heterodimerDG = calculateHeterodimerDG(fwd.seq, rev.seq, temperature);

        return [
          createPrimer({
            ...fwd,
            scoring: createScoring({
              ...fwd.scoring,
              equilibriumEfficiency: null,
              equilibriumScore: null,
              equilibriumLosses: { heterodimerDG: Math.round(heterodimerDG * 100) / 100 },
            }),
          }),
          createPrimer({
            ...rev,
            scoring: createScoring({
              ...rev.scoring,
              equilibriumEfficiency: null,
              equilibriumScore: null,
              equilibriumLosses: { heterodimerDG: Math.round(heterodimerDG * 100) / 100 },
            }),
          }),
        ];
      }

      const fwdBindingSite = template.slice(0, fwd.seq.length);
      const revBindingSite = rcSeq(template.slice(-rev.seq.length));

      const fwdInput = {
        seq: fwd.seq,
        bindingSite: fwdBindingSite,
      };
      const revInput = {
        seq: rev.seq,
        bindingSite: revBindingSite,
      };

      const equilibrium = calculateEquilibriumEfficiency(fwdInput, revInput, template, {
        temperature,
        includeOffTarget: true,
      });

      const equilibriumScore = efficiencyToScore(equilibrium.efficiency);

      const newFwd = createPrimer({
        ...fwd,
        scoring: createScoring({
          ...fwd.scoring,
          equilibriumEfficiency: Math.round(equilibrium.efficiencyFwd * 1000) / 1000,
          equilibriumScore: Math.round(equilibriumScore * 10) / 10,
          equilibriumLosses: {
            hairpin: Math.round(equilibrium.losses.fwd.hairpin * 1000) / 1000,
            homodimer: Math.round(equilibrium.losses.fwd.homodimer * 1000) / 1000,
            heterodimer: Math.round(equilibrium.losses.fwd.heterodimer * 1000) / 1000,
            offTarget: Math.round(equilibrium.losses.fwd.offTarget * 1000) / 1000,
            free: Math.round(equilibrium.losses.fwd.free * 1000) / 1000,
          },
        }),
      });

      const newRev = createPrimer({
        ...rev,
        scoring: createScoring({
          ...rev.scoring,
          equilibriumEfficiency: Math.round(equilibrium.efficiencyRev * 1000) / 1000,
          equilibriumScore: Math.round(equilibriumScore * 10) / 10,
          equilibriumLosses: {
            hairpin: Math.round(equilibrium.losses.rev.hairpin * 1000) / 1000,
            homodimer: Math.round(equilibrium.losses.rev.homodimer * 1000) / 1000,
            heterodimer: Math.round(equilibrium.losses.rev.heterodimer * 1000) / 1000,
            offTarget: Math.round(equilibrium.losses.rev.offTarget * 1000) / 1000,
            free: Math.round(equilibrium.losses.rev.free * 1000) / 1000,
          },
        }),
      });

      return [newFwd, newRev];
    },

    /**
     * Add piecewise logistic scoring to a primer pair.
     */
    addPiecewiseScoring(
      fwd: Primer,
      rev: Primer | null = null,
      temperature: number = 55,
      template: string | null = null
    ): [Primer, Primer | null] {
      const preset = ANALYSIS_PRESETS.amplification;

      const hairpinDGFwd = calculateHairpinDG(fwd.seq, temperature);
      const homodimerDGFwd = calculateHomodimerDG(fwd.seq, temperature);

      let offTargetAnalysis: any = null;
      if (template && rev) {
        try {
          offTargetAnalysis = analyzeOffTargetsPair(fwd.seq, rev.seq, template, { temperature });
        } catch (e) {
          offTargetAnalysis = null;
        }
      }

      const fwdG4Analysis = analyzeGQuadruplex(fwd.seq);
      const fwdTerminalDG = calculate3primeTerminalDG(fwd.seq).dG;

      const fwdScores: PiecewiseScores = {
        tm: scoreTm(fwd.tm, preset.tmOptions),
        gc: scoreGc(fwd.gc, preset.gcOptions),
        length: scoreLength(fwd.seq.length, preset.lengthOptions),
        gcClamp: scoreGcClamp(fwd.seq),
        homopolymer: scoreHomopolymer(fwd.seq),
        hairpin: scoreHairpin(hairpinDGFwd, { threshold: preset.hairpinThreshold }),
        homodimer: scoreHomodimer(homodimerDGFwd, { threshold: preset.homodimerThreshold }),
        offTarget: offTargetAnalysis
          ? scoreOffTargetClassification(offTargetAnalysis.fwd)
          : scoreOffTarget(fwd.offTargetCount),
        terminal3DG: scoreTerminal3DG(fwdTerminalDG),
        gQuadruplex: fwdG4Analysis.score,
        threePrimeComp: score3PrimeComposition(fwd.seq, fwdTerminalDG),
      };

      const fwdComposite = calculateCompositeScore({
        tmFwd: fwdScores.tm,
        gcFwd: fwdScores.gc,
        lengthFwd: fwdScores.length,
        gcClampFwd: fwdScores.gcClamp,
        homopolymerFwd: fwdScores.homopolymer,
        hairpinFwd: fwdScores.hairpin,
        selfDimerFwd: fwdScores.homodimer,
        offTarget: fwdScores.offTarget,
        terminal3DG: fwdScores.terminal3DG,
        gQuadruplexFwd: fwdScores.gQuadruplex,
        threePrimeCompFwd: fwdScores.threePrimeComp,
      });

      const newFwd = createPrimer({
        ...fwd,
        scoring: createScoring({
          ...fwd.scoring,
          piecewiseScores: fwdScores,
          compositeScore: fwdComposite.score,
          qualityTier: classifyQuality(fwdComposite.score).tier,
          gQuadruplex: fwdG4Analysis,
        }),
      });

      if (!rev) {
        return [newFwd, null];
      }

      const hairpinDGRev = calculateHairpinDG(rev.seq, temperature);
      const homodimerDGRev = calculateHomodimerDG(rev.seq, temperature);
      const heterodimerDG = calculateHeterodimerDG(fwd.seq, rev.seq, temperature);

      const revG4Analysis = analyzeGQuadruplex(rev.seq);
      const revTerminalDG = calculate3primeTerminalDG(rev.seq).dG;

      const revScores: PiecewiseScores = {
        tm: scoreTm(rev.tm, preset.tmOptions),
        gc: scoreGc(rev.gc, preset.gcOptions),
        length: scoreLength(rev.seq.length, preset.lengthOptions),
        gcClamp: scoreGcClamp(rev.seq),
        homopolymer: scoreHomopolymer(rev.seq),
        hairpin: scoreHairpin(hairpinDGRev, { threshold: preset.hairpinThreshold }),
        homodimer: scoreHomodimer(homodimerDGRev, { threshold: preset.homodimerThreshold }),
        offTarget: offTargetAnalysis
          ? scoreOffTargetClassification(offTargetAnalysis.rev)
          : scoreOffTarget(rev.offTargetCount),
        terminal3DG: scoreTerminal3DG(revTerminalDG),
        heterodimer: scoreHeterodimer(heterodimerDG, { threshold: preset.heterodimerThreshold }),
        tmDiff: scoreTmDiff(fwd.tm, rev.tm),
        gQuadruplex: revG4Analysis.score,
        threePrimeComp: score3PrimeComposition(rev.seq, revTerminalDG),
      };

      fwdScores.heterodimer = scoreHeterodimer(heterodimerDG, { threshold: preset.heterodimerThreshold });
      fwdScores.tmDiff = scoreTmDiff(fwd.tm, rev.tm);

      const pairComposite = calculateCompositeScore({
        tmFwd: fwdScores.tm,
        tmRev: revScores.tm,
        gcFwd: fwdScores.gc,
        gcRev: revScores.gc,
        lengthFwd: fwdScores.length,
        lengthRev: revScores.length,
        gcClampFwd: fwdScores.gcClamp,
        gcClampRev: revScores.gcClamp,
        homopolymerFwd: fwdScores.homopolymer,
        homopolymerRev: revScores.homopolymer,
        hairpinFwd: fwdScores.hairpin,
        hairpinRev: revScores.hairpin,
        selfDimerFwd: fwdScores.homodimer,
        selfDimerRev: revScores.homodimer,
        heterodimer: revScores.heterodimer,
        tmDiff: revScores.tmDiff,
        offTarget: Math.min(fwdScores.offTarget!, revScores.offTarget!),
        terminal3DG: Math.min(fwdScores.terminal3DG!, revScores.terminal3DG!),
        gQuadruplexFwd: fwdScores.gQuadruplex,
        gQuadruplexRev: revScores.gQuadruplex,
        threePrimeCompFwd: fwdScores.threePrimeComp,
        threePrimeCompRev: revScores.threePrimeComp,
      });

      const criticalWarnings: string[] = [];
      const lengthLow = preset.lengthOptions?.optimalLow ?? 18;
      const lengthHigh = preset.lengthOptions?.optimalHigh ?? 24;

      if (fwd.seq.length > lengthHigh + 15) {
        criticalWarnings.push(`Forward: Primer extremely long (${fwd.seq.length}bp)`);
      }
      let fwdMaxHomopolymer = 1, fwdCurrentRun = 1;
      for (let i = 1; i < fwd.seq.length; i++) {
        if (fwd.seq[i].toUpperCase() === fwd.seq[i - 1].toUpperCase()) {
          fwdCurrentRun++;
          fwdMaxHomopolymer = Math.max(fwdMaxHomopolymer, fwdCurrentRun);
        } else {
          fwdCurrentRun = 1;
        }
      }
      if (fwdMaxHomopolymer >= 6) {
        criticalWarnings.push(`Forward: Severe homopolymer run (${fwdMaxHomopolymer} bases)`);
      }

      if (rev.seq.length > lengthHigh + 15) {
        criticalWarnings.push(`Reverse: Primer extremely long (${rev.seq.length}bp)`);
      }
      let revMaxHomopolymer = 1, revCurrentRun = 1;
      for (let i = 1; i < rev.seq.length; i++) {
        if (rev.seq[i].toUpperCase() === rev.seq[i - 1].toUpperCase()) {
          revCurrentRun++;
          revMaxHomopolymer = Math.max(revMaxHomopolymer, revCurrentRun);
        } else {
          revCurrentRun = 1;
        }
      }
      if (revMaxHomopolymer >= 6) {
        criticalWarnings.push(`Reverse: Severe homopolymer run (${revMaxHomopolymer} bases)`);
      }

      const tmDiff = Math.abs(fwd.tm - rev.tm);
      if (tmDiff > 8) {
        criticalWarnings.push(`Tm difference: ${tmDiff.toFixed(1)}°C`);
      }

      const effectiveScore = Math.max(0, pairComposite.score - criticalWarnings.length * 20);
      const quality = classifyQuality(effectiveScore);

      const finalFwd = createPrimer({
        ...fwd,
        scoring: createScoring({
          ...fwd.scoring,
          piecewiseScores: fwdScores,
          compositeScore: pairComposite.score,
          effectiveScore,
          criticalWarnings: criticalWarnings.length,
          qualityTier: quality.tier,
          gQuadruplex: fwdG4Analysis,
        }),
      });

      const newRev = createPrimer({
        ...rev,
        scoring: createScoring({
          ...rev.scoring,
          piecewiseScores: revScores,
          compositeScore: pairComposite.score,
          effectiveScore,
          criticalWarnings: criticalWarnings.length,
          qualityTier: quality.tier,
          gQuadruplex: revG4Analysis,
        }),
      });

      return [finalFwd, newRev];
    },
  };
}

/**
 * Score primers from their sequence.
 */
export function score(
  fwd: string,
  rev: string = "",
  seq: string = "",
  offtargetCheck: string = "",
  {
    optimalTm = 62.0,
    optimalGc = 0.5,
    optimalLen = 22,
    penaltyTm = 1.0,
    penaltyGc = 0.2,
    penaltyLen = 0.5,
    penaltyTmDiff = 1.0,
    penaltyDg = 2.0,
    penaltyOffTarget = 20.0,
    includeEquilibrium = true,
    annealingTemperature = 55,
  }: ScoreOptions = {}
): [Primer, Primer | null] {
  if (fwd.length < LEN_MIN) {
    throw new Error(`\`fwd\` primer is too short: ${fwd.length} < ${LEN_MIN}`);
  }
  if (rev && rev.length < LEN_MIN) {
    throw new Error(`\`rev\` primer is too short: ${rev.length} < ${LEN_MIN}`);
  }

  fwd = fwd.toUpperCase();
  rev = rev.toUpperCase();
  const [parsedSeq, parsedOfftargetCheck] = parse(seq, offtargetCheck);
  const [addFwd, , addRev] = bindingSeq(fwd, rev, parsedSeq);

  const factory = createPrimerFactory({
    optimalTm,
    optimalGc,
    optimalLen,
    penaltyTm,
    penaltyGc,
    penaltyLen,
    penaltyTmDiff,
    penaltyDg,
    penaltyOffTarget,
  });

  const tmCacheFwd = tmCache(fwd);
  const dgCacheFwd = dgCache(fwd);
  const gcCacheFwd = gcCache(fwd);
  const otFwd = offTargets(fwd, parsedOfftargetCheck);

  const fwdPrimer = factory.build(
    fwd,
    tmCacheFwd[addFwd][fwd.length - 1],
    tmCacheFwd[0][fwd.length - 1],
    gcCacheFwd[0][fwd.length - 1],
    dgCacheFwd[0][fwd.length - 1],
    true,
    otFwd[fwd.length - 1]
  );

  if (!rev) {
    const [fwdWithScoring] = factory.addPiecewiseScoring(fwdPrimer, null, annealingTemperature);
    return [fwdWithScoring, null];
  }

  const tmCacheRev = tmCache(rev);
  const dgCacheRev = dgCache(rev);
  const gcCacheRev = gcCache(rev);
  const otRev = offTargets(rev, parsedOfftargetCheck);

  const revPrimer = factory.build(
    rev,
    tmCacheRev[addRev][rev.length - 1],
    tmCacheRev[0][rev.length - 1],
    gcCacheRev[0][rev.length - 1],
    dgCacheRev[0][rev.length - 1],
    false,
    otRev[rev.length - 1]
  );

  const [fwdWithTmDiff, revWithTmDiff] = factory.buildPair(fwdPrimer, revPrimer);

  let resultFwd = fwdWithTmDiff;
  let resultRev = revWithTmDiff;

  if (includeEquilibrium) {
    [resultFwd, resultRev] = factory.addEquilibriumScoring(
      resultFwd,
      resultRev,
      parsedSeq,
      annealingTemperature
    );
  }

  [resultFwd, resultRev] = factory.addPiecewiseScoring(
    resultFwd,
    resultRev,
    annealingTemperature,
    parsedOfftargetCheck || parsedSeq
  ) as [Primer, Primer];

  return [resultFwd, resultRev];
}

/**
 * Attempt to find the binding region between the fwd and rev primers in seq
 */
function bindingSeq(fwd: string, rev: string = "", seq: string = ""): [number, string, number] {
  if (!seq) {
    return [0, "", 0];
  }

  fwd = fwd.toUpperCase();
  rev = rev.toUpperCase();
  seq = seq.toUpperCase();
  seq = seq + seq + seq;

  let addFwd = 0;
  let addRev = 0;

  try {
    let startIndexFwd = -10;
    let startIndexSeq = seq.indexOf(fwd.slice(startIndexFwd));

    if (startIndexSeq === -1) {
      throw new Error("Not found");
    }

    while (
      startIndexSeq &&
      fwd.length + startIndexFwd > 0 &&
      seq[startIndexSeq] === fwd[fwd.length + startIndexFwd]
    ) {
      startIndexFwd -= 1;
      startIndexSeq -= 1;
    }

    addFwd = fwd.length + startIndexFwd;
    seq = seq.slice(startIndexSeq);
  } catch (err: any) {
    throw new Error(`failed to find \`fwd\` binding site in \`seq\`: ${err.message}`);
  }

  if (!rev) {
    return [addFwd, seq, addRev];
  }

  try {
    let endIndexRev = 10;
    rev = reverseComplement(rev);
    let endIndexSeq = seq.indexOf(rev.slice(0, endIndexRev));

    if (endIndexSeq === -1) {
      throw new Error("Not found");
    }

    while (
      endIndexRev < rev.length - 1 &&
      seq.slice(endIndexSeq).includes(rev.slice(0, endIndexRev + 2))
    ) {
      endIndexRev += 1;
    }

    addRev = rev.length - endIndexRev - 1;
    seq = seq.slice(0, endIndexSeq + endIndexRev + 1);
  } catch (err: any) {
    throw new Error(`failed to find \`rev\` binding site in \`seq\`: ${err.message}`);
  }

  return [addFwd, seq, addRev];
}

/**
 * Create primers for PCR amplification of the sequence.
 */
export function primers(
  seq: string,
  {
    addFwd = "",
    addRev = "",
    addFwdLen = [-1, -1] as [number, number],
    addRevLen = [-1, -1] as [number, number],
    offtargetCheck = "",
    optimalTm = 62.0,
    optimalGc = 0.5,
    optimalLen = 22,
    penaltyTm = 1.0,
    penaltyGc = 0.2,
    penaltyLen = 0.5,
    penaltyTmDiff = 1.0,
    penaltyDg = 2.0,
    penaltyOffTarget = 20.0,
    useCompositeScore = false,
    useSmartDesign = false,
    smartDesignTargetScore = 75,
    allowExtendedPrimers = true,
    exhaustiveSearch = true,
  }: PrimersOptions = {}
): [Primer, Primer] {
  const [parsedSeq, parsedOfftargetCheck] = parse(seq, offtargetCheck);
  const [parsedAddFwd] = parse(addFwd, "");
  const [parsedAddRev] = parse(addRev, "");

  const factory = createPrimerFactory({
    optimalTm,
    optimalGc,
    optimalLen,
    penaltyTm,
    penaltyGc,
    penaltyLen,
    penaltyTmDiff,
    penaltyDg,
    penaltyOffTarget,
  });

  const [addFwdMin, addFwdMax] = parseAddLen(parsedAddFwd, addFwdLen);
  const [addRevMin, addRevMax] = parseAddLen(parsedAddRev, addRevLen);

  const effectiveMaxLen = allowExtendedPrimers ? LEN_MAX_EXTENDED : LEN_MAX;

  const trimmedAddFwd = parsedAddFwd.slice(-addFwdMax);
  const trimmedAddRev = reverseComplement(parsedAddRev).slice(0, addRevMax);
  const seqFull = trimmedAddFwd + parsedSeq + trimmedAddRev;

  if (seqFull.length < LEN_MAX) {
    throw new Error(
      `Template sequence length is too short: ${seqFull.length}bp < ${LEN_MAX}bp`
    );
  }

  const optimalFwdLen = Math.round(optimalLen + (addFwdMin + addFwdMax) / 2);
  const fwdSeq = seqFull.slice(0, addFwdMax + effectiveMaxLen);
  const fwdPrimers = createPrimers(
    factory.withOptimalLen(optimalFwdLen),
    fwdSeq,
    parsedOfftargetCheck,
    range(0, addFwdMax - addFwdMin + 1),
    range(LEN_MIN - 1 + addFwdMax, effectiveMaxLen + addFwdMax),
    true,
    addFwdMax
  );

  const optimalRevLen = Math.round(optimalLen + (addRevMin + addRevMax) / 2);
  let revSeq = reverseComplement(seqFull);
  revSeq = revSeq.slice(0, addRevMax + effectiveMaxLen);
  const revPrimers = createPrimers(
    factory.withOptimalLen(optimalRevLen),
    revSeq,
    parsedOfftargetCheck,
    range(0, addRevMax - addRevMin + 1),
    range(LEN_MIN - 1 + addRevMax, effectiveMaxLen + addRevMin),
    false,
    addRevMax
  );

  const template = parsedOfftargetCheck || parsedSeq;
  const chooseBestResult = chooseBest(factory, fwdPrimers, revPrimers, {
    useCompositeScore,
    annealingTemp: optimalTm,
    template,
    exhaustiveSearch,
    returnCandidates: exhaustiveSearch,
  });

  let minFwd: Primer, minRev: Primer, allCandidates: PairCandidate[] | null;
  if (Array.isArray(chooseBestResult)) {
    [minFwd, minRev] = chooseBestResult;
    allCandidates = null;
  } else {
    [minFwd, minRev] = chooseBestResult.best;
    allCandidates = chooseBestResult.allCandidates;
  }

  let [fwdWithScoring, revWithScoring] = factory.addPiecewiseScoring(
    minFwd,
    minRev,
    optimalTm,
    template
  );

  if (useSmartDesign) {
    const smartResult = applySmartDesign(
      fwdWithScoring,
      revWithScoring!,
      seqFull,
      template,
      {
        targetScore: smartDesignTargetScore,
        temperature: optimalTm,
        addFwdMax,
        addRevMax,
      }
    );

    fwdWithScoring = smartResult.forward;
    revWithScoring = smartResult.reverse;

    if (smartResult.alternativeOptimized) {
      (fwdWithScoring as any).alternativeOptimized = smartResult.alternativeOptimized;
    }
  }

  try {
    if (allCandidates && allCandidates.length > 0) {
      const currentBest = { fwd: fwdWithScoring, rev: revWithScoring! };
      const diverseResult = generateDiverseAlternatives(allCandidates, currentBest);

      if (diverseResult.alternatives.length > 0) {
        fwdWithScoring.alternatives = diverseResult.alternatives;
        fwdWithScoring.alternativeCategories = diverseResult.categories;
        console.log(`[primers] Generated ${diverseResult.alternatives.length} diverse alternatives from exhaustive search`);
      }
    } else {
      const alternatives = generateAlternativesInternal(
        fwdPrimers,
        revPrimers,
        factory,
        template,
        optimalTm,
        { numAlternatives: 5, currentFwd: fwdWithScoring, currentRev: revWithScoring! }
      );

      if (alternatives.length > 0) {
        fwdWithScoring.alternatives = alternatives;
        console.log(`[primers] Generated ${alternatives.length} alternative pairs`);
      }
    }
  } catch (err: any) {
    console.warn(`[primers] Failed to generate alternatives: ${err.message}`);
  }

  return [fwdWithScoring, revWithScoring!];
}

/**
 * Apply smart design optimization to improve primer pair quality.
 */
function applySmartDesign(
  fwd: Primer,
  rev: Primer,
  fullSeq: string,
  template: string,
  options: {
    targetScore?: number;
    temperature?: number;
    maxIterations?: number;
    addFwdMax?: number;
    addRevMax?: number;
    verbose?: boolean;
  } = {}
): OptimizationResult {
  const {
    targetScore = 75,
    temperature = 55,
    maxIterations = 3,
    addFwdMax = 0,
    addRevMax = 0,
    verbose = true,
  } = options;

  const currentScore = fwd.scoring.compositeScore!;

  if (verbose) {
    console.log('\n[Smart Design] Analyzing primer pair...');
    console.log(`  Current score: ${currentScore}/100 (target: ${targetScore})`);
  }

  const fwdAnalysis = analyze3PrimeEnd(fwd.seq);
  const revAnalysis = analyze3PrimeEnd(rev.seq);

  const fwd3PrimeNeedsWork = fwdAnalysis.quality === 'poor' || fwdAnalysis.quality === 'marginal' || fwdAnalysis.gcCounts.last2 === 0;
  const rev3PrimeNeedsWork = revAnalysis.quality === 'poor' || revAnalysis.quality === 'marginal' || revAnalysis.gcCounts.last2 === 0;
  const has3PrimeIssues = fwd3PrimeNeedsWork || rev3PrimeNeedsWork;

  if (currentScore >= targetScore && !has3PrimeIssues) {
    if (verbose) {
      console.log('  ✓ Already meets target score with good 3\' ends - no optimization needed');
      console.log('  Forward primer 3\' end:');
      console.log(`    ...${fwd.seq.slice(-10)} | GC clamp: ${fwdAnalysis.gcCounts.last2}/2, Quality: ${fwdAnalysis.quality}`);
      console.log('  Reverse primer 3\' end:');
      console.log(`    ...${rev.seq.slice(-10)} | GC clamp: ${revAnalysis.gcCounts.last2}/2, Quality: ${revAnalysis.quality}\n`);
    }
    return {
      forward: addSmartDesignMetadata(fwd, { optimized: false, reason: 'Already meets target score with good 3\' ends' }),
      reverse: addSmartDesignMetadata(rev, { optimized: false, reason: 'Already meets target score with good 3\' ends' }),
      pairScore: currentScore,
      optimized: false,
    };
  }

  if (verbose && currentScore >= targetScore && has3PrimeIssues) {
    console.log(`  ⚠ Score ${currentScore} meets target, but 3' end quality needs improvement`);
  }

  const fwdNeedsOptim = fwdAnalysis.quality !== 'excellent' || fwdAnalysis.gcCounts.last2 === 0;
  const revNeedsOptim = revAnalysis.quality !== 'excellent' || revAnalysis.gcCounts.last2 === 0;

  if (verbose) {
    console.log('  Forward primer analysis:');
    console.log(`    Sequence: ${fwd.seq.slice(-10)}... (last 10bp)`);
    console.log(`    GC clamp: ${fwdAnalysis.gcCounts.last2}/2, Quality: ${fwdAnalysis.quality}`);
    console.log(`    Needs optimization: ${fwdNeedsOptim ? 'Yes' : 'No'}`);
    console.log('  Reverse primer analysis:');
    console.log(`    Sequence: ${rev.seq.slice(-10)}... (last 10bp)`);
    console.log(`    GC clamp: ${revAnalysis.gcCounts.last2}/2, Quality: ${revAnalysis.quality}`);
    console.log(`    Needs optimization: ${revNeedsOptim ? 'Yes' : 'No'}`);
  }

  let optimizedFwd = fwd;
  let optimizedRev = rev;
  const improvements: string[] = [];

  if (fwdNeedsOptim) {
    if (verbose) console.log('  Optimizing forward primer...');
    const fwdStartPos = addFwdMax;
    const fwdResult = optimizePrimer(fullSeq, fwdStartPos, fwd.seq.length, true, {
      targetGcClamp: true,
      target3PrimeDG: true,
      maxLengthChange: 3,
    });

    if (fwdResult.optimized && fwdResult.seq !== fwd.seq) {
      optimizedFwd = createOptimizedPrimer(fwdResult, fwd, true, fullSeq, template, temperature);
      improvements.push(`FWD: ${fwdResult.reason}`);
      if (verbose) {
        console.log(`    ✓ Optimized: ${fwd.seq.length}bp → ${fwdResult.seq.length}bp`);
        console.log(`    Reason: ${fwdResult.reason}`);
      }
    } else if (verbose) {
      console.log('    No improvement found');
    }
  }

  if (revNeedsOptim) {
    if (verbose) console.log('  Optimizing reverse primer...');
    const revStartPos = fullSeq.length - addRevMax - rev.seq.length;
    const revFullSeq = reverseComplement(fullSeq);
    const revResult = optimizePrimer(revFullSeq, addRevMax, rev.seq.length, false, {
      targetGcClamp: true,
      target3PrimeDG: true,
      maxLengthChange: 3,
    });

    if (revResult.optimized && revResult.seq !== rev.seq) {
      optimizedRev = createOptimizedPrimer(revResult, rev, false, fullSeq, template, temperature);
      improvements.push(`REV: ${revResult.reason}`);
      if (verbose) {
        console.log(`    ✓ Optimized: ${rev.seq.length}bp → ${revResult.seq.length}bp`);
        console.log(`    Reason: ${revResult.reason}`);
      }
    } else if (verbose) {
      console.log('    No improvement found');
    }
  }

  const originalFwdG4 = analyzeGQuadruplex(fwd.seq);
  const originalRevG4 = analyzeGQuadruplex(rev.seq);
  const optimizedFwdG4 = analyzeGQuadruplex(optimizedFwd.seq);
  const optimizedRevG4 = analyzeGQuadruplex(optimizedRev.seq);

  const fwdG4Worsened = !originalFwdG4.hasG4Motif && optimizedFwdG4.hasG4Motif;
  const revG4Worsened = !originalRevG4.hasG4Motif && optimizedRevG4.hasG4Motif;
  const g4Worsened = fwdG4Worsened || revG4Worsened;

  const originalFwdHairpin = calculateHairpinDG(fwd.seq, temperature);
  const originalRevHairpin = calculateHairpinDG(rev.seq, temperature);
  const optimizedFwdHairpin = calculateHairpinDG(optimizedFwd.seq, temperature);
  const optimizedRevHairpin = calculateHairpinDG(optimizedRev.seq, temperature);

  const hairpinWorsened = (
    (optimizedFwdHairpin < originalFwdHairpin - 3) ||
    (optimizedRevHairpin < originalRevHairpin - 3) ||
    (optimizedFwdHairpin < -6 && originalFwdHairpin >= -6) ||
    (optimizedRevHairpin < -6 && originalRevHairpin >= -6)
  );

  if (g4Worsened || hairpinWorsened) {
    if (verbose) {
      if (g4Worsened) {
        console.log('  ⚠ Rejecting optimization: Created G-quadruplex motif');
      }
      if (hairpinWorsened) {
        console.log('  ⚠ Rejecting optimization: Hairpin stability worsened significantly');
      }
    }
    optimizedFwd = fwd;
    optimizedRev = rev;
    improvements.length = 0;
  }

  const newPairScores = calculatePairCompositeScore(optimizedFwd, optimizedRev, temperature);
  const newScore = newPairScores.score;
  const scoreImproved = newScore > currentScore;
  const scoreWorsened = newScore < currentScore;

  const quality = classifyQuality(newScore);
  optimizedFwd.scoring.compositeScore = newScore;
  optimizedFwd.scoring.qualityTier = quality.tier;
  optimizedRev.scoring.compositeScore = newScore;
  optimizedRev.scoring.qualityTier = quality.tier;

  const originalTmDiff = Math.abs(fwd.tm - rev.tm);
  const newTmDiff = Math.abs(optimizedFwd.tm - optimizedRev.tm);
  const tmDiffWorsened = newTmDiff > originalTmDiff + 2;

  if (verbose) {
    console.log(`  Tm diff: ${originalTmDiff.toFixed(1)}°C → ${newTmDiff.toFixed(1)}°C`);
  }

  if ((tmDiffWorsened || scoreWorsened) && improvements.length > 0) {
    const revertReason = scoreWorsened
      ? `Score dropped from ${currentScore} to ${newScore}`
      : `Tm diff worsened from ${originalTmDiff.toFixed(1)}°C to ${newTmDiff.toFixed(1)}°C`;
    if (verbose) {
      console.log(`  ⚠ Reverting optimizations: ${revertReason}`);
      if (scoreWorsened) {
        console.log('  Score preservation takes priority over 3\' end optimization');
      } else {
        console.log('  Tm matching takes priority over GC clamp optimization');
      }
      console.log(`  Original primers retained (score=${currentScore})\n`);
    }
    return {
      forward: addSmartDesignMetadata(fwd, { optimized: false, reason: `Reverted - ${scoreWorsened ? 'Score drop' : 'Tm matching priority'}` }),
      reverse: addSmartDesignMetadata(rev, { optimized: false, reason: `Reverted - ${scoreWorsened ? 'Score drop' : 'Tm matching priority'}` }),
      pairScore: currentScore,
      originalScore: currentScore,
      optimized: false,
      scoreImproved: false,
      improvements: [],
      tmDiffPreserved: !tmDiffWorsened,
      scorePreserved: !scoreWorsened,
      alternativeOptimized: {
        forward: addSmartDesignMetadata(optimizedFwd, {
          optimized: true,
          reason: scoreWorsened ? 'Lower score, better 3\' ends' : 'Higher score, worse Tm match',
        }),
        reverse: addSmartDesignMetadata(optimizedRev, {
          optimized: true,
          reason: scoreWorsened ? 'Lower score, better 3\' ends' : 'Higher score, worse Tm match',
        }),
        pairScore: newScore,
        tmDiff: newTmDiff,
        tradeoff: scoreWorsened
          ? `Score ${newScore} (was ${currentScore}) with improved 3' ends`
          : `Score ${newScore} (+${newScore - currentScore}) but Tm diff ${newTmDiff.toFixed(1)}°C (was ${originalTmDiff.toFixed(1)}°C)`,
      },
    };
  }

  if (verbose) {
    if (improvements.length > 0) {
      console.log(`  Score: ${currentScore} → ${newScore} (${scoreImproved ? '↑' : '→'})`);
      console.log(`  Improvements: ${improvements.join(', ')}\n`);
    } else {
      console.log('  No optimizations applied\n');
    }
  }

  return {
    forward: addSmartDesignMetadata(optimizedFwd, {
      optimized: fwdNeedsOptim && improvements.some(i => i.startsWith('FWD')),
      originalAnalysis: fwdAnalysis,
      newAnalysis: analyze3PrimeEnd(optimizedFwd.seq),
    }),
    reverse: addSmartDesignMetadata(optimizedRev, {
      optimized: revNeedsOptim && improvements.some(i => i.startsWith('REV')),
      originalAnalysis: revAnalysis,
      newAnalysis: analyze3PrimeEnd(optimizedRev.seq),
    }),
    pairScore: newScore,
    originalScore: currentScore,
    optimized: improvements.length > 0,
    scoreImproved,
    improvements,
  };
}

/**
 * Calculate the proper PAIR composite score for a primer pair
 */
function calculatePairCompositeScore(
  fwd: Primer,
  rev: Primer,
  temperature: number = 55
): { score: number; breakdown: PiecewiseScores } {
  const sanitize = (val: number | undefined, fallback: number = 0): number => {
    return Number.isFinite(val) ? val! : fallback;
  };

  const heterodimerDG = sanitize(calculateHeterodimerDG(fwd.seq, rev.seq, temperature), 0);
  const hairpinDGFwd = sanitize(calculateHairpinDG(fwd.seq, temperature), 0);
  const hairpinDGRev = sanitize(calculateHairpinDG(rev.seq, temperature), 0);
  const homodimerDGFwd = sanitize(calculateHomodimerDG(fwd.seq, temperature), 0);
  const homodimerDGRev = sanitize(calculateHomodimerDG(rev.seq, temperature), 0);

  const fwdG4 = analyzeGQuadruplex(fwd.seq);
  const revG4 = analyzeGQuadruplex(rev.seq);

  let fwdTerminalDG = typeof fwd.dg === 'number' ? fwd.dg : null;
  let revTerminalDG = typeof rev.dg === 'number' ? rev.dg : null;

  if (fwdTerminalDG === null) {
    const calc = calculate3primeTerminalDG(fwd.seq);
    fwdTerminalDG = calc?.dG ?? -8;
  }
  if (revTerminalDG === null) {
    const calc = calculate3primeTerminalDG(rev.seq);
    revTerminalDG = calc?.dG ?? -8;
  }

  const threePrimeCompFwd = score3PrimeComposition(fwd.seq, fwdTerminalDG);
  const threePrimeCompRev = score3PrimeComposition(rev.seq, revTerminalDG);

  const pairScores: PiecewiseScores = {
    tm: sanitize(scoreTm(fwd.tm), 0.5),
    gc: sanitize(scoreGc(fwd.gc), 0.5),
    length: sanitize(scoreLength(fwd.seq.length), 0.5),
    gcClamp: sanitize(scoreGcClamp(fwd.seq), 0.5),
    homopolymer: sanitize(scoreHomopolymer(fwd.seq), 0.5),
    hairpin: sanitize(scoreHairpin(hairpinDGFwd), 0.5),
    homodimer: sanitize(scoreHomodimer(homodimerDGFwd), 0.5),
    heterodimer: sanitize(scoreHeterodimer(heterodimerDG), 0.5),
    tmDiff: sanitize(scoreTmDiff(fwd.tm, rev.tm), 0.5),
    offTarget: Math.min(
      sanitize(scoreOffTarget(fwd.offTargetCount || 0), 1),
      sanitize(scoreOffTarget(rev.offTargetCount || 0), 1)
    ),
    terminal3DG: Math.min(
      sanitize(scoreTerminal3DG(fwdTerminalDG), 0.5),
      sanitize(scoreTerminal3DG(revTerminalDG), 0.5)
    ),
    gQuadruplex: sanitize(fwdG4.score, 1),
    threePrimeComp: sanitize(threePrimeCompFwd, 0.5),
  };

  const result = calculateCompositeScore({
    tmFwd: pairScores.tm,
    tmRev: sanitize(scoreTm(rev.tm), 0.5),
    gcFwd: pairScores.gc,
    gcRev: sanitize(scoreGc(rev.gc), 0.5),
    lengthFwd: pairScores.length,
    lengthRev: sanitize(scoreLength(rev.seq.length), 0.5),
    gcClampFwd: pairScores.gcClamp,
    gcClampRev: sanitize(scoreGcClamp(rev.seq), 0.5),
    homopolymerFwd: pairScores.homopolymer,
    homopolymerRev: sanitize(scoreHomopolymer(rev.seq), 0.5),
    hairpinFwd: pairScores.hairpin,
    hairpinRev: sanitize(scoreHairpin(hairpinDGRev), 0.5),
    selfDimerFwd: pairScores.homodimer,
    selfDimerRev: sanitize(scoreHomodimer(homodimerDGRev), 0.5),
    heterodimer: pairScores.heterodimer,
    tmDiff: pairScores.tmDiff,
    offTarget: pairScores.offTarget,
    terminal3DG: pairScores.terminal3DG,
    gQuadruplexFwd: pairScores.gQuadruplex,
    gQuadruplexRev: sanitize(revG4.score, 1),
    threePrimeCompFwd: pairScores.threePrimeComp,
    threePrimeCompRev: sanitize(threePrimeCompRev, 0.5),
  });

  return {
    score: result.score,
    breakdown: pairScores,
  };
}

/**
 * Create an optimized primer object from optimization result
 */
function createOptimizedPrimer(
  optimResult: any,
  originalPrimer: Primer,
  isFwd: boolean,
  fullSeq: string,
  template: string,
  temperature: number
): Primer {
  const newScoring = scorePrimerVariant(optimResult.seq, template, isFwd, temperature);

  return createPrimer({
    seq: optimResult.seq,
    len: optimResult.seq.length,
    tm: newScoring.tm,
    tmTotal: newScoring.tm,
    gc: newScoring.gc,
    dg: newScoring.terminalDG,
    fwd: isFwd,
    offTargetCount: originalPrimer.offTargetCount,
    scoring: createScoring({
      ...originalPrimer.scoring,
      piecewiseScores: newScoring.scores,
      compositeScore: newScoring.compositeScore,
      qualityTier: newScoring.qualityTier,
    }),
  });
}

/**
 * Add smart design metadata to a primer
 */
function addSmartDesignMetadata(primer: Primer, metadata: Partial<SmartDesignMetadata>): Primer {
  return {
    ...primer,
    smartDesign: {
      optimized: metadata.optimized || false,
      reason: metadata.reason || (metadata.optimized ? 'Improved 3\' end composition' : 'No optimization needed'),
      originalAnalysis: metadata.originalAnalysis || null,
      newAnalysis: metadata.newAnalysis || null,
    },
  };
}

/**
 * Pair-Aware Optimization Algorithm
 */
export function optimizePairHolistically(
  seq: string,
  currentFwd: Primer,
  currentRev: Primer,
  options: {
    maxLengthDelta?: number;
    maxSinglePrimerDrop?: number;
    verbose?: boolean;
  } = {}
): any {
  const {
    maxLengthDelta = 3,
    maxSinglePrimerDrop = 5,
    verbose = false,
  } = options;

  const currentPairScore = (currentFwd.scoring?.compositeScore || 70) +
                           (currentRev.scoring?.compositeScore || 70);
  const currentTmDiff = Math.abs(currentFwd.tm - currentRev.tm);

  if (verbose) {
    console.log('\n[Pair-Aware Optimization]');
    console.log(`  Current pair score: ${currentPairScore} (Fwd: ${currentFwd.scoring?.compositeScore}, Rev: ${currentRev.scoring?.compositeScore})`);
    console.log(`  Current Tm diff: ${currentTmDiff.toFixed(1)}°C`);
  }

  const fwdLengths: number[] = [];
  const revLengths: number[] = [];

  for (let delta = -maxLengthDelta; delta <= maxLengthDelta; delta++) {
    const fwdLen = currentFwd.seq.length + delta;
    const revLen = currentRev.seq.length + delta;
    if (fwdLen >= LEN_MIN && fwdLen <= LEN_MAX) fwdLengths.push(fwdLen);
    if (revLen >= LEN_MIN && revLen <= LEN_MAX) revLengths.push(revLen);
  }

  const candidates: any[] = [];

  for (const fwdLen of fwdLengths) {
    for (const revLen of revLengths) {
      const fwdSeq = seq.slice(0, fwdLen);
      const revSeq = reverseComplement(seq).slice(0, revLen);

      const fwdTm = calcTm(fwdSeq);
      const revTm = calcTm(revSeq);
      const tmDiff = Math.abs(fwdTm - revTm);

      const fwd3Prime = analyze3PrimeEnd(fwdSeq);
      const rev3Prime = analyze3PrimeEnd(revSeq);

      const tmDiffPenalty = tmDiff > 5 ? (tmDiff - 5) * 2 : 0;
      const qualityBalance = Math.abs(
        (fwd3Prime.quality === 'excellent' ? 3 : fwd3Prime.quality === 'good' ? 2 : 1) -
        (rev3Prime.quality === 'excellent' ? 3 : rev3Prime.quality === 'good' ? 2 : 1)
      );

      const estimatedPairScore = 140 - tmDiffPenalty - qualityBalance * 3;

      candidates.push({
        fwdLen,
        revLen,
        tmDiff,
        fwd3PrimeQuality: fwd3Prime.quality,
        rev3PrimeQuality: rev3Prime.quality,
        estimatedPairScore,
        balanceScore: 4 - qualityBalance,
      });
    }
  }

  candidates.sort((a, b) => b.estimatedPairScore - a.estimatedPairScore);

  const paretoOptimal = candidates.filter((c, idx) => {
    return !candidates.some((other, otherIdx) =>
      otherIdx !== idx &&
      other.estimatedPairScore >= c.estimatedPairScore &&
      other.tmDiff <= c.tmDiff &&
      other.balanceScore >= c.balanceScore &&
      (other.estimatedPairScore > c.estimatedPairScore ||
       other.tmDiff < c.tmDiff ||
       other.balanceScore > c.balanceScore)
    );
  });

  if (verbose) {
    console.log(`  Evaluated ${candidates.length} length combinations`);
    console.log(`  Found ${paretoOptimal.length} Pareto-optimal solutions`);
    paretoOptimal.slice(0, 3).forEach((c, i) => {
      console.log(`    ${i + 1}. Fwd: ${c.fwdLen}bp, Rev: ${c.revLen}bp | ` +
        `TmΔ: ${c.tmDiff.toFixed(1)}° | 3' Quality: ${c.fwd3PrimeQuality}/${c.rev3PrimeQuality}`);
    });
  }

  return {
    candidates: paretoOptimal.slice(0, 5),
    bestCandidate: paretoOptimal[0] || null,
    currentPairScore,
    currentTmDiff,
  };
}

/**
 * Return a matrix of primers for (i, j) where (i, j) are the start/end indexes
 */
function createPrimers(
  factory: PrimerFactory,
  seq: string,
  offtargetCheck: string,
  startRange: number[],
  endRange: number[],
  fwd: boolean,
  addLen: number
): (Primer | null)[][] {
  const gc = gcCache(seq);
  const tm = tmCache(seq);
  const dg = dgCache(seq);
  const ot = offTargets(seq, offtargetCheck);

  const ps: (Primer | null)[][] = [];
  for (let i = 0; i < gc.length; i++) {
    ps.push(new Array(gc.length).fill(null));
  }

  for (const s of startRange) {
    for (const e of endRange) {
      const pSeq = seq.slice(s, e + 1);
      const pTm = tm[s + addLen][e];
      const pTmTotal = tm[s][e];
      const pGc = gc[s][e];
      const pDg = dg[s][e];
      const pOt = ot[e];

      ps[s][e] = factory.build(pSeq, pTm, pTmTotal, pGc, pDg, fwd, pOt);
    }
  }

  return ps;
}

/**
 * Assign primer pair candidates to quality tiers based on Tm difference and structure.
 */
function assignPrimerTiers(
  candidates: PairCandidate[],
  options: { minScoreForTier1?: number } = {}
): PrimerTiers {
  const { minScoreForTier1 = 70 } = options;

  const tiers: PrimerTiers = {
    tier1: [],
    tier2: [],
    tier3: [],
    tier4: [],
  };

  for (const c of candidates) {
    const tmDiff = Math.abs(c.fwd.tm - c.rev.tm);
    const fwdDg = c.fwd.dg || 0;
    const revDg = c.rev.dg || 0;
    const worstDg = Math.min(fwdDg, revDg);

    const hasCompositeScore = c.score && c.score > 0;
    const scoreTooLow = hasCompositeScore && c.score < minScoreForTier1;

    if (tmDiff <= 2 && worstDg > -3 && !scoreTooLow) {
      tiers.tier1.push(c);
    } else if (tmDiff <= 5 && worstDg > -5) {
      tiers.tier2.push(c);
    } else if (tmDiff <= 8) {
      tiers.tier3.push(c);
    } else {
      tiers.tier4.push(c);
    }
  }

  return tiers;
}

/**
 * Choose the best combo of primers.
 */
function chooseBest(
  factory: PrimerFactory,
  fwdPrimers: (Primer | null)[][],
  revPrimers: (Primer | null)[][],
  options: ChooseBestOptions = {}
): [Primer, Primer] | ChooseBestResult {
  const {
    useCompositeScore = false,
    annealingTemp = 55,
    template = null,
    exhaustiveSearch = false,
    returnCandidates = false,
  } = options;

  const CANDIDATES_PER_GROUP = exhaustiveSearch ? 100 : 10;
  const TOP_N_INITIAL = exhaustiveSearch ? 50 : 10;

  const allFwd: Primer[] = [];
  const allRev: Primer[] = [];

  for (const row of fwdPrimers) {
    for (const p of row) {
      if (!p) continue;
      allFwd.push(p);
    }
  }

  for (const row of revPrimers) {
    for (const p of row) {
      if (!p) continue;
      allRev.push(p);
    }
  }

  if (allFwd.length === 0) {
    throw new Error("Failed to create any primers in the FWD direction");
  }

  if (allRev.length === 0) {
    throw new Error("Failed to create any primers in the REV direction");
  }

  allFwd.sort((a, b) => a.scoring.penalty - b.scoring.penalty);
  allRev.sort((a, b) => a.scoring.penalty - b.scoring.penalty);

  const topNFwd = allFwd.slice(0, TOP_N_INITIAL);
  const topNRev = allRev.slice(0, TOP_N_INITIAL);

  const buildCandidates = (fwdList: Primer[], revList: Primer[], withScoring: boolean): PairCandidate[] => {
    const candidates: PairCandidate[] = [];
    for (const fwd of fwdList) {
      for (const rev of revList) {
        if (withScoring) {
          const [scoredFwd, scoredRev] = factory.addPiecewiseScoring(
            fwd, rev, annealingTemp, template
          );
          candidates.push({
            fwd: scoredFwd,
            rev: scoredRev!,
            score: scoredFwd.scoring.compositeScore!,
            penalty: scoredFwd.scoring.penalty + scoredRev!.scoring.penalty,
          });
        } else {
          const [newFwd, newRev] = factory.buildPair(fwd, rev);
          candidates.push({
            fwd: newFwd,
            rev: newRev,
            score: 0,
            penalty: newFwd.scoring.penalty + newRev.scoring.penalty,
          });
        }
      }
    }
    return candidates;
  };

  const selectFromTiers = (tiers: PrimerTiers, byScore: boolean): { best: PairCandidate | null; tier: number } => {
    if (byScore) {
      for (const tier of Object.values(tiers)) {
        tier.sort((a: any, b: any) => b.score - a.score);
      }
    } else {
      for (const tier of Object.values(tiers)) {
        tier.sort((a: any, b: any) => a.penalty - b.penalty);
      }
    }

    if (tiers.tier1.length > 0) return { best: tiers.tier1[0], tier: 1 };
    if (tiers.tier2.length > 0) return { best: tiers.tier2[0], tier: 2 };
    if (tiers.tier3.length > 0) return { best: tiers.tier3[0], tier: 3 };
    if (tiers.tier4.length > 0) return { best: tiers.tier4[0], tier: 4 };
    return { best: null, tier: 0 };
  };

  const initialCandidates = buildCandidates(topNFwd, topNRev, useCompositeScore);
  const initialTiers = assignPrimerTiers(initialCandidates);
  const initialResult = selectFromTiers(initialTiers, useCompositeScore);

  const debugTmDiff = initialResult.best
    ? Math.abs(initialResult.best.fwd.tm - initialResult.best.rev.tm).toFixed(1)
    : 'N/A';
  const debugScore = initialResult.best?.score || 'N/A';
  const modeLabel = exhaustiveSearch ? '[EXHAUSTIVE]' : '[standard]';
  console.log(`[chooseBest] ${modeLabel} Initial: tier=${initialResult.tier}, TmDiff=${debugTmDiff}°C, score=${debugScore}, candidates=${initialCandidates.length}, tiers: t1=${initialTiers.tier1.length}, t2=${initialTiers.tier2.length}, t3=${initialTiers.tier3.length}, t4=${initialTiers.tier4.length}`);

  const initialTmDiff = initialResult.best
    ? Math.abs(initialResult.best.fwd.tm - initialResult.best.rev.tm)
    : Infinity;

  if (initialResult.best && initialResult.tier === 1 && initialTmDiff <= 2 && !exhaustiveSearch) {
    console.log(`[chooseBest] Returning excellent tier 1 result (score=${debugScore})`);
    if (returnCandidates) {
      return {
        best: [initialResult.best.fwd, initialResult.best.rev],
        allCandidates: initialCandidates,
      };
    }
    return [initialResult.best.fwd, initialResult.best.rev];
  }

  if (exhaustiveSearch && initialResult.best && initialResult.tier === 1) {
    console.log(`[chooseBest] Exhaustive: continuing search for candidates with fewer issues...`);
  }

  const groupByTm = (primers: Primer[]): Map<number, Primer[]> => {
    const groups = new Map<number, Primer[]>();
    for (const p of primers) {
      const tmBucket = Math.round(p.tm / 2) * 2;
      if (!groups.has(tmBucket)) {
        groups.set(tmBucket, []);
      }
      groups.get(tmBucket)!.push(p);
    }
    return groups;
  };

  const fwdByTm = groupByTm(allFwd);
  const revByTm = groupByTm(allRev);

  const fwdTms = [...fwdByTm.keys()].sort((a, b) => a - b);
  const revTms = [...revByTm.keys()].sort((a, b) => a - b);

  const tmMatchedCandidates: PairCandidate[] = [];

  const calculateTmMatchPenalty = (fwd: Primer, rev: Primer): { penalty: number; tmDiff: number } => {
    let penalty = 0;
    const tmDiff = Math.abs(fwd.tm - rev.tm);

    penalty += 2.0 * Math.pow(tmDiff, 2);

    const optimalLen = 22;
    const lenPenalty = 0.05;
    penalty += Math.abs(fwd.seq.length - optimalLen) * lenPenalty;
    penalty += Math.abs(rev.seq.length - optimalLen) * lenPenalty;

    if (fwd.dg < -2) penalty += Math.abs(fwd.dg) * 0.5;
    if (rev.dg < -2) penalty += Math.abs(rev.dg) * 0.5;

    return { penalty, tmDiff };
  };

  let pairsEvaluated = 0;

  const selectTopByPenalty = (arr: Primer[], limit: number): Primer[] => {
    if (arr.length <= limit) return arr;
    const sorted = [...arr].sort((a, b) => a.scoring.penalty - b.scoring.penalty);
    return sorted.slice(0, limit);
  };

  for (const fwdTm of fwdTms) {
    for (const revTm of revTms) {
      if (Math.abs(fwdTm - revTm) <= 8) {
        const fwdGroup = fwdByTm.get(fwdTm)!;
        const revGroup = revByTm.get(revTm)!;

        const fwdSelected = selectTopByPenalty(fwdGroup, CANDIDATES_PER_GROUP);
        const revSelected = selectTopByPenalty(revGroup, CANDIDATES_PER_GROUP);

        for (const fwd of fwdSelected) {
          for (const rev of revSelected) {
            const { penalty: tmMatchPenalty, tmDiff } = calculateTmMatchPenalty(fwd, rev);

            if (tmDiff > 5) continue;

            pairsEvaluated++;

            if (useCompositeScore) {
              const [scoredFwd, scoredRev] = factory.addPiecewiseScoring(
                fwd, rev, annealingTemp, template
              );
              const candidate: PairCandidate = {
                fwd: scoredFwd,
                rev: scoredRev!,
                score: scoredFwd.scoring.compositeScore!,
                penalty: tmMatchPenalty,
                tmDiff,
              };
              tmMatchedCandidates.push(candidate);
            } else {
              const [newFwd, newRev] = factory.buildPair(fwd, rev);
              tmMatchedCandidates.push({
                fwd: newFwd,
                rev: newRev,
                score: 0,
                penalty: tmMatchPenalty,
                tmDiff,
              });
            }
          }
        }
      }
    }
  }

  if (exhaustiveSearch) {
    console.log(`[chooseBest] Exhaustive search evaluated ${pairsEvaluated} pairs`);
  }

  if (tmMatchedCandidates.length > 0) {
    const tmMatchedTiers = assignPrimerTiers(tmMatchedCandidates);
    const tmMatchedResult = selectFromTiers(tmMatchedTiers, false);

    if (tmMatchedResult.best) {
      const tmMatchedTmDiff = Math.abs(tmMatchedResult.best.fwd.tm - tmMatchedResult.best.rev.tm);

      const tmMatchedScore = tmMatchedResult.best.score || 0;
      console.log(`[chooseBest] Joint Tm opt: tier=${tmMatchedResult.tier}, TmDiff=${tmMatchedTmDiff.toFixed(1)}°C, score=${tmMatchedScore}`);

      let bestScoreCandidate = tmMatchedResult.best;

      if (exhaustiveSearch && useCompositeScore) {
        const countIssuesComprehensive = (c: PairCandidate): number => {
          let issues = 0;
          const fwdSeq = c.fwd.seq;
          const revSeq = c.rev.seq;
          const fwdScores = c.fwd.scoring?.piecewiseScores || {};
          const revScores = c.rev.scoring?.piecewiseScores || {};

          const fwdLast = fwdSeq.slice(-1).toUpperCase();
          const revLast = revSeq.slice(-1).toUpperCase();
          if (fwdLast !== 'G' && fwdLast !== 'C') issues += 1;
          if (revLast !== 'G' && revLast !== 'C') issues += 1;

          const fwdDg = c.fwd.dg || 0;
          const revDg = c.rev.dg || 0;
          if (fwdDg < -3) issues += 1;
          if (revDg < -3) issues += 1;

          if ((c.fwd.offTargetCount || 0) > 0) issues += 2;
          if ((c.rev.offTargetCount || 0) > 0) issues += 2;

          if ((fwdScores.terminal3DG || 1) < 0.7) issues += 1;
          if ((revScores.terminal3DG || 1) < 0.7) issues += 1;

          if ((fwdScores.heterodimer || 1) < 0.7) issues += 1;

          if ((fwdScores.selfDimerFwd || 1) < 0.5) issues += 1;
          if ((revScores.selfDimerRev || 1) < 0.5) issues += 1;

          if ((revScores.gQuadruplexRev || 1) < 0.5) issues += 2;
          if ((fwdScores.gQuadruplexFwd || 1) < 0.5) issues += 1;

          return issues;
        };

        const tier1Candidates = tmMatchedCandidates.filter(c => {
          const cTmDiff = Math.abs(c.fwd.tm - c.rev.tm);
          return cTmDiff <= 2 && c.score > 0;
        });

        if (tier1Candidates.length > 0) {
          tier1Candidates.sort((a, b) => {
            if (b.score !== a.score) return b.score - a.score;
            return countIssuesComprehensive(a) - countIssuesComprehensive(b);
          });

          const bestCandidate = tier1Candidates[0];

          const bestIssues = countIssuesComprehensive(bestCandidate);
          const currentIssues = countIssuesComprehensive(bestScoreCandidate);
          const bestTmDiff = Math.abs(bestCandidate.fwd.tm - bestCandidate.rev.tm);

          console.log(`[chooseBest] Exhaustive: evaluated ${tier1Candidates.length} tier-1 candidates`);
          console.log(`[chooseBest] Exhaustive: best score=${bestCandidate.score}, issues=${bestIssues}, TmDiff=${bestTmDiff.toFixed(1)}°C`);

          bestScoreCandidate = bestCandidate;
        } else {
          console.log(`[chooseBest] Exhaustive: no tier-1 candidates, using best from lower tiers`);
        }
      }

      const jointScore = bestScoreCandidate.score || 0;
      const initialScore = initialResult.best?.score || 0;
      const jointTmDiff = Math.abs(bestScoreCandidate.fwd.tm - bestScoreCandidate.rev.tm);

      const shouldUseJointResult = (
        tmMatchedResult.tier < initialResult.tier ||
        (tmMatchedResult.tier === initialResult.tier && jointTmDiff < initialTmDiff - 1) ||
        (tmMatchedResult.tier === initialResult.tier && jointTmDiff <= 2 && initialTmDiff <= 2 && jointScore > initialScore) ||
        exhaustiveSearch
      );

      if (shouldUseJointResult) {
        console.log(`[chooseBest] Using joint result: score=${jointScore}, TmDiff=${jointTmDiff.toFixed(1)}°C`);
        if (returnCandidates) {
          return {
            best: [bestScoreCandidate.fwd, bestScoreCandidate.rev],
            allCandidates: tmMatchedCandidates,
          };
        }
        return [bestScoreCandidate.fwd, bestScoreCandidate.rev];
      }
    }
  }

  if (initialResult.best) {
    console.log(`[chooseBest] Using initial result (TmDiff=${initialTmDiff.toFixed(1)}°C, score=${debugScore})`);
    if (returnCandidates) {
      return {
        best: [initialResult.best.fwd, initialResult.best.rev],
        allCandidates: initialCandidates,
      };
    }
    return [initialResult.best.fwd, initialResult.best.rev];
  }

  throw new Error("Failed to create a pair of PCR primers");
}

/**
 * Generate alternative primer pairs from existing primer matrices.
 */
function generateAlternativesInternal(
  fwdPrimers: (Primer | null)[][],
  revPrimers: (Primer | null)[][],
  factory: PrimerFactory,
  template: string,
  annealingTemp: number,
  options: { numAlternatives?: number; currentFwd?: Primer | null; currentRev?: Primer | null } = {}
): AlternativePair[] {
  const { numAlternatives = 5, currentFwd = null, currentRev = null } = options;

  const allFwd: Primer[] = [];
  const allRev: Primer[] = [];

  for (const row of fwdPrimers) {
    for (const p of row) {
      if (p) allFwd.push(p);
    }
  }

  for (const row of revPrimers) {
    for (const p of row) {
      if (p) allRev.push(p);
    }
  }

  if (allFwd.length === 0 || allRev.length === 0) {
    return [];
  }

  const groupByTm = (primers: Primer[]): Map<number, Primer[]> => {
    const groups = new Map<number, Primer[]>();
    for (const p of primers) {
      const tmBucket = Math.round(p.tm / 2) * 2;
      if (!groups.has(tmBucket)) groups.set(tmBucket, []);
      groups.get(tmBucket)!.push(p);
    }
    return groups;
  };

  const fwdByTm = groupByTm(allFwd);
  const revByTm = groupByTm(allRev);
  const fwdTms = [...fwdByTm.keys()].sort((a, b) => a - b);
  const revTms = [...revByTm.keys()].sort((a, b) => a - b);

  const calculateTmMatchPenalty = (fwd: Primer, rev: Primer): { penalty: number; tmDiff: number } => {
    let penalty = 0;
    const tmDiff = Math.abs(fwd.tm - rev.tm);
    penalty += 2.0 * Math.pow(tmDiff, 2);
    const lenPenalty = 0.05;
    penalty += Math.abs(fwd.seq.length - 22) * lenPenalty;
    penalty += Math.abs(rev.seq.length - 22) * lenPenalty;
    if (fwd.dg < -2) penalty += Math.abs(fwd.dg) * 0.5;
    if (rev.dg < -2) penalty += Math.abs(rev.dg) * 0.5;
    return { penalty, tmDiff };
  };

  const selectDiverse = (arr: Primer[], limit: number): Primer[] => {
    if (arr.length <= limit) return arr;
    const sorted = [...arr].sort((a, b) => a.seq.length - b.seq.length);
    const result: Primer[] = [];
    const step = (sorted.length - 1) / (limit - 1);
    for (let i = 0; i < limit; i++) {
      result.push(sorted[Math.round(i * step)]);
    }
    return result;
  };

  const allPairs: any[] = [];

  for (const fwdTm of fwdTms) {
    for (const revTm of revTms) {
      if (Math.abs(fwdTm - revTm) <= 10) {
        const fwdGroup = fwdByTm.get(fwdTm)!;
        const revGroup = revByTm.get(revTm)!;
        const fwdSelected = selectDiverse(fwdGroup, 6);
        const revSelected = selectDiverse(revGroup, 6);

        for (const fwd of fwdSelected) {
          for (const rev of revSelected) {
            const { penalty: tmMatchPenalty, tmDiff } = calculateTmMatchPenalty(fwd, rev);

            if (tmDiff > 8) continue;

            const [scoredFwd, scoredRev] = factory.addPiecewiseScoring(fwd, rev, annealingTemp, template);

            allPairs.push({
              fwd: scoredFwd,
              rev: scoredRev,
              tmDiff,
              penalty: tmMatchPenalty,
              compositeScore: scoredFwd.scoring.compositeScore!,
            });
          }
        }
      }
    }
  }

  allPairs.sort((a, b) => {
    const tierA = a.tmDiff <= 2 ? 1 : a.tmDiff <= 5 ? 2 : a.tmDiff <= 8 ? 3 : 4;
    const tierB = b.tmDiff <= 2 ? 1 : b.tmDiff <= 5 ? 2 : b.tmDiff <= 8 ? 3 : 4;
    if (tierA !== tierB) return tierA - tierB;
    return b.compositeScore - a.compositeScore;
  });

  const currentKey = currentFwd && currentRev
    ? `${currentFwd.seq.length}-${currentRev.seq.length}`
    : null;

  const seen = new Set<string>();
  if (currentKey) seen.add(currentKey);

  const alternatives: AlternativePair[] = [];

  const MIN_SCORE = 50;
  const MAX_PRIMER_TM = 72;

  for (const pair of allPairs) {
    const key = `${pair.fwd.seq.length}-${pair.rev.seq.length}`;
    if (seen.has(key)) continue;

    if (pair.compositeScore < MIN_SCORE) continue;
    if (pair.fwd.tm > MAX_PRIMER_TM || pair.rev.tm > MAX_PRIMER_TM) continue;

    seen.add(key);

    let tierLabel = 'poor';
    if (pair.compositeScore >= 90) tierLabel = 'excellent';
    else if (pair.compositeScore >= 75) tierLabel = 'good';
    else if (pair.compositeScore >= 60) tierLabel = 'acceptable';

    alternatives.push({
      forward: {
        sequence: pair.fwd.seq,
        length: pair.fwd.seq.length,
        tm: pair.fwd.tm,
        gc: pair.fwd.gc,
        dg: pair.fwd.dg,
        gcPercent: `${(pair.fwd.gc * 100).toFixed(1)}%`,
        hasGCClamp: /[GC]$/.test(pair.fwd.seq),
      },
      reverse: {
        sequence: pair.rev.seq,
        length: pair.rev.seq.length,
        tm: pair.rev.tm,
        gc: pair.rev.gc,
        dg: pair.rev.dg,
        gcPercent: `${(pair.rev.gc * 100).toFixed(1)}%`,
        hasGCClamp: /[GC]$/.test(pair.rev.seq),
      },
      tmDiff: Math.round(pair.tmDiff * 10) / 10,
      compositeScore: pair.compositeScore,
      qualityTier: tierLabel,
    });

    if (alternatives.length >= numAlternatives) break;
  }

  return alternatives;
}

/**
 * Generate diverse alternatives from pre-evaluated candidates.
 */
function generateDiverseAlternatives(
  candidates: PairCandidate[],
  currentBest: { fwd: Primer; rev: Primer },
  options: { maxPerCategory?: number } = {}
): { alternatives: AlternativePair[]; categories: Record<string, AlternativePair[]> } {
  const { maxPerCategory = 3 } = options;

  if (!candidates || candidates.length === 0) {
    return { alternatives: [], categories: {} };
  }

  const MIN_SCORE = 50;
  const MAX_PRIMER_TM = 72;
  const MAX_TM_DIFF = 5;

  const getQualityTier = (score: number): string => {
    if (score >= 90) return 'excellent';
    if (score >= 75) return 'good';
    if (score >= 60) return 'acceptable';
    return 'poor';
  };

  const formatCandidate = (c: PairCandidate, category: string): AlternativePair => {
    const tmDiff = Math.abs(c.fwd.tm - c.rev.tm);
    const score = c.score || c.compositeScore || 0;

    return {
      forward: {
        sequence: c.fwd.seq,
        length: c.fwd.seq.length,
        tm: c.fwd.tm,
        gc: c.fwd.gc,
        dg: c.fwd.dg,
        gcPercent: `${(c.fwd.gc * 100).toFixed(1)}%`,
        hasGCClamp: /[GC]$/.test(c.fwd.seq),
      },
      reverse: {
        sequence: c.rev.seq,
        length: c.rev.seq.length,
        tm: c.rev.tm,
        gc: c.rev.gc,
        dg: c.rev.dg,
        gcPercent: `${(c.rev.gc * 100).toFixed(1)}%`,
        hasGCClamp: /[GC]$/.test(c.rev.seq),
      },
      tmDiff: Math.round(tmDiff * 10) / 10,
      compositeScore: score,
      qualityTier: getQualityTier(score),
      category,
    };
  };

  const currentKey = currentBest
    ? `${currentBest.fwd.seq}-${currentBest.rev.seq}`
    : null;

  const filtered = candidates.filter(c => {
    const key = `${c.fwd.seq}-${c.rev.seq}`;
    if (key === currentKey) return false;

    const score = c.score || c.compositeScore || 0;
    const hasData = score > 0 || c.penalty !== undefined;
    if (!hasData) return false;

    if (score < MIN_SCORE) return false;
    if (c.fwd.tm > MAX_PRIMER_TM || c.rev.tm > MAX_PRIMER_TM) return false;
    if (Math.abs(c.fwd.tm - c.rev.tm) > MAX_TM_DIFF) return false;

    const heterodimerDG = c.heterodimerDG || 0;
    if (heterodimerDG < -10) return false;

    return true;
  });

  const countIssues = (c: PairCandidate): number => {
    let issues = 0;
    const fwdSeq = c.fwd.seq;
    const revSeq = c.rev.seq;

    if (!/[GC]$/.test(fwdSeq)) issues += 1;
    if (!/[GC]$/.test(revSeq)) issues += 1;

    if ((c.fwd.dg || 0) < -3) issues += 1;
    if ((c.rev.dg || 0) < -3) issues += 1;

    if ((c.fwd.offTargetCount || 0) > 0) issues += 1;
    if ((c.rev.offTargetCount || 0) > 0) issues += 1;

    if (c.fwd.tm > 68) issues += 1;
    if (c.rev.tm > 68) issues += 1;

    if (c.fwd.seq.length > 30) issues += 1;
    if (c.rev.seq.length > 30) issues += 1;

    const heterodimerDG = c.heterodimerDG || 0;
    if (heterodimerDG < -12) issues += 2;
    else if (heterodimerDG < -9) issues += 1;

    return issues;
  };

  const categories: Record<string, AlternativePair[]> = {};
  const seenKeys = new Set<string>();
  if (currentKey) seenKeys.add(currentKey);

  const addToCategory = (name: string, sortedCandidates: PairCandidate[]): void => {
    categories[name] = [];
    for (const c of sortedCandidates) {
      const key = `${c.fwd.seq}-${c.rev.seq}`;
      if (seenKeys.has(key)) continue;
      seenKeys.add(key);
      categories[name].push(formatCandidate(c, name));
      if (categories[name].length >= maxPerCategory) break;
    }
  };

  const byScore = [...filtered].sort((a, b) =>
    (b.score || b.compositeScore || 0) - (a.score || a.compositeScore || 0)
  );
  addToCategory('highestScore', byScore);

  const byShortestLength = [...filtered].sort((a, b) => {
    const lenA = a.fwd.seq.length + a.rev.seq.length;
    const lenB = b.fwd.seq.length + b.rev.seq.length;
    if (lenA !== lenB) return lenA - lenB;
    return (b.score || 0) - (a.score || 0);
  });
  addToCategory('shortestPrimers', byShortestLength);

  const byLongestLength = [...filtered].sort((a, b) => {
    const lenA = a.fwd.seq.length + a.rev.seq.length;
    const lenB = b.fwd.seq.length + b.rev.seq.length;
    if (lenA !== lenB) return lenB - lenA;
    return (b.score || 0) - (a.score || 0);
  });
  addToCategory('longestPrimers', byLongestLength);

  const byTmMatch = [...filtered].sort((a, b) => {
    const tmDiffA = Math.abs(a.fwd.tm - a.rev.tm);
    const tmDiffB = Math.abs(b.fwd.tm - b.rev.tm);
    if (Math.abs(tmDiffA - tmDiffB) > 0.1) return tmDiffA - tmDiffB;
    return (b.score || 0) - (a.score || 0);
  });
  addToCategory('bestTmMatch', byTmMatch);

  const byFewestIssues = [...filtered].sort((a, b) => {
    const issuesA = countIssues(a);
    const issuesB = countIssues(b);
    if (issuesA !== issuesB) return issuesA - issuesB;
    return (b.score || 0) - (a.score || 0);
  });
  addToCategory('fewestIssues', byFewestIssues);

  const allAlternatives: AlternativePair[] = [];
  const finalSeen = new Set<string>();
  if (currentKey) finalSeen.add(currentKey);

  for (let i = 0; i < maxPerCategory; i++) {
    for (const catName of ['highestScore', 'fewestIssues', 'bestTmMatch', 'shortestPrimers', 'longestPrimers']) {
      if (categories[catName] && categories[catName][i]) {
        const alt = categories[catName][i];
        const key = `${alt.forward.sequence}-${alt.reverse.sequence}`;
        if (!finalSeen.has(key)) {
          finalSeen.add(key);
          allAlternatives.push(alt);
        }
      }
    }
  }

  console.log(`[generateDiverseAlternatives] Generated ${allAlternatives.length} diverse alternatives from ${filtered.length} candidates`);

  return {
    alternatives: allAlternatives.slice(0, 10),
    categories,
  };
}

/**
 * Validate and parse the input sequence.
 */
function parse(seq: string, offtargetCheck: string): [string, string] {
  seq = seq.toUpperCase();
  offtargetCheck = (offtargetCheck || seq).toUpperCase();

  const diff = new Set([...seq].filter((c) => !"ATGC".includes(c)));

  if (diff.size > 0) {
    const desc = [...diff].join(",");
    throw new Error(`Invalid non-DNA bases found: ${desc}`);
  }

  return [seq, offtargetCheck];
}

/**
 * Parse the additional length range
 */
function parseAddLen(add: string, addLen: [number, number]): [number, number] {
  if (!add) {
    return [0, 0];
  }

  const [addMin, addMax] = addLen;

  if (addMin === -1 && addMax === -1) {
    const addAssume = Math.min(40, add.length);
    if (addAssume !== add.length) {
      console.warn(
        `${add.length}bp additional sequence added, but no \`addDirLen\` argument provided:\n\tAdding between ${addAssume - 10} and ${addAssume}bp`
      );
      return [addAssume - 10, addAssume];
    }
    return [addAssume, addAssume];
  }

  if (addMin < 0 && addMax >= 0) {
    return [0, Math.min(add.length, addMax)];
  }
  if (addMin >= 0 && addMax < 0) {
    return [Math.min(add.length, addMin), add.length];
  }

  if (addMin > -1 && addMax > -1 && addMin > addMax) {
    throw new Error(`addDirLen range has a min > max: ${addMin} > ${addMax}`);
  }

  return [Math.max(addMin, 0), Math.min(addMax, add.length)];
}

/**
 * Return the reverse complement of a DNA sequence.
 */
function reverseComplement(seq: string): string {
  const rc: Record<string, string> = { A: "T", T: "A", G: "C", C: "G" };
  return seq
    .split("")
    .reverse()
    .map((c) => rc[c])
    .join("");
}

/**
 * Create a range array
 */
function range(start: number, end: number): number[] {
  const result: number[] = [];
  for (let i = start; i < end; i++) {
    result.push(i);
  }
  return result;
}

/**
 * Calculate Tm of a sequence using the cached Tm calculator
 */
function calcTm(seq: string): number {
  if (!seq || seq.length < 2) return 0;
  const tm = tmCache(seq);
  return tm[0][seq.length - 1];
}

/**
 * Generate alternative primer pairs with different length combinations.
 */
export function generateAlternatives(
  seq: string,
  {
    addFwd = "",
    addRev = "",
    addFwdLen = [-1, -1] as [number, number],
    addRevLen = [-1, -1] as [number, number],
    offtargetCheck = "",
    optimalTm = 62.0,
    optimalGc = 0.5,
    optimalLen = 22,
    penaltyTm = 1.0,
    penaltyGc = 0.2,
    penaltyLen = 0.5,
    penaltyTmDiff = 1.0,
    penaltyDg = 2.0,
    penaltyOffTarget = 20.0,
    numAlternatives = 3,
  }: PrimersOptions & { numAlternatives?: number } = {}
): AlternativePair[] {
  const [parsedSeq, parsedOfftargetCheck] = parse(seq, offtargetCheck);
  const [parsedAddFwd] = parse(addFwd, "");
  const [parsedAddRev] = parse(addRev, "");

  const factory = createPrimerFactory({
    optimalTm,
    optimalGc,
    optimalLen,
    penaltyTm,
    penaltyGc,
    penaltyLen,
    penaltyTmDiff,
    penaltyDg,
    penaltyOffTarget,
  });

  const [addFwdMin, addFwdMax] = parseAddLen(parsedAddFwd, addFwdLen);
  const [addRevMin, addRevMax] = parseAddLen(parsedAddRev, addRevLen);

  const trimmedAddFwd = parsedAddFwd.slice(-addFwdMax);
  const trimmedAddRev = reverseComplement(parsedAddRev).slice(0, addRevMax);
  const seqFull = trimmedAddFwd + parsedSeq + trimmedAddRev;

  if (seqFull.length < LEN_MAX) {
    throw new Error(
      `Template sequence length is too short: ${seqFull.length}bp < ${LEN_MAX}bp`
    );
  }

  const optimalFwdLen = Math.round(optimalLen + (addFwdMin + addFwdMax) / 2);
  const fwdSeq = seqFull.slice(0, addFwdMax + LEN_MAX);
  const fwdPrimers = createPrimers(
    factory.withOptimalLen(optimalFwdLen),
    fwdSeq,
    parsedOfftargetCheck,
    range(0, addFwdMax - addFwdMin + 1),
    range(LEN_MIN - 1 + addFwdMax, LEN_MAX + addFwdMax),
    true,
    addFwdMax
  );

  const optimalRevLen = Math.round(optimalLen + (addRevMin + addRevMax) / 2);
  let revSeq = reverseComplement(seqFull);
  revSeq = revSeq.slice(0, addRevMax + LEN_MAX);
  const revPrimers = createPrimers(
    factory.withOptimalLen(optimalRevLen),
    revSeq,
    parsedOfftargetCheck,
    range(0, addRevMax - addRevMin + 1),
    range(LEN_MIN - 1 + addRevMax, LEN_MAX + addRevMin),
    false,
    addRevMax
  );

  const template = parsedOfftargetCheck || parsedSeq;
  const allPairs: any[] = [];

  const flatFwd: Primer[] = [];
  const flatRev: Primer[] = [];

  for (const row of fwdPrimers) {
    for (const p of row) {
      if (p) flatFwd.push(p);
    }
  }

  for (const row of revPrimers) {
    for (const p of row) {
      if (p) flatRev.push(p);
    }
  }

  const groupByTm = (primers: Primer[]): Map<number, Primer[]> => {
    const groups = new Map<number, Primer[]>();
    for (const p of primers) {
      const tmBucket = Math.round(p.tm / 2) * 2;
      if (!groups.has(tmBucket)) groups.set(tmBucket, []);
      groups.get(tmBucket)!.push(p);
    }
    return groups;
  };

  const fwdByTm = groupByTm(flatFwd);
  const revByTm = groupByTm(flatRev);
  const fwdTms = [...fwdByTm.keys()].sort((a, b) => a - b);
  const revTms = [...revByTm.keys()].sort((a, b) => a - b);

  const calculateTmMatchPenalty = (fwd: Primer, rev: Primer): { penalty: number; tmDiff: number } => {
    let penalty = 0;
    const tmDiff = Math.abs(fwd.tm - rev.tm);
    penalty += 2.0 * Math.pow(tmDiff, 2);
    const lenPenalty = 0.05;
    penalty += Math.abs(fwd.seq.length - 22) * lenPenalty;
    penalty += Math.abs(rev.seq.length - 22) * lenPenalty;
    if (fwd.dg < -2) penalty += Math.abs(fwd.dg) * 0.5;
    if (rev.dg < -2) penalty += Math.abs(rev.dg) * 0.5;
    return { penalty, tmDiff };
  };

  const selectDiverse = (arr: Primer[], limit: number): Primer[] => {
    if (arr.length <= limit) return arr;
    const sorted = [...arr].sort((a, b) => a.seq.length - b.seq.length);
    const result: Primer[] = [];
    const step = (sorted.length - 1) / (limit - 1);
    for (let i = 0; i < limit; i++) {
      result.push(sorted[Math.round(i * step)]);
    }
    return result;
  };

  for (const fwdTm of fwdTms) {
    for (const revTm of revTms) {
      if (Math.abs(fwdTm - revTm) <= 10) {
        const fwdGroup = fwdByTm.get(fwdTm)!;
        const revGroup = revByTm.get(revTm)!;
        const fwdSelected = selectDiverse(fwdGroup, 8);
        const revSelected = selectDiverse(revGroup, 8);

        for (const fwd of fwdSelected) {
          for (const rev of revSelected) {
            const { penalty: tmMatchPenalty, tmDiff } = calculateTmMatchPenalty(fwd, rev);
            const [scoredFwd, scoredRev] = factory.addPiecewiseScoring(fwd, rev, optimalTm, template);
            allPairs.push({
              forward: scoredFwd,
              reverse: scoredRev,
              compositeScore: scoredFwd.scoring.compositeScore!,
              qualityTier: scoredFwd.scoring.qualityTier,
              penalty: tmMatchPenalty,
              tmDiff: tmDiff,
            });
          }
        }
      }
    }
  }

  const tieredCandidates: PairCandidate[] = allPairs.map(pair => ({
    fwd: pair.forward,
    rev: pair.reverse!,
    score: pair.compositeScore,
    penalty: pair.penalty,
    qualityTier: pair.qualityTier,
    tmDiff: pair.tmDiff,
  }));

  const tiers = assignPrimerTiers(tieredCandidates);

  for (const tier of Object.values(tiers)) {
    tier.sort((a: any, b: any) => b.score - a.score);
  }

  const sortedPairs = [
    ...tiers.tier1,
    ...tiers.tier2,
    ...tiers.tier3,
    ...tiers.tier4,
  ];

  const seen = new Set<string>();
  const alternatives: AlternativePair[] = [];

  for (const pair of sortedPairs) {
    const key = `${pair.fwd.seq.length}-${pair.rev.seq.length}`;
    if (!seen.has(key)) {
      seen.add(key);
      const tmDiff = Math.abs(pair.fwd.tm - pair.rev.tm);
      let tierLabel = 'poor';
      if (tmDiff <= 2) tierLabel = 'excellent';
      else if (tmDiff <= 5) tierLabel = 'good';
      else if (tmDiff <= 8) tierLabel = 'acceptable';

      alternatives.push({
        forward: {
          sequence: pair.fwd.seq,
          length: pair.fwd.seq.length,
          tm: pair.fwd.tm,
          gc: pair.fwd.gc,
          dg: pair.fwd.dg,
          gcPercent: `${(pair.fwd.gc * 100).toFixed(1)}%`,
          hasGCClamp: pair.fwd.seq[pair.fwd.seq.length - 1] === 'G' ||
                      pair.fwd.seq[pair.fwd.seq.length - 1] === 'C',
        },
        reverse: {
          sequence: pair.rev.seq,
          length: pair.rev.seq.length,
          tm: pair.rev.tm,
          gc: pair.rev.gc,
          dg: pair.rev.dg,
          gcPercent: `${(pair.rev.gc * 100).toFixed(1)}%`,
          hasGCClamp: pair.rev.seq[pair.rev.seq.length - 1] === 'G' ||
                      pair.rev.seq[pair.rev.seq.length - 1] === 'C',
        },
        tmDiff: Math.round(tmDiff * 10) / 10,
        compositeScore: pair.score,
        qualityTier: tierLabel,
        penalty: Math.round(pair.penalty * 10) / 10,
      });

      if (alternatives.length >= numAlternatives) break;
    }
  }

  return alternatives;
}

// Alias for primers function
export const create = primers;
