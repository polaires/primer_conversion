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

/**
 * Create a Scoring object
 * @param {Object} params - Scoring parameters
 * @returns {Object} Scoring object
 */
function createScoring({
  penalty = 0,
  penaltyTm = 0,
  penaltyTmDiff = 0,
  penaltyGc = 0,
  penaltyLen = 0,
  penaltyDg = 0,
  penaltyOffTarget = 0,
  // Equilibrium efficiency metrics (Pythia model)
  equilibriumEfficiency = null,
  equilibriumScore = null,
  equilibriumLosses = null,
  // Piecewise logistic scores (0-1 scale, higher = better)
  piecewiseScores = null,
  compositeScore = null,
  qualityTier = null,
}) {
  return {
    penalty,
    penaltyTm,
    penaltyTmDiff,
    penaltyGc,
    penaltyLen,
    penaltyDg,
    penaltyOffTarget,
    // Equilibrium metrics
    equilibriumEfficiency,
    equilibriumScore,
    equilibriumLosses,
    // Piecewise logistic scores
    piecewiseScores,
    compositeScore,
    qualityTier,
  };
}

/**
 * Create a Primer object
 * @param {Object} params - Primer parameters
 * @returns {Object} Primer object
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
}) {
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
    get penalty() {
      return this.scoring.penalty;
    },
    dict() {
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
          // Equilibrium efficiency metrics (Pythia model)
          equilibrium_efficiency: this.scoring.equilibriumEfficiency,
          equilibrium_score: this.scoring.equilibriumScore,
          equilibrium_losses: this.scoring.equilibriumLosses,
          // Piecewise logistic scores (0-1 scale)
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
 * @param {Object} params - Factory parameters
 * @returns {Object} PrimerFactory object
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
}) {
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
    build(seq, tm, tmTotal, gc, dg, fwd, offTargetCount) {
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
    buildPair(fwd, rev) {
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
    withOptimalLen(newOptimalLen) {
      return createPrimerFactory({
        ...this,
        optimalLen: newOptimalLen,
      });
    },

    /**
     * Add equilibrium efficiency scoring to a primer pair.
     *
     * This uses the Pythia thermodynamic model to calculate the fraction
     * of primer that will bind to the intended target vs competing species:
     * - Hairpin (self-folding)
     * - Homodimer (self-dimer)
     * - Heterodimer (cross-dimer between fwd and rev)
     * - Off-target binding
     *
     * @param {Object} fwd - Forward primer object
     * @param {Object} rev - Reverse primer object
     * @param {string} template - Template sequence for binding site context
     * @param {number} temperature - Annealing temperature (default: 55°C)
     * @returns {[Object, Object]} Primer pair with equilibrium metrics added
     */
    addEquilibriumScoring(fwd, rev, template, temperature = 55) {
      if (!template) {
        // Without template, we can only calculate partial equilibrium
        // (hairpin, homodimer, heterodimer but no target/off-target)
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

      // Create primer objects for equilibrium calculation
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

      // Calculate full equilibrium efficiency
      const equilibrium = calculateEquilibriumEfficiency(fwdInput, revInput, template, {
        temperature,
        includeOffTarget: true,
      });

      const equilibriumScore = efficiencyToScore(equilibrium.efficiency);

      // Create updated primers with equilibrium metrics
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
     *
     * This computes individual feature scores using biologically-meaningful
     * piecewise logistic functions, then combines them into a composite score.
     *
     * @param {Object} fwd - Forward primer object
     * @param {Object} rev - Reverse primer object (optional)
     * @param {number} temperature - Annealing temperature (default: 55°C)
     * @param {string} template - Template sequence for advanced off-target classification (optional)
     * @returns {[Object, Object|null]} Primer pair with piecewise scores added
     */
    addPiecewiseScoring(fwd, rev = null, temperature = 55, template = null) {
      // Use amplification preset for consistent thresholds
      const preset = ANALYSIS_PRESETS.amplification;

      // Calculate individual scores for forward primer
      const hairpinDGFwd = calculateHairpinDG(fwd.seq, temperature);
      const homodimerDGFwd = calculateHomodimerDG(fwd.seq, temperature);

      // Use advanced off-target classification if template is provided
      let offTargetAnalysis = null;
      if (template && rev) {
        try {
          offTargetAnalysis = analyzeOffTargetsPair(fwd.seq, rev.seq, template, { temperature });
        } catch (e) {
          // Fall back to simple count-based scoring
          offTargetAnalysis = null;
        }
      }

      // Analyze G-Quadruplex risk for forward primer
      const fwdG4Analysis = analyzeGQuadruplex(fwd.seq);
      const fwdTerminalDG = calculate3primeTerminalDG(fwd.seq).dG;

      const fwdScores = {
        tm: scoreTm(fwd.tm, preset.tmOptions),
        gc: scoreGc(fwd.gc, preset.gcOptions),
        length: scoreLength(fwd.seq.length, preset.lengthOptions),
        gcClamp: scoreGcClamp(fwd.seq),
        homopolymer: scoreHomopolymer(fwd.seq),
        hairpin: scoreHairpin(hairpinDGFwd, { threshold: preset.hairpinThreshold }),
        homodimer: scoreHomodimer(homodimerDGFwd, { threshold: preset.homodimerThreshold }),
        // Use classification-based scoring if available, otherwise fall back to count-based
        offTarget: offTargetAnalysis
          ? scoreOffTargetClassification(offTargetAnalysis.fwd)
          : scoreOffTarget(fwd.offTargetCount),
        terminal3DG: scoreTerminal3DG(fwdTerminalDG),
        // G-Quadruplex scoring - critical for Q5 polymerase
        gQuadruplex: fwdG4Analysis.score,
        // 3' end composition score (4% weight - was missing!)
        threePrimeComp: score3PrimeComposition(fwd.seq, fwdTerminalDG),
      };

      // Calculate composite score for forward
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
          gQuadruplex: fwdG4Analysis,  // Include detailed G4 analysis for UI
        }),
      });

      if (!rev) {
        return [newFwd, null];
      }

      // Calculate individual scores for reverse primer
      const hairpinDGRev = calculateHairpinDG(rev.seq, temperature);
      const homodimerDGRev = calculateHomodimerDG(rev.seq, temperature);
      const heterodimerDG = calculateHeterodimerDG(fwd.seq, rev.seq, temperature);

      // Analyze G-Quadruplex risk for reverse primer
      const revG4Analysis = analyzeGQuadruplex(rev.seq);
      const revTerminalDG = calculate3primeTerminalDG(rev.seq).dG;

      const revScores = {
        tm: scoreTm(rev.tm, preset.tmOptions),
        gc: scoreGc(rev.gc, preset.gcOptions),
        length: scoreLength(rev.seq.length, preset.lengthOptions),
        gcClamp: scoreGcClamp(rev.seq),
        homopolymer: scoreHomopolymer(rev.seq),
        hairpin: scoreHairpin(hairpinDGRev, { threshold: preset.hairpinThreshold }),
        homodimer: scoreHomodimer(homodimerDGRev, { threshold: preset.homodimerThreshold }),
        // Use classification-based scoring if available
        offTarget: offTargetAnalysis
          ? scoreOffTargetClassification(offTargetAnalysis.rev)
          : scoreOffTarget(rev.offTargetCount),
        terminal3DG: scoreTerminal3DG(revTerminalDG),
        heterodimer: scoreHeterodimer(heterodimerDG, { threshold: preset.heterodimerThreshold }),
        tmDiff: scoreTmDiff(fwd.tm, rev.tm),
        // G-Quadruplex scoring - critical for Q5 polymerase
        gQuadruplex: revG4Analysis.score,
        // 3' end composition score (4% weight - was missing!)
        threePrimeComp: score3PrimeComposition(rev.seq, revTerminalDG),
      };

      // Add heterodimer and tmDiff to forward scores for completeness
      fwdScores.heterodimer = scoreHeterodimer(heterodimerDG, { threshold: preset.heterodimerThreshold });
      fwdScores.tmDiff = scoreTmDiff(fwd.tm, rev.tm);

      // Calculate pair composite score (using both primers)
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
        offTarget: Math.min(fwdScores.offTarget, revScores.offTarget), // Use worst case
        terminal3DG: Math.min(fwdScores.terminal3DG, revScores.terminal3DG),
        // G-Quadruplex scoring - critical for Q5 polymerase
        gQuadruplexFwd: fwdScores.gQuadruplex,
        gQuadruplexRev: revScores.gQuadruplex,
        // 3' end composition scores (8% total weight - was missing!)
        threePrimeCompFwd: fwdScores.threePrimeComp,
        threePrimeCompRev: revScores.threePrimeComp,
      });

      // Detect critical warnings for effective score calculation
      const criticalWarnings = [];
      const lengthLow = preset.lengthOptions?.optimalLow ?? 18;
      const lengthHigh = preset.lengthOptions?.optimalHigh ?? 24;

      // Check forward primer for critical issues
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

      // Check reverse primer for critical issues
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

      // Check Tm difference (critical if >8°C)
      const tmDiff = Math.abs(fwd.tm - rev.tm);
      if (tmDiff > 8) {
        criticalWarnings.push(`Tm difference: ${tmDiff.toFixed(1)}°C`);
      }

      // Calculate effective score with critical warning penalties (-20 per warning)
      const effectiveScore = Math.max(0, pairComposite.score - criticalWarnings.length * 20);
      const quality = classifyQuality(effectiveScore);

      // Update forward primer with pair-level composite
      const finalFwd = createPrimer({
        ...fwd,
        scoring: createScoring({
          ...fwd.scoring,
          piecewiseScores: fwdScores,
          compositeScore: pairComposite.score,
          effectiveScore,
          criticalWarnings: criticalWarnings.length,
          qualityTier: quality.tier,
          gQuadruplex: fwdG4Analysis,  // Include detailed G4 analysis for UI
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
          gQuadruplex: revG4Analysis,  // Include detailed G4 analysis for UI
        }),
      });

      return [finalFwd, newRev];
    },
  };
}

/**
 * Score primers from their sequence.
 *
 * Use-case: you already have primers and want to know their characteristics and scoring.
 *
 * @param {string} fwd - The sequence of the first primer
 * @param {string} rev - Optional sequence of the second primer
 * @param {string} seq - The sequence that the primers anneal to/amplify
 * @param {string} offtargetCheck - The sequence to check for offtarget binding sites
 * @param {Object} options - Scoring options
 * @returns {[Object, Object|null]} Primers for scoring
 */
export function score(
  fwd,
  rev = "",
  seq = "",
  offtargetCheck = "",
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
    // Equilibrium efficiency options (Pythia model)
    includeEquilibrium = true,
    annealingTemperature = 55,
  } = {}
) {
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
    // Add piecewise scoring for single primer
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

  // Add equilibrium efficiency scoring if enabled
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

  // Add piecewise logistic scoring (always enabled for comprehensive scoring)
  // Pass template for advanced off-target classification
  [resultFwd, resultRev] = factory.addPiecewiseScoring(
    resultFwd,
    resultRev,
    annealingTemperature,
    parsedOfftargetCheck || parsedSeq  // Template for off-target classification
  );

  return [resultFwd, resultRev];
}

/**
 * Attempt to find the binding region between the fwd and rev primers in seq
 */
function bindingSeq(fwd, rev = "", seq = "") {
  if (!seq) {
    return [0, "", 0];
  }

  fwd = fwd.toUpperCase();
  rev = rev.toUpperCase();
  seq = seq.toUpperCase();
  seq = seq + seq + seq; // account for amplifications across the zero-index

  let addFwd = 0;
  let addRev = 0;

  // remove start of seq prior to binding site
  try {
    let startIndexFwd = -10; // index on primer
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
  } catch (err) {
    throw new Error(`failed to find \`fwd\` binding site in \`seq\`: ${err.message}`);
  }

  if (!rev) {
    return [addFwd, seq, addRev];
  }

  // remove end of seq after the reverse primer binding site
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
  } catch (err) {
    throw new Error(`failed to find \`rev\` binding site in \`seq\`: ${err.message}`);
  }

  return [addFwd, seq, addRev];
}

/**
 * Create primers for PCR amplification of the sequence.
 *
 * @param {string} seq - The DNA sequence to amplify
 * @param {Object} options - Primer creation options
 * @param {boolean} options.useCompositeScore - Use calibrated composite score for ranking (slower but more accurate)
 * @param {boolean} options.useSmartDesign - Enable iterative smart design optimization (GC clamp, 3' end balancing)
 * @param {number} options.smartDesignTargetScore - Target score for smart design optimization (default: 75)
 * @param {boolean} options.allowExtendedPrimers - Allow primers up to 60bp for better Tm matching (default: true)
 * @returns {[Object, Object]} Primers for PCR amplification
 */
export function primers(
  seq,
  {
    addFwd = "",
    addRev = "",
    addFwdLen = [-1, -1],
    addRevLen = [-1, -1],
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
    useCompositeScore = false,  // Use calibrated composite score for ranking
    useSmartDesign = false,     // Enable iterative smart design optimization
    smartDesignTargetScore = 75, // Target score for smart design
    allowExtendedPrimers = true, // Allow primers up to 60bp for better Tm matching
    exhaustiveSearch = true,    // Exhaustive search by default: find optimal result (slightly slower but better quality)
  } = {}
) {
  // parse input
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

  // set min/max if additional sequence was provided at FWD/REV
  const [addFwdMin, addFwdMax] = parseAddLen(parsedAddFwd, addFwdLen);
  const [addRevMin, addRevMax] = parseAddLen(parsedAddRev, addRevLen);

  // Determine max primer length based on allowExtendedPrimers option
  // Extended mode (default) allows up to 60bp for better Tm matching in challenging sequences
  const effectiveMaxLen = allowExtendedPrimers ? LEN_MAX_EXTENDED : LEN_MAX;

  // create the template sequence
  const trimmedAddFwd = parsedAddFwd.slice(-addFwdMax);
  const trimmedAddRev = reverseComplement(parsedAddRev).slice(0, addRevMax);
  const seqFull = trimmedAddFwd + parsedSeq + trimmedAddRev;

  if (seqFull.length < LEN_MAX) {
    throw new Error(
      `Template sequence length is too short: ${seqFull.length}bp < ${LEN_MAX}bp`
    );
  }

  // create two 2D arrays of primers in the FWD and REV directions
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

  // Choose the best pair (by penalty or composite score)
  // With exhaustive search, also collect all candidates for diverse alternatives
  const template = parsedOfftargetCheck || parsedSeq;
  const chooseBestResult = chooseBest(factory, fwdPrimers, revPrimers, {
    useCompositeScore,
    annealingTemp: optimalTm,
    template,
    exhaustiveSearch,
    returnCandidates: exhaustiveSearch,  // Collect candidates for alternatives in exhaustive mode
  });

  // Handle both return formats (array or object with candidates)
  let minFwd, minRev, allCandidates;
  if (Array.isArray(chooseBestResult)) {
    [minFwd, minRev] = chooseBestResult;
    allCandidates = null;
  } else {
    [minFwd, minRev] = chooseBestResult.best;
    allCandidates = chooseBestResult.allCandidates;
  }

  // Add piecewise logistic scoring to design output (calibrated weights from Döring dataset)
  // Pass template for advanced off-target classification (Type A-F analysis)
  let [fwdWithScoring, revWithScoring] = factory.addPiecewiseScoring(
    minFwd,
    minRev,
    optimalTm,  // Use as annealing temperature approximation
    template    // Template for off-target classification
  );

  // Apply smart design optimization if enabled
  if (useSmartDesign) {
    const smartResult = applySmartDesign(
      fwdWithScoring,
      revWithScoring,
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

    // If Smart Design reverted due to Tm, attach alternative with higher score
    if (smartResult.alternativeOptimized) {
      fwdWithScoring.alternativeOptimized = smartResult.alternativeOptimized;
    }
  }

  // Generate alternatives - use diverse alternatives from exhaustive search if available
  try {
    if (allCandidates && allCandidates.length > 0) {
      // Use pre-evaluated candidates from exhaustive search (more efficient)
      const currentBest = { fwd: fwdWithScoring, rev: revWithScoring };
      const diverseResult = generateDiverseAlternatives(allCandidates, currentBest);

      if (diverseResult.alternatives.length > 0) {
        fwdWithScoring.alternatives = diverseResult.alternatives;
        fwdWithScoring.alternativeCategories = diverseResult.categories;
        console.log(`[primers] Generated ${diverseResult.alternatives.length} diverse alternatives from exhaustive search`);
      }
    } else {
      // Fallback: generate alternatives using separate pass (standard mode)
      const alternatives = generateAlternativesInternal(
        fwdPrimers,
        revPrimers,
        factory,
        template,
        optimalTm,
        { numAlternatives: 5, currentFwd: fwdWithScoring, currentRev: revWithScoring }
      );

      if (alternatives.length > 0) {
        fwdWithScoring.alternatives = alternatives;
        console.log(`[primers] Generated ${alternatives.length} alternative pairs`);
      }
    }
  } catch (err) {
    console.warn(`[primers] Failed to generate alternatives: ${err.message}`);
  }

  return [fwdWithScoring, revWithScoring];
}

/**
 * Apply smart design optimization to improve primer pair quality.
 *
 * This function iteratively adjusts primer length to achieve:
 * - GC clamp (G/C in last 2 bases)
 * - Optimal 3' end stability (ΔG in -6 to -11 kcal/mol range)
 * - Balanced 3' end composition
 *
 * @param {Object} fwd - Initial forward primer
 * @param {Object} rev - Initial reverse primer
 * @param {string} fullSeq - Full template sequence
 * @param {string} template - Off-target check template
 * @param {Object} options - Optimization options
 * @returns {Object} Optimized primer pair with design metadata
 */
function applySmartDesign(fwd, rev, fullSeq, template, options = {}) {
  const {
    targetScore = 75,
    temperature = 55,
    maxIterations = 3,
    addFwdMax = 0,
    addRevMax = 0,
    verbose = true,  // Enable console feedback by default
  } = options;

  // Calculate current pair score
  const currentScore = fwd.scoring.compositeScore;

  if (verbose) {
    console.log('\n[Smart Design] Analyzing primer pair...');
    console.log(`  Current score: ${currentScore}/100 (target: ${targetScore})`);
  }

  // Analyze current 3' ends (always analyze for feedback)
  const fwdAnalysis = analyze3PrimeEnd(fwd.seq);
  const revAnalysis = analyze3PrimeEnd(rev.seq);

  // Determine if 3' ends need optimization - poor/marginal quality should be optimized
  const fwd3PrimeNeedsWork = fwdAnalysis.quality === 'poor' || fwdAnalysis.quality === 'marginal' || fwdAnalysis.gcCounts.last2 === 0;
  const rev3PrimeNeedsWork = revAnalysis.quality === 'poor' || revAnalysis.quality === 'marginal' || revAnalysis.gcCounts.last2 === 0;
  const has3PrimeIssues = fwd3PrimeNeedsWork || rev3PrimeNeedsWork;

  // Only skip optimization if score meets target AND 3' ends are acceptable
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

  // Continue optimization if 3' quality is poor, even if score meets target
  if (verbose && currentScore >= targetScore && has3PrimeIssues) {
    console.log(`  ⚠ Score ${currentScore} meets target, but 3' end quality needs improvement`);
  }

  // Determine which primers need optimization (includes good quality that could be excellent)
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
  const improvements = [];

  // Optimize forward primer if needed
  if (fwdNeedsOptim) {
    if (verbose) console.log('  Optimizing forward primer...');
    const fwdStartPos = addFwdMax;  // Forward primer starts at position addFwdMax
    const fwdResult = optimizePrimer(fullSeq, fwdStartPos, fwd.seq.length, true, {
      targetGcClamp: true,
      target3PrimeDG: true,
      maxLengthChange: 3,
    });

    if (fwdResult.optimized && fwdResult.seq !== fwd.seq) {
      // Create updated primer with new sequence
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

  // Optimize reverse primer if needed
  if (revNeedsOptim) {
    if (verbose) console.log('  Optimizing reverse primer...');
    // Reverse primer is on reverse complement, calculate its position
    const revStartPos = fullSeq.length - addRevMax - rev.seq.length;
    const revFullSeq = reverseComplement(fullSeq);
    const revResult = optimizePrimer(revFullSeq, addRevMax, rev.seq.length, false, {
      targetGcClamp: true,
      target3PrimeDG: true,
      maxLengthChange: 3,
    });

    if (revResult.optimized && revResult.seq !== rev.seq) {
      // Create updated primer with new sequence
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

  // SAFETY CHECK: Reject optimization if it creates G-quadruplex or worse hairpin
  // This prevents optimizations that improve 3' ends but create secondary structure issues
  const originalFwdG4 = analyzeGQuadruplex(fwd.seq);
  const originalRevG4 = analyzeGQuadruplex(rev.seq);
  const optimizedFwdG4 = analyzeGQuadruplex(optimizedFwd.seq);
  const optimizedRevG4 = analyzeGQuadruplex(optimizedRev.seq);

  // Check if G-quadruplex risk increased
  const fwdG4Worsened = !originalFwdG4.hasG4Motif && optimizedFwdG4.hasG4Motif;
  const revG4Worsened = !originalRevG4.hasG4Motif && optimizedRevG4.hasG4Motif;
  const g4Worsened = fwdG4Worsened || revG4Worsened;

  // Check hairpin stability - reject if it got significantly worse (>3 kcal/mol more negative)
  const originalFwdHairpin = calculateHairpinDG(fwd.seq, temperature);
  const originalRevHairpin = calculateHairpinDG(rev.seq, temperature);
  const optimizedFwdHairpin = calculateHairpinDG(optimizedFwd.seq, temperature);
  const optimizedRevHairpin = calculateHairpinDG(optimizedRev.seq, temperature);

  const hairpinWorsened = (
    (optimizedFwdHairpin < originalFwdHairpin - 3) ||
    (optimizedRevHairpin < originalRevHairpin - 3) ||
    (optimizedFwdHairpin < -6 && originalFwdHairpin >= -6) ||  // Crossed critical threshold
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
    // Revert to original primers
    optimizedFwd = fwd;
    optimizedRev = rev;
    improvements.length = 0;  // Clear improvements
  }

  // Calculate new PAIR composite score (not just average of individual scores!)
  // This properly accounts for heterodimer, tmDiff, and other pair-level features
  const newPairScores = calculatePairCompositeScore(optimizedFwd, optimizedRev, temperature);
  const newScore = newPairScores.score;
  const scoreImproved = newScore > currentScore;
  // NEVER accept a score drop - the composite score already accounts for GC clamp, hairpin, etc.
  // If optimization lowers the score, it means it introduced other issues (worse heterodimer, etc.)
  const scoreWorsened = newScore < currentScore;

  // CRITICAL FIX: Update primers' compositeScore to use the PAIR score, not individual scores
  // This ensures the UI displays the correct pair score (e.g., 79) instead of individual scores (e.g., 39)
  const quality = classifyQuality(newScore);
  optimizedFwd.scoring.compositeScore = newScore;
  optimizedFwd.scoring.qualityTier = quality.tier;
  optimizedRev.scoring.compositeScore = newScore;
  optimizedRev.scoring.qualityTier = quality.tier;

  // CRITICAL: Check if Tm difference got worse - don't sacrifice Tm matching for GC clamp!
  const originalTmDiff = Math.abs(fwd.tm - rev.tm);
  const newTmDiff = Math.abs(optimizedFwd.tm - optimizedRev.tm);
  const tmDiffWorsened = newTmDiff > originalTmDiff + 2; // Allow up to 2°C degradation

  if (verbose) {
    console.log(`  Tm diff: ${originalTmDiff.toFixed(1)}°C → ${newTmDiff.toFixed(1)}°C`);
  }

  // If score got significantly worse OR Tm difference got significantly worse, revert
  // BUT provide the optimized version as an alternative for user to choose
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
      // Provide the optimized version as alternative so UI can show both options
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
 * This includes pair-level features like heterodimer, tmDiff, and cross-primer interactions
 * that are missing from individual primer scores.
 *
 * @param {Object} fwd - Forward primer object
 * @param {Object} rev - Reverse primer object
 * @param {number} temperature - Annealing temperature (default: 55°C)
 * @returns {Object} { score: number, breakdown: Object }
 */
function calculatePairCompositeScore(fwd, rev, temperature = 55) {
  // Helper to validate and sanitize numeric values
  const sanitize = (val, fallback = 0) => {
    return Number.isFinite(val) ? val : fallback;
  };

  // Calculate pair-level thermodynamics
  const heterodimerDG = sanitize(calculateHeterodimerDG(fwd.seq, rev.seq, temperature), 0);
  const hairpinDGFwd = sanitize(calculateHairpinDG(fwd.seq, temperature), 0);
  const hairpinDGRev = sanitize(calculateHairpinDG(rev.seq, temperature), 0);
  const homodimerDGFwd = sanitize(calculateHomodimerDG(fwd.seq, temperature), 0);
  const homodimerDGRev = sanitize(calculateHomodimerDG(rev.seq, temperature), 0);

  // Analyze G-Quadruplex risk
  const fwdG4 = analyzeGQuadruplex(fwd.seq);
  const revG4 = analyzeGQuadruplex(rev.seq);

  // Get terminal DG values - properly calculate from sequence if not available
  // Note: fwd.dg is the 3' terminal ΔG in kcal/mol, NOT a 0-1 score
  let fwdTerminalDG = typeof fwd.dg === 'number' ? fwd.dg : null;
  let revTerminalDG = typeof rev.dg === 'number' ? rev.dg : null;

  // Calculate if not available (don't use the score value as fallback!)
  if (fwdTerminalDG === null) {
    const calc = calculate3primeTerminalDG(fwd.seq);
    fwdTerminalDG = calc?.dG ?? -8;
  }
  if (revTerminalDG === null) {
    const calc = calculate3primeTerminalDG(rev.seq);
    revTerminalDG = calc?.dG ?? -8;
  }

  // Calculate 3' end composition scores (8% total weight in DEFAULT_WEIGHTS)
  const threePrimeCompFwd = score3PrimeComposition(fwd.seq, fwdTerminalDG);
  const threePrimeCompRev = score3PrimeComposition(rev.seq, revTerminalDG);

  // Build the complete pair scores object
  const pairScores = {
    tmFwd: sanitize(scoreTm(fwd.tm), 0.5),
    tmRev: sanitize(scoreTm(rev.tm), 0.5),
    gcFwd: sanitize(scoreGc(fwd.gc), 0.5),
    gcRev: sanitize(scoreGc(rev.gc), 0.5),
    lengthFwd: sanitize(scoreLength(fwd.seq.length), 0.5),
    lengthRev: sanitize(scoreLength(rev.seq.length), 0.5),
    gcClampFwd: sanitize(scoreGcClamp(fwd.seq), 0.5),
    gcClampRev: sanitize(scoreGcClamp(rev.seq), 0.5),
    homopolymerFwd: sanitize(scoreHomopolymer(fwd.seq), 0.5),
    homopolymerRev: sanitize(scoreHomopolymer(rev.seq), 0.5),
    hairpinFwd: sanitize(scoreHairpin(hairpinDGFwd), 0.5),
    hairpinRev: sanitize(scoreHairpin(hairpinDGRev), 0.5),
    selfDimerFwd: sanitize(scoreHomodimer(homodimerDGFwd), 0.5),
    selfDimerRev: sanitize(scoreHomodimer(homodimerDGRev), 0.5),
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
    gQuadruplexFwd: sanitize(fwdG4.score, 1),
    gQuadruplexRev: sanitize(revG4.score, 1),
    // 3' end composition scores (8% total weight - was missing!)
    threePrimeCompFwd: sanitize(threePrimeCompFwd, 0.5),
    threePrimeCompRev: sanitize(threePrimeCompRev, 0.5),
  };

  const result = calculateCompositeScore(pairScores);
  return {
    score: result.score,
    breakdown: pairScores,
  };
}

/**
 * Create an optimized primer object from optimization result
 */
function createOptimizedPrimer(optimResult, originalPrimer, isFwd, fullSeq, template, temperature) {
  // Score the new primer
  const newScoring = scorePrimerVariant(optimResult.seq, template, isFwd, temperature);

  return createPrimer({
    seq: optimResult.seq,
    len: optimResult.seq.length,
    tm: newScoring.tm,
    tmTotal: newScoring.tm,
    gc: newScoring.gc,
    dg: newScoring.terminalDG,
    fwd: isFwd,
    offTargetCount: originalPrimer.offTargetCount,  // Keep original (could re-check)
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
function addSmartDesignMetadata(primer, metadata) {
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
 *
 * This algorithm optimizes primer pairs holistically rather than individually.
 * It considers:
 * 1. Overall pair composite score
 * 2. Tm difference between primers (critical for PCR success)
 * 3. Balance of 3' end quality between primers
 * 4. Maximum acceptable sacrifice on any single primer
 *
 * The algorithm uses a grid search over length combinations and applies
 * Pareto optimization to find non-dominated solutions.
 *
 * @param {string} seq - Template sequence
 * @param {Object} currentFwd - Current forward primer
 * @param {Object} currentRev - Current reverse primer
 * @param {Object} options - Optimization options
 * @returns {Object} Best primer pair with trade-off analysis
 */
export function optimizePairHolistically(seq, currentFwd, currentRev, options = {}) {
  const {
    maxLengthDelta = 3,      // Max bp to add/remove from each primer
    maxSinglePrimerDrop = 5, // Max score drop allowed for single primer (trade-off bound)
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

  // Generate candidate lengths
  const fwdLengths = [];
  const revLengths = [];

  for (let delta = -maxLengthDelta; delta <= maxLengthDelta; delta++) {
    const fwdLen = currentFwd.seq.length + delta;
    const revLen = currentRev.seq.length + delta;
    if (fwdLen >= LEN_MIN && fwdLen <= LEN_MAX) fwdLengths.push(fwdLen);
    if (revLen >= LEN_MIN && revLen <= LEN_MAX) revLengths.push(revLen);
  }

  // Evaluate all combinations
  const candidates = [];

  for (const fwdLen of fwdLengths) {
    for (const revLen of revLengths) {
      // Generate primer at this length (simplified - actual implementation would
      // extract from sequence at correct position)
      const fwdSeq = seq.slice(0, fwdLen);
      const revSeq = reverseComplement(seq).slice(0, revLen);

      // Score each primer
      const fwdTm = calcTm(fwdSeq);
      const revTm = calcTm(revSeq);
      const tmDiff = Math.abs(fwdTm - revTm);

      // Calculate 3' end quality
      const fwd3Prime = analyze3PrimeEnd(fwdSeq);
      const rev3Prime = analyze3PrimeEnd(revSeq);

      // Pair-level scoring
      const tmDiffPenalty = tmDiff > 5 ? (tmDiff - 5) * 2 : 0;
      const qualityBalance = Math.abs(
        (fwd3Prime.quality === 'excellent' ? 3 : fwd3Prime.quality === 'good' ? 2 : 1) -
        (rev3Prime.quality === 'excellent' ? 3 : rev3Prime.quality === 'good' ? 2 : 1)
      );

      // Estimate pair score (simplified)
      const estimatedPairScore = 140 - tmDiffPenalty - qualityBalance * 3;

      candidates.push({
        fwdLen,
        revLen,
        tmDiff,
        fwd3PrimeQuality: fwd3Prime.quality,
        rev3PrimeQuality: rev3Prime.quality,
        estimatedPairScore,
        balanceScore: 4 - qualityBalance, // Higher is better balanced
      });
    }
  }

  // Sort by estimated pair score
  candidates.sort((a, b) => b.estimatedPairScore - a.estimatedPairScore);

  // Find Pareto-optimal candidates (non-dominated solutions)
  const paretoOptimal = candidates.filter((c, idx) => {
    // A candidate is Pareto-optimal if no other candidate is better in ALL objectives
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
  factory,
  seq,
  offtargetCheck,
  startRange,
  endRange,
  fwd,
  addLen
) {
  const gc = gcCache(seq);
  const tm = tmCache(seq);
  const dg = dgCache(seq);
  const ot = offTargets(seq, offtargetCheck);

  const ps = [];
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
 * Same tier system as mutagenesis for consistency.
 *
 * Tier 1 (Excellent): Tm diff ≤2°C AND structure ΔG > -3
 * Tier 2 (Good): Tm diff ≤5°C AND structure ΔG > -5
 * Tier 3 (Acceptable): Tm diff ≤8°C
 * Tier 4 (Poor): Everything else
 *
 * @param {Array} candidates - Array of { fwd, rev, score, penalty } objects
 * @returns {Object} { tier1: [], tier2: [], tier3: [], tier4: [] }
 */
function assignPrimerTiers(candidates, options = {}) {
  const {
    // Minimum composite score for tier 1 (0-100 scale)
    // If score is below this, downgrade from tier 1 to tier 2
    minScoreForTier1 = 70,
  } = options;

  const tiers = {
    tier1: [], // Excellent
    tier2: [], // Good
    tier3: [], // Acceptable
    tier4: [], // Poor
  };

  for (const c of candidates) {
    const tmDiff = Math.abs(c.fwd.tm - c.rev.tm);
    const fwdDg = c.fwd.dg || 0;
    const revDg = c.rev.dg || 0;
    const worstDg = Math.min(fwdDg, revDg);

    // Check composite score if available (score > 0 means it was calculated)
    const hasCompositeScore = c.score && c.score > 0;
    const scoreTooLow = hasCompositeScore && c.score < minScoreForTier1;

    if (tmDiff <= 2 && worstDg > -3 && !scoreTooLow) {
      tiers.tier1.push(c);
    } else if (tmDiff <= 5 && worstDg > -5) {
      // Candidates with low score but good Tm/dG go to tier 2
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
 * Choose the best combo of primers. One in FWD direction and one in the REV direction.
 * Uses tier-based selection to prioritize Tm difference (same as mutagenesis).
 *
 * JOINT TM OPTIMIZATION: When initial top-10 selection yields only tier 4 (poor) pairs,
 * expands search to find Tm-matched pairs by considering ALL candidates grouped by Tm.
 *
 * @param {Object} factory - Primer factory
 * @param {Array} fwdPrimers - 2D array of forward primers
 * @param {Array} revPrimers - 2D array of reverse primers
 * @param {Object} options - Selection options
 * @param {boolean} options.useCompositeScore - Rank by composite score instead of penalty
 * @param {number} options.annealingTemp - Annealing temperature for composite score calculation
 * @param {string} options.template - Template sequence for off-target classification
 * @param {boolean} options.returnCandidates - Return all evaluated candidates for alternatives
 * @returns {[Object, Object]|Object} Best primer pair, or object with best and candidates
 */
function chooseBest(factory, fwdPrimers, revPrimers, options = {}) {
  const {
    useCompositeScore = false,
    annealingTemp = 55,
    template = null,
    exhaustiveSearch = false,
    returnCandidates = false,  // New option: return all candidates for alternatives
  } = options;

  // Exhaustive search limit: how many candidates per group to consider
  // Standard: 10 per group (~100 pairs per Tm bucket)
  // Exhaustive: 100 per group (~10,000 pairs per Tm bucket, but more accurate)
  const CANDIDATES_PER_GROUP = exhaustiveSearch ? 100 : 10;
  const TOP_N_INITIAL = exhaustiveSearch ? 50 : 10;

  // Collect ALL primers from both directions
  const allFwd = [];
  const allRev = [];

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

  // Sort by penalty for initial selection
  allFwd.sort((a, b) => a.scoring.penalty - b.scoring.penalty);
  allRev.sort((a, b) => a.scoring.penalty - b.scoring.penalty);

  // Get top N from each for initial attempt (more for exhaustive search)
  const topNFwd = allFwd.slice(0, TOP_N_INITIAL);
  const topNRev = allRev.slice(0, TOP_N_INITIAL);

  // Helper to build candidates from two lists
  const buildCandidates = (fwdList, revList, withScoring) => {
    const candidates = [];
    for (const fwd of fwdList) {
      for (const rev of revList) {
        if (withScoring) {
          const [scoredFwd, scoredRev] = factory.addPiecewiseScoring(
            fwd, rev, annealingTemp, template
          );
          candidates.push({
            fwd: scoredFwd,
            rev: scoredRev,
            score: scoredFwd.scoring.compositeScore,
            penalty: scoredFwd.scoring.penalty + scoredRev.scoring.penalty,
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

  // Helper to select best from tiers
  const selectFromTiers = (tiers, byScore) => {
    if (byScore) {
      for (const tier of Object.values(tiers)) {
        tier.sort((a, b) => b.score - a.score);
      }
    } else {
      for (const tier of Object.values(tiers)) {
        tier.sort((a, b) => a.penalty - b.penalty);
      }
    }

    if (tiers.tier1.length > 0) return { best: tiers.tier1[0], tier: 1 };
    if (tiers.tier2.length > 0) return { best: tiers.tier2[0], tier: 2 };
    if (tiers.tier3.length > 0) return { best: tiers.tier3[0], tier: 3 };
    if (tiers.tier4.length > 0) return { best: tiers.tier4[0], tier: 4 };
    return { best: null, tier: 0 };
  };

  // First attempt: standard top-10 selection
  const initialCandidates = buildCandidates(topNFwd, topNRev, useCompositeScore);
  const initialTiers = assignPrimerTiers(initialCandidates);
  const initialResult = selectFromTiers(initialTiers, useCompositeScore);

  // Debug logging
  const debugTmDiff = initialResult.best
    ? Math.abs(initialResult.best.fwd.tm - initialResult.best.rev.tm).toFixed(1)
    : 'N/A';
  const debugScore = initialResult.best?.score || 'N/A';
  const modeLabel = exhaustiveSearch ? '[EXHAUSTIVE]' : '[standard]';
  console.log(`[chooseBest] ${modeLabel} Initial: tier=${initialResult.tier}, TmDiff=${debugTmDiff}°C, score=${debugScore}, candidates=${initialCandidates.length}, tiers: t1=${initialTiers.tier1.length}, t2=${initialTiers.tier2.length}, t3=${initialTiers.tier3.length}, t4=${initialTiers.tier4.length}`);

  // Calculate actual Tm diff for decision making
  const initialTmDiff = initialResult.best
    ? Math.abs(initialResult.best.fwd.tm - initialResult.best.rev.tm)
    : Infinity;

  // Only return early for excellent matches (tier 1 with Tm diff ≤2°C)
  // EXCEPT in exhaustive mode, continue searching for candidates with fewer issues
  // Tier 1 now requires score >= 70 (if composite scoring is enabled)
  // For tier 2-4 or any Tm diff > 2°C, continue to joint Tm optimization
  if (initialResult.best && initialResult.tier === 1 && initialTmDiff <= 2 && !exhaustiveSearch) {
    console.log(`[chooseBest] Returning excellent tier 1 result (score=${debugScore})`);
    if (returnCandidates) {
      return {
        best: [initialResult.best.fwd, initialResult.best.rev],
        allCandidates: initialCandidates,  // Return initial candidates for alternatives
      };
    }
    return [initialResult.best.fwd, initialResult.best.rev];
  }

  // For exhaustive search with tier-1 result, still continue to find candidates with fewer structural issues
  if (exhaustiveSearch && initialResult.best && initialResult.tier === 1) {
    console.log(`[chooseBest] Exhaustive: continuing search for candidates with fewer issues...`);
  }

  // JOINT TM OPTIMIZATION: Search ALL candidates to find better Tm-matched pairs
  // This is triggered when:
  // - No tier 1 match found, OR
  // - Tm diff > 2°C (even for tier 2-3)
  // Goal: Find primers with overlapping Tm ranges that top-10 selection missed

  // Group primers by Tm ranges (2°C buckets) for efficient matching
  const groupByTm = (primers) => {
    const groups = new Map();
    for (const p of primers) {
      // Round Tm to nearest 2°C for grouping
      const tmBucket = Math.round(p.tm / 2) * 2;
      if (!groups.has(tmBucket)) {
        groups.set(tmBucket, []);
      }
      groups.get(tmBucket).push(p);
    }
    return groups;
  };

  const fwdByTm = groupByTm(allFwd);
  const revByTm = groupByTm(allRev);

  // Find overlapping Tm ranges between forward and reverse
  const fwdTms = [...fwdByTm.keys()].sort((a, b) => a - b);
  const revTms = [...revByTm.keys()].sort((a, b) => a - b);

  // Collect candidates from overlapping and nearby Tm buckets
  // MUTAGENESIS-STYLE APPROACH: Don't filter by individual penalty (which penalizes
  // deviation from 62°C). Instead, take ALL primers from overlapping Tm buckets
  // and use Tm-matching penalty to rank pairs.
  const tmMatchedCandidates = [];

  // Calculate Tm-matching penalty (mutagenesis style):
  // - Quadratic penalty for Tm difference (PRIMARY criterion)
  // - Very light penalty for extreme lengths (secondary tiebreaker only)
  // - No penalty for absolute Tm value (allows finding optimal matches in GC-rich regions)
  const calculateTmMatchPenalty = (fwd, rev) => {
    let penalty = 0;
    const tmDiff = Math.abs(fwd.tm - rev.tm);

    // Tm difference: quadratic penalty - THIS IS THE PRIMARY FACTOR
    // No free zone - every degree of mismatch counts
    // 1°C → 2, 2°C → 8, 3°C → 18, 5°C → 50
    penalty += 2.0 * Math.pow(tmDiff, 2);

    // Length penalty: very light, only as tiebreaker
    // 0.05 per bp deviation means even 10bp deviation only adds 0.5 to penalty
    // Compare to Tm: 0.5°C Tm diff adds 0.5 to penalty
    const optimalLen = 22;
    const lenPenalty = 0.05;
    penalty += Math.abs(fwd.seq.length - optimalLen) * lenPenalty;
    penalty += Math.abs(rev.seq.length - optimalLen) * lenPenalty;

    // Secondary structure penalty (keep existing dG penalty)
    if (fwd.dg < -2) penalty += Math.abs(fwd.dg) * 0.5;
    if (rev.dg < -2) penalty += Math.abs(rev.dg) * 0.5;

    return { penalty, tmDiff };
  };

  // Track evaluation statistics
  let pairsEvaluated = 0;

  // For exhaustive search: select top candidates by penalty (quality), not by length diversity
  // This ensures we evaluate the best candidates, not just a spread of lengths
  const selectTopByPenalty = (arr, limit) => {
    if (arr.length <= limit) return arr;
    // Sort by penalty (lower is better) to get the highest quality primers
    const sorted = [...arr].sort((a, b) => a.scoring.penalty - b.scoring.penalty);
    return sorted.slice(0, limit);
  };

  for (const fwdTm of fwdTms) {
    // Look for reverse primers within 8°C (4 buckets) of this forward Tm
    for (const revTm of revTms) {
      if (Math.abs(fwdTm - revTm) <= 8) {
        const fwdGroup = fwdByTm.get(fwdTm);
        const revGroup = revByTm.get(revTm);

        // For exhaustive search: select top candidates by quality (penalty)
        // This ensures we evaluate the best primers, not just length diversity
        const fwdSelected = selectTopByPenalty(fwdGroup, CANDIDATES_PER_GROUP);
        const revSelected = selectTopByPenalty(revGroup, CANDIDATES_PER_GROUP);

        for (const fwd of fwdSelected) {
          for (const rev of revSelected) {

            const { penalty: tmMatchPenalty, tmDiff } = calculateTmMatchPenalty(fwd, rev);

            // PRUNING: Skip pairs with Tm diff > 5°C (unlikely to be good)
            if (tmDiff > 5) continue;

            pairsEvaluated++;

            if (useCompositeScore) {
              const [scoredFwd, scoredRev] = factory.addPiecewiseScoring(
                fwd, rev, annealingTemp, template
              );
              const candidate = {
                fwd: scoredFwd,
                rev: scoredRev,
                score: scoredFwd.scoring.compositeScore,
                penalty: tmMatchPenalty,
                tmDiff,
              };
              tmMatchedCandidates.push(candidate);
              // No early termination in exhaustive mode - evaluate ALL candidates
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

  // If we found Tm-matched candidates, evaluate them
  // For joint Tm optimization, ALWAYS sort by Tm-matching penalty, not composite score
  // This ensures we find the best Tm match regardless of other scoring factors
  if (tmMatchedCandidates.length > 0) {
    const tmMatchedTiers = assignPrimerTiers(tmMatchedCandidates);
    // Always use penalty-based selection (false) for joint Tm optimization
    // The penalty is our Tm-matching penalty, not the original PCR penalty
    const tmMatchedResult = selectFromTiers(tmMatchedTiers, false);

    if (tmMatchedResult.best) {
      const tmMatchedTmDiff = Math.abs(tmMatchedResult.best.fwd.tm - tmMatchedResult.best.rev.tm);

      const tmMatchedScore = tmMatchedResult.best.score || 0;
      console.log(`[chooseBest] Joint Tm opt: tier=${tmMatchedResult.tier}, TmDiff=${tmMatchedTmDiff.toFixed(1)}°C, score=${tmMatchedScore}`);

      // For exhaustive search, find the absolute BEST candidate among ALL evaluated
      // This is the key improvement: evaluate ALL candidates and pick the highest score
      let bestScoreCandidate = tmMatchedResult.best;

      if (exhaustiveSearch && useCompositeScore) {
        // Comprehensive issue counting that aligns with calibrated weights
        // Weights from calibration: offTarget=0.25, terminal3DG=0.20, gQuadruplexRev=0.15
        const countIssuesComprehensive = (c) => {
          let issues = 0;
          const fwdSeq = c.fwd.seq;
          const revSeq = c.rev.seq;
          const fwdScores = c.fwd.scoring?.piecewiseScores || {};
          const revScores = c.rev.scoring?.piecewiseScores || {};

          // GC clamp check (want G/C at 3' end)
          const fwdLast = fwdSeq.slice(-1).toUpperCase();
          const revLast = revSeq.slice(-1).toUpperCase();
          if (fwdLast !== 'G' && fwdLast !== 'C') issues += 1;
          if (revLast !== 'G' && revLast !== 'C') issues += 1;

          // Hairpin ΔG check (want > -3 kcal/mol)
          const fwdDg = c.fwd.dg || 0;
          const revDg = c.rev.dg || 0;
          if (fwdDg < -3) issues += 1;
          if (revDg < -3) issues += 1;

          // Off-target check (weight 0.25 - most important!)
          if ((c.fwd.offTargetCount || 0) > 0) issues += 2; // Double weight for off-targets
          if ((c.rev.offTargetCount || 0) > 0) issues += 2;

          // Terminal 3' ΔG quality (weight 0.20 - second most important)
          // Ideal range: -6 to -11 kcal/mol, score < 0.7 indicates issue
          if ((fwdScores.terminal3DG || 1) < 0.7) issues += 1;
          if ((revScores.terminal3DG || 1) < 0.7) issues += 1;

          // Heterodimer check (weight 0.06)
          if ((fwdScores.heterodimer || 1) < 0.7) issues += 1;

          // Homodimer checks (weight 0.04 each)
          if ((fwdScores.selfDimerFwd || 1) < 0.5) issues += 1;
          if ((revScores.selfDimerRev || 1) < 0.5) issues += 1;

          // G-Quadruplex check (weight 0.15 for reverse - very important!)
          if ((revScores.gQuadruplexRev || 1) < 0.5) issues += 2; // Double weight for G4 in reverse
          if ((fwdScores.gQuadruplexFwd || 1) < 0.5) issues += 1;

          return issues;
        };

        // Find ALL candidates with acceptable Tm diff (Tier 1: ≤ 2°C)
        const tier1Candidates = tmMatchedCandidates.filter(c => {
          const cTmDiff = Math.abs(c.fwd.tm - c.rev.tm);
          return cTmDiff <= 2 && c.score > 0;
        });

        if (tier1Candidates.length > 0) {
          // EXHAUSTIVE SEARCH: Simply pick the candidate with highest composite score
          // The composite score already incorporates all the calibrated weights
          // This is the simplest and most correct approach
          tier1Candidates.sort((a, b) => {
            // Primary: highest composite score wins
            if (b.score !== a.score) return b.score - a.score;
            // Tiebreaker: fewer issues
            return countIssuesComprehensive(a) - countIssuesComprehensive(b);
          });

          const bestCandidate = tier1Candidates[0];

          // Log the selection
          const bestIssues = countIssuesComprehensive(bestCandidate);
          const currentIssues = countIssuesComprehensive(bestScoreCandidate);
          const bestTmDiff = Math.abs(bestCandidate.fwd.tm - bestCandidate.rev.tm);

          console.log(`[chooseBest] Exhaustive: evaluated ${tier1Candidates.length} tier-1 candidates`);
          console.log(`[chooseBest] Exhaustive: best score=${bestCandidate.score}, issues=${bestIssues}, TmDiff=${bestTmDiff.toFixed(1)}°C`);

          // Always use the highest-scoring tier-1 candidate in exhaustive mode
          bestScoreCandidate = bestCandidate;
        } else {
          // No tier-1 candidates, fall back to best from tier-2+
          console.log(`[chooseBest] Exhaustive: no tier-1 candidates, using best from lower tiers`);
        }
      }

      // Compare exhaustive search result with initial result
      const jointScore = bestScoreCandidate.score || 0;
      const initialScore = initialResult.best?.score || 0;
      const jointTmDiff = Math.abs(bestScoreCandidate.fwd.tm - bestScoreCandidate.rev.tm);

      // For exhaustive search: always prefer the higher-scoring candidate
      // Remove arbitrary thresholds - just pick the best
      const shouldUseJointResult = (
        // Better tier always wins
        tmMatchedResult.tier < initialResult.tier ||
        // Same tier: prefer better Tm diff (if significant improvement)
        (tmMatchedResult.tier === initialResult.tier && jointTmDiff < initialTmDiff - 1) ||
        // Same tier, comparable Tm diff: prefer higher score (no arbitrary threshold)
        (tmMatchedResult.tier === initialResult.tier && jointTmDiff <= 2 && initialTmDiff <= 2 && jointScore > initialScore) ||
        // Exhaustive mode: always use the result from comprehensive evaluation
        exhaustiveSearch
      );

      if (shouldUseJointResult) {
        console.log(`[chooseBest] Using joint result: score=${jointScore}, TmDiff=${jointTmDiff.toFixed(1)}°C`);
        if (returnCandidates) {
          return {
            best: [bestScoreCandidate.fwd, bestScoreCandidate.rev],
            allCandidates: tmMatchedCandidates,  // All evaluated candidates for alternatives
          };
        }
        return [bestScoreCandidate.fwd, bestScoreCandidate.rev];
      }
    }
  }

  // Fallback to initial result if no improvement found
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
 * Uses the same Tm-matching algorithm as chooseBest() for consistency.
 *
 * This is called internally by primers() to always provide alternatives,
 * matching how mutagenesis always returns all candidates.
 *
 * @param {Array} fwdPrimers - 2D array of forward primers
 * @param {Array} revPrimers - 2D array of reverse primers
 * @param {Object} factory - Primer factory
 * @param {string} template - Template sequence
 * @param {number} annealingTemp - Annealing temperature
 * @param {Object} options - Options including numAlternatives and current primers
 * @returns {Array} Array of alternative primer pairs
 */
function generateAlternativesInternal(fwdPrimers, revPrimers, factory, template, annealingTemp, options = {}) {
  const { numAlternatives = 5, currentFwd = null, currentRev = null } = options;

  // Flatten primer arrays
  const allFwd = [];
  const allRev = [];

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

  // Group by Tm for efficient Tm-matching (same as chooseBest)
  const groupByTm = (primers) => {
    const groups = new Map();
    for (const p of primers) {
      const tmBucket = Math.round(p.tm / 2) * 2;
      if (!groups.has(tmBucket)) groups.set(tmBucket, []);
      groups.get(tmBucket).push(p);
    }
    return groups;
  };

  const fwdByTm = groupByTm(allFwd);
  const revByTm = groupByTm(allRev);
  const fwdTms = [...fwdByTm.keys()].sort((a, b) => a - b);
  const revTms = [...revByTm.keys()].sort((a, b) => a - b);

  // Calculate Tm-matching penalty (same formula as chooseBest)
  const calculateTmMatchPenalty = (fwd, rev) => {
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

  // Select diverse primers from each Tm bucket
  const selectDiverse = (arr, limit) => {
    if (arr.length <= limit) return arr;
    const sorted = [...arr].sort((a, b) => a.seq.length - b.seq.length);
    const result = [];
    const step = (sorted.length - 1) / (limit - 1);
    for (let i = 0; i < limit; i++) {
      result.push(sorted[Math.round(i * step)]);
    }
    return result;
  };

  // Collect all Tm-matched pairs from overlapping buckets
  const allPairs = [];

  for (const fwdTm of fwdTms) {
    for (const revTm of revTms) {
      if (Math.abs(fwdTm - revTm) <= 10) {
        const fwdGroup = fwdByTm.get(fwdTm);
        const revGroup = revByTm.get(revTm);
        const fwdSelected = selectDiverse(fwdGroup, 6);
        const revSelected = selectDiverse(revGroup, 6);

        for (const fwd of fwdSelected) {
          for (const rev of revSelected) {
            const { penalty: tmMatchPenalty, tmDiff } = calculateTmMatchPenalty(fwd, rev);

            // Skip pairs with very high Tm diff
            if (tmDiff > 8) continue;

            const [scoredFwd, scoredRev] = factory.addPiecewiseScoring(fwd, rev, annealingTemp, template);

            allPairs.push({
              fwd: scoredFwd,
              rev: scoredRev,
              tmDiff,
              penalty: tmMatchPenalty,
              compositeScore: scoredFwd.scoring.compositeScore,
            });
          }
        }
      }
    }
  }

  // Sort by tier first (Tm diff), then by composite score
  allPairs.sort((a, b) => {
    // Tier comparison: tier 1 (≤2) < tier 2 (≤5) < tier 3 (≤8) < tier 4
    const tierA = a.tmDiff <= 2 ? 1 : a.tmDiff <= 5 ? 2 : a.tmDiff <= 8 ? 3 : 4;
    const tierB = b.tmDiff <= 2 ? 1 : b.tmDiff <= 5 ? 2 : b.tmDiff <= 8 ? 3 : 4;
    if (tierA !== tierB) return tierA - tierB;
    // Within same tier, prefer higher composite score
    return b.compositeScore - a.compositeScore;
  });

  // Filter out the current best pair and select unique alternatives by length
  const currentKey = currentFwd && currentRev
    ? `${currentFwd.seq.length}-${currentRev.seq.length}`
    : null;

  const seen = new Set();
  if (currentKey) seen.add(currentKey);

  const alternatives = [];

  // Quality thresholds - only include usable alternatives
  const MIN_SCORE = 50;
  const MAX_PRIMER_TM = 72;

  for (const pair of allPairs) {
    const key = `${pair.fwd.seq.length}-${pair.rev.seq.length}`;
    if (seen.has(key)) continue;

    // Skip alternatives that don't meet quality thresholds
    if (pair.compositeScore < MIN_SCORE) continue;
    if (pair.fwd.tm > MAX_PRIMER_TM || pair.rev.tm > MAX_PRIMER_TM) continue;

    seen.add(key);

    // Determine tier label based on SCORE (not Tm diff)
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
 * This is more efficient than generateAlternativesInternal when we already have
 * scored candidates from exhaustive search.
 *
 * Provides alternatives in multiple categories:
 * - Highest Score: Best overall quality
 * - Shortest Primers: Lower synthesis cost
 * - Longest Primers: Higher specificity
 * - Best Tm Match: Optimal PCR conditions
 * - Fewest Issues: Most reliable (no warnings)
 *
 * @param {Array} candidates - Pre-evaluated candidates from chooseBest
 * @param {Object} currentBest - The selected best pair {fwd, rev}
 * @param {Object} options - Options
 * @returns {Object} Diverse alternatives organized by category
 */
function generateDiverseAlternatives(candidates, currentBest, options = {}) {
  const { maxPerCategory = 3 } = options;

  if (!candidates || candidates.length === 0) {
    return { alternatives: [], categories: {} };
  }

  // Quality thresholds - alternatives must meet these to be useful
  const MIN_SCORE = 50;  // Minimum acceptable composite score
  const MAX_PRIMER_TM = 72;  // Maximum individual primer Tm for standard PCR
  const MAX_TM_DIFF = 5;  // Maximum Tm difference between primers

  // Helper to classify quality tier based on score (consistent with scoring.js)
  const getQualityTier = (score) => {
    if (score >= 90) return 'excellent';
    if (score >= 75) return 'good';
    if (score >= 60) return 'acceptable';
    return 'poor';
  };

  // Helper to format a candidate for output
  const formatCandidate = (c, category) => {
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
      qualityTier: getQualityTier(score),  // Use score-based tier, not Tm-based
      category,
    };
  };

  // Exclude the current best from alternatives
  const currentKey = currentBest
    ? `${currentBest.fwd.seq}-${currentBest.rev.seq}`
    : null;

  // Filter candidates: must have data AND meet quality thresholds
  const filtered = candidates.filter(c => {
    const key = `${c.fwd.seq}-${c.rev.seq}`;
    if (key === currentKey) return false;

    // Must have score data
    const score = c.score || c.compositeScore || 0;
    const hasData = score > 0 || c.penalty !== undefined;
    if (!hasData) return false;

    // Quality thresholds: alternatives must be usable
    if (score < MIN_SCORE) return false;
    if (c.fwd.tm > MAX_PRIMER_TM || c.rev.tm > MAX_PRIMER_TM) return false;
    if (Math.abs(c.fwd.tm - c.rev.tm) > MAX_TM_DIFF) return false;

    // Reject pairs with severe heterodimer risk (ΔG < -10 kcal/mol)
    // This is the quality gate for user-facing alternatives
    const heterodimerDG = c.heterodimerDG || 0;
    if (heterodimerDG < -10) return false;

    return true;
  });

  // Count issues for sorting
  const countIssues = (c) => {
    let issues = 0;
    const fwdSeq = c.fwd.seq;
    const revSeq = c.rev.seq;

    // GC clamp
    if (!/[GC]$/.test(fwdSeq)) issues += 1;
    if (!/[GC]$/.test(revSeq)) issues += 1;

    // Hairpin
    if ((c.fwd.dg || 0) < -3) issues += 1;
    if ((c.rev.dg || 0) < -3) issues += 1;

    // Off-target
    if ((c.fwd.offTargetCount || 0) > 0) issues += 1;
    if ((c.rev.offTargetCount || 0) > 0) issues += 1;

    // High Tm (even within 72°C limit, prefer lower Tm)
    if (c.fwd.tm > 68) issues += 1;
    if (c.rev.tm > 68) issues += 1;

    // Very long primers (synthesis concern)
    if (c.fwd.seq.length > 30) issues += 1;
    if (c.rev.seq.length > 30) issues += 1;

    // Heterodimer risk (graduated penalty based on severity)
    // ΔG < -9: moderate risk, ΔG < -12: severe risk
    const heterodimerDG = c.heterodimerDG || 0;
    if (heterodimerDG < -12) issues += 2;       // Severe heterodimer
    else if (heterodimerDG < -9) issues += 1;  // Moderate heterodimer

    return issues;
  };

  const categories = {};
  const seenKeys = new Set();
  if (currentKey) seenKeys.add(currentKey);

  // Helper to add unique candidates to a category
  const addToCategory = (name, sortedCandidates) => {
    categories[name] = [];
    for (const c of sortedCandidates) {
      const key = `${c.fwd.seq}-${c.rev.seq}`;
      if (seenKeys.has(key)) continue;
      seenKeys.add(key);
      categories[name].push(formatCandidate(c, name));
      if (categories[name].length >= maxPerCategory) break;
    }
  };

  // 1. Highest Score - sorted by composite score descending
  const byScore = [...filtered].sort((a, b) =>
    (b.score || b.compositeScore || 0) - (a.score || a.compositeScore || 0)
  );
  addToCategory('highestScore', byScore);

  // 2. Shortest Primers - for lower synthesis cost
  const byShortestLength = [...filtered].sort((a, b) => {
    const lenA = a.fwd.seq.length + a.rev.seq.length;
    const lenB = b.fwd.seq.length + b.rev.seq.length;
    if (lenA !== lenB) return lenA - lenB;
    return (b.score || 0) - (a.score || 0);  // Tiebreaker: score
  });
  addToCategory('shortestPrimers', byShortestLength);

  // 3. Longest Primers - for higher specificity
  const byLongestLength = [...filtered].sort((a, b) => {
    const lenA = a.fwd.seq.length + a.rev.seq.length;
    const lenB = b.fwd.seq.length + b.rev.seq.length;
    if (lenA !== lenB) return lenB - lenA;
    return (b.score || 0) - (a.score || 0);
  });
  addToCategory('longestPrimers', byLongestLength);

  // 4. Best Tm Match - lowest Tm difference
  const byTmMatch = [...filtered].sort((a, b) => {
    const tmDiffA = Math.abs(a.fwd.tm - a.rev.tm);
    const tmDiffB = Math.abs(b.fwd.tm - b.rev.tm);
    if (Math.abs(tmDiffA - tmDiffB) > 0.1) return tmDiffA - tmDiffB;
    return (b.score || 0) - (a.score || 0);
  });
  addToCategory('bestTmMatch', byTmMatch);

  // 5. Fewest Issues - most reliable
  const byFewestIssues = [...filtered].sort((a, b) => {
    const issuesA = countIssues(a);
    const issuesB = countIssues(b);
    if (issuesA !== issuesB) return issuesA - issuesB;
    return (b.score || 0) - (a.score || 0);
  });
  addToCategory('fewestIssues', byFewestIssues);

  // Combine all unique alternatives (flattened, deduplicated)
  const allAlternatives = [];
  const finalSeen = new Set();
  if (currentKey) finalSeen.add(currentKey);

  // Interleave from categories for variety
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
    alternatives: allAlternatives.slice(0, 10),  // Max 10 total alternatives
    categories,
  };
}

/**
 * Validate and parse the input sequence.
 */
function parse(seq, offtargetCheck) {
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
function parseAddLen(add, addLen) {
  if (!add) {
    // we're not adding anything
    return [0, 0];
  }

  const [addMin, addMax] = addLen;

  if (addMin === -1 && addMax === -1) {
    // just add the whole sequence if <40bp
    // if greater than 40 bp, log a warning and add (30, 40) bp
    const addAssume = Math.min(40, add.length);
    if (addAssume !== add.length) {
      console.warn(
        `${add.length}bp additional sequence added, but no \`addDirLen\` argument provided:\n\tAdding between ${addAssume - 10} and ${addAssume}bp`
      );
      return [addAssume - 10, addAssume];
    }
    return [addAssume, addAssume];
  }

  // one or both of them were set
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
function reverseComplement(seq) {
  const rc = { A: "T", T: "A", G: "C", C: "G" };
  return seq
    .split("")
    .reverse()
    .map((c) => rc[c])
    .join("");
}

/**
 * Create a range array
 */
function range(start, end) {
  const result = [];
  for (let i = start; i < end; i++) {
    result.push(i);
  }
  return result;
}

/**
 * Generate alternative primer pairs with different length combinations.
 *
 * Returns the top N alternative primer pairs ranked by composite score.
 *
 * @param {string} seq - The DNA sequence to amplify
 * @param {Object} options - Same options as primers() plus:
 * @param {number} options.numAlternatives - Number of alternatives to return (default: 3)
 * @returns {Array} Array of alternative primer pairs with scoring
 */
export function generateAlternatives(
  seq,
  {
    addFwd = "",
    addRev = "",
    addFwdLen = [-1, -1],
    addRevLen = [-1, -1],
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
  } = {}
) {
  // parse input
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

  // set min/max if additional sequence was provided at FWD/REV
  const [addFwdMin, addFwdMax] = parseAddLen(parsedAddFwd, addFwdLen);
  const [addRevMin, addRevMax] = parseAddLen(parsedAddRev, addRevLen);

  // create the template sequence
  const trimmedAddFwd = parsedAddFwd.slice(-addFwdMax);
  const trimmedAddRev = reverseComplement(parsedAddRev).slice(0, addRevMax);
  const seqFull = trimmedAddFwd + parsedSeq + trimmedAddRev;

  if (seqFull.length < LEN_MAX) {
    throw new Error(
      `Template sequence length is too short: ${seqFull.length}bp < ${LEN_MAX}bp`
    );
  }

  // create two 2D arrays of primers in the FWD and REV directions
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

  // Collect all valid primer pairs
  const template = parsedOfftargetCheck || parsedSeq;
  const allPairs = [];

  // Flatten primer arrays
  const flatFwd = [];
  const flatRev = [];

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

  // FULL CARTESIAN PRODUCT: Evaluate ALL pairs like mutagenesis
  // This ensures we find optimal Tm matches even in GC-rich regions
  // Group by Tm for efficient matching
  const groupByTm = (primers) => {
    const groups = new Map();
    for (const p of primers) {
      const tmBucket = Math.round(p.tm / 2) * 2;
      if (!groups.has(tmBucket)) groups.set(tmBucket, []);
      groups.get(tmBucket).push(p);
    }
    return groups;
  };

  const fwdByTm = groupByTm(flatFwd);
  const revByTm = groupByTm(flatRev);
  const fwdTms = [...fwdByTm.keys()].sort((a, b) => a - b);
  const revTms = [...revByTm.keys()].sort((a, b) => a - b);

  // Calculate Tm-matching penalty (same as chooseBest)
  const calculateTmMatchPenalty = (fwd, rev) => {
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

  // Select diverse primers from each Tm bucket
  const selectDiverse = (arr, limit) => {
    if (arr.length <= limit) return arr;
    const sorted = [...arr].sort((a, b) => a.seq.length - b.seq.length);
    const result = [];
    const step = (sorted.length - 1) / (limit - 1);
    for (let i = 0; i < limit; i++) {
      result.push(sorted[Math.round(i * step)]);
    }
    return result;
  };

  // Evaluate pairs from overlapping Tm buckets (full search)
  for (const fwdTm of fwdTms) {
    for (const revTm of revTms) {
      if (Math.abs(fwdTm - revTm) <= 10) { // Wider range for alternatives
        const fwdGroup = fwdByTm.get(fwdTm);
        const revGroup = revByTm.get(revTm);
        const fwdSelected = selectDiverse(fwdGroup, 8);
        const revSelected = selectDiverse(revGroup, 8);

        for (const fwd of fwdSelected) {
          for (const rev of revSelected) {
            const { penalty: tmMatchPenalty, tmDiff } = calculateTmMatchPenalty(fwd, rev);
            const [scoredFwd, scoredRev] = factory.addPiecewiseScoring(fwd, rev, optimalTm, template);
            allPairs.push({
              forward: scoredFwd,
              reverse: scoredRev,
              compositeScore: scoredFwd.scoring.compositeScore,
              qualityTier: scoredFwd.scoring.qualityTier,
              penalty: tmMatchPenalty,
              tmDiff: tmDiff,
            });
          }
        }
      }
    }
  }

  // Use tier-based selection to prioritize Tm difference
  const tieredCandidates = allPairs.map(pair => ({
    fwd: pair.forward,
    rev: pair.reverse,
    score: pair.compositeScore,
    penalty: pair.penalty,
    qualityTier: pair.qualityTier,
    tmDiff: pair.tmDiff,
  }));

  const tiers = assignPrimerTiers(tieredCandidates);

  // Sort each tier by composite score (highest first)
  for (const tier of Object.values(tiers)) {
    tier.sort((a, b) => b.score - a.score);
  }

  // Combine tiers in order: tier1 first, then tier2, etc.
  const sortedPairs = [
    ...tiers.tier1,
    ...tiers.tier2,
    ...tiers.tier3,
    ...tiers.tier4,
  ];

  // Return unique alternatives (different lengths)
  const seen = new Set();
  const alternatives = [];

  for (const pair of sortedPairs) {
    const key = `${pair.fwd.seq.length}-${pair.rev.seq.length}`;
    if (!seen.has(key)) {
      seen.add(key);
      const tmDiff = Math.abs(pair.fwd.tm - pair.rev.tm);
      // Determine quality tier based on Tm diff
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

  // Return all alternatives - caller will filter based on current design
  return alternatives;
}

// Alias for primers function
export const create = primers;
