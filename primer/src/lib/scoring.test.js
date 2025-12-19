/**
 * Tests for Piecewise Logistic Scoring Functions
 */

import { expect, describe, it } from "vitest";
import {
  piecewiseLogistic,
  scoreTm,
  scoreGc,
  scoreTerminal3DG,
  scoreTmDiff,
  scoreHairpin,
  scoreHomodimer,
  scoreHeterodimer,
  scoreOffTarget,
  scoreOffTargetClassification,
  scoreLength,
  scoreGcClamp,
  scoreHomopolymer,
  scoreAmpliconLength,
  scoreDistanceToROI,
  scoreAmpliconStructure,
  calculateCompositeScore,
  classifyQuality,
} from "./scoring.js";

describe("piecewiseLogistic (generic)", () => {
  const opts = {
    optimalLow: 50,
    optimalHigh: 60,
    acceptableLow: 40,
    acceptableHigh: 70,
    steepness: 0.5,
  };

  it("should return 1.0 in optimal range", () => {
    expect(piecewiseLogistic(55, opts)).toBe(1.0);
    expect(piecewiseLogistic(50, opts)).toBe(1.0);
    expect(piecewiseLogistic(60, opts)).toBe(1.0);
  });

  it("should return 0.7-1.0 in acceptable range (below optimal)", () => {
    const score = piecewiseLogistic(45, opts);
    expect(score).toBeGreaterThan(0.7);
    expect(score).toBeLessThan(1.0);

    // At acceptable boundary, should be ~0.7
    const boundaryScore = piecewiseLogistic(40, opts);
    expect(boundaryScore).toBeCloseTo(0.7, 1);
  });

  it("should return 0.7-1.0 in acceptable range (above optimal)", () => {
    const score = piecewiseLogistic(65, opts);
    expect(score).toBeGreaterThan(0.7);
    expect(score).toBeLessThan(1.0);

    const boundaryScore = piecewiseLogistic(70, opts);
    expect(boundaryScore).toBeCloseTo(0.7, 1);
  });

  it("should decay below acceptable range", () => {
    const score = piecewiseLogistic(30, opts);
    expect(score).toBeLessThan(0.7);
    expect(score).toBeGreaterThan(0);
  });

  it("should decay above acceptable range", () => {
    const score = piecewiseLogistic(80, opts);
    expect(score).toBeLessThan(0.7);
    expect(score).toBeGreaterThan(0);
  });
});

describe("scoreTm", () => {
  it("should return 1.0 for optimal Tm range (55-60°C)", () => {
    expect(scoreTm(55)).toBe(1.0);
    expect(scoreTm(57.5)).toBe(1.0);
    expect(scoreTm(60)).toBe(1.0);
  });

  it("should return high score for acceptable Tm (50-55°C, 60-65°C)", () => {
    const lowAcceptable = scoreTm(52);
    expect(lowAcceptable).toBeGreaterThan(0.7);
    expect(lowAcceptable).toBeLessThan(1.0);

    const highAcceptable = scoreTm(63);
    expect(highAcceptable).toBeGreaterThan(0.7);
    expect(highAcceptable).toBeLessThan(1.0);
  });

  it("should penalize Tm below 50°C", () => {
    const score = scoreTm(45);
    expect(score).toBeLessThan(0.7);
  });

  it("should penalize Tm above 65°C", () => {
    const score = scoreTm(70);
    expect(score).toBeLessThan(0.7);
  });

  it("should accept custom parameters", () => {
    // RT-qPCR typically uses higher optimal Tm
    const rtqpcrScore = scoreTm(62, {
      optimalLow: 58,
      optimalHigh: 62,
    });
    expect(rtqpcrScore).toBe(1.0);
  });
});

describe("scoreGc", () => {
  it("should return 1.0 for optimal GC (40-60%)", () => {
    expect(scoreGc(0.5)).toBe(1.0);  // 50% as fraction
    expect(scoreGc(50)).toBe(1.0);   // 50% as percentage
    expect(scoreGc(0.4)).toBe(1.0);  // 40%
    expect(scoreGc(0.6)).toBe(1.0);  // 60%
  });

  it("should return high score for acceptable GC (30-40%, 60-70%)", () => {
    expect(scoreGc(0.35)).toBeGreaterThan(0.7);
    expect(scoreGc(0.35)).toBeLessThan(1.0);

    expect(scoreGc(0.65)).toBeGreaterThan(0.7);
    expect(scoreGc(0.65)).toBeLessThan(1.0);
  });

  it("should penalize extreme GC", () => {
    expect(scoreGc(0.25)).toBeLessThan(0.7);  // Low GC
    expect(scoreGc(0.75)).toBeLessThan(0.7);  // High GC
  });

  it("should handle both fraction and percentage input", () => {
    const fractionScore = scoreGc(0.45);
    const percentageScore = scoreGc(45);
    expect(fractionScore).toBeCloseTo(percentageScore, 2);
  });
});

describe("scoreTerminal3DG", () => {
  it("should return 1.0 for optimal 3' ΔG (-6 to -11 kcal/mol)", () => {
    expect(scoreTerminal3DG(-8)).toBe(1.0);
    expect(scoreTerminal3DG(-6)).toBe(1.0);
    expect(scoreTerminal3DG(-11)).toBe(1.0);
  });

  it("should penalize too-loose binding (ΔG > -6)", () => {
    const score = scoreTerminal3DG(-3);  // Too weak
    expect(score).toBeLessThan(1.0);
    expect(score).toBeGreaterThan(0);
  });

  it("should penalize too-tight binding (ΔG < -11)", () => {
    const score = scoreTerminal3DG(-15);  // Too strong
    expect(score).toBeLessThan(1.0);
    expect(score).toBeGreaterThan(0);
  });

  it("should penalize too-loose more harshly than too-tight", () => {
    const tooLoose = scoreTerminal3DG(-3);  // 3 kcal/mol above optimal
    const tooTight = scoreTerminal3DG(-14); // 3 kcal/mol below optimal

    // Too loose should have harsher penalty (lower score)
    expect(tooLoose).toBeLessThan(tooTight);
  });

  it("should handle positive ΔG (very weak)", () => {
    const score = scoreTerminal3DG(0);
    expect(score).toBeLessThan(0.5);
  });
});

describe("scoreTmDiff", () => {
  it("should return 1.0 for matched Tm (0-3°C difference)", () => {
    expect(scoreTmDiff(60, 60)).toBe(1.0);
    expect(scoreTmDiff(60, 62)).toBe(1.0);
    expect(scoreTmDiff(58, 61)).toBe(1.0);
  });

  it("should have mild penalty for 3-5°C difference", () => {
    const score = scoreTmDiff(60, 64);
    expect(score).toBeGreaterThan(0.8);
    expect(score).toBeLessThan(1.0);
  });

  it("should have moderate penalty for 5-8°C difference", () => {
    const score = scoreTmDiff(60, 66);
    expect(score).toBeGreaterThan(0.5);
    expect(score).toBeLessThan(0.8);
  });

  it("should have steep penalty for >8°C difference", () => {
    const score = scoreTmDiff(60, 70);
    expect(score).toBeLessThan(0.5);
  });

  it("should be symmetric (order doesn't matter)", () => {
    expect(scoreTmDiff(60, 65)).toBe(scoreTmDiff(65, 60));
  });
});

describe("scoreHairpin", () => {
  it("should return 1.0 for weak hairpin (ΔG >= -3)", () => {
    expect(scoreHairpin(0)).toBe(1.0);
    expect(scoreHairpin(-2)).toBe(1.0);
    expect(scoreHairpin(-3)).toBe(1.0);
  });

  it("should penalize stable hairpin (ΔG < -3)", () => {
    expect(scoreHairpin(-5)).toBeLessThan(1.0);
    expect(scoreHairpin(-8)).toBeLessThan(scoreHairpin(-5));
  });
});

describe("scoreHomodimer", () => {
  it("should return 1.0 for weak homodimer (ΔG >= -6)", () => {
    expect(scoreHomodimer(0)).toBe(1.0);
    expect(scoreHomodimer(-5)).toBe(1.0);
    expect(scoreHomodimer(-6)).toBe(1.0);
  });

  it("should penalize stable homodimer (ΔG < -6)", () => {
    expect(scoreHomodimer(-8)).toBeLessThan(1.0);
    expect(scoreHomodimer(-12)).toBeLessThan(scoreHomodimer(-8));
  });
});

describe("scoreHeterodimer", () => {
  it("should return 1.0 for weak heterodimer (ΔG >= -6)", () => {
    expect(scoreHeterodimer(0)).toBe(1.0);
    expect(scoreHeterodimer(-6)).toBe(1.0);
  });

  it("should penalize stable heterodimer (ΔG < -6)", () => {
    expect(scoreHeterodimer(-10)).toBeLessThan(1.0);
  });
});

describe("scoreOffTarget", () => {
  it("should return 1.0 for no off-targets", () => {
    expect(scoreOffTarget(0)).toBe(1.0);
  });

  it("should penalize off-targets exponentially", () => {
    const one = scoreOffTarget(1);
    const two = scoreOffTarget(2);

    expect(one).toBeCloseTo(0.7, 1);
    expect(two).toBeLessThan(one);
  });

  it("should disqualify 3+ off-targets", () => {
    expect(scoreOffTarget(3)).toBe(0.0);
    expect(scoreOffTarget(5)).toBe(0.0);
  });
});

describe("scoreOffTargetClassification", () => {
  it("should return 1.0 for no off-targets", () => {
    const classification = {
      counts: { highRisk: 0, mediumRisk: 0, lowRisk: 0 }
    };
    expect(scoreOffTargetClassification(classification)).toBe(1.0);
  });

  it("should return 0.3 for 1 high-risk site", () => {
    const classification = {
      counts: { highRisk: 1, mediumRisk: 0, lowRisk: 0 }
    };
    expect(scoreOffTargetClassification(classification)).toBe(0.3);
  });

  it("should return 0.1 for 2 high-risk sites", () => {
    const classification = {
      counts: { highRisk: 2, mediumRisk: 0, lowRisk: 0 }
    };
    expect(scoreOffTargetClassification(classification)).toBe(0.1);
  });

  it("should return 0.0 for 3+ high-risk sites", () => {
    const classification = {
      counts: { highRisk: 3, mediumRisk: 0, lowRisk: 0 }
    };
    expect(scoreOffTargetClassification(classification)).toBe(0.0);
  });

  it("should reduce score for medium-risk sites", () => {
    const noRisk = { counts: { highRisk: 0, mediumRisk: 0, lowRisk: 0 } };
    const mediumRisk = { counts: { highRisk: 0, mediumRisk: 2, lowRisk: 0 } };

    expect(scoreOffTargetClassification(mediumRisk)).toBeLessThan(
      scoreOffTargetClassification(noRisk)
    );
    expect(scoreOffTargetClassification(mediumRisk)).toBeCloseTo(0.7, 1);
  });

  it("should minimally reduce score for low-risk sites", () => {
    const noRisk = { counts: { highRisk: 0, mediumRisk: 0, lowRisk: 0 } };
    const lowRisk = { counts: { highRisk: 0, mediumRisk: 0, lowRisk: 5 } };

    const diff = scoreOffTargetClassification(noRisk) - scoreOffTargetClassification(lowRisk);
    expect(diff).toBeLessThan(0.2);  // Low risk has minimal impact
  });
});

describe("scoreLength", () => {
  it("should return 1.0 for optimal length (18-24bp)", () => {
    expect(scoreLength(20)).toBe(1.0);
    expect(scoreLength(18)).toBe(1.0);
    expect(scoreLength(24)).toBe(1.0);
  });

  it("should penalize short primers (<18bp)", () => {
    expect(scoreLength(16)).toBeLessThan(1.0);
    expect(scoreLength(16)).toBeGreaterThan(0.7);
  });

  it("should penalize long primers (>24bp)", () => {
    expect(scoreLength(28)).toBeLessThan(1.0);
    expect(scoreLength(28)).toBeGreaterThan(0.7);
  });

  it("should penalize very short primers harshly", () => {
    expect(scoreLength(12)).toBeLessThan(0.7);
  });
});

describe("scoreGcClamp", () => {
  it("should return 1.0 for ideal GC clamp (1 G/C in last 2 bases)", () => {
    expect(scoreGcClamp("ATCGATCGTG")).toBe(1.0);  // ends in TG (1 G/C)
    expect(scoreGcClamp("ATCGATCGAC")).toBe(1.0);  // ends in AC (1 G/C)
    expect(scoreGcClamp("ATCGATCGCT")).toBe(1.0);  // ends in CT (1 G/C)
    expect(scoreGcClamp("ATCGATCGAG")).toBe(1.0);  // ends in AG (1 G/C)
  });

  it("should return 0.85 for strong clamp (2 G/C in last 2 bases)", () => {
    expect(scoreGcClamp("ATCGATCGGC")).toBe(0.85);  // ends in GC
    expect(scoreGcClamp("ATCGATCGCG")).toBe(0.85);  // ends in CG
  });

  it("should return 0.5 for weak clamp (0 G/C in last 2 bases)", () => {
    expect(scoreGcClamp("ATCGATCGAA")).toBe(0.5);  // ends in AA
    expect(scoreGcClamp("ATCGATCGTT")).toBe(0.5);  // ends in TT
  });
});

describe("scoreHomopolymer", () => {
  it("should return 1.0 for no runs >3", () => {
    expect(scoreHomopolymer("ATCGATCGATCG")).toBe(1.0);
    expect(scoreHomopolymer("AAATTTCCCGGG")).toBe(1.0);  // runs of 3 OK
  });

  it("should penalize runs of 4", () => {
    expect(scoreHomopolymer("ATCGAAAATCG")).toBeLessThan(1.0);
  });

  it("should penalize longer runs more harshly", () => {
    const run4 = scoreHomopolymer("ATCGAAAATCG");   // 4 A's
    const run5 = scoreHomopolymer("ATCGAAAAATCG");  // 5 A's
    const run6 = scoreHomopolymer("ATCGAAAAAATCG"); // 6 A's

    expect(run5).toBeLessThan(run4);
    expect(run6).toBeLessThan(run5);
  });
});

describe("scoreAmpliconLength", () => {
  it("should return 1.0 for optimal length (400-800bp)", () => {
    expect(scoreAmpliconLength(500)).toBe(1.0);
    expect(scoreAmpliconLength(400)).toBe(1.0);
    expect(scoreAmpliconLength(800)).toBe(1.0);
  });

  it("should penalize short amplicons (<400bp)", () => {
    expect(scoreAmpliconLength(300)).toBeLessThan(1.0);
    expect(scoreAmpliconLength(300)).toBeGreaterThan(0.7);
  });

  it("should penalize long amplicons (>800bp)", () => {
    expect(scoreAmpliconLength(1000)).toBeLessThan(1.0);
    expect(scoreAmpliconLength(1000)).toBeGreaterThan(0.7);
  });

  it("should penalize very short amplicons harshly", () => {
    expect(scoreAmpliconLength(100)).toBeLessThan(0.7);
  });
});

describe("scoreDistanceToROI (Sanger-specific)", () => {
  it("should return 1.0 for optimal distance (100-500bp)", () => {
    expect(scoreDistanceToROI(200)).toBe(1.0);
    expect(scoreDistanceToROI(100)).toBe(1.0);
    expect(scoreDistanceToROI(500)).toBe(1.0);
  });

  it("should penalize too-close distances (<100bp)", () => {
    expect(scoreDistanceToROI(75)).toBeLessThan(1.0);
    expect(scoreDistanceToROI(75)).toBeGreaterThan(0.7);
  });

  it("should penalize too-far distances (>500bp)", () => {
    expect(scoreDistanceToROI(600)).toBeLessThan(1.0);
    expect(scoreDistanceToROI(600)).toBeGreaterThan(0.7);
  });

  it("should return 0.0 for negative distance (wrong direction)", () => {
    expect(scoreDistanceToROI(-50)).toBe(0.0);
    expect(scoreDistanceToROI(-1)).toBe(0.0);
  });

  it("should penalize very close distances harshly", () => {
    expect(scoreDistanceToROI(20)).toBeLessThan(0.7);
  });

  it("should penalize very far distances harshly", () => {
    expect(scoreDistanceToROI(900)).toBeLessThan(0.7);
  });
});

describe("scoreAmpliconStructure (Sanger-specific)", () => {
  it("should return 1.0 for normal GC content", () => {
    // 50% GC, no long runs
    const normalSeq = 'ATGCATGCATGCATGCATGCATGCATGCATGC';
    expect(scoreAmpliconStructure(normalSeq)).toBe(1.0);
  });

  it("should penalize GC-rich sequences", () => {
    // 80% GC in a window
    const gcRichSeq = 'GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC';
    const score = scoreAmpliconStructure(gcRichSeq);
    expect(score).toBeLessThan(1.0);
  });

  it("should penalize long homopolymer runs", () => {
    // Contains AAAAAA (6 A's)
    const homopolymerSeq = 'ATGCATGCAAAAAATGCATGCATGCATGCATGC';
    const score = scoreAmpliconStructure(homopolymerSeq);
    expect(score).toBeLessThan(1.0);
  });

  it("should return 1.0 for short sequences", () => {
    // Sequence shorter than window size
    expect(scoreAmpliconStructure('ATGC')).toBe(1.0);
    expect(scoreAmpliconStructure('')).toBe(1.0);
  });

  it("should penalize sequences with both GC-rich regions and homopolymers", () => {
    const badSeq = 'GCGCGCGCGCGCGCGCGCGCAAAAAAGCGCGCGCGCGCGCGCGCGCGCGC';
    const score = scoreAmpliconStructure(badSeq);
    expect(score).toBeLessThan(0.6);  // Multiple penalties
  });

  it("should handle AT-rich sequences well", () => {
    const atRichSeq = 'ATATATATATATATATATATATATATATATAT';
    expect(scoreAmpliconStructure(atRichSeq)).toBe(1.0);
  });
});

describe("calculateCompositeScore", () => {
  it("should calculate weighted score from individual scores", () => {
    const scores = {
      tmFwd: 1.0,
      tmRev: 1.0,
      gcFwd: 0.8,
      gcRev: 0.8,
    };

    const result = calculateCompositeScore(scores);

    expect(result.score).toBeGreaterThan(0);
    expect(result.score).toBeLessThanOrEqual(100);
    expect(result.breakdown).toHaveProperty('tmFwd');
    expect(result.breakdown).toHaveProperty('gcFwd');
  });

  it("should return perfect score when all scores are 1.0", () => {
    const scores = {
      offTarget: 1.0,
      terminal3DG: 1.0,
      tmFwd: 1.0,
      tmRev: 1.0,
      gcFwd: 1.0,
      gcRev: 1.0,
      hairpinFwd: 1.0,
      hairpinRev: 1.0,
      selfDimerFwd: 1.0,
      selfDimerRev: 1.0,
      heterodimer: 1.0,
      gcClampFwd: 1.0,
      gcClampRev: 1.0,
      homopolymerFwd: 1.0,
      homopolymerRev: 1.0,
      tmDiff: 1.0,
      lengthFwd: 1.0,
      lengthRev: 1.0,
    };

    const result = calculateCompositeScore(scores);
    expect(result.score).toBe(100);
  });

  it("should weight off-target heavily", () => {
    const goodOffTarget = {
      offTarget: 1.0,
      tmFwd: 0.5,
    };

    const badOffTarget = {
      offTarget: 0.5,
      tmFwd: 1.0,
    };

    const goodResult = calculateCompositeScore(goodOffTarget);
    const badResult = calculateCompositeScore(badOffTarget);

    // Off-target weight (0.18) > Tm weight (0.05), so good off-target should win
    expect(goodResult.score).toBeGreaterThan(badResult.score);
  });

  it("should accept custom weights", () => {
    const scores = { tmFwd: 1.0, gcFwd: 0.5 };

    const customWeights = {
      tmFwd: 0.9,  // Heavy weight on Tm
      gcFwd: 0.1,  // Light weight on GC
    };

    const result = calculateCompositeScore(scores, customWeights);

    // With these weights, should be closer to tmFwd score
    expect(result.rawScore).toBeGreaterThan(0.9);
  });
});

describe("classifyQuality", () => {
  it("should classify 90+ as excellent", () => {
    const result = classifyQuality(95);
    expect(result.tier).toBe('excellent');
    expect(result.color).toBe('green');
  });

  it("should classify 75-89 as good", () => {
    const result = classifyQuality(80);
    expect(result.tier).toBe('good');
    expect(result.color).toBe('blue');
  });

  it("should classify 60-74 as acceptable", () => {
    const result = classifyQuality(65);
    expect(result.tier).toBe('acceptable');
    expect(result.color).toBe('yellow');
  });

  it("should classify 40-59 as marginal", () => {
    const result = classifyQuality(50);
    expect(result.tier).toBe('marginal');
    expect(result.color).toBe('orange');
  });

  it("should classify <40 as poor", () => {
    const result = classifyQuality(30);
    expect(result.tier).toBe('poor');
    expect(result.color).toBe('red');
  });
});

describe("integration: realistic primer scoring", () => {
  it("should score a typical good primer pair highly", () => {
    // Typical good primers: ~20bp, 50% GC, ~58°C Tm
    const scores = {
      tmFwd: scoreTm(58),
      tmRev: scoreTm(57),
      gcFwd: scoreGc(0.5),
      gcRev: scoreGc(0.48),
      tmDiff: scoreTmDiff(58, 57),
      offTarget: scoreOffTarget(0),
      hairpinFwd: scoreHairpin(-1),
      hairpinRev: scoreHairpin(-0.5),
      selfDimerFwd: scoreHomodimer(-4),
      selfDimerRev: scoreHomodimer(-3),
      heterodimer: scoreHeterodimer(-5),
      gcClampFwd: scoreGcClamp("ATCGATCGTG"),
      gcClampRev: scoreGcClamp("ATCGATCGAC"),
      lengthFwd: scoreLength(20),
      lengthRev: scoreLength(21),
      homopolymerFwd: scoreHomopolymer("ATCGATCGATCGATCGATCG"),
      homopolymerRev: scoreHomopolymer("ATCGATCGATCGATCGATCGA"),
    };

    const result = calculateCompositeScore(scores);

    expect(result.score).toBeGreaterThan(90);
    expect(classifyQuality(result.score).tier).toBe('excellent');
  });

  it("should score a problematic primer pair poorly", () => {
    // Problematic: high GC, mismatched Tm, off-targets, hairpin
    const scores = {
      tmFwd: scoreTm(68),  // Too high
      tmRev: scoreTm(55),  // OK
      gcFwd: scoreGc(0.75),  // Too high
      gcRev: scoreGc(0.5),
      tmDiff: scoreTmDiff(68, 55),  // Large difference
      offTarget: scoreOffTarget(2),  // Multiple off-targets
      hairpinFwd: scoreHairpin(-6),  // Stable hairpin
      hairpinRev: scoreHairpin(-1),
      selfDimerFwd: scoreHomodimer(-4),
      selfDimerRev: scoreHomodimer(-3),
      heterodimer: scoreHeterodimer(-5),
      gcClampFwd: scoreGcClamp("ATCGATCGAA"),  // Weak clamp
      gcClampRev: scoreGcClamp("ATCGATCGTG"),
      lengthFwd: scoreLength(32),  // Too long
      lengthRev: scoreLength(20),
      homopolymerFwd: scoreHomopolymer("ATCGAAAAATCGATCGATCGATCGATCGATCG"),  // 5-run
      homopolymerRev: scoreHomopolymer("ATCGATCGATCGATCGATCG"),
    };

    const result = calculateCompositeScore(scores);

    expect(result.score).toBeLessThan(70);
    expect(classifyQuality(result.score).tier).not.toBe('excellent');
  });
});
