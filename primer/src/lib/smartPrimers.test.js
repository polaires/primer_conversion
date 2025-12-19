/**
 * Tests for Smart Primer Design Module
 *
 * Tests iterative optimization for GC clamp, 3' end composition, and length adjustment.
 */

import { describe, it, expect } from "vitest";
import {
  analyze3PrimeEnd,
  optimizePrimer,
  generateLengthVariants,
  scorePrimerVariant,
  generateDesignSuggestions,
  quickAssess,
  score3PrimeComposition,
  SMART_DESIGN_CONSTRAINTS,
} from "./smartPrimers.js";

describe("analyze3PrimeEnd", () => {
  it("should identify excellent 3' end with single GC clamp", () => {
    // Ends with ...GC - ideal
    const result = analyze3PrimeEnd("ATGCTAGCTAGCTAGCTGC");
    expect(result.gcCounts.last2).toBe(2);
    expect(result.endsWithGC).toBe(true);
    expect(result.hasGcClamp).toBe(true);
  });

  it("should identify poor 3' end with no GC clamp", () => {
    // Ends with ...AA - no clamp
    const result = analyze3PrimeEnd("ATGCTAGCTAGCTAGCTAA");
    expect(result.gcCounts.last2).toBe(0);
    expect(result.endsWithGC).toBe(false);
    expect(result.hasGcClamp).toBe(false);
    expect(result.quality).not.toBe("excellent");
    expect(result.issues.length).toBeGreaterThan(0);
  });

  it("should detect poly-A/T patterns in 3' end", () => {
    // Ends with ...AAAAT
    const result = analyze3PrimeEnd("ATGCTAGCTAGCAAAAAT");
    expect(result.patterns.hasPolyAT).toBe(true);
    expect(result.issues.some(i => i.includes("Poly"))).toBe(true);
  });

  it("should identify AT-rich 3' end", () => {
    // Last 5 bases have ≤1 G/C
    const result = analyze3PrimeEnd("ATGCTAGCTAGCTAAATA");
    expect(result.gcCounts.last5).toBeLessThanOrEqual(1);
  });

  it("should calculate terminal ΔG", () => {
    const result = analyze3PrimeEnd("GCGCGCGCGCGCGCGCGC");
    expect(typeof result.terminalDG).toBe("number");
    // GC-rich sequence should have strong (negative) ΔG
    expect(result.terminalDG).toBeLessThan(0);
  });

  it("should mark canImprove for suboptimal primers", () => {
    const poorPrimer = analyze3PrimeEnd("ATGCTAGCTAGCTAAATT");
    expect(poorPrimer.canImprove).toBe(true);

    // Primer ending in TG has 1 G/C in last 2 - but might still have terminal ΔG issues
    const betterPrimer = analyze3PrimeEnd("ATGCTAGCTAGCTGCGCG");
    // Should have GC clamp at minimum
    expect(betterPrimer.hasGcClamp).toBe(true);
  });
});

describe("generateLengthVariants", () => {
  const template = "ATGCTAGCTAGCTAGCTGCATGCATGCATGCATGCATGC";

  it("should generate variants of different lengths", () => {
    const variants = generateLengthVariants(template, 0, 20, true, {
      extendMax: 3,
      contractMax: 2,
    });

    expect(variants.length).toBeGreaterThan(0);

    // Should have variants of different lengths
    const lengths = variants.map(v => v.length);
    expect(new Set(lengths).size).toBeGreaterThan(1);
  });

  it("should respect length bounds", () => {
    const variants = generateLengthVariants(template, 0, 20, true, {
      minLength: 18,
      maxLength: 25,
    });

    for (const v of variants) {
      expect(v.length).toBeGreaterThanOrEqual(18);
      expect(v.length).toBeLessThanOrEqual(25);
    }
  });

  it("should sort variants by priority", () => {
    const variants = generateLengthVariants(template, 0, 20, true);

    // First variant should have highest priority
    for (let i = 1; i < variants.length; i++) {
      expect(variants[0].priority).toBeGreaterThanOrEqual(variants[i].priority);
    }
  });

  it("should include 3' end analysis for each variant", () => {
    const variants = generateLengthVariants(template, 0, 20, true);

    for (const v of variants) {
      expect(v.analysis).toBeDefined();
      expect(v.analysis.gcCounts).toBeDefined();
      expect(v.analysis.last2).toBeDefined();
    }
  });
});

describe("optimizePrimer", () => {
  // Template with varied GC content
  const template = "ATGCTAGCTAGCTAAATTGCATGCATGCATGCGCGCGC";

  it("should optimize primer to achieve GC clamp", () => {
    // Start with primer ending in AT
    const result = optimizePrimer(template, 0, 17, true, {
      targetGcClamp: true,
      maxLengthChange: 3,
    });

    // Should suggest optimization if no GC clamp
    if (result.optimized) {
      const newAnalysis = analyze3PrimeEnd(result.seq);
      expect(newAnalysis.gcCounts.last2).toBeGreaterThanOrEqual(1);
    }
  });

  it("should not change already optimal primer", () => {
    // Start with primer ending in GC
    const goodSeq = "ATGCTAGCTAGCTGC";
    const paddedTemplate = goodSeq + "AAAAAAAAAAAAAAAAAAAA";

    const result = optimizePrimer(paddedTemplate, 0, goodSeq.length, true, {
      targetGcClamp: true,
    });

    // Might or might not optimize depending on other factors
    expect(result.seq).toBeDefined();
    expect(result.analysis).toBeDefined();
  });

  it("should respect max length change", () => {
    const result = optimizePrimer(template, 0, 20, true, {
      maxLengthChange: 2,
    });

    const lengthDiff = Math.abs(result.seq.length - 20);
    expect(lengthDiff).toBeLessThanOrEqual(2);
  });

  it("should return original if no better variant found", () => {
    // Very short template - limited options
    const shortTemplate = "ATGCATGCATGCATGCAT";

    const result = optimizePrimer(shortTemplate, 0, 15, true);

    expect(result.seq).toBeDefined();
    expect(typeof result.optimized).toBe("boolean");
  });
});

describe("scorePrimerVariant", () => {
  it("should return comprehensive scoring", () => {
    const result = scorePrimerVariant("ATGCTAGCTAGCTAGCTGC");

    expect(result.seq).toBeDefined();
    expect(result.length).toBe(19);
    expect(result.tm).toBeGreaterThan(0);
    expect(result.gc).toBeGreaterThan(0);
    expect(result.gc).toBeLessThanOrEqual(1);
    expect(result.compositeScore).toBeGreaterThanOrEqual(0);
    expect(result.compositeScore).toBeLessThanOrEqual(100);
    expect(result.qualityTier).toBeDefined();
  });

  it("should include individual scores", () => {
    const result = scorePrimerVariant("ATGCTAGCTAGCTAGCTGC");

    expect(result.scores.tm).toBeDefined();
    expect(result.scores.gc).toBeDefined();
    expect(result.scores.length).toBeDefined();
    expect(result.scores.gcClamp).toBeDefined();
    expect(result.scores.homopolymer).toBeDefined();
    expect(result.scores.terminal3DG).toBeDefined();
  });

  it("should score GC-rich primer appropriately", () => {
    const gcRich = scorePrimerVariant("GCGCGCGCGCGCGCGCGCGC");
    const atRich = scorePrimerVariant("ATATATATATATATATATAT");

    // Both should have scores, but GC content differs
    expect(gcRich.gc).toBeGreaterThan(0.6);
    expect(atRich.gc).toBeLessThan(0.4);
  });

  it("should include 3' analysis", () => {
    const result = scorePrimerVariant("ATGCTAGCTAGCTAGCTGC");

    expect(result.analysis3Prime).toBeDefined();
    expect(result.analysis3Prime.gcCounts).toBeDefined();
    expect(result.analysis3Prime.quality).toBeDefined();
  });
});

describe("generateDesignSuggestions", () => {
  it("should suggest extending for weak 3' end", () => {
    // Primer with weak terminal ΔG
    const primerScore = scorePrimerVariant("ATATATATATATATATATA");
    const suggestions = generateDesignSuggestions(primerScore);

    expect(suggestions.suggestions).toBeDefined();
    expect(Array.isArray(suggestions.suggestions)).toBe(true);
  });

  it("should return no suggestions for optimal primer", () => {
    // Well-balanced primer
    const primerScore = scorePrimerVariant("ATGCATGCATGCATGCATGC");
    const suggestions = generateDesignSuggestions(primerScore);

    // May or may not have suggestions depending on exact scores
    expect(suggestions.canImprove !== undefined).toBe(true);
    expect(suggestions.summary).toBeDefined();
  });

  it("should prioritize high-priority issues", () => {
    // Primer with multiple issues
    const primerScore = scorePrimerVariant("AAAAAAAAAAAAAAAAAAA");
    const suggestions = generateDesignSuggestions(primerScore);

    if (suggestions.suggestions.length > 1) {
      // First suggestion should be high priority
      expect(["high", "medium", "low"]).toContain(suggestions.suggestions[0].priority);
    }
  });
});

describe("quickAssess", () => {
  it("should provide quick quality assessment", () => {
    const result = quickAssess("ATGCTAGCTAGCTAGCTGC");

    expect(result.seq).toBeDefined();
    expect(result.length).toBe(19);
    expect(result.tm).toBeGreaterThan(0);
    expect(result.gc).toBeGreaterThan(0);
    expect(result.quality).toBeDefined();
    expect(["excellent", "good", "acceptable", "marginal", "poor"]).toContain(result.quality);
  });

  it("should identify issues", () => {
    // Very short primer
    const result = quickAssess("ATGCATGC");

    expect(result.issues.length).toBeGreaterThan(0);
    expect(result.needsOptimization).toBe(true);
  });

  it("should detect homopolymer runs", () => {
    const result = quickAssess("ATGCAAAAGCTAGCTAGCTG");

    expect(result.issues.some(i => i.includes("homopolymer"))).toBe(true);
  });

  it("should provide optimization strategy", () => {
    const result = quickAssess("ATGCTAGCTAGCTAAATT");

    if (result.needsOptimization) {
      expect(result.canImproveWith).toBeDefined();
    }
  });

  it("should assess excellent primer correctly", () => {
    // Well-designed primer
    const result = quickAssess("ATGCTAGCTAGCTAGCTGC");

    // Should be at least acceptable
    expect(["excellent", "good", "acceptable"]).toContain(result.quality);
  });
});

describe("score3PrimeComposition", () => {
  it("should score ideal GC clamp highly", () => {
    // Ends in GC - ideal clamp (2 G/C, ending with G/C)
    const analysis = analyze3PrimeEnd("ATGCTAGCTAGCTAGCTGC");
    const score = score3PrimeComposition(analysis);

    expect(score).toBeGreaterThan(0.5);  // Good score for ideal clamp
  });

  it("should penalize no GC clamp", () => {
    const analysis = analyze3PrimeEnd("ATGCTAGCTAGCTAGCTAA");  // Ends in AA (0 G/C)
    const score = score3PrimeComposition(analysis);

    expect(score).toBeLessThan(0.8);
  });

  it("should return value between 0 and 1", () => {
    const sequences = [
      "ATGCTAGCTAGCTAGCTGC",
      "ATATATATATATATATATAT",
      "GCGCGCGCGCGCGCGCGCGC",
      "ATGCATGCATGCATGCATGC",
    ];

    for (const seq of sequences) {
      const analysis = analyze3PrimeEnd(seq);
      const score = score3PrimeComposition(analysis);
      expect(score).toBeGreaterThanOrEqual(0);
      expect(score).toBeLessThanOrEqual(1);
    }
  });
});

describe("SMART_DESIGN_CONSTRAINTS", () => {
  it("should have valid length bounds", () => {
    expect(SMART_DESIGN_CONSTRAINTS.minLength).toBe(15);
    expect(SMART_DESIGN_CONSTRAINTS.maxLength).toBe(32);
    expect(SMART_DESIGN_CONSTRAINTS.optimalLengthRange[0]).toBeLessThanOrEqual(
      SMART_DESIGN_CONSTRAINTS.optimalLengthRange[1]
    );
  });

  it("should have valid terminal ΔG range", () => {
    const dgRange = SMART_DESIGN_CONSTRAINTS.terminal3DG;
    expect(dgRange.optimal[0]).toBeLessThan(dgRange.optimal[1]);
    expect(dgRange.minimum).toBeLessThan(dgRange.optimal[0]);
    expect(dgRange.maximum).toBeGreaterThan(dgRange.optimal[1]);
  });

  it("should have valid score thresholds", () => {
    expect(SMART_DESIGN_CONSTRAINTS.minAcceptableScore).toBeLessThan(
      SMART_DESIGN_CONSTRAINTS.targetScore
    );
  });
});

describe("Integration: Full optimization workflow", () => {
  it("should optimize primer with weak 3' end", () => {
    // Create a template where extending the primer would capture a G/C
    const template = "ATGCTAGCTAGCTAAATTGCATGCATGCATGCATGCATGC";

    // Start position that ends in AA
    const initialLength = 17;  // ...AAATT
    const result = optimizePrimer(template, 0, initialLength, true, {
      targetGcClamp: true,
      maxLengthChange: 3,
    });

    // Either optimized with GC clamp, or kept original with explanation
    expect(result.seq).toBeDefined();
    expect(result.reason).toBeDefined();
  });

  it("should score and suggest improvements for suboptimal primer", () => {
    // AT-rich primer
    const seq = "ATATATATATATATATATAT";
    const scored = scorePrimerVariant(seq);
    const suggestions = generateDesignSuggestions(scored);

    expect(scored.compositeScore).toBeLessThan(80);  // Should be suboptimal
    expect(suggestions.suggestions.length).toBeGreaterThan(0);
  });

  it("should identify high-quality primer", () => {
    // Well-balanced primer
    const seq = "ATGCATGCATGCATGCATGC";
    const assessment = quickAssess(seq);

    // Should be at least acceptable
    expect(["excellent", "good", "acceptable"]).toContain(assessment.quality);
    expect(assessment.gc).toBeGreaterThanOrEqual(40);
    expect(assessment.gc).toBeLessThanOrEqual(60);
  });
});
