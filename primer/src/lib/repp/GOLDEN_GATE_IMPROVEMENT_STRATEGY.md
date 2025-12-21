# Golden Gate Algorithm Improvement Strategy

## Executive Summary

This document outlines a comprehensive improvement strategy for the Golden Gate junction and primer selection algorithms, focusing on **quality over efficiency** and **code unification** for maintainability.

---

## 1. Unified Scoring Architecture

### Current State
- `goldengate-primer-optimizer.ts` has its own `scoringWeights` (lines 600-619)
- `scoring.ts` has calibrated `DEFAULT_WEIGHTS` (F1=88.4%, AUC=0.929)
- `weightCalibration.ts` has `ASSEMBLY_WEIGHTS` for assembly-specific scoring
- Duplication causes maintenance burden and inconsistent behavior

### Proposed Solution: Layered Scoring

```typescript
// New file: lib/repp/goldengate-weights.ts

import { ASSEMBLY_WEIGHTS } from '../weightCalibration.js';

/**
 * Golden Gate-specific weight modifiers
 * Applied as ADDITIONS to ASSEMBLY_WEIGHTS, not replacements
 *
 * Philosophy: Base scoring captures primer quality,
 * GG modifiers capture assembly-specific concerns
 */
export const GG_WEIGHT_MODIFIERS = {
  // GG-specific additions (not in base scoring)
  overhangFidelity: 0.15,       // Ligation matrix fidelity
  gtMismatchRisk: 0.08,         // G:T wobble cross-ligation
  flankingQuality: 0.05,        // NEB flanking sequence
  recognitionSiteCreation: 0.10, // Avoid creating new sites

  // Overrides for assembly context (optional)
  overrides: {
    // In GG, heterodimer between primers is less critical
    // because fragments are mixed after PCR
    heterodimer: 0.04,  // Lower from 0.10
  },
};

/**
 * Compute unified GG primer score
 * Uses base ASSEMBLY_WEIGHTS + GG modifiers
 */
export function computeGGPrimerScore(
  primer: PrimerAnalysis,
  ggContext: GGContext
): UnifiedScore {
  // Step 1: Get base score from unified scoring
  const baseScore = calculateCompositeScore(primer.scores, ASSEMBLY_WEIGHTS);

  // Step 2: Apply GG-specific modifiers (additive layer)
  const ggModifiers = computeGGModifiers(primer, ggContext);

  // Step 3: Combine with configurable blend
  // Default: 70% base, 30% GG-specific
  const blendRatio = ggContext.options?.ggBlendRatio ?? 0.30;

  return {
    composite: baseScore.score * (1 - blendRatio) + ggModifiers.score * blendRatio,
    base: baseScore,
    ggSpecific: ggModifiers,
    quality: classifyQuality(/* combined */),
  };
}
```

### Benefits
- Single source of truth for base weights
- GG-specific concerns are clearly separated
- Easy to tune blend ratio for different use cases
- No duplication of scoring logic

---

## 2. PCR Annealing Temperature Strategy

### Current State
- No explicit temperature optimization
- Primers designed without considering thermal cycling constraints

### Proposed Solution: Temperature-Aware Primer Selection

```typescript
// In goldengate-primer-optimizer.ts

export interface TemperaturePreferences {
  // Primary target (NEB recommendation: 55-60°C for GG)
  targetTa: number;              // Default: 58

  // User preference: prioritize lower Ta when possible
  preferLowerTa: boolean;        // Default: true

  // Hard ceiling (for heat-sensitive applications)
  maxAcceptableTa: number;       // Default: 72

  // Acceptable range
  acceptableRange: {
    min: number;                 // Default: 50
    max: number;                 // Default: 72
  };
}

/**
 * Temperature-aware homology selection
 *
 * Strategy:
 * 1. Find all valid homology regions (18-25bp, acceptable Tm)
 * 2. If preferLowerTa: sort by Tm ascending
 * 3. Among candidates with Tm < maxAcceptableTa, choose best quality
 * 4. If no candidate < maxAcceptableTa, return best available with warning
 */
export function selectOptimalHomology(
  template: string,
  position: number,
  direction: 'forward' | 'reverse',
  options: HomologyOptions
): HomologySelection {
  const tempPrefs = options.temperature ?? getDefaultTemperaturePrefs();

  // Generate candidates with different lengths/shifts
  const candidates = generateHomologyCandidates(template, position, direction, {
    minLength: 18,
    maxLength: 28,
    shiftRange: [-5, 5],
  });

  // Score each candidate
  const scored = candidates.map(c => ({
    ...c,
    tm: calculateTmQ5(c.sequence),
    quality: scoreHomologyQuality(c),
  }));

  // Partition by temperature preference
  const belowThreshold = scored.filter(c => c.tm <= tempPrefs.maxAcceptableTa);
  const aboveThreshold = scored.filter(c => c.tm > tempPrefs.maxAcceptableTa);

  if (belowThreshold.length > 0) {
    // Sort by preference
    if (tempPrefs.preferLowerTa) {
      // Among those below threshold, prefer lower Tm with acceptable quality
      belowThreshold.sort((a, b) => {
        // Primary: quality score (must be acceptable)
        if (a.quality.score >= 70 && b.quality.score < 70) return -1;
        if (b.quality.score >= 70 && a.quality.score < 70) return 1;

        // Secondary: prefer lower Tm
        return a.tm - b.tm;
      });
    } else {
      // Sort by quality only
      belowThreshold.sort((a, b) => b.quality.score - a.quality.score);
    }

    return {
      selected: belowThreshold[0],
      alternatives: belowThreshold.slice(1, 4),
      warnings: [],
      temperatureOptimized: true,
    };
  }

  // Fallback: use best available, add warning
  aboveThreshold.sort((a, b) => b.quality.score - a.quality.score);

  return {
    selected: aboveThreshold[0],
    alternatives: aboveThreshold.slice(1, 4),
    warnings: [
      `No primer found with Tm ≤ ${tempPrefs.maxAcceptableTa}°C. ` +
      `Selected primer has Tm = ${aboveThreshold[0].tm.toFixed(1)}°C. ` +
      `Consider adjusting junction position or using shorter homology.`
    ],
    temperatureOptimized: false,
  };
}
```

### Advanced Options UI Integration

```typescript
// In GoldenGateDesigner.tsx - Advanced Options section

interface AdvancedPrimerOptions {
  temperature: {
    preferLowerTa: boolean;           // Checkbox: "Prefer lower annealing temp"
    maxAcceptableTa: number;          // Slider: 55-72°C, default 72
    targetTa: number;                 // Input: target, default 58
  };

  homology: {
    minLength: number;                // 15-20, default 18
    maxLength: number;                // 25-35, default 28
    allowShifting: boolean;           // Allow junction position adjustment
  };

  scoring: {
    ggBlendRatio: number;             // 0-1, default 0.30
    strictQuality: boolean;           // Reject marginal primers
  };
}
```

---

## 3. Overhang Optimization Improvements

### Issue 1: Monte Carlo Temperature Finding

**Current**: Fixed 1000 iterations, fixed target AR (5%)

**Proposed**: Adaptive with quality focus

```typescript
/**
 * Improved temperature calibration
 *
 * Changes:
 * 1. Use more iterations for better calibration (2000)
 * 2. Adaptive AR based on problem size
 * 3. Track quality of solutions found during calibration
 */
function findAnnealingTemperature(
  pool: string[],
  size: number,
  matrix: LigationMatrix,
  options: TemperatureOptions = {}
): TemperatureResult {
  const {
    calibrationIterations = 2000,
    // Adaptive AR: larger sets need lower AR to converge
    targetAR = Math.max(0.02, 0.10 - size * 0.005),
  } = options;

  // ... existing binary search logic ...

  return {
    temperature: exp,
    actualAR: ar,
    targetAR,
    calibrationQuality: bestFidelityDuringCalibration,
  };
}
```

### Issue 2: Greedy Fallback Without Improvement

**Current**: Greedy selects once, no refinement

**Proposed**: Add 2-opt local search

```typescript
/**
 * Enhanced greedy with 2-opt improvement
 *
 * After initial greedy selection, try swapping pairs
 * to find local improvements
 */
function findOptimalOverhangSetGreedy(
  numJunctions: number,
  enzyme: string,
  options: OverhangSearchOptions = {}
): OptimalOverhangSearchResult {
  // Step 1: Initial greedy selection (existing code)
  let resultSet = greedyInitialSelection(numJunctions, enzyme, options);
  let bestFidelity = calculateExperimentalFidelity(resultSet, enzyme).assemblyFidelity;

  // Step 2: 2-opt improvement phase
  const maxImprovementIterations = 50;
  let improved = true;
  let iterations = 0;

  while (improved && iterations < maxImprovementIterations) {
    improved = false;
    iterations++;

    // Try swapping each overhang with candidates from pool
    for (let i = 0; i < resultSet.length; i++) {
      if (fixedIndices.has(i)) continue; // Don't swap required overhangs

      const original = resultSet[i];

      for (const candidate of candidates) {
        if (resultSet.includes(candidate)) continue;
        if (hasRCConflict(resultSet, candidate, i)) continue;

        // Try swap
        const testSet = [...resultSet];
        testSet[i] = candidate;

        const testFidelity = calculateExperimentalFidelity(testSet, enzyme).assemblyFidelity;

        if (testFidelity > bestFidelity) {
          resultSet = testSet;
          bestFidelity = testFidelity;
          improved = true;
          break; // Restart from beginning of set
        }
      }

      if (improved) break;
    }
  }

  // Step 3: Random restarts for escaping local optima
  if (options.enableRandomRestarts !== false) {
    const numRestarts = 3;
    for (let r = 0; r < numRestarts; r++) {
      const randomStart = generateRandomValidSet(numJunctions, candidates, requiredOverhangs);
      const refined = localOptimize(randomStart);

      if (refined.fidelity > bestFidelity) {
        resultSet = refined.set;
        bestFidelity = refined.fidelity;
      }
    }
  }

  return {
    overhangs: resultSet,
    fidelity: bestFidelity,
    // ...
  };
}
```

### Issue 3: Site Recreation Check

**Current**: Simplified check (lines 774-796)

**Proposed**: Full context check

```typescript
/**
 * Comprehensive site recreation check
 *
 * Checks if assembly scar would recreate the recognition site
 * by examining the full junction context
 */
function wouldRecreateRecognitionSite(
  upstreamSeq: string,
  downstreamSeq: string,
  overhang: string,
  enzyme: GoldenGateEnzyme
): SiteRecreationResult {
  const recognition = enzyme.recognition;
  const recognitionRC = reverseComplement(recognition);
  const contextSize = recognition.length + 4; // Extra buffer

  // Get sequences around the junction
  const upstream = upstreamSeq.slice(-contextSize).toUpperCase();
  const downstream = downstreamSeq.slice(0, contextSize).toUpperCase();

  // The scar region: upstream-end + overhang + downstream-start
  // But in GG assembly, overhang is shared (comes from upstream OR downstream)
  // So the actual scar is: upstream[:-overhangLen] + overhang + downstream

  const overhangLen = overhang.length;
  const scarContext = upstream.slice(0, -overhangLen) + overhang + downstream;

  // Check for recognition site in scar region
  const fwdMatch = scarContext.indexOf(recognition);
  const revMatch = scarContext.indexOf(recognitionRC);

  if (fwdMatch >= 0 || revMatch >= 0) {
    return {
      recreatesSite: true,
      direction: fwdMatch >= 0 ? 'forward' : 'reverse',
      position: fwdMatch >= 0 ? fwdMatch : revMatch,
      context: scarContext,
      recommendation: 'Adjust junction position to avoid site recreation',
    };
  }

  return {
    recreatesSite: false,
    context: scarContext,
  };
}
```

---

## 4. Code Unification Plan

### 4.1 Unified Fidelity Calculation

Create single source for fidelity calculation:

```typescript
// New file: lib/repp/fidelity-core.ts

/**
 * Unified fidelity calculation module
 *
 * Consolidates:
 * - overhang-optimizer.ts::calculateLigationFidelity
 * - goldengate.ts::calculateExperimentalFidelity
 * - auto-domestication-optimizer.ts::validateDomesticationOverhangs
 */

export interface FidelityResult {
  assemblyFidelity: number;
  junctionFidelities: JunctionFidelity[];
  source: 'matrix' | 'static';
  warnings: string[];
  gtRisks?: GTMismatchRisk[];
}

/**
 * Calculate assembly fidelity with unified interface
 */
export function calculateFidelity(
  overhangs: string[],
  enzyme: string,
  options: FidelityOptions = {}
): FidelityResult {
  const {
    includeGTRisks = true,
    useMatrixData = true,
    fallbackFidelity = 0.85,
  } = options;

  // Single implementation used everywhere
  const enzymeData = getEnzymeLigationData(enzyme);

  if (!enzymeData || !useMatrixData) {
    return calculateStaticFidelity(overhangs, fallbackFidelity);
  }

  return calculateMatrixFidelity(overhangs, enzymeData, { includeGTRisks });
}
```

### 4.2 Unified Overhang Validation

```typescript
// New file: lib/repp/overhang-validation.ts

/**
 * Consolidated overhang validation
 *
 * Replaces scattered validation in:
 * - overhang-optimizer.ts::filterOverhangs
 * - auto-domestication-optimizer.ts::validateDomesticationOverhangs
 * - goldengate.ts::findOptimalOverhangSetExhaustive (inline checks)
 */

export interface OverhangValidation {
  isValid: boolean;
  issues: ValidationIssue[];
  score: number;  // 0-100 quality score
}

export function validateOverhang(
  overhang: string,
  context: ValidationContext
): OverhangValidation {
  const issues: ValidationIssue[] = [];
  let score = 100;

  // Check 1: Palindrome
  if (isPalindrome(overhang)) {
    issues.push({ type: 'palindrome', severity: 'error', message: 'Self-complementary' });
    score -= 100; // Invalid
  }

  // Check 2: Homopolymer
  if (isHomopolymer(overhang)) {
    issues.push({ type: 'homopolymer', severity: 'warning', message: 'All same base' });
    score -= 30;
  }

  // Check 3: GC content
  const gc = countGC(overhang) / overhang.length;
  if (gc === 0 || gc === 1) {
    issues.push({ type: 'extreme_gc', severity: 'warning', message: `${gc * 100}% GC` });
    score -= 15;
  }

  // Check 4: RC conflict with existing set
  if (context.existingOverhangs) {
    const rc = reverseComplement(overhang);
    for (const existing of context.existingOverhangs) {
      if (existing === overhang || existing === rc) {
        issues.push({ type: 'duplicate', severity: 'error', message: 'Duplicate/RC in set' });
        score -= 100;
      }
    }
  }

  // Check 5: Cross-ligation potential
  if (context.checkCrossLigation && context.matrix) {
    // Use matrix data to check cross-ligation
    const crossLigationScore = checkCrossLigation(overhang, context);
    if (crossLigationScore > 0.05) {
      issues.push({ type: 'cross_ligation', severity: 'warning',
                    message: `${(crossLigationScore * 100).toFixed(1)}% cross-ligation risk` });
      score -= crossLigationScore * 50;
    }
  }

  return {
    isValid: !issues.some(i => i.severity === 'error'),
    issues,
    score: Math.max(0, score),
  };
}

export function validateOverhangSet(
  overhangs: string[],
  enzyme: string
): SetValidation {
  // Unified set validation
  const individual = overhangs.map(oh => validateOverhang(oh, {
    existingOverhangs: overhangs.filter(o => o !== oh),
    checkCrossLigation: true,
    matrix: getEnzymeLigationData(enzyme)?.matrix,
  }));

  const fidelity = calculateFidelity(overhangs, enzyme);

  return {
    overhangs,
    individual,
    setFidelity: fidelity.assemblyFidelity,
    isValid: individual.every(v => v.isValid),
    recommendation: generateRecommendation(individual, fidelity),
  };
}
```

### 4.3 Unified Primer Design Pipeline

```typescript
// Refactored goldengate-primer-optimizer.ts

import { calculateCompositeScore, classifyQuality } from '../scoring.js';
import { ASSEMBLY_WEIGHTS } from '../weightCalibration.js';
import { calculateFidelity } from './fidelity-core.js';
import { validateOverhangSet } from './overhang-validation.js';
import { GG_WEIGHT_MODIFIERS, computeGGModifiers } from './goldengate-weights.js';

/**
 * Unified Golden Gate primer design pipeline
 */
export function designGoldenGatePrimers(
  parts: Part[],
  options: GGDesignOptions
): GGDesignResult {
  // Step 1: Validate/optimize overhangs
  const overhangResult = optimizeOverhangs(parts, options);

  // Step 2: Design primers for each part using unified scoring
  const primers = parts.map((part, i) => {
    const leftOH = overhangResult.overhangs[i];
    const rightOH = overhangResult.overhangs[i + 1];

    return designPartPrimers(part, leftOH, rightOH, {
      ...options,
      // Use unified scoring with GG modifiers
      scoringWeights: ASSEMBLY_WEIGHTS,
      ggModifiers: GG_WEIGHT_MODIFIERS,
      temperature: options.temperature,
    });
  });

  // Step 3: Calculate unified quality metrics
  const quality = calculateUnifiedQuality(primers, overhangResult);

  return {
    enzyme: options.enzyme,
    parts: primers,
    overhangs: overhangResult,
    quality,
    recommendations: generateRecommendations(primers, overhangResult, quality),
  };
}
```

---

## 5. Implementation Priority

### Phase 1: Foundation (High Priority)
1. Create `fidelity-core.ts` - unified fidelity calculation
2. Create `overhang-validation.ts` - consolidated validation
3. Create `goldengate-weights.ts` - GG-specific scoring modifiers

### Phase 2: Algorithm Improvements (Medium Priority)
1. Add 2-opt local search to greedy optimization
2. Implement comprehensive site recreation check
3. Add random restarts for escaping local optima

### Phase 3: Temperature Optimization (Medium Priority)
1. Implement temperature-aware homology selection
2. Add advanced options UI for temperature preferences
3. Add warnings for high-temperature primers

### Phase 4: Scoring Unification (High Priority)
1. Refactor `goldengate-primer-optimizer.ts` to use unified scoring
2. Remove duplicate weight definitions
3. Implement layered scoring (base + GG modifiers)

### Phase 5: Testing & Validation
1. Create comprehensive test suite for fidelity calculations
2. Validate against known good/bad assemblies
3. Performance benchmarks (quality metrics, not speed)

---

## 6. File Changes Summary

| File | Action | Description |
|------|--------|-------------|
| `lib/repp/fidelity-core.ts` | **NEW** | Unified fidelity calculation |
| `lib/repp/overhang-validation.ts` | **NEW** | Consolidated validation |
| `lib/repp/goldengate-weights.ts` | **NEW** | GG-specific weight modifiers |
| `lib/repp/overhang-optimizer.ts` | **MODIFY** | Add 2-opt, random restarts |
| `lib/repp/goldengate-primer-optimizer.ts` | **REFACTOR** | Use unified scoring |
| `lib/repp/goldengate.ts` | **REFACTOR** | Use fidelity-core |
| `lib/repp/auto-domestication-optimizer.ts` | **REFACTOR** | Use overhang-validation |

---

## 7. API Changes (Backward Compatible)

All existing exports will remain, with deprecation warnings for functions that are replaced:

```typescript
// goldengate.ts
/** @deprecated Use calculateFidelity from fidelity-core.ts */
export function calculateExperimentalFidelity(...) {
  console.warn('calculateExperimentalFidelity is deprecated, use calculateFidelity');
  return calculateFidelity(...);
}
```

---

## 8. Success Metrics

| Metric | Current | Target |
|--------|---------|--------|
| Code duplication | ~15% | <5% |
| Fidelity calculation consistency | Multiple implementations | Single source |
| Quality score correlation with success | Unknown | Validate with dataset |
| Temperature optimization | None | Available with <72°C preference |
| Local optima escape | None | Random restarts + 2-opt |

---

## Appendix: Your Feedback Integration

### A. Annealing Temperature (Your Point)
✅ **Agreed**: Implemented as temperature-aware homology selection with:
- `preferLowerTa: boolean` - prioritize options below threshold
- `maxAcceptableTa: number` - default 72°C
- Advanced options UI for user control

### B. Unified Scoring (Your Point)
✅ **Agreed**: Implemented as layered architecture:
- Base: `ASSEMBLY_WEIGHTS` from `weightCalibration.ts`
- Layer: `GG_WEIGHT_MODIFIERS` for assembly-specific concerns
- Configurable blend ratio (default 70% base, 30% GG)

### C. Code Unification (Your Point)
✅ **Agreed**: Three new modules to consolidate:
- `fidelity-core.ts` - single fidelity calculation
- `overhang-validation.ts` - unified validation
- `goldengate-weights.ts` - GG-specific weights only

### D. Quality Over Efficiency (Your Point)
✅ **Agreed**: All improvements prioritize quality:
- 2-opt local search (better solutions, more compute)
- Random restarts (escape local optima)
- Temperature optimization (better primers, not faster design)
- No mention of performance optimization
