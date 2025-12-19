# Empirically-Calibrated Primer Scoring for Plasmid PCR and Sanger Sequencing

## Manuscript Strategy Document

**Working Title**: "Empirically-Calibrated Primer Scoring: A Piecewise Logistic Model for Plasmid PCR and Sanger Sequencing"

**Target Journal**: Nucleic Acids Research (Methods), Bioinformatics, or PLOS Computational Biology

**Last Updated**: 2025-12-07

---

## 0. Quick Status Summary

### Current State (Dec 2025)
```
┌─────────────────────────────────────────────────────────────────┐
│                    IMPLEMENTATION STATUS                         │
├─────────────────────────────────────────────────────────────────┤
│  ✅ DONE: Equilibrium efficiency calculation (Pythia model)     │
│  ✅ DONE: DNA24 thermodynamic parameters integrated             │
│  ✅ DONE: All chemical species modeled (hairpin, dimer, etc.)   │
│  ✅ DONE: Integration with primers.js scoring system            │
│  ✅ DONE: Piecewise logistic scoring functions (scoring.js)     │
│  ✅ DONE: Composite scoring with empirical weights              │
│  ✅ DONE: Quality tier classification                           │
│  ✅ DONE: Off-target Type A-F classification (offTargetClass.)  │
│  ✅ DONE: Sanger-specific features (ROI distance, structure)    │
│  ✅ DONE: Weight calibration framework (weightCalibration.js)   │
│  ✅ DONE: Calibration on Döring dataset (829 pairs, F1=81.9%)   │
│  ✅ DONE: Cross-validation on Kayama dataset (2232 pairs, AUC=0.64)│
│  ✅ DONE: Full integration with primers.js design workflow      │
│  ✅ DONE: Off-target Type A-F classification integrated         │
│  ✅ DONE: 242 tests passing (all modules)                       │
│  ✅ DONE: Preliminary figures generated (ROC, feature importance)│
│                                                                  │
│  🔲 TODO: Wet-lab validation with partner sequencing core       │
│  🔄 IN PROGRESS: Manuscript preparation (figures drafted)       │
└─────────────────────────────────────────────────────────────────┘
```

### Immediate Next Steps (Priority Order)

| Priority | Task | Effort | Impact |
|----------|------|--------|--------|
| ✅ | **Piecewise logistic functions** - COMPLETED | - | - |
| ✅ | **Off-target Type A-F** - COMPLETED (offTargetClassification.js) | - | - |
| ✅ | **Sanger-specific features** - COMPLETED (ROI distance, structure) | - | - |
| ✅ | **Weight calibration framework** - COMPLETED (weightCalibration.js) | - | - |
| ✅ | **Calibration with real data** - COMPLETED (Döring dataset, 829 pairs) | - | - |
| ✅ | **Cross-domain validation** - COMPLETED (Kayama dataset, 2232 pairs, AUC=0.64) | - | - |
| ✅ | **Preliminary figures** - COMPLETED (ROC, feature importance) | - | - |
| ✅ | **Design workflow integration** - COMPLETED (primers.js, off-target classification) | - | - |
| 1️⃣ | **Wet-lab validation** - Partner with sequencing core | High | Critical |
| 2️⃣ | **Manuscript preparation** - Methods, results sections | Medium | High |

### Files to Modify for Next Phase

```
src/lib/
├── equilibrium.js            ✅ Completed - equilibrium efficiency model
├── scoring.js                ✅ Completed - piecewise logistic functions + calibrated weights
├── offTargetClassification.js ✅ Completed - Type A-F classification
├── weightCalibration.js      ✅ Completed - weight optimization framework
├── primers.js                ✅ Completed - full integration with scoring system
│   ├── primers() now outputs composite scores and quality tiers
│   ├── Off-target Type A-F classification integrated
│   ├── useCompositeScore option for calibrated ranking
│   └── 242 tests passing
└── offTargets.js             ← Simple count-based (kept for backward compatibility)
```

---

## 1. Problem Statement

### 1.1 Current State of Primer Scoring

Existing primer design tools (Primer3, PrimerBank, etc.) use **heuristic penalty weights** that are:
- Based on "accumulated knowledge" rather than experimental data
- Not systematically validated against wet-lab success/failure
- Designed for genome-scale applications, not plasmid-specific contexts

### 1.2 Gap in Literature

| Existing Tool | Limitation for Plasmid/Sanger |
|---------------|-------------------------------|
| Primer3 | Weights not empirically calibrated; genome-focused |
| PrimerBank | RT-qPCR only; fixed 60°C annealing |
| PrimerScore2 | Multiplex NGS; depth-ratio optimization not binary success |
| GM1 Model | Genome-scale off-target; overkill for plasmid |
| Pythia | Genome-scale; computationally expensive; SVM pre-filter needed |

### 1.3 Our Contribution

A primer scoring system that:
1. Uses **empirically-derived weights** from meta-analysis of published validation studies
2. Is **optimized for plasmid PCR and Sanger sequencing** context
3. Employs **piecewise logistic functions** for biologically meaningful scoring
4. Validates against **binary success/failure metrics** (not depth ratios)

---

## 2. Critical Clarifications

### 2.1 Issue: R² = 0.935 Context Mismatch

**Problem Identified**:
The PrimerScore2 R² = 0.935 correlation was for predicting **relative amplification efficiency ratios in multiplex competition**:

> "The depth ratios of the products were linearly correlated with the predicted efficiencies"

This metric is specific to:
- **Multiplex context**: Primers compete for resources
- **NGS depth ratios**: Relative abundance between amplicons
- **Continuous outcome**: Efficiency as a ratio

**For Monoplex Sanger, this is NOT applicable because**:
- No primer competition (single amplicon)
- Binary outcome: band/no band, readable sequence/failed
- Absolute success, not relative efficiency

**Resolution**:
Use the **mathematical framework** (piecewise logistic) from PrimerScore2, but validate using **binary success metrics**:

| Metric | Definition | Target |
|--------|------------|--------|
| **Accuracy** | (TP + TN) / Total | > 85% |
| **Precision** | TP / (TP + FP) | > 80% |
| **Recall** | TP / (TP + FN) | > 90% |
| **AUC-ROC** | Area under ROC curve | > 0.85 |
| **F1 Score** | Harmonic mean of P & R | > 0.85 |

**Reference Benchmarks**:
- ML-PCR Study (Preprint): 81% accuracy (binary)
- PrimerBank: 82.6% success rate (baseline)
- GM1 Model: Reduced failure 17% → 6% (binary)
- **Pythia: 81% precision, 97% recall** (binary, thermodynamic model)

---

### 2.3 Pythia Equilibrium Efficiency Model - Detailed Evaluation

#### 2.3.1 Overview

Pythia (Mann et al., 2009, Nucleic Acids Research) uses a **thermodynamic equilibrium approach** rather than heuristic penalty weights. It is the most theoretically grounded of all validation studies.

**Core Principle**:
> "If a primer pair works under the equilibrium model, then it will work in PCR conditions."

#### 2.3.2 Equilibrium Efficiency Calculation

```
Equilibrium Efficiency = min(η_left, η_right)

Where:
  η_left  = [Primer_L bound to target] / [Total Primer_L]
  η_right = [Primer_R bound to target] / [Total Primer_R]

The minimum is taken because PCR is only as efficient as its
least efficient priming reaction.
```

**Chemical Species Considered**:
| Species | Description |
|---------|-------------|
| Primer (unfolded) | Free primer in solution |
| Primer (folded) | Hairpin structure |
| Primer-Primer (homodimer) | Self-dimer |
| Primer-Template (target) | Desired binding |
| Primer-Template (off-target) | Mispriming |
| Primer_L-Primer_R (heterodimer) | Cross-dimer |

**Gibbs Free Energy Calculation**:
```
ΔG = -RT ln(Ka)

Where:
  R = gas constant (1.987 cal/mol·K)
  T = temperature (K)
  Ka = association constant

Equilibrium concentrations solved by minimizing total Gibbs energy G
using gradient descent optimization.
```

#### 2.3.3 Algorithm Steps

```
Step 1: CANDIDATE GENERATION
  - Enumerate all primer pairs satisfying constraints:
    * Primer Tm within range
    * Primer length within range
    * Amplicon length within range
  - Sort by |Tm_actual - Tm_desired|

Step 2: SVM PRE-FILTER (Speed Optimization)
  - SVM classifier predicts feasibility
  - If feasible → proceed to full equilibrium analysis
  - Reduces candidates from ~10,000 to ~900 in 10 min

Step 3: EQUILIBRIUM ANALYSIS
  - Compute ΔG for all duplex and folded forms
  - Solve for equilibrium concentrations
  - Calculate equilibrium efficiency metric
  - If metric > threshold → proceed

Step 4: SPECIFICITY CHECK
  - Use precomputed genomic index
  - Check for off-target binding sites
  - If specific → output primer pair
```

#### 2.3.4 Validation Results

| Metric | Pythia | Primer3 (P316) | Improvement |
|--------|--------|----------------|-------------|
| **Precision** | 81% | 81% | Equal |
| **Recall** | 97% | 48% | **+49%** |
| **Coverage (RepeatMasked)** | 89% | 51% | **+38%** |
| **Coverage (Interspersed repeats)** | 80% | 50% | **+30%** |

**Key Finding**:
> "About 95% of Pythia primers were scored as unacceptable by Primer3's scoring function, with only three of the unacceptable primer pairs resulting in failed PCR."

This demonstrates that **Primer3's heuristic rules are overly conservative** and reject many viable primers.

#### 2.3.5 Strengths for Our Use Case

| Strength | Relevance to Plasmid/Sanger |
|----------|----------------------------|
| **Thermodynamic basis** | Physically meaningful, not heuristic |
| **All species modeled** | Hairpin, dimer, heterodimer, off-target |
| **Binary validation** | Matches our success/failure metric |
| **High recall (97%)** | Won't miss good primers |
| **Handles difficult regions** | GC-rich, repeats common in plasmids |

#### 2.3.6 Limitations / Adaptations Needed

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| **Genome-scale indexing** | Overkill for plasmid | Use input-scoped search |
| **SVM pre-filter required** | Adds complexity | May not need for plasmid (smaller search space) |
| **~20 sec per design** | Too slow for real-time | Acceptable for batch; can optimize |
| **C++ implementation** | Harder to integrate | Port equilibrium logic to JS |
| **2009 parameters** | May be outdated | Update to DNA24 (Greenleaf 2024) |

#### 2.3.7 What to Adopt from Pythia

**ADOPT** (High Value):
1. **Equilibrium efficiency as core metric** - replaces arbitrary penalty sum
2. **All chemical species in model** - comprehensive thermodynamics
3. **Min(η_left, η_right)** - correct bottleneck identification
4. **Binary success validation** - matches our use case

**ADAPT** (Modify for Context):
1. **Off-target search** - use input-scoped, not genome index
2. **SVM pre-filter** - may not need for plasmid (smaller space)
3. **Thermodynamic parameters** - update to DNA24

**SKIP** (Not Needed):
1. Genome indexing infrastructure
2. Repeat-masked sequence handling (plasmids are designed)

---

### 2.4 Issue: Off-Target Definition for Plasmid Context

**Problem Identified**:
"Off-target (input-scoped)" is ambiguous. For plasmid work, we need precise definitions.

#### 2.2.1 Off-Target Scenario Classification

| Scenario | Description | Risk Level | Detection Method |
|----------|-------------|------------|------------------|
| **Type A: Alternative binding site** | Primer binds elsewhere on same plasmid (sense strand) | 🔴 HIGH - wrong product | Exact/near-exact match search |
| **Type B: Antisense binding** | Primer binds to reverse complement of template | 🔴 HIGH - primer-template hybrid | RevComp alignment |
| **Type C: Partial 3' homology** | 3' end (≥8bp) matches elsewhere | 🟡 MEDIUM - mispriming | 3'-anchored alignment |
| **Type D: Internal homology** | Homology not at 3' end | 🟢 LOW - usually tolerable | Full alignment |
| **Type E: Self-complementarity** | Primer binds to itself | 🟡 MEDIUM - primer dimer | Fold/dimer analysis |
| **Type F: Primer-primer** | Forward binds to reverse | 🟡 MEDIUM - primer dimer | Heterodimer analysis |

#### 2.2.2 Proposed Off-Target Algorithm

```
For each primer P against template T:

1. EXACT MATCH SEARCH (Type A)
   - Find all exact matches of P in T
   - If count > 1: HIGH penalty (alternative amplification)

2. NEAR-EXACT SEARCH (Type A variant)
   - Find matches with ≤2 mismatches (not at 3' end)
   - Weight by mismatch count and position

3. REVERSE COMPLEMENT CHECK (Type B)
   - Align P against revcomp(T)
   - If 3' end matches ≥10bp: HIGH penalty

4. 3' ANCHORED SEARCH (Type C)
   - Extract last 8-12bp of primer (3' end)
   - Search for exact matches elsewhere in T
   - If found: MEDIUM penalty (position-dependent)

5. STABILITY FILTER
   - For any off-target site, calculate binding ΔG
   - Only penalize if ΔG < -11 kcal/mol (stable enough to prime)
```

#### 2.2.3 Off-Target Scoring Formula

```javascript
/**
 * Off-target penalty for plasmid context
 *
 * Key insight from GM1: off-target is "dominant failure factor"
 * But for plasmid, we only search within input sequence
 */
function calculateOffTargetPenalty(primer, template, options = {}) {
  const {
    // Thresholds
    exactMatchPenalty = 50,        // Per additional exact match
    nearMatchPenalty = 20,         // Per near-match (≤2 mismatches)
    partialHomologyPenalty = 10,   // Per 3'-anchored partial match

    // 3' end parameters
    anchorLength = 8,              // Minimum 3' anchor for concern
    stabilityThreshold = -11,      // ΔG threshold for viable priming

    // Mismatch tolerance
    maxMismatches = 2,
    terminal3Mismatches = 0,       // No mismatches in last 3bp
  } = options;

  let penalty = 0;
  const sites = [];

  // Type A: Exact matches (excluding intended site)
  const exactMatches = findExactMatches(primer, template);
  if (exactMatches.length > 1) {
    penalty += (exactMatches.length - 1) * exactMatchPenalty;
    sites.push(...exactMatches.slice(1).map(m => ({ type: 'exact', ...m })));
  }

  // Type B: Reverse complement binding
  const rcMatches = findReverseComplementMatches(primer, template, anchorLength);
  for (const match of rcMatches) {
    if (match.bindingDG < stabilityThreshold) {
      penalty += exactMatchPenalty;  // As severe as exact match
      sites.push({ type: 'antisense', ...match });
    }
  }

  // Type C: Partial 3' homology
  const anchor = primer.slice(-anchorLength);
  const partialMatches = findPartialMatches(anchor, template, primer);
  for (const match of partialMatches) {
    if (match.is3primeAnchored && match.bindingDG < stabilityThreshold) {
      penalty += partialHomologyPenalty;
      sites.push({ type: 'partial_3prime', ...match });
    }
  }

  return {
    penalty,
    sites,
    count: sites.length,
    hasCritical: sites.some(s => s.type === 'exact' || s.type === 'antisense'),
  };
}
```

#### 2.2.4 Off-Target Scoring Weight

Based on GM1 finding that off-target is the "dominant failure factor":

| Off-Target Count | Current System | Proposed System |
|------------------|----------------|-----------------|
| 0 sites | 0 penalty | 0 penalty |
| 1 site | +20 (linear) | +30 (exponential base) |
| 2 sites | +40 (linear) | +90 (3× previous) |
| 3+ sites | +60+ (linear) | **Disqualify** (flag as critical) |

---

## 3. Validation Framework

### 3.1 Metric Selection Rationale

**Why Binary Success Metrics Over Depth Ratios**:

| Metric Type | Applicable Context | Our Context |
|-------------|-------------------|-------------|
| Depth ratio (R²) | Multiplex NGS | ❌ Not applicable |
| Relative efficiency | Competition assays | ❌ Not applicable |
| **Binary success** | Monoplex PCR, Sanger | ✅ Applicable |
| **Read quality score** | Sanger sequencing | ✅ Applicable |

### 3.2 Success Definition for Plasmid PCR/Sanger

```
SUCCESS = (PCR_SUCCESS) AND (SANGER_SUCCESS)

Where:
  PCR_SUCCESS:
    - Single band on gel at expected size (±10%)
    - Band intensity above threshold
    - No significant off-target bands

  SANGER_SUCCESS:
    - Phred Q20 read length > 400bp
    - Clean chromatogram (no mixed peaks)
    - Sequence matches expected reference
```

### 3.3 Validation Dataset Requirements

| Requirement | Minimum | Ideal |
|-------------|---------|-------|
| Total primer pairs | 200 | 500+ |
| Unique templates | 20 | 50+ |
| Template size range | 2-15 kb | 1-20 kb |
| GC content range | 30-70% | 20-80% |
| Success/failure ratio | ~80/20 | Natural distribution |
| Technical replicates | 2 | 3 |

### 3.4 Comparison Baselines

| Tool | How to Compare |
|------|----------------|
| Primer3 | Design primers with Primer3, test same templates |
| Current system | A/B test against current scoring |
| Random baseline | Random primers within constraints |
| "Perfect" primers | Manually optimized by expert |

---

## 4. Mathematical Framework

### 4.1 Why Piecewise Logistic (Not Linear)

**Linear penalties** have problems:
```
penalty = |value - optimal| × weight

Issues:
- No "acceptable range" (every deviation penalized equally)
- Extreme values not penalized enough
- Doesn't match biological reality (thresholds matter)
```

**Piecewise logistic** is better:
```
score = 1 / (1 + exp(k × (|value - optimal| - threshold)))

Benefits:
- "Free zone" near optimal (flat region)
- Sharp penalty at biological thresholds
- Bounded output (interpretable 0-1 range)
- Matches observed PCR behavior
```

### 4.2 Scoring Function Specification

```javascript
/**
 * Master scoring function for primer pair
 *
 * Output: Score from 0-100 (higher = better)
 * All individual scores are 0-1, then weighted and summed
 */
function scorePrimerPair(fwd, rev, template, options = {}) {
  const weights = {
    // Tier 1: Critical (~50% of total)
    offTarget: 0.18,
    terminal3DG: 0.12,
    tmFwd: 0.05,
    tmRev: 0.05,
    gcFwd: 0.04,
    gcRev: 0.04,

    // Tier 2: Important (~35% of total)
    hairpinFwd: 0.05,
    hairpinRev: 0.05,
    selfDimerFwd: 0.04,
    selfDimerRev: 0.04,
    heterodimer: 0.06,
    gcClampFwd: 0.03,
    gcClampRev: 0.03,
    homopolymerFwd: 0.02,
    homopolymerRev: 0.02,

    // Tier 3: Minor (~10% of total)
    tmDiff: 0.03,
    lengthFwd: 0.01,
    lengthRev: 0.01,
    terminalBaseFwd: 0.01,
    terminalBaseRev: 0.01,
    lengthDiff: 0.01,

    // Sanger-specific (~5% of total)
    ampliconLength: 0.02,
    ampliconStructure: 0.02,
    distanceToROI: 0.01,
  };

  // Calculate individual scores (each 0-1)
  const scores = {
    offTarget: scoreOffTarget(fwd, rev, template),
    terminal3DG: scoreTerminal3DG(fwd, rev),
    tmFwd: scoreTm(fwd.tm, options),
    tmRev: scoreTm(rev.tm, options),
    // ... etc
  };

  // Weighted sum
  let totalScore = 0;
  for (const [key, weight] of Object.entries(weights)) {
    totalScore += scores[key] * weight;
  }

  // Convert to 0-100 scale
  return Math.round(totalScore * 100);
}
```

### 4.3 Individual Scoring Functions

#### 4.3.1 Tm Scoring (Piecewise Logistic)

```javascript
/**
 * Tm scoring for Sanger sequencing
 *
 * Optimal: 55-60°C (lower than RT-qPCR due to different read requirements)
 * Acceptable: 50-65°C
 * Problematic: <50°C (won't anneal) or >65°C (secondary structure)
 */
function scoreTm(tm, options = {}) {
  const {
    optimalLow = 55,
    optimalHigh = 60,
    acceptableLow = 50,
    acceptableHigh = 65,
    steepness = 0.5,
  } = options;

  // In optimal range
  if (tm >= optimalLow && tm <= optimalHigh) {
    return 1.0;
  }

  // Below optimal
  if (tm < optimalLow) {
    if (tm >= acceptableLow) {
      // Linear decay in acceptable range
      return 0.7 + 0.3 * (tm - acceptableLow) / (optimalLow - acceptableLow);
    }
    // Logistic decay below acceptable
    const excess = acceptableLow - tm;
    return 0.7 / (1 + Math.exp(steepness * excess));
  }

  // Above optimal
  if (tm > optimalHigh) {
    if (tm <= acceptableHigh) {
      return 0.7 + 0.3 * (acceptableHigh - tm) / (acceptableHigh - optimalHigh);
    }
    const excess = tm - acceptableHigh;
    return 0.7 / (1 + Math.exp(steepness * excess));
  }
}
```

#### 4.3.2 3' Terminal ΔG Scoring

```javascript
/**
 * 3' terminal ΔG scoring
 *
 * Based on experimental finding:
 * "priming was detectable when 3'-terminal portion formed duplex more stable than -11 kcal/mol"
 *
 * Optimal range: -6 to -11 kcal/mol
 * Too loose (> -6): won't initiate efficiently
 * Too tight (< -12): mispriming risk
 */
function scoreTerminal3DG(dG) {
  const optimalLow = -11;   // Most stable acceptable
  const optimalHigh = -6;   // Least stable acceptable

  if (dG >= optimalLow && dG <= optimalHigh) {
    return 1.0;
  }

  if (dG > optimalHigh) {
    // Too loose - exponential penalty
    const excess = dG - optimalHigh;
    return Math.exp(-0.3 * excess);
  }

  if (dG < optimalLow) {
    // Too tight - milder penalty (still works, just less specific)
    const excess = optimalLow - dG;
    return Math.exp(-0.15 * excess);
  }
}
```

#### 4.3.3 Tm Difference Scoring (Reduced Weight)

```javascript
/**
 * Tm difference scoring
 *
 * Key finding from PrimerScore2: "little effect" on amplification
 * We use a generous free zone and mild penalties
 */
function scoreTmDiff(tmFwd, tmRev) {
  const diff = Math.abs(tmFwd - tmRev);

  // 0-3°C: No penalty (expanded from typical 1-2°C)
  if (diff <= 3) {
    return 1.0;
  }

  // 3-5°C: Mild penalty
  if (diff <= 5) {
    return 0.9 - 0.1 * (diff - 3) / 2;
  }

  // 5-8°C: Moderate penalty
  if (diff <= 8) {
    return 0.7 - 0.2 * (diff - 5) / 3;
  }

  // >8°C: Steep penalty
  return 0.5 * Math.exp(-0.2 * (diff - 8));
}
```

---

## 5. Implementation Status & Roadmap

### ✅ COMPLETED: Phase 1a - Equilibrium Efficiency Core (Dec 2025)

**Files Created:**
- `src/lib/equilibrium.js` - Core thermodynamic calculation module
- `src/lib/equilibrium.test.js` - 34 comprehensive unit tests

**Functions Implemented:**
| Function | Purpose | Status |
|----------|---------|--------|
| `calculateEquilibriumEfficiency()` | Main entry point for primer pair evaluation | ✅ Done |
| `calculateAllSpeciesDG()` | Computes ΔG for all chemical species | ✅ Done |
| `calculateHairpinDG()` | Self-folding using Zuker algorithm | ✅ Done |
| `calculateHomodimerDG()` | Self-dimer with antiparallel alignment | ✅ Done |
| `calculateHeterodimerDG()` | Cross-dimer between primer pairs | ✅ Done |
| `calculateDuplexDG()` | Target binding with NN model | ✅ Done |
| `calculateOffTargetDG()` | Off-target binding detection | ✅ Done |
| `solveEquilibrium()` | Partition function solver | ✅ Done |
| `efficiencyToScore()` | Piecewise logistic conversion | ✅ Done |

**Key Technical Decisions:**
- Uses **DNA24 parameters** (Greenleaf Lab 2024) - 50% more accurate than SantaLucia 1998
- Implements **antiparallel alignment** correctly for dimer calculations
- Returns **loss breakdown** by species (hairpin, homodimer, heterodimer, offTarget, free)

### ✅ COMPLETED: Phase 1b - Integration with Primer Scoring (Dec 2025)

**Files Modified:**
- `src/lib/primers.js` - Added equilibrium integration
- `src/lib/primers.test.js` - 3 new integration tests

**Integration Points:**
| Feature | Implementation | Status |
|---------|---------------|--------|
| Import equilibrium module | Added to primers.js | ✅ Done |
| Extend Scoring object | Added equilibriumEfficiency, equilibriumScore, equilibriumLosses | ✅ Done |
| `addEquilibriumScoring()` method | New PrimerFactory method | ✅ Done |
| `score()` function options | Added `includeEquilibrium`, `annealingTemperature` | ✅ Done |
| dict() output | Extended with equilibrium metrics | ✅ Done |

**Usage Example:**
```javascript
const [fwd, rev] = score(fwdSeq, revSeq, template, template, {
  includeEquilibrium: true,
  annealingTemperature: 55,
});

console.log(fwd.scoring.equilibriumEfficiency);  // 0.85
console.log(fwd.scoring.equilibriumScore);       // 82.5
console.log(fwd.scoring.equilibriumLosses);      // { hairpin: 0.001, ... }
```

---

### ✅ COMPLETED: Phase 2 - Piecewise Logistic Refinement (Dec 2025)

**Files Created:**
- `src/lib/scoring.js` - Piecewise logistic scoring functions module
- `src/lib/scoring.test.js` - 58 comprehensive unit tests

**Functions Implemented:**
| Function | Purpose | Status |
|----------|---------|--------|
| `piecewiseLogistic()` | Generic piecewise logistic curve | ✅ Done |
| `scoreTm()` | Tm scoring (55-60°C optimal) | ✅ Done |
| `scoreGc()` | GC content scoring (40-60% optimal) | ✅ Done |
| `scoreTerminal3DG()` | 3' terminal ΔG (-6 to -11 kcal/mol) | ✅ Done |
| `scoreTmDiff()` | Tm difference with reduced weight | ✅ Done |
| `scoreHairpin()` | Hairpin penalty score | ✅ Done |
| `scoreHomodimer()` | Homodimer penalty score | ✅ Done |
| `scoreHeterodimer()` | Heterodimer penalty score | ✅ Done |
| `scoreOffTarget()` | Off-target exponential penalty | ✅ Done |
| `scoreLength()` | Primer length scoring | ✅ Done |
| `scoreGcClamp()` | GC clamp evaluation | ✅ Done |
| `scoreHomopolymer()` | Homopolymer run penalty | ✅ Done |
| `scoreAmpliconLength()` | Sanger amplicon scoring | ✅ Done |
| `calculateCompositeScore()` | Weighted combination | ✅ Done |
| `classifyQuality()` | Quality tier classification | ✅ Done |

**Integration with primers.js:**
- Added `addPiecewiseScoring()` method to PrimerFactory
- Extended `score()` function to include piecewise scores
- Extended `dict()` output with `piecewise_scores`, `composite_score`, `quality_tier`

**Sanger-specific features** (completed in Phase 3):
- [x] Distance to region of interest (ROI) - `scoreDistanceToROI()`
- [x] Amplicon secondary structure penalty - `scoreAmpliconStructure()`

### ✅ COMPLETED: Phase 3 - Off-Target & Sanger Features (Dec 2025)

**Files Created:**
- `src/lib/offTargetClassification.js` - Type A-F classification
- `src/lib/offTargetClassification.test.js` - 31 comprehensive tests

**Functions Implemented:**
| Function | Purpose | Status |
|----------|---------|--------|
| `classifyOffTarget()` | Type A-F classification | ✅ Done |
| `findTypeA()` | Exact match detection | ✅ Done |
| `findTypeB()` | Antisense binding detection | ✅ Done |
| `findTypeC()` | Partial 3' homology detection | ✅ Done |
| `findTypeD()` | Internal homology detection | ✅ Done |
| `classifyAllOffTargets()` | Full classification with penalties | ✅ Done |
| `scoreDistanceToROI()` | Sanger ROI distance scoring | ✅ Done |
| `scoreAmpliconStructure()` | Sanger amplicon structure scoring | ✅ Done |

### ✅ COMPLETED: Phase 4a - Weight Calibration Framework (Dec 2025)

**Files Created:**
- `src/lib/weightCalibration.js` - Complete calibration module
- `src/lib/weightCalibration.test.js` - 34 comprehensive tests

**Functions Implemented:**
| Function | Purpose | Status |
|----------|---------|--------|
| `calculateMetrics()` | Accuracy, precision, recall, F1 | ✅ Done |
| `calculateAUC()` | AUC-ROC calculation | ✅ Done |
| `crossValidate()` | K-fold cross-validation | ✅ Done |
| `gridSearch()` | Exhaustive weight optimization | ✅ Done |
| `coordinateDescent()` | Iterative weight tuning | ✅ Done |
| `compareWeights()` | A/B weight comparison | ✅ Done |
| `generateCalibrationReport()` | Comprehensive report | ✅ Done |

### ✅ COMPLETED: Phase 4b - Calibration with Döring/openPrimeR Dataset (Dec 2025)

**Dataset:** Döring et al. immunoglobulin PCR evaluation (openPrimeR)
- Source: Figshare, CC-BY license
- 829 primer-template pairs with binary outcomes (Amplified/Unamplified)
- 365 amplified (44%), 464 unamplified (56%)
- Features: Tm, GC%, ΔG, mismatches, 3' position, hairpin, self-dimer

**Calibration Script:** `scripts/calibrateDoring.js`

**Results Summary:**

| Metric | Current Weights | Optimized Weights | Change |
|--------|-----------------|-------------------|--------|
| F1 Score | 77.1% | 81.9% | +4.8% |
| AUC-ROC | 0.831 | 0.848 | +0.017 |
| Accuracy | 80.2% | - | - |
| Precision | 78.6% | - | - |
| Recall | 75.6% | - | - |
| Optimal Threshold | 80 | 80 | - |

**Cross-Validation (5-fold):**
- Mean F1: 81.7% ± 4.0%
- Mean AUC: 0.840 ± 0.037

**Most Discriminative Features:**
1. `offTarget` (diff=0.515) - Mismatch/specificity is critical
2. `terminal3DG` (diff=0.487) - 3' binding strength matters
3. `tmFwd` (diff=-0.185) - Higher Tm correlated with failure

**Optimized Weight Changes:**
```
heterodimer:  0.10 → 0.04  (reduced - less critical in this context)
terminal3DG:  0.05 → 0.10  (doubled - more important than expected)
selfDimerFwd: 0.04 → 0.08  (doubled)
```

**Key Finding:** Off-target/mismatch features are the primary predictors of PCR success in the immunoglobulin primer context, not primer-intrinsic properties like Tm or GC content (which are constant per primer across templates)

### ✅ COMPLETED: Phase 4c - Cross-Domain Validation on Kayama Dataset (Dec 2025)

**Dataset:** Kayama et al. 16S rRNA PCR evaluation (Scientific Reports, 2021)
- Source: Paper Tables 1, 2, 4 - "Prediction of PCR amplification using RNN"
- 72 primer sets × 31 templates = 2,232 primer-template pairs
- 506 amplified (22.7%), 1,726 not amplified (77.3%)
- **Key difference from Döring:** Raw sequences only, no pre-computed features

**Validation Script:** `scripts/validateKayama.js`

**Approach:**
- Computed all primer features using our scoring algorithms:
  - Tm via `calculateTmQ5()` (NEB Q5 algorithm)
  - GC content, length, GC clamp, homopolymer via sequence analysis
  - 3' terminal ΔG via `calculate3primeTerminalDG()`
  - Hairpin/homodimer ΔG via heuristic estimation
- Applied Döring-calibrated weights without re-training
- This tests true generalization to an independent dataset

**Cross-Domain Validation Results:**

| Metric | Aggregated (72 pairs) | Expanded (2,232 pairs) |
|--------|----------------------|------------------------|
| Accuracy | 45.8% | 46.1% |
| Precision | 23.5% | 28.0% |
| Recall | 100.0% | 87.4% |
| F1 Score | 38.1% | 42.4% |
| **AUC-ROC** | **0.682** | **0.643** |
| 5-fold CV F1 | 35.1% ± 21.3% | 42.3% ± 1.0% |

**Feature Discrimination (Kayama dataset):**
1. `hairpinFwd` (diff=0.379) - Strongest discriminator
2. `gcFwd` (diff=0.081) - GC content matters
3. `tmFwd` (diff=0.046) - Modest Tm effect

**Interpretation:**
- AUC-ROC > 0.5 indicates **some cross-domain predictive power**
- Low F1/accuracy is largely due to **heavy class imbalance** (77% failures)
- The Kayama dataset includes intentionally poor primers for RNN training
- Template-specific interactions are not captured by primer-only features

**Key Insights:**
1. Döring weights **partially transfer** to 16S rRNA context (AUC=0.64)
2. Hairpin prediction is strong discriminator in both datasets
3. Off-target/mismatch data (missing in Kayama) is critical for high accuracy
4. Primer scoring alone cannot predict template-specific amplification

**Recommendation:** For robust cross-domain validation, include off-target analysis via BLAST or similar. Pure sequence-based scoring provides ~64% AUC baseline.

### 🔲 TODO: Phase 5 - Validation

- [ ] In-silico validation on PrimerBank subset
- [ ] Partner with sequencing core for wet-lab validation
- [ ] Calculate precision, recall, AUC-ROC
- [ ] Target: >85% accuracy, >90% recall

### 🔄 IN PROGRESS: Phase 6 - Manuscript Preparation

**Figure Data Generated:** `paper/figures/` (Dec 2025)

| File | Description | Status |
|------|-------------|--------|
| `roc_doring.csv` | ROC curve data (829 points) | ✅ Generated |
| `feature_importance.csv` | Feature discrimination data | ✅ Generated |
| `figures_summary.md` | Preliminary figures summary | ✅ Generated |

**Key Figure Metrics (Döring dataset):**
- AUC-ROC: **0.824**
- Top discriminative features:
  1. `offTarget` (+0.515 diff) - Mismatch/specificity is critical
  2. `terminal3DG` (+0.487 diff) - 3' binding strength matters
  3. `tmFwd/tmRev` (-0.185 diff) - Higher Tm correlated with failure

**Manuscript Preparation Checklist:**
- [x] Generate ROC curve data
- [x] Generate feature importance data
- [x] Create preliminary figures summary
- [ ] Write methods section
- [ ] Write results and discussion
- [ ] Prepare supplementary materials

---

## 6. Manuscript Outline

### Title
"Empirically-Calibrated Primer Scoring for Plasmid PCR and Sanger Sequencing: A Piecewise Logistic Approach"

### Abstract (150-250 words)
- Problem: Current primer scoring uses heuristic weights
- Method: Meta-analysis of 4 validation studies → piecewise logistic model
- Results: X% accuracy (vs Y% baseline)
- Conclusion: First empirically-calibrated tool for plasmid/Sanger context

### 1. Introduction
- 1.1 Importance of primer design
- 1.2 Limitations of current tools
- 1.3 Gap: No plasmid/Sanger-specific calibration
- 1.4 Our contribution

### 2. Methods
- 2.1 Literature meta-analysis (4 studies)
- 2.2 Feature selection and importance ranking
- 2.3 Piecewise logistic scoring framework
- 2.4 Off-target definition for plasmid context
- 2.5 Validation dataset and metrics

### 3. Results
- 3.1 Feature importance from meta-analysis
- 3.2 Optimal weight distribution
- 3.3 In-silico validation
- 3.4 Wet-lab validation
- 3.5 Comparison with Primer3

### 4. Discussion
- 4.1 Why Tm difference was over-weighted
- 4.2 Why off-target was under-weighted
- 4.3 Limitations
- 4.4 Future work

### 5. Conclusion

### Supplementary Materials
- S1: Complete scoring formulas
- S2: Validation dataset
- S3: Source code availability

---

## 7. Figures and Tables

### Figure 1: Scoring Framework Overview
- Flowchart showing: Input → Feature extraction → Piecewise scoring → Weighted sum → Quality tier

### Figure 2: Piecewise Logistic vs Linear
- Side-by-side comparison showing biological threshold handling

### Figure 3: Feature Importance Ranking
- Bar chart from meta-analysis of 4 studies
- Highlight: Off-target >> Tm difference

### Figure 4: ROC Curves
- Our model vs Primer3 vs random baseline

### Figure 5: Validation Results
- Confusion matrix
- Success rate by quality tier

### Table 1: Summary of Validation Studies
- PrimerBank, GM1, PrimerScore2, ML-PCR
- Dataset size, context, key findings

### Table 2: Off-Target Scenario Classification
- Type A-F with risk levels and detection methods

### Table 3: Empirical Weight Distribution
- Feature, evidence source, weight, rationale

### Table 4: Validation Metrics
- Accuracy, precision, recall, F1, AUC

---

## 8. Differentiation from Existing Work

### Novel Contributions

| Aspect | Existing Tools | Our Contribution |
|--------|---------------|------------------|
| Weight source | Heuristic | Meta-analysis of 4 studies |
| Scoring function | Linear | Piecewise logistic |
| Off-target scope | Genome-wide | Input-scoped (plasmid) |
| Off-target types | Exact match only | 6 types (A-F) classified |
| Validation metric | Depth ratio / ad-hoc | Binary success (precision/recall) |
| Sanger-specific | None | Amplicon analysis, ROI distance |
| Heterodimer | Post-hoc check | Integrated in primary score |

### Why This Matters

1. **Reproducibility**: Weights are derived from data, not intuition
2. **Appropriate context**: Optimized for plasmid/Sanger, not genome-scale
3. **Interpretability**: Piecewise functions have biological meaning
4. **Transparency**: Full methodology published with code

---

## 9. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Insufficient validation data | Medium | High | Partner with core facility early |
| No improvement over Primer3 | Low | High | Pre-test with in-silico validation |
| Overfitting to validation set | Medium | Medium | Use cross-validation, holdout set |
| Reviewer skepticism on meta-analysis | Medium | Low | Document literature review methodology |
| Wet-lab validation delays | High | Medium | Start partnerships in Phase 1 |

---

## 10. Comprehensive Comparison: All 5 Validation Studies

### 10.1 Head-to-Head Comparison

| Aspect | PrimerBank | GM1 Model | PrimerScore2 | ML-PCR | **Pythia** |
|--------|------------|-----------|--------------|--------|------------|
| **Year** | 2010 | 2008 | 2022 | 2021 | 2009 |
| **Journal** | NAR | NAR | Sci Rep | bioRxiv | NAR |
| **Dataset Size** | 26,855 pairs | 1,314 pairs | 57-plex | 290 reactions | ~1,000 loci |
| **Application** | RT-qPCR | Genotyping | NGS multiplex | General PCR | DNase-seq |
| **Validation Metric** | Success rate | Failure rate | Depth ratio | Binary | Binary |
| **Precision** | ~82.6%* | N/A | N/A | 81% | **81%** |
| **Recall** | N/A | N/A | N/A | N/A | **97%** |
| **Accuracy** | 82.6% | 83-94%** | R²=0.935 | 81% | ~89%*** |
| **Scoring Method** | Heuristic | Statistical | Piecewise logistic | Random Forest | **Thermodynamic** |
| **Off-target** | Transcriptome | Genome (dominant) | Amplicon panel | N/A | Genome index |
| **Open Source** | No | No | Yes (Perl) | No | Yes (C++) |

*Success rate ≈ precision when most primers work
**Failure rate reduced from 17% to 6%
***Coverage in difficult regions

### 10.2 Suitability for Plasmid/Sanger Context

| Criterion | PrimerBank | GM1 | PrimerScore2 | ML-PCR | **Pythia** |
|-----------|------------|-----|--------------|--------|------------|
| Binary success metric | ⚠️ Success rate | ✅ Failure rate | ❌ Depth ratio | ✅ Binary | ✅ Binary |
| Monoplex context | ❌ qPCR | ✅ Single | ❌ Multiplex | ✅ General | ✅ Single |
| Off-target scope | ❌ Large | ❌ Genome | ✅ Input | N/A | ⚠️ Genome |
| Thermodynamic basis | ❌ Heuristic | ⚠️ Partial | ⚠️ Features | ❌ ML black box | ✅ Full |
| Heterodimer modeled | ❌ No | ❌ No | ⚠️ Partial | ❌ No | ✅ Yes |
| Sanger-specific | ❌ No | ❌ No | ❌ No | ❌ No | ❌ No |
| Implementation ease | ❌ Closed | ❌ Closed | ⚠️ Perl | ❌ Closed | ⚠️ C++ |

**Legend**: ✅ Good fit | ⚠️ Partial fit | ❌ Poor fit

### 10.3 Key Insights by Study

| Study | Key Contribution | Key Limitation |
|-------|------------------|----------------|
| **PrimerBank** | Largest validated dataset (82.6% baseline) | RT-qPCR specific, no scoring formula |
| **GM1** | Off-target is dominant failure factor | Genome-scale, overkill for plasmid |
| **PrimerScore2** | Piecewise logistic math, feature importance | Multiplex depth ratio, not binary |
| **ML-PCR** | Real-world data, feature ranking | Small dataset (290), preprint |
| **Pythia** | Thermodynamic equilibrium, 97% recall | Genome-scale, computationally heavy |

### 10.4 Recommended Integration Strategy

Based on the comprehensive analysis, here is the recommended hybrid approach:

```
┌─────────────────────────────────────────────────────────────┐
│                    INTEGRATED MODEL                         │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  FROM PYTHIA (Core Framework):                              │
│  ├── Equilibrium efficiency as primary metric               │
│  ├── All chemical species modeled                           │
│  ├── min(η_left, η_right) bottleneck principle              │
│  └── Thermodynamic basis (update to DNA24 parameters)       │
│                                                             │
│  FROM PRIMERSCORE2 (Scoring Function):                      │
│  ├── Piecewise logistic curves (not linear penalties)       │
│  ├── Feature importance rankings                            │
│  └── Weighted sum normalization                             │
│                                                             │
│  FROM GM1 (Off-target):                                     │
│  ├── Off-target as dominant factor (highest weight)         │
│  └── Exponential penalty scaling                            │
│                                                             │
│  FROM ML-PCR (Feature Insights):                            │
│  ├── Tm and GC clamp as top features                        │
│  ├── Tm difference has "little effect"                      │
│  └── Real-world validation approach                         │
│                                                             │
│  NOVEL (Our Contribution):                                  │
│  ├── Input-scoped off-target (not genome)                   │
│  ├── Sanger-specific features (amplicon, ROI distance)      │
│  ├── JavaScript/web implementation                          │
│  └── DNA24 thermodynamic parameters (2024)                  │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### 10.5 Why Pythia Equilibrium Model is Best Foundation

| Reason | Explanation |
|--------|-------------|
| **Physically grounded** | Based on thermodynamics, not heuristics |
| **All interactions modeled** | Hairpin, homodimer, heterodimer, off-target in one framework |
| **Conservative approach** | "If it works at equilibrium, it works in PCR" |
| **Best recall (97%)** | Won't reject good primers (unlike Primer3's 48%) |
| **Binary validation** | Matches our success/failure use case |
| **Proven in difficult regions** | 89% coverage vs 51% for Primer3 |

**However**, Pythia needs adaptation:
1. Replace genome indexing with input-scoped search
2. Update from 2009 SantaLucia parameters to DNA24 (2024)
3. Port from C++ to JavaScript for web deployment
4. Add Sanger-specific features not in original

### 10.6 Proposed Hybrid: "Equilibrium Score"

Combine Pythia's thermodynamic model with PrimerScore2's scoring:

```javascript
/**
 * Equilibrium-based primer scoring
 *
 * Combines:
 * - Pythia: Equilibrium efficiency calculation
 * - PrimerScore2: Piecewise logistic normalization
 * - GM1: Off-target dominance
 */
function calculateEquilibriumScore(fwd, rev, template, options = {}) {
  // Step 1: Calculate equilibrium concentrations (Pythia-style)
  const species = calculateAllSpecies(fwd, rev, template, options.temperature);

  // Step 2: Calculate equilibrium efficiency (Pythia formula)
  const η_fwd = species.fwdBoundTarget / species.fwdTotal;
  const η_rev = species.revBoundTarget / species.revTotal;
  const equilibriumEfficiency = Math.min(η_fwd, η_rev);

  // Step 3: Convert to score using piecewise logistic (PrimerScore2-style)
  const efficiencyScore = piecewiseLogistic(equilibriumEfficiency, {
    optimal: 0.95,      // 95% bound to target
    acceptable: 0.80,   // 80% acceptable
    steepness: 5.0,
  });

  // Step 4: Additional feature scores (with empirical weights)
  const featureScores = {
    offTarget: scoreOffTarget(fwd, rev, template),      // GM1: dominant
    tmMatch: scoreTmMatch(fwd, rev),                    // PrimerScore2
    gcContent: scoreGcContent(fwd, rev),                // PrimerScore2
    ampliconLength: scoreAmpliconLength(fwd, rev),      // Sanger-specific
  };

  // Step 5: Weighted combination
  // Equilibrium efficiency is primary (50%), features are secondary (50%)
  const weights = {
    equilibrium: 0.50,
    offTarget: 0.20,
    tmMatch: 0.10,
    gcContent: 0.08,
    ampliconLength: 0.05,
    // ... other features: 0.07
  };

  const totalScore =
    weights.equilibrium * efficiencyScore +
    weights.offTarget * featureScores.offTarget +
    weights.tmMatch * featureScores.tmMatch +
    weights.gcContent * featureScores.gcContent +
    weights.ampliconLength * featureScores.ampliconLength;

  return {
    score: Math.round(totalScore * 100),
    equilibriumEfficiency,
    species,
    featureScores,
    quality: classifyQuality(totalScore),
  };
}

/**
 * Calculate all chemical species at equilibrium
 *
 * Species:
 * - Primer_unfolded: Free primer
 * - Primer_hairpin: Self-folded
 * - Primer_homodimer: Self-dimer
 * - Primer_target: Bound to intended site
 * - Primer_offtarget: Bound elsewhere
 * - Fwd_Rev_heterodimer: Cross-dimer
 */
function calculateAllSpecies(fwd, rev, template, temperature = 55) {
  // Initial concentrations (typical PCR)
  const C_primer = 0.5e-6;   // 500 nM each primer
  const C_template = 1e-9;   // 1 nM template (late-stage PCR)

  // Calculate ΔG for each interaction
  const dG = {
    fwd_hairpin: calculateHairpinDG(fwd.sequence, temperature),
    fwd_homodimer: calculateHomodimerDG(fwd.sequence, temperature),
    fwd_target: calculateDuplexDG(fwd.sequence, fwd.bindingSite, temperature),
    fwd_offtarget: calculateOffTargetDG(fwd.sequence, template, temperature),
    rev_hairpin: calculateHairpinDG(rev.sequence, temperature),
    rev_homodimer: calculateHomodimerDG(rev.sequence, temperature),
    rev_target: calculateDuplexDG(rev.sequence, rev.bindingSite, temperature),
    rev_offtarget: calculateOffTargetDG(rev.sequence, template, temperature),
    heterodimer: calculateHeterodimerDG(fwd.sequence, rev.sequence, temperature),
  };

  // Solve equilibrium using gradient descent on Gibbs energy
  // (simplified version - full implementation would use proper optimization)
  const K = {};
  for (const [key, value] of Object.entries(dG)) {
    K[key] = Math.exp(-value / (R * (temperature + 273.15)));
  }

  // Return equilibrium concentrations
  // (actual implementation requires solving coupled equilibria)
  return solveEquilibrium(C_primer, C_template, K);
}
```

---

## 11. References (Key Papers)

1. **Pythia** - Mann et al. (2009) - Nucleic Acids Research
   - "A thermodynamic approach to PCR primer design"
   - DOI: 10.1093/nar/gkp443
   - **Key contribution**: Equilibrium efficiency model, 97% recall

2. PrimerScore2 (2022) - Nature Scientific Reports
   - DOI: 10.1038/s41598-022-25561-z
   - **Key contribution**: Piecewise logistic scoring, feature importance

3. GM1 Model (2008) - Nucleic Acids Research
   - DOI: 10.1093/nar/gkn290
   - **Key contribution**: Off-target as dominant failure factor

4. ML-PCR Optimization (2021) - bioRxiv
   - DOI: 10.1101/2021.08.12.455589
   - **Key contribution**: Real-world feature importance

5. PrimerBank (2010) - Nucleic Acids Research
   - DOI: 10.1093/nar/gkp1014
   - **Key contribution**: 82.6% success rate baseline

6. RNN PCR Prediction (2021) - Scientific Reports
   - DOI: 10.1038/s41598-021-86357-1

7. DNA24 Parameters - Greenleaf Lab (2024)
   - **Key contribution**: Updated thermodynamic parameters

8. IDT OligoAnalyzer Guidelines
   - https://www.idtdna.com/pages/support/faqs/

9. PRIMEval (2019) - Scientific Reports
   - DOI: 10.1038/s41598-019-55883-4

---

## 11. Appendix: Current System Analysis

### Current Penalty Weights (from codebase analysis)

```javascript
// primers.js - PCR
penaltyTm = 1.0
penaltyGc = 0.2 × 100 = 20 (normalized)
penaltyLen = 0.5
penaltyDg = 2.0
penaltyOffTarget = 20.0
penaltyTmDiff = 1.0

// mutagenesis.js - SDM
Tm difference: 2.0 × (tmDiff - 1.0)²
Tm floor: 10.0 × (minTm - tm)²
GC: (asymmetric) 150× or 80× for high/low
Length: 1.5 per bp over 28

// sequencing.js
weights.offTarget = 15.0
weights.hairpin = 4.0
weights.selfDimer = 5.0
weights.gc = 1.5
weights.tm = 2.0
```

### Issues Identified

1. **Tm difference over-weighted**: Literature shows "little effect"
2. **Off-target under-weighted**: Should be dominant (GM1)
3. **No heterodimer in primary scoring**: Only post-hoc check
4. **Linear penalties**: Don't match biological thresholds
5. **No Sanger-specific features**: Amplicon length, ROI distance

---

*Document prepared for internal use. Will be refined during implementation.*
