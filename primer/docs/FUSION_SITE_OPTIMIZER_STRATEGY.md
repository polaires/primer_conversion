# Fusion Site Optimizer Strategy Document

## Executive Summary

This document outlines the strategy for implementing a **SplitSet-like Fusion Site Optimizer** that will search within DNA sequences to find optimal junction positions for Golden Gate assembly. The goal is to maximize user success rate by jointly optimizing:

1. Overhang fidelity (set-level, not just individual)
2. Primer quality at each junction
3. Fragment balance and practical constraints
4. Biological context (coding frames, domains)

---

## Table of Contents

1. [Implementation Status](#1-implementation-status)
2. [Two Optimization Approaches](#2-two-optimization-approaches)
3. [Architecture Overview](#3-architecture-overview)
4. [Implemented Components](#4-implemented-components)
5. [Integration with Existing Tools](#5-integration-with-existing-tools)
6. [Natural Integration with Primer Design](#6-natural-integration-with-primer-design)
7. [Primer Boundary Optimizer (Edge Cases)](#7-primer-boundary-optimizer-edge-cases)
8. [UI Implementation Strategy](#8-ui-implementation-strategy)
9. [API Reference](#9-api-reference)
10. [Testing Coverage](#10-testing-coverage)
11. [Future Enhancements](#11-future-enhancements)

---

## 1. Implementation Status

### Core Backend - COMPLETE ✅

| Component | File | Status | Tests |
|-----------|------|--------|-------|
| **Overhang Efficiency** | `overhang-efficiency.js` | ✅ Complete | ✅ |
| **Site Creation Check** | `site-creation-check.js` | ✅ Complete | ✅ |
| **Fusion Site Scanner** | `fusion-site-scanner.js` | ✅ Complete | ✅ |
| **Fusion Site Scorer** | `fusion-site-scorer.js` | ✅ Complete | ✅ |
| **Scar Preferences** | `scar-preferences.js` | ✅ Complete | ✅ |
| **Failure Prediction** | `failure-prediction.js` | ✅ Complete | ✅ |
| **Main Optimizer** | `fusion-site-optimizer.js` | ✅ Complete | ✅ |
| **Primer Boundary Optimizer** | `primer-boundary-optimizer.js` | ✅ Complete | 🔄 Pending |

### Integration with Existing Tools - COMPLETE ✅

| Integration | Status | Notes |
|-------------|--------|-------|
| `getOverhangFidelityExperimental()` | ✅ Integrated | From goldengate.js |
| `calculateTmQ5()` | ✅ Integrated | NEB Q5 nearest-neighbor algorithm |
| `calculateHairpinDG()` | ✅ Integrated | From equilibrium.js |
| `calculateHomodimerDG()` | ✅ Integrated | From equilibrium.js |
| `calculate3primeTerminalDG()` | ✅ Integrated | From tmQ5.js |
| Scoring functions | ✅ Integrated | From scoring.js |

### UI Implementation - COMPLETE ✅

| Component | Status | Priority |
|-----------|--------|----------|
| Fusion Site Optimizer Panel | ✅ Built | High |
| Junction Candidate Cards | ✅ Built | High |
| Failure Prediction Display | ✅ Built | High |
| Scoring Weights Editor | ✅ Built | Medium |
| **Backend Connection** | ✅ Complete | **Critical** |
| **Primer Design Integration** | ✅ Complete | **Critical** |
| **Natural Split-in-Place** | ✅ Complete | **Critical** |
| Export/Download | ❌ Pending | Low |

**Two Integration Modes:**
1. **Natural Integration (Recommended)**: Click the "Split" button on any large sequence (≥500bp) in the parts list to optimize junction positions and split into fragments inline
2. **Advanced Mode**: Toggle "Fusion Site Optimizer" for full control with algorithm selection, scoring weight customization, and detailed failure prediction

---

## 2. Two Optimization Approaches

This system provides **two distinct optimization tools** for different use cases:

### 2.1 Fusion Site Optimizer (Finding Junction Positions)

**Use case**: You have a **single long sequence** and need to determine WHERE to split it into fragments.

```
INPUT: Single sequence (e.g., 5000bp gene)
       ─────────────────────────────────────────────────────

OPTIMIZER FINDS: Optimal junction positions from scratch
       ────────────|──────────|──────────|──────────────────
                  J1         J2         J3
                  ↑          ↑          ↑
           Positions chosen by algorithm for:
           • High-fidelity overhangs
           • Good primer binding regions
           • Balanced fragment sizes

OUTPUT: Multiple fragments with optimized junctions
       [Fragment 1] [Fragment 2] [Fragment 3] [Fragment 4]
```

**Files**: `fusion-site-optimizer.js`, `fusion-site-scanner.js`, `fusion-site-scorer.js`

### 2.2 Primer Boundary Optimizer (Adjusting User-Defined Boundaries)

**Use case**: You have **multiple fragments with pre-defined boundaries** and the primers at those boundaries have poor quality.

```
INPUT: User-defined fragments with fixed boundaries
       [Fragment A]|[Fragment B]|[Fragment C]
                   ↑            ↑
            User's boundaries (might have poor primer regions)

PROBLEM: Primer at junction has poor Tm, hairpins, or other issues
       ...NNNNN[AAAA] | [GCGC]NNNNN...    ← Poor: poly-A, low Tm
              ↑ primer binding region

OPTIMIZER SHIFTS: Boundary left or right to get better primer
       ...NN[AGCTGC] | [TGCA]NNNNNNNN...  ← Better: balanced GC
              ↑ improved binding region

OUTPUT: Adjusted boundaries with better primers
       [Fragment A'] [Fragment B'] [Fragment C']
       (total sequence unchanged, boundaries shifted for primer quality)
```

**Files**: `primer-boundary-optimizer.js` (NEW)

### 2.3 Key Differences

| Aspect | Fusion Site Optimizer | Primer Boundary Optimizer |
|--------|----------------------|--------------------------|
| **Input** | Single sequence | Multiple user-defined fragments |
| **Goal** | Find WHERE to split | Adjust EXISTING boundaries |
| **User control** | Algorithm chooses positions | User defines positions, optimizer fine-tunes |
| **Boundary changes** | Creates new boundaries | Shifts existing boundaries ±50bp |
| **Use when** | Splitting a gene into fragments | User-defined fragments have primer issues |

### 2.4 When to Use Which

**Use Fusion Site Optimizer when:**
- You have a single long sequence (gene, plasmid region)
- You want the algorithm to find optimal split points
- You don't have pre-defined fragment boundaries

**Use Primer Boundary Optimizer when:**
- You have multiple fragments already defined
- The primer binding regions at junctions are poor quality
- You want to keep approximate boundaries but improve primers
- You're working with edge cases where fragments meet

---

## 3. Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           USER INTERFACE LAYER                              │
│                      (FusionSiteDesigner.jsx - NEW)                         │
│                                                                             │
│   ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  │
│   │ Sequence     │  │ Junction     │  │ Results      │  │ Export       │  │
│   │ Input        │  │ Visualizer   │  │ Panel        │  │ Options      │  │
│   └──────────────┘  └──────────────┘  └──────────────┘  └──────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         ORCHESTRATION LAYER                                 │
│                                                                             │
│   optimizeFusionSites()  ─────────────────────────────────────────────────  │
│     • Input validation                                                      │
│     • Algorithm selection (auto/greedy/monteCarlo/branchBound/hybrid)       │
│     • Result formatting with failure prediction                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        ▼                           ▼                           ▼
┌─────────────────┐       ┌─────────────────┐       ┌─────────────────────┐
│   CANDIDATE     │       │   SCORING       │       │   ANALYSIS          │
│   GENERATOR     │       │   ENGINE        │       │   ENGINE            │
│                 │       │                 │       │                     │
│ scanForFusion-  │       │ scoreFusionSite │       │ predictFailure-     │
│ Sites()         │       │ Composite()     │       │ Modes()             │
│                 │       │                 │       │                     │
│ generateTarget- │       │ quickScoreFu-   │       │ quickRisk-          │
│ Positions()     │       │ sionSite()      │       │ Assessment()        │
└─────────────────┘       └─────────────────┘       └─────────────────────┘
        │                           │                           │
        └───────────────────────────┼───────────────────────────┘
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          OPTIMIZER LAYER                                    │
│                                                                             │
│   Strategy Pattern - auto-selected based on problem size:                   │
│                                                                             │
│   ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────────────┐    │
│   │ Greedy          │  │ Branch & Bound  │  │ Monte Carlo             │    │
│   │ (any size)      │  │ (≤5 junctions)  │  │ (>10 junctions)         │    │
│   │                 │  │                 │  │                         │    │
│   │ Fast baseline   │  │ Guarantees      │  │ Simulated annealing     │    │
│   │                 │  │ global optimum  │  │ for large assemblies    │    │
│   └─────────────────┘  └─────────────────┘  └─────────────────────────┘    │
│                                    │                                        │
│                          ┌─────────────────┐                                │
│                          │ Hybrid          │                                │
│                          │ (5-10 junctions)│                                │
│                          │                 │                                │
│                          │ Combines all    │                                │
│                          │ approaches      │                                │
│                          └─────────────────┘                                │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            DATA LAYER                                       │
│                                                                             │
│   EXISTING (Reused):                    │   NEW (Implemented):              │
│   • ligation-data.json (matrices)       │   • EFFICIENCY_PENALTIES         │
│   • GOLDEN_GATE_ENZYMES                 │   • SCAR_PREFERENCES             │
│   • getOverhangFidelityExperimental()   │   • FAILURE_MODES                │
│   • calculateTmQ5()                     │   • DEFAULT_FUSION_WEIGHTS       │
│   • Scoring functions from scoring.js   │   • OPTIMIZER_DEFAULTS           │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Implemented Components

### 3.1 overhang-efficiency.js

Implements TNNA soft penalties and efficiency scoring.

**Key Exports:**
- `calculateEfficiency(overhang)` - Calculate efficiency factor with warnings
- `calculateSetEfficiency(overhangs)` - Batch efficiency for overhang set
- `isPalindrome(overhang)` - Self-ligation risk check
- `isHomopolymer(overhang)` - Poor fidelity check
- `categorizeEfficiency(efficiency)` - Quality categorization

**Technical Notes:**
- TNNA pattern: ~70% efficiency (soft penalty, NOT hard filter)
- Pattern penalties compound with specific penalties (reduced overlap at 30%)
- Integrates with experimental data from Pryor et al. 2020

### 3.2 site-creation-check.js

Detects if junction positions CREATE new restriction sites.

**Key Exports:**
- `checkSiteCreation(sequence, position, enzyme)` - Single junction check
- `checkMultipleSiteCreation(sequence, positions, enzyme)` - Batch check
- `findSafeJunctionNear(sequence, targetPos, enzyme)` - Find alternative
- `validateJunctionSet(sequence, junctions, enzyme)` - Full validation

**Technical Notes:**
- Uses enzyme-specific spacer lengths (not hardcoded)
- Checks both forward and reverse primer tails
- Distinguishes designed sites from problematic sites

### 3.3 fusion-site-scanner.js

Scans sequences for candidate junction positions.

**Key Exports:**
- `scanForFusionSites(sequence, options)` - Main scanning function
- `scanAndRankFusionSites(sequence, options)` - Scan with scoring
- `generateTargetPositions(seqLength, numFragments)` - Target calculation
- `filterByDistance(candidates, minDistance)` - Spacing filter
- `assessFeasibility(seqLength, numFragments)` - Feasibility check

**Technical Notes:**
- Uses centralized `getOverhangFidelityExperimental()` from goldengate.js
- Hard filters: palindromes, homopolymers
- Soft scoring: TNNA, GC content, fidelity

### 3.4 fusion-site-scorer.js

Composite scoring function combining all factors.

**Key Exports:**
- `scoreFusionSiteComposite(seq, pos, enzyme, options)` - Full scoring
- `scoreMultipleFusionSites(seq, positions, enzyme)` - Batch scoring
- `quickScoreFusionSite(seq, pos, enzyme)` - Fast scoring
- `DEFAULT_FUSION_WEIGHTS` - Default scoring weights

**Default Weights:**
```javascript
{
  overhangQuality: 0.20,   // Fidelity + efficiency
  forwardPrimer: 0.20,     // Primer quality downstream
  reversePrimer: 0.20,     // Primer quality upstream
  riskFactors: 0.25,       // Site creation, mispriming
  biologicalContext: 0.15, // Codons, domains, scars
}
```

**Technical Notes:**
- Uses `calculateTmQ5()` for accurate Tm calculation
- Uses `calculateHairpinDG()` and `calculateHomodimerDG()` from equilibrium.js
- Custom weights now properly respected (fixed from initial implementation)

### 3.5 scar-preferences.js

Scores overhangs based on biological context.

**Key Exports:**
- `scoreScarSequence(overhang, context)` - Context-aware scoring
- `scoreScarSet(overhangs, context)` - Batch scoring
- `checkStopCodons(overhang, contextBefore)` - Stop codon detection
- `translateOverhang(overhang, contextBefore, contextAfter)` - Translation
- `SCAR_PREFERENCES` - Context preferences

**Supported Contexts:**
- `coding` - Prefer flexible amino acids, avoid stop codons
- `linker` - Prefer Gly/Ser-rich sequences
- `nonCoding` - Focus on assembly efficiency
- `startCodon` - Contains ATG
- `stopCodon` - After stop

**Technical Notes:**
- Checks all 3 reading frames for stop codons (fixed from initial implementation)
- Frame 2 requires upstream context (2bp before overhang)

### 3.6 failure-prediction.js

Predicts failure modes and calculates expected success rates.

**Key Exports:**
- `predictFailureModes(overhangs, enzyme, options)` - Comprehensive prediction
- `predictCrossLigation(overhangs, enzyme)` - Cross-ligation analysis
- `predictSelfLigation(overhangs)` - Self-ligation analysis
- `predictEfficiencyIssues(overhangs)` - Efficiency analysis
- `quickRiskAssessment(overhangs, enzyme)` - Quick summary
- `FAILURE_MODES` - Failure mode definitions

**Failure Modes:**
- `CROSS_LIGATION` - Wrong fragments ligating
- `SELF_LIGATION` - Fragment circularizing
- `INCOMPLETE_ASSEMBLY` - Missing fragments
- `LOW_EFFICIENCY` - Poor ligation kinetics
- `G_T_MISMATCH` - G:T wobble mis-ligation
- `INTERNAL_SITE` - Recognition site in fragment
- `HAIRPIN_FORMATION` - Strong secondary structure
- `MISPRIMING` - Off-target annealing

**Technical Notes:**
- Uses additive penalty model with diminishing returns (not multiplicative)
- Provides penalty breakdown for debugging
- Calculates colonies needed for 3 correct on average

### 3.7 fusion-site-optimizer.js

Main optimization engine with multiple algorithms.

**Key Exports:**
- `optimizeFusionSites(sequence, numFragments, options)` - Main entry point
- `optimizeGreedy(sequence, numFragments, enzyme, options)` - Fast algorithm
- `optimizeMonteCarlo(sequence, numFragments, enzyme, options)` - SA algorithm
- `optimizeBranchBound(sequence, numFragments, enzyme, options)` - Exact algorithm
- `optimizeHybrid(sequence, numFragments, enzyme, options)` - Combined approach
- `quickOptimize(sequence, numFragments, enzyme)` - Quick optimization
- `OPTIMIZER_DEFAULTS` - Default configuration

**Algorithm Selection (auto):**
- ≤5 junctions → Branch & Bound (optimal)
- 5-10 junctions → Hybrid (greedy + B&B + MC)
- >10 junctions → Monte Carlo (scalable)

**Technical Notes:**
- Monte Carlo initialization ensures unique overhangs (fixed)
- B&B uses improved upper bound estimation (fixed)
- All algorithms fall back to greedy if no valid solution found

---

## 4. Integration with Existing Tools

### From tmQ5.js:
```javascript
import { calculateTmQ5, calculate3primeTerminalDG } from '../tmQ5.js';
```
- NEB Q5 Tm calculator with 100% validation accuracy
- Nearest-neighbor thermodynamics with Owczarzy Mg²⁺ correction

### From equilibrium.js:
```javascript
import { calculateHairpinDG, calculateHomodimerDG } from '../equilibrium.js';
```
- Pythia equilibrium model for secondary structure
- Free energy calculation for hairpins and dimers

### From scoring.js:
```javascript
import {
  scoreTm, scoreGc, scoreHairpin, scoreHomodimer,
  scoreLength, scoreGcClamp, scoreHomopolymer,
  scoreTerminal3DG, scoreGQuadruplex,
} from '../scoring.js';
```
- Piecewise logistic scoring with biological thresholds
- Calibrated from PrimerScore2 meta-analysis

### From goldengate.js:
```javascript
import {
  getOverhangFidelityExperimental,
  calculateExperimentalFidelity,
  getEnzymeLigationData,
  findInternalSites,
} from './goldengate.js';
```
- Experimental ligation data from Pryor et al. 2020
- Full n×n matrix fidelity calculation

---

## 5. Natural Integration with Primer Design

### 5.1 Problem Statement (SOLVED)

The Fusion Site Optimizer was initially designed as a **separate mode** in the GoldenGateDesigner. This has been solved with two integration approaches:

**Previous Workflow (Disconnected):**
```
Path A: Multiple Parts → Design Primers → Get primers
Path B: Single Sequence → Toggle Fusion Mode → Optimizer → (Dead end - no primers!)
```

**NEW: Natural Integration (Implemented):**
```
Single Large Sequence → Click "Split" button → Auto-optimize → Fragments created → Design Primers
```

### 5.1.1 Natural Split-in-Place Feature

The **Split into Fragments** button appears on any part card with a sequence ≥500bp when using Golden Gate assembly. Clicking it:

1. Opens an inline panel to select number of fragments (2-6)
2. Runs `optimizeFusionSites()` to find optimal junction positions
3. Replaces the single part with optimized fragment parts automatically
4. User can then click "Design Primers" as normal

This eliminates the need to switch modes for the common use case.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  PART CARD (sequence ≥500bp)                                      [Split]  │
├─────────────────────────────────────────────────────────────────────────────┤
│  ┌─ Split into Optimized Fragments ──────────────────────────────────────┐  │
│  │  This will find optimal junction positions to split your 3,500 bp     │  │
│  │  sequence into fragments with high-fidelity overhangs for BsaI.       │  │
│  │                                                                        │  │
│  │  Number of fragments: [2] [3] [4] [5] [6]    (3 junctions, ~875bp ea)  │  │
│  │                                                                        │  │
│  │                                        [Cancel]  [Split & Optimize]   │  │
│  └────────────────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│  PARTS LIST (after split)                                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│  [Fragment_1] 891 bp   ├─ AGCT ─┤                                           │
│  [Fragment_2] 873 bp   ├─ TGCA ─┤                                           │
│  [Fragment_3] 869 bp   ├─ GCTA ─┤                                           │
│  [Fragment_4] 867 bp                                                        │
│                                                                             │
│                              [Design Primers]                               │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 5.2 Natural Integration Goal

The optimizer should be seamlessly integrated into the primer design workflow, so when a user wants to split a sequence into fragments, they get both optimal junction positions AND designed primers in one unified flow.

**Target Workflow (Unified):**
```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    UNIFIED PRIMER DESIGN WORKFLOW                           │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   INPUT: Single Sequence + Number of Fragments                              │
│                     ↓                                                       │
│   STEP 1: Fusion Site Optimizer                                             │
│           • Scan for candidate junction positions                           │
│           • Run optimization algorithm (auto-selected)                      │
│           • Calculate set fidelity and failure predictions                  │
│                     ↓                                                       │
│   STEP 2: Junction Review (Optional)                                        │
│           • Display junction candidates with quality scores                 │
│           • Allow user to accept/refine positions                           │
│           • Show fidelity gauge and risk assessment                         │
│                     ↓                                                       │
│   STEP 3: Convert Junctions → Parts                                         │
│           • Split sequence at optimized junction positions                  │
│           • Create parts array with fragment sequences                      │
│           • Preserve overhang information for each junction                 │
│                     ↓                                                       │
│   STEP 4: Design Primers                                                    │
│           • Call designOptimizedGoldenGateAssemblyForUI()                   │
│           • Use pre-computed overhangs from optimizer                       │
│           • Generate forward/reverse primers for each fragment              │
│                     ↓                                                       │
│   STEP 5: Unified Results                                                   │
│           • Show junction quality + primer sequences together               │
│           • Display assembly fidelity and success prediction                │
│           • Export primers with full context                                │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 5.3 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         GoldenGateDesigner.jsx                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   ┌─────────────────────────────────────────────────────────────────────┐   │
│   │                    Fusion Mode (Single Sequence)                     │   │
│   ├─────────────────────────────────────────────────────────────────────┤   │
│   │                                                                     │   │
│   │   Input: fusionSequence, fusionNumFragments, enzyme                 │   │
│   │                          ↓                                          │   │
│   │   handleFusionOptimize() ─────────────────────────────────────────► │   │
│   │   │                                                                 │   │
│   │   │  ┌──────────────────────────────────────────────────────────┐  │   │
│   │   │  │  optimizeFusionSites(sequence, numFragments, options)     │  │   │
│   │   │  │  └─► Returns: { junctions, overhangs, score, ... }       │  │   │
│   │   │  └──────────────────────────────────────────────────────────┘  │   │
│   │   │                          ↓                                      │   │
│   │   │  setFusionResult(result)                                        │   │
│   │   │  setFusionCandidates(candidates)                                │   │
│   │                          ↓                                          │   │
│   │   FusionSiteOptimizerPanel                                          │   │
│   │   │  • Shows junction candidates                                    │   │
│   │   │  • Displays fidelity gauge                                      │   │
│   │   │  • Failure prediction                                           │   │
│   │                          ↓                                          │   │
│   │   [Accept & Design Primers] Button ──────────────────────────────►  │   │
│   │                          ↓                                          │   │
│   │   convertJunctionsToParts(fusionResult) ─────────────────────────►  │   │
│   │   │                                                                 │   │
│   │   │  ┌──────────────────────────────────────────────────────────┐  │   │
│   │   │  │  Creates parts array from optimized junction positions:   │  │   │
│   │   │  │                                                           │  │   │
│   │   │  │  parts = [                                                │  │   │
│   │   │  │    { id: 'Fragment_1', seq: seq[0..j1], leftOH, rightOH } │  │   │
│   │   │  │    { id: 'Fragment_2', seq: seq[j1..j2], leftOH, rightOH }│  │   │
│   │   │  │    ...                                                    │  │   │
│   │   │  │  ]                                                        │  │   │
│   │   │  └──────────────────────────────────────────────────────────┘  │   │
│   │                          ↓                                          │   │
│   │   designOptimizedGoldenGateAssemblyForUI(parts, { enzyme, ... })    │   │
│   │                          ↓                                          │   │
│   │   setResult(assemblyResult)                                         │   │
│   │                          ↓                                          │   │
│   │   PrimerResults Component                                           │   │
│   │   │  • Primer sequences for each fragment                           │   │
│   │   │  • Junction quality indicators                                  │   │
│   │   │  • Assembly fidelity summary                                    │   │
│   │                                                                     │   │
│   └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 5.4 Key Functions for Integration

#### 5.4.1 handleFusionOptimize (Updated)

Replace mock implementation with real optimizer:

```javascript
const handleFusionOptimize = useCallback(async (options) => {
  if (!fusionSequence || fusionSequence.length < 500) {
    setError('Sequence must be at least 500bp');
    return;
  }

  setIsOptimizing(true);
  setError(null);

  try {
    // Import and call real optimizer
    const { optimizeFusionSites } = await import('../lib/repp/fusion-site-optimizer');
    const { scanForFusionSites } = await import('../lib/repp/fusion-site-scanner');

    // Get all candidates for display
    const candidates = scanForFusionSites(fusionSequence, {
      enzyme,
      minDistanceFromEnds: options.constraints?.minDistanceFromEnds || 50,
    });

    // Run optimization
    const result = optimizeFusionSites(fusionSequence, fusionNumFragments, {
      enzyme,
      algorithm: options.algorithm || 'auto',
      minFragmentSize: options.constraints?.minFragmentSize || 200,
      maxFragmentSize: options.constraints?.maxFragmentSize || 3000,
      minDistanceFromEnds: options.constraints?.minDistanceFromEnds || 50,
      codingFrame: options.bioContext?.codingFrame,
      scarContext: options.bioContext?.scarPreferences || 'nonCoding',
    });

    // Format candidates with composite scores for UI
    const formattedCandidates = candidates.map(c => ({
      ...c,
      composite: c.fidelity * 100 || 70,
    }));

    setFusionCandidates(formattedCandidates);
    setFusionResult(result);

  } catch (err) {
    setError(`Optimization failed: ${err.message}`);
    console.error('Fusion optimization error:', err);
  } finally {
    setIsOptimizing(false);
  }
}, [fusionSequence, fusionNumFragments, enzyme]);
```

#### 5.4.2 convertJunctionsToParts (New)

Convert optimizer results to parts array for primer design:

```javascript
/**
 * Convert optimized junction positions to parts array for primer design
 */
function convertJunctionsToParts(sequence, fusionResult) {
  if (!fusionResult?.junctions || fusionResult.junctions.length === 0) {
    return [];
  }

  const junctions = fusionResult.junctions.sort((a, b) => a.position - b.position);
  const parts = [];

  // Create parts from junction positions
  let startPos = 0;

  for (let i = 0; i <= junctions.length; i++) {
    const isFirst = i === 0;
    const isLast = i === junctions.length;

    // Calculate fragment boundaries
    const fragStart = isFirst ? 0 : junctions[i - 1].position + 4; // +4 for overhang
    const fragEnd = isLast ? sequence.length : junctions[i].position + 4;

    // Extract fragment sequence
    const fragSeq = sequence.slice(fragStart, fragEnd);

    // Get overhangs
    const leftOverhang = isFirst ? null : junctions[i - 1].overhang;
    const rightOverhang = isLast ? null : junctions[i].overhang;

    parts.push({
      id: `Fragment_${i + 1}`,
      seq: fragSeq,
      type: 'fragment',
      leftOverhang,
      rightOverhang,
      // Store junction info for reference
      _junctionInfo: {
        startPos: fragStart,
        endPos: fragEnd,
        leftJunction: isFirst ? null : junctions[i - 1],
        rightJunction: isLast ? null : junctions[i],
      },
    });
  }

  return parts;
}
```

#### 5.4.3 handleAcceptAndDesignPrimers (New)

Complete the flow from junctions to primers:

```javascript
const handleAcceptAndDesignPrimers = useCallback(async () => {
  if (!fusionResult?.success) {
    setError('No valid optimization result to use');
    return;
  }

  setIsOptimizing(true);

  try {
    // Convert junctions to parts
    const parts = convertJunctionsToParts(fusionSequence, fusionResult);

    if (parts.length < 2) {
      setError('Need at least 2 fragments for assembly');
      return;
    }

    // Design primers using existing Golden Gate primer designer
    const { designOptimizedGoldenGateAssemblyForUI } = await import(
      '../lib/repp/goldengate-primer-optimizer'
    );

    const assemblyResult = designOptimizedGoldenGateAssemblyForUI(parts, {
      enzyme,
      circular: true,
      // Pass pre-computed overhangs to use
      precomputedOverhangs: fusionResult.overhangs,
    });

    // Merge fusion result metadata with assembly result
    const mergedResult = {
      ...assemblyResult,
      fusionOptimization: {
        algorithm: fusionResult.algorithm,
        setFidelity: fusionResult.score?.fidelity,
        setEfficiency: fusionResult.score?.efficiency,
        failurePrediction: fusionResult.failurePrediction,
        junctions: fusionResult.junctions,
      },
    };

    // Exit fusion mode and show results
    setFusionMode(false);
    setParts(parts);
    setResult(mergedResult);

  } catch (err) {
    setError(`Primer design failed: ${err.message}`);
    console.error('Primer design error:', err);
  } finally {
    setIsOptimizing(false);
  }
}, [fusionSequence, fusionResult, enzyme]);
```

### 5.5 UI Flow Enhancements

#### Add "Accept & Design Primers" Button

After optimization results are shown in FusionSiteOptimizerPanel, add a prominent button:

```jsx
{fusionResult?.success && (
  <div className="fusion-action-bar">
    <div className="fusion-summary">
      <span className="fidelity-badge">
        {(fusionResult.score?.fidelity * 100).toFixed(1)}% Fidelity
      </span>
      <span className="fragment-count">
        {fusionResult.numFragments} Fragments
      </span>
      <span className="risk-badge" data-risk={fusionResult.failurePrediction?.overallRisk}>
        {fusionResult.failurePrediction?.overallRisk || 'low'} risk
      </span>
    </div>

    <button
      className="accept-design-btn primary"
      onClick={handleAcceptAndDesignPrimers}
      disabled={isOptimizing}
    >
      <svg>...</svg>
      Accept & Design Primers
    </button>
  </div>
)}
```

### 5.6 Data Flow Summary

```
User Input                    Backend Processing              UI Updates
──────────────────────────────────────────────────────────────────────────────
fusionSequence        ──►     optimizeFusionSites()    ──►   fusionResult
fusionNumFragments                                            fusionCandidates
enzyme                                                         │
options                                                        ▼
                                                         [Show Candidates]
                                                         [Show Fidelity]
                                                         [Show Risk]
                                                               │
                             ◄── [Accept & Design] ◄──────────┘
                                       │
                      convertJunctionsToParts()
                                       │
                                       ▼
                      designOptimizedGoldenGateAssemblyForUI()
                                       │
                                       ▼
                                                         [Show Primers]
                                                         [Show Assembly]
                                                         [Export Options]
```

### 5.7 Benefits of Natural Integration

1. **Single Workflow**: Users don't need to understand two separate modes
2. **Optimized Results**: Every fragment assembly gets optimal junction positions
3. **Unified Quality Metrics**: See junction fidelity AND primer quality together
4. **Fewer Clicks**: Go from sequence to primers in one flow
5. **Better Success Rates**: Optimized overhangs lead to higher assembly success

---

## 7. Primer Boundary Optimizer (Edge Cases)

This section describes the **Primer Boundary Optimizer** - a separate tool for optimizing primer design when users have pre-defined fragment boundaries.

### 7.1 Problem Statement

When users provide multiple fragments with fixed boundaries, the primer binding regions at each junction may have poor quality:

```
User's fragments:
[Fragment A: ...NNNNNAAAAA] | [Fragment B: GCGCGNNNNN...]
                    ↑                    ↑
         Reverse primer of A      Forward primer of B
         binds here (poly-A!)     binds here (good)

Problem: The reverse primer for Fragment A has to bind to "AAAAA" which:
- Has very low Tm
- May form weak hairpins
- Has no GC clamp
```

### 7.2 Solution: Boundary Shifting

The Primer Boundary Optimizer shifts the boundary left or right to find better primer binding regions:

```
BEFORE (user boundary):
Fragment A: ...NNNNNAAAAA | GCGCGNNNNN...  Fragment B
                    └──┘ poor primer region

AFTER (shifted left 5bp):
Fragment A: ...NNNNN | AAAAAGCGCGNNNNN...  Fragment B
                   └──┘ now Fragment B's forward primer
                        starts from better region

Or AFTER (shifted right 5bp):
Fragment A: ...NNNNNAAAAAGCGCG | NNNNN...  Fragment B
                    └────────┘ Fragment A's reverse primer
                              now includes the GC-rich region
```

### 7.3 Implementation

**File**: `src/lib/repp/primer-boundary-optimizer.js`

#### Key Functions

```javascript
import {
  optimizeJunctionBoundary,
  optimizeAssemblyBoundaries,
  assessBoundaryOptimizationPotential,
} from './primer-boundary-optimizer.js';

// Optimize a single junction between two fragments
const result = optimizeJunctionBoundary(leftSeq, rightSeq, {
  maxShift: 50,           // Search ±50bp from original boundary
  minFragmentSize: 100,   // Don't shrink fragments below this
  targetTm: 60,           // Target primer Tm
});

// result = {
//   shift: -5,           // Shifted 5bp left
//   direction: 'left',
//   improvement: 0.15,   // 15% score improvement
//   beforePrimers: { left: {...}, right: {...} },
//   afterPrimers: { left: {...}, right: {...} },
//   reason: "Shift boundary 5bp left: primer quality improved from 65.2 to 80.1"
// }

// Optimize all boundaries in an assembly
const fullResult = optimizeAssemblyBoundaries(fragments, options);
// Returns optimized fragments with adjusted sequences
```

#### Algorithm

1. **Extract primer binding regions** at the junction (15-30bp each side)
2. **Score current primers** for Tm, GC, hairpins, homodimers, etc.
3. **Search left** (boundary into left fragment, ±1bp to ±maxShift)
   - For each position, score new primer binding regions
   - Check that resulting overhang is valid (no palindromes, no homopolymers)
4. **Search right** (boundary into right fragment)
   - Same scoring process
5. **Select best position** considering:
   - Primer quality improvement
   - Shift penalty (prefer smaller shifts)
   - Overhang validity
6. **Return optimized boundary** with before/after comparison

### 7.4 Configuration

```javascript
const BOUNDARY_OPTIMIZER_DEFAULTS = {
  maxShift: 50,              // Max bp to shift from original boundary
  minShiftThreshold: 3,      // Ignore shifts smaller than this
  minHomologyLength: 15,     // Minimum primer binding region
  maxHomologyLength: 30,     // Maximum primer binding region
  targetTm: 60,              // Target Tm for primers
  tmTolerance: 5,            // ±5°C acceptable
  minFragmentSize: 100,      // Don't shrink fragments below this

  weights: {
    primerQuality: 0.40,     // Weight for overall primer score
    tmMatch: 0.25,           // Weight for Tm accuracy
    balancedShift: 0.15,     // Penalty for large shifts
    overhangFidelity: 0.20,  // Weight for overhang quality
  },

  avoidPalindromes: true,    // Skip palindromic overhangs
  avoidHomopolymers: true,   // Skip homopolymer overhangs
};
```

### 7.5 UI Integration

The Primer Boundary Optimizer can be integrated into the GoldenGateDesigner when:
1. User has multiple parts defined
2. After primer design, some primers show poor quality
3. User clicks "Optimize Boundaries" to improve problem junctions

```jsx
// In GoldenGateDesigner.jsx
const handleOptimizeBoundaries = async () => {
  const assessment = assessBoundaryOptimizationPotential(parts);

  if (assessment.needsOptimization) {
    const optimized = await optimizeAssemblyBoundaries(parts, {
      enzyme,
      maxShift: 50,
    });

    // Update parts with optimized sequences
    setParts(optimized.optimizedFragments);

    // Re-design primers with new boundaries
    await handleDesignPrimers();
  }
};
```

### 7.6 Example Use Case

**Scenario**: User has 3 gene fragments from a paper that must be assembled in a specific order.

```javascript
const fragments = [
  { id: 'promoter', seq: 'ATGC...AAAA' },   // Ends with poly-A
  { id: 'coding', seq: 'GCGC...TTTT' },     // Starts with GC-rich
  { id: 'terminator', seq: 'GGGG...ATCG' }, // Ends with poly-G
];

// Check if boundaries need optimization
const assessment = assessBoundaryOptimizationPotential(fragments);
// assessment.issues = [
//   { junction: 0, side: 'left', quality: 'poor', issues: ['Low Tm: 45°C'] },
//   { junction: 1, side: 'left', quality: 'acceptable', issues: ['Weak GC clamp'] },
// ]

// Optimize boundaries
const optimized = optimizeAssemblyBoundaries(fragments);
// optimized.boundaries[0].shift = -8  (shifted left 8bp)
// optimized.boundaries[1].shift = 5   (shifted right 5bp)

// Result: promoter is 8bp shorter, coding is 13bp longer (gained 8 from promoter, lost 5 to terminator)
// But now all primers have good quality!
```

---

## 8. UI Implementation Strategy

### 6.1 Proposed Component Structure

```
src/components/golden-gate/
├── FusionSiteDesigner.jsx        # Main container component
├── FusionSiteDesigner.css        # Styles
├── components/
│   ├── SequenceInput.jsx         # Sequence input with validation
│   ├── ParameterPanel.jsx        # Configuration options
│   ├── JunctionVisualizer.jsx    # Visual sequence map
│   ├── ResultsPanel.jsx          # Optimization results
│   ├── JunctionDetails.jsx       # Individual junction info
│   ├── FailurePrediction.jsx     # Risk assessment display
│   ├── PrimerTable.jsx           # Primer sequences table
│   └── ExportPanel.jsx           # Download options
└── hooks/
    ├── useFusionOptimizer.js     # Optimization state management
    └── useJunctionSelection.js   # Interactive junction editing
```

### 6.2 Main Component: FusionSiteDesigner.jsx

```jsx
import React, { useState, useCallback } from 'react';
import { optimizeFusionSites } from '../../lib/repp/fusion-site-optimizer';
import { predictFailureModes } from '../../lib/repp/failure-prediction';

export function FusionSiteDesigner() {
  const [sequence, setSequence] = useState('');
  const [numFragments, setNumFragments] = useState(4);
  const [enzyme, setEnzyme] = useState('BsaI');
  const [options, setOptions] = useState({
    algorithm: 'auto',
    minFragmentSize: 200,
    maxFragmentSize: 5000,
    minDistanceFromEnds: 50,
  });
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);

  const handleOptimize = useCallback(async () => {
    setLoading(true);
    try {
      const optimizationResult = await optimizeFusionSites(
        sequence,
        numFragments,
        { enzyme, ...options }
      );
      setResult(optimizationResult);
    } catch (error) {
      console.error('Optimization failed:', error);
    } finally {
      setLoading(false);
    }
  }, [sequence, numFragments, enzyme, options]);

  return (
    <div className="fusion-site-designer">
      <div className="designer-header">
        <h2>Fusion Site Optimizer</h2>
        <p>Find optimal junction positions for Golden Gate assembly</p>
      </div>

      <div className="designer-layout">
        {/* Left Panel: Input */}
        <div className="input-panel">
          <SequenceInput
            value={sequence}
            onChange={setSequence}
          />
          <ParameterPanel
            numFragments={numFragments}
            onNumFragmentsChange={setNumFragments}
            enzyme={enzyme}
            onEnzymeChange={setEnzyme}
            options={options}
            onOptionsChange={setOptions}
          />
          <button
            onClick={handleOptimize}
            disabled={loading || !sequence}
            className="optimize-button"
          >
            {loading ? 'Optimizing...' : 'Optimize Junctions'}
          </button>
        </div>

        {/* Center Panel: Visualization */}
        <div className="visualization-panel">
          <JunctionVisualizer
            sequence={sequence}
            junctions={result?.junctions}
            enzyme={enzyme}
          />
        </div>

        {/* Right Panel: Results */}
        <div className="results-panel">
          {result && (
            <>
              <ResultsPanel result={result} />
              <FailurePrediction prediction={result.failurePrediction} />
              <PrimerTable primers={result.detailedJunctions} />
              <ExportPanel result={result} />
            </>
          )}
        </div>
      </div>
    </div>
  );
}
```

### 6.3 Junction Visualizer Component

```jsx
export function JunctionVisualizer({ sequence, junctions, enzyme }) {
  const seqLength = sequence?.length || 0;
  const scale = seqLength > 0 ? 800 / seqLength : 1;

  return (
    <div className="junction-visualizer">
      <svg width="100%" height="200" viewBox="0 0 850 200">
        {/* Sequence backbone */}
        <rect x="25" y="90" width="800" height="20" fill="#e0e0e0" rx="5" />

        {/* Scale markers */}
        {[0, 0.25, 0.5, 0.75, 1].map((frac, i) => (
          <g key={i}>
            <line
              x1={25 + frac * 800}
              y1="80"
              x2={25 + frac * 800}
              y2="120"
              stroke="#888"
              strokeWidth="1"
            />
            <text x={25 + frac * 800} y="75" textAnchor="middle" fontSize="10">
              {Math.round(frac * seqLength)}
            </text>
          </g>
        ))}

        {/* Junction markers */}
        {junctions?.map((junction, idx) => {
          const x = 25 + junction.position * scale;
          const quality = junction.score?.composite || 70;
          const color = quality >= 85 ? '#4caf50' :
                        quality >= 70 ? '#ff9800' :
                        quality >= 55 ? '#ff5722' : '#f44336';

          return (
            <g key={idx} className="junction-marker">
              {/* Junction line */}
              <line
                x1={x}
                y1="50"
                x2={x}
                y2="150"
                stroke={color}
                strokeWidth="3"
              />
              {/* Overhang label */}
              <text x={x} y="40" textAnchor="middle" fontSize="10" fill={color}>
                {junction.overhang}
              </text>
              {/* Fragment number */}
              <text x={x} y="170" textAnchor="middle" fontSize="10">
                F{idx + 1}|F{idx + 2}
              </text>
            </g>
          );
        })}

        {/* Fragment size annotations */}
        {junctions?.length > 0 && (
          <g className="fragment-sizes">
            {[0, ...junctions.map(j => j.position), seqLength].slice(0, -1).map((start, i, arr) => {
              const end = arr[i + 1] || seqLength;
              const size = end - start;
              const midX = 25 + ((start + end) / 2) * scale;
              return (
                <text key={i} x={midX} y="135" textAnchor="middle" fontSize="9" fill="#666">
                  {size}bp
                </text>
              );
            })}
          </g>
        )}
      </svg>
    </div>
  );
}
```

### 6.4 Parameter Panel Component

```jsx
export function ParameterPanel({
  numFragments, onNumFragmentsChange,
  enzyme, onEnzymeChange,
  options, onOptionsChange,
}) {
  const enzymes = ['BsaI', 'BsmBI', 'BbsI', 'Esp3I', 'SapI'];
  const algorithms = [
    { value: 'auto', label: 'Auto (recommended)' },
    { value: 'greedy', label: 'Greedy (fastest)' },
    { value: 'branchBound', label: 'Branch & Bound (optimal, ≤6 junctions)' },
    { value: 'monteCarlo', label: 'Monte Carlo (large assemblies)' },
    { value: 'hybrid', label: 'Hybrid (balanced)' },
  ];

  return (
    <div className="parameter-panel">
      <h3>Parameters</h3>

      <div className="param-group">
        <label>Number of Fragments</label>
        <input
          type="number"
          min="2"
          max="20"
          value={numFragments}
          onChange={(e) => onNumFragmentsChange(parseInt(e.target.value))}
        />
        <span className="hint">{numFragments - 1} junction(s)</span>
      </div>

      <div className="param-group">
        <label>Enzyme</label>
        <select value={enzyme} onChange={(e) => onEnzymeChange(e.target.value)}>
          {enzymes.map(enz => (
            <option key={enz} value={enz}>{enz}</option>
          ))}
        </select>
      </div>

      <div className="param-group">
        <label>Algorithm</label>
        <select
          value={options.algorithm}
          onChange={(e) => onOptionsChange({ ...options, algorithm: e.target.value })}
        >
          {algorithms.map(alg => (
            <option key={alg.value} value={alg.value}>{alg.label}</option>
          ))}
        </select>
      </div>

      <details className="advanced-options">
        <summary>Advanced Options</summary>

        <div className="param-group">
          <label>Min Fragment Size (bp)</label>
          <input
            type="number"
            value={options.minFragmentSize}
            onChange={(e) => onOptionsChange({
              ...options,
              minFragmentSize: parseInt(e.target.value)
            })}
          />
        </div>

        <div className="param-group">
          <label>Max Fragment Size (bp)</label>
          <input
            type="number"
            value={options.maxFragmentSize}
            onChange={(e) => onOptionsChange({
              ...options,
              maxFragmentSize: parseInt(e.target.value)
            })}
          />
        </div>

        <div className="param-group">
          <label>Min Distance from Ends (bp)</label>
          <input
            type="number"
            value={options.minDistanceFromEnds}
            onChange={(e) => onOptionsChange({
              ...options,
              minDistanceFromEnds: parseInt(e.target.value)
            })}
          />
        </div>
      </details>
    </div>
  );
}
```

### 6.5 Results Panel Component

```jsx
export function ResultsPanel({ result }) {
  if (!result) return null;

  const { success, algorithm, score, junctions, fragmentSizes, summary } = result;

  return (
    <div className={`results-panel ${success ? 'success' : 'failure'}`}>
      <div className="results-header">
        <h3>{success ? '✓ Optimization Complete' : '⚠ Optimization Issues'}</h3>
        <span className="algorithm-badge">{algorithm}</span>
      </div>

      <div className="score-summary">
        <div className="score-item">
          <span className="label">Set Fidelity</span>
          <span className="value">{(summary?.fidelity * 100 || 0).toFixed(1)}%</span>
        </div>
        <div className="score-item">
          <span className="label">Efficiency</span>
          <span className="value">{(summary?.efficiency * 100 || 0).toFixed(1)}%</span>
        </div>
        <div className="score-item">
          <span className="label">Expected Success</span>
          <span className="value">{summary?.expectedSuccessRate || 'N/A'}</span>
        </div>
      </div>

      <div className="junction-list">
        <h4>Junctions ({junctions?.length || 0})</h4>
        <table>
          <thead>
            <tr>
              <th>Position</th>
              <th>Overhang</th>
              <th>Quality</th>
            </tr>
          </thead>
          <tbody>
            {junctions?.map((j, i) => (
              <tr key={i}>
                <td>{j.position}</td>
                <td><code>{j.overhang}</code></td>
                <td>
                  <span className={`quality-badge quality-${j.quality || 'good'}`}>
                    {(j.score?.composite || j.composite || 70).toFixed(0)}
                  </span>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <div className="fragment-summary">
        <h4>Fragments ({fragmentSizes?.length || 0})</h4>
        <div className="fragment-sizes">
          {fragmentSizes?.map((size, i) => (
            <span key={i} className="fragment-size">
              F{i + 1}: {size}bp
            </span>
          ))}
        </div>
      </div>
    </div>
  );
}
```

### 6.6 Failure Prediction Component

```jsx
export function FailurePrediction({ prediction }) {
  if (!prediction) return null;

  const { overallRisk, expectedSuccessRate, recommendations, riskSummary } = prediction;

  const riskColors = {
    minimal: '#4caf50',
    low: '#8bc34a',
    medium: '#ff9800',
    high: '#f44336',
    critical: '#9c27b0',
  };

  return (
    <div className="failure-prediction">
      <h4>Risk Assessment</h4>

      <div className="risk-overview" style={{ borderLeftColor: riskColors[overallRisk] }}>
        <div className="risk-level">
          <span className="label">Overall Risk</span>
          <span className={`value risk-${overallRisk}`}>{overallRisk}</span>
        </div>
        <div className="success-rate">
          <span className="label">Expected Success Rate</span>
          <span className="value">{expectedSuccessRate?.ratePercent}</span>
        </div>
        <div className="colonies">
          <span className="label">Colonies to Screen</span>
          <span className="value">{expectedSuccessRate?.coloniesNeeded}</span>
        </div>
      </div>

      {riskSummary && (
        <div className="risk-breakdown">
          {riskSummary.critical > 0 && (
            <span className="risk-count critical">{riskSummary.critical} critical</span>
          )}
          {riskSummary.high > 0 && (
            <span className="risk-count high">{riskSummary.high} high</span>
          )}
          {riskSummary.medium > 0 && (
            <span className="risk-count medium">{riskSummary.medium} medium</span>
          )}
        </div>
      )}

      {recommendations?.length > 0 && (
        <div className="recommendations">
          <h5>Recommendations</h5>
          <ul>
            {recommendations.map((rec, i) => (
              <li key={i}>
                <strong>{rec.issue}:</strong> {rec.recommendation}
              </li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
}
```

### 6.7 UI Implementation Phases

#### Phase 1: Core UI (1-2 weeks)
1. Create `FusionSiteDesigner.jsx` main component
2. Implement `SequenceInput` with validation
3. Implement `ParameterPanel` with all options
4. Basic `ResultsPanel` showing junctions and scores
5. Add to application routing

#### Phase 2: Visualization (1 week)
1. Implement `JunctionVisualizer` SVG component
2. Add interactive junction selection
3. Add zoom/pan for long sequences
4. Color coding based on quality scores

#### Phase 3: Analysis Display (1 week)
1. Implement `FailurePrediction` component
2. Add `JunctionDetails` expandable panels
3. Show primer sequences with `PrimerTable`
4. Add warnings and recommendations display

#### Phase 4: Export & Polish (1 week)
1. Implement `ExportPanel` with formats:
   - CSV (primers table)
   - JSON (full result)
   - IDT order format
   - Twist order format
2. Add print-friendly view
3. Add keyboard shortcuts
4. Performance optimization for large sequences

---

## 7. API Reference

### Main Entry Point

```javascript
import { optimizeFusionSites } from './lib/repp/fusion-site-optimizer';

const result = await optimizeFusionSites(sequence, numFragments, {
  // Basic
  enzyme: 'BsaI',                    // Enzyme to use
  algorithm: 'auto',                 // 'auto' | 'greedy' | 'branchBound' | 'monteCarlo' | 'hybrid'

  // Constraints
  minFragmentSize: 200,              // Minimum fragment size in bp
  maxFragmentSize: 5000,             // Maximum fragment size in bp
  minDistanceFromEnds: 50,           // Minimum distance from sequence ends
  searchRadius: 50,                  // Search radius around target positions
  forbiddenRegions: [],              // [{start, end}, ...] - excluded regions

  // Biological context
  codingFrame: null,                 // 0, 1, or 2 for reading frame
  proteinDomains: [],                // [{start, end, name}, ...]
  scarContext: 'nonCoding',          // 'coding' | 'linker' | 'nonCoding'

  // Monte Carlo options
  mcIterations: 2000,                // Number of iterations
  mcTemperature: 1.0,                // Initial temperature
  mcCoolingRate: 0.995,              // Cooling rate

  // Branch & Bound options
  maxBranchDepth: 8,                 // Maximum search depth
  pruningThreshold: 0.7,             // Score threshold for pruning

  // Output
  verbose: false,                    // Log progress
});
```

### Result Structure

```javascript
{
  success: true,                     // Whether optimization succeeded
  algorithm: 'branchBound',          // Algorithm used
  enzyme: 'BsaI',                    // Enzyme used
  sequenceLength: 5000,              // Input sequence length
  numFragments: 5,                   // Number of fragments
  numJunctions: 4,                   // Number of junctions

  // Core results
  junctions: [                       // Junction details
    {
      position: 1200,
      overhang: 'AGCT',
      targetIndex: 0,
      idealPosition: 1250,
      deviation: -50,
      score: { composite: 85, ... },
    },
    // ...
  ],
  overhangs: ['AGCT', 'TGCA', ...], // Overhang sequences
  score: {                           // Set-level scores
    composite: 0.82,
    fidelity: 0.95,
    efficiency: 0.88,
    primerQuality: 0.75,
    positionQuality: 0.90,
  },

  // Detailed analysis
  detailedJunctions: [...],          // Full scoring for each junction
  failurePrediction: {               // Risk assessment
    overallRisk: 'low',
    predictions: [...],
    expectedSuccessRate: { rate: 0.85, tier: 'good' },
    recommendations: [...],
  },

  // Fragment info
  fragmentSizes: [1200, 1100, ...],  // Fragment sizes in bp
  minFragmentSize: 1100,
  maxFragmentSize: 1300,
  avgFragmentSize: 1200,

  // Warnings
  internalSites: null,               // Internal restriction sites if any
  summary: {
    fidelity: 0.95,
    efficiency: 0.88,
    expectedSuccessRate: '85%',
    riskLevel: 'low',
    recommendations: [...],
  },
}
```

---

## 8. Testing Coverage

### Test File: `fusion-site-optimizer.test.js`

**63 tests covering:**

1. **Overhang Efficiency Module** (8 tests)
   - isPalindrome detection
   - isHomopolymer detection
   - calculateEfficiency scoring
   - calculateSetEfficiency batch processing

2. **Site Creation Check Module** (6 tests)
   - checkSiteCreation single junction
   - checkMultipleSiteCreation batch
   - findSafeJunctionNear alternative finding

3. **Fusion Site Scanner Module** (8 tests)
   - scanForFusionSites candidate generation
   - generateTargetPositions calculation
   - filterByDistance spacing validation
   - assessFeasibility checks

4. **Fusion Site Scorer Module** (8 tests)
   - scoreFusionSiteComposite full scoring
   - quickScoreFusionSite fast path
   - Invalid position handling

5. **Scar Preferences Module** (6 tests)
   - scoreScarSequence context scoring
   - checkStopCodons frame detection
   - scoreScarSet batch processing

6. **Failure Prediction Module** (8 tests)
   - predictFailureModes comprehensive
   - predictCrossLigation analysis
   - predictSelfLigation palindrome detection
   - quickRiskAssessment summary

7. **Optimizer Module** (12 tests)
   - optimizeGreedy algorithm
   - optimizeMonteCarlo algorithm
   - optimizeBranchBound algorithm
   - optimizeHybrid combined approach
   - optimizeFusionSites main entry

8. **Integration Tests** (4 tests)
   - End-to-end optimization
   - Algorithm comparison
   - Edge cases

9. **Enzyme Support Tests** (3 tests)
   - BsaI, BsmBI, BbsI support

---

## 9. Future Enhancements

### Short-term (1-2 months)
- [ ] Interactive junction adjustment in UI
- [ ] Undo/redo for optimization
- [ ] Save/load optimization sessions
- [ ] Batch sequence processing

### Medium-term (3-6 months)
- [ ] Vector integration (pUC19, pET, etc.)
- [ ] Protocol generation (NEB format)
- [ ] Order file generation (IDT, Twist)
- [ ] Primer dimer analysis between fragments

### Long-term (6+ months)
- [ ] Machine learning for efficiency prediction
- [ ] Integration with external assembly planning tools
- [ ] Collaborative design features
- [ ] Assembly simulation and validation

---

## References

1. Potapov V, et al. (2018) "Comprehensive Profiling of Four Base Overhang Ligation Fidelity by T4 DNA Ligase" ACS Synth Biol. DOI: 10.1021/acssynbio.8b00333

2. Pryor JM, et al. (2020) "Enabling one-pot Golden Gate assemblies of unprecedented complexity" PLOS ONE. DOI: 10.1371/journal.pone.0238592

3. NEB Technical Guide: "Tips for Optimizing Golden Gate Assembly Reactions"

4. NEB NEBridge Tools: https://ligasefidelity.neb.com

---

*Document Version: 3.1*
*Last Updated: 2025-12-16*
*Status: Backend Complete, UI Complete, Natural Integration Complete*
*Author: Claude Code Assistant*

### Change Log

| Version | Date | Changes |
|---------|------|---------|
| 3.1 | 2025-12-16 | Added Natural Split-in-Place feature for seamless fragment assembly |
| 3.0 | 2025-12-16 | Added Section 5: Natural Integration with Primer Design workflow |
| 2.0 | 2025-12-16 | Complete backend implementation with all algorithms |
| 1.0 | 2025-12-15 | Initial strategy document |
