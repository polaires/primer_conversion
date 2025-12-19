# Gibson & Golden Gate Assembly Evaluation
## Building a State-of-the-Art NEB Replacement Tool

---

## Executive Summary

This codebase contains a **production-grade implementation** of DNA assembly and primer design algorithms with ~28,000 lines of specialized library code. It has significant potential to replace NEB tools (NEBChanger, HiFi Assembly Tool, Golden Gate Tool) by leveraging:

1. **Experimental Ligation Fidelity Data** - 473MB of Pryor et al. (2020) enzyme-specific data
2. **Advanced Thermodynamic Models** - Pythia equilibrium + DNA24 nearest-neighbor parameters
3. **Monte Carlo Overhang Optimization** - Direct port of NEB's algorithm
4. **Calibrated Scoring System** - Empirically validated on 829-pair Döring dataset (F1=81.9%)

---

## 1. Current Golden Gate Assembly Capabilities

### 1.1 Strengths ✅

| Feature | Implementation | NEB Tool Comparison |
|---------|---------------|---------------------|
| **Type IIS Enzymes** | BsaI-HFv2, BbsI-HF, BsmBI-v2, Esp3I, SapI | ✅ Matches NEB |
| **Ligation Fidelity Data** | Pryor et al. 2020 experimental matrices | ✅ **Superior** - NEB tools don't expose raw data |
| **Overhang Optimization** | Monte Carlo simulated annealing | ✅ Port of NEB's ggtools_optimize.pl |
| **High-Fidelity Sets** | Potapov et al. 2018 validated sets (2-30 parts) | ✅ Matches NEB |
| **Internal Site Detection** | Full domestication suggestions with silent mutations | ✅ Better than NEBChanger |
| **Multi-Enzyme Comparison** | Cross-enzyme fidelity analysis | ❌ NEB doesn't offer this |
| **Per-Junction Fidelity** | Detailed cross-ligation risk analysis | ❌ NEB shows only overall |

### 1.2 Key Algorithms

```
Golden Gate Assembly Flow:
┌─────────────────────────────────────────────────────────────┐
│  1. Input Parts                                              │
│     ↓                                                        │
│  2. Check Internal Sites → findInternalSites()              │
│     ↓                                                        │
│  3. Suggest Domestication → suggestDomestication()          │
│     ↓                                                        │
│  4. Assign Overhangs → getRecommendedOverhangSet()          │
│         OR → optimizeOverhangSet() (Monte Carlo)            │
│     ↓                                                        │
│  5. Calculate Fidelity → calculateExperimentalFidelity()    │
│     ↓                                                        │
│  6. Design Primers → designGoldenGatePrimers()              │
│     ↓                                                        │
│  7. Generate Protocol                                        │
└─────────────────────────────────────────────────────────────┘
```

### 1.3 Gaps for NEB Replacement 🔴

| Missing Feature | Priority | Implementation Effort |
|-----------------|----------|----------------------|
| Vector selection/design | High | Medium |
| Part library management | High | Medium |
| Batch assembly planning | Medium | Low |
| Colony screening suggestions | Medium | Low |
| Integration with inventory systems | Low | High |

---

## 2. Current Gibson Assembly Capabilities

### 2.1 Strengths ✅

| Feature | Implementation | Notes |
|---------|---------------|-------|
| **Dynamic Programming Assembly** | DAG traversal for optimal fragment selection | Optimal cost-based planning |
| **Cost Optimization** | IDT gBlocks pricing, PCR costs, assembly costs | Real-world pricing model |
| **Primer Design with Homology Tails** | Automatic homology region addition | 15-120bp overlap handling |
| **Synthetic Fragment Bridging** | Auto-generates gBlocks for gaps | Up to 3000bp synthetic |
| **Junction Validation** | Hairpin Tm checking at junctions | Max 47°C hairpin Tm |

### 2.2 Key Algorithms

```
Gibson Assembly Planning:
┌─────────────────────────────────────────────────────────────┐
│  1. Input Target + Fragment Database                        │
│     ↓                                                        │
│  2. Find Matches → findMatches()                            │
│     ↓                                                        │
│  3. Build Assembly DAG → createAssemblies()                 │
│     - Dynamic programming fragment selection                 │
│     - Cost estimation per junction                          │
│     ↓                                                        │
│  4. Fill Assemblies → fillAssemblies()                      │
│     - Design primers with homology tails                    │
│     - Add synthetic bridges                                 │
│     ↓                                                        │
│  5. Validate Junctions → validateJunctions()                │
│     ↓                                                        │
│  6. Return Ranked Solutions (by cost)                       │
└─────────────────────────────────────────────────────────────┘
```

### 2.3 Gaps for NEB Replacement 🔴

| Missing Feature | Priority | Implementation Effort |
|-----------------|----------|----------------------|
| **NEBuilder HiFi specific Tm calculation** | High | Low (already have Q5 Tm) |
| **Overlap optimizer** | High | Medium |
| **Fragment order optimization** | Medium | Medium |
| **Assembly simulation/visualization** | Medium | Medium |
| **Protocol generation (NEBuilder format)** | Low | Low |

---

## 3. Primer Design Algorithms - Major Competitive Advantages

### 3.1 Thermodynamic Foundation

The codebase implements **state-of-the-art** thermodynamic calculations:

#### Nearest-Neighbor Parameters
| Parameter Set | Source | Accuracy |
|--------------|--------|----------|
| SantaLucia 1998 | Classic reference | Baseline |
| **DNA24 (Greenleaf 2024)** | Latest research | **50% better mismatch accuracy** |

```javascript
// DNA24 has 25,000+ parameters vs SantaLucia's ~500
// Includes:
// - 1,062 tetraloop entries (vs 130)
// - Context-dependent G-T wobble parameters
// - TERMINAL_MM, HAIRPIN_MM, INTERNAL_MM tables
```

#### Pythia Equilibrium Model
```javascript
// Models 7 competing chemical species:
Species = {
  PrimerFree,        // Unfolded primer
  PrimerHairpin,     // Self-folded
  PrimerHomodimer,   // Self-dimer
  PrimerTemplate,    // Desired binding
  PrimerOffTarget,   // Mispriming
  PrimerHeterodimer, // Cross-dimer (fwd-rev)
  Solvent
}

// Calculates: η = [Primer·Target] / [Primer_total]
// This predicts actual PCR efficiency, not just Tm
```

### 3.2 Scoring System

The **piecewise logistic scoring** system provides biologically-meaningful scores:

```
Score
  1.0 ─────────────────┬───────────────────────
                      ╱ │ ╲  Optimal Zone
  0.7 ──────────────╱───┴────╲──────────────
                   ╱   Linear  ╲
  0.0 ────────────╱──────────────╲──────────
             Logistic         Logistic
```

#### Calibrated Weights (Grid-search optimized)

| Feature | Weight | Discrimination Power |
|---------|--------|---------------------|
| offTarget | 0.10 | +0.515 (highest) |
| heterodimer | 0.10 | High |
| terminal3DG | 0.05 | +0.487 |
| tmRev | 0.05 | |
| gcFwd | 0.06 | |
| hairpinRev | 0.05 | |

**Validation Performance:**
- F1 Score: 81.9%
- AUC: 0.848
- Dataset: 829 primer pairs (Döring)

### 3.3 Advanced Analysis Features

| Feature | Algorithm | Competitive Advantage |
|---------|-----------|----------------------|
| **G-Quadruplex Detection** | Regex + thermodynamic | Critical for Q5 polymerase |
| **Off-Target Classification** | Type A-F risk categories | More granular than Primer3 |
| **3' End Composition** | Multi-factor scoring | GC clamp + ΔG + pattern |
| **Smart Primer Optimization** | Iterative 3' end adjustment | Auto-improves weak primers |

---

## 4. Leveraging Algorithms for NEB Replacement

### 4.1 Direct Algorithm Reuse

```
Current Algorithm         →  NEB Tool Feature
────────────────────────────────────────────────────────
calculateTmQ5()           →  NEBuilder Tm Calculator
designGoldenGatePrimers() →  Golden Gate Primer Designer
optimizeOverhangSet()     →  Overhang Selection Tool
calculateExperimentalFidelity() → Assembly Fidelity Predictor
suggestDomestication()    →  NEBChanger Alternative
calculateHeterodimerDG()  →  Primer Pair Compatibility Check
```

### 4.2 Unique Capabilities Beyond NEB

| Feature | Our Implementation | NEB Status |
|---------|-------------------|------------|
| **Multi-enzyme fidelity comparison** | compareEnzymeFidelity() | Not available |
| **Batch overhang scoring** | batchScoreRandomSets() | Not available |
| **Cross-ligation risk matrix** | findProblematicPairs() | Hidden |
| **DNA24 mismatch accuracy** | dna24.js | Not used |
| **G-Quadruplex warnings** | analyzeGQuadruplex() | Not available |

### 4.3 Integration Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                  Unified Assembly Designer                   │
├──────────────────┬──────────────────┬───────────────────────┤
│   Golden Gate    │     Gibson       │    NEBuilder HiFi     │
│   Module         │     Module       │    Module             │
├──────────────────┴──────────────────┴───────────────────────┤
│                 Shared Primer Design Core                    │
│  ┌─────────────┬──────────────┬─────────────┬─────────────┐ │
│  │ Tm Calc     │ ΔG Calc      │ Off-Target  │ Equilibrium │ │
│  │ (Q5/DNA24)  │ (Zuker/NN)   │ Analysis    │ Efficiency  │ │
│  └─────────────┴──────────────┴─────────────┴─────────────┘ │
├─────────────────────────────────────────────────────────────┤
│              Experimental Ligation Data (473MB)              │
│              BsaI-HFv2 | BbsI-HF | BsmBI-v2 | Esp3I | SapI  │
└─────────────────────────────────────────────────────────────┘
```

---

## 5. Recommended Enhancements

### 5.1 High Priority (NEB Parity)

#### 1. NEBuilder HiFi Assembly Module
```javascript
// Extend Gibson module with NEB-specific features:
export function designNEBuilderAssembly(fragments, options) {
  // Use Q5 Tm calculation (already implemented)
  // Optimize overlaps for NEBuilder conditions:
  //   - 15-30bp overlaps (NEBuilder sweet spot)
  //   - Tm 48-65°C for overlaps
  //   - GC 40-60% in overlap region
  return {
    fragments,
    overlaps,
    protocol: generateNEBuilderProtocol(),
    estimatedEfficiency
  };
}
```

#### 2. Overlap Optimization Algorithm
```javascript
// Add overlap region optimization
export function optimizeOverlaps(assembly, options) {
  const { targetTm = 55, minLength = 15, maxLength = 30 } = options;

  for (const junction of assembly.junctions) {
    // Find overlap that maximizes:
    // 1. Tm near target
    // 2. Low secondary structure
    // 3. Balanced GC content
    junction.optimizedOverlap = findBestOverlap(junction, targetTm);
  }
}
```

#### 3. Assembly Simulation
```javascript
// Predict assembly outcome
export function simulateAssembly(design) {
  return {
    expectedYield: calculateExpectedYield(design),
    potentialIssues: identifyRisks(design),
    recommendations: generateRecommendations(design),
    alternativeDesigns: suggestAlternatives(design)
  };
}
```

### 5.2 Medium Priority (Competitive Advantage)

#### 4. Part Library System
```javascript
// MoClo-compatible part management
export const PartLibrary = {
  addPart(part, options) {},
  findCompatible(leftOverhang, rightOverhang) {},
  validateMoCloCompliance(part) {},
  exportToGenBank(parts) {}
};
```

#### 5. Batch Assembly Planning
```javascript
// Multi-construct optimization
export function planBatchAssembly(constructs, options) {
  // Identify shared parts
  // Optimize primer ordering (minimize synthesis)
  // Schedule assembly reactions
  return {
    sharedParts,
    uniqueParts,
    primerOrder,
    assemblySchedule,
    totalCost
  };
}
```

### 5.3 Low Priority (Nice to Have)

- Colony screening PCR primer suggestions
- Restriction digest verification primers
- Integration with Benchling/Geneious APIs
- Inventory management hooks

---

## 6. Competitive Analysis: Our Tool vs NEB

| Feature | NEB Tools | Our Implementation | Advantage |
|---------|-----------|-------------------|-----------|
| **Tm Calculation** | Proprietary | Q5-calibrated + DNA24 | ✅ We expose methodology |
| **Fidelity Data** | Black box | 473MB experimental data | ✅ Full transparency |
| **Overhang Selection** | Curated sets | Monte Carlo optimization | ✅ Custom optimization |
| **Domestication** | Manual lookup | Auto-suggestion with silent mutations | ✅ Automated |
| **Multi-enzyme** | One at a time | Side-by-side comparison | ✅ More efficient |
| **Primer Design** | Basic | Equilibrium + piecewise scoring | ✅ More accurate |
| **Cost Estimation** | None | Full cost model | ✅ Better planning |
| **Open Source** | No | Yes | ✅ Customizable |

---

## 7. Domestication Algorithm: Silent Mutation vs Junction-Based

### 7.1 Problem Statement

Internal restriction sites within coding sequences interfere with Golden Gate assembly. Traditional approaches have limitations:

| Approach | Method | Issue |
|----------|--------|-------|
| **Manual NEBChanger** | User finds silent mutations | Time-consuming, error-prone |
| **Junction-Based** | Split sequence at internal site | Recreates site in assembled product |

**Critical Insight**: Junction-based domestication is **incompatible with one-pot Golden Gate assembly**.

### 7.2 The One-Pot Problem with Junction-Based Domestication

```
ONE-POT GOLDEN GATE REACTION
════════════════════════════════════════════════════════════

Standard Protocol:
  Fragments + BsaI + T4 Ligase + Buffer
  Thermocycling: 37°C (cut) ↔ 16°C (ligate)

WHAT HAPPENS WITH JUNCTION-BASED DOMESTICATION:
────────────────────────────────────────────────────────────

1. Fragments ligate correctly at 16°C
   Fragment 1: ...ATCGGG─TC────┐
                         TCNN   │ Ligate
   Fragment 2:       ────TCNN───┘

2. Result: Internal site is RECREATED
   ...ATCGGG[TCTC]NNNN... = GGTCTC (BsaI site!)

3. At 37°C: BsaI cuts the recreated site AGAIN

4. Endless cycle: ligate → cut → ligate → cut...

RESULT: Reduced efficiency, incomplete assemblies
```

### 7.3 Solution: Silent Mutation-Based Domestication

The new `silent-mutation-domesticator.js` module implements a scientifically correct approach:

```
SILENT MUTATION APPROACH
════════════════════════════════════════════════════════════

Original:     ATG AAA GGT CTC AAA TGA
                      ↑↑↑ ↑↑↑
                      BsaI site (GGTCTC)

Step 1: Identify codons overlapping the site
        GGT = Gly (G)    CTC = Leu (L)

Step 2: Find synonymous codon that breaks site
        GGT → GGC (still Gly)    OR
        CTC → CTG (still Leu)

Step 3: Apply minimal mutation
        ATG AAA GGC CTC AAA TGA  ← Site broken!
                  ↑
             Single nucleotide change

Result:
  ✓ Same protein sequence
  ✓ No internal BsaI site
  ✓ Same number of fragments
  ✓ Compatible with one-pot assembly
```

### 7.4 Algorithm Implementation

```javascript
// New unified domestication optimizer
import { optimizeDomestication } from './domestication-optimizer.js';

const result = optimizeDomestication(sequence, 'BsaI', {
  frame: 0,                    // Reading frame
  isCodingSequence: true,      // Enable silent mutations
  organism: 'ecoli',           // Codon usage table
  allowJunctionFallback: true, // Fallback for non-coding
});

// Result includes:
// - domesticatedSequence: Modified sequence with sites removed
// - strategy: 'silent_mutation' | 'junction_based' | 'hybrid'
// - mutations: Array of applied mutations with details
// - warnings: Any compatibility warnings
```

### 7.5 Strategy Priority

The optimizer uses this priority order (ALL are one-pot compatible!):

| Priority | Strategy | When to Use | Fragment Count | One-Pot |
|----------|----------|-------------|----------------|---------|
| 1 (Best) | **Direct Primer Mutation** | Site near existing junction | Same | ✅ Yes |
| 2 | **Mutagenic Junction** | Site in middle of fragment | +1 per site | ✅ Yes |
| 3 | **Alternative Enzyme** | If available without sites | Same | ✅ Yes |
| ~~4~~ | ~~Legacy Junction~~ | ~~DEPRECATED~~ | ~~+1 per site~~ | ❌ No |

### 7.5.1 Mutagenic Junction Strategy (NEW)

For internal sites in the **middle of fragments** (too far from existing junctions), we use a hybrid approach:

```
MUTAGENIC JUNCTION STRATEGY
════════════════════════════════════════════════════════════════

Problem: Site is 500bp from nearest junction - primer can't reach it

Solution: Place new junction AT the site, but use mutagenic primers

Step 1: Place junction within/near the internal site
    ...ATGAAA[GGT|CTC]AAATGA...
               ↑
          Junction here

Step 2: Design primers with silent mutations in homology region

Fragment 1 Reverse Primer:
    5'-[Flank]-[BsaI]-N-[NNNN]-[...ATGAAAGCC]-3'
                       free     ↑ Mutation in
                       overhang   homology region
                       choice!    (GGT→GGC, Gly)

Fragment 2 Forward Primer:
    5'-[Flank]-[BsaI]-N-[NNNN]-[CTCAAATGA...]-3'
                       same
                       overhang

Step 3: PCR incorporates the mutation into fragments

Step 4: Golden Gate assembly

Result: ...ATGAAAGCC[NNNN]CTCAAATGA...
               ↑↑↑
           GCC ≠ GGT → Site permanently removed!
           Overhang chosen for optimal fidelity!

KEY ADVANTAGES:
✅ One-pot compatible (site doesn't exist in product)
✅ Overhang freely chosen (not constrained by original)
✅ Seamless workflow (mutation during same PCR)
✅ Protein sequence preserved (silent mutation)
```

### 7.6 Mutation Scoring System

When multiple silent mutations are available, the algorithm scores them:

```
MUTATION SCORING (0-100)
════════════════════════════════════════════════════════════

Score Component              Weight    Description
─────────────────────────────────────────────────────────────
Site Breaking                Required  Must break recognition site
No New Sites Created         0.95      Critical - avoid new problems
Codon Frequency              0.30      Prefer common codons (not priority)
Position Preference          0.10      Middle of site slightly better
Wobble Position              0.20      3rd codon position preferred

Penalties:
  - Creates new restriction site:    -80 points
  - Introduces rare codon (<10/1000): -20 points

Bonuses:
  - Wobble position mutation:         +5 points
  - Maintains codon frequency:       +10 points
```

### 7.7 Codon Usage Tables

The module includes organism-specific codon usage:

```javascript
// E. coli optimized (default)
const result = optimizeDomestication(seq, 'BsaI', { organism: 'ecoli' });

// Yeast optimized
const result = optimizeDomestication(seq, 'BsaI', { organism: 'yeast' });

// Common codons are preferred to maintain expression levels
// Example: E. coli Leucine codons
//   CTG: 52.6/1000 (preferred)
//   CTA:  3.9/1000 (avoided)
```

### 7.8 Validation & Verification

The algorithm includes comprehensive verification:

```javascript
import { validateDomestication, verifyProteinSequence } from './domestication-optimizer.js';

// Verify protein sequence preserved
const proteinCheck = verifyProteinSequence(original, domesticated, frame);
// Returns: { identical: true/false, differences: [...] }

// Full validation
const validation = validateDomestication(result, 'BsaI', { frame: 0 });
// Checks:
//   1. Protein sequence identical
//   2. No internal sites remain
//   3. No new sites created (optional: all enzymes)
//   4. Sequence length unchanged
```

### 7.9 Fallback Strategy

For non-coding regions or when silent mutations aren't possible:

```javascript
const result = optimizeDomestication(sequence, 'BsaI', {
  isCodingSequence: false,  // Non-coding region
  allowJunctionFallback: true,
});

// result.strategy === 'junction_based'
// result.warnings includes:
//   "Junction-based domestication recreates internal sites in assembled product.
//    NOT compatible with one-pot Golden Gate reactions.
//    Use sequential assembly: (1) digest, (2) heat-kill enzyme, (3) ligate."
```

### 7.10 Files & API

| File | Purpose |
|------|---------|
| `silent-mutation-domesticator.js` | Core silent mutation algorithm |
| `mutagenic-junction-domesticator.js` | **NEW**: Mutagenic junction splitting |
| `domestication-optimizer.js` | Unified optimizer (recommended entry point) |
| `auto-domestication-optimizer.js` | Legacy junction-based (DEPRECATED) |

Key exports:
```javascript
// Main entry point (auto-selects best strategy)
export { optimizeDomestication } from './domestication-optimizer.js';

// Strategy selection and analysis
export {
  selectDomesticationStrategy,
  analyzeDomesticationOptions,
  generateDomesticationReport,
} from './domestication-optimizer.js';

// Mutagenic junction functions (for sites in middle of fragments)
export {
  designMutagenicJunction,
  designAllMutagenicJunctions,
} from './domestication-optimizer.js';

// Low-level silent mutation functions
export {
  findAllSilentMutationCandidates,
  scoreMutationCandidates,
  verifyProteinSequence,
} from './silent-mutation-domesticator.js';
```

### 7.11 Comparison: All Strategies

| Feature | Legacy Junction | Direct Primer | Mutagenic Junction |
|---------|-----------------|---------------|-------------------|
| One-pot compatible | ❌ No | ✅ Yes | ✅ Yes |
| Fragment increase | +1 per site | 0 | +1 per site |
| Sequence change | None | 1-2 bp | 1-2 bp |
| Protein change | None | None | None |
| Site in final product | ✅ Recreated | ❌ Removed | ❌ Removed |
| Overhang flexibility | Constrained | N/A | **Free choice** |
| Site location | Any | Near junction | **Any (middle)** |
| Recommended | ❌ No | ✅ Yes | ✅ Yes |

---

## 8. Implementation Roadmap

### Phase 1: Core Parity
- [x] **Silent mutation domestication** - Implemented in `silent-mutation-domesticator.js`
- [x] **Mutagenic junction domestication** - Implemented in `mutagenic-junction-domesticator.js`
- [ ] Unify Golden Gate and Gibson UI
- [ ] Add NEBuilder-specific Tm targets
- [ ] Implement overlap optimization
- [ ] Generate NEB-format protocols

### Phase 2: Differentiation
- [ ] Multi-enzyme comparison dashboard
- [ ] Batch assembly planner
- [ ] Part library with MoClo compliance
- [ ] Assembly simulation/visualization

### Phase 3: Polish (1-2 weeks effort)
- [ ] Export to GenBank/Benchling
- [ ] Save/load assembly projects
- [ ] Comprehensive documentation
- [ ] Tutorial workflows

---

## 9. Conclusion

This codebase has **strong foundations** to build a state-of-the-art NEB replacement tool:

1. **Thermodynamic algorithms** are already superior (DNA24, Pythia equilibrium)
2. **Ligation fidelity data** matches NEB's experimental results
3. **Monte Carlo optimization** is a direct port of NEB's algorithm
4. **Scoring system** is empirically calibrated
5. **Silent mutation domestication** - One-pot compatible, maintains assembly efficiency

**Key differentiators** over NEB:
- Open-source and customizable
- Transparent methodology
- Multi-enzyme comparison
- Advanced off-target analysis
- Cost optimization
- **Silent mutation domestication** (scientifically correct, one-pot compatible)

**Recommended next steps:**
1. Unify the assembly modules under a single interface
2. Add overlap optimization for Gibson/NEBuilder
3. Implement assembly simulation
4. Create comprehensive documentation and tutorials

The existing algorithms provide 80%+ of what's needed - the remaining work is primarily **integration and UX polish**.
