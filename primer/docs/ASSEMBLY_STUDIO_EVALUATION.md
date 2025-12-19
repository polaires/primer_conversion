# Assembly Studio Evaluation & Improvement Roadmap

## Executive Summary

This document evaluates the current Primer Assembly Studio implementation and provides a comprehensive roadmap to build a **state-of-the-art primer designer** that can compete with and surpass NEB HiFi Assembly and Golden Gate Assembly tools.

**Current Status**: Production-ready core with 36,000+ lines of code across 28 modules. The system already exceeds NEB tools in transparency and algorithmic sophistication, but requires strategic enhancements to achieve market leadership.

**Target Use Case**: Plasmid assembly from pure DNA fragments. This simplifies off-target concerns significantly - we only need to check against the plasmid sequences in the current assembly, not entire genomes.

---

## Current Implementation Assessment

### Architecture Overview

| Component | Files | Lines | Status |
|-----------|-------|-------|--------|
| Core Thermodynamics | `tm.js`, `tmQ5.js`, `dna24.js`, `equilibrium.js` | ~4,000 | Production |
| Scoring System | `scoring.js`, `weightCalibration.js`, `presets.js` | ~2,500 | Production |
| Primer Design | `primers.js`, `smartPrimers.js`, `unifiedPrimerDesign.js` | ~5,000 | Production |
| Assembly Modules | `nebuilder.js`, `goldengate.js`, `assemblyCore.js` | ~4,500 | Production |
| UI Components | 17 React components | ~8,000 | Production |
| Tests | 15 test files | ~8,000 | Good coverage |

### Scoring System Analysis

**Strengths:**
- Piecewise logistic functions with biologically-meaningful thresholds
- Calibrated weights from 829-pair Döring validation dataset
- Performance: **F1=88.7%**, **AUC=0.930** (cross-validated)
- Domain-specific presets for amplification, mutagenesis, sequencing, assembly

**Current Weight Distribution:**
```javascript
// Top discriminative features (from calibration)
offTarget:      0.25  // Most important - off-target binding
terminal3DG:    0.20  // 3' end stability
gQuadruplexRev: 0.15  // G-quadruplex risk
hairpinFwd:     0.10  // Secondary structure
tmFwd:          0.08  // Melting temperature
```

**Gaps Identified:**
1. No mismatch tolerance modeling at binding sites
2. G-quadruplex detection misses some non-canonical motifs
3. Internal mispriming detection within plasmid could be more robust

### NEBuilder HiFi Module Analysis

**Strengths:**
- Complete overlap optimization with sliding window
- Comprehensive scoring: Tm, GC, hairpin, patterns, GC-clamp
- Automatic primer design with homology tails
- Full protocol generation matching NEB format

**Current Capabilities:**
```javascript
NEBUILDER_PARAMS = {
  overlap: { min: 15, max: 35, optimal: 20, optimalTm: 55 },
  incubation: { '2-3': 15min, '4+': 60min },
  fragmentAmounts: { '2-3': 0.03-0.2 pmol, '4-6': 0.2-0.5 pmol }
}
```

**Gaps Identified:**
1. No assembly yield prediction model
2. Missing recombination frequency analysis
3. No synthetic bridge/linker generation for difficult junctions

### Golden Gate Assembly Module Analysis

**Strengths:**
- Complete Type IIS enzyme support (BsaI, BbsI, BsmBI, Esp3I, SapI)
- Experimental ligation data integration (Pryor et al. 2020)
- High-fidelity overhang sets (Potapov et al. 2018) up to 30 parts
- Domestication suggestions for internal sites
- Monte Carlo overhang optimization

**Current Overhang Sets:**
| Set | Parts | Fidelity | Source |
|-----|-------|----------|--------|
| MoClo Standard | 4 | 96% | Weber et al. 2011 |
| NEB Level 1 | 8 | 93% | NEB Application Note |
| Potapov Set 2 | 20 | 98% | ACS Synth Biol 2018 |
| Extended | 30 | 90% | Potapov et al. 2018 |

**Gaps Identified:**
1. 473MB ligation matrix not fully integrated (only ~1.1MB loaded)
2. No real-time cross-reactivity visualization
3. Missing assembly simulation/yield prediction
4. No support for hierarchical (Level 0→1→2) assembly planning

---

## Competitive Analysis: NEB vs Current System

### Feature Comparison

| Feature | NEB Tools | Our System | Advantage |
|---------|-----------|------------|-----------|
| **Transparency** | Black box | Open algorithms | **Ours** |
| **Ligation Matrices** | Hidden | Exposed raw data | **Ours** |
| **Internal Site Detection** | Basic | Comprehensive + domestication | **Ours** |
| **Domestication** | Manual | Automated with suggestions | **Ours** |
| **Scoring** | Proprietary | Published & calibrated | **Ours** |
| **Thermodynamics** | SantaLucia | DNA24 + Pythia equilibrium | **Ours** |
| **Multi-enzyme Comparison** | No | Cross-enzyme fidelity | **Ours** |
| **Assembly Simulation** | Yes | Not implemented | NEB |
| **Batch Processing** | Yes | Single assembly | NEB |
| **API Access** | Commercial | Not available | NEB |
| **User Accounts** | Yes | Not implemented | NEB |

### What NEB Does Better (Current Gaps)

1. **Assembly yield prediction** - Experimental correlation data
2. **Batch processing** - Design 100s of primers at once
3. **Commercial-grade API** - REST endpoints for automation
4. **User experience** - Polished, well-tested UI

> **Note on Off-Target**: For plasmid assembly from pure DNA, genomic off-target databases (BLAST) are unnecessary. We only need to check primer specificity against the plasmid sequences in the current assembly. The existing `findInternalSites()` and primer-dimer checking already covers the critical cases.

---

## Improvement Roadmap

### Phase 1: Critical Infrastructure (High Priority)

#### 1.1 Full Ligation Matrix Integration ✅ COMPLETE
**Status:** Already integrated - 473KB with full 256×256 matrices for all 5 enzymes
**Extended:** Pre-computed optimal sets now available for 2-30 parts

**What's Available:**
```javascript
// Full ligation data for all enzymes
import ligationData from './ligation-data.json';
// BsaI-HFv2, BsmBI-v2, Esp3I, BbsI-HF, SapI

// Optimal sets from 2-30 parts
ligationData.enzymes['BsaI-HFv2'].optimalSets
// { 2: {fidelity: 98.6%}, ..., 30: {fidelity: 13.3%} }

// Cross-reactivity analysis (NEW)
import { generateHeatmapData, calculateJunctionFidelities }
  from './repp/cross-reactivity.js';
```

**Added Capabilities:**
- Pre-computed optimal sets up to 30 parts
- Cross-reactivity heatmap generation
- Per-junction fidelity analysis
- Better alternative overhang suggestions

#### 1.2 DNA24 Parameter Validation ✅ COMPLETE
**Status:** Comprehensive validation suite created and passing
**Run:** `node scripts/validate-dna24.js`

**Validation Results (20/20 tests passing):**
```
Parameter Completeness:
  ✓ All 16 Watson-Crick NN pairs present
  ✓ 576+ context-dependent internal mismatch entries
  ✓ Expanded tetraloops dataset (vs SantaLucia 130)
  ✓ Hairpin loops 3-30 defined

Thermodynamic Values:
  ✓ dH values in range (-15 to +10 kcal/mol)
  ✓ dS values in range (-50 to +30 cal/mol/K)
  ✓ CG/GC more stable than AT/TA (expected)

DNA24 vs SantaLucia:
  ✓ Watson-Crick pairs agree within 30%
  AA/TT: DNA24 dH=-7.9, SantaLucia dH=-7.6 (ratio 1.04)

Performance:
  ✓ Q5 Tm calc < 1ms (0.015ms average)
```

### Phase 2: Quality & Accuracy Improvements (Interactive Design Focus)

> **Note:** This tool targets single-user interactive design, not large-scale automation.
> Batch processing, REST API, and collaboration features are deprioritized.

#### 2.1 Assembly Yield Prediction Model
**Goal:** Help users understand if their design will work before ordering

```javascript
function predictAssemblyYield(assembly) {
  return {
    overallFidelity: 0.85,           // Product of junction fidelities
    expectedCorrectRate: '85%',       // % of colonies with correct assembly
    coloniesToScreen: 12,             // Recommended colonies to pick
    riskFactors: [
      { type: 'low_fidelity_junction', junction: 3, fidelity: 0.91 },
      { type: 'large_fragment', fragment: 2, size: 4500 }
    ],
    confidence: 'high'                // Based on data quality
  };
}
```

#### 2.2 Synthetic Bridge Generation
**Goal:** Automatically solve difficult junctions with synthetic oligos

```javascript
function designSyntheticBridge(junction) {
  // For junctions with:
  // - Poor overlap Tm (<50°C)
  // - High secondary structure (ΔG < -5 kcal/mol)
  // - Repeat sequences causing mispriming

  return {
    bridgeSequence: 'ATGC...80bp...GCTA',
    leftOverlap: { tm: 62, dg: -2.1 },
    rightOverlap: { tm: 60, dg: -1.8 },
    orderAs: 'gBlock or Ultramer',
    estimatedCost: '$15-45'
  };
}
```

#### 2.3 Hierarchical Assembly (MoClo/Modular Cloning)
**Goal:** Support multi-level Golden Gate workflows

```javascript
// Level 0: Basic parts (promoters, RBS, CDS, terminators)
// Level 1: Transcription units (assembled from Level 0)
// Level 2: Multigene constructs (assembled from Level 1)

function planHierarchicalAssembly(parts) {
  return {
    levels: [
      { name: 'Level 0', enzyme: 'BsaI', parts: basicParts },
      { name: 'Level 1', enzyme: 'BsmBI', parts: transcriptionUnits },
      { name: 'Level 2', enzyme: 'BsaI', parts: multiGene }
    ],
    standardOverhangs: MOCLO_STANDARD_SET,
    compatibilityCheck: validateMoCloCompatibility(parts)
  };
}
```

#### 2.4 Enhanced Primer Quality Scoring
**Goal:** Improve primer success prediction beyond current F1=88.7%

- Add mismatch tolerance modeling at 3' end
- Improve G-quadruplex detection (non-canonical motifs)
- Better secondary structure prediction for long primers
- Context-dependent scoring for assembly vs amplification

### Phase 3: Competitive Differentiation (Innovation)

#### 3.1 Machine Learning Primer Scoring
**Goal:** Train model on experimental PCR success data

```javascript
// Collect training data
const trainingData = {
  döringDataset: 829,  // Already integrated
  additionalSources: [
    'IDT primer QC data',
    'Published primer validation studies',
    'User-submitted results'
  ]
};

// Features for ML model
const features = [
  'tm', 'gc', 'length', 'hairpinDG', 'homodimerDG',
  'terminal3DG', 'gcClamp', 'homopolymerRuns',
  'gQuadruplexScore', 'offTargetCount', 'templateGC',
  'ampliconLength', 'templateSecondaryStructure'
];

// Model architecture
// 1. Gradient boosting for interpretability
// 2. Neural network for maximum accuracy
// 3. Ensemble for production
```

#### 3.2 Synthetic Bridge Generation
**Goal:** Auto-design synthetic bridges for difficult junctions

```javascript
function designSyntheticBridge(junction, options = {}) {
  // For junctions with:
  // - Poor overlap Tm
  // - High secondary structure
  // - Repeat sequences

  return {
    bridgeSequence: optimizedBridge,
    leftOverlap: { seq, tm, dg },
    rightOverlap: { seq, tm, dg },
    synthesisComplexity: 'low|medium|high',
    orderingInfo: { vendor: 'IDT', product: 'gBlock' }
  };
}
```

#### 3.3 Hierarchical Assembly Planning
**Goal:** Multi-level Golden Gate (MoClo, Modular Cloning)

```javascript
function planHierarchicalAssembly(parts, options = {}) {
  // Level 0: Basic parts (promoters, RBS, CDS, terminators)
  // Level 1: Transcription units (assembled from Level 0)
  // Level 2: Multigene constructs (assembled from Level 1)

  return {
    levels: [
      { parts: level0Parts, enzyme: 'BsaI', overhangs: mocoStandard },
      { parts: level1TUs, enzyme: 'BsmBI', overhangs: level1Set },
      { parts: level2Construct, enzyme: 'BsaI', overhangs: level2Set }
    ],
    timeline: estimateAssemblyTimeline(levels),
    cost: estimateSynthesisCost(levels)
  };
}
```

#### 3.4 Real-time Collaboration
**Goal:** Multiple users designing assembly simultaneously

```javascript
// WebSocket-based collaboration
const assemblySession = {
  id: 'session_123',
  users: ['user1', 'user2'],
  assembly: sharedAssemblyState,
  cursor: { user1: 'fragment_2', user2: 'junction_4' },
  chat: [],
  history: []  // Undo/redo stack
};
```

### Phase 4: Platform & Infrastructure

#### 4.1 REST API Layer
```javascript
// API endpoints
POST   /api/primers/design          // Design primers for template
POST   /api/primers/batch           // Batch primer design
POST   /api/assembly/golden-gate    // Plan Golden Gate assembly
POST   /api/assembly/nebuilder      // Plan NEBuilder assembly
GET    /api/assembly/:id            // Get saved assembly
POST   /api/score                   // Score primer pair
GET    /api/enzymes                 // List supported enzymes
GET    /api/overhangs/:enzyme       // Get overhang fidelity data
```

#### 4.2 Database Schema
```sql
-- Core entities
CREATE TABLE users (id, email, name, organization, created_at);
CREATE TABLE projects (id, user_id, name, description, created_at);
CREATE TABLE assemblies (id, project_id, type, config, status, created_at);
CREATE TABLE fragments (id, assembly_id, sequence, position, primers);
CREATE TABLE primers (id, fragment_id, type, sequence, tm, score);

-- Analytics
CREATE TABLE primer_results (id, primer_id, success, notes, created_at);
CREATE TABLE assembly_results (id, assembly_id, colonies, correct_rate, created_at);
```

#### 4.3 Export Formats
```javascript
const exportFormats = {
  'csv': standardPrimerTable,
  'tsv': idtOrderFormat,
  'xlsx': multiSheetWorkbook,
  'genbank': annotatedSequence,
  'fasta': sequenceOnly,
  'json': fullAssemblyData,
  'snapgene': snapGeneFormat,  // .dna file
  'benchling': benchlingJSON   // API compatible
};
```

---

## Performance Benchmarks

### Current Performance
| Operation | Time | Target |
|-----------|------|--------|
| Primer design (20bp) | <100ms | Maintain |
| Hairpin ΔG calculation | 5-20ms | Maintain |
| Golden Gate (5 parts) | 50-100ms | Maintain |
| NEBuilder (6 parts) | 100-200ms | Maintain |
| Overhang optimization (20 parts) | 500-2000ms | <500ms |

### Optimization Opportunities
1. **Web Workers** for parallel primer scoring
2. **Wasm** for thermodynamic calculations (10x speedup potential)
3. **IndexedDB** caching for ligation matrices
4. **Service Worker** for offline support

---

## Implementation Priority Matrix

| Feature | Impact | Effort | Priority | Status |
|---------|--------|--------|----------|--------|
| Full ligation data | High | Medium | P1 | ✅ Complete |
| DNA24 parameter validation | High | Low | P1 | ✅ Complete |
| Cross-reactivity visualization | High | Medium | P1 | ✅ Complete |
| **Assembly yield prediction** | High | Medium | **P2** | 🔄 In Progress |
| **Synthetic bridge generation** | High | Medium | **P2** | Pending |
| **Hierarchical assembly (MoClo)** | Medium | Medium | **P2** | Pending |
| **Enhanced primer scoring** | Medium | Low | **P2** | Pending |
| Codon optimization | Medium | High | P3 | Pending |
| ML scoring model | Medium | High | P3 | Pending |
| ~~Batch processing~~ | ~~High~~ | ~~Medium~~ | ~~P4~~ | Deprioritized |
| ~~REST API~~ | ~~High~~ | ~~Medium~~ | ~~P4~~ | Deprioritized |
| ~~Real-time collaboration~~ | ~~Low~~ | ~~High~~ | ~~P4~~ | Deprioritized |

---

## Success Metrics

### Technical Metrics
- **Primer Success Rate:** Target 95%+ (vs 85% current benchmark)
- **Scoring Accuracy:** Target F1 > 92% (vs 88.7% current)
- **Assembly Fidelity Prediction:** Target ±5% of actual (not measured currently)

### User Metrics
- **Design Time:** <30 seconds for standard assembly
- **Error Rate:** <1% incorrect primer designs
- **User Satisfaction:** NPS > 50

### Business Metrics
- **Active Users:** Track monthly active designers
- **Designs Per User:** Average assemblies per month
- **API Calls:** For integration adoption

---

## Conclusion

The current Assembly Studio has a **strong foundation** that already exceeds NEB tools in algorithm transparency and theoretical sophistication. For **plasmid assembly from pure DNA**, genomic off-target concerns are minimal - we only need internal plasmid checking which is already implemented.

### Phase 1 Complete ✅

All Phase 1 items have been implemented:

1. **Ligation Matrix Integration** ✅
   - Full 256×256 matrices for all 5 enzymes (BsaI, BsmBI, Esp3I, BbsI, SapI)
   - Pre-computed optimal sets extended to 30 parts
   - New cross-reactivity analysis module

2. **DNA24 Thermodynamic Validation** ✅
   - 20/20 validation tests passing
   - Parameter completeness verified
   - Performance validated (<1ms per Tm calc)

3. **Cross-Reactivity Visualization** ✅
   - Heatmap data generation
   - Per-junction fidelity analysis
   - Better alternative suggestions

### Phase 2 Focus: Quality & Accuracy (Interactive Design)

**New priorities** (batch/API/collaboration deprioritized):

1. **Assembly yield prediction** - Help users know if design will work
2. **Synthetic bridge generation** - Solve difficult junctions automatically
3. **Hierarchical assembly (MoClo)** - Multi-level Golden Gate workflows
4. **Enhanced primer scoring** - Improve beyond F1=88.7%

**Recommended Next Steps:**
1. Implement assembly yield prediction model (Phase 2.1)
2. Add synthetic bridge generation for difficult junctions (Phase 2.2)
3. Support MoClo/modular cloning workflows (Phase 2.3)
