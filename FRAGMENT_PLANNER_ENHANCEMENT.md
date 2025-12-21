# Fragment Planner Enhancement Plan

## P0 Features

### 1. Design History + Undo/Redo
### 2. Constraint-Based Overhang Solver
### 3. Modern UI Overhaul

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         ENHANCED FRAGMENT PLANNER                        │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────────┐  │
│  │  Design Store   │◄──►│  History Stack  │    │  Constraint Engine  │  │
│  │  (Zustand)      │    │  (Undo/Redo)    │    │  (SAT Solver)       │  │
│  └────────┬────────┘    └─────────────────┘    └──────────┬──────────┘  │
│           │                                                │             │
│           ▼                                                ▼             │
│  ┌────────────────────────────────────────────────────────────────────┐ │
│  │                         UI LAYER (React + Tailwind)                 │ │
│  │                                                                      │ │
│  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────────────────┐   │ │
│  │  │  Timeline    │  │  Fragment    │  │  Assembly Visualization  │   │ │
│  │  │  Navigator   │  │  Editor      │  │  (Circular + Linear)     │   │ │
│  │  └──────────────┘  └──────────────┘  └──────────────────────────┘   │ │
│  │                                                                      │ │
│  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────────────────┐   │ │
│  │  │  Constraint  │  │  Overhang    │  │  Fidelity Dashboard      │   │ │
│  │  │  Panel       │  │  Optimizer   │  │  (Live Metrics)          │   │ │
│  │  └──────────────┘  └──────────────┘  └──────────────────────────┘   │ │
│  └────────────────────────────────────────────────────────────────────┘ │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 1. State Management Architecture

### Design Store (Zustand)

```typescript
// src/stores/designStore.ts

interface Fragment {
  id: string;
  name: string;
  sequence: string;
  type: 'promoter' | 'rbs' | 'cds' | 'terminator' | 'vector' | 'insert' | 'other';
  leftOverhang?: string;
  rightOverhang?: string;
  rangeStart?: number;
  rangeEnd?: number;
  // Computed
  length: number;
  gcContent: number;
}

interface AssemblyConstraints {
  minFidelity: number;           // 0-1, e.g., 0.95 for 95%
  maxFragments: number;          // e.g., 8
  enzyme: string;                // e.g., 'BsaI'
  targetOrganism: 'ecoli' | 'yeast' | 'mammalian';
  avoidSequences: string[];      // e.g., ['GGTCTC'] - sites to avoid
  preferHighFidelityOverhangs: boolean;
}

interface DesignState {
  // Core data
  fragments: Fragment[];
  constraints: AssemblyConstraints;
  assemblyMethod: 'golden_gate' | 'gibson' | 'nebuilder';

  // Computed results
  overhangs: string[];
  fidelity: number;
  warnings: string[];

  // UI state
  selectedFragmentId: string | null;
  viewMode: 'circular' | 'linear';
}

interface HistoryEntry {
  id: string;
  timestamp: number;
  label: string;            // e.g., "Added GFP fragment"
  state: DesignState;
}

interface DesignStore extends DesignState {
  // History
  history: HistoryEntry[];
  historyIndex: number;

  // Actions
  addFragment: (fragment: Omit<Fragment, 'id' | 'length' | 'gcContent'>) => void;
  updateFragment: (id: string, updates: Partial<Fragment>) => void;
  removeFragment: (id: string) => void;
  reorderFragments: (fromIndex: number, toIndex: number) => void;

  // Constraints
  updateConstraints: (updates: Partial<AssemblyConstraints>) => void;

  // History navigation
  undo: () => void;
  redo: () => void;
  canUndo: () => boolean;
  canRedo: () => boolean;

  // Optimization
  optimizeOverhangs: () => Promise<void>;

  // Snapshots
  createSnapshot: (label: string) => void;
  restoreSnapshot: (historyId: string) => void;
}
```

### History Implementation

```typescript
// Middleware for automatic history tracking
const historyMiddleware = (config) => (set, get, api) =>
  config(
    (args) => {
      const prevState = get();
      set(args);
      const nextState = get();

      // Auto-snapshot on significant changes
      if (shouldSnapshot(prevState, nextState)) {
        const label = generateLabel(prevState, nextState);
        get().createSnapshot(label);
      }
    },
    get,
    api
  );

// Determine what constitutes a "significant change"
function shouldSnapshot(prev: DesignState, next: DesignState): boolean {
  return (
    prev.fragments.length !== next.fragments.length ||
    prev.fragments.some((f, i) => f.sequence !== next.fragments[i]?.sequence) ||
    prev.assemblyMethod !== next.assemblyMethod
  );
}

// Generate human-readable labels
function generateLabel(prev: DesignState, next: DesignState): string {
  if (next.fragments.length > prev.fragments.length) {
    const added = next.fragments.find(f => !prev.fragments.find(p => p.id === f.id));
    return `Added ${added?.name || 'fragment'}`;
  }
  if (next.fragments.length < prev.fragments.length) {
    return 'Removed fragment';
  }
  return 'Modified design';
}
```

---

## 2. Constraint-Based Overhang Solver

### Solver Architecture

```typescript
// src/lib/constraintSolver.ts

interface OverhangConstraint {
  type: 'min_fidelity' | 'no_cross_ligation' | 'avoid_sequence' | 'gc_balanced';
  params: Record<string, any>;
}

interface SolverResult {
  success: boolean;
  overhangs: string[];
  fidelity: number;
  crossLigationRisks: Array<{
    overhang1: string;
    overhang2: string;
    frequency: number;
  }>;
  suggestions?: string[];
}

class OverhangConstraintSolver {
  private ligationData: LigationMatrix;
  private enzyme: string;

  constructor(enzyme: string) {
    this.enzyme = enzyme;
    this.ligationData = loadLigationData(enzyme);
  }

  /**
   * Find optimal overhang set for N fragments
   * Uses greedy algorithm with backtracking
   */
  solve(
    fragmentCount: number,
    constraints: OverhangConstraint[]
  ): SolverResult {
    const candidates = this.generateCandidates();
    const validSets = this.findValidSets(candidates, fragmentCount, constraints);

    if (validSets.length === 0) {
      return this.suggestRelaxation(constraints);
    }

    // Rank by composite score
    const ranked = validSets.sort((a, b) =>
      this.scoreSet(b, constraints) - this.scoreSet(a, constraints)
    );

    return {
      success: true,
      overhangs: ranked[0],
      fidelity: this.calculateFidelity(ranked[0]),
      crossLigationRisks: this.findCrossLigations(ranked[0]),
    };
  }

  /**
   * Check if a proposed overhang set satisfies all constraints
   */
  validate(overhangs: string[], constraints: OverhangConstraint[]): ValidationResult {
    const violations: Violation[] = [];

    for (const constraint of constraints) {
      if (!this.checkConstraint(overhangs, constraint)) {
        violations.push({
          constraint,
          message: this.getViolationMessage(constraint, overhangs),
        });
      }
    }

    return {
      valid: violations.length === 0,
      violations,
      fidelity: this.calculateFidelity(overhangs),
    };
  }

  /**
   * Suggest alternative enzyme if constraints can't be satisfied
   */
  suggestAlternativeEnzyme(
    fragmentCount: number,
    constraints: OverhangConstraint[]
  ): EnzymeSuggestion[] {
    const enzymes = ['BsaI', 'BsmBI', 'BbsI', 'Esp3I', 'SapI'];
    const suggestions: EnzymeSuggestion[] = [];

    for (const enzyme of enzymes) {
      if (enzyme === this.enzyme) continue;

      const solver = new OverhangConstraintSolver(enzyme);
      const result = solver.solve(fragmentCount, constraints);

      if (result.success) {
        suggestions.push({
          enzyme,
          fidelity: result.fidelity,
          overhangs: result.overhangs,
        });
      }
    }

    return suggestions.sort((a, b) => b.fidelity - a.fidelity);
  }
}
```

### High-Fidelity Overhang Sets (NEB Data)

```typescript
// Pre-computed high-fidelity sets from Pryor et al. 2020
const HIGH_FIDELITY_SETS = {
  BsaI: {
    4: ['GGAG', 'TACT', 'AATG', 'AGGT', 'GCTT'],  // 4 fragments + vector
    5: ['GGAG', 'TACT', 'AATG', 'AGGT', 'TTCG', 'GCTT'],
    6: ['GGAG', 'TACT', 'AATG', 'AGGT', 'TTCG', 'GGTA', 'GCTT'],
    // ... more validated sets
  },
  BsmBI: {
    // Similar validated sets
  }
};
```

---

## 3. Modern UI Components

### Component Hierarchy

```
FragmentPlanner/
├── FragmentPlannerShell.tsx       # Main container with layout
├── components/
│   ├── DesignTimeline.tsx         # Undo/redo + history navigation
│   ├── FragmentList.tsx           # Drag-drop fragment management
│   ├── FragmentCard.tsx           # Individual fragment editor
│   ├── ConstraintPanel.tsx        # Constraint configuration
│   ├── OverhangOptimizer.tsx      # Solver UI + results
│   ├── AssemblyViz/
│   │   ├── CircularPlasmid.tsx    # Circular visualization
│   │   ├── LinearAssembly.tsx     # Linear fragment view
│   │   └── JunctionDetail.tsx     # Overhang detail popup
│   ├── FidelityDashboard.tsx      # Live fidelity metrics
│   └── ActionBar.tsx              # Run/export actions
├── hooks/
│   ├── useDesignStore.ts          # Store hook
│   ├── useHistory.ts              # Undo/redo hook
│   └── useConstraintSolver.ts     # Async solver hook
└── stores/
    └── designStore.ts             # Zustand store
```

### UI Wireframe

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  Fragment Planner                                      [Undo] [Redo] [Save]  │
├───────────────────────┬─────────────────────────────────────────────────────┤
│                       │                                                      │
│  FRAGMENTS            │           ASSEMBLY VIEW                              │
│  ┌─────────────────┐  │  ┌────────────────────────────────────────────────┐ │
│  │ 1. pTrc Promoter│  │  │                                                │ │
│  │    ●━━━━━━━━━━━━┿──┼──┼───► GGAG                                       │ │
│  └─────────────────┘  │  │         ╲                                      │ │
│  ┌─────────────────┐  │  │          ╲   ┌───────────────────┐             │ │
│  │ 2. GFP CDS      │  │  │           ╲──│    Circular       │             │ │
│  │    ━━━━━━━━━━━━●┿──┼──┼───► AATG    │    Plasmid View   │             │ │
│  └─────────────────┘  │  │             └───────────────────┘             │ │
│  ┌─────────────────┐  │  │                                                │ │
│  │ 3. T7 Term      │  │  │                                                │ │
│  │    ━━━━━━━━━●   │  │  └────────────────────────────────────────────────┘ │
│  └─────────────────┘  │                                                      │
│                       │  ┌────────────────────────────────────────────────┐ │
│  [+ Add Fragment]     │  │  FIDELITY DASHBOARD                            │ │
│                       │  │  ████████████████████░░  96.2%                 │ │
├───────────────────────┤  │                                                │ │
│  CONSTRAINTS          │  │  Overhangs: GGAG → AATG → TTCG → GCTT          │ │
│  ┌─────────────────┐  │  │  Cross-ligation: None detected                 │ │
│  │ Min Fidelity    │  │  └────────────────────────────────────────────────┘ │
│  │ [████████░] 95% │  │                                                      │
│  └─────────────────┘  │  ┌────────────────────────────────────────────────┐ │
│  ┌─────────────────┐  │  │  HISTORY                                       │ │
│  │ Enzyme: BsaI ▼  │  │  │  ○ Added T7 Terminator           2 min ago    │ │
│  └─────────────────┘  │  │  ● Added GFP CDS                  5 min ago    │ │
│  ┌─────────────────┐  │  │  ○ Added pTrc Promoter           10 min ago   │ │
│  │ Max Fragments   │  │  │  ○ Initial design                15 min ago   │ │
│  │ [═══●═════] 5   │  │  └────────────────────────────────────────────────┘ │
│  └─────────────────┘  │                                                      │
│                       │                                                      │
│  [Optimize Overhangs] │           [Design Primers] [Export GenBank]         │
└───────────────────────┴─────────────────────────────────────────────────────┘
```

---

## 4. File Structure

```
src/
├── components/
│   └── FragmentPlanner/
│       ├── index.tsx                    # Main export
│       ├── FragmentPlannerShell.tsx     # Layout container
│       ├── DesignTimeline.tsx           # History UI
│       ├── FragmentList.tsx             # Fragment management
│       ├── FragmentCard.tsx             # Individual fragment
│       ├── ConstraintPanel.tsx          # Constraint config
│       ├── OverhangOptimizer.tsx        # Solver UI
│       ├── FidelityDashboard.tsx        # Metrics display
│       ├── AssemblyViz/
│       │   ├── index.tsx
│       │   ├── CircularPlasmid.tsx
│       │   ├── LinearAssembly.tsx
│       │   └── JunctionDetail.tsx
│       └── styles.css                   # SeqViz styles only
├── stores/
│   └── designStore.ts                   # Zustand store + history
├── lib/
│   ├── constraintSolver.ts              # Overhang solver
│   ├── ligationData.ts                  # NEB fidelity data
│   └── fragmentUtils.ts                 # Helper functions
├── hooks/
│   ├── useDesignStore.ts
│   ├── useHistory.ts
│   └── useConstraintSolver.ts
└── types/
    └── fragmentPlanner.ts               # Type definitions
```

---

## 5. Implementation Phases

### Phase 1: Foundation (Current Focus)
- [x] Create architecture document
- [ ] Set up project structure
- [ ] Create Zustand store with history middleware
- [ ] Implement basic undo/redo

### Phase 2: Solver Engine
- [ ] Port ligation data from existing code
- [ ] Implement constraint solver
- [ ] Add validation layer
- [ ] Create solver hook

### Phase 3: UI Shell
- [ ] Build FragmentPlannerShell layout
- [ ] Create DesignTimeline component
- [ ] Build FragmentList with drag-drop
- [ ] Add ConstraintPanel

### Phase 4: Visualization
- [ ] Port CircularPlasmid from existing code
- [ ] Create LinearAssembly view
- [ ] Build FidelityDashboard
- [ ] Add JunctionDetail popups

### Phase 5: Integration
- [ ] Connect solver to UI
- [ ] Add keyboard shortcuts (Cmd+Z, Cmd+Shift+Z)
- [ ] Implement export functionality
- [ ] Polish animations and transitions

---

## 6. Key Design Decisions

### Why Zustand over Redux?
- Simpler API, less boilerplate
- Built-in middleware support for history
- Better TypeScript integration
- Smaller bundle size

### Why custom solver over external SAT library?
- Domain-specific optimizations
- Pre-computed high-fidelity sets
- No external dependencies
- Faster for small problem sizes (< 20 overhangs)

### Why separate history from state?
- Clean separation of concerns
- Configurable snapshot granularity
- Efficient memory usage (only store diffs)
- Easy to implement branching later

---

## 7. Tailwind Design Tokens

```javascript
// tailwind.config.js additions
module.exports = {
  theme: {
    extend: {
      colors: {
        // Fragment type colors
        'frag-promoter': '#ef4444',
        'frag-rbs': '#f97316',
        'frag-cds': '#22c55e',
        'frag-terminator': '#3b82f6',
        'frag-other': '#8b5cf6',

        // Fidelity colors
        'fidelity-excellent': '#22c55e',
        'fidelity-good': '#84cc16',
        'fidelity-marginal': '#f59e0b',
        'fidelity-poor': '#ef4444',

        // UI accents
        'primary': {
          50: '#eff6ff',
          // ... full scale
          600: '#2563eb',
          700: '#1d4ed8',
        }
      }
    }
  }
}
```
