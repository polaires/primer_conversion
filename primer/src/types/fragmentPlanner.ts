/**
 * Fragment Planner Type Definitions
 * State-of-the-art plasmid assembly tool
 */

// ============================================================================
// Core Fragment Types
// ============================================================================

export type FragmentType =
  | 'promoter'
  | 'rbs'
  | 'cds'
  | 'terminator'
  | 'vector'
  | 'insert'
  | 'other';

export interface Fragment {
  id: string;
  name: string;
  sequence: string;
  type: FragmentType;

  // Overhangs (assigned by solver or manually)
  leftOverhang?: string;
  rightOverhang?: string;

  // Range selection (for partial fragment use)
  rangeStart?: number;
  rangeEnd?: number;

  // Computed properties
  length: number;
  gcContent: number;

  // Metadata
  source?: 'addgene' | 'igem' | 'dnasu' | 'custom' | 'synthetic';
  url?: string;
  cost?: number;

  // Domestication state
  domestication?: {
    needsDomestication: boolean;
    strategy?: 'silent_mutation' | 'mutagenic_junction' | 'alternative_enzyme';
    internalSites?: Array<{ position: number; enzyme: string }>;
    approved?: boolean;
  };
}

// ============================================================================
// Assembly Method Types
// ============================================================================

export type AssemblyMethod = 'golden_gate' | 'gibson' | 'nebuilder';

export type GoldenGateEnzyme = 'BsaI' | 'BsmBI' | 'BbsI' | 'Esp3I' | 'SapI';

export interface EnzymeInfo {
  name: string;
  recognition: string;
  cutDistance: number;
  overhangLength: number;
  optimalFlankingSequence: string;
  supplier: string;
}

// ============================================================================
// Constraint Types
// ============================================================================

export interface AssemblyConstraints {
  // Fidelity requirements
  minFidelity: number; // 0-1, e.g., 0.95 for 95%

  // Fragment limits
  maxFragments: number;

  // Enzyme selection
  enzyme: GoldenGateEnzyme;

  // Optimization preferences
  targetOrganism: 'ecoli' | 'yeast' | 'mammalian';
  preferHighFidelityOverhangs: boolean;
  allowDomestication: boolean;

  // Sequence constraints
  avoidSequences: string[]; // Sites to avoid in overhangs

  // Advanced
  gcRange?: { min: number; max: number };
  tmRange?: { min: number; max: number };
}

export const DEFAULT_CONSTRAINTS: AssemblyConstraints = {
  minFidelity: 0.95,
  maxFragments: 5,
  enzyme: 'BsaI',
  targetOrganism: 'ecoli',
  preferHighFidelityOverhangs: true,
  allowDomestication: true,
  avoidSequences: [],
};

// ============================================================================
// Solver Types
// ============================================================================

export interface OverhangSet {
  overhangs: string[];
  fidelity: number;
  source: 'neb_validated' | 'computed' | 'user_defined';
}

export interface CrossLigationRisk {
  overhang1: string;
  overhang2: string;
  frequency: number; // 0-1
  severity: 'low' | 'medium' | 'high';
}

export interface SolverResult {
  success: boolean;
  overhangs: string[];
  fidelity: number;
  crossLigationRisks: CrossLigationRisk[];
  alternativeEnzymes?: Array<{
    enzyme: GoldenGateEnzyme;
    fidelity: number;
    overhangs: string[];
  }>;
  suggestions?: string[];
  computationTime?: number;
}

export interface ValidationResult {
  valid: boolean;
  violations: Array<{
    type: string;
    message: string;
    severity: 'error' | 'warning';
  }>;
  fidelity: number;
}

// ============================================================================
// History Types
// ============================================================================

export interface HistoryEntry {
  id: string;
  timestamp: number;
  label: string;
  type: 'auto' | 'manual' | 'restore';
  snapshot: DesignSnapshot;
}

export interface DesignSnapshot {
  fragments: Fragment[];
  constraints: AssemblyConstraints;
  assemblyMethod: AssemblyMethod;
  overhangs: string[];
}

// ============================================================================
// Store Types
// ============================================================================

export interface DesignState {
  // Core data
  fragments: Fragment[];
  constraints: AssemblyConstraints;
  assemblyMethod: AssemblyMethod;

  // Computed/solver results
  overhangs: string[];
  fidelity: number;
  solverResult: SolverResult | null;

  // Validation
  warnings: string[];
  errors: string[];

  // UI state
  selectedFragmentId: string | null;
  viewMode: 'circular' | 'linear';
  isOptimizing: boolean;

  // History
  history: HistoryEntry[];
  historyIndex: number;
}

export interface DesignActions {
  // Fragment management
  addFragment: (fragment: Omit<Fragment, 'id' | 'length' | 'gcContent'>) => void;
  updateFragment: (id: string, updates: Partial<Fragment>) => void;
  removeFragment: (id: string) => void;
  reorderFragments: (fromIndex: number, toIndex: number) => void;
  setFragments: (fragments: Fragment[]) => void;

  // Constraint management
  updateConstraints: (updates: Partial<AssemblyConstraints>) => void;

  // Method selection
  setAssemblyMethod: (method: AssemblyMethod) => void;

  // Overhang management
  setOverhangs: (overhangs: string[]) => void;
  optimizeOverhangs: () => Promise<SolverResult>;

  // History
  undo: () => void;
  redo: () => void;
  canUndo: () => boolean;
  canRedo: () => boolean;
  createSnapshot: (label: string) => void;
  restoreSnapshot: (historyId: string) => void;
  clearHistory: () => void;

  // UI
  selectFragment: (id: string | null) => void;
  setViewMode: (mode: 'circular' | 'linear') => void;

  // Import/Export
  exportDesign: () => DesignSnapshot;
  importDesign: (snapshot: DesignSnapshot) => void;
  reset: () => void;
}

export type DesignStore = DesignState & DesignActions;

// ============================================================================
// UI Component Types
// ============================================================================

export interface FragmentCardProps {
  fragment: Fragment;
  index: number;
  isSelected: boolean;
  isDragging: boolean;
  onSelect: () => void;
  onUpdate: (updates: Partial<Fragment>) => void;
  onRemove: () => void;
  onDragStart: (e: React.DragEvent) => void;
  onDragOver: (e: React.DragEvent) => void;
  onDrop: (e: React.DragEvent) => void;
}

export interface FidelityLevel {
  label: string;
  color: string;
  minValue: number;
}

export const FIDELITY_LEVELS: FidelityLevel[] = [
  { label: 'Excellent', color: '#22c55e', minValue: 0.95 },
  { label: 'Good', color: '#84cc16', minValue: 0.85 },
  { label: 'Marginal', color: '#f59e0b', minValue: 0.70 },
  { label: 'Poor', color: '#ef4444', minValue: 0 },
];

export function getFidelityLevel(fidelity: number): FidelityLevel {
  for (const level of FIDELITY_LEVELS) {
    if (fidelity >= level.minValue) {
      return level;
    }
  }
  return FIDELITY_LEVELS[FIDELITY_LEVELS.length - 1];
}

// ============================================================================
// Fragment Type Metadata
// ============================================================================

export interface FragmentTypeMeta {
  name: string;
  color: string;
  icon: string;
  description: string;
}

export const FRAGMENT_TYPES: Record<FragmentType, FragmentTypeMeta> = {
  promoter: {
    name: 'Promoter',
    color: '#ef4444',
    icon: 'P',
    description: 'Transcription initiation',
  },
  rbs: {
    name: "RBS/5'UTR",
    color: '#f97316',
    icon: 'R',
    description: 'Ribosome binding site',
  },
  cds: {
    name: 'CDS',
    color: '#22c55e',
    icon: 'C',
    description: 'Coding sequence',
  },
  terminator: {
    name: 'Terminator',
    color: '#3b82f6',
    icon: 'T',
    description: 'Transcription termination',
  },
  vector: {
    name: 'Vector',
    color: '#6366f1',
    icon: 'V',
    description: 'Plasmid backbone',
  },
  insert: {
    name: 'Insert',
    color: '#8b5cf6',
    icon: 'I',
    description: 'Generic insert',
  },
  other: {
    name: 'Other',
    color: '#64748b',
    icon: 'O',
    description: 'Other sequence',
  },
};

// ============================================================================
// Enzyme Data
// ============================================================================

export const GOLDEN_GATE_ENZYMES: Record<GoldenGateEnzyme, EnzymeInfo> = {
  BsaI: {
    name: 'BsaI-HFv2',
    recognition: 'GGTCTC',
    cutDistance: 1,
    overhangLength: 4,
    optimalFlankingSequence: 'CGTCTC',
    supplier: 'NEB R3733',
  },
  BsmBI: {
    name: 'BsmBI-v2',
    recognition: 'CGTCTC',
    cutDistance: 1,
    overhangLength: 4,
    optimalFlankingSequence: 'CGTCTC',
    supplier: 'NEB R0739',
  },
  BbsI: {
    name: 'BbsI-HF',
    recognition: 'GAAGAC',
    cutDistance: 2,
    overhangLength: 4,
    optimalFlankingSequence: 'GAAGAC',
    supplier: 'NEB R3539',
  },
  Esp3I: {
    name: 'Esp3I',
    recognition: 'CGTCTC',
    cutDistance: 1,
    overhangLength: 4,
    optimalFlankingSequence: 'CGTCTC',
    supplier: 'NEB R0734',
  },
  SapI: {
    name: 'SapI',
    recognition: 'GCTCTTC',
    cutDistance: 1,
    overhangLength: 3,
    optimalFlankingSequence: 'GCTCTTC',
    supplier: 'NEB R0569',
  },
};
