/**
 * Enhanced Fragment Planner
 * State-of-the-art plasmid assembly tool
 *
 * Features:
 * - Design history with undo/redo (Ctrl+Z / Ctrl+Shift+Z)
 * - Constraint-based overhang solver (NEB experimental data)
 * - Drag-and-drop fragment management
 * - Circular and linear visualization
 * - Live fidelity metrics
 */

export { FragmentPlannerShell as FragmentPlanner } from './FragmentPlannerShell';
export { DesignTimeline } from './DesignTimeline';
export { FragmentList } from './FragmentList';
export { ConstraintPanel } from './ConstraintPanel';
export { FidelityDashboard } from './FidelityDashboard';
export { ActionBar } from './ActionBar';
export { CircularPlasmidView, LinearAssemblyView } from './AssemblyViz';

// Re-export store hooks
export { useDesignStore } from '../../stores/designStore';

// Re-export types
export type {
  Fragment,
  FragmentType,
  AssemblyConstraints,
  AssemblyMethod,
  GoldenGateEnzyme,
  SolverResult,
  HistoryEntry,
  DesignSnapshot,
} from '../../types/fragmentPlanner';
