/**
 * Design Store - Zustand state management with undo/redo
 * Manages fragment planner state with automatic history tracking
 */

import { create } from 'zustand';
import { subscribeWithSelector } from 'zustand/middleware';
import { immer } from 'zustand/middleware/immer';
import {
  DesignStore,
  DesignState,
  Fragment,
  AssemblyConstraints,
  AssemblyMethod,
  HistoryEntry,
  DesignSnapshot,
  SolverResult,
  DEFAULT_CONSTRAINTS,
} from '../types/fragmentPlanner';
import { solveOverhangs } from '../lib/constraintSolver';

// ============================================================================
// Utility Functions
// ============================================================================

function generateId(): string {
  return `frag_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
}

function generateHistoryId(): string {
  return `hist_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
}

function calculateGcContent(sequence: string): number {
  if (!sequence) return 0;
  const gc = (sequence.toUpperCase().match(/[GC]/g) || []).length;
  return (gc / sequence.length) * 100;
}

function createFragment(
  input: Omit<Fragment, 'id' | 'length' | 'gcContent'>
): Fragment {
  return {
    ...input,
    id: generateId(),
    length: input.sequence.length,
    gcContent: calculateGcContent(input.sequence),
  };
}

function createSnapshot(state: DesignState): DesignSnapshot {
  return {
    fragments: JSON.parse(JSON.stringify(state.fragments)),
    constraints: JSON.parse(JSON.stringify(state.constraints)),
    assemblyMethod: state.assemblyMethod,
    overhangs: [...state.overhangs],
  };
}

function generateHistoryLabel(
  prevState: DesignState,
  nextState: DesignState
): string {
  const prevCount = prevState.fragments.length;
  const nextCount = nextState.fragments.length;

  if (nextCount > prevCount) {
    const added = nextState.fragments.find(
      (f) => !prevState.fragments.find((p) => p.id === f.id)
    );
    return `Added ${added?.name || 'fragment'}`;
  }

  if (nextCount < prevCount) {
    const removed = prevState.fragments.find(
      (f) => !nextState.fragments.find((n) => n.id === f.id)
    );
    return `Removed ${removed?.name || 'fragment'}`;
  }

  // Check for reorder
  const orderChanged = prevState.fragments.some(
    (f, i) => nextState.fragments[i]?.id !== f.id
  );
  if (orderChanged) {
    return 'Reordered fragments';
  }

  // Check for sequence changes
  const seqChanged = prevState.fragments.find((f) => {
    const next = nextState.fragments.find((n) => n.id === f.id);
    return next && next.sequence !== f.sequence;
  });
  if (seqChanged) {
    return `Modified ${seqChanged.name}`;
  }

  // Check for constraint changes
  if (
    JSON.stringify(prevState.constraints) !==
    JSON.stringify(nextState.constraints)
  ) {
    return 'Updated constraints';
  }

  // Check for method changes
  if (prevState.assemblyMethod !== nextState.assemblyMethod) {
    return `Switched to ${nextState.assemblyMethod.replace('_', ' ')}`;
  }

  // Check for overhang changes
  if (
    JSON.stringify(prevState.overhangs) !== JSON.stringify(nextState.overhangs)
  ) {
    return 'Updated overhangs';
  }

  return 'Modified design';
}

function shouldCreateSnapshot(
  prevState: DesignState,
  nextState: DesignState
): boolean {
  // Always snapshot fragment count changes
  if (prevState.fragments.length !== nextState.fragments.length) {
    return true;
  }

  // Snapshot sequence changes
  const seqChanged = prevState.fragments.some((f) => {
    const next = nextState.fragments.find((n) => n.id === f.id);
    return next && next.sequence !== f.sequence;
  });
  if (seqChanged) return true;

  // Snapshot reorder
  const orderChanged = prevState.fragments.some(
    (f, i) => nextState.fragments[i]?.id !== f.id
  );
  if (orderChanged) return true;

  // Snapshot method changes
  if (prevState.assemblyMethod !== nextState.assemblyMethod) {
    return true;
  }

  // Snapshot constraint changes (debounced by the UI)
  if (
    JSON.stringify(prevState.constraints) !==
    JSON.stringify(nextState.constraints)
  ) {
    return true;
  }

  return false;
}

// ============================================================================
// Initial State
// ============================================================================

const initialState: DesignState = {
  fragments: [],
  constraints: DEFAULT_CONSTRAINTS,
  assemblyMethod: 'golden_gate',
  overhangs: [],
  fidelity: 0,
  solverResult: null,
  warnings: [],
  errors: [],
  selectedFragmentId: null,
  viewMode: 'circular',
  isOptimizing: false,
  history: [],
  historyIndex: -1,
};

// ============================================================================
// Store Creation
// ============================================================================

export const useDesignStore = create<DesignStore>()(
  subscribeWithSelector(
    immer((set: (fn: (draft: DesignState) => void) => void, get: () => DesignState) => {
      // Helper to create history entry
      const pushHistory = (label: string, type: 'auto' | 'manual' = 'auto') => {
        const state = get();
        const entry: HistoryEntry = {
          id: generateHistoryId(),
          timestamp: Date.now(),
          label,
          type,
          snapshot: createSnapshot(state),
        };

        set((draft: DesignState) => {
          // If we're not at the end of history, truncate forward history
          if (draft.historyIndex < draft.history.length - 1) {
            draft.history = draft.history.slice(0, draft.historyIndex + 1);
          }
          draft.history.push(entry);
          draft.historyIndex = draft.history.length - 1;

          // Limit history size to 50 entries
          if (draft.history.length > 50) {
            draft.history.shift();
            draft.historyIndex--;
          }
        });
      };

      return {
        ...initialState,

        // ====================================================================
        // Fragment Management
        // ====================================================================

        addFragment: (input: Omit<Fragment, 'id' | 'length' | 'gcContent'>) => {
          const prevState = get();
          const fragment = createFragment(input);

          set((draft: DesignState) => {
            draft.fragments.push(fragment);
            draft.selectedFragmentId = fragment.id;
          });

          const nextState = get();
          if (shouldCreateSnapshot(prevState, nextState)) {
            pushHistory(generateHistoryLabel(prevState, nextState));
          }
        },

        updateFragment: (id: string, updates: Partial<Fragment>) => {
          const prevState = get();

          set((draft: DesignState) => {
            const index = draft.fragments.findIndex((f: Fragment) => f.id === id);
            if (index !== -1) {
              const fragment = draft.fragments[index];
              Object.assign(fragment, updates);

              // Recalculate computed fields if sequence changed
              if (updates.sequence !== undefined) {
                fragment.length = updates.sequence.length;
                fragment.gcContent = calculateGcContent(updates.sequence);
              }
            }
          });

          const nextState = get();
          if (shouldCreateSnapshot(prevState, nextState)) {
            pushHistory(generateHistoryLabel(prevState, nextState));
          }
        },

        removeFragment: (id: string) => {
          const prevState = get();

          set((draft: DesignState) => {
            const index = draft.fragments.findIndex((f: Fragment) => f.id === id);
            if (index !== -1) {
              draft.fragments.splice(index, 1);
              if (draft.selectedFragmentId === id) {
                draft.selectedFragmentId = null;
              }
            }
          });

          const nextState = get();
          if (shouldCreateSnapshot(prevState, nextState)) {
            pushHistory(generateHistoryLabel(prevState, nextState));
          }
        },

        reorderFragments: (fromIndex: number, toIndex: number) => {
          const prevState = get();

          set((draft: DesignState) => {
            const [fragment] = draft.fragments.splice(fromIndex, 1);
            draft.fragments.splice(toIndex, 0, fragment);
          });

          const nextState = get();
          if (shouldCreateSnapshot(prevState, nextState)) {
            pushHistory(generateHistoryLabel(prevState, nextState));
          }
        },

        setFragments: (fragments: Fragment[]) => {
          const prevState = get();

          set((draft: DesignState) => {
            draft.fragments = fragments.map((f: Fragment) => ({
              ...f,
              length: f.sequence.length,
              gcContent: calculateGcContent(f.sequence),
            }));
          });

          const nextState = get();
          if (shouldCreateSnapshot(prevState, nextState)) {
            pushHistory('Loaded fragments');
          }
        },

        // ====================================================================
        // Constraint Management
        // ====================================================================

        updateConstraints: (updates: Partial<AssemblyConstraints>) => {
          const prevState = get();

          set((draft: DesignState) => {
            Object.assign(draft.constraints, updates);
          });

          const nextState = get();
          if (shouldCreateSnapshot(prevState, nextState)) {
            pushHistory(generateHistoryLabel(prevState, nextState));
          }
        },

        // ====================================================================
        // Method Selection
        // ====================================================================

        setAssemblyMethod: (method: AssemblyMethod) => {
          const prevState = get();

          set((draft: DesignState) => {
            draft.assemblyMethod = method;
          });

          const nextState = get();
          if (shouldCreateSnapshot(prevState, nextState)) {
            pushHistory(generateHistoryLabel(prevState, nextState));
          }
        },

        // ====================================================================
        // Overhang Management
        // ====================================================================

        setOverhangs: (overhangs: string[]) => {
          set((draft: DesignState) => {
            draft.overhangs = overhangs;
          });
        },

        optimizeOverhangs: async () => {
          set((draft: DesignState) => {
            draft.isOptimizing = true;
            draft.errors = [];
          });

          const state = get();

          try {
            const result = await solveOverhangs(
              state.fragments.length + 1, // +1 for circular assembly junction
              state.constraints
            );

            set((draft: DesignState) => {
              draft.solverResult = result;
              draft.overhangs = result.overhangs;
              draft.fidelity = result.fidelity;
              draft.isOptimizing = false;

              if (!result.success) {
                draft.warnings = result.suggestions || [
                  'Could not find optimal overhangs',
                ];
              } else {
                draft.warnings = result.crossLigationRisks
                  .filter((r) => r.severity !== 'low')
                  .map(
                    (r) =>
                      `Cross-ligation risk: ${r.overhang1} / ${r.overhang2} (${(r.frequency * 100).toFixed(1)}%)`
                  );
              }
            });

            return result;
          } catch (error) {
            const message =
              error instanceof Error ? error.message : 'Optimization failed';

            set((draft: DesignState) => {
              draft.isOptimizing = false;
              draft.errors = [message];
            });

            return {
              success: false,
              overhangs: [],
              fidelity: 0,
              crossLigationRisks: [],
              suggestions: [message],
            } as SolverResult;
          }
        },

        // ====================================================================
        // History Navigation
        // ====================================================================

        undo: () => {
          const state = get();
          if (state.historyIndex <= 0) return;

          const prevEntry = state.history[state.historyIndex - 1];
          if (!prevEntry) return;

          set((draft: DesignState) => {
            draft.historyIndex--;
            draft.fragments = JSON.parse(
              JSON.stringify(prevEntry.snapshot.fragments)
            );
            draft.constraints = JSON.parse(
              JSON.stringify(prevEntry.snapshot.constraints)
            );
            draft.assemblyMethod = prevEntry.snapshot.assemblyMethod;
            draft.overhangs = [...prevEntry.snapshot.overhangs];
          });
        },

        redo: () => {
          const state = get();
          if (state.historyIndex >= state.history.length - 1) return;

          const nextEntry = state.history[state.historyIndex + 1];
          if (!nextEntry) return;

          set((draft: DesignState) => {
            draft.historyIndex++;
            draft.fragments = JSON.parse(
              JSON.stringify(nextEntry.snapshot.fragments)
            );
            draft.constraints = JSON.parse(
              JSON.stringify(nextEntry.snapshot.constraints)
            );
            draft.assemblyMethod = nextEntry.snapshot.assemblyMethod;
            draft.overhangs = [...nextEntry.snapshot.overhangs];
          });
        },

        canUndo: () => {
          const state = get();
          return state.historyIndex > 0;
        },

        canRedo: () => {
          const state = get();
          return state.historyIndex < state.history.length - 1;
        },

        createSnapshot: (label: string) => {
          pushHistory(label, 'manual');
        },

        restoreSnapshot: (historyId: string) => {
          const state = get();
          const index = state.history.findIndex((h: HistoryEntry) => h.id === historyId);
          if (index === -1) return;

          const entry = state.history[index];

          set((draft: DesignState) => {
            draft.historyIndex = index;
            draft.fragments = JSON.parse(
              JSON.stringify(entry.snapshot.fragments)
            );
            draft.constraints = JSON.parse(
              JSON.stringify(entry.snapshot.constraints)
            );
            draft.assemblyMethod = entry.snapshot.assemblyMethod;
            draft.overhangs = [...entry.snapshot.overhangs];
          });

          // Push a new entry for the restore action
          pushHistory(`Restored: ${entry.label}`, 'manual');
        },

        clearHistory: () => {
          set((draft: DesignState) => {
            draft.history = [];
            draft.historyIndex = -1;
          });
        },

        // ====================================================================
        // UI State
        // ====================================================================

        selectFragment: (id: string | null) => {
          set((draft: DesignState) => {
            draft.selectedFragmentId = id;
          });
        },

        setViewMode: (mode: 'circular' | 'linear') => {
          set((draft: DesignState) => {
            draft.viewMode = mode;
          });
        },

        // ====================================================================
        // Import/Export
        // ====================================================================

        exportDesign: () => {
          return createSnapshot(get());
        },

        importDesign: (snapshot: DesignSnapshot) => {
          const prevState = get();

          set((draft: DesignState) => {
            draft.fragments = snapshot.fragments.map((f: Fragment) => ({
              ...f,
              length: f.sequence.length,
              gcContent: calculateGcContent(f.sequence),
            }));
            draft.constraints = snapshot.constraints;
            draft.assemblyMethod = snapshot.assemblyMethod;
            draft.overhangs = snapshot.overhangs;
          });

          pushHistory('Imported design', 'manual');
        },

        reset: () => {
          set((draft: DesignState) => {
            Object.assign(draft, initialState);
          });
        },
      };
    })
  )
);

// ============================================================================
// Keyboard Shortcuts Hook
// ============================================================================

export function useHistoryKeyboardShortcuts() {
  const { undo, redo, canUndo, canRedo } = useDesignStore();

  if (typeof window === 'undefined') return;

  const handleKeyDown = (e: KeyboardEvent) => {
    const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
    const modKey = isMac ? e.metaKey : e.ctrlKey;

    if (modKey && e.key === 'z' && !e.shiftKey) {
      e.preventDefault();
      if (canUndo()) undo();
    }

    if (modKey && e.key === 'z' && e.shiftKey) {
      e.preventDefault();
      if (canRedo()) redo();
    }

    // Also support Ctrl+Y for redo on Windows
    if (!isMac && e.ctrlKey && e.key === 'y') {
      e.preventDefault();
      if (canRedo()) redo();
    }
  };

  window.addEventListener('keydown', handleKeyDown);
  return () => window.removeEventListener('keydown', handleKeyDown);
}

// ============================================================================
// Selectors (for performance optimization)
// ============================================================================

export const selectFragments = (state: DesignStore) => state.fragments;
export const selectConstraints = (state: DesignStore) => state.constraints;
export const selectOverhangs = (state: DesignStore) => state.overhangs;
export const selectFidelity = (state: DesignStore) => state.fidelity;
export const selectHistory = (state: DesignStore) => state.history;
export const selectHistoryIndex = (state: DesignStore) => state.historyIndex;
export const selectIsOptimizing = (state: DesignStore) => state.isOptimizing;
export const selectViewMode = (state: DesignStore) => state.viewMode;
