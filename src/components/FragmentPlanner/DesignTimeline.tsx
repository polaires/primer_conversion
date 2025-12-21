/**
 * Design Timeline
 * Shows design history with visual timeline and restore functionality
 */

import { useDesignStore } from '../../stores/designStore';
import { HistoryEntry } from '../../types/fragmentPlanner';

export function DesignTimeline() {
  const { history, historyIndex, restoreSnapshot, createSnapshot } = useDesignStore();

  const formatTime = (timestamp: number): string => {
    const now = Date.now();
    const diff = now - timestamp;

    if (diff < 60000) {
      return 'Just now';
    }
    if (diff < 3600000) {
      const mins = Math.floor(diff / 60000);
      return `${mins} min ago`;
    }
    if (diff < 86400000) {
      const hours = Math.floor(diff / 3600000);
      return `${hours} hr ago`;
    }

    return new Date(timestamp).toLocaleDateString();
  };

  const getIcon = (entry: HistoryEntry) => {
    if (entry.label.includes('Added')) {
      return (
        <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
        </svg>
      );
    }
    if (entry.label.includes('Removed')) {
      return (
        <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M20 12H4" />
        </svg>
      );
    }
    if (entry.label.includes('Reordered')) {
      return (
        <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 16V4m0 0L3 8m4-4l4 4m6 0v12m0 0l4-4m-4 4l-4-4" />
        </svg>
      );
    }
    if (entry.label.includes('constraint') || entry.label.includes('Constraint')) {
      return (
        <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10.325 4.317c.426-1.756 2.924-1.756 3.35 0a1.724 1.724 0 002.573 1.066c1.543-.94 3.31.826 2.37 2.37a1.724 1.724 0 001.065 2.572c1.756.426 1.756 2.924 0 3.35a1.724 1.724 0 00-1.066 2.573c.94 1.543-.826 3.31-2.37 2.37a1.724 1.724 0 00-2.572 1.065c-.426 1.756-2.924 1.756-3.35 0a1.724 1.724 0 00-2.573-1.066c-1.543.94-3.31-.826-2.37-2.37a1.724 1.724 0 00-1.065-2.572c-1.756-.426-1.756-2.924 0-3.35a1.724 1.724 0 001.066-2.573c-.94-1.543.826-3.31 2.37-2.37.996.608 2.296.07 2.572-1.065z" />
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
        </svg>
      );
    }
    if (entry.label.includes('Restored')) {
      return (
        <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
        </svg>
      );
    }
    // Default: edit icon
    return (
      <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M11 5H6a2 2 0 00-2 2v11a2 2 0 002 2h11a2 2 0 002-2v-5m-1.414-9.414a2 2 0 112.828 2.828L11.828 15H9v-2.828l8.586-8.586z" />
      </svg>
    );
  };

  const handleSaveCheckpoint = () => {
    const label = prompt('Enter a name for this checkpoint:');
    if (label) {
      createSnapshot(label);
    }
  };

  return (
    <div className="h-full flex flex-col">
      {/* Header */}
      <div className="p-4 border-b border-slate-200 dark:border-slate-700">
        <div className="flex items-center justify-between">
          <h2 className="text-sm font-semibold text-slate-900 dark:text-white uppercase tracking-wide">
            History
          </h2>
          <button
            onClick={handleSaveCheckpoint}
            className="text-xs px-2 py-1 rounded bg-slate-100 dark:bg-slate-700 text-slate-600 dark:text-slate-300 hover:bg-slate-200 dark:hover:bg-slate-600 transition-colors"
            title="Save a named checkpoint"
          >
            + Checkpoint
          </button>
        </div>
        <p className="text-xs text-slate-500 dark:text-slate-400 mt-1">
          {history.length} entries | Click to restore
        </p>
      </div>

      {/* Timeline */}
      <div className="flex-1 overflow-y-auto">
        {history.length === 0 ? (
          <div className="p-4 text-center text-slate-500 dark:text-slate-400 text-sm">
            <p>No history yet</p>
            <p className="text-xs mt-1">Changes will appear here</p>
          </div>
        ) : (
          <div className="p-2">
            {history
              .slice()
              .reverse()
              .map((entry, reverseIdx) => {
                const idx = history.length - 1 - reverseIdx;
                const isActive = idx === historyIndex;
                const isFuture = idx > historyIndex;

                return (
                  <button
                    key={entry.id}
                    onClick={() => restoreSnapshot(entry.id)}
                    className={`w-full flex items-start gap-3 p-3 rounded-lg text-left transition-all group ${
                      isActive
                        ? 'bg-amber-50 dark:bg-amber-900/20 border border-amber-200 dark:border-amber-800'
                        : isFuture
                          ? 'opacity-50 hover:opacity-75'
                          : 'hover:bg-slate-50 dark:hover:bg-slate-700/50'
                    }`}
                  >
                    {/* Timeline connector */}
                    <div className="flex flex-col items-center">
                      <div
                        className={`w-8 h-8 rounded-full flex items-center justify-center ${
                          isActive
                            ? 'bg-amber-500 text-white'
                            : isFuture
                              ? 'bg-slate-200 dark:bg-slate-600 text-slate-400 dark:text-slate-500'
                              : 'bg-slate-100 dark:bg-slate-700 text-slate-500 dark:text-slate-400 group-hover:bg-slate-200 dark:group-hover:bg-slate-600'
                        }`}
                      >
                        {getIcon(entry)}
                      </div>
                      {reverseIdx < history.length - 1 && (
                        <div
                          className={`w-0.5 h-4 mt-1 ${
                            isActive || isFuture
                              ? 'bg-slate-200 dark:bg-slate-700'
                              : 'bg-slate-300 dark:bg-slate-600'
                          }`}
                        />
                      )}
                    </div>

                    {/* Content */}
                    <div className="flex-1 min-w-0">
                      <div className="flex items-center gap-2">
                        <span
                          className={`text-sm font-medium truncate ${
                            isActive
                              ? 'text-amber-700 dark:text-amber-400'
                              : isFuture
                                ? 'text-slate-400 dark:text-slate-500'
                                : 'text-slate-700 dark:text-slate-300'
                          }`}
                        >
                          {entry.label}
                        </span>
                        {entry.type === 'manual' && (
                          <span className="flex-shrink-0 px-1.5 py-0.5 text-[10px] font-medium bg-blue-100 dark:bg-blue-900/30 text-blue-600 dark:text-blue-400 rounded">
                            saved
                          </span>
                        )}
                        {isActive && (
                          <span className="flex-shrink-0 px-1.5 py-0.5 text-[10px] font-medium bg-amber-100 dark:bg-amber-900/30 text-amber-600 dark:text-amber-400 rounded">
                            current
                          </span>
                        )}
                      </div>
                      <div className="flex items-center gap-2 mt-0.5">
                        <span className="text-xs text-slate-500 dark:text-slate-400">
                          {formatTime(entry.timestamp)}
                        </span>
                        <span className="text-xs text-slate-400 dark:text-slate-500">
                          {entry.snapshot.fragments.length} fragments
                        </span>
                      </div>
                    </div>
                  </button>
                );
              })}
          </div>
        )}
      </div>

      {/* Footer with keyboard hints */}
      <div className="p-3 border-t border-slate-200 dark:border-slate-700 bg-slate-50 dark:bg-slate-800/50">
        <div className="flex items-center justify-center gap-4 text-xs text-slate-500 dark:text-slate-400">
          <span className="flex items-center gap-1">
            <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-[10px] font-mono">
              Ctrl
            </kbd>
            <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-[10px] font-mono">
              Z
            </kbd>
            <span>Undo</span>
          </span>
          <span className="flex items-center gap-1">
            <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-[10px] font-mono">
              Ctrl
            </kbd>
            <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-[10px] font-mono">
              Shift
            </kbd>
            <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-[10px] font-mono">
              Z
            </kbd>
            <span>Redo</span>
          </span>
        </div>
      </div>
    </div>
  );
}

export default DesignTimeline;
