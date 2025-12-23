/**
 * Fragment List
 * Drag-and-drop fragment management with compact card view
 */

import React, { useState, useCallback } from 'react';
import { useDesignStore } from '../../stores/designStore';
import { Fragment, FRAGMENT_TYPES, FragmentType } from '../../types/fragmentPlanner';

export function FragmentList() {
  const {
    fragments,
    selectedFragmentId,
    addFragment,
    updateFragment,
    removeFragment,
    reorderFragments,
    selectFragment,
  } = useDesignStore();

  const [draggedIndex, setDraggedIndex] = useState<number | null>(null);
  const [dragOverIndex, setDragOverIndex] = useState<number | null>(null);
  const [showAddForm, setShowAddForm] = useState(false);

  // Drag handlers
  const handleDragStart = useCallback((e: React.DragEvent, index: number) => {
    setDraggedIndex(index);
    e.dataTransfer.effectAllowed = 'move';
    e.dataTransfer.setData('text/plain', index.toString());
  }, []);

  const handleDragOver = useCallback((e: React.DragEvent, index: number) => {
    e.preventDefault();
    if (draggedIndex !== null && draggedIndex !== index) {
      setDragOverIndex(index);
    }
  }, [draggedIndex]);

  const handleDragLeave = useCallback(() => {
    setDragOverIndex(null);
  }, []);

  const handleDrop = useCallback(
    (e: React.DragEvent, toIndex: number) => {
      e.preventDefault();
      if (draggedIndex !== null && draggedIndex !== toIndex) {
        reorderFragments(draggedIndex, toIndex);
      }
      setDraggedIndex(null);
      setDragOverIndex(null);
    },
    [draggedIndex, reorderFragments]
  );

  const handleDragEnd = useCallback(() => {
    setDraggedIndex(null);
    setDragOverIndex(null);
  }, []);

  return (
    <div className="p-2">
      {/* Fragment Cards */}
      <div className="space-y-2">
        {fragments.map((fragment: Fragment, index: number) => (
          <FragmentCard
            key={fragment.id}
            fragment={fragment}
            index={index}
            isSelected={selectedFragmentId === fragment.id}
            isDragging={draggedIndex === index}
            isDragOver={dragOverIndex === index}
            onSelect={() => selectFragment(fragment.id)}
            onUpdate={(updates) => updateFragment(fragment.id, updates)}
            onRemove={() => removeFragment(fragment.id)}
            onDragStart={(e) => handleDragStart(e, index)}
            onDragOver={(e) => handleDragOver(e, index)}
            onDragLeave={handleDragLeave}
            onDrop={(e) => handleDrop(e, index)}
            onDragEnd={handleDragEnd}
          />
        ))}
      </div>

      {/* Add Fragment Button/Form */}
      {showAddForm ? (
        <AddFragmentForm
          onAdd={(fragment) => {
            addFragment(fragment);
            setShowAddForm(false);
          }}
          onCancel={() => setShowAddForm(false)}
        />
      ) : (
        <button
          onClick={() => setShowAddForm(true)}
          className="w-full mt-3 p-3 border-2 border-dashed border-slate-200 dark:border-slate-700 rounded-lg text-slate-500 dark:text-slate-400 hover:border-amber-400 hover:text-amber-600 dark:hover:border-amber-500 dark:hover:text-amber-400 transition-colors flex items-center justify-center gap-2"
        >
          <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
          </svg>
          <span className="font-medium">Add Fragment</span>
        </button>
      )}
    </div>
  );
}

interface FragmentCardProps {
  fragment: Fragment;
  index: number;
  isSelected: boolean;
  isDragging: boolean;
  isDragOver: boolean;
  onSelect: () => void;
  onUpdate: (updates: Partial<Fragment>) => void;
  onRemove: () => void;
  onDragStart: (e: React.DragEvent) => void;
  onDragOver: (e: React.DragEvent) => void;
  onDragLeave: () => void;
  onDrop: (e: React.DragEvent) => void;
  onDragEnd: () => void;
}

function FragmentCard({
  fragment,
  index,
  isSelected,
  isDragging,
  isDragOver,
  onSelect,
  onUpdate,
  onRemove,
  onDragStart,
  onDragOver,
  onDragLeave,
  onDrop,
  onDragEnd,
}: FragmentCardProps) {
  const [isExpanded, setIsExpanded] = useState(false);
  const [isEditing, setIsEditing] = useState(false);

  const typeInfo = FRAGMENT_TYPES[fragment.type];

  return (
    <div
      draggable
      onDragStart={onDragStart}
      onDragOver={onDragOver}
      onDragLeave={onDragLeave}
      onDrop={onDrop}
      onDragEnd={onDragEnd}
      onClick={onSelect}
      className={`group relative rounded-lg border transition-all cursor-pointer ${
        isDragging
          ? 'opacity-50 border-dashed border-slate-400 dark:border-slate-500'
          : isDragOver
            ? 'border-amber-400 dark:border-amber-500 bg-amber-50 dark:bg-amber-900/20'
            : isSelected
              ? 'border-amber-500 dark:border-amber-400 bg-amber-50/50 dark:bg-amber-900/10 shadow-sm'
              : 'border-slate-200 dark:border-slate-700 hover:border-slate-300 dark:hover:border-slate-600 bg-white dark:bg-slate-800'
      }`}
    >
      {/* Main Row */}
      <div className="flex items-center gap-2 p-2">
        {/* Drag Handle */}
        <div className="flex-shrink-0 cursor-grab active:cursor-grabbing text-slate-400 dark:text-slate-500 hover:text-slate-600 dark:hover:text-slate-300">
          <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 24 24">
            <path d="M11 18c0 1.1-.9 2-2 2s-2-.9-2-2 .9-2 2-2 2 .9 2 2zm-2-8c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm0-6c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm6 4c1.1 0 2-.9 2-2s-.9-2-2-2-2 .9-2 2 .9 2 2 2zm0 2c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm0 6c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2z" />
          </svg>
        </div>

        {/* Index Badge */}
        <div className="flex-shrink-0 w-6 h-6 rounded-full bg-slate-100 dark:bg-slate-700 flex items-center justify-center text-xs font-medium text-slate-600 dark:text-slate-300">
          {index + 1}
        </div>

        {/* Type Tag */}
        <div
          className="flex-shrink-0 w-6 h-6 rounded flex items-center justify-center text-xs font-bold text-white"
          style={{ backgroundColor: typeInfo.color }}
          title={typeInfo.name}
        >
          {typeInfo.icon}
        </div>

        {/* Name & Info */}
        <div className="flex-1 min-w-0">
          <div className="flex items-center gap-2">
            <span className="text-sm font-medium text-slate-900 dark:text-white truncate">
              {fragment.name}
            </span>
          </div>
          <div className="flex items-center gap-2 text-xs text-slate-500 dark:text-slate-400">
            <span>{fragment.length.toLocaleString()} bp</span>
            <span>|</span>
            <span>GC: {fragment.gcContent.toFixed(1)}%</span>
          </div>
        </div>

        {/* Overhangs */}
        {(fragment.leftOverhang || fragment.rightOverhang) && (
          <div className="flex-shrink-0 flex items-center gap-1 text-xs font-mono">
            {fragment.leftOverhang && (
              <span className="px-1.5 py-0.5 bg-blue-100 dark:bg-blue-900/30 text-blue-700 dark:text-blue-400 rounded">
                {fragment.leftOverhang}
              </span>
            )}
            <span className="text-slate-400">-</span>
            {fragment.rightOverhang && (
              <span className="px-1.5 py-0.5 bg-violet-100 dark:bg-violet-900/30 text-violet-700 dark:text-violet-400 rounded">
                {fragment.rightOverhang}
              </span>
            )}
          </div>
        )}

        {/* Expand/Actions */}
        <div className="flex-shrink-0 flex items-center gap-1">
          <button
            onClick={(e) => {
              e.stopPropagation();
              setIsExpanded(!isExpanded);
            }}
            className="p-1 rounded hover:bg-slate-100 dark:hover:bg-slate-700 text-slate-400 hover:text-slate-600 dark:hover:text-slate-300"
          >
            <svg
              className={`w-4 h-4 transition-transform ${isExpanded ? 'rotate-180' : ''}`}
              fill="none"
              viewBox="0 0 24 24"
              stroke="currentColor"
            >
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
            </svg>
          </button>
          <button
            onClick={(e) => {
              e.stopPropagation();
              onRemove();
            }}
            className="p-1 rounded hover:bg-red-100 dark:hover:bg-red-900/30 text-slate-400 hover:text-red-600 dark:hover:text-red-400 opacity-0 group-hover:opacity-100 transition-opacity"
            title="Remove fragment"
          >
            <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>
      </div>

      {/* Expanded View */}
      {isExpanded && (
        <div className="border-t border-slate-200 dark:border-slate-700 p-3 space-y-3">
          {/* Name Edit */}
          <div>
            <label className="block text-xs font-medium text-slate-500 dark:text-slate-400 mb-1">
              Name
            </label>
            <input
              type="text"
              value={fragment.name}
              onChange={(e) => onUpdate({ name: e.target.value })}
              className="w-full px-2 py-1.5 text-sm border border-slate-200 dark:border-slate-600 rounded-md bg-white dark:bg-slate-700 text-slate-900 dark:text-white focus:ring-2 focus:ring-amber-500 focus:border-transparent"
            />
          </div>

          {/* Type Select */}
          <div>
            <label className="block text-xs font-medium text-slate-500 dark:text-slate-400 mb-1">
              Type
            </label>
            <select
              value={fragment.type}
              onChange={(e) => onUpdate({ type: e.target.value as FragmentType })}
              className="w-full px-2 py-1.5 text-sm border border-slate-200 dark:border-slate-600 rounded-md bg-white dark:bg-slate-700 text-slate-900 dark:text-white focus:ring-2 focus:ring-amber-500 focus:border-transparent"
            >
              {Object.entries(FRAGMENT_TYPES).map(([key, info]) => (
                <option key={key} value={key}>
                  {info.name}
                </option>
              ))}
            </select>
          </div>

          {/* Sequence Preview */}
          <div>
            <label className="block text-xs font-medium text-slate-500 dark:text-slate-400 mb-1">
              Sequence ({fragment.length.toLocaleString()} bp)
            </label>
            <div className="relative">
              <textarea
                value={fragment.sequence}
                onChange={(e) => {
                  const seq = e.target.value.toUpperCase().replace(/[^ATGCNRYSWKMBDHV]/g, '');
                  onUpdate({ sequence: seq });
                }}
                className="w-full px-2 py-1.5 text-xs font-mono border border-slate-200 dark:border-slate-600 rounded-md bg-white dark:bg-slate-700 text-slate-900 dark:text-white focus:ring-2 focus:ring-amber-500 focus:border-transparent resize-none"
                rows={3}
                spellCheck={false}
              />
              <button
                onClick={() => navigator.clipboard.writeText(fragment.sequence)}
                className="absolute top-1 right-1 p-1 rounded bg-slate-100 dark:bg-slate-600 text-slate-500 dark:text-slate-300 hover:bg-slate-200 dark:hover:bg-slate-500 text-xs"
                title="Copy sequence"
              >
                <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 16H6a2 2 0 01-2-2V6a2 2 0 012-2h8a2 2 0 012 2v2m-6 12h8a2 2 0 002-2v-8a2 2 0 00-2-2h-8a2 2 0 00-2 2v8a2 2 0 002 2z" />
                </svg>
              </button>
            </div>
          </div>

          {/* Range Selection (for partial fragments) */}
          {fragment.length > 100 && (
            <div>
              <label className="block text-xs font-medium text-slate-500 dark:text-slate-400 mb-1">
                Use Range (optional)
              </label>
              <div className="flex items-center gap-2">
                <input
                  type="number"
                  min={1}
                  max={fragment.length}
                  value={fragment.rangeStart || 1}
                  onChange={(e) => onUpdate({ rangeStart: parseInt(e.target.value) || 1 })}
                  className="w-20 px-2 py-1 text-sm border border-slate-200 dark:border-slate-600 rounded-md bg-white dark:bg-slate-700 text-slate-900 dark:text-white"
                  placeholder="Start"
                />
                <span className="text-slate-400">to</span>
                <input
                  type="number"
                  min={1}
                  max={fragment.length}
                  value={fragment.rangeEnd || fragment.length}
                  onChange={(e) => onUpdate({ rangeEnd: parseInt(e.target.value) || fragment.length })}
                  className="w-20 px-2 py-1 text-sm border border-slate-200 dark:border-slate-600 rounded-md bg-white dark:bg-slate-700 text-slate-900 dark:text-white"
                  placeholder="End"
                />
                <button
                  onClick={() => onUpdate({ rangeStart: 1, rangeEnd: fragment.length })}
                  className="text-xs text-amber-600 dark:text-amber-400 hover:underline"
                >
                  Full length
                </button>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

interface AddFragmentFormProps {
  onAdd: (fragment: Omit<Fragment, 'id' | 'length' | 'gcContent'>) => void;
  onCancel: () => void;
}

function AddFragmentForm({ onAdd, onCancel }: AddFragmentFormProps) {
  const [name, setName] = useState('');
  const [sequence, setSequence] = useState('');
  const [type, setType] = useState<FragmentType>('other');

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (!name.trim() || !sequence.trim()) return;

    onAdd({
      name: name.trim(),
      sequence: sequence.toUpperCase().replace(/[^ATGCNRYSWKMBDHV]/g, ''),
      type,
    });
  };

  return (
    <form onSubmit={handleSubmit} className="mt-3 p-3 border border-slate-200 dark:border-slate-700 rounded-lg bg-slate-50 dark:bg-slate-800/50">
      <div className="space-y-3">
        <div>
          <label className="block text-xs font-medium text-slate-600 dark:text-slate-300 mb-1">
            Name
          </label>
          <input
            type="text"
            value={name}
            onChange={(e) => setName(e.target.value)}
            placeholder="e.g., GFP"
            className="w-full px-2 py-1.5 text-sm border border-slate-200 dark:border-slate-600 rounded-md bg-white dark:bg-slate-700 text-slate-900 dark:text-white focus:ring-2 focus:ring-amber-500 focus:border-transparent"
            autoFocus
          />
        </div>

        <div>
          <label className="block text-xs font-medium text-slate-600 dark:text-slate-300 mb-1">
            Type
          </label>
          <select
            value={type}
            onChange={(e) => setType(e.target.value as FragmentType)}
            className="w-full px-2 py-1.5 text-sm border border-slate-200 dark:border-slate-600 rounded-md bg-white dark:bg-slate-700 text-slate-900 dark:text-white focus:ring-2 focus:ring-amber-500 focus:border-transparent"
          >
            {Object.entries(FRAGMENT_TYPES).map(([key, info]) => (
              <option key={key} value={key}>
                {info.name}
              </option>
            ))}
          </select>
        </div>

        <div>
          <label className="block text-xs font-medium text-slate-600 dark:text-slate-300 mb-1">
            Sequence
          </label>
          <textarea
            value={sequence}
            onChange={(e) => setSequence(e.target.value)}
            placeholder="Paste DNA sequence (ATGC...)"
            className="w-full px-2 py-1.5 text-sm font-mono border border-slate-200 dark:border-slate-600 rounded-md bg-white dark:bg-slate-700 text-slate-900 dark:text-white focus:ring-2 focus:ring-amber-500 focus:border-transparent resize-none"
            rows={4}
            spellCheck={false}
          />
          {sequence && (
            <p className="text-xs text-slate-500 dark:text-slate-400 mt-1">
              {sequence.replace(/[^ATGCNRYSWKMBDHV]/gi, '').length} bp
            </p>
          )}
        </div>
      </div>

      <div className="flex items-center justify-end gap-2 mt-4">
        <button
          type="button"
          onClick={onCancel}
          className="px-3 py-1.5 text-sm text-slate-600 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 rounded-md transition-colors"
        >
          Cancel
        </button>
        <button
          type="submit"
          disabled={!name.trim() || !sequence.trim()}
          className="px-3 py-1.5 text-sm bg-amber-500 hover:bg-amber-600 text-white font-medium rounded-md disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
        >
          Add Fragment
        </button>
      </div>
    </form>
  );
}

export default FragmentList;
