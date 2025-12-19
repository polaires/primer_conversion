# Primer Repository Conversion Guide

A comprehensive guide for converting the `polaires/primer` repository from JavaScript + Custom CSS to TypeScript + Tailwind CSS, matching the patterns used in Protein_engineering_tools.

---

## Table of Contents

1. [Project Setup](#1-project-setup)
2. [Directory Structure](#2-directory-structure)
3. [TypeScript Configuration](#3-typescript-configuration)
4. [Tailwind Configuration](#4-tailwind-configuration)
5. [CSS to Tailwind Mapping](#5-css-to-tailwind-mapping)
6. [Component Conversion Patterns](#6-component-conversion-patterns)
7. [Type Definitions](#7-type-definitions)
8. [State Management Patterns](#8-state-management-patterns)
9. [Import/Export Conventions](#9-importexport-conventions)
10. [Dark Mode Implementation](#10-dark-mode-implementation)
11. [Component-by-Component Conversion Order](#11-component-by-component-conversion-order)
12. [Testing & Validation](#12-testing--validation)
13. [Common Pitfalls](#13-common-pitfalls)

---

## 1. Project Setup

### 1.1 Install Dependencies

```bash
# In primer repo root
npm install -D typescript @types/react @types/react-dom
npm install -D tailwindcss postcss autoprefixer
npm install -D @vitejs/plugin-react
npm install clsx lucide-react

# Initialize Tailwind
npx tailwindcss init -p
```

### 1.2 Update package.json

```json
{
  "name": "primers",
  "version": "1.0.0",
  "type": "module",
  "scripts": {
    "dev": "vite",
    "build": "tsc && vite build",
    "preview": "vite preview",
    "test": "vitest",
    "lint": "eslint . --ext ts,tsx --report-unused-disable-directives --max-warnings 0"
  },
  "dependencies": {
    "clsx": "^2.0.0",
    "lucide-react": "^0.294.0",
    "react": "^18.2.0",
    "react-dom": "^18.2.0",
    "seqviz": "^3.10.10"
  },
  "devDependencies": {
    "@types/react": "^18.2.43",
    "@types/react-dom": "^18.2.17",
    "@vitejs/plugin-react": "^4.2.1",
    "autoprefixer": "^10.4.16",
    "postcss": "^8.4.32",
    "tailwindcss": "^3.3.6",
    "typescript": "^5.2.2",
    "vite": "^5.0.8",
    "vitest": "^1.0.4"
  }
}
```

### 1.3 Create/Update vite.config.ts

```typescript
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import path from 'path';

export default defineConfig({
  plugins: [react()],

  server: {
    port: 3000,
  },

  build: {
    target: ['es2020', 'chrome87', 'firefox78', 'safari14', 'edge88'],
    sourcemap: true,
    chunkSizeWarningLimit: 5000,
  },

  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
});
```

---

## 2. Directory Structure

### 2.1 Target Structure (After Conversion)

```
primer/
├── src/
│   ├── main.tsx                    # Entry point
│   ├── App.tsx                     # Main app component
│   ├── components/
│   │   ├── primers/                # Primer sub-components
│   │   │   ├── AlternativesPanel.tsx
│   │   │   ├── EnhancedAnalysisSection.tsx
│   │   │   ├── ScoreBreakdownPopup.tsx
│   │   │   ├── SummaryStatusPanel.tsx
│   │   │   └── index.ts
│   │   ├── AssemblyDesigner.tsx
│   │   ├── GoldenGateDesigner.tsx
│   │   ├── TmCalculator.tsx
│   │   ├── UnifiedPrimerDesigner.tsx
│   │   ├── SequencingDesigner.tsx
│   │   ├── EnhancedScorer.tsx
│   │   ├── StandaloneViewer.tsx
│   │   └── ... (other components)
│   ├── lib/                        # Keep as JS (gradual migration)
│   │   ├── index.js                # Main exports
│   │   ├── primers.js
│   │   ├── scoring.js
│   │   ├── mutagenesis.js
│   │   ├── repp/
│   │   └── ...
│   ├── types/
│   │   └── index.ts                # All type definitions
│   ├── styles/
│   │   └── index.css               # Tailwind directives + custom
│   └── utils/                      # Utility functions (if needed)
├── public/
├── index.html
├── package.json
├── tsconfig.json
├── tailwind.config.js
├── postcss.config.js
└── vite.config.ts
```

### 2.2 File Renaming Convention

```
# JSX → TSX
App.jsx → App.tsx
TmCalculator.jsx → TmCalculator.tsx

# Keep lib as JS (can add .d.ts files for types)
lib/primers.js → lib/primers.js (keep)
lib/primers.d.ts → (create type definitions)
```

---

## 3. TypeScript Configuration

### 3.1 tsconfig.json

```json
{
  "compilerOptions": {
    "target": "ES2020",
    "useDefineForClassFields": true,
    "lib": ["ES2020", "DOM", "DOM.Iterable"],
    "module": "ESNext",
    "skipLibCheck": true,

    "moduleResolution": "bundler",
    "allowImportingTsExtensions": true,
    "resolveJsonModule": true,
    "isolatedModules": true,
    "noEmit": true,
    "jsx": "react-jsx",

    "strict": true,
    "noUnusedLocals": true,
    "noUnusedParameters": true,
    "noFallthroughCasesInSwitch": true,

    "baseUrl": ".",
    "paths": {
      "@/*": ["./src/*"]
    },

    "allowJs": true
  },
  "include": ["src"],
  "references": [{ "path": "./tsconfig.node.json" }]
}
```

### 3.2 tsconfig.node.json

```json
{
  "compilerOptions": {
    "composite": true,
    "skipLibCheck": true,
    "module": "ESNext",
    "moduleResolution": "bundler",
    "allowSyntheticDefaultImports": true
  },
  "include": ["vite.config.ts"]
}
```

---

## 4. Tailwind Configuration

### 4.1 tailwind.config.js

```javascript
/** @type {import('tailwindcss').Config} */
export default {
  darkMode: 'class',
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        // Primary colors (sky blue - matches Protein_engineering_tools)
        primary: {
          50: '#f0f9ff',
          100: '#e0f2fe',
          200: '#bae6fd',
          300: '#7dd3fc',
          400: '#38bdf8',
          500: '#0ea5e9',
          600: '#0284c7',
          700: '#0369a1',
          800: '#075985',
          900: '#0c4a6e',
        },
        // Primer-specific colors (for forward/reverse primers)
        primer: {
          fwd: '#3b82f6',      // Blue for forward
          rev: '#8b5cf6',      // Purple for reverse
          template: '#10b981', // Green for template
          mutation: '#f59e0b', // Amber for mutations
        },
        // Science accent colors
        science: {
          beaker: '#4ade80',
          flask: '#818cf8',
          lab: '#f59e0b',
        }
      },
      fontFamily: {
        sans: ['Inter', 'system-ui', 'sans-serif'],
        mono: ['JetBrains Mono', 'Monaco', 'Menlo', 'monospace'],
      },
      boxShadow: {
        'glass': '0 8px 32px 0 rgba(31, 38, 135, 0.37)',
      },
    },
  },
  plugins: [],
}
```

### 4.2 postcss.config.js

```javascript
export default {
  plugins: {
    tailwindcss: {},
    autoprefixer: {},
  },
}
```

---

## 5. CSS to Tailwind Mapping

### 5.1 Core CSS Variables → Tailwind

| Primer CSS Variable | Tailwind Equivalent |
|---------------------|---------------------|
| `--primary-color: #2563eb` | `primary-600` (we use `#0284c7`) |
| `--primary-dark: #1d4ed8` | `primary-700` |
| `--secondary-color: #10b981` | `emerald-500` or `science-beaker` |
| `--error-color: #ef4444` | `red-500` |
| `--text-color: #1f2937` | `slate-800` |
| `--text-light: #6b7280` | `slate-500` |
| `--bg-color: #f9fafb` | `slate-50` |
| `--card-bg: #ffffff` | `white` / `dark:slate-800` |
| `--border-color: #e5e7eb` | `slate-200` / `dark:slate-700` |
| `--fwd-color: #3b82f6` | `primer-fwd` or `blue-500` |
| `--rev-color: #8b5cf6` | `primer-rev` or `violet-500` |

### 5.2 Common Class Mappings

#### Layout & Spacing

| Primer CSS | Tailwind |
|------------|----------|
| `padding: 1rem` | `p-4` |
| `padding: 1.5rem` | `p-6` |
| `padding: 2rem` | `p-8` |
| `margin-bottom: 1rem` | `mb-4` |
| `gap: 0.5rem` | `gap-2` |
| `gap: 1rem` | `gap-4` |

#### Typography

| Primer CSS | Tailwind |
|------------|----------|
| `font-size: 2.5rem` | `text-4xl` |
| `font-size: 1.1rem` | `text-lg` |
| `font-size: 1rem` | `text-base` |
| `font-size: 0.9rem` | `text-sm` |
| `font-size: 0.875rem` | `text-sm` |
| `font-size: 0.8rem` | `text-xs` |
| `font-weight: 700` | `font-bold` |
| `font-weight: 600` | `font-semibold` |
| `font-weight: 500` | `font-medium` |

#### Borders & Radius

| Primer CSS | Tailwind |
|------------|----------|
| `border-radius: 8px` | `rounded-lg` |
| `border-radius: 12px` | `rounded-xl` |
| `border-radius: 6px` | `rounded-md` |
| `border: 2px solid var(--border-color)` | `border-2 border-slate-200 dark:border-slate-700` |
| `border: 1px solid` | `border` |

#### Flexbox & Grid

| Primer CSS | Tailwind |
|------------|----------|
| `display: flex` | `flex` |
| `flex-direction: column` | `flex-col` |
| `justify-content: space-between` | `justify-between` |
| `align-items: center` | `items-center` |
| `flex: 1` | `flex-1` |
| `display: grid` | `grid` |
| `grid-template-columns: 1fr 1fr` | `grid-cols-2` |
| `grid-template-columns: repeat(auto-fit, minmax(150px, 1fr))` | `grid-cols-[repeat(auto-fit,minmax(150px,1fr))]` |

### 5.3 Component Class Mappings

#### .app-header → Header Component
```css
/* Before (Primer) */
.app-header {
  background: linear-gradient(135deg, var(--primary-color), var(--primary-dark));
  color: white;
  padding: 2rem;
  text-align: center;
}

/* After (Tailwind) */
className="bg-gradient-to-br from-primary-600 to-primary-700 text-white p-8 text-center"
```

#### .primer-form → Form Component
```css
/* Before */
.primer-form {
  background: var(--card-bg);
  padding: 1.5rem;
  border-radius: 12px;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

/* After - use the card class */
className="card"
/* or explicitly: */
className="glass-card p-6"
```

#### .mode-tabs button → Tab Button
```css
/* Before */
.mode-tabs button {
  flex: 1;
  padding: 1rem;
  border: 2px solid var(--border-color);
  background: var(--card-bg);
  color: var(--text-color);
  font-size: 1rem;
  font-weight: 500;
  cursor: pointer;
  border-radius: 8px;
  transition: all 0.2s;
}

.mode-tabs button.active {
  background: var(--primary-color);
  border-color: var(--primary-color);
  color: white;
}

/* After - use calc-mode-tab class */
className={`calc-mode-tab ${isActive ? 'active' : ''}`}
```

#### .submit-btn → Primary Button
```css
/* Before */
.submit-btn {
  width: 100%;
  padding: 1rem;
  background: var(--primary-color);
  color: white;
  border: none;
  border-radius: 8px;
  font-size: 1rem;
  font-weight: 600;
  cursor: pointer;
  transition: background-color 0.2s;
}

/* After */
className="btn-primary w-full"
```

#### .form-group → Form Field
```css
/* Before */
.form-group {
  margin-bottom: 1rem;
}
.form-group label {
  display: block;
  font-weight: 500;
  margin-bottom: 0.5rem;
  color: var(--text-color);
}
.form-group input {
  width: 100%;
  padding: 0.75rem;
  border: 2px solid var(--border-color);
  border-radius: 8px;
  font-family: monospace;
}

/* After */
<div className="mb-4">
  <label className="input-label">Label</label>
  <input className="input-field font-mono" />
</div>
```

#### .results-container → Results Card
```css
/* Before */
.results-container {
  background: var(--card-bg);
  padding: 1.5rem;
  border-radius: 12px;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

/* After */
className="card"
```

#### .loading / .spinner → Loading State
```css
/* Before */
.spinner {
  width: 40px;
  height: 40px;
  border: 3px solid var(--border-color);
  border-top-color: var(--primary-color);
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

/* After - use existing spinner class */
<div className="spinner" />
/* or with size: */
<div className="spinner w-10 h-10" />
```

---

## 6. Component Conversion Patterns

### 6.1 Basic Component Template

```tsx
// Before (JSX)
import React, { useState, useCallback } from 'react';
import './ComponentStyles.css';

export default function MyComponent({ prop1, prop2 }) {
  const [value, setValue] = useState('');

  return (
    <div className="my-component">
      <h2 className="title">{prop1}</h2>
      <input
        type="text"
        value={value}
        onChange={(e) => setValue(e.target.value)}
        className="input"
      />
    </div>
  );
}

// After (TSX)
import { useState, useCallback } from 'react';

interface MyComponentProps {
  prop1: string;
  prop2?: number;
}

export default function MyComponent({ prop1, prop2 }: MyComponentProps) {
  const [value, setValue] = useState<string>('');

  return (
    <div className="card">
      <h2 className="section-title">{prop1}</h2>
      <input
        type="text"
        value={value}
        onChange={(e) => setValue(e.target.value)}
        className="input-field"
      />
    </div>
  );
}
```

### 6.2 Event Handler Typing

```tsx
// Before
const handleChange = (e) => {
  setValue(e.target.value);
};

const handleSubmit = (e) => {
  e.preventDefault();
  // ...
};

// After
const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
  setValue(e.target.value);
};

const handleSubmit = (e: React.FormEvent<HTMLFormElement>) => {
  e.preventDefault();
  // ...
};

// Other common types
const handleClick = (e: React.MouseEvent<HTMLButtonElement>) => {};
const handleKeyDown = (e: React.KeyboardEvent<HTMLInputElement>) => {};
const handleSelect = (e: React.ChangeEvent<HTMLSelectElement>) => {};
```

### 6.3 Ref Typing

```tsx
// Before
const inputRef = useRef(null);
const fileInputRef = useRef(null);

// After
const inputRef = useRef<HTMLInputElement>(null);
const fileInputRef = useRef<HTMLInputElement>(null);
const divRef = useRef<HTMLDivElement>(null);
const textareaRef = useRef<HTMLTextAreaElement>(null);
```

### 6.4 State Typing with Union Types

```tsx
// Before
const [mode, setMode] = useState('primer-designer');
const [activeTab, setActiveTab] = useState('design');

// After
type Mode = 'primer-designer' | 'sequencing' | 'score' | 'tm' | 'viewer' | 'assembly-studio' | 'fragment-planner';
type Tab = 'design' | 'analysis' | 'assembly';

const [mode, setMode] = useState<Mode>('primer-designer');
const [activeTab, setActiveTab] = useState<Tab>('design');
```

### 6.5 Callback Props Pattern

```tsx
// Before
function ChildComponent({ onSubmit, onCancel }) {
  return (
    <button onClick={() => onSubmit({ data: 'value' })}>Submit</button>
  );
}

// After
interface SubmitData {
  data: string;
  timestamp?: number;
}

interface ChildComponentProps {
  onSubmit: (data: SubmitData) => void;
  onCancel?: () => void;
}

function ChildComponent({ onSubmit, onCancel }: ChildComponentProps) {
  return (
    <button
      onClick={() => onSubmit({ data: 'value' })}
      className="btn-primary"
    >
      Submit
    </button>
  );
}
```

---

## 7. Type Definitions

### 7.1 src/types/index.ts - Core Types

```typescript
/**
 * Primer Design Types
 */

// ============================================================================
// Primer Types
// ============================================================================

export interface Primer {
  sequence: string;
  tm: number;
  gc: number;
  length: number;
  direction: 'forward' | 'reverse';
  score?: number;
  hairpin?: HairpinResult;
  homodimer?: DimerResult;
}

export interface PrimerPair {
  forward: Primer;
  reverse: Primer;
  tmDifference: number;
  heterodimer?: DimerResult;
  quality: QualityTier;
}

export type QualityTier = 'excellent' | 'good' | 'acceptable' | 'poor';

// ============================================================================
// Thermodynamics Types
// ============================================================================

export interface TmResult {
  tm: number;
  method: 'Q5' | 'SantaLucia' | 'DNA24';
  parameters?: TmParameters;
}

export interface TmParameters {
  primerConc: number;    // nM
  saltConc: number;      // mM
  mgConc: number;        // mM
  dNTPConc: number;      // mM
}

export interface HairpinResult {
  dG: number;
  tm: number;
  structure?: string;
  position?: number;
}

export interface DimerResult {
  dG: number;
  tm: number;
  position?: number;
  severity: 'none' | 'low' | 'moderate' | 'high';
}

// ============================================================================
// Assembly Types
// ============================================================================

export type AssemblyMethod = 'golden_gate' | 'gibson' | 'nebuilder_hifi';

export interface AssemblyFragment {
  id: string;
  name: string;
  sequence: string;
  length: number;
  type?: FragmentType;
}

export type FragmentType = 'promoter' | 'rbs' | 'cds' | 'terminator' | 'other';

export interface GoldenGateEnzyme {
  name: string;
  recognition: string;
  overhang: number;
  temperature: number;
}

export interface AssemblyResult {
  success: boolean;
  fragments: AssemblyFragment[];
  primers: PrimerPair[];
  overhangs?: string[];
  fidelity?: number;
  warnings?: string[];
}

// ============================================================================
// Mutagenesis Types
// ============================================================================

export type MutationType = 'substitution' | 'insertion' | 'deletion' | 'codon_change';

export interface MutationSpec {
  type: MutationType;
  position: number;
  original?: string;
  replacement?: string;
  codon?: string;
  aminoAcid?: string;
}

export interface MutagenesisResult {
  primers: PrimerPair;
  mutantSequence: string;
  mutationType: MutationType;
  verified: boolean;
}

// ============================================================================
// Sequencing Types
// ============================================================================

export interface SequencingPrimer extends Primer {
  position: number;
  coverage: [number, number];  // Start and end positions covered
}

export interface SequencingDesign {
  primers: SequencingPrimer[];
  coverage: number;  // Percentage
  gaps?: [number, number][];
}

// ============================================================================
// UI State Types
// ============================================================================

export type ToolMode =
  | 'primer-designer'
  | 'sequencing'
  | 'score'
  | 'tm'
  | 'viewer'
  | 'assembly-studio'
  | 'fragment-planner';

export type ToolCategory = 'design' | 'analysis' | 'assembly';

export interface Tool {
  id: ToolMode;
  label: string;
  category: ToolCategory;
  description: string;
}
```

### 7.2 Library Type Declarations (src/lib/primers.d.ts)

```typescript
/**
 * Type declarations for primers.js library
 */

import type { Primer, PrimerPair, TmParameters } from '../types';

export const LEN_MIN: number;
export const LEN_MAX: number;
export const LEN_MAX_EXTENDED: number;

export interface PrimerOptions {
  minTm?: number;
  maxTm?: number;
  optimalTm?: number;
  minGC?: number;
  maxGC?: number;
  minLength?: number;
  maxLength?: number;
  circular?: boolean;
  offtargetCheck?: boolean;
  tmParams?: TmParameters;
}

export function primers(
  sequence: string,
  options?: PrimerOptions
): [Primer, Primer];

export function create(
  sequence: string,
  options?: PrimerOptions
): [Primer, Primer];

export function score(
  fwdPrimer: string,
  revPrimer: string,
  sequence: string,
  offtargetCheck?: boolean,
  options?: PrimerOptions
): [Primer, Primer];

export function generateAlternatives(
  sequence: string,
  options?: PrimerOptions
): PrimerPair[];
```

---

## 8. State Management Patterns

### 8.1 Local Component State

```tsx
// Simple state
const [value, setValue] = useState<string>('');
const [count, setCount] = useState<number>(0);
const [isOpen, setIsOpen] = useState<boolean>(false);

// Complex state with interface
interface FormState {
  sequence: string;
  options: PrimerOptions;
  results: PrimerPair | null;
}

const [formState, setFormState] = useState<FormState>({
  sequence: '',
  options: { minTm: 55, maxTm: 65 },
  results: null,
});

// Update nested state
setFormState(prev => ({
  ...prev,
  options: { ...prev.options, minTm: 60 }
}));
```

### 8.2 Reducer Pattern for Complex State

```tsx
interface DesignerState {
  template: string;
  mode: 'amplify' | 'mutagenesis';
  results: PrimerPair | null;
  loading: boolean;
  error: string | null;
}

type DesignerAction =
  | { type: 'SET_TEMPLATE'; payload: string }
  | { type: 'SET_MODE'; payload: 'amplify' | 'mutagenesis' }
  | { type: 'SET_LOADING'; payload: boolean }
  | { type: 'SET_RESULTS'; payload: PrimerPair }
  | { type: 'SET_ERROR'; payload: string }
  | { type: 'RESET' };

function designerReducer(state: DesignerState, action: DesignerAction): DesignerState {
  switch (action.type) {
    case 'SET_TEMPLATE':
      return { ...state, template: action.payload };
    case 'SET_MODE':
      return { ...state, mode: action.payload };
    case 'SET_LOADING':
      return { ...state, loading: action.payload, error: null };
    case 'SET_RESULTS':
      return { ...state, results: action.payload, loading: false };
    case 'SET_ERROR':
      return { ...state, error: action.payload, loading: false };
    case 'RESET':
      return initialState;
    default:
      return state;
  }
}

// Usage
const [state, dispatch] = useReducer(designerReducer, initialState);
```

### 8.3 Computed Values with useMemo

```tsx
// Before (computed in render)
const gcContent = sequence.split('').filter(b => b === 'G' || b === 'C').length / sequence.length * 100;

// After (memoized)
const gcContent = useMemo(() => {
  if (!sequence) return 0;
  const gcCount = sequence.split('').filter(b => b === 'G' || b === 'C').length;
  return (gcCount / sequence.length) * 100;
}, [sequence]);
```

---

## 9. Import/Export Conventions

### 9.1 Component Imports

```tsx
// React hooks (named imports, no React prefix needed in React 18)
import { useState, useEffect, useMemo, useCallback, useRef } from 'react';

// Icons from lucide-react
import {
  Dna,
  FlaskConical,
  Play,
  Copy,
  Download,
  AlertCircle,
  CheckCircle,
  Info
} from 'lucide-react';

// Internal components (absolute paths with @/)
import TmCalculator from '@/components/TmCalculator';
import { AlternativesPanel, EnhancedAnalysisSection } from '@/components/primers';

// Types
import type { Primer, PrimerPair, ToolMode } from '@/types';

// Library functions (keep as JS imports)
import { primers, score, generateAlternatives } from '@/lib/primers.js';
import { calculateTmQ5, calculateGC } from '@/lib/tmQ5.js';

// Utilities
import { clsx } from 'clsx';
```

### 9.2 Component Exports

```tsx
// Default export for main components
export default function TmCalculator() { ... }

// Named exports for sub-components
export function PrimerCard({ primer }: PrimerCardProps) { ... }
export function ResultsTable({ results }: ResultsTableProps) { ... }

// Barrel exports in index.ts
// src/components/primers/index.ts
export { default as AlternativesPanel } from './AlternativesPanel';
export { default as EnhancedAnalysisSection } from './EnhancedAnalysisSection';
export { default as ScoreBreakdownPopup } from './ScoreBreakdownPopup';
export { default as SummaryStatusPanel } from './SummaryStatusPanel';
```

---

## 10. Dark Mode Implementation

### 10.1 Pattern for All Colors

Always provide both light and dark variants:

```tsx
// Text colors
className="text-slate-900 dark:text-slate-100"     // Primary text
className="text-slate-600 dark:text-slate-400"     // Secondary text
className="text-slate-500 dark:text-slate-500"     // Muted text

// Background colors
className="bg-white dark:bg-slate-800"              // Card background
className="bg-slate-50 dark:bg-slate-900"           // Page background
className="bg-slate-100 dark:bg-slate-700"          // Subtle background

// Border colors
className="border-slate-200 dark:border-slate-700"  // Standard border
className="border-slate-300 dark:border-slate-600"  // Input border

// Interactive elements
className="hover:bg-slate-100 dark:hover:bg-slate-700"
className="focus:ring-primary-500 dark:focus:ring-primary-400"
```

### 10.2 Primer-Specific Colors with Dark Mode

```tsx
// Forward primer (blue)
className="text-primer-fwd dark:text-blue-400"
className="bg-blue-50 dark:bg-blue-900/20"

// Reverse primer (purple)
className="text-primer-rev dark:text-violet-400"
className="bg-violet-50 dark:bg-violet-900/20"

// Template (green)
className="text-primer-template dark:text-emerald-400"

// Mutations (amber)
className="text-primer-mutation dark:text-amber-400"
className="bg-amber-50 dark:bg-amber-900/20"
```

### 10.3 Quality Indicators

```tsx
// Excellent (green)
className="text-green-600 dark:text-green-400"
className="bg-green-100 dark:bg-green-900/30"

// Good (blue)
className="text-blue-600 dark:text-blue-400"
className="bg-blue-100 dark:bg-blue-900/30"

// Warning (amber)
className="text-amber-600 dark:text-amber-400"
className="bg-amber-100 dark:bg-amber-900/30"

// Error/Poor (red)
className="text-red-600 dark:text-red-400"
className="bg-red-100 dark:bg-red-900/30"
```

---

## 11. Component-by-Component Conversion Order

### Priority Order (Recommended)

Convert in this order based on dependencies and complexity:

#### Phase 1: Foundation (Days 1-2)
1. **main.tsx** - Entry point
2. **App.tsx** - Main app shell, navigation
3. **Types** - Create all type definitions

#### Phase 2: Simple Components (Days 3-5)
4. **TmCalculator.tsx** - Simple, standalone
5. **PrimerForm.tsx** - Basic form component
6. **PrimerResults.tsx** - Results display
7. **StandaloneViewer.tsx** - Sequence viewer wrapper

#### Phase 3: Medium Components (Days 6-9)
8. **EnhancedScorer.tsx** - Primer scoring
9. **SequenceViewer.tsx** - Core visualization
10. **HairpinDiagram.tsx** - Structure visualization
11. **FornaViewer.tsx** - RNA structure
12. **SecondaryStructureViewer.tsx**
13. **primers/SummaryStatusPanel.tsx**
14. **primers/ScoreBreakdownPopup.tsx**

#### Phase 4: Complex Components (Days 10-14)
15. **SequencingDesigner.tsx** - Sequencing tool
16. **PrimerStructureViewer.tsx** - Detailed structure
17. **PrimerOnTemplateViewer.tsx** - Template viz
18. **primers/AlternativesPanel.tsx**
19. **primers/EnhancedAnalysisSection.tsx**
20. **UnifiedPrimerDesigner.tsx** - Main designer

#### Phase 5: Assembly Components (Days 15-18)
21. **AssemblyDesigner.tsx** - Fragment planner
22. **IsothermalAssemblyPanel.tsx** - Gibson/HiFi
23. **CrossLigationHeatmap.tsx**
24. **SequenceConflictMap.tsx**
25. **FusionSiteOptimizerPanel.tsx**
26. **EnhancedDomesticationPanel.tsx**
27. **DomesticationWorkflowGuide.tsx**
28. **GoldenGateDesigner.tsx** - Most complex

---

## 12. Testing & Validation

### 12.1 Visual Comparison Checklist

For each converted component, verify:

- [ ] Colors match original (light mode)
- [ ] Colors work in dark mode
- [ ] Spacing/padding matches
- [ ] Font sizes match
- [ ] Border radius matches
- [ ] Shadows match
- [ ] Hover states work
- [ ] Focus states work
- [ ] Transitions are smooth
- [ ] Responsive breakpoints work
- [ ] Loading states match
- [ ] Error states match

### 12.2 Functional Testing

- [ ] All calculations produce same results
- [ ] Form validation works
- [ ] File upload/download works
- [ ] Copy to clipboard works
- [ ] All tooltips display
- [ ] Keyboard navigation works
- [ ] Screen reader compatibility

### 12.3 Screenshot Comparison Tool

```bash
# Install comparison tools
npm install -D puppeteer pixelmatch pngjs

# Create comparison script (scripts/compare-screenshots.js)
```

```javascript
// scripts/compare-screenshots.js
const puppeteer = require('puppeteer');
const { PNG } = require('pngjs');
const pixelmatch = require('pixelmatch');
const fs = require('fs');

async function comparePages(url1, url2, outputPath) {
  const browser = await puppeteer.launch();

  // Screenshot original
  const page1 = await browser.newPage();
  await page1.goto(url1);
  await page1.screenshot({ path: 'original.png', fullPage: true });

  // Screenshot converted
  const page2 = await browser.newPage();
  await page2.goto(url2);
  await page2.screenshot({ path: 'converted.png', fullPage: true });

  // Compare
  const img1 = PNG.sync.read(fs.readFileSync('original.png'));
  const img2 = PNG.sync.read(fs.readFileSync('converted.png'));
  const diff = new PNG({ width: img1.width, height: img1.height });

  const numDiffPixels = pixelmatch(
    img1.data, img2.data, diff.data,
    img1.width, img1.height,
    { threshold: 0.1 }
  );

  fs.writeFileSync(outputPath, PNG.sync.write(diff));
  console.log(`Diff pixels: ${numDiffPixels}`);

  await browser.close();
}

comparePages(
  'http://localhost:3000',  // Original
  'http://localhost:3001',  // Converted
  'diff.png'
);
```

---

## 13. Common Pitfalls

### 13.1 TypeScript Errors

```tsx
// ❌ Error: Property 'value' does not exist on type 'EventTarget'
const handleChange = (e) => setValue(e.target.value);

// ✅ Fix: Type the event
const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
  setValue(e.target.value);
};

// ❌ Error: Object is possibly 'null'
const element = document.getElementById('root');
element.innerHTML = ''; // Error!

// ✅ Fix: Add null check or assertion
const element = document.getElementById('root');
if (element) element.innerHTML = '';
// or
const element = document.getElementById('root')!;
```

### 13.2 Dark Mode Forgetting

```tsx
// ❌ Missing dark mode
className="bg-white text-black"

// ✅ With dark mode
className="bg-white dark:bg-slate-800 text-slate-900 dark:text-slate-100"
```

### 13.3 Tailwind Class Conflicts

```tsx
// ❌ Conflicting classes (last wins, unpredictable)
className="p-4 p-6"

// ✅ Use clsx for conditional classes
import { clsx } from 'clsx';
className={clsx('p-4', isLarge && 'p-6')}
```

### 13.4 Missing Transition Classes

```tsx
// ❌ Abrupt state changes
className="bg-white hover:bg-slate-100"

// ✅ Smooth transitions
className="bg-white hover:bg-slate-100 transition-colors duration-200"
```

### 13.5 Z-Index Issues

```tsx
// Use consistent z-index scale:
// z-10: Dropdowns, tooltips
// z-20: Modals backdrop
// z-30: Modal content
// z-40: Header/navigation
// z-50: Toast notifications
```

### 13.6 Font Mono for Sequences

```tsx
// ❌ Default font for DNA sequences
<span>{sequence}</span>

// ✅ Monospace for sequences
<span className="font-mono">{sequence}</span>
```

---

## Appendix A: Custom CSS Classes to Add

Add these to `src/styles/index.css`:

```css
@tailwind base;
@tailwind components;
@tailwind utilities;

@layer components {
  /* Primer-specific sequence display */
  .sequence-display {
    @apply font-mono text-sm tracking-wide;
    @apply bg-slate-50 dark:bg-slate-900;
    @apply p-3 rounded-lg;
    @apply overflow-x-auto;
  }

  /* Primer forward/reverse indicators */
  .primer-fwd-bg {
    @apply bg-blue-50 dark:bg-blue-900/20;
    @apply border-l-4 border-blue-500;
  }

  .primer-rev-bg {
    @apply bg-violet-50 dark:bg-violet-900/20;
    @apply border-l-4 border-violet-500;
  }

  /* Score indicators */
  .score-excellent {
    @apply text-green-600 dark:text-green-400;
  }

  .score-good {
    @apply text-blue-600 dark:text-blue-400;
  }

  .score-acceptable {
    @apply text-amber-600 dark:text-amber-400;
  }

  .score-poor {
    @apply text-red-600 dark:text-red-400;
  }

  /* Tm display */
  .tm-value {
    @apply font-mono font-semibold;
    @apply text-lg;
  }

  /* Tool navigation sections */
  .tool-nav-section {
    @apply flex items-center gap-2;
  }

  .tool-nav-label {
    @apply text-xs font-medium uppercase tracking-wider;
    @apply text-slate-500 dark:text-slate-400;
  }
}
```

---

## Appendix B: Checklist for Each Component

Use this checklist when converting each component:

```markdown
## Component: [ComponentName]

### Setup
- [ ] Rename .jsx to .tsx
- [ ] Add TypeScript interface for props
- [ ] Type all useState hooks
- [ ] Type all useRef hooks
- [ ] Type all event handlers

### Styling
- [ ] Remove CSS import
- [ ] Convert all className strings to Tailwind
- [ ] Add dark mode variants
- [ ] Verify responsive breakpoints
- [ ] Add transition classes

### Testing
- [ ] Component renders without errors
- [ ] All interactive elements work
- [ ] Visual match with original (light mode)
- [ ] Visual match with original (dark mode)
- [ ] TypeScript compiles without errors

### Documentation
- [ ] Add JSDoc comments for complex functions
- [ ] Document any props
```

---

## Quick Reference Card

### Most Common Conversions

| Primer CSS | Tailwind |
|------------|----------|
| `background: var(--card-bg)` | `bg-white dark:bg-slate-800` |
| `color: var(--text-color)` | `text-slate-900 dark:text-slate-100` |
| `color: var(--text-light)` | `text-slate-500 dark:text-slate-400` |
| `border: 1px solid var(--border-color)` | `border border-slate-200 dark:border-slate-700` |
| `border-radius: 8px` | `rounded-lg` |
| `padding: 1rem` | `p-4` |
| `margin-bottom: 1rem` | `mb-4` |
| `font-weight: 600` | `font-semibold` |
| `font-size: 0.875rem` | `text-sm` |
| `transition: all 0.2s` | `transition-all duration-200` |

### Common Component Classes

| Use Case | Class |
|----------|-------|
| Card container | `card` or `glass-card p-6` |
| Primary button | `btn-primary` |
| Secondary button | `btn-secondary` |
| Icon button | `btn-icon` |
| Text input | `input-field` |
| Select dropdown | `select-field` |
| Label | `input-label` |
| Section heading | `section-title` |
| Tab button | `calc-mode-tab` |
| Loading spinner | `spinner` |
| Result highlight | `result-display` |
