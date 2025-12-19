# Primer Conversion Project

## Project Goal
Convert `polaires/primer` from **JavaScript + Custom CSS** to **TypeScript + Tailwind CSS** while maintaining **100% identical UI/UX and functionality**.

## Critical Constraints

### DO NOT CHANGE
- Component layouts or visual appearance
- Calculation results (Tm, scoring, folding, etc.)
- File import/export formats (GenBank, FASTA, CSV, XLSX)
- Any user-facing behavior or workflow

### ALLOWED
- JS → TypeScript conversion
- CSS classes → Tailwind classes
- Consolidating duplicate code
- Adding proper types
- Internal refactoring (same external behavior)
- Adding dark mode support

## Key Files

| File | Purpose |
|------|---------|
| `PRIMER_CONVERSION_GUIDE.md` | Detailed conversion guide with mappings |
| `primer/` | Cloned source repository (gitignored) |

## Repository Structure

```
primer/src/
├── components/     # 25 React components (largest: GoldenGateDesigner 222KB)
├── lib/           # 60+ algorithm files (thermodynamics, scoring, assembly)
└── App.css        # 451 KB monolithic CSS (24,756 lines)
```

## Conversion Order

1. **Foundation**: main.tsx, App.tsx, types/index.ts, configs
2. **Small components** (< 20 KB): PrimerForm, SummaryStatusPanel, etc.
3. **Medium components** (20-50 KB): TmCalculator, HairpinDiagram, etc.
4. **Large components** (50-100 KB): SequencingDesigner, IsothermalAssemblyPanel
5. **Very large** (> 100 KB): UnifiedPrimerDesigner (109 KB), GoldenGateDesigner (222 KB)
6. **Library files**: Parallel with component work

## Code Patterns

### TypeScript Props
```tsx
interface ComponentProps {
  sequence: string;
  onResult: (result: PrimerResult) => void;
  options?: DesignOptions;
}

export default function Component({ sequence, onResult, options }: ComponentProps) {
```

### Event Handlers
```tsx
const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
  setValue(e.target.value);
};
```

### Tailwind Classes (from CSS)
| CSS Variable | Tailwind |
|--------------|----------|
| `--primary-color` | `primary-600` |
| `--card-bg` | `bg-white dark:bg-slate-800` |
| `--border-color` | `border-slate-200 dark:border-slate-700` |
| `--fwd-color` | `text-blue-500` (primer-fwd) |
| `--rev-color` | `text-violet-500` (primer-rev) |

## Known Issues to Fix

1. **Code duplication**: IUPAC tables in SequenceViewer AND sequenceUtils
2. **Inline styles**: 851+ occurrences to convert to Tailwind
3. **No error boundaries**: Add React error boundaries
4. **Missing accessibility**: Add ARIA labels, focus indicators

## Testing Requirements

For each converted component:
- [ ] TypeScript compiles without errors
- [ ] Visual comparison < 1% pixel diff
- [ ] All calculations produce identical results
- [ ] Existing tests still pass

## Commands

```bash
# Type check
npx tsc --noEmit

# Run dev server
npm run dev

# Run tests
npm run test

# Build
npm run build
```

## SeqViz Integration

The SeqViz library uses `.la-vz-*` CSS classes. Keep these in a separate CSS file - do NOT convert to Tailwind.

## Don't Forget

- Always read the original component before converting
- Run visual comparison after each component
- Verify calculations with test sequences
- Keep lib files as .js initially, add .d.ts type declarations
