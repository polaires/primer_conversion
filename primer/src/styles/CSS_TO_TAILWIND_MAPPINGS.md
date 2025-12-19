# CSS to Tailwind Class Mappings

This document maps the legacy CSS classes in `App.css` to their Tailwind equivalents in `index.css`.

## Quick Reference

### Colors (CSS Variables → Tailwind)

| CSS Variable | Tailwind Class | Usage |
|--------------|----------------|-------|
| `--primary-color` (#2563eb) | `primary-600` | Primary buttons, links |
| `--primary-dark` (#1d4ed8) | `primary-700` | Hover states |
| `--secondary-color` (#10b981) | `emerald-500` | Success states |
| `--error-color` (#ef4444) | `red-500` | Error states |
| `--text-color` (#1f2937) | `slate-800` | Body text |
| `--text-light` (#6b7280) | `slate-500` | Secondary text |
| `--bg-color` (#f9fafb) | `slate-50` | Page background |
| `--card-bg` (#ffffff) | `white` / `bg-white` | Card backgrounds |
| `--border-color` (#e5e7eb) | `slate-200` | Borders |
| `--fwd-color` (#3b82f6) | `blue-500` | Forward primer |
| `--rev-color` (#8b5cf6) | `violet-500` | Reverse primer |

### Common Layout Patterns

| Legacy CSS | Tailwind Classes |
|------------|------------------|
| `display: flex` | `flex` |
| `flex-direction: column` | `flex-col` |
| `justify-content: center` | `justify-center` |
| `justify-content: space-between` | `justify-between` |
| `align-items: center` | `items-center` |
| `gap: 8px` | `gap-2` |
| `gap: 16px` | `gap-4` |
| `gap: 24px` | `gap-6` |

### Spacing (margin/padding)

| Legacy CSS | Tailwind |
|------------|----------|
| `4px` | `1` |
| `8px` | `2` |
| `12px` | `3` |
| `16px` | `4` |
| `20px` | `5` |
| `24px` | `6` |
| `32px` | `8` |

### Component Class Mappings

#### Cards
```
.card, .primer-form, .results-container → card
.glass-card → glass-card
```

#### Buttons
```
.btn-primary, .submit-btn → btn-primary
.btn-secondary, .export-btn → btn-secondary
.btn-icon → btn-icon
```

#### Form Inputs
```
.form-input, input[type="text"] → input-field
.form-label → input-label
.form-select → select-field
```

#### Mode/Tab Navigation
```
.mode-tabs button → calc-mode-tab
.mode-tabs button.active → calc-mode-tab active
```

#### Status Badges
```
.status-excellent → status-badge excellent
.status-good → status-badge good
.status-warning → status-badge warning
.status-poor → status-badge poor
```

#### Primer-Specific
```
.primer-fwd → primer-fwd-bg (background) or primer-label fwd (text)
.primer-rev → primer-rev-bg (background) or primer-label rev (text)
.sequence-display → sequence-display
.tm-value → tm-value
```

#### Score Indicators
```
.score-excellent, .text-green → score-excellent
.score-good, .text-blue → score-good
.score-acceptable, .text-amber → score-acceptable
.score-poor, .text-red → score-poor
```

### Inline Style Conversions

| Inline Style | Tailwind Class |
|--------------|----------------|
| `style={{ marginBottom: '16px' }}` | `className="mb-4"` |
| `style={{ padding: '12px' }}` | `className="p-3"` |
| `style={{ display: 'flex', gap: '8px' }}` | `className="flex gap-2"` |
| `style={{ fontSize: '14px' }}` | `className="text-sm"` |
| `style={{ fontWeight: 'bold' }}` | `className="font-bold"` |
| `style={{ color: '#ef4444' }}` | `className="text-red-500"` |
| `style={{ backgroundColor: '#f9fafb' }}` | `className="bg-slate-50"` |
| `style={{ borderRadius: '8px' }}` | `className="rounded-lg"` |
| `style={{ border: '1px solid #e5e7eb' }}` | `className="border border-slate-200"` |

### Font Sizes

| Size | Tailwind |
|------|----------|
| 12px | `text-xs` |
| 14px | `text-sm` |
| 16px | `text-base` |
| 18px | `text-lg` |
| 20px | `text-xl` |
| 24px | `text-2xl` |

### Border Radius

| Size | Tailwind |
|------|----------|
| 4px | `rounded` |
| 6px | `rounded-md` |
| 8px | `rounded-lg` |
| 12px | `rounded-xl` |
| 9999px (full) | `rounded-full` |

## Migration Strategy

1. **Keep SeqViz styles separate** - All `.la-vz-*` classes remain in `seqviz.css`
2. **Use Tailwind component classes** - Prefer classes defined in `index.css` `@layer components`
3. **Convert inline styles** - Replace `style={{}}` with Tailwind utility classes
4. **Dark mode support** - All component classes include `dark:` variants

## Files Structure

```
src/styles/
├── index.css       # Tailwind config + component classes
├── seqviz.css      # SeqViz library styles (DO NOT CONVERT)
└── CSS_TO_TAILWIND_MAPPINGS.md  # This file
```

## Notes

- App.css is gradually being deprecated as styles move to Tailwind
- Components should import from `index.css` only
- Test visual appearance after converting each component
