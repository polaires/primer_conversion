/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  darkMode: 'class',
  theme: {
    extend: {
      colors: {
        // Fragment type colors
        'frag-promoter': '#ef4444',
        'frag-rbs': '#f97316',
        'frag-cds': '#22c55e',
        'frag-terminator': '#3b82f6',
        'frag-vector': '#6366f1',
        'frag-insert': '#8b5cf6',
        'frag-other': '#64748b',

        // Fidelity colors
        'fidelity-excellent': '#22c55e',
        'fidelity-good': '#84cc16',
        'fidelity-marginal': '#f59e0b',
        'fidelity-poor': '#ef4444',
      },
      animation: {
        'spin-slow': 'spin 3s linear infinite',
      },
    },
  },
  plugins: [],
}
