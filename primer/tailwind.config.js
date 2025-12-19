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
        // Primary colors (matches original --primary-color)
        primary: {
          50: '#eff6ff',
          100: '#dbeafe',
          200: '#bfdbfe',
          300: '#93c5fd',
          400: '#60a5fa',
          500: '#3b82f6',
          600: '#2563eb',  // --primary-color
          700: '#1d4ed8',  // --primary-dark
          800: '#1e40af',
          900: '#1e3a8a',
        },
        // Primer-specific colors (for forward/reverse primers)
        primer: {
          fwd: '#3b82f6',      // Blue for forward --fwd-color
          rev: '#8b5cf6',      // Purple for reverse --rev-color
          template: '#10b981', // Green for template
          mutation: '#f59e0b', // Amber for mutations
        },
      },
      fontFamily: {
        sans: ['Inter', 'system-ui', '-apple-system', 'BlinkMacSystemFont', 'Segoe UI', 'Roboto', 'sans-serif'],
        mono: ['JetBrains Mono', 'Monaco', 'Consolas', 'Liberation Mono', 'Courier New', 'monospace'],
      },
      boxShadow: {
        'glass': '0 8px 32px 0 rgba(31, 38, 135, 0.37)',
      },
    },
  },
  plugins: [],
}
