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
    outDir: 'dist',
    sourcemap: true,
    chunkSizeWarningLimit: 5000,
  },

  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
});
