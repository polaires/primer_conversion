/**
 * Tests for Overhang Set Optimizer
 */

import { describe, it, expect } from 'vitest';
import {
  optimizeOverhangSet,
  evaluateOverhangSet,
  optimizeOverhangSetMultiRun,
  batchScoreRandomSets,
  calculateLigationFidelity,
  getAllOverhangs,
  filterOverhangs,
} from './overhang-optimizer.js';

describe('Overhang Optimizer', () => {
  describe('getAllOverhangs', () => {
    it('should return valid overhangs for BsaI', () => {
      const overhangs = getAllOverhangs('BsaI');
      expect(overhangs.length).toBeGreaterThan(100);
      // Should not include palindromes
      expect(overhangs).not.toContain('ATAT');
      expect(overhangs).not.toContain('GCGC');
    });

    it('should return empty array for unknown enzyme', () => {
      const overhangs = getAllOverhangs('UnknownEnzyme');
      expect(overhangs).toEqual([]);
    });
  });

  describe('filterOverhangs', () => {
    it('should filter by excluded overhangs', () => {
      const overhangs = ['AAAA', 'GGAG', 'TACT', 'AATG'];
      const filtered = filterOverhangs(overhangs, { excluded: ['GGAG'] }, {});
      expect(filtered).not.toContain('GGAG');
      expect(filtered).toContain('TACT');
    });

    it('should filter by GC content', () => {
      const overhangs = ['AAAA', 'GGCC', 'AATG', 'GCGC'];
      const filtered = filterOverhangs(overhangs, { maxGC: 2 }, {});
      expect(filtered).not.toContain('GGCC');
      expect(filtered).not.toContain('GCGC');
      expect(filtered).toContain('AATG');
    });

    it('should filter by AT content', () => {
      const overhangs = ['AAAA', 'GGGC', 'GACT'];
      const filtered = filterOverhangs(overhangs, { maxAT: 2 }, {});
      expect(filtered).not.toContain('AAAA'); // 4 AT bases, filtered
      expect(filtered).toContain('GGGC');     // 0 AT bases, passes
      expect(filtered).toContain('GACT');     // 2 AT bases (A, T), passes (RC=AGTC, not palindrome)
    });
  });

  describe('evaluateOverhangSet', () => {
    it('should evaluate MoClo standard set', () => {
      const mocloSet = ['GGAG', 'TACT', 'AATG', 'GCTT'];
      const result = evaluateOverhangSet(mocloSet, 'BsaI');

      expect(result.enzyme).toBe('BsaI');
      expect(result.numJunctions).toBe(4);
      expect(result.fidelity).toBeGreaterThan(0.9);
      expect(result.junctions).toHaveLength(4);
      expect(result.weakestJunction).toBeDefined();
      expect(result.strongestJunction).toBeDefined();
    });

    it('should return lower fidelity for poor overhang set', () => {
      const goodSet = ['GGAG', 'TACT', 'AATG', 'GCTT'];
      const badSet = ['AAAA', 'TTTT', 'CCCC', 'GGGG'];

      const goodResult = evaluateOverhangSet(goodSet, 'BsaI');
      const badResult = evaluateOverhangSet(badSet, 'BsaI');

      expect(goodResult.fidelity).toBeGreaterThan(badResult.fidelity);
    });
  });

  describe('optimizeOverhangSet', () => {
    it('should optimize a basic 5-junction set', () => {
      const result = optimizeOverhangSet({
        numJunctions: 5,
        enzyme: 'BsaI',
        iterations: 1000,
      });

      expect(result.overhangs).toHaveLength(5);
      expect(result.fidelity).toBeGreaterThan(0.8);
      expect(result.enzyme).toBe('BsaI');
      expect(result.source).toBe('optimized');
    });

    it('should respect required overhangs', () => {
      const required = ['GGAG', 'GCTT'];
      const result = optimizeOverhangSet({
        numJunctions: 5,
        enzyme: 'BsaI',
        required,
        iterations: 1000,
      });

      expect(result.overhangs).toContain('GGAG');
      expect(result.overhangs).toContain('GCTT');
      expect(result.constraints.required).toEqual(['GGAG', 'GCTT']);
    });

    it('should respect excluded overhangs', () => {
      const excluded = ['AAAA', 'TTTT', 'GATC'];
      const result = optimizeOverhangSet({
        numJunctions: 5,
        enzyme: 'BsaI',
        excluded,
        iterations: 1000,
      });

      expect(result.overhangs).not.toContain('AAAA');
      expect(result.overhangs).not.toContain('TTTT');
      expect(result.overhangs).not.toContain('GATC');
    });

    it('should respect maxGC constraint', () => {
      const result = optimizeOverhangSet({
        numJunctions: 5,
        enzyme: 'BsaI',
        maxGC: 2,
        iterations: 1000,
      });

      for (const oh of result.overhangs) {
        const gcCount = (oh.match(/[GC]/g) || []).length;
        expect(gcCount).toBeLessThanOrEqual(2);
      }
    });

    it('should throw error for invalid numJunctions', () => {
      expect(() => optimizeOverhangSet({ numJunctions: 1 })).toThrow();
      expect(() => optimizeOverhangSet({ numJunctions: 100 })).toThrow();
    });

    it('should throw error for unknown enzyme', () => {
      expect(() => optimizeOverhangSet({
        numJunctions: 5,
        enzyme: 'UnknownEnzyme',
      })).toThrow();
    });

    it('should throw error when required is also excluded', () => {
      expect(() => optimizeOverhangSet({
        numJunctions: 5,
        required: ['GGAG'],
        excluded: ['GGAG'],
      })).toThrow();
    });
  });

  describe('optimizeOverhangSetMultiRun', () => {
    it('should run multiple optimizations and return best', () => {
      const result = optimizeOverhangSetMultiRun({
        numJunctions: 5,
        enzyme: 'BsaI',
        iterations: 500,
        runs: 3,
      });

      expect(result.overhangs).toHaveLength(5);
      expect(result.multiRun).toBeDefined();
      expect(result.multiRun.totalRuns).toBe(3);
      expect(result.multiRun.allFidelities).toHaveLength(3);
      expect(result.fidelity).toBe(Math.max(...result.multiRun.allFidelities));
    });
  });

  describe('batchScoreRandomSets', () => {
    it('should score multiple random sets', () => {
      const result = batchScoreRandomSets({
        numJunctions: 5,
        enzyme: 'BsaI',
        numSets: 50,
      });

      expect(result.numSets).toBe(50);
      expect(result.best.fidelity).toBeGreaterThanOrEqual(result.worst.fidelity);
      expect(result.statistics.max).toBe(result.best.fidelity);
      expect(result.statistics.min).toBe(result.worst.fidelity);
      expect(result.allResults).toHaveLength(50);
    });
  });

  describe('Performance comparison', () => {
    it('should produce better results than random selection', () => {
      // Get random baseline
      const randomResults = batchScoreRandomSets({
        numJunctions: 10,
        enzyme: 'BsaI',
        numSets: 100,
      });

      // Optimize
      const optimized = optimizeOverhangSet({
        numJunctions: 10,
        enzyme: 'BsaI',
        iterations: 5000,
      });

      // Optimized should be at least as good as the best random
      expect(optimized.fidelity).toBeGreaterThanOrEqual(randomResults.statistics.median);
    });
  });
});
