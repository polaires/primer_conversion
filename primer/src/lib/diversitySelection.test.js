/**
 * Tests for Diversity Selection Module
 */

import { describe, it, expect } from 'vitest';
import {
  calculateDistance,
  selectDiverseCandidates,
  selectParetoOptimal,
  selectHybrid,
  identifyStrengths,
  generateLabel,
  generateExplanation,
} from './diversitySelection.js';

describe('calculateDistance', () => {
  const defaultRanges = {
    maxPositionDiff: 100,
    maxLengthDiff: 10,
    maxTmDiff: 15,
  };

  it('returns 0 for identical candidates', () => {
    const a = { startPos: 10, length: 20, tm: 55 };
    const b = { startPos: 10, length: 20, tm: 55 };
    expect(calculateDistance(a, b, defaultRanges)).toBe(0);
  });

  it('returns positive distance for different candidates', () => {
    const a = { startPos: 10, length: 20, tm: 55 };
    const b = { startPos: 50, length: 25, tm: 60 };
    const distance = calculateDistance(a, b, defaultRanges);
    expect(distance).toBeGreaterThan(0);
    expect(distance).toBeLessThanOrEqual(1);
  });

  it('caps distance at 1.0', () => {
    const a = { startPos: 0, length: 15, tm: 50 };
    const b = { startPos: 200, length: 30, tm: 70 };
    const distance = calculateDistance(a, b, defaultRanges);
    expect(distance).toBeLessThanOrEqual(1);
  });

  it('handles missing startPos gracefully', () => {
    const a = { length: 20, tm: 55 };
    const b = { length: 25, tm: 60 };
    const distance = calculateDistance(a, b, defaultRanges);
    expect(distance).toBeGreaterThan(0);
  });
});

describe('selectDiverseCandidates', () => {
  it('returns empty array for empty input', () => {
    const result = selectDiverseCandidates([], 5);
    expect(result).toEqual([]);
  });

  it('returns all candidates if fewer than numToSelect', () => {
    const candidates = [
      { score: 90, startPos: 10, length: 20, tm: 55 },
      { score: 85, startPos: 20, length: 22, tm: 57 },
    ];
    const result = selectDiverseCandidates(candidates, 5);
    expect(result.length).toBe(2);
  });

  it('selects best by score first', () => {
    const candidates = [
      { score: 70, startPos: 10, length: 20, tm: 55 },
      { score: 90, startPos: 50, length: 22, tm: 57 },
      { score: 80, startPos: 30, length: 21, tm: 56 },
    ];
    const result = selectDiverseCandidates(candidates, 2);
    expect(result[0].score).toBe(90);
  });

  it('favors diverse candidates over similar high-scoring ones', () => {
    const candidates = [
      { score: 95, startPos: 10, length: 20, tm: 55, id: 'A' },
      { score: 94, startPos: 11, length: 20, tm: 55.1, id: 'B' }, // Very similar to A
      { score: 93, startPos: 12, length: 20, tm: 55.2, id: 'C' }, // Very similar to A
      { score: 85, startPos: 80, length: 25, tm: 62, id: 'D' },   // Very different
    ];

    const result = selectDiverseCandidates(candidates, 2);

    // First should be best (A)
    expect(result[0].id).toBe('A');

    // Second should favor diversity (D) over similar B or C
    // D is much more different from A despite lower score
    expect(result[1].id).toBe('D');
  });

  it('respects minScore threshold', () => {
    const candidates = [
      { score: 90, startPos: 10, length: 20, tm: 55 },
      { score: 30, startPos: 50, length: 25, tm: 60 }, // Below threshold
      { score: 85, startPos: 30, length: 22, tm: 57 },
    ];

    const result = selectDiverseCandidates(candidates, 3, { minScore: 50 });
    expect(result.length).toBe(2);
    expect(result.every(c => c.score >= 50)).toBe(true);
  });
});

describe('selectParetoOptimal', () => {
  it('returns empty array for empty input', () => {
    const result = selectParetoOptimal([], []);
    expect(result).toEqual([]);
  });

  it('identifies single Pareto-optimal candidate', () => {
    const candidates = [
      { score: 95, tmDiff: 1, dimerDG: -2 },
      { score: 90, tmDiff: 2, dimerDG: -3 },
      { score: 85, tmDiff: 3, dimerDG: -4 },
    ];

    const objectives = [
      c => c.score,      // Higher is better
      c => -c.tmDiff,    // Lower tmDiff is better (negate)
      c => c.dimerDG,    // Higher (less negative) is better
    ];

    const result = selectParetoOptimal(candidates, objectives);

    // First candidate dominates all others
    expect(result.length).toBe(1);
    expect(result[0].score).toBe(95);
  });

  it('identifies multiple Pareto-optimal candidates with trade-offs', () => {
    const candidates = [
      { id: 'A', score: 92, tmDiff: 1.5, dimerDG: -4 },  // Best score
      { id: 'B', score: 88, tmDiff: 0.5, dimerDG: -3 },  // Best tmDiff
      { id: 'C', score: 85, tmDiff: 2.0, dimerDG: -1.5 },// Best dimer
      { id: 'D', score: 86, tmDiff: 1.8, dimerDG: -4 },  // Dominated by A
    ];

    const objectives = [
      c => c.score,
      c => -c.tmDiff,
      c => c.dimerDG,
    ];

    const result = selectParetoOptimal(candidates, objectives);

    // A, B, C are Pareto-optimal; D is dominated by A
    expect(result.length).toBe(3);
    expect(result.map(c => c.id).sort()).toEqual(['A', 'B', 'C']);
  });
});

describe('selectHybrid', () => {
  it('returns all candidates with labels if fewer than numToSelect', () => {
    const candidates = [
      { score: 90, startPos: 10, length: 20, tm: 55 },
    ];

    const result = selectHybrid(candidates, 5);
    expect(result.length).toBe(1);
    expect(result[0].selectionReason).toBe('bestOverall');
  });

  it('selects top by score then adds diverse alternatives', () => {
    const candidates = [
      { score: 95, startPos: 10, length: 20, tm: 55, id: 'A' },
      { score: 94, startPos: 11, length: 20, tm: 55, id: 'B' },
      { score: 93, startPos: 12, length: 20, tm: 55, id: 'C' },
      { score: 85, startPos: 80, length: 26, tm: 62, id: 'D' },
      { score: 80, startPos: 60, length: 24, tm: 60, id: 'E' },
    ];

    const result = selectHybrid(candidates, 4, {
      numByScore: 2,
    });

    expect(result.length).toBe(4);

    // First two should be by score
    expect(result[0].selectionReason).toBe('bestOverall');
    expect(result[1].selectionReason).toBe('topByScore');

    // Rest should be diverse
    expect(['diverse', 'paretoOptimal']).toContain(result[2].selectionReason);
    expect(['diverse', 'paretoOptimal']).toContain(result[3].selectionReason);
  });
});

describe('identifyStrengths', () => {
  it('identifies bestOverall when clearly best (not tied)', () => {
    const candidate = { compositeScore: 95 };
    const allCandidates = [
      { compositeScore: 95 },
      { compositeScore: 90 }, // 5 points lower - significant gap
      { compositeScore: 85 },
    ];

    const strengths = identifyStrengths(candidate, allCandidates);
    expect(strengths).toContain('bestOverall');
  });

  it('does NOT identify bestOverall when scores are similar', () => {
    const candidate = { compositeScore: 95 };
    const allCandidates = [
      { compositeScore: 95 },
      { compositeScore: 95 }, // Tied
      { compositeScore: 94 }, // Only 1 point lower
    ];

    const strengths = identifyStrengths(candidate, allCandidates);
    expect(strengths).not.toContain('bestOverall');
  });

  it('identifies bestTmMatch when significantly better', () => {
    const candidate = {
      compositeScore: 85,
      forward: { tm: 58 },
      reverse: { tm: 58.2 }, // 0.2°C diff
    };
    const allCandidates = [
      { compositeScore: 95, forward: { tm: 60 }, reverse: { tm: 55 } }, // 5°C diff
      candidate,
    ];

    const strengths = identifyStrengths(candidate, allCandidates);
    expect(strengths).toContain('bestTmMatch');
  });

  it('identifies safestDimer when significantly safer', () => {
    const candidate = {
      compositeScore: 80,
      heterodimerDG: -2,
      forward: { tm: 58 },
      reverse: { tm: 58 }, // Same Tm diff as others - no advantage
    };
    const allCandidates = [
      { compositeScore: 95, heterodimerDG: -6, forward: { tm: 58 }, reverse: { tm: 58 } },
      candidate,
    ];

    const strengths = identifyStrengths(candidate, allCandidates);
    expect(strengths).toContain('safestDimer');
  });

  it('does NOT identify safestDimer when all similar', () => {
    const candidate = {
      compositeScore: 80,
      heterodimerDG: -2.5,
      forward: { tm: 58 },
      reverse: { tm: 58 },
    };
    const allCandidates = [
      { compositeScore: 95, heterodimerDG: -2.6, forward: { tm: 58 }, reverse: { tm: 58 } },
      candidate,
    ];

    const strengths = identifyStrengths(candidate, allCandidates);
    expect(strengths).not.toContain('safestDimer');
  });

  it('identifies shortestAmplicon when significantly shorter', () => {
    const candidate = {
      compositeScore: 80,
      ampliconLength: 150,
      heterodimerDG: -5, // Same dimer as others
      forward: { tm: 58 },
      reverse: { tm: 58 },
    };
    const allCandidates = [
      { compositeScore: 95, ampliconLength: 500, heterodimerDG: -5, forward: { tm: 58 }, reverse: { tm: 58 } },
      candidate,
    ];

    const strengths = identifyStrengths(candidate, allCandidates);
    expect(strengths).toContain('shortestAmplicon');
  });
});

describe('generateLabel', () => {
  it('generates icon label for bestOverall', () => {
    const label = generateLabel(['bestOverall'], 'excellent');
    expect(label).toBe('⭐ Best');
  });

  it('generates icon label for bestTmMatch', () => {
    const label = generateLabel(['bestTmMatch'], 'good');
    expect(label).toBe('🎯 Tm');
  });

  it('generates icon label for safestDimer', () => {
    const label = generateLabel(['safestDimer']);
    expect(label).toBe('🔗 Safe');
  });

  it('returns null with no strengths', () => {
    const label = generateLabel([], 'acceptable');
    expect(label).toBeNull();
  });

  it('prioritizes bestOverall over other strengths', () => {
    const label = generateLabel(['bestOverall', 'bestTmMatch', 'safestDimer']);
    expect(label).toBe('⭐ Best');
  });
});

describe('generateExplanation', () => {
  it('explains bestOverall with practical benefit', () => {
    const explanation = generateExplanation({}, ['bestOverall']);
    expect(explanation).toContain('composite score');
  });

  it('explains bestTmMatch with practical benefit', () => {
    const explanation = generateExplanation({}, ['bestTmMatch']);
    expect(explanation).toContain('annealing');
  });

  it('explains safestDimer with practical benefit', () => {
    const explanation = generateExplanation({}, ['safestDimer']);
    expect(explanation).toContain('primer-primer binding');
  });

  it('explains compact with practical benefit', () => {
    const explanation = generateExplanation({}, ['compact']);
    expect(explanation).toContain('synthesis cost');
  });

  it('returns null for unlabeled alternatives', () => {
    const explanation = generateExplanation({}, []);
    expect(explanation).toBeNull();
  });
});

describe('Integration: Diversity improves alternative variety', () => {
  it('selects more spread-out alternatives than pure ranking', () => {
    // Simulate clustered candidates (common in GC-rich regions)
    const candidates = [];

    // Cluster 1: positions 10-20 (best scores)
    for (let i = 0; i < 10; i++) {
      candidates.push({
        id: `cluster1-${i}`,
        score: 95 - i * 0.5,
        startPos: 10 + i,
        length: 20,
        tm: 55 + i * 0.1,
      });
    }

    // Cluster 2: positions 50-60 (medium scores)
    for (let i = 0; i < 5; i++) {
      candidates.push({
        id: `cluster2-${i}`,
        score: 85 - i * 0.5,
        startPos: 50 + i * 2,
        length: 22,
        tm: 58 + i * 0.2,
      });
    }

    // Cluster 3: positions 80-90 (lower scores)
    for (let i = 0; i < 5; i++) {
      candidates.push({
        id: `cluster3-${i}`,
        score: 75 - i * 0.5,
        startPos: 80 + i * 2,
        length: 24,
        tm: 60 + i * 0.2,
      });
    }

    // Pure ranking: would select all from cluster 1
    const byRanking = [...candidates]
      .sort((a, b) => b.score - a.score)
      .slice(0, 5);

    const rankingPositions = byRanking.map(c => c.startPos);
    const rankingRange = Math.max(...rankingPositions) - Math.min(...rankingPositions);

    // Diversity selection: should spread across clusters
    const byDiversity = selectDiverseCandidates(candidates, 5);

    const diversityPositions = byDiversity.map(c => c.startPos);
    const diversityRange = Math.max(...diversityPositions) - Math.min(...diversityPositions);

    // Diversity selection should have much larger position spread
    expect(diversityRange).toBeGreaterThan(rankingRange);

    // Should include candidates from multiple clusters
    const clusters = new Set(byDiversity.map(c => c.id.split('-')[0]));
    expect(clusters.size).toBeGreaterThan(1);
  });
});
