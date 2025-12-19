/**
 * Site Creation Check for Golden Gate Assembly
 *
 * Checks if a junction position + primer design would CREATE new restriction sites.
 * This is different from findInternalSites() which finds EXISTING sites.
 *
 * Critical for preventing:
 * 1. New recognition sites at junction boundaries
 * 2. Sites created by primer annealing to template
 * 3. Sites in the final assembled product
 */

import { reverseComplement } from './enzymes.js';
import { GOLDEN_GATE_ENZYMES } from './goldengate.js';
import { OPTIMAL_FLANKING_SEQUENCES } from './goldengate-primer-optimizer.js';

/**
 * Check if a junction position + primers would CREATE new restriction sites
 *
 * @param {string} sequence - Template sequence
 * @param {number} position - Junction position (start of 4bp overhang)
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Additional options
 * @returns {Object} Site creation risk assessment
 */
export function checkSiteCreation(sequence, position, enzyme, options = {}) {
  const {
    flankingSequence = null,  // If known, use it; otherwise use default
    checkUpstream = true,     // Check upstream fragment
    checkDownstream = true,   // Check downstream fragment
  } = options;

  const enzData = GOLDEN_GATE_ENZYMES[enzyme];
  if (!enzData) {
    return {
      hasRisk: false,
      risks: [],
      error: `Unknown enzyme: ${enzyme}`,
    };
  }

  const recognition = enzData.recognition;
  const recognitionRC = reverseComplement(recognition);
  const recognitionLen = recognition.length;

  const seq = sequence.toUpperCase();
  const risks = [];

  // Get flanking sequence (what will be added to primer)
  const flanking = flankingSequence ||
    OPTIMAL_FLANKING_SEQUENCES[enzyme]?.default ||
    'AAGAAG';

  const overhang = seq.slice(position, position + 4);
  if (overhang.length !== 4) {
    return {
      hasRisk: false,
      risks: [],
      error: 'Position too close to sequence end',
    };
  }

  // ═══════════════════════════════════════════════════════════════════
  // Check 1: Junction context - sites spanning the junction
  // ═══════════════════════════════════════════════════════════════════
  // Look at sequence around the overhang to see if the junction itself
  // creates or sits within a recognition site

  const contextStart = Math.max(0, position - recognitionLen + 1);
  const contextEnd = Math.min(seq.length, position + 4 + recognitionLen - 1);
  const junctionContext = seq.slice(contextStart, contextEnd);

  // Check for recognition sites in this context
  for (let i = 0; i <= junctionContext.length - recognitionLen; i++) {
    const window = junctionContext.slice(i, i + recognitionLen);
    const absolutePos = contextStart + i;

    if (window === recognition || window === recognitionRC) {
      // Calculate where this site is relative to the junction
      const siteStart = absolutePos;
      const siteEnd = absolutePos + recognitionLen;

      // Is this site at the junction boundary? (would be part of the design)
      const overlapsOverhang = siteStart < position + 4 && siteEnd > position;

      // If the site spans the junction, it's expected behavior for GG
      // But if it's near but not at the junction, it's a problem
      if (!overlapsOverhang) {
        risks.push({
          type: 'junction_context',
          position: siteStart,
          site: window,
          orientation: window === recognition ? 'forward' : 'reverse',
          severity: 'high',
          message: `Recognition site ${window} found near junction at position ${siteStart}`,
          recommendation: 'Choose different junction position',
        });
      }
    }
  }

  // ═══════════════════════════════════════════════════════════════════
  // Check 2: Primer creates site when annealed
  // ═══════════════════════════════════════════════════════════════════
  // The primer adds flanking + recognition + spacer + overhang
  // Check if this combined sequence creates internal sites

  // Get spacer from enzyme data (varies by enzyme)
  // BsaI/BsmBI: 1 base spacer, BbsI: no spacer, SapI: 3 base spacer
  const spacerLength = enzData.spacer || 1;
  const spacer = 'N'.repeat(spacerLength); // Placeholder for spacer bases
  const fwdPrimerTail = flanking + recognition + spacer + overhang;

  // Check for internal sites in primer tail (excluding the designed site)
  const designedSiteStart = flanking.length; // Where our recognition site starts
  for (let i = 0; i < fwdPrimerTail.length - recognitionLen; i++) {
    // Skip the designed recognition site
    if (i >= designedSiteStart && i < designedSiteStart + recognitionLen) continue;

    const window = fwdPrimerTail.slice(i, i + recognitionLen);
    if (window === recognition || window === recognitionRC) {
      risks.push({
        type: 'primer_internal_site',
        position: i,
        site: window,
        primerRegion: 'forward_tail',
        severity: 'medium',
        message: `Flanking+overhang combination creates internal ${enzyme} site`,
        recommendation: 'Use different flanking sequence',
      });
    }
  }

  // ═══════════════════════════════════════════════════════════════════
  // Check 3: Assembled product creates site
  // ═══════════════════════════════════════════════════════════════════
  // After assembly, the scar (overhang) region should not contain sites
  // Check if overhang + surrounding sequence in product could form a site

  // The assembled product at this junction will have:
  // [upstream sequence][overhang][downstream sequence]
  // Check the region spanning the scar

  if (checkUpstream && position >= recognitionLen) {
    const upstreamContext = seq.slice(position - recognitionLen + 4, position + 4);
    for (let i = 0; i <= upstreamContext.length - recognitionLen; i++) {
      const window = upstreamContext.slice(i, i + recognitionLen);
      if (window === recognition || window === recognitionRC) {
        risks.push({
          type: 'upstream_product_site',
          position: position - recognitionLen + 4 + i,
          site: window,
          severity: 'high',
          message: `Assembly product will contain ${enzyme} site upstream of junction`,
          recommendation: 'Choose different junction position',
        });
      }
    }
  }

  if (checkDownstream && position + 4 + recognitionLen <= seq.length) {
    const downstreamContext = seq.slice(position, position + recognitionLen);
    for (let i = 0; i <= downstreamContext.length - recognitionLen; i++) {
      const window = downstreamContext.slice(i, i + recognitionLen);
      if (window === recognition || window === recognitionRC) {
        // Skip if this is exactly at the overhang start (expected for some designs)
        if (i !== 0) {
          risks.push({
            type: 'downstream_product_site',
            position: position + i,
            site: window,
            severity: 'high',
            message: `Assembly product will contain ${enzyme} site at junction`,
            recommendation: 'Choose different junction position',
          });
        }
      }
    }
  }

  // Determine overall severity
  const highRisks = risks.filter(r => r.severity === 'high');
  const mediumRisks = risks.filter(r => r.severity === 'medium');

  let overallSeverity = 'safe';
  if (highRisks.length > 0) {
    overallSeverity = 'critical';
  } else if (mediumRisks.length > 0) {
    overallSeverity = 'warning';
  }

  return {
    enzyme,
    position,
    overhang,
    hasRisk: risks.length > 0,
    risks,
    riskCount: risks.length,
    highRiskCount: highRisks.length,
    mediumRiskCount: mediumRisks.length,
    severity: overallSeverity,
    recommendation: risks.length > 0
      ? 'Consider alternative junction position to avoid unintended restriction sites'
      : null,
    isSafe: risks.length === 0,
  };
}

/**
 * Check multiple positions for site creation risks
 *
 * @param {string} sequence - Template sequence
 * @param {number[]} positions - Array of junction positions
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Additional options
 * @returns {Object} Batch analysis results
 */
export function checkMultipleSiteCreation(sequence, positions, enzyme, options = {}) {
  const results = positions.map(pos => ({
    position: pos,
    ...checkSiteCreation(sequence, pos, enzyme, options),
  }));

  const safePositions = results.filter(r => r.isSafe);
  const riskyPositions = results.filter(r => !r.isSafe);

  return {
    enzyme,
    totalPositions: positions.length,
    results,
    safePositions: safePositions.map(r => r.position),
    riskyPositions: riskyPositions.map(r => ({
      position: r.position,
      riskCount: r.riskCount,
      severity: r.severity,
    })),
    summary: {
      safe: safePositions.length,
      risky: riskyPositions.length,
      critical: results.filter(r => r.severity === 'critical').length,
      warning: results.filter(r => r.severity === 'warning').length,
    },
  };
}

/**
 * Find safe junction positions near a target position
 *
 * @param {string} sequence - Template sequence
 * @param {number} targetPosition - Desired junction position
 * @param {string} enzyme - Enzyme name
 * @param {Object} options - Search options
 * @returns {Object} Safe alternatives analysis
 */
export function findSafeJunctionNear(sequence, targetPosition, enzyme, options = {}) {
  const {
    searchRadius = 20,  // How far to search from target
    maxAlternatives = 5, // Maximum alternatives to return
  } = options;

  const alternatives = [];

  // Check positions from closest to furthest
  for (let offset = 0; offset <= searchRadius; offset++) {
    for (const dir of [0, 1, -1]) {
      if (offset === 0 && dir !== 0) continue; // Only check 0 once

      const pos = targetPosition + (offset * (dir || 1));

      // Bounds check
      if (pos < 0 || pos + 4 > sequence.length) continue;

      const check = checkSiteCreation(sequence, pos, enzyme, options);

      if (check.isSafe) {
        alternatives.push({
          position: pos,
          offset: pos - targetPosition,
          overhang: check.overhang,
          check,
        });

        if (alternatives.length >= maxAlternatives) {
          break;
        }
      }
    }

    if (alternatives.length >= maxAlternatives) break;
  }

  // Sort by absolute offset (closest first)
  alternatives.sort((a, b) => Math.abs(a.offset) - Math.abs(b.offset));

  const targetCheck = checkSiteCreation(sequence, targetPosition, enzyme, options);

  return {
    targetPosition,
    targetIsSafe: targetCheck.isSafe,
    targetCheck,
    alternatives: alternatives.slice(0, maxAlternatives),
    foundSafeAlternative: alternatives.length > 0,
    closestSafe: alternatives[0] || null,
    recommendation: targetCheck.isSafe
      ? 'Target position is safe'
      : alternatives.length > 0
        ? `Consider position ${alternatives[0].position} (offset ${alternatives[0].offset > 0 ? '+' : ''}${alternatives[0].offset})`
        : `No safe positions found within ${searchRadius}bp of target`,
  };
}

/**
 * Validate an entire junction set for site creation issues
 *
 * @param {string} sequence - Template sequence
 * @param {Array} junctions - Array of {position, overhang} objects
 * @param {string} enzyme - Enzyme name
 * @returns {Object} Full validation results
 */
export function validateJunctionSet(sequence, junctions, enzyme) {
  const positions = junctions.map(j => j.position);
  const batchResults = checkMultipleSiteCreation(sequence, positions, enzyme);

  // Additional check: verify no sites created between junctions in final product
  const interJunctionRisks = [];

  for (let i = 0; i < junctions.length - 1; i++) {
    const start = junctions[i].position + 4; // After first overhang
    const end = junctions[i + 1].position;   // Before second overhang

    if (end > start) {
      const segment = sequence.slice(start, end).toUpperCase();
      const enzData = GOLDEN_GATE_ENZYMES[enzyme];
      const recognition = enzData.recognition;
      const recognitionRC = reverseComplement(recognition);

      let idx = segment.indexOf(recognition);
      while (idx !== -1) {
        interJunctionRisks.push({
          type: 'inter_junction_site',
          fragmentIndex: i,
          position: start + idx,
          site: recognition,
          severity: 'high',
        });
        idx = segment.indexOf(recognition, idx + 1);
      }

      idx = segment.indexOf(recognitionRC);
      while (idx !== -1) {
        interJunctionRisks.push({
          type: 'inter_junction_site',
          fragmentIndex: i,
          position: start + idx,
          site: recognitionRC,
          severity: 'high',
        });
        idx = segment.indexOf(recognitionRC, idx + 1);
      }
    }
  }

  const allRisks = [
    ...batchResults.results.flatMap(r => r.risks),
    ...interJunctionRisks,
  ];

  return {
    enzyme,
    junctions,
    junctionResults: batchResults.results,
    interJunctionRisks,
    allRisks,
    isValid: allRisks.length === 0,
    summary: {
      totalJunctions: junctions.length,
      junctionsWithRisks: batchResults.riskyPositions.length,
      interJunctionSites: interJunctionRisks.length,
      totalRisks: allRisks.length,
    },
    recommendation: allRisks.length === 0
      ? 'All junctions validated - no site creation risks'
      : `Found ${allRisks.length} potential issue(s) - review before proceeding`,
  };
}
