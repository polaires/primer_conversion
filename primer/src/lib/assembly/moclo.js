/**
 * MoClo (Modular Cloning) / Hierarchical Assembly Module
 *
 * Supports multi-level Golden Gate assembly workflows:
 * - Level 0: Basic parts (promoters, RBS, CDS, terminators)
 * - Level 1: Transcription units (assembled from Level 0)
 * - Level 2: Multigene constructs (assembled from Level 1)
 *
 * References:
 * - Weber et al. 2011 PLOS ONE (original MoClo)
 * - Werner et al. 2012 PLOS ONE (plant MoClo)
 * - Lee et al. 2015 ACS Synth Biol (MoClo toolkit)
 */

import { reverseComplement } from '../repp/enzymes.js';

/**
 * Standard MoClo part types and their fusion sites
 * Based on Weber et al. 2011 / Phytobricks standard
 */
export const MOCLO_PART_TYPES = {
  // Level 0 parts (basic genetic elements)
  promoter: {
    code: 'A',
    description: 'Promoter with 5\' UTR',
    level: 0,
    upstream: 'GGAG',   // Entry fusion site
    downstream: 'TACT', // Exit fusion site
  },
  rbs: {
    code: 'B',
    description: 'Ribosome binding site',
    level: 0,
    upstream: 'TACT',
    downstream: 'AATG',
  },
  cds: {
    code: 'C',
    description: 'Coding sequence (start to stop-1)',
    level: 0,
    upstream: 'AATG',   // ATG start codon embedded
    downstream: 'GCTT', // Before stop codon
  },
  cds_ns: {
    code: 'C*',
    description: 'CDS without stop (for C-terminal fusions)',
    level: 0,
    upstream: 'AATG',
    downstream: 'CGCT', // Alternative for fusions
  },
  terminator: {
    code: 'D',
    description: 'Terminator with 3\' UTR',
    level: 0,
    upstream: 'GCTT',
    downstream: 'GGTA',
  },
  tag_n: {
    code: 'N',
    description: 'N-terminal tag',
    level: 0,
    upstream: 'AATG',
    downstream: 'AGGT',
  },
  tag_c: {
    code: 'CT',
    description: 'C-terminal tag',
    level: 0,
    upstream: 'CGCT',
    downstream: 'GCTT',
  },
  linker: {
    code: 'L',
    description: 'Linker sequence',
    level: 0,
    upstream: 'AGGT',
    downstream: 'AATG',
  },
};

/**
 * Standard Level 1 destination vectors
 */
export const LEVEL1_POSITIONS = {
  position1: {
    upstream: 'GGAG',
    downstream: 'TGCC',
    enzyme: 'BsaI',
  },
  position2: {
    upstream: 'TGCC',
    downstream: 'CAGA',
    enzyme: 'BsaI',
  },
  position3: {
    upstream: 'CAGA',
    downstream: 'GCTT',
    enzyme: 'BsaI',
  },
  position4: {
    upstream: 'GCTT',
    downstream: 'GGTA',
    enzyme: 'BsaI',
  },
  position5: {
    upstream: 'GGTA',
    downstream: 'CGCT',
    enzyme: 'BsaI',
  },
  position6: {
    upstream: 'CGCT',
    downstream: 'TACT',
    enzyme: 'BsaI',
  },
};

/**
 * Standard Level 2 positions for multi-TU assemblies
 */
export const LEVEL2_POSITIONS = {
  tu1: {
    upstream: 'GGAG',
    downstream: 'AATG',
    enzyme: 'BsmBI',
  },
  tu2: {
    upstream: 'AATG',
    downstream: 'TAGG',
    enzyme: 'BsmBI',
  },
  tu3: {
    upstream: 'TAGG',
    downstream: 'CGAA',
    enzyme: 'BsmBI',
  },
  tu4: {
    upstream: 'CGAA',
    downstream: 'ATCC',
    enzyme: 'BsmBI',
  },
  tu5: {
    upstream: 'ATCC',
    downstream: 'TTAC',
    enzyme: 'BsmBI',
  },
  tu6: {
    upstream: 'TTAC',
    downstream: 'GCTT',
    enzyme: 'BsmBI',
  },
};

/**
 * Validate that parts are compatible for assembly
 */
function validatePartCompatibility(parts) {
  const issues = [];

  for (let i = 0; i < parts.length - 1; i++) {
    const current = parts[i];
    const next = parts[i + 1];

    // Check fusion site compatibility
    if (current.downstream !== next.upstream) {
      issues.push({
        type: 'incompatible_fusion',
        position: i,
        parts: [current.name || current.type, next.name || next.type],
        expected: current.downstream,
        found: next.upstream,
        message: `Part ${i + 1} downstream (${current.downstream}) doesn't match part ${i + 2} upstream (${next.upstream})`,
      });
    }
  }

  return {
    compatible: issues.length === 0,
    issues,
  };
}

/**
 * Plan a Level 0 to Level 1 assembly
 * Assembles basic parts into a transcription unit
 *
 * @param {Object[]} parts - Array of Level 0 parts
 * @param {Object} options - Assembly options
 * @returns {Object} Assembly plan
 */
export function planLevel1Assembly(parts, options = {}) {
  const {
    vectorUpstream = 'GGAG',
    vectorDownstream = 'GGTA',
    enzyme = 'BsaI',
    vector = 'pICH47742', // Standard MoClo Level 1 acceptor
  } = options;

  // Add vector fusion sites
  const assemblyParts = [
    { name: 'Vector (5\')', upstream: null, downstream: vectorUpstream, isVector: true },
    ...parts,
    { name: 'Vector (3\')', upstream: vectorDownstream, downstream: null, isVector: true },
  ];

  // Validate compatibility
  const validation = validatePartCompatibility(assemblyParts.slice(1, -1));

  // Check first and last parts match vector
  if (parts.length > 0) {
    if (parts[0].upstream !== vectorUpstream) {
      validation.issues.push({
        type: 'vector_mismatch_5',
        expected: vectorUpstream,
        found: parts[0].upstream,
        message: `First part upstream (${parts[0].upstream}) doesn't match vector (${vectorUpstream})`,
      });
      validation.compatible = false;
    }

    if (parts[parts.length - 1].downstream !== vectorDownstream) {
      validation.issues.push({
        type: 'vector_mismatch_3',
        expected: vectorDownstream,
        found: parts[parts.length - 1].downstream,
        message: `Last part downstream (${parts[parts.length - 1].downstream}) doesn't match vector (${vectorDownstream})`,
      });
      validation.compatible = false;
    }
  }

  // Collect all overhangs for fidelity calculation
  const overhangs = [vectorUpstream];
  for (const part of parts) {
    if (part.downstream && part.downstream !== vectorDownstream) {
      overhangs.push(part.downstream);
    }
  }
  overhangs.push(vectorDownstream);

  return {
    level: 1,
    enzyme,
    vector,
    parts: parts.map((p, i) => ({
      position: i + 1,
      name: p.name || `Part ${i + 1}`,
      type: p.type,
      upstream: p.upstream,
      downstream: p.downstream,
      sequence: p.sequence,
    })),
    overhangs,
    numJunctions: overhangs.length,
    validation,
    protocol: generateMoCloProtocol(parts.length, enzyme, 1),
  };
}

/**
 * Plan a Level 1 to Level 2 assembly
 * Assembles multiple transcription units into a multigene construct
 *
 * @param {Object[]} transcriptionUnits - Array of Level 1 parts
 * @param {Object} options - Assembly options
 * @returns {Object} Assembly plan
 */
export function planLevel2Assembly(transcriptionUnits, options = {}) {
  const {
    enzyme = 'BsmBI',
    vector = 'pAGM4673', // Standard MoClo Level 2 acceptor
  } = options;

  // Assign positions to TUs
  const positions = Object.entries(LEVEL2_POSITIONS);
  if (transcriptionUnits.length > positions.length) {
    throw new Error(`Maximum ${positions.length} transcription units supported in Level 2`);
  }

  const assemblyParts = transcriptionUnits.map((tu, i) => {
    const [posName, posData] = positions[i];
    return {
      position: i + 1,
      positionName: posName,
      name: tu.name || `TU${i + 1}`,
      upstream: posData.upstream,
      downstream: posData.downstream,
      level1Parts: tu.parts || [],
      sequence: tu.sequence,
    };
  });

  // Collect overhangs
  const overhangs = assemblyParts.map(p => p.upstream);
  overhangs.push(assemblyParts[assemblyParts.length - 1].downstream);

  return {
    level: 2,
    enzyme,
    vector,
    transcriptionUnits: assemblyParts,
    overhangs,
    numJunctions: overhangs.length,
    totalParts: transcriptionUnits.reduce((sum, tu) => sum + (tu.parts?.length || 1), 0),
    protocol: generateMoCloProtocol(transcriptionUnits.length, enzyme, 2),
  };
}

/**
 * Plan a complete hierarchical assembly from basic parts
 *
 * @param {Object} design - Design specification
 * @param {Object[]} design.transcriptionUnits - Array of TU specifications
 * @returns {Object} Complete assembly plan
 */
export function planHierarchicalAssembly(design) {
  const { transcriptionUnits } = design;

  const level1Assemblies = [];

  // Plan each Level 1 assembly
  for (let i = 0; i < transcriptionUnits.length; i++) {
    const tu = transcriptionUnits[i];
    const l1 = planLevel1Assembly(tu.parts, {
      vectorUpstream: LEVEL2_POSITIONS[`tu${i + 1}`]?.upstream || 'GGAG',
      vectorDownstream: LEVEL2_POSITIONS[`tu${i + 1}`]?.downstream || 'GGTA',
    });

    level1Assemblies.push({
      name: tu.name || `TU${i + 1}`,
      ...l1,
    });
  }

  // Plan Level 2 assembly
  const level2 = planLevel2Assembly(transcriptionUnits.map((tu, i) => ({
    name: tu.name || `TU${i + 1}`,
    parts: tu.parts,
  })));

  // Calculate total parts and reactions
  const totalLevel0Parts = level1Assemblies.reduce(
    (sum, l1) => sum + l1.parts.length, 0
  );

  return {
    summary: {
      totalTranscriptionUnits: transcriptionUnits.length,
      totalLevel0Parts,
      totalReactions: level1Assemblies.length + 1, // L1s + L2
    },

    level1: {
      enzyme: 'BsaI',
      assemblies: level1Assemblies,
      estimatedTime: '2-3 hours per assembly',
    },

    level2: {
      enzyme: 'BsmBI',
      assembly: level2,
      estimatedTime: '2-3 hours',
    },

    timeline: {
      day1: 'Level 1 assemblies (can run in parallel)',
      day2: 'Transform Level 1, miniprep',
      day3: 'Level 2 assembly',
      day4: 'Transform Level 2, verify',
      totalDays: '4-5 days',
    },

    costEstimate: {
      level0Parts: `${totalLevel0Parts} parts × $0.10/bp average`,
      enzymes: 'BsaI-HFv2 + BsmBI-v2',
      vectors: `${level1Assemblies.length} Level 1 + 1 Level 2`,
    },
  };
}

/**
 * Generate MoClo protocol for assembly
 */
function generateMoCloProtocol(numParts, enzyme, level) {
  const enzymeVolumes = {
    BsaI: { enzyme: 0.5, ligase: 0.5 },
    BsmBI: { enzyme: 0.75, ligase: 0.5 },
  };

  const volumes = enzymeVolumes[enzyme] || enzymeVolumes.BsaI;

  return {
    reaction: {
      partVolume: `${(40 / numParts).toFixed(1)} fmol each`,
      vectorVolume: '40 fmol',
      buffer: '2 µL 10× T4 ligase buffer',
      enzyme: `${volumes.enzyme} µL ${enzyme}-HFv2`,
      ligase: `${volumes.ligase} µL T4 DNA ligase`,
      water: 'to 20 µL',
    },
    cycling: {
      steps: [
        { temp: 37, time: '5 min', description: 'Digestion' },
        { temp: 16, time: '5 min', description: 'Ligation' },
      ],
      cycles: 30,
      finalDigest: { temp: 55, time: '10 min' },
      heatInactivation: { temp: 80, time: '10 min' },
    },
    transformation: {
      cells: 'NEB 10-beta or similar',
      volume: '2-5 µL reaction',
      recovery: '1 hour at 37°C',
      plating: 'LB + appropriate antibiotic',
    },
    screening: {
      method: 'Colony PCR or restriction digest',
      coloniesPerAssembly: level === 1 ? 3 : 6,
      expectedEfficiency: level === 1 ? '>90%' : '>80%',
    },
  };
}

/**
 * Suggest part arrangements for common expression constructs
 */
export function suggestPartArrangement(constructType, options = {}) {
  const arrangements = {
    // Simple expression
    simple_expression: {
      description: 'Basic protein expression cassette',
      parts: [
        { type: 'promoter', suggestion: 'T7, lac, or constitutive' },
        { type: 'rbs', suggestion: 'Standard or optimized RBS' },
        { type: 'cds', suggestion: 'Your gene of interest' },
        { type: 'terminator', suggestion: 'T7 or rrnB terminator' },
      ],
    },

    // N-terminal tagged
    nterm_tagged: {
      description: 'N-terminal fusion protein',
      parts: [
        { type: 'promoter', suggestion: 'T7 or lac' },
        { type: 'rbs', suggestion: 'Standard RBS' },
        { type: 'tag_n', suggestion: 'His6, GST, MBP, etc.' },
        { type: 'linker', suggestion: 'GS linker optional' },
        { type: 'cds', suggestion: 'Your gene (no start codon)' },
        { type: 'terminator', suggestion: 'Standard terminator' },
      ],
    },

    // C-terminal tagged
    cterm_tagged: {
      description: 'C-terminal fusion protein',
      parts: [
        { type: 'promoter', suggestion: 'T7 or lac' },
        { type: 'rbs', suggestion: 'Standard RBS' },
        { type: 'cds_ns', suggestion: 'Your gene (no stop codon)' },
        { type: 'linker', suggestion: 'GS linker optional' },
        { type: 'tag_c', suggestion: 'His6, FLAG, etc.' },
        { type: 'terminator', suggestion: 'Standard terminator' },
      ],
    },

    // Operon
    operon: {
      description: 'Multi-gene operon',
      parts: [
        { type: 'promoter', suggestion: 'Strong promoter' },
        { type: 'rbs', suggestion: 'RBS 1' },
        { type: 'cds', suggestion: 'Gene 1' },
        { type: 'rbs', suggestion: 'RBS 2 (intercistronic)' },
        { type: 'cds', suggestion: 'Gene 2' },
        // ... can extend
        { type: 'terminator', suggestion: 'Strong terminator' },
      ],
    },
  };

  return arrangements[constructType] || arrangements.simple_expression;
}

/**
 * Check if a sequence is compatible with MoClo
 * (i.e., doesn't contain internal BsaI/BsmBI sites)
 */
export function checkMoCloCompatibility(sequence, options = {}) {
  const { level = 0 } = options;

  const sites = {
    BsaI: /GGTCTC|GAGACC/gi,
    BsmBI: /CGTCTC|GAGACG/gi,
  };

  const issues = [];

  // Level 0 parts should not have BsaI sites
  // Level 1 parts should not have BsmBI sites
  const enzymeToCheck = level === 0 ? 'BsaI' : 'BsmBI';
  const regex = sites[enzymeToCheck];

  let match;
  while ((match = regex.exec(sequence)) !== null) {
    issues.push({
      enzyme: enzymeToCheck,
      site: match[0],
      position: match.index,
      message: `Internal ${enzymeToCheck} site at position ${match.index}`,
    });
  }

  return {
    compatible: issues.length === 0,
    issues,
    suggestion: issues.length > 0
      ? `Remove internal ${enzymeToCheck} sites by silent mutation`
      : null,
  };
}

/**
 * Generate domestication primers for MoClo
 * Adds appropriate fusion sites to a part
 */
export function generateDomesticationPrimers(part, options = {}) {
  const {
    partType,
    sequence,
    name = 'Part',
  } = part;

  const {
    enzyme = 'BsaI',
    addBuffer = true,
  } = options;

  const partSpec = MOCLO_PART_TYPES[partType];
  if (!partSpec) {
    throw new Error(`Unknown part type: ${partType}`);
  }

  // BsaI recognition: GGTCTC (N)
  // Cut produces 4-base 5' overhang
  const enzymeSite = enzyme === 'BsaI' ? 'GGTCTC' : 'CGTCTC';
  const enzymeBuffer = 'AAGAAG'; // Buffer for efficient cutting

  // Forward primer: buffer + enzyme site + spacer + fusion site + annealing
  const fwdPrefix = addBuffer
    ? enzymeBuffer + enzymeSite + 'N' + partSpec.upstream
    : enzymeSite + 'N' + partSpec.upstream;

  // Reverse primer: buffer + enzyme site (RC) + spacer + fusion site (RC) + annealing
  const revEnzymeSite = reverseComplement(enzymeSite);
  const revFusionSite = reverseComplement(partSpec.downstream);
  const revPrefix = addBuffer
    ? enzymeBuffer + revEnzymeSite + 'N' + revFusionSite
    : revEnzymeSite + 'N' + revFusionSite;

  // Get annealing regions from sequence
  const fwdAnnealing = sequence.substring(0, 20);
  const revAnnealing = reverseComplement(sequence.substring(sequence.length - 20));

  return {
    partType,
    partName: name,
    fusionSites: {
      upstream: partSpec.upstream,
      downstream: partSpec.downstream,
    },

    forward: {
      name: `${name}_F`,
      sequence: fwdPrefix + fwdAnnealing,
      length: fwdPrefix.length + 20,
      annealingRegion: fwdAnnealing,
    },

    reverse: {
      name: `${name}_R`,
      sequence: revPrefix + revAnnealing,
      length: revPrefix.length + 20,
      annealingRegion: revAnnealing,
    },

    notes: [
      `Adds ${enzyme} sites for MoClo assembly`,
      `Fusion sites: ${partSpec.upstream} / ${partSpec.downstream}`,
      'Verify sequence matches your gene exactly',
    ],
  };
}
