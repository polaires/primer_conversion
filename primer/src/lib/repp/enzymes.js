/**
 * Restriction Enzymes Module
 * Provides enzyme database and digestion logic for DNA assembly
 * Based on: https://github.com/Lattice-Automation/repp
 */

/**
 * Restriction enzyme database
 * Format: "NAME": "RECOGNITION_SEQUENCE"
 * ^ marks the cut on the top strand (5' to 3')
 * _ marks the cut on the bottom strand (3' to 5')
 */
export const ENZYMES = {
  // Common cloning enzymes
  "EcoRI": "G^AATT_C",
  "BamHI": "G^GATC_C",
  "HindIII": "A^AGCT_T",
  "XhoI": "C^TCGA_G",
  "XbaI": "T^CTAG_A",
  "SpeI": "A^CTAG_T",
  "PstI": "C_TGCA^G",
  "SalI": "G^TCGA_C",
  "NotI": "GC^GGCC_GC",
  "NcoI": "C^CATG_G",
  "NdeI": "CA^TA_TG",
  "BglII": "A^GATC_T",
  "KpnI": "G_GTAC^C",
  "SacI": "G_AGCT^C",
  "SmaI": "CCC^_GGG",

  // Golden Gate enzymes
  "BsaI": "GGTCTCN^NNNN_N",
  "BsmBI": "CGTCTCN^NNNN_N",
  "BbsI": "GAAGACNN^NNNN_N",
  "Esp3I": "CGTCTCN^NNNN_N",
  "SapI": "GCTCTTCN^NNN_N",
  "BspQI": "GCTCTTCN^NNN_N",

  // Rare cutters
  "AscI": "GG^CGCG_CC",
  "FseI": "GG_CCGG^CC",
  "PacI": "TTA_AT^TAA",
  "PmeI": "GTTT^_AAAC",
  "SfiI": "GGCCN_NNN^NGGCC",
  "SwaI": "ATTT^_AAAT",
  "SrfI": "GCCC^_GGGC",

  // Other common enzymes
  "AatII": "G_ACGT^C",
  "Acc65I": "G^GTAC_C",
  "AclI": "AA^CG_TT",
  "AfeI": "AGC^_GCT",
  "AflII": "C^TTAA_G",
  "AgeI": "A^CCGG_T",
  "ApaI": "G_GGCC^C",
  "ApaLI": "G^TGCA_C",
  "AseI": "AT^TA_AT",
  "AvaI": "C^YCGR_G",
  "AvrII": "C^CTAG_G",
  "BclI": "T^GATC_A",
  "BglI": "GCCN_NNN^NGGC",
  "BlpI": "GC^TNA_GC",
  "BsiWI": "C^GTAC_G",
  "BspEI": "T^CCGG_A",
  "BspHI": "T^CATG_A",
  "BsrGI": "T^GTAC_A",
  "BssHII": "G^CGCG_C",
  "BstBI": "TT^CG_AA",
  "BstEII": "G^GTNAC_C",
  "BstXI": "CCAN_NNNN^NTGG",
  "ClaI": "AT^CG_AT",
  "DraI": "TTT^_AAA",
  "DraIII": "CAC_NNN^GTG",
  "EagI": "C^GGCC_G",
  "EcoNI": "CCTNN^N_NNAGG",
  "EcoRV": "GAT^_ATC",
  "FspI": "TGC^_GCA",
  "HaeII": "R_GCGC^Y",
  "HincII": "GTY^_RAC",
  "HpaI": "GTT^_AAC",
  "KasI": "G^GCGC_C",
  "MfeI": "C^AATT_G",
  "MluI": "A^CGCG_T",
  "MscI": "TGG^_CCA",
  "NaeI": "GCC^_GGC",
  "NarI": "GG^CG_CC",
  "NgoMIV": "G^CCGG_C",
  "NheI": "G^CTAG_C",
  "NruI": "TCG^_CGA",
  "NsiI": "A_TGCA^T",
  "PciI": "A^CATG_T",
  "PmlI": "CAC^_GTG",
  "PsiI": "TTA^_TAA",
  "PspOMI": "G^GGCC_C",
  "PvuI": "CG_AT^CG",
  "PvuII": "CAG^_CTG",
  "SacII": "CC_GC^GG",
  "ScaI": "AGT^_ACT",
  "SfoI": "GGC^_GCC",
  "SnaBI": "TAC^_GTA",
  "SphI": "G_CATG^C",
  "SspI": "AAT^_ATT",
  "StuI": "AGG^_CCT",
  "XcmI": "CCANNNN_N^NNNNTGG",
  "XmaI": "C^CCGG_G",
  "XmnI": "GAANN^_NNTTC",
  "ZraI": "GAC^_GTC",

  // Frequent cutters
  "AluI": "AG^_CT",
  "BfaI": "C^TA_G",
  "DpnII": "^GATC_N",
  "HaeIII": "GG^_CC",
  "HhaI": "G_CG^C",
  "HinfI": "G^ANT_C",
  "HpaII": "C^CG_G",
  "MboI": "^GATC_N",
  "MseI": "T^TA_A",
  "MspI": "C^CG_G",
  "RsaI": "GT^_AC",
  "Sau3AI": "^GATC_N",
  "TaqI": "T^CG_A",
};

/**
 * IUPAC ambiguity codes for nucleotides
 */
const IUPAC_CODES = {
  'A': 'A',
  'C': 'C',
  'G': 'G',
  'T': 'T',
  'M': '[AC]',
  'R': '[AG]',
  'W': '[AT]',
  'Y': '[CT]',
  'S': '[CG]',
  'K': '[GT]',
  'H': '[ACT]',
  'D': '[AGT]',
  'V': '[ACG]',
  'B': '[CGT]',
  'N': '[ACGT]',
  'X': '[ACGT]',
};

/**
 * Parse an enzyme recognition sequence
 * @param {string} name - Enzyme name
 * @param {string} recogSeq - Recognition sequence with ^ and _ marks
 * @returns {Object} Parsed enzyme object
 */
export function parseEnzyme(name, recogSeq) {
  let cutIndex = recogSeq.indexOf('^');
  let hangIndex = recogSeq.indexOf('_');

  // Adjust indices after removing markers
  if (cutIndex < hangIndex) {
    hangIndex--;
  } else {
    cutIndex--;
  }

  const cleanSeq = recogSeq.replace(/[\^_]/g, '');

  return {
    name,
    recog: cleanSeq,
    seqCutIndex: cutIndex,
    compCutIndex: hangIndex,
    overhang: cutIndex - hangIndex,
  };
}

/**
 * Convert recognition sequence to regex pattern
 * @param {string} recog - Recognition sequence (without ^ and _)
 * @returns {string} Regex pattern
 */
export function recogToRegex(recog) {
  return recog.split('').map(c => IUPAC_CODES[c] || c).join('');
}

/**
 * Get reverse complement of a DNA sequence
 * @param {string} seq - DNA sequence
 * @returns {string} Reverse complement
 */
export function reverseComplement(seq) {
  const complement = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
    'W': 'W', 'S': 'S', 'H': 'D', 'D': 'H',
    'V': 'B', 'B': 'V', 'N': 'N', 'X': 'X',
  };

  return seq.split('').reverse().map(c => complement[c] || c).join('');
}

/**
 * Find all cut sites for enzymes in a sequence
 * @param {string} seq - DNA sequence
 * @param {Array} enzymes - Array of parsed enzyme objects
 * @returns {Object} { cuts: Array, lengths: Array }
 */
export function findCutSites(seq, enzymes) {
  const cuts = [];
  const seqUpper = seq.toUpperCase();
  const rcSeq = reverseComplement(seqUpper);

  for (const enzyme of enzymes) {
    const regex = new RegExp(recogToRegex(enzyme.recog), 'g');

    // Search forward strand
    let match;
    while ((match = regex.exec(seqUpper)) !== null) {
      cuts.push({
        index: match.index,
        enzyme,
        strand: true,  // forward
      });
    }

    // Check if palindrome - if so, don't search reverse
    const rcRecog = reverseComplement(enzyme.recog);
    if (rcRecog === enzyme.recog) {
      continue;
    }

    // Search reverse complement
    const rcRegex = new RegExp(recogToRegex(enzyme.recog), 'g');
    while ((match = rcRegex.exec(rcSeq)) !== null) {
      const fwdIndex = seq.length - match.index - enzyme.recog.length;
      cuts.push({
        index: fwdIndex,
        enzyme,
        strand: false,  // reverse
      });
    }
  }

  // Sort by index
  cuts.sort((a, b) => a.index - b.index);

  // Calculate band lengths
  const lengths = cuts.map((cut, i) => {
    const next = (i + 1) % cuts.length;
    return (cuts[next].index - cut.index + seq.length) % seq.length;
  });

  return { cuts, lengths };
}

/**
 * Digest a circular sequence with enzymes
 * @param {string} seq - DNA sequence
 * @param {Array<string>} enzymeNames - Array of enzyme names
 * @returns {Object} { linearized: string, backbone: Object } or null if no cuts
 */
export function digest(seq, enzymeNames) {
  // Parse enzymes
  const enzymes = enzymeNames.map(name => {
    const recogSeq = ENZYMES[name];
    if (!recogSeq) {
      throw new Error(`Unknown enzyme: ${name}`);
    }
    return parseEnzyme(name, recogSeq);
  });

  // Handle doubled sequences (circular fragments in databases)
  const seqUpper = seq.toUpperCase();
  const firstHalf = seqUpper.slice(0, seqUpper.length / 2);
  const secondHalf = seqUpper.slice(seqUpper.length / 2);
  const cleanSeq = firstHalf === secondHalf ? firstHalf : seqUpper;

  // Find cut sites
  const { cuts, lengths } = findCutSites(cleanSeq, enzymes);

  if (cuts.length === 0) {
    return null;  // No cut sites found
  }

  // Single cut site
  if (cuts.length === 1) {
    const cut = cuts[0];
    const overhang = cut.enzyme.seqCutIndex - cut.enzyme.compCutIndex;
    let linearized;

    if (overhang >= 0) {
      const cutIdx = (cut.index + cut.enzyme.seqCutIndex) % cleanSeq.length;
      linearized = cleanSeq.slice(cutIdx) + cleanSeq.slice(0, cutIdx);
    } else {
      const bottomIdx = (cut.index + cut.enzyme.seqCutIndex) % cleanSeq.length;
      const topIdx = (cut.index + cut.enzyme.compCutIndex) % cleanSeq.length;
      linearized = cleanSeq.slice(topIdx) + cleanSeq.slice(0, bottomIdx);
    }

    return {
      linearized,
      backbone: {
        seq: cleanSeq,
        enzymes: [cut.enzyme.name],
        cutSites: [cut.index],
        strands: [cut.strand],
      },
    };
  }

  // Multiple cut sites - find largest band
  let largestBand = 0;
  for (let i = 0; i < lengths.length; i++) {
    if (lengths[i] > lengths[largestBand]) {
      largestBand = i;
    }
  }

  const cut1 = cuts[largestBand];
  const cut2 = cuts[(largestBand + 1) % cuts.length];
  const doubled = cleanSeq + cleanSeq;

  let cut1Index = cut1.index + cut1.enzyme.compCutIndex;
  if (cut1.enzyme.seqCutIndex - cut1.enzyme.compCutIndex < 0) {
    cut1Index = cut1.index + cut1.enzyme.seqCutIndex;
  }

  let cut2Index = cut2.index + cut2.enzyme.compCutIndex;
  if (cut2.enzyme.seqCutIndex - cut2.enzyme.compCutIndex < 0) {
    cut2Index = cut2.index + cut2.enzyme.seqCutIndex;
  }

  if (cut2Index < cut1Index) {
    cut2Index += cleanSeq.length;
  }

  const linearized = doubled.slice(cut1Index, cut2Index);

  return {
    linearized,
    backbone: {
      seq: cleanSeq,
      enzymes: [cut1.enzyme.name, cut2.enzyme.name],
      cutSites: [cut1Index, cut2Index % cleanSeq.length],
      strands: [cut1.strand, cut2.strand],
    },
  };
}

/**
 * Get list of common enzyme names for UI selection
 * @returns {Array<string>} Sorted enzyme names
 */
export function getEnzymeNames() {
  return Object.keys(ENZYMES).sort();
}

/**
 * Check if an enzyme exists in the database
 * @param {string} name - Enzyme name
 * @returns {boolean}
 */
export function hasEnzyme(name) {
  return name in ENZYMES;
}

/**
 * Get enzyme recognition sequence
 * @param {string} name - Enzyme name
 * @returns {string|null} Recognition sequence or null
 */
export function getEnzymeRecognition(name) {
  return ENZYMES[name] || null;
}
