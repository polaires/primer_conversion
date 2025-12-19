import { primers, LEN_MIN, LEN_MAX } from '../src/lib/primers.js';
import { calculateTmQ5, calculateGC } from '../src/lib/tmQ5.js';
import { reverseComplement } from '../src/lib/sequenceUtils.js';

// pUC19 sequence
const puc19 = `TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTC`;

// Extract region 1137-1639 (0-indexed: 1136-1638)
const start = 1136;
const end = 1638;
const regionSeq = puc19.slice(start, end);

console.log("=== pUC19 1137-1639 ANALYSIS ===\n");
console.log(`Region length: ${regionSeq.length}bp`);
console.log(`Region GC: ${(calculateGC(regionSeq) * 100).toFixed(1)}%\n`);

// Analyze 5' end (forward primer region)
console.log("=== FORWARD PRIMER REGION (5' end) ===");
const fwd5End = regionSeq.slice(0, LEN_MAX + 5);
console.log(`First ${fwd5End.length}bp: ${fwd5End}`);
console.log(`GC: ${(calculateGC(fwd5End) * 100).toFixed(1)}%\n`);

// Check forward primer Tm at different lengths
console.log("Forward primer Tm by length:");
for (let len = LEN_MIN; len <= LEN_MAX; len++) {
  const seq = regionSeq.slice(0, len);
  const tm = calculateTmQ5(seq);
  const gc = calculateGC(seq);
  console.log(`  ${len}bp: Tm=${tm.toFixed(1)}C, GC=${(gc*100).toFixed(0)}%`);
}

// Analyze 3' end (reverse primer region)
console.log("\n=== REVERSE PRIMER REGION (3' end) ===");
const rev3End = regionSeq.slice(-LEN_MAX - 5);
const revCompSeq = reverseComplement(rev3End);
console.log(`Last ${rev3End.length}bp: ${rev3End}`);
console.log(`Reverse complement: ${revCompSeq}`);
console.log(`GC: ${(calculateGC(rev3End) * 100).toFixed(1)}%\n`);

// Check reverse primer Tm at different lengths
console.log("Reverse primer Tm by length (from 3' end, reverse complemented):");
for (let len = LEN_MIN; len <= LEN_MAX; len++) {
  const templateRegion = regionSeq.slice(-len);
  const seq = reverseComplement(templateRegion);
  const tm = calculateTmQ5(seq);
  const gc = calculateGC(seq);
  console.log(`  ${len}bp: Tm=${tm.toFixed(1)}C, GC=${(gc*100).toFixed(0)}%`);
}

// Now run primers() to see what it actually picks
console.log("\n=== ACTUAL PRIMERS SELECTED ===");
const [fwd, rev] = primers(regionSeq, {
  optimalTm: 62,
  useCompositeScore: true,
});

console.log(`Forward: ${fwd.seq.length}bp, Tm=${fwd.tm.toFixed(1)}C, GC=${(fwd.gc*100).toFixed(0)}%`);
console.log(`Reverse: ${rev.seq.length}bp, Tm=${rev.tm.toFixed(1)}C, GC=${(rev.gc*100).toFixed(0)}%`);
console.log(`Tm diff: ${Math.abs(fwd.tm - rev.tm).toFixed(1)}C`);

// Show what would happen if we could extend the reverse primer beyond 32bp
console.log("\n=== WHAT IF WE EXTENDED REVERSE PRIMER BEYOND 32bp? ===");
for (let len = LEN_MAX; len <= 45; len++) {
  const templateRegion = regionSeq.slice(-len);
  const seq = reverseComplement(templateRegion);
  const tm = calculateTmQ5(seq);
  const gc = calculateGC(seq);
  const tmDiff = Math.abs(fwd.tm - tm);
  const marker = tmDiff <= 2 ? ' <-- GOOD Tm match!' : '';
  console.log(`  ${len}bp: Tm=${tm.toFixed(1)}C, GC=${(gc*100).toFixed(0)}%, TmDiff=${tmDiff.toFixed(1)}C${marker}`);
}

// KEY INSIGHT: What if we look at DIFFERENT positions on the 3' end?
console.log("\n=== KEY: SLIDING WINDOW ON 3' END ===");
console.log("Looking for reverse primer positions with better Tm...\n");

// Try sliding window positions at 3' end (move back from exact end)
for (let offset = 0; offset <= 30; offset += 5) {
  console.log(`--- Offset ${offset}bp from 3' end ---`);
  for (let len = 18; len <= LEN_MAX; len += 4) {
    const templateStart = regionSeq.length - len - offset;
    if (templateStart < 0) continue;
    const templateRegion = regionSeq.slice(templateStart, templateStart + len);
    const seq = reverseComplement(templateRegion);
    const tm = calculateTmQ5(seq);
    const gc = calculateGC(seq);
    const tmDiff = Math.abs(fwd.tm - tm);
    const marker = tmDiff <= 2 ? ' <-- GOOD!' : (tmDiff <= 5 ? ' ok' : '');
    console.log(`  ${len}bp @ pos ${templateStart}: Tm=${tm.toFixed(1)}C, GC=${(gc*100).toFixed(0)}%, TmDiff=${tmDiff.toFixed(1)}C${marker}`);
  }
}

// What's at the 3' end region?
console.log("\n=== DETAILED 3' END SEQUENCE ANALYSIS ===");
const last50 = regionSeq.slice(-50);
console.log(`Last 50bp: ${last50}`);
console.log(`GC: ${(calculateGC(last50) * 100).toFixed(1)}%`);

// Check each 20bp window
console.log("\nSliding 22bp window GC content:");
for (let i = last50.length - 22; i >= 0; i -= 5) {
  const window = last50.slice(i, i + 22);
  const windowGC = calculateGC(window);
  console.log(`  pos ${i}: ${window} GC=${(windowGC*100).toFixed(0)}%`);
}
