/**
 * Find off-target binding sites.
 * Ported from primers Python library
 */

/**
 * Return a list of off-target counts for primers that end at that index.
 *
 * For example, offTargets[20] -> returns the number of offtarget binding
 * sites whose last bp ends in the 20th index of `seq`
 *
 * @param seq - The template sequence being checked
 * @param checkSeq - The sequence being checked for offtarget binding sites
 * @returns A list of binding site counts for primers whose last bp is within that index
 */
export function offTargets(seq: string, checkSeq: string): number[] {
  const checkMap: Record<string, number> = {};
  const mutateMap: Record<string, string> = { A: "TGC", T: "AGC", G: "ATC", C: "ATG" };

  /**
   * Generate all single-base mutations of a sequence
   * @param s - Input sequence
   * @returns Array of mutated sequences including original
   */
  function mutate(s: string): string[] {
    const mutatedSites = [s];
    for (let i = 0; i < s.length; i++) {
      const c = s[i];
      if (mutateMap[c]) {
        for (const m of mutateMap[c]) {
          mutatedSites.push(s.slice(0, i) + m + s.slice(i + 1));
        }
      }
    }
    return mutatedSites;
  }

  // Build the check map with all possible binding sites (with 0 or 1 mismatch)
  for (let s = 0; s <= checkSeq.length - 10; s++) {
    const subseq = checkSeq.slice(s, s + 10);
    for (const m of mutate(subseq)) {
      checkMap[m] = (checkMap[m] || 0) + 1;
    }
  }

  // Create the cache
  const cache: number[] = new Array(seq.length).fill(0);

  for (let e = 10; e <= seq.length; e++) {
    const ss = seq.slice(e - 10, e);

    // assumed to bind at least once
    const count = (checkMap[ss] || 0) + (checkMap[reverseComplement(ss)] || 0) - 1;
    cache[e - 1] = Math.max(0, count);
  }

  return cache;
}

/**
 * Return the reverse complement of a DNA sequence.
 *
 * @param seq - The template sequence
 * @returns The reverse complement of the template sequence
 */
function reverseComplement(seq: string): string {
  const rc: Record<string, string> = { A: "T", T: "A", G: "C", C: "G" };
  return seq
    .split("")
    .reverse()
    .map((c) => rc[c])
    .join("");
}

export { reverseComplement };
