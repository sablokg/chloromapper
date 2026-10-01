//! Find the pair of inverted repeats (IRa / IRb) within a single sequence.
//!
//! Approach: a chloroplast genome `S` (length `n`) that has already been
//! correctly assembled into one circular contig contains two regions,
//! `region1 = S[a..b)` and `region2 = S[c..d)` (with `a < b <= c < d <= n`),
//! such that `region2` is approximately the reverse complement of `region1`.
//!
//! Let `T = revcomp(S)`. Then `T[n-b .. n-a)` is exactly `revcomp(region1)`,
//! so finding long near-exact matches between `T` and `S` is equivalent to
//! finding inverted repeats within `S`. We do this with classic seed
//! (k-mer index) + diagonal chaining + X-drop extension, which is linear-ish
//! in practice and needs no external alignment crate.

use crate::seq::{bases_match, norm, revcomp};
use std::collections::HashMap;

/*
Gaurav Sablok
gsablok@proton.me
*/

/// A candidate inverted-repeat pair, in original `S` coordinates (0-based,
/// half-open), with `region1` occurring before `region2`.
#[derive(Debug, Clone, Copy)]
pub struct IrPair {
    pub region1_start: usize,
    pub region1_end: usize,
    pub region2_start: usize,
    pub region2_end: usize,
    pub matches: usize,
    pub mismatches: usize,
}

impl IrPair {
    pub fn region1_len(&self) -> usize {
        self.region1_end - self.region1_start
    }
    pub fn region2_len(&self) -> usize {
        self.region2_end - self.region2_start
    }
    pub fn identity(&self) -> f64 {
        let total = self.matches + self.mismatches;
        if total == 0 {
            0.0
        } else {
            self.matches as f64 / total as f64
        }
    }
}

pub struct IrFinderConfig {
    /// Seed k-mer length. Must be <= 31 (packed into a u64, 2 bits/base).
    pub kmer: usize,
    /// Max gap (in T-coordinates) allowed between consecutive seeds on the
    /// same diagonal before starting a new chain.
    pub max_seed_gap: usize,
    /// X-drop threshold for greedy boundary extension.
    pub x_drop: i32,
    /// Minimum final (extended) match length to be reported as a candidate.
    pub min_len: usize,
}

impl Default for IrFinderConfig {
    fn default() -> Self {
        IrFinderConfig {
            kmer: 21,
            max_seed_gap: 60,
            x_drop: 40,
            min_len: 1000,
        }
    }
}

/// 2-bit-pack a k-mer starting at `pos` in `seq`. Returns `None` if the
/// window contains a non-ACGT base (case-insensitive) or runs off the end.
fn pack_kmer(seq: &[u8], pos: usize, k: usize) -> Option<u64> {
    if pos + k > seq.len() {
        return None;
    }
    let mut code: u64 = 0;
    for &b in &seq[pos..pos + k] {
        let bits = match norm(b) {
            b'A' => 0u64,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        code = (code << 2) | bits;
    }
    Some(code)
}

struct ChainState {
    i_start: usize,
    j_start: usize,
    i_end: usize,
    j_end: usize,
}

/// Find candidate inverted-repeat pairs in `s`, sorted by descending
/// extended match length. Only non-self-overlapping pairs are returned.
pub fn find_inverted_repeats(s: &[u8], cfg: &IrFinderConfig) -> Vec<IrPair> {
    let n = s.len();
    if n < cfg.kmer * 4 {
        return Vec::new();
    }
    let t = revcomp(s);

    // Index every ACGT k-mer of s -> positions.
    let mut index: HashMap<u64, Vec<u32>> = HashMap::with_capacity(n);
    for i in 0..=n - cfg.kmer {
        if let Some(code) = pack_kmer(s, i, cfg.kmer) {
            index.entry(code).or_default().push(i as u32);
        }
    }

    // Seed + chain along diagonals (diag = i - j) in a single left-to-right
    // pass over t. HashMap key is the diagonal.
    let mut active: HashMap<i64, ChainState> = HashMap::new();
    let mut completed: Vec<ChainState> = Vec::new();

    for j in 0..=n - cfg.kmer {
        let Some(code) = pack_kmer(&t, j, cfg.kmer) else {
            continue;
        };
        let Some(positions) = index.get(&code) else {
            continue;
        };
        for &i_u32 in positions {
            let i = i_u32 as usize;
            let diag = i as i64 - j as i64;
            match active.get_mut(&diag) {
                Some(state) if j.saturating_sub(state.j_end) <= cfg.max_seed_gap => {
                    state.i_end = i + cfg.kmer;
                    state.j_end = j + cfg.kmer;
                }
                _ => {
                    if let Some(prev) = active.remove(&diag) {
                        completed.push(prev);
                    }
                    active.insert(
                        diag,
                        ChainState {
                            i_start: i,
                            j_start: j,
                            i_end: i + cfg.kmer,
                            j_end: j + cfg.kmer,
                        },
                    );
                }
            }
        }
    }
    completed.extend(active.into_values());

    // Extend each chain with X-drop, dropping degenerate seed chains early.
    let mut candidates: Vec<IrPair> = Vec::new();
    for chain in completed {
        if chain.j_end - chain.j_start < cfg.kmer {
            continue;
        }
        let (i_start, i_end, j_start, j_end) = extend_chain(
            s,
            &t,
            chain.i_start,
            chain.i_end,
            chain.j_start,
            chain.j_end,
            cfg.x_drop,
        );
        let len = i_end - i_start;
        if len < cfg.min_len || len != j_end - j_start {
            continue;
        }
        // Map the t-region back to s-coordinates: region1 = revcomp source.
        let region1_start = n - j_end;
        let region1_end = n - j_start;
        let region2_start = i_start;
        let region2_end = i_end;
        if region1_end > region2_start && region2_end > region1_start {
            // Overlaps itself in original coordinates -> not a valid IR pair.
            continue;
        }
        let (r1s, r1e, r2s, r2e) = if region1_start <= region2_start {
            (region1_start, region1_end, region2_start, region2_end)
        } else {
            (region2_start, region2_end, region1_start, region1_end)
        };
        let (matches, mismatches) = count_identity(s, &t, i_start, j_start, len);
        candidates.push(IrPair {
            region1_start: r1s,
            region1_end: r1e,
            region2_start: r2s,
            region2_end: r2e,
            matches,
            mismatches,
        });
    }

    candidates.sort_by(|a, b| b.region1_len().cmp(&a.region1_len()));
    dedup_overlapping(candidates)
}

/// Greedy X-drop extension of a seed-chain match `s[i_start..i_end)` vs
/// `t[j_start..j_end)` in both directions.
fn extend_chain(
    s: &[u8],
    t: &[u8],
    mut i_start: usize,
    mut i_end: usize,
    mut j_start: usize,
    mut j_end: usize,
    x_drop: i32,
) -> (usize, usize, usize, usize) {
    // Extend right.
    {
        let mut score: i32 = 0;
        let mut best_score: i32 = 0;
        let mut best_off: usize = 0;
        let mut off: usize = 0;
        while i_end + off < s.len() && j_end + off < t.len() {
            score += if bases_match(s[i_end + off], t[j_end + off]) {
                1
            } else {
                -2
            };
            off += 1;
            if score > best_score {
                best_score = score;
                best_off = off;
            }
            if best_score - score > x_drop {
                break;
            }
        }
        i_end += best_off;
        j_end += best_off;
    }
    // Extend left.
    {
        let mut score: i32 = 0;
        let mut best_score: i32 = 0;
        let mut best_off: usize = 0;
        let mut off: usize = 0;
        while i_start > off && j_start > off {
            let (si, tj) = (s[i_start - off - 1], t[j_start - off - 1]);
            score += if bases_match(si, tj) { 1 } else { -2 };
            off += 1;
            if score > best_score {
                best_score = score;
                best_off = off;
            }
            if best_score - score > x_drop {
                break;
            }
        }
        i_start -= best_off;
        j_start -= best_off;
    }
    (i_start, i_end, j_start, j_end)
}

fn count_identity(
    s: &[u8],
    t: &[u8],
    i_start: usize,
    j_start: usize,
    len: usize,
) -> (usize, usize) {
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    for k in 0..len {
        if bases_match(s[i_start + k], t[j_start + k]) {
            matches += 1;
        } else {
            mismatches += 1;
        }
    }
    (matches, mismatches)
}

/// Keep candidates in descending length order, dropping any whose region1 or
/// region2 interval overlaps an already-kept candidate's intervals.
fn dedup_overlapping(candidates: Vec<IrPair>) -> Vec<IrPair> {
    let mut kept: Vec<IrPair> = Vec::new();
    'outer: for c in candidates {
        for k in &kept {
            let overlaps = |a0: usize, a1: usize, b0: usize, b1: usize| a0 < b1 && b0 < a1;
            if overlaps(
                c.region1_start,
                c.region1_end,
                k.region1_start,
                k.region1_end,
            ) || overlaps(
                c.region1_start,
                c.region1_end,
                k.region2_start,
                k.region2_end,
            ) || overlaps(
                c.region2_start,
                c.region2_end,
                k.region1_start,
                k.region1_end,
            ) || overlaps(
                c.region2_start,
                c.region2_end,
                k.region2_start,
                k.region2_end,
            ) {
                continue 'outer;
            }
        }
        kept.push(c);
    }
    kept
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::revcomp as rc;

    fn synthetic_genome(lsc_len: usize, ir_len: usize, ssc_len: usize, seed: u64) -> Vec<u8> {
        let mut rng = seed;
        let mut next_base = || {
            // xorshift PRNG, deterministic, dependency-free
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            match rng % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }
        };
        let lsc: Vec<u8> = (0..lsc_len).map(|_| next_base()).collect();
        let ir: Vec<u8> = (0..ir_len).map(|_| next_base()).collect();
        let ssc: Vec<u8> = (0..ssc_len).map(|_| next_base()).collect();
        let ir_rc = rc(&ir);
        let mut genome = Vec::new();
        genome.extend_from_slice(&lsc);
        genome.extend_from_slice(&ir);
        genome.extend_from_slice(&ssc);
        genome.extend_from_slice(&ir_rc);
        genome
    }

    #[test]
    fn finds_ir_in_synthetic_genome() {
        let genome = synthetic_genome(8000, 2000, 3000, 12345);
        let cfg = IrFinderConfig {
            min_len: 500,
            ..Default::default()
        };
        let hits = find_inverted_repeats(&genome, &cfg);
        assert!(!hits.is_empty(), "should find at least one IR candidate");
        let best = &hits[0];
        assert!(best.region1_len() >= 1900 && best.region1_len() <= 2000 + 5);
        assert!(best.identity() > 0.99);
        // region1 should start right after LSC (position 8000), region2 after SSC.
        assert!((best.region1_start as i64 - 8000).abs() <= 5);
    }

    #[test]
    fn tolerates_mismatches_between_copies() {
        let mut genome = synthetic_genome(5000, 4000, 2500, 999);
        // Introduce a few point mutations into the second IR copy.
        let n = genome.len();
        let ir2_start = n - 4000;
        for k in [10usize, 500, 1500, 3000] {
            let pos = ir2_start + k;
            genome[pos] = match genome[pos] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        let cfg = IrFinderConfig {
            min_len: 500,
            ..Default::default()
        };
        let hits = find_inverted_repeats(&genome, &cfg);
        assert!(!hits.is_empty());
        assert!(hits[0].region1_len() > 3900);
        assert!(hits[0].identity() > 0.98);
    }
}
