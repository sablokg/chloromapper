//! Heuristics for classifying multiple input contigs (as opposed to one
//! already-complete circular contig) into LSC / IR / SSC roles.
//!
//! Two supported scenarios:
//!  1. The IR was assembled as two separate, near-mirror-image contigs.
//!     Detected by pairwise seed-matching every contig against every other
//!     contig's reverse complement.
//!  2. The assembly graph collapsed the repeat into a single IR contig
//!     (very common in practice). With exactly 3 contigs left we fall back
//!     to a length-based heuristic (LSC longest, SSC shortest, IR the
//!     middle one) which holds for the great majority of land-plant
//!     plastomes but is not guaranteed — callers should let users override
//!     it explicitly via CLI flags.

use crate::ir_finder::IrFinderConfig;
use std::collections::HashMap;

/*
Gaurav Sablok
gsablok@proton.me
*/

pub struct Classified<'a> {
    pub lsc: &'a [u8],
    pub lsc_id: &'a str,
    pub ir: &'a [u8],
    pub ir_id: &'a str,
    pub ssc: &'a [u8],
    pub ssc_id: &'a str,
    /// Second IR contig, if the repeat was assembled as two separate
    /// contigs (rather than collapsed into one by the assembly graph).
    pub ir_second: Option<Vec<u8>>,
    pub ir_second_id: Option<&'a str>,
    pub used_length_heuristic: bool,
}

/// 2-bit-pack an ACGT k-mer; `None` on any ambiguous base or out-of-range window.
fn pack(seq: &[u8], pos: usize, k: usize) -> Option<u64> {
    if pos + k > seq.len() {
        return None;
    }
    let mut code = 0u64;
    for &b in &seq[pos..pos + k] {
        let bits = match b.to_ascii_uppercase() {
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

/// Longest approximately-collinear seed run between `a` and `b` (both read
/// forward), found by k-mer seeding + diagonal grouping with a small allowed
/// gap between consecutive seeds on the same diagonal. This is a coarse,
/// gap-tolerant "how much of a and b look like the same sequence" measure —
/// good enough to decide whether two contigs are the two IR copies.
fn best_collinear_run(a: &[u8], b: &[u8], k: usize, max_gap: usize) -> usize {
    if a.len() < k || b.len() < k {
        return 0;
    }
    let mut index: HashMap<u64, Vec<u32>> = HashMap::new();
    for i in 0..=a.len() - k {
        if let Some(code) = pack(a, i, k) {
            index.entry(code).or_default().push(i as u32);
        }
    }

    // diag -> (run_start_j, run_end_j)
    let mut runs: HashMap<i64, (usize, usize)> = HashMap::new();
    let mut best = 0usize;
    for j in 0..=b.len() - k {
        let Some(code) = pack(b, j, k) else { continue };
        let Some(positions) = index.get(&code) else {
            continue;
        };
        for &i_u32 in positions {
            let i = i_u32 as usize;
            let diag = i as i64 - j as i64;
            let run = runs.entry(diag).or_insert((j, j + k));
            if j > run.1 + max_gap {
                // Gap too large: start a fresh run on this diagonal.
                *run = (j, j + k);
            } else {
                run.1 = (j + k).max(run.1);
            }
            best = best.max(run.1 - run.0);
        }
    }
    best
}

/// Try to find two contigs that are (near-)reverse-complements of each
/// other across most of their length: the "IR assembled as two separate
/// contigs" case.
fn find_matching_ir_pair(
    contigs: &[(&str, &[u8])],
    cfg: &IrFinderConfig,
) -> Option<(usize, usize)> {
    use crate::seq::revcomp;
    for i in 0..contigs.len() {
        for j in (i + 1)..contigs.len() {
            let (_, seq_i) = contigs[i];
            let (_, seq_j) = contigs[j];
            let rc_j = revcomp(seq_j);
            let run = best_collinear_run(seq_i, &rc_j, cfg.kmer.min(17), 200);
            let shorter = seq_i.len().min(rc_j.len());
            if shorter > 0 && run as f64 >= 0.9 * shorter as f64 && run >= cfg.min_len.min(shorter)
            {
                return Some((i, j));
            }
        }
    }
    None
}

/// Classify a set of named contigs. `contigs` must have length >= 2.
pub fn classify<'a>(
    contigs: &[(&'a str, &'a [u8])],
    cfg: &IrFinderConfig,
) -> Result<Classified<'a>, String> {
    if contigs.len() < 2 {
        return Err("classify() requires at least 2 contigs".to_string());
    }

    if let Some((i, j)) = find_matching_ir_pair(contigs, cfg) {
        let mut rest: Vec<usize> = (0..contigs.len()).filter(|&k| k != i && k != j).collect();
        if rest.len() < 2 {
            return Err(format!(
                "found an IR pair ({} / {}) but fewer than 2 remaining contigs for LSC/SSC",
                contigs[i].0, contigs[j].0
            ));
        }
        rest.sort_by_key(|&k| std::cmp::Reverse(contigs[k].1.len()));
        let (lsc_idx, ssc_idx) = (rest[0], rest[1]);
        return Ok(Classified {
            lsc: contigs[lsc_idx].1,
            lsc_id: contigs[lsc_idx].0,
            ir: contigs[i].1,
            ir_id: contigs[i].0,
            ssc: contigs[ssc_idx].1,
            ssc_id: contigs[ssc_idx].0,
            ir_second: Some(contigs[j].1.to_vec()),
            ir_second_id: Some(contigs[j].0),
            used_length_heuristic: false,
        });
    }

    // Fallback: no matching IR pair found -> assume the repeat was
    // collapsed into a single contig. Requires exactly 3 contigs to apply
    // the length heuristic unambiguously.
    if contigs.len() != 3 {
        return Err(format!(
            "no inverted-repeat contig pair found, and {} contigs were given (expected exactly 3 \
             for the collapsed-repeat heuristic: LSC, IR, SSC). Re-run with --lsc/--ir/--ssc to \
             specify roles explicitly.",
            contigs.len()
        ));
    }
    let mut idx: Vec<usize> = (0..3).collect();
    idx.sort_by_key(|&k| std::cmp::Reverse(contigs[k].1.len()));
    let (lsc_idx, ir_idx, ssc_idx) = (idx[0], idx[1], idx[2]);
    Ok(Classified {
        lsc: contigs[lsc_idx].1,
        lsc_id: contigs[lsc_idx].0,
        ir: contigs[ir_idx].1,
        ir_id: contigs[ir_idx].0,
        ssc: contigs[ssc_idx].1,
        ssc_id: contigs[ssc_idx].0,
        ir_second: None,
        ir_second_id: None,
        used_length_heuristic: true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::revcomp;

    fn dna(seed: u64, len: usize) -> Vec<u8> {
        let mut rng = seed;
        (0..len)
            .map(|_| {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                match rng % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                }
            })
            .collect()
    }

    #[test]
    fn classifies_two_separate_ir_contigs() {
        let lsc = dna(1, 6000);
        let ir = dna(2, 2200);
        let ssc = dna(3, 3000);
        let ir2 = revcomp(&ir);
        let contigs: Vec<(&str, &[u8])> = vec![
            ("lsc_contig", &lsc),
            ("ssc_contig", &ssc),
            ("irA", &ir),
            ("irB", &ir2),
        ];
        let cfg = IrFinderConfig {
            min_len: 500,
            ..Default::default()
        };
        let c = classify(&contigs, &cfg).expect("should classify");
        assert!(!c.used_length_heuristic);
        assert_eq!(c.lsc_id, "lsc_contig");
        assert_eq!(c.ssc_id, "ssc_contig");
        assert!(c.ir_id == "irA" || c.ir_id == "irB");
        assert!(c.ir_second.is_some());
    }

    #[test]
    fn falls_back_to_length_heuristic_for_collapsed_repeat() {
        let lsc = dna(10, 9000);
        let ir = dna(11, 2500);
        let ssc = dna(12, 1800);
        let contigs: Vec<(&str, &[u8])> = vec![("c1", &ssc), ("c2", &lsc), ("c3", &ir)];
        let cfg = IrFinderConfig {
            min_len: 500,
            ..Default::default()
        };
        let c = classify(&contigs, &cfg).expect("should classify via heuristic");
        assert!(c.used_length_heuristic);
        assert_eq!(c.lsc_id, "c2");
        assert_eq!(c.ir_id, "c3");
        assert_eq!(c.ssc_id, "c1");
    }
}
