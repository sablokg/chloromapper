//! Turn detected regions into a single, canonically-ordered, end-to-end
//! chloroplast genome sequence (LSC - IR - SSC - IR), filling `N` at any
//! junction whose exact length is genuinely unknown (only arises when
//! stitching separately-assembled contigs together).

use crate::ir_finder::IrPair;
use crate::seq::{revcomp, rotate};

/*
Gaurav Sablok
gsablok@proton.me
*/

/// Coordinates (0-based, half-open, into the *output* sequence) of each of
/// the four canonical regions, for reporting/annotation.
#[derive(Debug, Clone, Copy)]
pub struct RegionMap {
    pub lsc: (usize, usize),
    pub ir_b: (usize, usize),
    pub ssc: (usize, usize),
    pub ir_a: (usize, usize),
}

pub struct ArrangedGenome {
    pub seq: Vec<u8>,
    pub regions: RegionMap,
    /// True if `ir_a` was synthesized as the exact reverse complement of
    /// `ir_b` rather than taken from directly-observed sequence (only
    /// happens in multi-contig mode when the second IR copy is missing).
    pub ir_a_synthesized: bool,
}

/// Single-contig case: `s` is one already-circular-assembled sequence
/// containing both IR copies. We only need to *rotate* it (never
/// reverse-complement any piece) so it reads LSC -> IR(b) -> SSC -> IR(a),
/// because a correctly assembled single contig is already self-consistent
/// on one strand.
pub fn arrange_single_contig(s: &[u8], ir: &IrPair) -> ArrangedGenome {
    let n = s.len();
    let (a, b, c, d) = (
        ir.region1_start,
        ir.region1_end,
        ir.region2_start,
        ir.region2_end,
    );
    let gap_mid = c - b; // S[b..c): the "inner" single-copy segment
    let gap_wrap = (n - d) + a; // S[d..n) + S[0..a): the "wrap" single-copy segment

    let (rotation_start, lsc_len, ir_len_first): (usize, usize, usize);
    if gap_mid >= gap_wrap {
        // LSC is the inner segment; rotate so LSC starts at b.
        rotation_start = b;
        lsc_len = gap_mid;
        ir_len_first = d - c;
    } else {
        // LSC is the wrap-around segment; rotate so LSC starts at d.
        rotation_start = d;
        lsc_len = gap_wrap;
        ir_len_first = b - a;
    }

    let rotated = rotate(s, rotation_start);
    let ssc_len = n - lsc_len - (d - c) - (b - a);

    let lsc = (0, lsc_len);
    let ir_b = (lsc_len, lsc_len + ir_len_first);
    let ssc = (ir_b.1, ir_b.1 + ssc_len);
    let ir_a = (ssc.1, n);

    debug_assert_eq!(ir_a.1 - ir_a.0, n - lsc_len - ir_len_first - ssc_len);

    ArrangedGenome {
        seq: rotated,
        regions: RegionMap {
            lsc,
            ir_b,
            ssc,
            ir_a,
        },
        ir_a_synthesized: false,
    }
}

/// Multi-contig case: contigs have already been *classified* by the caller
/// into lsc / ir / ssc (and optionally a second, independently-assembled IR
/// copy). Any contig the caller marks as needing orientation flipped should
/// already have been reverse-complemented before calling this function.
///
/// `ir_second`: `Some(seq)` if a second IR contig was found and matched as
/// the inverted-repeat partner of `ir`; `None` if the assembly only
/// contains one IR copy (common when the graph collapses the repeat), in
/// which case it is synthesized as `revcomp(ir)`.
///
/// `gap_n`: number of `N` inserted at each junction between independently
/// assembled contigs, since the exact junction sequence/length is not
/// known from the contigs alone. Pass `0` to abut contigs directly.
pub fn arrange_from_contigs(
    lsc: &[u8],
    ir: &[u8],
    ssc: &[u8],
    ir_second: Option<&[u8]>,
    gap_n: usize,
) -> ArrangedGenome {
    let (ir_a_seq, synthesized): (Vec<u8>, bool) = match ir_second {
        Some(seq) => (seq.to_vec(), false),
        None => (revcomp(ir), true),
    };
    let gap = vec![b'N'; gap_n];

    let mut out = Vec::with_capacity(lsc.len() + ir.len() + ssc.len() + ir_a_seq.len() + gap_n * 4);
    let lsc_start = out.len();
    out.extend_from_slice(lsc);
    if gap_n > 0 {
        out.extend_from_slice(&gap);
    }
    let ir_b_start = out.len();
    out.extend_from_slice(ir);
    if gap_n > 0 {
        out.extend_from_slice(&gap);
    }
    let ssc_start = out.len();
    out.extend_from_slice(ssc);
    if gap_n > 0 {
        out.extend_from_slice(&gap);
    }
    let ir_a_start = out.len();
    out.extend_from_slice(&ir_a_seq);
    let end = out.len();

    ArrangedGenome {
        seq: out,
        regions: RegionMap {
            lsc: (
                lsc_start,
                ir_b_start.saturating_sub(if gap_n > 0 { gap_n } else { 0 }),
            ),
            ir_b: (
                ir_b_start,
                ssc_start.saturating_sub(if gap_n > 0 { gap_n } else { 0 }),
            ),
            ssc: (
                ssc_start,
                ir_a_start.saturating_sub(if gap_n > 0 { gap_n } else { 0 }),
            ),
            ir_a: (ir_a_start, end),
        },
        ir_a_synthesized: synthesized,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir_finder::{find_inverted_repeats, IrFinderConfig};

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
    fn rotation_reproduces_canonical_layout() {
        // Build LSC+IR+SSC+revcomp(IR), then rotate the *test input* itself
        // to simulate an assembler that started the contig mid-genome, and
        // check we recover a consistent canonical layout.
        let lsc = dna(42, 600);
        let ir = dna(43, 300);
        let ssc = dna(44, 400);
        let mut genome = Vec::new();
        genome.extend_from_slice(&lsc);
        genome.extend_from_slice(&ir);
        genome.extend_from_slice(&ssc);
        genome.extend_from_slice(&revcomp(&ir));

        let shifted = rotate(&genome, 137); // simulate arbitrary assembler start point

        let cfg = IrFinderConfig {
            min_len: 200,
            kmer: 15,
            ..Default::default()
        };
        let hits = find_inverted_repeats(&shifted, &cfg);
        assert!(!hits.is_empty());
        let arranged = arrange_single_contig(&shifted, &hits[0]);

        assert_eq!(arranged.regions.lsc.1 - arranged.regions.lsc.0, lsc.len());
        assert_eq!(arranged.regions.ssc.1 - arranged.regions.ssc.0, ssc.len());
        assert_eq!(
            &arranged.seq[arranged.regions.ir_a.0..arranged.regions.ir_a.1],
            &revcomp(&arranged.seq[arranged.regions.ir_b.0..arranged.regions.ir_b.1])[..]
        );
    }

    #[test]
    fn multi_contig_synthesizes_missing_ir_copy() {
        let lsc = b"AAAACCCCGGGGTTTT".to_vec();
        let ir = b"ACGTACGTACGTACGTGGGGCCCCAAAATTTT".to_vec();
        let ssc = b"TTTTGGGGCCCCAAAA".to_vec();
        let arranged = arrange_from_contigs(&lsc, &ir, &ssc, None, 5);
        assert!(arranged.ir_a_synthesized);
        assert_eq!(
            &arranged.seq[arranged.regions.ir_a.0..arranged.regions.ir_a.1],
            &revcomp(&ir)[..]
        );
        // 3 junctions * 5 N's inserted (before ir, before ssc, before ir_a)
        let n_count = arranged.seq.iter().filter(|&&b| b == b'N').count();
        assert_eq!(n_count, 15);
    }
}
