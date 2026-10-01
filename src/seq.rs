//! Minimal DNA sequence helpers (no external dependencies).

/*
Gaurav Sablok
gsablok@proton.me
*/

/// Complement a single IUPAC nucleotide byte, preserving case.
#[inline]
pub fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'U' => b'A',
        b'u' => b'a',
        b'N' => b'N',
        b'n' => b'n',
        // IUPAC ambiguity codes
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'r' => b'y',
        b'y' => b'r',
        b's' => b's',
        b'w' => b'w',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        other => other,
    }
}

/// Reverse-complement a sequence.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// Rotate a sequence so that it starts at `offset` (circular rotation).
/// `offset` must be < seq.len() (or seq is returned unchanged if seq is empty).
pub fn rotate(seq: &[u8], offset: usize) -> Vec<u8> {
    if seq.is_empty() {
        return Vec::new();
    }
    let offset = offset % seq.len();
    let mut out = Vec::with_capacity(seq.len());
    out.extend_from_slice(&seq[offset..]);
    out.extend_from_slice(&seq[..offset]);
    out
}

/// Uppercase-normalize a base for identity comparisons.
#[inline]
pub fn norm(b: u8) -> u8 {
    b.to_ascii_uppercase()
}

/// True if two normalized bases should count as a "match" for scoring purposes.
/// Ambiguous/N bases never count as a match (conservative).
#[inline]
pub fn bases_match(a: u8, b: u8) -> bool {
    let (a, b) = (norm(a), norm(b));
    matches!(a, b'A' | b'C' | b'G' | b'T') && a == b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT".to_vec()); // palindrome
        assert_eq!(revcomp(b"AACCGGTT"), b"AACCGGTT".to_vec());
        assert_eq!(revcomp(b"ATGC"), b"GCAT".to_vec());
        assert_eq!(revcomp(b"NNNACGT"), b"ACGTNNN".to_vec());
    }

    #[test]
    fn test_rotate() {
        assert_eq!(rotate(b"ABCDEF", 2), b"CDEFAB".to_vec());
        assert_eq!(rotate(b"ABCDEF", 0), b"ABCDEF".to_vec());
        assert_eq!(rotate(b"ABCDEF", 6), b"ABCDEF".to_vec());
    }
}
