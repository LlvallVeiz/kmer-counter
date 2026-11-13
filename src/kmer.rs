//! k-mer encoding, canonicalisation and reverse complements.
//!
//! A k-mer is stored as a 2-bit packed `u64`: A = 00, C = 01, G = 10, T = 11.
//! This supports k up to 32 (2 * 32 = 64 bits). The encoding is the same
//! idea used by production counters such as Jellyfish (Marçais & Kingsford,
//! 2011), which packs k-mers into integers for fast hashing.

use crate::errors::KmerError;

/// Largest k that fits in a `u64` with 2 bits per base.
pub const MAX_K: usize = 32;

/// Validate a k-mer size: must be in 1..=32.
pub fn validate_k(k: usize) -> Result<(), KmerError> {
    if k == 0 {
        return Err(KmerError::InvalidK("k must be at least 1".to_string()));
    }
    if k > MAX_K {
        return Err(KmerError::InvalidK(format!(
            "k = {k} exceeds the maximum of {MAX_K} (2-bit packing in u64)"
        )));
    }
    Ok(())
}

/// Encode one base (case-insensitive). Returns `None` for any non-ACGT
/// character (N and other IUPAC ambiguity codes included).
pub fn encode_base(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Decode one 2-bit base value back to an ASCII base.
pub fn decode_base(value: u64) -> u8 {
    match value & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

/// Bit mask keeping only the lowest `2 * k` bits.
pub fn kmer_mask(k: usize) -> u64 {
    if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    }
}

/// Encode a full k-mer string into its 2-bit code.
/// Returns `None` if the length is not `k` or any base is ambiguous.
pub fn encode_kmer(seq: &[u8], k: usize) -> Option<u64> {
    if seq.len() != k {
        return None;
    }
    let mut code: u64 = 0;
    for &b in seq {
        code = (code << 2) | encode_base(b)?;
    }
    Some(code)
}

/// Decode a 2-bit k-mer code back to its string form.
pub fn decode_kmer(code: u64, k: usize) -> String {
    let mut out = Vec::with_capacity(k);
    for i in (0..k).rev() {
        out.push(decode_base(code >> (2 * i)));
    }
    // SAFETY: only ASCII ACGT bytes are produced.
    String::from_utf8(out).expect("decoded k-mer is always valid ASCII")
}

/// Reverse complement of a packed k-mer: complement each base (xor 0b11)
/// and reverse the order.
pub fn reverse_complement_code(code: u64, k: usize) -> u64 {
    let mut rc: u64 = 0;
    for i in 0..k {
        let base = (code >> (2 * i)) & 0b11;
        rc = (rc << 2) | (base ^ 0b11);
    }
    rc
}

/// Canonical k-mer: the lexicographically smaller of a k-mer and its
/// reverse complement (compared as packed integers, which preserves
/// lexicographic order for the A<C<G<T encoding).
pub fn canonical_code(code: u64, k: usize) -> u64 {
    let rc = reverse_complement_code(code, k);
    code.min(rc)
}

/// Reverse complement of a k-mer string (convenience wrapper).
pub fn reverse_complement_str(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'T' | 't' => 'A',
            other => other,
        })
        .collect()
}

/// Fraction of G/C bases in a packed k-mer, in [0, 1].
pub fn gc_fraction(code: u64, k: usize) -> f64 {
    if k == 0 {
        return 0.0;
    }
    let mut gc = 0usize;
    for i in 0..k {
        let base = (code >> (2 * i)) & 0b11;
        if base == 1 || base == 2 {
            gc += 1;
        }
    }
    gc as f64 / k as f64
}

/// All k-mers of a sequence, as packed codes, skipping every window that
/// contains an ambiguous base. Also returns the number of skipped windows
/// (total possible windows minus emitted windows).
pub fn sliding_kmers(seq: &[u8], k: usize) -> (Vec<u64>, usize) {
    let mask = kmer_mask(k);
    let mut out = Vec::new();
    let mut code: u64 = 0;
    let mut run = 0usize; // length of current unambiguous run
    for &b in seq {
        match encode_base(b) {
            Some(v) => {
                code = ((code << 2) | v) & mask;
                run += 1;
                if run >= k {
                    out.push(code);
                }
            }
            None => {
                // An ambiguous base breaks the window: every k-mer that
                // would span it is skipped, and the rolling code resets.
                code = 0;
                run = 0;
            }
        }
    }
    let total_windows = seq.len().saturating_sub(k - 1);
    let skipped = total_windows.saturating_sub(out.len());
    (out, skipped)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_decode_roundtrip() {
        for (s, k) in [("ACGT", 4), ("AAAA", 4), ("TTTT", 4), ("GATTACA", 7)] {
            let code = encode_kmer(s.as_bytes(), k).unwrap();
            assert_eq!(decode_kmer(code, k), s);
        }
    }

    #[test]
    fn encode_rejects_ambiguous() {
        assert_eq!(encode_kmer(b"ACN", 3), None);
        assert_eq!(encode_kmer(b"ACG", 4), None); // wrong length
        assert_eq!(encode_kmer(b"ACX", 3), None);
    }

    #[test]
    fn reverse_complement_known() {
        // rc(ACGT) = ACGT (palindromic); rc(AAA) = TTT.
        assert_eq!(decode_kmer(reverse_complement_code(encode_kmer(b"ACGT", 4).unwrap(), 4), 4), "ACGT");
        assert_eq!(decode_kmer(reverse_complement_code(encode_kmer(b"AAA", 3).unwrap(), 3), 3), "TTT");
        assert_eq!(decode_kmer(reverse_complement_code(encode_kmer(b"CGT", 3).unwrap(), 3), 3), "ACG");
        assert_eq!(reverse_complement_str("GATTACA"), "TGTAATC");
    }

    #[test]
    fn canonical_picks_smaller() {
        // CGT and ACG are reverse complements; canonical is ACG in both cases.
        let acg = encode_kmer(b"ACG", 3).unwrap();
        let cgt = encode_kmer(b"CGT", 3).unwrap();
        assert_eq!(canonical_code(acg, 3), acg);
        assert_eq!(canonical_code(cgt, 3), acg);
        // GTA and TAC: canonical is GTA.
        let gta = encode_kmer(b"GTA", 3).unwrap();
        let tac = encode_kmer(b"TAC", 3).unwrap();
        assert_eq!(canonical_code(tac, 3), gta);
    }

    #[test]
    fn gc_fraction_known() {
        assert!((gc_fraction(encode_kmer(b"ACG", 3).unwrap(), 3) - 2.0 / 3.0).abs() < 1e-12);
        assert!((gc_fraction(encode_kmer(b"AAA", 3).unwrap(), 3) - 0.0).abs() < 1e-12);
        assert!((gc_fraction(encode_kmer(b"GCGC", 4).unwrap(), 4) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn validate_k_bounds() {
        assert!(validate_k(1).is_ok());
        assert!(validate_k(32).is_ok());
        assert!(validate_k(0).is_err());
        assert!(validate_k(33).is_err());
    }

    #[test]
    fn sliding_kmers_skips_ambiguous_windows() {
        // "ACNTAC", k=3: windows ACN, CNT, NTA skipped; TAC counted.
        let (kmers, skipped) = sliding_kmers(b"ACNTAC", 3);
        assert_eq!(kmers.len(), 1);
        assert_eq!(decode_kmer(kmers[0], 3), "TAC");
        assert_eq!(skipped, 3);
    }

    #[test]
    fn sliding_kmers_toy_sequence() {
        // "ACGTACGT", k=3 -> ACG CGT GTA TAC ACG CGT.
        let (kmers, skipped) = sliding_kmers(b"ACGTACGT", 3);
        assert_eq!(skipped, 0);
        let decoded: Vec<String> = kmers.iter().map(|&c| decode_kmer(c, 3)).collect();
        assert_eq!(decoded, ["ACG", "CGT", "GTA", "TAC", "ACG", "CGT"]);
    }
}
