//! Integration tests for the k-mer encoding layer (public API).

use kmer_counter::kmer::*;
use std::collections::HashSet;

#[test]
fn all_3mers_encode_decode() {
    // All 4^3 = 64 3-mers round-trip through the packed encoding.
    let bases = *b"ACGT";
    let mut seen = HashSet::new();
    for &a in &bases {
        for &b in &bases {
            for &c in &bases {
                let s = [a, b, c];
                let code = encode_kmer(&s, 3).unwrap();
                assert!(seen.insert(code), "duplicate code for {:?}", s);
                assert_eq!(decode_kmer(code, 3).as_bytes(), s);
            }
        }
    }
    assert_eq!(seen.len(), 64);
}

#[test]
fn canonical_pairs_collapse() {
    // For every 3-mer, canonical(kmer) == canonical(reverse complement).
    let bases = *b"ACGT";
    for &a in &bases {
        for &b in &bases {
            for &c in &bases {
                let s = [a, b, c];
                let code = encode_kmer(&s, 3).unwrap();
                let rc = reverse_complement_code(code, 3);
                assert_eq!(canonical_code(code, 3), canonical_code(rc, 3));
            }
        }
    }
}

#[test]
fn k32_boundary_works() {
    // A 32-mer of all A is code 0; all T is u64::MAX; rc swaps them.
    let all_a = vec![b'A'; 32];
    let all_t = vec![b'T'; 32];
    let ca = encode_kmer(&all_a, 32).unwrap();
    let ct = encode_kmer(&all_t, 32).unwrap();
    assert_eq!(ca, 0);
    assert_eq!(ct, u64::MAX);
    assert_eq!(reverse_complement_code(ca, 32), ct);
    assert_eq!(reverse_complement_code(ct, 32), ca);
    assert!(validate_k(32).is_ok());
}
