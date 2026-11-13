//! Integration tests for sample comparison (public API).

use kmer_counter::compare::{bray_curtis, compare, cosine, jaccard};
use kmer_counter::counter::KmerCounter;

fn count(seq: &str, k: usize, canonical: bool) -> KmerCounter {
    let mut c = KmerCounter::new(k, canonical).unwrap();
    c.count_sequence(seq.as_bytes());
    c
}

#[test]
fn jaccard_symmetry_and_bounds() {
    let a = count("ACGTACGT", 3, true);
    let b = count("ACGTACGA", 3, true);
    assert!((jaccard(&a, &b) - jaccard(&b, &a)).abs() < 1e-12);
    assert!((0.0..=1.0).contains(&jaccard(&a, &b)));
    assert!((0.0..=1.0).contains(&bray_curtis(&a, &b)));
    assert!((0.0..=1.0).contains(&cosine(&a, &b)));
}

#[test]
fn empty_vs_empty_is_identical() {
    let a = KmerCounter::new(3, true).unwrap();
    let b = KmerCounter::new(3, true).unwrap();
    let cmp = compare(&a, &b).unwrap();
    assert!((cmp.jaccard - 1.0).abs() < 1e-12);
    assert!((cmp.bray_curtis - 0.0).abs() < 1e-12);
}

#[test]
fn nested_sample_has_containment_one() {
    // B's k-mers are a subset of A's: containment of B in A is 1.
    let a = count("ACGTACGT", 3, true);
    let b = count("ACGT", 3, true);
    let cmp = compare(&a, &b).unwrap();
    assert_eq!(cmp.distinct_b, 1); // only ACG (canonical) in "ACGT"
    assert!((cmp.containment_b - 1.0).abs() < 1e-12);
    assert!(cmp.containment_a < 1.0);
}
