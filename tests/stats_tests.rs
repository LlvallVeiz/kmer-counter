//! Integration tests for summary statistics (public API).

use kmer_counter::counter::KmerCounter;
use kmer_counter::stats::{shannon_entropy, summarize};

#[test]
fn summary_of_empty_counter_is_zeroed() {
    let c = KmerCounter::new(3, true).unwrap();
    let s = summarize(&c, 5);
    assert_eq!(s.total_kmers, 0);
    assert_eq!(s.distinct_kmers, 0);
    assert_eq!(s.entropy_bits, 0.0);
    assert!(s.top.is_empty());
    assert_eq!(s.genome_size_estimate, None);
}

#[test]
fn entropy_scales_with_alphabet_size() {
    // Uniform over 4 outcomes: 2 bits; uniform over 16: 4 bits.
    assert!((shannon_entropy(&[1, 1, 1, 1]) - 2.0).abs() < 1e-12);
    assert!((shannon_entropy(&[1; 16]) - 4.0).abs() < 1e-12);
}

#[test]
fn top_n_is_respected() {
    let mut c = KmerCounter::new(3, false).unwrap();
    c.count_sequence(b"ACGTACGT");
    let s = summarize(&c, 2);
    assert_eq!(s.top.len(), 2);
    let s0 = summarize(&c, 0);
    assert!(s0.top.is_empty());
}
