//! Integration tests for spectrum + stats + filter (public API).

use kmer_counter::counter::KmerCounter;
use kmer_counter::filter::{filter_min, retain_top, solid_kmers};
use kmer_counter::spectrum::Spectrum;
use kmer_counter::stats::summarize;

fn toy() -> KmerCounter {
    let mut c = KmerCounter::new(3, false).unwrap();
    c.count_sequence(b"ACGTACGT");
    c
}

#[test]
fn spectrum_and_filter_agree() {
    let s = Spectrum::from_counter(&toy());
    assert_eq!(s.singletons(), 2);
    let solid = solid_kmers(&toy());
    assert_eq!(solid.distinct(), 2);
    assert_eq!(Spectrum::from_counter(&solid).singletons(), 0);
    assert_eq!(Spectrum::from_counter(&solid).total_kmers(), 4);
}

#[test]
fn filter_min_beyond_max_empties_table() {
    let f = filter_min(&toy(), 99);
    assert_eq!(f.distinct(), 0);
    assert_eq!(f.total, 0);
    let s = summarize(&f, 5);
    assert_eq!(s.entropy_bits, 0.0);
    assert!(s.top.is_empty());
}

#[test]
fn retain_top_zero_empties_table() {
    assert_eq!(retain_top(&toy(), 0).distinct(), 0);
}

#[test]
fn periodic_spectrum_is_point_mass() {
    // "ACGT" x 30, k=5: 4 distinct 5-mers, each x29 -> spectrum {29: 4}.
    let seq = "ACGT".repeat(30);
    let mut c = KmerCounter::new(5, false).unwrap();
    c.count_sequence(seq.as_bytes());
    let s = Spectrum::from_counter(&c);
    assert_eq!(s.histogram.len(), 1);
    assert_eq!(s.histogram.get(&29), Some(&4));
    assert_eq!(s.mode_frequency(), 29);
    assert_eq!(s.mean_frequency(), 29.0);
    // Naive genome-size estimate: 116 / 29 = 4 (the repeat unit, as expected
    // for a perfectly periodic toy sequence).
    assert_eq!(s.estimate_genome_size(), Some(4));
}
