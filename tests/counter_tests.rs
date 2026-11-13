//! Integration tests for counting + spectrum + stats on hand-verified
//! sequences (public API).

use kmer_counter::counter::KmerCounter;
use kmer_counter::fasta::parse_fasta;
use kmer_counter::spectrum::Spectrum;
use kmer_counter::stats::summarize;

fn count(seq: &str, k: usize, canonical: bool) -> KmerCounter {
    let mut c = KmerCounter::new(k, canonical).unwrap();
    c.count_sequence(seq.as_bytes());
    c
}

#[test]
fn toy_sequence_full_pipeline() {
    let c = count("ACGTACGT", 3, false);
    let s = Spectrum::from_counter(&c);
    assert_eq!(s.total_kmers(), 6);
    assert_eq!(s.singletons(), 2);
    let summary = summarize(&c, 10);
    assert_eq!(summary.distinct_kmers, 4);
    assert!((summary.entropy_bits - 1.9182958340544896).abs() < 1e-9);
    assert!((summary.gc.mean - 0.5).abs() < 1e-12);
}

#[test]
fn gattaca_repeat_k5_hand_counted() {
    // "GATTACA" x 17 + "G" = 120 bp, k=5, non-canonical:
    // 7 distinct 5-mers; GATTA/ATTAC/TTACA/TACAG x17, ACAGA/CAGAT/AGATT x16.
    let seq = "GATTACA".repeat(17) + "G";
    assert_eq!(seq.len(), 120);
    let c = count(&seq, 5, false);
    assert_eq!(c.total, 116);
    assert_eq!(c.distinct(), 7);
    for kmer in ["GATTA", "ATTAC", "TTACA", "TACAG"] {
        assert_eq!(c.get_str(kmer), 17, "{kmer}");
    }
    for kmer in ["ACAGA", "CAGAT", "AGATT"] {
        assert_eq!(c.get_str(kmer), 16, "{kmer}");
    }
    // Canonical: the 7 k-mers pair up with reverse complements outside the
    // set, so the table still has 7 entries, stored under the smaller code
    // of each pair: GATTA and ATTAC and ACAGA stay; TACAG -> CTGTA,
    // TTACA -> TGTAA, CAGAT -> ATCTG, AGATT -> AATCT. `get_str` canonicalises
    // the query, so lookups by the forward form still hit the right bucket.
    let cc = count(&seq, 5, true);
    assert_eq!(cc.distinct(), 7);
    assert_eq!(cc.total, 116);
    assert_eq!(cc.get_str("TACAG"), 17); // same bucket as CTGTA
    assert_eq!(cc.get_str("CTGTA"), 17);
    assert_eq!(cc.get_str("GATTA"), 17);
    assert_eq!(cc.get_str("TTACA"), 17); // same bucket as TGTAA
    assert_eq!(cc.get_str("CAGAT"), 16); // same bucket as ATCTG
}

#[test]
fn multi_record_total_matches_sum_of_parts() {
    let recs = parse_fasta(">a\nACGTACGT\n>b\nACGTACGT\n").unwrap();
    let mut c = KmerCounter::new(3, false).unwrap();
    c.count_records(&recs);
    assert_eq!(c.total, 12);
    assert_eq!(c.get_str("ACG"), 4);
    assert_eq!(c.records, 2);
}

#[test]
fn counting_is_case_insensitive() {
    let upper = count("ACGTACGT", 3, false);
    let lower = count("acgtacgt", 3, false);
    assert_eq!(upper.counts, lower.counts);
}
