//! Frequency filtering of a k-mer count table.
//!
//! The standard QC step before assembly: drop low-frequency k-mers, which
//! in real read data are overwhelmingly sequencing errors (the frequency-1
//! peak of the k-mer spectrum; see `spectrum.rs`). Jellyfish and KMC
//! expose the same idea as a minimum-count cutoff (`jellyfish count -L`,
//! `kmc -ci`).

use crate::counter::KmerCounter;
use std::collections::HashMap;

/// Keep only k-mers with `count >= min_count`. Metadata (total, skipped,
/// per-record summaries) is adjusted to the filtered table.
pub fn filter_min(counter: &KmerCounter, min_count: u64) -> KmerCounter {
    let counts: HashMap<u64, u64> = counter
        .counts
        .iter()
        .filter(|(_, &n)| n >= min_count)
        .map(|(&code, &n)| (code, n))
        .collect();
    let total: u64 = counts.values().sum();
    KmerCounter {
        k: counter.k,
        canonical: counter.canonical,
        counts,
        total,
        skipped_windows: counter.skipped_windows,
        total_bases: counter.total_bases,
        records: counter.records,
        per_record: counter.per_record.clone(),
    }
}

/// Keep only the `n` most frequent k-mers (ties broken by k-mer code,
/// same order as [`KmerCounter::sorted_counts`]).
pub fn retain_top(counter: &KmerCounter, n: usize) -> KmerCounter {
    let counts: HashMap<u64, u64> = counter
        .sorted_counts()
        .into_iter()
        .take(n)
        .collect();
    let total: u64 = counts.values().sum();
    KmerCounter {
        k: counter.k,
        canonical: counter.canonical,
        counts,
        total,
        skipped_windows: counter.skipped_windows,
        total_bases: counter.total_bases,
        records: counter.records,
        per_record: counter.per_record.clone(),
    }
}

/// "Solid" k-mers: those seen at least twice. In error-correction
/// pipelines these are the k-mers trusted as real sequence.
pub fn solid_kmers(counter: &KmerCounter) -> KmerCounter {
    filter_min(counter, 2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counter::KmerCounter;

    fn toy() -> KmerCounter {
        let mut c = KmerCounter::new(3, false).unwrap();
        c.count_sequence(b"ACGTACGT");
        c
    }

    #[test]
    fn filter_min_drops_singletons() {
        // Toy counts 2,2,1,1 -> min 2 keeps ACG, CGT only.
        let f = filter_min(&toy(), 2);
        assert_eq!(f.distinct(), 2);
        assert_eq!(f.total, 4);
        assert_eq!(f.get_str("ACG"), 2);
        assert_eq!(f.get_str("GTA"), 0);
    }

    #[test]
    fn filter_min_one_keeps_everything() {
        assert_eq!(filter_min(&toy(), 1).distinct(), 4);
    }

    #[test]
    fn solid_kmers_equals_min_two() {
        assert_eq!(solid_kmers(&toy()).counts, filter_min(&toy(), 2).counts);
    }

    #[test]
    fn retain_top_keeps_most_frequent() {
        let f = retain_top(&toy(), 2);
        assert_eq!(f.distinct(), 2);
        assert_eq!(f.total, 4);
        // Ties at count 2 resolve by code ascending: ACG, CGT.
        assert_eq!(f.get_str("ACG"), 2);
        assert_eq!(f.get_str("CGT"), 2);
    }
}
