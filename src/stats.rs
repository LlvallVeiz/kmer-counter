//! Summary statistics over a k-mer count table: Shannon entropy of the
//! count distribution, most frequent k-mers, and GC-content statistics
//! across distinct k-mers.

use crate::counter::KmerCounter;
use crate::kmer::{decode_kmer, gc_fraction};
use crate::spectrum::Spectrum;

/// One entry of the "most frequent k-mers" list.
#[derive(Debug, Clone, PartialEq)]
pub struct TopKmer {
    /// k-mer string.
    pub kmer: String,
    /// Occurrence count.
    pub count: u64,
    /// Share of all k-mer occurrences, in [0, 1].
    pub fraction: f64,
}

/// GC-content statistics across distinct k-mers.
#[derive(Debug, Clone, PartialEq)]
pub struct GcStats {
    /// Mean GC fraction over distinct k-mers.
    pub mean: f64,
    /// Minimum GC fraction.
    pub min: f64,
    /// Maximum GC fraction.
    pub max: f64,
}

/// Full summary of a counting run.
#[derive(Debug, Clone)]
pub struct Summary {
    /// k-mer size.
    pub k: usize,
    /// Canonical mode in effect.
    pub canonical: bool,
    /// Records processed.
    pub records: usize,
    /// Total bases processed.
    pub total_bases: u64,
    /// Total k-mer windows counted.
    pub total_kmers: u64,
    /// Distinct k-mers.
    pub distinct_kmers: usize,
    /// Singleton k-mers (frequency 1).
    pub singletons: u64,
    /// Windows skipped due to ambiguous bases.
    pub skipped_windows: u64,
    /// Mean frequency (total / distinct).
    pub mean_frequency: f64,
    /// Median frequency.
    pub median_frequency: f64,
    /// Modal frequency (peak of the spectrum).
    pub mode_frequency: u64,
    /// Shannon entropy of the count distribution, in bits.
    pub entropy_bits: f64,
    /// GC statistics over distinct k-mers.
    pub gc: GcStats,
    /// Most frequent k-mers (up to `top_n`).
    pub top: Vec<TopKmer>,
    /// Naive genome-size estimate (total / modal frequency), if defined.
    pub genome_size_estimate: Option<u64>,
}

/// Shannon entropy (base 2) of a count distribution.
/// `p_i = count_i / total`; H = -sum(p_i * log2(p_i)).
pub fn shannon_entropy(counts: &[u64]) -> f64 {
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return 0.0;
    }
    counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / total as f64;
            -p * p.log2()
        })
        .sum()
}

/// GC statistics over a set of packed k-mers.
pub fn gc_stats(codes: &[u64], k: usize) -> GcStats {
    if codes.is_empty() {
        return GcStats {
            mean: 0.0,
            min: 0.0,
            max: 0.0,
        };
    }
    let mut sum = 0.0;
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &code in codes {
        let g = gc_fraction(code, k);
        sum += g;
        if g < min {
            min = g;
        }
        if g > max {
            max = g;
        }
    }
    GcStats {
        mean: sum / codes.len() as f64,
        min,
        max,
    }
}

/// Build the full summary for a counter (top `top_n` k-mers included).
pub fn summarize(counter: &KmerCounter, top_n: usize) -> Summary {
    let spectrum = Spectrum::from_counter(counter);
    let codes: Vec<u64> = counter.counts.keys().copied().collect();
    let counts: Vec<u64> = counter.counts.values().copied().collect();
    let top: Vec<TopKmer> = counter
        .sorted_counts()
        .into_iter()
        .take(top_n)
        .map(|(code, count)| TopKmer {
            kmer: decode_kmer(code, counter.k),
            count,
            fraction: if counter.total == 0 {
                0.0
            } else {
                count as f64 / counter.total as f64
            },
        })
        .collect();
    Summary {
        k: counter.k,
        canonical: counter.canonical,
        records: counter.records,
        total_bases: counter.total_bases,
        total_kmers: counter.total,
        distinct_kmers: counter.distinct(),
        singletons: spectrum.singletons(),
        skipped_windows: counter.skipped_windows,
        mean_frequency: spectrum.mean_frequency(),
        median_frequency: spectrum.median_frequency(),
        mode_frequency: spectrum.mode_frequency(),
        entropy_bits: shannon_entropy(&counts),
        gc: gc_stats(&codes, counter.k),
        top,
        genome_size_estimate: spectrum.estimate_genome_size(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counter::KmerCounter;

    #[test]
    fn entropy_uniform_is_log2_n() {
        // Four equally frequent k-mers -> H = log2(4) = 2 bits.
        assert!((shannon_entropy(&[1, 1, 1, 1]) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn entropy_single_outcome_is_zero() {
        assert!((shannon_entropy(&[10]) - 0.0).abs() < 1e-12);
        assert!((shannon_entropy(&[]) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn toy_entropy_hand_counted() {
        // Non-canonical toy counts 2,2,1,1 over total 6:
        // H = 2*(1/3)*log2(3) + 2*(1/6)*log2(6) = 1.9182958340544896 bits.
        let mut c = KmerCounter::new(3, false).unwrap();
        c.count_sequence(b"ACGTACGT");
        let s = summarize(&c, 10);
        assert!((s.entropy_bits - 1.9182958340544896).abs() < 1e-9);
        assert_eq!(s.total_kmers, 6);
        assert_eq!(s.distinct_kmers, 4);
        assert_eq!(s.singletons, 2);
        assert!((s.mean_frequency - 1.5).abs() < 1e-12);
        // Distinct k-mers ACG(2/3), CGT(2/3), GTA(1/3), TAC(1/3) -> mean GC 0.5.
        assert!((s.gc.mean - 0.5).abs() < 1e-12);
        assert!((s.gc.min - 1.0 / 3.0).abs() < 1e-12);
        assert!((s.gc.max - 2.0 / 3.0).abs() < 1e-12);
        // Top k-mers: ACG and CGT tied at 2.
        assert_eq!(s.top.len(), 4);
        assert_eq!(s.top[0].count, 2);
        assert!((s.top[0].fraction - 2.0 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn canonical_toy_entropy_hand_counted() {
        // Canonical counts 4,2 over total 6: H = 0.9182958340544896 bits.
        let mut c = KmerCounter::new(3, true).unwrap();
        c.count_sequence(b"ACGTACGT");
        let s = summarize(&c, 10);
        assert!((s.entropy_bits - 0.9182958340544896).abs() < 1e-9);
        assert_eq!(s.distinct_kmers, 2);
        assert_eq!(s.top[0].kmer, "ACG");
        assert_eq!(s.top[0].count, 4);
        assert_eq!(s.top[1].kmer, "GTA");
        assert_eq!(s.top[1].count, 2);
    }

    #[test]
    fn gc_stats_empty() {
        let g = gc_stats(&[], 3);
        assert_eq!(g, GcStats { mean: 0.0, min: 0.0, max: 0.0 });
    }
}
