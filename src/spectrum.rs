//! The k-mer frequency spectrum: how many distinct k-mers occur exactly
//! `f` times, for each frequency `f`.
//!
//! The spectrum is the classic diagnostic plot of genome projects: in a
//! shotgun dataset the true genomic k-mers form a peak around the sequencing
//! coverage, while sequencing errors pile up at frequency 1 (see
//! Vurture et al., GenomeScope, 2017, and the UC Davis genome-assembly
//! workshop notes). From the spectrum this module also derives a naive
//! genome-size estimate: total k-mers divided by the modal (peak)
//! frequency, following the UConn genome-size estimation tutorial.

use crate::counter::KmerCounter;
use std::collections::BTreeMap;

/// Frequency spectrum: frequency -> number of distinct k-mers with it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Spectrum {
    /// Sorted by frequency ascending.
    pub histogram: BTreeMap<u64, u64>,
}

impl Spectrum {
    /// Build the spectrum from a counter's counts.
    pub fn from_counter(counter: &KmerCounter) -> Self {
        let mut histogram: BTreeMap<u64, u64> = BTreeMap::new();
        for &count in counter.counts.values() {
            *histogram.entry(count).or_insert(0) += 1;
        }
        Self { histogram }
    }

    /// Total distinct k-mers (sum of the histogram).
    pub fn total_distinct(&self) -> u64 {
        self.histogram.values().sum()
    }

    /// Total k-mer occurrences (sum of frequency * count).
    pub fn total_kmers(&self) -> u64 {
        self.histogram.iter().map(|(&f, &n)| f * n).sum()
    }

    /// Number of singleton k-mers (frequency 1) — the classic error proxy.
    pub fn singletons(&self) -> u64 {
        self.histogram.get(&1).copied().unwrap_or(0)
    }

    /// Highest observed frequency.
    pub fn max_frequency(&self) -> u64 {
        self.histogram.keys().next_back().copied().unwrap_or(0)
    }

    /// Modal frequency: the frequency with the most distinct k-mers.
    /// Ties resolve to the smaller frequency.
    pub fn mode_frequency(&self) -> u64 {
        self.histogram
            .iter()
            .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
            .map(|(&f, _)| f)
            .unwrap_or(0)
    }

    /// Mean frequency (total k-mers / distinct k-mers).
    pub fn mean_frequency(&self) -> f64 {
        let d = self.total_distinct();
        if d == 0 {
            0.0
        } else {
            self.total_kmers() as f64 / d as f64
        }
    }

    /// Median frequency (over distinct k-mers).
    pub fn median_frequency(&self) -> f64 {
        let d = self.total_distinct();
        if d == 0 {
            return 0.0;
        }
        let mut freqs: Vec<u64> = Vec::with_capacity(d as usize);
        for (&f, &n) in &self.histogram {
            freqs.extend(std::iter::repeat_n(f, n as usize));
        }
        let mid = freqs.len() / 2;
        if freqs.len() % 2 == 1 {
            freqs[mid] as f64
        } else {
            (freqs[mid - 1] as f64 + freqs[mid] as f64) / 2.0
        }
    }

    /// Naive genome-size estimate: total k-mers / peak frequency.
    /// Returns `None` when the peak frequency is 0.
    pub fn estimate_genome_size(&self) -> Option<u64> {
        self.total_kmers().checked_div(self.mode_frequency())
    }

    /// Render an ASCII histogram: `frequency  distinct  bar`.
    /// Bars scale to at most `max_width` characters.
    pub fn ascii_histogram(&self, max_width: usize) -> String {
        let max_count = self.histogram.values().max().copied().unwrap_or(0);
        let mut lines = Vec::new();
        for (&freq, &n) in &self.histogram {
            let width = if max_count == 0 {
                0
            } else {
                ((n as usize * max_width) / max_count as usize).max(1)
            };
            lines.push(format!("{freq:>8} {n:>10} {}", "#".repeat(width)));
        }
        lines.join("\n")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counter::KmerCounter;

    fn toy_counter() -> KmerCounter {
        let mut c = KmerCounter::new(3, false).unwrap();
        c.count_sequence(b"ACGTACGT");
        c
    }

    #[test]
    fn toy_spectrum_hand_counted() {
        // Counts 2,2,1,1 -> spectrum: freq 1 -> 2 distinct, freq 2 -> 2.
        let s = Spectrum::from_counter(&toy_counter());
        assert_eq!(s.histogram.get(&1), Some(&2));
        assert_eq!(s.histogram.get(&2), Some(&2));
        assert_eq!(s.total_distinct(), 4);
        assert_eq!(s.total_kmers(), 6);
        assert_eq!(s.singletons(), 2);
        assert_eq!(s.max_frequency(), 2);
        assert!((s.mean_frequency() - 1.5).abs() < 1e-12);
        assert!((s.median_frequency() - 1.5).abs() < 1e-12);
    }

    #[test]
    fn canonical_toy_spectrum() {
        // Canonical counts 4,2 -> spectrum: freq 2 -> 1, freq 4 -> 1.
        let mut c = KmerCounter::new(3, true).unwrap();
        c.count_sequence(b"ACGTACGT");
        let s = Spectrum::from_counter(&c);
        assert_eq!(s.histogram.get(&2), Some(&1));
        assert_eq!(s.histogram.get(&4), Some(&1));
        assert_eq!(s.mode_frequency(), 2); // tie 1 vs 1 -> smaller frequency
        assert_eq!(s.estimate_genome_size(), Some(6 / 2));
    }

    #[test]
    fn empty_spectrum_is_zero() {
        let c = KmerCounter::new(3, false).unwrap();
        let s = Spectrum::from_counter(&c);
        assert_eq!(s.total_distinct(), 0);
        assert_eq!(s.mode_frequency(), 0);
        assert_eq!(s.estimate_genome_size(), None);
        assert_eq!(s.ascii_histogram(40), "");
    }

    #[test]
    fn ascii_histogram_scales() {
        let s = Spectrum::from_counter(&toy_counter());
        let h = s.ascii_histogram(10);
        assert!(h.contains("       1          2 ##########"));
        assert!(h.contains("       2          2 ##########"));
    }
}
