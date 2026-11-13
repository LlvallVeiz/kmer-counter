//! Comparing two samples by their k-mer sets and count profiles.
//!
//! Metrics (all computed on canonical k-mer sets by default, since the
//! counters compared are usually canonical):
//! - **Jaccard similarity** |A ∩ B| / |A ∪ B| over distinct k-mers.
//! - **Containment** |A ∩ B| / |A| (and symmetrically for B).
//! - **Bray–Curtis dissimilarity** Σ|aᵢ − bᵢ| / Σ(aᵢ + bᵢ) over the union.
//! - **Cosine similarity** of the two count vectors.
//!
//! These are the simple set/profile distances used e.g. by alignment-free
//! sequence comparison; Mash (Ondov et al., 2016) builds on the same
//! Jaccard idea with MinHash sketches for whole genomes.

use crate::counter::KmerCounter;

/// Result of comparing two k-mer counters.
#[derive(Debug, Clone, PartialEq)]
pub struct Comparison {
    /// Distinct k-mers in sample A.
    pub distinct_a: usize,
    /// Distinct k-mers in sample B.
    pub distinct_b: usize,
    /// Distinct k-mers present in both.
    pub shared: usize,
    /// Distinct k-mers only in A.
    pub only_a: usize,
    /// Distinct k-mers only in B.
    pub only_b: usize,
    /// Jaccard similarity in [0, 1].
    pub jaccard: f64,
    /// Containment of A in B: |A ∩ B| / |A|.
    pub containment_a: f64,
    /// Containment of B in A: |A ∩ B| / |B|.
    pub containment_b: f64,
    /// Bray–Curtis dissimilarity in [0, 1] (0 = identical profiles).
    pub bray_curtis: f64,
    /// Cosine similarity of the count vectors in [0, 1].
    pub cosine: f64,
}

/// Jaccard similarity of two distinct-k-mer sets.
pub fn jaccard(a: &KmerCounter, b: &KmerCounter) -> f64 {
    let shared = a.counts.keys().filter(|k| b.counts.contains_key(*k)).count();
    let union = a.distinct() + b.distinct() - shared;
    if union == 0 {
        1.0 // both empty: identical by convention
    } else {
        shared as f64 / union as f64
    }
}

/// Bray–Curtis dissimilarity of two count profiles.
pub fn bray_curtis(a: &KmerCounter, b: &KmerCounter) -> f64 {
    let mut num: u64 = 0;
    let mut den: u64 = 0;
    for (&code, &ca) in &a.counts {
        let cb = b.counts.get(&code).copied().unwrap_or(0);
        num += ca.abs_diff(cb);
        den += ca + cb;
    }
    for (&code, &cb) in &b.counts {
        if !a.counts.contains_key(&code) {
            num += cb;
            den += cb;
        }
    }
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}

/// Cosine similarity of two count vectors.
pub fn cosine(a: &KmerCounter, b: &KmerCounter) -> f64 {
    let mut dot: u64 = 0;
    let mut na: u64 = 0;
    let mut nb: u64 = 0;
    for (&code, &ca) in &a.counts {
        na += ca * ca;
        if let Some(&cb) = b.counts.get(&code) {
            dot += ca * cb;
        }
    }
    for &cb in b.counts.values() {
        nb += cb * cb;
    }
    let denom = (na as f64).sqrt() * (nb as f64).sqrt();
    if denom == 0.0 {
        1.0
    } else {
        dot as f64 / denom
    }
}

/// Compare two counters on all metrics. Requires equal k (canonical mode
/// may differ, but the result is only meaningful for equal settings).
pub fn compare(a: &KmerCounter, b: &KmerCounter) -> Result<Comparison, crate::errors::KmerError> {
    if a.k != b.k {
        return Err(crate::errors::KmerError::Usage(format!(
            "cannot compare counters with different k ({} vs {})",
            a.k, b.k
        )));
    }
    let shared = a.counts.keys().filter(|k| b.counts.contains_key(*k)).count();
    let only_a = a.distinct() - shared;
    let only_b = b.distinct() - shared;
    let union = shared + only_a + only_b;
    Ok(Comparison {
        distinct_a: a.distinct(),
        distinct_b: b.distinct(),
        shared,
        only_a,
        only_b,
        jaccard: if union == 0 { 1.0 } else { shared as f64 / union as f64 },
        containment_a: if a.distinct() == 0 {
            1.0
        } else {
            shared as f64 / a.distinct() as f64
        },
        containment_b: if b.distinct() == 0 {
            1.0
        } else {
            shared as f64 / b.distinct() as f64
        },
        bray_curtis: bray_curtis(a, b),
        cosine: cosine(a, b),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counter::KmerCounter;

    fn counter_of(seq: &str, k: usize, canonical: bool) -> KmerCounter {
        let mut c = KmerCounter::new(k, canonical).unwrap();
        c.count_sequence(seq.as_bytes());
        c
    }

    #[test]
    fn compare_identical_samples() {
        let a = counter_of("ACGTACGT", 3, true);
        let b = counter_of("ACGTACGT", 3, true);
        let cmp = compare(&a, &b).unwrap();
        assert_eq!(cmp.shared, 2);
        assert_eq!(cmp.only_a, 0);
        assert_eq!(cmp.only_b, 0);
        assert!((cmp.jaccard - 1.0).abs() < 1e-12);
        assert!((cmp.bray_curtis - 0.0).abs() < 1e-12);
        assert!((cmp.cosine - 1.0).abs() < 1e-12);
    }

    #[test]
    fn compare_toy_pair_hand_counted() {
        // A = "ACGTACGT" canonical set {ACG, GTA}; B = "ACGTACGA" -> {ACG, GTA, CGA}
        // (rc(CGA) = TCG, but CGA < TCG as a packed code, so canonical is CGA).
        // Jaccard = 2/3; containment A = 1.0; containment B = 2/3.
        let a = counter_of("ACGTACGT", 3, true);
        let b = counter_of("ACGTACGA", 3, true);
        assert_eq!(b.distinct(), 3);
        let cmp = compare(&a, &b).unwrap();
        assert_eq!(cmp.shared, 2);
        assert_eq!(cmp.only_a, 0);
        assert_eq!(cmp.only_b, 1);
        assert!((cmp.jaccard - 2.0 / 3.0).abs() < 1e-12);
        assert!((cmp.containment_a - 1.0).abs() < 1e-12);
        assert!((cmp.containment_b - 2.0 / 3.0).abs() < 1e-12);
        // Bray-Curtis: counts A {ACG:4, GTA:2}, B {ACG:3, GTA:2, CGA:1}:
        // |4-3|+|2-2|+|0-1| = 2; total = (4+3)+(2+2)+(0+1) = 12 -> 2/12.
        assert!((cmp.bray_curtis - 2.0 / 12.0).abs() < 1e-12);
        // Cosine: dot = 4*3 + 2*2 = 16; |A| = sqrt(20), |B| = sqrt(14).
        let expected = 16.0 / (20.0f64.sqrt() * 14.0f64.sqrt());
        assert!((cmp.cosine - expected).abs() < 1e-12);
    }

    #[test]
    fn compare_disjoint_samples() {
        let a = counter_of("AAAA", 3, false);
        let b = counter_of("CCCC", 3, false);
        let cmp = compare(&a, &b).unwrap();
        assert_eq!(cmp.shared, 0);
        assert!((cmp.jaccard - 0.0).abs() < 1e-12);
        assert!((cmp.bray_curtis - 1.0).abs() < 1e-12);
        assert!((cmp.cosine - 0.0).abs() < 1e-12);
    }

    #[test]
    fn compare_rejects_different_k() {
        let a = counter_of("ACGT", 3, false);
        let b = counter_of("ACGT", 4, false);
        assert!(compare(&a, &b).is_err());
    }
}
