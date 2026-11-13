//! k-mer counting over one or many FASTA records.
//!
//! Counting uses a rolling 2-bit hash (`HashMap<u64, u64>`). With canonical
//! mode on (the default, as in Jellyfish's `-C` flag), each k-mer is stored
//! under `min(kmer, reverse_complement(kmer))`, so a k-mer and its reverse
//! complement are counted together — the correct behaviour for
//! double-stranded DNA, where the sequenced strand is unknown.

use crate::errors::KmerError;
use crate::fasta::Record;
use crate::kmer::{canonical_code, kmer_mask, sliding_kmers, validate_k};
use std::collections::HashMap;

/// Per-record counting summary (batch mode).
#[derive(Debug, Clone, PartialEq)]
pub struct RecordSummary {
    /// Record id.
    pub id: String,
    /// Sequence length in bases.
    pub length: usize,
    /// Number of k-mer windows counted (excludes skipped windows).
    pub total_kmers: u64,
    /// Distinct k-mers contributed by this record alone.
    pub distinct_kmers: usize,
    /// Windows skipped because they contain an ambiguous base.
    pub skipped_windows: u64,
}

/// The k-mer counter: configuration + accumulated counts.
#[derive(Debug, Clone)]
pub struct KmerCounter {
    /// k-mer size.
    pub k: usize,
    /// Whether k-mers are stored in canonical form.
    pub canonical: bool,
    /// Packed k-mer code -> occurrence count.
    pub counts: HashMap<u64, u64>,
    /// Total windows counted (sum of all counts).
    pub total: u64,
    /// Windows skipped due to ambiguous bases.
    pub skipped_windows: u64,
    /// Total bases processed.
    pub total_bases: u64,
    /// Number of records processed.
    pub records: usize,
    /// Per-record summaries (batch mode).
    pub per_record: Vec<RecordSummary>,
}

impl KmerCounter {
    /// Create a counter for k-mer size `k` (1..=32).
    pub fn new(k: usize, canonical: bool) -> Result<Self, KmerError> {
        validate_k(k)?;
        Ok(Self {
            k,
            canonical,
            counts: HashMap::new(),
            total: 0,
            skipped_windows: 0,
            total_bases: 0,
            records: 0,
            per_record: Vec::new(),
        })
    }

    /// Count one sequence (raw bytes, any case). Ambiguous bases break the
    /// window; windows spanning them are skipped and reported.
    pub fn count_sequence(&mut self, seq: &[u8]) {
        let (kmers, skipped) = sliding_kmers(seq, self.k);
        self.total_bases += seq.len() as u64;
        self.skipped_windows += skipped as u64;
        for code in kmers {
            let key = if self.canonical {
                canonical_code(code, self.k)
            } else {
                code
            };
            *self.counts.entry(key).or_insert(0) += 1;
            self.total += 1;
        }
    }

    /// Count one FASTA record and store its per-record summary.
    pub fn count_record(&mut self, record: &Record) {
        let before_total = self.total;
        let before_skipped = self.skipped_windows;
        self.count_sequence(record.seq.as_bytes());
        // Distinct k-mers contributed by this record alone: recount it in
        // isolation (cheap for the record sizes this tool targets).
        let mut alone = KmerCounter::new(self.k, self.canonical)
            .expect("k already validated");
        alone.count_sequence(record.seq.as_bytes());
        self.per_record.push(RecordSummary {
            id: record.id.clone(),
            length: record.len(),
            total_kmers: self.total - before_total,
            distinct_kmers: alone.counts.len(),
            skipped_windows: self.skipped_windows - before_skipped,
        });
        self.records += 1;
    }

    /// Count many records (batch mode).
    pub fn count_records(&mut self, records: &[Record]) {
        for rec in records {
            self.count_record(rec);
        }
    }

    /// Number of distinct k-mers observed.
    pub fn distinct(&self) -> usize {
        self.counts.len()
    }

    /// Count of one k-mer given as a string (0 if absent or ambiguous).
    pub fn get_str(&self, kmer: &str) -> u64 {
        use crate::kmer::encode_kmer;
        let Some(code) = encode_kmer(kmer.as_bytes(), self.k) else {
            return 0;
        };
        let key = if self.canonical {
            canonical_code(code, self.k)
        } else {
            code
        };
        self.counts.get(&key).copied().unwrap_or(0)
    }

    /// All (code, count) pairs sorted by count descending, then code
    /// ascending (deterministic order for reports).
    pub fn sorted_counts(&self) -> Vec<(u64, u64)> {
        let mut v: Vec<(u64, u64)> = self.counts.iter().map(|(&c, &n)| (c, n)).collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        v
    }

    /// Merge another counter with the same k and canonical setting.
    pub fn merge(&mut self, other: &KmerCounter) -> Result<(), KmerError> {
        if other.k != self.k || other.canonical != self.canonical {
            return Err(KmerError::Usage(format!(
                "cannot merge counters with different settings (k={} canonical={} vs k={} canonical={})",
                self.k, self.canonical, other.k, other.canonical
            )));
        }
        for (&code, &n) in &other.counts {
            *self.counts.entry(code).or_insert(0) += n;
        }
        self.total += other.total;
        self.skipped_windows += other.skipped_windows;
        self.total_bases += other.total_bases;
        self.records += other.records;
        self.per_record.extend(other.per_record.iter().cloned());
        Ok(())
    }

    /// The 2-bit mask for this k (used by callers that roll their own scan).
    pub fn mask(&self) -> u64 {
        kmer_mask(self.k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::parse_fasta;

    fn counter_of(seq: &str, k: usize, canonical: bool) -> KmerCounter {
        let mut c = KmerCounter::new(k, canonical).unwrap();
        c.count_sequence(seq.as_bytes());
        c
    }

    #[test]
    fn toy_sequence_non_canonical_hand_counted() {
        // "ACGTACGT", k=3: ACG x2, CGT x2, GTA x1, TAC x1 (6 total, 4 distinct).
        let c = counter_of("ACGTACGT", 3, false);
        assert_eq!(c.total, 6);
        assert_eq!(c.distinct(), 4);
        assert_eq!(c.get_str("ACG"), 2);
        assert_eq!(c.get_str("CGT"), 2);
        assert_eq!(c.get_str("GTA"), 1);
        assert_eq!(c.get_str("TAC"), 1);
        assert_eq!(c.get_str("AAA"), 0);
        assert_eq!(c.skipped_windows, 0);
    }

    #[test]
    fn toy_sequence_canonical_hand_counted() {
        // Canonical k=3: ACG absorbs CGT (rc), GTA absorbs TAC (rc):
        // ACG x4, GTA x2 (6 total, 2 distinct).
        let c = counter_of("ACGTACGT", 3, true);
        assert_eq!(c.total, 6);
        assert_eq!(c.distinct(), 2);
        assert_eq!(c.get_str("ACG"), 4);
        assert_eq!(c.get_str("CGT"), 4); // same canonical bucket
        assert_eq!(c.get_str("GTA"), 2);
        assert_eq!(c.get_str("TAC"), 2);
    }

    #[test]
    fn ambiguous_bases_skip_windows() {
        // "ACNTAC", k=3: only TAC survives; 3 windows skipped.
        let c = counter_of("ACNTAC", 3, false);
        assert_eq!(c.total, 1);
        assert_eq!(c.distinct(), 1);
        assert_eq!(c.skipped_windows, 3);
    }

    #[test]
    fn short_sequence_yields_no_kmers() {
        let c = counter_of("AC", 3, false);
        assert_eq!(c.total, 0);
        assert_eq!(c.distinct(), 0);
    }

    #[test]
    fn batch_mode_per_record_summaries() {
        let recs = parse_fasta(">r1\nACGTACGT\n>r2\nAAAA\n").unwrap();
        let mut c = KmerCounter::new(3, false).unwrap();
        c.count_records(&recs);
        assert_eq!(c.records, 2);
        assert_eq!(c.total, 6 + 2);
        assert_eq!(c.per_record[0].id, "r1");
        assert_eq!(c.per_record[0].total_kmers, 6);
        assert_eq!(c.per_record[0].distinct_kmers, 4);
        assert_eq!(c.per_record[1].id, "r2");
        assert_eq!(c.per_record[1].total_kmers, 2);
        assert_eq!(c.per_record[1].distinct_kmers, 1); // AAA only
    }

    #[test]
    fn merge_adds_counts() {
        let mut a = counter_of("ACGT", 3, false);
        let b = counter_of("ACGT", 3, false);
        a.merge(&b).unwrap();
        assert_eq!(a.total, 4);
        assert_eq!(a.get_str("ACG"), 2);
    }

    #[test]
    fn merge_rejects_mismatched_settings() {
        let mut a = counter_of("ACGT", 3, false);
        let b = counter_of("ACGT", 4, false);
        assert!(a.merge(&b).is_err());
    }

    #[test]
    fn periodic_repeat_hand_counted() {
        // "ACGT" x 30 = 120 bp, k=5: 116 windows, 4 distinct 5-mers x29 each.
        let seq = "ACGT".repeat(30);
        let c = counter_of(&seq, 5, false);
        assert_eq!(c.total, 116);
        assert_eq!(c.distinct(), 4);
        for kmer in ["ACGTA", "CGTAC", "GTACG", "TACGT"] {
            assert_eq!(c.get_str(kmer), 29, "k-mer {kmer}");
        }
        // Canonical: ACGTA absorbs TACGT, CGTAC absorbs GTACG -> 2 x58.
        let cc = counter_of(&seq, 5, true);
        assert_eq!(cc.distinct(), 2);
        assert_eq!(cc.get_str("ACGTA"), 58);
        assert_eq!(cc.get_str("CGTAC"), 58);
    }
}
