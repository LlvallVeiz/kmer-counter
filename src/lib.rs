//! kmer-counter — k-mer counting, frequency spectra, statistics and
//! sample comparison for DNA FASTA files.
//!
//! A **k-mer** is a substring of length *k* of a DNA sequence. Counting how
//! often each k-mer occurs is a central step in genome assembly, error
//! correction, genome-size estimation and alignment-free comparison
//! (Marçais & Kingsford, 2011). This crate implements that step in pure
//! Rust (std-only, no external dependencies):
//!
//! - 2-bit packed k-mer encoding (A=00, C=01, G=10, T=11), k up to 32;
//! - canonical k-mers (a k-mer and its reverse complement counted together);
//! - multi-record FASTA parsing with batch mode;
//! - frequency spectra (histograms) and naive genome-size estimates;
//! - summary statistics: Shannon entropy, GC per k-mer, top k-mers;
//! - frequency filtering (`--min-count`) for error-k-mer removal;
//! - two-sample comparison: Jaccard, containment, Bray–Curtis, cosine.
//!
//! See `docs/RESEARCH.md` for the background and references, and the
//! crate README for worked examples with hand-verified numbers.

pub mod cli;
pub mod compare;
pub mod counter;
pub mod errors;
pub mod fasta;
pub mod filter;
pub mod kmer;
pub mod output;
pub mod spectrum;
pub mod stats;
