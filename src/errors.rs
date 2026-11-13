//! Error type for the whole crate.

use std::fmt;
use std::io;

/// All recoverable errors produced by kmer-counter.
#[derive(Debug)]
pub enum KmerError {
    /// Invalid k-mer size (0 or > 32).
    InvalidK(String),
    /// A FASTA record has no sequence lines.
    EmptySequence(String),
    /// FASTA structure is broken (sequence before a header, empty file, ...).
    FastaParse(String),
    /// A file could not be read.
    Io(io::Error),
    /// Command-line usage error.
    Usage(String),
    /// Output could not be written.
    Output(String),
}

impl fmt::Display for KmerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            KmerError::InvalidK(msg) => write!(f, "invalid k-mer size: {msg}"),
            KmerError::EmptySequence(id) => {
                write!(f, "FASTA record '{id}' has an empty sequence")
            }
            KmerError::FastaParse(msg) => write!(f, "FASTA parse error: {msg}"),
            KmerError::Io(err) => write!(f, "I/O error: {err}"),
            KmerError::Usage(msg) => write!(f, "{msg}"),
            KmerError::Output(msg) => write!(f, "output error: {msg}"),
        }
    }
}

impl std::error::Error for KmerError {}

impl From<io::Error> for KmerError {
    fn from(err: io::Error) -> Self {
        KmerError::Io(err)
    }
}
