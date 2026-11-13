//! Minimal FASTA parsing: multi-record files, descriptions, blank lines.
//!
//! A record starts with a `>` header line; the record id is the first
//! whitespace-separated token after `>`, the rest of the line is the
//! description. Sequence lines are concatenated, upper-cased and stripped
//! of whitespace. Bases are kept as-is (including `N`); the counter decides
//! how ambiguous bases are handled.

use crate::errors::KmerError;
use std::fs;
use std::path::Path;

/// One FASTA record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// First token after `>`.
    pub id: String,
    /// Remainder of the header line (may be empty).
    pub description: String,
    /// Sequence, upper-cased, whitespace removed.
    pub seq: String,
}

impl Record {
    /// Length of the sequence in bases.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// True when the sequence is empty (never happens for records returned
    /// by [`parse_fasta`], which rejects them).
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// GC fraction of the whole record (G+C)/length; 0.0 for empty input.
    pub fn gc_content(&self) -> f64 {
        if self.seq.is_empty() {
            return 0.0;
        }
        let gc = self
            .seq
            .bytes()
            .filter(|&b| b == b'G' || b == b'C')
            .count();
        gc as f64 / self.seq.len() as f64
    }
}

/// Parse FASTA text into records.
///
/// Errors: empty input, sequence lines before any header, a record with an
/// empty sequence (explicit validation — a silent empty record would make
/// per-record statistics misleading).
pub fn parse_fasta(text: &str) -> Result<Vec<Record>, KmerError> {
    let mut records: Vec<Record> = Vec::new();
    let mut current: Option<Record> = None;

    for (lineno, raw_line) in text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(rec) = current.take() {
                if rec.seq.is_empty() {
                    return Err(KmerError::EmptySequence(rec.id));
                }
                records.push(rec);
            }
            let mut parts = header.splitn(2, char::is_whitespace);
            let id = parts.next().unwrap_or("").to_string();
            if id.is_empty() {
                return Err(KmerError::FastaParse(format!(
                    "line {}: header '>' with no record id",
                    lineno + 1
                )));
            }
            let description = parts.next().unwrap_or("").trim().to_string();
            current = Some(Record {
                id,
                description,
                seq: String::new(),
            });
        } else {
            match current.as_mut() {
                Some(rec) => {
                    // Strip internal whitespace, upper-case.
                    rec.seq.extend(line.chars().filter(|c| !c.is_whitespace()).map(|c| c.to_ascii_uppercase()));
                }
                None => {
                    return Err(KmerError::FastaParse(format!(
                        "line {}: sequence data before any '>' header",
                        lineno + 1
                    )));
                }
            }
        }
    }

    if let Some(rec) = current.take() {
        if rec.seq.is_empty() {
            return Err(KmerError::EmptySequence(rec.id));
        }
        records.push(rec);
    }
    if records.is_empty() {
        return Err(KmerError::FastaParse(
            "no records found (empty file or no '>' header)".to_string(),
        ));
    }
    Ok(records)
}

/// Read and parse a FASTA file.
pub fn read_fasta_file(path: &Path) -> Result<Vec<Record>, KmerError> {
    let text = fs::read_to_string(path).map_err(|e| {
        KmerError::Io(std::io::Error::new(
            e.kind(),
            format!("cannot read '{}': {e}", path.display()),
        ))
    })?;
    parse_fasta(&text)
}

/// Read and parse several FASTA files (batch mode), concatenating records.
pub fn read_fasta_files(paths: &[&Path]) -> Result<Vec<Record>, KmerError> {
    let mut all = Vec::new();
    for path in paths {
        all.extend(read_fasta_file(path)?);
    }
    Ok(all)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = ">seq1 first record\nACGTACGT\n>seq2\nACNT\nAC\n";

    #[test]
    fn parses_multi_record() {
        let recs = parse_fasta(SAMPLE).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "seq1");
        assert_eq!(recs[0].description, "first record");
        assert_eq!(recs[0].seq, "ACGTACGT");
        assert_eq!(recs[1].id, "seq2");
        assert_eq!(recs[1].seq, "ACNTAC");
    }

    #[test]
    fn lower_case_and_spaces_normalised() {
        let recs = parse_fasta(">x\nac gt\nac\n").unwrap();
        assert_eq!(recs[0].seq, "ACGTAC");
    }

    #[test]
    fn rejects_sequence_before_header() {
        assert!(parse_fasta("ACGT\n>x\nACGT\n").is_err());
    }

    #[test]
    fn rejects_empty_file_and_empty_record() {
        assert!(parse_fasta("").is_err());
        assert!(parse_fasta(">only\n").is_err());
        assert!(parse_fasta(">a\nACGT\n>b\n").is_err());
    }

    #[test]
    fn rejects_bare_header() {
        assert!(parse_fasta(">\nACGT\n").is_err());
    }

    #[test]
    fn gc_content_known() {
        let rec = Record {
            id: "t".into(),
            description: String::new(),
            seq: "ACGT".into(),
        };
        assert!((rec.gc_content() - 0.5).abs() < 1e-12);
    }
}
