//! Integration tests for the FASTA parser (public API).

use kmer_counter::fasta::parse_fasta;

#[test]
fn parses_realistic_multi_record_file() {
    let text = ">chr1 synthetic chromosome 1\nACGTACGTAC\nGTAC\n>chr2\nNNNNACGT\n";
    let recs = parse_fasta(text).unwrap();
    assert_eq!(recs.len(), 2);
    assert_eq!(recs[0].id, "chr1");
    assert_eq!(recs[0].description, "synthetic chromosome 1");
    assert_eq!(recs[0].seq, "ACGTACGTACGTAC");
    assert_eq!(recs[1].seq, "NNNNACGT");
}

#[test]
fn windows_line_endings_ok() {
    let text = ">r1\r\nACGT\r\nACGT\r\n";
    let recs = parse_fasta(text).unwrap();
    assert_eq!(recs[0].seq, "ACGTACGT");
}

#[test]
fn parse_errors_are_descriptive() {
    let err = parse_fasta("ACGT\n").unwrap_err();
    assert!(err.to_string().contains("before any"));

    let err = parse_fasta(">empty\n").unwrap_err();
    assert!(err.to_string().contains("empty sequence"));
}
