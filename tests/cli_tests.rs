//! End-to-end CLI tests: run the real binary on fixture FASTA files.

use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kmer-counter"))
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join(name)
}

fn run_cli(args: &[&str]) -> (String, String, i32) {
    let out = Command::new(bin())
        .args(args)
        .output()
        .expect("failed to run kmer-counter binary");
    (
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
        out.status.code().unwrap_or(-1),
    )
}

#[test]
fn cli_count_toy_matches_hand_count() {
    let toy = fixture("toy.fasta");
    let (stdout, _, code) = run_cli(&[
        "count",
        "-k",
        "3",
        "--no-canonical",
        toy.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "stderr: {stdout}");
    assert!(stdout.contains("ACG"));
    assert!(stdout.contains("2"));
    assert!(stdout.contains("4 distinct k-mers, 6 total"));
}

#[test]
fn cli_count_canonical_toy() {
    let toy = fixture("toy.fasta");
    let (stdout, _, code) = run_cli(&["count", "-k", "3", toy.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(stdout.contains("2 distinct k-mers, 6 total"));
}

#[test]
fn cli_count_tsv_format() {
    let toy = fixture("toy.fasta");
    let (stdout, _, code) = run_cli(&[
        "count",
        "-k",
        "3",
        "--no-canonical",
        "--format",
        "tsv",
        toy.to_str().unwrap(),
    ]);
    assert_eq!(code, 0);
    assert!(stdout.contains("ACG\t2\n"));
}

#[test]
fn cli_spectrum_toy() {
    let toy = fixture("toy.fasta");
    let (stdout, _, code) = run_cli(&[
        "spectrum",
        "-k",
        "3",
        "--no-canonical",
        toy.to_str().unwrap(),
    ]);
    assert_eq!(code, 0);
    assert!(stdout.contains("singletons: 2"));
    assert!(stdout.contains("mode frequency: 1"));
}

#[test]
fn cli_stats_toy_reports_entropy() {
    let toy = fixture("toy.fasta");
    let (stdout, _, code) = run_cli(&[
        "stats",
        "-k",
        "3",
        "--no-canonical",
        toy.to_str().unwrap(),
    ]);
    assert_eq!(code, 0);
    assert!(stdout.contains("total k-mers:          6"));
    assert!(stdout.contains("entropy:               1.9183 bits"));
    assert!(stdout.contains("most frequent k-mers:"));
}

#[test]
fn cli_compare_toy_pair() {
    let a = fixture("compare_a.fasta");
    let b = fixture("compare_b.fasta");
    let (stdout, _, code) = run_cli(&[
        "compare",
        "-k",
        "3",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0);
    assert!(stdout.contains("Jaccard similarity:    0.6667"));
    assert!(stdout.contains("shared:                2"));
    assert!(stdout.contains("only in B:             1"));
}

#[test]
fn cli_top_limits_rows() {
    let toy = fixture("toy.fasta");
    let (stdout, _, code) = run_cli(&[
        "top",
        "-k",
        "3",
        "--no-canonical",
        "--top",
        "2",
        toy.to_str().unwrap(),
    ]);
    assert_eq!(code, 0);
    assert!(stdout.contains("ACG"));
    assert!(!stdout.contains("GTA"));
}

#[test]
fn cli_errors_are_nonzero() {
    // Missing file.
    let (_, stderr, code) = run_cli(&["count", "-k", "3", "no-such-file.fasta"]);
    assert_ne!(code, 0);
    assert!(stderr.contains("error"));
    // Invalid k.
    let (_, _, code) = run_cli(&["count", "-k", "99", "x.fasta"]);
    assert_ne!(code, 0);
    // No subcommand.
    let (_, _, code) = run_cli(&[]);
    assert_ne!(code, 0);
}

#[test]
fn cli_version_and_help() {
    let (stdout, _, code) = run_cli(&["--version"]);
    assert_eq!(code, 0);
    assert!(stdout.contains("kmer-counter 1.0.0"));
    let (stdout, _, code) = run_cli(&["--help"]);
    assert_eq!(code, 0);
    assert!(stdout.contains("SUBCOMMANDS"));
}
