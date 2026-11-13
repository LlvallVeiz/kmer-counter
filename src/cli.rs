//! Command-line interface: subcommands, flags, and dispatch.
//!
//! The CLI is parsed by hand (std-only, no clap) and has five subcommands:
//!
//! ```text
//! kmer-counter count     <files...>  [-k N] [--canonical|--no-canonical]
//!                                  [--min-count N] [--format table|tsv|csv|json]
//!                                  [-o FILE] [--top N]
//! kmer-counter spectrum  <files...>  [-k N] [--canonical|--no-canonical]
//!                                  [--format table|json]
//! kmer-counter stats     <files...>  [-k N] [--canonical|--no-canonical]
//!                                  [--format table|json] [--top N]
//! kmer-counter compare   <fileA> <fileB> [-k N] [--canonical|--no-canonical]
//! kmer-counter top       <files...>  [-k N] [--canonical|--no-canonical] [--top N]
//! ```

use crate::compare;
use crate::counter::KmerCounter;
use crate::errors::KmerError;
use crate::fasta::read_fasta_files;
use crate::filter;
use crate::kmer::validate_k;
use crate::output::{self, Format};
use crate::spectrum::Spectrum;
use crate::stats;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// k-mer size default: 21, the value recommended by GenomeScope
/// (Vurture et al., 2017) for short-read data — long enough that most
/// k-mers are unique, short enough to tolerate sequencing errors.
pub const DEFAULT_K: usize = 21;
pub const DEFAULT_TOP: usize = 10;
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub const USAGE: &str = "kmer-counter 1.0.0 — k-mer counting, spectra and comparison for DNA FASTA files

USAGE:
    kmer-counter <subcommand> [options] <input.fasta> [...]

SUBCOMMANDS:
    count      Count k-mers and print the count table
    spectrum   Print the k-mer frequency spectrum (histogram)
    stats      Print summary statistics (entropy, GC, top k-mers, ...)
    compare    Compare two samples (Jaccard, Bray-Curtis, cosine)
    top        Print only the N most frequent k-mers

OPTIONS:
    -k, --kmer-size N        k-mer size, 1..=32 (default: 21)
        --canonical          Count k-mers canonically (default)
        --no-canonical       Count forward-strand k-mers only
        --min-count N        Drop k-mers with count < N (count/stats only)
        --format FMT         table (default), tsv, csv, json
    -o, --output FILE        Write output to FILE instead of stdout
        --top N              Number of top k-mers to show (default: 10)
    -h, --help               Show this help
    -V, --version            Show version

EXAMPLES:
    kmer-counter count -k 3 examples/toy.fasta
    kmer-counter spectrum -k 5 examples/synthetic.fasta
    kmer-counter stats -k 5 --format json examples/synthetic.fasta
    kmer-counter compare -k 3 examples/compare_a.fasta examples/compare_b.fasta
";

/// Options shared by most subcommands.
#[derive(Debug, Clone)]
pub struct Options {
    pub k: usize,
    pub canonical: bool,
    pub min_count: u64,
    pub format: Format,
    pub output: Option<PathBuf>,
    pub top_n: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            k: DEFAULT_K,
            canonical: true,
            min_count: 1,
            format: Format::Table,
            output: None,
            top_n: DEFAULT_TOP,
        }
    }
}

/// A parsed subcommand invocation.
#[derive(Debug)]
pub enum Command {
    Count { files: Vec<PathBuf>, opts: Options },
    Spectrum { files: Vec<PathBuf>, opts: Options },
    Stats { files: Vec<PathBuf>, opts: Options },
    Compare { file_a: PathBuf, file_b: PathBuf, opts: Options },
    Top { files: Vec<PathBuf>, opts: Options },
    Help,
    Version,
}

/// Parse `std::env::args()` (skipping argv[0]) into a [`Command`].
pub fn parse_args(args: &[String]) -> Result<Command, KmerError> {
    if args.is_empty() {
        return Err(KmerError::Usage(USAGE.to_string()));
    }
    match args[0].as_str() {
        "-h" | "--help" | "help" => return Ok(Command::Help),
        "-V" | "--version" | "version" => return Ok(Command::Version),
        "count" | "spectrum" | "stats" | "compare" | "top" => {}
        other => {
            return Err(KmerError::Usage(format!(
                "unknown subcommand '{other}'.\n{USAGE}"
            )))
        }
    }
    let sub = args[0].as_str();
    let mut opts = Options::default();
    let mut files: Vec<PathBuf> = Vec::new();
    let mut i = 1;
    while i < args.len() {
        let arg = &args[i];
        match arg.as_str() {
            "-h" | "--help" => return Ok(Command::Help),
            "-k" | "--kmer-size" => {
                i += 1;
                let v = args.get(i).ok_or_else(|| {
                    KmerError::Usage("missing value for --kmer-size".to_string())
                })?;
                opts.k = v.parse().map_err(|_| {
                    KmerError::Usage(format!("invalid k-mer size '{v}' (expected an integer)"))
                })?;
                validate_k(opts.k)?;
            }
            "--canonical" => opts.canonical = true,
            "--no-canonical" => opts.canonical = false,
            "--min-count" => {
                i += 1;
                let v = args.get(i).ok_or_else(|| {
                    KmerError::Usage("missing value for --min-count".to_string())
                })?;
                opts.min_count = v.parse().map_err(|_| {
                    KmerError::Usage(format!("invalid --min-count '{v}' (expected an integer)"))
                })?;
            }
            "--format" => {
                i += 1;
                let v = args.get(i).ok_or_else(|| {
                    KmerError::Usage("missing value for --format".to_string())
                })?;
                opts.format = Format::parse(v).ok_or_else(|| {
                    KmerError::Usage(format!(
                        "unknown format '{v}' (expected table, tsv, csv, or json)"
                    ))
                })?;
            }
            "-o" | "--output" => {
                i += 1;
                let v = args.get(i).ok_or_else(|| {
                    KmerError::Usage("missing value for --output".to_string())
                })?;
                opts.output = Some(PathBuf::from(v));
            }
            "--top" => {
                i += 1;
                let v = args.get(i).ok_or_else(|| {
                    KmerError::Usage("missing value for --top".to_string())
                })?;
                opts.top_n = v.parse().map_err(|_| {
                    KmerError::Usage(format!("invalid --top '{v}' (expected an integer)"))
                })?;
            }
            _ if arg.starts_with('-') => {
                return Err(KmerError::Usage(format!("unknown option '{arg}'.\n{USAGE}")));
            }
            _ => files.push(PathBuf::from(arg)),
        }
        i += 1;
    }
    if sub == "compare" {
        if files.len() != 2 {
            return Err(KmerError::Usage(
                "compare needs exactly two input files.\n".to_string() + USAGE,
            ));
        }
        let mut it = files.into_iter();
        return Ok(Command::Compare {
            file_a: it.next().unwrap(),
            file_b: it.next().unwrap(),
            opts,
        });
    }
    if files.is_empty() {
        return Err(KmerError::Usage(format!(
            "'{sub}' needs at least one input FASTA file.\n{USAGE}"
        )));
    }
    Ok(match sub {
        "count" => Command::Count { files, opts },
        "spectrum" => Command::Spectrum { files, opts },
        "stats" => Command::Stats { files, opts },
        "top" => Command::Top { files, opts },
        _ => unreachable!(),
    })
}

/// Open the output destination: a file when `--output` is given, stdout
/// otherwise. Returns a boxed writer.
fn open_output(path: &Option<PathBuf>) -> Result<Box<dyn Write>, KmerError> {
    match path {
        Some(p) => Ok(Box::new(BufWriter::new(File::create(p).map_err(|e| {
            KmerError::Io(io::Error::new(
                e.kind(),
                format!("cannot create output file '{}': {e}", p.display()),
            ))
        })?))),
        None => Ok(Box::new(io::stdout().lock())),
    }
}

fn load_counter(files: &[PathBuf], opts: &Options) -> Result<KmerCounter, KmerError> {
    let paths: Vec<&Path> = files.iter().map(|p| p.as_path()).collect();
    let records = read_fasta_files(&paths)?;
    let mut counter = KmerCounter::new(opts.k, opts.canonical)?;
    counter.count_records(&records);
    Ok(counter)
}

/// Run a parsed command. Returns the process exit code (0 on success).
pub fn run(cmd: Command) -> Result<(), KmerError> {
    match cmd {
        Command::Help => {
            print!("{USAGE}");
            Ok(())
        }
        Command::Version => {
            println!("kmer-counter {VERSION}");
            Ok(())
        }
        Command::Count { files, opts } => {
            let counter = load_counter(&files, &opts)?;
            let counter = if opts.min_count > 1 {
                filter::filter_min(&counter, opts.min_count)
            } else {
                counter
            };
            let mut out = open_output(&opts.output)?;
            // `count` prints the full table; `--top` only applies to the
            // `top` and `stats` subcommands.
            output::write_counts(&mut out, &counter, opts.format, 0)
                .map_err(|e| KmerError::Output(e.to_string()))?;
            let _ = out.flush();
            Ok(())
        }
        Command::Spectrum { files, opts } => {
            let counter = load_counter(&files, &opts)?;
            let spectrum = Spectrum::from_counter(&counter);
            let mut out = open_output(&opts.output)?;
            output::write_spectrum(&mut out, &spectrum, opts.format)
                .map_err(|e| KmerError::Output(e.to_string()))?;
            let _ = out.flush();
            Ok(())
        }
        Command::Stats { files, opts } => {
            let counter = load_counter(&files, &opts)?;
            let counter = if opts.min_count > 1 {
                filter::filter_min(&counter, opts.min_count)
            } else {
                counter
            };
            let summary = stats::summarize(&counter, opts.top_n);
            let mut out = open_output(&opts.output)?;
            output::write_summary(&mut out, &summary, opts.format)
                .map_err(|e| KmerError::Output(e.to_string()))?;
            let _ = out.flush();
            Ok(())
        }
        Command::Top { files, opts } => {
            let counter = load_counter(&files, &opts)?;
            let summary = stats::summarize(&counter, opts.top_n);
            let mut out = open_output(&opts.output)?;
            output::write_top(&mut out, &summary.top, counter.k)
                .map_err(|e| KmerError::Output(e.to_string()))?;
            let _ = out.flush();
            Ok(())
        }
        Command::Compare { file_a, file_b, opts } => {
            let paths_a: Vec<&Path> = vec![file_a.as_path()];
            let paths_b: Vec<&Path> = vec![file_b.as_path()];
            let mut ca = KmerCounter::new(opts.k, opts.canonical)?;
            ca.count_records(&read_fasta_files(&paths_a)?);
            let mut cb = KmerCounter::new(opts.k, opts.canonical)?;
            cb.count_records(&read_fasta_files(&paths_b)?);
            let cmp = compare::compare(&ca, &cb)?;
            let mut out = open_output(&opts.output)?;
            output::write_comparison(&mut out, &cmp)
                .map_err(|e| KmerError::Output(e.to_string()))?;
            let _ = out.flush();
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_count_with_flags() {
        let cmd = parse_args(&args(&["count", "-k", "5", "--no-canonical", "a.fasta", "b.fasta"])).unwrap();
        match cmd {
            Command::Count { files, opts } => {
                assert_eq!(opts.k, 5);
                assert!(!opts.canonical);
                assert_eq!(files.len(), 2);
            }
            _ => panic!("wrong command"),
        }
    }

    #[test]
    fn defaults_are_k21_canonical_table() {
        let cmd = parse_args(&args(&["stats", "a.fasta"])).unwrap();
        match cmd {
            Command::Stats { opts, .. } => {
                assert_eq!(opts.k, DEFAULT_K);
                assert!(opts.canonical);
                assert_eq!(opts.format, Format::Table);
            }
            _ => panic!("wrong command"),
        }
    }

    #[test]
    fn rejects_bad_k_and_unknown_options() {
        assert!(parse_args(&args(&["count", "-k", "0", "a.fasta"])).is_err());
        assert!(parse_args(&args(&["count", "-k", "33", "a.fasta"])).is_err());
        assert!(parse_args(&args(&["count", "-k", "abc", "a.fasta"])).is_err());
        assert!(parse_args(&args(&["count", "--bogus", "a.fasta"])).is_err());
        assert!(parse_args(&args(&["count"])).is_err());
        assert!(parse_args(&args(&["frobnicate", "a.fasta"])).is_err());
    }

    #[test]
    fn compare_needs_two_files() {
        assert!(parse_args(&args(&["compare", "a.fasta"])).is_err());
        assert!(parse_args(&args(&["compare", "a.fasta", "b.fasta", "c.fasta"])).is_err());
        assert!(parse_args(&args(&["compare", "a.fasta", "b.fasta"])).is_ok());
    }

    #[test]
    fn rejects_bad_format() {
        assert!(parse_args(&args(&["count", "--format", "xml", "a.fasta"])).is_err());
    }
}
