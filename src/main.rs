//! Binary entry point: parse argv, dispatch, report errors on stderr.

use kmer_counter::cli::{parse_args, run, USAGE};
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(&args) {
        Ok(cmd) => match run(cmd) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("kmer-counter: error: {e}");
                ExitCode::FAILURE
            }
        },
        Err(e) => {
            eprintln!("{e}");
            if !e.to_string().contains(USAGE.trim()) {
                eprintln!("{USAGE}");
            }
            ExitCode::from(2)
        }
    }
}
