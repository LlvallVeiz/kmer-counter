#!/usr/bin/env bash
# Regenerate every example output in examples/ from the real binary.
# Run from the repository root after `cargo build --release`.
set -euo pipefail

B="${B:-target/release/kmer-counter}"
if [ ! -x "$B" ]; then
  echo "binary not found: $B (run 'cargo build --release' first)" >&2
  exit 1
fi

"$B" count -k 3 --no-canonical examples/toy.fasta > examples/toy.k3.output.txt
"$B" count -k 3 examples/toy.fasta > examples/toy.k3.canonical.output.txt
"$B" count -k 3 --no-canonical --format csv examples/toy.fasta > examples/toy.k3.output.csv
"$B" count -k 3 --no-canonical --format json examples/toy.fasta > examples/toy.k3.output.json
"$B" spectrum -k 3 --no-canonical examples/toy.fasta > examples/toy.k3.spectrum.txt
"$B" stats -k 3 --no-canonical examples/toy.fasta > examples/toy.k3.stats.txt
"$B" count -k 5 examples/synthetic.fasta > examples/synthetic.k5.output.txt
"$B" spectrum -k 5 examples/synthetic.fasta > examples/synthetic.k5.spectrum.txt
"$B" stats -k 5 examples/synthetic.fasta > examples/synthetic.k5.stats.txt
"$B" stats -k 5 --format json examples/synthetic.fasta > examples/synthetic.k5.stats.json
"$B" compare -k 3 examples/compare_a.fasta examples/compare_b.fasta > examples/compare.output.txt
"$B" top -k 5 --top 5 examples/synthetic.fasta > examples/synthetic.k5.top.txt
"$B" count -k 5 --min-count 17 examples/synthetic.fasta > examples/synthetic.k5.mincount17.output.txt

echo "All example outputs regenerated."
