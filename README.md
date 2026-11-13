# kmer-counter

Fast k-mer counting, frequency spectra, statistics and sample comparison
for DNA FASTA files — in pure Rust (std-only, zero dependencies).

A **k-mer** is a substring of length *k* of a DNA sequence. Counting how
often each k-mer occurs is the first step of genome assembly, read error
correction, reference-free genome profiling and alignment-free sequence
comparison (Marçais & Kingsford, 2011). This tool implements that step
for small-to-medium FASTA files, with every number hand-verifiable on
the toy examples in `examples/`.

## Features

- **k-mer counting** from multi-record FASTA files (batch mode), k = 1..=32,
  2-bit packed encoding (A=00, C=01, G=10, T=11) with a rolling window.
- **Canonical k-mers** (default): a k-mer and its reverse complement are
  counted together, as in Jellyfish `-C` / KMC. `--no-canonical` counts
  forward-strand k-mers only.
- **Frequency spectrum** (histogram): how many distinct k-mers occur
  exactly *f* times, with a naive genome-size estimate (`total / modal
  frequency`).
- **Summary statistics**: total/distinct k-mers, singletons, mean/median/
  modal frequency, Shannon entropy of the count distribution, GC content
  per distinct k-mer, most frequent k-mers.
- **Frequency filtering** (`--min-count`): drop low-frequency k-mers —
  the standard error-k-mer removal step before assembly.
- **Sample comparison**: Jaccard similarity, containment, Bray–Curtis
  dissimilarity and cosine similarity between two samples.
- **Output formats**: human-readable table (default), TSV, CSV, JSON.
- **Validation**: ambiguous bases (`N`, …) break the k-mer window and
  are reported as skipped windows instead of being silently counted;
  empty sequences and `k > sequence length` are handled explicitly.

## Installation

Requires a Rust toolchain (cargo ≥ 1.70; tested with rustc 1.99.1).

```bash
cargo build --release
# binary: target/release/kmer-counter
```

No external crates are needed — the project builds fully offline.

## Usage

```text
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
```

## Worked example — hand-verified

Input `examples/toy.fasta` is the 8 bp sequence `ACGTACGT`.
For k = 3 it contains 6 windows: ACG, CGT, GTA, TAC, ACG, CGT, so by
hand: **ACG ×2, CGT ×2, GTA ×1, TAC ×1 — 6 total, 4 distinct**.
Canonical counting merges each k-mer with its reverse complement
(rc(CGT) = ACG, rc(TAC) = GTA), giving **ACG ×4, GTA ×2 — 2 distinct**.

```bash
$ kmer-counter count -k 3 --no-canonical examples/toy.fasta
k-mer       count
---------------
ACG           2
CGT           2
GTA           1
TAC           1
(4 distinct k-mers, 6 total)

$ kmer-counter count -k 3 examples/toy.fasta   # canonical (default)
k-mer       count
---------------
ACG           4
GTA           2
(2 distinct k-mers, 6 total)

$ kmer-counter stats -k 3 --no-canonical examples/toy.fasta
k-mer summary (k = 3, canonical = false)
records processed:     1
total bases:           8
total k-mers:          6
distinct k-mers:       4
singletons (freq 1):   2
skipped windows (N):   0
mean frequency:        1.50
median frequency:      1.50
mode frequency:        1
entropy:               1.9183 bits
GC per distinct k-mer: mean 0.5000  min 0.3333  max 0.6667
genome-size estimate:    6 bp (total / mode)

most frequent k-mers:
k-mer       count  fraction
-------------------------
ACG           2    33.33%
CGT           2    33.33%
GTA           1    16.67%
TAC           1    16.67%
```

The entropy is checkable by hand: counts 2, 2, 1, 1 over 6 give
H = −(2·⅓·log₂⅓ + 2·⅙·log₂⅙) = **1.9183 bits**; canonical counts 4, 2
give H = **0.9183 bits**. Both are asserted in the test suite.

## Larger example — synthetic multi-record file

`examples/synthetic.fasta` has 3 synthetic records (300 bp total):
`ACGT` ×30, `GATTACA` ×17 + `G`, and a 60 bp mixed sequence.

```bash
$ kmer-counter stats -k 5 examples/synthetic.fasta
k-mer summary (k = 5, canonical = true)
records processed:     3
total bases:           300
total k-mers:          288
distinct k-mers:       32
singletons (freq 1):   11
skipped windows (N):   0
mean frequency:        9.00
median frequency:      2.00
mode frequency:        1
entropy:               3.8467 bits
GC per distinct k-mer: mean 0.4562  min 0.2000  max 0.8000
genome-size estimate:    288 bp (total / mode)

most frequent k-mers:
k-mer       count  fraction
---------------------------
ACGTA          58    20.14%
CGTAC          58    20.14%
ATTAC          17     5.90%
CTGTA          17     5.90%
GATTA          17     5.90%
```

Hand-check: the `ACGT` repeat yields 116 windows over 4 distinct 5-mers
(29 each); canonical merging leaves 2 buckets of 58 — exactly the top
two rows above. The `GATTACA` repeat yields 116 windows over 7 distinct
5-mers (four ×17, three ×16), visible as the 17/16 rows. The spectrum
(`examples/synthetic.k5.spectrum.txt`) shows the same structure:
frequency 58 → 2 distinct k-mers, 17 → 4, 16 → 3, and 11 singletons
from the mixed record.

Filtering to solid k-mers keeps only the repeats:

```bash
$ kmer-counter count -k 5 --min-count 17 examples/synthetic.fasta
k-mer       count
-----------------
ACGTA          58
CGTAC          58
ATTAC          17
CTGTA          17
GATTA          17
TGTAA          17
(6 distinct k-mers, 184 total)
```

## Comparing two samples

`examples/compare_a.fasta` (`ACGTACGT`) vs `compare_b.fasta`
(`ACGTACGA`) at k = 3, canonical: A = {ACG, GTA}, B = {ACG, GTA, CGA}.

```bash
$ kmer-counter compare -k 3 examples/compare_a.fasta examples/compare_b.fasta
sample comparison (distinct k-mers)
distinct in A:         2
distinct in B:         3
shared:                2
only in A:             0
only in B:             1
Jaccard similarity:    0.6667
containment A in B:    1.0000
containment B in A:    0.6667
Bray-Curtis distance:  0.1667
cosine similarity:     0.9562
```

Jaccard = 2/3 by hand (|A ∩ B| / |A ∪ B|); Bray–Curtis = 2/12 = 0.1667
from counts {ACG:4, GTA:2} vs {ACG:3, GTA:2, CGA:1}.

## Project structure

```text
kmer-counter/
├── Cargo.toml            # std-only crate, no dependencies
├── README.md
├── LICENSE               # MIT
├── docs/
│   └── RESEARCH.md       # background + references (Jellyfish, KMC 3, DSK, GenomeScope, ...)
├── src/
│   ├── main.rs           # binary entry point
│   ├── lib.rs            # crate docs
│   ├── cli.rs            # subcommands, flag parsing, dispatch
│   ├── kmer.rs           # 2-bit encoding, reverse complement, canonical form
│   ├── fasta.rs          # multi-record FASTA parser
│   ├── counter.rs        # counting, batch mode, merge
│   ├── spectrum.rs       # frequency histogram + genome-size estimate
│   ├── stats.rs          # entropy, GC stats, top k-mers, summary
│   ├── filter.rs         # --min-count filtering, solid k-mers, top-N
│   ├── compare.rs        # Jaccard / containment / Bray-Curtis / cosine
│   ├── output.rs         # table / TSV / CSV / JSON writers
│   └── errors.rs         # error type
├── tests/                # integration tests (public API + CLI end-to-end)
│   ├── kmer_tests.rs  fasta_tests.rs  counter_tests.rs
│   ├── spectrum_tests.rs  stats_tests.rs  compare_tests.rs
│   └── cli_tests.rs
└── examples/
    ├── toy.fasta / compare_a.fasta / compare_b.fasta / synthetic.fasta
    ├── *.output.txt / *.output.csv / *.output.json   # real outputs (see run_examples.sh)
    └── run_examples.sh   # regenerates every example output
```

## Tests

```bash
cargo test        # 78 tests: 49 unit + 29 integration (CLI end-to-end included)
cargo clippy --all-targets   # 0 warnings
```

The suite includes hand-computed values asserted exactly: the toy
3-mer counts (canonical and not), the periodic-repeat 5-mer counts
(29/58 and 17/16), the toy entropy (1.9183 / 0.9183 bits), the
spectrum shape, and the comparison metrics (Jaccard 2/3, Bray–Curtis
2/12, cosine 16/√(20·14)).

## Design decisions

- **2-bit packing in `u64`** (k ≤ 32), the same encoding Jellyfish
  uses; rolling-window counting is O(1) per base.
- **Canonical by default**, matching KMC and Jellyfish `-C`.
- **Ambiguous bases break the window** and are reported
  (`skipped_windows`); nothing is silently invented.
- **Default k = 21**, the value GenomeScope recommends for short reads.
- **Std-only, zero dependencies**: builds offline, single static binary.
- **Deterministic output**: count tables sort by count desc, then k-mer
  asc, so diffs between runs are meaningful.

## Limitations

- FASTA only (no FASTQ/quality scores); counts are quality-blind.
- In-memory hash table: memory grows with distinct k-mers, so a
  mammalian genome at k = 21 is out of scope — use KMC 3 or Jellyfish 2.
- Single-threaded; no disk-based binning (DSK/KMC territory).
- The genome-size estimate is naive (`total / modal frequency`); for
  real data use GenomeScope 2.0. See `docs/RESEARCH.md` §6.

## References

See `docs/RESEARCH.md` for the full background and reference list
(Jellyfish, KMC 3, KMC 2, DSK, GenomeScope, Mash, …).
