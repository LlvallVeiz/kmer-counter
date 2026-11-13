# Research Notes — k-mer Counting

Background research for `kmer-counter`: what a k-mer is, why counting
them matters, how the established tools do it, and which design choices
this project took from them. Every claim below is traceable to the
references at the end; nothing here is invented.

## 1. What is a k-mer?

A **k-mer** is a substring of length *k* of a DNA sequence. A sequence of
length *L* contains exactly *L − k + 1* (overlapping) k-mers. For example,
`ACGTACGT` (L = 8) contains six 3-mers: ACG, CGT, GTA, TAC, ACG, CGT.

Because DNA is double-stranded and the sequenced strand is usually
unknown, most tools count **canonical k-mers**: a k-mer and its reverse
complement are treated as the same object, represented by the
lexicographically smaller of the two. Jellyfish exposes this as the `-C`
flag; KMC counts canonical k-mers by default. This project follows the
same convention (canonical on by default, `--no-canonical` to disable).

## 2. Why count k-mers?

Marçais & Kingsford (2011) list the central applications, and they have
not changed since:

- **Genome assembly.** De Bruijn-graph assemblers are built on k-mers:
  reads are decomposed into k-mers, and the assembly graph's nodes/edges
  are k-mers and their overlaps. Counting is the first step of the
  pipeline.
- **Error correction of sequencing reads.** In a shotgun dataset, true
  genomic k-mers appear many times (around the sequencing coverage),
  while k-mers created by sequencing errors appear once or twice. The
  k-mer frequency spectrum therefore separates signal from error: the
  low-frequency peak is (mostly) error. Correctors such as Quake and
  SOAPec exploit exactly this.
- **Genome profiling without a reference.** The k-mer spectrum of raw
  reads is enough to estimate genome size, heterozygosity and repeat
  content — the approach of GenomeScope (Vurture et al., 2017), which
  fits a mixture of negative-binomial distributions to the spectrum.
  GenomeScope recommends *k = 21* for short reads: long enough that most
  k-mers are unique, short enough to stay robust to sequencing errors.
  This project therefore defaults to *k = 21*.
- **Alignment-free comparison.** Two samples can be compared by their
  k-mer sets alone (Jaccard similarity, containment) or by their count
  profiles (Bray–Curtis, cosine). Mash (Ondov et al., 2016) scales the
  Jaccard idea to whole genomes with MinHash sketches.

## 3. The k-mer frequency spectrum

Plotting "how many distinct k-mers occur exactly *f* times" against *f*
gives the spectrum. In a real dataset it has a characteristic shape
(UC Davis Genome Assembly Workshop, 2020; UConn genome-size tutorial):

- a sharp peak at frequency 1 — sequencing errors;
- a main peak around the mean coverage — true single-copy genomic k-mers;
- smaller peaks at 2×, 3×, … coverage — repetitive sequence;
- for a heterozygous diploid, an extra peak at half coverage.

Two numbers from the spectrum are used directly in this project:

- **Singletons** (frequency 1) — a proxy for the error load.
- **Naive genome-size estimate**: total k-mers divided by the modal
  (peak) frequency, `G ≈ total / coverage`, as in the UConn tutorial.
  This is a rough estimate only (it assumes a homozygous, repeat-free
  genome); GenomeScope's model fitting is the rigorous version.

## 4. How the established tools count

| Tool | Approach | Key facts (from the cited papers) |
|------|----------|-----------------------------------|
| **Jellyfish** (Marçais & Kingsford, 2011) | Multithreaded, lock-free hash table; 2-bit packed k-mers up to k = 31 | C++, `jellyfish count -m 21 -C`; binary database + `histo`/`dump` subcommands |
| **KMC 3** (Kokot et al., 2017) | Disk-based, two-stage: reads split into bins by k-mer "signatures", each bin sorted (MSD radix sort) | Counts the 28-mers of a 106 GB human read set in ~20 min on a 6-core PC with ≤ 12 GB RAM |
| **DSK** (Rizk et al., 2013) | Disk streaming: fixed, user-defined memory + disk; k-mers partitioned by hash across passes | Counted all 27-mers of a human dataset with 4 GB RAM (17.9 h) |
| **KAnalyze** (Audano & Vannberg, 2014) | Split-and-merge: partition k-mer space, sort in memory, merge | Pipelined toolkit for k-mer counting and filtering |
| **GenomeScope** (Vurture et al., 2017) | Not a counter — consumes a k-mer histogram (e.g. from Jellyfish) and fits a statistical model | Genome-size estimates 99.7% accurate on 16 real datasets; recommends k = 21 |

A benchmark of these approaches (Manekar & Sathe, *GigaScience* 2018)
compares them on runtime, memory, disk usage and accuracy, and confirms
the basic trade-off: in-memory hash tables are fastest but memory-hungry;
disk-based sorters trade I/O for a small memory footprint.

## 5. Design choices in this project

This is a teaching-scale counter, not a competitor to Jellyfish/KMC: it
holds the whole count table in memory (`HashMap<u64, u64>`), which is
fine for genomes and read sets up to a few megabases, and it prioritises
transparency over speed.

- **2-bit packing** (A=00, C=01, G=10, T=11) in a `u64`, exactly the
  encoding Jellyfish uses — k up to 32, rolling-window counting in O(1)
  per base.
- **Canonical by default**, matching KMC/Jellyfish `-C`.
- **Ambiguous bases (N, …) break the window.** Any k-mer that would span
  an `N` is skipped and counted in `skipped_windows`, instead of
  silently inventing a base. This mirrors how production pipelines
  treat Ns (they cannot be trusted), and it is reported, not hidden.
- **`--min-count` filtering** is the same cutoff Jellyfish (`-L`) and
  KMC (`-ci`) offer to drop error k-mers before assembly.
- **Spectrum + naive genome-size estimate** follow the UConn tutorial's
  `total / peak` formula, labelled "naive" because the rigorous version
  needs GenomeScope's model fit.
- **Std-only Rust, zero dependencies.** The whole tool builds offline
  with `cargo build --release`; JSON output is written by hand.

## 6. What this tool does NOT do (limitations)

- No FASTQ/quality-score support (counts are quality-blind, like
  Jellyfish 1.x; KMC 3 reads gzipped FASTQ).
- No multithreading, no disk-based binning — memory grows with the
  number of *distinct* k-mers (≈ 4^k worst case), so a mammalian genome
  at k = 21 does not fit comfortably; use KMC 3 or Jellyfish 2 there.
- No Bloom-filter prefiltering (khmer/BFCounter territory).
- The genome-size estimate is naive (see §3); for real data use
  GenomeScope 2.0.
- k is capped at 32 by the `u64` packing (Jellyfish 1.x capped at 31
  for the same reason).

## References

1. Marçais, G. & Kingsford, C. (2011). A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. *Bioinformatics* 27(6):764–770. doi:10.1093/bioinformatics/btr011 — https://pubmed.ncbi.nlm.nih.gov/21217122/
2. Kokot, M., Długosz, M. & Deorowicz, S. (2017). KMC 3: counting and manipulating k-mer statistics. *Bioinformatics* 33(17):2759–2761. doi:10.1093/bioinformatics/btx304 — https://pubmed.ncbi.nlm.nih.gov/28472236/
3. Deorowicz, S., Kokot, M., Grabowski, S. & Debudaj-Grabysz, A. (2015). KMC 2: fast and resource-frugal k-mer counting. *Bioinformatics* 31(10):1569–1576. doi:10.1093/bioinformatics/btv022
4. Rizk, G., Lavenier, D. & Chikhi, R. (2013). DSK: k-mer counting with very low memory usage. *Bioinformatics* 29(5):652–653. doi:10.1093/bioinformatics/btt020
5. Vurture, G. W. et al. (2017). GenomeScope: fast reference-free genome profiling from short reads. *Bioinformatics* 33(14):2202–2204. doi:10.1093/bioinformatics/btx153 — https://pmc.ncbi.nlm.nih.gov/articles/PMC5870704/
6. Ondov, B. D. et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology* 17:132. doi:10.1186/s13059-016-0997-x
7. Manekar, S. C. & Sathe, S. R. (2018). A benchmark study of k-mer counting methods for high-throughput sequencing. *GigaScience* 7(12):giy125. doi:10.1093/gigascience/giy125 — https://pmc.ncbi.nlm.nih.gov/articles/PMC6280066/
8. Audano, P. & Vannberg, F. (2014). KAnalyze: a fast versatile pipelined Kmer toolkit. *Bioinformatics* 30(14):2070–2072.
9. UC Davis Bioinformatics Core (2020). Genome Assembly Workshop — K-mer frequency plots. https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/kmers/kmers
10. University of Connecticut, Computational Biology Core. Genome Size Estimation Tutorial. https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
11. Chikhi, R. (2018). De novo assembly & kmers (workshop slides). http://rayan.chikhi.name/pdf/2018-evomics-assembly-kmers.pdf
12. Ranallo-Benavidez, A. R. et al. (2020). GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature Communications* 11:1432.

*All references were checked on 2026-11-12 (PubMed/PMC/publisher pages
or the cited lecture slides).*
