//! Output writers: human-readable tables, TSV/CSV, and JSON.
//!
//! JSON is written by hand (no serde dependency — this crate is std-only)
//! and is safe because k-mer strings contain only ACGT characters.

use crate::compare::Comparison;
use crate::counter::KmerCounter;
use crate::kmer::decode_kmer;
use crate::spectrum::Spectrum;
use crate::stats::{Summary, TopKmer};
use std::io::Write;

/// Output format for count tables.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// Human-readable aligned table.
    Table,
    /// Tab-separated values: `kmer<TAB>count`.
    Tsv,
    /// Comma-separated values: `kmer,count`.
    Csv,
    /// JSON object/array.
    Json,
}

impl Format {
    /// Parse a format name from the CLI (`table`, `tsv`, `csv`, `json`).
    pub fn parse(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "table" => Some(Format::Table),
            "tsv" => Some(Format::Tsv),
            "csv" => Some(Format::Csv),
            "json" => Some(Format::Json),
            _ => None,
        }
    }
}

/// Write the k-mer count table in the requested format.
/// `limit` caps the number of rows (0 = no limit).
pub fn write_counts(
    w: &mut dyn Write,
    counter: &KmerCounter,
    format: Format,
    limit: usize,
) -> std::io::Result<()> {
    let rows = counter.sorted_counts();
    let rows: &[(u64, u64)] = if limit == 0 {
        &rows
    } else {
        &rows[..rows.len().min(limit)]
    };
    match format {
        Format::Table => {
            writeln!(w, "{:>width$}  {:>10}", "k-mer", "count", width = counter.k)?;
            writeln!(w, "{}", "-".repeat(counter.k + 12))?;
            for &(code, n) in rows {
                writeln!(w, "{:>width$}  {:>10}", decode_kmer(code, counter.k), n, width = counter.k)?;
            }
            writeln!(w, "({} distinct k-mers, {} total)", counter.distinct(), counter.total)?;
        }
        Format::Tsv => {
            for &(code, n) in rows {
                writeln!(w, "{}\t{}", decode_kmer(code, counter.k), n)?;
            }
        }
        Format::Csv => {
            writeln!(w, "kmer,count")?;
            for &(code, n) in rows {
                writeln!(w, "{},{}", decode_kmer(code, counter.k), n)?;
            }
        }
        Format::Json => {
            writeln!(w, "{{")?;
            writeln!(w, "  \"k\": {},", counter.k)?;
            writeln!(w, "  \"canonical\": {},", counter.canonical)?;
            writeln!(w, "  \"total\": {},", counter.total)?;
            writeln!(w, "  \"distinct\": {},", counter.distinct())?;
            writeln!(w, "  \"counts\": {{")?;
            let last = rows.len().saturating_sub(1);
            for (i, &(code, n)) in rows.iter().enumerate() {
                let comma = if i == last { "" } else { "," };
                writeln!(w, "    \"{}\": {}{}", decode_kmer(code, counter.k), n, comma)?;
            }
            writeln!(w, "  }}")?;
            writeln!(w, "}}")?;
        }
    }
    Ok(())
}

/// Write the frequency spectrum (histogram) in table or JSON form.
pub fn write_spectrum(
    w: &mut dyn Write,
    spectrum: &Spectrum,
    format: Format,
) -> std::io::Result<()> {
    match format {
        Format::Json => {
            writeln!(w, "{{")?;
            writeln!(w, "  \"total_distinct\": {},", spectrum.total_distinct())?;
            writeln!(w, "  \"total_kmers\": {},", spectrum.total_kmers())?;
            writeln!(w, "  \"singletons\": {},", spectrum.singletons())?;
            writeln!(w, "  \"mode_frequency\": {},", spectrum.mode_frequency())?;
            writeln!(w, "  \"histogram\": {{")?;
            let n = spectrum.histogram.len();
            for (i, (&freq, &count)) in spectrum.histogram.iter().enumerate() {
                let comma = if i + 1 == n { "" } else { "," };
                writeln!(w, "    \"{freq}\": {count}{comma}")?;
            }
            writeln!(w, "  }}")?;
            writeln!(w, "}}")?;
        }
        _ => {
            writeln!(w, "frequency  distinct_kmers  histogram")?;
            writeln!(w, "{}", spectrum.ascii_histogram(50))?;
            writeln!(
                w,
                "singletons: {}  mode frequency: {}  mean frequency: {:.2}",
                spectrum.singletons(),
                spectrum.mode_frequency(),
                spectrum.mean_frequency()
            )?;
            if let Some(g) = spectrum.estimate_genome_size() {
                writeln!(w, "naive genome-size estimate (total / mode): {g} bp")?;
            }
        }
    }
    Ok(())
}

/// Write the top-N k-mers as a small table.
pub fn write_top(w: &mut dyn Write, top: &[TopKmer], k: usize) -> std::io::Result<()> {
    writeln!(w, "{:>width$}  {:>10}  {:>8}", "k-mer", "count", "fraction", width = k)?;
    writeln!(w, "{}", "-".repeat(k + 22))?;
    for t in top {
        writeln!(w, "{:>width$}  {:>10}  {:>7.2}%", t.kmer, t.count, t.fraction * 100.0, width = k)?;
    }
    Ok(())
}

/// Write the full summary report (table or JSON).
pub fn write_summary(w: &mut dyn Write, summary: &Summary, format: Format) -> std::io::Result<()> {
    match format {
        Format::Json => {
            writeln!(w, "{{")?;
            writeln!(w, "  \"k\": {},", summary.k)?;
            writeln!(w, "  \"canonical\": {},", summary.canonical)?;
            writeln!(w, "  \"records\": {},", summary.records)?;
            writeln!(w, "  \"total_bases\": {},", summary.total_bases)?;
            writeln!(w, "  \"total_kmers\": {},", summary.total_kmers)?;
            writeln!(w, "  \"distinct_kmers\": {},", summary.distinct_kmers)?;
            writeln!(w, "  \"singletons\": {},", summary.singletons)?;
            writeln!(w, "  \"skipped_windows\": {},", summary.skipped_windows)?;
            writeln!(w, "  \"mean_frequency\": {:.6},", summary.mean_frequency)?;
            writeln!(w, "  \"median_frequency\": {:.6},", summary.median_frequency)?;
            writeln!(w, "  \"mode_frequency\": {},", summary.mode_frequency)?;
            writeln!(w, "  \"entropy_bits\": {:.6},", summary.entropy_bits)?;
            writeln!(w, "  \"gc_mean\": {:.6},", summary.gc.mean)?;
            writeln!(w, "  \"gc_min\": {:.6},", summary.gc.min)?;
            writeln!(w, "  \"gc_max\": {:.6},", summary.gc.max)?;
            match summary.genome_size_estimate {
                Some(g) => writeln!(w, "  \"genome_size_estimate\": {g},")?,
                None => writeln!(w, "  \"genome_size_estimate\": null,")?,
            }
            writeln!(w, "  \"top\": [")?;
            let n = summary.top.len();
            for (i, t) in summary.top.iter().enumerate() {
                let comma = if i + 1 == n { "" } else { "," };
                writeln!(
                    w,
                    "    {{\"kmer\": \"{}\", \"count\": {}, \"fraction\": {:.6}}}{}",
                    t.kmer, t.count, t.fraction, comma
                )?;
            }
            writeln!(w, "  ]")?;
            writeln!(w, "}}")?;
        }
        _ => {
            writeln!(w, "k-mer summary (k = {}, canonical = {})", summary.k, summary.canonical)?;
            writeln!(w, "records processed:     {}", summary.records)?;
            writeln!(w, "total bases:           {}", summary.total_bases)?;
            writeln!(w, "total k-mers:          {}", summary.total_kmers)?;
            writeln!(w, "distinct k-mers:       {}", summary.distinct_kmers)?;
            writeln!(w, "singletons (freq 1):   {}", summary.singletons)?;
            writeln!(w, "skipped windows (N):   {}", summary.skipped_windows)?;
            writeln!(w, "mean frequency:        {:.2}", summary.mean_frequency)?;
            writeln!(w, "median frequency:      {:.2}", summary.median_frequency)?;
            writeln!(w, "mode frequency:        {}", summary.mode_frequency)?;
            writeln!(w, "entropy:               {:.4} bits", summary.entropy_bits)?;
            writeln!(
                w,
                "GC per distinct k-mer: mean {:.4}  min {:.4}  max {:.4}",
                summary.gc.mean, summary.gc.min, summary.gc.max
            )?;
            match summary.genome_size_estimate {
                Some(g) => writeln!(w, "genome-size estimate:    {g} bp (total / mode)")?,
                None => writeln!(w, "genome-size estimate:    n/a")?,
            }
            if !summary.top.is_empty() {
                writeln!(w)?;
                writeln!(w, "most frequent k-mers:")?;
                write_top(w, &summary.top, summary.k)?;
            }
        }
    }
    Ok(())
}

/// Write a two-sample comparison report.
pub fn write_comparison(w: &mut dyn Write, cmp: &Comparison) -> std::io::Result<()> {
    writeln!(w, "sample comparison (distinct k-mers)")?;
    writeln!(w, "distinct in A:         {}", cmp.distinct_a)?;
    writeln!(w, "distinct in B:         {}", cmp.distinct_b)?;
    writeln!(w, "shared:                {}", cmp.shared)?;
    writeln!(w, "only in A:             {}", cmp.only_a)?;
    writeln!(w, "only in B:             {}", cmp.only_b)?;
    writeln!(w, "Jaccard similarity:    {:.4}", cmp.jaccard)?;
    writeln!(w, "containment A in B:    {:.4}", cmp.containment_a)?;
    writeln!(w, "containment B in A:    {:.4}", cmp.containment_b)?;
    writeln!(w, "Bray-Curtis distance:  {:.4}", cmp.bray_curtis)?;
    writeln!(w, "cosine similarity:     {:.4}", cmp.cosine)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counter::KmerCounter;

    fn toy() -> KmerCounter {
        let mut c = KmerCounter::new(3, false).unwrap();
        c.count_sequence(b"ACGTACGT");
        c
    }

    #[test]
    fn tsv_output_is_sorted_by_count() {
        let mut buf = Vec::new();
        write_counts(&mut buf, &toy(), Format::Tsv, 0).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0], "ACG\t2");
        assert_eq!(lines[1], "CGT\t2");
        assert!(lines.contains(&"GTA\t1"));
    }

    #[test]
    fn csv_has_header() {
        let mut buf = Vec::new();
        write_counts(&mut buf, &toy(), Format::Csv, 0).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.starts_with("kmer,count\n"));
    }

    #[test]
    fn json_output_is_well_formed_text() {
        let mut buf = Vec::new();
        write_counts(&mut buf, &toy(), Format::Json, 0).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("\"k\": 3"));
        assert!(text.contains("\"ACG\": 2"));
        assert!(text.contains("\"distinct\": 4"));
    }

    #[test]
    fn table_limit_applies() {
        let mut buf = Vec::new();
        write_counts(&mut buf, &toy(), Format::Table, 2).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("ACG"));
        assert!(!text.contains("GTA"));
    }

    #[test]
    fn summary_table_mentions_key_numbers() {
        let s = crate::stats::summarize(&toy(), 2);
        let mut buf = Vec::new();
        write_summary(&mut buf, &s, Format::Table).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("total k-mers:          6"));
        assert!(text.contains("distinct k-mers:       4"));
        assert!(text.contains("entropy:               1.9183 bits"));
    }
}
