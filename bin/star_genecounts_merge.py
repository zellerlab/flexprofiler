#!/usr/bin/env python3
"""Merge the STAR reads per gene (--quantMode GeneCounts) of the samples into tables of genes by samples.

usage: star_genecounts_merge.py <prefix> <samples.tsv>

<samples.tsv> has one line per sample, without a header: the sample name and the path of its ReadsPerGene.out.tab,
separated by a tab.

A ReadsPerGene.out.tab has four summary lines (N_unmapped, N_multimapping, N_noFeature, N_ambiguous), then one line per
gene of the STAR index: the gene ID and its counts for three strandednesses of the library, unstranded, forward (the
first read on the strand of the gene, htseq-count -s yes) and reverse (the second read on the strand of the gene,
htseq-count -s reverse, e.g. Illumina TruSeq stranded).

Writes:
- <prefix>.{unstranded,forward,reverse}.tsv: the counts of each gene (rows, in the order of the index) in each sample
  (columns, sorted)
- <prefix>.summary.tsv: for each strandedness, the four summary counts and the reads assigned to a gene (N_assigned),
  in each sample

All the samples must have the same genes, in the same order (the same STAR index).
"""
import csv
import sys

import pandas as pd

STRANDEDNESS = ["unstranded", "forward", "reverse"]
SUMMARY = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]


def read_tsv(path, columns):
    """A headerless TSV of strings, checking that it has the given number of columns."""
    table = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
    assert table.shape[1] == columns, f"{path}: {table.shape[1]} columns, expected {columns}"
    return table


def read_counts(path):
    """The summary counts and the gene counts of a ReadsPerGene.out.tab, each a table indexed by ID, by strandedness.

    The summary has the four summary lines and the reads assigned to a gene (N_assigned), the sum of the gene counts.
    """
    table = read_tsv(path, 1 + len(STRANDEDNESS)).set_index(0).astype("int64")
    table.columns = STRANDEDNESS
    summary = table.iloc[: len(SUMMARY)].copy()
    genes = table.iloc[len(SUMMARY) :]
    summary_labels = list(summary.index)
    assert summary_labels == SUMMARY, f"{path}: summary lines {summary_labels}, expected {SUMMARY}"
    summary.loc["N_assigned"] = genes.sum()
    return summary, genes


def main(prefix, sheet):
    pairs = read_tsv(sheet, 2).set_axis(["sample", "path"], axis=1)
    assert not pairs.empty, f"{sheet}: no samples"
    duplicates = sorted(pairs.loc[pairs["sample"].duplicated(), "sample"].unique())
    assert not duplicates, f"duplicate sample names: {', '.join(duplicates)}"
    pairs = pairs.sort_values("sample")

    summaries, genes = {}, {}
    for sample, path in zip(pairs["sample"], pairs["path"]):
        summaries[sample], genes[sample] = read_counts(path)
    first = pairs["sample"].iloc[0]
    for sample, counts in genes.items():
        assert counts.index.equals(genes[first].index), (
            f"{sample}: the genes differ from those of {first}: were they aligned with the same STAR index?"
        )

    # columns (sample, strandedness), in the order of the samples
    counts = pd.concat(genes, axis=1).rename_axis("gene_id")
    for strandedness in STRANDEDNESS:
        counts.xs(strandedness, axis=1, level=1).to_csv(f"{prefix}.{strandedness}.tsv", sep="\t", lineterminator="\n")

    # rows (strandedness, metric), columns the samples
    summary = pd.concat({sample: table.T.stack() for sample, table in summaries.items()}, axis=1)
    summary = summary.rename_axis(["strandedness", "metric"])
    summary.to_csv(f"{prefix}.summary.tsv", sep="\t", lineterminator="\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
