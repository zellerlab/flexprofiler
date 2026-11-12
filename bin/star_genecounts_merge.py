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
  in each sample. The strandedness of a library is the one with the most reads assigned

All the samples must have the same genes, in the same order (the same STAR index).
"""
import csv
import sys

STRANDEDNESS = ["unstranded", "forward", "reverse"]
SUMMARY = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]


def read_counts(path):
    """The summary counts and the gene IDs and counts of a ReadsPerGene.out.tab, each a list of [id, c1, c2, c3]."""
    with open(path, newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t")]
    for i, row in enumerate(rows, 1):
        assert len(row) == 4, f"{path}:{i}: {len(row)} columns, expected 4"
    ids = [row[0] for row in rows[: len(SUMMARY)]]
    assert ids == SUMMARY, f"{path}: summary lines {ids}, expected {SUMMARY}"
    counts = [[row[0]] + [int(value) for value in row[1:]] for row in rows]
    return counts[: len(SUMMARY)], counts[len(SUMMARY) :]


def read_sheet(path):
    """The [sample, ReadsPerGene.out.tab] pairs of the sample sheet."""
    with open(path, newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if row]
    for i, row in enumerate(rows, 1):
        assert len(row) == 2, f"{path}:{i}: {len(row)} columns, expected 2 (sample name, ReadsPerGene.out.tab)"
    assert rows, f"{path}: no samples"
    return rows


def main(prefix, sheet):
    pairs = read_sheet(sheet)
    samples = [sample for sample, _ in pairs]
    duplicates = sorted({sample for sample in samples if samples.count(sample) > 1})
    assert not duplicates, f"duplicate sample names: {', '.join(duplicates)}"

    tables = {sample: read_counts(path) for sample, path in pairs}
    samples = sorted(tables)
    genes = [row[0] for row in tables[samples[0]][1]]
    for sample in samples[1:]:
        assert [row[0] for row in tables[sample][1]] == genes, (
            f"{sample}: the genes differ from those of {samples[0]}: were they aligned with the same STAR index?"
        )

    for column, strandedness in enumerate(STRANDEDNESS, 1):
        with open(f"{prefix}.{strandedness}.tsv", "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id"] + samples)
            for i, gene in enumerate(genes):
                writer.writerow([gene] + [tables[sample][1][i][column] for sample in samples])

    with open(f"{prefix}.summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["strandedness", "metric"] + samples)
        for column, strandedness in enumerate(STRANDEDNESS, 1):
            for i, metric in enumerate(SUMMARY):
                writer.writerow([strandedness, metric] + [tables[sample][0][i][column] for sample in samples])
            assigned = [sum(row[column] for row in tables[sample][1]) for sample in samples]
            writer.writerow([strandedness, "N_assigned"] + assigned)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
