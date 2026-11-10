#!/usr/bin/env python3
"""Merge the Kraken2 reports of the samples of a database into tables of taxa by samples.

usage: kraken2_merge.py <manifest.tsv> <reports_dir> <read_counts.tsv|-> <prefix>

The manifest has a header and one row per report, with the columns file (its name in <reports_dir>), sample (the
column of the sample in the tables), id and run_accession (the row of the sample in the read counts), and mates: 2 for a report of read pairs, 1 for one of single
reads (single-end reads, or the reads of pairs left without their mate by the PathSeq-T2T filtering). Kraken2 counts a
pair as one fragment: the counts of a report are multiplied by its mates, so that the tables count reads, each mate
counted, and the reports of a sample are summed. This differs from upstream PathSeq-T2T, which multiplies the counts
of the single reads by 2 as well (scripts/pst2t_summarize.py:253-254).

Writes:
- <prefix>.reads_clade.tsv: the reads assigned to each taxon or a taxon below it
- <prefix>.reads_taxon.tsv: the reads assigned to each taxon itself
- <prefix>.taxonomy.tsv: the name, rank and lineage (the taxids from the root) of each taxon
- <prefix>.rpm_clade.tsv, if <read_counts.tsv> is given: the reads of reads_clade per million raw reads of the sample
  (column raw_reads of the read counts, e.g. the primary reads aligned to GRCh38 with --pathseq_t2t), as in
  Dohlman et al., Cell 2026 (https://doi.org/10.1016/j.cell.2026.04.015)

The rows are the taxa (taxid), in the order of the taxonomy, with the unclassified reads (taxid 0) first, and the
columns the samples, sorted. A sample with no reads (an empty report) has a column of zeros.
"""
import pathlib
import sys

import pandas as pd

COLUMNS = ["percent", "reads_clade", "reads_taxon", "minimizers", "distinct_minimizers", "rank", "taxid", "name"]


def read_report(path):
    """The rows of a Kraken2 report (with or without --report-minimizer-data), with the lineage of each taxon."""
    rows = []
    lineage = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) == 6:
            # without --report-minimizer-data
            fields = fields[:3] + ["", ""] + fields[3:]
        assert len(fields) == 8, f"Not a Kraken2 report line in {path}: {line}"
        name = fields[7]
        # the depth of a taxon is the indentation of its name, two spaces per level
        depth = (len(name) - len(name.lstrip(" "))) // 2
        taxid = int(fields[6])
        if taxid == 0:
            taxon_lineage = "0"
        else:
            lineage = lineage[:depth] + [str(taxid)]
            taxon_lineage = ";".join(lineage)
        rows.append(
            {
                "taxid": taxid,
                "name": name.strip(),
                "rank": fields[5],
                "lineage": taxon_lineage,
                "reads_clade": int(fields[1]),
                "reads_taxon": int(fields[2]),
            }
        )
    return pd.DataFrame(rows, columns=["taxid", "name", "rank", "lineage", "reads_clade", "reads_taxon"])


def main(manifest, reports_dir, read_counts, prefix):
    entries = pd.read_csv(
        manifest, sep="\t", dtype={"file": str, "sample": str, "id": str, "run_accession": str, "mates": int}
    )
    assert set(entries["mates"]) <= {1, 2}, f"mates must be 1 or 2: {sorted(set(entries['mates']))}"
    samples = sorted(entries["sample"].unique())

    reports = []
    for entry in entries.itertuples():
        report = read_report(pathlib.Path(reports_dir) / entry.file)
        report[["reads_clade", "reads_taxon"]] *= entry.mates
        report["sample"] = entry.sample
        reports.append(report)
    # an empty report has columns of no type
    reports = pd.concat(reports, ignore_index=True).astype({"taxid": int, "reads_clade": int, "reads_taxon": int})

    # the taxonomy, from all the reports (the same database), in the order of the lineages
    taxonomy = reports.drop_duplicates("taxid")[["taxid", "name", "rank", "lineage"]]
    conflicting = taxonomy["taxid"].duplicated()
    assert not conflicting.any(), "Taxa with different names or ranks in the reports"
    taxonomy = taxonomy.assign(order=taxonomy["taxid"].ne(0)).sort_values(["order", "lineage"]).drop(columns="order")
    taxonomy.to_csv(f"{prefix}.taxonomy.tsv", sep="\t", index=False)

    tables = {}
    for column in ["reads_clade", "reads_taxon"]:
        table = (
            reports.pivot_table(index="taxid", columns="sample", values=column, aggfunc="sum", fill_value=0)
            .reindex(index=taxonomy["taxid"], columns=samples, fill_value=0)
            .astype(int)
        )
        table.columns.name = None
        table.to_csv(f"{prefix}.{column}.tsv", sep="\t")
        tables[column] = table

    if read_counts != "-":
        counts = pd.read_csv(read_counts, sep="\t", dtype={"sample": str, "run_accession": str})
        assert "raw_reads" in counts.columns, (
            f"No raw_reads in {read_counts}: the reads per million need the raw reads of the samples"
        )
        # the raw reads of each column, from the row of its sample and run in the read counts
        raw = (
            entries.drop_duplicates(["sample", "id", "run_accession"])
            .merge(counts.rename(columns={"sample": "id"})[["id", "run_accession", "raw_reads"]], how="left")
        )
        denominators = raw.groupby("sample")["raw_reads"].sum(min_count=1).reindex(samples)
        missing = denominators[denominators.isna() | (denominators == 0)]
        assert missing.empty, f"No raw reads in {read_counts} for {', '.join(missing.index)}"
        (tables["reads_clade"] * 1e6 / denominators).to_csv(f"{prefix}.rpm_clade.tsv", sep="\t", float_format="%.6g")


if __name__ == "__main__":
    assert len(sys.argv) == 5, __doc__
    main(*sys.argv[1:])
