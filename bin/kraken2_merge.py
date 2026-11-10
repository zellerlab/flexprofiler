#!/usr/bin/env python3
"""Merge the Kraken2 reports of the samples of a database into tables of taxa by samples.

usage: kraken2_merge.py <db_name> <pairs_dir> <single_reads_dir> <prefix>

<pairs_dir> has the reports of read pairs, <single_reads_dir> those of single reads (single-end reads, and the reads of
pairs left without their mate by the PathSeq-T2T filtering), named <sample>_<db_name>[_singletons].kraken2.report.txt.
Kraken2 counts a pair as one fragment: the table of the pairs is counted twice and added to that of the single reads,
so that the tables count reads, each mate counted. This differs from upstream PathSeq-T2T, which counts the single
reads twice as well (scripts/pst2t_summarize.py:253-254).

Writes:
- <prefix>.reads_clade.tsv: the reads assigned to each taxon or a taxon below it
- <prefix>.reads_taxon.tsv: the reads assigned to each taxon itself
- <prefix>.taxonomy.tsv: the name, rank and lineage (the taxids from the root) of each taxon

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


def sample_name(path, db_name):
    """The sample of a report, from its name: <sample>_<db_name>[_singletons].kraken2.report.txt."""
    return path.name.removesuffix(".kraken2.report.txt").removesuffix("_singletons").removesuffix(f"_{db_name}")


def main(db_name, pairs_dir, single_reads_dir, prefix):
    reports = []
    for reports_dir, mates in [(pairs_dir, 2), (single_reads_dir, 1)]:
        for path in pathlib.Path(reports_dir).iterdir():
            report = read_report(path)
            report[["reads_clade", "reads_taxon"]] *= mates
            report["sample"] = sample_name(path, db_name)
            reports.append(report)
    # from the file names, as an empty report (a sample with no reads) has no rows
    samples = sorted({sample_name(path, db_name) for d in [pairs_dir, single_reads_dir] for path in pathlib.Path(d).iterdir()})
    # an empty report has columns of no type
    reports = pd.concat(reports, ignore_index=True).astype({"taxid": int, "reads_clade": int, "reads_taxon": int})

    # the taxonomy, from all the reports (the same database), in the order of the lineages
    taxonomy = reports.drop_duplicates("taxid")[["taxid", "name", "rank", "lineage"]]
    conflicting = taxonomy["taxid"].duplicated()
    assert not conflicting.any(), "Taxa with different names or ranks in the reports"
    taxonomy = taxonomy.assign(order=taxonomy["taxid"].ne(0)).sort_values(["order", "lineage"]).drop(columns="order")
    taxonomy.to_csv(f"{prefix}.taxonomy.tsv", sep="\t", index=False)

    for column in ["reads_clade", "reads_taxon"]:
        table = (
            reports.pivot_table(index="taxid", columns="sample", values=column, aggfunc="sum", fill_value=0)
            .reindex(index=taxonomy["taxid"], columns=samples, fill_value=0)
            .astype(int)
        )
        table.columns.name = None
        table.to_csv(f"{prefix}.{column}.tsv", sep="\t")


if __name__ == "__main__":
    assert len(sys.argv) == 5, __doc__
    main(*sys.argv[1:])
