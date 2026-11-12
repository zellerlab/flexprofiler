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
import csv
import pathlib
import sys

import pandas as pd

# the columns of a Kraken2 report, with and without --report-minimizer-data
COLUMNS_MINIMIZERS = ["percent_fragments", "fragments_clade", "fragments_taxon", "minimizers", "distinct_minimizers", "rank", "taxid", "name"]
COLUMNS = ["percent_fragments", "fragments_clade", "fragments_taxon", "rank", "taxid", "name"]
OUT_FIELDS = ["taxid", "name", "rank", "lineage", "fragments_clade", "fragments_taxon"]


def read_report(path):
    """The rows of a Kraken2 report (with or without --report-minimizer-data), with the lineage of each taxon."""
    try:
        # names can have quotes and look like missing values
        report = pd.read_csv(path, sep="\t", header=None, quoting=csv.QUOTE_NONE, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=OUT_FIELDS)
    n_columns = len(report.columns)
    assert n_columns == len(COLUMNS_MINIMIZERS) or n_columns == len(COLUMNS), f"{path}: {n_columns} columns, not a Kraken2 report"
    report.columns = COLUMNS_MINIMIZERS if n_columns == len(COLUMNS_MINIMIZERS) else COLUMNS
    return report.assign(name=report["name"].str.strip(), lineage=lineages(report))[OUT_FIELDS]


def lineages(report):
    """The lineage of each taxon of a report (the taxids from the root, joined by ";"), from the tree of the report.

    The report lists the taxa depth-first, with the depth of a taxon as the indentation of its name, two spaces per
    level. The unclassified reads (taxid 0) are at depth 0, like the root, so their lineage is "0".
    """
    depth = report["name"].str.extract(r"^( *)")[0].str.len() // 2
    lineage, result = [], []
    for taxid, d in zip(report["taxid"], depth):
        # the ancestors are the first d taxa of the lineage of the row above
        lineage = lineage[:d] + [str(taxid)]
        result.append(";".join(lineage))
    return result


def sample_name(path, db_name):
    """The sample of a report, from its name: <sample>_<db_name>[_singletons].kraken2.report.txt."""
    return path.name.removesuffix(".kraken2.report.txt").removesuffix("_singletons").removesuffix(f"_{db_name}")


def main(db_name, pairs_dir, single_reads_dir, prefix):
    reports = []
    for reports_dir, mates in [(pairs_dir, 2), (single_reads_dir, 1)]:
        for path in pathlib.Path(reports_dir).iterdir():
            report = read_report(path)
            # the fragments (pairs or single reads) are converted to reads
            for level in ["clade", "taxon"]:
                report[f"reads_{level}"] = report.pop(f"fragments_{level}") * mates
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
            # the pairs and the singletons of a sample are summed, and a taxon absent from a sample counts 0
            reports.pivot_table(index="taxid", columns="sample", values=column, aggfunc="sum", fill_value=0)
            # the rows in the order of the taxonomy (the same taxa), and a column of zeros for the samples with no
            # reads, which have no rows in the reports
            .reindex(index=taxonomy["taxid"], columns=samples, fill_value=0)
            .astype(int)
        )
        table.to_csv(f"{prefix}.{column}.tsv", sep="\t")


if __name__ == "__main__":
    assert len(sys.argv) == 5, __doc__
    main(*sys.argv[1:])
