#!/usr/bin/env python3
"""Merge the Bracken tables of the samples of a Kraken2 database into a table of taxa by samples.

usage: bracken_merge.py <db_name> <pairs_dir> <single_reads_dir> <prefix>

<pairs_dir> has the Bracken tables of read pairs, <single_reads_dir> those of single reads (single-end reads, and the
reads of pairs left without their mate by the PathSeq-T2T filtering), named
<sample>_<db_name>[_singletons].bracken.tsv. Bracken re-estimates the fragments of a Kraken2 report, in which a pair is
one fragment: as for Kraken2 (kraken2_merge.py), the table of the pairs is counted twice and added to that of the single
reads, so that the table counts reads, each mate counted.

Writes:
- <prefix>.reads.tsv: the reads estimated by Bracken for each taxon of its rank (new_est_reads)
- <prefix>.taxonomy.tsv: the name and rank of each taxon

The rows are the taxa (taxid), sorted, and the columns the samples, sorted. A sample with no taxa (no reads at the rank of
the estimation) has a column of zeros. The reads that Bracken does not assign to a taxon of its rank (unclassified, or
classified above it and not redistributed, or of taxa below its threshold) are not in the table.
"""
import pathlib
import sys

import pandas as pd

SUFFIX = ".bracken.tsv"


def sample_name(path, db_name):
    """The sample of a table, from its name: <sample>_<db_name>[_singletons].bracken.tsv."""
    return path.name.removesuffix(SUFFIX).removesuffix("_singletons").removesuffix(f"_{db_name}")


def main(db_name, pairs_dir, single_reads_dir, prefix):
    tables = []
    paths = []
    for tables_dir, mates in [(pairs_dir, 2), (single_reads_dir, 1)]:
        for path in pathlib.Path(tables_dir).glob(f"*{SUFFIX}"):
            paths.append(path)
            # names can look like missing values
            table = pd.read_csv(path, sep="\t", keep_default_na=False, dtype={"taxonomy_id": int, "new_est_reads": int})
            # the fragments (pairs or single reads) are converted to reads
            tables.append(
                pd.DataFrame(
                    {
                        "taxid": table["taxonomy_id"],
                        "name": table["name"],
                        "rank": table["taxonomy_lvl"],
                        "reads": table["new_est_reads"] * mates,
                        "sample": sample_name(path, db_name),
                    }
                )
            )
    # from the file names, as a sample with no taxa has no rows
    samples = sorted({sample_name(path, db_name) for path in paths})
    reads = pd.concat(tables, ignore_index=True).astype({"taxid": int, "reads": int})

    taxonomy = reads.drop_duplicates(["taxid", "name", "rank"])[["taxid", "name", "rank"]].sort_values("taxid")
    assert not taxonomy["taxid"].duplicated().any(), "Taxa with different names or ranks in the Bracken tables"
    taxonomy.to_csv(f"{prefix}.taxonomy.tsv", sep="\t", index=False)

    if reads.empty:
        # no sample has taxa: a table of the samples with no rows
        pd.DataFrame(columns=samples, index=pd.Index([], name="taxid")).to_csv(f"{prefix}.reads.tsv", sep="\t")
        return
    (
        # the pairs and the singletons of a sample are summed, and a taxon absent from a sample counts 0
        reads.pivot_table(index="taxid", columns="sample", values="reads", aggfunc="sum", fill_value=0)
        # the rows of the taxonomy, and a column of zeros for the samples with no taxa
        .reindex(index=taxonomy["taxid"], columns=samples, fill_value=0)
        .astype(int)
        .to_csv(f"{prefix}.reads.tsv", sep="\t")
    )


if __name__ == "__main__":
    assert len(sys.argv) == 5, __doc__
    main(*sys.argv[1:])
