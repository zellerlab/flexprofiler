#!/usr/bin/env python3
"""Merge the IDTAXA profiles (of IDTAXA_CLASSIFY) of the samples into count tables.

usage: idtaxa_merge.py <prefix> <suffix> <profiles.txt> [<ranks>]

<profiles.txt> lists the profiles to merge (<profile.tsv>), one path per line. <ranks> is a comma-separated list of
the ranks to write a table for, all the ranks of the lineages if not given.

The sample name of a profile is its file name without <suffix>. Writes, with the samples as columns (sorted by name):
- <prefix>.lineage.tsv: the counts of each lineage, as assigned by IDTAXA (from domain to the deepest rank assigned)
- <prefix>.<rank>.tsv, for each rank of the lineages: the counts of the lineages up to that rank
Both with the rows unassigned (the sequences not assigned down to the rank of the table) and total_reads last.
"""
import collections
import pathlib
import sys

import pandas as pd


def cut(lineage, lineage_ranks, rank):
    # the lineage up to rank, found by name, as a taxonomy can skip a rank in some lineages; unassigned if it does not
    # reach it (as the unassigned row of the profiles, which has no ranks)
    lineage_ranks = lineage_ranks.split(";")
    return ";".join(lineage.split(";")[: lineage_ranks.index(rank) + 1]) if rank in lineage_ranks else "unassigned"


def main(prefix, suffix, paths, ranks=None):
    names = [pathlib.Path(path).name for path in paths]
    # the suffix is set by the module from the meta, the file names by the ext.prefix of IDTAXA_CLASSIFY: an error here
    # means that the two do not agree anymore
    wrong = [name for name in names if not name.endswith(suffix) or name == suffix]
    if wrong:
        raise ValueError(f"The profiles {', '.join(wrong)} do not end with {suffix}, the sample names cannot be taken from them")
    samples = [name.removesuffix(suffix) for name in names]
    duplicated = sorted(sample for sample, n in collections.Counter(samples).items() if n > 1)
    if duplicated:
        raise ValueError(f"Duplicated samples: {', '.join(duplicated)}")

    rows = pd.concat(
        pd.read_csv(path, sep="\t", keep_default_na=False).assign(sample=sample) for path, sample in zip(paths, samples)
    )

    # the total of each sample, the same in all the tables
    totals = rows[rows["lineage"] == "total_reads"].set_index("sample")["count"]
    rows = rows[rows["lineage"] != "total_reads"]
    # the ranks to write a table for, by default all of them, in their order in the lineages
    ranks = ranks or list(dict.fromkeys(rank for lineage_ranks in rows["ranks"] if lineage_ranks for rank in lineage_ranks.split(";")))

    # the whole lineages, then the lineages cut at each rank
    tables = {"lineage": rows}
    for rank in ranks:
        tables[rank] = rows.assign(lineage=[cut(lineage, lineage_ranks, rank) for lineage, lineage_ranks in zip(rows["lineage"], rows["ranks"])])

    # one row per lineage, one column per sample
    for name, table_rows in tables.items():
        table = table_rows.pivot_table(index="lineage", columns="sample", values="count", aggfunc="sum", fill_value=0)
        # the samples without any lineage in the table are added, with 0
        table = table.reindex(columns=sorted(samples), fill_value=0)
        # the lineages sorted, then unassigned and total_reads
        table = table.loc[sorted(table.index.drop("unassigned")) + ["unassigned"]]
        table.loc["total_reads"] = totals
        table.to_csv(f"{prefix}.{name}.tsv", sep="\t", index_label="lineage")


if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3]).read_text().splitlines(), sys.argv[4].split(",") if len(sys.argv) == 5 else None)
