#!/usr/bin/env python3
"""Merge the gqlib (and cayman) profiles of the samples, for one category, into one table per count column.

usage: gqlib_merge.py <prefix> <suffix> <profiles.txt>

<profiles.txt> lists the profiles to merge (<profile.txt.gz>), one path per line.

The sample name of a profile is its file name without <suffix>. Writes <prefix>.<column>.txt.gz for each count
column of the profiles (uniq_raw, uniq_lnorm, ..., then combined_raw, ...), with the samples as columns.

The tables are as those of collate_counts (gqlib 2.14.4, gqlib/bin/collate_counts.py), which this replaces: the
features of all the samples as rows, total_reads, filtered_reads, category and unannotated first and the others
sorted, the samples sorted by name, and the same pandas reading and writing; but 0 for a feature missing in the
profile of a sample (no reads in it), where collate_counts writes NA.
gqlib writes the combined_ columns (the counts including the ambiguous reads) only for the samples with ambiguous
reads, as they would be the same as the uniq_ ones otherwise: for the other samples, the uniq_ column is used.
collate_counts falls back to the uniq_ column too, but then keeps using it for all the following samples.
"""

import collections
import pathlib
import sys

import pandas as pd

HEADER_ROWS = ["total_reads", "filtered_reads", "category", "unannotated"]


def main(prefix, suffix, paths):
    names = [pathlib.Path(path).name for path in paths]
    # the suffix is set by the module from the meta, the file names by the ext.prefix of the profilers: an error here
    # means that the two do not agree anymore
    wrong = [name for name in names if not name.endswith(suffix) or name == suffix]
    if wrong:
        raise ValueError(f"The profiles {', '.join(wrong)} do not end with {suffix}, the sample names cannot be taken from them")
    samples = [name.removesuffix(suffix) for name in names]
    duplicated = sorted(sample for sample, n in collections.Counter(samples).items() if n > 1)
    if duplicated:
        raise ValueError(f"Duplicated samples: {', '.join(duplicated)}")

    # sample -> profile, one row per feature (gene for gene_counts, feature for the other categories)
    profiles = {sample: pd.read_csv(path, sep="\t", index_col=0) for sample, path in sorted(zip(samples, paths))}
    # the features of all the samples as rows, (sample, count column) as columns
    merged = pd.concat(profiles, axis=1)
    features = set(merged.index)
    index = pd.Index([row for row in HEADER_ROWS if row in features] + sorted(features - set(HEADER_ROWS)), name=merged.index.name)

    # the samples without the combined_ columns have the same counts in the uniq_ ones: the uniq_ columns, renamed
    # combined_, fill them in for the combined_ columns of the other samples (combine_first only fills what is missing)
    names = merged.columns.get_level_values(1)
    as_combined = merged.loc[:, names.str.startswith("uniq_")].rename(columns=lambda column: column.replace("uniq_", "combined_", 1), level=1)
    merged = merged.combine_first(as_combined.loc[:, as_combined.columns.get_level_values(1).isin(names)])

    # a feature missing in the profile of a sample has no reads in it: 0
    for column in merged.columns.unique(level=1):
        merged.xs(column, axis=1, level=1).reindex(index=index, columns=list(profiles)).fillna(0).to_csv(f"{prefix}.{column}.txt.gz", sep="\t")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3]).read_text().splitlines())
