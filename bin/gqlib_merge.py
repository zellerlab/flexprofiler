#!/usr/bin/env python3
"""Merge the gqlib (and cayman) profiles of the samples, for one category, into one table per count column.

usage: gqlib_merge.py <prefix> <suffix> <profiles.txt>

<profiles.txt> lists the profiles to merge (<profile.txt.gz>), one path per line.

The sample name of a profile is its file name without <suffix>. Writes <prefix>.<column>.txt.gz for each count
column of the profiles (uniq_raw, uniq_lnorm, ..., then combined_raw, ...), with the samples as columns.

The tables are as those of collate_counts (gqlib 2.14.4, gqlib/bin/collate_counts.py), which this replaces: the
features of all the samples as rows, total_reads, filtered_reads, category and unannotated first and the others
sorted, the samples sorted by name, NA for a feature missing in a sample, and the same pandas reading and writing.
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
    feature_label = next(iter(profiles.values())).index.name

    # the count columns, the uniq_ ones then the combined_ ones
    columns = list(dict.fromkeys(column for profile in profiles.values() for column in profile.columns))
    columns = sorted(columns, key=lambda column: column.startswith("combined_"))
    features = set().union(*(profile.index for profile in profiles.values()))
    index = [row for row in HEADER_ROWS if row in features] + sorted(features.difference(HEADER_ROWS))

    for column in columns:
        # a sample without the combined_ column has the same counts in the uniq_ one
        table = pd.concat(
            {sample: profile[column if column in profile else column.replace("combined_", "uniq_", 1)] for sample, profile in profiles.items()},
            axis=1,
        )
        table.reindex(index).to_csv(f"{prefix}.{column}.txt.gz", sep="\t", na_rep="NA", index_label=feature_label)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3]).read_text().splitlines())
