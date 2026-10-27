#!/usr/bin/env python3
"""Give the same count columns to the gqlib profiles of all the samples, before merging them.

usage: gqlib_merge_columns.py <profile.txt.gz>...

gqlib reports the combined_ columns (the counts including the ambiguous reads) only for the
samples with ambiguous reads, as they would be the same as the uniq_ columns otherwise. The
missing combined_ columns are added here as copies of the uniq_ ones, and the profiles are
rewritten in place (a symlink is replaced by a file, its target is not changed).

collate_counts (gqlib 2.14.4) falls back to the uniq_ column for a sample without the combined_
one, but then keeps using the uniq_ column for all the following samples, so it cannot be used
on profiles with different columns.

Prints the count columns of all the profiles (uniq_ then combined_, as gqlib writes them).
"""

import gzip
import os
import sys


def read_header(path):
    with gzip.open(path, "rt") as profile:
        return profile.readline().rstrip("\n").split("\t")


def uniq_column(column):
    return column.replace("combined_", "uniq_", 1)


def main(paths):
    headers = {path: read_header(path) for path in paths}

    columns = []
    for header in headers.values():
        columns += [column for column in header[1:] if column not in columns]
    columns = [c for c in columns if not c.startswith("combined_")] + [c for c in columns if c.startswith("combined_")]

    for path, header in headers.items():
        if header[1:] == columns:
            continue
        missing = [column for column in columns if column not in header[1:]]
        unfillable = [column for column in missing if not column.startswith("combined_") or uniq_column(column) not in header]
        if unfillable:
            sys.exit(f"{path} does not have the columns {', '.join(unfillable)}, which cannot be filled from its uniq_ columns")

        tmp_path = f"{path}.tmp"
        with gzip.open(path, "rt") as profile, gzip.open(tmp_path, "wt") as out:
            next(profile)
            print(header[0], *columns, sep="\t", file=out)
            for line in profile:
                fields = line.rstrip("\n").split("\t")
                values = dict(zip(header[1:], fields[1:]))
                print(fields[0], *(values.get(column, values.get(uniq_column(column))) for column in columns), sep="\t", file=out)
        os.replace(tmp_path, path)

    print(*columns)


if __name__ == "__main__":
    main(sys.argv[1:])
