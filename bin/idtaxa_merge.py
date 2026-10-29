#!/usr/bin/env python3
"""Merge the IDTAXA profiles (of IDTAXA_CLASSIFY) of the samples into count tables.

usage: idtaxa_merge.py <prefix> <suffix> <profile.tsv>...

The sample name of a profile is its file name without <suffix>. Writes, with the samples as columns (sorted by name):
- <prefix>.lineage.tsv: the counts of each lineage, as assigned by IDTAXA (from domain to the deepest rank assigned)
- <prefix>.<rank>.tsv, for each rank of the lineages: the counts of the lineages up to that rank, the sequences not
  assigned down to it being in the row unassigned
Both with the rows unassigned and total_reads last.
"""
import collections
import csv
import pathlib
import sys

SPECIAL = ["unassigned", "total_reads"]


def read_profile(path):
    with open(path) as handle:
        return [(row["lineage"], row["ranks"], int(row["count"])) for row in csv.DictReader(handle, delimiter="\t")]


def write_table(path, counts, samples):
    rows = sorted(lineage for lineage in counts if lineage not in SPECIAL) + SPECIAL
    with open(path, "w") as handle:
        handle.write("\t".join(["lineage"] + samples) + "\n")
        for lineage in rows:
            handle.write("\t".join([lineage] + [str(counts[lineage].get(s, 0)) for s in samples]) + "\n")


def main(prefix, suffix, paths):
    profiles = {}
    for path in map(pathlib.Path, paths):
        # the suffix is set by the module from the meta, the file names by the ext.prefix of IDTAXA_CLASSIFY: an error
        # here means that the two do not agree anymore
        if not path.name.endswith(suffix) or path.name == suffix:
            sys.exit(f"The profile {path.name} does not end with {suffix}, the sample name cannot be taken from it")
        sample = path.name.removesuffix(suffix)
        if sample in profiles:
            sys.exit(f"Duplicate sample {sample}")
        profiles[sample] = read_profile(path)
    samples = sorted(profiles)

    # the ranks in their order in the lineages
    ranks = []
    for profile in profiles.values():
        for _, lineage_ranks, _ in profile:
            for rank in filter(None, lineage_ranks.split(";")):
                if rank not in ranks:
                    ranks.append(rank)

    lineage_counts = collections.defaultdict(dict)
    rank_counts = {rank: collections.defaultdict(dict) for rank in ranks}
    for sample, profile in profiles.items():
        for lineage, lineage_ranks, count in profile:
            lineage_counts[lineage][sample] = count
            for rank in ranks:
                counts = rank_counts[rank]
                if lineage in SPECIAL:
                    key = lineage
                elif rank in lineage_ranks.split(";"):
                    depth = lineage_ranks.split(";").index(rank) + 1
                    key = ";".join(lineage.split(";")[:depth])
                else:
                    key = "unassigned"
                counts[key][sample] = counts[key].get(sample, 0) + count
    for counts in [lineage_counts] + list(rank_counts.values()):
        for special in SPECIAL:
            counts.setdefault(special, {})

    write_table(f"{prefix}.lineage.tsv", lineage_counts, samples)
    for rank, counts in rank_counts.items():
        write_table(f"{prefix}.{rank}.tsv", counts, samples)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
