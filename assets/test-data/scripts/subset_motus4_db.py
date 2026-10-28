"""Subset a mOTUs v4 database (db_mOTU folder) to some marker genes, keeping the format of all its files.

usage: subset_motus4_db.py <db_mOTU> <marker_genes.txt> <out_db_mOTU>

The bwa index of the subset is not built: run `bwa index <out_db_mOTU>/mOTUsv4.0.db.fna.gz` afterwards.
"""
import csv
import gzip
import hashlib
import pathlib
import sys

src, mgs_file, out = pathlib.Path(sys.argv[1]), sys.argv[2], pathlib.Path(sys.argv[3])
keep_mgs = {line.strip() for line in open(mgs_file) if line.strip()}
out.mkdir(parents=True, exist_ok=True)


def read_gz(name):
    return gzip.open(src / name, "rt")


def write_gz(name):
    # no file name and timestamp in the header, so that the output is reproducible
    return gzip.GzipFile(filename="", mode="wb", fileobj=open(out / name, "wb"), mtime=0)


# marker genes (MG), with their marker gene cluster (MGC) and mOTU
keep_mgcs, keep_motus = set(), set()
with read_gz("mOTUsv4.0.map.tsv.gz") as handle, write_gz("mOTUsv4.0.map.tsv.gz") as out_handle:
    out_handle.write(handle.readline().encode())
    for line in handle:
        motu, mgc, _cog, mg = line.split("\t", 4)[:4]
        if mg in keep_mgs:
            out_handle.write(line.encode())
            keep_mgcs.add(mgc)
            keep_motus.add(motu)
print(f"marker genes: {len(keep_mgs)}, marker gene clusters: {len(keep_mgcs)}, mOTUs: {len(keep_motus)}")
if not any("unassigned" in motu for motu in keep_motus):
    sys.exit("The marker genes must include some of the unassigned mOTU, or motus considers the database corrupted")

# marker gene sequences
n_seqs = 0
with read_gz("mOTUsv4.0.db.fna.gz") as handle, write_gz("mOTUsv4.0.db.fna.gz") as out_handle:
    write = False
    for line in handle:
        if line.startswith(">"):
            write = line[1:].split()[0] in keep_mgs
            n_seqs += write
        if write:
            out_handle.write(line.encode())
if n_seqs != len(keep_mgs):
    sys.exit(f"Found {n_seqs} sequences for {len(keep_mgs)} marker genes")

with read_gz("mOTUsv4.0.db.blocklist.gz") as handle, write_gz("mOTUsv4.0.db.blocklist.gz") as out_handle:
    for line in handle:
        if line.strip() in keep_mgs:
            out_handle.write(line.encode())

# tables with one row per mOTU (or genome), with a header
for name, motu_column in [
    ("mOTUsv4.0.gtdb.taxonomy.rep.tsv.gz", "MOTU"),
    ("mOTUsv4.0.gtdb.taxonomy.80mv.tsv.gz", "MOTU"),
    ("mOTUsv4.0.genomes.tsv.gz", "MOTU4"),
]:
    with read_gz(name) as handle, write_gz(name) as out_handle:
        header = handle.readline()
        index = header.rstrip("\n").split("\t").index(motu_column)
        out_handle.write(header.encode())
        for line in handle:
            if line.rstrip("\n").split("\t")[index] in keep_motus:
                out_handle.write(line.encode())

# map of the mOTUs v3.1 to the mOTUs v4 (<motu_v3> <motu_v4>=<fraction>,...), no header
with read_gz("mOTUs3.1_2_mOTUs4.0.map.tsv.gz") as handle, write_gz("mOTUs3.1_2_mOTUs4.0.map.tsv.gz") as out_handle:
    for line in handle:
        targets = {target.split("=")[0] for target in line.rstrip("\n").split("\t")[1].split(",")}
        if targets & keep_motus:
            out_handle.write(line.encode())

# versions file: only the version and date (first two lines) are read by motus, the rest is updated to describe the subset
with open(src / "mOTUsv4.0.db") as handle:
    version, date = handle.readline(), handle.readline()
with gzip.open(out / "mOTUsv4.0.map.tsv.gz", "rt") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
with gzip.open(out / "mOTUsv4.0.genomes.tsv.gz", "rt") as handle:
    n_genomes = sum(1 for _ in handle) - 1
unassigned = [row for row in rows if "unassigned" in row["#MOTU"]]
assigned = [row for row in rows if "unassigned" not in row["#MOTU"]]
counts = {
    "mgs": len(rows),
    "mgcs": len({row["MGC"] for row in rows}),
    "motus": len({row["#MOTU"] for row in assigned}),
    "assigned_mgs": len(assigned),
    "assigned_mgcs": len({row["MGC"] for row in assigned}),
    "assigned_motus": len({row["#MOTU"] for row in assigned}),
    "unassigned_mgs": len(unassigned),
    "unassigned_mgcs": len({row["MGC"] for row in unassigned}),
    "unassigned_motus": len({row["#MOTU"] for row in unassigned}),
    "genomes": n_genomes,
}
with open(out / "mOTUsv4.0.db", "w") as handle:
    handle.write(version + date)
    for key, value in counts.items():
        handle.write(f"{key}: {value}\n")
    for name in [
        "mOTUsv4.0.db.fna.gz",
        "mOTUsv4.0.genomes.tsv.gz",
        "mOTUsv4.0.map.tsv.gz",
        "mOTUsv4.0.db.blocklist.gz",
        "mOTUsv4.0.gtdb.taxonomy.rep.tsv.gz",
        "mOTUsv4.0.gtdb.taxonomy.80mv.tsv.gz",
    ]:
        handle.write(f"{name}: {hashlib.md5((out / name).read_bytes()).hexdigest()}\n")

# marks the database as downloaded for motus
(out / "db_mOTU.downloaded").touch()
print(counts)
