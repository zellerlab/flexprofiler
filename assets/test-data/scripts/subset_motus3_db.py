"""Subset a mOTUs v3 database (db_mOTU folder) to some marker genes, keeping the format of all its files.

usage: subset_motus3_db.py <db_mOTU> <marker_genes.txt> <out_db_mOTU>

The marker genes are those of the NR database (db_mOTU_DB_NR.fasta). The CEN database (one centroid per marker gene
cluster, used for SNV calling) is restricted to the marker gene clusters of the kept marker genes. The bwa indices
are not built: run `bwa index` on the two fastas afterwards. The test data and public profiles are not copied.
"""
import pathlib
import shutil
import sys

src, mgs_file, out = pathlib.Path(sys.argv[1]), sys.argv[2], pathlib.Path(sys.argv[3])
keep_genes = {line.strip() for line in open(mgs_file) if line.strip()}
out.mkdir(parents=True, exist_ok=True)


def filter_table(name, key, header_lines=0, sep="\t"):
    """Copy the rows of a table for which key(fields) is true, and the header lines."""
    with open(src / name) as handle, open(out / name, "w") as out_handle:
        for _ in range(header_lines):
            out_handle.write(handle.readline())
        for line in handle:
            if key(line.rstrip("\n").split(sep)):
                out_handle.write(line)


def filter_fasta(name, keep):
    n = 0
    with open(src / name) as handle, open(out / name, "w") as out_handle:
        write = False
        for line in handle:
            if line.startswith(">"):
                write = line[1:].split()[0] in keep
                n += write
            if write:
                out_handle.write(line)
    return n


# marker genes (NR) and their marker gene clusters (MGC), named <mOTU>.<COG>
keep_mgcs = set()
with open(src / "db_mOTU_MAP_genes_to_MGCs.tsv") as handle:
    for line in handle:
        gene, _length, _cog, mgc = line.rstrip("\n").split("\t")
        if gene in keep_genes:
            keep_mgcs.add(mgc)
keep_motus = set()
with open(src / "db_mOTU_MAP_MGCs_to_mOTUs.tsv") as handle:
    for line in handle:
        mgc, motu = line.rstrip("\n").split("\t")
        if mgc in keep_mgcs:
            keep_motus.add(motu)
print(f"marker genes: {len(keep_genes)}, marker gene clusters: {len(keep_mgcs)}, mOTUs: {len(keep_motus)}")
if "unassigned" not in keep_motus:
    sys.exit("The marker genes must include some of unassigned marker gene clusters")

# NR database, used for profiling
n = filter_fasta("db_mOTU_DB_NR.fasta", keep_genes)
if n != len(keep_genes):
    sys.exit(f"Found {n} sequences for {len(keep_genes)} marker genes")
filter_table("db_mOTU_bam_header_NR", lambda f: f[1].removeprefix("SN:") in keep_genes)
filter_table("db_mOTU_genes_length_NR", lambda f: f[0] in keep_genes)
filter_table("db_mOTU_padding_coordinates_NR.tsv", lambda f: f[0] in keep_genes)
filter_table("db_mOTU_MAP_genes_to_MGCs.tsv", lambda f: f[0] in keep_genes)
filter_table("db_mOTU_MAP_MGCs_to_mOTUs.tsv", lambda f: f[0] in keep_mgcs)
# <mOTU> <MGC>;<MGC>;..., restricted to the kept marker gene clusters, as motus expects to find all of them in the other
# files (with all the marker genes of a mOTU kept, all its marker gene clusters are kept)
with open(src / "db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv") as handle, open(out / "db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv", "w") as out_handle:
    for line in handle:
        motu, mgcs = line.rstrip("\n").split("\t")
        if motu in keep_motus:
            out_handle.write(motu + "\t" + ";".join(mgc for mgc in mgcs.split(";") if mgc in keep_mgcs) + "\n")

# CEN database, one centroid per marker gene cluster, named <mOTU>.<COG>_<gene>
def cen_mgc(cen_gene):
    return cen_gene.rsplit("_", 1)[0]


keep_cen = set()
with open(src / "db_mOTU_bam_header_CEN") as handle:
    for line in handle:
        cen_gene = line.split("\t")[1].removeprefix("SN:")
        if cen_mgc(cen_gene) in keep_mgcs:
            keep_cen.add(cen_gene)
filter_fasta("db_mOTU_DB_CEN.fasta", keep_cen)
filter_table("db_mOTU_DB_CEN.fasta.annotations", lambda f: f[1] in keep_cen)
filter_table("db_mOTU_bam_header_CEN", lambda f: f[1].removeprefix("SN:") in keep_cen)
filter_table("db_mOTU_padding_coordinates_CEN.tsv", lambda f: f[0] in keep_cen)

# taxonomy of the kept mOTUs
filter_table("db_mOTU_taxonomy_ref-mOTUs.tsv", lambda f: f[1] in keep_motus, header_lines=1)
filter_table("db_mOTU_taxonomy_ref-mOTUs_short_names.tsv", lambda f: f[0] in keep_motus, header_lines=1)
filter_table("db_mOTU_taxonomy_meta-mOTUs.tsv", lambda f: f[0] in keep_motus, header_lines=1)
# the first line is a date
filter_table("db_mOTU_taxonomy_CAMI.tsv", lambda f: f[0] in keep_motus, header_lines=1)

for name in ["db_mOTU_versions", "README"]:
    shutil.copy(src / name, out / name)
print(f"CEN marker genes: {len(keep_cen)}")
