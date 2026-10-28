"""Subset a gqlib annotation sqlite to a list of genes, keeping schema and ids unchanged."""
import os
import sqlite3
import sys

src, genes_file, out = sys.argv[1:4]
genes = sorted({line.strip() for line in open(genes_file) if line.strip()})
if os.path.exists(out):
    os.remove(out)

con = sqlite3.connect(out)
con.execute("ATTACH DATABASE ? AS src", (src,))

# same tables and indexes as the source, in creation order
for (sql,) in con.execute("SELECT sql FROM src.sqlite_master WHERE sql IS NOT NULL ORDER BY rowid"):
    con.execute(sql)

con.execute("CREATE TEMP TABLE keep_genes (seqid VARCHAR PRIMARY KEY)")
con.executemany("INSERT INTO keep_genes VALUES (?)", ((g,) for g in genes))
con.execute("INSERT INTO main.annotatedsequence SELECT a.* FROM src.annotatedsequence a JOIN keep_genes k ON a.seqid = k.seqid")

# annotation_str is <category_id>=<feature_id>,<feature_id>;...
feature_ids = set()
for (annot,) in con.execute("SELECT annotation_str FROM main.annotatedsequence"):
    for block in (annot or "").split(";"):
        if "=" in block:
            feature_ids.update(int(f) for f in block.split("=", 1)[1].split(",") if f)
con.execute("CREATE TEMP TABLE keep_features (id INTEGER PRIMARY KEY)")
con.executemany("INSERT INTO keep_features VALUES (?)", ((f,) for f in feature_ids))
con.execute("INSERT INTO main.feature SELECT f.* FROM src.feature f JOIN keep_features k ON f.id = k.id")
con.execute("INSERT INTO main.category SELECT * FROM src.category")
con.commit()
con.execute("DETACH DATABASE src")
con.execute("VACUUM")
# the source is in WAL mode (set by build_gene_database), the subset is not, so that reading it does not create -wal and -shm files next to it
con.execute("PRAGMA journal_mode=delete")

for table in ("annotatedsequence", "feature", "category", "sample", "observation"):
    print(table, con.execute(f"SELECT count(*) FROM {table}").fetchone()[0])
print("genes requested", len(genes), "| features referenced", len(feature_ids))
con.close()
