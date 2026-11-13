#!/bin/bash
# Build the Bracken k-mer distribution of the Kraken2 database of the WGS test (test_kraken2, with --run_bracken),
# data/kraken2_wgs/database100mers.kmer_distrib (100: the default read length of Bracken), from the same genomes and
# taxonomy as make_kraken2_wgs_db.sh. The Kraken2 database itself (its hash.k2d, opts.k2d and taxo.k2d, not byte-identical
# from one build to the next) is not rebuilt: the genomes are classified with it, as bracken-build does.
# Run from assets/test-data, after make_kraken2_wgs_db.sh, with docker, curl and unzip.
set -euo pipefail
OUT=data/kraken2_wgs
READ_LEN=100
BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT
mkdir -p "$BUILD/db/library" "$BUILD/db/taxonomy"
cp "$OUT"/{hash,opts,taxo}.k2d "$BUILD/db/"

# the genomes, as in make_kraken2_wgs_db.sh, and the taxid of each sequence (seqid2taxid.map, written by kraken2-build)
for genome in GCF_000016525.1:2173 GCF_010509075.1:74426; do
    accession=${genome%%:*}
    taxid=${genome##*:}
    curl -sSfL -o "$BUILD/$accession.zip" \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$accession/download?include_annotation_type=GENOME_FASTA"
    unzip -p "$BUILD/$accession.zip" "ncbi_dataset/data/$accession/*_genomic.fna" \
        | sed "s/^>\([^ ]*\).*/>\1|kraken:taxid|$taxid/" >> "$BUILD/db/library/genomes.fna"
done
grep '^>' "$BUILD/db/library/genomes.fna" | sed 's/^>//; s/|kraken:taxid|\(.*\)$/&\t\1/' > "$BUILD/db/seqid2taxid.map"

# the NCBI taxonomy, kept only for the lineages of the two species, as in make_kraken2_wgs_db.sh
curl -sSfL -o "$BUILD/taxdump.tar.gz" https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf "$BUILD/taxdump.tar.gz" -C "$BUILD" nodes.dmp names.dmp
awk -F'\t[|]\t' '
    NR == FNR { parent[$1] = $2; next }
    FNR == 1 {
        split("2173 74426", species, " ")
        for (s in species) for (t = species[s]; !(t in keep); t = parent[t]) keep[t] = 1
    }
    $1 in keep
' "$BUILD/nodes.dmp" "$BUILD/nodes.dmp" > "$BUILD/db/taxonomy/nodes.dmp"
awk -F'\t[|]\t' 'NR == FNR { keep[$1] = 1; next } ($1 in keep)' \
    "$BUILD/db/taxonomy/nodes.dmp" "$BUILD/names.dmp" > "$BUILD/db/taxonomy/names.dmp"

# the container of the BRACKEN_BRACKEN module (it also has kraken2); the k-mer length is that of the database (35)
docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$BUILD:/w" -w /w \
    community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01 \
    bracken-build -d db -t 2 -k 35 -l "$READ_LEN"
test -s "$BUILD/db/database${READ_LEN}mers.kmer_distrib"
cp "$BUILD/db/database${READ_LEN}mers.kmer_distrib" "$OUT/"
