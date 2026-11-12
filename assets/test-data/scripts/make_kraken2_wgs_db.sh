#!/bin/bash
# Build the Kraken2 database of the WGS test (test_kraken2), data/kraken2_wgs/, from the NCBI RefSeq reference genomes of
# two species found in the WGS test reads (data/fastq/ERX547493*), with their NCBI taxonomy:
#   Methanobrevibacter smithii ATCC 35061 (GCF_000016525.1), species taxid 2173
#   Collinsella aerofaciens ATCC 25986 (GCF_010509075.1), species taxid 74426
# The taxonomy is the NCBI one (taxdump) cut down to the lineages of the two species, so that the files used are small.
# Run from assets/test-data, with docker, curl and unzip. The database is not byte-identical from one build to the next.
set -euo pipefail
OUT=data/kraken2_wgs
BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT
mkdir -p "$BUILD/db/taxonomy" "$OUT"

# the genomes, each sequence labelled with the species taxid (kraken:taxid)
for genome in GCF_000016525.1:2173 GCF_010509075.1:74426; do
    accession=${genome%%:*}
    taxid=${genome##*:}
    curl -sSfL -o "$BUILD/$accession.zip" \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$accession/download?include_annotation_type=GENOME_FASTA"
    unzip -p "$BUILD/$accession.zip" "ncbi_dataset/data/$accession/*_genomic.fna" \
        | sed "s/^>\([^ ]*\).*/>\1|kraken:taxid|$taxid/" >> "$BUILD/genomes.fna"
done

# the NCBI taxonomy, kept only for the lineages of the two species (nodes.dmp: taxid, parent, rank, ...)
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

docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$BUILD:/w" -w /w \
    community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71 sh -c '
        kraken2-build --add-to-library genomes.fna --db db --no-masking &&
        kraken2-build --build --db db --threads 2'
cp "$BUILD"/db/{hash,opts,taxo}.k2d "$OUT/"
