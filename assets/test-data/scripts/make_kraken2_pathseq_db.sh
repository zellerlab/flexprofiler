#!/bin/bash
# Build the Kraken2 database of the PathSeq-T2T test, data/kraken2_pathseq/, from the genomes of make_pathseq_t2t_data.py
# (data/pathseq_t2t/microbe_genomes.fna.gz), with a taxonomy of a few made-up taxa:
#   1 root > 2 Bacteria (superkingdom) > 900001 Testella (genus) > 900002 Testella microbica (species): microbe_genome
#   1 root > 10239 Viruses (superkingdom) > 900003 Decoyvirus testi (species): decoy_virus
# Run from assets/test-data, with docker. The database is not byte-identical from one build to the next.
set -euo pipefail
OUT=data/kraken2_pathseq
BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT
mkdir -p "$BUILD/db/taxonomy" "$OUT"

# nodes.dmp and names.dmp, in the NCBI format: fields separated by "\t|\t", lines ending with "\t|"
d=$'\t|\t'
cat > "$BUILD/db/taxonomy/nodes.dmp" <<END
1${d}1${d}no rank${d}	|
2${d}1${d}superkingdom${d}	|
10239${d}1${d}superkingdom${d}	|
900001${d}2${d}genus${d}	|
900002${d}900001${d}species${d}	|
900003${d}10239${d}species${d}	|
END
cat > "$BUILD/db/taxonomy/names.dmp" <<END
1${d}root${d}${d}scientific name	|
2${d}Bacteria${d}${d}scientific name	|
10239${d}Viruses${d}${d}scientific name	|
900001${d}Testella${d}${d}scientific name	|
900002${d}Testella microbica${d}${d}scientific name	|
900003${d}Decoyvirus testi${d}${d}scientific name	|
END
# the taxon of each genome, in the header (kraken:taxid)
gzip -dc data/pathseq_t2t/microbe_genomes.fna.gz \
    | sed 's/^>microbe_genome$/>microbe_genome|kraken:taxid|900002/; s/^>decoy_virus$/>decoy_virus|kraken:taxid|900003/' \
    > "$BUILD/genomes.fna"

docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$BUILD:/w" -w /w \
    community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71 sh -c '
        kraken2-build --add-to-library genomes.fna --db db --no-masking &&
        kraken2-build --build --db db --threads 2 --kmer-len 25 --minimizer-len 21 --minimizer-spaces 3'
cp "$BUILD"/db/{hash,opts,taxo}.k2d "$OUT/"
