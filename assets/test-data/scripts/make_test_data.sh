#!/bin/bash
# Regenerate the test data of the pipeline, calling the scripts of this folder in order (see ../README.md for each step).
#
# usage: make_test_data.sh [options] [step|group ...]
#
# The steps run in the order below, whatever the order given. Without steps, the `local` group runs.
#
# groups:
#   local      suffixed unordered mixed hostremoval adapters amplicon idtaxa
#              made only from files in the repo (needs python3, and Rscript with DECIPHER for idtaxa)
#   download   hostref ref16s
#              downloaded from NCBI, which may change over time and then change the files made from them
#   databases  motus4 motus3 gmgc10 cayman gqlib
#              subsets of the full databases, given with the options below (needs bwa, and samtools for motus4)
#   all        all of the above
#
# steps:
#   suffixed     data/fastq/ERX5474937_ERR5766181_suffixed_{1,2}.fastq.gz
#   unordered    data/fastq/ERX5474937_ERR5766181_unordered_2.fastq.gz
#   mixed        symbolic links in data/fastq/mixed/
#   hostref      data/hostremoval/NC_012920.1.fna.gz
#   hostremoval  data/fastq/hostremoval/ (from hostref)
#   adapters     data/fastq/adapters/adapters_*.fastq.gz (from hostref)
#   ref16s       data/idtaxa/reference_16s.fasta.gz
#   amplicon     data/fastq/amplicon/ (from ref16s)
#   idtaxa       data/idtaxa/trainingset.RData (from ref16s)
#   motus4       data/motus4/marker_genes.txt and data/motus4/db_mOTU/ (--motus4-db)
#   motus3       data/motus3/marker_genes.txt and data/motus3/db_mOTU/ (--motus3-db)
#   gmgc10       data/gmgc10/ fasta and bwa index, from data/gmgc10/genes.txt (--gmgc10-fna)
#   cayman       data/cayman/ (--cayman-csv)
#   gqlib        data/gqlib/ (--gqlib-db)
#
# options:
#   --motus4-db <dir>    full mOTUs v4.0 db_mOTU folder, with its bwa index
#   --motus3-db <dir>    full mOTUs v3.1.0 db_mOTU folder
#   --gmgc10-fna <file>  GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz
#   --cayman-csv <file>  GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv
#   --gqlib-db <file>    GMGC10.human-gut.95nr.0.5pv.20230422.2101.sqlite3
#
# The original reads (data/fastq/ERX*, ERR3201952) and fastas (data/fasta/) come from nf-core/taxprofiler and are not
# regenerated, nor are the hand-written files (data/fastq/adapters/adapters.fasta, data/gmgc10/genes.txt).
set -euo pipefail

STEPS=(suffixed unordered mixed hostref hostremoval adapters ref16s amplicon idtaxa motus4 motus3 gmgc10 cayman gqlib)
LOCAL=(suffixed unordered mixed hostremoval adapters amplicon idtaxa)
DOWNLOAD=(hostref ref16s)
DATABASES=(motus4 motus3 gmgc10 cayman gqlib)

usage() {
    sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'
    exit "${1:-0}"
}

absolute() {
    [[ -e $1 ]] || { echo "Not found: $1" >&2; exit 1; }
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

MOTUS4_DB="" MOTUS3_DB="" GMGC10_FNA="" CAYMAN_CSV="" GQLIB_DB=""
requested=()
while (( $# )); do
    case $1 in
        -h|--help) usage ;;
        --motus4-db) MOTUS4_DB=$(absolute "$2"); shift 2 ;;
        --motus3-db) MOTUS3_DB=$(absolute "$2"); shift 2 ;;
        --gmgc10-fna) GMGC10_FNA=$(absolute "$2"); shift 2 ;;
        --cayman-csv) CAYMAN_CSV=$(absolute "$2"); shift 2 ;;
        --gqlib-db) GQLIB_DB=$(absolute "$2"); shift 2 ;;
        local) requested+=("${LOCAL[@]}"); shift ;;
        download) requested+=("${DOWNLOAD[@]}"); shift ;;
        databases) requested+=("${DATABASES[@]}"); shift ;;
        all) requested+=("${STEPS[@]}"); shift ;;
        -*) echo "Unknown option: $1" >&2; usage 1 ;;
        *)
            [[ " ${STEPS[*]} " == *" $1 "* ]] || { echo "Unknown step: $1" >&2; usage 1; }
            requested+=("$1"); shift ;;
    esac
done
(( ${#requested[@]} )) || requested=("${LOCAL[@]}")

# the database steps need their full database, checked before running anything
for step in "${requested[@]}"; do
    case $step in
        motus4) [[ -n $MOTUS4_DB ]] || { echo "motus4 needs --motus4-db" >&2; exit 1; } ;;
        motus3) [[ -n $MOTUS3_DB ]] || { echo "motus3 needs --motus3-db" >&2; exit 1; } ;;
        gmgc10) [[ -n $GMGC10_FNA ]] || { echo "gmgc10 needs --gmgc10-fna" >&2; exit 1; } ;;
        cayman) [[ -n $CAYMAN_CSV ]] || { echo "cayman needs --cayman-csv" >&2; exit 1; } ;;
        gqlib) [[ -n $GQLIB_DB ]] || { echo "gqlib needs --gqlib-db" >&2; exit 1; } ;;
    esac
done

cd "$(dirname "$0")/.."
SCRATCH=$(mktemp -d)
trap 'rm -rf "$SCRATCH"' EXIT

step_suffixed() {
    for i in 1 2; do
        gzip -cd data/fastq/ERX5474937_ERR5766181_$i.fastq.gz \
            | awk -v s="/$i" 'NR % 4 == 1 { n = split($0, f, " "); name = f[1] s; comment = substr($0, length(f[1]) + 2); sub("/[12]$", "", comment); print (comment == "" ? name : name " " comment); next } { print }' \
            | gzip -n > data/fastq/ERX5474937_ERR5766181_suffixed_$i.fastq.gz
    done
}

step_unordered() {
    gzip -cd data/fastq/ERX5474937_ERR5766181_2.fastq.gz \
        | paste - - - - \
        | LC_ALL=C sort -t $'\t' -k2,2 -k1,1 \
        | tr '\t' '\n' \
        | gzip -n > data/fastq/ERX5474937_ERR5766181_unordered_2.fastq.gz
}

step_mixed() {
    mkdir -p data/fastq/mixed
    while read -r link target; do
        ln -sfn "../$target" "data/fastq/mixed/$link"
    done <<'EOF'
2612_mixed_pe_1.fastq.gz ERX5474932_ERR5766176_1.fastq.gz
2612_mixed_pe_2.fastq.gz ERX5474932_ERR5766176_2.fastq.gz
2612_mixed_se_1.fastq.gz ERX5474936_ERR5766180_1.fastq.gz
2613_mixed_suffixed_1.fastq.gz ERX5474937_ERR5766181_suffixed_1.fastq.gz
2613_mixed_suffixed_2.fastq.gz ERX5474937_ERR5766181_suffixed_2.fastq.gz
2613_mixed_unordered_1.fastq.gz ERX5474937_ERR5766181_1.fastq.gz
2613_mixed_unordered_2.fastq.gz ERX5474937_ERR5766181_unordered_2.fastq.gz
2613_unordered_1.fastq.gz ERX5474937_ERR5766181_1.fastq.gz
2613_unordered_2.fastq.gz ERX5474937_ERR5766181_unordered_2.fastq.gz
EOF
}

step_hostref() {
    mkdir -p data/hostremoval
    curl -sf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" \
        | gzip -n > data/hostremoval/NC_012920.1.fna.gz
}

step_hostremoval() {
    python3 scripts/make_host_contaminated_reads.py \
        data/fastq/ERX5474937_ERR5766181_1.fastq.gz data/fastq/ERX5474937_ERR5766181_2.fastq.gz \
        data/hostremoval/NC_012920.1.fna.gz data/fastq/hostremoval
}

step_adapters() {
    python3 scripts/make_adapter_reads.py data/hostremoval/NC_012920.1.fna.gz data/fastq/adapters
}

step_ref16s() {
    mkdir -p data/idtaxa
    python3 scripts/fetch_16s_reference.py data/idtaxa/reference_16s.fasta.gz
}

step_amplicon() {
    python3 scripts/make_amplicon_reads.py data/idtaxa/reference_16s.fasta.gz data/fastq/amplicon
}

step_idtaxa() {
    Rscript scripts/make_idtaxa_trainingset.R data/idtaxa/reference_16s.fasta.gz data/idtaxa/trainingset.RData
}

step_motus4() {
    # marker genes hit by the reads of the samples profiled by motus (2612, 2612_se, 2613)
    for reads in data/fastq/ERX547493{2,6,7}_*.fastq.gz; do
        scripts/motus_marker_genes_hit.sh "$MOTUS4_DB/mOTUsv4.0.db.fna.gz" "$SCRATCH/$(basename "$reads" .fastq.gz).txt" "$reads"
    done
    cat "$SCRATCH"/*.txt | sort -u > "$SCRATCH/marker_genes_hit.txt"
    gzip -dc "$MOTUS4_DB/mOTUsv4.0.map.tsv.gz" \
        | awk -F'\t' 'NR==FNR {hit[$1]=1; next} ($4 in hit)' "$SCRATCH/marker_genes_hit.txt" - > "$SCRATCH/map_hit.tsv"
    {
        awk -F'\t' '$1 ~ /^mOTUv4.0_(014686|000924|001912|014998)$/ {print $4}' "$SCRATCH/map_hit.tsv"
        awk -F'\t' '$1 == "mOTUv4.0_unassigned" {print $4}' "$SCRATCH/map_hit.tsv" | sort | head -50
    } | sort -u > data/motus4/marker_genes.txt
    python3 scripts/subset_motus4_db.py "$MOTUS4_DB" data/motus4/marker_genes.txt data/motus4/db_mOTU
    bwa index data/motus4/db_mOTU/mOTUsv4.0.db.fna.gz
}

step_motus3() {
    {
        awk -F'\t' '$4 ~ /^(ref_mOTU_v31_03626|ref_mOTU_v31_03695|ext_mOTU_v31_17631|ref_mOTU_v31_03379)\./ {print $1}' "$MOTUS3_DB/db_mOTU_MAP_genes_to_MGCs.tsv"
        awk -F'\t' '$4 ~ /^unassigned\./ {print $1}' "$MOTUS3_DB/db_mOTU_MAP_genes_to_MGCs.tsv" | head -50
    } > data/motus3/marker_genes.txt
    python3 scripts/subset_motus3_db.py "$MOTUS3_DB" data/motus3/marker_genes.txt data/motus3/db_mOTU
    bwa index data/motus3/db_mOTU/db_mOTU_DB_NR.fasta
    bwa index data/motus3/db_mOTU/db_mOTU_DB_CEN.fasta
}

step_gmgc10() {
    gzip -dc "$GMGC10_FNA" \
        | awk 'NR==FNR {keep[">"$1]=1; next} /^>/ {p = ($1 in keep)} p' data/gmgc10/genes.txt - \
        | gzip -n > data/gmgc10/GMGC10.human-gut.95nr.0.5.percent.prevalence.test_subset.fna.gz
    bwa index data/gmgc10/GMGC10.human-gut.95nr.0.5.percent.prevalence.test_subset.fna.gz
}

step_cayman() {
    awk -F, 'NR==FNR {keep[$1]=1; next} FNR==1 || ($1 in keep)' data/gmgc10/genes.txt "$CAYMAN_CSV" \
        > data/cayman/GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.test_subset.csv
}

step_gqlib() {
    python3 scripts/subset_gqlib_db.py "$GQLIB_DB" data/gmgc10/genes.txt \
        data/gqlib/GMGC10.human-gut.95nr.0.5pv.20230422.2101.test_subset.sqlite3
}

for step in "${STEPS[@]}"; do
    [[ " ${requested[*]} " == *" $step "* ]] || continue
    echo "== $step" >&2
    "step_$step"
done
