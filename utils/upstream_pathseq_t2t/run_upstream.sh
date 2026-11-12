#!/bin/bash
# Run upstream PathSeq-T2T (https://github.com/abdohlman/pathseq-t2t, commit 51d74430, with t2tfilter.patch) on one
# BAM aligned to GRCh38 with bwa: filter (prefilter, qcfilter, t2tfilter), classify with Kraken2, and summarize, with
# the options of the paper (Dohlman et al., Cell 2026, https://doi.org/10.1016/j.cell.2026.04.015). See README.md.
#
# usage: run_upstream.sh <upstream_dir> <picard.jar> <grch38.bam> <sample> <decoys.bed|None> <pathseq_host_dir> \
#            <t2t.fa> <kraken2_db> <outdir> [threads] [ram_gb]
#
# Run it with the tools of environment.yml on the PATH, and JAVA_HOME set to their Java 17. <t2t.fa> needs its bwa
# index alongside (bwa index <t2t.fa>). Outputs in <outdir>, as upstream names them: bams/, filter_stats/,
# classification_stats/ and results/<sample>.{summary,kraken}.tsv|txt, with the logs of each command.
set -euo pipefail

if [[ $# -lt 9 || $# -gt 11 ]]; then
    sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'
    exit 1
fi
UPSTREAM=$1 PICARD=$2 BAM=$3 SAMPLE=$4 DECOYS=$5 HOST_DIR=$6 T2T=$7 KRAKEN2_DB=$8 OUT=$9
THREADS=${10:-4}
RAM_GB=${11:-48}

export PATH="$UPSTREAM/src:$PATH"
mkdir -p "$OUT/tmp"

# the three host filters, run in sequence by upstream; --min-clipped-read-length is upstream's default, 60
pathseq-t2t filter \
    --input-bam "$BAM" --aligner bwa --decoys-to-mask "$DECOYS" --sample-id "$SAMPLE" \
    --hostdir "$HOST_DIR" --reference "$T2T" --outdir "$OUT" --threads "$THREADS" \
    --qcfilter-args "--ram-gb $RAM_GB --tmpdir $OUT/tmp --picard-jar $PICARD" \
    --t2tfilter-args "--picard-jar $PICARD" \
    > "$OUT/filter.log" 2>&1
# Kraken2 on the pairs and on the singletons, with --confidence 0.15 --report-minimizer-data (fixed in upstream)
pathseq-t2t classify \
    --sample-id "$SAMPLE" --outdir "$OUT" --classifiers kraken --kraken-index "$KRAKEN2_DB" \
    --threads "$THREADS" --picard-jar "$PICARD" --keep-intermediate \
    > "$OUT/classify.log" 2>&1
# the counts of each step and the reads per million of each taxon
pathseq-t2t summarize --sample-id "$SAMPLE" --outdir "$OUT" --results-dir "$OUT/results" -v \
    > "$OUT/summarize.log" 2>&1
