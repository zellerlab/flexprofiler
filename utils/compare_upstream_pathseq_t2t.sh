#!/bin/bash
# Compare the PathSeq-T2T mode of flexprofiler (--pathseq_t2t) with upstream PathSeq-T2T
# (https://github.com/abdohlman/pathseq-t2t, commit 51d74430), read by read, on the same alignment to GRCh38.
#
# usage: compare_upstream_pathseq_t2t.sh <flexprofiler_outdir> <sample> <t2t_fasta> <pathseq_host_dir> <decoys_bed> \
#            <kraken2_db> <workdir>
#
# <flexprofiler_outdir> is the output of a flexprofiler run with --pathseq_t2t --save_pathseq_t2t_bams
# --save_analysis_ready_fastqs --run_kraken2 --run_profile_standardisation, and <sample> the name of the sample in its
# files (<sample>_<run_accession> without run merging). Upstream starts from the GRCh38 BAM of flexprofiler
# (pathseq_t2t/grch38/<sample>.grch38.bam), so that only the filtering and the classification are compared.
# <t2t_fasta> must be uncompressed, with its bwa index alongside (built if missing).
#
# Upstream is run as saved in utils/upstream_pathseq_t2t/ (see its README.md): cloned at the commit, with
# t2tfilter.patch applied so that it runs at all, in the conda environment of environment.yml (the tool versions of
# flexprofiler), built in <workdir> with docker and micromamba, with the commands of run_upstream.sh.
#
# Prints a table of the reads after each step (same reads?) and of the Kraken2 counts. Paired-end samples only:
# upstream stops on single-end reads (no paired output of PathSeqFilterSpark, qcfilter.sh:304-314).
set -euo pipefail

if [[ $# -ne 7 ]]; then
    sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'
    exit 1
fi
OUTDIR=$(cd "$1" && pwd)
SAMPLE=$2
T2T_FASTA=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
HOST_DIR=$(cd "$4" && pwd)
DECOYS=$(cd "$(dirname "$5")" && pwd)/$(basename "$5")
KRAKEN2_DB=$(cd "$6" && pwd)
mkdir -p "$7"
WORK=$(cd "$7" && pwd)
COMMIT=51d74430b6a4b34073f2d878612ac24b6a1d1e80
IMAGE=mambaorg/micromamba:1.5.10
RECIPE=$(cd "$(dirname "$0")/upstream_pathseq_t2t" && pwd)

# docker mounts: the work directory, and each input directory read-only
MOUNTS=(-v "$WORK:/c" -v "$RECIPE:/recipe:ro" -v "$OUTDIR:/flex:ro" -v "$(dirname "$T2T_FASTA"):/t2t" -v "$HOST_DIR:/host:ro"
    -v "$(dirname "$DECOYS"):/decoys:ro" -v "$KRAKEN2_DB:/kraken2:ro")
# as root: Spark (Hadoop) needs the user to be known to the container, which the user of the host is not. The files
# are given back to the user of the host at the end of each command
OWNER="$(id -u):$(id -g)"
run() {
    docker run --rm --platform linux/amd64 -u 0:0 "${MOUNTS[@]}" "$IMAGE" \
        bash -c "trap 'chown -R $OWNER /c /t2t/$(basename "$T2T_FASTA").* 2>/dev/null || true' EXIT; $1"
}

if [[ ! -d "$WORK/upstream" ]]; then
    git clone -q https://github.com/abdohlman/pathseq-t2t.git "$WORK/upstream"
    git -C "$WORK/upstream" checkout -q "$COMMIT"
    git -C "$WORK/upstream" apply "$RECIPE/t2tfilter.patch"
fi
if [[ ! -x "$WORK/env/bin/gatk" ]]; then
    run "micromamba create -y -q -p /c/env -f /recipe/environment.yml"
fi

PICARD=/c/env/share/picard-3.4.0-0/picard.jar
T2T=/t2t/$(basename "$T2T_FASTA")
run "export PATH=/c/env/bin:\$PATH JAVA_HOME=/c/env/lib/jvm
    set -e
    [ -e $T2T.bwt ] || bwa index $T2T
    rm -rf /c/out_$SAMPLE
    /recipe/run_upstream.sh /c/upstream $PICARD /flex/pathseq_t2t/grch38/$SAMPLE.grch38.bam $SAMPLE \
        /decoys/$(basename "$DECOYS") /host $T2T /kraken2 /c/out_$SAMPLE 4 8"

cat > "$WORK/compare.py" <<'EOF'
import csv, glob, gzip, subprocess, sys

S, U, O = sys.argv[1], f"/c/out_{sys.argv[1]}", "/flex"


def bam_reads(*paths):
    # the primary records: name, mate (flags 0x40/0x80) and sequence as sequenced
    reads = set()
    for path in paths:
        lines = subprocess.run(["samtools", "view", "-F", "0x900", path], capture_output=True, text=True, check=True).stdout
        for line in lines.splitlines():
            f = line.split("\t")
            flag, seq = int(f[1]), f[9]
            reads.add((f[0], flag & 0xC0, seq[::-1].translate(str.maketrans("ACGTN", "TGCAN")) if flag & 0x10 else seq))
    return reads


def fastq_reads(path, mate):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    # Picard SamToFastq (upstream) adds /1 and /2 to the names
    return {(lines[i][1:].split()[0].removesuffix("/1").removesuffix("/2"), mate, lines[i + 1]) for i in range(0, len(lines), 4)}


def report(step, upstream, ours):
    same = upstream == ours
    detail = "yes" if same else f"no: {len(upstream - ours)} only upstream, {len(ours - upstream)} only flexprofiler"
    print(f"| {step} | {len(upstream)} | {len(ours)} | {detail} |")
    return same


print("| Step | Upstream reads | flexprofiler reads | Same reads? |")
print("|---|---|---|---|")
ok = report("Not aligned to GRCh38 (prefilter)", bam_reads(f"{U}/bams/{S}.prefilter.unaligned.bam"), bam_reads(f"{O}/pathseq_t2t/prefilter/{S}.prefilter.unaligned.bam"))
# flexprofiler writes no decoy BAM if there are no decoy reads
ours_decoys = glob.glob(f"{O}/pathseq_t2t/prefilter/{S}.prefilter.decoys.bam")
ok &= report("Aligned to the decoy contigs (prefilter)", bam_reads(f"{U}/bams/{S}.prefilter.decoys.bam"), bam_reads(*ours_decoys))
for kind in ("paired", "unpaired"):
    ours = [f"{O}/pathseq_t2t/qcfilter/{S}.{part}.qcfilt_{kind}.bam" for part in ("unaligned", "decoys")]
    ours = [path for path in ours if glob.glob(path)]
    ok &= report(f"PathSeqFilterSpark, {kind}", bam_reads(f"{U}/bams/{S}.qcfilt_{kind}.bam"), bam_reads(*ours))
ready = f"{O}/analysis_ready_fastqs/{S}.t2tfilt"
ok &= report("Not aligned to T2T-CHM13, pairs", fastq_reads(f"{U}/bams/{S}.classify.R1.fq.gz", 64) | fastq_reads(f"{U}/bams/{S}.classify.R2.fq.gz", 128), fastq_reads(f"{ready}_1.fastq.gz", 64) | fastq_reads(f"{ready}_2.fastq.gz", 128))
ok &= report("Not aligned to T2T-CHM13, singletons", fastq_reads(f"{U}/bams/{S}.classify.FU.fq.gz", 0), fastq_reads(f"{ready}_unpaired.fastq.gz", 0))
print("\nAll the read sets are the same." if ok else "\nSOME READ SETS DIFFER.")

# Kraken2: the reports of the same reads, then the counts of the merged tables (of the first kraken2 database)
merged = sorted(glob.glob(f"{O}/kraken2/kraken2_*.reads_clade.tsv"))[0]
db = merged.split("/kraken2_", 1)[1].removesuffix(".reads_clade.tsv")
for kind, suffix in (("paired", ""), ("unpaired", "_singletons")):
    up_report = open(f"{U}/classification_stats/{S}.{kind}.kraken.report.txt").read()
    our_report = open(f"{O}/kraken2/{db}/{S}_{db}{suffix}.kraken2.report.txt").read()
    print(f"Kraken2 report of the {kind} reads: {'same' if up_report == our_report else 'DIFFERENT'}")
up = {r["tax_id"]: r for r in csv.DictReader(open(f"{U}/results/{S}.kraken.txt"), delimiter="\t")}
clade = {r["taxid"]: r for r in csv.DictReader(open(merged), delimiter="\t")}
rpm = {r["taxid"]: r for r in csv.DictReader(open(f"{O}/kraken2/kraken2_{db}.rpm_clade.tsv"), delimiter="\t")}
print("\n| Taxon | Upstream reads | flexprofiler reads | Upstream RPM | flexprofiler RPM |")
print("|---|---|---|---|---|")
for taxid, row in sorted(up.items(), key=lambda item: -int(item[1]["reads_clade"]))[:30]:
    print(f"| {row['name']} ({taxid}) | {row['reads_clade']} | {clade[taxid][S]} | {float(row['reads_clade_per_million']):.4g} | {float(rpm[taxid][S]):.4g} |")
EOF
run "export PATH=/c/env/bin:\$PATH; python /c/compare.py $SAMPLE"
