#!/bin/bash
#SBATCH --job-name=flexprofiler
#SBATCH --partition=cpu-skylake
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=7-00:00:00
# Run flexprofiler on the ALICE cluster with the inputs and results on SHARK: the same params file as for run_shark.sh
# (with SHARK paths for input and outdir) is used. As one job on ALICE, it:
#
#   1. copies the params file and the samplesheet from SHARK, and the input files of the samplesheet (fastq_1,
#      fastq_2) to RUN_DIR/inputs, keeping their SHARK paths under it (a path P on SHARK is RUN_DIR/inputs/P)
#   2. runs the pipeline on them, with the results in RUN_DIR/results
#   3. if the run succeeded, copies the results to outdir on SHARK
#
# Nothing is deleted on ALICE: check with utils/check_transfer_alice.sh that the results are all on SHARK, then delete
# RUN_DIR by hand.
#
# The files are copied with JOBS rsync processes at once (one stream through the ALICE gateway reaches ~40-75 MB/s).
# Files already copied are skipped, so after a failure or a change the job can simply be submitted again (with
# -resume, the tasks already run are reused). The databases must have been copied to ALICE first, with
# utils/sync_databases_alice.sh on SHARK.
#
# Submit it on ALICE, or from SHARK without logging in (the log, slurm-<job id>.out, goes to the directory given to cd):
#
#   ssh alice 'cd /projects/pi-zeller/<YOUR USER NAME> && sbatch' < run_alice.sh
#
# Stop it with: scancel --batch --signal=TERM <job id> (on ALICE). The script passes the signal on to Nextflow, which
# cancels its tasks; check afterwards that `squeue -u $USER` shows no nf-* jobs left.
set -euo pipefail

# set this to the params file on SHARK, as for run_shark.sh: input and outdir are SHARK paths
PARAMS="/exports/archive/lucid-grpzeller-primary/<YOUR USER NAME>/projects/<PROJECT NAME>/params.yaml"
# set this to the directory of the run on ALICE: it holds the staged inputs, the work directory and the results
RUN_DIR="/projects/pi-zeller/<YOUR USER NAME>/<PROJECT NAME>/flexprofiler"
# set this to the Nextflow profiles you want to enable: zellerlab_alice (which also runs the containers with apptainer)
# and usually either zellerlab_metagenomics or zellerlab_amplicon
PROFILE="zellerlab_alice,zellerlab_metagenomics"
# be sure that this points to the pipeline version you want to run
VERSION="0.1.1"

# you should not need to edit this
SHARK_HOST="shark"
JOBS=4
PIPELINE="zellerlab/flexprofiler"
# the Nextflow module of ALICE (25.10.2) is older than the version the pipeline requires, and it is a self-contained
# build that ignores NXF_VER: the standard launcher shared by the lab is used instead, which downloads (once, to
# ~/.nextflow/framework) and runs the version set by NXF_VER
NEXTFLOW="/projects/pi-zeller/SHARED/flexprofiler/nextflow"
export NXF_VER="26.10.5"
module purge
module load ALICE/default Java/21

INPUTS_DIR="$RUN_DIR/inputs"
RESULTS_DIR="$RUN_DIR/results"
WORKDIR="$RUN_DIR/work"
FILES_DIR="$RUN_DIR/run_files"
mkdir -p "$INPUTS_DIR" "$RESULTS_DIR" "$WORKDIR" "$FILES_DIR"

if [ ! -x "$NEXTFLOW" ]; then
    echo "The Nextflow launcher is not in $NEXTFLOW: install it there with" >&2
    echo "  cd $(dirname "$NEXTFLOW") && curl -fsSL https://get.nextflow.io | bash" >&2
    exit 1
fi

# copy the files of LIST (one path per line, relative to SRC) from SRC to DEST with JOBS rsync processes, then check
# with a dry run that nothing is left to copy. SRC or DEST can be on SHARK (host:path); the other arguments are added to
# the rsync options
parallel_rsync() {
    local list="$1" src="$2" dest="$3"
    shift 3
    local opts=(-a --partial "$@")
    local parts_dir part pid failed=0 left
    local pids=()
    parts_dir=$(mktemp -d -p "$FILES_DIR")
    split -n "r/$JOBS" "$list" "$parts_dir/part_"
    for part in "$parts_dir"/part_*; do
        if [ -s "$part" ]; then
            rsync "${opts[@]}" --files-from="$part" "$src" "$dest" &
            pids+=("$!")
        fi
    done
    for pid in "${pids[@]}"; do
        wait "$pid" || failed=1
    done
    rm -rf "$parts_dir"
    if [ "$failed" -ne 0 ]; then
        echo "Some rsync processes failed copying $src to $dest" >&2
        return 1
    fi
    left=$(rsync "${opts[@]}" --dry-run --itemize-changes --files-from="$list" "$src" "$dest" | grep -c '^[<>c]' || true)
    if [ "$left" -ne 0 ]; then
        echo "$left files still differ between $src and $dest" >&2
        return 1
    fi
}

# the value of a top-level key of the params file (quotes and comments removed)
yaml_value() {
    sed -n "s/^$1:[[:space:]]*//p" "$2" | sed 's/[[:space:]]*#.*$//; s/^["'\'']//; s/["'\'']$//'
}

# 1. params file, samplesheet and inputs
rm -f "$FILES_DIR/results_copied"
rsync -a "$SHARK_HOST:$PARAMS" "$FILES_DIR/params.shark.yaml"
SHARK_INPUT=$(yaml_value input "$FILES_DIR/params.shark.yaml")
SHARK_OUTDIR=$(yaml_value outdir "$FILES_DIR/params.shark.yaml")
if [[ "$SHARK_INPUT" != /* || "$SHARK_OUTDIR" != /* ]]; then
    echo "input and outdir must be set to absolute SHARK paths in $PARAMS" >&2
    exit 1
fi
rsync -a "$SHARK_HOST:$SHARK_INPUT" "$FILES_DIR/samplesheet.shark.csv"

# the input columns of the samplesheet, found by name
INPUT_COLUMNS='$i == "fastq_1" || $i == "fastq_2"'
tr -d '\r' < "$FILES_DIR/samplesheet.shark.csv" \
    | awk -F, "NR == 1 { for (i = 1; i <= NF; i++) if ($INPUT_COLUMNS) c[i] = 1; next }
               { for (i in c) if (\$i != \"\") print \$i }" \
    | sort -u > "$FILES_DIR/inputs.txt"
if grep -qv '^/' "$FILES_DIR/inputs.txt"; then
    echo "The input files of the samplesheet must be absolute SHARK paths" >&2
    exit 1
fi
tr -d '\r' < "$FILES_DIR/samplesheet.shark.csv" \
    | awk -F, -v OFS=, -v root="$INPUTS_DIR" "NR == 1 { for (i = 1; i <= NF; i++) if ($INPUT_COLUMNS) c[i] = 1; print; next }
               { for (i in c) if (\$i != \"\") \$i = root \$i; print }" \
    > "$FILES_DIR/samplesheet.csv"
sed -e "s#^input:.*#input: \"$FILES_DIR/samplesheet.csv\"#" \
    -e "s#^outdir:.*#outdir: \"$RESULTS_DIR\"#" \
    "$FILES_DIR/params.shark.yaml" > "$FILES_DIR/params.yaml"
if grep -q '/exports/' "$FILES_DIR/params.yaml"; then
    echo "Warning: the params file has other SHARK paths, which do not exist on ALICE:" >&2
    grep '/exports/' "$FILES_DIR/params.yaml" >&2
fi

echo "Staging $(wc -l < "$FILES_DIR/inputs.txt") input files from $SHARK_HOST to $INPUTS_DIR"
# the inputs are often symlinks (e.g. in SHARED/DATA), which would be broken on ALICE: the files they point to are copied
parallel_rsync "$FILES_DIR/inputs.txt" "$SHARK_HOST:/" "$INPUTS_DIR/" --copy-links

# 2. the pipeline, in the background so that a SIGTERM to this script (scancel --batch --signal=TERM) is passed on to
# Nextflow, which then cancels its tasks; the script waits for it to finish
cd "$RUN_DIR"
# sbatch is run through the wrapper of the pipeline (utils/sbatch_wrapper), which keeps a task from being submitted
# twice when sbatch times out although the job was submitted; versions without it run with the plain sbatch
"$NEXTFLOW" pull "$PIPELINE" -r "$VERSION"
SBATCH_WRAPPER_DIR="${NXF_ASSETS:-${NXF_HOME:-$HOME/.nextflow}/assets}/$PIPELINE/utils/sbatch_wrapper"
if [ -x "$SBATCH_WRAPPER_DIR/sbatch" ]; then
    export PATH="$SBATCH_WRAPPER_DIR:$PATH"
fi
"$NEXTFLOW" run "$PIPELINE" -r "$VERSION" -resume -params-file "$FILES_DIR/params.yaml" -profile "$PROFILE" -w "$WORKDIR" &
NF_PID=$!
trap 'kill -TERM "$NF_PID" 2> /dev/null || true' TERM
NF_STATUS=0
wait "$NF_PID" || NF_STATUS=$?
# wait returns early when the trap runs: wait again until Nextflow has exited
while kill -0 "$NF_PID" 2> /dev/null; do
    NF_STATUS=0
    wait "$NF_PID" || NF_STATUS=$?
done
trap - TERM
if [ "$NF_STATUS" -ne 0 ]; then
    echo "Nextflow exited with status $NF_STATUS: the results are not copied to SHARK" >&2
    exit "$NF_STATUS"
fi

# 3. the results, to SHARK
(cd "$RESULTS_DIR" && find . \( -type f -o -type l \) -printf '%P\n') > "$FILES_DIR/results.txt"
echo "Copying $(wc -l < "$FILES_DIR/results.txt") result files to $SHARK_HOST:$SHARK_OUTDIR"
ssh "$SHARK_HOST" mkdir -p "$SHARK_OUTDIR"
parallel_rsync "$FILES_DIR/results.txt" "$RESULTS_DIR/" "$SHARK_HOST:$SHARK_OUTDIR/"

# the results are on SHARK: utils/check_transfer_alice.sh passes only if this file is there
printf '%s\n' "$SHARK_OUTDIR" > "$FILES_DIR/results_copied"
echo "Done: the results are in $SHARK_HOST:$SHARK_OUTDIR"
echo "Before deleting $RUN_DIR, check the copy with: utils/check_transfer_alice.sh $RUN_DIR"
