#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=30-00:00:00

# Run flexprofiler on a large dataset, one processing batch at a time (the profile zellerlab_large_scale, see
# docs/zellerlab.md "Large datasets processed in batches"). Every run of the samplesheet must have a processing_batch.
# Each Nextflow run processes the next batch not done (the smallest first), deleting at its start the task directories
# of the batches done before; the run after the last batch makes the tables of all the samples. This job runs the
# pipeline again and again until that is done (OUTDIR/aggregation_inputs/aggregated.json). A run that fails stops the
# job: submitting it again resumes it, the batches done being skipped.
#
# Submit it from the launch directory: sbatch run_shark_large_scale.sh. Stop it with
# `scancel --batch --signal=TERM <job id>`: the signal is passed on to Nextflow, which cancels its tasks (nf-* jobs)
# before exiting. Check afterwards that `squeue -u $USER` shows no nf-* job of the run left.

# set this to the correct path
PARAMS="params.yaml"
# set this to the outdir of the params file
OUTDIR="<THE OUTDIR OF THE PARAMS FILE>"
# set this to the correct path, ideally on the lucid-grpzeller-work folder
WORKDIR="/exports/lucid-grpzeller-work/<YOUR USER NAME>/<PROJECT NAME>/flexprofiler"
# set this to the Nextflow profiles you want to enable: zellerlab_large_scale with those of run_shark.sh (zellerlab_shark,
# which also runs the containers with apptainer, and zellerlab_metagenomics or zellerlab_amplicon)
PROFILE="zellerlab_shark,zellerlab_metagenomics,zellerlab_large_scale"
# be sure that this points to the pipeline version you want to run
VERSION="0.1.1"

# you should not need to edit this
module purge
module load library/java/OpenJDK-20.0.1
module load container/apptainer/1.5.3/gcc-8.5.0
# the newest Nextflow module of SHARK is older than the one the pipeline requires: its launcher downloads (once, to
# ~/.nextflow/framework) and runs the version set by NXF_VER
module load bioinformatics/tools/Nextflow/25.04.6
export NXF_VER="26.10.5"

PIPELINE="zellerlab/flexprofiler"
mkdir -p $WORKDIR
# sbatch is run through the wrapper of the pipeline (utils/sbatch_wrapper), which keeps a task from being submitted
# twice when sbatch times out although the job was submitted; versions without it run with the plain sbatch
nextflow pull $PIPELINE -r $VERSION
SBATCH_WRAPPER_DIR="${NXF_ASSETS:-${NXF_HOME:-$HOME/.nextflow}/assets}/$PIPELINE/utils/sbatch_wrapper"
if [ -x "$SBATCH_WRAPPER_DIR/sbatch" ]; then
    export PATH="$SBATCH_WRAPPER_DIR:$PATH"
fi

# Nextflow runs in the background, so that the SIGTERM of scancel --batch is passed on to it (bash runs a trap only once
# the foreground command exits)
NF_PID=""
trap '[ -n "$NF_PID" ] && kill -TERM "$NF_PID" && wait "$NF_PID"; exit 143' TERM
while true; do
    nextflow run $PIPELINE -r $VERSION -resume -params-file $PARAMS -profile $PROFILE -w $WORKDIR &
    NF_PID=$!
    wait "$NF_PID" || exit $?
    NF_PID=""
    [ -e "$OUTDIR/aggregation_inputs/aggregated.json" ] && break
done
