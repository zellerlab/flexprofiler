#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=30-00:00:00

# Run flexprofiler on a large dataset, one processing batch at a time (the profile zellerlab_large_scale, see
# docs/zellerlab.md "Large datasets processed in batches"). Every run of the samplesheet must have a processing_batch.
# Each Nextflow run processes the next batch not done (the smallest first), deleting at its start the task directories
# of the batches done before; the run after the last batch makes the tables of all the samples. One run per job: each
# run that processes a batch submits this script again, as a job that starts once this one has ended successfully
# (--large_scale_slurm_resubmit), until the tables are made. A run that fails ends the chain: submitting the script again
# resumes it, the batches done being skipped.
#
# Submit it from the launch directory: sbatch run_shark_large_scale.sh (options for sbatch go in the #SBATCH lines:
# those of the command line are not passed on to the next jobs). Stop it with `scancel --batch --signal=TERM <job id>`
# of the running job: Nextflow cancels its tasks (nf-* jobs) before exiting, and the next job, if already submitted, is
# removed by SLURM. Check afterwards that `squeue -u $USER` shows no nf-* job of the run left, nor a pending job of
# this script.

# set this to the correct path
PARAMS="params.yaml"
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

# exec: Nextflow gets the SIGTERM of scancel --batch, and cancels its tasks before exiting. Each run that processes a
# batch submits this script again, as a job that starts once this one has ended successfully (the profile
# zellerlab_shark sets --large_scale_slurm_resubmit, used with zellerlab_large_scale only)
exec nextflow run $PIPELINE -r $VERSION -resume -params-file $PARAMS -profile $PROFILE -w $WORKDIR
