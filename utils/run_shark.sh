#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=3-00:00:00

# set this to the correct path
PARAMS="params.yaml"
# set this to the correct path, ideally on the lucid-grpzeller-work folder
WORKDIR="/exports/lucid-grpzeller-work/<YOUR USER NAME>/<PROJECT NAME>/flexprofiler"
# set this to the Nextflow profiles you want to enable
# you usually need zellerlab_shark (which also runs the containers with apptainer) and either zellerlab_metagenomics or
# zellerlab_amplicon
PROFILE="zellerlab_shark,zellerlab_metagenomics"
# be sure that this points to the pipeline version you want to run
VERSION="0.1.0"

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
nextflow run $PIPELINE -r $VERSION -resume -params-file $PARAMS -profile $PROFILE -w $WORKDIR
