#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=3-00:00:00

# set this to the correct path
PARAMS="params.yaml"
# set this to the correct path, ideally on the lucid-grpzeller-work folder
WORKDIR="/exports/lucid-grpzeller-work/<YOUR USER NAME>/<PROJECT NAME>/flexprofiler"
# set this to the Nextflow profiles you want to enable
# you usually need singularity,zellerlab_shark and either zellerlab_metagenomics or zellerlab_16s
PROFILE="singularity,zellerlab_shark,zellerlab_metagenomics"
# be sure that this points to the pipeline version you want to run
VERSION="0.0.1"

# you should not need to edit this
module purge
module load library/java/OpenJDK-20.0.1
module load container/singularity/3.11.4/gcc.8.5.0
module load bioinformatics/tools/Nextflow/25.04.6

PIPELINE="zellerlab/flexprofiler"
mkdir -p $WORKDIR
nextflow run $PIPELINE -r $VERSION -resume -params-file $PARAMS -profile $PROFILE -w $WORKDIR
