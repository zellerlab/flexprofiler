#!/bin/bash
#SBATCH --job-name=fetch_flexprofiler_dbs
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=fetch_flexprofiler_dbs_%j.log

# Downloads the database that is missing on SHARK (Kraken2) to the shared database directory with
# utils/fetch_databases.sh. Submit from any directory: sbatch utils/fetch_missing_databases_shark.sh
# The bwa indexes and host files of PathSeq-T2T are already there (FETCH_PATHSEQ_T2T=1, job 25996578). The STAR index
# of GRCh38.d1.vd1 (FETCH_PATHSEQ_T2T_STAR=1, for --pathseq_t2t_grch38_aligner star) is not built yet: add it to the
# command below when needed, with --cpus-per-task=8 and --mem=64G (STAR needs ~40 GB). The PathSeq-T2T paths are
# already set in the zellerlab_shark profile of conf/zellerlab.config. Afterwards, add the kraken2 row of flexprofiler_databases.csv (written in WORKDIR) to
# assets/zellerlab_shark_databases.csv: it cannot be added before, as the paths of the database sheet must exist.

set -euo pipefail

DBPATH=/exports/archive/lucid-grpzeller-primary/SHARED/DATA/flexprofiler_databases
# the temporary files (the Kraken2 archive alone is ~70 GB) go here and not to the /tmp of the node
WORKDIR=/exports/lucid-grpzeller-work/${USER}/fetch_flexprofiler_databases

# micromamba is needed by fetch_databases.sh: set these if it is not in the PATH of the job
export PATH=${HOME}/.local/bin:${PATH}
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-${HOME}/micromamba}
command -v micromamba > /dev/null || { echo "micromamba not found" >&2 ; exit 1 ; }

# fetch_databases.sh is next to this script: under sbatch $0 is a copy in the spool directory of slurm, so the path of
# the submitted script is taken from slurm
if [[ -n ${SLURM_JOB_ID:-} ]] ; then
    SELF=$(scontrol show job ${SLURM_JOB_ID} | grep -oP '^\s*Command=\K\S+')
else
    SELF=$0
fi
SCRIPT=$(realpath "$(dirname "${SELF}")/fetch_databases.sh")
mkdir -p ${WORKDIR}/tmp
export TMPDIR=${WORKDIR}/tmp
cd ${WORKDIR}

# fetch_databases.sh is run from a copy: bash reads a script as it runs it, so a change of the repository during the
# job (e.g. a git pull) would otherwise change the rest of the run (this stopped job 25996578 after PathSeq-T2T)
cp ${SCRIPT} ${WORKDIR}/fetch_databases.${SLURM_JOB_ID:-$$}.sh

FETCH_KRAKEN2=1 bash ${WORKDIR}/fetch_databases.${SLURM_JOB_ID:-$$}.sh ${DBPATH} ${SLURM_CPUS_PER_TASK:-1}
rm ${WORKDIR}/fetch_databases.${SLURM_JOB_ID:-$$}.sh

echo "Database sheet rows in ${WORKDIR}/flexprofiler_databases.csv:"
cat ${WORKDIR}/flexprofiler_databases.csv
