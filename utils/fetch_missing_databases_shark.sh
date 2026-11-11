#!/bin/bash
#SBATCH --job-name=fetch_flexprofiler_dbs
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=fetch_flexprofiler_dbs_%j.log

# Downloads the databases that are missing on SHARK (Kraken2 and the PathSeq-T2T references) to the shared database
# directory with utils/fetch_databases.sh. Submit from any directory: sbatch utils/fetch_missing_databases_shark.sh
# The PathSeq-T2T paths are already set in the zellerlab_shark profile of conf/zellerlab.config. Afterwards, add the
# kraken2 row of flexprofiler_databases.csv (written in WORKDIR) to assets/zellerlab_shark_databases.csv: it cannot be
# added before, as the paths of the database sheet must exist.

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

# bwa index is single-threaded, so 1 cpu is enough
# the T2T-CHM13v2.0 FASTA of the host removal is the same file as the one of PathSeq-T2T (same md5): it is linked, and
# only its bwa index is built
T2T_FASTA=${DBPATH}/host_genome/t2t_chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz \
FETCH_KRAKEN2=1 FETCH_PATHSEQ_T2T=1 bash ${SCRIPT} ${DBPATH} ${SLURM_CPUS_PER_TASK:-1}

echo "Database sheet rows in ${WORKDIR}/flexprofiler_databases.csv:"
cat ${WORKDIR}/flexprofiler_databases.csv
