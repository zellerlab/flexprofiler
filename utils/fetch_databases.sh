#!/bin/bash

# fetch_databases.sh <out_database_path> [cpus]
# Author: Saul Pierotti
# WARINING: many part of this script were not explicitly tested, they are mostly for reference
# 17/02/2026

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]] ; then
    echo "This script downloads the databases used by flexprofiler and saves them to a specified directory."
    echo "It also generates a CSV file in the local directory with the paths to the downloaded databases, which can be used as input for flexprofiler."
    echo "The gqlib database is built by annotating the gene catalogue with eggNOG-mapper, which takes long (days with few cpus)."
    echo "Set SKIP_GQLIB=1 to skip it. The gqlib part of this script has not been tested."
    echo "Requires micromamba, wget and unzip."
    echo "Usage: $0 <out_database_path> [cpus (default: 8)]"
    exit 0
fi

mkdir -p "$1"
DBPATH=$(realpath "$1")
CPUS=${2:-8}
SKIP_GQLIB=${SKIP_GQLIB:-0}
TMPDIR=$(mktemp -d)
trap 'rm -rf "${TMPDIR}"' EXIT

# mOTUs 4.0.4
MOTUSDIR=${TMPDIR}/mOTUs4.0.4
micromamba create -y \
    -p ${MOTUSDIR} \
    -c bioconda \
    -c conda-forge \
    pip \
    'bwa>=0.7.19' \
    'vsearch>=2.30.4' \
    'python>=3.12'
micromamba run -p ${MOTUSDIR} pip install motus-tool==4.0.4
micromamba run -p ${MOTUSDIR} motus downloadMGDB
MOTUSDB_FROM="$(dirname $(micromamba run -p ${MOTUSDIR} python -c 'import motus.mutils as m; print(m.__file__)'))/db_mOTU"
MOTUS4_DB_TO=${DBPATH}/motus_v4.0.4
mkdir -p ${MOTUS4_DB_TO}
cp -r ${MOTUSDB_FROM} ${MOTUS4_DB_TO}

# mOTUs 3.1
MOTUSDIR=${TMPDIR}/mOTUs3.1
micromamba create -y \
    -p ${MOTUSDIR} \
    -c bioconda \
    -c conda-forge \
    motus=3.1
micromamba run -p ${MOTUSDIR} motus downloadDB
MOTUSDB_FROM=$(micromamba run -p ${MOTUSDIR} python -c 'import motus.motus as m; print(m.DATABASE)')
MOTUS3_DB_TO=${DBPATH}/motus_v3.1
mkdir -p ${MOTUS3_DB_TO}
cp -r ${MOTUSDB_FROM} ${MOTUS3_DB_TO}

# GMGC10 human gut gene catalogue, used by cayman and gqlib, with its bwa index and the cayman annotation
# (the downloaded archive contains all the cayman gene catalogues, ~19GB)
TOOLSDIR=${TMPDIR}/tools
micromamba create -y \
    -p ${TOOLSDIR} \
    -c bioconda \
    -c conda-forge \
    'bwa=0.7.19' \
    'gqlib=2.14.4'
GMGC_NAME=GMGC10.human-gut.95nr.0.5.percent.prevalence
GMGC_TO=${DBPATH}/gmgc10/${GMGC_NAME}
mkdir -p ${GMGC_TO} ${TMPDIR}/cayman
cd ${TMPDIR}/cayman
wget https://zenodo.org/records/10473258/files/gene_catalogue_annotations.zip
wget https://zenodo.org/records/10473258/files/gene_catalogues.zip
unzip gene_catalogue_annotations.zip
unzip gene_catalogues.zip
mv gene_catalogues/${GMGC_NAME}.fna.gz ${GMGC_TO}
mv annots/GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv ${GMGC_TO}
micromamba run -p ${TOOLSDIR} bwa index -b 90000000000 ${GMGC_TO}/${GMGC_NAME}.fna.gz
cd - > /dev/null
rm -rf ${TMPDIR}/cayman

# gqlib annotation database of the gene catalogue (NOT TESTED)
# The database used in the Zeller lab was built from the eggNOG-mapper 2.1.7 annotation (diamond mode) of
# GMGC10.human-gut.95nr.protein.0.5.percent.prevalence.faa, the proteins of the gene catalogue: here eggNOG-mapper
# translates the genes itself (--itype CDS), which gives the same proteins (not tested). The eggNOG database is version 5.0.2
# (inferred from the eggNOG-mapper version, 2.1.7, which downloads and expects it).
# See assets/test-data/README.md in the pipeline for the details of the database.
GQLIB_DB=${GMGC_TO}/GMGC10.human-gut.95nr.0.5pv.emapper_v2.1.7.sqlite3
if [[ ${SKIP_GQLIB} != 1 ]] ; then
    EMAPPERDIR=${TMPDIR}/emapper
    micromamba create -y \
        -p ${EMAPPERDIR} \
        -c bioconda \
        -c conda-forge \
        'eggnog-mapper=2.1.7'
    mkdir -p ${TMPDIR}/eggnog_db ${TMPDIR}/emapper_tmp
    micromamba run -p ${EMAPPERDIR} download_eggnog_data.py -y --data_dir ${TMPDIR}/eggnog_db
    gzip -cd ${GMGC_TO}/${GMGC_NAME}.fna.gz > ${TMPDIR}/${GMGC_NAME}.fna
    micromamba run -p ${EMAPPERDIR} emapper.py \
        -i ${TMPDIR}/${GMGC_NAME}.fna \
        --itype CDS \
        --data_dir ${TMPDIR}/eggnog_db \
        --output emapper_output \
        --output_dir ${TMPDIR} \
        -m diamond \
        --cpu ${CPUS} \
        --dbmem \
        --temp_dir ${TMPDIR}/emapper_tmp
    # --emapper_version v2.1.2 is the parser for the output of eggNOG-mapper >= 2.1.2
    micromamba run -p ${TOOLSDIR} python -m gqlib.bin.build_gene_database \
        ${GQLIB_DB} \
        ${TMPDIR}/emapper_output.emapper.annotations \
        --initialise_db \
        --emapper_version v2.1.2
fi

# the database sheet, with the paths to the downloaded databases
cat <<-EOT > flexprofiler_databases.csv
tool,db_name,db_params,db_path,annot_db_name,annot_db_path
motus4,db_mOTUv4,,${MOTUS4_DB_TO}/db_mOTU,,
motus3,db_mOTUv3,,${MOTUS3_DB_TO}/db_mOTU,,
cayman,${GMGC_NAME},,${GMGC_TO}/${GMGC_NAME}.fna.gz,cayman_v3,${GMGC_TO}/GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv
EOT
if [[ ${SKIP_GQLIB} != 1 ]] ; then
    echo "gqlib,${GMGC_NAME},,${GMGC_TO}/${GMGC_NAME}.fna.gz,emapper_v2.1.7,${GQLIB_DB}" >> flexprofiler_databases.csv
fi
