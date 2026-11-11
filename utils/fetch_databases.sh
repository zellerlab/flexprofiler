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
    echo "All the databases are opt-in, set these variables to 1 to download them:"
    echo "  FETCH_MOTUS4=1: the mOTUs 4.0.4 database."
    echo "  FETCH_MOTUS3=1: the mOTUs 3.1 database."
    echo "  FETCH_CAYMAN=1: the GMGC10 human gut gene catalogue with its bwa index and the cayman annotation."
    echo "  FETCH_GQLIB=1: the gqlib database of the GMGC10 human gut gene catalogue (not tested, also downloads the catalogue)."
    echo "  FETCH_PATHSEQ_T2T=1: the references of the PathSeq-T2T host filtering (--pathseq_t2t, ~200 GB)."
    echo "    Set T2T_FASTA=<path> to link an existing GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz instead of downloading it (its md5 is checked)."
    echo "  FETCH_KRAKEN2=1: the Kraken2 PlusPF 2024-12-28 database (~75 GB)."
    echo "Requires micromamba, wget and unzip."
    echo "Usage: [FETCH_<DB>=1 ...] $0 <out_database_path> [cpus (default: 8)]"
    exit 0
fi

FETCH_MOTUS4=${FETCH_MOTUS4:-0}
FETCH_MOTUS3=${FETCH_MOTUS3:-0}
FETCH_CAYMAN=${FETCH_CAYMAN:-0}
FETCH_GQLIB=${FETCH_GQLIB:-0}
FETCH_PATHSEQ_T2T=${FETCH_PATHSEQ_T2T:-0}
FETCH_KRAKEN2=${FETCH_KRAKEN2:-0}
if [[ ${FETCH_MOTUS4} != 1 && ${FETCH_MOTUS3} != 1 && ${FETCH_CAYMAN} != 1 && ${FETCH_GQLIB} != 1 && ${FETCH_PATHSEQ_T2T} != 1 && ${FETCH_KRAKEN2} != 1 ]] ; then
    echo "No database selected: set at least one of FETCH_MOTUS4, FETCH_MOTUS3, FETCH_CAYMAN, FETCH_GQLIB, FETCH_PATHSEQ_T2T, FETCH_KRAKEN2 to 1 (run $0 without arguments for the help)." >&2
    exit 1
fi

mkdir -p "$1"
DBPATH=$(realpath "$1")
CPUS=${2:-8}
TMPDIR=$(mktemp -d)
trap 'rm -rf "${TMPDIR}"' EXIT

# bwa and gqlib, used for the gene catalogue, gqlib and PathSeq-T2T
TOOLSDIR=${TMPDIR}/tools
if [[ ${FETCH_CAYMAN} == 1 || ${FETCH_GQLIB} == 1 || ${FETCH_PATHSEQ_T2T} == 1 ]] ; then
    micromamba create -y \
        -p ${TOOLSDIR} \
        -c bioconda \
        -c conda-forge \
        'bwa=0.7.19' \
        'gqlib=2.14.4'
fi

# mOTUs 4.0.4, only with FETCH_MOTUS4=1
MOTUS4_DB_TO=${DBPATH}/motus_v4.0.4
if [[ ${FETCH_MOTUS4} == 1 ]] ; then
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
    mkdir -p ${MOTUS4_DB_TO}
    cp -r ${MOTUSDB_FROM} ${MOTUS4_DB_TO}
fi

# mOTUs 3.1, only with FETCH_MOTUS3=1
MOTUS3_DB_TO=${DBPATH}/motus_v3.1
if [[ ${FETCH_MOTUS3} == 1 ]] ; then
    MOTUSDIR=${TMPDIR}/mOTUs3.1
    micromamba create -y \
        -p ${MOTUSDIR} \
        -c bioconda \
        -c conda-forge \
        motus=3.1
    micromamba run -p ${MOTUSDIR} motus downloadDB
    MOTUSDB_FROM=$(micromamba run -p ${MOTUSDIR} python -c 'import motus.motus as m; print(m.DATABASE)')
    mkdir -p ${MOTUS3_DB_TO}
    cp -r ${MOTUSDB_FROM} ${MOTUS3_DB_TO}
fi

# GMGC10 human gut gene catalogue, used by cayman and gqlib, with its bwa index and the cayman annotation
# (the downloaded archive contains all the cayman gene catalogues, ~19GB)
# only with FETCH_CAYMAN=1 or FETCH_GQLIB=1 (the catalogue is needed by both)
GMGC_NAME=GMGC10.human-gut.95nr.0.5.percent.prevalence
GMGC_TO=${DBPATH}/gmgc10/${GMGC_NAME}
if [[ ${FETCH_CAYMAN} == 1 || ${FETCH_GQLIB} == 1 ]] ; then
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
fi

# gqlib annotation database of the gene catalogue (NOT TESTED), only with FETCH_GQLIB=1
# The database used in the Zeller lab was built from the eggNOG-mapper 2.1.7 annotation (diamond mode) of
# GMGC10.human-gut.95nr.protein.0.5.percent.prevalence.faa, the proteins of the gene catalogue: here eggNOG-mapper
# translates the genes itself (--itype CDS), which gives the same proteins (not tested). The eggNOG database is version 5.0.2
# (inferred from the eggNOG-mapper version, 2.1.7, which downloads and expects it).
# See assets/test-data/README.md in the pipeline for the details of the database.
GQLIB_DB=${GMGC_TO}/GMGC10.human-gut.95nr.0.5pv.emapper_v2.1.7.sqlite3
if [[ ${FETCH_GQLIB} == 1 ]] ; then
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

# PathSeq-T2T host filtering (--pathseq_t2t, see docs/pathseq_t2t.md), the references of Dohlman et al., Cell 2026
# (https://doi.org/10.1016/j.cell.2026.04.015), only with FETCH_PATHSEQ_T2T=1 (~200 GB of disk, the bwa indexes
# take hours). NOT TESTED as a whole: the URLs and checksums are those of the upstream scripts.
if [[ ${FETCH_PATHSEQ_T2T} == 1 ]] ; then
    PATHSEQ_TO=${DBPATH}/pathseq_t2t
    mkdir -p ${PATHSEQ_TO}/grch38_d1_vd1/bwa ${PATHSEQ_TO}/t2t_chm13v2/bwa ${PATHSEQ_TO}/pathseq_host
    # GRCh38 of the NCI Genomic Data Commons, with its viral decoy contigs (GRCh38.d1.vd1)
    wget -O ${TMPDIR}/GRCh38.d1.vd1.fa.tar.gz https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834
    tar -xzf ${TMPDIR}/GRCh38.d1.vd1.fa.tar.gz -C ${PATHSEQ_TO}/grch38_d1_vd1
    rm ${TMPDIR}/GRCh38.d1.vd1.fa.tar.gz
    micromamba run -p ${TOOLSDIR} bwa index -p ${PATHSEQ_TO}/grch38_d1_vd1/bwa/GRCh38.d1.vd1 ${PATHSEQ_TO}/grch38_d1_vd1/GRCh38.d1.vd1.fa
    # T2T-CHM13v2.0, linked from T2T_FASTA if set (e.g. the one of the host removal) instead of downloaded, checked either way
    if [[ -n ${T2T_FASTA:-} ]] ; then
        ln -sfn $(realpath ${T2T_FASTA}) ${PATHSEQ_TO}/t2t_chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
    else
        wget -O ${PATHSEQ_TO}/t2t_chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
    fi
    echo "9e6bf6b586bc8954208d1cc1d5f2fc99  ${PATHSEQ_TO}/t2t_chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" | md5sum -c -
    micromamba run -p ${TOOLSDIR} bwa index -p ${PATHSEQ_TO}/t2t_chm13v2/bwa/T2T-CHM13v2.0 ${PATHSEQ_TO}/t2t_chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
    # the host files of GATK PathSeq (pathseq_host.bfi, pathseq_host.fa.img)
    wget -O ${TMPDIR}/pathseq_host.tar.gz ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_host.tar.gz
    tar -xzf ${TMPDIR}/pathseq_host.tar.gz -C ${PATHSEQ_TO}/pathseq_host
    rm ${TMPDIR}/pathseq_host.tar.gz
    echo "PathSeq-T2T references in ${PATHSEQ_TO}: use them with"
    echo "  --pathseq_t2t_grch38_bwa_index ${PATHSEQ_TO}/grch38_d1_vd1/bwa --pathseq_t2t_chm13_bwa_index ${PATHSEQ_TO}/t2t_chm13v2/bwa --pathseq_host_dir <the directory of pathseq_host.bfi in ${PATHSEQ_TO}/pathseq_host>"
fi

# Kraken2 PlusPF 2024-12-28 (~75 GB), the latest release (Dohlman et al. used the one of 2024-06-05), only with FETCH_KRAKEN2=1
KRAKEN2_TO=${DBPATH}/kraken2/k2_pluspf_20241228
if [[ ${FETCH_KRAKEN2} == 1 ]] ; then
    mkdir -p ${KRAKEN2_TO}
    wget -O ${TMPDIR}/k2_pluspf_20241228.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20241228.tar.gz
    tar -xzf ${TMPDIR}/k2_pluspf_20241228.tar.gz -C ${KRAKEN2_TO}
    rm ${TMPDIR}/k2_pluspf_20241228.tar.gz
fi

# the database sheet, with the paths to the downloaded databases
echo "tool,db_name,db_params,db_path,annot_db_name,annot_db_path" > flexprofiler_databases.csv
if [[ ${FETCH_MOTUS4} == 1 ]] ; then
    echo "motus4,db_mOTUv4,,${MOTUS4_DB_TO}/db_mOTU,," >> flexprofiler_databases.csv
fi
if [[ ${FETCH_MOTUS3} == 1 ]] ; then
    echo "motus3,db_mOTUv3,,${MOTUS3_DB_TO}/db_mOTU,," >> flexprofiler_databases.csv
fi
if [[ ${FETCH_CAYMAN} == 1 ]] ; then
    echo "cayman,${GMGC_NAME},,${GMGC_TO}/${GMGC_NAME}.fna.gz,cayman_v3,${GMGC_TO}/GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv" >> flexprofiler_databases.csv
fi
if [[ ${FETCH_KRAKEN2} == 1 ]] ; then
    echo "kraken2,k2_pluspf_20241228,,${KRAKEN2_TO},," >> flexprofiler_databases.csv
fi
if [[ ${FETCH_GQLIB} == 1 ]] ; then
    echo "gqlib,${GMGC_NAME},,${GMGC_TO}/${GMGC_NAME}.fna.gz,emapper_v2.1.7,${GQLIB_DB}" >> flexprofiler_databases.csv
fi
