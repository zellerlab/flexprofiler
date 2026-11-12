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
    echo "  FETCH_PATHSEQ_T2T_STAR=1: the STAR index of GRCh38.d1.vd1 with the GENCODE v36 genes, for --pathseq_t2t_grch38_aligner star (~30 GB, ~1 h, needs ~40 GB of memory)."
    echo "    The GRCh38.d1.vd1 FASTA of FETCH_PATHSEQ_T2T is used if present, downloaded otherwise. Set STAR_RAM=<bytes> to change the memory given to STAR (default: 90% of the memory of the slurm job, or 40 GB)."
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
FETCH_PATHSEQ_T2T_STAR=${FETCH_PATHSEQ_T2T_STAR:-0}
FETCH_KRAKEN2=${FETCH_KRAKEN2:-0}
if [[ ${FETCH_MOTUS4} != 1 && ${FETCH_MOTUS3} != 1 && ${FETCH_CAYMAN} != 1 && ${FETCH_GQLIB} != 1 && ${FETCH_PATHSEQ_T2T} != 1 && ${FETCH_PATHSEQ_T2T_STAR} != 1 && ${FETCH_KRAKEN2} != 1 ]] ; then
    echo "No database selected: set at least one of FETCH_MOTUS4, FETCH_MOTUS3, FETCH_CAYMAN, FETCH_GQLIB, FETCH_PATHSEQ_T2T, FETCH_PATHSEQ_T2T_STAR, FETCH_KRAKEN2 to 1 (run $0 without arguments for the help)." >&2
    exit 1
fi

mkdir -p "$1"
DBPATH=$(realpath "$1")
CPUS=${2:-8}
TMPDIR=$(mktemp -d)
trap 'rm -rf "${TMPDIR}"' EXIT

# bash reads a script as it runs it: if this file is changed (or rewritten by git) while it runs, the rest is read from
# the changed file. A git pull during a run stopped it silently after the PathSeq-T2T block: edit a copy, or wait.

# bwa, for the indexes of the gene catalogue and of the PathSeq-T2T references, and gqlib, for its database: in
# separate environments, as bwa 0.7.19 and gqlib 2.14.4 cannot be installed together (zlib versions)
BWADIR=${TMPDIR}/bwa
if [[ ${FETCH_CAYMAN} == 1 || ${FETCH_GQLIB} == 1 || ${FETCH_PATHSEQ_T2T} == 1 ]] ; then
    micromamba create -y \
        -p ${BWADIR} \
        -c bioconda \
        -c conda-forge \
        'bwa=0.7.19'
fi
GQLIBDIR=${TMPDIR}/gqlib
if [[ ${FETCH_GQLIB} == 1 ]] ; then
    micromamba create -y \
        -p ${GQLIBDIR} \
        -c bioconda \
        -c conda-forge \
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
    micromamba run -p ${BWADIR} bwa index -b 90000000000 ${GMGC_TO}/${GMGC_NAME}.fna.gz
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
    micromamba run -p ${GQLIBDIR} python -m gqlib.bin.build_gene_database \
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
    micromamba run -p ${BWADIR} bwa index -p ${PATHSEQ_TO}/grch38_d1_vd1/bwa/GRCh38.d1.vd1 ${PATHSEQ_TO}/grch38_d1_vd1/GRCh38.d1.vd1.fa
    # T2T-CHM13v2.0, linked from T2T_FASTA if set (e.g. the one of the host removal) instead of downloaded, checked either way
    if [[ -n ${T2T_FASTA:-} ]] ; then
        ln -sfn $(realpath ${T2T_FASTA}) ${PATHSEQ_TO}/t2t_chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
    else
        wget -O ${PATHSEQ_TO}/t2t_chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
    fi
    echo "9e6bf6b586bc8954208d1cc1d5f2fc99  ${PATHSEQ_TO}/t2t_chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" | md5sum -c -
    micromamba run -p ${BWADIR} bwa index -p ${PATHSEQ_TO}/t2t_chm13v2/bwa/T2T-CHM13v2.0 ${PATHSEQ_TO}/t2t_chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
    # the host files of GATK PathSeq (pathseq_host.bfi, pathseq_host.fa.img)
    wget -O ${TMPDIR}/pathseq_host.tar.gz ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_host.tar.gz
    mkdir -p ${TMPDIR}/pathseq_host
    tar -xzf ${TMPDIR}/pathseq_host.tar.gz -C ${TMPDIR}/pathseq_host
    rm ${TMPDIR}/pathseq_host.tar.gz
    # the files are moved to pathseq_host/ itself, whatever the directory layout of the archive
    find ${TMPDIR}/pathseq_host -name 'pathseq_host.bfi' -exec mv {} ${PATHSEQ_TO}/pathseq_host/ \;
    find ${TMPDIR}/pathseq_host -name 'pathseq_host.fa.img' -exec mv {} ${PATHSEQ_TO}/pathseq_host/ \;
    rm -rf ${TMPDIR}/pathseq_host
    ls ${PATHSEQ_TO}/pathseq_host/pathseq_host.bfi ${PATHSEQ_TO}/pathseq_host/pathseq_host.fa.img > /dev/null
    echo "PathSeq-T2T references in ${PATHSEQ_TO}: use them with"
    echo "  --pathseq_t2t_grch38_bwa_index ${PATHSEQ_TO}/grch38_d1_vd1/bwa --pathseq_t2t_chm13_bwa_index ${PATHSEQ_TO}/t2t_chm13v2/bwa --pathseq_host_dir ${PATHSEQ_TO}/pathseq_host"
fi

# STAR index of GRCh38.d1.vd1 for PathSeq-T2T on RNA-seq (--pathseq_t2t_grch38_aligner star, see docs/pathseq_t2t.md),
# only with FETCH_PATHSEQ_T2T_STAR=1. Built as STAR_GENOMEGENERATE_GRCH38 of the pipeline would: the same STAR version
# (2.7.11b), the GTF (needed for --quantMode GeneCounts), and the defaults of --genomeSAindexNbases (14, what the
# pipeline computes for GRCh38: min(14, log2(genome length)/2 - 1)) and --sjdbOverhang (100). The GTF is GENCODE v36
# of the GDC reference files, the annotation of GRCh38.d1.vd1 used by the GDC
# (https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files).
PATHSEQ_STAR_TO=${DBPATH}/pathseq_t2t/grch38_d1_vd1/star
if [[ ${FETCH_PATHSEQ_T2T_STAR} == 1 ]] ; then
    GRCH38_DIR=${DBPATH}/pathseq_t2t/grch38_d1_vd1
    mkdir -p ${GRCH38_DIR} ${PATHSEQ_STAR_TO}
    # the FASTA of FETCH_PATHSEQ_T2T, downloaded here if missing
    if [[ ! -s ${GRCH38_DIR}/GRCh38.d1.vd1.fa ]] ; then
        wget -O ${TMPDIR}/GRCh38.d1.vd1.fa.tar.gz https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834
        echo "3ffbcfe2d05d43206f57f81ebb251dc9  ${TMPDIR}/GRCh38.d1.vd1.fa.tar.gz" | md5sum -c -
        tar -xzf ${TMPDIR}/GRCh38.d1.vd1.fa.tar.gz -C ${GRCH38_DIR}
        rm ${TMPDIR}/GRCh38.d1.vd1.fa.tar.gz
    fi
    # the GTF is kept next to the index (gzipped), STAR reads it uncompressed
    wget -O ${GRCH38_DIR}/gencode.v36.annotation.gtf.gz https://api.gdc.cancer.gov/data/be002a2c-3b27-43f3-9e0f-fd47db92a6b5
    echo "c03931958d4572148650d62eb6dec41a  ${GRCH38_DIR}/gencode.v36.annotation.gtf.gz" | md5sum -c -
    gzip -cd ${GRCH38_DIR}/gencode.v36.annotation.gtf.gz > ${TMPDIR}/gencode.v36.annotation.gtf
    STARDIR=${TMPDIR}/star
    micromamba create -y \
        -p ${STARDIR} \
        -c bioconda \
        -c conda-forge \
        'star=2.7.11b'
    # 90% of the memory of the slurm job (SLURM_MEM_PER_NODE is in MB) if not set, as the pipeline does
    if [[ -z ${STAR_RAM:-} ]] ; then
        if [[ -n ${SLURM_MEM_PER_NODE:-} ]] ; then
            STAR_RAM=$(( SLURM_MEM_PER_NODE * 1024 * 1024 * 9 / 10 ))
        else
            STAR_RAM=40000000000
        fi
    fi
    # the temporary files and the logs of STAR go to TMPDIR (--outTmpDir must not exist), the index to PATHSEQ_STAR_TO
    micromamba run -p ${STARDIR} STAR \
        --runMode genomeGenerate \
        --genomeDir ${PATHSEQ_STAR_TO} \
        --genomeFastaFiles ${GRCH38_DIR}/GRCh38.d1.vd1.fa \
        --sjdbGTFfile ${TMPDIR}/gencode.v36.annotation.gtf \
        --runThreadN ${CPUS} \
        --limitGenomeGenerateRAM ${STAR_RAM} \
        --outTmpDir ${TMPDIR}/star_tmp \
        --outFileNamePrefix ${TMPDIR}/star_
    cp ${TMPDIR}/star_Log.out ${PATHSEQ_STAR_TO}/genomeGenerate.Log.out
    rm ${TMPDIR}/gencode.v36.annotation.gtf
    # the pipeline checks this file: the index must have been built with a GTF
    ls ${PATHSEQ_STAR_TO}/SA ${PATHSEQ_STAR_TO}/exonGeTrInfo.tab > /dev/null
    echo "STAR index of GRCh38.d1.vd1 (GENCODE v36) in ${PATHSEQ_STAR_TO}: use it with"
    echo "  --pathseq_t2t_grch38_aligner star --pathseq_t2t_grch38_star_index ${PATHSEQ_STAR_TO}"
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
