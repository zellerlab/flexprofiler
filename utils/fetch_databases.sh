#!/bin/bash

# fetch_databases.sh <out_database_path>
# Author: Saul Pierotti
# 17/02/2026

if [[ $# -ne 1 ]] ; then
    echo "This script downloads the databases used by flexprofiler and saves them to a specified directory."
    echo "It also generates a CSV file in the local directory with the paths to the downloaded databases, which can be used as input for flexprofiler."
    echo "Usage: $0 <out_database_path>"
    exit 0
fi

DBPATH=$1
TMPDIR=$(mktemp -d)

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
MOTUSDB_TO=${DBPATH}/motus_v4.0.4
mkdir -p ${MOTUSDB_TO}
cp -r ${MOTUSDB_FROM} ${MOTUSDB_TO}

# cayman
CAYMANDIR=${TMPDIR}/cayman
cd $CAYMANDIR
wget https://zenodo.org/records/10473258/files/gene_catalogue_annotations.zip
wget https://zenodo.org/records/10473258/files/gene_catalogues.zip
unzip gene_catalogue_annotations.zip
unzip gene_catalogues.zip
bwa index -b 90000000000 gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz
CAYMAN_TO=${DBPATH}/cayman
cp annots/GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv $CAYMAN_TO
cp gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz $CAYMAN_TO
cp gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz.* $CAYMAN_TO

# heredoc to generate the CSV file with the paths to the downloaded databases
cat <<-EOT > flexprofiler_databases.csv
tool,db_name,db_params,db_path
motus,db_mOTU,,${MOTUSDB_TO}/db_mOTU
cayman,annot,,${CAYMAN_TO}/annot/GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv
cayman,gene_cat,,${CAYMAN_TO}/gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz
cayman,bwa_amb,,${CAYMAN_TO}/gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz.amb
cayman,bwa_ann,,${CAYMAN_TO}/gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz.ann
cayman,bwa_bwt,,${CAYMAN_TO}/gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz.bwt
cayman,bwa_pac,,${CAYMAN_TO}/gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz.pac
cayman,bwa_sa,,${CAYMAN_TO}/gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz.sa
EOT

# cleanup
rm -rf ${TMPDIR}