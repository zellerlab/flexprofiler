#!/bin/bash
# Samplesheet generator - make a copy of this file and edit it for your project
# Autor: Saul Pierotti

# Edit these as required
INDIR="<SET INPUT DIR HERE>" # where are the fastq files
SAMPLESHEET_OUT="samplesheet.csv" # where to save the samplesheet

# you should not need to edit these
HEADER_STR="sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta"
TMPDIR=$(mktemp -d)
PLATFORM="${TMPDIR}/platform.txt"
SAMPLE="${TMPDIR}/sample.txt"
RUN="${TMPDIR}/run.txt"
FQ1="${TMPDIR}/fastq1.txt"
FQ2="${TMPDIR}/fastq2.txt"
FASTA="${TMPDIR}/fasta.txt"
HEAD_F="${TMPDIR}/header.csv"
BODY_F="${TMPDIR}/body.csv"
echo $HEADER_STR > $HEAD_F

# Edit this as required to fit the filenames of your samples
# Pay attention that often files with different naming structure are present in the same experiment
# you can use the -e flag in sed to pass multiple regex patterns (see example below)
# to leave a column empty just use the following to generate a file with empty lines but of the right length:
# cat $FQ1|sed 's/^.*$//'
#
# FQ1: full path to the fastq1 file for pair end reads or path to the fastq file for single ended reads
# FQ2: full path to the fastq2 file for pair end reads, leave empty for single ended reads
# SAMPLE: unique ID for each sample (there can be multiple RUN with the same SAMPLE if the same actual biological sample was sequenced more than once and you want to collate the reads). This can also be set to the same value as RUN.
# RUN: unique ID for each sequencing file (or pair of files for PE). Normally this should be the name of the fastq file without extension and without 1/2 for PE reads
# PLATFORM: a literal string indicating the sequencing technology, from the ENA controlled vocabulary (see https://www.ebi.ac.uk/ena/portal/api/controlledVocab?field=instrument_platform)
# FASTA: normally leave empty, set to a literal "true" without quotes if the files are fasta and not fastq 
find $INDIR -name "*_1.fastq.gz" > $FQ1
cat $FQ1|sed -e 's/_1.fastq.gz/_2.fastq.gz/' > $FQ2
cat $FQ1|sed -e "s%.*/%%" -e "s/_1\\.fastq\\.gz//" -e "s/_R1_001\\.fastq\\.gz//" > $RUN
cat $RUN|sed -E  "s/^.*lane1(MG)?/MG/" > $SAMPLE
cat $FQ1|sed 's/^.*$/ILLUMINA/' > $PLATFORM
cat $FQ1|sed 's/^.*$//' > $FASTA

# you should not need to edit these
paste -d ',' $SAMPLE $RUN $PLATFORM $FQ1 $FQ2 $FASTA > $BODY_F
cat $HEAD_F $BODY_F > $SAMPLESHEET_OUT
rm -rf $TMPDIR