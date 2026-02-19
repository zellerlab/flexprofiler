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
FQ1="${TMPDIR}/fastq1.txt"
FQ2="${TMPDIR}/fastq2.txt"
FASTA="${TMPDIR}/fasta.txt"
HEAD_F="${TMPDIR}/header.csv"
BODY_F="${TMPDIR}/body.csv"
echo $HEADER_STR > $HEAD_F

# Edit this as required to fit the filenames of your samples
find $INDIR -name "*_1.fastq.gz" > $FQ1
cat $FQ1|sed -e 's/_1.fastq.gz/_2.fastq.gz/' > $FQ2
cat $FQ1|sed -e 's%.*/%%' -e 's/_1.fastq.gz//' > $SAMPLE
cat $FQ1|sed 's/^.*$/ILLUMINA/' > $PLATFORM
cat $FQ1|sed 's/^.*$//' > $FASTA

# you should not need to edit these
paste -d ',' $SAMPLE $SAMPLE $PLATFORM $FQ1 $FQ2 $FASTA > $BODY_F
cat $HEAD_F $BODY_F > $SAMPLESHEET_OUT
rm -rf $TMPDIR