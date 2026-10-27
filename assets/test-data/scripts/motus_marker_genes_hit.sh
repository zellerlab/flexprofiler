#!/bin/bash
# list the mOTUs marker genes with alignments of the reads of one read file that pass the motus filters,
# aligned as motus does (single end, bwa mem -a): alignment length (M+I+D) >= 75, identity ((length - NM) / length) >= 97%
set -euo pipefail
IDX=$1; OUT=$2; shift 2
bwa mem -t 8 -a -v 1 $IDX "$@" \
    | samtools view -F 4 \
    | awk -F'\t' '{
        len = 0; cigar = $6
        while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
            op = substr(cigar, RLENGTH, 1); n = substr(cigar, 1, RLENGTH - 1) + 0
            if (op == "M" || op == "I" || op == "D") len += n
            cigar = substr(cigar, RLENGTH + 1)
        }
        nm = 0; for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) nm = substr($i, 6) + 0
        if (len >= 75 && (len - nm) / len * 100 >= 97) print $3
    }' \
    | sort -u > $OUT
