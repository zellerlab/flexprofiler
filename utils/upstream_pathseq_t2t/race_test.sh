#!/bin/bash
# Test the process-substitution race of upstream PathSeq-T2T (issue 6 of README.md): run the prefilter command of
# upstream (src/commands/prefilter.sh:150-155, aligner bwa) and measure how long after it returns its output BAM is
# complete; then the extraction and merge of the T2T-CHM13 filter with a decoy BED (src/commands/t2tfilter.sh:251-264,
# as patched), and check that the merged BAM has all the reads. Run it in the directory to test (e.g. a network file
# system), with samtools 1.22 on the PATH (environment.yml).
#
# usage: race_test.sh <grch38.bam> <decoys.bed> [runs] [threads]
set -uo pipefail
if [[ $# -lt 2 ]]; then
    sed -n '2,/^set -uo/p' "$0" | sed '$d; s/^# \{0,1\}//'
    exit 1
fi
BAM=$1 DECOYS=$2 RUNS=${3:-5} threads=${4:-4}

lags=""
for run in $(seq "$RUNS"); do
    rm -f race_unaligned.bam
    samtools view -@ "${threads}" -b "$BAM" \
      | samtools view -@ "${threads}" -bh - -f 3 -e '[AS]>35' \
          -U >(samtools view -@ "${threads}" -bh -F 2048 -x SA -x OQ -x MD -o race_unaligned.bam) -o /dev/null
    start=$(date +%s%N)
    until samtools quickcheck race_unaligned.bam 2>/dev/null; do sleep 0.01; done
    lags="$lags $(( ($(date +%s%N) - start) / 1000000 ))"
done
echo "prefilter: output BAM complete after (ms, upstream waits 1000):$lags"

bad=0
for run in $(seq "$RUNS"); do
    rm -f race_out.bam race_decoys.bam race_merged.bam
    samtools view -@ "${threads}" -bh "$BAM" -f 3 -e '[AS]>35' \
      -U >(samtools view -@ "${threads}" -bh -F 2048 -x SA -x OQ -x MD -o race_out.bam) -o /dev/null
    samtools view -@ "${threads}" -bh "$BAM" -f 3 -e '[AS]>35' -L "$DECOYS" \
      | samtools view -@ "${threads}" -bh -F 2048 -x SA -x OQ -x MD -o race_decoys.bam
    samtools cat -o race_merged.bam race_out.bam race_decoys.bam 2>/dev/null
    got=$(samtools view -c race_merged.bam 2>/dev/null || echo error)
    want=$(( $(samtools view -c race_out.bam) + $(samtools view -c race_decoys.bam) ))
    [[ "$got" == "$want" ]] || bad=$((bad + 1))
done
echo "t2tfilter: merged BAMs with missing reads: $bad of $RUNS"
rm -f race_unaligned.bam race_out.bam race_decoys.bam race_merged.bam
