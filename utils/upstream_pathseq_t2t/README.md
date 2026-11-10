# Running upstream PathSeq-T2T, and its known issues

The PathSeq-T2T mode of flexprofiler (`--pathseq_t2t`, see [`docs/pathseq_t2t.md`](../../docs/pathseq_t2t.md))
reimplements the host filtering of upstream PathSeq-T2T:

> Dohlman A.B., Mjelle R., Wood H.M., et al. Biodiversity and biogeography of the multi-kingdom cancer microbiome.
> Cell 189, 3817-3837 (2026). doi: [10.1016/j.cell.2026.04.015](https://doi.org/10.1016/j.cell.2026.04.015)

Code: [abdohlman/pathseq-t2t](https://github.com/abdohlman/pathseq-t2t), v0.3.0, commit
`51d74430b6a4b34073f2d878612ac24b6a1d1e80` (the latest on 2026-11-10), MIT license.

This folder saves how upstream was run to compare it with flexprofiler, and the issues found in it.

## Recipe

| File               | What it is                                                                                                                                                                                        |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `environment.yml`  | The conda environment, with the tool versions of flexprofiler: GATK 4.6.2.0 on Java 17, Picard 3.4.0, samtools 1.22.1, bwa 0.7.19, Kraken2 2.1.6. Upstream's own `envs/main.yml` leaves them open. |
| `t2tfilter.patch`  | The two changes to `src/commands/t2tfilter.sh` without which upstream cannot run its T2T-CHM13 filter (issues 1 and 2 below). They do not change which reads are selected.                        |
| `run_upstream.sh`  | The commands: `pathseq-t2t filter`, `classify` (Kraken2) and `summarize` on one BAM aligned to GRCh38, with the options of the paper.                                                             |
| `race_test.sh`     | The test of the process-substitution race (issue 6), to run on the file system to test.                                                                                                          |

By hand, on Linux (x86-64):

```bash
git clone https://github.com/abdohlman/pathseq-t2t.git
git -C pathseq-t2t checkout 51d74430b6a4b34073f2d878612ac24b6a1d1e80
git -C pathseq-t2t apply /path/to/flexprofiler/utils/upstream_pathseq_t2t/t2tfilter.patch

micromamba create -p ./env -f /path/to/flexprofiler/utils/upstream_pathseq_t2t/environment.yml
export PATH=$PWD/env/bin:$PATH JAVA_HOME=$PWD/env/lib/jvm

bwa index t2t.fa   # once
/path/to/flexprofiler/utils/upstream_pathseq_t2t/run_upstream.sh \
    $PWD/pathseq-t2t $PWD/env/share/picard-3.4.0-0/picard.jar \
    sample.grch38.bam sample non_human_decoys.bed pathseq_host/ t2t.fa k2_pluspf_20240605/ out_sample/ 16 48
```

The input is a BAM of the reads aligned to GRCh38 with `bwa mem -T 0`, with a read group (flexprofiler writes it as
`pathseq_t2t/grch38/<sample>.grch38.bam` with `--save_pathseq_t2t_bams`). GATK PathSeqFilterSpark (Spark, Hadoop)
needs the user running it to be known to the system: in a container, run as a user that exists in it.

[`utils/compare_upstream_pathseq_t2t.sh`](../compare_upstream_pathseq_t2t.sh) does all of this with docker, on the
GRCh38 BAM of a flexprofiler run, and compares the reads after each step and the Kraken2 results. The results on the
test data are in [`docs/pathseq_t2t.md`](../../docs/pathseq_t2t.md#comparison-with-upstream).

## Issues of upstream PathSeq-T2T

At commit `51d74430`. "Reproduced" means seen when running upstream as above. flexprofiler is affected by none of
them.

### 1. The T2T-CHM13 filter always fails: `samtools view -o` without a file

Reproduced. `src/commands/t2tfilter.sh:251-253` (and `280-282` for the unpaired reads):

```bash
time samtools view -@ "${threads}" -bh "${bam_aligned_paired}" -f 3 -e '[AS]>35' \
  -U >(samtools view ... -o "${output_paired}") \
  -o >/dev/null #(samtools flagstat --output-fmt tsv - > "${flagstat_aligned_paired}")
```

Bash takes `>/dev/null` as a redirection of the output of the command, so `samtools view` gets `-o` without its file
name, and stops with `samtools view: invalid option -- 'o'`. `pathseq-t2t t2tfilter` (and `filter`) fails on any
input. Fix: `-o /dev/null` (`t2tfilter.patch`). The flagstat that the comment would write is not written either way.

### 2. The decoy reads cannot be merged with current samtools: `samtools cat -@`

Reproduced with samtools 1.22.1. `src/commands/t2tfilter.sh:263` (and `292`):

```bash
time samtools cat -@ "${threads}" -o "${tmp_merge_p}" "${output_paired}" "${bam_decoys_paired}"
```

samtools 1.22.1 `cat` has no `-@` option (`cat: invalid option -- '@'`); 1.16.1 and 1.20 have it. Upstream asks for
samtools >= 1.16 (`lib/requirements.sh`), so a current samtools passes its check and then fails here, whenever a decoy
BED is given (`--decoys-to-mask`): the decoy BAM is never empty, as it has a header (`[[ -s ... ]]`). Fix: drop `-@`
(`t2tfilter.patch`), or require samtools < 1.22.

### 3. qcfilter fails when no paired or no unpaired read is left

Reproduced with single-end reads. `src/commands/qcfilter.sh:304-314`: the outputs of PathSeqFilterSpark are required
to exist (`Expected qcfilter output missing .sbi: ...paired.bam.sbi`), but PathSeqFilterSpark writes no BAM for a
category with no reads (`No paired reads to write - BAM will not be written`). This happens for single-end input
(no paired reads), and for samples whose reads are all filtered out, e.g. a sample with no microbial reads, which is
expected in low-biomass data. The `FINAL_PAIRED_READS` and `FINAL_UNPAIRED_READS` of the filter metrics tell whether
a BAM is expected.

Single-end input has a second problem: the prefilter keeps the reads that are not in a proper pair with an alignment
score above 35 (`samtools view -f 3 -e '[AS]>35' -U`), so it keeps all the single-end reads, which are never in a
proper pair (PathSeqFilterSpark then removes those aligned to the host, `--is-host-aligned true`). flexprofiler selects
single-end reads by their alignment score alone, as upstream does for the singletons in `t2tfilter.sh`.

### 4. The error message of qcfilter is hidden by its cleanup

Reproduced. When qcfilter fails, its `EXIT` trap (`src/commands/qcfilter.sh:205-217`) runs after `qcfilter` has
returned, when its `local` variables no longer exist, and fails in turn under `set -u`:
`qcfilter.sh: line 210: PART_PAIRED_DIR: unbound variable`. This is the last line of the log, after the real error.

### 5. The singletons are counted twice in the Kraken2 results

Reproduced. `scripts/pst2t_summarize.py:248-254` multiplies the counts of the Kraken2 report of the unpaired reads by
2, as for the pairs (which Kraken2 counts as one fragment). An unpaired read is one read: the reads and reads per
million of each taxon are too high by the number of its unpaired reads. On the test data, 929 reads are left after the
filters and `FINAL_READ_COUNT_INCLUDING_MATES` is 1108. This is also discussed in `docs/read-counting-review.md` of
[fabispringer/pathseq-t2t-nextflow](https://github.com/fabispringer/pathseq-t2t-nextflow). flexprofiler counts each
read once.

### 6. BAMs written by background processes, then checked after a delay

The race is reproduced; no wrong result was. `src/commands/prefilter.sh:142-155` and
`src/commands/t2tfilter.sh:251-253, 280-281` write BAMs through process substitutions (`>(samtools view ... -o
out.bam)`), which bash does not wait for: the command returns while the BAM is still being written. `lib/io.sh` sleeps
one second before checking a BAM (`samtools quickcheck`). Reported as a race by the
[fabispringer/pathseq-t2t-nextflow](https://github.com/fabispringer/pathseq-t2t-nextflow) wrapper, which replaces
these scripts (`overrides/`).

Tested with the prefilter command of upstream on the GRCh38 BAM of the test sample repeated to 2.2 and 11 million
reads, on a local disk (docker on macOS, 4 threads):

- right after the command returned, the unaligned BAM was not complete in 19 of 20 runs
- it was complete 28-39 ms later, for both sizes (the lag is the end of the writing, not proportional to the size):
  well within the one-second wait of upstream
- in the T2T-CHM13 filter with a decoy BED, the output is merged (`samtools cat`) after the decoy extraction, which
  reads the whole alignment again: the 10 merges tested all had the expected reads

So on a local disk upstream waits long enough. If the end of the writing took more than a second (e.g. on a network
file system or on a loaded node; not tested), the check of `lib/io.sh` would fail and stop the run, rather than give
wrong results; the merge in the T2T-CHM13 filter has no check, and could read an incomplete BAM. The test is
`race_test.sh` (run it in a directory of the file system to test, e.g. on the cluster:
`race_test.sh sample.grch38.bam non_human_decoys.bed`).

### Not an issue: `-VALIDATION_STRINGENCY` with one dash

`src/commands/t2tfilter.sh:231, 236` pass `-VALIDATION_STRINGENCY SILENT` to Picard `SamToFastq`, where the other calls
use two dashes. Picard 3.4.0 accepts it (checked in its log).
