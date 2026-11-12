# zellerlab/flexprofiler: PathSeq-T2T host filtering

`--pathseq_t2t` filters the host reads of low-biomass human samples (e.g. tumour or blood whole-genome sequencing)
with PathSeq-T2T, the host subtraction of:

> Dohlman A.B., Mjelle R., Wood H.M., et al. Biodiversity and biogeography of the multi-kingdom cancer microbiome.
> Cell 189, 3817-3837 (2026). doi: [10.1016/j.cell.2026.04.015](https://doi.org/10.1016/j.cell.2026.04.015)

It reimplements, as Nextflow modules, the filtering commands of the upstream code,
[abdohlman/pathseq-t2t](https://github.com/abdohlman/pathseq-t2t) v0.3.0 (commit
`51d74430b6a4b34073f2d878612ac24b6a1d1e80`, MIT license). The [differences](#differences-from-upstream-pathseq-t2t)
are listed below. Please cite the paper, and GATK PathSeq (Walker et al. 2018), when you use this mode.

## What it does

The reads go through the pipeline as they were sequenced, and are filtered in three steps, as in the paper:

1. **Alignment to GRCh38 and prefilter** (`BWA_GRCH38`, `PATHSEQT2T_PREFILTER`): the reads are aligned to GRCh38 with
   `bwa mem -T 0`, the parameters of the NCI Genomic Data Commons. The reads that are not in a proper pair with an
   alignment score above 35 (`samtools view -f 3 -e '[AS]>35'`) are kept. So are the proper pairs aligned to the viral
   decoy contigs of the reference (`--pathseq_decoys_bed`), which are not host reads. The primary reads of the
   alignment (`samtools flagstat`) are the raw reads of the read counts. With `--pathseq_t2t_grch38_aligner star`, the
   reads are aligned with STAR instead, which also counts the reads per host gene (see [RNA-seq](#rna-seq-star)).
2. **PathSeqFilterSpark** (`PATHSEQT2T_QCFILTER`): GATK PathSeqFilterSpark removes the reads still aligned to GRCh38,
   those of low quality (Phred < 15) or low complexity, those with k-mers of the host or aligned to it (the GRCh38,
   HLA, UniVec, Gencode and breakpoint sequences of the GATK PathSeq host bundle, `--pathseq_host_dir`), the
   duplicates, and the reads shorter than 60 bases after clipping
   (`--pathseq_min_clipped_read_length`). The reads aligned to the decoy contigs are reverted to unaligned reads first.
   It writes the reads left in pairs, and those left without their mate (singletons).
3. **T2T-CHM13 filter** (`PATHSEQT2T_T2TFILTER`): the pairs and the singletons are aligned separately to the complete
   human genome T2T-CHM13 with `bwa mem -T 0`, and the proper pairs and the singletons with an alignment score above
   35 are removed.

The pairs left go to the profilers; the singletons are kept apart (`analysis_ready_fastqs/<sample>.t2tfilt_unpaired.fastq.gz`),
and are classified by Kraken2 (`--run_kraken2`), on their own, as in the paper. mOTUs (`--run_motus4`, `--run_motus3`)
gets them too, as unpaired reads (`-s`) in the same run as the pairs, so that its profile includes them; mOTUs counts
each of them as one insert, like a pair. For single-end reads, all the reads are singletons, and go to all the profilers.
The standard preprocessing (without `--pathseq_t2t`) drops the reads left without their mate.

With `--run_kraken2` and `--run_profile_standardisation`, the Kraken2 reports are merged into tables of reads per taxon
(see [the output documentation](output.md#kraken2)). The reads per million of the paper are not computed: they are the
reads of `reads_clade` times 10^6 divided by the `raw_reads` of the read counts.

FastQC (or Falco) is run on the raw reads, and on the filtered reads, as the reads given to the profilers (`fastqc/processed/`, with the reads left without their mate as `*_singletons_*`). The reads left after
each step are in the [read counts](output.md#read-counts), as in Figure 1C of the paper.

The steps are run after [run merging](usage.md#run-merging): with `--perform_runmerging`, the runs of a sample are
filtered together, so that the duplicates are removed over the whole sample, as the paper does with one BAM per sample.
Without it, each run is filtered, deduplicated and profiled on its own.

## Usage

```bash
nextflow run zellerlab/flexprofiler \
    -profile docker \
    --input samplesheet.csv \
    --databases databases.csv \
    --outdir results \
    --pathseq_t2t \
    --pathseq_t2t_grch38_bwa_index /refs/grch38_d1_vd1/bwa \
    --pathseq_t2t_chm13_bwa_index /refs/t2t_chm13v2/bwa \
    --pathseq_host_dir /refs/pathseq_host \
    --perform_runmerging \
    --run_kraken2 \
    --run_profile_standardisation
```

with a `kraken2` row in the database sheet (see [the database documentation](databases.md#kraken2)).

Only Illumina FASTQ of `WGS` runs are accepted: other rows of the samplesheet are an error. The steps that change the
reads before the filtering are an error too, as PathSeqFilterSpark does the quality and complexity filtering itself:
`--perform_shortread_qc` and `--perform_shortread_hostremoval`. The repair of the
order of the read pairs (`--perform_repair_pe_order`) can be used. Cayman and gqlib (`--run_cayman`, `--run_gqlib`)
are not supported, and are an error too.

### RNA-seq (STAR)

With `--pathseq_t2t_grch38_aligner star` (`STAR_GRCH38`), the reads are aligned to GRCh38 with STAR, a splice-aware
aligner, in place of `bwa mem`, for RNA-seq. This is not part of the paper or of upstream, which takes BAMs of DRAGEN
or bwa. STAR writes an unsorted BAM with all the reads, the unaligned ones included (`--outSAMunmapped Within
KeepPairs`), one alignment per read (`--outSAMmultNmax 1`), and counts the reads per gene of its index
(`--quantMode GeneCounts`). As in the STAR route of
[pathseq-t2t-nextflow](#differences-from-the-star-route-of-pathseq-t2t-nextflow), the reads are aligned in two passes
(`--twopassMode Basic`): the splice junctions found in the first pass, also those not in the GTF, are added to the
index, and the reads are aligned again, so that the reads with only a few bases across a new junction are aligned
spliced. The prefilter then takes as host reads the proper pairs (single-end: the aligned reads), with no alignment
score threshold, as upstream does for DRAGEN (`--aligner dragen`): STAR reports only the
alignments that pass its own filters (by default, at least 66% of the bases of the pair matched), and its alignment
score is of the pair. The reads aligned to more than 10 loci (`--outFilterMultimapNmax`) are left unaligned by STAR,
and are left to the host k-mer filter of PathSeqFilterSpark. The rest of the filtering (PathSeqFilterSpark,
T2T-CHM13) is unchanged; note that PathSeqFilterSpark removes the duplicates, which in RNA-seq are also those of the
highly expressed transcripts.

The STAR index (`--pathseq_t2t_grch38_star_index`) must have been built with a GTF of the genes, whose reads are
counted. If it is not given, it is built from the FASTA (`--pathseq_t2t_grch38_reference`) and the GTF
(`--pathseq_t2t_grch38_gtf`; for GRCh38.d1.vd1, the GENCODE v36 annotation of the
[GDC reference files](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)), decompressed first
if gzipped (about 1 h, 48 GB), and saved with `--save_hostremoval_index`. The bwa index of GRCh38 is not needed then.
`STAR_GRCH38` asks for 12 CPUs, 40 GB and 24 h. The runs must still be `WGS` in the samplesheet.

#### Differences from the STAR route of pathseq-t2t-nextflow

[fabispringer/pathseq-t2t-nextflow](https://github.com/fabispringer/pathseq-t2t-nextflow) (commit `9ef48f45`), a
Nextflow wrapper of upstream, also aligns RNA-seq reads to GRCh38 with STAR (`--aligner star`, its default; process
`ALIGN_STAR` of its `main.nf`). The STAR step of flexprofiler differs from it as follows:

| pathseq-t2t-nextflow                                                                                                                                                                                                     | flexprofiler                                                                                                                                                                                                                      | Why, and can it change results?                                                                                                                                                                                                    |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| The reads not aligned to GRCh38 are the unmapped reads written by STAR (`--outReadsUnmapped Fastx`: the pairs not aligned as a pair, both mates), converted to a BAM with Picard `FastqToSam`, which adds the read group | STAR writes all the reads in its BAM, the unaligned ones included (`--outSAMunmapped Within KeepPairs`), with the read group (`--outSAMattrRGline`), and the prefilter keeps those not in a proper pair (single-end: not aligned) | Same reads: STAR flags as a proper pair every pair aligned with both mates, and with its default filters (66% of the bases of the pair matched) a pair is rarely aligned by one mate alone. One BAM for bwa and STAR, no Picard    |
| The reads given to PathSeqFilterSpark (`--is-host-aligned true`) are all unaligned (`FastqToSam`), so its pre-aligned host filter removes none; in a pair with one mate aligned to GRCh38 (in the unmapped reads of STAR with its mate), the aligned mate goes on to the host k-mer filter | The reads keep their alignment, so the aligned mate of such a pair is removed by the pre-aligned host filter, as in the bwa route of both | Almost always the same reads: aligned to GRCh38, the mate is removed by the host k-mer filter in pathseq-t2t-nextflow. It is counted in `prealigned_host_filter_reads` instead of `host_kmer_filter_reads`. Such pairs need one mate with 66% of the bases of the pair: rare without trimming, more frequent after the BBDuk trimming of pathseq-t2t-nextflow |
| No rescue of the reads aligned to the viral decoy contigs: aligned to GRCh38, they are not in the unmapped reads, and are removed as host reads (its `docs/star-decoy-handling.md`)                                      | The proper pairs (single-end: aligned reads) that overlap the decoy BED go to PathSeqFilterSpark, as with bwa                                                                                                                     | More viral reads kept (HPV, EBV, HBV, HCV, HIV, ...), as in the bwa route of both. Its rescue at the T2T-CHM13 step does not make up for it, as T2T-CHM13 has no decoy contigs                                                     |
| All the alignments of the multi-mapping reads (STAR default)                                                                                                                                                             | One alignment per read (`--outSAMmultNmax 1`)                                                                                                                                                                                     | Same reads: a read aligned to 2-10 loci is a host read either way. One record per read, as with bwa, so that a read is rescued only if its best alignment is on a decoy contig, and is not duplicated in the BAMs of the prefilter |
| BAM sorted by coordinate (`--outSAMtype BAM SortedByCoordinate`) and indexed                                                                                                                                             | Unsorted BAM (the default of the nf-core module), the mates together                                                                                                                                                              | Nothing reads the coordinate order; the T2T-CHM13 step groups the mates (`samtools collate`). Less memory and time. The BAM is not published                                                                                       |
| Transcript alignments (`--quantMode TranscriptomeSAM`, `Aligned.toTranscriptome.out.bam`), not used or published                                                                                                         | Not written                                                                                                                                                                                                                       | None: the gene counts (`GeneCounts`) are the same                                                                                                                                                                                  |
| Raw reads: the input reads of `Log.final.out` (pairs times two), by `scripts/write_star_primary_flagstat.sh`                                                                                                             | Raw reads: the primary reads of `samtools flagstat` of the STAR BAM, as with bwa                                                                                                                                                  | Same number: the BAM has each input read once                                                                                                                                                                                      |
| Reads trimmed by BBDuk before STAR; the pairs and the reads left without their mate are aligned by two STAR runs, whose BAMs and gene counts are merged                                                                  | Reads not trimmed (PathSeqFilterSpark does the quality filtering, as for bwa); one STAR run, of the pairs or of the single-end reads                                                                                              | The trimming is in both routes of pathseq-t2t-nextflow, and is not part of upstream PathSeq-T2T                                                                                                                                    |
| Gene counts collated into cohort matrices (unstranded, forward, reverse), and tables of the summary lines | Gene counts published per sample (`pathseq_t2t/grch38/`), and with `--pathseq_t2t_merge_gene_counts` merged into tables (unstranded, forward, reverse), with one summary table that also has the reads assigned to genes; the gene counts are also in MultiQC | Same counts |
| Read group `ID` and `SM`: the sample; `PL`: `params.read_group_platform` (`ILLUMINA`)                                                                                                                                    | `ID`: the sample, or `<sample>_<run>` without run merging; `SM`: the sample; `PL`: the `instrument_platform` of the samplesheet                                                                                                   | Same as the bwa route of flexprofiler                                                                                                                                                                                              |
| The STAR index must be given                                                                                                                                                                                             | Built from the FASTA and the GTF if not given                                                                                                                                                                                     | Convenience only                                                                                                                                                                                                                   |

Both use two-pass alignment (`--twopassMode Basic`), count the reads per gene (`--quantMode GeneCounts`), and keep
STAR's default alignment filters.

### References

| Parameter                                                         | Reference                                                                                                                                                                                                           |
| ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--pathseq_t2t_grch38_bwa_index` or `--pathseq_t2t_grch38_reference` | GRCh38 of the NCI Genomic Data Commons, with its decoy contigs (GRCh38.d1.vd1, [GRCh38.d1.vd1.fa.tar.gz](https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834)). The bwa index is built from the FASTA if not given (about 1 h, 24 GB). |
| `--pathseq_t2t_chm13_bwa_index` or `--pathseq_t2t_chm13_reference` | T2T-CHM13v2.0 ([GCF_009914755.1](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz)). |
| `--pathseq_host_dir`                                              | The directory of `pathseq_host.bfi` and `pathseq_host.fa.img`, from the GATK PathSeq resource bundle (`ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_host.tar.gz`, 12 GB).            |
| `--pathseq_decoys_bed`                                            | The viral decoy contigs of GRCh38.d1.vd1, by default `assets/pathseq_t2t/non_human_decoys.bed` (from upstream). Empty to remove the reads aligned to them as host reads.                                          |
| `--pathseq_t2t_grch38_star_index` or `--pathseq_t2t_grch38_reference` and `--pathseq_t2t_grch38_gtf` | With `--pathseq_t2t_grch38_aligner star`, in place of the bwa index of GRCh38: see [RNA-seq](#rna-seq-star). |

The indexes built by the pipeline are saved with `--save_hostremoval_index`.

### Resources

PathSeqFilterSpark loads the BWA index image of the host bundle in memory: `PATHSEQT2T_QCFILTER` asks for 16 CPUs and
64 GB, of which 48 GB are given to Java, as upstream. The alignment to GRCh38 of a 30x genome takes several hours.

## Output

- `pathseq_t2t/prefilter/<sample>.flagstat.tsv`: `samtools flagstat` of the alignment to GRCh38
- with `--pathseq_t2t_grch38_aligner star`, `pathseq_t2t/grch38/<sample>.grch38.ReadsPerGene.out.tab`: the reads per
  host gene (STAR `GeneCounts`: after four summary lines, one line per gene, with the counts of unstranded, forward and
  reverse stranded libraries; also in MultiQC), and `<sample>.grch38.Log.final.out`, the summary of the alignment (also
  in MultiQC)
- with `--pathseq_t2t_grch38_aligner star` and `--pathseq_t2t_merge_gene_counts`,
  `pathseq_t2t/grch38/grch38_gene_counts.{unstranded,forward,reverse}.tsv`:
  the reads per host gene of all the samples merged (`bin/star_genecounts_merge.py`), one table per strandedness of the
  library, with the genes (rows, in the order of the STAR index) by samples (columns); and
  `grch38_gene_counts.summary.tsv`, for each strandedness, the four summary counts of STAR (`N_unmapped`,
  `N_multimapping`, `N_noFeature`, `N_ambiguous`) and the reads assigned to a gene (`N_assigned`) of each sample: the
  strandedness of a library is the one with the most reads assigned
- `pathseq_t2t/qcfilter/<sample>.prefilter.{unaligned,decoys}.filter_metrics.txt`: the reads after each filter of
  PathSeqFilterSpark, for the reads not aligned to GRCh38 and for those aligned to its decoy contigs
- `pathseq_t2t/t2tfilter/<sample>.t2tfilter.log`: the reads before and after the T2T-CHM13 filter
- `analysis_ready_fastqs/<sample>.t2tfilt_{1,2,unpaired}.fastq.gz` (`--save_analysis_ready_fastqs`, or in
  `pathseq_t2t/t2tfilter/` with `--save_hostremoval_unmapped`): the filtered reads
- with `--save_pathseq_t2t_bams`, the BAMs of each step

Without run merging, `<sample>` is `<sample>_<run_accession>`.

The read counts have, for each sample, the columns `raw_reads` (the primary reads aligned to GRCh38),
`grch38_unaligned_reads`, `prealigned_host_filter_reads`, `quality_complexity_filter_reads`, `host_kmer_filter_reads`,
`deduplication_filter_reads` and `t2t_filter_reads`, each mate counted.

## Differences from upstream PathSeq-T2T

### Fixed bugs

These are in the upstream code at commit `51d74430`, still its latest when this mode was written; flexprofiler is
affected by none of them. They are described in detail, with how they were found, in
[`utils/upstream_pathseq_t2t/README.md`](../utils/upstream_pathseq_t2t/README.md#issues-of-upstream-pathseq-t2t).
"Reproduced": seen when running upstream (see [the comparison](#comparison-with-upstream)).

| Upstream                                                                                                                                                                                                                  | flexprofiler                                                                                                     | Effect in upstream                                                                                                       |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| The T2T-CHM13 filter passes `-o >/dev/null` to `samtools view` (`src/commands/t2tfilter.sh:253, 282`): bash takes `>/dev/null` as a redirection, and `samtools view` stops with `invalid option -- 'o'`. Reproduced.       | Each output is written to a file                                                                                 | The T2T-CHM13 filter fails on any input.                                                                                 |
| The decoy reads are merged back with `samtools cat -@` (`src/commands/t2tfilter.sh:263, 292`), an option that samtools 1.22.1 no longer has (1.16.1 and 1.20 have it), although upstream asks for samtools >= 1.16. Reproduced. | No decoy rescue at the T2T-CHM13 step (see the deliberate differences)                                                                                      | The T2T-CHM13 filter fails with current samtools whenever a decoy BED is given.                                         |
| The step stops if PathSeqFilterSpark writes no paired or no unpaired reads (`src/commands/qcfilter.sh:304-314`), which it does when none is left. Reproduced with single-end reads.                                        | A missing output is an error only if the filter metrics of its input (the unaligned or the decoy reads) count reads for it (`FINAL_PAIRED_READS`, `FINAL_UNPAIRED_READS`), and the next steps take empty inputs | Fails on single-end samples, and on samples with no paired or no unpaired reads left (e.g. no microbial reads).         |
| When qcfilter fails, its cleanup (`EXIT` trap, `src/commands/qcfilter.sh:205-217`) fails in turn on an unbound variable, `PART_PAIRED_DIR`. Reproduced.                                                                    | Not applicable (Nextflow keeps the work directory of a failed task)                                              | The real error is followed by a misleading one.                                                                          |
| The singletons are counted twice in the Kraken2 results, like the pairs (`scripts/pst2t_summarize.py:248-254`). Reproduced: 1108 reads counted for 929 reads left.                                                         | Each read is counted once: the pairs twice (two reads), the singletons once (`bin/kraken2_merge.py`)             | The reads and reads per million of each taxon are too high by its classified singletons.                                |
| The BAMs are written by background process substitutions, and checked after a `sleep 1` (`src/commands/prefilter.sh:142-155`, `src/commands/t2tfilter.sh:251-253, 280-281`, `lib/io.sh`). Reproduced: the BAM is usually incomplete when the command returns, and complete 30-40 ms later on a local disk. | The prefilter also reads its input once, through named pipes, but waits for each reader and checks its exit status; the T2T-CHM13 filter writes each BAM by a command that has finished before the next one starts | None seen on a local disk, where the one-second wait is enough. If writing took longer (network file system, loaded node; not tested), the run would stop, or the decoy merge of the T2T-CHM13 filter, which has no check, could read an incomplete BAM. |

Picard `SamToFastq` is called with the single-dash option `-VALIDATION_STRINGENCY` (`src/commands/t2tfilter.sh:231,
236`): Picard 3.4.0 accepts it, so it is not an issue. flexprofiler uses `samtools fastq` for both conversions.

### Deliberate differences

| Upstream / paper                                                                                                                                       | flexprofiler                                                                                                                                                                                         | Why, and can it change results?                                                                                                                                                                         |
| ------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Input: a BAM already aligned to GRCh38 (DRAGEN for the cohorts of the paper, `bwa mem -T 0` for its in silico benchmark)                               | Input: FASTQ, aligned by the pipeline with `bwa mem -T 0` to GRCh38                                                                                                                                  | flexprofiler takes FASTQ. As in the benchmark of the paper; DRAGEN alignments can differ slightly.                                                                                                      |
| Single-end reads are not supported: not documented either way, but the code assumes paired-end reads throughout. The prefilter keeps the reads not in a proper pair (`-f 3`), i.e. all of them, and the run then stops at qcfilter (see the fixed bugs) | Single-end reads are kept if they have no alignment score above 35 to GRCh38, as upstream does for the singletons at the T2T-CHM13 step                                                               | An addition, not a change: upstream does not remove any host read from single-end samples at the prefilter, and fails on them at the next step.                                                       |
| One BAM per sample                                                                                                                                     | Per sample with `--perform_runmerging`, else per run                                                                                                                                                 | Without run merging, the duplicates across runs are not removed.                                                                                                                                         |
| No pair repair                                                                                                                                         | `--perform_repair_pe_order` can be used before the alignment                                                                                                                                         | Reads without their mate in the input are dropped by it, and are not in the raw reads.                                                                                                          |
| The decoy reads are merged with the other reads after PathSeqFilterSpark (`samtools merge`), and the T2T-CHM13 filter writes BAMs                       | The BAMs of the unaligned and of the decoy reads are not merged: the T2T-CHM13 step converts each to FASTQ and concatenates them (pairs and singletons apart) before aligning them, and writes FASTQ | Same reads. A read of a proper pair left without its mate by the filter is dropped, as upstream does when converting the BAMs to FASTQ for the classification (`classify.sh`, Picard `SamToFastq`). |
| The BAMs are converted to FASTQ with Picard `SamToFastq` (`--VALIDATION_STRINGENCY LENIENT`; `src/commands/t2tfilter.sh:231, 236`, `classify.sh`) | `samtools fastq -n -F 0xB00`, after `samtools collate` for the pairs; the reads without their mate are written apart (`-s`) and dropped, and counted in the log of the T2T-CHM13 step | Picard is not in the bwa/samtools container, and samtools streams with several threads. Same reads: `-F 0xB00` leaves out the secondary, supplementary and QC-failed (`0x200`) reads, as `SamToFastq` does by default (`INCLUDE_NON_PF_READS=false`), and `LENIENT` makes `SamToFastq` drop the reads without their mate with a warning. |
| The T2T-CHM13 filter keeps the reads aligned to it that overlap the decoy contigs (`--decoys-to-mask`, the same BED as the prefilter; `src/commands/t2tfilter.sh:255-266, 284-295`) | No decoy rescue at the T2T-CHM13 step: the BED is used only by the prefilter | None: the BED names viral contigs of GRCh38.d1.vd1, which T2T-CHM13 does not have, so upstream rescues no read (and warns that the BED chromosomes are not in the reference, `src/commands/t2tfilter.sh:114-119`). |
| The `pathseq-t2t` command line tool, with its conda environment                                                                                        | Three modules, with the containers of bwa 0.7.19 / samtools 1.22.1 and GATK 4.6.2.0 (the paper: BWA-MEM, samtools 1.16.1, GATK 4)                                                                    | One container per tool, as in the rest of the pipeline. GATK 4.6.2.0 needs Java 17 for Spark.                                                                                                           |
| Classification with Kraken2 (PlusPF, `--confidence 0.15`), MetaPhlAn4 and Sylph                                                                          | Kraken2 (`--run_kraken2`), with the same options (`db_params` of kraken2 databases rejected), pairs and singletons classified separately as upstream; the pairs and the singletons also go to mOTUs, in the same run (the singletons as unpaired reads, `-s`); Cayman and gqlib cannot be used. MetaPhlAn4 and Sylph are not included | Kraken2 2.1.6 instead of 2.1.3 (the version of the nf-core module).                                                                                                                                      |
| Reads per million, per sample (`<sample>.kraken.txt`), with `reads_clade_per_million` and `reads_taxon_per_million`                                       | Merged tables of all the samples (`kraken2/kraken2_<db>.{reads_clade,reads_taxon,taxonomy}.tsv`), no reads per million (from `raw_reads`)                                                              | Same numbers except for the singleton counts (see the fixed bugs).                                                                                                                                       |
| Metagenome assembly and binning, and decontamination across cancer types (pan-cancer equiprevalence score)                                              | Not included                                                                                                                                                                                         | The decontamination is a cohort analysis, to be done on the merged tables; the paper takes 0.1 reads per million as the background.                                                                      |
| Fixed value: `--min-clipped-read-length 60`                                                                                       | The same default, as a parameter (`--pathseq_min_clipped_read_length`)                                                                                                | Same results with the defaults.                                                                                                                                                                          |
| FastQC on the raw reads                                                                                                                                | Also on the filtered reads                                                                                                                                                                           | Reporting only.                                                                                                                                                                                          |
| `summary.tsv` per sample                                                                                                                               | The read counts table (`read_counts/read_counts.tsv`): `raw_reads` is `PRIMARY_READS` of the prefilter flagstat, the PathSeqFilterSpark columns are the `READS_AFTER_*` of its metrics (summed over the unaligned and decoy reads), and `t2t_filter_reads` is the reads not aligned to T2T-CHM13 | Same numbers, in one table for all the samples.                                                                                                                                                          |

## Comparison with upstream

`utils/compare_upstream_pathseq_t2t.sh` runs upstream PathSeq-T2T (commit `51d74430`) on the GRCh38 alignment of a
flexprofiler run, with the same tool versions, as saved in [`utils/upstream_pathseq_t2t/`](../utils/upstream_pathseq_t2t/README.md)
(environment, patch and commands, which can also be used by hand), and compares the reads left after each step, read by read (name, mate
and sequence), and the Kraken2 results. Upstream needs four of its lines changed to run at all (the `samtools view -o`
and `samtools cat -@` bugs above, `t2tfilter.patch`); none changes which reads are selected. Run flexprofiler with
`--pathseq_t2t --save_pathseq_t2t_bams --save_analysis_ready_fastqs --run_kraken2 --run_profile_standardisation`, then:

```bash
utils/compare_upstream_pathseq_t2t.sh <outdir> <sample> <t2t.fa> <pathseq_host_dir> <decoys.bed> <kraken2_db> <workdir>
```

On the paired-end sample of the test data (`test_pathseq_t2t` profile, two runs merged), both give the same reads at
every step, and the same Kraken2 reports:

| Step                                   | Upstream reads | flexprofiler reads | Same reads? |
| -------------------------------------- | -------------- | ------------------ | ----------- |
| Not aligned to GRCh38 (prefilter)      | 1580           | 1580               | yes         |
| Aligned to the decoy contigs           | 200            | 200                | yes         |
| PathSeqFilterSpark, paired             | 1034           | 1034               | yes         |
| PathSeqFilterSpark, unpaired           | 278            | 278                | yes         |
| Not aligned to T2T-CHM13, pairs        | 750            | 750                | yes         |
| Not aligned to T2T-CHM13, singletons   | 179            | 179                | yes         |

The read counts of each filter (`summary.tsv` of upstream, `read_counts.tsv` here) are the same. The Kraken2 reports
of the pairs and of the singletons are identical; the merged counts differ only by the counting of the singletons:

| Taxon                         | Pairs | Singletons | Upstream reads (2 × pairs + 2 × singletons) | flexprofiler reads (2 × pairs + singletons) |
| ----------------------------- | ----- | ---------- | ------------------------------------------- | ------------------------------------------- |
| Testella microbica (bacteria) | 284   | 166        | 900                                         | 734                                         |
| Decoyvirus testi (virus)      | 87    | 13         | 200                                         | 187                                         |
| Unclassified                  | 4     | 0          | 8                                           | 8                                           |

Upstream counts 1108 reads (`FINAL_READ_COUNT_INCLUDING_MATES`) for the 929 reads left after the filters. The
single-end sample of the test data cannot be compared: upstream stops at PathSeqFilterSpark, which writes no paired
reads for it (the `qcfilter.sh:304-314` bug above). The comparison has not been run on real data yet.

### Known behaviour of PathSeqFilterSpark

On the simulated reads of the tests, PathSeqFilterSpark removes about 10% of the error-free microbial reads in its
quality and complexity filters, and leaves some reads without their mate. The paper reports a similar loss (79.8%
of the microbial reads retained in its in silico mixtures, most of the loss at the quality and complexity filter,
Figure S1C).
