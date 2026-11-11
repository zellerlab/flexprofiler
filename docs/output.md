# zellerlab/flexprofiler: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [UNTAR](#untar) - Optionally saved decompressed input databases
- [FastQC](#fastqc) - Raw read QC
- [falco](#fastqc) - Alternative to FastQC for raw read QC
- [fastp](#fastp) - Adapter trimming for Illumina data
- [Bowtie2](#bowtie2) - Host removal for Illumina reads
- [SAMtools stats](#samtools-stats) - Statistics from host removal
- [PathSeq-T2T](pathseq_t2t.md#output) - Host filtering of low-biomass samples (`--pathseq_t2t`)
- [Analysis Ready Reads](#analysis-ready-reads) - Optional results directory containing the final processed reads used as input for classification/profiling.
- [mOTUs](#motus) - Tool for marker gene-based OTU (mOTU) profiling (v4 and v3).
- [bwa](#bwa) - Optional alignments of the reads to the gene catalogues of cayman and gqlib
- [cayman](#cayman) - Profiling of carbohydrate-active enzymes from a gene catalogue
- [gqlib](#gqlib) - Functional profiling from the annotation of a gene catalogue
- [Kraken2](#kraken2) - Taxonomic classification of the reads, with merged tables of reads
- [Read counts](#read-counts) - Per-run read counts after QC, host removal and alignment to each profiling database
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### untar

untar is used in zellerlab/flexprofiler to decompress various input files ending in `.tar.gz`. This process is mainly used for decompressing input database archive files.

<details markdown="1">
<summary>Output files</summary>

- `untar/`
  - `database/`
    - `<database_file_name>`: directory containing contents of the decompressed archive

</details>

This directory will only be present if `--save_untarred_databases` is supplied. The contained directories can be useful for moving the decompressed directories to a central 'cache' location allowing users to re-use the same databases. This is useful to save unnecessary computational time of decompressing the archives on every run.

### FastQC or Falco

<details markdown="1">
<summary>Output files</summary>

- `{fastqc,falco}/`
  - `raw/`: the reads as given in the samplesheet
  - `processed/`: the reads given to the profilers, after all the preprocessing (and run merging, and PathSeq-T2T, with its `*_singletons_*` reports of the paired-end reads left without their mate)
    - `*html`: FastQC or Falco report containing quality metrics in HTML format.
    - `*.txt`: FastQC or Falco report containing quality metrics in TXT format.
    - `*.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images (FastQC only).

The reports are named `<sample>_<run_accession>[_singletons]_<stage>[_1|_2]`, with the stage (`raw` or `processed`) and the mate of the paired-end reads, one report per file, for both tools (e.g. `<prefix>_1_fastqc.zip` with FastQC, `<prefix>_1_fastqc_data.txt`, `<prefix>_1_summary.txt` and `<prefix>_1_report.html` with Falco). In MultiQC, the reports are split by stage into two sections, and the samples are named without the stage (`<sample>_<run_accession>[_singletons][_1|_2]`).

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

zellerlab/flexprofiler runs FastQC/Falco on the raw reads (unless `--skip_preprocessing_qc`) and on the reads given to the profilers, after all the preprocessing steps (unless `--skip_postprocessing_qc`), to allow evaluation of the preprocessing as a whole; the effect of each step is in the [read counts](#read-counts). The files without reads (e.g. a sample with all its reads removed) are not checked. Note in the General Stats table, the columns of these two instances of FastQC/Falco are placed next to each other to make it easier to evaluate. However, the columns of the actual preprocessing steps (e.g. fastp) will be displayed _after_ the two FastQC/Falco columns, even if they were run 'between' the two FastQC/Falco jobs in the pipeline itself.

:::info
Falco produces identical output to FastQC but in the `falco/` directory.
:::

### fastp

[fastp](https://github.com/OpenGene/fastp) is a FASTQ pre-processing tool for quality control, trimmming of adapters, quality filtering and other features.

It is used in zellerlab/flexprofiler for adapter trimming of short-reads.

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `<sample_id>.fastp.fastq.gz`: File with the trimmed unmerged fastq reads.
  - `<sample_id>.merged.fastq.gz`: File with the reads that were successfully merged.
  - `<sample_id>.*{log,html,json}`: Log files in different formats.

</details>

By default zellerlab/flexprofiler will only provide the `<sample_id>.fastp.fastq.gz` file. The file `<sample_id>.merged.fastq.gz` will be available in the output folder if you provide the argument ` --shortread_qc_mergepairs` (optionally retaining un-merged pairs when in combination with `--shortread_qc_includeunmerged`).

You will only find the `.fastq` files in the results directory if you provide ` --save_preprocessed_reads`. If this is selected, you may receive different combinations of `.fastq` files for each sample depending on the input types - e.g. whether you have merged or not, or if you're supplying both single- and paired-end reads. Alternatively, if you wish only to have the 'final' reads that go into classification/profiling (i.e., that may have additional processing), do not specify this flag but rather specify `--save_analysis_ready_fastqs`, in which case the reads will be in the folder `analysis_ready_reads`.

:::warning
The resulting `.fastq` files may _not_ always be the 'final' reads that go into taxprofiling, if you also run other steps such as host removal, run merging etc..
:::

### Bowtie2

[Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes.

It is used with zellerlab/flexprofiler to allow removal of 'host' (e.g. human) and/or other possible contaminant reads (e.g. Phi X) from short-read `.fastq` files prior to profiling.

The two mates of paired-end reads are aligned separately, as single-end reads, and a pair is kept only if neither mate aligns to the host, as recommended by [Forbes et al. 2025](https://doi.org/10.1016/j.crmeth.2025.101218). The unaligned reads of the two mates are then paired again with [BBTools repair.sh](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/repair-guide/), which drops the reads whose mate aligned to the host and keeps the order of the reads.

<details markdown="1">
<summary>Output files</summary>

- `bowtie2/`
  - `build/`
    - `*.bt2`: Bowtie2 indicies of reference genome, only if `--save_hostremoval_index` supplied.
  - `align/`
    - `<sample_id>.bam`: BAM file containing reads that aligned against the user-supplied reference genome as well as unmapped reads. For paired-end reads, one per mate (`<sample_id>_1.bam`, `<sample_id>_2.bam`).
    - `<sample_id>.bowtie2.log`: log file about the mapped reads. For paired-end reads, one per mate (`<sample_id>_1.bowtie2.log`, `<sample_id>_2.bowtie2.log`).
    - `<sample_id>.unmapped.fastq.gz`: single-end reads, the off-target reads from the mapping that are used in downstream steps.
    - `<sample_id>.unmapped_{1,2}_repaired.fastq.gz`: paired-end reads, the pairs with neither mate aligned to the host, that are used in downstream steps.
    - `<sample_id>.unmapped.repair.sh.log`: paired-end reads, log of repair.sh, with the number of pairs kept and of reads dropped because their mate aligned to the host (`Singletons`).

</details>

By default zellerlab/flexprofiler will only provide the `.log` file if host removal is turned on. You will only have a `.bam` file if you specify `--save_hostremoval_bam`. This will contain _both_ mapped and unmapped reads. You will only get FASTQ files if you specify to save `--save_hostremoval_unmapped` - these contain only unmapped reads. Alternatively, if you wish only to have the 'final' reads that go into classification/profiling (i.e., that may have additional processing), do not specify this flag but rather specify `--save_analysis_ready_fastqs`, in which case the reads will be in the folder `analysis_ready_reads`.

:::info
Unmapped reads in FASTQ are only found in this directory for short-reads, for long-reads see [`samtools/fastq/`](#samtools-fastq).
:::

:::info
The resulting `.fastq` files may _not_ always be the 'final' reads that go into taxprofiling, if you also run other steps such as run merging etc..
:::

:::info
While there is a dedicated section in the MultiQC HTML for Bowtie2, these values are not displayed by default in the General Stats table. Rather, alignment statistics to host genome is reported via samtools stats module in MultiQC report (see below).
:::

### Amplicon primers

<details markdown="1">
<summary>Output files</summary>

- `amplicon_primers/`
  - `primers.tsv`: the primers found by PIMENTO, one row per run and mate, with the region of the primer library and the proportion of the reads with the primer (not written if the primers are given with `--amplicon_primer_fwd` and `--amplicon_primer_rev`).
  - `pimento/`: the outputs of PIMENTO (`<sample>_<run>_<mate>_std_primers.fasta`, `<sample>_<run>_<mate>_std_primer_out.txt`).
  - `vsearch_orient/<sample>_<run>.vsearch_orient.log`: the orientation of the long reads (vsearch), with the reads oriented and not oriented.
  - `cutadapt/`: the cutadapt logs, of the primers (`*_primers.cutadapt.log`) and read-through (`*_readthrough.cutadapt.log`) trimming of the short reads, and of the linked primer trimming of the long reads (`*_linked.cutadapt.log`).

</details>

### PacBio QC

<details markdown="1">
<summary>Output files</summary>

- `vsearch_fastqfilter/`
  - `<sample>_<run>.vsearch_fastqfilter.log`: the PacBio reads kept and discarded by the expected-errors and N filter (`--perform_longread_qc`).
  - `<sample>_<run>.filtered.fastq.gz`: the filtered reads, with `--save_preprocessed_reads`.
- `bbduk_pacbio/`
  - `<sample>_<run>_bbduk.bbduk.log`: the PacBio reads kept and discarded by the minimum base quality filter, after vsearch.
  - `<sample>_<run>_bbduk.fastq.gz`: the filtered reads, with `--save_preprocessed_reads`.

</details>

### Analysis Ready Reads

:::info
This optional results directory will only be present in the pipeline results when supplying `--save_analysis_ready_fastqs`.
:::

<details markdown="1">
<summary>Output files</summary>

- `analysis_ready_fastqs/`
  - `<sample_id>_{fq,fastq}.gz`: Final reads that underwent preprocessing and were sent for classification/profiling.
  - `motus4/<sample_id>_read{1,2}_renamed_reads.fastq.gz` (without `--perform_runmerging`, `<sample_id>_<run_accession>_read{1,2}_...`): With `--motus4_rename_reads`, the paired-end reads given to mOTUs v4, with the `/1` and `/2` suffixes removed from the read names (the other profilers get the reads above).

</details>

The results directory will contain the 'final' processed reads used as input for classification/profiling. It will _only_ include the output of the _last_ step of any combinations of preprocessing steps that may have been specified in the run configuration. For example, if you perform the read QC and host-removal preprocessing steps, the final reads that are sent to classification/profiling are the host-removed FASTQ files - those will be the ones present in this directory.

The reads of the `AMPLICON` runs are those at the end of the amplicon preprocessing (primer trimming and pair handling), named `<sample_id>_<run_accession>.fastq.gz`, or `<sample_id>_<run_accession>_{1,2}.fastq.gz` for the pairs kept with `--amplicon_pair_handling none`; with `--perform_runmerging`, those of the samples with more than one run are the merged reads.

:::warning
If you turn off all preprocessing steps, then no results will be present in this directory. This happens independently for short- and long-reads. I.e. you will only have FASTQ files for short reads in this directory if you skip all long-read preprocessing.
:::

### SAMtools stats

[SAMtools stats](http://www.htslib.org/doc/samtools-stats.html) collects statistics from a `.sam`, `.bam`, or `.cram` alignment file and outputs in a text format.

<details markdown="1">
<summary>Output files</summary>

- `samtools/stats/`
  - `<sample_id>.stats`: File containing samtools stats output.

</details>

In most cases you do not need to check this file, as it is rendered in the MultiQC run report.

### Run Merging

zellerlab/flexprofiler offers the option to merge FASTQ files of multiple sequencing runs or libraries that derive from the same sample, as specified in the input samplesheet.

This is the last possible preprocessing step, so if you have multiple runs or libraries (and run merging turned on), this will represent the final reads that will go into classification/profiling steps.

<details markdown="1">
<summary>Output files</summary>

- `run_merging/`
  - `*.fastq.gz`: Concatenated FASTQ files on a per-sample basis, named `<sample>_merged_<library_strategy>` (e.g. `2612_merged_WGS`)

</details>

Note that you will only find samples that went through the run merging step in this directory. For samples that had a single run or library will not go through this step of the pipeline and thus will not be present in this directory.

This directory and its FASTQ files will only be present if you supply `--save_runmerged_reads`. Alternatively, if you wish only to have the 'final' reads that go into classification/profiling (i.e., that may have additional processing), do not specify this flag but rather specify `--save_analysis_ready_fastqs`, in which case the reads will be in the folder `analysis_ready_reads`.

### mOTUs

[mOTUs](https://github.com/motu-tool/mOTUs) is a taxonomic profiler that maps reads to a unique marker specific database and estimates the relative abundance of known and unknown species. mOTUs v4 (`--run_motus4`) and v3 (`--run_motus3`) have the same outputs, in `motus4/` and `motus3/`.

<details markdown="1">
<summary>Output files</summary>

- `motus4/` (or `motus3/`)
  - `profile/<db_name>/`
    - `<sample_id>_<db_name>.log`: A log file that contains summary statistics
    - `<sample_id>_<db_name>.out`: A classification file that summarises taxonomic identifiers and their read counts in the profiled sample
  - `merge/`
    - `motus4_<db_name>_combined_reports.txt` (or `motus3_...`): A combined profile of all samples aligned to a given database

</details>

Normally `*_combined_reports.txt` is the most useful file for downstream analyses, but the per sample `.out` file can provide additional more specific information. Without `--perform_runmerging`, the files are per run, `<sample_id>_<run_accession>_<db_name>.*`.

### bwa

The reads are aligned with [bwa mem](https://github.com/lh3/bwa) to each gene catalogue used by cayman and gqlib, once for all the cayman and gqlib databases with the same gene catalogue.

<details markdown="1">
<summary>Output files</summary>

- `bwa/<db_name>/`
  - `<sample_id>_<db_name>.bam`: All the alignments of each read (`bwa mem -a`), not sorted

</details>

This directory is only present if you supply `--save_bwa_alignments`, as the alignments can be large.

### cayman

[cayman](https://github.com/zellerlab/cayman) profiles the carbohydrate-active enzymes (CAZymes) from the reads aligned to a gene catalogue annotated with CAZy domains. The pipeline runs the profiler used by cayman, gqlib, directly on the bwa alignment, with the same outputs as `cayman profile`.

<details markdown="1">
<summary>Output files</summary>

- `cayman/`
  - `profile/<db_name>/<annot_db_name>/`
    - `<sample_id>_<db_name>_<annot_db_name>.cazy.txt.gz`: The counts of each CAZy family, raw and in RPKM
    - `<sample_id>_<db_name>_<annot_db_name>.gene_counts.txt.gz`: The counts of each gene of the catalogue
    - `<sample_id>_<db_name>_<annot_db_name>.aln_stats.txt.gz`: Alignment statistics
    - `<sample_id>_<db_name>_<annot_db_name>.cayman.log`: The log, with the read counts
  - `merge/<db_name>/<annot_db_name>/`
    - `<db_name>_<annot_db_name>.cazy.uniq_raw.txt.gz`, `<db_name>_<annot_db_name>.cazy.uniq_rpkm.txt.gz`: The raw and RPKM counts of each CAZy family (rows) in all the samples (columns)

</details>

The `total_reads` and `filtered_reads` rows of the tables are the reads aligned to the gene catalogue and those passing the alignment filters (identity >= 97% and length >= 45). Without `--perform_runmerging`, the files are per run.

### gqlib

[gqlib](https://github.com/cschu/gqlib) counts the reads aligned to each gene of a gene catalogue, and summarises the counts for each functional category of the annotation database of the catalogue (e.g. eggNOG orthologous groups, KEGG orthologs, pathways and modules, CAZy, COG categories, GO terms, PFAM domains).

<details markdown="1">
<summary>Output files</summary>

- `gqlib/`
  - `profile_genes/<db_name>/<annot_db_name>/`
    - `<sample_id>_<db_name>_<annot_db_name>.<category>.txt.gz`: The counts of each feature of a functional category (e.g. `KEGG_ko`), raw, length-normalised (`lnorm`), scaled to the raw counts (`scaled`) and in RPKM
    - `<sample_id>_<db_name>_<annot_db_name>.gene_counts.txt.gz`: The counts of each gene of the catalogue
    - `<sample_id>_<db_name>_<annot_db_name>.aln_stats.txt.gz`: Alignment statistics
    - `<sample_id>_<db_name>_<annot_db_name>.gqlib.log`: The log, with the read counts
  - `merge/<db_name>/<annot_db_name>/`
    - `<db_name>_<annot_db_name>.<category>.<uniq|combined>_<metric>.txt.gz`: The counts of each feature of a category (rows) in all the samples (columns), for each of the metrics `raw`, `lnorm`, `scaled` and `rpkm`, including the gene counts (`gene_counts`)

</details>

The `total_reads` and `filtered_reads` rows of the category tables are the reads aligned to the gene catalogue and those passing the alignment filters (identity >= 97% and length >= 45), and the `category` row the reads assigned to any feature of the category. The gene count tables do not have these rows. Without `--perform_runmerging`, the files are per run.

#### What gqlib does

gqlib is the library of [gffquant](https://github.com/cschu/gff_quantifier), a functional profiler: from the alignments of the reads to a set of reference sequences, it counts the reads aligned to each annotated region of the references, and passes the counts on to the functional features the region is annotated with (e.g. the KEGG orthologs, pathways and GO terms of a gene). The result is a count table per functional category (e.g. `KEGG_ko`), with the abundance of each feature of the category in the sample.

gffquant has three modes, depending on the references and on how they are annotated:

- **`genes`** (used by gqlib in this pipeline): the references are the genes of a gene catalogue (e.g. GMGC), and the annotation gives the functional features of each gene. As every aligned read falls in an annotated gene, there is no need to look for overlaps between alignments and annotated regions. In gffquant the annotation is a gzipped table following the GMGC convention (the gene id in the first column, the functional categories from the seventh); gqlib reads it from a sqlite3 database built from that table (see the [database documentation](databases.md#gqlib)).
- **`domain`** (used by cayman): the references are again the genes of a gene catalogue, but the annotated regions are domains within the genes (the CAZy domains, with their coordinates in the gene, in the cayman `.csv` annotation), so only the alignments overlapping a domain count for it.
- **`genome`** (not used): the references are genome contigs, annotated with the coordinates of their genes in a gff3 file, and the alignments are assigned to the genes they overlap.

gffquant also aligns the reads itself (with bwa or minimap2), or takes an existing alignment. In this pipeline the reads are aligned separately with `bwa mem`, with the same settings, and gqlib and cayman take the alignment as input.

#### How gqlib and cayman count the reads

gqlib counts the reads as gffquant, with the same settings as the gffquant profiler of [nevermore](https://github.com/cschu/nevermore) (`--mode genes --ambig_mode 1overN`, and the default alignment filters). cayman uses the same counting, on the CAZy domains of the genes instead of the whole genes.

- **Alignment filters**: an alignment is counted only if its identity (matches over aligned length) is at least 97% and its aligned length at least 45 bp. `aln_stats.txt.gz` reports the alignments (not the reads) in total (`Total`), passing the filters (`Passed`), and removed by the identity (`Seqid`) and length (`Length`) filters.
- **Reads and pairs**: all the alignments of a read are processed together, so the alignment must be grouped by read name (the bwa output is, and it is not sorted by position). The two mates of a pair have the same name, and count as one read: each mate contributes half of the count of the pair.
- **Ambiguous reads**: a read with a single alignment contributes 1 to its gene. A read with `N` alignments passing the filters (`bwa mem -a` reports all the alignments of a read, not only the best one) contributes `1/N` to the gene of each of them (`1overN`, as the `dist1` mode of NGLess), counted separately for each mate of a pair. `N` includes the alignments to genes without annotation in a category, whose share is not counted in that category.
- **From genes to features**: the counts of a gene are added to each feature it is annotated with in a category (e.g. all its KEGG orthologs), so a gene with several features counts once for each of them, and the counts of the features of a category can sum to more than the counts of the genes.

The count tables have one row per feature, and one column per metric:

| Column    | Value                                                                                                                                                                         |
| --------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `*_raw`    | The counts of the reads                                                                                                                                                        |
| `*_lnorm`  | The counts normalised by the length of each gene (`count / length`), summed over the genes of the feature                                                                      |
| `*_scaled` | `lnorm * sum(raw) / sum(lnorm)`, with the sums over the genes annotated in the category (over all the genes for `gene_counts`): counts corrected for gene length, on the scale of the raw counts |
| `*_rpkm`   | `lnorm * 10^9 / filtered_reads`: reads per kilobase of gene per million reads passing the filters                                                                             |

The `uniq_` columns count only the reads with a single alignment. The `combined_` columns count also the ambiguous reads, with their `1/N` contributions, and are present only if the sample has ambiguous reads (otherwise they would be the same as the `uniq_` columns). cayman reports only the `raw` and `rpkm` metrics.

For example, take a gene catalogue with three genes, and their KEGG orthologs:

| Gene | Length   | `KEGG_ko`     |
| ---- | -------- | ------------- |
| A    | 1,000 bp | `K1` and `K2` |
| B    | 2,000 bp | `K1`          |
| C    | 500 bp   | none          |

and 10 single-end reads passing the alignment filters (`filtered_reads` = 10): 6 aligned only to A, 2 only to B, and 2 aligned both to A and to C.

1. **Genes** (`gene_counts`): the 8 reads with a single alignment count 1 each, so `uniq_raw` is 6 for A and 2 for B. The 2 ambiguous reads count `1/2` for A and `1/2` for C each, so `combined_raw` is 6 + 1 = 7 for A, 2 for B and 1 for C.
2. **Features** (`KEGG_ko`): the counts of each gene go to all its orthologs. `uniq_raw` is 6 + 2 = 8 for `K1` (from A and B) and 6 for `K2` (from A), and `combined_raw` 7 + 2 = 9 for `K1` and 7 for `K2`. The share of the ambiguous reads going to C, which has no ortholog, is not counted in `KEGG_ko`, and A is counted for both its orthologs, so `K1` and `K2` sum to 14 `uniq_raw` counts from 8 reads.
3. **Length normalisation**: `uniq_lnorm` is 6 / 1,000 = 0.006 for A and 2 / 2,000 = 0.001 for B, so 0.007 for `K1` and 0.006 for `K2`. B has a third of the reads of A, but is twice as long, so its share of `K1` is smaller in `lnorm` than in `raw`.
4. **Scaling**: the genes with a KEGG ortholog (A and B) sum to 8 in `uniq_raw` and 0.007 in `uniq_lnorm`, so `uniq_scaled` is `lnorm * 8 / 0.007`: 8 for `K1` and 6.86 for `K2` (compared to 8 and 6 raw). In `gene_counts`, the sums are over all the genes, and `uniq_scaled` is 6.86 for A and 1.14 for B: they sum to the 8 reads, but split according to the length-normalised counts.
5. **RPKM**: `uniq_rpkm` is `lnorm * 10^9 / 10`: 700,000 for `K1` and 600,000 for `K2`. These are comparable across samples with different sequencing depths, while `raw` and `scaled` are not.

With paired-end reads, the same holds for each pair: a pair with both mates aligned only to A counts `1/2` + `1/2` = 1 for A.

### Kraken2

[Kraken2](https://ccb.jhu.edu/software/kraken2/) assigns each read (or read pair) to the lowest common ancestor of the genomes of the database that share its k-mers (`--run_kraken2`).

<details markdown="1">
<summary>Output files</summary>

- `kraken2/`
  - `<db_name>/<sample>_<db_name>[_singletons].kraken2.report.txt`: the Kraken2 report of a sample (`_singletons`: of the reads left without their mate by `--pathseq_t2t`, classified on their own). Kraken2 counts a pair as one fragment.
  - `kraken2_<db_name>.reads_clade.tsv` (`--run_profile_standardisation`): the reads assigned to each taxon or to a taxon below it, by sample, each mate counted (the pairs counted twice, the singletons once), with the unclassified reads (taxid 0) first
  - `kraken2_<db_name>.reads_taxon.tsv`: the reads assigned to each taxon itself
  - `kraken2_<db_name>.taxonomy.tsv`: the name, rank and lineage (taxids from the root) of each taxon of the tables

</details>

`<sample>` is the sample, or `<sample>_<run_accession>` without run merging. The reports are also in the MultiQC report. The reads per million of PathSeq-T2T are not computed: with `--pathseq_t2t` they are `reads_clade` times 10^6 divided by the `raw_reads` of the [read counts](#read-counts), except that here the singletons are counted once instead of twice (see [the PathSeq-T2T page](pathseq_t2t.md#differences-from-upstream-pathseq-t2t)); the paper takes 0.1 reads per million as the background level of contamination.

### Read counts

<details markdown="1">
<summary>Output files</summary>

- `read_counts/`
  - `read_counts.tsv`: one row per run, with the reads left after each step.

</details>

The counts are parsed from the logs of the tools. There is one row per sample and run (columns `sample` and `run_accession`): with `--perform_runmerging`, the runs of a sample with more than one run are summed into its merged run (`merged_WGS` or `merged_AMPLICON`), and a sample with a single run keeps its `run_accession`. The table is also in the MultiQC report, in the section Read counts. A column is present only if its step was run: the raw and QC reads only with `--shortread_qc_tool fastp`, the pairs retained by host removal only with `--perform_shortread_hostremoval`, and one column per profiler and database.

| Column                      | Unit                                     | Source                                                                                                     |
| --------------------------- | ---------------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `raw_reads`                 | reads, each mate counted                 | fastp, reads before filtering                                                                              |
| `qc_reads`                  | reads, each mate counted                 | fastp, reads after filtering                                                                               |
| `hostremoval_retained_pairs` | read pairs (single-end: reads)           | repair.sh, pairs with neither mate aligned to the host genome (single-end: bowtie2, reads not aligned), i.e. the reads given to the profilers |
| `grch38_unaligned_reads`    | reads, each mate counted                 | PathSeq-T2T: reads not aligned to GRCh38 and reads of its viral decoy contigs, given to PathSeqFilterSpark (`PRIMARY_READS` of its metrics); with `--pathseq_t2t`, `raw_reads` is the primary reads aligned to GRCh38 (`samtools flagstat`) |
| `prealigned_host_filter_reads` | reads, each mate counted              | PathSeq-T2T: PathSeqFilterSpark, `READS_AFTER_PREALIGNED_HOST_FILTER`                                     |
| `quality_complexity_filter_reads` | reads, each mate counted           | PathSeq-T2T: PathSeqFilterSpark, `READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER`                              |
| `host_kmer_filter_reads`    | reads, each mate counted                 | PathSeq-T2T: PathSeqFilterSpark, `READS_AFTER_HOST_FILTER`                                                 |
| `deduplication_filter_reads` | reads, each mate counted                | PathSeq-T2T: PathSeqFilterSpark, `READS_AFTER_DEDUPLICATION`                                               |
| `t2t_filter_reads`          | reads, each mate counted                 | PathSeq-T2T: reads not aligned to T2T-CHM13 (pairs and singletons), i.e. the reads given to the profilers |
| `motus4_<db_name>_pairs`    | read pairs                               | mOTUs v4, inserts aligned to the database and passing the mOTUs filters (`Read N aligned inserts`); with `--pathseq_t2t`, each singleton counts as one insert |
| `motus3_<db_name>_pairs`    | read pairs                               | mOTUs v3, inserts aligned to the database and passing the mOTUs filters (`Total number of inserts`); with `--pathseq_t2t`, each singleton counts as one insert |
| `cayman_<db_name>_pairs`    | read pairs                               | cayman, reads aligned to the gene catalogue and passing the filters (`Reads passing filters`)              |
| `gqlib_<db_name>_pairs`     | read pairs                               | gqlib, reads aligned to the gene catalogue and passing the filters (`Reads passing filters`)               |
| `kraken2_<db_name>_reads`   | reads, each mate counted                 | Kraken2, reads classified (the `root` line of the report; pairs counted twice, with the singletons of `--pathseq_t2t`) |

For cayman and gqlib, the alignment and its filters are the same for all the annotations of a gene catalogue, so there is one column per gene catalogue. A count that could not be found in a log is reported as `NA`, with a warning in the log of the task.

:::note
With `--shortread_qc_mergepairs`, fastp merges the overlapping pairs into single reads, so the columns after QC count fragments (merged pairs and unmerged reads) rather than pairs. Long reads have no QC and host removal counts.
:::

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

All tools in flexprofiler supported by MultiQC will have a dedicated section showing summary statistics of each tool based on information stored in log files.

You can expect in the MultiQC reports either sections and/or general stats columns for the following tools:

- fastqc (or falco)
- fastp
- bowtie2
- samtools (stats)
- kraken
- motus

:::info
The 'General Stats' table by default will only show statistics referring to pre-processing steps, and will not display possible values from each classifier/profiler, unless turned on by the user within the 'Configure Columns' menu or via a custom MultiQC config file (`--multiqc_config`).

To show the columns of a profiler in the general stats table, copy the zellerlab/flexprofiler [MultiQC config](https://github.com/zellerlab/flexprofiler/blob/master/assets/multiqc_config.yml) and change its entry in `table_columns_visible:` (e.g. `Kraken`) to True.
:::

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
