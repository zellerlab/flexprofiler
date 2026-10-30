# zellerlab/flexprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## v0.0.5 [30-10-2026]

### Added
- `library_strategy` column in the samplesheet, `WGS` (shotgun metagenomics) or `AMPLICON` (e.g. 16S), from the ENA controlled vocabulary; `WGS` if empty or missing. `utils/make_samplesheet.sh` writes it
- Taxonomic classification of the `AMPLICON` samples with IDTAXA (DECIPHER) (`--run_idtaxa`), against the training sets of the `idtaxa` rows of the database sheet (`--idtaxa_db_names` to choose among them), with the confidence threshold `--idtaxa_threshold` (default 40). The profiles are merged across samples into a table of the whole lineages and one table per rank (`--idtaxa_ranks` to choose the ranks)
- The read pairs of paired-end amplicon reads are merged before IDTAXA, and the pairs that do not overlap are handled as set by `--idtaxa_pair_handling`: joined with a padding of Ns by vsearch (`merged_plus_joined`, default), or kept as two sequences by fastp (`merged_plus_unmerged`)
- Primer detection with PIMENTO and primer and read-through trimming with cutadapt for the `AMPLICON` runs (`--perform_amplicon_primertrimming`, on by default), with the reads trimmed in the read counts (`primer_trimmed_reads`, `readthrough_trimmed_reads`)
- `zellerlab_amplicon` profile and `test_idtaxa` profile
- The reads given to mOTUs v4 with `--motus4_rename_reads` are saved with `--save_analysis_ready_fastqs`, in `analysis_ready_fastqs/motus4`
- Tests for IDTAXA, primer removal, adapter removal, host removal, and a mix of `WGS` and `AMPLICON` samples in one run

### Changed
- Short-read host removal aligns the two mates of paired-end reads separately with bowtie2, as single-end reads, and keeps a pair only if neither mate aligns to the host, as recommended in https://doi.org/10.1016/j.crmeth.2025.101218. Before, the pairs not aligned concordantly were kept (`bowtie2 --un-conc`), including those with one mate aligned to the host. The unaligned mates are paired again with `repair.sh`, whose log gives the `hostremoval_retained_pairs` read count (before `host_removed_pairs`)
- The `AMPLICON` runs are not host-removed, and Nonpareil is run only on the `WGS` runs
- The `WGS` and `AMPLICON` runs of a sample are merged separately with `--perform_runmerging`: a merged run is named `merged_WGS` or `merged_AMPLICON` (also in the names of the concatenated FASTQ files), and a sample with a single run of a strategy keeps its run accession
- The read counts have one row per run (columns `sample` and `run_accession`), or per merged run with `--perform_runmerging`, instead of one per sample
- The gqlib and cayman profiles are merged by the pipeline instead of `collate_counts` of gqlib, with the same tables, but with 0 instead of NA for the features without reads in a sample
- The profiles merged by gqlib and mOTUs v4 are listed in a file, so that the number of samples is not limited by the length of a command line

## v0.0.4 [27-10-2026]

### Added
- gqlib functional profiling of gene catalogues (`--run_gqlib`), with the profiles merged across samples for each functional category (KEGG, CAZy, COG, GO, ...)
- The reads are aligned to each gene catalogue once with `bwa mem`, and the alignment is shared by cayman and gqlib and by all the annotations of the catalogue (`--save_bwa_alignments` to keep it)
- Columns `annot_db_name` and `annot_db_path` in the database sheet, for the annotations of the cayman and gqlib gene catalogues
- Per-sample read counts after QC, host removal and alignment to each profiling database, in `read_counts/read_counts.tsv` and in the MultiQC report
- `zellerlab_minimal` profile, with the new parameters `--skip_multiqc`, `--skip_postprocessing_qc` and `--skip_hostremoval_qc`
- `zellerlab_alice` profile for the ALICE cluster
- Tests for gqlib, the minimal profile, read names with /1 and /2 suffixes, paired-end reads in a different order in the two files, a mix of these in one run, and samples with paired-end and single-end runs, and `test_all` profile with all the profilers
- Count checks of the profiles and of the read counts in the tests
- `utils/fetch_databases.sh` builds the gqlib database from the eggNOG-mapper annotation of the gene catalogue (not tested; `SKIP_GQLIB=1` to skip it)
- Documentation of the outputs of mOTUs v3, cayman, gqlib and bwa, of the gqlib database, and of how gqlib and cayman count the reads (adapted from the gffquant documentation)

### Changed
- In the database sheet, the `db_path` of cayman is the gene catalogue fasta (indexed with `bwa index`) instead of the directory with the gene catalogue and its annotation, and the annotation `.csv` is given in the new `annot_db_path` column, with a name in `annot_db_name`
- cayman runs the gqlib domain profiler, which `cayman profile` (0.13.0) uses, directly on the shared alignment, with the same results: only the `.csv` annotations are supported, and the software versions report gqlib
- The tests use small subsets of the databases and of the host genome in `assets/test-data`, with nothing to download, and a single database sheet for all the profilers
- `zellerlab_shark` excludes the node res-hpc-mem02, where motus fails
- Lower resource requirements

### Fixed
- motus4 profiling of a single sample, which failed when merging the profiles
- The motus4 profiles change slightly between runs, because mOTUs 4 assigns at random, without a seed, the reads aligned equally well to several genes of a marker gene cluster: it is now seeded
- The runs of a sample are concatenated in the order of their run accession, instead of the order in which their processing finishes, so that the results of mOTUs do not depend on it
- bowtie2 host removal keeps the reads in the input order (`--reorder`), instead of the order in which its threads process them, for the same reason
- The samples are in the same order (by name) in the merged profiles, instead of the order in which their profiling finishes
- Merging of the gqlib and cayman profiles when only some samples have ambiguous reads (and so the `combined_` columns): the merged `combined_` tables were missing if the first sample had no ambiguous reads, and had the `uniq_` counts of the samples following one without ambiguous reads otherwise (a bug in `collate_counts` of gqlib 2.14.4)
- The settings of cayman and gqlib (ambiguous reads and alignment filters) are given explicitly, instead of relying on the defaults of gqlib
- `utils/fetch_databases.sh` writes the database sheet in the current format, saves the mOTUs v3 database to its own directory instead of over the mOTUs v4 one, and stops at the first error
- Outdated and wrong paths, profiles and descriptions in the documentation
- Minor bug fixes

### Removed
- Unused kraken2 and krona local modules

## v0.0.3 [31-05-2026]
- Addition of mOTUs v3 profiling
- nf-test tests for mOTUs v4, mOTUs v3 and cayman
- Support for Nextflow 26 and update of the nf-core modules
- Profiling is disabled for non-Illumina reads, and runs of different technologies cannot be merged
- BBmap repair passes single-end reads through
- MultiQC uses a scratch directory outside `/dev/shm`

## v0.0.2 [06-03-2026]
- Addition of cayman functional profiler
- Addition of BBmap repair to fix mismatched PE reads upstream of profiling
- Addition of the MOTUS_REPAIR process to fix read tags before mOTUs (the tool expects identical tag names from PE reads)

## v0.0.1 [24-02-2026]
- Full rebranding to zellerlab/flexprofiler
- Dropped support for all profilers until testing
- Implementation of zellerlab-specific nextflow profiles
- Addition of utility scripts in utils
- Implementation of mOTUs v4.0.4

## v0.0.0 [06-02-2026]
- Initial fork from **nf-core/taxprofiler** v1.2.5
