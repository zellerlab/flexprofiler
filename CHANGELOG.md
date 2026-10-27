# zellerlab/flexprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


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
- Documentation of the outputs of mOTUs v3, cayman, gqlib and bwa, and of the gqlib database

### Changed
- In the database sheet, the `db_path` of cayman is the gene catalogue fasta (indexed with `bwa index`) instead of the directory with the gene catalogue and its annotation, and the annotation `.csv` is given in the new `annot_db_path` column, with a name in `annot_db_name`
- cayman runs the gqlib domain profiler, which `cayman profile` (0.10.2) uses, directly on the shared alignment, with the same results: only the `.csv` annotations are supported, and the software versions report gqlib
- The tests use small subsets of the databases and of the host genome in `assets/test-data`, with nothing to download, and a single database sheet for all the profilers
- `zellerlab_shark` excludes the node res-hpc-mem02, where motus fails
- Lower resource requirements

### Fixed
- motus4 profiling of a single sample, which failed when merging the profiles
- The motus4 profiles change slightly between runs, because mOTUs 4 assigns at random, without a seed, the reads aligned equally well to several genes of a marker gene cluster: it is now seeded
- The runs of a sample are concatenated in the order of their run accession, instead of the order in which their processing finishes, so that the results of mOTUs do not depend on it
- bowtie2 host removal keeps the reads in the input order (`--reorder`), instead of the order in which its threads process them, for the same reason
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
