# zellerlab/flexprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## Unreleased

### Added
- PathSeq-T2T host filtering for low-biomass human samples (`--pathseq_t2t`), after Dohlman et al., Cell 2026 (https://doi.org/10.1016/j.cell.2026.04.015) and https://github.com/abdohlman/pathseq-t2t (commit `51d74430`): alignment to GRCh38 with `bwa mem -T 0` and selection of the reads not aligned to it (with the reads of its viral decoy contigs), GATK PathSeqFilterSpark (quality, complexity, host k-mer and duplicate filters), and removal of the reads aligned to T2T-CHM13. It replaces the short-read QC and host removal, accepts only Illumina `WGS` FASTQ, cannot be used with Cayman and gqlib, and runs after run merging. New parameters `--pathseq_t2t_grch38_reference`, `--pathseq_t2t_grch38_bwa_index`, `--pathseq_t2t_chm13_reference`, `--pathseq_t2t_chm13_bwa_index`, `--pathseq_host_dir`, `--pathseq_decoys_bed`, `--pathseq_min_clipped_read_length` and `--save_pathseq_t2t_bams`; the reads after each filter in the read counts; FastQC of the filtered reads; `test_pathseq_t2t` profile and `tests/test_pathseq_t2t.nf.test`. The differences from upstream, including the upstream bugs that are fixed, are listed in `docs/pathseq_t2t.md`
- Kraken2 (`--run_kraken2`, `kraken2` rows in the database sheet), with `--report-minimizer-data` by default, and with `--pathseq_t2t` the confidence threshold of PathSeq-T2T (`--confidence 0.15`, the `db_params` of the kraken2 databases being rejected): the reads left without their mate by `--pathseq_t2t` are classified too (they also go to mOTUs, as unpaired reads in the same run as the pairs), and the reports are merged into tables of reads (each mate counted, the singletons once) per taxon, with the classified reads in the read counts (`kraken2_<db_name>_reads`) and the reports in MultiQC
- `utils/compare_upstream_pathseq_t2t.sh`, to compare the PathSeq-T2T mode with upstream PathSeq-T2T read by read (on the test data: the same reads after every step and the same Kraken2 reports, see `docs/pathseq_t2t.md`), and `utils/upstream_pathseq_t2t/`, the recipe to run upstream (environment, patch, commands) and its known issues
- PacBio (`PACBIO_SMRT`) FASTQ input for the `AMPLICON` runs, adapted from the best practice for full-length 16S (Callahan et al. 2019, PacBio pb-16S-nf): quality filters with vsearch (`--fastq_maxee 10 --fastq_maxns 0`, relaxed from the maxEE 2 of Callahan et al. 2019) and BBDuk (`mbq=3`), with `--perform_longread_qc`, orientation with `vsearch --orient` against `--amplicon_orient_reference` before PIMENTO, and primer trimming as a pair, both primers required (cutadapt linked adapter, the reads without both primers discarded), in the orientation of the forward primer, with a length window of 1000-1600 bp. IDTAXA classifies only the top strand of the long reads, oriented as the gene (both strands of the short reads)
- `--amplicon_pimento_primers_pacbio`: primers added to the PIMENTO standard library for the PacBio AMPLICON runs, by default the full-length 16S primers of PacBio (Kinnex 16S: 27F `AGRGTTYGATYMTGGCTCAG`, 1492R `RGYTACCTTGTTACGACTT`), not in the PIMENTO library
- `--amplicon_primer_fwd` and `--amplicon_primer_rev`, to give the primers of the `AMPLICON` runs instead of finding them with PIMENTO
- The read counts of the long amplicon reads: `raw_reads` (vsearch) and `qc_reads` (BBDuk), `oriented_reads`, and `primer_trimmed_reads` (with both primers, in the length window)
- `--amplicon_pair_handling none`, to keep the read pairs of the amplicon reads (IDTAXA cannot be run then)
- PacBio amplicon reads in the test data, and a test of them (`tests/test_pacbio.nf.test`)

### Changed
- The preprocessing is split in two stages: of the read type (`SHORTREAD_PREPROCESSING`: QC; `LONGREAD_PREPROCESSING`: PacBio quality filtering), then of the library strategy (`WGS_PREPROCESSING`: host removal; `AMPLICON_PREPROCESSING`: primer trimming and pair handling), then run merging (`RUN_MERGING`); the read counts and the MultiQC report (with the software versions) are made in their own subworkflows (`READ_COUNTING`, `MULTIQC_REPORTING`), from the read count logs emitted by each subworkflow
- The primers of the amplicon reads are trimmed after the read QC (fastp), instead of before
- `--idtaxa_pair_handling` is renamed `--amplicon_pair_handling`, and the read pairs are handled in the amplicon preprocessing, before run merging, for all the amplicon profilers: the analysis-ready reads of the `AMPLICON` runs are the merged (or joined) pairs. `FASTP_IDTAXA_MERGE` is renamed `FASTP_AMPLICON_MERGE` (output in `fastp_amplicon_merge/`), and the vsearch and fastp merge outputs are named per run
- The analysis-ready reads of the `AMPLICON` runs are named `<sample_id>_<run_accession>[_{1,2}].fastq.gz`
- The read counts columns are in the order of the steps: `raw_reads`, `qc_reads`, then `primer_trimmed_reads` and `readthrough_trimmed_reads`
- Unsupported inputs fail the samplesheet validation, before the pipeline starts, instead of their runs being skipped with a warning: `instrument_platform` accepts only `ILLUMINA` and `PACBIO_SMRT`, and the `PACBIO_SMRT` runs must be `AMPLICON` (not `WGS`, the default `library_strategy`) and single-end
- The database sheet is checked against the samplesheet and the enabled profilers before the pipeline starts, as the reads and the databases that do not match are dropped silently: it is an error if an enabled profiler (`--run_<tool>`) has no database in the sheet, if a run matches no database of the enabled profilers (by read type, from `instrument_platform` and `db_type`, and by library strategy), and if an enabled profiler has no run to match, which replaces the warning given when there were no runs of its library strategy. Without any enabled profiler, the reads are only preprocessed, as before
- FastQC/Falco runs on the raw reads and on the reads given to the profilers, after all the preprocessing (and run merging and PathSeq-T2T), in `fastqc/{raw,processed}/` (or `falco/`), instead of on the short reads after adapter trimming (`processed`) and on the reads after PathSeq-T2T (`pathseq_t2t`): the processed reads are checked also without `--perform_shortread_qc`, and for the long reads. The reports are named `<sample>_<run_accession>_<stage>`, also after run merging, and the files without reads are left out
- `LICENSE`: added the copyright notice of zellerlab/flexprofiler (Saul Pierotti, Zeller Lab, LUMC), above the notice of nf-core/taxprofiler, which is kept as the MIT license requires and now names the version of the fork (v1.2.5)
- The results are published with the workflow outputs of Nextflow (the `output` block of `main.nf`) instead of `publishDir` in `conf/modules.config`, which Nextflow plans to deprecate: each subworkflow chooses the files of its steps to publish, and their directory. The layout of the output directory is unchanged. `--publish_dir_mode` sets `workflow.output.mode`
- The types of the parameters are declared in a `params` block in `main.nf`, so that the values given on the command line are converted to them: with the strict syntax of Nextflow (the default since 26.04), `--<flag> false` was the string `'false'`, which is true, and passed the parameter validation, so a flag could not be turned off from the command line (`-params-file` and config files were not affected). An invalid value (e.g. `--perform_shortread_qc maybe`) is now an error
- Removed the Nextflow options deprecated in 25.10-26.10 (`nextflow.enable.strict`, `nextflow.enable.configProcessNamesValidation`, `manifest.defaultBranch`), the `when` sections of the local modules, and the implicit closure parameters (`it`) and `Channel` factories of the local code
- The Zeller lab setup uses Apptainer instead of Singularity: the `zellerlab_shark` profile runs the containers with Apptainer (`-profile apptainer` is not needed with it), with the images in a cache shared by the lab (`/exports/lucid-grpzeller-work/SOFTWARE_CACHE/flexprofiler/apptainer`; `zellerlab_alice` is marked as incomplete), `utils/run_shark.sh` loads the Apptainer module of SHARK, and sets `NXF_VER` to the Nextflow version the pipeline requires, which SHARK has no module for (the launcher of its Nextflow module downloads it); `utils/user.config` sets the Apptainer cache. The local modules use the Galaxy depot images with Apptainer too, as the nf-core ones do, instead of converting the Docker images
- All the local modules use Seqera Containers images (`community.wave.seqera.io`, and their Singularity images by https), built from the same conda packages, instead of the BioContainers and Galaxy depot ones: IDTAXA no longer needs the Docker image converted (the Galaxy depot has no DECIPHER 3.6.0 image), and the GATK image of PathSeq-T2T pins Java 17 itself
- `utils/pull_images.sh`, to download all the container images of the pipeline into the apptainer cache shared by the lab on SHARK
- The local modules report their software versions to the `versions` topic channel (`emit: versions_<tool>, topic: versions`, the version obtained with `eval`), like the nf-core modules, instead of writing a `versions.yml`; the `versions` channels are removed from the subworkflows. The mOTUs v4 merge reports its version as `motus` also in stub runs (was `motus4`)

- Falco is the default read QC tool (`--preprocessing_qc_tool falco`, was `fastqc`): it is a drop-in replacement of FastQC, faster and with much less memory (FastQC was given 2 CPUs and 15 GB), and its reports are read by MultiQC like those of FastQC. The reports are in `falco/{raw,processed}/` instead of `fastqc/`; `--preprocessing_qc_tool fastqc` gives the previous behaviour. Falco runs with 1 CPU and 4 GB (times the attempt)

### Removed
- Oxford Nanopore input (`OXFORD_NANOPORE` is no longer an `instrument_platform` of the samplesheet) and its preprocessing: Porechop, Porechop_ABI, Filtlong and nanoq, with `--longread_adapterremoval_tool`, `--longread_qc_skipadaptertrim`, `--longread_qc_skipqualityfilter`, `--longread_filter_tool` and `--longread_qc_qualityfilter_*`. `--perform_longread_qc` now only turns on the PacBio quality filtering
- FASTA input: the `fasta` column of the samplesheet (FASTA runs were already skipped with a warning). `fastq_1` is now required, and samplesheets that still have a `fasta` column should drop it
- The long-read host removal with minimap2 (`--perform_longread_hostremoval`, `--longread_hostremoval_index`)
- The redundancy estimation with Nonpareil (`--perform_shortread_redundancyestimation`, `--shortread_redundancyestimation_mode`)
- The completion notifications to Microsoft Teams or Slack (`--hook_url`, `HOOK_URL`)
- AdapterRemoval as short-read preprocessing tool (`--shortread_qc_tool`, which is removed): the short reads are preprocessed with fastp. Adapter lists (`--shortread_qc_adapterlist`) are FASTA only
- The short-read complexity filtering with BBDuk or fastp (`--perform_shortread_complexityfilter`, `--shortread_complexityfilter_tool`, `--shortread_complexityfilter_entropy`, `--shortread_complexityfilter_bbduk_*`, `--shortread_complexityfilter_fastp_threshold`, `--save_complexityfiltered_reads`) and the PRINSEQ++ one (`--shortread_complexityfilter_tool prinseqplusplus`, `--shortread_complexityfilter_prinseqplusplus_*`). PathSeq-T2T still filters the low-complexity reads itself
- Leftovers of nf-core/taxprofiler: the institutional `pipeline/taxprofiler.config` of nf-core/configs is no longer loaded, and the MultiQC config, `tower.yml`, citations and methods text no longer list the profilers of nf-core/taxprofiler that flexprofiler does not run (Bracken, KrakenUniq, Centrifuge, Kaiju, DIAMOND, MALT, MetaPhlAn, ganon, KMCP, Krona, TAXPASTA)

### Fixed
- In MultiQC, the post-trimming FastQC/Falco section (and the General Stats) kept only one mate of each paired-end sample, and not always the same one: the `_processed` name cleaning cut the names at the stage, dropping `_1`/`_2`, so that both mates (and the PathSeq-T2T singletons) had the same name. The stage is now removed with a regex that keeps the mates, and only from the FastQC/Falco sample names: the samples of both sections are named `<sample>_<run_accession>[_singletons][_1|_2]`
- The FastQC/Falco reports were split into the MultiQC sections by `raw` or `processed` anywhere in their path, so a sample or run name containing either (e.g. `straw`) put its reports in the wrong section or in none. The stage now ends the report names (`<sample>_<run_accession>[_singletons]_<stage>`, was `<sample>_<run_accession>_<stage>[_singletons]`), and the sections match it only right before the endings written by the tools
- Falco reports of paired-end reads were named after the input files instead of the prefix, so they were left out of both MultiQC sections, and MultiQC named all Falco samples after the input files. The nf-core Falco module is patched to link each file to `<prefix>` or `<prefix>_{1,2}` (as the FastQC module does) and to run Falco once per file, writing `<prefix>[_{1,2}]_{fastqc_data,summary}.txt` and `<prefix>[_{1,2}]_report.html`; the `_falco` suffix of the Falco reports is dropped (they are in `falco/`). The stub of the module is fixed too
- The reverse primer found by PIMENTO in single-end reads (`--merged`), alone or in a pair, was complemented before trimming, so that cutadapt never found it: PIMENTO writes it complemented only for the paired-end mates

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
