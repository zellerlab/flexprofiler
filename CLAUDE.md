# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Working rules

- Do not run tests or update snapshots (`nf-test test`, `--update-snapshot`) unless explicitly asked. Changes are often made many in a row, so make the requested change and stop.

## Overview

Nextflow DSL2 pipeline (zellerlab/flexprofiler), a stripped-down fork of nf-core/taxprofiler, for taxonomic and functional profiling of shotgun (WGS) and 16S (AMPLICON) data. It follows nf-core conventions: `modules/nf-core/` and `subworkflows/nf-core/` are installed from nf-core/modules (tracked in `modules.json`), and local code lives in `modules/local/` and `subworkflows/local/`.

## Commands

```bash
# all pipeline tests (each test is tagged "pipeline")
nf-test test --tag pipeline --profile +docker
# a single test, by file or by its tag (tags are named like the file, e.g. test_motus4)
nf-test test tests/test_motus4.nf.test --profile +docker
nf-test test --tag test_motus4 --profile +docker
# update snapshots (only when asked)
nf-test test tests/test_motus4.nf.test --profile +docker --update-snapshot

# run with a test profile directly
nextflow run . -profile test_motus4,docker --outdir results
```

- Tests are pipeline-level only (`tests/*.nf.test` + `.snap`), each using a `test_*` profile from `conf/test_*.config`. `tests/nextflow.config` caps resources (8 CPUs, 20 GB). Test data and the test samplesheet/database sheet are in `assets/test-data/`.
- Pre-commit runs prettier, trailing-whitespace and end-of-file fixers (excluding nf-core components and `.snap` files).
- `docs/parameters.md` is the rendered `nextflow_schema.json`. New params go in `nextflow.config` (defaults), the `params` block of `main.nf` (types: `Boolean`, `Integer`, `String?`, so that command-line values are converted), `nextflow_schema.json`, and the docs.
- `CHANGELOG.md` (Keep a Changelog, `Unreleased` section) is updated with detailed entries for user-facing changes.
- At every version bump (`manifest.version` in `nextflow.config`), also update `VERSION` in `utils/run_shark.sh` (the example SHARK launcher in `docs/zellerlab.md`) and in `utils/pull_images.sh` (pulls all the images into the shared apptainer cache) to the new release, and, if `nextflowVersion` was raised, their `NXF_VER` (the launcher of SHARK's Nextflow module downloads that version; SHARK has no module newer than 25.04.6), or its Nextflow module if SHARK gets one recent enough.

## Architecture

Entry: `main.nf` → `PIPELINE_INITIALISATION` (`subworkflows/local/utils_zellerlab_flexprofiler_pipeline`, parses and validates `--input` against `assets/schema_input.json` and `--databases` against `assets/schema_database.json`) → `FLEXPROFILER` (`workflows/flexprofiler.nf`) → `PIPELINE_COMPLETION`.

Data flow in `workflows/flexprofiler.nf`:

1. The samplesheet is split by read type into `meta.type` `short`/`long`. Only FASTQ is accepted, of Illumina and of PacBio AMPLICON runs: the samplesheet schema allows only these platforms, and `validatePacbio` in the utils subworkflow rejects PacBio runs that are not AMPLICON or are paired.
2. Optional `BBMAP_REPAIR`. `READS_QC` (FastQC/Falco) is called once, just before profiling, on the raw reads and on the reads given to the profilers (`meta.qc_stage` `raw`/`processed`).
3. Preprocessing in two layers: first by read type (`SHORTREAD_PREPROCESSING`, `LONGREAD_PREPROCESSING`), then by `meta.library_strategy` (`WGS_PREPROCESSING`: host removal; `AMPLICON_PREPROCESSING`: PIMENTO primer detection, cutadapt, orientation, pair handling).
4. `RUN_MERGING` concatenates runs per sample (`run_accession` becomes `merged_<library_strategy>` for the samples with more than one run).
5. `PATHSEQ_T2T` (`--pathseq_t2t`): stringent human host filtering. It is a passthrough when disabled and also emits singletons (orphaned mates), which go to Kraken2/mOTUs. See `docs/pathseq_t2t.md`. Related pipelines:
   - Upstream [abdohlman/pathseq-t2t](https://github.com/abdohlman/pathseq-t2t) (commit `51d7443`), shell scripts that take a BAM already aligned to GRCh38. `PATHSEQ_T2T` reimplements its filtering commands (`prefilter`, `qcfilter`, `t2tfilter`) as Nextflow modules. `utils/upstream_pathseq_t2t/` records how upstream was run to compare it with flexprofiler, and the issues found in it.
   - Fabian's pipeline, [fabispringer/pathseq-t2t-nextflow](https://github.com/fabispringer/pathseq-t2t-nextflow) (commit `9ef48f45`), a Nextflow wrapper of upstream that starts from FASTQ. The STAR route (`--pathseq_t2t_grch38_aligner star`, two-pass) was modelled on its `ALIGN_STAR`. It is a reference implementation for comparison, not a dependency. It trims and filters the reads with BBDuk before the alignment (flexprofiler only trims the adapters, with fastp, in the `zellerlab_pathseq_t2t` profiles), and it unaligns the reads before PathSeqFilterSpark. The differences are listed in `docs/pathseq_t2t.md` ("Differences from the STAR route of pathseq-t2t-nextflow").
6. `PROFILING` (`subworkflows/local/profiling.nf`): reads are combined with database rows by `type` and by the profiler's library strategy (`toolLibraryStrategy()`). Profilers are mOTUs v3/v4, Cayman, gqlib and Kraken2 (WGS), and IDTAXA (AMPLICON). A tool runs only when `params.run_<tool>` is set and the database sheet has a row for it. Mismatches (an enabled tool without databases or without reads, reads without any database) are errors at startup (`validateDatabaseMatching`), unless no tool is enabled (preprocessing only). Cayman and gqlib share one BWA alignment per gene catalogue.
7. `READ_COUNTING`: every step emits `read_count_logs` (`[meta + tool, log]`), which `bin/read_counts.py` collates into `read_counts/read_counts.tsv`.
8. `STANDARDISATION_PROFILES` merges per-sample profiles into tables (`bin/*_merge.py`), then `MULTIQC_REPORTING`.

Conventions:

- Subworkflows emit `reads`, `mqc`, `read_count_logs` and `published`, which the main workflow mixes together. Software versions are not passed through channels: every module (nf-core and local) emits `tuple val("${task.process}"), val(<tool>), eval(<cmd>), emit: versions_<tool>, topic: versions`, collected by `MULTIQC_REPORTING` from `channel.topic("versions")`.
- All process options (`ext.args`, `ext.prefix`, resources) are in `conf/modules.config` via `withName`, not in the modules.
- Results are published with workflow outputs, not `publishDir`: each subworkflow emits `published`, a channel of `[file, target]` built with `publishFiles(channel, directory, pattern)` (utils subworkflow) next to the step that makes the files, `FLEXPROFILER` mixes them into one `published` channel, and the single `results` output of `main.nf` publishes it (`source >> target`) to `outputDir` (`--outdir`). Only the entry workflow can publish, so a new published channel must be emitted up to `main.nf`. Files outside the work directory are not published: the amplicon analysis-ready reads, which can be the input FASTQs, are emitted apart (`analysis_ready`) to their own output, the only one with `includeInputs true`. The layout is in `docs/output.md`.
- Databases ending in `.tar.gz` are untarred once per file, only if their tool is enabled.
- `conf/zellerlab.config` and `assets/zellerlab_*_databases.csv` hold the lab's institutional profiles (e.g. `zellerlab_shark`). `utils/` has helper scripts (database download, SLURM launcher, samplesheet maker).
- `TODO.md` lists known limitations and differences from upstream tools.
