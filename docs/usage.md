# zellerlab/flexprofiler: Usage

## Introduction

zellerlab/flexprofiler is a pipeline for highly-parallelised taxonomic anf functional classification and profiling of shotgun metagenomic and 16S data across multiple tools simultaneously.
In addition to multiple classification and profiling tools, at the same time it allows you to performing taxonomic classification and profiling across multiple databases and settings per tool, as well as produces standardised output tables to allow immediate cross comparison of results between tools.

## General Usage

To run zellerlab/flexprofiler, at a minimum you require two inputs:

- a sequencing read samplesheet
- a database samplesheet

Both contain metadata and paths to the data of your input samples and databases.
When running zellerlab/flexprofiler, every step and tool is 'opt in'.
To run a given classifier or profiler you must make sure to supply both a database in your `<database>.csv` and supply `--run_<profiler>` flag to your command.
Omitting either will result in the profiling tool not executing.
zellerlab/flexprofiler also includes optional pre-processing (adapter clipping, merge running etc.) or post-processing (visualisation) steps.
These are also opt in with a `--perform_<step>` flag.
In some cases, the pre- and post-processing steps may also require additional files. Please check the parameters tab of this documentation for more information.
Please see the rest of this page for information about how to prepare input samplesheets and databases and how to run Nextflow pipelines.
See the [parameters](docs/parameters.md) documentation for more information about specific options the pipeline also offers.

## Samplesheet inputs

zellerlab/flexprofiler can accept as input raw or preprocessed single- or paired-end short-read (e.g. Illumina) FASTQ files, long-read FASTQ files (e.g. Oxford Nanopore), or FASTA sequences (available for a subset of profilers).
For PACBIO_SMRT data, please convert BAM files to FASTQ format before input.
You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline.
Use this parameter to specify its location.
It has to be a comma-separated file with 6 columns, and a header row as shown in the examples below.
Furthermore, zellerlab/flexprofiler also requires a second comma-separated file of 3 columns with a header row as in the examples below.
This samplesheet is then specified on the command line as follows:

```console
--input '[path to samplesheet file]' --databases '[path to database sheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth.
The pipeline will concatenate different runs FASTQ files of the same sample before performing profiling, when `--perform_runmerging` is supplied.
Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
2612,lane1,ILLUMINA,2612_lane1_R1.fq.gz,2612_lane1_R2.fq.gz,
2612,lane2,ILLUMINA,2612_lane2_R1.fq.gz,2612_lane2_R2.fq.gz,
2612,lane3,ILLUMINA,2612_lane3_R1.fq.gz,,
```

::: info
Please note that the column name `run_accession` follows the definition of an ENA 'run'.
A 'run' corresponds to a single or paired-end set of demultiplexed FASTQs.
Given that demultiplexing of a given library happens per lane, each sequencing pair from each lane is a 'run'.
Therefore, for each sample, you may get multiple 'runs' consisting of _both_ lanes (of the same library) _and_ sequencing libraries.
Therefore ensure that each `run_accession` ID is unique, even if from the same sample!
:::

:::warning
Runs of the same sample sequenced on Illumina platforms with a combination of single and paired-end data will **not** be run-wise concatenated, unless pair-merging is specified. In the example above, `run3` will be profiled independently of `run1` and `run2` if pairs are not merged.
:::

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet.
The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 6 columns to match those defined in the table below.
A final samplesheet file consisting of both single- and paired-end data, as well as long-read FASTA files may look something like the one below.
This is for 6 samples, where `2612` has been sequenced twice.

```csv title="samplesheet.csv"
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
2611,ERR5766174,ILLUMINA,,,/<path>/<to>/fasta/ERX5474930_ERR5766174_1.fa.gz
2612,ERR5766176,ILLUMINA,/<path>/<to>/fastq/ERX5474932_ERR5766176_1.fastq.gz,/<path>/<to>/fastq/ERX5474932_ERR5766176_2.fastq.gz,
2612,ERR5766180,ILLUMINA,/<path>/<to>/fastq/ERX5474936_ERR5766180_1.fastq.gz,,
2613,ERR5766181,ILLUMINA,/<path>/<to>/fastq/ERX5474937_ERR5766181_1.fastq.gz,/<path>/<to>/fastq/ERX5474937_ERR5766181_2.fastq.gz,
ERR3201952,ERR3201952,OXFORD_NANOPORE,/<path>/<to>/fastq/ERR3201952.fastq.gz,,
```

:::warning
Input FASTQ and FASTA files _must_ be gzipped.
:::

:::warning
While one can include both short-read and long-read data in one run, we recommend that you split these across _two_ pipeline runs and database sheets (see below).
This will allow classification optimisation for each data type, and make MultiQC run-reports more readable (due to run statistics having vary large number differences).
:::

| Column                | Description                                                                                                                                                                                               |
| --------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`              | Unique sample name [required].                                                                                                                                                                            |
| `run_accession`       | Run ID or name unique for each (pairs of) file(s) .Can also supply sample name again here, if only a single run was generated [required].                                                                 |
| `instrument_platform` | Sequencing platform reads generated on, selected from the EBI ENA [controlled vocabulary](https://www.ebi.ac.uk/ena/portal/api/controlledVocab?field=instrument_platform) [required].                     |
| `fastq_1`             | Path or URL to sequencing reads or for Illumina R1 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if data in FASTA is specified. Cannot be combined with `fasta`. |
| `fastq_2`             | Path or URL to Illumina R2 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if single end data. Cannot be combined with `fasta`.                                    |
| `fasta`               | Path or URL to long-reads or contigs in FASTA format. GZipped compressed files accepted. Can be left empty if data in FASTA is specified. Cannot be combined with `fastq_1` or `fastq_2`.                 |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

:::warning
FASTA input will not go through any preprocessing steps, and will go directly to profiling.
:::

### Full database sheet

zellerlab/flexprofiler supports multiple databases being classified/profiled against in parallel for each tool.
Databases can be supplied either in the form of a compressed `.tar.gz` archive of a directory containing all relevant database files or the path to a directory on the filesystem.

:::warning
zellerlab/flexprofiler does not provide any databases by default, nor does it currently generate them for you.
This must be performed manually by the user.
See the [databases documentation](databases.md) for details.
:::

The pipeline takes the paths and specific classification/profiling parameters of the tool of these databases as input via a four (or five) column comma-separated sheet.
The optional `db_type` column allows to use specific database/parameters against specific data types.
By specifying if a database is for short-or long-reads (or even both), the samples sequenced with Illumina are combined with the short-read databases and the samples sequenced with Nanopore are combined with long-read databases.
If `db_type` is not provided, it is assumed the database and parameters are applicable for both short- and long-read data.

:::warning
To allow user freedom, zellerlab/flexprofiler does not check for mandatory or the validity of non-file database parameters for correct execution of the tool - excluding options offered via pipeline level parameters!
:::

An example database sheet can look as follows.
Since the `db_type` column is missing, it is therefore assumed that the database and parameters are suitable for both short- and long-read data.
In the second example database sheet, the `db_type` column has been provided.
The valid options are `short`, `long` and `short;long`.

```csv
tool,db_name,db_params,db_path
motus4,db_mOTU,,/path/to/db_mOTU
cayman,GMGC10.human-gut.95nr.0.5.percent.prevalence,,/path/to/GMGC10.human-gut.95nr.0.5.percent.prevalence
```

```csv
tool,db_name,db_params,db_type,db_path
motus4,db_mOTU,,short,/path/to/db_mOTU
cayman,GMGC10.human-gut.95nr.0.5.percent.prevalence,,long,/path/to/GMGC10.human-gut.95nr.0.5.percent.prevalence
```

Column specifications are as follows:

| Column      | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `tool`      | Profiling tool (supported by zellerlab/flexprofiler) that the database has been indexed for [required].                                                                                                                                                                                                                                                                                                                                                                                                              |
| `db_name`   | A unique name per tool for the particular database [required].                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `db_params` | Any parameters of the given classifier/profiler that you wish to specify that the classifier/profiling tool should use when profiling against this specific database. Can be empty to use classifier/profiler defaults. Must not be surrounded by quotes [required]. We generally do not recommend specifying parameters here that turn on/off saving of output files or specifying particular file extensions - this should be already addressed via pipeline parameters. |
| `db_type`   | An optional column to distinguish between short- and long-read databases. If the column is empty, the pipeline will assume all databases (and their settings specified in `db_params`!) will be applicable for both short and long read data. Possible values: `long`, `short`, or `short;long`. If the `db_type` column is missing from the database.csv, it will take the default value short;long                                                                                                                                                                                                        |
| `db_path`   | Path to the database. Can either be a path to a directory containing the database index files or a `.tar.gz` file which contains the compressed database directory with the same name as the tar archive, minus `.tar.gz` [required].                                                                                                                                                                                                                                                                                                                                                                       |

:::tip
You can also specify the same database directory/file twice (ensuring unique `db_name`s) and specify different parameters for each database to compare the effect of different parameters during classification/profiling.
:::

zellerlab/flexprofiler will automatically decompress and extract any compressed archives for you.
The optional `db_type` column enables the use of specific databases or parameters for different data types.
By specifying if a database is for short-reads, long-reads, or both, Illumina samples are combined with short-read databases, while Nanopore samples are combined with long-read databases.

:::tip
Click the links in the list below for short quick-reference tutorials how to generate download 'pre-made' and/or custom databases for each tool.
:::

The database paths (`db_path`) for each tool are expected to contain:

- [**mOTUs v4**:](databases.md#motus4) the directory `db_mOTU/` that is downloaded via `motus4 downloadMGDB`.
  - Note that you must use `motus downloadMGDB` and if installed via `conda`, will be placed in a specific `site-package` directory in the conda environment. For more details see the [mOTUs database documentation](databases.md#motus4).

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run zellerlab/flexprofiler --input samplesheet.csv --databases databases.csv --outdir <OUTDIR> -profile docker --run_<TOOL1> --run_<TOOL2>
```

This will launch the pipeline with the `docker` configuration profile.
See below for more information about profiles.

When running zellerlab/flexprofiler, every step and tool is 'opt in'.
To run a given classifier/profiler you must make sure to supply both a database in your `<database>.csv` and supply `--run_<profiler>` flag to your command.
Omitting either will result in the classification/profiling tool not executing.
If you wish to perform pre-processing (adapter clipping, merge running etc.) or post-processing (visualisation) steps, these are also opt in with a `--perform_<step>` flag.
In some cases, the pre- and post-processing steps may also require additional files. Please check the parameters tab of this documentation for more information.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.
Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run zellerlab/flexprofiler -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

### Sequencing quality control

[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads.
It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.
zellerlab/flexprofiler offers [`falco`](https://github.com/smithlabcode/falco) as an drop-in replacement, with supposedly better improvement particularly for long reads.

### Preprocessing Steps

zellerlab/flexprofiler offers four main preprocessing steps for preprocessing raw sequencing reads:

- [**Paired-end read repair**](#paired-end-repair): correctly pair read order in PE files that have become unpaired.
- [**Read processing**](#read-processing): adapter clipping and pair-merging.
- [**Redundancy estimation**](#redundancy-estimation): short-read metagenome coverage estimation.
- [**Complexity filtering**](#complexity-filtering): removal of low-sequence complexity reads.
- [**Host read-removal**](#host-read-removal): removal of reads aligning to reference genome(s) of a host.
- [**Run merging**](#run-merging): concatenation of multiple FASTQ chunks/sequencing runs/libraries of a sample.

:::info
You can save the 'final' reads used for classification/profiling from any combination of these steps with `--save_analysis_ready_fastqs`.
:::

#### Pared-read repair
If sequencing files are processed with tools that are not PE-aware, reads may become unpaired (that is, not in the same order) across PE files.
This breaks almost every PE tool, so PE repair is done using [BBDUK's repair.sh](https://github.com/BioInfoTools/BBMap/blob/master/sh/repair.sh) BEFORE any other step when `--perform_repair_pe_order` is set.
This is also before FASTQC raw is run, because with incorrectly paired reads the FASTQC output is nonsensical.

#### Read Processing

Raw sequencing read processing in the form of adapter clipping and paired-end read merging can be activated via the `--perform_shortread_qc` or `--perform_longread_qc` flags.

It is highly recommended to run this on raw reads to remove artifacts from sequencing that can cause false positive identification of taxa (e.g. contaminated reference genomes) and/or skews in taxonomic abundance profiles.
If you have public data, normally these should have been corrected for, however you should still check that these steps have indeed been already performed.

There are currently two options for short-read preprocessing: [`fastp`](https://github.com/OpenGene/fastp) or [`adapterremoval`](https://github.com/MikkelSchubert/adapterremoval).

For adapter clipping, you can either rely on the tool's default adapter sequences, or supply your own adapters (`--shortread_qc_adapter1` and `--shortread_qc_adapter2`)
By default, paired-end merging is not activated.
In this case paired-end 'alignment' against the reference databases is performed where supported, and if not, supported pairs will be independently classified/profiled.
If paired-end merging is activated you can also specify whether to include unmerged reads in the reads sent for classification/profiling (`--shortread_qc_mergepairs` and `--shortread_qc_includeunmerged`).
You can also turn off clipping and only perform paired-end merging, if requested.
This can be useful when processing data downloaded from the ENA, SRA, or DDBJ (`--shortread_qc_skipadaptertrim`).
Both tools support length filtering of reads and can be tuned with `--shortread_qc_minlength`.
Performing length filtering can be useful to remove short (often low sequencing complexity) sequences that result in unspecific classification and therefore slow down runtime during classification/profiling, with minimal gain.

There are currently two options for long-read Oxford Nanopore processing: [`porechop`](https://github.com/rrwick/Porechop), [`porechop_abi`](https://github.com/bonsai-team/Porechop_ABI).

For both short-read and long-read preprocessing, you can optionally save the resulting processed reads with `--save_preprocessed_reads`.

#### Redundancy Estimation

Metagenome 'coverage' or sequencing complexity estimations of short-read datasets can be activated with `--perform_shortread_redundancyestimation`.

This turns on checking of read redundancy in a sequencing library using [Nonpareil](https://nonpareil.readthedocs.io/en/latest/), to provide an estimation of whether you have sequenced enough to capture all possible genomes present in your metagenomic sample (with the assumption that once you've sequenced enough, you will keep sequencing PCR amplicons rather than unique reads).

This is only suitable for short-read, and in zellerlab/flexprofiler specifically, FASTQ files.

Nonpareil is performed on processed reads (i.e. after fastp or AdapterRemoval). This will run on either the first read of each read pair (as recommended by the authors), or on merged reads.

Before using this tool please note the following caveats:

:::warning

- It is not recommended to keep unmerged (`--shortread_qc_includeunmerged`) reads when using the calculation.
- Your shortest reads _after_ processing should not go below 24bp
  If the 'kmer' value is not correct, make sure your shortest reads _after_ processing is not less than 24bp.

      If this is the case you will need to specify in a custom config

      ```nextflow
      process {
        withName: NONPAREIL_NONPAREIL {
          ext.args = { "-k <NUMBER>" }
          }
      }
      ```

      Where `<NUMBER>` should be at least the shortest read in your library

  :::

#### Complexity Filtering

Complexity filtering can be activated via the `--perform_shortread_complexityfilter` flag.

Complexity filtering is primarily a run-time optimisation step.
It is not necessary for accurate taxonomic classification/profiling, however it can speed up run-time of each tool by removing reads with low-diversity of nucleotides (e.g. with mono-nucleotide - `AAAAAAAA`, or di-nucleotide repeats `GAGAGAGAGAGAGAG`) that have a low-chance of giving an informative taxonomic ID as they can be associated with many different taxa.
Removing these reads therefore saves computational time and resources.

There are currently three options for short-read complexity filtering: [`bbduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/), [`prinseq++`](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus), and [`fastp`](https://github.com/OpenGene/fastp#low-complexity-filter).

There are two options for long-read quality filtering: [`Filtlong`](https://github.com/rrwick/Filtlong) and [`nanoq`](https://github.com/esteinig/nanoq), with `nanoq` being the default option.

The tools offer different algorithms and parameters for removing low complexity reads and quality filtering.
We therefore recommend reviewing the pipeline's [parameter documentation](docs/parameters.md) and the documentation of the tools (see links above) to decide on optimal methods and parameters for your dataset.

You can optionally save the FASTQ output of the run merging with the `--save_complexityfiltered_reads`.
If running with `fastp`, complexity filtering happens inclusively within the earlier shortread preprocessing step.
Therefore there will not be an independent pipeline step for complexity filtering, and no independent FASTQ file (i.e. `--save_complexityfiltered_reads` will be ignored) - your complexity filtered reads will also be in the `fastp/` folder in the same file(s) as the preprocessed read.

:::warning
For nanopore data: we do not recommend performing any read preprocessing or complexity filtering if you are using ONTs Guppy toolkit for basecalling and post-processing.
:::

#### Host-Read Removal

Removal of possible-host reads from FASTQ files prior classification/profiling can be activated with `--perform_shortread_hostremoval` or `--perform_longread_hostremoval`.

Similarly to complexity filtering, host-removal can be useful for runtime optimisation and reduction in misclassified reads.
It is not always necessary to report classification of reads from a host when you already know the host of the sample, therefore you can gain a run-time and computational advantage by removing these prior typically resource-heavy classification/profiling with more efficient methods.
Furthermore, particularly with human samples, you can reduce the number of false positives during classification/profiling that occur due to host-sequence contamination in reference genomes on public databases.

zellerlab/flexprofiler currently offers host-removal via alignment against a reference genome with Bowtie2 for short reads and minimap2 for long reads, and the use of the unaligned reads for downstream classification/profiling.

You can supply your reference genome in FASTA format with `--hostremoval_reference`.
You can also optionally supply a directory containing pre-indexed Bowtie2 index files with `--shortread_hostremoval_index` or a minimap2 `.mmi` file for `--longread_hostremoval_index`, however zellerlab/flexprofiler will generate these for you if necessary.
Pre-supplying the index directory or files can greatly speed up the process, and these can be re-used.

:::tip
If you have multiple taxa or sequences you wish to remove (e.g., the host genome and then also PhiX - common quality-control reagent during sequencing) you can simply concatenate the FASTAs of each taxa or sequences into a single reference file.
:::

#### Run Merging

For samples that may have been sequenced over multiple runs, or for FASTQ files split into multiple chunks, you can activate the ability to merge across all runs or chunks with `--perform_runmerging`.
For more information how to set up your input samplesheet, see [Multiple runs of the same sample](#multiple-runs-of-the-same-sample).
Activating this functionality will concatenate the FASTQ files with the same sample name _after_ the optional preprocessing steps and _before_ classification/profiling.
Note that libraries with runs of different pairing types will **not** be merged and this will be indicated on output files with a `_se` or `_pe` suffix to the sample name accordingly.

You can optionally save the FASTQ output of the run merging with the `--save_runmerged_reads`.

#### Classification and Profiling

The following sections provide tips and suggestions for running the different classification and profiling tools _within the pipeline_.
For advice and/or guidance whether you should run a particular tool on your specific data, please see the documentation of each tool!
An important distinction between the different tools in included in the pipeline is classification versus profiling.
For advice as to which tool to run in your context, please see the documentation of each tool.

<!-- TODO: add section
##### MapSeq

##### IDTaxa

##### mOTUs v3
-->

##### mOTUs v4

mOTUs v4 has a number of breaking changes compared to v3, in particular:
- The taxonomy is based on GDBM and not NCBI
- Many command line options have been removed

To run mOTUs v4, you need to specify the parameter `--run_motus4`.
You also need to provide appropriate databases (see the [database documentation](databases.md) for details).
The parameter `--motus4_rename_reads` is useful when working with paired-end reads with a `/1` and `/2` suffix in the read names.
mOTUs v4 expects the reads in a pair to have identical names, so this steps renames all the reads accordingly removing the suffixes.

##### Cayman

To run mOTUs v4, you need to specify the parameter `--run_cayman`.
You also need to provide appropriate databases (see the [database documentation](databases.md) for details).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version.
When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since.
To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull zellerlab/flexprofiler
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data.
This ensures that a specific version of the pipeline code and software are used when you run your pipeline.
If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [zellerlab/flexprofiler releases page](https://github.com/zellerlab/flexprofiler/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.
Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.
For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time.
For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`.
This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test_<tool_name>`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline.
Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
For input to be considered the same, not only the names must be identical but the files' contents as well.
For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command).
See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests.
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original).
If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool.
By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects.
However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter.
You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session.
The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
