# zellerlab/flexprofiler

## Introduction

**zellerlab/flexprofiler** is a bioinformatics best-practice analysis pipeline for (microbiome) taxonomic classification and profiling of shotgun and 16S metagenomic data. It allows for in-parallel taxonomic identification of reads or taxonomic abundance estimation with multiple classification and profiling tools against multiple databases, and produces standardised output tables for facilitating results comparison between different tools and databases.
It aims to enable microbiome profiling not only of single datasets but also of collections of studies for meta-analisis.

**zellerlab/flexprofiler** is an adaptation of [**nf-core/taxprofiler**](https://github.com/nf-core/taxprofiler) for the needs of the Zeller lab.
Many features of **nf-core/taxprofiler** that are not needed for the Zeller lab have been removed for simplicity, and support for 16S sequencing data and meta-analyses has been added.
This is not strictly a fork of **nf-core/taxprofiler**, as the two pipelines are planned to diverge consistently with their different intended use.
We plan to contribute to **nf-core** modules or other components of more general interest, and to implement relevent **nf-core** developments in this pipeline on an irregular schedule.

## Pipeline simplified overview

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Read pre-processing
   - Adapter clipping and merging: [fastp](https://github.com/OpenGene/fastp)
   - Host-read removal: [BowTie2](http://bowtie-bio.sourceforge.net/bowtie2/)
3. Performs taxonomic classification and/or profiling using one or more of:
   - [mOTUs](https://motus-tool.org/)
4. Standardises output tables (tool-specific and [`Taxpasta`](https://taxpasta.readthedocs.io))
5. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

```csv title="samplesheet.csv"
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
2612,run1,ILLUMINA,2612_run1_R1.fq.gz,,
2612,run2,ILLUMINA,2612_run2_R1.fq.gz,,
2612,run3,ILLUMINA,2612_run3_R1.fq.gz,2612_run3_R2.fq.gz,
```

Each row represents a fastq file (single-end), a pair of fastq files (paired end), or a fasta (with long reads).

Additionally, you will need a database sheet that looks as follows:

```csv title="databases.csv"
tool,db_name,db_params,db_path
kraken2,db2,--quick,/<path>/<to>/kraken2/testdb-kraken2.tar.gz
metaphlan,db1,,/<path>/<to>/metaphlan/metaphlan_database/
```

That includes directories or `.tar.gz` archives containing databases for the tools you wish to run the pipeline against.

Now, you can run the pipeline using:

```bash
nextflow run zellerlab/flexprofiler \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --databases databases.csv \
   --outdir <OUTDIR>  \
   --run_motus
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](docs/usage.md) and the [parameter documentation](docs/parameters.md).


### Database download

To use **zellerlab/flexprofiler** you need the databases that the profiling tools that you want to run require.
These can be downloaded using [this script](utils/fetch_databases.sh), which requires as argument the a path where the databases will be downloaded.
The script also generates a `flexprofiler_databases.csv` file in the local directory, which can be used in input for the pipeline.

```bash
fetch_databases.sh <out_database_path>
```

> [!IMPORTANT]
> If you are a member of the Zeller lab, most likely you do not need to download the databases as this has already been done for you in a centralised location.
> Databases have been already set up for the following compute environments with corresponding profiles:
> - LUMC SHARK cluster (nextflow profile name: `zellerlab_shark`)
>
> You just need to use the correct profile in the nextflow run so that they are loaded correctly, for example:
>```bash
>nextflow run zellerlab/flexprofiler -profile zellerlab_shark <YOUR OTHER NEXTFLOW ARGUMENTS>
>```


## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/taxprofiler/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/taxprofiler/output).

## Credits

**zellerlab/flexprofiler** is developped and maintained by Saul Pierotti.
**nf-core/taxprofiler** was originally written by James A. Fellows Yates, Sofia Stamouli, Moritz E. Beber, Lili Andersson-Li, and the nf-core/taxprofiler team.
We refer the reader to [the pipeline page](https://github.com/nf-core/taxprofiler) for further credits.

### Team

- [Saul Pierotti](https://github.com/saulpierotti)

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
