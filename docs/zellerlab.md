# zellerlab/flexprofiler: Internal zellerlab pre-configured options and scripts

## Running the pipeline on the SHARK cluster
To run zellerlab/flexprofiler on the SHARK cluster, copy [this example run file](../utils/run_shark.sh) and edit as required.
You can then run it with `sbatch run_shark.sh`, or start an interactive shell on a worker (not login!) node and just execute it with `./run_shark.sh`

## Params file
To run the pipeline you need to set up a `params.yaml` file with the required parameters.
An example file that you can copy and edit to suit your needs is [here](../utils/params.yaml).

## Profiles
A set of Nextflow profiles has been developed to standardise microbiome profiling in the Zellerlab.
These profiles are defined [here](../conf/zellerlab.config).

You can load them by appending `-profile <profile name>` to your nextflow command.
Multiple profiles can be specified with a comma (`-profile <profile1>,<profile2>`).

For reproducibility, always use apptainer to run the pipeline (by specifying the nextflow profile `-profile apptainer`; `zellerlab_shark` sets it already).
Note that you need to have apptainer installed for this to work (`-profile singularity` also works, with Singularity).

The following profiles are available:
- `zellerlab_shark`: for use on the LUMC shark cluster, sets paths to databases etc., and runs the containers with apptainer (no need for `-profile apptainer`), with the images in the cache shared by the lab (`/exports/lucid-grpzeller-work/SOFTWARE_CACHE/flexprofiler/apptainer`)
- `zellerlab_metagenomics`: for short read metagenomics profiling, sets defaults for preprocessing and collation
- `zellerlab_amplicon`: for 16S amplicon profiling with IDTAXA, sets defaults for preprocessing (with primer trimming), run merging, collation and IDTAXA (turn it on with `--run_idtaxa`): the pairs not merged are joined with Ns (`--amplicon_pair_handling merged_plus_joined`), and classified against the `gtdb_r220_motus4.0` training set of the database CSV (GTDB R220, the taxonomy of mOTUs 4.0 used by `--run_motus4`, so that 16S and shotgun profiles share the taxon names; change with `--idtaxa_db_names`), with a merged table for the genus only besides the one of the whole lineages (`--idtaxa_ranks genus`)
- `zellerlab_pathseq_t2t`: for low-biomass human samples (e.g. tumour WGS), the [PathSeq-T2T host filtering](pathseq_t2t.md) with Kraken2, run merging and collation; use it instead of `zellerlab_metagenomics`. On SHARK, `zellerlab_shark` sets the paths of its references (downloaded by `utils/fetch_missing_databases_shark.sh`); on other clusters give them with `--pathseq_t2t_grch38_bwa_index`, `--pathseq_t2t_chm13_bwa_index` and `--pathseq_host_dir`
- `zellerlab_pathseq_t2t_rna`: as `zellerlab_pathseq_t2t`, for RNA-seq: the reads are aligned to GRCh38 with STAR (`--pathseq_t2t_grch38_aligner star`), which also counts the reads per host gene (see [RNA-seq](pathseq_t2t.md#rna-seq-star)). On SHARK, `zellerlab_shark` sets the T2T-CHM13 index and the PathSeq host files, but not the STAR index of GRCh38, which is not on the cluster yet: give it with `--pathseq_t2t_grch38_star_index`, or the FASTA and GTF to build it with `--pathseq_t2t_grch38_reference` and `--pathseq_t2t_grch38_gtf`
- `zellerlab_minimal`: skip all non-essential QC steps, useful for large meta-analysis runs
- `zellerlab_alice`: **incomplete**, for testing only, on the ALICE cluster, sets paths to databases etc.: only the mOTUs v4 database is set, there is no reference to orient the long amplicon reads (`--amplicon_orient_reference`), and the container engine must be given (`-profile apptainer`)

## Containers and environments
All the pipeline steps are set up to use appropriate Apptainer containers, so you need to have the apptainer profile enabled (see above; `zellerlab_shark` enables it) and the apptainer program should be available in your environment (for the SHARK cluster if you use [this script](../utils/run_shark.sh) the apptainer module is automatically loaded so you don't need to do anything).
The newest Nextflow module of SHARK (25.04.6) is older than the version the pipeline requires: the script sets `NXF_VER`, so that the Nextflow launcher downloads the right version the first time (to `~/.nextflow/framework`).

The images are kept in the apptainer cache shared by the lab (see `zellerlab_shark`). To download all of them at once, e.g. after a new version of the pipeline, so that the runs do not have to, submit [this script](../utils/pull_images.sh) on SHARK (`sbatch utils/pull_images.sh`, with `VERSION` set as in `run_shark.sh`): it skips the images already in the cache.

All the steps run in containers, so Conda is not needed (`zellerlab_shark` turns it off).

Off SHARK, it is good to set an apptainer cache directory so that containers are not downloaded again every time you run the pipeline (on SHARK, `zellerlab_shark` already uses the apptainer cache shared by the lab).
If you use [this config file](../utils/user.config) this is also set up for you, you only need to set the right paths.

## Nextflow Tower / Seqera cloud
To monitor more easily the progress of your run use [Sequera Cloud](https://cloud.seqera.io), which can be reached also outside of the LUMC VPN.
Create an account on the website and generate a Token.
Then, copy the content of [this file](../utils/user.config) to `~/.nextflow/config` (generate it if not present, this is loaded automatically by nextflow).
Make sure to copy your own token generated in Seqera Cloud to the config!

## Database generation
The required databases have been already generated for you and should be automatically loaded on the SHARK cluster.
If not, contact Saul.