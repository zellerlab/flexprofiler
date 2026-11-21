# zellerlab/flexprofiler: Internal zellerlab pre-configured options and scripts

## Running the pipeline on the SHARK cluster
To run zellerlab/flexprofiler on the SHARK cluster, copy [this example run file](../utils/run_shark.sh) and edit as required.
You can then run it with `sbatch run_shark.sh`, or start an interactive shell on a worker (not login!) node and just execute it with `./run_shark.sh`.
To stop a run, use `scancel --batch --signal=TERM <job id>`: Nextflow (which the script runs with `exec`) then cancels its tasks before exiting. A plain `scancel <job id>` kills Nextflow before it can, leaving its tasks (the `nf-*` jobs) running and queued: check afterwards with `squeue -u $USER` that none of them is left, and cancel those that are

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
- `zellerlab_shark`: for use on the LUMC shark cluster, sets paths to databases etc., and runs the containers with apptainer (no need for `-profile apptainer`), with the images in the cache shared by the lab (`/exports/lucid-grpzeller-work/SOFTWARE_CACHE/flexprofiler/apptainer`). The PathSeq-T2T QC filter (`PATHSEQT2T_QCFILTER`) runs in the node-local `/tmp` (Nextflow `scratch`), many times faster than in the NFS work directory, as its Spark temporary files are written to its task directory
- `zellerlab_metagenomics`: for short read metagenomics profiling, sets defaults for preprocessing and collation
- `zellerlab_amplicon`: for 16S amplicon profiling with IDTAXA, sets defaults for preprocessing (with primer trimming), run merging, collation and IDTAXA (turn it on with `--run_idtaxa`): the pairs not merged are joined with Ns (`--amplicon_pair_handling merged_plus_joined`), and classified against the `gtdb_r220_motus4.0` training set of the database CSV (GTDB R220, the taxonomy of mOTUs 4.0 used by `--run_motus4`, so that 16S and shotgun profiles share the taxon names; change with `--idtaxa_db_names`), with a merged table for the genus only besides the one of the whole lineages (`--idtaxa_ranks genus`)
- `zellerlab_pathseq_t2t`: for low-biomass human samples (e.g. tumour WGS), the [PathSeq-T2T host filtering](pathseq_t2t.md) with Kraken2, run merging and collation; use it instead of `zellerlab_metagenomics`. On SHARK, `zellerlab_shark` sets the paths of its references (downloaded by `utils/fetch_missing_databases_shark.sh`); on other clusters give them with `--pathseq_t2t_grch38_bwa_index`, `--pathseq_t2t_chm13_bwa_index` and `--pathseq_host_dir`
- `zellerlab_pathseq_t2t_rna`: as `zellerlab_pathseq_t2t`, for RNA-seq: the reads are aligned to GRCh38 with STAR (`--pathseq_t2t_grch38_aligner star`), which also counts the reads per host gene (see [RNA-seq](pathseq_t2t.md#rna-seq-star)). On SHARK, `zellerlab_shark` sets the STAR index of GRCh38 (GRCh38.d1.vd1 with GENCODE v36), the T2T-CHM13 index and the PathSeq host files
- `zellerlab_minimal`: skip all non-essential QC steps, useful for large meta-analysis runs
- `zellerlab_large_scale`: for large datasets, processed one batch at a time, see [below](#large-datasets-processed-in-batches)
- `zellerlab_alice`: for the ALICE cluster (Leiden University), sets the paths to the databases and runs the containers with apptainer. The databases are copies of those of SHARK, made with `utils/sync_databases_alice.sh` (run on SHARK: `sbatch utils/sync_databases_alice.sh [<database directory> ...]`), which also writes the database sheet with only the databases copied so far. To run on ALICE with the inputs and results on SHARK, use [this script](../utils/run_alice.sh): set it up as `run_shark.sh` (it takes the same params file, with SHARK paths), and submit it on ALICE, or from SHARK with `ssh alice 'cd /projects/pi-zeller/<YOUR USER NAME> && sbatch' < run_alice.sh`. It copies the inputs from SHARK, runs the pipeline and copies the results back to `outdir` on SHARK. It deletes nothing on ALICE: when the run is done, check with `utils/check_transfer_alice.sh <RUN_DIR>` that the results are all on SHARK, then delete the run directory. It runs the Nextflow launcher shared by the lab (`/projects/pi-zeller/SHARED/flexprofiler/nextflow`), which downloads the version the pipeline requires the first time: the Nextflow module of ALICE is too old

## Large datasets processed in batches
A run of thousands of samples needs a work directory as large as all the intermediate files of all the samples (the trimmed and host-filtered reads, the alignments), kept until the end of the run. The profile `zellerlab_large_scale` (`--large_scale`) keeps it to one batch of samples (e.g. one study) at a time:

1. Add a `processing_batch` column to the samplesheet, with the batch of each run (all the runs of a sample in the same batch; letters, digits, `.`, `_` and `-`). The column is not part of the meta map of the runs, so adding it to a samplesheet that was already run does not change what `-resume` reuses.
2. Add `zellerlab_large_scale` to the profiles (e.g. `-profile zellerlab_shark,zellerlab_metagenomics,zellerlab_large_scale`), and run the pipeline as usual, with `-resume`, again and again with the same command, until it reports that nothing is left to do. Each run:
   - deletes the task directories of the batches done before (once per batch, `--large_scale_cleanup`), from the traces of the runs that processed them (the cached tasks too, e.g. those of an earlier run of the whole samplesheet in the same work directory). The task directories of the last batch are deleted by the run after it.
   - processes the next batch not done (the one with the fewest runs; `--processing_batch <batch>` to choose), and publishes the profiles of its samples, with their manifests (`aggregation_inputs/`, see [Aggregation inputs](output.md#aggregation-inputs)), but no table across the samples. If it succeeds, the batch is recorded as done (`aggregation_inputs/batches/<batch>.json`); if not, the next run processes it again, resuming it.
   - when all the batches are done, processes no sample, and makes the tables of all the samples instead, as a run of all of them would: the merged profiles (`<tool>/merge/...`) and `read_counts/read_counts.tsv`. Then `aggregation_inputs/aggregated.json` is written.

   All the runs must use the same `--outdir` (an error otherwise) and the same profilers. If the samples of a batch change in the samplesheet, the batch is processed again. The runs are on the same launch directory, one at a time.

   On SHARK, copy [this script](../utils/run_shark_large_scale.sh), set it up as `run_shark.sh`, and submit it with `sbatch run_shark_large_scale.sh`. It runs the pipeline once; `zellerlab_shark` sets `--large_scale_slurm_resubmit` (which does nothing without `--large_scale`), so each run that processes a batch submits the script of its job again (found with `scontrol show job`) (`sbatch --dependency=afterok:<its job> --kill-on-invalid-dep=yes`), so the next run starts once this one has ended successfully, until the tables are made. A run that fails ends the chain; submitted again, it resumes. Stop it with `scancel --batch --signal=TERM <job id>` of the running job: Nextflow cancels its tasks, and the next job, if already submitted, is removed by SLURM. Options for `sbatch` go in the `#SBATCH` lines of the script, as those of the command line are not passed on. No next run is submitted when `--processing_batch` is given, as each run would process that batch again. Without SLURM, run the same command again until the pipeline reports that nothing is left to do.

With `zellerlab_shark`, a `zellerlab_large_scale` run must be in a SLURM job submitted with `sbatch <script>` (`--large_scale_slurm_resubmit`; add `--large_scale_slurm_resubmit false` to run the batches by hand). `zellerlab_large_scale` turns on `--run_profile_standardisation` (with `--large_scale` it is an error to turn it off, as the run after the last batch always merges the profiles), skips the QC reports and MultiQC (as `zellerlab_minimal`), and adds the work directory of each task to the trace, which the cleanup reads. The cleanup keeps the tasks that prepare a database or reference, shared by the batches: UNTAR, and the indexes built when not given (`BOWTIE2_BUILD`, `BWA_INDEX_GRCH38`, `BWA_INDEX_T2T`). A database given as a `.tar.gz` is untarred again in each batch (its UNTAR task is not deleted, but cached only within the batch): give the databases untarred.

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