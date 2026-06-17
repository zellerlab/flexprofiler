# zellerlab/flexprofiler: Internal zellerlab pre-configured options and scripts

## Running the pipeline on the SHARK cluster
To run zellerlab/flexprofiler on the SHARK cluster, copy [this example run file](../utils/run_shark.sh) and edit as required.
You can then run it with `sbatch run_shark.sh`, or start an interactive shell on a worker (not login!) node and just execute it with `./run_shark.sh`

## Params file
To run the pipeline you need to set up a `params.yaml` file with the required parameters.
An example file that you can copy and edit to suit your needs is [here](../utils/params.yaml).

## Profiles
A set of Nextflow profiles has been developped to standardise microbiome profiling in the Zellerlab.
These profiles are defined [here](../conf/zellerlab.config).

You can load them by appending `-p <profile name>` to your nextflow command.
Multiple profiles can be specified with a comma (`-p <profle1>,<profile2>`).

For reproducibility, always use singularity to run the pipeline (by specifying the nextflow profile `-p singularity`).
Note that you need to have singularity installed for this to work.

The following profiles are available:
- `zellerlab_shark`: for use on the LUMC shark cluster, sets paths to databases etc.
- `zellerlab_metagenomics`: for short read metagenomics profiling, sets defaults for preprocessing and collation
- `zellerlab_minimal`: skip all non-essential QC steps, useful for large meta-analysis runs

## Containers and environments
All the pipeline steps are set up to use appropriate Singularity containers, so you need to have the singularity profile enabled (see above) and the singularity program should be available in your environment (for the SHARK cluster if you use [this script](../utils/run_shark.sh) the singularity module is automatically loaded so you don't need to do anything).

For some steps Conda/mamba/micromamba is also required.
This need to be installed in your system and available in the `$PATH`.
If you use micromamba or mamba instead than conda (as you should), you need to set the relevant flag in `~/.nextflow/config`.
Copy [this config file](../utils/user.config) to `~/.nextflow/config` and adapt it as required to do it.

It is good also if you set an appropriate cache directory for singularity and conda so that containers and environments are not downloaded again every time you run the pipeline.
If you use [this config file](../utils/user.config) this is also set up for you, you only need to set the right paths.

## Nextflow Tower / Seqera cloud
To monitor more easily the progress of your run use [Sequera Cloud](https://cloud.seqera.io), which can be reached also outside of the LUMC VPN.
Create an account on the website and generate a Token.
Then, copy the content of [this file](../utils/user.config) to `~/.nextflow/config` (generate it if not present, this is loaded automatically by nextflow).
Make sure to copy your own token generated in Seqera Cloud to the config!

## Database generation
The required databases have been already generated for you and should be automatically loaded on the SHARK cluster.
If not, contact Saul.