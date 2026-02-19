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

## Nextflow Tower / Seqera cloud
To monitor more easily the progress of your run use [Sequera Cloud](https://cloud.seqera.io), which can be reached also outside of the LUMC VPN.
Create an account on the website and generate a Token.
Then, copy the content of [this file](../utils/tower.config) to `~/.nextflow/config` (generate it if not present, this is loaded automatically by nextflow).
Make sure to copy your own token generated in Seqera Cloud to the config!

## Database generation
The required databases have been already generated for you and should be automatically loaded on the SHARK cluster.
If not, contact Saul.