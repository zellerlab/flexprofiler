# zellerlab/flexprofiler: Internal zellerlab pre-configured options
A set of Nextflow profiles has been developped to standardise microbiome profiling in the Zellerlab.
These profiles are defined (here)[conf/zellerlab.config].

You can load them by appending `-p <profile name>` to your nextflow command.
Multiple profiles can be specified with a comma (`-p <profle1>,<profile2>`).

For reproducibility, always use singularity to run the pipeline (by specifying the nextflow profile `-p singularity`).
Note that you need to have singularity installed for this to work.

The following profiles are available:
- `zellerlab_shark`: for use on the LUMC shark cluster, sets paths to databases etc.
- `zellerlab_metagenomics`: for short read metagenomics profiling, sets defaults for preprocessing and collation