process IDTAXA_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84e799e45917de0293324e738f0c0b7d1fd0f85ac0536d39d1161d7ef372fbaa/data'
        : 'community.wave.seqera.io/library/pandas:2.2.1--2558de69969e58b9'}"

    input:
    // the profiles of all the samples for one training set (meta.id), staged in a folder, so that the command lists
    // them with a glob, and not one by one
    tuple val(meta), path(input, stageAs: 'profiles/*')

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the comma-separated ranks to write a table for, all of them if not set
    def ranks = task.ext.ranks ?: ''
    // the profiles are named <sample>_<training set>.idtaxa.profile.tsv
    """
    # the list of the profiles is written by printf, a shell builtin, so that it is not limited by the length of a
    # command line
    printf '%s\\n' profiles/* > profiles.txt
    idtaxa_merge.py ${prefix} _${meta.id}.idtaxa.profile.tsv profiles.txt ${ranks}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lineage.tsv
    """
}
