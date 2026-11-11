process GQLIB_MERGE {
    tag "${meta.id}|${meta.category}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b3/b36790ad56052b28ef374a2d230862b0d7dd9d0d599b820d2e6bca34d0be57c4/data'
        : 'community.wave.seqera.io/library/gqlib:2.14.4--1a47e08f7db4071d'}"

    input:
    // the profiles of all the samples for one database and annotation (meta.id) and one category (meta.category), staged
    // in a folder, so that the command lists them with a glob, and not one by one
    tuple val(meta), path(input, stageAs: 'profiles/*')

    output:
    tuple val(meta), path("*.txt.gz"), emit: txt
    tuple val("${task.process}"), val('gqlib'), eval('python -c "import gqlib; print(gqlib.__version__)"'), emit: versions_gqlib, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the profiles are named <sample>_<database>_<annotation>.<category>.txt.gz, meta.id being <database>_<annotation>;
    // one table per count column (uniq_raw, combined_raw, ...)
    """
    # the list of the profiles is written by printf, a shell builtin, so that it is not limited by the length of a
    # command line
    printf '%s\\n' profiles/* > profiles.txt
    gqlib_merge.py ${prefix}.${meta.category} _${meta.id}.${meta.category}.txt.gz profiles.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.${meta.category}.uniq_raw.txt.gz"
    """

}
