process GQLIB_MERGE {
    tag "${meta.id}|${meta.category}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gqlib:2.14.4--pyhdfd78af_0'
        : 'biocontainers/gqlib:2.14.4--pyhdfd78af_0'}"

    input:
    // the profiles of all the samples for one database and annotation (meta.id) and one category (meta.category), staged
    // in a folder, so that the command lists them with a glob, and not one by one
    tuple val(meta), path(input, stageAs: 'profiles/*')

    output:
    tuple val(meta), path("*.txt.gz"), emit: txt
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the profiles are named <sample>_<database>_<annotation>.<category>.txt.gz, meta.id being <database>_<annotation>;
    // one table per count column (uniq_raw, combined_raw, ...)
    """
    # the list of the profiles is written by printf, a shell builtin, so that it is not limited by the length of a
    # command line
    printf '%s\\n' profiles/* > profiles.txt
    gqlib_merge.py ${prefix}.${meta.category} _${meta.id}.${meta.category}.txt.gz profiles.txt

    VERSION=\$(python -c 'import gqlib; print(gqlib.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gqlib: \$VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.${meta.category}.uniq_raw.txt.gz"

    VERSION=\$(python -c 'import gqlib; print(gqlib.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gqlib: \$VERSION
    END_VERSIONS
    """

}
