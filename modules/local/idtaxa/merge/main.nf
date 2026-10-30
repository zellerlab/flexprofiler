process IDTAXA_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'biocontainers/pandas:2.2.1'}"

    input:
    // the profiles of all the samples for one training set (meta.id), staged in a folder, so that the command lists
    // them with a glob, and not one by one
    tuple val(meta), path(input, stageAs: 'profiles/*')

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lineage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
