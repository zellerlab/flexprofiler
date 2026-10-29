process READ_COUNTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'biocontainers/python:3.9--1'}"

    input:
    // the labels of each file ([sample, run, kind, db_name, annot_db_name]) and the files, in the same order. Each file
    // is staged in its own folder, inputs/<index>/, as files of different steps can have the same name
    tuple val(labels), path(files, stageAs: 'inputs/*/*')

    output:
    path "read_counts.tsv"    , emit: tsv
    path "read_counts_mqc.tsv", emit: mqc
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // one row per file, single-quoted for the shell, on one line (a multi-line value breaks the indentation of the script)
    def manifest = labels.withIndex().collect { label, i -> "'" + ([i + 1] + label.values()).join('\t').replace("'", "'\\''") + "'" }.join(' ')
    """
    printf '%s\\n' ${manifest} > manifest.tsv

    read_counts.py \\
        manifest.tsv \\
        inputs \\
        read_counts.tsv \\
        read_counts_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch read_counts.tsv
    touch read_counts_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
