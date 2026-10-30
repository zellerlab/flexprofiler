process READ_COUNTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'biocontainers/pandas:2.2.1'}"

    input:
    // each file in its own folder, inputs/<index>/, as files of different steps can have the same name
    path files, stageAs: 'inputs/*/*'
    // with a header, one row per file, in the order of the files: sample run_accession kind db_name annot_db_name
    path manifest

    output:
    path "read_counts.tsv"    , emit: tsv
    path "read_counts_mqc.tsv", emit: mqc
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    read_counts.py \\
        ${manifest} \\
        inputs \\
        read_counts.tsv \\
        read_counts_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch read_counts.tsv
    touch read_counts_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
