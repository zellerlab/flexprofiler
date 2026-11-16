process READ_COUNTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bdb4d8c79866b73925b2a5b70b17fdfb04d873b078eb1f030c987c4ef932e035/data'
        : 'community.wave.seqera.io/library/pandas_matplotlib:2ef9b2c7a64400ba'}"

    input:
    // named uniquely, as in the manifest
    path files, stageAs: 'inputs/*'
    // with a header, one row per file: file sample run_accession kind db_name
    path manifest

    output:
    path "read_counts.tsv"    , emit: tsv
    path "read_counts{_mqc.tsv,_plot_mqc.json}", emit: mqc
    path "read_counts.{png,pdf}", emit: plot
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions

    script:
    """
    read_counts.py \\
        ${manifest} \\
        inputs \\
        read_counts.tsv \\
        read_counts_mqc.tsv \\
        read_counts_plot_mqc.json \\
        read_counts
    """

    stub:
    """
    touch read_counts.tsv
    touch read_counts_mqc.tsv
    touch read_counts_plot_mqc.json
    touch read_counts.png
    touch read_counts.pdf
    """
}
