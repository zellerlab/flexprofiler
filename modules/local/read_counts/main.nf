process READ_COUNTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84e799e45917de0293324e738f0c0b7d1fd0f85ac0536d39d1161d7ef372fbaa/data'
        : 'community.wave.seqera.io/library/pandas:2.2.1--2558de69969e58b9'}"

    input:
    // named uniquely, as in the manifest
    path files, stageAs: 'inputs/*'
    // with a header, one row per file: file sample run_accession kind db_name
    path manifest

    output:
    path "read_counts.tsv"    , emit: tsv
    path "read_counts_mqc.tsv", emit: mqc
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions

    script:
    """
    read_counts.py \\
        ${manifest} \\
        inputs \\
        read_counts.tsv \\
        read_counts_mqc.tsv
    """

    stub:
    """
    touch read_counts.tsv
    touch read_counts_mqc.tsv
    """
}
