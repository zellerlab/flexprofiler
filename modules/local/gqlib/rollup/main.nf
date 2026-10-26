process GQLIB_ROLLUP {
    tag "$meta.id"
    label 'process_single'

    // gqlib_rollup.py needs only python, the same container of GQLIB_COUNT is used
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cayman:0.10.2--pyh7e72e81_0'
        : 'biocontainers/cayman:0.10.2--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(gene_counts), path(read_counts)
    path db

    output:
    tuple val(meta), path("*.txt.gz"), emit: out
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def annot_file = task.ext.annot_file ?: "${db}/*.annotation.tsv*"
    """
    gqlib_rollup.py \\
        --gene_counts $gene_counts \\
        --read_counts $read_counts \\
        --annotation $annot_file \\
        --out_prefix $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.category.txt.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

}
