process GQLIB_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gqlib:2.14.4--pyhdfd78af_0'
        : 'biocontainers/gqlib:2.14.4--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.gene_counts.txt.gz"), emit: gene_counts
    tuple val(meta), path("*.read_counts.txt"   ), emit: read_counts
    tuple val(meta), path("*.aln_stats.txt.gz"  ), emit: aln_stats
    tuple val(meta), path("*.gqlib.log"         ), emit: log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # gqlib_patch.py runs only the gene counting, the functional roll-up is done by GQLIB_ROLLUP
    gqlib_patch.py \\
        --mode genes \\
        --bam $bam \\
        --out_prefix $prefix \\
        $args \\
        &> ${prefix}.gqlib.log

    VERSION=\$(python -c 'import gqlib; print(gqlib.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gqlib: \$VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.gene_counts.txt.gz"
    touch "${prefix}.read_counts.txt"
    touch "${prefix}.aln_stats.txt.gz"
    touch "${prefix}.gqlib.log"

    VERSION=\$(python -c 'import gqlib; print(gqlib.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gqlib: \$VERSION
    END_VERSIONS
    """

}
