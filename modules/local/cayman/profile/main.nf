process CAYMAN_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gqlib:2.14.4--pyhdfd78af_0'
        : 'biocontainers/gqlib:2.14.4--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(bam)
    path annot

    output:
    tuple val(meta), path("*.cazy.txt.gz"       ), emit: out
    tuple val(meta), path("*.aln_stats.txt.gz"  ), emit: aln_stats
    tuple val(meta), path("*.gene_counts.txt.gz"), emit: gene_counts
    tuple val(meta), path("*.cayman.log"        ), emit: log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # `cayman profile` runs bwa mem itself, then the gqlib domain profiler: the alignment is done
    # separately by BWA_ALIGN, and the gqlib domain profiler is run directly on it, as cayman 0.10.2
    # does. The cayman annotation csv is in the hmmer format.
    python -m gqlib \\
        --mode domain \\
        --db_format hmmer \\
        --bam $bam \\
        --db $annot \\
        --out_prefix $prefix \\
        $args \\
        &> ${prefix}.cayman.log

    # gqlib names the category of the annotation `feature`, cayman names it `cazy`
    mv ${prefix}.feature.txt.gz ${prefix}.cazy.txt.gz

    VERSION=\$(python -c 'import gqlib; print(gqlib.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gqlib: \$VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.aln_stats.txt.gz"
    touch "${prefix}.cazy.txt.gz"
    touch "${prefix}.gene_counts.txt.gz"
    touch "${prefix}.cayman.log"

    VERSION=\$(python -c 'import gqlib; print(gqlib.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gqlib: \$VERSION
    END_VERSIONS
    """

}
