process CAYMAN_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cayman:0.10.2--pyh7e72e81_0'
        : 'biocontainers/cayman:0.10.2--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(bam)
    path db

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
    def annot_file = task.ext.annot_file ?: "${db}/*.csv"
    """
    # alignment is done separately by BWA_MEM, cayman_patch.py makes cayman read
    # the alignment file instead of running bwa mem (the bwa index is not used)
    cayman_patch.py profile \\
        --singles $bam \\
        --out_prefix $prefix \\
        $annot_file \\
        no_bwa_index \\
        $args \\
        &> ${prefix}.cayman.log

    VERSION=\$(cayman --version|sed 's/cayman //')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cayman: \$VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.aln_stats.txt.gz"
    touch "${prefix}.cazy.txt.gz"
    touch "${prefix}.gene_counts.txt.gz"
    touch "${prefix}.cayman.log"

    VERSION=\$(cayman --version|sed 's/cayman //')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cayman: \$VERSION
    END_VERSIONS
    """

}
