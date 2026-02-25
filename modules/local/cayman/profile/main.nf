process CAYMAN_PROFILE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cayman:0.10.2--pyh7e72e81_0'
        : 'biocontainers/cayman:0.10.2--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(annot), path(gmgc10_fa), path(gmgc10_amb), path(gmgc10_ann), path(gmgc10_bwt), path(gmgc10_pac), path(gmgc10_sa)

    output:
    tuple val(meta), path("${meta.id}.cazy.txt.gz"       ), emit: out
    tuple val(meta), path("${meta.id}.aln_stats.txt.gz"  ), emit: aln_stats
    tuple val(meta), path("${meta.id}.gene_counts.txt.gz"), emit: gene_counts
    tuple val(meta), path("${meta.id}.cayman.log"        ), emit: log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = meta.single_end ? "-singles $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    cayman profile \\
        $inputs \\
        --out_prefix $prefix \\
        $annot \\
        $gmgc10_fa \\
        $args \\
        &> ${prefix}.cayman.log
        
    VERSION=\$(cayman --version|sed 's/cayman //')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cayman: \$VERSION
    END_VERSIONS
    """

    stub:
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
