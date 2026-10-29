process PIMENTO_STD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mi-pimento:1.3.3--pyhdfd78af_0'
        : 'biocontainers/mi-pimento:1.3.3--pyhdfd78af_0'}"

    input:
    // one mate of paired-end reads, or single-end reads (searched with --merged, set in ext.args)
    tuple val(meta), path(reads)

    output:
    // empty if no primer was found
    tuple val(meta), path("*_std_primers.fasta"), emit: primers
    tuple val(meta), path("*_std_primer_out.txt"), emit: out
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pimento std -i ${reads} -o ${prefix} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pimento: \$(pimento --version 2>&1 | sed -n 's/.*version //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_std_primers.fasta ${prefix}_std_primer_out.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pimento: \$(pimento --version 2>&1 | sed -n 's/.*version //p')
    END_VERSIONS
    """
}
