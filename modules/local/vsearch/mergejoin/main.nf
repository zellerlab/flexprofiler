process VSEARCH_MERGEJOIN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vsearch:2.32.0--hd2be7a0_0'
        : 'biocontainers/vsearch:2.32.0--hd2be7a0_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.mergejoin.fastq.gz"), emit: reads
    tuple val(meta), path("*.vsearch.log")       , emit: log
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''   // vsearch --fastq_mergepairs
    def args2 = task.ext.args2 ?: '' // vsearch --fastq_join
    def prefix = task.ext.prefix ?: "${meta.id}"
    // One sequence per read pair: the pairs that overlap are merged, the others are joined, R1 followed by a padding
    // of Ns and the reverse complement of R2
    """
    vsearch \\
        --fastq_mergepairs ${reads[0]} \\
        --reverse ${reads[1]} \\
        --fastqout merged.fastq \\
        --fastqout_notmerged_fwd notmerged_1.fastq \\
        --fastqout_notmerged_rev notmerged_2.fastq \\
        --threads $task.cpus \\
        $args \\
        2> ${prefix}.vsearch.log

    vsearch \\
        --fastq_join notmerged_1.fastq \\
        --reverse notmerged_2.fastq \\
        --fastqout joined.fastq \\
        $args2 \\
        2>> ${prefix}.vsearch.log

    cat merged.fastq joined.fastq | gzip -n > ${prefix}.mergejoin.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | sed -n '1s/^vsearch v\\([0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.mergejoin.fastq.gz
    touch ${prefix}.vsearch.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | sed -n '1s/^vsearch v\\([0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """
}
