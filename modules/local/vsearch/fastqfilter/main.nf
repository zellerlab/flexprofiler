process VSEARCH_FASTQFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vsearch:2.32.0--hd2be7a0_0'
        : 'biocontainers/vsearch:2.32.0--hd2be7a0_0'}"

    input:
    // single-end reads
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.filtered.fastq.gz")      , emit: reads
    tuple val(meta), path("*.vsearch_fastqfilter.log"), emit: log
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // the filters, e.g. --fastq_maxee 2 --fastq_maxns 0 (--fastq_qmax 93 for the PacBio HiFi base qualities)
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vsearch \\
        --fastq_filter ${reads} \\
        --fastqout ${prefix}.filtered.fastq \\
        --threads $task.cpus \\
        $args \\
        2> ${prefix}.vsearch_fastqfilter.log

    gzip -n ${prefix}.filtered.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | sed -n '1s/^vsearch v\\([0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.filtered.fastq.gz
    touch ${prefix}.vsearch_fastqfilter.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | sed -n '1s/^vsearch v\\([0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """
}
