process VSEARCH_ORIENT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vsearch:2.32.0--hd2be7a0_0'
        : 'biocontainers/vsearch:2.32.0--hd2be7a0_0'}"

    input:
    // single-end reads
    tuple val(meta), path(reads)
    // the reference sequences, of the same gene as the reads, in its orientation (FASTA)
    path reference

    output:
    // the reads oriented as the reference, reverse-complemented if needed
    tuple val(meta), path("*.oriented.fastq.gz")   , emit: reads
    // the reads that could not be oriented
    tuple val(meta), path("*.notmatched.fastq.gz") , emit: notmatched
    tuple val(meta), path("*.vsearch_orient.log")  , emit: log
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g;s/,.*//g;s/^v//;s/_.*//"'), emit: versions_vsearch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vsearch \\
        --orient ${reads} \\
        --db ${reference} \\
        --fastqout ${prefix}.oriented.fastq \\
        --notmatched ${prefix}.notmatched.fastq \\
        --threads $task.cpus \\
        $args \\
        2> ${prefix}.vsearch_orient.log

    gzip -n ${prefix}.oriented.fastq ${prefix}.notmatched.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.oriented.fastq.gz
    echo "" | gzip > ${prefix}.notmatched.fastq.gz
    touch ${prefix}.vsearch_orient.log
    """
}
