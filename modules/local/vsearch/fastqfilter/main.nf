process VSEARCH_FASTQFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de673c52a3310baca71b9cd912cd62fbfc24685899c4cafcb84f996303dc/data'
        : 'community.wave.seqera.io/library/vsearch:2.32.0--bee2e68f4ccce807'}"

    input:
    // single-end reads
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.filtered.fastq.gz")      , emit: reads
    tuple val(meta), path("*.vsearch_fastqfilter.log"), emit: log
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g;s/,.*//g;s/^v//;s/_.*//"'), emit: versions_vsearch, topic: versions

    script:
    // the filters, e.g. --fastq_maxee 10 --fastq_maxns 0 (--fastq_qmax 93 for the PacBio HiFi base qualities)
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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.filtered.fastq.gz
    touch ${prefix}.vsearch_fastqfilter.log
    """
}
