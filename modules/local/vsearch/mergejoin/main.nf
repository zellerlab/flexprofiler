process VSEARCH_MERGEJOIN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de673c52a3310baca71b9cd912cd62fbfc24685899c4cafcb84f996303dc/data'
        : 'community.wave.seqera.io/library/vsearch:2.32.0--bee2e68f4ccce807'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.mergejoin.fastq.gz"), emit: reads
    tuple val(meta), path("*.vsearch.log")       , emit: log
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g;s/,.*//g;s/^v//;s/_.*//"'), emit: versions_vsearch, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.mergejoin.fastq.gz
    touch ${prefix}.vsearch.log
    """
}
