process CAYMAN_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b3/b36790ad56052b28ef374a2d230862b0d7dd9d0d599b820d2e6bca34d0be57c4/data'
        : 'community.wave.seqera.io/library/gqlib:2.14.4--1a47e08f7db4071d'}"

    input:
    tuple val(meta), path(bam)
    path annot

    output:
    tuple val(meta), path("*.cazy.txt.gz"       ), emit: out
    tuple val(meta), path("*.aln_stats.txt.gz"  ), emit: aln_stats
    tuple val(meta), path("*.gene_counts.txt.gz"), emit: gene_counts
    tuple val(meta), path("*.cayman.log"        ), emit: log
    tuple val("${task.process}"), val('gqlib'), eval('python -c "import gqlib; print(gqlib.__version__)"'), emit: versions_gqlib, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # `cayman profile` runs bwa mem itself, then the gqlib domain profiler: the alignment is done
    # separately by BWA_ALIGN, and the gqlib domain profiler is run directly on it, as cayman 0.13.0
    # does. The cayman annotation csv is in the hmmer format.
    python -m gqlib \\
        --mode domain \\
        --db_format hmmer \\
        --bam $bam \\
        --db $annot \\
        --out_prefix $prefix \\
        $args \\
        &> ${prefix}.cayman.log

    # gqlib names the category of the annotation `feature`, cayman names it `cazy`, and cayman does
    # not report the `category` row (the counts of the reads assigned to any CAZy domain)
    gzip -cd ${prefix}.feature.txt.gz \\
        | awk -F '\\t' '\$1 != "category"' \\
        | gzip -n > ${prefix}.cazy.txt.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.aln_stats.txt.gz"
    touch "${prefix}.cazy.txt.gz"
    touch "${prefix}.gene_counts.txt.gz"
    touch "${prefix}.cayman.log"
    """

}
