process GQLIB_PROFILE_GENES {
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
    tuple val(meta), path("*.gene_counts.txt.gz"   ), emit: gene_counts
    tuple val(meta), path("aln_stats/*.aln_stats.txt.gz"), emit: aln_stats
    // all the profiles: one per functional category of the annotation database, plus the gene counts
    tuple val(meta), path("*.txt.gz"               ), emit: profiles
    tuple val(meta), path("*.gqlib.log"            ), emit: log
    tuple val("${task.process}"), val('gqlib'), eval('python -c "import gqlib; print(gqlib.__version__)"'), emit: versions_gqlib, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # alignment is done separately by BWA_MEM, gqlib takes the alignment file as input
    python -m gqlib \\
        --mode genes \\
        --bam $bam \\
        --db $annot \\
        --out_prefix $prefix \\
        $args \\
        &> ${prefix}.gqlib.log

    # moved to a subdirectory so that it is not captured with the profiles
    mkdir aln_stats
    mv ${prefix}.aln_stats.txt.gz aln_stats/
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.gene_counts.txt.gz"
    mkdir aln_stats
    touch "aln_stats/${prefix}.aln_stats.txt.gz"
    touch "${prefix}.gqlib.log"
    """

}
