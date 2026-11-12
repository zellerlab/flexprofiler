process IDTAXA_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/58/58630872068a0e2afac569fe03330f6f1b83209da57c5b9e01d38f2d8dbd51c7/data'
        : 'community.wave.seqera.io/library/bioconductor-decipher:3.6.0--297c39bd2d9527b3'}"

    input:
    // one sequence per read pair (merged or joined) or single-end read
    tuple val(meta), path(reads)
    path training_set

    output:
    // the number of sequences of each lineage (the taxa assigned above the threshold, ";"-separated), with its ranks, and
    // the rows unassigned and total_reads
    tuple val(meta), path("*.idtaxa.profile.tsv"), emit: profile
    // the number of sequences and of assigned sequences, for the read counts
    tuple val(meta), path("*.idtaxa.log")        , emit: log
    tuple val("${task.process}"), val('decipher'), eval('Rscript -e "cat(as.character(packageVersion(\'DECIPHER\')))"'), emit: versions_decipher, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the IdTaxa options, with its defaults
    def threshold = task.ext.threshold ?: 60
    def strand = task.ext.strand ?: 'both'
    // the R script is in bin (it is not the script of the process, as Nextflow allows eval outputs, the version, only in
    // Bash scripts)
    """
    idtaxa_classify.R \\
        ${training_set} \\
        ${reads} \\
        ${prefix} \\
        ${strand} \\
        ${threshold} \\
        ${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idtaxa.profile.tsv ${prefix}.idtaxa.log
    """
}
