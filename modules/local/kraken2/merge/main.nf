// Merge the Kraken2 reports of the samples of a database into tables of taxa by samples, counting reads (each mate
// counted) (bin/kraken2_merge.py)
process KRAKEN2_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'biocontainers/pandas:2.2.1'}"

    input:
    // meta.id: the database; the reports of the read pairs, and those of the single reads (single-end reads, and the
    // reads left without their mate by --pathseq_t2t)
    tuple val(meta), path(pairs, stageAs: 'pairs/*'), path(single_reads, stageAs: 'single_reads/*')

    output:
    tuple val(meta), path("*.reads_clade.tsv"), emit: reads_clade
    tuple val(meta), path("*.reads_taxon.tsv"), emit: reads_taxon
    tuple val(meta), path("*.taxonomy.tsv")   , emit: taxonomy
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "kraken2_${meta.id}"
    """
    # a directory is not staged if it has no reports
    mkdir -p pairs single_reads
    kraken2_merge.py \\
        ${meta.id} \\
        pairs \\
        single_reads \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "kraken2_${meta.id}"
    """
    touch ${prefix}.reads_clade.tsv ${prefix}.reads_taxon.tsv ${prefix}.taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
