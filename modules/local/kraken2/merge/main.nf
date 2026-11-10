// Merge the Kraken2 reports of the samples of a database into tables of taxa by samples, counting reads (each mate
// counted), and the reads per million raw reads (bin/kraken2_merge.py)
process KRAKEN2_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'biocontainers/pandas:2.2.1'}"

    input:
    // meta.id: the database
    tuple val(meta), path(reports, stageAs: 'reports/*')
    // with a header, one row per report: file sample id run_accession mates
    path manifest
    // the read counts, for the reads per million, [] for none
    path read_counts

    output:
    tuple val(meta), path("*.reads_clade.tsv"), emit: reads_clade
    tuple val(meta), path("*.reads_taxon.tsv"), emit: reads_taxon
    tuple val(meta), path("*.rpm_clade.tsv")  , emit: rpm_clade, optional: true
    tuple val(meta), path("*.taxonomy.tsv")   , emit: taxonomy
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "kraken2_${meta.id}"
    """
    kraken2_merge.py \\
        ${manifest} \\
        reports \\
        ${read_counts ?: '-'} \\
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
