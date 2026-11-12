// Merge the Kraken2 reports of the samples of a database into tables of taxa by samples, counting reads (each mate
// counted) (bin/kraken2_merge.py)
process KRAKEN2_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84e799e45917de0293324e738f0c0b7d1fd0f85ac0536d39d1161d7ef372fbaa/data'
        : 'community.wave.seqera.io/library/pandas:2.2.1--2558de69969e58b9'}"

    input:
    // meta.id: the database; the reports of the read pairs, and those of the single reads (single-end reads, and the
    // reads left without their mate by --pathseq_t2t)
    tuple val(meta), path(pairs, stageAs: 'pairs/*'), path(single_reads, stageAs: 'single_reads/*')

    output:
    tuple val(meta), path("*.reads_clade.tsv"), emit: reads_clade
    tuple val(meta), path("*.reads_taxon.tsv"), emit: reads_taxon
    tuple val(meta), path("*.taxonomy.tsv")   , emit: taxonomy
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "kraken2_${meta.id}"
    """
    touch ${prefix}.reads_clade.tsv ${prefix}.reads_taxon.tsv ${prefix}.taxonomy.tsv
    """
}
