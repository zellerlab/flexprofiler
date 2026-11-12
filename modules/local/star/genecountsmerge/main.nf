// Merge the reads per gene of STAR (--quantMode GeneCounts) of the samples into tables of genes by samples, one per
// strandedness of the library, and a summary table (bin/star_genecounts_merge.py)
process STAR_GENECOUNTSMERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84e799e45917de0293324e738f0c0b7d1fd0f85ac0536d39d1161d7ef372fbaa/data'
        : 'community.wave.seqera.io/library/pandas:2.2.1--2558de69969e58b9'}"

    input:
    // a sheet of the sample names and the names of their ReadsPerGene.out.tab, one per line, separated by a tab
    // (bin/star_genecounts_merge.py), and the ReadsPerGene.out.tab
    tuple val(meta), path(sheet), path(counts)

    output:
    tuple val(meta), path("*.{unstranded,forward,reverse}.tsv"), emit: counts
    tuple val(meta), path("*.summary.tsv")                      , emit: summary
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    star_genecounts_merge.py \\
        ${prefix} \\
        ${sheet}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unstranded.tsv ${prefix}.forward.tsv ${prefix}.reverse.tsv ${prefix}.summary.tsv
    """
}
