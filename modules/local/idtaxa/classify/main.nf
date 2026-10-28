process IDTAXA_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-decipher:3.6.0--r45h01b2380_0'
        : 'biocontainers/bioconductor-decipher:3.6.0--r45h01b2380_0'}"

    input:
    // one sequence per read pair (merged or joined) or single-end read
    tuple val(meta), path(reads)
    path training_set

    output:
    tuple val(meta), path("*.idtaxa.classifications.tsv.gz"), emit: classifications
    tuple val(meta), path("*.idtaxa.profile.tsv")            , emit: profile
    tuple val(meta), path("*.idtaxa.log")                    , emit: log
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    idtaxa_classify.R \\
        $reads \\
        $training_set \\
        $prefix \\
        --processors $task.cpus \\
        $args \\
        > ${prefix}.idtaxa.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        decipher: \$(Rscript -e 'cat(as.character(packageVersion("DECIPHER")))')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.idtaxa.classifications.tsv.gz
    touch ${prefix}.idtaxa.profile.tsv ${prefix}.idtaxa.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        decipher: \$(Rscript -e 'cat(as.character(packageVersion("DECIPHER")))')
    END_VERSIONS
    """
}
