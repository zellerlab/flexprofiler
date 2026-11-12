// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MOTUS_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.out"), emit: out
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.mgc"), emit: mgc, optional: true
    tuple val(meta), path("*.log"), emit: log
    // WARN: Version information not provided by tool on CLI.  Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('motus'), val("3.1.0"), topic: versions, emit: versions_motus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the paired-end reads can be followed by the reads left without their mate, given as unpaired reads, or be only
    // those if no pair is left (a single file is staged as a path, not a list)
    def fastqs = reads instanceof List ? reads : [reads]
    def inputs = reads[0].getExtension() == 'bam' ? "-i ${reads}" :
                reads[0].getExtension() == 'mgc' ? "-m $reads" :
                    meta.single_end || fastqs.size() == 1 ? "-s $reads" :
                    "-f ${fastqs[0]} -r ${fastqs[1]}" + (fastqs.size() > 2 ? " -s ${fastqs[2]}" : "")
    def refdb = db ? "-db ${db}" : ""
    """
    motus profile \\
        ${args} \\
        ${inputs} \\
        ${refdb} \\
        -t ${task.cpus} \\
        -n ${prefix} \\
        -o ${prefix}.out \\
        2>| >(tee ${prefix}.log >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.out
    touch ${prefix}.log
    """
}
