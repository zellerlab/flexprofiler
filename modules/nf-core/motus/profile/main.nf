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
                    "-f ${fastqs[0]} -r ${fastqs[1]}"
    // the unpaired reads after the pairs (the singletons of --pathseq_t2t), given with -s only if the file has reads:
    // it is empty when no read was left without its mate, and mOTUs rejects an empty file ("is not a fastq file")
    def unpaired = !meta.single_end && fastqs.size() > 2 ? fastqs[2] : ''
    def refdb = db ? "-db ${db}" : ""
    """
    # whether a FASTQ (gzipped or not) has any read: only its first byte is read, not the whole file
    has_reads() {
        python3 - "\$1" <<'END_PYTHON'
    import gzip, sys
    with open(sys.argv[1], 'rb') as f:
        gzipped = f.read(2) == b'\\x1f\\x8b'
        f.seek(0)
        content = gzip.open(f) if gzipped else f
        sys.exit(0 if content.read(1) else 1)
    END_PYTHON
    }

    unpaired=""
    if [ -n "${unpaired}" ] && has_reads "${unpaired}"; then
        unpaired="-s ${unpaired}"
    fi

    motus profile \\
        ${args} \\
        ${inputs} \\
        \$unpaired \\
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
