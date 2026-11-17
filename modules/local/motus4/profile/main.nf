process MOTUS4_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/3930809d11a1016e5b92efcdb77a68678b63ec2a8ee0a1e0cf5656cc7efdce1a/data'
        : 'community.wave.seqera.io/library/motus:4.0.4--79a63faa164322d2'}"

    input:
    tuple val(meta), path(reads)
    path db, stageAs: 'db_mOTU' // motus4_patch.py expects the database at ./db_mOTU

    output:
    tuple val(meta), path("*.out"), emit: out
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.mgc"), optional: true, emit: mgc
    tuple val(meta), path("*.log")                , emit: log
    tuple val("${task.process}"), val('motus4'), eval('motus 2>&1 | sed -n "s/^.*Version: //p"'), emit: versions_motus4, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // mOTUs 4 does not support anymore additional inputs formats except for FASTA,FASTQ. The paired-end reads can be
    // followed by the reads left without their mate by --pathseq_t2t, given as unpaired reads, or be only those if no
    // pair is left (a single file is staged as a path, not a list)
    def fastqs = reads instanceof List ? reads : [reads]
    def inputs = meta.single_end || fastqs.size() == 1 ? "-s ${fastqs.join(' ')}" : "-f ${fastqs[0]} -r ${fastqs[1]}"
    // the unpaired reads after the pairs (the singletons of --pathseq_t2t), given with -s only if the file has reads:
    // it is empty when no read was left without its mate, and mOTUs fails on an empty file (ZeroDivisionError)
    def unpaired = !meta.single_end && fastqs.size() > 2 ? fastqs[2] : ''
    """
    # whether a FASTQ (gzipped or not) has any read: only its first byte is read, not the whole file. With python3, as
    # the container has no gzip
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

    motus4_patch.py profile \\
        $args \\
        $inputs \\
        \$unpaired \\
        -t $task.cpus \\
        -n $prefix \\
        -o ${prefix}.out \\
        2> ${prefix}.log || {
        exit_code=\$?
        # show the log in .command.err, which is included in the error report
        cat ${prefix}.log >&2
        exit \$exit_code
    }
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.out
    touch ${prefix}.log
    """
}
