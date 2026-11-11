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
    def inputs = meta.single_end || fastqs.size() == 1 ? "-s ${fastqs.join(' ')}" : "-f ${fastqs[0]} -r ${fastqs[1]}" + (fastqs.size() > 2 ? " -s ${fastqs[2]}" : "")
    """
    motus4_patch.py profile \\
        $args \\
        $inputs \\
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
