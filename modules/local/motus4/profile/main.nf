process MOTUS4_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:4.0.4--pyhdfd78af_0':
        'biocontainers/motus:4.0.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db, stageAs: 'db_mOTU' // motus4_patch.py expects the database at ./db_mOTU

    output:
    tuple val(meta), path("*.out"), emit: out
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.mgc"), optional: true, emit: mgc
    tuple val(meta), path("*.log")                , emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

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
    
    # no way to get the version number without triggering an exit code 2
    VERSION=\$(motus 2>&1 | grep "Version" | sed "s%^.*Version: %%" || true)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus4: \$VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.out
    touch ${prefix}.log

    # no way to get the version number without triggering an exit code 2
    VERSION=\$(motus 2>&1 | grep "Version" | sed "s%^.*Version: %%" || true)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus4: \$VERSION
    END_VERSIONS
    """
}
