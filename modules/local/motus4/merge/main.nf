process MOTUS4_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:4.0.4--pyhdfd78af_0':
        'biocontainers/motus:4.0.4--pyhdfd78af_0' }"

    input:
    // the profiles of all the samples for one database, staged in a folder, so that the command lists them with a glob,
    // and not one by one
    tuple val(meta), path(input, stageAs: 'profiles/*')
    path db, stageAs: 'db_mOTU' // motus4_patch.py expects the database at ./db_mOTU

    output:
    tuple val(meta), path("*.txt") , optional: true, emit: txt
    tuple val("${task.process}"), val('motus'), eval('motus 2>&1 | sed -n "s/^.*Version: //p"'), emit: versions_motus, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // motus merge needs at least 2 profiles, and a single profile has already the format of a merged one
    def merge_cmd = input instanceof Path
        ? "cp profiles/* ${prefix}.txt"
        : "motus4_patch.py merge -i profile_list.txt $args -o ${prefix}.txt"
    """
    # motus merge treats a single -i argument as a file listing the profiles, so a list is passed
    printf '%s\\n' profiles/* > profile_list.txt

    $merge_cmd
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """

}
