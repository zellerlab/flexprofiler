process MOTUS4_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:4.0.4--pyhdfd78af_0':
        'biocontainers/motus:4.0.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input)
    path db, stageAs: 'db_mOTU' // motus4_patch.py expects the database at ./db_mOTU

    output:
    tuple val(meta), path("*.txt") , optional: true, emit: txt
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # motus merge treats a single -i argument as a file listing the profiles,
    # so always pass a list to also support a single input profile
    printf '%s\\n' ${input} > profile_list.txt

    motus4_patch.py \\
        merge \\
        -i profile_list.txt \\
        $args \\
        -o ${prefix}.txt

    # no way to get the version number without triggering an exit code 2
    VERSION=\$(motus 2>&1 | grep "Version" | sed "s%^.*Version: %%" || true)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus: \$VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    # no way to get the version number without triggering an exit code 2
    version=\$(motus 2>&1 | grep "version" | sed "s%^.*version: %%" || true)

    cat <<-end_versions > versions.yml
    "${task.process}":
        motus4: \$version
    end_versions
    """

}
