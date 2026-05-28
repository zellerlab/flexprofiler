process MOTUS_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input)
    path db

    output:
    tuple val(meta), path("*.txt") , optional: true, emit: txt
    tuple val(meta), path("*.biom"), optional: true, emit: biom
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmd_input = input.size() > 1 ? "-i ${input.join(',')}" : input.isDirectory() ? "-d ${input}" : "-i ${input}"
    def suffix = task.ext.args?.contains("-B") ? "biom" : "txt"
    """
    motus \\
        merge \\
        -db $db \\
        ${cmd_input} \\
        $args \\
        -o ${prefix}.${suffix}

    ## Take version from the mOTUs/profile module output, as cannot reconstruct
    ## version without having database staged in this directory.
    VERSION=\$(cat ${profile_version_yml} | grep '/*motus:.*' | sed 's/.*otus: //g')

    ## mOTUs version number is not available from command line.
    ## mOTUs save the version number in index database folder.
    ## mOTUs will check the database version is same version as exec version.
    if [ "$db" == "" ]; then
    else
        VERSION=\$(grep motus $db/db_mOTU_versions | sed 's/motus\\t//g')
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus3: \$VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.args?.contains("-B") ? "biom" : "txt"

    """
    touch ${prefix}.${suffix}

    if [ "$db" == "" ]; then
        VERSION=\$(echo \$(motus -h 2>&1) | sed 's/^.*Version: //; s/References.*\$//')
    else
        VERSION=\$(grep motus $db/db_mOTU_versions | sed 's/motus\\t//g')
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus3: \$VERSION
    END_VERSIONS
    """

}
