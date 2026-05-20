process MOTUS4_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:4.0.4--pyhdfd78af_0':
        'biocontainers/motus:4.0.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input)
    path db // to stop docker saying it can't find it... would have to have the module in upstream steps anyway

    output:
    tuple val(meta), path("*.txt") , optional: true, emit: txt
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # mOTUs 4 does not support anymore the -db flag to pass a custom database
    # the database path is hardcoded in the module, and this snippet overrides it to
    # point to the local work directory
    rm -f motus_patch.py
    cat <<-'PATCH_EOF' > motus_patch.py
    import sys
    import pathlib
    from motus.motus import main
    import motus.mutils as m

    m.DEFAULT_MOTUS_MGDB_PARENT_LOCATION = pathlib.Path(__file__).resolve().parent
    m.DEFAULT_MOTUS_MGDB_LOCATION = m.DEFAULT_MOTUS_MGDB_PARENT_LOCATION.joinpath('db_mOTU')
    m.DEFAULT_MOTUS_MGDB_LOCATION_MARKER = m.DEFAULT_MOTUS_MGDB_LOCATION.joinpath('db_mOTU.downloaded')
    m.DEFAULT_MOTUS_ANNODB_LOCATION = m.DEFAULT_MOTUS_MGDB_LOCATION.joinpath('mOTUsv4.0.annotation.db')
    m.DEFAULT_MOTUS_ANNODB_LOCATION_MARKER = m.DEFAULT_MOTUS_MGDB_LOCATION.joinpath('mOTUsv4.0.annotation.db.downloaded')

    if __name__ == '__main__':
        sys.argv[0] = sys.argv[0].removesuffix('.exe')
        sys.exit(main())
    PATCH_EOF

    python motus_patch.py \\
        merge \\
        -i ${input} \\
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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    # no way to get the version number without triggering an exit code 2
    VERSION=\$(motus 2>&1 | grep "Version" | sed "s%^.*Version: %%" || true)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus4: \$VERSION
    END_VERSIONS
    """

}
