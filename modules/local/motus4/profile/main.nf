process MOTUS4_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:4.0.4--pyhdfd78af_0':
        'biocontainers/motus:4.0.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

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
    // mOTUs 4 does not support anymore additional inputs formats except for FASTA,FASTQ
    def inputs = meta.single_end ? "-s $reads" : "-f ${reads[0]} -r ${reads[1]}"
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

    python motus_patch.py profile \\
        $args \\
        $inputs \\
        -t $task.cpus \\
        -n $prefix \\
        -o ${prefix}.out \\
        2>| >(tee ${prefix}.log >&2)
    
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
    def inputs = meta.single_end ? "-s $reads" : "-f ${reads[0]} -r ${reads[1]}"
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
