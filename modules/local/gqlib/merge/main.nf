process GQLIB_MERGE {
    tag "${meta.id}|${meta.category}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gqlib:2.14.4--pyhdfd78af_0'
        : 'biocontainers/gqlib:2.14.4--pyhdfd78af_0'}"

    input:
    // the profiles of all the samples for one database and annotation (meta.id) and one category (meta.category)
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.txt.gz"), emit: txt
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def files = input instanceof List ? input : [input]
    // collate_counts takes the sample name from the file name, so the database and annotation names are removed from it
    def suffix = "_${meta.id}.${meta.category}.txt.gz"
    def links = files.collect { "ln -s ../${it} counts/${it.name - suffix}.${meta.category}.txt.gz" }.join('\n    ')
    """
    mkdir counts
    ${links}

    # one table per count column (uniq_raw, combined_raw, ...), as they are in the profiles
    columns=\$(python -c 'import gzip, sys; print(*gzip.open(sys.argv[1], "rt").readline().split()[1:])' ${files[0]})
    for column in \$columns; do
        python -m gqlib.bin.collate_counts counts \\
            --out_prefix ${prefix} \\
            --column \$column \\
            $args
    done

    VERSION=\$(python -c 'import gqlib; print(gqlib.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gqlib: \$VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.${meta.category}.uniq_raw.txt.gz"

    VERSION=\$(python -c 'import gqlib; print(gqlib.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gqlib: \$VERSION
    END_VERSIONS
    """

}
