process PIMENTO_STD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mi-pimento:1.3.3--pyhdfd78af_0'
        : 'biocontainers/mi-pimento:1.3.3--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    // one per mate (single-end: one), named <prefix>_<mate>_std_primers.fasta, empty if no primer was found
    tuple val(meta), path("*_std_primers.fasta"), emit: primers
    tuple val(meta), path("*_std_primer_out.txt"), emit: out
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Paired-end mates are searched separately, as unmerged reads: all the primers of the library, forward and
    // reverse, are searched at the start of the reads, so each mate gets the primer it starts with. Single-end reads
    // are searched in --merged mode, for the forward primers at the start and the reverse ones at the end.
    // PIMENTO writes all_standard_primer_proportions.txt in the working directory, so each search runs in its own.
    def searches = meta.single_end
        ? [[1, reads instanceof List ? reads[0] : reads, '--merged']]
        : [[1, reads[0], ''], [2, reads[1], '']]
    def commands = searches.collect { mate, fastq, mode ->
        """
    mkdir mate${mate}
    (cd mate${mate} && pimento std -i ../${fastq} -o ${prefix}_${mate} ${mode} $args > ../${prefix}_${mate}.pimento.log 2>&1)
    mv mate${mate}/${prefix}_${mate}_std_primers.fasta mate${mate}/${prefix}_${mate}_std_primer_out.txt .
    mv mate${mate}/all_standard_primer_proportions.txt ${prefix}_${mate}_all_standard_primer_proportions.txt"""
    }.join('\n')
    """
    ${commands}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pimento: \$(pimento --version 2>&1 | sed -n 's/.*version //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_1_std_primers.fasta ${prefix}_1_std_primer_out.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pimento: \$(pimento --version 2>&1 | sed -n 's/.*version //p')
    END_VERSIONS
    """
}
