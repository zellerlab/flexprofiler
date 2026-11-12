// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FALCO {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/falco:1.2.5--h077b44d_0'
        : 'quay.io/biocontainers/falco:1.2.5--h077b44d_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('falco'), eval("falco --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_falco

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Each file is linked to a name made from the prefix (<prefix>, or <prefix>_1 and <prefix>_2 for paired-end reads,
    // as in the FastQC module) and run on its own, so that the outputs and the Filename line of the report, which
    // MultiQC uses as sample name, are named after the prefix. Falco detects the format from the extension.
    def files = reads instanceof Path ? [reads] : reads as List
    def names = files.size() == 1 ? [prefix] : (1..files.size()).collect { index -> "${prefix}_${index}" }
    def commands = [files, names].transpose().collect { read, name ->
        def link = "${name}.${read.name.endsWith('.gz') ? 'fastq.gz' : 'fastq'}"
        "[ -e ${link} ] || ln -s ${read} ${link}\n" +
            "falco ${args} --threads ${task.cpus} ${link} -D ${name}_fastqc_data.txt -S ${name}_summary.txt -R ${name}_report.html"
    }
    """
    ${commands.join('\n')}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def files = reads instanceof Path ? [reads] : reads as List
    def names = files.size() == 1 ? [prefix] : (1..files.size()).collect { index -> "${prefix}_${index}" }
    """
    ${names.collect { name -> "touch ${name}_fastqc_data.txt ${name}_summary.txt ${name}_report.html" }.join('\n')}
    """
}
