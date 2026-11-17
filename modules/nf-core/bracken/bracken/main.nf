// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BRACKEN_BRACKEN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f30aa99d8d4f6ff1104f56dbacac95c1dc0905578fb250c80f145b6e80703bd1/data':
        'community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01' }"

    input:
    tuple val(meta), path(kraken_report)
    path database

    output:
    tuple val(meta), path(bracken_report)        , emit: reports
    tuple val(meta), path(bracken_kraken_style_report), emit: txt
    tuple val("${task.process}"), val('bracken'), eval('bracken -v | cut -f2 -d"v"'), topic: versions, emit: versions_bracken

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    bracken_report = "${prefix}.tsv"
    bracken_kraken_style_report = "${prefix}.kraken2.report_bracken.txt"
    // A report with no reads at the rank of the estimation (an empty report, none classified, or all the taxa of the
    // rank below the threshold, which Bracken fails on with "no reads found") gives a table with only the header and an
    // empty report, so that the sample is kept. The errors of the bracken script, which exits with status 0 (e.g. no
    // k-mer distribution for the read length), fail the task.
    """
    if [ ! -s '${kraken_report}' ]; then
        echo "Empty Kraken2 report: no taxa"
    else
        status=0
        bracken \\
            ${args} \\
            -d '${database}' \\
            -i '${kraken_report}' \\
            -o '${bracken_report}' \\
            -w '${bracken_kraken_style_report}' \\
            > bracken.log 2>&1 \\
            || status=\$?
        cat bracken.log
        if grep -q ' ERROR:' bracken.log; then
            exit 1
        fi
        if [ \$status -ne 0 ]; then
            grep -q '^Error: no reads found' bracken.log || exit \$status
        fi
    fi
    if [ ! -s '${bracken_report}' ]; then
        printf 'name\\ttaxonomy_id\\ttaxonomy_lvl\\tkraken_assigned_reads\\tadded_reads\\tnew_est_reads\\tfraction_total_reads\\n' > '${bracken_report}'
        : > '${bracken_kraken_style_report}'
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    bracken_report = "${prefix}.tsv"
    bracken_kraken_style_report = "${prefix}.kraken2.report_bracken.txt"
    """
    touch ${prefix}.tsv
    touch ${bracken_kraken_style_report}
    """
}
