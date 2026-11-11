// Second step of the PathSeq-T2T host filtering (Dohlman et al., Cell 2026, https://doi.org/10.1016/j.cell.2026.04.015):
// GATK PathSeqFilterSpark on the reads not aligned to GRCh38, and on those aligned to its viral decoy contigs. It
// removes the reads aligned to the host, of low quality or complexity, with host k-mers, and the duplicates.
//
// Reimplements src/commands/qcfilter.sh of https://github.com/abdohlman/pathseq-t2t (v0.3.0, commit 51d74430), with the
// same commands. Differences (see docs/pathseq_t2t.md):
// - upstream stops if PathSeqFilterSpark writes no paired or no unpaired reads (qcfilter.sh:304-314), which it does
//   when none is left: here a missing BAM is an error only if the filter metrics count reads for it
// - the BAMs of the two inputs are not merged here, but when they are converted to FASTQ (PATHSEQT2T_T2TFILTER)
process PATHSEQT2T_QCFILTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // with Java 17 (environment.yml): the Spark tools of GATK 4.6 do not run on Java 24
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33fdebc95eb14966914fe7294717b9e852dcd94beaed459d832f3693875a8c17/data'
        : 'community.wave.seqera.io/library/gatk4_openjdk:cb0a04a09c5bc584'}"
    // Spark (Hadoop) needs the user to be known to the system, and the user set by docker.runOptions (-u) is not known to
    // the container: with Docker the task runs as root, the last -u being the one used, and gives its files back to the
    // owner of the task directory when it ends
    containerOptions "${workflow.containerEngine == 'docker' ? '-u 0:0' : ''}"

    input:
    // the reads not aligned to the host, and those aligned to its decoy contigs ([] if none), of PATHSEQT2T_PREFILTER
    tuple val(meta), path(unaligned), path(decoys)
    // the host files of PathSeq: pathseq_host.bfi (k-mers) and pathseq_host.fa.img (BWA index image)
    path host_dir

    output:
    // the reads passing the filters with their mate, and without it; of each input, none if no read is left
    tuple val(meta), path("*.qcfilt_paired.bam")  , emit: paired  , optional: true
    tuple val(meta), path("*.qcfilt_unpaired.bam"), emit: unpaired, optional: true
    // the reads after each filter, of each input
    tuple val(meta), path("*.filter_metrics.txt") , emit: metrics
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), emit: versions_gatk4, topic: versions

    script:
    // e.g. --min-clipped-read-length
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the BWA index image is loaded outside of the Java heap, as in upstream (48 of 64 GB)
    def heap_mb = task.memory ? (task.memory.toMega() * 0.75).intValue() : 12288
    def decoys_bam = decoys ? "${decoys}" : ''
    """
    set -euo pipefail

    # with Docker the task runs as root (see containerOptions)
    give_back() {
        if [ "\$(id -u)" = 0 ]; then
            chown -R "\$(stat -c '%u:%g' .)" . || true
        fi
    }
    trap give_back EXIT

    mkdir -p tmp

    # PathSeqFilterSpark on input_bam, writing the outputs named by name, with the BAM partition size partition_size
    pathseq_filter() {
        local input_bam=\$1 name=\$2 partition_size=\$3
        local kind count bam
        gatk --java-options "-Xmx${heap_mb}M" PathSeqFilterSpark \\
            --input "\$input_bam" \\
            --tmp-dir tmp \\
            --spark-master "local[${task.cpus}]" \\
            --bam-partition-size "\$partition_size" \\
            --is-host-aligned true \\
            --kmer-file ${host_dir}/pathseq_host.bfi \\
            --filter-bwa-image ${host_dir}/pathseq_host.fa.img \\
            --filter-metrics ${prefix}.prefilter.\$name.filter_metrics.txt \\
            --paired-output ${prefix}.\$name.qcfilt_paired.bam \\
            --unpaired-output ${prefix}.\$name.qcfilt_unpaired.bam \\
            ${args}

        # PathSeqFilterSpark writes no BAM if it has no reads for it: a BAM is expected only if the metrics count reads
        for kind in PAIRED UNPAIRED; do
            count=\$(awk -F '\\t' -v column=FINAL_\${kind}_READS '
                header { print \$index_of[column]; found = 1; exit }
                \$1 == "PRIMARY_READS" { for (i = 1; i <= NF; i++) index_of[\$i] = i; header = 1 }
                END { if (!found) exit 1 }
            ' ${prefix}.prefilter.\$name.filter_metrics.txt)
            bam=${prefix}.\$name.qcfilt_\$(echo \$kind | tr '[:upper:]' '[:lower:]').bam
            if [ "\$count" -gt 0 ] && [ ! -s "\$bam" ]; then
                echo "PathSeqFilterSpark counts \$count reads (FINAL_\${kind}_READS) but did not write \$bam" >&2
                exit 1
            fi
        done
    }

    pathseq_filter ${unaligned} unaligned 4000000

    if [ -n "${decoys_bam}" ]; then
        # the alignments to the decoy contigs are removed, or the reads would be filtered as aligned to the host
        gatk --java-options "-Xmx${heap_mb}M" RevertSam \\
            --INPUT "${decoys_bam}" \\
            --OUTPUT decoys.reverted.tmp.bam \\
            --REMOVE_DUPLICATE_INFORMATION false \\
            --RESTORE_ORIGINAL_QUALITIES false \\
            --TMP_DIR tmp \\
            --VERBOSITY ERROR
        pathseq_filter decoys.reverted.tmp.bam decoys 0
    fi

    rm -rf tmp decoys.reverted.tmp.bam *.sbi *.bai
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unaligned.qcfilt_paired.bam ${prefix}.unaligned.qcfilt_unpaired.bam
    touch ${prefix}.prefilter.unaligned.filter_metrics.txt
    """
}
