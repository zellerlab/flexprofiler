// First step of the PathSeq-T2T host filtering (Dohlman et al., Cell 2026, https://doi.org/10.1016/j.cell.2026.04.015):
// from the reads aligned to GRCh38, select those not aligned to it, and those aligned to its viral decoy contigs.
//
// Reimplements src/commands/prefilter.sh of https://github.com/abdohlman/pathseq-t2t (v0.3.0, commit 51d74430), with the
// same read selection. Differences (see docs/pathseq_t2t.md):
// - upstream writes its BAMs through background process substitutions (prefilter.sh:142-155) and validates them after
//   a sleep, so a BAM can be read while still being written: here each BAM is written to a file by a command that
//   has finished before the next one starts
// - single-end reads are selected by their alignment score alone, as upstream does for the unpaired reads in
//   t2tfilter.sh: upstream requires a proper pair here, which no single-end read is
// - the decoy BAM is written only if it has reads
process PATHSEQT2T_PREFILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    // the reads aligned to GRCh38 with bwa mem -T 0, with a read group (PathSeqFilterSpark drops the reads without one)
    tuple val(meta), path(bam)
    // the viral decoy contigs of the reference, [] to not rescue the reads aligned to them
    path decoys_bed

    output:
    // the reads not aligned to the host: all but the proper pairs (single-end: reads) with an alignment score above 35
    tuple val(meta), path("*.prefilter.unaligned.bam"), emit: unaligned
    // the reads aligned to the host that overlap the decoy contigs, only if there are any
    tuple val(meta), path("*.prefilter.decoys.bam")   , emit: decoys, optional: true
    // samtools flagstat of the input: its primary reads are the raw reads of the read counts
    tuple val(meta), path("*.flagstat.tsv")           , emit: flagstat
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // paired-end: properly paired (prefilter.sh, aligner bwa); single-end: aligned (the reads not aligned have no AS)
    def host_aligned = meta.single_end ? "-e '[AS]>35'" : "-f 3 -e '[AS]>35'"
    // no supplementary alignments, and none of the tags that depend on the alignment
    def clean = "-F 2048 -x SA -x OQ -x MD"
    def bed = decoys_bed ? "${decoys_bed}" : ''
    """
    set -euo pipefail

    samtools flagstat --output-fmt tsv -@ ${task.cpus} ${bam} > ${prefix}.flagstat.tsv

    samtools view -@ ${task.cpus} -b -h ${host_aligned} -U not_host_aligned.tmp.bam -o host_aligned.tmp.bam ${bam}
    samtools view -@ ${task.cpus} -b -h ${clean} -o ${prefix}.prefilter.unaligned.bam not_host_aligned.tmp.bam
    if [ -n "${bed}" ]; then
        samtools view -@ ${task.cpus} -b -h ${clean} -L "${bed}" -o decoys.tmp.bam host_aligned.tmp.bam
        if [ "\$(samtools view -c decoys.tmp.bam)" -gt 0 ]; then
            mv decoys.tmp.bam ${prefix}.prefilter.decoys.bam
        fi
    fi
    rm -f not_host_aligned.tmp.bam host_aligned.tmp.bam decoys.tmp.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.prefilter.unaligned.bam ${prefix}.flagstat.tsv
    """
}
