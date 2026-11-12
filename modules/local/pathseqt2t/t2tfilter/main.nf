// Third step of the PathSeq-T2T host filtering (Dohlman et al., Cell 2026, https://doi.org/10.1016/j.cell.2026.04.015):
// the reads left by PathSeqFilterSpark are aligned to the complete human genome T2T-CHM13 with bwa mem -T 0, the
// paired and the unpaired reads separately, and those aligned to it are removed.
//
// Reimplements src/commands/t2tfilter.sh of https://github.com/abdohlman/pathseq-t2t (v0.3.0, commit 51d74430), with
// the same read selection. Differences (see docs/pathseq_t2t.md):
// - upstream writes its BAMs through background process substitutions (t2tfilter.sh:251-253, 280-281): here each BAM
//   is written to a file by a command that has finished before the next one starts
// - the BAMs are converted to FASTQ with samtools fastq instead of Picard SamToFastq (t2tfilter.sh:231, 236), which is
//   not in the container, leaving out the same reads as Picard by default (-F 0xB00: secondary, supplementary and
//   failing the vendor quality checks)
// - no reads of a kind (paired or unpaired) give an empty FASTQ instead of an error
// - the reads are written as FASTQ, and a read of a proper pair kept without its mate (the mate being aligned to
//   T2T-CHM13) is dropped here: upstream keeps it in the BAM of the paired reads, and drops it when converting that
//   BAM to FASTQ for the classification (classify.sh, Picard SamToFastq with a lenient validation)
// - the reads aligned to T2T-CHM13 that overlap the viral decoy contigs of GRCh38 are not rescued (t2tfilter.sh:255-266,
//   284-295): T2T-CHM13 has none of these contigs, so upstream rescues no read
process PATHSEQT2T_T2TFILTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    // the paired and the unpaired reads of PATHSEQT2T_QCFILTER, any number of BAMs of each, [] if none
    tuple val(meta), path(paired, stageAs: 'paired/*'), path(unpaired, stageAs: 'unpaired/*')
    // the bwa index of T2T-CHM13
    tuple val(meta2), path(index)

    output:
    // the reads not aligned to T2T-CHM13: the pairs (paired-end reads only), and the reads without their mate (for
    // single-end reads, all the reads); empty if no read is left
    tuple val(meta), path("*.t2tfilt_{1,2}.fastq.gz")    , emit: paired  , optional: true
    tuple val(meta), path("*.t2tfilt_unpaired.fastq.gz") , emit: unpaired
    // the alignments to T2T-CHM13, of the pairs (paired-end reads only) and of the reads without their mate
    tuple val(meta), path("*.t2t.{paired,unpaired}.bam") , emit: bam
    // the reads before and after the filter, for the read counts
    tuple val(meta), path("*.t2tfilter.log")             , emit: log
    tuple val("${task.process}"), val('bwa'), eval('bwa 2>&1 | sed -n "s/^Version: //p"'), emit: versions_bwa, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    script:
    // bwa mem, e.g. -T 0
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : 'true'
    // no supplementary alignments, and none of the tags that depend on the alignment
    def clean = "-F 2048 -x SA -x OQ -x MD"
    """
    set -euo pipefail

    INDEX=\$(find -L "${index}" -name "*.amb" | sed 's/\\.amb\$//')

    # the number of reads of a FASTQ
    count_fastq() {
        local fastq=\$1
        bgzip -dc "\$fastq" | awk 'END { print NR / 4 }'
    }

    # the reads of input_bam that are not aligned to T2T-CHM13 (flags: the flags of those aligned, unquoted so that
    # none can be given) are written to output_bam
    not_aligned() {
        local input_bam=\$1 flags=\$2 output_bam=\$3
        samtools view -@ ${task.cpus} -b -h \$flags -e '[AS]>35' -U not_aligned.tmp.bam -o /dev/null "\$input_bam"
        samtools view -@ ${task.cpus} -b -h ${clean} -o "\$output_bam" not_aligned.tmp.bam
        rm -f not_aligned.tmp.bam
    }

    input_pairs=0
    retained_pairs=0
    dropped_mates=0
    if [ -n "${paired_end}" ]; then
        : > in_1.tmp.fastq
        : > in_2.tmp.fastq
        for bam in paired/*.bam; do
            [ -e "\$bam" ] || continue
            samtools collate -@ ${task.cpus} -u -O "\$bam" collate.tmp \\
                | samtools fastq -n -F 0xB00 -1 one_1.tmp.fastq -2 one_2.tmp.fastq -0 /dev/null -s /dev/null -
            cat one_1.tmp.fastq >> in_1.tmp.fastq
            cat one_2.tmp.fastq >> in_2.tmp.fastq
        done
        input_pairs=\$(awk 'END { print NR / 4 }' in_1.tmp.fastq)
        # without reads, bwa mem writes only the header, and all the outputs are empty
        bwa mem ${args} -t ${task.cpus} \$INDEX in_1.tmp.fastq in_2.tmp.fastq \\
            | samtools view -@ ${task.cpus} -b -o ${prefix}.t2t.paired.bam -
        # properly paired, with an alignment score above 35
        not_aligned ${prefix}.t2t.paired.bam "-f 3" kept.tmp.bam
        # a read kept without its mate (-s) is dropped
        samtools collate -@ ${task.cpus} -u -O kept.tmp.bam collate.tmp \\
            | samtools fastq -n -F 0xB00 -1 ${prefix}.t2tfilt_1.fastq.gz -2 ${prefix}.t2tfilt_2.fastq.gz -0 /dev/null -s dropped.tmp.fastq -
        dropped_mates=\$(awk 'END { print NR / 4 }' dropped.tmp.fastq)
        retained_pairs=\$(count_fastq ${prefix}.t2tfilt_1.fastq.gz)
        if [ "\$retained_pairs" != "\$(count_fastq ${prefix}.t2tfilt_2.fastq.gz)" ]; then
            echo "The mates of the pairs not aligned to T2T-CHM13 are not the same number" >&2
            exit 1
        fi
    fi

    : > in_unpaired.tmp.fastq
    for bam in unpaired/*.bam; do
        [ -e "\$bam" ] || continue
        samtools fastq -n -F 0xB00 "\$bam" >> in_unpaired.tmp.fastq
    done
    input_unpaired=\$(awk 'END { print NR / 4 }' in_unpaired.tmp.fastq)
    bwa mem ${args} -t ${task.cpus} \$INDEX in_unpaired.tmp.fastq \\
        | samtools view -@ ${task.cpus} -b -o ${prefix}.t2t.unpaired.bam -
    # with an alignment score above 35 (the reads not aligned have none)
    not_aligned ${prefix}.t2t.unpaired.bam "" kept.tmp.bam
    # aligned as single reads, so none is flagged as a mate (READ1 or READ2): all are written to -0
    samtools fastq -@ ${task.cpus} -n -F 0xB00 -0 ${prefix}.t2tfilt_unpaired.fastq.gz kept.tmp.bam
    retained_unpaired=\$(count_fastq ${prefix}.t2tfilt_unpaired.fastq.gz)

    rm -f *.tmp.fastq *.tmp.bam

    cat <<-END_LOG > ${prefix}.t2tfilter.log
    Input pairs: \$input_pairs
    Input unpaired reads: \$input_unpaired
    Retained pairs: \$retained_pairs
    Dropped reads kept without their mate: \$dropped_mates
    Retained unpaired reads: \$retained_unpaired
    END_LOG
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : 'true'
    """
    if [ -n "${paired_end}" ]; then
        echo "" | gzip > ${prefix}.t2tfilt_1.fastq.gz
        echo "" | gzip > ${prefix}.t2tfilt_2.fastq.gz
        touch ${prefix}.t2t.paired.bam
    fi
    echo "" | gzip > ${prefix}.t2tfilt_unpaired.fastq.gz
    touch ${prefix}.t2t.unpaired.bam
    touch ${prefix}.t2tfilter.log
    """
}
