//
// Remove host reads via alignment and export off-target reads
//

include { BOWTIE2_BUILD                            } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                            } from '../../modules/nf-core/bowtie2/align/main'
include { BBMAP_REPAIR as BBMAP_REPAIR_HOSTREMOVAL } from '../../modules/local/bbmap/repair/main'
include { SAMTOOLS_INDEX                           } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                           } from '../../modules/nf-core/samtools/stats/main'
include { publishFiles                             } from './utils_zellerlab_flexprofiler_pipeline'

workflow SHORTREAD_HOSTREMOVAL {
    take:
    reads     // [ [ meta ], [ reads ] ]
    reference // /path/to/fasta
    index     // /path/to/index

    main:
    ch_multiqc_files = channel.empty()

    ch_published = channel.empty()

    if (!params.shortread_hostremoval_index) {
        ch_bowtie2_index = BOWTIE2_BUILD([[], reference]).index
        if (params.save_hostremoval_index) {
            ch_published = ch_published.mix(publishFiles(BOWTIE2_BUILD.out.index, 'bowtie2/build', 'bowtie2'))
        }
    }
    else {
        ch_bowtie2_index = index.first()
    }

    // The mates of the paired-end reads are aligned separately, as single-end reads, and a pair is kept only if
    // neither mate aligns to the host, following https://doi.org/10.1016/j.crmeth.2025.101218. Aligning the pairs as
    // pairs would keep the pairs with one mate aligned to the host, as bowtie2 --un-conc writes all the pairs not
    // aligned concordantly. meta.mate (1 or 2) tells the mates apart, and the single-end reads have none.
    ch_reads_for_bowtie2 = reads.flatMap { meta, fastqs ->
        meta.single_end
            ? [[meta, fastqs]]
            : [[meta + [single_end: true, mate: 1], fastqs[0]], [meta + [single_end: true, mate: 2], fastqs[1]]]
    }

    // Map, generate BAM with all reads and unmapped reads in FASTQ for downstream
    BOWTIE2_ALIGN(ch_reads_for_bowtie2, ch_bowtie2_index, [[], reference], true, true)
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log)

    ch_unmapped = BOWTIE2_ALIGN.out.fastq.branch { meta, _fastq ->
        mate: meta.mate
        single_end: true
    }

    // Pair the unmapped reads of the two mates again: the reads whose mate aligned to the host are the singletons,
    // which are dropped. The bowtie2 output is in the input order (--reorder), and so is that of repair.sh.
    ch_unmapped_mates = ch_unmapped.mate
        .map { meta, fastq -> [meta.findAll { key, _value -> key != 'mate' } + [single_end: false], meta.mate, fastq] }
        // grouped by the whole meta, which is the same for the two mates once mate is removed. An incomplete group (a
        // field of the meta differing between the mates) is an error, rather than dropped silently by groupTuple
        .groupTuple(by: 0, size: 2, remainder: true)
        .map { meta, mates, fastqs ->
            // the mates are in the order in which they arrived
            if (mates.toSet() != [1, 2].toSet()) {
                error("[zellerlab/flexprofiler] ERROR: host removal: the two mates of ${meta.id} ${meta.run_accession} could not be paired again, found mates ${mates}")
            }
            [meta, [fastqs[mates.indexOf(1)], fastqs[mates.indexOf(2)]]]
        }

    BBMAP_REPAIR_HOSTREMOVAL(ch_unmapped_mates, false)

    // Indexing whole BAM for host removal statistics
    if (!params.skip_hostremoval_qc) {
        SAMTOOLS_INDEX(BOWTIE2_ALIGN.out.bam)
        bam_bai = BOWTIE2_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.index, remainder: true)

        SAMTOOLS_STATS(bam_bai, [[], reference, []])
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats)
        ch_published = ch_published.mix(publishFiles(SAMTOOLS_STATS.out.stats, 'samtools/stats', '*stats'))
    }

    // the logs, the alignments (--save_hostremoval_bam), and the reads not aligned to the host
    // (--save_hostremoval_unmapped): those of the mates of the paired-end reads after BBMAP_REPAIR_HOSTREMOVAL pairs them
    // again, those of the single-end reads from bowtie2
    ch_published = ch_published
        .mix(publishFiles(BOWTIE2_ALIGN.out.log.mix(BBMAP_REPAIR_HOSTREMOVAL.out.log), 'bowtie2/align', '*.log'))
        .mix(params.save_hostremoval_bam ? publishFiles(BOWTIE2_ALIGN.out.bam, 'bowtie2/align', '*.bam') : channel.empty())
    if (params.save_hostremoval_unmapped) {
        ch_published = ch_published
            .mix(publishFiles(ch_unmapped.single_end, 'bowtie2/align', '*.fastq.gz'))
            .mix(publishFiles(BBMAP_REPAIR_HOSTREMOVAL.out.repaired, 'bowtie2/align', '*_repaired.fastq.gz'))
    }
    // the reads given to the profilers (--save_analysis_ready_fastqs), unless their sample has more than one run and
    // the runs are merged (published by the run merging)
    if (params.save_analysis_ready_fastqs) {
        def is_analysis_ready = { meta, _files -> !(params.perform_runmerging && meta.is_multirun) }
        ch_published = ch_published
            .mix(publishFiles(ch_unmapped.single_end.filter(is_analysis_ready), 'analysis_ready_fastqs', '*.fastq.gz'))
            .mix(publishFiles(BBMAP_REPAIR_HOSTREMOVAL.out.repaired.filter(is_analysis_ready), 'analysis_ready_fastqs', '*_repaired.fastq.gz'))
    }

    // the reads left after host removal are counted from the repair.sh log for the paired-end reads, as the
    // bowtie2 logs count each mate separately, and from the bowtie2 log for the single-end reads
    ch_read_count_logs = BOWTIE2_ALIGN.out.log
        .filter { meta, _log_file -> !meta.mate }
        .map { meta, log_file -> [meta + [tool: 'bowtie2'], log_file] }
        .mix(BBMAP_REPAIR_HOSTREMOVAL.out.log.map { meta, log_file -> [meta + [tool: 'bbmap_repair'], log_file] })

    emit:
    reads           = BBMAP_REPAIR_HOSTREMOVAL.out.repaired.mix(ch_unmapped.single_end) // channel: [ val(meta), [ reads ] ]
    read_count_logs = ch_read_count_logs // channel: [ val(meta), log ], with the tool in meta.tool
    mqc             = ch_multiqc_files
    published       = ch_published // channel: [ file, target ], the files to publish (publishFiles)
}
