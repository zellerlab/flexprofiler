//
// PathSeq-T2T host filtering of the short shotgun reads (--pathseq_t2t), for low-biomass samples, after
// Dohlman et al., Cell 2026 (https://doi.org/10.1016/j.cell.2026.04.015) and https://github.com/abdohlman/pathseq-t2t
// (v0.3.0, commit 51d74430): alignment to GRCh38, selection of the reads not aligned to it (prefilter), quality,
// complexity, host k-mer and duplicate filtering with GATK PathSeqFilterSpark (qcfilter), and removal of the reads
// aligned to T2T-CHM13 (t2tfilter). The differences from upstream are listed in docs/pathseq_t2t.md.
//
// The two stages are kept apart for a later RNA mode, which would change only the first:
// - host pre-alignment: the reads not aligned to GRCh38, those aligned to its decoy contigs, and the primary reads
// - host filtering: PathSeqFilterSpark, then T2T-CHM13
//

include { BWA_INDEX as BWA_INDEX_GRCH38     } from '../../modules/nf-core/bwa/index/main'
include { BWA_INDEX as BWA_INDEX_T2T        } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_GRCH38             } from '../../modules/nf-core/bwa/mem/main'
include { PATHSEQT2T_PREFILTER              } from '../../modules/local/pathseqt2t/prefilter/main'
include { PATHSEQT2T_QCFILTER               } from '../../modules/local/pathseqt2t/qcfilter/main'
include { PATHSEQT2T_T2TFILTER              } from '../../modules/local/pathseqt2t/t2tfilter/main'
include { READS_QC                          } from './reads_qc'

workflow PATHSEQ_T2T {
    take:
    reads           // [ [ meta ], [ reads ] ], short shotgun reads (passed through without --pathseq_t2t), as they were sequenced, per run or sample (run merging)
    grch38_fasta    // /path/to/fasta of GRCh38, [] if its index is given
    grch38_index    // /path/to/directory of the bwa index of GRCh38, [] to build it
    t2t_fasta       // /path/to/fasta of T2T-CHM13, [] if its index is given
    t2t_index       // /path/to/directory of the bwa index of T2T-CHM13, [] to build it
    host_dir        // /path/to/directory of the PathSeq host files (pathseq_host.bfi, pathseq_host.fa.img)
    decoys_bed      // /path/to/bed of the viral decoy contigs of GRCh38, [] to not rescue the reads aligned to them

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // without --pathseq_t2t the reads are passed through unchanged
    ch_reads_filtered = reads
    ch_singletons = channel.empty()
    ch_read_count_logs = channel.empty()

    if (params.pathseq_t2t) {
        if (grch38_index) {
            ch_grch38_index = channel.value([[id: 'grch38'], grch38_index])
        }
        else {
            ch_grch38_index = BWA_INDEX_GRCH38([[id: 'grch38'], grch38_fasta]).index.first()
        }
        if (t2t_index) {
            ch_t2t_index = channel.value([[id: 't2t_chm13'], t2t_index])
        }
        else {
            ch_t2t_index = BWA_INDEX_T2T([[id: 't2t_chm13'], t2t_fasta]).index.first()
        }

        /*
            HOST PRE-ALIGNMENT
        */

        // not sorted: the mates stay together, and nothing downstream needs the coordinate order
        BWA_GRCH38(reads, ch_grch38_index, [[], []], false)

        PATHSEQT2T_PREFILTER(BWA_GRCH38.out.bam, decoys_bed)
        ch_versions = ch_versions.mix(PATHSEQT2T_PREFILTER.out.versions.first())

        /*
            HOST FILTERING
        */

        // the reads aligned to the decoy contigs are given only if there are any
        ch_reads_for_qcfilter = PATHSEQT2T_PREFILTER.out.unaligned
            .join(PATHSEQT2T_PREFILTER.out.decoys, remainder: true)
            .map { meta, unaligned, decoys -> [meta, unaligned, decoys ?: []] }

        PATHSEQT2T_QCFILTER(ch_reads_for_qcfilter, host_dir)
        ch_versions = ch_versions.mix(PATHSEQT2T_QCFILTER.out.versions.first())

        // PathSeqFilterSpark writes no paired or unpaired reads if none is left: the metrics are always written
        ch_reads_for_t2tfilter = PATHSEQT2T_QCFILTER.out.metrics
            .join(PATHSEQT2T_QCFILTER.out.paired, remainder: true)
            .join(PATHSEQT2T_QCFILTER.out.unpaired, remainder: true)
            .map { meta, _metrics, paired, unpaired -> [meta, paired ?: [], unpaired ?: []] }

        PATHSEQT2T_T2TFILTER(ch_reads_for_t2tfilter, ch_t2t_index, decoys_bed)
        ch_versions = ch_versions.mix(PATHSEQT2T_T2TFILTER.out.versions.first())

        // paired-end reads: the pairs, the reads left without their mate (singletons) being kept apart for the profilers
        // that can use them; single-end reads: all the reads are unpaired
        ch_unpaired = PATHSEQT2T_T2TFILTER.out.unpaired.branch { meta, _fastq ->
            single_end: meta.single_end
            singletons: true
        }
        ch_reads_filtered = PATHSEQT2T_T2TFILTER.out.paired
            .mix(ch_unpaired.single_end)
            .map { meta, fastqs -> [meta, [fastqs].flatten()] }

        /*
            QC OF THE FILTERED READS
        */

        // the singletons as single-end reads, named apart; a sample can be left with no reads, which FastQC cannot read
        ch_reads_for_qc = ch_reads_filtered
            .mix(ch_unpaired.singletons.map { meta, fastq -> [meta + [single_end: true, singletons: true], [fastq]] })
            .filter { _meta, fastqs -> fastqs[0].countFastq() > 0 }

        READS_QC(ch_reads_for_qc, 'pathseq_t2t')
        ch_versions = ch_versions.mix(READS_QC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(READS_QC.out.mqc)

        // the reads after each step, for the read counts: the primary reads given to the filtering (flagstat), the reads
        // after each filter of PathSeqFilterSpark (the metrics of the reads not aligned to the host and of those aligned to
        // its decoy contigs, which are summed), and the reads not aligned to T2T-CHM13
        ch_read_count_logs = PATHSEQT2T_PREFILTER.out.flagstat
            .map { meta, flagstat -> [meta + [tool: 'pathseq_flagstat'], flagstat] }
            .mix(
                PATHSEQT2T_QCFILTER.out.metrics
                    .transpose()
                    .map { meta, metrics -> [meta + [tool: 'pathseq_filter'], metrics] }
            )
            .mix(PATHSEQT2T_T2TFILTER.out.log.map { meta, log_file -> [meta + [tool: 'pathseq_t2tfilter'], log_file] })
        ch_singletons = ch_unpaired.singletons
    }

    emit:
    reads           = ch_reads_filtered     // channel: [ val(meta), [ reads ] ], the pairs, or the single-end reads
    singletons      = ch_singletons         // channel: [ val(meta), fastq ], paired-end reads left without their mate
    read_count_logs = ch_read_count_logs    // channel: [ val(meta), log ], with the tool in meta.tool
    versions        = ch_versions           // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
