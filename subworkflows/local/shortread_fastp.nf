//
// Process short raw reads with FastP
//

include { FASTP as FASTP_SINGLE } from '../../modules/nf-core/fastp/main'
include { FASTP as FASTP_PAIRED } from '../../modules/nf-core/fastp/main'
include { publishFiles          } from './utils_zellerlab_flexprofiler_pipeline'

workflow SHORTREAD_FASTP {
    take:
    reads       // [[meta], [reads]]
    adapterlist

    main:
    ch_multiqc_files = channel.empty()

    ch_input_for_fastp = reads.branch { meta, _input_reads ->
        single: meta.single_end == true
        paired: meta.single_end == false
    }

    // the adapter list goes with the reads ([] if none)
    FASTP_SINGLE(ch_input_for_fastp.single.map { meta, input_reads -> [meta, input_reads, adapterlist] }, false, false, false)
    // Last parameter here turns on merging of PE data
    FASTP_PAIRED(ch_input_for_fastp.paired.map { meta, input_reads -> [meta, input_reads, adapterlist] }, false, false, params.shortread_qc_mergepairs)

    if (params.shortread_qc_mergepairs) {
        ch_fastp_reads_prepped_pe = FASTP_PAIRED.out.reads_merged.map { meta, merged_reads ->
            [meta + [single_end: true], [merged_reads].flatten()]
        }

        ch_fastp_reads_prepped = ch_fastp_reads_prepped_pe.mix(FASTP_SINGLE.out.reads)
    }
    else {
        ch_fastp_reads_prepped = FASTP_PAIRED.out.reads.mix(FASTP_SINGLE.out.reads)
    }
    ch_processed_reads = ch_fastp_reads_prepped

    ch_multiqc_files = ch_multiqc_files.mix(FASTP_SINGLE.out.json)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP_PAIRED.out.json)

    // the reports, and the reads with --save_preprocessed_reads
    ch_fastp_reports = FASTP_SINGLE.out.log.mix(FASTP_SINGLE.out.html, FASTP_SINGLE.out.json, FASTP_PAIRED.out.log, FASTP_PAIRED.out.html, FASTP_PAIRED.out.json)
    ch_fastp_reads = FASTP_SINGLE.out.reads.mix(FASTP_SINGLE.out.reads_fail, FASTP_SINGLE.out.reads_merged, FASTP_PAIRED.out.reads, FASTP_PAIRED.out.reads_fail, FASTP_PAIRED.out.reads_merged)
    ch_published = publishFiles(ch_fastp_reports, 'fastp')
        .mix(params.save_preprocessed_reads ? publishFiles(ch_fastp_reads, 'fastp', '*.fastq.gz') : channel.empty())

    // the reads given to the profilers (--save_analysis_ready_fastqs) are those of fastp for the WGS runs without host
    // removal, unless their sample has more than one run and the runs are merged (published by the run merging)
    if (params.save_analysis_ready_fastqs && !params.perform_shortread_hostremoval) {
        def is_analysis_ready = { meta, _files -> meta.library_strategy == 'WGS' && !(params.perform_runmerging && meta.is_multirun) }
        ch_paired_analysis_ready = params.shortread_qc_mergepairs ? FASTP_PAIRED.out.reads_merged : FASTP_PAIRED.out.reads
        ch_published = ch_published
            .mix(publishFiles(FASTP_SINGLE.out.reads.mix(FASTP_SINGLE.out.reads_fail, FASTP_SINGLE.out.reads_merged).filter(is_analysis_ready), 'analysis_ready_fastqs', '*.fastq.gz'))
            .mix(publishFiles(ch_paired_analysis_ready.filter(is_analysis_ready), 'analysis_ready_fastqs', params.shortread_qc_mergepairs ? '*merged.fastq.gz' : '*.fastp.fastq.gz'))
    }

    emit:
    reads     = ch_processed_reads // channel: [ val(meta), [ reads ] ]
    json      = FASTP_SINGLE.out.json.mix(FASTP_PAIRED.out.json) // channel: [ val(meta), json ]
    mqc       = ch_multiqc_files
    published = ch_published       // channel: [ file, target ], the files to publish (publishFiles)
}
