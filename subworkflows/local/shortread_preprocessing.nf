//
// Preprocessing of the short reads, whatever their library strategy: read trimming and merging
//


include { SHORTREAD_FASTP               } from './shortread_fastp'

workflow SHORTREAD_PREPROCESSING {
    take:
    reads //  [ [ meta ], [ reads ] ]
    adapterlist // file

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_fastp_json = channel.empty()

    if (params.perform_shortread_qc) {
        SHORTREAD_FASTP(reads, adapterlist)
        ch_processed_reads = SHORTREAD_FASTP.out.reads
        ch_versions = ch_versions.mix(SHORTREAD_FASTP.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_FASTP.out.mqc)
        ch_fastp_json = SHORTREAD_FASTP.out.json
    }
    else {
        ch_processed_reads = reads
    }

    emit:
    reads           = ch_processed_reads // channel: [ val(meta), [ reads ] ]
    json            = ch_fastp_json // channel: [ val(meta), json ], only with fastp
    read_count_logs = ch_fastp_json.map { meta, json -> [meta + [tool: 'fastp'], json] } // channel: [ val(meta), log ], with the tool in meta.tool
    versions        = ch_versions // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
