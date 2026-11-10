//
// Preprocessing of the short reads, whatever their library strategy: read trimming and merging, and complexity
// filtering
//


include { SHORTREAD_FASTP               } from './shortread_fastp'
include { SHORTREAD_ADAPTERREMOVAL      } from './shortread_adapterremoval'
include { SHORTREAD_COMPLEXITYFILTERING } from './shortread_complexityfiltering'
include { READS_QC                      } from './reads_qc'

workflow SHORTREAD_PREPROCESSING {
    take:
    reads //  [ [ meta ], [ reads ] ]
    adapterlist // file

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_fastp_json = channel.empty()

    if (params.perform_shortread_qc) {
        if (params.shortread_qc_tool == "fastp") {
            SHORTREAD_FASTP(reads, adapterlist)
            ch_processed_reads = SHORTREAD_FASTP.out.reads
            ch_versions = ch_versions.mix(SHORTREAD_FASTP.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_FASTP.out.mqc)
            ch_fastp_json = SHORTREAD_FASTP.out.json
        }
        else if (params.shortread_qc_tool == "adapterremoval") {
            SHORTREAD_ADAPTERREMOVAL(reads, adapterlist)
            ch_processed_reads = SHORTREAD_ADAPTERREMOVAL.out.reads
            ch_versions = ch_versions.mix(SHORTREAD_ADAPTERREMOVAL.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_ADAPTERREMOVAL.out.mqc)
        }
        else {
            ch_processed_reads = reads
        }

        READS_QC(ch_processed_reads, 'processed')
        ch_versions = ch_versions.mix(READS_QC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(READS_QC.out.mqc)
    }
    else {
        ch_processed_reads = reads
    }

    // fastp complexity filtering is activated via modules.conf in SHORTREAD_FASTP
    if (params.perform_shortread_complexityfilter && params.shortread_complexityfilter_tool != 'fastp') {
        SHORTREAD_COMPLEXITYFILTERING(ch_processed_reads)
        ch_filtered_reads = SHORTREAD_COMPLEXITYFILTERING.out.reads
        ch_versions = ch_versions.mix(SHORTREAD_COMPLEXITYFILTERING.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_COMPLEXITYFILTERING.out.mqc)
    }
    else {
        ch_filtered_reads = ch_processed_reads
    }

    emit:
    reads           = ch_filtered_reads // channel: [ val(meta), [ reads ] ]
    json            = ch_fastp_json // channel: [ val(meta), json ], only with fastp
    read_count_logs = ch_fastp_json.map { meta, json -> [meta + [tool: 'fastp'], json] } // channel: [ val(meta), log ], with the tool in meta.tool
    versions        = ch_versions // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
