//
// Preprocessing specific to shotgun reads (library_strategy WGS), after the preprocessing of their read type
// (short/long): redundancy estimation and host removal
//

include { NONPAREIL             } from './nonpareil'
include { SHORTREAD_HOSTREMOVAL } from './shortread_hostremoval'
include { LONGREAD_HOSTREMOVAL  } from './longread_hostremoval'

workflow WGS_PREPROCESSING {
    take:
    reads                  // [ [ meta ], [ reads ] ], WGS runs, short or long reads (meta.type)
    reference              // /path/to/fasta, the host reference, [] if not given
    shortread_index        // the bowtie2 index of the host reference, [] if not given
    longread_index         // the minimap2 index of the host reference, [] if not given

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    // logs of the steps, for the read counts, with the tool in meta.tool: [ meta, log ]
    ch_read_count_logs = channel.empty()

    ch_reads_by_type = reads.branch { meta, _input_reads ->
        shortreads: meta.type == 'short'
        longreads: meta.type == 'long'
    }

    /*
        SUBWORKFLOW: REDUNDANCY ESTIMATION
    */

    if (params.perform_shortread_redundancyestimation) {
        NONPAREIL(ch_reads_by_type.shortreads)
        ch_versions = ch_versions.mix(NONPAREIL.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(NONPAREIL.out.mqc)
    }

    /*
        SUBWORKFLOW: HOST REMOVAL
    */

    if (params.perform_shortread_hostremoval) {
        SHORTREAD_HOSTREMOVAL(ch_reads_by_type.shortreads, reference, shortread_index)
        ch_shortreads_hostremoved = SHORTREAD_HOSTREMOVAL.out.reads
        ch_versions = ch_versions.mix(SHORTREAD_HOSTREMOVAL.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_HOSTREMOVAL.out.mqc)
        ch_read_count_logs = ch_read_count_logs.mix(SHORTREAD_HOSTREMOVAL.out.read_count_logs)
    }
    else {
        ch_shortreads_hostremoved = ch_reads_by_type.shortreads
    }

    if (params.perform_longread_hostremoval) {
        LONGREAD_HOSTREMOVAL(ch_reads_by_type.longreads, reference, longread_index)
        ch_longreads_hostremoved = LONGREAD_HOSTREMOVAL.out.reads
        ch_versions = ch_versions.mix(LONGREAD_HOSTREMOVAL.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_HOSTREMOVAL.out.mqc)
    }
    else {
        ch_longreads_hostremoved = ch_reads_by_type.longreads
    }

    emit:
    reads           = ch_shortreads_hostremoved.mix(ch_longreads_hostremoved) // channel: [ val(meta), [ reads ] ]
    read_count_logs = ch_read_count_logs // channel: [ val(meta), log ], with the tool in meta.tool
    versions        = ch_versions        // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
