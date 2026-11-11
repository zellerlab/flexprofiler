//
// Preprocessing specific to shotgun reads (library_strategy WGS), after the preprocessing of their read type: host
// removal. The WGS runs are all short reads, as PacBio is accepted only for AMPLICON runs
//

include { SHORTREAD_HOSTREMOVAL } from './shortread_hostremoval'

workflow WGS_PREPROCESSING {
    take:
    reads                  // [ [ meta ], [ reads ] ], WGS runs, short reads
    reference              // /path/to/fasta, the host reference, [] if not given
    shortread_index        // the bowtie2 index of the host reference, [] if not given

    main:
    ch_multiqc_files = channel.empty()
    // logs of the steps, for the read counts, with the tool in meta.tool: [ meta, log ]
    ch_read_count_logs = channel.empty()

    /*
        SUBWORKFLOW: HOST REMOVAL
    */

    if (params.perform_shortread_hostremoval) {
        SHORTREAD_HOSTREMOVAL(reads, reference, shortread_index)
        ch_reads_hostremoved = SHORTREAD_HOSTREMOVAL.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_HOSTREMOVAL.out.mqc)
        ch_read_count_logs = ch_read_count_logs.mix(SHORTREAD_HOSTREMOVAL.out.read_count_logs)
    }
    else {
        ch_reads_hostremoved = reads
    }

    emit:
    reads           = ch_reads_hostremoved // channel: [ val(meta), [ reads ] ]
    read_count_logs = ch_read_count_logs   // channel: [ val(meta), log ], with the tool in meta.tool
    mqc             = ch_multiqc_files
}
