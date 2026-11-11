//
// Preprocessing of the long reads (PacBio HiFi), whatever their library strategy: quality filtering
//

include { VSEARCH_FASTQFILTER        } from '../../modules/local/vsearch/fastqfilter/main'
include { BBMAP_BBDUK as BBDUK_PACBIO } from '../../modules/nf-core/bbmap/bbduk/main'

workflow LONGREAD_PREPROCESSING {
    take:
    reads // [ [ meta ], [ reads ] ], long reads, PacBio

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    // logs of the steps, for the read counts, with the tool in meta.tool: [ meta, log ]
    ch_read_count_logs = channel.empty()
    // the reads filtered by each platform, mixed in by its block below: [ meta, [ reads ] ]
    ch_processed_reads = channel.empty()

    // the reads are filtered by platform, so that other long-read platforms can be given their own filtering (a branch
    // here and a block below that mixes its reads into ch_processed_reads); the reads without a filtering for their
    // platform, or all of them without --perform_longread_qc, are left as they are
    ch_reads_by_platform = reads.branch { meta, _input_reads ->
        pacbio: params.perform_longread_qc && meta.instrument_platform == 'PACBIO_SMRT'
        other: true
    }

    /*
        PACBIO
    */

    // PacBio HiFi reads need neither adapter nor quality trimming, but, as in the best practice for full-length 16S
    // (filterAndTrim of DADA2 in Callahan et al. 2019), the reads with more than 10 expected errors (relaxed from 2, see
    // conf/modules.config) or with Ns are removed (vsearch), and then those with a base of quality below 3 (BBDuk, as
    // vsearch has no such filter) (ext.args)
    VSEARCH_FASTQFILTER(ch_reads_by_platform.pacbio)
    ch_versions = ch_versions.mix(VSEARCH_FASTQFILTER.out.versions.first())
    ch_read_count_logs = ch_read_count_logs.mix(VSEARCH_FASTQFILTER.out.log.map { meta, log_file -> [meta + [tool: 'vsearch_fastqfilter'], log_file] })

    BBDUK_PACBIO(VSEARCH_FASTQFILTER.out.reads, [])
    ch_read_count_logs = ch_read_count_logs.mix(BBDUK_PACBIO.out.log.map { meta, log_file -> [meta + [tool: 'bbduk_pacbio'], log_file] })
    // the processed reads are a single file, put in a list as the other reads
    ch_processed_reads = ch_processed_reads.mix(BBDUK_PACBIO.out.reads.map { meta, processed_reads -> [meta, [processed_reads].flatten()] })

    emit:
    reads           = ch_processed_reads.mix(ch_reads_by_platform.other) // channel: [ val(meta), [ reads ] ]
    read_count_logs = ch_read_count_logs // channel: [ val(meta), log ], with the tool in meta.tool
    versions        = ch_versions // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
