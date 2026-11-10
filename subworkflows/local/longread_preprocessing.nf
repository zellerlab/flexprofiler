//
// Preprocessing of the long reads, whatever their library strategy: read trimming and filtering of the Oxford Nanopore
// reads, and quality filtering of the PacBio (HiFi) reads
//

include { READS_QC                   } from './reads_qc'

include { LONGREAD_ADAPTERREMOVAL    } from './longread_adapterremoval.nf'
include { LONGREAD_FILTERING         } from './longread_filtering.nf'
include { VSEARCH_FASTQFILTER        } from '../../modules/local/vsearch/fastqfilter/main'
include { BBMAP_BBDUK as BBDUK_PACBIO } from '../../modules/nf-core/bbmap/bbduk/main'

workflow LONGREAD_PREPROCESSING {
    take:
    reads // [ [ meta ], [ reads ] ], long reads, Oxford Nanopore or PacBio (meta.instrument_platform)

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    // logs of the steps, for the read counts, with the tool in meta.tool: [ meta, log ]
    ch_read_count_logs = channel.empty()

    ch_reads_by_platform = reads.branch { meta, _input_reads ->
        nanopore: meta.instrument_platform == 'OXFORD_NANOPORE' && params.perform_longread_qc
        pacbio: meta.instrument_platform == 'PACBIO_SMRT' && params.perform_longread_qc
        other: true
    }

    // PacBio HiFi reads need neither adapter nor quality trimming, but, as in the best practice for full-length 16S
    // (filterAndTrim of DADA2 in Callahan et al. 2019), the reads with more than 10 expected errors (relaxed from 2, see
    // conf/modules.config) or with Ns are removed (vsearch), and then those with a base of quality below 3 (BBDuk, as
    // vsearch has no such filter) (ext.args)
    VSEARCH_FASTQFILTER(ch_reads_by_platform.pacbio)
    ch_versions = ch_versions.mix(VSEARCH_FASTQFILTER.out.versions.first())
    ch_read_count_logs = ch_read_count_logs.mix(VSEARCH_FASTQFILTER.out.log.map { meta, log_file -> [meta + [tool: 'vsearch_fastqfilter'], log_file] })

    BBDUK_PACBIO(VSEARCH_FASTQFILTER.out.reads, [])
    ch_versions = ch_versions.mix(BBDUK_PACBIO.out.versions.first())
    ch_read_count_logs = ch_read_count_logs.mix(BBDUK_PACBIO.out.log.map { meta, log_file -> [meta + [tool: 'bbduk_pacbio'], log_file] })

    if (!params.perform_longread_qc) {
        ch_processed_reads = channel.empty()
    }
    else if (!params.longread_qc_skipadaptertrim && params.longread_qc_skipqualityfilter) {
        LONGREAD_ADAPTERREMOVAL(ch_reads_by_platform.nanopore)
        ch_processed_reads = LONGREAD_ADAPTERREMOVAL.out.reads
        ch_versions = ch_versions.mix(LONGREAD_ADAPTERREMOVAL.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_ADAPTERREMOVAL.out.mqc)
    }
    else if (params.longread_qc_skipadaptertrim && !params.longread_qc_skipqualityfilter) {
        LONGREAD_FILTERING(ch_reads_by_platform.nanopore)
        ch_processed_reads = LONGREAD_FILTERING.out.reads
        ch_versions = ch_versions.mix(LONGREAD_FILTERING.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_FILTERING.out.mqc)
    }
    else {
        LONGREAD_ADAPTERREMOVAL(ch_reads_by_platform.nanopore)
        ch_clipped_reads = LONGREAD_ADAPTERREMOVAL.out.reads.map { meta, clipped_long_reads -> [meta + [single_end: true], clipped_long_reads] }
        LONGREAD_FILTERING(ch_clipped_reads)
        ch_processed_reads = LONGREAD_FILTERING.out.reads
        ch_versions = ch_versions.mix(LONGREAD_ADAPTERREMOVAL.out.versions.first())
        ch_versions = ch_versions.mix(LONGREAD_FILTERING.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_ADAPTERREMOVAL.out.mqc)
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_FILTERING.out.mqc)
    }

    // no reads without --perform_longread_qc
    READS_QC(ch_processed_reads, 'processed')
    ch_versions = ch_versions.mix(READS_QC.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(READS_QC.out.mqc)

    emit:
    // the processed reads are a single file, put in a list as the other reads
    reads           = ch_processed_reads
        .mix(BBDUK_PACBIO.out.reads)
        .map { meta, processed_reads -> [meta, [processed_reads].flatten()] }
        .mix(ch_reads_by_platform.other) // channel: [ val(meta), [ reads ] ]
    read_count_logs = ch_read_count_logs // channel: [ val(meta), log ], with the tool in meta.tool
    versions        = ch_versions // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
