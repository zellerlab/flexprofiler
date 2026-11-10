//
// Quality control of the reads with the tool of --preprocessing_qc_tool (FastQC or Falco), at any stage of the
// pipeline: the stage names the outputs and their directory (see conf/modules.config)
//

include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { FALCO  } from '../../modules/nf-core/falco/main'

workflow READS_QC {
    take:
    reads   // [ [ meta ], [ reads ] ]
    stage   // val: 'raw', 'processed' or 'pathseq_t2t'

    main:
    ch_reads = reads.map { meta, fastqs -> [meta + [qc_stage: stage], fastqs] }
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // the raw reads are skipped with --skip_preprocessing_qc, those of the later stages with --skip_postprocessing_qc;
    // only the files read by MultiQC are emitted
    if (stage == 'raw' ? params.skip_preprocessing_qc : params.skip_postprocessing_qc) {
        // no QC
    }
    else if (params.preprocessing_qc_tool == 'falco') {
        FALCO(ch_reads)
        ch_versions = FALCO.out.versions.first()
        ch_multiqc_files = FALCO.out.txt
            .transpose()
            .filter { _meta, report -> report.name.endsWith('_data.txt') }
    }
    else {
        FASTQC(ch_reads)
        ch_versions = FASTQC.out.versions.first()
        ch_multiqc_files = FASTQC.out.zip
            .transpose()
            .filter { _meta, report -> report.name.endsWith('.zip') }
    }

    emit:
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files     // channel: [ val(meta), report ], one report per channel element
}
