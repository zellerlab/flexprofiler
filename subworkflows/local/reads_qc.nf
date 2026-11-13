//
// Quality control of the reads with the tool of --preprocessing_qc_tool (FastQC or Falco), of the raw reads and of the
// reads given to the profilers: the stage (meta.qc_stage, raw or processed) names the outputs and their directory (see
// conf/modules.config)
//

include { FASTQC       } from '../../modules/nf-core/fastqc/main'
include { FALCO        } from '../../modules/nf-core/falco/main'
include { publishFiles } from './utils_zellerlab_flexprofiler_pipeline'

workflow READS_QC {
    take:
    raw_reads       // [ [ meta ], [ reads ] ], the reads as given in the samplesheet (before the pair repair)
    processed_reads // [ [ meta ], [ reads ] ], the reads given to the profilers

    main:
    ch_multiqc_files = channel.empty()

    // the raw reads are skipped with --skip_preprocessing_qc, the processed ones with --skip_postprocessing_qc. The files
    // without reads are left out, as FastQC cannot read them (e.g. a sample with all its reads removed by a filter)
    ch_reads = (params.skip_preprocessing_qc ? channel.empty() : raw_reads.map { meta, fastqs -> [meta + [qc_stage: 'raw'], fastqs] })
        .mix(params.skip_postprocessing_qc ? channel.empty() : processed_reads.map { meta, fastqs -> [meta + [qc_stage: 'processed'], fastqs] })
        // a single-end sample has one file, not a list (indexing a path would give its first name element)
        .filter { _meta, fastqs -> hasReads(fastqs instanceof List ? fastqs[0] : fastqs) }

    // only the files read by MultiQC are emitted
    if (params.preprocessing_qc_tool == 'falco') {
        FALCO(ch_reads)
        ch_multiqc_files = FALCO.out.txt
            .transpose()
            .filter { _meta, report -> report.name.endsWith('_data.txt') }
        ch_published = publishFiles(FALCO.out.html.mix(FALCO.out.txt), { meta -> "falco/${meta.qc_stage}" }, '*.{html,txt,zip}')
        // the reads of each file (Total Sequences), for the read counts: falco_raw and falco_processed
        ch_read_count_logs = ch_multiqc_files.map { meta, report -> [meta + [tool: "falco_${meta.qc_stage}"], report] }
    }
    else {
        FASTQC(ch_reads)
        ch_multiqc_files = FASTQC.out.zip
            .transpose()
            .filter { _meta, report -> report.name.endsWith('.zip') }
        ch_published = publishFiles(FASTQC.out.html.mix(FASTQC.out.zip), { meta -> "fastqc/${meta.qc_stage}" }, '*.{html,zip}')
        // the reads of each file (Total Sequences, in the zip), for the read counts: fastqc_raw and fastqc_processed
        ch_read_count_logs = ch_multiqc_files.map { meta, report -> [meta + [tool: "fastqc_${meta.qc_stage}"], report] }
    }

    emit:
    mqc             = ch_multiqc_files   // channel: [ val(meta), report ], one report per channel element; the versions go to the versions topic
    read_count_logs = ch_read_count_logs // channel: [ val(meta), report ], with the tool in meta.tool
    published       = ch_published       // channel: [ file, target ], the files to publish (publishFiles)
}

//
// Whether a gzipped FASTQ has any read: only its first byte is decompressed, so that large files are not read
//
def hasReads(fastq) {
    return new java.util.zip.GZIPInputStream(fastq.newInputStream()).withCloseable { stream -> stream.read() != -1 }
}
