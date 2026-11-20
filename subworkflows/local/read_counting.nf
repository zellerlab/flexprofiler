//
// The read counts after each step, parsed from the logs of the steps, per run or, with run merging, per merged run
//

include { READ_COUNTS       } from '../../modules/local/read_counts/main'
include { publishFiles      } from './utils_zellerlab_flexprofiler_pipeline'
include { readCountColumns  } from './utils_zellerlab_flexprofiler_pipeline'
include { readCountLabel    } from './utils_zellerlab_flexprofiler_pipeline'

workflow READ_COUNTING {
    take:
    logs // [ [ meta ], log ], the logs of the steps, with the tool in meta.tool (the kind of log for read_counts.py)

    main:
    // The [ meta, log ] of each step are split into the logs, under their unique names, and the rows of the read
    // counts manifest (readCountLabel)
    def columns = readCountColumns()
    ch_read_count = logs
        .multiMap { meta, log_file ->
            def label = readCountLabel(meta, log_file)
            files: [label.file, log_file]
            rows: columns.collect { column -> label[column] == null ? '' : label[column] }.join('\t')
        }
    ch_read_count_manifest = ch_read_count.rows
        .collectFile(name: 'read_counts_manifest.tsv', seed: columns.join('\t'), newLine: true)
    ch_read_count_files = ch_read_count.files
        .collectFile { file_name, log_file -> [file_name, log_file] }
        .collect()

    READ_COUNTS(ch_read_count_files, ch_read_count_manifest)

    emit:
    tsv       = READ_COUNTS.out.tsv // channel: path(read_counts.tsv)
    mqc       = READ_COUNTS.out.mqc // channel: [ path(read_counts_mqc.tsv), path(read_counts_plot_mqc.json) ]
    // the MultiQC files go only to MultiQC
    published = publishFiles(
        READ_COUNTS.out.tsv.mix(READ_COUNTS.out.plot),
        'read_counts',
        'read_counts.{tsv,png,pdf}',
    ) // channel: [ file, target ], the files to publish (publishFiles)
}
