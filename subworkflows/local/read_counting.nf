//
// The read counts after each step, parsed from the logs of the steps, per run or, with run merging, per merged run
//

include { READ_COUNTS  } from '../../modules/local/read_counts/main'
include { publishFiles } from './utils_zellerlab_flexprofiler_pipeline'

workflow READ_COUNTING {
    take:
    logs // [ [ meta ], log ], the logs of the steps, with the tool in meta.tool (the kind of log for read_counts.py)

    main:
    // The [ meta, log ] of each step are split into the logs and the rows of the read counts manifest. Logs of
    // different steps can have the same name, so each is renamed, the manifest giving the new name: a hash of its path,
    // unique, and stable with -resume, the work directories being kept.
    // The database names are defined only for the profiler logs, and are left empty for the other steps
    def required_columns = ['sample', 'run_accession', 'kind']
    def optional_columns = ['db_name']
    def columns = ['file'] + required_columns + optional_columns
    ch_read_count = logs
        .multiMap { meta, log_file ->
            // with run merging, the steps before it (QC, host removal, primer trimming, pair handling) still have the
            // run_accession of each run: the counts of the runs of a multi-run sample go to its merged run, named as by
            // the run merging
            def run = params.perform_runmerging && meta.is_multirun ? "merged_${meta.library_strategy}" : meta.run_accession
            def file_name = "${log_file.toString().md5()}_${log_file.name}"
            def label = [file: file_name, sample: meta.id, run_accession: run, kind: meta.tool, db_name: meta.db_name]
            // only null or blank values count as missing, not falsy ones such as 0
            def missing = required_columns.findAll { column ->
                def value = label[column]
                value == null || value.toString().trim().isEmpty()
            }
            if (missing) {
                error("Read counts: ${missing.join(', ')} not defined for ${log_file}")
            }
            files: [file_name, log_file]
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
    mqc       = READ_COUNTS.out.mqc // channel: path(read_counts_mqc.tsv)
    // read_counts_mqc.tsv goes only to MultiQC
    published = publishFiles(READ_COUNTS.out.tsv, 'read_counts', 'read_counts.tsv') // channel: [ file, target ], the files to publish (publishFiles)
}
