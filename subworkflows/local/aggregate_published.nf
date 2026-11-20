//
// --large_scale, all the processing batches done: the tables of all the samples of the samplesheet (merged profiles and
// read counts), made from the per-sample outputs published by the runs of the batches, without processing any sample.
// The manifests of each sample, in <outdir>/aggregation_inputs/, list its profiles (published by PROFILING) and its read
// count logs: the profiles are merged by STANDARDISATION_PROFILES and the logs by READ_COUNTS, as in a run of all the
// samples. The runs of the batches must have had the same --outdir as this one, which their manifests record (checked
// by readAggregationManifests), so that the tables are published with the per-sample outputs they are made from
//

include { UNTAR                    } from '../../modules/nf-core/untar/main'
include { READ_COUNTS              } from '../../modules/local/read_counts/main'
include { STANDARDISATION_PROFILES } from './standardisation_profiles'
include { MULTIQC_REPORTING        } from './multiqc_reporting'
include { publishFiles             } from './utils_zellerlab_flexprofiler_pipeline'
include { profilerEnabled          } from './utils_zellerlab_flexprofiler_pipeline'
include { readAggregationManifests } from './utils_zellerlab_flexprofiler_pipeline'
include { aggregationOutdir        } from './utils_zellerlab_flexprofiler_pipeline'
include { readCountColumns         } from './utils_zellerlab_flexprofiler_pipeline'

workflow AGGREGATE_PUBLISHED {
    take:
    samplesheet // [ meta, run_accession, instrument_platform, fastq_1, fastq_2 ], the runs of the samples to aggregate
    databases   // [ db_meta, db_path ], from --databases

    main:
    def source = aggregationOutdir()

    // the samples of the samplesheet, each once, sorted
    ch_ids = samplesheet
        .map { meta, _run_accession, _instrument_platform, _fastq_1, _fastq_2 -> meta.id }
        .unique()
        .toSortedList { a, b -> a.toString() <=> b.toString() }

    // the profiles of the enabled profilers, as emitted by PROFILING: [ meta, profile ], or [ meta, [ profiles ] ]
    ch_profiles = ch_ids.flatMap { ids ->
        readAggregationManifests('profiles', ids)
            .findAll { entry -> profilerEnabled(entry.meta.tool) }
            .collect { entry ->
                def profiles = entry.paths.collect { profile -> file("${source}/${profile}", checkIfExists: true) }
                [new LinkedHashMap(entry.meta), entry.list ? profiles : profiles.first()]
            }
    }

    // the mOTUs merges need the database, untarred if it is a .tar.gz, as in FLEXPROFILER
    ch_dbs = databases
        .filter { db_meta, _db -> db_meta.tool in ['motus4', 'motus3'] && profilerEnabled(db_meta.tool) }
        .branch { _db_meta, db ->
            untar: db.name.endsWith('.tar.gz')
            skip: true
        }
    UNTAR(ch_dbs.untar.map { db_meta, db -> [[id: db.baseName, meta: db_meta], db] })
    ch_final_dbs = ch_dbs.skip.mix(UNTAR.out.untar.map { meta, db -> [meta.meta, db] })

    STANDARDISATION_PROFILES(ch_profiles, ch_final_dbs)

    // the read counts of all the samples, from the logs and the manifests, which give the file names of the logs
    def columns = readCountColumns()
    ch_read_count_entries = ch_ids.flatMap { ids -> readAggregationManifests('read_counts', ids) }
    ch_read_count_manifest = ch_read_count_entries
        .map { entry -> columns.collect { column -> entry[column] == null ? '' : entry[column].toString() }.join('\t') }
        // sorted, so that the task input is the same in every run, for -resume
        .collectFile(name: 'read_counts_manifest.tsv', seed: columns.join('\t'), newLine: true, sort: true)
    ch_read_count_files = ch_read_count_entries
        .map { entry -> file("${source}/aggregation_inputs/read_counts/logs/${entry.file}", checkIfExists: true) }
        .toSortedList { a, b -> a.name <=> b.name }

    READ_COUNTS(ch_read_count_files, ch_read_count_manifest)

    MULTIQC_REPORTING(
        READ_COUNTS.out.mqc.mix(STANDARDISATION_PROFILES.out.mqc.collect { _meta, report -> report }.ifEmpty([]))
    )

    emit:
    multiqc_report = MULTIQC_REPORTING.out.report // channel: /path/to/multiqc_report.html
    published      = STANDARDISATION_PROFILES.out.published.mix(
        publishFiles(READ_COUNTS.out.tsv.mix(READ_COUNTS.out.plot), 'read_counts', 'read_counts.{tsv,png,pdf}'),
        MULTIQC_REPORTING.out.published,
    ) // channel: [ file, target ], the files to publish (publishFiles)
}
