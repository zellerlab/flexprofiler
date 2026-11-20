/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BBMAP_REPAIR                  } from '../modules/local/bbmap/repair/main'
include { READS_QC                      } from '../subworkflows/local/reads_qc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { SHORTREAD_PREPROCESSING       } from '../subworkflows/local/shortread_preprocessing'
include { LONGREAD_PREPROCESSING        } from '../subworkflows/local/longread_preprocessing'
include { WGS_PREPROCESSING             } from '../subworkflows/local/wgs_preprocessing'
include { AMPLICON_PREPROCESSING        } from '../subworkflows/local/amplicon_preprocessing'
include { RUN_MERGING                   } from '../subworkflows/local/run_merging'
include { PATHSEQ_T2T                   } from '../subworkflows/local/pathseq_t2t'
include { PROFILING                     } from '../subworkflows/local/profiling'
include { STANDARDISATION_PROFILES      } from '../subworkflows/local/standardisation_profiles'
include { READ_COUNTING                 } from '../subworkflows/local/read_counting'
include { MULTIQC_REPORTING             } from '../subworkflows/local/multiqc_reporting'
include { publishFiles                  } from '../subworkflows/local/utils_zellerlab_flexprofiler_pipeline'
include { profileDirectory              } from '../subworkflows/local/utils_zellerlab_flexprofiler_pipeline'
include { toJsonValue                   } from '../subworkflows/local/utils_zellerlab_flexprofiler_pipeline'
include { aggregationOutdir             } from '../subworkflows/local/utils_zellerlab_flexprofiler_pipeline'
include { readCountLabel                } from '../subworkflows/local/utils_zellerlab_flexprofiler_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { UNTAR                         } from '../modules/nf-core/untar/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLEXPROFILER {
    take:
    samplesheet // channel: samplesheet read in from --input
    databases // channel: databases from --databases

    main:

    ch_multiqc_files = channel.empty()

    // Preprocessing auxiliary file input channel preperation
    adapterlist = params.shortread_qc_adapterlist ? file(params.shortread_qc_adapterlist, checkIfExists: true) : []
    ch_reference = params.hostremoval_reference ? file(params.hostremoval_reference, checkIfExists: true) : []
    ch_amplicon_orient_reference = params.amplicon_orient_reference ? file(params.amplicon_orient_reference, checkIfExists: true) : []
    if (params.shortread_hostremoval_index) {
        ch_shortread_reference_index = channel.fromPath(params.shortread_hostremoval_index, checkIfExists: true).map { index -> [[], index] }
    }
    else {
        ch_shortread_reference_index = []
    }

    // PathSeq-T2T host filtering (--pathseq_t2t): a reference (and for STAR its GTF) is needed only to build its index,
    // of the aligner to GRCh38 (--pathseq_t2t_grch38_aligner)
    def pathseq_grch38_index_param = params.pathseq_t2t_grch38_aligner == 'star' ? params.pathseq_t2t_grch38_star_index : params.pathseq_t2t_grch38_bwa_index
    ch_pathseq_grch38_reference = params.pathseq_t2t && !pathseq_grch38_index_param ? file(params.pathseq_t2t_grch38_reference, checkIfExists: true) : []
    ch_pathseq_grch38_index = params.pathseq_t2t && pathseq_grch38_index_param ? file(pathseq_grch38_index_param, checkIfExists: true) : []
    ch_pathseq_grch38_gtf = params.pathseq_t2t && params.pathseq_t2t_grch38_aligner == 'star' && !pathseq_grch38_index_param ? file(params.pathseq_t2t_grch38_gtf, checkIfExists: true) : []
    ch_pathseq_chm13_reference = params.pathseq_t2t && !params.pathseq_t2t_chm13_bwa_index ? file(params.pathseq_t2t_chm13_reference, checkIfExists: true) : []
    ch_pathseq_chm13_index = params.pathseq_t2t && params.pathseq_t2t_chm13_bwa_index ? file(params.pathseq_t2t_chm13_bwa_index, checkIfExists: true) : []
    ch_pathseq_host_dir = params.pathseq_t2t ? file(params.pathseq_host_dir, checkIfExists: true) : []
    ch_pathseq_decoys_bed = params.pathseq_t2t && params.pathseq_decoys_bed ? file(params.pathseq_decoys_bed, checkIfExists: true) : []

    // Separate channels for paired-end and single-end short-read FASTQ, and PacBio FASTQ (the samplesheet is validated in
    // PIPELINE_INITIALISATION: only ILLUMINA, and PACBIO_SMRT single-end AMPLICON runs)
    ch_input = samplesheet
        .map { meta, run_accession, instrument_platform, fastq_1, fastq_2 ->
            meta.run_accession = run_accession
            meta.instrument_platform = instrument_platform
            meta.single_end = !fastq_2
            return [meta, instrument_platform, fastq_1, fastq_2]
        }
        .branch { meta, instrument_platform, fastq_1, fastq_2 ->
            pacbio: instrument_platform == 'PACBIO_SMRT'
            return [meta + [type: "long"], [fastq_1]]
            fastq_pe: !meta.single_end
            return [meta + [type: "short"], [fastq_1, fastq_2]]
            fastq_se: true
            return [meta + [type: "short"], [fastq_1]]
        }

    /*
        MODULE: Run BBmap repair
    */

    if (params.perform_repair_pe_order) {
        BBMAP_REPAIR(ch_input.fastq_pe, channel.value(false))
        ch_fastq_repaired = BBMAP_REPAIR.out.repaired.mix(ch_input.fastq_se)
        ch_repair_published = publishFiles(BBMAP_REPAIR.out.repaired.mix(BBMAP_REPAIR.out.singleton, BBMAP_REPAIR.out.log), 'bbmap')
        ch_repair_read_count_logs = BBMAP_REPAIR.out.log.map { meta, log_file -> [meta + [tool: 'bbmap_repair_pe_order'], log_file] }
    }
    else {
        ch_fastq_repaired = ch_input.fastq_pe.mix(ch_input.fastq_se)
        ch_repair_published = channel.empty()
        ch_repair_read_count_logs = channel.empty()
    }
    
    // the short and long reads as given in the samplesheet (before the pair repair), for the QC of the raw reads
    ch_input_for_fastqc = ch_input.fastq_pe.mix(ch_input.fastq_se, ch_input.pacbio)

    // Validate and decompress databases
    ch_dbs_for_untar = databases.branch { db_meta, db_path ->
        if (!db_meta.db_type) {
            db_meta = db_meta + [db_type: "short;long"]
        }
        untar: db_path.name.endsWith(".tar.gz")
        skip: true
    }
    // Filter the channel to untar only those databases for tools that are selected to be run by the user.
    // Also, to ensure only untar once per file, group together all databases of one file
    ch_inputdb_untar = ch_dbs_for_untar.untar
        .filter { db_meta, _db_path ->
            params["run_${db_meta.tool}"]
        }
        .groupTuple(by: 1)
        .map { meta, dbfile ->
            def new_meta = ['id': dbfile.baseName] + ['meta': meta]
            [new_meta, dbfile]
        }

    // Untar the databases
    UNTAR(ch_inputdb_untar)
    ch_databases_published = params.save_untarred_databases ? publishFiles(UNTAR.out.untar, 'untar/databases') : channel.empty()
    // Spread out the untarred and shared databases
    ch_outputdb_from_untar = UNTAR.out.untar
        .map { meta, db ->
            [meta.meta, db]
        }
        .transpose(by: 0)

    ch_final_dbs = ch_dbs_for_untar.skip
        .mix(ch_outputdb_from_untar)
        .map { db_meta, db ->
            def corrected_db_params = db_meta.db_params ? [db_params: db_meta.db_params] : [db_params: '']
            [db_meta + corrected_db_params, db]
        }


    /*
        SUBWORKFLOW: READ TYPE PREPROCESSING
    */

    // the preprocessing of each read type (short/long), whatever the library strategy
    SHORTREAD_PREPROCESSING(ch_fastq_repaired, adapterlist)

    LONGREAD_PREPROCESSING(ch_input.pacbio)

    /*
        SUBWORKFLOW: LIBRARY STRATEGY PREPROCESSING
    */

    // then the preprocessing of each library strategy (WGS/AMPLICON), whatever the read type
    ch_reads_by_strategy = SHORTREAD_PREPROCESSING.out.reads
        .mix(LONGREAD_PREPROCESSING.out.reads)
        .branch { meta, _reads ->
            wgs: meta.library_strategy == 'WGS'
            amplicon: meta.library_strategy == 'AMPLICON'
            other: true
        }
    // the samplesheet schema allows only WGS and AMPLICON, this guards against the two going out of sync
    ch_reads_by_strategy.other.subscribe { meta, _reads ->
        error("[zellerlab/flexprofiler] ERROR: unexpected library_strategy '${meta.library_strategy}' for ${meta.id}, it must be WGS or AMPLICON")
    }

    WGS_PREPROCESSING(ch_reads_by_strategy.wgs, ch_reference, ch_shortread_reference_index)

    AMPLICON_PREPROCESSING(ch_reads_by_strategy.amplicon, ch_amplicon_orient_reference)

    ch_reads_preprocessed = WGS_PREPROCESSING.out.reads.mix(AMPLICON_PREPROCESSING.out.reads)

    /*
        SUBWORKFLOW: RUN MERGING
    */

    RUN_MERGING(ch_reads_preprocessed)

    /*
        SUBWORKFLOW: PATHSEQ-T2T HOST FILTERING
    */

    PATHSEQ_T2T(
        RUN_MERGING.out.reads,
        ch_pathseq_grch38_reference,
        ch_pathseq_grch38_index,
        ch_pathseq_grch38_gtf,
        ch_pathseq_chm13_reference,
        ch_pathseq_chm13_index,
        ch_pathseq_host_dir,
        ch_pathseq_decoys_bed,
    )
    ch_multiqc_files = ch_multiqc_files.mix(PATHSEQ_T2T.out.mqc.collect { _meta, report -> report }.ifEmpty([]))
    // the paired-end reads left without their mate, for the profilers that can use them
    ch_singletons = PATHSEQ_T2T.out.singletons

    ch_reads_runmerged = PATHSEQ_T2T.out.reads

    /*
        SUBWORKFLOW: READS QC
    */

    // the raw reads, and the reads given to the profilers, with the singletons of PathSeq-T2T as single-end reads,
    // named apart
    READS_QC(
        ch_input_for_fastqc,
        ch_reads_runmerged.mix(ch_singletons.map { meta, fastq -> [meta + [single_end: true, singletons: true], [fastq]] }),
    )
    ch_multiqc_files = ch_multiqc_files.mix(READS_QC.out.mqc.collect { _meta, report -> report }.ifEmpty([]))

    /*
        SUBWORKFLOW: PROFILING
    */

    PROFILING(ch_reads_runmerged, ch_singletons, ch_final_dbs)

    /*
        SUBWORKFLOW: READ COUNTS
    */

    // the logs of the steps, with the tool in meta.tool
    ch_read_count_logs = ch_repair_read_count_logs
        .mix(READS_QC.out.read_count_logs)
        .mix(SHORTREAD_PREPROCESSING.out.read_count_logs)
        .mix(LONGREAD_PREPROCESSING.out.read_count_logs)
        .mix(WGS_PREPROCESSING.out.read_count_logs)
        .mix(AMPLICON_PREPROCESSING.out.read_count_logs)
        .mix(PATHSEQ_T2T.out.read_count_logs)
        .mix(PROFILING.out.read_count_logs)

    // the read counts after each step: not in the --large_scale runs, each of one processing batch, whose logs are
    // published with their manifests instead (AGGREGATION INPUTS below): the run after the last batch makes the table of
    // all the batches (AGGREGATE_PUBLISHED)
    if (!params.large_scale) {
        READ_COUNTING(ch_read_count_logs)
        ch_read_counts_mqc = READ_COUNTING.out.mqc
        ch_read_counts_published = READ_COUNTING.out.published
    }
    else {
        ch_read_counts_mqc = channel.empty()
        ch_read_counts_published = channel.empty()
    }

    /*
        SUBWORKFLOW: PROFILING STANDARDISATION
    */
    // not in the --large_scale runs, each of one processing batch: the run after the last batch merges the profiles of
    // all the batches (AGGREGATE_PUBLISHED)
    if (params.run_profile_standardisation && !params.large_scale) {
        STANDARDISATION_PROFILES(PROFILING.out.profiles, ch_final_dbs)
        ch_merged_profiles_published = STANDARDISATION_PROFILES.out.published
    }
    else {
        ch_merged_profiles_published = channel.empty()
    }

    /*
        AGGREGATION INPUTS
    */

    // --large_scale: per sample, the manifests from which AGGREGATE_PUBLISHED makes the tables of all the samples, as
    // STANDARDISATION_PROFILES and READ_COUNTING do in a run of all of them, one JSON object per line, each with the
    // output directory, which AGGREGATE_PUBLISHED checks against its own:
    // - of its profiles: their meta map, their paths under the output directory (where PROFILING publishes them), and
    //   whether the profiler emits a list of files
    // - of its read count logs: the row of the read counts manifest (readCountLabel); the logs are published under their
    //   unique names
    if (params.large_scale) {
        def outdir = aggregationOutdir()
        ch_profile_manifests_published = PROFILING.out.profiles
            .map { meta, profiles ->
                def files = profiles instanceof List ? profiles : [profiles]
                def entry = [
                    meta: toJsonValue(meta),
                    paths: files.collect { profile -> "${profileDirectory(meta)}/${profile.name}".toString() },
                    list: profiles instanceof List,
                    outdir: outdir,
                ]
                [meta.id, groovy.json.JsonOutput.toJson(entry)]
            }
            // sorted, so that the manifest of a sample is the same in every run
            .collectFile(newLine: true, sort: true) { sample, line -> ["${sample}.jsonl", line] }
            .map { manifest -> [manifest, 'aggregation_inputs/profiles/'] }
        ch_read_count_labels = ch_read_count_logs.map { meta, log_file -> [readCountLabel(meta, log_file), log_file] }
        ch_read_count_manifests_published = ch_read_count_labels
            .map { label, log_file -> [log_file, "aggregation_inputs/read_counts/logs/${label.file}"] }
            .mix(
                ch_read_count_labels
                    .map { label, _log_file -> [label.sample, groovy.json.JsonOutput.toJson(toJsonValue(label + [outdir: outdir]))] }
                    // sorted, so that the manifest of a sample is the same in every run
                    .collectFile(newLine: true, sort: true) { sample, line -> ["${sample}.jsonl", line] }
                    .map { manifest -> [manifest, 'aggregation_inputs/read_counts/'] }
            )
    }
    else {
        ch_profile_manifests_published = channel.empty()
        ch_read_count_manifests_published = channel.empty()
    }

    /*
        SUBWORKFLOW: MULTIQC
    */

    // the reports of the steps for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_PREPROCESSING.out.mqc.collect { _meta, report -> report }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_PREPROCESSING.out.mqc.collect { _meta, report -> report }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(WGS_PREPROCESSING.out.mqc.collect { _meta, report -> report }.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(PROFILING.out.mqc.collect { _meta, report -> report }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_read_counts_mqc)

    if (params.run_profile_standardisation && !params.large_scale) {
        ch_multiqc_files = ch_multiqc_files.mix(STANDARDISATION_PROFILES.out.mqc.collect { _meta, report -> report }.ifEmpty([]))
    }

    MULTIQC_REPORTING(ch_multiqc_files)

    emit:
    multiqc_report = MULTIQC_REPORTING.out.report // channel: /path/to/multiqc_report.html
    // the files to publish of all the stages, for the output block of main.nf: [ file, target ] (publishFiles)
    published      = ch_databases_published.mix(
        ch_repair_published,
        SHORTREAD_PREPROCESSING.out.published,
        LONGREAD_PREPROCESSING.out.published,
        WGS_PREPROCESSING.out.published,
        AMPLICON_PREPROCESSING.out.published,
        RUN_MERGING.out.published,
        PATHSEQ_T2T.out.published,
        READS_QC.out.published,
        PROFILING.out.published,
        ch_merged_profiles_published,
        ch_profile_manifests_published,
        ch_read_count_manifests_published,
        ch_read_counts_published,
        MULTIQC_REPORTING.out.published,
    )
    // the analysis-ready reads of the amplicon runs: [ file, target ], which can be input files
    amplicon_analysis_ready = AMPLICON_PREPROCESSING.out.analysis_ready
}

