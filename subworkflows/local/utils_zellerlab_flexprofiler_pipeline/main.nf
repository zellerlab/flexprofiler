//
// Subworkflow with functionality specific to the nf-core/pipeline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message
    databases         //  string: Path to databases
    processing_batch  //  string: the processing_batch of the samplesheet to process (all the runs if null): --processing_batch, or the batch chosen by --large_scale (largeScalePlan)

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m-----------------------------------------------\033[0m-
\033[0;34m  _____    _ _           _       _         \033[0m
\033[0;34m |__  /___| | | ___ _ __| | __ _| |__      \033[0m
\033[0;34m   / // _ \\ | |/ _ \\ '__| |/ _` | '_ \\     \033[0m
\033[0;34m  / /|  __/ | |  __/ |  | | (_| | |_) |    \033[0m
\033[0;34m /____\\___|_|_|\\___|_|  |_|\\__,_|_.__/     \033[0m
                                       
\033[0;35m  zellerlab/flexprofiler ${workflow.manifest.version}\033[0m
\033[0;35m  commit ${getPipelineCommit()}\033[0m
-\033[2m-----------------------------------------------\033[0m-
"""
    after_text = """
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/zellerlab/flexprofiler/blob/master/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null,       // cli_typecast: the default of nf-schema
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Validate and create channel from input file provided through params.input
    //

    // only the runs of the processing batch, if given; the processing_batch column is then dropped, so that it is not
    // in the meta map (which would change the hash of every task, and so what -resume reuses)
    def samplesheet_rows = selectProcessingBatch(samplesheetToList(input, "assets/schema_input.json"), processing_batch)
    validatePacbio(samplesheet_rows)
    if (params.perform_runmerging) {
        validateRunMerging(samplesheet_rows)
    }
    // first, as it rejects the inputs that the other checks are about
    if (params.pathseq_t2t) {
        validatePathseqT2t(samplesheet_rows)
    }
    // the reads per host gene are counted only by STAR, in PathSeq-T2T
    if (params.pathseq_t2t_merge_gene_counts && !(params.pathseq_t2t && params.pathseq_t2t_grch38_aligner == 'star')) {
        error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t_merge_gene_counts merges the reads per host gene counted by STAR, and needs --pathseq_t2t with --pathseq_t2t_grch38_aligner star: please turn it off, or set them")
    }
    validateAmpliconPrimers(samplesheet_rows)
    // Samples with more than one run are published as analysis-ready reads only after run merging. The runs of a
    // sample are merged separately for each library strategy, so they are counted per sample and strategy
    def runs_per_sample = samplesheet_rows.countBy { row -> [row[0].id, row[0].library_strategy] }
    ch_samplesheet = channel.fromList(samplesheet_rows)
        .map { meta, run_accession, instrument_platform, fastq_1, fastq_2 ->
            [meta + [is_multirun: runs_per_sample[[meta.id, meta.library_strategy]] > 1], run_accession, instrument_platform, fastq_1, fastq_2]
        }

    //
    // Validate and create channel from databases file provided through params.databases
    //
    def database_rows = selectDatabases(
        selectDatabases(samplesheetToList(databases, "assets/schema_database.json"), 'idtaxa', params.idtaxa_db_names),
        'kraken2',
        params.kraken2_db_names,
    )
    validateDatabaseSheet(database_rows)
    validateDatabaseMatching(samplesheet_rows, database_rows)
    ch_databases = channel.fromList(database_rows)

    //
    // Validate parameter inputs
    //

    // Bracken re-estimates the Kraken2 reports, with the k-mer distribution in the kraken2 databases. It follows
    // --run_kraken2 when --run_bracken is not set (runBracken), so only an explicit --run_bracken needs Kraken2
    if (params.run_bracken && !params.run_kraken2) {
        error("[zellerlab/flexprofiler] ERROR: --run_bracken re-estimates the Kraken2 reports and needs --run_kraken2: please set it, or turn Bracken off")
    }

    // Preprocessing
    if (params.run_idtaxa && params.amplicon_pair_handling == 'none') {
        error("[zellerlab/flexprofiler] ERROR: IDTAXA classifies single sequences, so --run_idtaxa needs the read pairs of the amplicon reads made into single sequences: please set --amplicon_pair_handling to merged_plus_joined or merged_plus_unmerged")
    }

    if (params.shortread_qc_includeunmerged && !params.shortread_qc_mergepairs) {
        error("ERROR: [zellerlab/flexprofiler] cannot include unmerged reads when merging is not turned on. Please specify --shortread_qc_mergepairs")
    }

    if (params.shortread_qc_adapterlist) {
        def adapterlist = file(params.shortread_qc_adapterlist, checkIfExists: true)
        if (!adapterlist.extension.matches(".*(fa|fasta|fna|fas)")) {
            error("[zellerlab/flexprofiler] ERROR: fastp adapter list requires a `.fasta` format and extension (or fa, fas, fna). Check input: --shortread_qc_adapterlist ${params.shortread_qc_adapterlist}")
        }
    }

    // The reference is needed also when an index is supplied, e.g. for host removal statistics
    if (params.perform_shortread_hostremoval && !params.hostremoval_reference) {
        error("ERROR: [zellerlab/flexprofiler] --perform_shortread_hostremoval requested but no --hostremoval_reference FASTA supplied (required also when providing --shortread_hostremoval_index). Check input.")
    }

    emit:
    samplesheet = ch_samplesheet
    databases   = ch_databases
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report
    large_scale_plan //    map: the plan of a --large_scale run (largeScalePlan), null otherwise

    main:
    summary_params = paramsSummaryMap(parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

        // --large_scale: the batch processed, or the aggregation, is recorded as done only if the run succeeded (checked
        // in recordLargeScaleRun: `workflow` is null in the onComplete closure, its delegate being the script binding)
        if (large_scale_plan) {
            recordLargeScaleRun(large_scale_plan)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// The files to publish, for the output block of main.nf: a channel of [ file, target ], the target being the directory
// of the file under the output directory, with a trailing slash (the file keeps its name). Each subworkflow emits the
// files of its steps (`published`), so that what is published, and where, is decided next to the step that makes it.
//
// ch: a channel of files, or of tuples of a meta map and files (or lists of files, or directories); directory: a
// string, or a closure of the meta map; pattern: a glob matched against the file names (as the publishDir pattern)
//
def publishFiles(ch, directory, pattern = '*') {
    def matcher = java.nio.file.FileSystems.getDefault().getPathMatcher("glob:${pattern}")
    return ch.flatMap { value ->
        def items = value instanceof List ? value : [value]
        def meta = items && items[0] instanceof Map ? items[0] : [:]
        def target = directory instanceof Closure ? directory.call(meta) : directory
        items.flatten()
            .findAll { item -> item instanceof Path && matcher.matches(item.fileName) }
            .collect { file_path -> [file_path, "${target}/"] }
    }
}

//
// The directory, under the output directory, where PROFILING publishes the profiles of a profiler (meta.tool), as
// recorded by the profile manifests of --large_scale
//
def profileDirectory(meta) {
    def directories = [
        motus4:  "motus4/profile/${meta.db_name}",
        motus3:  "motus3/profile/${meta.db_name}",
        cayman:  "cayman/profile/${meta.db_name}/${meta.annot_db_name}",
        gqlib:   "gqlib/profile_genes/${meta.db_name}/${meta.annot_db_name}",
        idtaxa:  "idtaxa/${meta.db_name}",
        kraken2: "kraken2/${meta.db_name}",
        bracken: "bracken/${meta.db_name}",
    ]
    if (!directories.containsKey(meta.tool)) {
        error("[zellerlab/flexprofiler] ERROR: no output directory known for the profiles of '${meta.tool}'")
    }
    return directories[meta.tool].toString()
}

//
// Whether the profiles of a profiler (meta.tool of the profiles) are made in this run: run_<tool>, and for Bracken
// runBracken()
//
def profilerEnabled(tool) {
    return tool == 'bracken' ? runBracken() : params["run_${tool}"] as Boolean
}

//
// A value made fit for JSON: the paths and strings (GStrings) as strings, in maps and lists too
//
def toJsonValue(value) {
    if (value instanceof Map) {
        return value.collectEntries { key, item -> [key.toString(), toJsonValue(item)] }
    }
    if (value instanceof Path || value instanceof CharSequence) {
        return value.toString()
    }
    if (value instanceof Collection) {
        return value.collect { item -> toJsonValue(item) }
    }
    return value
}

//
// The columns of the read counts manifest (read_counts.py), also those of the read count manifests of --large_scale
//
def readCountColumns() {
    return ['file', 'sample', 'run_accession', 'kind', 'db_name']
}

//
// The row of the read counts manifest of a log of a step ([ meta, log ], with the tool in meta.tool), as a map of its
// columns (readCountColumns). Logs of different steps can have the same name, so each is renamed (file): a hash of its
// path, unique, and stable with -resume, the work directories being kept. With run merging, the steps before it (QC,
// host removal, primer trimming, pair handling) still have the run_accession of each run: the counts of the runs of a
// multi-run sample go to its merged run, named as by the run merging. The database names are defined only for the
// profiler logs, and are left empty for the other steps
//
def readCountLabel(meta, log_file) {
    // a step emits one log per [ meta ]: more than one file means that the glob of its output matched another file in
    // its task directory (e.g. a file written there by hand into a cached task, found again by -resume)
    if (!(log_file instanceof Path)) {
        error("[zellerlab/flexprofiler] ERROR: read counts: the ${meta.tool} log of ${meta.id} ${meta.run_accession} is not one file but ${log_file}: remove from the task directory the files that are not outputs of the step, and run again")
    }
    def run = params.perform_runmerging && meta.is_multirun ? "merged_${meta.library_strategy}" : meta.run_accession
    def label = [
        file: "${log_file.toString().md5()}_${log_file.name}".toString(),
        sample: meta.id,
        run_accession: run,
        kind: meta.tool,
        db_name: meta.db_name,
    ]
    // only null or blank values count as missing, not falsy ones such as 0
    def missing = ['sample', 'run_accession', 'kind'].findAll { column ->
        def value = label[column]
        value == null || value.toString().trim().isEmpty()
    }
    if (missing) {
        error("Read counts: ${missing.join(', ')} not defined for ${log_file}")
    }
    return label
}

//
// The output directory (--outdir) as an absolute path, without a trailing slash: recorded in the manifests of
// --large_scale, and checked by AGGREGATE_PUBLISHED against its own
//
def aggregationOutdir() {
    return file(params.outdir).toAbsolutePath().normalize().toString()
}

//
// Whether two output directories are the same: the same absolute path, or the same directory through a symlink
//
def sameOutdir(a, b) {
    if (a == b) {
        return true
    }
    def path_a = file(a)
    def path_b = file(b)
    return path_a.exists() && path_b.exists() && path_a.toRealPath() == path_b.toRealPath()
}

//
// The entries of the per-sample manifests written by the --large_scale runs of the batches, of the samples with these ids:
// <outdir>/aggregation_inputs/<kind>/<sample>.jsonl, one JSON object per line. Errors, each listing all the samples
// concerned:
// - a sample without its manifest (its batch not processed with --large_scale and this --outdir)
// - a sample whose manifest records another output directory than the --outdir of this run (outdir of each entry):
//   the batches and the aggregation must all use the same --outdir, so that the tables are published with the
//   per-sample outputs they are made from
//
def readAggregationManifests(kind, ids) {
    def outdir = aggregationOutdir()
    def directory = "${outdir}/aggregation_inputs/${kind}"
    def shorten = { list -> list.take(20).join(', ') + (list.size() > 20 ? ", ... (${list.size()} in all)" : '') }
    def missing = ids.findAll { id -> !file("${directory}/${id}.jsonl").exists() }
    if (missing) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale: no ${kind} manifest in ${directory}/ for ${missing.size()} of the ${ids.size()} samples of the samplesheet, although their batches are recorded as done: ${shorten.call(missing)}")
    }
    def slurper = new groovy.json.JsonSlurper()
    def entries = ids.collectMany { id ->
        file("${directory}/${id}.jsonl").readLines()
            .findAll { line -> line.trim() }
            .collect { line -> [sample: id, manifest: slurper.parseText(line)] }
    }
    // the output directories recorded, each checked once
    def recorded = entries.collect { entry -> entry.manifest.outdir?.toString() }.unique()
    def other = recorded.findAll { recorded_outdir -> !recorded_outdir || !sameOutdir(recorded_outdir, outdir) }
    if (other) {
        def samples = entries.findAll { entry -> entry.manifest.outdir?.toString() in other }.collect { entry -> entry.sample }.unique()
        def outdirs = other.collect { recorded_outdir -> recorded_outdir ?: 'none recorded' }
        error("[zellerlab/flexprofiler] ERROR: --large_scale: the ${kind} manifests of ${samples.size()} samples were written by runs with another --outdir (${outdirs.join(', ')}) than this one (${outdir}): the runs of the processing batches and the aggregation must use the same --outdir. Samples: ${shorten.call(samples)}")
    }
    return entries.collect { entry -> entry.manifest }
}

//
// --large_scale: the plan of this run, made before anything else (main.nf). The samplesheet is processed one
// processing_batch at a time, one batch per run, so that the task directories of a batch can be deleted once it is done:
// - a batch is done if a run processing it succeeded, which wrote its record (aggregation_inputs/batches/<batch>.json,
//   recordLargeScaleRun) with the samples it had: if its samples changed since, it is processed again. A record written
//   with another --outdir is an error
// - the task directories of the batches done are deleted (--large_scale_cleanup), once per batch: those listed in the
//   traces of all the runs that processed it (aggregation_inputs/batches/<batch>.traces, recorded when each run
//   starts, so also those of the runs that failed or were killed). The traces of the earlier runs are complete, unlike
//   that of the current run, which is why this is done here and not when a run ends. The tasks that prepare a
//   database or reference, shared by the batches, are kept (largeScaleSharedTask)
// - the batch processed is --processing_batch if given, else the batch not done with the fewest runs (then by name);
//   when all the batches are done, the run makes the tables of all the samples instead (AGGREGATE_PUBLISHED)
// - with --large_scale_slurm_resubmit, the run that succeeds submits the sbatch script of its job again, for the next
//   batch (resubmitLargeScaleRun): the script is found here (largeScaleJobScript), before the batch is processed
// Returns [ aggregate: whether to aggregate, batch: the batch to process (null to aggregate), samples: its samples,
//   forced: whether the batch is --processing_batch, resubmit: the script to submit again (--large_scale_slurm_resubmit) ]
//
def largeScalePlan(input, processing_batch) {
    def resubmit = params.large_scale_slurm_resubmit ? largeScaleJobScript() : null
    if (!input || !params.outdir) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale needs --input and --outdir")
    }
    // the run after the last batch merges the profiles of all the batches: checked first, before anything is deleted
    if (!params.run_profile_standardisation) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale merges the profiles of all the samples once all the batches are done, and needs --run_profile_standardisation: please set it (the profile zellerlab_large_scale does, unless --run_profile_standardisation false is given)")
    }
    def rows = samplesheetToList(input, "assets/schema_input.json").collect { row -> samplesheetRow(row) }
    def without_batch = rows.findAll { row -> row.processing_batch == null }
    if (without_batch) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale processes the samplesheet one processing_batch at a time: ${without_batch.size()} runs have no processing_batch, e.g. ${without_batch.take(5).collect { row -> "${row.meta.id} ${row.run_accession}" }.join(', ')}")
    }
    def bad_names = rows.collect { row -> row.processing_batch }.unique().findAll { b -> !(b ==~ /[A-Za-z0-9._-]+/) }
    if (bad_names) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale: the processing_batch names are file names, and can have only letters, digits, '.', '_' and '-': ${bad_names.join(', ')}")
    }
    // the samples of each batch, sorted (all the runs of a sample are in one batch, checked by selectProcessingBatch)
    def batch_rows = rows.groupBy { row -> row.processing_batch }
    def batch_samples = batch_rows.collectEntries { b, b_rows -> [b, b_rows.collect { row -> row.meta.id.toString() }.unique().sort()] }
    if (processing_batch != null && !batch_samples.containsKey(processing_batch.toString())) {
        error("[zellerlab/flexprofiler] ERROR: --processing_batch ${processing_batch}: no run of the samplesheet is in this processing_batch (batches: ${batch_samples.keySet().sort().join(', ')})")
    }
    def forced = processing_batch?.toString()

    def outdir = aggregationOutdir()
    def batches_dir = "${outdir}/aggregation_inputs/batches"
    def done = batch_samples.keySet().findAll { b -> largeScaleBatchDone(batches_dir, b, batch_samples[b], outdir) }.sort()
    // the batch given with --processing_batch is processed again, so its task directories are kept for -resume
    if (params.large_scale_cleanup) {
        done.findAll { b -> b != forced }.each { b -> cleanupLargeScaleBatch(batches_dir, b) }
    }
    def pending = batch_samples.keySet().findAll { b -> !(b in done) }
        .sort { a, b -> batch_rows[a].size() <=> batch_rows[b].size() ?: a <=> b }
    def batch = forced ?: (pending ? pending.first() : null)

    log.info("--large_scale: ${batch_samples.size()} processing batches, ${done.size()} done${done ? " (${done.join(', ')})" : ''}, ${pending.size()} to do${pending ? " (${pending.join(', ')})" : ''}")
    if (batch == null) {
        log.info("--large_scale: all the batches are done: making the tables of all the samples")
        return [aggregate: true, batch: null, samples: batch_samples.values().flatten().sort(), forced: false, resubmit: resubmit]
    }
    log.info("--large_scale: processing batch ${batch} (${batch_samples[batch].size()} samples, ${batch_rows[batch].size()} runs)")
    // the trace of this run, whose task directories are deleted once the batch is done (again, if it was done before);
    // the tables of all the samples, if made before, are out of date
    file(batches_dir).mkdirs()
    file("${batches_dir}/${batch}.traces").append("${largeScaleTraceFile()}\t${workflow.workDir.toAbsolutePath().normalize()}\n")
    def cleaned = file("${batches_dir}/${batch}.cleaned")
    if (cleaned.exists()) {
        cleaned.delete()
    }
    def aggregated = file("${outdir}/aggregation_inputs/aggregated.json")
    if (aggregated.exists()) {
        aggregated.delete()
    }
    return [aggregate: false, batch: batch, samples: batch_samples[batch], forced: forced != null, resubmit: resubmit]
}

//
// The trace file of this run, as set in nextflow.config
//
def largeScaleTraceFile() {
    return "${aggregationOutdir()}/pipeline_info/execution_trace_${params.trace_report_suffix}.txt".toString()
}

//
// Whether a batch is done (--large_scale): its record exists, was written with this --outdir (an error otherwise), and
// has the samples that the batch has now
//
def largeScaleBatchDone(batches_dir, batch, samples, outdir) {
    def record_file = file("${batches_dir}/${batch}.json")
    if (!record_file.exists()) {
        return false
    }
    def record = new groovy.json.JsonSlurper().parseText(record_file.text)
    if (!record.outdir || !sameOutdir(record.outdir.toString(), outdir)) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale: batch ${batch} was processed with another --outdir (${record.outdir}) than this one (${outdir}), according to ${record_file}: all the runs must use the same --outdir")
    }
    if (record.samples.collect { sample -> sample.toString() }.sort() != samples) {
        log.warn("--large_scale: the samples of batch ${batch} changed since it was processed: it is processed again")
        return false
    }
    return true
}

//
// The tasks of a --large_scale run kept by the cleanup: those that prepare a database or reference, the same in all the
// batches, which the next batches reuse from their task directories (-resume): the databases untarred (UNTAR), and the
// indexes built when not given, of the host genome (BOWTIE2_BUILD, --shortread_hostremoval_index) and of the
// PathSeq-T2T references (BWA_INDEX_GRCH38, BWA_INDEX_T2T). The task name is as in the trace: the process, with its
// workflow path, and its tag in parentheses
//
def largeScaleSharedTask(name) {
    return name ==~ /(.*:)?(UNTAR|BOWTIE2_BUILD|BWA_INDEX_GRCH38|BWA_INDEX_T2T)( .*)?/
}

//
// The task directories listed in a trace (its workdir column), of the form <work_dir>/xx/<hash>, but those of the tasks
// shared by the batches (largeScaleSharedTask)
//
def taskDirsOfTrace(trace_path, work_dir) {
    def trace = file(trace_path)
    if (!trace.exists()) {
        log.warn("--large_scale: trace ${trace_path} not found, its task directories are not deleted")
        return []
    }
    def lines = trace.readLines()
    def header = lines ? lines.first().split('\t', -1).toList() : []
    if (!('workdir' in header)) {
        log.warn("--large_scale: trace ${trace_path} has no workdir column, its task directories are not deleted")
        return []
    }
    def prefix = "${work_dir}/".toString()
    return lines.drop(1)
        // each task as a map of its columns
        .collect { line -> [header, line.split('\t', -1).toList()].transpose().collectEntries() }
        .findAll { task -> task.workdir && !largeScaleSharedTask(task.name) }
        .collect { task -> task.workdir }
        .findAll { dir -> dir.startsWith(prefix) && dir.substring(prefix.length()) ==~ /[0-9a-f]{2}\/[0-9a-f]{30}/ }
}

//
// Delete the task directories of a batch done (--large_scale), once: those of the traces of all the runs that processed
// it. The file <batch>.cleaned records that it was done
//
def cleanupLargeScaleBatch(batches_dir, batch) {
    def marker = file("${batches_dir}/${batch}.cleaned")
    def traces = file("${batches_dir}/${batch}.traces")
    if (marker.exists() || !traces.exists()) {
        return
    }
    def dirs = traces.readLines()
        .findAll { line -> line.trim() }
        .collectMany { line ->
            def (trace_path, work_dir) = line.split('\t', -1).toList()
            work_dir ? taskDirsOfTrace(trace_path, work_dir) : []
        }
        .unique()
        .collect { dir -> file(dir) }
        .findAll { dir -> dir.exists() }
    log.info("--large_scale: batch ${batch} is done: deleting its ${dirs.size()} task directories")
    def deleted = dirs.count { dir -> deleteTaskDir(dir) }
    if (deleted < dirs.size()) {
        log.warn("--large_scale: batch ${batch}: ${dirs.size() - deleted} of its ${dirs.size()} task directories were not deleted (see the warnings above)")
    }
    marker.text = "${deleted} task directories deleted on ${java.time.OffsetDateTime.now()}\n"
}

//
// Delete a task directory, without following the symlinks in it, if it really is a task directory of the work
// directory of this run (the real path of -w, however the trace got the task directory). Returns whether it was
// deleted. The checks are on the
// filesystem, not on the path as given (from a trace), so that nothing outside the work directory can be deleted:
// - its real path (all symlinks resolved) is <work_root>/xx/<30 hex characters>: none of its components is a symlink
//   (e.g. an xx directory linked elsewhere), and it is not in another work directory
// - it is a directory with a .command.run, as every Nextflow task directory
// The inputs of a task are symlinks, e.g. to a database or index directory outside the work directory, which must not
// be touched: Path.deleteDir() follows the symlinks to directories and deletes what is in them (it emptied the shared
// host index), `rm -rf` removes the symlinks themselves
//
def deleteTaskDir(dir) {
    def work_root = java.nio.file.Paths.get(workflow.workDir.toString()).toRealPath()
    def path = java.nio.file.Paths.get(dir.toString()).toAbsolutePath().normalize()
    def skip = { reason ->
        log.warn("--large_scale: ${path} not deleted: ${reason}")
        return false
    }
    if (!java.nio.file.Files.isDirectory(path, java.nio.file.LinkOption.NOFOLLOW_LINKS)) {
        return skip.call('not a directory, or a symlink')
    }
    def real = path.toRealPath()
    if (real != path) {
        return skip.call("its path goes through a symlink (real path: ${real})")
    }
    if (real.nameCount < 2 || real.parent.parent != work_root
        || !(real.parent.fileName.toString() ==~ /[0-9a-f]{2}/) || !(real.fileName.toString() ==~ /[0-9a-f]{30}/)) {
        return skip.call("not a task directory of the work directory of this run (${work_root})")
    }
    if (!java.nio.file.Files.isRegularFile(real.resolve('.command.run'), java.nio.file.LinkOption.NOFOLLOW_LINKS)) {
        return skip.call('no .command.run: not a Nextflow task directory')
    }
    def command = ['rm', '-rf', '--', real.toString()]
    def process = new ProcessBuilder(command).redirectErrorStream(true).start()
    def output = process.inputStream.text.trim()
    if (process.waitFor() != 0) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale: `${command.join(' ')}` failed (exit status ${process.exitValue()}): ${output}")
    }
    return true
}

//
// Record the run of a --large_scale plan that succeeded (PIPELINE_COMPLETION): the batch processed is done
// (aggregation_inputs/batches/<batch>.json, with the samples it had), or the tables of all the samples are made
// (aggregation_inputs/aggregated.json). Then, with --large_scale_slurm_resubmit, the next run is submitted if there is
// one: the next batch, or the tables of all the samples. The record of a batch is read back first as the next run will
// (largeScaleBatchDone): if it would not be taken as done, the next run would process the same batch again, and the
// chain would never end, so none is submitted. An error here does not make the run fail (it is in onComplete): it is
// logged, with what to do, and the chain ends
//
def recordLargeScaleRun(plan) {
    if (!workflow.success) {
        return
    }
    def what = plan.aggregate ? 'the tables of all the samples' : "batch ${plan.batch}"
    try {
        def outdir = aggregationOutdir()
        def record = [
            outdir: outdir,
            samples: plan.samples,
            run_name: workflow.runName,
            session_id: workflow.sessionId.toString(),
            completed: java.time.OffsetDateTime.now().toString(),
        ]
        if (plan.aggregate) {
            file("${outdir}/aggregation_inputs/aggregated.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(record)) + '\n'
            log.info("--large_scale: the tables of all the samples are made: nothing left to do")
            return
        }
        // the batch is done only if all the files the aggregation will need are published: once it is recorded, the
        // next run deletes its task directories, and they could no longer be published again
        def unpublished = largeScaleUnpublished(plan.samples)
        if (unpublished) {
            log.error("--large_scale: batch ${plan.batch} is NOT recorded as done, and no next run was submitted: ${unpublished.size()} files listed by its manifests are not in ${outdir}, e.g. ${unpublished.take(5).join(', ')}. Its task directories are kept: fix the publishing, and run the pipeline again (it resumes, from the cache)")
            return
        }
        def batches_dir = "${outdir}/aggregation_inputs/batches"
        file(batches_dir).mkdirs()
        file("${batches_dir}/${plan.batch}.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([batch: plan.batch] + record)) + '\n'
        if (!largeScaleBatchDone(batches_dir, plan.batch, plan.samples, outdir)) {
            log.error("--large_scale: batch ${plan.batch} was recorded (${batches_dir}/${plan.batch}.json), but its record is not read back as done: the next run would process it again. No next run submitted: please report this as a bug")
            return
        }
        log.info("--large_scale: batch ${plan.batch} is done: run the pipeline again for the next batch, or the tables of all the samples once all the batches are done")
        // a batch given with --processing_batch is processed in every run: running again would not move on
        if (plan.resubmit && !plan.forced) {
            resubmitLargeScaleRun(plan.resubmit)
        }
    }
    catch (Exception e) {
        log.error("--large_scale: the run succeeded, but ${what} could NOT be recorded as done (${e.message ?: e}), and no next run was submitted: run the pipeline again (it resumes, from the cache). See .nextflow.log for the stack trace")
        log.debug("--large_scale: recording the run failed", e)
    }
}

//
// --large_scale: the files listed by the manifests of these samples (aggregation_inputs/profiles/ and read_counts/,
// read by readAggregationManifests, as the aggregation will) that are not in the output directory: the profiles (paths
// relative to the output directory, where PROFILING publishes them) and the read count logs (published under
// aggregation_inputs/read_counts/logs/). A sample without its manifests is an error (readAggregationManifests)
//
def largeScaleUnpublished(samples) {
    def outdir = aggregationOutdir()
    def profiles = readAggregationManifests('profiles', samples)
        .collectMany { entry -> entry.paths.collect { path -> "${outdir}/${path}".toString() } }
    def logs = readAggregationManifests('read_counts', samples)
        .collect { entry -> "${outdir}/aggregation_inputs/read_counts/logs/${entry.file}".toString() }
    return (profiles + logs).unique().findAll { path -> !file(path).isFile() }
}

//
// --large_scale_slurm_resubmit: the sbatch script of the SLURM job running the pipeline (the Command of the job, an absolute
// path, as given to sbatch). An error if the executor is not slurm, if not run in a SLURM job, or if the script is not
// found
//
def largeScaleJobScript() {
    def executor = workflow.session.config.navigate('process.executor') ?: 'local'
    if (executor != 'slurm') {
        error("[zellerlab/flexprofiler] ERROR: --large_scale_slurm_resubmit submits the next run as a SLURM job, and needs the slurm executor (process.executor is ${executor}): use a SLURM profile (e.g. zellerlab_shark), or --large_scale_slurm_resubmit false")
    }
    def job_id = System.getenv('SLURM_JOB_ID')
    if (!job_id) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale_slurm_resubmit submits the sbatch script of the job running the pipeline again: the pipeline must be run in a SLURM job (no SLURM_JOB_ID)")
    }
    def command = ['scontrol', 'show', 'job', '-o', job_id]
    def output = ''
    try {
        def process = new ProcessBuilder(command).redirectErrorStream(true).start()
        output = process.inputStream.text
        process.waitFor()
    }
    catch (Exception e) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale_slurm_resubmit: `${command.join(' ')}` failed: ${e.message}")
    }
    def matcher = output =~ /(?:^|\s)Command=(\S+)/
    def script = matcher.find() ? new File(matcher.group(1)) : null
    if (!script || !script.isAbsolute() || !script.isFile()) {
        error("[zellerlab/flexprofiler] ERROR: --large_scale_slurm_resubmit: no sbatch script found for job ${job_id} (`${command.join(' ')}`: ${output.trim()}): submit the pipeline with `sbatch <script>`")
    }
    log.info("--large_scale_slurm_resubmit: the next run will be submitted with ${script}")
    return script.toString()
}

//
// --large_scale_slurm_resubmit: submit the next run, from the launch directory, as a SLURM job that starts once this one has
// ended with exit status 0 (afterok: the Nextflow session is closed, and its lock on the launch directory released). If
// this job is cancelled or fails instead, the next one never starts, and SLURM removes it (--kill-on-invalid-dep). A
// failed submission ends the chain: submit the script again to resume
//
def resubmitLargeScaleRun(script) {
    def job_id = System.getenv('SLURM_JOB_ID')
    def command = ['sbatch', "--dependency=afterok:${job_id}".toString(), '--kill-on-invalid-dep=yes', script.toString()]
    try {
        def process = new ProcessBuilder(command).directory(workflow.launchDir.toFile()).redirectErrorStream(true).start()
        def output = process.inputStream.text.trim()
        if (process.waitFor() == 0) {
            log.info("--large_scale_slurm_resubmit: next run submitted, to start after job ${job_id}: ${output}")
        }
        else {
            log.error("--large_scale_slurm_resubmit: `${command.join(' ')}` failed (exit status ${process.exitValue()}): ${output}. Submit it again to process the next batch")
        }
    }
    catch (Exception e) {
        log.error("--large_scale_slurm_resubmit: `${command.join(' ')}` failed: ${e.message}. Submit it again to process the next batch")
    }
}

//
// A row of the samplesheet (samplesheetToList: [ meta, run_accession, instrument_platform, fastq_1, fastq_2,
// processing_batch ], in the order of assets/schema_input.json) as a map. The empty values are [] (nf-schema): the
// processing_batch is then null, and a string otherwise (the integers are converted)
//
def samplesheetRow(row) {
    def (meta, run_accession, instrument_platform, fastq_1, fastq_2, processing_batch) = row
    return [
        meta: meta,
        run_accession: run_accession,
        instrument_platform: instrument_platform,
        fastq_1: fastq_1,
        fastq_2: fastq_2,
        processing_batch: processing_batch == null || processing_batch == [] ? null : processing_batch.toString(),
    ]
}

//
// Select the runs of a processing batch (the processing_batch column of the samplesheet), all if batch is null, and drop
// the column: rows [ meta, run_accession, instrument_platform, fastq_1, fastq_2, processing_batch ] ->
// [ meta, run_accession, instrument_platform, fastq_1, fastq_2 ]. All the runs of a sample must be in the same batch
// (they are merged, and the tables of a sample are made from all of them)
//
def selectProcessingBatch(rows, batch) {
    def named_rows = rows.collect { row -> samplesheetRow(row) }
    def split = named_rows.groupBy { row -> row.meta.id }
        .collectEntries { sample, sample_rows -> [sample, sample_rows.collect { row -> row.processing_batch }.unique()] }
        .findAll { _sample, batches -> batches.size() > 1 }
    if (split) {
        def shown = split.collect { sample, batches -> "${sample} (${batches.collect { b -> b ?: 'none' }.join(', ')})" }
        error("[zellerlab/flexprofiler] ERROR: all the runs of a sample must be in the same processing_batch of the samplesheet: ${shown.take(20).join('; ')}")
    }
    def selected = named_rows
    if (batch != null) {
        selected = named_rows.findAll { row -> row.processing_batch == batch.toString() }
        if (!selected) {
            def batches = named_rows.collect { row -> row.processing_batch }.unique().findAll { b -> b != null }.sort()
            error("[zellerlab/flexprofiler] ERROR: --processing_batch ${batch}: no run of the samplesheet is in this processing_batch (batches: ${batches ? batches.join(', ') : 'none, the column is empty or missing'})")
        }
    }
    return selected.collect { row -> [row.meta, row.run_accession, row.instrument_platform, row.fastq_1, row.fastq_2] }
}

//
// The commit of the pipeline code that is run: Nextflow knows it only when the pipeline is run from a repository
// (nextflow run zellerlab/flexprofiler -r ...), so for a local clone it is read with git, and marked as dirty if the
// clone has uncommitted changes
//
def getPipelineCommit() {
    if (workflow.commitId) {
        return workflow.commitId
    }
    def commit = runGit(['rev-parse', 'HEAD'])
    if (!commit) {
        return 'unknown'
    }
    def status = runGit(['status', '--porcelain', '--untracked-files=no'])
    return status ? "${commit}-dirty" : commit
}

// The output of a git command run on the pipeline code, or null if it fails (e.g. git is missing or it is not a clone)
def runGit(args) {
    try {
        def proc = (['git', '-C', projectDir.toString()] + args).execute()
        def out = proc.text.trim()
        proc.waitFor()
        return proc.exitValue() == 0 ? out : null
    }
    catch (Exception e) {
        return null
    }
}

// Validate that all runs of each sample can be merged together
//
def validateRunMerging(rows) {
    // rows: [ meta, run_accession, instrument_platform, fastq_1, fastq_2 ]
    // With pair merging, paired-end reads become single-end before run merging, so endedness can differ: with
    // --shortread_qc_mergepairs, and for amplicon reads unless --amplicon_pair_handling is none
    def pairs_merged = params.perform_shortread_qc && params.shortread_qc_mergepairs
    // the WGS and AMPLICON runs of a sample are merged separately
    rows.groupBy { row -> [row[0].id, row[0].library_strategy] }.each { sample_strategy, sample_rows ->
        def (sample, strategy) = sample_strategy
        def check_endedness = !pairs_merged && (strategy != 'AMPLICON' || params.amplicon_pair_handling == 'none')
        def datatypes = sample_rows.collect { row ->
            def datatype = [row[2]]
            if (check_endedness) {
                datatype << (row[4] ? 'paired-end' : 'single-end')
            }
            datatype
        }.unique()
        if (datatypes.size() > 1) {
            error(
                "Please check input samplesheet -> with --perform_runmerging, all runs of a sample must have the same " +
                "instrument platform and, unless pairs are merged (--shortread_qc_mergepairs, or for AMPLICON runs --amplicon_pair_handling), " +
                "endedness (single/paired-end).\n" +
                "Sample '${sample}' (${strategy}) has: ${datatypes.collect { datatype -> datatype.join(' ') }.join('; ')}"
            )
        }
    }
}

//
// Validate the PacBio runs: their reads are preprocessed only as amplicons, and are single-end
//
def validatePacbio(rows) {
    // rows: [ meta, run_accession, instrument_platform, fastq_1, fastq_2 ]
    rows.findAll { row -> row[2] == 'PACBIO_SMRT' }.each { row ->
        if (row[0].library_strategy != 'AMPLICON') {
            error("[zellerlab/flexprofiler] ERROR: PACBIO_SMRT runs are supported only with library_strategy AMPLICON: please check sample ${row[0].id}, run ${row[1]} (library_strategy ${row[0].library_strategy})")
        }
        if (row[4]) {
            error("[zellerlab/flexprofiler] ERROR: PACBIO_SMRT runs are single-end, fastq_2 must be empty: please check sample ${row[0].id}, run ${row[1]}")
        }
    }
}

//
// Validate the primer parameters: the primers are given both or none, and without them the long-read AMPLICON runs
// need the reference to orient the reads with
//
def validateAmpliconPrimers(rows) {
    // rows: [ meta, run_accession, instrument_platform, fastq_1, fastq_2 ]
    if (!params.amplicon_primer_fwd != !params.amplicon_primer_rev) {
        error("Please check the primer parameters -> --amplicon_primer_fwd and --amplicon_primer_rev must be given together.")
    }
    def long_amplicon = rows.any { row -> row[0].library_strategy == 'AMPLICON' && row[2] == 'PACBIO_SMRT' }
    if (params.perform_amplicon_primertrimming && long_amplicon && !params.amplicon_primer_fwd && !params.amplicon_orient_reference) {
        error(
            "Please check the primer parameters -> the long reads of the AMPLICON runs are oriented before finding their " +
            "primers with PIMENTO: give the reference sequences of their gene (e.g. 16S) with --amplicon_orient_reference, " +
            "or give the primers with --amplicon_primer_fwd and --amplicon_primer_rev."
        )
    }
}

//
// Keep only the rows of a tool in the database sheet named in --<tool>_db_names (comma-separated, db_names), all of them
// if not set; the rows of the other tools are kept
//
def selectDatabases(rows, tool, db_names) {
    // rows: [ meta, db_path ]
    if (!db_names) {
        return rows
    }
    def selected = db_names.tokenize(',')
    def available = rows.findAll { meta, _db_path -> meta.tool == tool }.collect { meta, _db_path -> meta.db_name }
    def missing = selected - available
    if (missing) {
        error("Please check --${tool}_db_names -> ${missing.join(', ')} not among the ${tool} rows of the database sheet. Available: ${available.join(', ') ?: 'none'}")
    }
    return rows.findAll { meta, _db_path -> meta.tool != tool || meta.db_name in selected }
}

//
// The library strategy of the reads a profiler is run on: AMPLICON for the amplicon profilers, WGS for the shotgun ones
//
def toolLibraryStrategy(tool) {
    return tool in ['idtaxa'] ? 'AMPLICON' : 'WGS'
}

//
// Whether Bracken is run: --run_bracken if it is set, otherwise with Kraken2 (--run_kraken2)
//
def runBracken() {
    return params.run_bracken != null ? params.run_bracken : params.run_kraken2
}

//
// Validate that the reads and the databases of the enabled profilers match each other. In PROFILING, the reads meet
// the databases with the same read type (short/long, from instrument_platform and db_type) and library strategy (of
// the reads, and of the profiler, toolLibraryStrategy), and what has no match is dropped silently, so it is an error
// here: an enabled profiler without databases, reads that no enabled profiler can run on, and an enabled profiler
// without reads to run on. Without any enabled profiler the pipeline only preprocesses the reads, which is allowed
//
def validateDatabaseMatching(samplesheet_rows, database_rows) {
    // samplesheet_rows: [ meta, run_accession, instrument_platform, fastq_1, fastq_2 ]; database_rows: [ meta, db_path ]
    def enabled = ['motus4', 'motus3', 'cayman', 'gqlib', 'idtaxa', 'kraken2'].findAll { tool -> params["run_${tool}"] }
    if (!enabled) {
        return
    }
    // the read type and library strategy that each enabled profiler can run on, from its databases
    def tool_keys = enabled.collectEntries { tool ->
        def keys = database_rows
            .findAll { meta, _db_path -> meta.tool == tool }
            .collectMany { meta, _db_path -> (meta.db_type ?: 'short;long').tokenize(';').collect { type -> [type, toolLibraryStrategy(tool)] } }
            .unique()
        [tool, keys]
    }

    def no_db = enabled.findAll { tool -> !tool_keys[tool] }
    if (no_db) {
        error("[zellerlab/flexprofiler] ERROR: ${no_db.collect { tool -> "--run_${tool}" }.join(', ')} set, but the database sheet has no ${no_db.join(', ')} database: add one, or turn the profiler off.")
    }

    def all_keys = tool_keys.values().sum()
    def unprofiled = samplesheet_rows.findAll { row -> !(readKey(row) in all_keys) }
    if (unprofiled) {
        def listed = unprofiled.collect { row -> "${row[0].id} ${row[1]} (${readKey(row).join(' ')})" }.join(', ')
        error(
            "[zellerlab/flexprofiler] ERROR: no enabled profiler has a database for these runs, which would be preprocessed but not profiled: ${listed}. " +
            "The ${enabled.join(', ')} databases match: ${all_keys.unique().collect { key -> key.join(' ') }.join(', ')} reads (read type from instrument_platform and the db_type of the database sheet, library strategy of the profiler)."
        )
    }

    def read_keys = samplesheet_rows.collect { row -> readKey(row) }.unique()
    def no_reads = enabled.findAll { tool -> !tool_keys[tool].any { key -> key in read_keys } }
    if (no_reads) {
        error(
            "[zellerlab/flexprofiler] ERROR: no runs of the samplesheet for ${no_reads.collect { tool -> "--run_${tool} (${tool_keys[tool].collect { key -> key.join(' ') }.join(', ')})" }.join(', ')}: " +
            "check the library_strategy and instrument_platform of the samplesheet and the db_type of the database sheet, or turn the profiler off."
        )
    }
}

//
// The read type (short/long, as meta.type and db_type) and library strategy of a samplesheet row
//
def readKey(row) {
    // row: [ meta, run_accession, instrument_platform, fastq_1, fastq_2 ]
    return [row[2] == 'PACBIO_SMRT' ? 'long' : 'short', row[0].library_strategy]
}

//
// Validate the annotation columns, used only by cayman and gqlib, and that the alignment can be shared
//
def validateDatabaseSheet(rows) {
    // rows: [ meta, db_path ]
    def gene_tools = ['cayman', 'gqlib']
    rows.each { meta, _db_path ->
        if (meta.tool in gene_tools && !(meta.annot_db_name && meta.annot_db_path)) {
            error("Please check database sheet -> ${meta.tool} database '${meta.db_name}' needs both annot_db_name and annot_db_path.")
        }
        if (!(meta.tool in gene_tools) && (meta.annot_db_name || meta.annot_db_path)) {
            error("Please check database sheet -> ${meta.tool} database '${meta.db_name}' has annot_db_name or annot_db_path, which are only used by ${gene_tools.join(' and ')}.")
        }
        // PathSeq-T2T classifies with its own Kraken2 options
        if (params.pathseq_t2t && params.run_kraken2 && meta.tool == 'kraken2' && meta.db_params) {
            error("Please check database sheet -> kraken2 database '${meta.db_name}' has db_params, which cannot be used with --pathseq_t2t (its options, `--confidence 0.15 --report-minimizer-data`, are used): please remove them.")
        }
    }
    // the reads are aligned once per db_name, shared by all the annotations and by cayman and gqlib
    rows.findAll { meta, _db_path -> meta.tool in gene_tools }
        .groupBy { meta, _db_path -> meta.db_name }
        .each { db_name, db_rows ->
            def db_paths = db_rows.collect { _meta, db_path -> db_path.toString() }.unique()
            if (db_paths.size() > 1) {
                error("Please check database sheet -> the ${gene_tools.join(' and ')} rows with db_name '${db_name}' must all have the same db_path. Found: ${db_paths.join(', ')}")
            }
        }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def text_seq_qc = [
        "Sequencing quality control with",
        params.preprocessing_qc_tool == "falco" ? "Falco (de Sena Brandine and Smith 2021)." : "FastQC (Andrews 2010).",
    ].join(' ').trim()

    def text_shortread_qc = [
        "Short read preprocessing was performed with fastp (Chen et al. 2018).",
    ].join(' ').trim()


    def text_longread_qc = [
        "The PacBio reads were filtered with VSEARCH (Rognes et al. 2016) and BBDuk (Bushnell 2022)."
    ].join(' ').trim()

    def text_shortreadhostremoval = [
        "Host read removal was performed for short reads with Bowtie2 (Langmead and Salzberg 2012) and SAMtools (Danecek et al. 2021)."
    ].join(' ').trim()

    def text_pathseq_t2t = [
        "Host reads were removed with PathSeq-T2T (Dohlman et al. 2026): the reads were aligned to GRCh38 with " + (params.pathseq_t2t_grch38_aligner == 'star' ? "STAR (Dobin et al. 2013), which also counted the reads per host gene" : "BWA-MEM (Li and Durbin 2009)") + ", the reads not aligned to it and those aligned to its viral decoy contigs selected with SAMtools (Danecek et al. 2021), filtered with GATK PathSeqFilterSpark (Walker et al. 2018), and the reads aligned to T2T-CHM13 (Nurk et al. 2022) removed."
    ].join(' ').trim()

    def text_amplicon_primertrimming = [
        params.amplicon_primer_fwd ? "The primers of amplicon reads were trimmed with cutadapt (Martin 2011)." : "The primers of amplicon reads were detected with PIMENTO (Atallah et al. 2026), the long reads after their orientation with VSEARCH (Rognes et al. 2016), and trimmed with cutadapt (Martin 2011)."
    ].join(' ').trim()

    def text_amplicon_pairmerging = [
        "The read pairs of amplicon reads were merged with VSEARCH (Rognes et al. 2016)."
    ].join(' ').trim()


    def text_classification = [
        "Taxonomic classification or profiling was carried out with:",
        params.run_motus4 ? "mOTUs 4 (Ruscheweyh et al. 2022)," : "",
        params.run_cayman ? "Cayman (Ducarmon et al. 2026)," : "",
        params.run_gqlib ? "gqlib (https://github.com/cschu/gqlib)," : "",
        params.run_idtaxa ? "IDTAXA (Murali et al. 2018) of DECIPHER (Wright 2016)," : "",
        params.run_kraken2 ? "Kraken2 (Wood et al. 2019)," : "",
        runBracken() ? "Bracken (Lu et al. 2017)," : "",
        ".",
    ].join(' ').trim()

    def citation_text = [
        "Tools used in the workflow included:",
        text_seq_qc,
        params.perform_shortread_qc ? text_shortread_qc : "",
        params.perform_longread_qc ? text_longread_qc : "",
        params.perform_shortread_hostremoval ? text_shortreadhostremoval : "",
        params.pathseq_t2t ? text_pathseq_t2t : "",
        params.perform_amplicon_primertrimming ? text_amplicon_primertrimming : "",
        params.amplicon_pair_handling == 'merged_plus_joined' ? text_amplicon_pairmerging : "",
        [
            params.run_motus4,
            params.run_cayman,
            params.run_gqlib,
            params.run_idtaxa,
            params.run_kraken2
        ].any()
            ? text_classification
            : "",
        "Pipeline results statistics were summarised with MultiQC (Ewels et al. 2016).",
    ].join(' ').trim().replaceAll("[,|.] +\\.", ".")

    return citation_text
}

def toolBibliographyText() {
    def text_seq_qc = [
        params.preprocessing_qc_tool == "falco" ? "<li>de Sena Brandine, G., & Smith, A. D. (2021). Falco: high-speed FastQC emulation for quality control of sequencing data. F1000Research, 8(1874), 1874.  <a href=\"https://doi.org/10.12688/f1000research.21142.2\">10.12688/f1000research.21142.2</li>" : "",
        params.preprocessing_qc_tool == "fastqc" ? "<li>Andrews S. (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data, URL: <a href=\"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/\">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a></li>" : "",
    ].join(' ').trim()



    def text_longread_qc = [
        // the PacBio reads (VSEARCH is cited with the amplicon pair merging)
        "<li>Bushnell B. (2022) BBMap, URL:  <a href=\"http://sourceforge.net/projects/bbmap/\">http://sourceforge.net/projects/bbmap/</a></li>",
    ].join(' ').trim()

    def text_shortreadhostremoval = [
        "<li>Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. <a href=\"https://doi.org/10.1038/nmeth.1923\">10.1038/nmeth.1923</a></li>"
    ].join(' ').trim()

    def text_amplicon_primertrimming = [
        params.amplicon_primer_fwd ? "" : "<li>Atallah, C., Richardson, L., Beracochea, M., & Finn, R. D. (2026). PIMENTO: a primer inference toolkit to facilitate large-scale calling of amplicon sequence variants. GigaScience, 15. <a href=\"https://doi.org/10.1093/gigascience/giag083\">10.1093/gigascience/giag083</a></li>",
        "<li>Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10–12. <a href=\"https://doi.org/10.14806/ej.17.1.200\">10.14806/ej.17.1.200</a></li>"
    ].join(' ').trim()

    def text_amplicon_pairmerging = [
        "<li>Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584. <a href=\"https://doi.org/10.7717/peerj.2584\">10.7717/peerj.2584</a></li>"
    ].join(' ').trim()

    def text_pathseq_t2t = [
        "<li>Dohlman, A. B., Mjelle, R., Wood, H. M., Jiang, K., Shumate, A., Lee, I., Piccinno, G., Serna, G., Yakubu, A.-R., Nuciforo, P., Quirke, P., Huttenhower, C., Segata, N., & Meyerson, M. (2026). Biodiversity and biogeography of the multi-kingdom cancer microbiome. Cell, 189, 3817–3837. <a href=\"https://doi.org/10.1016/j.cell.2026.04.015\">10.1016/j.cell.2026.04.015</a></li>",
        "<li>Walker, M. A., Pedamallu, C. S., Ojesina, A. I., Bullman, S., Sharpe, T., Whelan, C. W., & Meyerson, M. (2018). GATK PathSeq: a customizable computational tool for the discovery and identification of microbial sequences in libraries from eukaryotic hosts. Bioinformatics, 34(24), 4287–4289. <a href=\"https://doi.org/10.1093/bioinformatics/bty501\">10.1093/bioinformatics/bty501</a></li>",
        "<li>Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25(14), 1754–1760. <a href=\"https://doi.org/10.1093/bioinformatics/btp324\">10.1093/bioinformatics/btp324</a></li>",
        params.pathseq_t2t_grch38_aligner == 'star' ? "<li>Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15–21. <a href=\"https://doi.org/10.1093/bioinformatics/bts635\">10.1093/bioinformatics/bts635</a></li>" : "",
        "<li>Nurk, S., Koren, S., Rhie, A., et al. (2022). The complete sequence of a human genome. Science, 376(6588), 44–53. <a href=\"https://doi.org/10.1126/science.abj6987\">10.1126/science.abj6987</a></li>",
    ].join(' ').trim()


    def text_classification = [
        params.run_motus4 ? "<li>Ruscheweyh, H.-J., Milanese, A., Paoli, L., Karcher, N., Clayssen, Q., Keller, M. I., Wirbel, J., Bork, P., Mende, D. R., Zeller, G., & Sunagawa, S. (2022). Cultivation-independent genomes greatly expand taxonomic-profiling capabilities of mOTUs across various environments. Microbiome, 10(1), 212. <a href=\"https://doi.org/10.1186/s40168-022-01410-z\">10.1186/s40168-022-01410-z</a></li>" : "",
        params.run_cayman ? "<li>Ducarmon, Q.R., Karcher, N., Giri, S. et al. Cayman enables large-scale analysis of gut microbiome carbohydrate-active enzyme repertoires. Nat Microbiol (2026). <a href=\"https://doi.org/10.1038/s41564-026-02318-2\">https://doi.org/10.1038/s41564-026-02318-2</a></li>" : "",
        params.run_gqlib ? "<li>Schudoma, C. gqlib: gene quantification library (gffquant). <a href=\"https://github.com/cschu/gqlib\">https://github.com/cschu/gqlib</a></li>" : "",
        params.run_idtaxa ? "<li>Murali, A., Bhargava, A., & Wright, E. S. (2018). IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. Microbiome, 6(1), 140. <a href=\"https://doi.org/10.1186/s40168-018-0521-5\">10.1186/s40168-018-0521-5</a></li>" : "",
        params.run_kraken2 ? "<li>Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1), 257. <a href=\"https://doi.org/10.1186/s13059-019-1891-0\">10.1186/s13059-019-1891-0</a></li>" : "",
        runBracken() ? "<li>Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104. <a href=\"https://doi.org/10.7717/peerj-cs.104\">10.7717/peerj-cs.104</a></li>" : "",
        params.run_idtaxa ? "<li>Wright, E. S. (2016). Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R. The R Journal, 8(1), 352–359. <a href=\"https://doi.org/10.32614/RJ-2016-025\">10.32614/RJ-2016-025</a></li>" : "",
    ].join(' ').trim()

    def text_extras = [
        params.perform_shortread_qc ? "<li>Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics , 34(17), i884–i890. <a href=\"https://doi.org/10.1093/bioinformatics/bty560\">10.1093/bioinformatics/bty560</a></li>" : "",
        params.perform_shortread_hostremoval ? "<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2). <a href=\"https://doi.org/10.1093/gigascience/giab008\">10.1093/gigascience/giab008</a></li>" : "",
    ].join(' ').trim()

    def reference_text = [
        text_seq_qc,
        params.perform_longread_qc ? text_longread_qc : "",
        params.perform_shortread_hostremoval ? text_shortreadhostremoval : "",
        params.pathseq_t2t ? text_pathseq_t2t : "",
        params.perform_amplicon_primertrimming ? text_amplicon_primertrimming : "",
        // VSEARCH: pair merging, orientation of the long amplicon reads, and QC of the PacBio reads
        params.amplicon_pair_handling == 'merged_plus_joined' || (params.perform_amplicon_primertrimming && !params.amplicon_primer_fwd) || params.perform_longread_qc ? text_amplicon_pairmerging : "",
        text_extras,
        text_classification,
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. <a href=\"https:/doi.org/10.1093/bioinformatics/btw354\">10.1093/bioinformatics/btw354.</a></li>",
    ].join(' ').trim().replaceAll("[,|.] +\\.", ".")

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]

    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // meta["tool_citations"] = ""
    // meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// PathSeq-T2T host filtering (--pathseq_t2t): it is for Illumina shotgun reads, so the other inputs, and the steps that
// would do again what it does or break its pairs, are errors, rather than being skipped silently
//
def validatePathseqT2t(samplesheet_rows) {
    // rows: [ meta, run_accession, instrument_platform, fastq_1, fastq_2 ]
    def unsupported = samplesheet_rows.findAll { row -> row[2] != 'ILLUMINA' || row[0].library_strategy != 'WGS' }
    if (unsupported) {
        def listed = unsupported.collect { row -> "${row[0].id} ${row[1]} (${row[2]} ${row[0].library_strategy})" }.join(', ')
        error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t accepts only ILLUMINA FASTQ of WGS runs, please remove from the samplesheet: ${listed}")
    }

    // the steps that PathSeq-T2T does itself (host removal, deduplication over the whole sample), and the merging of the
    // pairs, which would leave no pairs for its proper-pair filters. The adapter trimming of fastp (--perform_shortread_qc)
    // can be used
    def conflicting = ['perform_shortread_hostremoval']
        .plus(params.perform_shortread_qc ? ['shortread_qc_dedup', 'shortread_qc_mergepairs'] : [])
        .findAll { name -> params[name] }
    if (conflicting) {
        error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t removes the host reads and the duplicates itself, and filters read pairs, so it cannot be used with ${conflicting.collect { name -> "--${name}" }.join(', ')}: please turn ${conflicting.size() > 1 ? 'them' : 'it'} off")
    }

    // the gene catalogue profilers are not supported on the reads filtered by PathSeq-T2T
    def unsupported_profilers = ['run_cayman', 'run_gqlib'].findAll { name -> params[name] }
    if (unsupported_profilers) {
        error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t cannot be used with ${unsupported_profilers.collect { name -> "--${name}" }.join(', ')}: please turn ${unsupported_profilers.size() > 1 ? 'them' : 'it'} off")
    }

    if (params.pathseq_t2t_grch38_aligner == 'star') {
        if (params.pathseq_t2t_grch38_star_index) {
            // the reads per gene need the genes in the index
            def star_index = file(params.pathseq_t2t_grch38_star_index, checkIfExists: true)
            if (!star_index.resolve('exonGeTrInfo.tab').exists()) {
                error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t_grch38_star_index ${params.pathseq_t2t_grch38_star_index} has no exonGeTrInfo.tab: the STAR index must be built with a GTF (--sjdbGTFfile) to count the reads per gene, please rebuild it, or give the FASTA and GTF of GRCh38 (--pathseq_t2t_grch38_reference, --pathseq_t2t_grch38_gtf) to build it")
            }
        }
        else if (!params.pathseq_t2t_grch38_reference || !params.pathseq_t2t_grch38_gtf) {
            error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t with --pathseq_t2t_grch38_aligner star needs GRCh38: please give its STAR index (--pathseq_t2t_grch38_star_index), or its FASTA and GTF to build it (--pathseq_t2t_grch38_reference and --pathseq_t2t_grch38_gtf)")
        }
    }
    else if (!params.pathseq_t2t_grch38_reference && !params.pathseq_t2t_grch38_bwa_index) {
        error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t needs GRCh38: please give its bwa index (--pathseq_t2t_grch38_bwa_index) or its FASTA to build it (--pathseq_t2t_grch38_reference)")
    }
    if (!params.pathseq_t2t_chm13_reference && !params.pathseq_t2t_chm13_bwa_index) {
        error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t needs T2T-CHM13: please give its bwa index (--pathseq_t2t_chm13_bwa_index) or its FASTA to build it (--pathseq_t2t_chm13_reference)")
    }
    if (!params.pathseq_host_dir) {
        error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t needs the PathSeq host files (pathseq_host.bfi and pathseq_host.fa.img): please give their directory (--pathseq_host_dir)")
    }
    def host_dir = file(params.pathseq_host_dir, checkIfExists: true)
    def missing = ['pathseq_host.bfi', 'pathseq_host.fa.img'].findAll { name -> !host_dir.resolve(name).exists() }
    if (missing) {
        error("[zellerlab/flexprofiler] ERROR: --pathseq_host_dir ${params.pathseq_host_dir} has no ${missing.join(' and ')}")
    }
}
