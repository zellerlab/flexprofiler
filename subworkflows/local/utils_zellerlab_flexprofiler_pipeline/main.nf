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

    def samplesheet_rows = samplesheetToList(input, "assets/schema_input.json")
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
// PathSeq-T2T host filtering (--pathseq_t2t): it is for Illumina shotgun reads as they were sequenced, so the other
// inputs and the steps that change the reads before it are errors, rather than being skipped silently
//
def validatePathseqT2t(samplesheet_rows) {
    // rows: [ meta, run_accession, instrument_platform, fastq_1, fastq_2 ]
    def unsupported = samplesheet_rows.findAll { row -> row[2] != 'ILLUMINA' || row[0].library_strategy != 'WGS' }
    if (unsupported) {
        def listed = unsupported.collect { row -> "${row[0].id} ${row[1]} (${row[2]} ${row[0].library_strategy})" }.join(', ')
        error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t accepts only ILLUMINA FASTQ of WGS runs, please remove from the samplesheet: ${listed}")
    }

    // the steps that PathSeq-T2T does itself (QC, host removal)
    def conflicting = ['perform_shortread_qc', 'perform_shortread_hostremoval'].findAll { name -> params[name] }
    if (conflicting) {
        error("[zellerlab/flexprofiler] ERROR: --pathseq_t2t filters the reads as they were sequenced, and cannot be used with ${conflicting.collect { name -> "--${name}" }.join(', ')}: please turn ${conflicting.size() > 1 ? 'them' : 'it'} off")
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
