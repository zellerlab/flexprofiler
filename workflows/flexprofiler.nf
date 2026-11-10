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

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Preprocessing auxiliary file input channel preperation
    adapterlist = params.shortread_qc_adapterlist ? file(params.shortread_qc_adapterlist, checkIfExists: true) : []
    ch_reference = params.hostremoval_reference ? file(params.hostremoval_reference, checkIfExists: true) : []
    ch_amplicon_orient_reference = params.amplicon_orient_reference ? file(params.amplicon_orient_reference, checkIfExists: true) : []
    if (params.shortread_hostremoval_index) {
        ch_shortread_reference_index = channel.fromPath(params.shortread_hostremoval_index, checkIfExists: true).map { [[], it] }
    }
    else {
        ch_shortread_reference_index = []
    }
    if (params.longread_hostremoval_index) {
        ch_longread_reference_index = file(params.longread_hostremoval_index, checkIfExists: true)
    }
    else {
        ch_longread_reference_index = []
    }

    // PathSeq-T2T host filtering (--pathseq_t2t): a reference is needed only to build its bwa index
    ch_pathseq_grch38_reference = params.pathseq_t2t && !params.pathseq_t2t_grch38_bwa_index ? file(params.pathseq_t2t_grch38_reference, checkIfExists: true) : []
    ch_pathseq_grch38_index = params.pathseq_t2t && params.pathseq_t2t_grch38_bwa_index ? file(params.pathseq_t2t_grch38_bwa_index, checkIfExists: true) : []
    ch_pathseq_chm13_reference = params.pathseq_t2t && !params.pathseq_t2t_chm13_bwa_index ? file(params.pathseq_t2t_chm13_reference, checkIfExists: true) : []
    ch_pathseq_chm13_index = params.pathseq_t2t && params.pathseq_t2t_chm13_bwa_index ? file(params.pathseq_t2t_chm13_bwa_index, checkIfExists: true) : []
    ch_pathseq_host_dir = params.pathseq_t2t ? file(params.pathseq_host_dir, checkIfExists: true) : []
    ch_pathseq_decoys_bed = params.pathseq_t2t && params.pathseq_decoys_bed ? file(params.pathseq_decoys_bed, checkIfExists: true) : []

    // Validate input files and create separate channels for FASTQ, FASTA, and Nanopore data
    ch_input = samplesheet
        .map { meta, run_accession, instrument_platform, fastq_1, fastq_2, fasta ->
            meta.run_accession = run_accession
            meta.instrument_platform = instrument_platform

            // Define single_end based on the conditions
            meta.single_end = (fastq_1 && !fastq_2 && instrument_platform != 'OXFORD_NANOPORE' && instrument_platform != 'PACBIO_SMRT')

            // Define is_fasta based on the presence of fasta
            meta.is_fasta = fasta ? true : false

            if (!meta.is_fasta && !fastq_1) {
                error("ERROR: Please check input samplesheet: entry `fastq_1` doesn't exist!")
            }
            if (meta.instrument_platform == 'OXFORD_NANOPORE' && fastq_2) {
                error("Error: Please check input samplesheet: for Oxford Nanopore reads entry `fastq_2` should be empty!")
            }
            if (meta.instrument_platform == 'PACBIO_SMRT' && fastq_2) {
                error("Error: Please check input samplesheet: for PacBio reads entry `fastq_2` should be empty!")
            }
            // this is probably not needed as the condition is impossible to reach (meta.single_end is defined as true only if fastq_2 is not present)
            if (meta.single_end && fastq_2) {
                error("Error: Please check input samplesheet: for single-end reads entry `fastq_2` should be empty")
            }
            return [meta, run_accession, instrument_platform, fastq_1, fastq_2, fasta]
        }
        .filter { meta, _run_accession, instrument_platform, _fastq_1, _fastq_2, _fasta ->
            // TODO: this is an initial way to simplify development and focus on adding profilers instead than supporting any data modality
            // Support for long reads and fasta will be enabled on a case-by-case basis depending on need
            def supported = !meta.is_fasta && (instrument_platform == 'ILLUMINA' || (instrument_platform == 'PACBIO_SMRT' && meta.library_strategy == 'AMPLICON'))
            if (!supported) {
                log.warn("[zellerlab/flexprofiler] currently accepts only ILLUMINA FASTQ, and PACBIO_SMRT FASTQ of AMPLICON runs, as input. Skipping sample ${meta.id}, run ${meta.run_accession}.")
            }
            supported
        }
        .branch { meta, _run_accession, instrument_platform, fastq_1, fastq_2, fasta ->
            fastq_pe: !meta.single_end && fastq_2
            return [meta + [type: "short"], [fastq_1, fastq_2]]
            fastq_se: meta.single_end && !fastq_2
            return [meta + [type: "short"], [fastq_1]]
            nanopore: instrument_platform == 'OXFORD_NANOPORE' && !meta.is_fasta
            meta.single_end = true
            return [meta + [type: "long"], [fastq_1]]
            pacbio: instrument_platform == 'PACBIO_SMRT' && !meta.is_fasta
            meta.single_end = true
            return [meta + [type: "long"], [fastq_1]]
            fasta_short: meta.is_fasta && instrument_platform == 'ILLUMINA'
            meta.single_end = true
            return [meta + [type: "short"], [fasta]]
            fasta_long: meta.is_fasta && (instrument_platform == 'OXFORD_NANOPORE' || instrument_platform == 'PACBIO_SMRT')
            meta.single_end = true
            return [meta + [type: "long"], [fasta]]
        }

    /*
        MODULE: Run BBmap repair
    */

    if (params.perform_repair_pe_order) {
        BBMAP_REPAIR(ch_input.fastq_pe, Channel.value(false))
        ch_fastq_repaired = BBMAP_REPAIR.out.repaired.mix(ch_input.fastq_se)
        ch_versions = ch_versions.mix(BBMAP_REPAIR.out.versions.first())
    }
    else {
        ch_fastq_repaired = ch_input.fastq_pe.mix(ch_input.fastq_se)
    }
    
    // Merge ch_fastq_repaired and ch_input.nanopore into a single channel
    ch_input_for_fastqc = ch_fastq_repaired.mix(ch_input.nanopore, ch_input.pacbio)

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
    ch_versions = ch_versions.mix(UNTAR.out.versions.first())

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
        SUBWORKFLOW: QC of the raw reads
    */

    READS_QC(ch_input_for_fastqc, 'raw')
    ch_versions = ch_versions.mix(READS_QC.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(READS_QC.out.mqc.collect { it[1] }.ifEmpty([]))

    /*
        SUBWORKFLOW: READ TYPE PREPROCESSING
    */

    // the preprocessing of each read type (short/long), whatever the library strategy
    SHORTREAD_PREPROCESSING(ch_fastq_repaired, adapterlist)
    ch_versions = ch_versions.mix(SHORTREAD_PREPROCESSING.out.versions)

    LONGREAD_PREPROCESSING(ch_input.nanopore.mix(ch_input.pacbio))
    ch_versions = ch_versions.mix(LONGREAD_PREPROCESSING.out.versions)

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

    WGS_PREPROCESSING(ch_reads_by_strategy.wgs, ch_reference, ch_shortread_reference_index, ch_longread_reference_index)
    ch_versions = ch_versions.mix(WGS_PREPROCESSING.out.versions)

    AMPLICON_PREPROCESSING(ch_reads_by_strategy.amplicon, ch_amplicon_orient_reference)
    ch_versions = ch_versions.mix(AMPLICON_PREPROCESSING.out.versions)

    ch_reads_preprocessed = WGS_PREPROCESSING.out.reads.mix(AMPLICON_PREPROCESSING.out.reads)

    /*
        SUBWORKFLOW: RUN MERGING
    */

    RUN_MERGING(ch_reads_preprocessed)
    ch_versions = ch_versions.mix(RUN_MERGING.out.versions)

    /*
        SUBWORKFLOW: PATHSEQ-T2T HOST FILTERING
    */

    PATHSEQ_T2T(
        RUN_MERGING.out.reads,
        ch_pathseq_grch38_reference,
        ch_pathseq_grch38_index,
        ch_pathseq_chm13_reference,
        ch_pathseq_chm13_index,
        ch_pathseq_host_dir,
        ch_pathseq_decoys_bed,
    )
    ch_versions = ch_versions.mix(PATHSEQ_T2T.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PATHSEQ_T2T.out.mqc.collect { it[1] }.ifEmpty([]))
    // the paired-end reads left without their mate, for the profilers that can use them
    ch_singletons = PATHSEQ_T2T.out.singletons

    // the FASTA inputs are not merged
    // NOTE: FASTA fails with pathseq_t2t, so the mix is only useful when PATHSEQ_T2T is a passthrough
    ch_reads_runmerged = PATHSEQ_T2T.out.reads.mix(ch_input.fasta_short, ch_input.fasta_long)

    /*
        SUBWORKFLOW: PROFILING
    */

    PROFILING(ch_reads_runmerged, ch_singletons, ch_final_dbs)
    ch_versions = ch_versions.mix(PROFILING.out.versions)

    /*
        SUBWORKFLOW: READ COUNTS
    */

    // the read counts after each step, from the logs of the steps, with the tool in meta.tool
    READ_COUNTING(
        SHORTREAD_PREPROCESSING.out.read_count_logs
            .mix(LONGREAD_PREPROCESSING.out.read_count_logs)
            .mix(WGS_PREPROCESSING.out.read_count_logs)
            .mix(AMPLICON_PREPROCESSING.out.read_count_logs)
            .mix(PATHSEQ_T2T.out.read_count_logs)
            .mix(PROFILING.out.read_count_logs)
    )
    ch_versions = ch_versions.mix(READ_COUNTING.out.versions)

    /*
        SUBWORKFLOW: PROFILING STANDARDISATION
    */
    if (params.run_profile_standardisation) {
        STANDARDISATION_PROFILES(PROFILING.out.profiles, ch_final_dbs)
        ch_versions = ch_versions.mix(STANDARDISATION_PROFILES.out.versions)
    }

    /*
        SUBWORKFLOW: MULTIQC
    */

    // the reports of the steps for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_PREPROCESSING.out.mqc.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_PREPROCESSING.out.mqc.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(WGS_PREPROCESSING.out.mqc.collect { it[1] }.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(PROFILING.out.mqc.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(READ_COUNTING.out.mqc)

    if (params.run_profile_standardisation) {
        ch_multiqc_files = ch_multiqc_files.mix(STANDARDISATION_PROFILES.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    MULTIQC_REPORTING(ch_versions, ch_multiqc_files)

    emit:
    multiqc_report = MULTIQC_REPORTING.out.report // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}

