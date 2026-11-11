#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    zellerlab/flexprofiler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/zellerlab/flexprofiler
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FLEXPROFILER            } from './workflows/flexprofiler'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_zellerlab_flexprofiler_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_zellerlab_flexprofiler_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The types of the parameters: Nextflow converts the values given on the command line to them (with the strict syntax,
// the default since Nextflow 26.04, `--perform_shortread_qc false` would otherwise be the string 'false', which is
// true). The defaults are in nextflow.config, which sets them over this block, and the descriptions and allowed values
// in nextflow_schema.json, validated by nf-schema. Every parameter of the schema is declared here, but `help`, which is
// a boolean or the name of a parameter
params {
    input:                           String?
    skip_multiqc:                    Boolean
    multiqc_config:                  String?
    multiqc_title:                   String?
    multiqc_logo:                    String?
    max_multiqc_email_size:          String?
    multiqc_methods_description:     String?
    outdir:                          String?
    publish_dir_mode:                String?
    email:                           String?
    email_on_fail:                   String?
    plaintext_email:                 Boolean
    monochrome_logs:                 Boolean
    help_full:                       Boolean
    show_hidden:                     Boolean
    version:                         Boolean
    pipelines_testdata_base_path:    String?
    trace_report_suffix:             String?
    config_profile_name:             String?
    config_profile_description:      String?
    custom_config_version:           String?
    custom_config_base:              String?
    config_profile_contact:          String?
    config_profile_url:              String?
    validate_params:                 Boolean
    databases:                       String?
    save_untarred_databases:         Boolean
    perform_repair_pe_order:         Boolean
    skip_preprocessing_qc:           Boolean
    skip_postprocessing_qc:          Boolean
    preprocessing_qc_tool:           String?
    perform_shortread_qc:            Boolean
    shortread_qc_skipadaptertrim:    Boolean
    shortread_qc_mergepairs:         Boolean
    shortread_qc_includeunmerged:    Boolean
    shortread_qc_adapter1:           String?
    shortread_qc_adapter2:           String?
    shortread_qc_adapterlist:        String?
    shortread_qc_minlength:          Integer
    shortread_qc_dedup:              Boolean
    perform_longread_qc:             Boolean
    save_preprocessed_reads:         Boolean
    perform_amplicon_primertrimming: Boolean
    amplicon_primer_fwd:             String?
    amplicon_primer_rev:             String?
    amplicon_orient_reference:       String?
    amplicon_pimento_primers_pacbio: String?
    amplicon_pair_handling:          String?
    perform_runmerging:              Boolean
    save_runmerged_reads:            Boolean
    skip_hostremoval_qc:             Boolean
    perform_shortread_hostremoval:   Boolean
    hostremoval_reference:           String?
    shortread_hostremoval_index:     String?
    save_hostremoval_index:          Boolean
    save_hostremoval_bam:            Boolean
    save_hostremoval_unmapped:       Boolean
    pathseq_t2t:                     Boolean
    pathseq_t2t_grch38_reference:    String?
    pathseq_t2t_grch38_bwa_index:    String?
    pathseq_t2t_chm13_reference:     String?
    pathseq_t2t_chm13_bwa_index:     String?
    pathseq_host_dir:                String?
    pathseq_decoys_bed:              String?
    pathseq_min_clipped_read_length: Integer
    pathseq_filter_duplicates:       Boolean
    save_pathseq_t2t_bams:           Boolean
    save_analysis_ready_fastqs:      Boolean
    run_motus4:                      Boolean
    run_motus3:                      Boolean
    motus4_rename_reads:             Boolean
    run_cayman:                      Boolean
    run_gqlib:                       Boolean
    save_bwa_alignments:             Boolean
    run_kraken2:                     Boolean
    run_idtaxa:                      Boolean
    idtaxa_db_names:                 String?
    idtaxa_threshold:                Integer
    idtaxa_ranks:                    String?
    run_profile_standardisation:     Boolean
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        params.databases,
    )

    //
    // WORKFLOW: Run main workflow
    //
    ZELLERLAB_FLEXPROFILER(
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.databases,
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        ZELLERLAB_FLEXPROFILER.out.multiqc_report,
    )

    // the files to publish, as [ file, target ] (publishFiles), see the output block below
    publish:
    results                        = ZELLERLAB_FLEXPROFILER.out.published
    amplicon_analysis_ready_fastqs = ZELLERLAB_FLEXPROFILER.out.amplicon_analysis_ready
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The results, in the output directory (--outdir, outputDir in nextflow.config), published with --publish_dir_mode
// (workflow.output.mode). Each subworkflow decides which files of its steps are published, and in which directory
// (publishFiles in subworkflows/local/utils_zellerlab_flexprofiler_pipeline), as a channel of [ file, target ], the
// target being a directory, with a trailing slash, or the path of the file. The layout is described in docs/output.md
output {
    results {
        path { source, target -> source >> target }
    }
    // the analysis-ready reads of the amplicon runs (--save_analysis_ready_fastqs): those of a run that no step changed
    // are the input reads, outside the work directory, which are published only with includeInputs
    amplicon_analysis_ready_fastqs {
        path { source, target -> source >> target }
        includeInputs true
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow ZELLERLAB_FLEXPROFILER {
    take:
    samplesheet // channel: samplesheet read in from --input
    databases   // channel: databases in from --databases

    main:

    //
    // WORKFLOW: Run pipeline
    //
    FLEXPROFILER(
        samplesheet,
        databases,
    )

    emit:
    multiqc_report          = FLEXPROFILER.out.multiqc_report          // channel: /path/to/multiqc_report.html
    published               = FLEXPROFILER.out.published               // channel: [ file, target ], the files to publish
    amplicon_analysis_ready = FLEXPROFILER.out.amplicon_analysis_ready // channel: [ file, target ], which can be input files
}
