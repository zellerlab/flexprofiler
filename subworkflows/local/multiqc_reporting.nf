//
// Collate the software versions and make the MultiQC report
//

include { MULTIQC                } from '../../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from './utils_zellerlab_flexprofiler_pipeline'
include { getPipelineCommit      } from './utils_zellerlab_flexprofiler_pipeline'
include { publishFiles           } from './utils_zellerlab_flexprofiler_pipeline'

workflow MULTIQC_REPORTING {
    take:
    multiqc_files // channel: [ path ], the reports of the steps for MultiQC

    main:

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    // The commit of the pipeline code is added to the Workflow section made by the nf-core utils (workflowVersionToYAML),
    // so that the versions file and the MultiQC report record which code was run, also for a local clone (where
    // workflow.commitId is not set and getWorkflowVersion has no commit). The nf-core utils indent the entries of this
    // section with 4 spaces, unlike the 2 of the other sections (dumped by the YAML library or made above): they are
    // re-indented with 2, so that the file is consistent and the added entry does not depend on the indentation of the
    // nf-core utils (a mismatched one would make the YAML invalid)
    def pipeline_commit = getPipelineCommit()
    ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .map { version ->
            if (!version.startsWith('Workflow:')) {
                return version
            }
            def lines = version.readLines()
            def entries = lines.drop(1).collect { line -> "  ${line.trim()}" }
            return ([lines[0]] + entries + ["  ${workflow.manifest.name}_commit: ${pipeline_commit}"]).join('\n')
        }
        .mix(topic_versions_string)
        .collectFile(
            name: 'zellerlab_' + 'flexprofiler_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    def summary_params  = paramsSummaryMap(parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = multiqc_files
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_collated_versions)
        .mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true,
            )
        )

    if (!params.skip_multiqc) {
        // a single input: the reports, the configs (the pipeline's, then the custom one), the logo, and no sample
        // renaming; each list is wrapped so that combine keeps it whole
        ch_multiqc_input = ch_multiqc_files.collect().map { files -> [files] }
            .combine(ch_multiqc_config.concat(ch_multiqc_custom_config).toList().map { configs -> [configs] })
            .combine(ch_multiqc_logo.toList().map { logo -> [logo] })
            .map { files, configs, logo -> [[id: 'multiqc'], files, configs, logo, [], []] }
        MULTIQC(ch_multiqc_input)
        multiqc_report = MULTIQC.out.report.map { _meta, report -> report }.toList()
        ch_multiqc_published = publishFiles(MULTIQC.out.report.mix(MULTIQC.out.data, MULTIQC.out.plots), 'multiqc')
    }
    else {
        multiqc_report = channel.empty()
        ch_multiqc_published = channel.empty()
    }

    emit:
    report    = multiqc_report // channel: /path/to/multiqc_report.html
    // the software versions, in pipeline_info, and the MultiQC report
    published = publishFiles(ch_collated_versions, 'pipeline_info').mix(ch_multiqc_published) // channel: [ file, target ], the files to publish (publishFiles)
}
