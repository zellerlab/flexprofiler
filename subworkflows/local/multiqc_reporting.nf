//
// Collate the software versions and make the MultiQC report
//

include { MULTIQC                } from '../../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from './utils_zellerlab_flexprofiler_pipeline'
include { getPipelineCommit      } from './utils_zellerlab_flexprofiler_pipeline'

workflow MULTIQC_REPORTING {
    take:
    versions      // channel: [ path(versions.yml) ], the versions not emitted to the versions topic
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
    ch_collated_versions = softwareVersionsToYAML(versions.mix(topic_versions.versions_file))
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
            storeDir: "${params.outdir}/pipeline_info",
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

    def summary_params  = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
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
        MULTIQC(
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            [],
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
    else {
        multiqc_report = channel.empty()
    }

    emit:
    report = multiqc_report // channel: /path/to/multiqc_report.html
}
