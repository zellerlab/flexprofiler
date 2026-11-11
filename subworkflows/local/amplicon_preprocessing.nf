//
// Preprocessing specific to amplicon reads (library_strategy AMPLICON), after the preprocessing of their read type
// (short/long): primer trimming, and the handling of the read pairs, so that each amplicon is a single sequence
//

include { AMPLICON_PRIMERTRIMMING             } from './amplicon_primertrimming'
include { VSEARCH_MERGEJOIN                   } from '../../modules/local/vsearch/mergejoin/main'
include { FASTP as FASTP_AMPLICON_MERGE       } from '../../modules/nf-core/fastp/main'

workflow AMPLICON_PREPROCESSING {
    take:
    reads            // [ [ meta ], [ reads ] ], AMPLICON runs, short (paired-end or single-end) or long reads
    orient_reference // the reference to orient the long reads with (FASTA), [] if not given

    main:
    // logs of the steps, for the read counts, with the tool in meta.tool: [ meta, log ]
    ch_read_count_logs = channel.empty()

    /*
        SUBWORKFLOW: PRIMER TRIMMING
    */

    if (params.perform_amplicon_primertrimming) {
        AMPLICON_PRIMERTRIMMING(reads, orient_reference)
        ch_reads_primertrimmed = AMPLICON_PRIMERTRIMMING.out.reads
        ch_read_count_logs = ch_read_count_logs.mix(AMPLICON_PRIMERTRIMMING.out.read_count_logs)
    }
    else {
        ch_reads_primertrimmed = reads
    }

    /*
        PAIR HANDLING
    */

    // The amplicon profilers classify single sequences: the pairs of paired-end reads are merged if they overlap, and
    // the others joined with a padding of Ns (merged_plus_joined, vsearch) or kept as two separate mates
    // (merged_plus_unmerged, fastp). With none, the pairs are kept as they are
    ch_reads_for_pairhandling = ch_reads_primertrimmed.branch { meta, _input_reads ->
        pe: !meta.single_end && params.amplicon_pair_handling != 'none'
        other: true
    }

    if (params.amplicon_pair_handling == 'merged_plus_joined') {
        VSEARCH_MERGEJOIN(ch_reads_for_pairhandling.pe)
        ch_read_count_logs = ch_read_count_logs.mix(VSEARCH_MERGEJOIN.out.log.map { meta, log_file -> [meta + [tool: 'vsearch'], log_file] })
        ch_reads_pairhandled = VSEARCH_MERGEJOIN.out.reads
    }
    else if (params.amplicon_pair_handling == 'merged_plus_unmerged') {
        FASTP_AMPLICON_MERGE(ch_reads_for_pairhandling.pe.map { meta, input_reads -> [meta, input_reads, []] }, false, false, true)
        ch_read_count_logs = ch_read_count_logs.mix(FASTP_AMPLICON_MERGE.out.log.map { meta, log_file -> [meta + [tool: 'fastp_merge'], log_file] })
        ch_reads_pairhandled = FASTP_AMPLICON_MERGE.out.reads_merged
    }
    else {
        ch_reads_pairhandled = channel.empty()
    }

    ch_reads_preprocessed = ch_reads_pairhandled
        .map { meta, merged_reads -> [meta + [single_end: true], [merged_reads]] }
        .mix(ch_reads_for_pairhandling.other)

    /*
        ANALYSIS-READY READS
    */

    // Which step writes the analysis-ready reads of an amplicon run is known only once the primers are searched (the
    // runs without primers are not trimmed), so they are copied here, at the end of the amplicon preprocessing, rather
    // than published by each step. The runs of the samples with more than one run are published after run merging
    if (params.save_analysis_ready_fastqs) {
        ch_reads_preprocessed
            .filter { meta, _input_reads -> !(params.perform_runmerging && meta.is_multirun) }
            .subscribe { meta, input_reads ->
                def fastqs = [input_reads].flatten()
                fastqs.eachWithIndex { fastq, i ->
                    def suffix = fastqs.size() > 1 ? "_${i + 1}" : ''
                    def target = file("${params.outdir}/analysis_ready_fastqs/${meta.id}_${meta.run_accession}${suffix}.fastq.gz")
                    target.parent.mkdirs()
                    fastq.copyTo(target)
                }
            }
    }

    emit:
    reads           = ch_reads_preprocessed // channel: [ val(meta), [ reads ] ], single-end unless --amplicon_pair_handling none
    read_count_logs = ch_read_count_logs    // channel: [ val(meta), log ], with the tool in meta.tool
}
