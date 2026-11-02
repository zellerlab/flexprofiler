//
// Concatenate the runs of each sample (--perform_runmerging), separately for each library strategy
//

include { CAT_FASTQ as MERGE_RUNS } from '../../modules/nf-core/cat/fastq/main'

workflow RUN_MERGING {
    take:
    reads // [ [ meta ], [ reads ] ], the preprocessed reads of each run

    main:
    ch_versions = channel.empty()

    if (params.perform_runmerging) {
        ch_reads_for_cat_branch = reads
            .map { meta, input_reads ->
                def meta_new = meta - meta.subMap('run_accession')
                [meta_new, meta.run_accession, input_reads]
            }
            .groupTuple()
            // groupTuple collects the runs in the order in which they finish, so they
            // are sorted to concatenate them always in the same order (the result of
            // mOTUs depends on the order of the reads). The runs are grouped by the meta without run_accession, which
            // contains library_strategy, so the WGS and AMPLICON runs of a sample are merged separately: the merged
            // runs are named merged_<strategy>, a single run keeps its run_accession
            .map { meta, run_accessions, run_reads ->
                def run_accession = run_accessions.size() > 1 ? "merged_${meta.library_strategy}" : run_accessions[0]
                [meta + [run_accession: run_accession], [run_accessions, run_reads].transpose().sort { run_reads_pair -> run_reads_pair[0].toString() }.collect { run_reads_pair -> run_reads_pair[1] }.flatten()]
            }
            .branch { meta, run_reads ->
                cat: (meta.single_end && run_reads.size() > 1) || (!meta.single_end && run_reads.size() > 2)
                skip: true
            }
        // we can't concatenate files if there is not a second run, so we branch
        // here to separate them out, and mix back in after for efficiency
        MERGE_RUNS(ch_reads_for_cat_branch.cat)
        ch_versions = ch_versions.mix(MERGE_RUNS.out.versions)

        ch_reads_runmerged = MERGE_RUNS.out.reads
            .mix(ch_reads_for_cat_branch.skip)
            .map { meta, merged_reads ->
                [meta, [merged_reads].flatten()]
            }
    }
    else {
        ch_reads_runmerged = reads
    }

    emit:
    reads    = ch_reads_runmerged // channel: [ val(meta), [ reads ] ], one per sample and library strategy with --perform_runmerging, else per run
    versions = ch_versions        // channel: [ versions.yml ]
}
