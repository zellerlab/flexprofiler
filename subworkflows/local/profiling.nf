//
// Run profiling
//

include { BIOAWK as MOTUS4_RENAME_READS } from '../../modules/nf-core/bioawk/main'
include { MOTUS4_PROFILE                } from '../../modules/local/motus4/profile/main'
include { MOTUS_PROFILE  as MOTUS3_PROFILE  } from '../../modules/nf-core/motus/profile'
include { CAYMAN_PROFILE                } from '../../modules/local/cayman/profile/main'
include { BWA_MEM as BWA_ALIGN         } from '../../modules/nf-core/bwa/mem/main'
include { GQLIB_PROFILE_GENES           } from '../../modules/local/gqlib/profile_genes/main'
include { VSEARCH_MERGEJOIN             } from '../../modules/local/vsearch/mergejoin/main'
include { FASTP as FASTP_IDTAXA_MERGE   } from '../../modules/nf-core/fastp/main'
include { IDTAXA_CLASSIFY               } from '../../modules/local/idtaxa/classify/main'

workflow PROFILING {
    take:
    reads     // [ [ meta ], [ reads ] ]
    databases // [ [ meta ], path ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    // per-read taxonomic assignments
    ch_raw_classifications = Channel.empty()
    // count tables
    ch_raw_profiles = Channel.empty()
    // logs of the profilers and of the steps preparing the reads for a profiler, for the read counts, with the tool in
    // meta.tool: [ meta, log ]
    ch_read_count_logs = Channel.empty()

    /*
        COMBINE READS WITH POSSIBLE DATABASES
    */

    // Separate default 'short;long' (when necessary) databases when short/long specified in database sheet. Each
    // database is also keyed by the library strategy of its profiler, so that WGS reads only meet the databases of
    // the shotgun profilers and AMPLICON reads those of the amplicon profilers
    ch_dbs = databases
        .map { meta_db, db ->
            [[meta_db.db_type.split(";")].flatten(), meta_db, db]
        }
        .transpose(by: 0)
        .map { type, meta_db, db ->
            [[type: type, library_strategy: toolLibraryStrategy(meta_db.tool)], meta_db.subMap(meta_db.keySet() - 'db_type') + [type: type], db]
        }

    // a profiler run on no reads is most likely a mistake in the samplesheet (library_strategy) or in the parameters
    ['motus4', 'motus3', 'cayman', 'gqlib', 'idtaxa'].findAll { tool -> params["run_${tool}"] }.each { tool ->
        def strategy = toolLibraryStrategy(tool)
        reads
            .filter { meta, _input_reads -> meta.library_strategy == strategy }
            .count()
            .subscribe { n ->
                if (n == 0) {
                    log.warn("[zellerlab/flexprofiler] --run_${tool} is set, but there are no ${strategy} samples (library_strategy column of the samplesheet) to run it on.")
                }
            }
    }

    // Join short and long reads with their corresponding short/long database
    // Note that for not-specified `short;long`, it will match with the database.
    // E.g. if there is no 'long' reads the above generated 'long' database channel element
    //  will have nothing to join to and will be discarded
    // Final output [DUMP: reads_plus_db] [['id':'2612', 'run_accession':'combined', 'instrument_platform':'ILLUMINA', 'single_end':false, 'is_fasta':false, 'type':'short'], <reads_path>/2612.merged.fastq.gz, ['tool':'malt', 'db_name':'malt95', 'db_params':'"-id 90"', 'type':'short'], <db_path>/malt95]

    // mOTUs v4 needs read names to be identical across PE mates, without /1 /2 suffix.
    // Reads are renamed once here, before being paired with the motus4 databases
    if (params.run_motus4 && params.motus4_rename_reads) {
        ch_reads_for_rename = reads
            .branch { meta, _input_reads ->
                pe: !meta.single_end && meta.library_strategy == 'WGS'
                other: true
            }

        ch_input_for_motus4_rename_reads = ch_reads_for_rename.pe
            .flatMap { meta, input_reads ->
                [
                    [ meta + [paired_idx: 1], input_reads[0] ],
                    [ meta + [paired_idx: 2], input_reads[1] ]
                ]
            }

        MOTUS4_RENAME_READS(ch_input_for_motus4_rename_reads)
        ch_versions = ch_versions.mix(MOTUS4_RENAME_READS.out.versions.first())

        ch_reads_for_motus4 = MOTUS4_RENAME_READS.out.output
            .map { meta, input_reads ->
                [ groupKey(meta.findAll { k, _v -> k != "paired_idx" }, 2), [meta.paired_idx, input_reads] ]
            }
            .groupTuple()
            .map { meta, read_list ->
                [ meta.getGroupTarget(), read_list.sort { idx_read -> idx_read[0] }.collect { idx_read -> idx_read[1] } ]
            }
            .mix(ch_reads_for_rename.other)
    } else {
        ch_reads_for_motus4 = reads
    }

    // cayman and gqlib databases are the gene catalogue fasta (db_path), with its bwa index alongside,
    // and the annotation file (annot_db_path). The alignment is shared by all the annotations of a
    // gene catalogue, and between cayman and gqlib: the databases with the same name and gene
    // catalogue fasta are grouped, so that the reads are aligned only once to each of them
    ch_dbs_bwa = ch_dbs
        .filter { _type, db_meta, _db -> db_meta.tool in ['cayman', 'gqlib'] && params["run_${db_meta.tool}"] }
        .map { type, db_meta, fasta -> [type, db_meta.db_name, fasta, db_meta] }
        .groupTuple(by: [0, 1, 2])
        .map { type, db_name, fasta, db_metas ->
            // the index files are named as created by `bwa index <fasta>`
            def index = ['amb', 'ann', 'bwt', 'pac', 'sa'].collect { ext -> file("${fasta}.${ext}", checkIfExists: true) }
            [type, [tool: 'bwa', db_name: db_name, type: type.type, profilers: db_metas], index]
        }

    // IDTAXA classifies single sequences: the pairs of paired-end amplicon reads are merged if they overlap, and the
    // others joined with a padding of Ns (merged_plus_joined, vsearch) or kept as two separate mates
    // (merged_plus_unmerged, fastp)
    if (params.run_idtaxa) {
        ch_reads_for_pairhandling = reads
            .filter { meta, _input_reads -> meta.library_strategy == 'AMPLICON' }
            .branch { meta, _input_reads ->
                pe: !meta.single_end
                se: true
            }
        if (params.idtaxa_pair_handling == 'merged_plus_joined') {
            VSEARCH_MERGEJOIN(ch_reads_for_pairhandling.pe)
            ch_versions = ch_versions.mix(VSEARCH_MERGEJOIN.out.versions.first())
            ch_read_count_logs = ch_read_count_logs.mix(VSEARCH_MERGEJOIN.out.log.map { meta, log_file -> [meta + [tool: 'vsearch'], log_file] })
            ch_reads_pairhandled = VSEARCH_MERGEJOIN.out.reads
        }
        else if (params.idtaxa_pair_handling == 'merged_plus_unmerged') {
            FASTP_IDTAXA_MERGE(ch_reads_for_pairhandling.pe, [], false, false, true)
            ch_versions = ch_versions.mix(FASTP_IDTAXA_MERGE.out.versions.first())
            ch_read_count_logs = ch_read_count_logs.mix(FASTP_IDTAXA_MERGE.out.log.map { meta, log_file -> [meta + [tool: 'fastp_merge'], log_file] })
            ch_reads_pairhandled = FASTP_IDTAXA_MERGE.out.reads_merged
        }
        ch_reads_for_idtaxa = ch_reads_pairhandled
            .map { meta, merged_reads -> [meta + [single_end: true], [merged_reads]] }
            .mix(ch_reads_for_pairhandling.se)
    }
    else {
        ch_reads_for_idtaxa = Channel.empty()
    }

    ch_input_for_profiling = combineReadsWithDatabases(reads, ch_dbs.filter { _type, db_meta, _db -> !(db_meta.tool in ['motus4', 'cayman', 'gqlib', 'idtaxa']) }.mix(ch_dbs_bwa))
        .mix(combineReadsWithDatabases(ch_reads_for_motus4, ch_dbs.filter { _type, db_meta, _db -> db_meta.tool == 'motus4' }))
        .mix(combineReadsWithDatabases(ch_reads_for_idtaxa, ch_dbs.filter { _type, db_meta, _db -> db_meta.tool == 'idtaxa' }))
        .branch { _meta, _input_reads, db_meta, _db ->
            motus4: db_meta.tool == 'motus4'
            motus3: db_meta.tool == 'motus3'
            bwa: db_meta.tool == 'bwa'
            idtaxa: db_meta.tool == 'idtaxa'
            unknown: true
        }

    /*
        PREPARE PROFILER INPUT CHANNELS & RUN PROFILING
    */

    // Each tool as a slightly different input structure and generally separate
    // input channels for reads vs databases. We restructure the channel tuple
    // for each tool and make liberal use of multiMap to keep reads/databases
    // channel element order in sync with each other

    if (params.run_motus4) {
        ch_input_for_motus4_profile = ch_input_for_profiling.motus4
            .multiMap { read_meta, input_reads, db_meta, db ->
                reads: [read_meta + db_meta, input_reads]
                db: db
            }

        MOTUS4_PROFILE(ch_input_for_motus4_profile.reads, ch_input_for_motus4_profile.db)
        ch_versions = ch_versions.mix(MOTUS4_PROFILE.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(MOTUS4_PROFILE.out.out)
        ch_multiqc_files = ch_multiqc_files.mix(MOTUS4_PROFILE.out.log)
        ch_read_count_logs = ch_read_count_logs.mix(MOTUS4_PROFILE.out.log)
    }

    if (params.run_motus3) {
        ch_input_for_motus3 = ch_input_for_profiling.motus3
            .multiMap { it ->
                reads: [it[0] + it[2], it[1]]
                db: it[3]
            }

        MOTUS3_PROFILE(ch_input_for_motus3.reads, ch_input_for_motus3.db)

        ch_versions = ch_versions.mix(MOTUS3_PROFILE.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(MOTUS3_PROFILE.out.out)
        ch_multiqc_files = ch_multiqc_files.mix(MOTUS3_PROFILE.out.log)
        ch_read_count_logs = ch_read_count_logs.mix(MOTUS3_PROFILE.out.log)
    }

    if (params.run_cayman || params.run_gqlib) {
        ch_input_for_align = ch_input_for_profiling.bwa
            .multiMap { read_meta, input_reads, db_meta, db ->
                reads: [read_meta + db_meta, input_reads]
                db: [read_meta + db_meta, db]
            }

        // the bam is not sorted, as cayman and gqlib expect the alignments grouped by read
        BWA_ALIGN(ch_input_for_align.reads, ch_input_for_align.db, [[], []], false)

        // spread out each alignment to the profilers using its gene catalogue, each with its annotation file
        ch_bams_for_profiling = BWA_ALIGN.out.bam
            .flatMap { meta, bam ->
                meta.profilers.collect { db_meta -> [meta.subMap(meta.keySet() - 'profilers') + db_meta, bam, file(db_meta.annot_db_path)] }
            }
            .branch { meta, _bam, _annot ->
                cayman: meta.tool == 'cayman'
                gqlib: meta.tool == 'gqlib'
            }
    }

    if (params.run_cayman) {
        ch_input_for_cayman = ch_bams_for_profiling.cayman
            .multiMap { meta, bam, annot ->
                bam: [meta, bam]
                annot: annot
            }

        CAYMAN_PROFILE(ch_input_for_cayman.bam, ch_input_for_cayman.annot)
        ch_versions = ch_versions.mix(CAYMAN_PROFILE.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(CAYMAN_PROFILE.out.out)
        ch_read_count_logs = ch_read_count_logs.mix(CAYMAN_PROFILE.out.log)
    }

    if (params.run_gqlib) {
        ch_input_for_gqlib = ch_bams_for_profiling.gqlib
            .multiMap { meta, bam, annot ->
                bam: [meta, bam]
                annot: annot
            }

        GQLIB_PROFILE_GENES(ch_input_for_gqlib.bam, ch_input_for_gqlib.annot)
        ch_versions = ch_versions.mix(GQLIB_PROFILE_GENES.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(GQLIB_PROFILE_GENES.out.profiles)
        ch_read_count_logs = ch_read_count_logs.mix(GQLIB_PROFILE_GENES.out.log)
    }

    if (params.run_idtaxa) {
        ch_input_for_idtaxa = ch_input_for_profiling.idtaxa
            .multiMap { read_meta, input_reads, db_meta, db ->
                reads: [read_meta + db_meta, input_reads]
                db: db
            }

        IDTAXA_CLASSIFY(ch_input_for_idtaxa.reads, ch_input_for_idtaxa.db)
        ch_versions = ch_versions.mix(IDTAXA_CLASSIFY.out.versions.first())
        ch_raw_classifications = ch_raw_classifications.mix(IDTAXA_CLASSIFY.out.classifications)
        ch_raw_profiles = ch_raw_profiles.mix(IDTAXA_CLASSIFY.out.profile)
        ch_read_count_logs = ch_read_count_logs.mix(IDTAXA_CLASSIFY.out.log)
    }

    emit:
    classifications = ch_raw_classifications
    profiles        = ch_raw_profiles // channel: [ val(meta), [ reads ] ] - should be text files or biom
    read_count_logs = ch_read_count_logs // channel: [ val(meta), log ], the logs of motus4, motus3, cayman, gqlib and idtaxa (with meta.db_name), and of vsearch or fastp (merging the pairs for idtaxa), with the tool in meta.tool
    versions        = ch_versions // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}


// Custom Functions

/**
* Pair reads with all databases of the matching read type (short/long) and library strategy (WGS/AMPLICON).
*
* @param ch_reads A channel of [ meta, reads ], where meta contains the keys `type` and `library_strategy`
* @param ch_databases A channel of [ [type: type, library_strategy: library_strategy], db_meta, db ]
* @return A channel of [ meta, reads, db_meta, db ]
*/
def combineReadsWithDatabases(ch_reads, ch_databases) {
    return ch_reads
        .map { meta, input_reads -> [[type: meta.type, library_strategy: meta.library_strategy], meta, input_reads] }
        .combine(ch_databases, by: 0)
        .map { _db_type, meta, input_reads, db_meta, db -> [meta, input_reads, db_meta, db] }
}

/**
* The library strategy of the reads a profiler is run on.
*
* @param tool The profiler, as in the tool column of the database sheet.
* @return AMPLICON for the amplicon profilers, WGS for the shotgun ones.
*/
def toolLibraryStrategy(tool) {
    return tool in ['idtaxa'] ? 'AMPLICON' : 'WGS'
}
