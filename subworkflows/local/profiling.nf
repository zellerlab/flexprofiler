//
// Run profiling
//

include { samplesheetToList             } from 'plugin/nf-schema'
include { BIOAWK as MOTUS4_RENAME_READS } from '../../modules/nf-core/bioawk/main'
include { MOTUS4_PROFILE                } from '../../modules/local/motus4/profile/main'
include { MOTUS_PROFILE  as MOTUS3_PROFILE  } from '../../modules/nf-core/motus/profile'
include { CAYMAN_PROFILE                } from '../../modules/local/cayman/profile/main'
include { BWA_MEM as BWA_ALIGN         } from '../../modules/nf-core/bwa/mem/main'
include { GQLIB_PROFILE_GENES           } from '../../modules/local/gqlib/profile_genes/main'
include { IDTAXA_CLASSIFY               } from '../../modules/local/idtaxa/classify/main'
include { KRAKEN2_KRAKEN2                } from '../../modules/nf-core/kraken2/kraken2/main'

workflow PROFILING {
    take:
    reads      // [ [ meta ], [ reads ] ]
    singletons // [ [ meta ], fastq ], the paired-end reads left without their mate by --pathseq_t2t, for Kraken2 and mOTUs
    databases  // [ [ meta ], path ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    // count tables
    ch_raw_profiles = Channel.empty()
    // logs of the profilers, for the read counts, with the tool in meta.tool: [ meta, log ]
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
    ['motus4', 'motus3', 'cayman', 'gqlib', 'idtaxa', 'kraken2'].findAll { tool -> params["run_${tool}"] }.each { tool ->
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

    // mOTUs takes the reads left without their mate by --pathseq_t2t in the same run as the pairs, as unpaired reads:
    // [ meta, [ fwd, rev, singletons ] ], or [ meta, [ singletons ] ] if no pair is left. The empty files are left out.
    // They are not renamed for mOTUs v4, as the unpaired reads are not matched by name
    ch_singletons_for_motus = singletons.filter { _meta, fastq -> fastq.countFastq() > 0 }
    ch_reads_for_motus4 = addSingletons(ch_reads_for_motus4, ch_singletons_for_motus)
    ch_reads_for_motus3 = addSingletons(reads, ch_singletons_for_motus)

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

    // the read counts of cayman and gqlib are those of the alignment, the same for all the annotations of a gene
    // catalogue: only the log of the annotation first in the database sheet is used for them. The sheet is read again,
    // as the order of its rows is lost in the database channel
    def first_annots = samplesheetToList(params.databases, "assets/schema_database.json")
        .collect { db_meta, _db_path -> db_meta }
        .findAll { db_meta -> db_meta.tool in ['cayman', 'gqlib'] }
        // the first row of each gene catalogue and profiler
        .unique { db_meta -> [db_meta.tool, db_meta.db_name] }
        .collect { db_meta -> [db_meta.tool, db_meta.db_name, db_meta.annot_db_name] }
    def is_first_annot_log = { meta, _log_file -> [meta.tool, meta.db_name, meta.annot_db_name] in first_annots }

    // the amplicon reads are single sequences (merged or joined read pairs, mates of unmerged pairs, or single-end
    // reads), as made by the amplicon preprocessing
    ch_input_for_profiling = combineReadsWithDatabases(reads, ch_dbs.filter { _type, db_meta, _db -> !(db_meta.tool in ['motus4', 'motus3', 'cayman', 'gqlib']) }.mix(ch_dbs_bwa))
        .mix(combineReadsWithDatabases(ch_reads_for_motus4, ch_dbs.filter { _type, db_meta, _db -> db_meta.tool == 'motus4' }))
        .mix(combineReadsWithDatabases(ch_reads_for_motus3, ch_dbs.filter { _type, db_meta, _db -> db_meta.tool == 'motus3' }))
        .branch { _meta, _input_reads, db_meta, _db ->
            motus4: db_meta.tool == 'motus4'
            motus3: db_meta.tool == 'motus3'
            bwa: db_meta.tool == 'bwa'
            idtaxa: db_meta.tool == 'idtaxa'
            kraken2: db_meta.tool == 'kraken2'
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
                meta.profilers.collect { db_meta -> [meta.subMap(meta.keySet() - 'profilers') + db_meta, bam, file(db_meta.annot_db_path, checkIfExists: true)] }
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
        ch_read_count_logs = ch_read_count_logs.mix(CAYMAN_PROFILE.out.log.filter(is_first_annot_log))
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
        ch_read_count_logs = ch_read_count_logs.mix(GQLIB_PROFILE_GENES.out.log.filter(is_first_annot_log))
    }

    if (params.run_idtaxa) {
        ch_input_for_idtaxa = ch_input_for_profiling.idtaxa
            .multiMap { read_meta, input_reads, db_meta, db ->
                reads: [read_meta + db_meta, input_reads]
                db: db
            }

        IDTAXA_CLASSIFY(ch_input_for_idtaxa.reads, ch_input_for_idtaxa.db)
        ch_versions = ch_versions.mix(IDTAXA_CLASSIFY.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(IDTAXA_CLASSIFY.out.profile)
        ch_read_count_logs = ch_read_count_logs.mix(IDTAXA_CLASSIFY.out.log)
    }

    if (params.run_kraken2) {
        // the singletons are classified on their own, as single-end reads, as in PathSeq-T2T: their counts are added to
        // those of the pairs by the merge (KRAKEN2_MERGE), counting each read once
        ch_singletons_for_kraken2 = combineReadsWithDatabases(
            singletons.map { meta, fastq -> [meta + [single_end: true, singletons: true], [fastq]] },
            ch_dbs.filter { _type, db_meta, _db -> db_meta.tool == 'kraken2' }
        )
        ch_input_for_kraken2 = ch_input_for_profiling.kraken2
            .mix(ch_singletons_for_kraken2)
            .multiMap { read_meta, input_reads, db_meta, db ->
                reads: [read_meta + db_meta, input_reads]
                db: db
            }

        KRAKEN2_KRAKEN2(ch_input_for_kraken2.reads, ch_input_for_kraken2.db, false, false)
        ch_raw_profiles = ch_raw_profiles.mix(KRAKEN2_KRAKEN2.out.report)
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report)
        // the reads of a pair are one fragment for Kraken2: counted twice in the read counts
        ch_read_count_logs = ch_read_count_logs.mix(
            KRAKEN2_KRAKEN2.out.report.map { meta, report -> [meta + [tool: meta.single_end ? 'kraken2' : 'kraken2_pairs'], report] }
        )
    }

    emit:
    profiles        = ch_raw_profiles // channel: [ val(meta), [ reads ] ] - should be text files or biom
    read_count_logs = ch_read_count_logs // channel: [ val(meta), log ], the logs of motus4, motus3, cayman, gqlib (of the first annotation of each gene catalogue) and idtaxa (with meta.db_name), with the tool in meta.tool
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
* Add to the reads of each sample its reads left without their mate, if any.
*
* @param ch_reads A channel of [ meta, [ reads ] ]
* @param ch_singletons A channel of [ meta, fastq ], with the same meta as the reads (or with no reads, if no pair is left)
* @return A channel of [ meta, [ reads, singletons ] ], [ meta, [ reads ] ] without singletons, [ meta, [ singletons ] ]
*     without reads
*/
def addSingletons(ch_reads, ch_singletons) {
    return ch_reads
        .join(ch_singletons, remainder: true)
        .map { meta, input_reads, singletons_fastq -> [meta, (input_reads ?: []) + (singletons_fastq ? [singletons_fastq] : [])] }
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
