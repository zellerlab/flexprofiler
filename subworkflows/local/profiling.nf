//
// Run profiling
//

include { BIOAWK as MOTUS4_RENAME_READS } from '../../modules/nf-core/bioawk/main'
include { MOTUS4_PROFILE                } from '../../modules/local/motus4/profile/main'
include { MOTUS_PROFILE  as MOTUS3_PROFILE  } from '../../modules/nf-core/motus/profile'
include { CAYMAN_PROFILE                } from '../../modules/local/cayman/profile/main'
include { BWA_MEM as BWA_ALIGN         } from '../../modules/nf-core/bwa/mem/main'
include { GQLIB_PROFILE_GENES           } from '../../modules/local/gqlib/profile_genes/main'

workflow PROFILING {
    take:
    reads     // [ [ meta ], [ reads ] ]
    databases // [ [ meta ], path ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_raw_classifications = Channel.empty()
    // These per-read ID taxonomic assingment
    ch_raw_profiles = Channel.empty()
    // These are count tables
    
    // TODO: add a step to produce a total read count here

    /*
        COMBINE READS WITH POSSIBLE DATABASES
    */

    // Separate default 'short;long' (when necessary) databases when short/long specified in database sheet
    ch_dbs = databases
        .map { meta_db, db ->
            [[meta_db.db_type.split(";")].flatten(), meta_db, db]
        }
        .transpose(by: 0)
        .map { type, meta_db, db ->
            [[type: type], meta_db.subMap(meta_db.keySet() - 'db_type') + [type: type], db]
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
                pe: !meta.single_end
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

    ch_input_for_profiling = combineReadsWithDatabases(reads, ch_dbs.filter { _type, db_meta, _db -> !(db_meta.tool in ['motus4', 'cayman', 'gqlib']) }.mix(ch_dbs_bwa))
        .mix(combineReadsWithDatabases(ch_reads_for_motus4, ch_dbs.filter { _type, db_meta, _db -> db_meta.tool == 'motus4' }))
        .branch { _meta, _input_reads, db_meta, _db ->
            motus4: db_meta.tool == 'motus4'
            motus3: db_meta.tool == 'motus3'
            bwa: db_meta.tool == 'bwa'
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
    }

    emit:
    classifications = ch_raw_classifications
    profiles        = ch_raw_profiles // channel: [ val(meta), [ reads ] ] - should be text files or biom
    versions        = ch_versions // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}


// Custom Functions

/**
* Pair reads with all databases of the matching read type (short/long).
*
* @param ch_reads A channel of [ meta, reads ], where meta contains the key `type`
* @param ch_databases A channel of [ [type: type], db_meta, db ]
* @return A channel of [ meta, reads, db_meta, db ]
*/
def combineReadsWithDatabases(ch_reads, ch_databases) {
    return ch_reads
        .map { meta, input_reads -> [[type: meta.type], meta, input_reads] }
        .combine(ch_databases, by: 0)
        .map { _db_type, meta, input_reads, db_meta, db -> [meta, input_reads, db_meta, db] }
}
