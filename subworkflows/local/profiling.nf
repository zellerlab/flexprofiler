//
// Run profiling
//

include { BIOAWK as MOTUS4_RENAME_READS } from '../../modules/nf-core/bioawk/main'
include { MOTUS4_PROFILE                } from '../../modules/local/motus4/profile/main'
include { CAYMAN_PROFILE               } from '../../modules/local/cayman/profile/main'

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

    ch_input_for_profiling = reads
        .map { meta, input_reads ->
            [[type: meta.type], meta, input_reads]
        }
        .combine(ch_dbs, by: 0)
        .map { _db_type, meta, input_reads, db_meta, db ->
            [meta, input_reads, db_meta, db]
        }
        .branch { _meta, _input_reads, db_meta, _db ->
            motus4: db_meta.tool == 'motus4'
            cayman: db_meta.tool == 'cayman'
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
        if (params.motus4_rename_reads) {
            // mOTUs needs the read names to be identical across PE reads, without /1 /2 suffix
            // this step takes care of this
            ch_input_for_motus4_raw = ch_input_for_profiling.motus4
                .branch { read_meta, reads, db_meta, db -> 
                    se: read_meta.single_end
                    pe: !read_meta.single_end
                }
            
            ch_input_for_motus4_rename_reads = ch_input_for_motus4_raw.pe
                .map { read_meta, reads, db_meta, db -> [read_meta, reads] }
                .flatMap { meta, reads ->
                    [
                        [ meta + [paired_idx: 1], reads[0] ],
                        [ meta + [paired_idx: 2], reads[1] ]
                    ]
                }

            MOTUS4_RENAME_READS(ch_input_for_motus4_rename_reads)
            ch_versions = ch_versions.mix(MOTUS4_RENAME_READS.out.versions.first())

            ch_input_for_motus4_profile = MOTUS4_RENAME_READS.out.output
                .map {
                    meta, read -> [ meta.findAll { k, v -> k != "paired_idx" }, [meta.paired_idx, read] ]
                }
                .groupTuple(by: 0)
                .map {
                    meta, read_list ->
                    def r1 = read_list.find { idx, read -> idx == 1 }[1]
                    def r2 = read_list.find { idx, read -> idx == 2 }[1]
                    [ meta, [r1, r2] ]
                }
                .join(ch_input_for_motus4_raw.pe.map { read_meta, reads, db_meta, db -> [read_meta, db_meta, db] }, failOnDuplicate: true, failOnMismatch: true)
                .mix(ch_input_for_motus4_raw.se)
                .multiMap { read_meta, reads, db_meta, db ->
                    reads: [read_meta + db_meta, reads]
                    db: db
                }
        } else {
            ch_input_for_motus4_profile = ch_input_for_profiling.motus4
                .multiMap { read_meta, reads, db_meta, db ->
                    reads: [read_meta + db_meta, reads]
                    db: db
                }
        }

        MOTUS4_PROFILE(ch_input_for_motus4_profile.reads, ch_input_for_motus4_profile.db)
        ch_versions = ch_versions.mix(MOTUS4_PROFILE.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(MOTUS4_PROFILE.out.out)
        ch_multiqc_files = ch_multiqc_files.mix(MOTUS4_PROFILE.out.log)
    }

    if (params.run_cayman) {
        ch_input_for_cayman = ch_input_for_profiling.cayman
            .multiMap { read_meta, reads, db_meta, db ->
                reads: [read_meta + db_meta, reads]
                db: db
            }

        CAYMAN_PROFILE(ch_input_for_cayman.reads, ch_input_for_cayman.db)
        ch_versions = ch_versions.mix(CAYMAN_PROFILE.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(CAYMAN_PROFILE.out.out)
        ch_multiqc_files = ch_multiqc_files.mix(CAYMAN_PROFILE.out.log)
    }

    emit:
    classifications = ch_raw_classifications
    profiles        = ch_raw_profiles // channel: [ val(meta), [ reads ] ] - should be text files or biom
    versions        = ch_versions // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}


// Custom Functions

/**
* Combine profiles with their original database, then separate into two channels.
*
* The channel elements are assumed to be tuples one of [ meta, profile ], and the
* database to be of [db_key, meta, database_file].
*
* @param ch_profile A channel containing a meta and the profilign report of a given profiler
* @param ch_database A channel containing a key, the database meta, and the database file/folders itself
* @return A multiMap'ed output channel with two sub channels, one with the profile and the other with the db
*/
def combineProfilesWithDatabase(ch_profile, ch_database) {

    return ch_profile
        .map { meta, profile -> [meta.db_name, meta, profile] }
        .combine(ch_database, by: 0)
        .multiMap { _key, meta, profile, _db_meta, db ->
            profile: [meta, profile]
            db: db
        }
}
