//
// Standardise output files e.g. aggregation
//

include { MOTUS4_MERGE                 } from '../../modules/local/motus4/merge/main'
include { MOTUS_MERGE  as MOTUS3_MERGE } from '../../../modules/nf-core/motus/merge'
include { CAYMAN_MERGE                 } from '../../modules/local/cayman/merge/main'

workflow STANDARDISATION_PROFILES {
    take:
    classifications
    profiles
    databases

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*
        Split profile results based on tool they come from
    */
    ch_input_profiles = profiles.branch {
        motus4: it[0]['tool'] == 'motus4'
        motus3: it[0]['tool'] == 'motus3'
        cayman: it[0]['tool'] == 'cayman'
        unknown: true
    }

    ch_input_classifications = classifications.branch {
        unknown: true
    }

    ch_input_databases = databases.branch {
        motus4: it[0]['tool'] == 'motus4'
        motus3: it[0]['tool'] == 'motus3'
        cayman: it[0]['tool'] == 'cayman'
        unknown: true
    }

    /*
        Standardise and aggregate
    */


    // mOTUs v4
    ch_profiles_for_motus4 = groupProfiles(ch_input_profiles.motus4)
    ch_input_for_motus4merge = combineProfilesWithDatabase(ch_profiles_for_motus4, ch_input_databases.motus4)
    MOTUS4_MERGE(ch_input_for_motus4merge.profile, ch_input_for_motus4merge.db)
    ch_versions = ch_versions.mix(MOTUS4_MERGE.out.versions)
    
    // mOTUs v3
    ch_profiles_for_motus3 = groupProfiles(ch_input_profiles.motus3)
    ch_input_for_motus3merge = combineProfilesWithDatabase(ch_profiles_for_motus3, ch_input_databases.motus3)
    MOTUS3_MERGE(ch_input_for_motus3merge.profile, ch_input_for_motus3merge.db)
    ch_versions = ch_versions.mix(MOTUS3_MERGE.out.versions)
    
    // cayman
    ch_profiles_for_cayman = groupProfiles(ch_input_profiles.cayman)
    ch_input_for_caymanmerge = combineProfilesWithDatabase(ch_profiles_for_cayman, ch_input_databases.cayman)
    CAYMAN_MERGE(ch_input_for_caymanmerge.profile, ch_input_for_caymanmerge.db)
    ch_versions = ch_versions.mix(CAYMAN_MERGE.out.versions)

    emit:
    versions = ch_versions
    mqc      = ch_multiqc_files
}

// Custom Functions

/**
* Group all profiles per reference database.
*
* @param ch_profiles A channel containing pairs of a meta map and the report of
*   a given profiler, where meta must contain a key `db_name`.
* @return A channel with one element per reference database. Each element is a
*   pair of a meta map with an `id` key and all corresponding profiles.
*/
def groupProfiles(ch_profiles, groupTupleOptions = [:]) {
    return ch_profiles
        .map { meta, profile -> [meta.db_name, profile] }
        .groupTuple(groupTupleOptions)
        .map { 
            db_name, profiles -> 
            def names = profiles.collect { it.getName() }
            def duplicates = names.findAll { name -> names.count(name) > 1 }.unique()
            if (duplicates) {
                error (
                    "Duplicate sample names identified in profile standardisation: ${duplicates.join(', ')}\n" +
                    "This usually happens if sample names collide and run merging is not activated.\n" +
                    "If run merging is active, runs obtained with different sequencing technology or a mix of paired and single end files cannot be merged.\n" +
                    "This can result in multiple merged samples with the same name.\n" +
                    "Please rename your sample column in the samplesheet accordingly."
                )
            }
            [[id: db_name], profiles] 
        }
}

/**
* Combine profiles with their corresponding reference database, then separate into two channels.
*
* The combined results are returned on multiple channels, where the element
* position for the profiles in one channel is the same as the position of the
* corresponding database element in the other channel.
*
* @param ch_profiles A channel containing pairs of a meta map with an `id` key
*   for a reference database, and all the corresponding profiling reports.
* @param ch_database A channel containing pairs of a database meta map and the
*   database itself.
* @return A multiMap'ed output channel with two sub channels, one with the
*   profiles (`profile`) and the other with the corresponding database (`db`).
*/
def combineProfilesWithDatabase(ch_profiles, ch_database) {
    return ch_profiles
        .map { meta, profile -> [meta.id, meta, profile] }
        .combine(ch_database.map { db_meta, db -> [db_meta.db_name, db] }, by: 0)
        .multiMap { _key, meta, profile, db ->
            profile: [meta, profile]
            db: db
        }
}
