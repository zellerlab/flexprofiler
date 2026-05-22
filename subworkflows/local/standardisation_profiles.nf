//
// Standardise output files e.g. aggregation
//

include { MOTUS4_MERGE } from '../../modules/local/motus4/merge/main'
include { CAYMAN_MERGE } from '../../modules/local/cayman/merge/main'

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
        cayman: it[0]['tool'] == 'cayman'
        unknown: true
//        bracken: it[0]['tool'] == 'bracken'
//        centrifuge: it[0]['tool'] == 'centrifuge'
//        ganon: it[0]['tool'] == 'ganon'
//        kmcp: it[0]['tool'] == 'kmcp'
//        kraken2: it[0]['tool'] == 'kraken2' || it[0]['tool'] == 'kraken2-bracken'
//        metaphlan: it[0]['tool'] == 'metaphlan'
    }

    ch_input_classifications = classifications.branch {
//        kaiju: it[0]['tool'] == 'kaiju'
        unknown: true
    }

    ch_input_databases = databases.branch {
        motus4: it[0]['tool'] == 'motus4'
        cayman: it[0]['tool'] == 'cayman'
//        kaiju: it[0]['tool'] == 'kaiju'
        unknown: true
    }

    /*
        Standardise and aggregate
    */


    // mOTUs

    ch_profiles_for_motus4 = groupProfiles(ch_input_profiles.motus4)
    ch_input_for_motus4merge = combineProfilesWithDatabase(ch_profiles_for_motus4, ch_input_databases.motus4)
    MOTUS4_MERGE(ch_input_for_motus4merge.profile, ch_input_for_motus4merge.db)
    ch_versions = ch_versions.mix(MOTUS4_MERGE.out.versions)
    
    // cayman
    // TODO: this is not yet tested
    ch_profiles_for_cayman = groupProfiles(ch_input_profiles.cayman)
    ch_input_for_caymanmerge = combineProfilesWithDatabase(ch_profiles_for_cayman, ch_input_databases.cayman)
    CAYMAN_MERGE(ch_input_for_caymanmerge.profile, ch_input_for_caymanmerge.db)
    ch_versions = ch_versions.mix(CAYMAN_MERGE.out.versions)

    emit:
    //taxpasta = TAXPASTA_MERGE.out.merged_profiles
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
        .map { db_name, profiles -> [[id: db_name], profiles] }
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
