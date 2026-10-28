//
// Standardise output files e.g. aggregation
//

include { MOTUS4_MERGE                 } from '../../modules/local/motus4/merge/main'
include { MOTUS_MERGE  as MOTUS3_MERGE } from '../../modules/nf-core/motus/merge/main'
include { GQLIB_MERGE                  } from '../../modules/local/gqlib/merge/main'
// cayman profiles have the same format as the gqlib ones
include { GQLIB_MERGE  as CAYMAN_MERGE } from '../../modules/local/gqlib/merge/main'
include { IDTAXA_MERGE                 } from '../../modules/local/idtaxa/merge/main'

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
        gqlib: it[0]['tool'] == 'gqlib'
        idtaxa: it[0]['tool'] == 'idtaxa'
        unknown: true
    }

    ch_input_classifications = classifications.branch {
        unknown: true
    }

    ch_input_databases = databases.branch {
        motus4: it[0]['tool'] == 'motus4'
        motus3: it[0]['tool'] == 'motus3'
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
    ch_profiles_for_cayman = groupProfilesByCategory(ch_input_profiles.cayman)
    CAYMAN_MERGE(ch_profiles_for_cayman)
    ch_versions = ch_versions.mix(CAYMAN_MERGE.out.versions.first())

    // gqlib, merged separately for each category (gene_counts, KEGG_ko, ...) of each database and annotation
    ch_profiles_for_gqlib = groupProfilesByCategory(ch_input_profiles.gqlib.transpose())
    GQLIB_MERGE(ch_profiles_for_gqlib)
    ch_versions = ch_versions.mix(GQLIB_MERGE.out.versions.first())

    // IDTAXA, merged separately for each training set, one table per rank
    ch_profiles_for_idtaxa = groupProfiles(ch_input_profiles.idtaxa)
    IDTAXA_MERGE(ch_profiles_for_idtaxa)
    ch_versions = ch_versions.mix(IDTAXA_MERGE.out.versions.first())

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
            // groupTuple collects the profiles in the order in which they finish, so they are
            // sorted by name to always have the samples in the same order in the merged profile
            [[id: db_name], profiles.sort { it.getName() }]
        }
}

/**
* Group gqlib (and cayman) profiles per reference database, annotation and category.
*
* @param ch_profiles A channel containing pairs of a meta map, which must contain
*   the keys `db_name` and `annot_db_name`, and a single profile.
* @return A channel with one element per reference database, annotation and
*   category. Each element is a pair of a meta map with `id` (database and
*   annotation names, as in the profile names), `db_name`, `annot_db_name` and
*   `category` keys, and all corresponding profiles.
*/
def groupProfilesByCategory(ch_profiles) {
    return ch_profiles
        .map { meta, profile ->
            def key = [id: "${meta.db_name}_${meta.annot_db_name}", db_name: meta.db_name, annot_db_name: meta.annot_db_name, category: gqlibCategory(profile)]
            [key, profile]
        }
        .groupTuple()
        // sorted by name, as in groupProfiles
        .map { key, profiles -> [key, profiles.sort { it.getName() }] }
}

/**
* Get the category of a gqlib profile from its file name, <prefix>.<category>.txt.gz.
*
* The prefix can contain dots (the database name usually does), the category cannot.
*/
def gqlibCategory(profile) {
    return (profile.getName() - ~/\.txt\.gz$/).tokenize('.')[-1]
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
