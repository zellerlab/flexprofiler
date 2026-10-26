//
// Run profiling
//

include { BIOAWK as MOTUS4_RENAME_READS } from '../../modules/nf-core/bioawk/main'
include { MOTUS4_PROFILE                } from '../../modules/local/motus4/profile/main'
include { MOTUS_PROFILE  as MOTUS3_PROFILE  } from '../../modules/nf-core/motus/profile'
include { CAYMAN_PROFILE                } from '../../modules/local/cayman/profile/main'
include { BWA_MEM as BWA_ALIGN         } from '../../modules/nf-core/bwa/mem/main'
include { GQLIB_COUNT                   } from '../../modules/local/gqlib/count/main'
include { GQLIB_ROLLUP                  } from '../../modules/local/gqlib/rollup/main'

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

    // cayman and gqlib both profile a bwa alignment to the database: the databases with the same
    // path are grouped, so that the reads are aligned only once to each of them
    ch_dbs_bwa = ch_dbs
        .filter { _type, db_meta, _db -> db_meta.tool in ['cayman', 'gqlib'] && params["run_${db_meta.tool}"] }
        .groupTuple(by: [0, 2])
        .map { type, db_metas, db ->
            def db_names = db_metas.collect { db_meta -> db_meta.db_name }.unique().sort()
            [type, [tool: 'bwa', db_name: db_names.join('-'), type: type.type, profilers: db_metas], db]
        }

    ch_input_for_profiling = combineReadsWithDatabases(reads, ch_dbs.filter { _type, db_meta, _db -> !(db_meta.tool in ['motus4', 'cayman', 'gqlib']) }.mix(ch_dbs_bwa))
        .mix(combineReadsWithDatabases(ch_reads_for_motus4, ch_dbs.filter { _type, db_meta, _db -> db_meta.tool == 'motus4' }))
        .filter{
            // TODO: this is an initial way to simplify development and focus on adding profilers instead than supporting any data modality
            // Support for long reads and fasta will be enabled on a case-by-case basis depending on need
            // I can filter pairs of tools/samples accordingly at this step
            read_meta, input_reads, db_meta, db ->
            if (read_meta.is_fasta || read_meta.instrument_platform != "ILLUMINA") {
                log.warn("[zellerlab/flexprofiler] currently accepts only ILLUMINA FASTQ as input. Skipping sample ${read_meta.id}.")
            }
            !read_meta.is_fasta && read_meta.instrument_platform == "ILLUMINA"
        }
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

        // spread out each alignment to the profilers using its database
        ch_bams_for_profiling = BWA_ALIGN.out.bam
            .join(ch_input_for_align.db)
            .flatMap { meta, bam, db ->
                meta.profilers.collect { db_meta -> [meta.subMap(meta.keySet() - 'profilers') + db_meta, bam, db] }
            }
            .branch { meta, _bam, _db ->
                cayman: meta.tool == 'cayman'
                gqlib: meta.tool == 'gqlib'
            }
    }

    if (params.run_cayman) {
        ch_input_for_cayman = ch_bams_for_profiling.cayman
            .multiMap { meta, bam, db ->
                bam: [meta, bam]
                db: db
            }

        CAYMAN_PROFILE(ch_input_for_cayman.bam, ch_input_for_cayman.db)
        ch_versions = ch_versions.mix(CAYMAN_PROFILE.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(CAYMAN_PROFILE.out.out)
        ch_multiqc_files = ch_multiqc_files.mix(CAYMAN_PROFILE.out.log)
    }

    if (params.run_gqlib) {
        GQLIB_COUNT(ch_bams_for_profiling.gqlib.map { meta, bam, _db -> [meta, bam] })
        ch_versions = ch_versions.mix(GQLIB_COUNT.out.versions.first())

        // functional roll-up of the gene counts, only for databases with an annotation table
        ch_input_for_gqlib_rollup = GQLIB_COUNT.out.gene_counts
            .join(GQLIB_COUNT.out.read_counts)
            .join(ch_bams_for_profiling.gqlib.map { meta, _bam, db -> [meta, db] })
            .filter { _meta, _gene_counts, _read_counts, db ->
                db.listFiles().any { f -> f.name ==~ /.*\.annotation\.tsv(\.gz)?/ }
            }
            .multiMap { meta, gene_counts, read_counts, db ->
                counts: [meta, gene_counts, read_counts]
                db: db
            }

        GQLIB_ROLLUP(ch_input_for_gqlib_rollup.counts, ch_input_for_gqlib_rollup.db)
        ch_versions = ch_versions.mix(GQLIB_ROLLUP.out.versions.first())
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
