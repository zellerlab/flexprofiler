//
// Find the primers of amplicon reads with PIMENTO, or take those given (--amplicon_primer_fwd, --amplicon_primer_rev),
// and trim them with cutadapt
//

include { VSEARCH_ORIENT                      } from '../../modules/local/vsearch/orient/main'
include { PIMENTO_STD                         } from '../../modules/local/pimento/std/main'
include { CUTADAPT as CUTADAPT_PRIMERS        } from '../../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_READTHROUGH    } from '../../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_LINKED         } from '../../modules/nf-core/cutadapt/main'

workflow AMPLICON_PRIMERTRIMMING {
    take:
    reads            // [ [ meta ], [ reads ] ], AMPLICON runs, short (meta.type) paired-end or single-end, or long reads
    orient_reference // the reference to orient the long reads with (FASTA), [] if not given

    main:
    ch_versions = channel.empty()
    ch_read_count_logs = channel.empty()
    ch_primers_tsv = channel.empty()

    if (params.amplicon_primer_fwd && params.amplicon_primer_rev) {
        // the primers given, for all the runs: the long reads are oriented by cutadapt (--revcomp)
        ch_reads_with_primers = reads.map { meta, input_reads ->
            [meta + primerMeta([params.amplicon_primer_fwd.toUpperCase()], [params.amplicon_primer_rev.toUpperCase()]), input_reads]
        }
    }
    else {
        ch_reads_by_type = reads.branch { meta, _input_reads ->
            longreads: meta.type == 'long'
            shortreads: true
        }

        // Long reads (PacBio) are sequenced from either end of the amplicon: they are oriented as the reference
        // sequences (of the same gene) first, so that, as single-end short reads, they start with the forward primer
        // and end with the reverse complement of the reverse primer. The reads that cannot be oriented (not of the gene
        // of the reference) are dropped
        VSEARCH_ORIENT(ch_reads_by_type.longreads, orient_reference)
        ch_versions = ch_versions.mix(VSEARCH_ORIENT.out.versions.first())
        ch_read_count_logs = ch_read_count_logs.mix(VSEARCH_ORIENT.out.log.map { meta, log_file -> [meta + [tool: 'vsearch_orient'], log_file] })
        ch_reads_for_detection = ch_reads_by_type.shortreads
            .mix(VSEARCH_ORIENT.out.reads.map { meta, oriented_reads -> [meta, [oriented_reads]] })

        // PIMENTO searches each mate of paired-end reads separately (meta.mate, single-end reads: 1): the primers found
        // for the two mates (normally the forward primer for R1 and the reverse one for R2) are the primers of the run
        ch_mates = ch_reads_for_detection.flatMap { meta, input_reads ->
            [input_reads].flatten().indexed(1).collect { mate, fastq -> [meta + [mate: mate], fastq] }
        }
        PIMENTO_STD(ch_mates)
        ch_versions = ch_versions.mix(PIMENTO_STD.out.versions.first())

        // the primer files of the mates of a run, grouped again, in the order of the mates (<prefix>_<mate>_std_primers.fasta),
        // so that the primers, and the cutadapt arguments, are always in the same order
        ch_pimento_primers = PIMENTO_STD.out.primers
            .map { meta, fasta -> [groupKey(meta - meta.subMap('mate'), meta.single_end ? 1 : 2), fasta] }
            .groupTuple(sort: { fasta -> fasta.name })
            .map { meta, fastas -> [meta.getGroupTarget(), fastas] }

        // the forward and reverse primers, by the last letter of their names in the PIMENTO library (F or R)
        ch_reads_with_primers = ch_reads_for_detection
            .join(ch_pimento_primers)
            .map { meta, input_reads, primer_fastas ->
                def primers = primer_fastas.collectMany { fasta -> readPimentoPrimers(fasta, meta.single_end) }.unique()
                def primers_fwd = primers.findAll { name, _seq -> !name.endsWith('R') }.collect { _name, seq -> seq }
                def primers_rev = primers.findAll { name, _seq -> name.endsWith('R') }.collect { _name, seq -> seq }
                [meta + primerMeta(primers_fwd, primers_rev), input_reads]
            }

        // the primers found, one row per run and mate (single-end: mate 1), with the primer (if any), its proportion
        // of reads and the region of the library it belongs to
        ch_primers_tsv = PIMENTO_STD.out.out
            .map { meta, out_file ->
                def lines = out_file.text.readLines().findAll { line -> line.trim() }
                def region = lines ? lines[0] : ''
                def found = lines.drop(1).collect { line -> line.split(': ')*.trim() }
                if (!found) {
                    found = [['', '']]
                }
                found.collect { name_prop -> [meta.id, meta.run_accession, meta.mate, region, name_prop[0], name_prop[1]].join('\t') }.join('\n')
            }
            .collectFile(
                name: 'primers.tsv',
                seed: ['sample', 'run_accession', 'mate', 'region', 'primer', 'proportion'].join('\t'),
                newLine: true,
                sort: true,
            )
            // copied to the output directory here and not with storeDir, which creates the directory also without
            // amplicon samples
            .map { primers_tsv -> primers_tsv.copyTo("${params.outdir}/amplicon_primers/primers.tsv") }
    }

    // The long reads are trimmed only with a pair of primers, both required (linked), as in the best practice for PacBio
    // full-length 16S: the reads without both primers, truncated or chimeric, are discarded
    ch_reads_for_trimming = ch_reads_with_primers.branch { meta, _input_reads ->
        linked: meta.type == 'long' && meta.primer_linked
        trim: meta.type != 'long' && meta.primers
        skip: true
    }
    ch_reads_for_trimming.skip
        .filter { meta, _input_reads -> meta.type == 'long' && meta.primers }
        .subscribe { meta, _input_reads ->
            log.warn("[zellerlab/flexprofiler] the primers found for the long reads of ${meta.id}, run ${meta.run_accession}, are not a pair of a forward and a reverse primer: the reads are not trimmed.")
        }

    // All the primers of the run are searched on both mates, so that the pairs sequenced in the other orientation
    // (R1 starting with the reverse primer) are trimmed too: first at the start of the reads, then, at their end, the
    // reverse complement of the primers, for the reads running past the end of the amplicon
    CUTADAPT_PRIMERS(ch_reads_for_trimming.trim)
    CUTADAPT_READTHROUGH(CUTADAPT_PRIMERS.out.reads)

    // the forward primer at the start and the reverse complement of the reverse primer at the end, both required, in
    // either orientation (--revcomp), with the length window of the amplicon (ext.args)
    CUTADAPT_LINKED(ch_reads_for_trimming.linked)

    // the primers are removed from the meta, so that the runs of a sample are merged together whatever their primers
    def primer_keys = ['primers', 'primers_revcomp', 'primer_linked']
    ch_trimmed_reads = CUTADAPT_READTHROUGH.out.reads
        .mix(CUTADAPT_LINKED.out.reads.map { meta, trimmed_reads -> [meta, [trimmed_reads].flatten()] })
        .mix(ch_reads_for_trimming.skip)
        .map { meta, trimmed_reads -> [meta - meta.subMap(primer_keys), trimmed_reads] }

    ch_read_count_logs = ch_read_count_logs
        .mix(CUTADAPT_PRIMERS.out.log.map { meta, log_file -> [meta - meta.subMap(primer_keys) + [tool: 'cutadapt_primers'], log_file] })
        .mix(CUTADAPT_READTHROUGH.out.log.map { meta, log_file -> [meta - meta.subMap(primer_keys) + [tool: 'cutadapt_readthrough'], log_file] })
        .mix(CUTADAPT_LINKED.out.log.map { meta, log_file -> [meta - meta.subMap(primer_keys) + [tool: 'cutadapt_linked'], log_file] })

    emit:
    reads           = ch_trimmed_reads   // channel: [ val(meta), [ reads ] ]
    read_count_logs = ch_read_count_logs // channel: [ val(meta), log ], with the tool in meta.tool
    primers         = ch_primers_tsv     // channel: path(primers.tsv), only with the primers found by PIMENTO
    versions        = ch_versions        // channel: [ versions.yml ]
}

// Custom Functions

/**
* The meta fields of the primers of a run, for cutadapt.
*
* @param primers_fwd The forward primers, as in the primer library.
* @param primers_rev The reverse primers, as in the primer library (5'-3' on the reverse strand).
* @return primers (all of them) and primers_revcomp (their reverse complements), for the short reads, and primer_linked,
* the linked adapter of a pair of primers (one forward and one reverse), for the long reads, null otherwise. Both
* primers of the linked adapter are required in full (o: minimum overlap), so that a partial match of a few bases at
* the end of a read does not count as the primer.
*/
def primerMeta(primers_fwd, primers_rev) {
    def primers = primers_fwd + primers_rev
    def primer_linked = null
    if (primers_fwd.size() == 1 && primers_rev.size() == 1) {
        def fwd = primers_fwd[0]
        def rev = primers_rev[0]
        primer_linked = "${fwd};required;o=${fwd.length()}...${reverseComplement(rev)};required;o=${rev.length()}"
    }
    return [primers: primers, primers_revcomp: primers.collect { seq -> reverseComplement(seq) }, primer_linked: primer_linked]
}

/**
* Read the primers found by PIMENTO (pimento std), as they are in the primer library.
*
* PIMENTO writes the reverse primer of a pair complemented, but in --merged mode (single-end reads): the complement is
* undone. In --merged mode PIMENTO searches the reverse primers in the reversed end of the reads, so it matches them
* against their complement (the reverse of their reverse complement), stored complemented in its primer dictionary
* (parse_std_primers); write_std_output complements the reverse primer of a pair, and in --merged mode also a single
* reverse primer, so that in --merged mode both are written as in the library, and the reverse primer of a pair is
* written complemented only without --merged (pimento/bin/standard_primer_matching.py, PIMENTO 1.3.3: to check again
* when PIMENTO is updated).
*
* @param fasta The <prefix>_std_primers.fasta file of pimento std, empty if no primer was found.
* @param merged Whether pimento std was run with --merged (single-end reads).
* @return A list of [ name, sequence ].
*/
def readPimentoPrimers(fasta, merged) {
    def records = fasta.splitFasta(record: [id: true, seqString: true])
    return records.withIndex().collect { record, i ->
        def complemented = !merged && records.size() == 2 && i == 1
        [record.id, complemented ? complement(record.seqString) : record.seqString]
    }
}

/**
* Complement of a nucleotide sequence, with the IUPAC ambiguity codes (Cornish-Bowden, Nucleic Acids Res. 1985,
* doi:10.1093/nar/13.9.3021): an unknown base becomes X, which never matches in cutadapt.
*/
def complement(seq) {
    def bases = [A: 'T', T: 'A', U: 'A', G: 'C', C: 'G', Y: 'R', R: 'Y', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', D: 'H', H: 'D', V: 'B', N: 'N']
    return seq.toUpperCase().collect { base -> bases[base] ?: 'X' }.join('')
}

/**
* Reverse complement of a nucleotide sequence, with IUPAC ambiguity codes.
*/
def reverseComplement(seq) {
    return complement(seq).reverse()
}
