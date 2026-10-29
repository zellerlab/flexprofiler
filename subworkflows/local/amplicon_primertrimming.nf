//
// Find the primers of amplicon reads with PIMENTO and trim them with cutadapt
//

include { PIMENTO_STD                         } from '../../modules/local/pimento/std/main'
include { CUTADAPT as CUTADAPT_PRIMERS        } from '../../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_READTHROUGH    } from '../../modules/nf-core/cutadapt/main'

workflow AMPLICON_PRIMERTRIMMING {
    take:
    reads // [ [ meta ], [ reads ] ], AMPLICON runs

    main:
    ch_versions = channel.empty()

    // PIMENTO searches each mate of paired-end reads separately: the primers found for the two mates (normally the
    // forward primer for R1 and the reverse one for R2) are the primers of the run
    PIMENTO_STD(reads)
    ch_versions = ch_versions.mix(PIMENTO_STD.out.versions.first())

    ch_reads_for_trimming = reads
        .join(PIMENTO_STD.out.primers)
        .map { meta, input_reads, primer_fastas ->
            def primers = [primer_fastas].flatten().collectMany { fasta -> readPimentoPrimers(fasta, meta.single_end) }
                .collect { _name, seq -> seq }
                .unique()
            [meta + [primers: primers, primers_revcomp: primers.collect { seq -> reverseComplement(seq) }], input_reads]
        }
        .branch { meta, _input_reads ->
            trim: meta.primers
            skip: true
        }

    // All the primers of the run are searched on both mates, so that the pairs sequenced in the other orientation
    // (R1 starting with the reverse primer) are trimmed too: first at the start of the reads, then, at their end, the
    // reverse complement of the primers, for the reads running past the end of the amplicon
    CUTADAPT_PRIMERS(ch_reads_for_trimming.trim)
    CUTADAPT_READTHROUGH(CUTADAPT_PRIMERS.out.reads)

    // the primers are removed from the meta, so that the runs of a sample are merged together whatever their primers
    ch_trimmed_reads = CUTADAPT_READTHROUGH.out.reads
        .mix(ch_reads_for_trimming.skip)
        .map { meta, trimmed_reads -> [meta - meta.subMap(['primers', 'primers_revcomp']), trimmed_reads] }

    ch_read_count_logs = CUTADAPT_PRIMERS.out.log
        .map { meta, log_file -> [meta - meta.subMap(['primers', 'primers_revcomp']), 'cutadapt_primers', log_file] }
        .mix(
            CUTADAPT_READTHROUGH.out.log.map { meta, log_file -> [meta - meta.subMap(['primers', 'primers_revcomp']), 'cutadapt_readthrough', log_file] }
        )

    // the primers found, one row per run and mate (single-end: mate 1), with the primer (if any), its proportion of
    // reads and the region of the library it belongs to
    ch_primers_tsv = PIMENTO_STD.out.out
        .flatMap { meta, out_files ->
            [out_files].flatten().collect { out_file ->
                def mate = (out_file.name =~ /_(\d)_std_primer_out\.txt$/)[0][1]
                def lines = out_file.text.readLines().findAll { line -> line.trim() }
                def region = lines ? lines[0] : ''
                def found = lines.drop(1).collect { line -> line.split(': ')*.trim() }
                if (!found) {
                    found = [['', '']]
                }
                found.collect { name_prop -> [meta.id, meta.run_accession, mate, region, name_prop[0], name_prop[1]].join('\t') }.join('\n')
            }
        }
        .collectFile(
            name: 'primers.tsv',
            seed: ['sample', 'run_accession', 'mate', 'region', 'primer', 'proportion'].join('\t'),
            newLine: true,
            sort: true,
        )
        // copied to the output directory here and not with storeDir, which creates the directory also without amplicon
        // samples
        .map { primers_tsv -> primers_tsv.copyTo("${params.outdir}/amplicon_primers/primers.tsv") }

    emit:
    reads           = ch_trimmed_reads   // channel: [ val(meta), [ reads ] ]
    read_count_logs = ch_read_count_logs // channel: [ val(meta), kind, log ]
    primers         = ch_primers_tsv     // channel: path(primers.tsv)
    versions        = ch_versions        // channel: [ versions.yml ]
}

// Custom Functions

/**
* Read the primers found by PIMENTO (pimento std), as they are in the primer library.
*
* PIMENTO writes the reverse primer of a pair complemented, and also a single reverse primer when searching single-end
* reads (--merged): the complement is undone. In --merged mode PIMENTO searches the reverse primers in the reversed end
* of the reads, so it matches them against their complement (the reverse of their reverse complement), and it writes
* the sequence it searched with (write_std_output in pimento/bin/standard_primer_matching.py, PIMENTO 1.3.3: to check
* again when PIMENTO is updated).
*
* @param fasta The <prefix>_std_primers.fasta file of pimento std, empty if no primer was found.
* @param merged Whether pimento std was run with --merged (single-end reads).
* @return A list of [ name, sequence ].
*/
def readPimentoPrimers(fasta, merged) {
    def records = fasta.splitFasta(record: [id: true, seqString: true])
    return records.withIndex().collect { record, i ->
        def complemented = (records.size() == 2 && i == 1) || (records.size() == 1 && merged && record.id.contains('R'))
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
