# Primers added to the PIMENTO library

The FASTA files given with `--amplicon_pimento_primers_pacbio` (default: `V1-V9_PacBio.fasta`) are added to the PIMENTO
standard library for the PacBio (`PACBIO_SMRT`) AMPLICON runs: `pimento std` is run on a copy of the standard library
with them (`-p`). As in the PIMENTO library, each file is a region (its name), and the names of the forward primers end
with `F`, those of the reverse primers with `R`.

`V1-V9_PacBio.fasta`: the full-length 16S primers of PacBio (also used for Kinnex 16S), without the barcodes and the
`GCATC` tail: forward `AGRGTTYGATYMTGGCTCAG` (27F) and reverse `RGYTACCTTGTTACGACTT` (1492R, 19 nt). From PacBio,
"Procedure & checklist – Amplification of bacterial full-length 16S rRNA gene with barcoded primers"
(PN 101-599-700 REV05, May 2022); the same as the defaults of PacBio's pb-16S-nf (`--front_p`, `--adapter_p`).

They are not in the PIMENTO library (1.3.3), whose full-length 16S primers (`V1-V9`) are 27F `AGAGTTTGATCMTGGCTCAG`,
1492R `TACGGYTACCTTGTTACGACTT` and 1492modR `TACGGYTACCTTGTTAYGACTT`. With fuzzy matching, PIMENTO finds those in
reads amplified with the PacBio primers, but 1492R and 1492modR are 3 bases longer than the end of the reads, and
cutadapt, which requires the primers of the long reads in full, discards nearly all of them.

The PacBio primers are in a region of their own, not in `V1-V9`: within a region PIMENTO picks, between primers found in
proportions within 0.03 of each other, the longest, so 1492modR (22 nt) would be picked over the PacBio 1492R (19 nt)
also in the reads that end with the latter. Between regions, PIMENTO picks the pair found in the highest mean
proportion of reads.
