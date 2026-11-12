#!/usr/bin/env python3
"""Make the GTF and the RNA-seq reads of the PathSeq-T2T test with STAR (--pathseq_t2t_grch38_aligner star), from the
references of the PathSeq-T2T test (make_pathseq_t2t_data.py), which are left unchanged.

usage: make_pathseq_t2t_rna_data.py <reference_dir> <out_fastq_dir>

Reads from <reference_dir> grch38_test.fna.gz (host_grch38, decoy_virus), t2t_test.fna.gz (host_t2t) and
microbe_genomes.fna.gz (microbe_genome), and writes:
- <reference_dir>/grch38_test.gtf: two genes on host_grch38 (1-based coordinates):
  - gene_spliced, + strand, two exons (1001-1600 and 3001-3600) and an intron between them
  - gene_unspliced, - strand, one exon (5001-6500)
- <out_fastq_dir>/lowbiomass_rna_{1,2}.fastq.gz, paired-end reads of a reverse stranded library (read 1 antisense, read
  2 sense, as TruSeq stranded), shuffled:
  - spliced_pair_<i> (100): fragments of the mature transcript of gene_spliced that span its exon junction, removed by
    the alignment to "GRCh38" (STAR aligns them across the intron), counted for gene_spliced
  - gene_pair_<i> (60): from the transcript of gene_unspliced, removed by the alignment to "GRCh38", counted for it
  - microbe_pair_<i> (200): from the microbe genome, kept as pairs
  - decoy_pair_<i> (50): from the decoy contig, aligned to "GRCh38" and kept as pairs (decoy rescue)
  - t2t_pair_<i> (100): from the host genome after host_grch38, removed by the alignment to "T2T-CHM13"
- <out_fastq_dir>/lowbiomass_rna_se.fastq.gz, single-end reads (antisense): spliced_single_<i> (100) and
  gene_single_<i> (60), counted for their gene and removed, microbe_single_<i> (200) and decoy_single_<i> (50), kept,
  t2t_single_<i> (100), removed

Error-free reads of 151 bp from fragments of 350 bp. The fragments of each class of a sample all start at different
positions, so that no read is a duplicate. The two genes have different numbers of reads, so that the counts per gene
show that the reads of each are counted for it.
"""
import gzip
import pathlib
import random
import sys

READ, FRAGMENT = 151, 350
T2T_START = 8100
SEED = 7
# 0-based, end-exclusive
SPLICED_EXONS = [(1000, 1600), (3000, 3600)]
UNSPLICED_EXON = (5000, 6500)


def read_fasta(path):
    contigs, name = {}, None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                contigs[name] = []
            else:
                contigs[name].append(line.upper())
    return {name: "".join(seqs) for name, seqs in contigs.items()}


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def fastq(name, seq, mate=None):
    return [f"@{name}{'' if mate is None else f' {mate}'}", seq, "+", "I" * len(seq)]


def write_gz(path, text):
    # mtime=0 so that the files are the same at each run
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as out:
        out.write(text.encode())


def write_fastq(path, records):
    write_gz(path, "".join(line + "\n" for record in records for line in record))


def gtf_line(feature, start, end, strand, gene, exon=None):
    attributes = f'gene_id "{gene}"; transcript_id "{gene}.1"; gene_name "{gene}";'
    if exon is not None:
        attributes += f" exon_number {exon};"
    return "\t".join(["host_grch38", "test", feature, str(start + 1), str(end), ".", strand, ".", attributes])


def main(reference_dir, fastq_dir):
    rng = random.Random(SEED)
    reference_dir, fastq_dir = pathlib.Path(reference_dir), pathlib.Path(fastq_dir)
    fastq_dir.mkdir(parents=True, exist_ok=True)

    grch38 = read_fasta(reference_dir / "grch38_test.fna.gz")
    host_grch38, decoy = grch38["host_grch38"], grch38["decoy_virus"]
    host_genome = read_fasta(reference_dir / "t2t_test.fna.gz")["host_t2t"]
    microbe = read_fasta(reference_dir / "microbe_genomes.fna.gz")["microbe_genome"]

    lines = [gtf_line("gene", SPLICED_EXONS[0][0], SPLICED_EXONS[-1][1], "+", "gene_spliced"),
             gtf_line("transcript", SPLICED_EXONS[0][0], SPLICED_EXONS[-1][1], "+", "gene_spliced")]
    lines += [gtf_line("exon", start, end, "+", "gene_spliced", i + 1) for i, (start, end) in enumerate(SPLICED_EXONS)]
    lines += [gtf_line(feature, *UNSPLICED_EXON, "-", "gene_unspliced") for feature in ("gene", "transcript")]
    lines += [gtf_line("exon", *UNSPLICED_EXON, "-", "gene_unspliced", 1)]
    (reference_dir / "grch38_test.gtf").write_text("".join(line + "\n" for line in lines))

    # the mature transcripts, in the sense orientation
    spliced = "".join(host_grch38[start:end] for start, end in SPLICED_EXONS)
    junction = SPLICED_EXONS[0][1] - SPLICED_EXONS[0][0]
    unspliced = revcomp(host_grch38[UNSPLICED_EXON[0]:UNSPLICED_EXON[1]])

    def sense_fragments(seq, starts, n):
        # n fragments in the sense orientation, from distinct starts
        return [seq[start:start + FRAGMENT] for start in rng.sample(starts, n)]

    def random_strand_fragments(seq, starts, n):
        return [frag if rng.random() < 0.5 else revcomp(frag) for frag in sense_fragments(seq, starts, n)]

    def classes(n_host, n_gene, n_microbe, n_decoy):
        # the fragments of each class, as [ (class, fragment) ], in the sense orientation for the genes
        # the fragments of the spliced transcript span the junction, at least 20 bases on each side
        spliced_starts = list(range(junction - FRAGMENT + 20, junction - 20))
        return (
            [("spliced", frag) for frag in sense_fragments(spliced, spliced_starts, n_host)]
            + [("gene", frag) for frag in sense_fragments(unspliced, range(len(unspliced) - FRAGMENT), n_gene)]
            + [("microbe", frag) for frag in random_strand_fragments(microbe, range(len(microbe) - FRAGMENT), n_microbe)]
            + [("decoy", frag) for frag in random_strand_fragments(decoy, range(len(decoy) - FRAGMENT), n_decoy)]
            + [("t2t", frag) for frag in random_strand_fragments(host_genome, range(T2T_START, len(host_genome) - FRAGMENT), n_host)]
        )

    def numbered(fragments, kind):
        counts, named = {}, []
        for cls, frag in fragments:
            i = counts.get(cls, 0)
            counts[cls] = i + 1
            named.append((f"{cls}_{kind}_{i}", frag))
        return named

    # paired-end, reverse stranded: read 1 is the antisense end of the fragment, read 2 its sense start
    pairs = numbered(classes(100, 60, 200, 50), "pair")
    rng.shuffle(pairs)
    write_fastq(fastq_dir / "lowbiomass_rna_1.fastq.gz", [fastq(name, revcomp(frag[-READ:]), 1) for name, frag in pairs])
    write_fastq(fastq_dir / "lowbiomass_rna_2.fastq.gz", [fastq(name, frag[:READ], 2) for name, frag in pairs])

    # single-end, reverse stranded: the antisense end of the fragment
    singles = numbered(classes(100, 60, 200, 50), "single")
    rng.shuffle(singles)
    write_fastq(fastq_dir / "lowbiomass_rna_se.fastq.gz", [fastq(name, revcomp(frag[-READ:])) for name, frag in singles])


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
