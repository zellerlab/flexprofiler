#!/usr/bin/env python3
"""Make the references and the reads of the PathSeq-T2T test, from the host genome (the human mitochondrial genome).

usage: make_pathseq_t2t_data.py <host.fna.gz> <out_reference_dir> <out_fastq_dir>

The three host references of PathSeq-T2T are parts of the host genome, so that each filter has reads that only it
removes. Writes to <out_reference_dir>:
- grch38_test.fna.gz: the "GRCh38" of the test, with the contigs host_grch38 (the first 8000 bp of the host genome) and
  decoy_virus (a random sequence of 5000 bp, standing for a viral decoy contig)
- decoys_test.bed: the decoy contig
- pathseq_host/pathseq_host.fa: the host of PathSeqFilterSpark, the contig host_grch38 alone. Its k-mers and BWA index
  image (pathseq_host.bfi, pathseq_host.fa.img) are made from it with GATK, see make_test_data.sh
- t2t_test.fna.gz: the "T2T-CHM13" of the test, the whole host genome (contig host_t2t): the reads from after its first
  8100 bp are removed only by the alignment to it
- microbe_genomes.fna.gz: the genomes the other reads are from, microbe_genome (a random sequence of 20000 bp) and
  decoy_virus, for the test database of the profilers

and to <out_fastq_dir> error-free reads of 151 bp from fragments of 300 to 400 bp (so that bwa estimates an insert size
distribution with some spread, as in real data, and calls the pairs proper), named by what they are, shuffled:
- lowbiomass_run1_{1,2}.fastq.gz, paired-end reads:
  - microbe_pair_<i> (300): from the microbe genome, kept as pairs
  - decoy_pair_<i> (100): from the decoy contig, aligned to "GRCh38" and kept as pairs (decoy rescue)
  - grch38_pair_<i> (200): from host_grch38, removed by the alignment to "GRCh38"
  - t2t_pair_<i> (200): from the host genome after host_grch38, removed by the alignment to "T2T-CHM13"
  - duplicate_pair_<g>_<k> (10 groups of 2): from the microbe genome, the two pairs of a group being the same: one is
    kept
  - lowcomplexity_pair_<i> (20): repeats of 1 to 3 bases, removed by the complexity filter
  - lowquality_pair_<i> (20): from the microbe genome, with base quality 2, removed by the quality filter
  - short_pair_<i> (20): from the microbe genome, reads of 50 bp, removed by the minimum read length (60)
  - chimera_grch38r1_microbe_<i>, chimera_grch38r2_microbe_<i> (50 each): one mate from host_grch38 (R1 or R2), removed,
    the other from the microbe genome, kept without its mate (singletons)
  - chimera_grch38_t2t_<i> (50): R1 from host_grch38, removed, R2 from the host genome after it, left without its mate
    and removed by the alignment to "T2T-CHM13"
- lowbiomass_run2_{1,2}.fastq.gz, a second run of the same sample:
  - microbe2_pair_<i> (50): from the microbe genome, kept as pairs
  - crossrun_duplicate_<i> (10): the same as microbe_pair_<i> of the first run, removed as duplicates only if the runs
    are merged
  - grch38_pair2_<i> (20): from host_grch38, removed
- lowbiomass_se.fastq.gz, single-end reads: microbe_single_<i> (300) and decoy_single_<i> (100), kept,
  grch38_single_<i> (200), t2t_single_<i> (200), lowcomplexity_single_<i> (20), removed, and duplicate_single_<g>_<k>
  (10 groups of 2), one kept per group

The fragments from the microbe genome and from the decoy contig all start at different positions, so that only the
reads made as duplicates are duplicates.
"""
import gzip
import pathlib
import random
import sys

GRCH38_END, T2T_START = 8000, 8100
MICROBE_LENGTH, DECOY_LENGTH = 20000, 5000
READ, SHORT_READ = 151, 50
MIN_FRAGMENT, MAX_FRAGMENT = 300, 400
SEED = 42


def read_fasta(path):
    with gzip.open(path, "rt") as handle:
        return "".join(line.strip() for line in handle if not line.startswith(">")).upper()


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def fragment_reads(genome, start, rng, length=READ):
    # the two mates of the fragment starting at start, of a random length, read from a random strand
    fragment = genome[start:start + rng.randint(MIN_FRAGMENT, MAX_FRAGMENT)]
    if rng.random() < 0.5:
        fragment = revcomp(fragment)
    return fragment[:length], revcomp(fragment[-length:])


def fastq(name, seq, mate=None, quality="I"):
    return [f"@{name}{'' if mate is None else f' {mate}'}", seq, "+", quality * len(seq)]


def pair(name, seqs, quality="I"):
    return fastq(name, seqs[0], 1, quality), fastq(name, seqs[1], 2, quality)


def write_gz(path, text):
    # mtime=0 so that the files are the same at each run
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as out:
        out.write(text.encode())


def write_fasta(path, contigs):
    text = "".join(f">{name}\n" + "".join(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80)) for name, seq in contigs)
    if str(path).endswith(".gz"):
        write_gz(path, text)
    else:
        pathlib.Path(path).write_text(text)


def write_fastq(path, records):
    write_gz(path, "".join(line + "\n" for record in records for line in record))


def main(host, reference_dir, fastq_dir):
    rng = random.Random(SEED)
    reference_dir, fastq_dir = pathlib.Path(reference_dir), pathlib.Path(fastq_dir)
    (reference_dir / "pathseq_host").mkdir(parents=True, exist_ok=True)
    fastq_dir.mkdir(parents=True, exist_ok=True)

    host_genome = read_fasta(host)
    microbe = "".join(rng.choice("ACGT") for _ in range(MICROBE_LENGTH))
    decoy = "".join(rng.choice("ACGT") for _ in range(DECOY_LENGTH))
    host_grch38 = host_genome[:GRCH38_END]

    write_fasta(reference_dir / "grch38_test.fna.gz", [("host_grch38", host_grch38), ("decoy_virus", decoy)])
    (reference_dir / "decoys_test.bed").write_text(f"decoy_virus\t0\t{DECOY_LENGTH}\n")
    write_fasta(reference_dir / "pathseq_host" / "pathseq_host.fa", [("host_grch38", host_grch38)])
    write_fasta(reference_dir / "t2t_test.fna.gz", [("host_t2t", host_genome)])
    write_fasta(reference_dir / "microbe_genomes.fna.gz", [("microbe_genome", microbe), ("decoy_virus", decoy)])

    # fragment starts: unique for the microbe genome and the decoy contig, taken in order from a shuffled list
    microbe_starts = iter(rng.sample(range(MICROBE_LENGTH - MAX_FRAGMENT), 1500))
    decoy_starts = iter(rng.sample(range(DECOY_LENGTH - MAX_FRAGMENT), 300))

    def from_microbe(length=READ):
        return fragment_reads(microbe, next(microbe_starts), rng, length)

    def from_decoy():
        return fragment_reads(decoy, next(decoy_starts), rng)

    def from_grch38():
        return fragment_reads(host_grch38, rng.randrange(GRCH38_END - MAX_FRAGMENT), rng)

    def from_t2t():
        return fragment_reads(host_genome, rng.randrange(T2T_START, len(host_genome) - MAX_FRAGMENT), rng)

    def low_complexity():
        unit = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 3)))
        return (unit * READ)[:READ]

    # the first run
    microbe_pairs = [from_microbe() for _ in range(300)]
    run1 = [pair(f"microbe_pair_{i}", seqs) for i, seqs in enumerate(microbe_pairs)]
    run1 += [pair(f"decoy_pair_{i}", from_decoy()) for i in range(100)]
    run1 += [pair(f"grch38_pair_{i}", from_grch38()) for i in range(200)]
    run1 += [pair(f"t2t_pair_{i}", from_t2t()) for i in range(200)]
    for group in range(10):
        seqs = from_microbe()
        run1 += [pair(f"duplicate_pair_{group}_{copy}", seqs) for copy in range(2)]
    run1 += [pair(f"lowcomplexity_pair_{i}", (low_complexity(), low_complexity())) for i in range(20)]
    run1 += [pair(f"lowquality_pair_{i}", from_microbe(), quality="#") for i in range(20)]
    run1 += [pair(f"short_pair_{i}", from_microbe(SHORT_READ)) for i in range(20)]
    run1 += [pair(f"chimera_grch38r1_microbe_{i}", (from_grch38()[0], from_microbe()[1])) for i in range(50)]
    run1 += [pair(f"chimera_grch38r2_microbe_{i}", (from_microbe()[0], from_grch38()[1])) for i in range(50)]
    run1 += [pair(f"chimera_grch38_t2t_{i}", (from_grch38()[0], from_t2t()[1])) for i in range(50)]
    rng.shuffle(run1)

    # the second run of the sample
    run2 = [pair(f"microbe2_pair_{i}", from_microbe()) for i in range(50)]
    run2 += [pair(f"crossrun_duplicate_{i}", microbe_pairs[i]) for i in range(10)]
    run2 += [pair(f"grch38_pair2_{i}", from_grch38()) for i in range(20)]
    rng.shuffle(run2)

    # the single-end reads
    single = [fastq(f"microbe_single_{i}", from_microbe()[0]) for i in range(300)]
    single += [fastq(f"decoy_single_{i}", from_decoy()[0]) for i in range(100)]
    single += [fastq(f"grch38_single_{i}", from_grch38()[0]) for i in range(200)]
    single += [fastq(f"t2t_single_{i}", from_t2t()[0]) for i in range(200)]
    single += [fastq(f"lowcomplexity_single_{i}", low_complexity()) for i in range(20)]
    for group in range(10):
        seq = from_microbe()[0]
        single += [fastq(f"duplicate_single_{group}_{copy}", seq) for copy in range(2)]
    rng.shuffle(single)

    for name, pairs in (("lowbiomass_run1", run1), ("lowbiomass_run2", run2)):
        write_fastq(fastq_dir / f"{name}_1.fastq.gz", [first for first, _ in pairs])
        write_fastq(fastq_dir / f"{name}_2.fastq.gz", [second for _, second in pairs])
    write_fastq(fastq_dir / "lowbiomass_se.fastq.gz", single)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
