#!/usr/bin/env python3
"""Make the adapter-contaminated shotgun reads of the adapter removal test, from the human mitochondrial genome.

usage: make_adapter_reads.py <genome.fna.gz> <out_dir>

The reads are error-free 2x150 pairs of fragments of the genome, from a random strand. For each library kit (TruSeq,
Nextera), writes to <out_dir>:
- adapters_<kit>_{1,2}.fastq.gz: N_CLEAN pairs of 300 bp fragments, longer than the reads (clean_<i>), and N_ADAPTER
  pairs of fragments shorter than the reads, between 40 and 140 bp, so that the reads continue into the adapter of the
  kit and then random bases (adapter_<insert length>_<i>), shuffled
- adapters_<kit>_se.fastq.gz: the R1 of the pairs
and adapters_truseq_se_large.fastq.gz, the R1 of N_CLEAN_LARGE and N_ADAPTER_LARGE TruSeq pairs made in the same way
(with another seed), as fastp detects the adapters of single-end reads only in files of at least 10000 reads.
After adapter removal, the reads must have no adapter left, and the reads of the adapter_ pairs must be as long as
their insert.
"""
import gzip
import pathlib
import random
import sys

SEED = 42
N_CLEAN, N_ADAPTER = 700, 300
N_CLEAN_LARGE, N_ADAPTER_LARGE = 8400, 3600
LENGTH = 150
CLEAN_FRAGMENT = 300
ADAPTERS = {
    "truseq": ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
    "nextera": ("CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", "CTGTCTCTTATACACATCTGACGCTGCCGACGAA"),
}


def read_fasta(path):
    with gzip.open(path, "rt") as handle:
        return "".join(line.strip() for line in handle if not line.startswith(">")).upper()


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def fragment(genome, length, rng):
    start = rng.randrange(len(genome) - length)
    seq = genome[start:start + length]
    return revcomp(seq) if rng.random() < 0.5 else seq


def read(insert, adapter, rng):
    return (insert + adapter + "".join(rng.choice("ACGT") for _ in range(LENGTH)))[:LENGTH]


def write(path, records):
    # no file name and timestamp in the header, so that the output is reproducible
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as handle:
        for name, seq in records:
            handle.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n".encode())


def make_pairs(genome, adapter_1, adapter_2, n_clean, n_adapter, rng):
    pairs = []
    for i in range(n_clean):
        insert = fragment(genome, CLEAN_FRAGMENT, rng)
        pairs.append((f"clean_{i}", insert[:LENGTH], revcomp(insert)[:LENGTH]))
    for i in range(n_adapter):
        insert = fragment(genome, rng.randint(40, 140), rng)
        pairs.append((f"adapter_{len(insert)}_{i}", read(insert, adapter_1, rng), read(revcomp(insert), adapter_2, rng)))
    rng.shuffle(pairs)
    return pairs


def main(genome_path, out_dir):
    rng = random.Random(SEED)
    genome = read_fasta(genome_path)
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for kit, (adapter_1, adapter_2) in ADAPTERS.items():
        pairs = make_pairs(genome, adapter_1, adapter_2, N_CLEAN, N_ADAPTER, rng)
        write(out_dir / f"adapters_{kit}_1.fastq.gz", [(f"{name}/1", r1) for name, r1, _ in pairs])
        write(out_dir / f"adapters_{kit}_2.fastq.gz", [(f"{name}/2", r2) for name, _, r2 in pairs])
        write(out_dir / f"adapters_{kit}_se.fastq.gz", [(name, r1) for name, r1, _ in pairs])
    pairs = make_pairs(genome, *ADAPTERS["truseq"], N_CLEAN_LARGE, N_ADAPTER_LARGE, random.Random(SEED + 1))
    write(out_dir / "adapters_truseq_se_large.fastq.gz", [(name, r1) for name, r1, _ in pairs])


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
