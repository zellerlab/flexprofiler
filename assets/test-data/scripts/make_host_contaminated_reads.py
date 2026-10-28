#!/usr/bin/env python3
"""Make the host-contaminated reads of the host removal test, from the reads of sample 2613 and the host genome.

usage: make_host_contaminated_reads.py <R1.fastq.gz> <R2.fastq.gz> <host.fna.gz> <out_dir>

Writes to <out_dir>:
- host_contaminated_1.fastq.gz, host_contaminated_2.fastq.gz: paired-end reads, made of
  - the first 2000 pairs of the input (non-host pairs, with their names)
  - 200 host pairs, both mates from the host genome (host_pair_<i>)
  - 100 chimeric pairs with R1 from the host genome and R2 from an input pair (chimera_r1host_<i>)
  - 100 chimeric pairs with R2 from the host genome and R1 from an input pair (chimera_r2host_<i>)
  shuffled. The input pairs used for the chimeras come after the first 2000.
- host_contaminated_se.fastq.gz: single-end reads, the R1 of the 2000 non-host pairs and 200 host reads
  (host_single_<i>), shuffled.

The host reads are error-free fragments of 350 bp of the host genome, from a random strand, read from both ends with
the read length of the input. All the host-derived reads have names starting with host_ or chimera_, which must not
be left after host removal, and all the other pairs should be kept.
"""
import gzip
import random
import sys
from itertools import islice

N_NONHOST, N_HOST, N_CHIMERA = 2000, 200, 100
FRAGMENT = 350
SEED = 42


def records(path):
    with gzip.open(path, "rt") as handle:
        while True:
            record = [line.rstrip("\n") for line in islice(handle, 4)]
            if not record:
                return
            yield record


def read_fasta(path):
    with gzip.open(path, "rt") as handle:
        return "".join(line.strip() for line in handle if not line.startswith(">")).upper()


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def host_pair(genome, length, rng):
    start = rng.randrange(len(genome) - FRAGMENT)
    fragment = genome[start:start + FRAGMENT]
    if rng.random() < 0.5:
        fragment = revcomp(fragment)
    return fragment[:length], revcomp(fragment[-length:])


def fastq(name, seq, comment=""):
    return [f"@{name}{comment}", seq, "+", "I" * len(seq)]


def main(r1, r2, host, out_dir):
    rng = random.Random(SEED)
    genome = read_fasta(host)
    pairs = list(islice(zip(records(r1), records(r2)), N_NONHOST + 2 * N_CHIMERA))
    length = len(pairs[0][0][1])

    pe = [(a, b) for a, b in pairs[:N_NONHOST]]
    for i in range(N_HOST):
        s1, s2 = host_pair(genome, length, rng)
        pe.append((fastq(f"host_pair_{i}", s1, " 1"), fastq(f"host_pair_{i}", s2, " 2")))
    for i, (a, b) in enumerate(pairs[N_NONHOST:N_NONHOST + N_CHIMERA]):
        s1, _ = host_pair(genome, length, rng)
        pe.append((fastq(f"chimera_r1host_{i}", s1, " 1"), [f"@chimera_r1host_{i} 2"] + b[1:]))
    for i, (a, b) in enumerate(pairs[N_NONHOST + N_CHIMERA:]):
        _, s2 = host_pair(genome, length, rng)
        pe.append(([f"@chimera_r2host_{i} 1"] + a[1:], fastq(f"chimera_r2host_{i}", s2, " 2")))
    rng.shuffle(pe)

    se = [a for a, _ in pairs[:N_NONHOST]]
    se += [fastq(f"host_single_{i}", host_pair(genome, length, rng)[0]) for i in range(N_HOST)]
    rng.shuffle(se)

    # mtime=0 so that the files are the same at each run
    def write(name, recs):
        with open(f"{out_dir}/{name}", "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as out:
            out.write("".join(line + "\n" for rec in recs for line in rec).encode())

    write("host_contaminated_1.fastq.gz", [a for a, _ in pe])
    write("host_contaminated_2.fastq.gz", [b for _, b in pe])
    write("host_contaminated_se.fastq.gz", se)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
