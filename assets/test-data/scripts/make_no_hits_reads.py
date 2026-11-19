#!/usr/bin/env python3
"""Make the reads of the no-hits tests: random sequences, which match no database of the profilers.

usage: make_no_hits_reads.py <out_dir>

Writes to <out_dir>:
- no_hits_{1,2}.fastq.gz: one 2x150 pair of random bases (random_pair, the same name for both mates)
- no_hits_se.fastq.gz: one random 150 bp read (random_single)
The reads are of good quality and without adapters, so that they pass the QC and are given to the profilers, which
must handle a sample of which no read aligns or is classified.
"""
import gzip
import pathlib
import random
import sys

SEED = 42
LENGTH = 150


def write(path, records):
    # no file name and timestamp in the header, so that the output is reproducible
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as handle:
        for name, seq in records:
            handle.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n".encode())


def main(out_dir):
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    r1, r2, se = ("".join(rng.choice("ACGT") for _ in range(LENGTH)) for _ in range(3))
    write(out_dir / "no_hits_1.fastq.gz", [("random_pair", r1)])
    write(out_dir / "no_hits_2.fastq.gz", [("random_pair", r2)])
    write(out_dir / "no_hits_se.fastq.gz", [("random_single", se)])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
