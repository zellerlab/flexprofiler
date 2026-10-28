#!/usr/bin/env python3
"""Make the 16S amplicon reads of the amplicon tests, from the reference sequences of fetch_16s_reference.py.

usage: make_amplicon_reads.py <reference_16s.fasta.gz> <out_dir>

The reads are error-free, simulated from the V3-V4 amplicon (341F-805R) of each reference sequence, with names
`<kind>_<genus>_<i>` telling what each read is (the genus being the truth for the classification). Writes to <out_dir>:
- amplicon_clean_{1,2}.fastq.gz: primers already removed (the insert is the region between the primers), mostly
  overlapping 2x250 pairs (ovl_), and some non-overlapping 2x150 pairs (noovl_); amplicon_clean_se.fastq.gz: their R1.
- amplicon_primers_{1,2}.fastq.gz: 2x250 pairs of the whole amplicon, starting with 341F (R1) and 805R (R2) (pr_), some
  after a 1-7 base spacer (spacer_), and some short-insert pairs reading through the end of the amplicon into the
  reverse complement of the other primer and the adapter (rt_); amplicon_primers_se.fastq.gz: their R1.
- amplicon_mixed_{1,2}.fastq.gz: as amplicon_primers, with a quarter of the pairs flipped (flip_: R1 starting with 805R
  and R2 with 341F), as in libraries ligated in both orientations.
"""
import gzip
import pathlib
import random
import re
import sys

SEED = 42
N_OVL, N_NOOVL, N_PRIMERS, N_SPACER, N_RT = 5, 2, 5, 2, 1  # pairs per reference sequence
LONG, SHORT = 250, 150
RT_INSERT = 120  # region length of the read-through pairs
IUPAC = {"N": "ACGT", "W": "AT", "H": "ACT", "V": "ACG", "R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "B": "CGT", "D": "AGT"}
PRIMER_341F = "CCTACGGGNGGCWGCAG"
PRIMER_805R = "GACTACHVGGGTATCTAATCC"
# the Illumina adapters read after the insert, in R1 and R2
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTNWHVRYMKSBD", "TGCANWDBYRKMSVH"))


def primer_regex(primer):
    return re.compile("".join(f"[{IUPAC[base]}]" if base in IUPAC else base for base in primer))


def read_fasta(path):
    with gzip.open(path, "rt") as handle:
        for entry in handle.read().split(">")[1:]:
            header, *lines = entry.splitlines()
            yield header, "".join(lines).upper()


def amplicons(path):
    """The primer sites and the region between them, of each reference: (genus, fwd site, region, rev site)."""
    fwd, rev = primer_regex(PRIMER_341F), primer_regex(revcomp(PRIMER_805R))
    for header, seq in read_fasta(path):
        genus = header.split(";g__")[1]
        f = fwd.search(seq)
        r = rev.search(seq, f.end())
        # the rev site is reported as the primer, i.e. reverse complemented
        yield genus, f.group(), seq[f.end():r.start()], revcomp(r.group())


def pair(insert_fwd, insert_rev, length, rng, adapters=True):
    """The two reads of an insert, given as read from either end, padded past its end with adapter and random bases."""
    pad = lambda adapter: (adapter if adapters else "") + "".join(rng.choice("ACGT") for _ in range(length))
    return (insert_fwd + pad(ADAPTER_R1))[:length], (insert_rev + pad(ADAPTER_R2))[:length]


def spacer(rng):
    return "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 7)))


def write(path, records):
    with gzip.open(path, "wt", compresslevel=9) as handle:
        for name, seq in records:
            handle.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def write_pairs(out_dir, prefix, pairs, se=True):
    write(out_dir / f"{prefix}_1.fastq.gz", [(f"{name}/1", r1) for name, r1, _ in pairs])
    write(out_dir / f"{prefix}_2.fastq.gz", [(f"{name}/2", r2) for name, _, r2 in pairs])
    if se:
        write(out_dir / f"{prefix}_se.fastq.gz", [(name, r1) for name, r1, _ in pairs])


def main(reference, out_dir):
    rng = random.Random(SEED)
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    clean, primers = [], []
    counter = 0
    for genus, fwd_site, region, rev_site in amplicons(reference):
        amplicon = fwd_site + region + revcomp(rev_site)
        for kind, n in (("ovl", N_OVL), ("noovl", N_NOOVL)):
            for _ in range(n):
                counter += 1
                r1, r2 = pair(region, revcomp(region), LONG if kind == "ovl" else SHORT, rng)
                clean.append((f"{kind}_{genus}_{counter}", r1, r2))
        for _ in range(N_PRIMERS):
            counter += 1
            primers.append((f"pr_{genus}_{counter}", *pair(amplicon, revcomp(amplicon), LONG, rng)))
        for _ in range(N_SPACER):
            counter += 1
            primers.append((f"spacer_{genus}_{counter}", *pair(spacer(rng) + amplicon, spacer(rng) + revcomp(amplicon), LONG, rng)))
        for _ in range(N_RT):
            counter += 1
            short = fwd_site + region[:RT_INSERT] + revcomp(rev_site)
            primers.append((f"rt_{genus}_{counter}", *pair(short, revcomp(short), LONG, rng)))
    rng.shuffle(clean)
    rng.shuffle(primers)
    write_pairs(out_dir, "amplicon_clean", clean)
    write_pairs(out_dir, "amplicon_primers", primers)
    mixed = [
        (f"flip_{name.split('_', 1)[1]}", r2, r1) if i % 4 == 0 else (name, r1, r2)
        for i, (name, r1, r2) in enumerate(primers)
    ]
    write_pairs(out_dir, "amplicon_mixed", mixed, se=False)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
