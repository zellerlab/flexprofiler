#!/usr/bin/env python3
"""Fetch the 16S reference sequences of the IDtaxa tests from the NCBI RefSeq targeted loci (16S, BioProject 33175).

usage: fetch_16s_reference.py <out.fasta.gz>

Writes the first N_PER_GENUS full-length (>= MIN_LENGTH bp) RefSeq 16S sequences of each genus in GENERA, with headers
`<accession> <lineage>`, the lineage being `Root;d__...;p__...;c__...;o__...;f__...;g__...` as in GENERA. The genera are
common gut bacteria from five phyla, two of them in the same phylum (Bacillota). Only sequences containing the 341F and
805R primer sites are kept, as the reads of the tests are simulated from the V3-V4 region between them.
"""
import gzip
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
N_PER_GENUS = 15
MIN_LENGTH = 1300
GENERA = {
    "Escherichia": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae",
    "Bacteroides": "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae",
    "Bifidobacterium": "d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Bifidobacteriales;f__Bifidobacteriaceae",
    "Akkermansia": "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiia;o__Verrucomicrobiales;f__Akkermansiaceae",
    "Streptococcus": "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae",
    "Clostridium": "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae",
}
IUPAC = {"N": "ACGT", "W": "AT", "H": "ACT", "V": "ACG", "R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "B": "CGT", "D": "AGT"}
PRIMER_341F = "CCTACGGGNGGCWGCAG"
PRIMER_805R = "GACTACHVGGGTATCTAATCC"


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTNWHVRYMKSBD", "TGCANWDBYRKMSVH"))


def primer_regex(primer):
    return "".join(f"[{IUPAC[base]}]" if base in IUPAC else base for base in primer)


def eutils(tool, **params):
    url = f"{EUTILS}/{tool}.fcgi?{urllib.parse.urlencode(params)}"
    time.sleep(0.4)  # at most 3 requests per second without an API key
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode()


def parse_fasta(text):
    records = []
    for entry in text.split(">")[1:]:
        header, *lines = entry.splitlines()
        records.append((header.split()[0], "".join(lines).upper()))
    return records


def main(out_path):
    fwd = re.compile(primer_regex(PRIMER_341F))
    rev = re.compile(primer_regex(revcomp(PRIMER_805R)))
    with gzip.open(out_path, "wt") as out:
        for genus, lineage in GENERA.items():
            ids = re.findall(
                r"<Id>(\d+)</Id>",
                eutils("esearch", db="nuccore", term=f"33175[BioProject] AND {genus}[Organism]", retmax=200),
            )
            records = parse_fasta(eutils("efetch", db="nuccore", id=",".join(ids), rettype="fasta", retmode="text"))
            kept = 0
            for accession, seq in sorted(records):
                if len(seq) >= MIN_LENGTH and fwd.search(seq) and rev.search(seq):
                    out.write(f">{accession} Root;{lineage};g__{genus}\n{seq}\n")
                    kept += 1
                    if kept == N_PER_GENUS:
                        break
            print(f"{genus}: {kept} sequences", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
