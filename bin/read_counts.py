#!/usr/bin/env python3
"""Per-run read counts after QC, host removal and alignment to each profiling database, parsed from the logs.

usage: read_counts.py <manifest.tsv> <inputs_dir> <out.tsv> <out_mqc.tsv>

The manifest has one row per input file, without header: <index> <sample> <run> <kind> <db_name> <annot_db_name>, the
file being the only one in <inputs_dir>/<index>/. The run is the run_accession, or, with run merging, the merged run of a
multi-run sample (merged_<library_strategy>). The kinds are fastp, bowtie2, bbmap_repair, motus4, motus3, cayman and
gqlib, and for amplicons cutadapt_primers, cutadapt_readthrough, vsearch or fastp_merge, and idtaxa. The counts of the
files of a sample and run (e.g. of the runs merged into it) are summed. A column is written only if some file has its
count.
"""
import collections
import json
import pathlib
import re
import sys

PROFILERS = ["motus4", "motus3", "cayman", "gqlib", "idtaxa"]
ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")


def search_int(pattern, text):
    match = re.search(pattern, text, re.MULTILINE)
    return int(match.group(1)) if match else None


def parse_fastp(text):
    # individual reads, each mate counted
    summary = json.loads(text)["summary"]
    return {
        "raw_reads": summary["before_filtering"]["total_reads"],
        "qc_reads": summary["after_filtering"]["total_reads"],
    }


def parse_bowtie2(text):
    # single-end reads, the reads not aligned (--un)
    return {"host_removed_pairs": search_int(r"^\s*(\d+) \([\d.]+%\) aligned 0 times", text)}


def parse_bbmap_repair(text):
    # paired-end reads, the pairs with neither mate aligned to the host, re-paired by repair.sh after aligning the
    # mates separately: repair.sh counts the reads of the pairs, i.e. both mates
    count = search_int(r"^Pairs:\s+(\d+) reads", text)
    return {"host_removed_pairs": None if count is None else count // 2}


def parse_cutadapt(text):
    # individual reads, each mate counted: the reads with a primer (adapter) found and trimmed. cutadapt writes the
    # numbers with thousands separators, and the single-end and paired-end reports differ
    text = text.replace(",", "")
    reads = search_int(r"^Reads with adapters:\s+(\d+)", text)
    if reads is None:
        r1 = search_int(r"^\s*Read 1 with adapter:\s+(\d+)", text)
        r2 = search_int(r"^\s*Read 2 with adapter:\s+(\d+)", text)
        reads = None if r1 is None or r2 is None else r1 + r2
    return reads


def parse_vsearch(text):
    # read pairs, merged or joined into one sequence
    return {
        "merged_pairs": search_int(r"^\s*(\d+)\s+Merged \(", text),
        "joined_pairs": search_int(r"^(\d+) pairs joined", text),
    }


def parse_fastp_merge(text):
    # read pairs merged into one sequence (the mates of the others are kept as separate sequences)
    return {"merged_pairs": search_int(r"^Read pairs merged: (\d+)", text)}


def parse_idtaxa(text):
    # the sequences (read pairs merged or joined, mates of unmerged pairs, or single-end reads) assigned to a taxon
    return search_int(r"^Assigned sequences: (\d+)", text)


def parse_motus4(text):
    # the inserts (read pairs) passing the motus filters
    return search_int(r"Read (\d+) aligned inserts", text)


def parse_motus3(text):
    # the inserts (read pairs) passing the motus filters
    return search_int(r"Total number of inserts: (\d+)", text)


def parse_gene_catalogue(text):
    # gqlib (also used by cayman) logs the counts only if there are alignments
    count = search_int(r"Reads passing filters: (\d+)", text)
    if count is None and re.search(r"Finished\.$", text, re.MULTILINE):
        count = 0
    return count


PARSERS = {
    "motus4": parse_motus4,
    "motus3": parse_motus3,
    "cayman": parse_gene_catalogue,
    "gqlib": parse_gene_catalogue,
    "idtaxa": parse_idtaxa,
}

DESCRIPTIONS = {
    "primer_trimmed_reads": (
        "Primer-trimmed reads",
        "Amplicon reads with a primer found by PIMENTO trimmed at their start (cutadapt), each mate counted",
    ),
    "readthrough_trimmed_reads": (
        "Read-through-trimmed reads",
        "Amplicon reads with the opposite primer trimmed at their end (cutadapt), each mate counted",
    ),
    "raw_reads": ("Raw reads", "Reads before QC (fastp), each mate counted"),
    "qc_reads": ("QC reads", "Reads after QC (fastp), each mate counted"),
    "host_removed_pairs": (
        "Host-removed pairs",
        "Read pairs with neither mate aligned to the host genome (single-end: reads not aligned), given to the"
        " profilers",
    ),
    "motus4": ("mOTUs v4 {db} pairs", "Read pairs aligned to the mOTUs v4 database {db} and passing the mOTUs filters"),
    "motus3": ("mOTUs v3 {db} pairs", "Read pairs aligned to the mOTUs v3 database {db} and passing the mOTUs filters"),
    "cayman": ("cayman {db} pairs", "Read pairs aligned to the gene catalogue {db} and passing the cayman filters"),
    "gqlib": ("gqlib {db} pairs", "Read pairs aligned to the gene catalogue {db} and passing the gqlib filters"),
    "merged_pairs": ("Merged pairs", "Amplicon read pairs merged into one sequence (vsearch or fastp), for IDTAXA"),
    "joined_pairs": (
        "Joined pairs",
        "Amplicon read pairs not merged, joined with a padding of Ns into one sequence (vsearch), for IDTAXA",
    ),
    "idtaxa": (
        "IDTAXA {db} sequences",
        "Amplicon sequences (merged or joined read pairs, mates of unmerged pairs, or single-end reads) assigned to a"
        " taxon by IDTAXA with the training set {db}",
    ),
}


def main(manifest, inputs_dir, out_tsv, out_mqc):
    entries = []
    with open(manifest) as handle:
        for line in handle:
            if line.strip():
                index, sample, run, kind, db_name, annot_db_name = line.rstrip("\n").split("\t")
                (path,) = (pathlib.Path(inputs_dir) / index).iterdir()
                entries.append(((sample, run), kind, db_name, annot_db_name, path))

    # cayman and gqlib write one log per annotation of a gene catalogue, with the same counts: only those of the first
    # annotation are used
    first_annot = {}
    for sample, kind, db_name, annot_db_name, _path in entries:
        key = (sample, kind, db_name)
        first_annot[key] = min(first_annot.get(key, annot_db_name), annot_db_name)

    counts = collections.defaultdict(dict)  # (sample, run) -> column -> count, None if not found in a log
    columns = {}  # column -> (title, description)
    for sample, kind, db_name, annot_db_name, path in entries:
        text = ANSI_ESCAPE.sub("", path.read_text())
        if kind == "fastp":
            values = parse_fastp(text)
        elif kind == "bowtie2":
            values = parse_bowtie2(text)
        elif kind == "bbmap_repair":
            values = parse_bbmap_repair(text)
        elif kind == "cutadapt_primers":
            values = {"primer_trimmed_reads": parse_cutadapt(text)}
        elif kind == "cutadapt_readthrough":
            values = {"readthrough_trimmed_reads": parse_cutadapt(text)}
        elif kind == "vsearch":
            values = parse_vsearch(text)
        elif kind == "fastp_merge":
            values = parse_fastp_merge(text)
        elif kind in PARSERS:
            if annot_db_name != first_annot[(sample, kind, db_name)]:
                continue
            column = f"{kind}_{db_name}_pairs"
            columns[column] = tuple(x.format(db=db_name) for x in DESCRIPTIONS[kind])
            values = {column: PARSERS[kind](text)}
        else:
            sys.exit(f"Unknown kind {kind} for {path}")
        for column, value in values.items():
            if value is None:
                print(
                    f"WARNING: count {column} not found in {path}, it is reported as NA for {sample[0]} {sample[1]}",
                    file=sys.stderr,
                )
            previous = counts[sample].get(column, 0)
            counts[sample][column] = None if previous is None or value is None else previous + value
            columns.setdefault(column, DESCRIPTIONS.get(column))

    # preprocessing first, in the order of the steps, then the profilers in a fixed order
    order = [
        "primer_trimmed_reads",
        "readthrough_trimmed_reads",
        "raw_reads",
        "qc_reads",
        "host_removed_pairs",
        "merged_pairs",
        "joined_pairs",
    ]
    profiler_columns = sorted(
        (c for c in columns if c not in order), key=lambda c: (PROFILERS.index(c.split("_")[0]), c)
    )
    header = [c for c in order if c in columns] + profiler_columns
    samples = sorted(counts)

    def row(sample, missing):
        values = [counts[sample].get(c) for c in header]
        return list(sample) + [missing if v is None else str(v) for v in values]

    with open(out_tsv, "w") as handle:
        handle.write("\t".join(["sample", "run_accession"] + header) + "\n")
        for sample in samples:
            handle.write("\t".join(row(sample, "NA")) + "\n")

    # MultiQC custom content, json strings being valid yaml
    with open(out_mqc, "w") as handle:
        handle.write(
            "# id: 'read_counts'\n"
            "# section_name: 'Read counts'\n"
            "# description: 'Reads after each step. The raw and QC reads count each mate, the following columns count"
            " read pairs (single-end: reads).'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            "#     id: 'read_counts_table'\n"
            "#     title: 'Read counts'\n"
            "#     col1_header: 'Sample'\n"
            "# headers:\n"
        )
        for column in header:
            title, description = columns[column]
            handle.write(
                f"#     {column}:\n"
                f"#         title: {json.dumps(title)}\n"
                f"#         description: {json.dumps(description)}\n"
                "#         format: '{:,.0f}'\n"
                "#         scale: false\n"
            )
        handle.write("\t".join(["sample"] + header) + "\n")
        for sample in samples:
            # MultiQC needs a single id column
            values = row(sample, "")
            handle.write("\t".join([f"{values[0]}_{values[1]}"] + values[2:]) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
