#!/usr/bin/env python3
"""Per-run read counts after QC, host removal and alignment to each profiling database, parsed from the logs.

usage: read_counts.py <manifest.tsv> <inputs_dir> <out.tsv> <out_mqc.tsv>

The manifest has a header and one row per input file, in the order of the files, with the columns sample, run_accession,
kind, db_name and annot_db_name, the file of row <i> (from 1) being the only one in <inputs_dir>/<i>/. The run_accession
is that of the run, or, with run merging, the merged run of a multi-run sample (merged_<library_strategy>). The kinds
are fastp, bowtie2, bbmap_repair, motus4, motus3, cayman and gqlib, and for amplicons cutadapt_primers,
cutadapt_readthrough, vsearch or fastp_merge, and idtaxa. The counts of the files of a sample and run (e.g. of the runs
merged into it) are summed. A column is written only if some file has its count.
"""
import json
import pathlib
import re
import sys

import pandas as pd

PROFILERS = ["motus4", "motus3", "cayman", "gqlib", "idtaxa"]
ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")


def search_int(pattern, text):
    match = re.search(pattern, text, re.MULTILINE)
    return int(match.group(1)) if match else None


def parse_fastp(text, kind, db):
    # individual reads, each mate counted
    summary = json.loads(text)["summary"]
    return {
        "raw_reads": summary["before_filtering"]["total_reads"],
        "qc_reads": summary["after_filtering"]["total_reads"],
    }


def parse_bowtie2(text, kind, db):
    # single-end reads, the reads not aligned (--un)
    return {"host_removed_pairs": search_int(r"^\s*(\d+) \([\d.]+%\) aligned 0 times", text)}


def parse_bbmap_repair(text, kind, db):
    # paired-end reads, the pairs with neither mate aligned to the host, re-paired by repair.sh after aligning the
    # mates separately: repair.sh counts the reads of the pairs, i.e. both mates
    count = search_int(r"^Pairs:\s+(\d+) reads", text)
    return {"host_removed_pairs": None if count is None else count // 2}


def parse_cutadapt(text, kind, db):
    # individual reads, each mate counted: the reads with a primer (adapter) found and trimmed. cutadapt writes the
    # numbers with thousands separators, and the single-end and paired-end reports differ
    text = text.replace(",", "")
    reads = search_int(r"^Reads with adapters:\s+(\d+)", text)
    if reads is None:
        r1 = search_int(r"^\s*Read 1 with adapter:\s+(\d+)", text)
        r2 = search_int(r"^\s*Read 2 with adapter:\s+(\d+)", text)
        reads = None if r1 is None or r2 is None else r1 + r2
    column = {"cutadapt_primers": "primer_trimmed_reads", "cutadapt_readthrough": "readthrough_trimmed_reads"}[kind]
    return {column: reads}


def parse_vsearch(text, kind, db):
    # read pairs, merged or joined into one sequence
    return {
        "merged_pairs": search_int(r"^\s*(\d+)\s+Merged \(", text),
        "joined_pairs": search_int(r"^(\d+) pairs joined", text),
    }


def parse_fastp_merge(text, kind, db):
    # read pairs merged into one sequence (the mates of the others are kept as separate sequences)
    return {"merged_pairs": search_int(r"^Read pairs merged: (\d+)", text)}


def parse_idtaxa(text, kind, db):
    # the sequences (read pairs merged or joined, mates of unmerged pairs, or single-end reads) assigned to a taxon
    return {f"{kind}_{db}_pairs": search_int(r"^Assigned sequences: (\d+)", text)}


def parse_motus4(text, kind, db):
    # the inserts (read pairs) passing the motus filters
    return {f"{kind}_{db}_pairs": search_int(r"Read (\d+) aligned inserts", text)}


def parse_motus3(text, kind, db):
    # the inserts (read pairs) passing the motus filters
    return {f"{kind}_{db}_pairs": search_int(r"Total number of inserts: (\d+)", text)}


def parse_gene_catalogue(text, kind, db):
    # gqlib (also used by cayman) logs the counts only if there are alignments
    count = search_int(r"Reads passing filters: (\d+)", text)
    if count is None and re.search(r"Finished\.$", text, re.MULTILINE):
        count = 0
    return {f"{kind}_{db}_pairs": count}


# kind -> parser, giving the counts by column from the log, the kind and the database (profilers) of the log
PARSERS = {
    "fastp": parse_fastp,
    "bowtie2": parse_bowtie2,
    "bbmap_repair": parse_bbmap_repair,
    "cutadapt_primers": parse_cutadapt,
    "cutadapt_readthrough": parse_cutadapt,
    "vsearch": parse_vsearch,
    "fastp_merge": parse_fastp_merge,
    "motus4": parse_motus4,
    "motus3": parse_motus3,
    "cayman": parse_gene_catalogue,
    "gqlib": parse_gene_catalogue,
    "idtaxa": parse_idtaxa,
}

# column, or, for the profilers, kind -> (title, description), {db} being the database
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
    entries = pd.read_csv(manifest, sep="\t", dtype=str, keep_default_na=False)

    # cayman and gqlib write one log per annotation of a gene catalogue, with the same counts: only those of the first
    # annotation are used
    first_annot = entries.groupby(["sample", "run_accession", "kind", "db_name"])["annot_db_name"].transform("min")
    entries = entries[entries["annot_db_name"] == first_annot]

    rows = []  # one per count: sample, run, column, count (None if not found in a log)
    descriptions = {}  # column -> (title, description)
    for entry in entries.itertuples():
        (path,) = (pathlib.Path(inputs_dir) / str(entry.Index + 1)).iterdir()
        if entry.kind not in PARSERS:
            sys.exit(f"Unknown kind {entry.kind} for {path}")
        values = PARSERS[entry.kind](ANSI_ESCAPE.sub("", path.read_text()), entry.kind, entry.db_name)
        for column, value in values.items():
            if value is None:
                print(
                    f"WARNING: count {column} not found in {path}, it is reported as NA for {entry.sample}"
                    f" {entry.run_accession}",
                    file=sys.stderr,
                )
            rows.append((entry.sample, entry.run_accession, column, value))
            description = DESCRIPTIONS.get(column) or DESCRIPTIONS[entry.kind]
            descriptions.setdefault(column, tuple(x.format(db=entry.db_name) for x in description))

    # the counts of a sample and run summed, NA if any is; NA for a column without counts for a sample and run
    counts = (
        pd.DataFrame(rows, columns=["sample", "run_accession", "column", "count"])
        .astype({"count": "Int64"})
        .groupby(["sample", "run_accession", "column"])["count"]
        .agg(lambda counts: counts.sum(skipna=False))
        .unstack("column")
    )

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
        (c for c in counts.columns if c not in order), key=lambda c: (PROFILERS.index(c.split("_")[0]), c)
    )
    counts = counts[[c for c in order if c in counts.columns] + profiler_columns].astype("Int64")
    counts.columns.name = None
    counts.to_csv(out_tsv, sep="\t", na_rep="NA")

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
        for column in counts.columns:
            title, description = descriptions[column]
            handle.write(
                f"#     {column}:\n"
                f"#         title: {json.dumps(title)}\n"
                f"#         description: {json.dumps(description)}\n"
                "#         format: '{:,.0f}'\n"
                "#         scale: false\n"
            )
        # MultiQC needs a single id column
        counts.index = counts.index.map("_".join)
        counts.to_csv(handle, sep="\t", na_rep="", index_label="sample")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
