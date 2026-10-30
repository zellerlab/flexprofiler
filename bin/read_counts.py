#!/usr/bin/env python3
"""Per-run read counts after QC, host removal and alignment to each profiling database, parsed from the logs.

usage: read_counts.py <manifest.tsv> <inputs_dir> <out.tsv> <out_mqc.tsv>

The manifest has a header and one row per input file, with the columns file (its name in <inputs_dir>), sample,
run_accession, kind and db_name. The run_accession is that of the run, or, with run merging, the merged run of a multi-
run sample (merged_<library_strategy>). The kinds are fastp, bowtie2, bbmap_repair, motus4, motus3, cayman and gqlib,
and for amplicons cutadapt_primers, cutadapt_readthrough, vsearch or fastp_merge, and idtaxa. The counts of the files of
a sample and run (e.g. of the runs merged into it) are summed. A column is written only if some file has its count.
"""
import json
import pathlib
import re
import sys
import textwrap

import pandas as pd

PROFILERS = ["motus4", "motus3", "cayman", "gqlib", "idtaxa"]
# the columns of the preprocessing steps, in their order
PREPROCESSING_COLUMNS = [
    "primer_trimmed_reads",
    "readthrough_trimmed_reads",
    "raw_reads",
    "qc_reads",
    "hostremoval_retained_pairs",
    "merged_pairs",
    "joined_pairs",
]
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
    return {"hostremoval_retained_pairs": search_int(r"^\s*(\d+) \([\d.]+%\) aligned 0 times", text)}


def parse_bbmap_repair(text, kind, db):
    # paired-end reads, the pairs with neither mate aligned to the host, re-paired by repair.sh after aligning the
    # mates separately: repair.sh counts the reads of the pairs, i.e. both mates
    count = search_int(r"^Pairs:\s+(\d+) reads", text)
    return {"hostremoval_retained_pairs": None if count is None else count // 2}


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
    "hostremoval_retained_pairs": (
        "Host removal retained pairs",
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


def read_log(path):
    return ANSI_ESCAPE.sub("", path.read_text())


def main(manifest, inputs_dir, out_tsv, out_mqc):
    entries = pd.read_csv(manifest, sep="\t", dtype=str, keep_default_na=False)

    # the counts are joined to the entries by file
    duplicated = entries.loc[entries["file"].duplicated(), "file"]
    assert duplicated.empty, f"Files in the manifest more than once: {', '.join(duplicated.unique())}"
    unknown = set(entries["kind"]) - set(PARSERS)
    assert not unknown, f"Unknown kinds {', '.join(sorted(unknown))}"

    entries["path"] = entries["file"].map(lambda file: pathlib.Path(inputs_dir) / file)
    absent = entries.loc[~entries["path"].map(pathlib.Path.is_file), "file"]
    assert absent.empty, f"Files of the manifest not in {inputs_dir}: {', '.join(absent)}"
    # the counts of the files, one row per count, joined to the entries
    counts = pd.DataFrame(
        [
            {"file": e.file, "column": column, "count": count}
            for e in entries.itertuples()
            for column, count in PARSERS[e.kind](read_log(e.path), e.kind, e.db_name).items()
        ]
    ).astype({"count": "Int64"})
    entries = entries.merge(counts, on="file")

    missing = entries[entries["count"].isna()]
    assert missing.empty, "Counts not found in their files:\n" + "\n".join(
        f"{e.column} in {e.path} ({e.sample} {e.run_accession})" for e in missing.itertuples()
    )

    # the counts of a sample and run summed; NA for a column without counts for a sample and run
    counts = entries.groupby(["sample", "run_accession", "column"])["count"].sum().unstack("column")

    # preprocessing first, in the order of the steps, then the profilers in the order of PROFILERS, each by database
    column_kinds = dict(zip(entries["column"], entries["kind"]))
    profiler_columns = sorted(
        (c for c in counts.columns if c not in PREPROCESSING_COLUMNS),
        key=lambda c: (PROFILERS.index(column_kinds[c]), c),
    )
    counts = counts[[c for c in PREPROCESSING_COLUMNS if c in counts.columns] + profiler_columns].astype("Int64")
    counts.columns.name = None
    counts.to_csv(out_tsv, sep="\t", na_rep="NA")

    # the MultiQC header of each column, the profiler descriptions being by kind
    headers = {}
    for e in entries.drop_duplicates("column").itertuples():
        title, description = DESCRIPTIONS.get(e.column) or DESCRIPTIONS[e.kind]
        headers[e.column] = {
            "title": title.format(db=e.db_name),
            "description": description.format(db=e.db_name),
            "format": "{:,.0f}",
            "scale": False,
        }

    # MultiQC custom content: the config as a commented header, json being valid yaml
    config = {
        "id": "read_counts",
        "section_name": "Read counts",
        "description": "Reads after each step. The raw and QC reads count each mate, the following columns count read"
        " pairs (single-end: reads).",
        "plot_type": "table",
        "pconfig": {"id": "read_counts_table", "title": "Read counts", "col1_header": "Sample"},
        "headers": {column: headers[column] for column in counts.columns},
    }
    with open(out_mqc, "w") as handle:
        print(textwrap.indent(json.dumps(config, indent=4), "# "), file=handle)
        # MultiQC needs a single id column
        counts.index = counts.index.map("_".join)
        counts.to_csv(handle, sep="\t", na_rep="", index_label="sample")


if __name__ == "__main__":
    assert len(sys.argv) == 5, __doc__
    main(*sys.argv[1:])
