#!/usr/bin/env python3
"""Per-run read counts after QC, host removal and alignment to each profiling database, parsed from the logs.

usage: read_counts.py <manifest.tsv> <inputs_dir> <out.tsv> <out_mqc.tsv>

The manifest has a header and one row per input file, with the columns file (its name in <inputs_dir>), sample,
run_accession, kind and db_name. The run_accession is that of the run, or, with run merging, the merged run of a multi-
run sample (merged_<library_strategy>). The kinds are falco_raw and falco_processed (or fastqc_raw and fastqc_processed:
the reads of the samplesheet and those given to the profilers), bbmap_repair_pe_order, fastp, bowtie2, bbmap_repair,
motus4, motus3, cayman and gqlib, and for amplicons cutadapt_primers, cutadapt_readthrough, vsearch or fastp_merge, and
idtaxa, and for long amplicon reads vsearch_fastqfilter and bbduk_pacbio (PacBio QC), vsearch_orient and
cutadapt_linked, and for the PathSeq-T2T host filtering pathseq_flagstat, pathseq_filter and pathseq_t2tfilter. The
counts of the files of a sample and run (e.g. of the runs merged into it) are summed. A column of a step is written only
if some file has its count; the summary columns (raw_reads, qc_reads, hostremoval_reads, profiler_input_reads) are always
written, computed from the counts of the steps.
"""
import json
import pathlib
import re
import sys
import textwrap
import zipfile

import pandas as pd

PROFILERS = ["motus4", "motus3", "cayman", "gqlib", "kraken2", "kraken2_pairs", "idtaxa"]
# the summary columns, each mate counted, computed from those of the steps (summed, NA if none of them was run):
# the raw reads, after the QC, after the host removal, and given to the profilers
SUMMARY_COLUMNS = {
    "raw_reads": ["falco_raw_reads", "fastqc_raw_reads"],
    "qc_reads": ["shortread_qc_reads", "longread_qc_reads"],
    "hostremoval_reads": ["hostremoval_retained_reads", "t2t_filter_reads"],
    "profiler_input_reads": ["falco_processed_reads", "fastqc_processed_reads"],
}
# the columns of the preprocessing steps, in their order
PREPROCESSING_COLUMNS = [
    *SUMMARY_COLUMNS,
    "falco_raw_reads",
    "fastqc_raw_reads",
    "repair_input_reads",
    "repair_paired_reads",
    "shortread_qc_input_reads",
    "shortread_qc_reads",
    "longread_qc_input_reads",
    "longread_qc_reads",
    "oriented_reads",
    "primer_trimmed_reads",
    "readthrough_trimmed_reads",
    "hostremoval_retained_reads",
    "hostremoval_retained_pairs",
    "pathseq_t2t_input_reads",
    "grch38_unaligned_reads",
    "prealigned_host_filter_reads",
    "quality_complexity_filter_reads",
    "host_kmer_filter_reads",
    "deduplication_filter_reads",
    "t2t_filter_reads",
    "merged_pairs",
    "joined_pairs",
    "falco_processed_reads",
    "fastqc_processed_reads",
]
ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")


def search_int(pattern, text):
    match = re.search(pattern, text, re.MULTILINE)
    return int(match.group(1)) if match else None


def parse_fastp(text, kind, db):
    # individual reads, each mate counted
    summary = json.loads(text)["summary"]
    return {
        "shortread_qc_input_reads": summary["before_filtering"]["total_reads"],
        "shortread_qc_reads": summary["after_filtering"]["total_reads"],
    }


def parse_bowtie2(text, kind, db):
    # single-end reads, the reads not aligned (--un)
    count = search_int(r"^\s*(\d+) \([\d.]+%\) aligned 0 times", text)
    return {"hostremoval_retained_reads": count, "hostremoval_retained_pairs": count}


def parse_bbmap_repair(text, kind, db):
    # paired-end reads, the pairs with neither mate aligned to the host, re-paired by repair.sh after aligning the
    # mates separately: repair.sh counts the reads of the pairs, i.e. both mates
    count = search_int(r"^Pairs:\s+(\d+) reads", text)
    return {"hostremoval_retained_reads": count, "hostremoval_retained_pairs": None if count is None else count // 2}


def parse_reads_qc(text, kind, db):
    # the reads of a FASTQ file (Total Sequences of the report of Falco or FastQC), of the raw reads or of those given
    # to the profilers (kind <tool>_<stage>): a paired-end run has a report per mate
    return {f"{kind}_reads": search_int(r"^Total Sequences\t(\d+)$", text)}


def parse_bbmap_repair_pe_order(text, kind, db):
    # paired-end reads, before the QC: the reads given to repair.sh, and those it pairs (the reads without their mate
    # are dropped), each mate counted
    return {
        "repair_input_reads": search_int(r"^Input:\s+(\d+) reads", text),
        "repair_paired_reads": search_int(r"^Pairs:\s+(\d+) reads", text),
    }


def parse_pathseq_flagstat(text, kind, db):
    # samtools flagstat (tsv) of the alignment to GRCh38, before any filter: the primary reads, aligned or not, i.e. the
    # reads given to PathSeq-T2T, each mate counted
    return {"pathseq_t2t_input_reads": search_int(r"^(\d+)\t\d+\tprimary$", text)}


# the columns of the filter metrics of GATK PathSeqFilterSpark -> those of the read counts, in the order of the filters
PATHSEQ_FILTER_COLUMNS = {
    "PRIMARY_READS": "grch38_unaligned_reads",
    "READS_AFTER_PREALIGNED_HOST_FILTER": "prealigned_host_filter_reads",
    "READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER": "quality_complexity_filter_reads",
    "READS_AFTER_HOST_FILTER": "host_kmer_filter_reads",
    "READS_AFTER_DEDUPLICATION": "deduplication_filter_reads",
}


def parse_pathseq_filter(text, kind, db):
    # the reads left after each filter of PathSeqFilterSpark, each mate counted: a table with one row, after the
    # header lines (#). There is a file for the reads not aligned to GRCh38 and one for those aligned to its decoy
    # contigs, whose counts are summed
    rows = [line.split("\t") for line in text.splitlines() if line.strip() and not line.startswith("#")]
    metrics = dict(zip(*rows)) if len(rows) == 2 and len(rows[0]) == len(rows[1]) else {}
    return {
        column: int(metrics[name]) if metrics.get(name, "").isdigit() else None
        for name, column in PATHSEQ_FILTER_COLUMNS.items()
    }


def parse_pathseq_t2tfilter(text, kind, db):
    # the reads not aligned to T2T-CHM13, given to the profilers, each mate counted: the pairs, and the reads left
    # without their mate by the filters (single-end: all the reads)
    pairs = search_int(r"^Retained pairs: (\d+)$", text)
    unpaired = search_int(r"^Retained unpaired reads: (\d+)$", text)
    return {"t2t_filter_reads": None if pairs is None or unpaired is None else 2 * pairs + unpaired}


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


def parse_cutadapt_linked(text, kind, db):
    # long reads with both primers (linked adapter, the others discarded) and in the length window: the reads written
    text = text.replace(",", "")
    return {"primer_trimmed_reads": search_int(r"^Reads written \(passing filters\):\s+(\d+)", text)}


def parse_vsearch_fastqfilter(text, kind, db):
    # PacBio reads, before the quality filters (vsearch, then BBDuk, which gives the reads after them)
    kept = search_int(r"^(\d+) sequences kept", text)
    discarded = search_int(r"(\d+) sequences discarded", text)
    return {"longread_qc_input_reads": None if kept is None or discarded is None else kept + discarded}


def parse_bbduk_pacbio(text, kind, db):
    # PacBio reads, after the quality filters (vsearch, then BBDuk)
    return {"longread_qc_reads": search_int(r"^Result:\s+(\d+) reads", text)}


def parse_vsearch_orient(text, kind, db):
    # long reads oriented as the reference (the others are dropped)
    return {"oriented_reads": search_int(r"^All oriented sequences:\s+(\d+)", text)}


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


def parse_kraken2(text, kind, db):
    # the reads classified by Kraken2: the fragments of the clade of the root (taxid 1, the second column), each mate
    # counted, as a pair is one fragment for Kraken2 (kraken2_pairs), counted twice. An empty report, or one with no
    # root (no fragment classified), is of no reads
    fragments = search_int(r"^[^\t\n]*\t *(\d+)\t(?:[^\t\n]*\t)*1\troot$", text) or 0
    return {f"kraken2_{db}_reads": fragments * (2 if kind == "kraken2_pairs" else 1)}


def parse_motus4(text, kind, db):
    # the inserts (read pairs, and the singletons of --pathseq_t2t) passing the motus filters
    return {f"{kind}_{db}_pairs": search_int(r"Read (\d+) aligned inserts", text)}


def parse_motus3(text, kind, db):
    # the inserts (read pairs, and the singletons of --pathseq_t2t) passing the motus filters
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
    "bbmap_repair_pe_order": parse_bbmap_repair_pe_order,
    "falco_raw": parse_reads_qc,
    "falco_processed": parse_reads_qc,
    "fastqc_raw": parse_reads_qc,
    "fastqc_processed": parse_reads_qc,
    "pathseq_flagstat": parse_pathseq_flagstat,
    "pathseq_filter": parse_pathseq_filter,
    "pathseq_t2tfilter": parse_pathseq_t2tfilter,
    "cutadapt_primers": parse_cutadapt,
    "cutadapt_readthrough": parse_cutadapt,
    "cutadapt_linked": parse_cutadapt_linked,
    "vsearch_fastqfilter": parse_vsearch_fastqfilter,
    "bbduk_pacbio": parse_bbduk_pacbio,
    "vsearch_orient": parse_vsearch_orient,
    "vsearch": parse_vsearch,
    "fastp_merge": parse_fastp_merge,
    "motus4": parse_motus4,
    "motus3": parse_motus3,
    "cayman": parse_gene_catalogue,
    "gqlib": parse_gene_catalogue,
    "idtaxa": parse_idtaxa,
    "kraken2": parse_kraken2,
    "kraken2_pairs": parse_kraken2,
}

# column, or, for the profilers, kind -> (title, description), {db} being the database
DESCRIPTIONS = {
    "primer_trimmed_reads": (
        "Primer-trimmed reads",
        "Amplicon reads with a primer (found by PIMENTO or given) trimmed at their start (cutadapt), each mate counted;"
        " long reads: with both primers trimmed, and in the length window",
    ),
    "readthrough_trimmed_reads": (
        "Read-through-trimmed reads",
        "Amplicon reads with the opposite primer trimmed at their end (cutadapt), each mate counted",
    ),
    "raw_reads": ("Raw reads", "Summary: the reads of the samplesheet (Falco/FastQC), each mate counted"),
    "qc_reads": ("QC reads", "Summary: the reads after the QC (fastp, or vsearch and BBDuk), each mate counted"),
    "hostremoval_reads": (
        "Host removal reads",
        "Summary: the reads after the host removal (bowtie2, or PathSeq-T2T), each mate counted",
    ),
    "profiler_input_reads": (
        "Profiler input reads",
        "Summary: the reads given to the profilers (Falco/FastQC), each mate counted; amplicons: the sequences after"
        " the pair handling",
    ),
    "falco_raw_reads": ("Falco raw", "Reads of the samplesheet (Falco, Total Sequences), each mate counted"),
    "fastqc_raw_reads": ("FastQC raw", "Reads of the samplesheet (FastQC, Total Sequences), each mate counted"),
    "falco_processed_reads": (
        "Falco processed",
        "Reads given to the profilers (Falco, Total Sequences), each mate counted, with the singletons of PathSeq-T2T",
    ),
    "fastqc_processed_reads": (
        "FastQC processed",
        "Reads given to the profilers (FastQC, Total Sequences), each mate counted, with the singletons of PathSeq-T2T",
    ),
    "repair_input_reads": ("Pair repair input", "Paired-end reads given to repair.sh, each mate counted"),
    "repair_paired_reads": (
        "Pair repair",
        "Paired-end reads paired by repair.sh (the reads without their mate are dropped), each mate counted",
    ),
    "shortread_qc_input_reads": ("Short-read QC input", "Short reads given to fastp, each mate counted"),
    "shortread_qc_reads": ("Short-read QC", "Short reads after fastp, each mate counted"),
    "longread_qc_input_reads": ("Long-read QC input", "PacBio reads given to the quality filter of vsearch"),
    "longread_qc_reads": ("Long-read QC", "PacBio reads after the quality filters (vsearch, then BBDuk)"),
    "oriented_reads": ("Oriented reads", "Long amplicon reads oriented as the reference sequences (vsearch), the others dropped"),
    "hostremoval_retained_reads": (
        "Host removal retained reads",
        "Reads of the pairs with neither mate aligned to the host genome (single-end: reads not aligned), each mate"
        " counted",
    ),
    "hostremoval_retained_pairs": (
        "Host removal retained pairs",
        "Read pairs with neither mate aligned to the host genome (single-end: reads not aligned), given to the"
        " profilers",
    ),
    "pathseq_t2t_input_reads": (
        "PathSeq-T2T input",
        "PathSeq-T2T: the reads given to it, each mate counted (the primary reads of the alignment to GRCh38, aligned or"
        " not, samtools flagstat)",
    ),
    "grch38_unaligned_reads": (
        "Not aligned to GRCh38",
        "PathSeq-T2T: reads not aligned to GRCh38 in a proper pair with an alignment score above 35 (single-end:"
        " not aligned with such a score), and reads aligned to its viral decoy contigs, each mate counted",
    ),
    "prealigned_host_filter_reads": (
        "Prealigned host filter",
        "PathSeq-T2T: reads left by PathSeqFilterSpark after removing those aligned to GRCh38, each mate counted",
    ),
    "quality_complexity_filter_reads": (
        "Quality/complexity filter",
        "PathSeq-T2T: reads left by PathSeqFilterSpark after the quality and complexity filters, each mate counted",
    ),
    "host_kmer_filter_reads": (
        "Host k-mer filter",
        "PathSeq-T2T: reads left by PathSeqFilterSpark after removing those with host k-mers or aligned to the host,"
        " each mate counted",
    ),
    "deduplication_filter_reads": (
        "Deduplication filter",
        "PathSeq-T2T: reads left by PathSeqFilterSpark after removing the duplicates, each mate counted",
    ),
    "t2t_filter_reads": (
        "T2T-CHM13 filter",
        "PathSeq-T2T: reads not aligned to T2T-CHM13, given to the profilers, each mate counted: the pairs, and the"
        " reads left without their mate",
    ),
    "motus4": (
        "mOTUs v4 {db} pairs",
        "Inserts (read pairs, and the reads left without their mate by PathSeq-T2T) aligned to the mOTUs v4 database {db}"
        " and passing the mOTUs filters",
    ),
    "motus3": (
        "mOTUs v3 {db} pairs",
        "Inserts (read pairs, and the reads left without their mate by PathSeq-T2T) aligned to the mOTUs v3 database {db}"
        " and passing the mOTUs filters",
    ),
    "cayman": ("cayman {db} pairs", "Read pairs aligned to the gene catalogue {db} and passing the cayman filters"),
    "gqlib": ("gqlib {db} pairs", "Read pairs aligned to the gene catalogue {db} and passing the gqlib filters"),
    "kraken2": ("Kraken2 {db} reads", "Reads classified by Kraken2 with the database {db}, each mate counted"),
    "kraken2_pairs": ("Kraken2 {db} reads", "Reads classified by Kraken2 with the database {db}, each mate counted"),
    "merged_pairs": ("Merged pairs", "Amplicon read pairs merged into one sequence (vsearch or fastp)"),
    "joined_pairs": (
        "Joined pairs",
        "Amplicon read pairs not merged, joined with a padding of Ns into one sequence (vsearch)",
    ),
    "idtaxa": (
        "IDTAXA {db} sequences",
        "Amplicon sequences (merged or joined read pairs, mates of unmerged pairs, or single-end reads) assigned to a"
        " taxon by IDTAXA with the training set {db}",
    ),
}


def read_log(path):
    # the report of FastQC is in its zip (<name>_fastqc/fastqc_data.txt)
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            name = next(name for name in archive.namelist() if name.endswith("/fastqc_data.txt"))
            return archive.read(name).decode()
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

    # the summary columns, always written: the sum of the columns of their steps, NA if none of them was run
    for summary, columns in SUMMARY_COLUMNS.items():
        present = [c for c in columns if c in counts.columns]
        counts[summary] = counts[present].sum(axis=1, min_count=1) if present else pd.NA

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
    # the columns of the logs, and raw_reads, computed from them
    columns = entries.drop_duplicates("column")[["column", "kind", "db_name"]].itertuples(index=False)
    for column, kind, db_name in [*columns, *((summary, None, "") for summary in SUMMARY_COLUMNS)]:
        title, description = DESCRIPTIONS.get(column) or DESCRIPTIONS[kind]
        headers[column] = {
            "title": title.format(db=db_name),
            "description": description.format(db=db_name),
            "format": "{:,.0f}",
            "scale": False,
        }

    # MultiQC custom content: the config as a commented header, json being valid yaml
    config = {
        "id": "read_counts",
        "section_name": "Read counts",
        "description": "Reads after each step. The columns ending in _pairs count read pairs (single-end: reads), the"
        " others count each mate.",
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
