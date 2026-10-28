#!/usr/bin/env Rscript
# Classify sequences with IDTAXA (DECIPHER::IdTaxa) and count them per lineage.
#
# usage: idtaxa_classify.R <reads.fastq.gz> <training_set.rds|.RData> <prefix> [--threshold N] [--strand S] [--processors N]
#
# The sequences are classified as they are: the pairs of paired-end reads are given merged, and the others joined
# (VSEARCH_MERGEJOIN) or as their two mates (FASTP_IDTAXA_MERGE), see --idtaxa_pair_handling. The options are those of IdTaxa, with its defaults (threshold 60, strand both). The
# training set must have ranks, as the ones distributed with DECIPHER.
#
# Writes:
# - <prefix>.idtaxa.classifications.tsv.gz: one row per sequence, with its lineage (the taxa assigned above the
#   threshold, Root excluded, ";"-separated), their ranks and their confidences; an empty lineage if no taxon was
#   assigned
# - <prefix>.idtaxa.profile.tsv: the number of sequences of each lineage (lineage, ranks, count), with the rows
#   unassigned (sequences with an empty lineage) and total_reads (all the sequences)
# - stdout: the number of sequences and of assigned sequences, for the read counts

suppressPackageStartupMessages(library(DECIPHER))

args <- commandArgs(trailingOnly = TRUE)
usage <- "usage: idtaxa_classify.R <reads.fastq.gz> <training_set.rds|.RData> <prefix> [--threshold N] [--strand S] [--processors N]"
if (length(args) < 3 || length(args) %% 2 != 1) stop(usage)
options <- list(threshold = "60", strand = "both", processors = "1")
for (i in seq_len((length(args) - 3) / 2) * 2 + 2) {
    name <- sub("^--", "", args[i])
    if (!name %in% names(options)) stop(paste("unknown option", args[i], "\n", usage))
    options[[name]] <- args[i + 1]
}
fastq <- args[1]
prefix <- args[3]

# the training sets distributed with DECIPHER are .RData files with a trainingSet object, saveRDS gives .rds files
training_set_path <- args[2]
if (grepl("\\.rds$", training_set_path, ignore.case = TRUE)) {
    training_set <- readRDS(training_set_path)
} else {
    env <- new.env()
    load(training_set_path, envir = env)
    is_taxa <- vapply(ls(env), function(x) inherits(get(x, envir = env), "Taxa"), logical(1))
    if (sum(is_taxa) != 1) stop(paste("expected one training set (class Taxa) in", training_set_path))
    training_set <- get(ls(env)[is_taxa], envir = env)
}
if (is.null(training_set$ranks)) stop(paste("the training set", training_set_path, "has no ranks"))

seqs <- readDNAStringSet(fastq, format = "fastq")
names(seqs) <- sub(" .*", "", names(seqs))
ids <- IdTaxa(
    seqs,
    training_set,
    strand = options$strand,
    threshold = as.numeric(options$threshold),
    processors = as.integer(options$processors),
    verbose = FALSE
)

# the assigned taxa, without Root and the unclassified_ placeholders of the levels below the threshold
assigned <- lapply(ids, function(x) {
    keep <- x$rank != "rootrank" & x$taxon != "Root" & !startsWith(x$taxon, "unclassified_")
    list(taxon = x$taxon[keep], rank = x$rank[keep], confidence = x$confidence[keep])
})
lineages <- vapply(assigned, function(x) paste(x$taxon, collapse = ";"), character(1))
ranks <- vapply(assigned, function(x) paste(x$rank, collapse = ";"), character(1))
confidences <- vapply(assigned, function(x) paste(sprintf("%.1f", x$confidence), collapse = ";"), character(1))

classifications <- data.frame(read = names(seqs), lineage = lineages, ranks = ranks, confidence = confidences)
handle <- gzfile(paste0(prefix, ".idtaxa.classifications.tsv.gz"), "w")
write.table(classifications, handle, sep = "\t", quote = FALSE, row.names = FALSE)
close(handle)

is_assigned <- lineages != ""
profile <- aggregate(
    list(count = rep(1L, sum(is_assigned))),
    by = list(lineage = lineages[is_assigned], ranks = ranks[is_assigned]),
    FUN = sum
)
profile <- profile[order(profile$lineage), , drop = FALSE]
profile <- rbind(
    profile,
    data.frame(lineage = c("unassigned", "total_reads"), ranks = "", count = c(sum(!is_assigned), length(seqs)))
)
write.table(profile, paste0(prefix, ".idtaxa.profile.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Sequences: %d\nAssigned sequences: %d\n", length(seqs), sum(is_assigned)))
