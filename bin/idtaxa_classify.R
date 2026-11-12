#!/usr/bin/env Rscript
# Classify amplicon sequences with IDTAXA (DECIPHER) and count the sequences of each lineage.
#
# usage: idtaxa_classify.R <training_set.RData> <reads.fastq[.gz]> <prefix> <strand> <threshold> <processors>
#
# Writes <prefix>.idtaxa.profile.tsv, the number of sequences of each lineage (the taxa assigned above the threshold,
# ";"-separated), with its ranks, and the rows unassigned and total_reads, and <prefix>.idtaxa.log, the number of
# sequences and of assigned sequences, for the read counts.

args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: idtaxa_classify.R <training_set> <reads> <prefix> <strand> <threshold> <processors>" = length(args) == 6)
training_set_file <- args[1]
reads <- args[2]
prefix <- args[3]
strand <- args[4]
threshold <- as.numeric(args[5])
processors <- as.integer(args[6])

library(DECIPHER)

# the training set, the only object of the .RData file, whatever its name (trainingSet in those of DECIPHER)
training_set <- get(load(training_set_file))
seqs <- readDNAStringSet(reads, format = "fastq")
# IdTaxa needs at least one sequence, and sapply below would return an empty list for none
if (length(seqs) == 0) stop("no sequences in ", reads)
# the progress bar of IdTaxa (verbose) is unreadable in a log file: the setTxtProgressBar that DECIPHER imports is
# replaced by one that draws nothing and logs the progress (with the time) at each whole percent instead. It is called
# at (almost) each step i of the loop of IdTaxa (DECIPHER 3.6.0) over I, read from its frame: the unique sequences (it
# dereplicates them), then the reverse complements of those classified on both strands; value is i/length(I).
# The next percent to log is kept in an environment, as it is modified by reference across the calls
imports <- parent.env(asNamespace("DECIPHER"))
progress <- new.env()
progress$percent <- 0
unlockBinding("setTxtProgressBar", imports)
assign("setTxtProgressBar", function(pb, value, ...) {
    if (100 * value >= progress$percent) {
        progress$percent <- floor(100 * value) + 1
        message(
            format(Sys.time(), "%F %T"), " IdTaxa: ", floor(100 * value), "% done (", parent.frame()$i, "/",
            length(parent.frame()$I), " unique sequences and reverse complements)"
        )
    }
}, envir = imports)
ids <- IdTaxa(seqs, training_set, strand = strand, threshold = threshold, processors = processors, verbose = TRUE)

# the lineage of each sequence and its ranks (";"-separated, tab-separated from each other), without the root and the
# unclassified_ placeholder of the rank where the assignment stopped; NA if no taxon was assigned
lineages <- sapply(ids, function(x) {
    keep <- x$rank != "rootrank" & !startsWith(x$taxon, "unclassified_")
    if (!any(keep)) return(NA_character_)
    paste(paste(x$taxon[keep], collapse = ";"), paste(x$rank[keep], collapse = ";"), sep = "\t")
})
stopifnot(
    "the lineages are not a character vector" = is.character(lineages),
    "not one lineage (or NA) per sequence" = length(lineages) == length(seqs)
)

# the number of sequences of each lineage (table leaves out the NA), of the unassigned sequences and of all of them
counts <- table(lineages)
writeLines(c(
    "lineage\tranks\tcount",
    paste(names(counts), counts, sep = "\t"),
    paste0("unassigned\t\t", sum(is.na(lineages))),
    paste0("total_reads\t\t", length(lineages))
), paste0(prefix, ".idtaxa.profile.tsv"))

writeLines(c(paste("Sequences:", length(lineages)), paste("Assigned sequences:", sum(!is.na(lineages)))), paste0(prefix, ".idtaxa.log"))
