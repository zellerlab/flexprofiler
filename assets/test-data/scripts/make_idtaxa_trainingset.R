#!/usr/bin/env Rscript
# Train the IDtaxa training set of the tests on the 16S reference sequences fetched by fetch_16s_reference.py.
#
# usage: make_idtaxa_trainingset.R <reference_16s.fasta.gz> <trainingset.RData>
#
# The reference headers are `<accession> <lineage>`, the lineage being `Root;d__...;...;g__...`. The training set has
# the ranks (rootrank, domain, phylum, class, order, family, genus) from the lineage prefixes, as the pre-trained
# DECIPHER training sets, through the taxid table of LearnTaxa.

suppressPackageStartupMessages(library(DECIPHER))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: make_idtaxa_trainingset.R <reference_16s.fasta.gz> <trainingset.RData>")

seqs <- readDNAStringSet(args[1])
lineages <- sub("^\\S+ ", "", names(seqs))
names(seqs) <- sub(" .*", "", names(seqs))

prefix_ranks <- c(d = "domain", p = "phylum", c = "class", o = "order", f = "family", g = "genus")

# taxid table: one row per taxon, its parent being the row of the lineage up to the rank above
taxa <- unique(unlist(lapply(strsplit(lineages, ";"), function(parts) {
    vapply(seq_along(parts), function(i) paste(parts[seq_len(i)], collapse = ";"), character(1))
})))
taxa <- taxa[order(lengths(strsplit(taxa, ";")), taxa)]
names_ <- vapply(strsplit(taxa, ";"), function(parts) parts[length(parts)], character(1))
parents <- vapply(strsplit(taxa, ";"), function(parts) {
    if (length(parts) == 1) "" else paste(parts[-length(parts)], collapse = ";")
}, character(1))
taxid <- data.frame(
    Index = seq_along(taxa) - 1L,
    Name = names_,
    Parent = ifelse(parents == "", 0L, match(parents, taxa) - 1L),
    Level = lengths(strsplit(taxa, ";")) - 1L,
    Rank = ifelse(names_ == "Root", "rootrank", prefix_ranks[substr(names_, 1, 1)]),
    stringsAsFactors = FALSE
)

set.seed(42)
# saved as the pre-trained DECIPHER training sets, an .RData file with a trainingSet object
trainingSet <- LearnTaxa(seqs, paste0(lineages, ";"), taxid, verbose = FALSE)
save(trainingSet, file = args[2])
print(trainingSet)
