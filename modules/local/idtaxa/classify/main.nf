process IDTAXA_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-decipher:3.6.0--r45h01b2380_0'
        : 'biocontainers/bioconductor-decipher:3.6.0--r45h01b2380_0'}"

    input:
    // one sequence per read pair (merged or joined) or single-end read
    tuple val(meta), path(reads)
    path training_set

    output:
    // the number of sequences of each lineage (the taxa assigned above the threshold, ";"-separated), with its ranks, and
    // the rows unassigned and total_reads
    tuple val(meta), path("*.idtaxa.profile.tsv"), emit: profile
    // the number of sequences and of assigned sequences, for the read counts
    tuple val(meta), path("*.idtaxa.log")        , emit: log
    tuple val("${task.process}"), val('decipher'), eval('Rscript -e "cat(as.character(packageVersion(\\"DECIPHER\\")))"'), emit: versions_decipher, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the IdTaxa options, with its defaults
    def threshold = task.ext.threshold ?: 60
    def strand = task.ext.strand ?: 'both'
    // R code: lists are indexed with [[ ]], as $ would be taken by Nextflow
    """
    #!/usr/bin/env Rscript
    library(DECIPHER)

    # the training set, the only object of the .RData file, whatever its name (trainingSet in those of DECIPHER)
    training_set <- get(load("${training_set}"))
    seqs <- readDNAStringSet("${reads}", format = "fastq")
    # IdTaxa needs at least one sequence, and sapply below would return an empty list for none
    stopifnot("no sequences in ${reads}" = length(seqs) > 0)
    # the progress bar of IdTaxa (verbose) is unreadable in a log file: the setTxtProgressBar that DECIPHER imports is
    # replaced by one that draws nothing and logs the progress (with the time) at each whole percent instead. It is called
    # at (almost) each step i of the loop of IdTaxa (DECIPHER 3.6.0) over I, read from its frame: the unique sequences (it
    # dereplicates them), then the reverse complements of those classified on both strands; value is i/length(I).
    # The next percent to log is kept in an environment, as it is modified by reference across the calls
    imports <- parent.env(asNamespace("DECIPHER"))
    progress <- new.env()
    progress[["percent"]] <- 0
    unlockBinding("setTxtProgressBar", imports)
    assign("setTxtProgressBar", function(pb, value, ...) {
        if (100 * value >= progress[["percent"]]) {
            progress[["percent"]] <- floor(100 * value) + 1
            message(
                format(Sys.time(), "%F %T"), " IdTaxa: ", floor(100 * value), "% done (", parent.frame()[["i"]], "/",
                length(parent.frame()[["I"]]), " unique sequences and reverse complements)"
            )
        }
    }, envir = imports)
    ids <- IdTaxa(seqs, training_set, strand = "${strand}", threshold = ${threshold}, processors = ${task.cpus}, verbose = TRUE)

    # the lineage of each sequence and its ranks (";"-separated, tab-separated from each other), without the root and the
    # unclassified_ placeholder of the rank where the assignment stopped; NA if no taxon was assigned
    lineages <- sapply(ids, function(x) {
        keep <- x[["rank"]] != "rootrank" & !startsWith(x[["taxon"]], "unclassified_")
        if (!any(keep)) return(NA_character_)
        paste(paste(x[["taxon"]][keep], collapse = ";"), paste(x[["rank"]][keep], collapse = ";"), sep = "\\t")
    })
    stopifnot(
        "the lineages are not a character vector" = is.character(lineages),
        "not one lineage (or NA) per sequence" = length(lineages) == length(seqs)
    )

    # the number of sequences of each lineage (table leaves out the NA), of the unassigned sequences and of all of them
    counts <- table(lineages)
    writeLines(c(
        "lineage\\tranks\\tcount",
        paste(names(counts), counts, sep = "\\t"),
        paste0("unassigned\\t\\t", sum(is.na(lineages))),
        paste0("total_reads\\t\\t", length(lineages))
    ), "${prefix}.idtaxa.profile.tsv")

    writeLines(c(paste("Sequences:", length(lineages)), paste("Assigned sequences:", sum(!is.na(lineages)))), "${prefix}.idtaxa.log")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idtaxa.profile.tsv ${prefix}.idtaxa.log
    """
}
