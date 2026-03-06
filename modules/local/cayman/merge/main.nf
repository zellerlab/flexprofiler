process CAYMAN_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/saulpierotti/r_datatable_tidyverse"

    input:
    tuple val(meta), path(input)
    path db

    output:
    tuple val(meta), path("*.txt.gz") , optional: true, emit: txt
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_rvec = input.collect { "\"${it}\"" }.join(', ')
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    
    library("data.table")

    inputs <- c(${input_rvec})

    read_f <- function(f) {
        df <- fread(f)
        df[, sample := gsub(".txt.gz", "", basename(f))]
        return(df)
    }

    df <- rbindlist(lapply(inputs, read_f))
    fwrite(df, file="${prefix}_cayman_merged.txt.gz", sep="\t")

    writeLines(sprintf("data.table: %s", packageVersion('data.table')), "versions.yml")
    writeLines(sprintf("r-base: %s", getRversion()), "versions.yml", append=TRUE)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cayman_merged.txt.gz

    Rscript -e "writeLines(sprintf("data.table: %s", packageVersion('data.table')), "versions.yml")
    Rscript -e "writeLines(sprintf("r-base: %s", getRversion()), "versions.yml", append=TRUE)"
    """

}
