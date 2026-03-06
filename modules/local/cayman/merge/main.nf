process CAYMAN_MERGE {
    tag "$meta.id"
    label 'process_single'

    // this should be set to something more stable in time
    conda "${moduleDir}/environment.yml"
    container "quay.io/saulpierotti/r_datatable_tidyverse:1.0.0"

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
        df[, sample := gsub("_${meta.id}.cazy.txt.gz", "", basename(f))]
        return(df)
    }

    df <- rbindlist(lapply(inputs, read_f))
    fwrite(df, file="${prefix}_cayman_merged.txt.gz", sep="\t")

    handle <- file("versions.yml", "a")
    writeLines("\\"${task.process}\\":", handle)
    writeLines(sprintf("   data.table: %s", packageVersion('data.table')), handle)
    writeLines(sprintf("   r-base: %s", getRversion()), handle)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cayman_merged.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        data.table: \$(Rscript -e "sprintf("%s", packageVersion('data.table'))")
        r-base: \$(Rscript -e "sprintf("%s", getRversion())")
    END_VERSIONS
    """

}
