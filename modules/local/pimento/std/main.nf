process PIMENTO_STD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/cc161c6256c24520e3d2cb0e16571b46a1d79116bd296b74e5747de61931bb51/data'
        : 'community.wave.seqera.io/library/mi-pimento:1.3.3--db3751b71e154452'}"

    input:
    // one mate of paired-end reads, or single-end reads (searched with --merged, set in ext.args), and FASTA files of
    // primers to add to the PIMENTO standard library (one region per file, named as the file), [] for none
    tuple val(meta), path(reads), path(extra_primers, stageAs: 'extra_primers/*')

    output:
    // empty if no primer was found
    tuple val(meta), path("*_std_primers.fasta"), emit: primers
    tuple val(meta), path("*_std_primer_out.txt"), emit: out
    tuple val("${task.process}"), val('pimento'), eval("pimento --version 2>&1 | sed -n 's/.*version //p'"), emit: versions_pimento, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the standard library of PIMENTO, copied from its installation, with the extra primers
    def library_setup = extra_primers ? """
    mkdir primer_library
    cp "\$(python -c 'import os, pimento.standard_primers as library; print(os.path.dirname(library.__file__))')"/*.fasta primer_library/
    cp extra_primers/* primer_library/
    """ : ''
    def library = extra_primers ? '-p primer_library' : ''
    """
    ${library_setup}
    pimento std -i ${reads} -o ${prefix} ${library} ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_std_primers.fasta ${prefix}_std_primer_out.txt
    """
}
