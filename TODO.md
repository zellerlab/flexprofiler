# TODOs
- Add support for 16s profilers: IDTaxa, MapSeq
- Add support for metagenomics profilers: kraken2, PathSeq/pythseq, sylph
- Host removal: run bowtie2 in single-end mode (`-U`) on R1 and R2 separately, and drop a pair if either mate maps to the host, as recommended in https://doi.org/10.1016/j.crmeth.2025.101218
