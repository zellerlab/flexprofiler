## Flexprofiler test specific information

> This test data and documentation have been adapted from the [nf-core/taxprofiler test data](https://github.com/nf-core/test-datasets/tree/taxprofiler)

### FASTQ

The main test data used for zellerlab/flexprofiler is from [Maixner et al. (2021) _Curr. Bio._](https://doi.org/10.1016/j.cub.2021.09.031), with ENA project accession ID: PRJEB44507.
The following selected libraries were all sequenced on an Illumina MiSeq, and were selected due to their small size (~1million reads, <100MB) and known mixture of (gut) bacteria, (ancient human) eukaryotes, and (yeast) fungi (according to the results of the paper).

- ERX5474937
- ERX5474932
- ERX5474930
- ERX5474936

Data was downloaded with nf-core/fetchNGS 1.5 (with Nextflow 21.10.06):

```bash
nextflow run nf-core/fetchngs --input maixner2021_acc_codes.txt --input_type sra
```

FASTQ files are stored under `data/fastq/`

Test data for long reads with ENA project accession ID: PRJEB29152. They were subsampled with seqtk 1.3-r106:

```bash
seqtk sample ERR3201952.fastq.gz 10000 > ERR3201952.fastq.gz
```

### FASTA

One of the files was converted to FASTA file with seqtk 1.3-r106

```bash
seqtk seq -a  ERX5474930_ERR5766174_1.fastq.gz > ERX5474930_ERR5766174_1.fa.gz
```

FASTA files are stored under `data/fasta/`