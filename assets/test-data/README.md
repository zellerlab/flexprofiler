## Flexprofiler test specific information

> This test data and documentation have been adapted from the [nf-core/taxprofiler test data](https://github.com/nf-core/test-datasets/tree/taxprofiler)

The samplesheets (`samplesheet*.csv`) and the database sheet (`databases.csv`, with the databases of all the profilers) used by the `test_*` profiles give the paths relative to the directory where nf-test launches the pipeline (`.nf-test/tests/<hash>/`), so they work only in the nf-test tests.

### Regenerating the test data

[`scripts/make_test_data.sh`](scripts/make_test_data.sh) runs the commands of the sections below in order, to regenerate the test data (`scripts/make_test_data.sh --help` lists its steps and the files they write):

```bash
scripts/make_test_data.sh            # the files made only from files in the repo (the `local` group)
scripts/make_test_data.sh download   # the files downloaded from NCBI (host reference, 16S reference)
scripts/make_test_data.sh \
    --motus4-db <db_mOTU v4> --motus3-db <db_mOTU v3> \
    --gmgc10-fna <catalogue.fna.gz> --cayman-csv <annotation.csv> --gqlib-db <annotation.sqlite3> \
    databases                        # the subsets of the full databases
```

The downloaded files can change over time, and with them the files made from them, so they are regenerated only when asked. The original reads and fastas from nf-core/taxprofiler and the hand-written files (`data/fastq/adapters/adapters.fasta`, `data/gmgc10/genes.txt`) are not regenerated. The reads are simulated with fixed seeds and compressed without timestamps, so regenerating them gives identical files.

### FASTQ

The main test data used for zellerlab/flexprofiler is from [Maixner et al. (2021) _Curr. Bio._](https://doi.org/10.1016/j.cub.2021.09.031), with ENA project accession ID: PRJEB44507.
The following selected libraries were all sequenced on an Illumina MiSeq, and were selected due to their small size (~1million reads, <100MB) and known mixture of (gut) bacteria, (ancient human) eukaryotes, and (yeast) fungi (according to the results of the paper).

- ERX5474937
- ERX5474932
- ERX5474930
- ERX5474936

Data was downloaded with nf-core/fetchNGS 1.5 (with Nextflow 21.10.06) (from the nf-core/taxprofiler documentation, not re-run):

```bash
nextflow run nf-core/fetchngs --input maixner2021_acc_codes.txt --input_type sra
```

FASTQ files are stored under `data/fastq/`

Test data for long reads with ENA project accession ID: PRJEB29152. They were subsampled with seqtk 1.3-r106 (from the nf-core/taxprofiler documentation, not re-run):

```bash
seqtk sample ERR3201952.fastq.gz 10000 > ERR3201952.fastq.gz
```

#### Read names with /1 and /2 suffixes

In the files above, the /1 and /2 of the mates are in the comment of the read names (`@ERR5766181.1 M00380:...:10763/1`), so the read names of the mates are identical. To test that the pipeline handles read names with the suffixes (`@ERR5766181.1/1 M00380:...:10763`), as mOTUs v4 requires identical names (`--motus4_rename_reads`), the suffix was moved to the read name in a copy of ERX5474937 (sample 2613 in the tests), `ERX5474937_ERR5766181_suffixed_{1,2}.fastq.gz`, used in `samplesheet_read_suffixes.csv`:

```bash
for i in 1 2; do
    gzip -cd ERX5474937_ERR5766181_$i.fastq.gz \
        | awk -v s="/$i" 'NR % 4 == 1 { n = split($0, f, " "); name = f[1] s; comment = substr($0, length(f[1]) + 2); sub("/[12]$", "", comment); print (comment == "" ? name : name " " comment); next } { print }' \
        | gzip -n > ERX5474937_ERR5766181_suffixed_$i.fastq.gz
done
```

#### Paired-end reads in a different order in the two files

To test the repair of the order of the pairs (`--perform_repair_pe_order`), the reads of the second file of ERX5474937 (sample 2613 in the tests) were sorted by sequence, `ERX5474937_ERR5766181_unordered_2.fastq.gz`. It is used with the original first file in `samplesheet_unordered_pairs.csv`:

```bash
gzip -cd ERX5474937_ERR5766181_2.fastq.gz \
    | paste - - - - \
    | LC_ALL=C sort -t $'\t' -k2,2 -k1,1 \
    | tr '\t' '\n' \
    | gzip -n > ERX5474937_ERR5766181_unordered_2.fastq.gz
```

`samplesheet_mixed.csv` combines in one run the reads of ERX5474937 as they are, with the suffixes and with the pairs in a different order, the last two also as two runs of one sample (`2613_mixed`), and the single-end reads of ERX5474936. `samplesheet_mixed_endedness.csv` has a sample with a paired-end run (ERX5474932) and a single-end run (ERX5474936), `2612_mixed`, which can be merged only with `--shortread_qc_mergepairs`. As the samplesheet does not allow the same file in two rows, the files used more than once are symbolic links in `data/fastq/mixed/`.

### FASTA

One of the files was converted to FASTA file with seqtk 1.3-r106 (from the nf-core/taxprofiler documentation, not re-run):

```bash
seqtk seq -a  ERX5474930_ERR5766174_1.fastq.gz > ERX5474930_ERR5766174_1.fa.gz
```

FASTA files are stored under `data/fasta/`

### Host removal reference

The tests use the human mitochondrial genome (revised Cambridge Reference Sequence, [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1)) as host reference instead of a full human genome, so that building its bowtie2 index (done by the pipeline) and the host removal are fast. The test data are not expected to have many human reads, and most tests only check that host removal runs.

#### Host-contaminated reads

To test that host removal removes the host reads and keeps the others (`tests/test_hostremoval.nf.test`, `samplesheet_hostremoval.csv`), reads simulated from the host reference were added to the first 2000 pairs of ERX5474937 (sample 2613 in the other tests), in `data/fastq/hostremoval/`: 200 pairs with both mates from the host, and 200 chimeric pairs with one mate from the host and the other from ERX5474937 (100 with the host mate as R1, 100 as R2), and, for single-end reads, 200 host reads added to the R1 of the 2000 pairs. The host-derived reads have names starting with `host_` or `chimera_`, and the reads are shuffled. They were made with:

```bash
python3 scripts/make_host_contaminated_reads.py \
    data/fastq/ERX5474937_ERR5766181_1.fastq.gz data/fastq/ERX5474937_ERR5766181_2.fastq.gz \
    data/hostremoval/NC_012920.1.fna.gz data/fastq/hostremoval
```

The T2T-CHM13v2.0 RefSeq assembly recommended for real runs (GCF_009914755.1) does not include the mitochondrial genome, so it was downloaded from NCBI:

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" \
    | gzip -n > data/hostremoval/NC_012920.1.fna.gz
```

The reference is stored under `data/hostremoval/`.

### Adapter-contaminated reads

To test that adapter removal removes the adapters and keeps the inserts (`tests/test_adapter_removal.nf.test`, `samplesheet_adapters.csv`), error-free 2x150 reads were simulated from the host reference (see above) for two library kits, TruSeq and Nextera, in `data/fastq/adapters/`. For each kit, `adapters_<kit>_{1,2}.fastq.gz` have 1000 pairs, shuffled:

- 700 pairs of 300 bp fragments, longer than the reads, so without adapter (`clean_<i>`)
- 300 pairs of fragments of 40-140 bp, shorter than the reads, which continue into the adapter of the kit and then random bases (`adapter_<insert length>_<i>`). After adapter removal, they must be as long as their insert.

`adapters_<kit>_se.fastq.gz` has the R1 of the pairs. fastp detects the adapters of single-end reads only in files of at least 10000 reads, so it detects none in these. To test the detection, `adapters_truseq_se_large.fastq.gz` (`samplesheet_adapters_se.csv`) has the R1 of 12000 TruSeq pairs made in the same way, with another seed (8400 `clean_` and 3600 `adapter_` pairs). They were made with [`scripts/make_adapter_reads.py`](scripts/make_adapter_reads.py) (Python 3 standard library only):

```bash
python3 scripts/make_adapter_reads.py data/hostremoval/NC_012920.1.fna.gz data/fastq/adapters
```

`data/fastq/adapters/adapters.fasta`, written by hand, has the TruSeq R1 and R2 adapters and the Nextera adapter, for the test with an adapter list (`--shortread_qc_adapterlist`).

### 16S amplicon reads and IDtaxa training set

The amplicon tests (`tests/test_primer_removal.nf.test`, `tests/test_idtaxa.nf.test`) use 16S V3-V4 reads simulated from a small set of 16S reference sequences, and an IDtaxa training set trained on the same sequences, so that the genus of each read is known.

#### 16S reference sequences

The reference sequences are full-length (>= 1300 bp) 16S sequences from the NCBI RefSeq targeted loci (BioProject 33175) of six gut bacteria genera from five phyla (two of them in Bacillota): Escherichia, Bacteroides, Bifidobacterium, Akkermansia, Streptococcus and Clostridium. Only the sequences with the 341F (`CCTACGGGNGGCWGCAG`) and 805R (`GACTACHVGGGTATCTAATCC`) primer sites are kept, and at most the first 15 (by accession) of each genus. The headers are `<accession> Root;d__...;p__...;c__...;o__...;f__...;g__<genus>`. They were downloaded with [`scripts/fetch_16s_reference.py`](scripts/fetch_16s_reference.py) (Python 3 standard library only):

```bash
python3 scripts/fetch_16s_reference.py data/idtaxa/reference_16s.fasta.gz
```

When downloaded (28 Oct 2026), fewer sequences passed the filters for Escherichia (7) and Akkermansia (2) than for the other genera (15). As RefSeq is updated, downloading them again can give different sequences, and then different reads and training set.

#### Amplicon reads

The reads were simulated, error-free, from the V3-V4 amplicon (341F-805R) of each reference sequence with [`scripts/make_amplicon_reads.py`](scripts/make_amplicon_reads.py) (Python 3 standard library only), with names `<kind>_<genus>_<i>` (the genus being the truth for the classification), in `data/fastq/amplicon/`:

```bash
python3 scripts/make_amplicon_reads.py data/idtaxa/reference_16s.fasta.gz data/fastq/amplicon
```

- `amplicon_clean_{1,2}.fastq.gz`: reads with the primers already removed (the insert is the region between the primers), 5 mostly overlapping 2x250 pairs (`ovl_`) and 2 non-overlapping 2x150 pairs (`noovl_`) per reference sequence. `amplicon_clean_se.fastq.gz` has their R1. Used in `samplesheet_amplicon.csv` (IDtaxa test).
- `amplicon_primers_{1,2}.fastq.gz`: 2x250 pairs of the whole amplicon, starting with 341F (R1) and 805R (R2), per reference sequence 5 starting with the primer (`pr_`), 2 after a 1-7 base spacer (`spacer_`), and 1 with a 120 bp region, reading through the end of the amplicon into the reverse complement of the other primer and the adapter (`rt_`). `amplicon_primers_se.fastq.gz` has their R1.
- `amplicon_mixed_{1,2}.fastq.gz`: as `amplicon_primers`, with a quarter of the pairs flipped (`flip_`: R1 starting with 805R and R2 with 341F), as in libraries ligated in both orientations.

The last two are used in `samplesheet_primers.csv` (primer removal test), with `amplicon_clean` as a sample without primers.

#### IDtaxa training set

The training set (`data/idtaxa/trainingset.RData`, the `trainingSet` object saved as the pre-trained DECIPHER training sets) was trained on the reference sequences with `LearnTaxa` of DECIPHER, with the ranks (rootrank, domain, phylum, class, order, family, genus) from the prefixes of the lineages, with [`scripts/make_idtaxa_trainingset.R`](scripts/make_idtaxa_trainingset.R):

```bash
Rscript scripts/make_idtaxa_trainingset.R data/idtaxa/reference_16s.fasta.gz data/idtaxa/trainingset.RData
```

The reference sequences and the training set are stored under `data/idtaxa/`.

### mOTUs v4 database

The tests use the mOTUs v4.0 database ([zenodo](https://zenodo.org/records/17668622), ~9GB with its bwa index) restricted to 450 marker genes, so that it is small (~1.5MB) and loading it is fast. It keeps the marker genes, with alignments of the test reads passing the motus filters, of four species found in the test data:

- `mOTUv4.0_014686` (Collinsella sp900768795)
- `mOTUv4.0_000924` (Methanobrevibacter_A smithii)
- `mOTUv4.0_001912` (Paratractidigestivibacter sp003150175)
- `mOTUv4.0_014998` (Absicoccus sp000434355)

and 50 of those of the unassigned mOTU (`mOTUv4.0_unassigned`), which must be in the database. The first two are checked by the tests, the others were chosen as detected in the test data with few marker genes.

The marker genes hit were obtained by aligning each short-read file of the samples profiled by motus (2612, 2612_se, 2613; raw reads, a superset of the reads left after QC and host removal) to the full database as motus does (single end, `bwa mem -a`), keeping the marker genes with alignments passing the motus filters (alignment length >= 75, the default `-l`, and identity >= 97%), with [`scripts/motus_marker_genes_hit.sh`](scripts/motus_marker_genes_hit.sh) (BWA 0.7.19-r1273, samtools 1.24):

```bash
for reads in data/fastq/ERX547493{2,6,7}_*.fastq.gz; do
    scripts/motus_marker_genes_hit.sh db_mOTU/mOTUsv4.0.db.fna.gz $(basename $reads .fastq.gz).txt $reads
done
cat *.txt | sort -u > marker_genes_hit.txt
```

Then, with `mOTUsv4.0.map.tsv.gz` (columns `#MOTU`, `MGC`, `COG`, `MG`) the marker genes of the four species and the first 50 (in sorted order) of the unassigned mOTU were selected ([`data/motus4/marker_genes.txt`](data/motus4/marker_genes.txt)):

```bash
gzip -dc db_mOTU/mOTUsv4.0.map.tsv.gz | awk -F'\t' 'NR==FNR {hit[$1]=1; next} ($4 in hit)' marker_genes_hit.txt - > map_hit.tsv
{
    awk -F'\t' '$1 ~ /^mOTUv4.0_(014686|000924|001912|014998)$/ {print $4}' map_hit.tsv
    awk -F'\t' '$1 == "mOTUv4.0_unassigned" {print $4}' map_hit.tsv | sort | head -50
} | sort -u > data/motus4/marker_genes.txt
```

The database was subset to these marker genes with [`scripts/subset_motus4_db.py`](scripts/subset_motus4_db.py) (Python 3 standard library only), which keeps the format of all its files and the rows of the kept marker genes and mOTUs, and rewrites the counts and md5s in `mOTUsv4.0.db` (motus only reads its version and date). It was then indexed with BWA 0.7.19-r1273:

```bash
python scripts/subset_motus4_db.py db_mOTU data/motus4/marker_genes.txt data/motus4/db_mOTU
bwa index data/motus4/db_mOTU/mOTUsv4.0.db.fna.gz
```

With motus 4.0.4, the four species are detected in the test data with counts close to those with the full database. The reads of the other species are either not aligned or aligned to the marker genes of the unassigned mOTU, whose counts are higher.

The database is stored under `data/motus4/db_mOTU/`.

### mOTUs v3 database

The tests use the mOTUs v3.1.0 database ([zenodo](https://zenodo.org/record/7778108), ~7GB with its bwa indices) restricted to all the marker genes of four mOTUs found in the test data, so that it is small (~7MB) and loading it is fast:

- `ref_mOTU_v31_03626` (Collinsella aerofaciens), checked by the test
- `ref_mOTU_v31_03695` (Methanobrevibacter smithii)
- `ext_mOTU_v31_17631` (Collinsella species incertae sedis)
- `ref_mOTU_v31_03379` (Halococcus morrhuae)

and the first 50 marker genes of the unassigned marker gene clusters (`unassigned` in `db_mOTU_MAP_MGCs_to_mOTUs.tsv`), which must be in the database ([`data/motus3/marker_genes.txt`](data/motus3/marker_genes.txt)):

```bash
{
    awk -F'\t' '$4 ~ /^(ref_mOTU_v31_03626|ref_mOTU_v31_03695|ext_mOTU_v31_17631|ref_mOTU_v31_03379)\./ {print $1}' db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv
    awk -F'\t' '$4 ~ /^unassigned\./ {print $1}' db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv | head -50
} > data/motus3/marker_genes.txt
```

The database was subset to these marker genes with [`scripts/subset_motus3_db.py`](scripts/subset_motus3_db.py) (Python >= 3.9 standard library only), which keeps the format of all its files and the rows of the kept marker genes, marker gene clusters and mOTUs. The CEN database (one centroid per marker gene cluster, used for SNV calling) is restricted to the kept marker gene clusters, and the test data and public profiles of the database are not included. Both fastas were then indexed with BWA 0.7.19-r1273:

```bash
python scripts/subset_motus3_db.py db_mOTU data/motus3/marker_genes.txt data/motus3/db_mOTU
bwa index data/motus3/db_mOTU/db_mOTU_DB_NR.fasta
bwa index data/motus3/db_mOTU/db_mOTU_DB_CEN.fasta
```

With motus 3.1.0, the four mOTUs are detected in the test data with counts close to those with the full database. No reads align to the 50 unassigned marker genes, so the unassigned counts are 0.

The database is stored under `data/motus3/db_mOTU/`.

### GMGC10 gene catalogue

cayman and gqlib align the reads to a gene catalogue with bwa mem, then annotate the alignments. The tests use the GMGC10 human gut gene catalogue (`GMGC10.human-gut.95nr.0.5.percent.prevalence`, from the [cayman databases](https://zenodo.org/records/10473258)), restricted to two genes so that the database is small and the alignment fast ([`data/gmgc10/genes.txt`](data/gmgc10/genes.txt)):

- `GMGC10.266_871_379.GLGB` (glycogen branching enzyme)
- `GMGC10.266_660_280.AMY` (alpha-amylase)

They were chosen among the genes with reads of the test data aligned to them in both samples 2612_se and 2613 (from the `gene_counts` profiles of a run of the cayman test with the full catalogue) as those with the most counts and annotated in the most gqlib categories (all but KEGG_TC). Both are annotated by cayman (GH13, and CBM48 or CBM34).

The sequences were extracted from the full catalogue and indexed with BWA 0.7.19-r1273:

```bash
gzip -dc GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz \
    | awk 'NR==FNR {keep[">"$1]=1; next} /^>/ {p = ($1 in keep)} p' data/gmgc10/genes.txt - \
    | gzip -n > data/gmgc10/GMGC10.human-gut.95nr.0.5.percent.prevalence.test_subset.fna.gz
bwa index data/gmgc10/GMGC10.human-gut.95nr.0.5.percent.prevalence.test_subset.fna.gz
```

The fasta and its index are stored under `data/gmgc10/`. All the reads that would align elsewhere in the full catalogue are either unaligned or aligned to these two genes (not tested), so the profiles of the test data differ from those with the full catalogue.

### cayman annotation

The cayman annotation of the catalogue (`GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv`, from the [cayman databases](https://zenodo.org/records/10473258)) was restricted to the rows of the two genes, keeping the header:

```bash
awk -F, 'NR==FNR {keep[$1]=1; next} FNR==1 || ($1 in keep)' \
    data/gmgc10/genes.txt \
    GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv \
    > data/cayman/GMGC10.human-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.test_subset.csv
```

The subset is stored under `data/cayman/`.

### gqlib annotation database

gqlib annotates the alignments to a gene catalogue with an sqlite3 database (`annot_db_path` in the database sheet). The test uses the database for the GMGC10 human gut gene catalogue, restricted to the two genes of the test catalogue (see above).

#### Building the full database

The full database (`GMGC10.human-gut.95nr.0.5pv.20230422.2101.sqlite3`, ~2GB, available on the lab shared storage) was built with gqlib from the eggNOG-mapper annotation table of the gene catalogue.

The annotation table was generated (13 Oct 2022) with eggNOG-mapper 2.1.7 in diamond mode, from the protein sequences of the gene catalogue (`GMGC10.human-gut.95nr.protein.0.5.percent.prevalence.faa`), with the eggNOG database version 5.0.2 (inferred from the eggNOG-mapper version, which downloads and expects it), as recorded in its header:

```bash
emapper.py \
    -i GMGC10.human-gut.95nr.protein.0.5.percent.prevalence.faa \
    --data_dir eggnog_db \
    --output emapper_output \
    -m diamond \
    --cpu 8 \
    --dbmem \
    --temp_dir tmp
```

The protein file (`GMGC10.human-gut.95nr.protein.0.5.percent.prevalence.faa`) is internal to the Bork group at EMBL and not publicly available. To regenerate the annotation table, run eggNOG-mapper on the nucleotide sequences of the gene catalogue (`GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz`), which it translates itself (`--itype CDS`, see the gqlib section of [the database documentation](../../docs/databases.md#gqlib)).

This is not the eggNOG-mapper annotation table available from the [GMGC downloads](https://gmgc.embl.de/download.cgi) (`GMGC10.emapper2.annotations.tsv.gz`), which was generated with eggNOG-mapper 1.0.3 and has a different format and set of columns (for example, it has no PFAMs).

The database was then built from the annotation table (`emapper_output.emapper.annotations`) with gqlib:

```bash
python -m gqlib.bin.build_gene_database \
    GMGC10.human-gut.95nr.0.5pv.20230422.2101.sqlite3 \
    emapper_output.emapper.annotations \
    --initialise_db \
    --emapper_version v2.1.2
```

- `--emapper_version v2.1.2` is the parser for the output of eggNOG-mapper >= 2.1.2 (the annotations are from eggNOG-mapper 2.1.7, `annot_db_name` `emapper_v2.1.7` in the test). It reads the columns eggNOG_OGs, COG_cat, GOs, EC, KEGG_ko, KEGG_Pathway, KEGG_Module, KEGG_Reaction, KEGG_rclass, BRITE, KEGG_TC, CAZy, BiGG_Reaction and PFAMs, which become the categories of the database and the profiles. The older formats (`v1`, `v2`) have different columns and category names.
- The table is read twice: the first pass encodes each category and feature as an integer (`category` and `feature` tables, also dumped to `<db>.code_map.json.gz`), the second stores the encoded annotations of each gene (`annotatedsequence` table, one row per annotated gene, `annotation_str` as `<category_id>=<feature_id>,<feature_id>;...`). Genes without annotations are not stored (from the gqlib 2.14.4 code, not tested).
- The database is left in WAL journal mode.

#### Building the test subset

To keep the test database small (~28KB), it is restricted to the two genes of the test gene catalogue.

The full database was subset to these genes with [`scripts/subset_gqlib_db.py`](scripts/subset_gqlib_db.py) (Python 3 standard library only):

```bash
python scripts/subset_gqlib_db.py \
    GMGC10.human-gut.95nr.0.5pv.20230422.2101.sqlite3 \
    data/gmgc10/genes.txt \
    data/gqlib/GMGC10.human-gut.95nr.0.5pv.20230422.2101.test_subset.sqlite3
```

The subset has the same schema and the same category and feature ids as the full database. It keeps the `annotatedsequence` rows of the two genes, the features they refer to (122) and all the categories. Unlike the full database, it is not in WAL mode, so that reading it does not create `-wal` and `-shm` files next to it.

The subset is stored under `data/gqlib/`.
