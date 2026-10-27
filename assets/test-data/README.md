## Flexprofiler test specific information

> This test data and documentation have been adapted from the [nf-core/taxprofiler test data](https://github.com/nf-core/test-datasets/tree/taxprofiler)

The samplesheets (`samplesheet*.csv`) and the database sheet (`databases.csv`, with the databases of all the profilers) used by the `test_*` profiles give the paths relative to the directory where nf-test launches the pipeline (`.nf-test/tests/<hash>/`), so they work only in the nf-test tests.

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

### FASTA

One of the files was converted to FASTA file with seqtk 1.3-r106

```bash
seqtk seq -a  ERX5474930_ERR5766174_1.fastq.gz > ERX5474930_ERR5766174_1.fa.gz
```

FASTA files are stored under `data/fasta/`

### Host removal reference

The tests use the human mitochondrial genome (revised Cambridge Reference Sequence, [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1)) as host reference instead of a full human genome, so that building its bowtie2 index (done by the pipeline) and the host removal are fast. The test data are not expected to have many human reads, and the tests only check that host removal runs.

The T2T-CHM13v2.0 RefSeq assembly recommended for real runs (GCF_009914755.1) does not include the mitochondrial genome, so it was downloaded from NCBI:

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" \
    | gzip -n > data/hostremoval/NC_012920.1.fna.gz
```

The reference is stored under `data/hostremoval/`.

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

The fasta and its index are stored under `data/gmgc10/`. All the reads that would align elsewhere in the full catalogue are either unaligned or aligned to these two genes, so the profiles of the test data differ from those with the full catalogue.

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

The annotation table was generated (13 Oct 2022) with eggNOG-mapper 2.1.7 in diamond mode, from the protein sequences of the gene catalogue (`GMGC10.human-gut.95nr.protein.0.5.percent.prevalence.faa`), as recorded in its header:

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
- The table is read twice: the first pass encodes each category and feature as an integer (`category` and `feature` tables, also dumped to `<db>.code_map.json.gz`), the second stores the encoded annotations of each gene (`annotatedsequence` table, one row per annotated gene, `annotation_str` as `<category_id>=<feature_id>,<feature_id>;...`). Genes without annotations are not stored.
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
