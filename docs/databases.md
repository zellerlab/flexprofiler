# zellerlab/flexprofiler: Databases

Flexprofiler expects a csv file with the relevant database paths.
Refer to the [usage documentation](usage.md) for details on how to set this up.
These instruction explain how databases for the different profiling tools can be generated.
The relevant tools are assumed to be installed and in your path for the following instructions.
The script [`utils/fetch_databases.sh`](../utils/fetch_databases.sh) runs them and writes the corresponding database sheet: each database is downloaded only when asked, with `FETCH_MOTUS4=1`, `FETCH_MOTUS3=1`, `FETCH_CAYMAN=1`, `FETCH_GQLIB=1`, `FETCH_PATHSEQ_T2T=1`, `FETCH_PATHSEQ_T2T_STAR=1` (the STAR index of GRCh38.d1.vd1 with the GENCODE v36 genes, for `--pathseq_t2t_grch38_aligner star`) and `FETCH_KRAKEN2=1`.

## mOTUs v4

```
# ensure that you installed the correct mOTUs version before proceeding!
# mOTUs databases are not exchangeable across versions

motus downloadMGDB
DB_PATH=$(dirname $(python -c 'import motus.mutils as m; print(m.__file__)'))/db_mOTU
cp -r $DB_PATH .
```

## mOTUs v3

```
# ensure that you installed the correct mOTUs version before proceeding!
# mOTUs databases are not exchangeable across versions

motus downloadDB
DB_PATH=$(python -c 'import motus.motus as m; print(m.DATABASE)')
cp -r ${DB_PATH} .
```

## Kraken2

PathSeq-T2T (`--pathseq_t2t`, [Dohlman et al. 2026](https://doi.org/10.1016/j.cell.2026.04.015)) used the PlusPF database of 2024-06-05 (bacteria, archaea, viruses, protozoa, fungi, human and UniVec), from the [Kraken2 index collection](https://benlangmead.github.io/aws-indexes/k2).
`utils/fetch_databases.sh` downloads the latest PlusPF release, of 2024-12-28 (replace the date with 20240605 to use the database of the paper, on which its `--confidence 0.15` was chosen):

```
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20241228.tar.gz
mkdir k2_pluspf_20241228
tar -xzf k2_pluspf_20241228.tar.gz -C k2_pluspf_20241228
```

Database sheet row: `kraken2,k2_pluspf_20241228,,<path>/k2_pluspf_20241228,,`

## Cayman

```
wget https://zenodo.org/records/10473258/files/gene_catalogue_annotations.zip
wget https://zenodo.org/records/10473258/files/gene_catalogues.zip
unzip gene_catalogue_annotations.zip
unzip gene_catalogues.zip

# change accordingly for other databases among the downloaded ones
CAYMAN_DB=GMGC10.human-gut.95nr.0.5.percent.prevalence
mkdir $CAYMAN_DB
mv gene_catalogues/$CAYMAN_DB* $CAYMAN_DB
mv annots/$CAYMAN_DB* $CAYMAN_DB
bwa index $CAYMAN_DB/${CAYMAN_DB}.fna.gz
```

:::note
The pipeline does not run `cayman profile`, which aligns the reads itself: the reads are aligned with `bwa mem` (shared with gqlib), and the profiling library used by cayman, gqlib, is run directly on the alignment, as `cayman profile` does. Only the cayman annotations in the `.csv` format of the release above are supported: the `.bed` annotations, which `cayman profile` also accepts, are not.
:::

## gqlib

gqlib uses the same gene catalogue as cayman, with its bwa index (see above), and an annotation database in sqlite3 format, built with gqlib from the [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) annotation of the genes of the catalogue.

:::warning
These instructions have not been tested.
:::

The database used in the Zeller lab (`GMGC10.human-gut.95nr.0.5pv.20230422.2101.sqlite3`) was built from the annotation with eggNOG-mapper 2.1.7 in diamond mode of `GMGC10.human-gut.95nr.protein.0.5.percent.prevalence.faa`, the proteins of the GMGC10 human gut gene catalogue (95% non-redundant, genes present in at least 0.5% of the samples), the same catalogue whose genes the reads are aligned to (see [the test data documentation](../assets/test-data/README.md#gqlib-annotation-database) for the details).
Here eggNOG-mapper translates the genes of the catalogue itself (`--itype CDS`), which gives the same proteins (not tested).
The eggNOG database is version 5.0.2 (inferred from the eggNOG-mapper version, 2.1.7, which downloads and expects it).
The annotation of the full catalogue takes long (days with a few cpus) and needs ~50GB of memory with `--dbmem`.

```
# eggNOG-mapper 2.1.7 and its database (eggNOG 5.0.2)
mkdir eggnog_db
download_eggnog_data.py -y --data_dir eggnog_db

# annotation of the gene catalogue (see the cayman section above for its download)
CAYMAN_DB=GMGC10.human-gut.95nr.0.5.percent.prevalence
gzip -cd $CAYMAN_DB/${CAYMAN_DB}.fna.gz > ${CAYMAN_DB}.fna
mkdir emapper_tmp
emapper.py \
    -i ${CAYMAN_DB}.fna \
    --itype CDS \
    --data_dir eggnog_db \
    --output emapper_output \
    -m diamond \
    --cpu 8 \
    --dbmem \
    --temp_dir emapper_tmp

# the gqlib database (gqlib 2.14.4); v2.1.2 is the parser for the output of eggNOG-mapper >= 2.1.2
python -m gqlib.bin.build_gene_database \
    $CAYMAN_DB/GMGC10.human-gut.95nr.0.5pv.emapper_v2.1.7.sqlite3 \
    emapper_output.emapper.annotations \
    --initialise_db \
    --emapper_version v2.1.2
```

In the database sheet, `db_path` is the gene catalogue fasta (as for cayman) and `annot_db_path` the sqlite3 database.
