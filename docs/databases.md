# zellerlab/flexprofiler: Databases

Flexprofiler expects a csv file with the relevant database paths.
Refer to the [usage documentation](usage.md) for details on how to set this up.
These instruction explain how databases for the different profiling tools can be generated.
The relevant tools are assumed to be installed and in your path for the following instructions.

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
DB_PATH=$(dirname $(python -c 'import motus.mutils as m; print(m.__file__)'))/db_mOTU
cp -r $DB_PATH .
```

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