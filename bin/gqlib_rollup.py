#!/usr/bin/env python3
# functional roll-up of gene counts from gqlib_patch.py, reproducing the gqlib (gffquant) one:
# - the counts of each gene are added in full to each of its features (a gene annotated with
#   more features counts towards all of them), and once to the category row
# - scaled = lnorm * sum(raw) / sum(lnorm), with the sums taken over the annotated genes
# - rpkm = lnorm * 1e9 / filtered_reads
# - unannotated genes are not reported, as in gqlib gene mode
# differently from gqlib, all the features of a gene are used also when the annotation has
# one row per feature (e.g. one row per domain), and not only the first one
# the annotation is a table with a header, with a gene id column and one column per category
# (multiple features in the same cell separated by --feature_sep), and the genes can be on more
# rows; one output table is written for each category, as <out_prefix>.<category>.txt.gz

import argparse
import gzip
import sys

METRICS = ('raw', 'lnorm', 'scaled', 'rpkm')


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gene_counts', required=True, help='gene_counts table from gqlib_patch.py')
    ap.add_argument('--read_counts', required=True, help='read_counts table from gqlib_patch.py')
    ap.add_argument('--annotation', required=True, help='table with one or more rows per gene (can be gzipped)')
    ap.add_argument('--gene_col', default=None, help='column of the annotation with the gene ids [default: first column]')
    ap.add_argument('--feature_cols', default=None, help='comma-separated columns of the annotation to use as categories [default: all except --gene_col]')
    ap.add_argument('--sep', default=None, help='annotation column separator [default: , for .csv(.gz) files, tab otherwise]')
    ap.add_argument('--feature_sep', default=',', help='separator for multiple features in the same cell [default: ,]')
    ap.add_argument('--restrict_metrics', default='raw,lnorm,scaled', help=f'comma-separated subset of {",".join(METRICS)}')
    ap.add_argument('--out_prefix', required=True)
    args = ap.parse_args()
    args.restrict_metrics = args.restrict_metrics.split(',')
    if not set(args.restrict_metrics) <= set(METRICS):
        ap.error(f'--restrict_metrics must be a subset of {",".join(METRICS)}')
    if args.sep is None:
        args.sep = ',' if args.annotation.removesuffix('.gz').endswith('.csv') else '\t'
    return args


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_gene_counts(path):
    # returns {gene: [uniq_raw, uniq_lnorm, combined_raw, combined_lnorm]}, has_ambig_counts
    with open_text(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        missing = {'uniq_raw', 'uniq_lnorm'} - set(header)
        if missing:
            sys.exit(f'Missing columns {missing} in {path}: run gqlib with --restrict_metrics including raw,lnorm')
        # gqlib writes the combined columns only when there are ambiguous counts
        has_ambig = 'combined_raw' in header
        prefixes = ('uniq', 'combined') if has_ambig else ('uniq', 'uniq')
        cols = [header.index(f'{p}_{m}') for p in prefixes for m in ('raw', 'lnorm')]
        counts = {}
        for line in f:
            row = line.rstrip('\n').split('\t')
            counts[row[0]] = [float(row[c]) for c in cols]
    return counts, has_ambig


def read_read_counts(path):
    with open(path) as f:
        return dict((k, int(v)) for k, v in (line.split() for line in f))


def read_annotation(path, gene_counts, gene_col, feature_cols, sep, feature_sep):
    # returns {category: {gene: set(features)}}, only for genes with counts
    with open_text(path) as f:
        lines = (line for line in f if line.strip() and not line.startswith('#'))
        header = next(lines).rstrip('\n').split(sep)
        gene_col = gene_col or header[0]
        feature_cols = feature_cols.split(',') if feature_cols else [c for c in header if c != gene_col]
        missing = set([gene_col] + feature_cols) - set(header)
        if missing:
            sys.exit(f'Missing columns {missing} in {path}')
        gi = header.index(gene_col)
        fis = {c: header.index(c) for c in feature_cols}
        annotation = {c: {} for c in feature_cols}
        for line in lines:
            row = line.rstrip('\n').split(sep)
            gene = row[gi]
            if gene not in gene_counts:
                continue
            for c, fi in fis.items():
                # eggnog-mapper uses '-' for missing annotations
                features = {x.strip() for x in row[fi].split(feature_sep)} - {'', '-'}
                if features:
                    annotation[c].setdefault(gene, set()).update(features)
    return annotation


def rollup(gene_counts, has_ambig, read_counts, annotation, restrict_metrics, out_file):
    category = [0.0] * 4
    features = {}
    for gene, gene_features in annotation.items():
        counts = gene_counts[gene]
        category = [a + b for a, b in zip(category, counts)]
        for feature in gene_features:
            features[feature] = [a + b for a, b in zip(features.get(feature, [0.0] * 4), counts)]

    def scaling_factor(raw, lnorm):
        return raw / lnorm if lnorm else 0
    uniq_factor = scaling_factor(category[0], category[1])
    ambig_factor = scaling_factor(category[2], category[3])
    rpkm_factor = 1e9 / read_counts['filtered_reads'] if read_counts['filtered_reads'] else 0

    def output_row(counts):
        row = []
        blocks = ((counts[0], counts[1], uniq_factor),)
        if has_ambig:
            blocks += ((counts[2], counts[3], ambig_factor),)
        for raw, lnorm, factor in blocks:
            values = dict(raw=raw, lnorm=lnorm, scaled=lnorm * factor, rpkm=lnorm * rpkm_factor)
            row += [values[m] for m in METRICS if m in restrict_metrics]
        return row

    metrics = [m for m in METRICS if m in restrict_metrics]
    header = [f'uniq_{m}' for m in metrics] + ([f'combined_{m}' for m in metrics] if has_ambig else [])

    def write_row(name, values, stream):
        print(name, *(f'{c:.5f}' for c in values), sep='\t', file=stream)

    with gzip.open(out_file, 'wt') as out:
        print('feature', *header, sep='\t', file=out)
        for name in ('total_reads', 'filtered_reads'):
            if read_counts[name]:
                write_row(name, [read_counts[name]] * len(header), out)
        if annotation:
            write_row('category', output_row(category), out)
        for feature in sorted(features):
            write_row(feature, output_row(features[feature]), out)


def main():
    args = parse_args()
    gene_counts, has_ambig = read_gene_counts(args.gene_counts)
    read_counts = read_read_counts(args.read_counts)
    annotation = read_annotation(args.annotation, gene_counts, args.gene_col, args.feature_cols, args.sep, args.feature_sep)
    for category, category_annotation in annotation.items():
        rollup(
            gene_counts, has_ambig, read_counts, category_annotation, args.restrict_metrics,
            f'{args.out_prefix}.{category}.txt.gz',
        )


if __name__ == '__main__':
    main()
