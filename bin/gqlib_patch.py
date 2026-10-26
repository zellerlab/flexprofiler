#!/usr/bin/env python3
# runs only the gene-level counting of gqlib (gffquant), the functional roll-up is done
# separately by gqlib_rollup.py
# - gqlib requires an annotation database, an empty one is used so that no roll-up is done
# - counts are written at full precision, so that the roll-up does not accumulate rounding errors
# - the read counts used for normalisation are written to <out_prefix>.read_counts.txt
# usage is the same as `gffquant --mode genes`, without --db

import logging
import os
import sys
import tempfile
from unittest.mock import patch
import gqlib.__main__
from gqlib.__main__ import main
from gqlib.annotation.count_writer import CountWriter
from gqlib.db.db_import import SmallDatabaseImporter
from gqlib.profilers.feature_quantifier import FeatureQuantifier

logger = logging.getLogger(__name__)
handle_args = gqlib.__main__.handle_args
finalise = FeatureQuantifier.finalise


def handle_args_empty_db(argv):
    if '--db' in argv:
        raise ValueError('--db is not used, the functional roll-up is done by gqlib_rollup.py')
    fd, empty_db = tempfile.mkstemp(suffix='.bed')
    os.close(fd)
    try:
        args = handle_args(argv + ['--db', empty_db])
        args.annotation_db = SmallDatabaseImporter(logger, empty_db, single_category='feature', db_format='bed')
    finally:
        os.remove(empty_db)
    if args.mode not in ('gene', 'genes'):
        raise ValueError('only --mode genes is supported')
    return args


def write_row_full_precision(header, data, stream=sys.stdout):
    print(header, *(f'{c:.12g}' for c in data), flush=True, sep='\t', file=stream)


def finalise_with_read_counts(self, *args, **kwargs):
    finalise(self, *args, **kwargs)
    # same names used by gqlib in the feature tables, plus all reads (including unaligned)
    with open(f'{self.out_prefix}.read_counts.txt', 'w') as f:
        print('total_reads', self.aln_counter['read_count'], sep='\t', file=f)
        print('filtered_reads', self.aln_counter['filtered_read_count'], sep='\t', file=f)
        print('all_reads', self.aln_counter['full_read_count'], sep='\t', file=f)


with (
    patch('gqlib.__main__.handle_args', new=handle_args_empty_db),
    patch.object(CountWriter, 'write_row', new=staticmethod(write_row_full_precision)),
    patch.object(FeatureQuantifier, 'finalise', new=finalise_with_read_counts),
):
    if __name__ == '__main__':
        sys.argv[0] = sys.argv[0].removesuffix('.exe')
        sys.exit(main())
