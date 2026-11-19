#!/usr/bin/env python3
# mOTUs 4 does not support anymore the -db flag to pass a custom database
# the database path is hardcoded in the module, and this patch overrides it to
# point to the local work directory

# mOTUs 4 assigns an insert that aligns equally well to several genes of the
# same marker gene cluster to one of them at random, without a seed, so the
# profiles change slightly between runs. The seed makes them reproducible.

# mOTUs 4 quits with exit code 1, without writing the profile, when no read
# aligns to the marker genes ("The alignmentfile ... has no valid alignments").
# The patch writes instead an empty table of marker gene cluster counts, from
# which mOTUs writes an empty profile (only its header lines), as mOTUs 3 does.

import sys,pathlib,random,logging
from unittest.mock import patch
import pysam
from motus import motus, mutils
from motus.mentities import MOTUS_PARAMETERS
from motus.motus import main

random.seed(0)

workdir = pathlib.Path.cwd()

calc_mgc = motus.calc_mgc

def calc_mgc_allowing_no_alignments():
    # mOTUs stops when the BAM written by its alignment step has no record (it
    # keeps only the aligned reads)
    with pysam.AlignmentFile(str(MOTUS_PARAMETERS.get_alignment_file()), 'r') as alignments:
        has_alignments = next(alignments, None) is not None
    if has_alignments:
        return calc_mgc()
    logging.warning('No read aligned to the marker genes: writing an empty profile')
    # as mOTUs logs it, for the read counts (bin/read_counts.py)
    logging.info('Read 0 aligned inserts of which 0.0% are multimappers')
    # the header of the table written by calc_mgc, without any marker gene cluster
    with open(MOTUS_PARAMETERS.get_mgc_file(), 'w') as handle:
        handle.write(f'{mutils.create_mgc_header_line()}\n')
        handle.write('MGC\tINSERT_RAW\tINSERT_NORM\tINSERT_SCALED\tBASE_RAW\tBASE_NORM\n')
    return None

with (
    patch('motus.mutils.DEFAULT_MOTUS_MGDB_PARENT_LOCATION', new=workdir),
    patch('motus.mutils.DEFAULT_MOTUS_MGDB_LOCATION', new=workdir.joinpath('db_mOTU')),
    patch('motus.mutils.DEFAULT_MOTUS_MGDB_LOCATION_MARKER', new=workdir.joinpath('db_mOTU/db_mOTU.downloaded')),
    patch('motus.motus.calc_mgc', new=calc_mgc_allowing_no_alignments)
):
    if __name__ == '__main__':
        sys.argv[0] = sys.argv[0].removesuffix('.exe')
        sys.exit(main())
