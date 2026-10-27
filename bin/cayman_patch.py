#!/usr/bin/env python3
# cayman profile always runs bwa mem itself and streams its output to the profiler
# this script runs the cayman profiler on pre-computed alignment files instead, so that
# the alignment can run in a separate process, mirroring cayman.functions.run_profile
# (cayman 0.10.2) minus the alignment
# alignment files (SAM/BAM) must be in the order emitted by bwa mem (grouped by read
# name, with unmapped reads included), and paired-end counting is applied based on
# the read group in the alignments (ID:2)

import argparse
import logging
import os
import pathlib
import sys
from unittest.mock import patch
from gqlib.db.db_import import SmallDatabaseImporter
from gqlib.profilers import RegionQuantifier
from cayman.handle_args import handle_args

# reuse the `cayman profile` arguments that are meaningful without the alignment: handle_args
# parses right away, so parse_args is patched to return the parser (it also sets up logging)
with patch.object(argparse.ArgumentParser, 'parse_args', new=lambda self, args=None: self):
    (subparsers,) = handle_args(sys.argv[1:])._subparsers._group_actions
parser = argparse.ArgumentParser()
for action in subparsers.choices['profile']._actions:
    if action.dest in ('annotation_db', 'out_prefix', 'min_identity', 'min_seqlen'):
        parser._add_action(action)
parser.add_argument('bam', nargs='+', help='Alignment files (SAM/BAM) produced by bwa mem.')
args = parser.parse_args()

logger = logging.getLogger()

# same input checks as cayman
for bam in args.bam:
    if not os.path.isfile(bam):
        raise ValueError(f"{bam} is not a valid alignment file")
if not os.path.exists(args.annotation_db):
    raise ValueError(f"{args.annotation_db} is not a valid annotation database")
if os.path.dirname(args.out_prefix):
    pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(exist_ok=True, parents=True)

db_format = None
with open(args.annotation_db) as db_in:
    for line in db_in:
        line = line.strip()
        if line and line[0] != '#':
            if ',' in line:
                db_format = 'hmmer'
                break
            if '\t' in line:
                db_format = 'bed'
                break
if db_format is None:
    logger.error("Cannot determine database format in %s.", args.annotation_db)
    raise ValueError(f"Cannot determine database format in {args.annotation_db}.")
logger.info("Identified database format as `%s`.", db_format)

db_importer = SmallDatabaseImporter(logger, args.annotation_db, single_category='cazy', db_format=db_format)
logger.info("Finished loading database.")

profiler = RegionQuantifier(
    db=db_importer,
    out_prefix=args.out_prefix,
    ambig_mode='1overN',
    reference_type='domain',
)
for bam in args.bam:
    logger.info("Reading alignments: %s", bam)
    try:
        # aln_format only sets the pysam open mode, and pysam autodetects SAM/BAM when reading
        with open(bam, 'rb') as aln_stream:
            profiler.count_alignments(aln_stream, aln_format='bam', min_identity=args.min_identity, min_seqlen=args.min_seqlen)
    except Exception as err:
        logger.error("Encountered problems digesting the alignment file %s:", bam)
        logger.error("%s", err)
        logger.error("Shutting down.")
        sys.exit(1)
profiler.finalise(restrict_reports=('raw', 'rpkm'))
