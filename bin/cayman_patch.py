#!/usr/bin/env python3
# cayman profile always runs bwa mem itself and streams its output to the profiler
# this patch replaces the aligner with a reader of pre-computed alignment files, so
# that the alignment can run in a separate process
# usage is the same as `cayman profile`, but the files passed as reads must be
# alignment files (SAM/BAM) in the order emitted by bwa mem (grouped by read name,
# with unmapped reads included) and the bwa index is not checked for existence

import sys
from types import SimpleNamespace
from unittest.mock import patch
from cayman.__main__ import main


class AlignmentFileReader:
    def __init__(self, *args, **kwargs):
        pass

    def run(self, input_files, single_end_reads=False):
        assert len(input_files) == 1, "pass a single alignment file per read set (e.g. --singles <bam>)"
        # cayman reads the alignments from proc.stdout, and pysam autodetects SAM/BAM
        return SimpleNamespace(stdout=open(input_files[0], 'rb')), f"reading {input_files[0]}"


with (
    patch('cayman.functions.BwaMemRunner', new=AlignmentFileReader),
    patch('cayman.functions.check_bwa_index', new=lambda prefix: True)
):
    if __name__ == '__main__':
        sys.argv[0] = sys.argv[0].removesuffix('.exe')
        sys.exit(main())
