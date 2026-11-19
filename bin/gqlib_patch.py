#!/usr/bin/env python3
# Runs gqlib (python -m gqlib, with the same arguments), writing the profiles also
# when no alignment passes its filters (e.g. a sample of which no read aligns to
# the gene catalogue). gqlib 2.14.4 then writes only the alignment stats
# (FeatureQuantifier.finalise skips process_counters when aln_count is 0), and
# the modules fail on the missing profiles. The patch runs process_counters
# anyway, so that gqlib writes the gene counts with only their header; the tables
# of the functional categories are written only for the categories with counts,
# as without the patch, so none.

import logging
import runpy
import sys

from gqlib.profilers.feature_quantifier import FeatureQuantifier

logger = logging.getLogger(__name__)

finalise = FeatureQuantifier.finalise


def finalise_allowing_no_alignments(self, restrict_reports=None, report_category=False, report_unannotated=False, dump_counters=False):
    finalise(
        self,
        restrict_reports=restrict_reports,
        report_category=report_category,
        report_unannotated=report_unannotated,
        dump_counters=dump_counters,
    )
    if not self.aln_counter.get("aln_count"):
        logger.warning("No alignment passed the filters: writing empty profiles")
        self.process_counters(
            restrict_reports=restrict_reports,
            report_category=report_category,
            report_unannotated=report_unannotated,
            dump_counters=dump_counters,
        )


FeatureQuantifier.finalise = finalise_allowing_no_alignments

if __name__ == "__main__":
    sys.argv[0] = "gqlib"
    runpy.run_module("gqlib", run_name="__main__", alter_sys=True)
