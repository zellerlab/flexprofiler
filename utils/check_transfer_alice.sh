#!/bin/bash
# Check that the results of a run of utils/run_alice.sh are all on SHARK, the same as on ALICE, before deleting the run
# directory (RUN_DIR) on ALICE by hand. Run it on ALICE, or from SHARK without logging in:
#
#   bash check_transfer_alice.sh <RUN_DIR>
#   ssh alice bash -s -- <RUN_DIR> < check_transfer_alice.sh
#
# It succeeds (exit status 0) only if run_alice.sh finished copying the results to SHARK (RUN_DIR/run_files/
# results_copied, removed when the run is submitted again, so not there while a run is going) and an rsync dry run finds
# no result file missing or different on SHARK. It deletes nothing.
set -euo pipefail

SHARK_HOST="shark"

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <RUN_DIR>" >&2
    exit 1
fi
RUN_DIR="${1%/}"
FILES_DIR="$RUN_DIR/run_files"
MARKER="$FILES_DIR/results_copied"

if [ ! -f "$FILES_DIR/params.shark.yaml" ]; then
    echo "NOT OK: $RUN_DIR is not the directory of a run of run_alice.sh (no run_files/params.shark.yaml)" >&2
    exit 1
fi
if [ ! -f "$MARKER" ]; then
    echo "NOT OK: the results of $RUN_DIR were not copied to SHARK (the run is still going, failed, or was submitted again)" >&2
    exit 1
fi
SHARK_OUTDIR=$(cat "$MARKER")

LEFT=$(rsync -a --dry-run --itemize-changes --files-from="$FILES_DIR/results.txt" \
    "$RUN_DIR/results/" "$SHARK_HOST:$SHARK_OUTDIR/" | grep -c '^[<>c]' || true)
if [ "$LEFT" -ne 0 ]; then
    echo "NOT OK: $LEFT result files are missing or differ in $SHARK_HOST:$SHARK_OUTDIR" >&2
    exit 1
fi
echo "OK: all the $(wc -l < "$FILES_DIR/results.txt") result files are in $SHARK_HOST:$SHARK_OUTDIR"
echo "The run directory can be deleted: rm -rf $RUN_DIR"
