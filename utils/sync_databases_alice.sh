#!/bin/bash
#SBATCH --job-name=sync_databases_alice
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
# Copy the lab's flexprofiler databases from SHARK to ALICE, with the same layout under the destination directory, so
# that a path on ALICE is the SHARK path with SRC replaced by DEST. Run it on SHARK, as a job (its compute nodes reach
# ALICE through the "alice" host of ~/.ssh/config):
#
#   sbatch utils/sync_databases_alice.sh                       # everything in SRC (~200 GB)
#   sbatch utils/sync_databases_alice.sh motus_v4.0.4 gmgc10   # only these subdirectories of SRC
#
# (or with bash on a login node). Files already on ALICE with the same size and modification time are skipped, so it can
# be run again to resume an interrupted copy or to add new databases. Nothing is deleted on ALICE.
#
# The files are copied with several rsync processes at once (JOBS), each with its own share of the files, balanced by
# size: one stream through the ALICE SSH gateway reaches ~40-75 MB/s, four ~165 MB/s.
#
# At the end, the database sheet of ALICE (assets/zellerlab_alice_databases.csv) is made from the one of SHARK
# (assets/zellerlab_shark_databases.csv), with SRC replaced by DEST in the paths, keeping only the rows whose files are
# all on ALICE (also those copied by earlier runs).
set -euo pipefail

SRC="${SRC:-/exports/archive/lucid-grpzeller-primary/SHARED/DATA/flexprofiler_databases}"
DEST="${DEST:-/projects/pi-zeller/SHARED/flexprofiler_databases}"
ALICE_HOST="${ALICE_HOST:-alice}"
JOBS="${JOBS:-4}"

# the pipeline directory, for the database sheets: sbatch runs a copy of the script in its spool directory, so under
# SLURM it is found from the path of the submitted script
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT=$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p')
else
    SCRIPT="${BASH_SOURCE[0]}"
fi
PIPELINE_DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
SHARK_SHEET="$PIPELINE_DIR/assets/zellerlab_shark_databases.csv"
ALICE_SHEET="$PIPELINE_DIR/assets/zellerlab_alice_databases.csv"
if [ ! -f "$SHARK_SHEET" ]; then
    echo "Database sheet not found: $SHARK_SHEET" >&2
    exit 1
fi

cd "$SRC"
if [ "$#" -gt 0 ]; then
    SUBDIRS=("$@")
else
    SUBDIRS=(*)
fi
for d in "${SUBDIRS[@]}"; do
    if [ ! -e "$d" ]; then
        echo "Not found in $SRC: $d" >&2
        exit 1
    fi
done

LIST_DIR=$(mktemp -d)
trap 'rm -rf "$LIST_DIR"' EXIT

# the files and symlinks to copy, relative to SRC, largest first, each given to the bucket with the fewest bytes so far
find "${SUBDIRS[@]}" \( -type f -o -type l \) -printf '%s\t%p\n' \
    | sort -t $'\t' -k1,1nr \
    | awk -F '\t' -v n="$JOBS" -v dir="$LIST_DIR" '
        { b = 0; for (i = 1; i < n; i++) if (tot[i] < tot[b]) b = i
          tot[b] += $1; print $2 > (dir "/bucket_" b) }
        END { for (i = 0; i < n; i++) printf "bucket %d: %.1f GB\n", i, tot[i] / 1e9 }'
cat "$LIST_DIR"/bucket_* > "$LIST_DIR/all"
echo "$(wc -l < "$LIST_DIR/all") files from $SRC to $ALICE_HOST:$DEST"

ssh "$ALICE_HOST" mkdir -p "$DEST"

# --files-from keeps the paths relative to SRC (creating the directories); the .fq.gz/.gz files are already compressed,
# so no -z. Directories are made group-writable and setgid, for the lab to share them. Symlinks pointing outside SRC or
# with an absolute target (e.g. pathseq_t2t/t2t_chm13v2/*.fna.gz -> host_genome/...) would be broken on ALICE: the file
# they point to is copied in their place
RSYNC_OPTS=(-a --partial --copy-unsafe-links --no-group --chmod=D2775,F664 --files-from)
pids=()
for bucket in "$LIST_DIR"/bucket_*; do
    rsync "${RSYNC_OPTS[@]}" "$bucket" "$SRC/" "$ALICE_HOST:$DEST/" &
    pids+=("$!")
done
failed=0
for pid in "${pids[@]}"; do
    wait "$pid" || failed=1
done
if [ "$failed" -ne 0 ]; then
    echo "Some rsync processes failed: run the script again to resume" >&2
    exit 1
fi

# check: a dry run over all the files must find nothing left to copy
LEFT=$(rsync "${RSYNC_OPTS[@]}" "$LIST_DIR/all" --dry-run --itemize-changes "$SRC/" "$ALICE_HOST:$DEST/" | grep -c '^[<>c]' || true)
if [ "$LEFT" -ne 0 ]; then
    echo "$LEFT files still differ on ALICE: run the script again" >&2
    exit 1
fi
echo "Done: all the files are on $ALICE_HOST:$DEST"

# the database sheet of ALICE: the paths (db_path, annot_db_path) under SRC moved to DEST, and only the rows whose paths
# all exist on ALICE. The paths are checked with one connection
MOVED_SHEET="$LIST_DIR/moved.csv"
sed "s#${SRC}/#${DEST}/#g" "$SHARK_SHEET" > "$MOVED_SHEET"
tail -n +2 "$MOVED_SHEET" | awk -F, '{ print $4; if ($6 != "") print $6 }' | sort -u \
    | ssh "$ALICE_HOST" 'while IFS= read -r p; do if [ -e "$p" ]; then printf "%s\n" "$p"; fi; done' \
    > "$LIST_DIR/present"
awk -F, -v present="$LIST_DIR/present" '
    BEGIN { while ((getline p < present) > 0) ok[p] = 1 }
    NR == 1 { print; next }
    ($4 in ok) && ($6 == "" || $6 in ok) { print; next }
    { print "Not on ALICE, left out of the sheet: " $1 " " $2 > "/dev/stderr" }
' "$MOVED_SHEET" > "$LIST_DIR/alice.csv"
cp "$LIST_DIR/alice.csv" "$ALICE_SHEET"
echo "Wrote $ALICE_SHEET ($(( $(wc -l < "$ALICE_SHEET") - 1 )) databases)"
