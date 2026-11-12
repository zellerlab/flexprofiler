#!/bin/bash
#SBATCH --job-name=pull_images
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
# Pull all the container images of the pipeline into the apptainer cache of the zellerlab_shark profile (shared by the
# lab), so that the runs find them there and do not download them. Run it once for each new version of the pipeline,
# as a job:
#
#   sbatch utils/pull_images.sh
#
# (or with bash on a login node). The images already in the cache are skipped. The list of the images comes from the
# pipeline itself (nextflow inspect), for all its processes, also those that a run does not use, and the cache from the
# zellerlab_shark profile.
set -euo pipefail

# be sure that this points to the pipeline version you want to run (as VERSION in run_shark.sh): a release tag, a
# branch or a commit
VERSION="0.1.0"
PROFILE="zellerlab_shark"

# as in run_shark.sh
module purge
module load library/java/OpenJDK-20.0.1
module load container/apptainer/1.5.3/gcc-8.5.0
module load bioinformatics/tools/Nextflow/25.04.6
export NXF_VER="26.10.5"

# the cache is shared by the lab: the images are made writable by the group
umask 002

# nextflow inspect writes a work directory and a log where it is run: keep them out of the submission directory
RUN_DIR=$(mktemp -d)
cd "$RUN_DIR"

PIPELINE="zellerlab/flexprofiler"
nextflow pull "$PIPELINE" -r "$VERSION"
CACHE_DIR=$(nextflow config "$PIPELINE" -r "$VERSION" -profile "$PROFILE" -value apptainer.cacheDir)
mkdir -p "$CACHE_DIR"

# the images of all the processes, as the pipeline gives them to apptainer
IMAGES=$(
    nextflow inspect "$PIPELINE" -r "$VERSION" -profile "$PROFILE" -format json \
        | grep -o '"container": *"[^"]*"' \
        | sed 's/^"container": *"//; s/"$//' \
        | sort -u
)
if [ -z "$IMAGES" ]; then
    echo "No container images found for $PIPELINE $VERSION" >&2
    exit 1
fi

# a partly downloaded image and the run directory are removed if the script stops, also when SLURM cancels the job or
# it reaches its time limit (SIGTERM)
TMP_IMAGE=""
cleanup() {
    if [ -n "$TMP_IMAGE" ]; then rm -f "$CACHE_DIR/$TMP_IMAGE"; fi
    rm -rf "$RUN_DIR"
}
trap cleanup EXIT
trap 'exit 143' TERM
trap 'exit 130' INT

pulled=0
cached=0
for image in $IMAGES; do
    # the name Nextflow gives to the image in the cache (SingularityCache.simpleName): the URL without its scheme, with
    # ':' and '/' replaced by '-', and the extension .img (.sif for the images named *.sif)
    name="${image#*://}"
    ext=".img"
    if [[ "$name" == *.sif:* ]]; then
        ext=".sif"
        name="${name/.sif:/-}"
    elif [[ "$name" == *.sif ]]; then
        ext=".sif"
        name="${name%.sif}"
    fi
    name="${name//:/-}"
    name="${name//\//-}$ext"
    if [ -s "$CACHE_DIR/$name" ]; then
        echo "cached:  $image"
        cached=$((cached + 1))
        continue
    fi
    echo "pulling: $image"
    # pulled to a temporary name and renamed when complete, as Nextflow does, so that a run never finds a partial image;
    # without the apptainer cache, which would keep a second copy of each image in the home
    TMP_IMAGE="$name.pulling.$$"
    (cd "$CACHE_DIR" && apptainer pull --disable-cache --name "$TMP_IMAGE" "$image" > /dev/null)
    mv "$CACHE_DIR/$TMP_IMAGE" "$CACHE_DIR/$name"
    TMP_IMAGE=""
    pulled=$((pulled + 1))
done

echo "Done: $pulled images pulled, $cached already in $CACHE_DIR"
