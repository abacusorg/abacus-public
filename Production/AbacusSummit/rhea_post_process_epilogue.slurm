#!/usr/bin/env bash
#SBATCH -A AST145
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -p batch
#SBATCH --mem 0  # all memory
#SBATCH -t 06:00:00  # 48 hours is max rhea wall clock

set -e  # quit on error
shopt -s nullglob

# This script echos success to the log file

SIM_SET=${SIM_SET:-AbacusSummit}
SIM_NAME=$1

WORKINGDIR=$ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME
pushd $WORKINGDIR > /dev/null  # work inside sim dir

echo "Executing epilogue on $WORKINGDIR"

# Zip and delete the logs in groups of 100 steps
echo "Starting zip of logs at $(date)"
pushd log/ > /dev/null
i=0
while : ; do
    start=$(printf 'step%02d??' $i)
    LOGSTEPDIRS=($start/)
    if [ "${#LOGSTEPDIRS[@]}" -eq 0 ]; then
        break
    fi
    name=$(printf 'log%02d00' $i)
    zip --move -qr ${name}.zip ${LOGSTEPDIRS[@]} &
    i=$((i + 1))
done

wait  # wait for all zip to finish
echo "Done zip at $(date)"

# Checksum the logs
$ABACUS/external/fast-cksum/bin/fast_cksum log*.zip > checksums.crc32

popd > /dev/null  # back to WORKINGDIR

# Now delete:
# - Log/last symlink
# - Wisdom file: no sense in carting around
# - Read directory (?): the final header is in the halos anyway, and the full state is in the logs
# - Any disbatch files, or maybe we'll want to be able to check on them later

echo "Deleting read"
rm -rf read/

echo "Deleting extra files"
rm -f fftw_*.wisdom
rm -f log/last  # symlink
popd > /dev/null  # back to submission dir
#rm tmp/${SIM_NAME}.*

# This is a special string the assembly line looks for to know that it's safe to start Globus and htar
echo "Post processing complete."
