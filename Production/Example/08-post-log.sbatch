#!/bin/bash

# Compress the log files from the simulation and post-processing into zip files.
# Depends on: all

#SBATCH -N 1
#SBATCH -t 0-0:10
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-log
#SBATCH --mem=0

sim=$1

#########

set -e

if [ -z "$sim" ]; then
    echo "Usage: sbatch $0 simpath"
    exit 1
fi

cd $sim/log

rm -f last

i=0
while : ; do
    stepdirs=$(find -maxdepth 1 -name "$(printf 'step%02d??' $i)")
    if [ -z "$stepdirs" ]; then
        break
    fi
    name=$(printf 'log%02d00' $i)
    zip --move -qr ${name}.zip $stepdirs &
    pids+=($!)
    i=$((i + 1))
done

for pid in ${pids[*]}; do
    wait $pid
done
unset pids

$ABACUS/external/fast-cksum/bin/fast_cksum *.zip > checksums.crc32

cd $sim/post_log

for d in */; do
    zip --move -qr $(basename $d).zip $d &
    pids+=($!)
done

for pid in ${pids[*]}; do
    wait $pid
done
unset pids

$ABACUS/external/fast-cksum/bin/fast_cksum *.zip > checksums.crc32
