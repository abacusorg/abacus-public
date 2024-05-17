#!/bin/bash

# Process the time slices, packing the packN binary into compressed ASDF
# Depends on: none

#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH -c 8
#SBATCH -t 0-1
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-slice
#SBATCH --mem=0

sim=$1
chunk=10

#########

set -e

if [ -z "$sim" ]; then
    echo "Usage: sbatch $0 simpath"
    exit 1
fi

slices=$sim/slice*.*
workdir=$sim/post_log/$SLURM_JOBID-slice
disBatch=$ABACUS/external/disBatch/disBatch
merge_cksum=$ABACUS/external/fast-cksum/bin/merge_checksum_files.py

mkdir -p $workdir
cd $workdir

python -m Abacus.post.timeslice_tasks --delete -c $chunk -l $workdir $slices > tasks
$disBatch -e --status-header tasks

for d in $sim/slices/z*/*; do
    cd $d
    $merge_cksum --delete *.crc32 > checksums
    mv checksums checksums.crc32
    cd -
done
