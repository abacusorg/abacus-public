#!/bin/bash

# Process the light cones, packing the raw binary into compressed ASDF
# Depends on: none

#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH -c 4
#SBATCH -t 0-4
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-lc
#SBATCH --mem=0

sim=$1

#########

set -e

if [ -z "$sim" ]; then
    echo "Usage: sbatch $0 simpath"
    exit 1
fi

workdir=$sim/post_log/$SLURM_JOBID-lc
disBatch=$ABACUS/external/disBatch/disBatch
merge_cksum=$ABACUS/external/fast-cksum/bin/merge_checksum_files.py

mkdir -p $workdir
cd $workdir

python -m Abacus.post.lc_tasks --delete -l $workdir $sim/lightcone > tasks
$disBatch -e --status-header tasks

for d in $sim/lightcone-post/*; do
    cd $d
    $merge_cksum --delete *.crc32 > checksums
    mv checksums checksums.crc32
    cd -
done

rm -rf $sim/lightcone
mv $sim/lightcone-post $sim/lightcones
