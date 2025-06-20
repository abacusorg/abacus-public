#!/bin/bash

# Process the light cones, packing the raw binary into compressed ASDF
# Depends on: none

#SBATCH -n 32
#SBATCH -c 4
#SBATCH --mem-per-cpu=64G
#SBATCH -t 1-0
#SBATCH -p scc
#SBATCH --exclusive=USER
#SBATCH -J post-lc

#########

set -eu

cosm=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
sim=$HOME/ceph/AbacusBacklight/AbacusBacklight_base_c${cosm}_ph${cosm}
target_width=30  # Mpc/h
bit_trunc=0
dtype=quantize
coarsify_trans_vel=2

workdir=$sim/post_log/$SLURM_JOBID-lc
merge_cksum=$ABACUS/external/fast-cksum/bin/merge_checksum_files.py

mkdir -p $workdir
python -m Abacus.post.lc_tasks \
    -s \
    -w $target_width \
    -b $bit_trunc \
    -d $dtype \
    -l $workdir \
    -a greedy \
    --coarsify-trans-vel $coarsify_trans_vel \
    --discard-partialsky \
    --delete \
    $sim/lightcone \
    > $workdir/tasks

cd $workdir
disBatch --fill -e --status-header tasks
cd -

for d in $sim/lightcone-post/*/*/; do
    cd $d
    $merge_cksum --delete *.crc32 | sponge checksums.crc32
    cd -
done

rm -rf $sim/lightcone
mv $sim/lightcone-post $sim/lightcone_healpix
