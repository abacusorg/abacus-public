#!/bin/bash

# Take the halo time slice catalog and merger trees and paint them onto the light cone.
# This is done serially by the build_mt.py script.
# Then make a reduced halo catalog that only has the halos that ended up in the light cone.
# This can be done in parallel on each redshift by the Abacus.post.halo_lc_gather module.

# Depends on: post-groups, post-associations, post-cleaning

#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 2-0
#SBATCH -p scc
#SBATCH --exclusive
#SBATCH -J post-halo-lc

#########

set -eu

cosm=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
simname=AbacusBacklight_base_c${cosm}_ph000
simpath=$HOME/ceph/AbacusBacklight/$simname

workdir=$simpath/post_log/$SLURM_JOBID-halo-lc
mkdir -p $workdir
cd $workdir

merge_cksum=$ABACUS/external/fast-cksum/bin/merge_checksum_files.py

outtmp=/dev/shm/$USER/halo_lc_out/
worktmp=/dev/shm/$USER/halo_lc_tmp/

srun -n 1 --cpu-bind=none \
    python -u $HOME/abacus_cookiecutter/build_mt.py $simpath \
    --output-parent=$outtmp \
    --tmpdir=$worktmp \
    --merger-dir=$simpath/merger &> $workdir/build_mt.log
srun -n 1 rm -rf $worktmp

echo "#DISBATCH PREFIX python -m Abacus.post.halo_lc_gather -d $outtmp/$simname/" > tasks
echo "#DISBATCH SUFFIX  $simpath &> halo_lc_gather_\${DISBATCH_TASKID}.txt" >> tasks
(cd $outtmp/$simname; find -maxdepth 1 -name 'z*' -type d) | sort >> tasks

disBatch -e --status-header tasks

for d in $simpath/lightcone_halos/z*/; do
    cd $d
    $merge_cksum --delete *.crc32 | sponge checksums.crc32
    cd -
done
