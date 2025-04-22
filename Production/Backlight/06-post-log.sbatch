#!/bin/bash

# Compress the log files from the simulation and post-processing into zip files.
# Depends on: post-lc

#SBATCH -n 1
#SBATCH -t 0-0:30
#SBATCH -p scc
#SBTACH --exclusive
#SBATCH -J post-log

#########

set -eu

cosm=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
sim=$HOME/ceph/AbacusBacklight/AbacusBacklight_base_c${cosm}_ph000

cd $sim

rm -f log/last

srun zip --move -qr log.zip log/
srun zip --move -qr post_log.zip post_log/

$ABACUS/external/fast-cksum/bin/fast_cksum *.zip >> checksums.crc32
