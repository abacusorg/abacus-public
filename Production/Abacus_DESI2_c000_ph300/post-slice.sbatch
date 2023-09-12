#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH -t 0-4
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-slice.Abacus_DESI2_c000_ph300
#SBATCH --mem=0

set -e

sim=$PSCRATCH/Abacus_DESI2_c000_ph300
chunk=10

cd $HOME/abacus
. env/environment.sh
cd Production

slices=$sim/slice*.*
workdir=$sim/post_log/$SLURM_JOBID-slice
disBatch=$HOME/abacus/external/disBatch/disBatch
merge_cksum=$HOME/abacus/external/fast-cksum/bin/merge_checksum_files.py

mkdir -p $workdir
cd $workdir

python -m Abacus.post.timeslice_tasks -c $chunk -l $workdir $slices > tasks
$disBatch -e --status-header tasks

for d in $sim/slices/z*/*; do
    cd $d
    $merge_cksum --delete *.crc32 > checksums
    mv checksums checksums.crc32
    cd -
done
