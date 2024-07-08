#!/bin/bash

# Compute the matter power spectrum from the particle data in the halo catalogs
# and/or time slices. Uses a grid-based FFT method.
# Depends on: post-groups, post-slice

#SBATCH -N 8
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH -t 0-4
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-power
#SBATCH --mem=0

sim=$1
nfft=2048
nthread=30
nreader=2

#########

set -e

if [ -z "$sim" ]; then
    echo "Usage: sbatch $0 simpath"
    exit 1
fi

inputs=$(echo $sim/{halos,slices}/z*)
workdir=$sim/post_log/$SLURM_JOBID-power
disBatch=$HOME/abacus/external/disBatch/disBatch

mkdir -p $workdir
cd $workdir

python -m Abacus.post.power_tasks --nfft $nfft -t $nthread -r $nreader -l $workdir $inputs > tasks
$disBatch -e --status-header tasks
