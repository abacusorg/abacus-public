#!/bin/bash
#SBATCH -N 8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH -t 0-4
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-power
#SBATCH --mem=0

set -e

sim=$PSCRATCH/Abacus_DESI2_c000_ph300
nfft=2048
nthread=30
nreader=2

cd $HOME/abacus
. env/environment.sh
cd Production

inputs=$(echo $sim/{halos,slices}/z*)
workdir=$sim/post_log/$SLURM_JOBID-groups
disBatch=$HOME/abacus/external/disBatch/disBatch

mkdir -p $workdir
cd $workdir

python -m Abacus.post.power_tasks --nfft $nfft -t $nthread -r $nreader -l $workdir $inputs > tasks
$disBatch -c $SLURM_CPUS_PER_TASK -e --status-header tasks
