#!/bin/bash
#SBATCH -N 2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH -t 0-1
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-groups
#SBATCH --mem=0

set -e

sim=$PSCRATCH/Abacus_DESI2_c000_ph300
chunk=50
workers=4

cd $HOME/abacus
. env/environment.sh
cd Production

groups=$sim/group/Step*
workdir=$sim/post_log/$SLURM_JOBID-groups
disBatch=$HOME/abacus/external/disBatch/disBatch

mkdir -p $workdir
cd $workdir

python -m Abacus.post.group_tasks -c $chunk -l $workdir -w $workers $groups > tasks
$disBatch -e --status-header tasks
