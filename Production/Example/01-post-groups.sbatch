#!/bin/bash

# Process the groups, packing the raw binary data into compressed ASDF
# Depends on: none

#SBATCH -N 2
#SBATCH --ntasks-per-node=8
#SBATCH -c 16
#SBATCH -t 0-1
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-groups
#SBATCH --mem=0

sim=$1
chunk=50
workers=4  # cpus-per-task / 4 is a good value

#########

set -e

if [ -z "$sim" ]; then
    echo "Usage: sbatch $0 simpath"
    exit 1
fi

groups=$sim/group/Step*
workdir=$sim/post_log/$SLURM_JOBID-groups

mkdir -p $workdir
python -m Abacus.post.group_tasks --delete -c $chunk -l $workdir -w $workers $groups > $workdir/tasks

cd $workdir
disBatch --fill -e --status-header tasks
