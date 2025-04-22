#!/bin/bash

# Process the groups, packing the raw binary data into compressed ASDF
# Depends on: none

#SBATCH -N 2
#SBATCH -c 4
#SBATCH -t 0-2
#SBATCH -p scc
#SBATCH --exclusive
#SBATCH -J post-groups

chunk=1000
workers=1  # max(cpus-per-task / 4, n_superslab) is a good value

#########

set -eu

cosm=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
sim=$HOME/ceph/AbacusBacklight/AbacusBacklight_base_c${cosm}_ph000

groups=$sim/group/Step*
workdir=$sim/post_log/$SLURM_JOBID-groups

mkdir -p $workdir
python -m Abacus.post.group_tasks --delete -c $chunk -l $workdir -w $workers $groups > $workdir/tasks

cd $workdir
disBatch --fill -e --status-header tasks
