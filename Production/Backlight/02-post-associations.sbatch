#!/bin/bash

# Create merger tree asssociations
# Depends on: post-groups

#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 0-12
#SBATCH -p scc
#SBATCH --exclusive
#SBATCH -J post-associations

#########

set -eu

cosm=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
simdir=$HOME/ceph/AbacusBacklight/AbacusBacklight_base_c${cosm}_ph${cosm}

logdir=$simdir/post_log/$SLURM_JOBID-merger
mkdir -p $logdir

srun python $ABACUS_MERGERTREE/create_associations_slabwise.py \
    -inputdir $(dirname $simdir) \
    -simname $(basename $simdir) \
    -num_cores $SLURM_CPUS_PER_TASK \
    -outputdir $simdir/merger/ \
    -num_epochs 100 \
    &> $logdir/associations.log
