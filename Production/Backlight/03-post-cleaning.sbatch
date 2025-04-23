#!/bin/bash

# Create the cleaned halo catalogs.
# Depends on: post-associations

#SBATCH -t 0-4
#SBATCH -N 4
#SBATCH --ntasks-per-node=10
#SBATCH -p scc
#SBATCH --exclusive
#SBATCH -J post-cleaning
#SBATCH --mem=0

#########

set -eu

cosm=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
simdir=$HOME/ceph/AbacusBacklight/AbacusBacklight_base_c${cosm}_ph${cosm}

workdir=$simdir/post_log/$SLURM_JOBID-cleaning
mkdir -p $workdir

cd $workdir

# write the disBatch tasks
echo "#DISBATCH PREFIX (python $ABACUS_MERGERTREE/postprocess_trees/extract_mergerhistory_slabwise_slim_mpi.py $simdir/merger $simdir " > tasks
echo "#DISBATCH SUFFIX  $simdir/mergerhistory/ ) &> clean-\$DISBATCH_TASKID.log" >> tasks

# sorted zdirs
zdirs=( $(ls -d $simdir/halos/z* | sort) )

for zdir in ${zdirs[@]}; do
   echo "$(basename ${zdir} | cut -c2-)" >> tasks
done

echo "#DISBATCH BARRIER CHECK" >> tasks

# build up the next task list
echo "#DISBATCH PREFIX (python $ABACUS_MERGERTREE/postprocess_trees/build_particle_lists.py $simdir $simdir 0 " >> tasks
echo "#DISBATCH SUFFIX ) &> particles-\$DISBATCH_TASKID.log" >> tasks

for zdir in ${zdirs[@]}; do
    echo "$(basename ${zdir} | cut -c2-)" >> tasks
done

disBatch -e --status-header tasks
