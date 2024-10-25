#!/bin/bash

# Create the cleaned halo catalogs.
# Depends on: post-associations

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 128
#SBATCH -t 0-8
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-cleaning
#SBATCH --mem=0

sim=$1

#########

set -e

if [ -z "$sim" ]; then
    echo "Usage: sbatch $0 simpath"
    exit 1
fi

for zdir in $sim/halos/z*; do
    python $HOME/abacus_mergertree/postprocess_trees/extract_mergerhistory_slabwise_slim_mpi.py \
        $(dirname $sim) \
        $(basename $sim) \
        $(basename ${zdir} | cut -c2-) \
        $(dirname $sim)/mergerhistory/ &
    pids+=($!)
done

for pid in ${pids[*]}; do
    wait $pid
done
unset pids

for zdir in $sim/halos/z*; do
    python $HOME/abacus_mergertree/postprocess_trees/build_particle_lists.py \
        $(basename $sim) \
        $(basename ${zdir} | cut -c2-) &
    pids+=($!)
done

for pid in ${pids[*]}; do
    wait $pid
done
