#!/bin/bash

# Create merger tree asssociations
# Depends on: post-groups

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 32
#SBATCH -t 0-8
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-associations
#SBATCH --mem=0

sim=$1

#########

set -e

if [ -z "$sim" ]; then
    echo "Usage: sbatch $0 simpath"
    exit 1
fi

$ABACUS_MERGERTREE/create_associations_slabwise.py \
    -inputdir $(dirname $sim) \
    -simname $(basename $sim) \
    -num_cores $SLURM_CPUS_PER_TASK \
    -outputdir $(dirname $sim)/merger/ \
    -num_epochs 100
