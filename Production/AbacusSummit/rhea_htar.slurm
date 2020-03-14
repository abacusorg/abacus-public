#!/usr/bin/env bash
#SBATCH -A AST145
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p batch
#SBATCH --mem 0  # all memory
#SBATCH -t 24:00:00  # 48 hours is max rhea wall clock

set -e  # quit on error

# This script runs htar to store files on tape (HPSS).
# Deletion of the originals will be taken care of later by the assembly line
# once Globus is confirmed to have completed too.

# While it's not obvious that one needs an allocation to run htar,
# it is multi-threaded and it thus seems prudent to avoid using a login node.
# We could imagine using the DTNs, but it's hard to know if a job is running or crashed

SIM_SET=${SIM_SET:-AbacusSummit}
SIM_NAME=$1

WORKINGDIR=$ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME

NTHREAD=${SLURM_CPUS_PER_TASK:-16}

echo "Executing htar on $WORKINGDIR using $NTHREAD threads"

htar -T $NTHREAD -cPvf "/hpss/prod/ast145/proj-shared/$SIM_SET/${SIM_NAME}.tar" $WORKINGDIR

# This is a magic string that the assembly line looks for
echo "Abacus htar complete."
