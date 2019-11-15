#!/bin/bash

##SBATCH -n 16 # Number of tasks/processes (should be a multiple of 8 for rockstar)
#SBATCH -c 1 # Number of cores per task
#SBATCH -N 1 # Number of nodes to divide it across
##SBATCH --contiguous # do we need this?
##SBATCH --ntasks-per-node 4
#SBATCH -t 0-24:00 # Runtime in D-HH:MM
#SBATCH -p shared # Partition to submit to
##SBATCH --mem=128000 # Memory pool per node in MB (see also --mem-per-cpu)
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lgarrison@flatironinstitute.org # Email to which notifications will be sent
##SBATCH --array=2-38
##SBATCH --exclusive
##SBATCH -x holy7c07301  # bad node?

source $HOME/.bashrc
ulimit -c 0

SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

module list

# Set the location of the time slices
SIM_NAME="DESI_L2000_N7000_Summit_June19"
BASEDIR=/n/scratchlfs/eisenstein_lab/lgarrison/

#SLICEPAT="slice*"
SLICEPAT="z0.100"

# Find all slices; figure out which one this slurm array task will do
CATALOGS=$BASEDIR/$SIM_NAME\_products/$SIM_NAME\_rockstar_halos/$SLICEPAT
NCAT=$(echo $CATALOGS | wc -w)
CAT=$(echo $CATALOGS | cut -d' ' -f $SLURM_ARRAY_TASK_ID)

#for SLICE in $BASEDIR/$SIM_NAME/$SLICEPAT; do

#echo "Starting find_parents_binary on $CAT"
# find_parents needs to load all the halos at once
#$ABACUS/Analysis/Rockstar/util/find_parents_binary -c $CAT/rockstar.cfg
#echo "Done with exit code $?"

echo "Starting convert_to_HDF5.py on $CAT"
# This could be done in parallel!
$ABACUS/Analysis/Rockstar/util/convert_to_HDF5.py --nproc=1 "$CAT"
echo "Done with exit code $?"

#done
