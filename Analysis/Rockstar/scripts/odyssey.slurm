#!/bin/bash

#SBATCH -n 1280 # Number of tasks/processes (should be a multiple of 8 for rockstar)
##SBATCH -c 8 # Number of cores per task
##SBATCH -N 1 # Number of nodes to divide it across
##SBATCH --contiguous # do we need this?
##SBATCH --ntasks-per-node 4
#SBATCH -t 0-12:00 # Runtime in D-HH:MM
#SBATCH -p shared # Partition to submit to
##SBATCH --mem=128000 # Memory pool per node in MB (see also --mem-per-cpu)
#SBATCH --mem-per-cpu=16G
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
SLICEPAT="slice0.100"

# Find all slices; figure out which one this slurm array task will do
SLICES=$BASEDIR/$SIM_NAME/$SLICEPAT
NSLICE=$(echo $SLICES | wc -w)
SLICE=$(echo $SLICES | cut -d' ' -f $SLURM_ARRAY_TASK_ID)

#for SLICE in $BASEDIR/$SIM_NAME/$SLICEPAT; do

ZSLICE=$(echo $SLICE | sed 's/slice/z/' | sed "s,$SIM_NAME,$SIM_NAME\_products/$SIM_NAME\_rockstar_halos,")
CLI_CFG=$ZSLICE/auto-rockstar.cfg
rm -f $CLI_CFG

echo "Starting slice $SLICE server"
$ABACUS/Analysis/Rockstar/rockstar.py --ncpu=$SLURM_NTASKS --nnode=1 $SLICE --SO &

echo "Waiting for config file $CLI_CFG"

# Wait for the server to generate the client config file
while [ ! -f $CLI_CFG ]; do
    sleep 1
done

echo "Starting client in $ZSLICE"
stdbuf -o0 -e0 srun -u $ABACUS/Analysis/Rockstar/rockstar -c $CLI_CFG

#done
