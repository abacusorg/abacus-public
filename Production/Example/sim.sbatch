#!/bin/bash

# Run the simulation.
# Depends on: ic, derivatives

#SBATCH -N 32
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=4
#SBATCH -c 128
#SBATCH -t 0-20
#SBATCH -C gpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J abacus-sim
#SBATCH --mem=0
#SBATCH --mail-type=ALL

#########

if [ $# -gt 0 ]; then
    echo "Usage: sbatch $0"
    exit 1
fi

sim=.

$ABACUS/Production/run_sim.py $sim
abacus_exit=$?

echo "Abacus exited with code $abacus_exit"

exit $abacus_exit
