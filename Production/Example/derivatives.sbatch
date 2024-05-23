#!/bin/bash

# Create derivatives.
# Depends on: none

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 64
#SBATCH -t 0-4
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J derivatives
#SBATCH --mem=0
##SBATCH --mail-type=ALL

#########

if [ $# -gt 0 ]; then
    echo "Usage: sbatch $0"
    exit 1
fi

sim=.

$ABACUS/Production/run_sim.py --just-params $sim
python -m Abacus.derivatives $sim/abacus.par
exit_code=$?

echo "Exited with code $exit_code"
