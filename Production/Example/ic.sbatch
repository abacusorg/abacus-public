#!/bin/bash

# Create ICs.
# Depends on: none

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 64
#SBATCH -t 1-0
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J abacus-ic
#SBATCH --mem=0
##SBATCH --mail-type=ALL

#########

if [ $# -gt 0 ]; then
    echo "Usage: sbatch $0"
    exit 1
fi

sim=.

$ABACUS/Production/run_sim.py --just-params $sim
python -m Abacus.zeldovich $sim/abacus.par --no-parallel
#python -m Abacus.transpose_ic $HOME/pscratch/AbacusPNG_hugebase_c030_ph000/ic 2304 rvzel --nthread 32
exit_code=$?

echo "Exited with code $exit_code"
