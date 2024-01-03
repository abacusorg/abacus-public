#!/bin/bash
#SBATCH -N 100
##SBATCH -N 4
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=128
#SBATCH -t 0-12
#SBATCH -C gpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J Abacus_DESI2_c000_ph300
#SBATCH --mem=0
#SBATCH --mail-type=ALL

cd $HOME/abacus
. env/environment.sh
cd Production

./run_sim.py Abacus_DESI2_c000_ph300
abacus_exit=$?

echo "Abacus exited with code $abacus_exit"

exit $abacus_exit
