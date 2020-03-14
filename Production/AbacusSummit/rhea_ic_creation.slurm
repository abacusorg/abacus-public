#!/usr/bin/env bash

## This script runs on rhea and calls the zeldovich module to create ics. 

#SBATCH -A AST145
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -p batch
#SBATCH --mem 0  # all memory
#SBATCH -o logs/%x.out
#SBATCH -t 18:00:00  #run for X hours.

# Exit if any command fails
set -e

PARFN=$1

echo "Starting ICs using parameter file: $PARFN"
echo "Using zeldovich module in: $ABACUS"

$ABACUS/Abacus/zeldovich.py --no-parallel $PARFN

# This is a special string that the assembly line will look for
echo -e "\n\nIC creation complete."
