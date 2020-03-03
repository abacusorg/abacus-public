#!/usr/bin/env bash
#SBATCH -A AST145
#SBATCH -N 1
#SBATCH -p batch
#SBATCH --mem 0  # all memory
#SBATCH -o logs/%x.out
#SBATCH -t 00:01:00  # 48 hours is max rhea wall clock

# This script echo success to the log file

SIM_NAME=$1
echo "Post processing complete." >> logs/${SIM_NAME}_PostProcess.out
