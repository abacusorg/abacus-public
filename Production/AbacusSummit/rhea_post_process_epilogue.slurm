#!/usr/bin/env bash
#SBATCH -A AST145
#SBATCH -N 1
#SBATCH -p batch
#SBATCH --mem 0  # all memory
#SBATCH -o logs/%x.out
#SBATCH -t 06:00:00  # 48 hours is max rhea wall clock

set -e  # quit on error

# This script echos success to the log file

SIM_SET=${SIM_SET:-AbacusSummit}
SIM_NAME=$1

WORKINGDIR=$ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME
cd $WORKINGDIR  # work inside sim dir

echo "Executing epilogue on $WORKINGDIR"

# Zip and delete the logs  # --move
echo "Zipping logs"
time zip -qr log.zip log/

# Now delete:
# - Group directory: should basically be empty, since the converter should remove files when it's finished
# - IC directory: at this point, we're reasonably confident we won't need them again
# - Wisdom file: no sense in carting around
# - Read directory (?): the final header is in the halos anyway, and the full state is in the logs

echo "Deleting group"
#rm -rf group/

echo "Deleting ic"
#rm -rf ic/

echo "Deleting read"
#rm -rf read/

echo "Deleting extra files"
#rm fftw_*.wisdom
cd -  # back to submission dir
#rm tmp/${SIM_NAME}.{groupdirs,ldirs}

echo "Post processing complete."
