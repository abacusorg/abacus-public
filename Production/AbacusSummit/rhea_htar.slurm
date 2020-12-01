#!/usr/bin/env bash
#SBATCH -A AST145
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH -p batch
#SBATCH --mem 0  # all memory
#SBATCH -t 24:00:00  # dtn wall clock unlimited

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
HPSSDIR=/hpss/prod/ast145/proj-shared/$SIM_SET/$SIM_NAME
cd $WORKINGDIR

NTHREAD=${SLURM_CPUS_PER_TASK:-16}

echo "Executing $NTHREAD parallel htar jobs on $WORKINGDIR"

#htar -T $NTHREAD -cPvf "/hpss/prod/ast145/proj-shared/$SIM_SET/${SIM_NAME}.tar" $WORKINGDIR

TARJOBS="$(find -maxdepth 1 -type d -name halos -o -name lightcones) $(find -maxdepth 2 -type d -path './slices/z*')"
echo $TARJOBS
hsi "mkdir $HPSSDIR; cd $HPSSDIR; mkdir -p $(dirname ${TARJOBS} | sort -u)"

for JOB in $TARJOBS; do
    htar -T $NTHREAD -cvf $HPSSDIR/${JOB}.tar $JOB
done

# One more quick job
TOPLEVEL="$(find -maxdepth 1 -type f) log/ info/"
htar -T $NTHREAD -cvf "$HPSSDIR/log_and_info.tar" $TOPLEVEL

# This is a magic string that the assembly line looks for
echo "Abacus htar complete."
