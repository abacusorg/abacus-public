#!/bin/bash

# Copy the output of a simulation to HPSS (tape)
# Depends on: all

#SBATCH -q xfer
#SBATCH -t 0-21
#SBATCH --mem=48G
#SBATCH -c 16
#SBATCH -J post-hpss

simdir=$1
simname=$(basename $simdir)
hpssdir=/nersc/projects/desi/cosmosim/Abacus/$simname

#########

set -e

if [ -z "$simdir" ]; then
    echo "Usage: sbatch $0 simpath"
    exit 1
fi

cd $simdir

for d in slices/z*; do
    htar -T 1 -cvPf $hpssdir/slices/${simname}-slices-$(basename $d).tar ${d} &
done
wait

for d in halos*; do
    htar -T 16 -cvPf $hpssdir/${simname}-$d.tar $d &
done
wait

htar -T 16 -cvPf $hpssdir/${simname}-lightcones.tar lightcones

# One more quick job
toplevel="$(find -maxdepth 1 -type f) log/ info/ post_log/"
htar -T $SLURM_CPUS_PER_TASK -cvPf "$hpssdir/${simname}-log_and_info.tar" $toplevel
