#!/bin/bash
#SBATCH -q xfer
#SBATCH -t 1-0
#SBATCH --mem=48G
#SBATCH -c 16

sim=Abacus_DESI2_c000_ph300
hpssdir=/nersc/projects/desi/cosmosim/Abacus/$sim

cd $PSCRATCH/$sim

for d in slices/z*; do
    htar -T 4 -cvPf $hpssdir/slices/${sim}-slices-$(basename $d).tar ${d} &
done
wait

for d in halos*; do
    htar -T 16 -cvPf $hpssdir/${sim}-$d.tar $d &
done
wait

# One more quick job
toplevel="$(find -maxdepth 1 -type f) log/ info/"
htar -T $SLURM_CPUS_PER_TASK -cvPf "$hpssdir/${sim}-log_and_info.tar" $toplevel
