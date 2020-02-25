#!/bin/sh

# Sets up modules for running Abacus on summit or rhea

module purge

module use $ABACUS/modulefiles
module load abacus

# Common modules
module load git
module load hsi
module load nano

if [ "$LMOD_SYSTEM_NAME" = "summit" ]; then
    module load cuda/10.1.243
    module load gcc/8.1.1
    module load spectrum-mpi
    module load lsf-tools
    module load screen
    module load valgrind

    if [ -z "$TBBROOT" ]; then
        #source $HOME/tbb/build/linux_ppc64le_gcc_cc6.4.0_libc2.17_kernel4.14.0_release/tbbvars.sh
        #source $HOME/tbb/build/linux_ppc64le_xl_cc4.8.5_libc2.17_kernel4.14.0_release/tbbvars.sh
        #source $HOME/tbb/build/linux_ppc64le_gcc_cc7.4.0_libc2.17_kernel4.14.0_release/tbbvars.sh
        source $HOME/tbb/build/linux_ppc64le_gcc_cc8.1.1_libc2.17_kernel4.14.0_release/tbbvars.sh
    fi
fi

if [ "$LMOD_SYSTEM_NAME" = "rhea" ]; then
    module load openmpi
    module load slurm
    module load gcc
    module load gsl
fi

# Needs to happen after compiler is loaded
module load fftw