#!/usr/bin/env bash

# Sets up modules for running Abacus on summit or rhea

# PROJWORK: large, fast, no backup. For simulation data.
export _PROJWORK=${_PROJWORK:-$PROJWORK}
export PROJWORK=$_PROJWORK/ast145  # tack on ast145 so we don't have to type it every time

# PROJHOME: small, slow, backed up.  For code and scripts.
export PROJHOME=${PROJHOME:-"/ccs/proj/ast145"}

module purge

#Abacus enviroment variables
export ABACUS_SSD=/dev/shm/$USER
export ABACUS_TMP=/dev/shm/$USER
export ABACUS=$HOME/abacus
export ABACUS_PERSIST=$PROJWORK/$USER

# Load the abacus module to set the include and python paths
module use $ABACUS/modulefiles
module load abacus

# Contains the directories with the .par2 files
export ABACUSSUMMIT_SPEC=$ABACUS/external/AbacusSummit/Simulations
# Contains the AbacusSummit directory with simulation data
export ABACUSSUMMIT_PERSIST=$PROJWORK

# Common modules
module load git
module load hsi
module load nano

if [ "$LMOD_SYSTEM_NAME" = "summit" ]; then
    module load gcc/8.1.1
    module load cuda/10.1.243
    module load spectrum-mpi
    module load lsf-tools
    module load screen
    module load valgrind

    if [ -z "$TBBROOT" ]; then
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
