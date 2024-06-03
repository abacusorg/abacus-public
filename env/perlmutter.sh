# Create a Python venv in /global/common (recommended for software environments)
# $ module load python
# $ python -m venv --system-site-packages /global/common/software/desi/users/$USER/venv-abacus
# $ ln -s /global/common/software/desi/users/$USER/venv-abacus venv

module load PrgEnv-gnu
module load cudatoolkit
module load cray-fftw
module load evp-patch
module load python
export MPICH_GPU_SUPPORT_ENABLED=0

# keep everything in the python env
export PYTHONNOUSERSITE=1

# in the past darshan has caused problems
# module unload darshan

export ABACUS=$HOME/abacus
export ABACUS_TMP=$SCRATCH
export ABACUS_PERSIST=$SCRATCH
export ABACUS_SSD=/dev/shm/$USER

export LD_LIBRARY_PATH=$ABACUS/clibs:$ABACUS/external/gperftools/install/lib64:$LD_LIBRARY_PATH
export    LIBRARY_PATH=$ABACUS/clibs:$ABACUS/external/gperftools/install/lib64:$LIBRARY_PATH
export CPATH=$ABACUS/clibs:$ABACUS/external/gperftools/install/include:$CPATH
export PYTHONPATH=$ABACUS:$PYTHONPATH

. $ABACUS/external/oneTBB/build/gnu_12.3_cxx11_64_relwithdebinfo/vars.sh
# need to find the installed oneTBB before the system TBB
export LDFLAGS="$LDFLAGS -L$ABACUS/external/oneTBB/build/gnu_12.3_cxx11_64_relwithdebinfo"

. $ABACUS/venv/bin/activate

export CXX=CC
export CC=cc

export SBATCH_ACCOUNT=desi
export SALLOC_ACCOUNT=desi
