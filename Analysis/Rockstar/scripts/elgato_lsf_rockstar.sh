#!/bin/bash
#BSUB -n 32  # Request number of cores
#BSUB -R "span[ptile=16]"
#BSUB -o lsf-%J_%I.out
#BSUB -q medium
##BSUB -R select[type==any]
#BSUB -x  # Request exclusive node
#BSUB -J "AbacusCosmosRockstar[1-40]"  # Inclusive; must start from 1
##BSUB -J "AbacusCosmos[1]"  # Inclusive; must start from 1
#BSUB -We 8:00  #Estimated runtime (non-binding)
#BSUB -R rusage[mem=200000]

# Don't run this script directly; submit it to the queue with
# bsub < lsf_script.sh
# Check this directory for a file named lsf-jobid.out, which should
# contain stderr and stdout

# Fail if any command returns a non-zero exit status
# If we intend to use requeueing, we'll have to be a little clever about
# capturing the requeue exit code.
set -e

# Set these if not running through LSF
#LSB_DJOB_NUMPROC=16
#LSB_JOBINDEX=6

export OMP_NUM_THREADS=$LSB_DJOB_NUMPROC

# Make a sim-specific tmp folder, in case we're running two runs on the same node
# We start folder numbers at 00, but LSF job arrays start at 01
#export SIM_DIR=$(printf "emulator_planck_PoP_test_00-%d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "emulator_planck_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "PLT_cosmo/PLT_cosmo_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "PLT_cosmo_08-%d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "2LPTic_cosmo_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "glass_cosmo/glass_cosmo_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "emulator_1100box_planck/emulator_1100box_planck_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "emulator_PoP_test_720box_planck_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR="emulator_1100box_planck_paired/emulator_1100box_planck_fixed_flipped_00-17"
#export SIM_DIR=$(printf "2LPTic_cosmo_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "emulator_720box_Neff3/emulator_720box_Neff3_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "AbacusCosmos_720box/AbacusCosmos_720box_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "AbacusCosmos_1100box/AbacusCosmos_1100box_%02d" $(expr $LSB_JOBINDEX - 1))
#export SIM_DIR=$(printf "AbacusCosmos_1100box_planck/AbacusCosmos_1100box_planck_00-%d" $(expr $LSB_JOBINDEX - 1))
export SIM_DIR=$(printf "AbacusCosmos_720box_planck/AbacusCosmos_720box_planck_00-%d" $(expr $LSB_JOBINDEX - 1))

# Split SIM_DIR into name and project
export SIM_NAME=$(echo "$SIM_DIR" | awk -F '/' -- '{print $2}' -)
export PROJECT=$(echo "$SIM_DIR" | awk -F '/' -- '{print $1}' -)

echo -e "* Source .bashrc:\n"
source ~/.bashrc
echo -e "\n\n\n\n"

echo -e "* Compute node:"
echo `hostname`
echo -e "\n\n\n\n"

rm -rf $ABACUS_TMP

# Note: this version of the script does not copy Abacus locally

echo -e "* Checking if we want to make high-resolution power spectra:\n"
if [ ]; then
  echo -e "Making power spectra."
  for SLICE in $ABACUS_PERSIST/$SIM_DIR/slice*/; do
    $ABACUS/Analysis/PowerSpectrum/run_PS.py --nfft 2048 $SLICE
  done  
else
  echo -e "No PS requested."
fi
echo -e "\n\n\n\n"


echo -e "* Checking if we need to run FoF:\n"
if [ 0 ]; then
  echo -e "Running FoF."
  $ABACUS/Analysis/FoF/FoF.py $ABACUS_PERSIST/$SIM_DIR/slice* --tar-mode TAR --tar-remove-source-files --boundary-slabs=4
else
  echo -e "No FoF requested."
fi
echo -e "\n\n\n\n"


echo -e "* Checking if we need to run Rockstar:\n"
if [ 0 ]; then
  echo -e "Running Rockstar."
  for SLICE in $ABACUS_PERSIST/$SIM_DIR/slice*; do
    ZSLICE=$ABACUS_PERSIST/$PROJECT/$SIM_NAME\_products/$SIM_NAME\_rockstar_halos/$(echo "$SLICE" | awk -F '/' '{print $NF}' - | sed 's/slice/z/')
    CLI_CFG=$ZSLICE/auto-rockstar.cfg
    rm -f $CLI_CFG  # Always remove an existing config file

    echo "Starting slice $SLICE server"
    $ABACUS/Analysis/Rockstar/rockstar.py --ncpu $LSB_DJOB_NUMPROC $SLICE --SO --tar-mode TAR --tar-remove-source-files &
    #$ABACUS/Analysis/Rockstar/rockstar.py --ncpu $LSB_DJOB_NUMPROC $SLICE --SO --tar-mode ONLY_TAR_INFER --tar-remove-source-files &

    # Wait for the server to generate the client config file
    while [ ! -f $CLI_CFG ]; do
        sleep 1
    done

    echo "Starting client in $ZSLICE"
    mpirun $ABACUS/Analysis/Rockstar/rockstar -c $CLI_CFG
  done
else
  echo -e "No Rockstar requested."
fi
echo -e "\n\n\n\n"


echo -e "* Checking if we need to tarball the outputs:\n"
if [ 0 ]; then
    $ABACUS/Abacus/archive_sim.py $ABACUS_PERSIST/$SIM_DIR --nthreads=$LSB_DJOB_NUMPROC --inplace --remove-source
else
  echo -e "No tar requested."
fi
echo -e "\n\n\n\n"


exit $ABACUS_EXIT
