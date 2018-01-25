#!/bin/bash

# Fail if any command returns a non-zero exit status
set -e

# Set these if not running through LSF
NTHREADS=24

export OMP_NUM_THREADS=$NTHREADS
export ABACUS_PERSIST="/mnt/gosling2/bigsims/AbacusCosmos_720box"

# Make a sim-specific tmp folder, in case we're running two runs on the same node
# We start folder numbers at 00, but LSF job arrays start at 01
export SIM_NAME="AbacusCosmos_720box_planck"

# Split SIM_DIR into name and project
echo -e "* Checking if we need to run Rockstar:\n"
if [ 0 ]; then
  echo -e "Running Rockstar."
  for SLICE in $ABACUS_PERSIST/$SIM_NAME/slice*; do
    ZSLICE=$ABACUS_PERSIST/$SIM_NAME\_products/$SIM_NAME\_rockstar_halos/$(echo "$SLICE" | awk -F '/' '{print $NF}' - | sed 's/slice/z/')
    CLI_CFG=$ZSLICE/auto-rockstar.cfg
    rm -f $CLI_CFG  # Always remove an existing config file

    echo "Starting slice $SLICE server"
    $ABACUS/Analysis/Rockstar/rockstar.py --ncpu $NTHREADS $SLICE --SO --tar-mode TAR --tar-remove-source-files &

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
if [ ]; then
    $ABACUS/python/Abacus/archive_sim.py $ABACUS_PERSIST/$SIM_DIR --nthreads=$NTHREADS --inplace --remove-source
else
  echo -e "No tar requested."
fi
echo -e "\n\n\n\n"


exit $ABACUS_EXIT
