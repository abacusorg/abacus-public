#!/bin/bash

# Fail if any command returns a non-zero exit status
set -e

export OMP_NUM_THREADS=24
export ABACUS_PERSIST="/mnt/gosling2/bigsims/emulator_1100box_planck/"
export SIM_NAME="emulator_1100box_planck_spline_smallsoft_00"

# Split SIM_DIR into name and project
echo -e "* Checking if we need to run Rockstar:\n"
if [ 0 ]; then
  echo -e "Running Rockstar."
  for SLICE in $ABACUS_PERSIST/$SIM_NAME/slice*; do
    ZSLICE=$ABACUS_PERSIST/$SIM_NAME\_products/$SIM_NAME\_rockstar_halos/$(echo "$SLICE" | awk -F '/' '{print $NF}' - | sed 's/slice/z/')
    CLI_CFG=$ZSLICE/auto-rockstar.cfg
    rm -f $CLI_CFG  # Always remove an existing config file

    echo "Starting slice $SLICE server"
    $ABACUS/Analysis/Rockstar/rockstar.py --ncpu $OMP_NUM_THREADS $SLICE --SO &  #--tar-mode TAR --tar-remove-source-files &

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
    $ABACUS/python/Abacus/archive_sim.py $ABACUS_PERSIST/$SIM_DIR --nthreads=$OMP_NUM_THREADS --inplace --remove-source
else
  echo -e "No tar requested."
fi
echo -e "\n\n\n\n"


exit $ABACUS_EXIT
