#!/bin/bash

# Fail if any command returns a non-zero exit status
set -e

export OMP_NUM_THREADS=24
export ABACUS_PERSIST="/mnt/gosling2/bigsims/AbacusCosmos_720box_planck/"
export SIM_NAMES=$(echo AbacusCosmos_720box_planck_00-{4,5})

# Split SIM_DIR into name and project
echo -e "* Checking if we need to run Rockstar:\n"
if [ 0 ]; then
  echo -e "Running Rockstar."
  for SIM_NAME in $SIM_NAMES; do 
      for SLICE in $ABACUS_PERSIST/$SIM_NAME/slice*; do
        ZSLICE=$ABACUS_PERSIST/$SIM_NAME\_products/$SIM_NAME\_rockstar_halos/$(echo "$SLICE" | awk -F '/' '{print $NF}' - | sed 's/slice/z/')
        CLI_CFG=$ZSLICE/auto-rockstar.cfg
        rm -f $CLI_CFG  # Always remove an existing config file

        echo "Starting slice $SLICE server"
        $ABACUS/Analysis/Rockstar/rockstar.py --ncpu $OMP_NUM_THREADS $SLICE --SO --tar-mode TAR --tar-remove-source-files &

        # Wait for the server to generate the client config file
        while [ ! -f $CLI_CFG ]; do
            sleep 1
        done

        echo "Starting client in $ZSLICE"
        mpirun -np $OMP_NUM_THREADS $ABACUS/Analysis/Rockstar/rockstar -c $CLI_CFG
      done
  done
else
  echo -e "No Rockstar requested."
fi
echo -e "\n\n\n\n"


echo -e "* Checking if we need to tarball the outputs:\n"
if [ ]; then
    for SIM_NAME in $SIM_NAMES; do
        $ABACUS/python/Abacus/archive_sim.py $ABACUS_PERSIST/$SIM_NAME --nthreads=$OMP_NUM_THREADS --inplace --remove-source
    done
else
  echo -e "No tar requested."
fi
echo -e "\n\n\n\n"


exit $ABACUS_EXIT
