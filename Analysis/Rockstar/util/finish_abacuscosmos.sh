#!/usr/bin/env bash

# A "RUN_ON_SUCCESS" script called by Rockstar
# First arg: Rockstar dir
# Second arg: Output dir

# Do postprocessing steps:
# Identify parent/child halos
# Convert binary files to HDF5

$1/util/find_parents_binary -c $2/rockstar.cfg
$1/util/convert_to_HDF5.py $2
