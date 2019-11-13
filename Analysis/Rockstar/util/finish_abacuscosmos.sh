#!/usr/bin/env bash

# A "RUN_ON_SUCCESS" script called by Rockstar
# Note this only runs on one rank
# First arg: Rockstar dir
# Second arg: Output dir

#UTIL_DIR="$(dirname "$(readlink -f "$0")")"
#OUTPUT_DIR="$(dirname "$(readlink -f "$5")")"

# Do postprocessing steps:
# Identify parent/child halos
# Convert binary files to HDF5

$1/util/find_parents_binary -c $2/rockstar.cfg
$1/util/convert_to_HDF5.py $2
