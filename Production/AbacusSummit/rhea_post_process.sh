#!/usr/bin/env bash

set -e  # exit on error
set -o pipefail

export SIM_SET=${SIM_SET:-AbacusSummit}
SIM_NAME=$1

mkdir -p tmp
TMPFN=tmp/${SIM_NAME}.groupdirs

GROUPDIRS=($ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME/group/*/)
cat <<< "${GROUPDIRS[@]}" > $TMPFN
NDIRS=${#GROUPDIRS[@]}
NDIRSm1=$(($NDIRS - 1))

# TODO: slurm dependency is pretty ugly, requires an extra script.
# Do we want to save the job id and query the scheduler? How long can we query old jobs?
JOBID=$(sbatch --parsable --job-name=${SIM_NAME}_PostProcess --array=0-$NDIRSm1 rhea_post_process.slurm $SIM_NAME)
sbatch --depend=afterok:$JOBID rhea_post_process_epilogue.slurm $SIM_NAME
