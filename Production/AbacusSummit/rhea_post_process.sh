#!/usr/bin/env bash

set -e  # exit on error
set -o pipefail

export SIM_SET=${SIM_SET:-AbacusSummit}
SIM_NAME=$1

##### queue up group processing #####

mkdir -p tmp
TMPFN=tmp/${SIM_NAME}.groupdirs

GROUPDIRS=($ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME/group/*/)
cat <<< "${GROUPDIRS[@]}" > $TMPFN
NDIRS=${#GROUPDIRS[@]}
NDIRSm1=$(($NDIRS - 1))

JOBID1=$(sbatch --parsable --job-name=${SIM_NAME}_PostProcessGroups --array=0-$NDIRSm1 rhea_post_process_groups.slurm $SIM_NAME)

##### queue up lightcone processing #####

TMPFN=tmp/${SIM_NAME}.lcdirs

LCDIRS=($ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME/lightcone/*/)
cat <<< "${LCDIRS[@]}" > $TMPFN
NDIRS=${#LCDIRS[@]}
NDIRSm1=$(($NDIRS - 1))

JOBID2=$(sbatch --parsable --job-name=${SIM_NAME}_PostProcessLC --array=0-$NDIRSm1 rhea_post_process_lightcones.slurm $SIM_NAME)

##### Queue up the epilogue ###
# This job will clean up and mark completion after all previous jobs have completed
# It carries the special job name that the assembly line will look for to determine if post-processing is queued
sbatch --job-name=${SIM_NAME}_PostProcessEpilogue --depend=afterok:$JOBID1:$JOBID2 --kill-on-invalid-dep=yes rhea_post_process_epilogue.slurm $SIM_NAME
#sbatch --job-name=${SIM_NAME}_PostProcessEpilogue --depend=afterok:$JOBID2 --kill-on-invalid-dep=yes rhea_post_process_epilogue.slurm $SIM_NAME