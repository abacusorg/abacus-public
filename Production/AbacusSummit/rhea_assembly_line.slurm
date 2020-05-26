#!/usr/bin/env bash

## This script runs on rhea and calls AbacusSummit_assembly_line.py once. 
## It requeues itself, with a one hour delay, when done. 


#SBATCH -A ast145
#SBATCH -J AbacusSummit_assembly_line
#SBATCH -N 1
#SBATCH -t 1:00:00 #run for X hours, just in case. Should only take a minute though. 
#SBATCH -o logs/assembly_line-%j.out

set -e  # exit on error

python $ABACUS/Production/scripts/AbacusSummit_assembly_line.py

echo "Requeuing this job, w/ 1 hr delay."

export jobid=<JOB ID>
scontrol requeue $jobid 
scontrol update JobID=$jobid StartTime=now+1hour

