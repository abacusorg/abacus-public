#!/bin/bash

# Run the simulation.
# Depends on: derivatives

#SBATCH -N 1
#SBATCH --gpus-per-node=4
#SBATCH --exclusive
#SBATCH -t 3-0
##SBATCH -p gpuxl
#SBATCH -p gpupreempt -q gpupreempt
#SBATCH --requeue
#SBATCH -C h100&genoa
#SBATCH --open-mode=append
#SBATCH --comment=FI_JOB_RESOURCES=tmpfs
#SBATCH -x workergpu204

#########

set -eu

cosm=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
sim=AbacusBacklight_base_c${cosm}_ph000

# TODO: checkpoint more often if using the preempt queue, less often if not

srun $ABACUS/Production/Backlight/run_sim.py ./param/${sim}.par2
abacus_exit=$?

echo "Abacus exited with code $abacus_exit"

exit $abacus_exit
