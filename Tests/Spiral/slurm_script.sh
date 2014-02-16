#!/bin/bash
#SBATCH --get-user-env #Do we need to do this if everthing is in .bashrc?
#SBATCH -n 64  #Number of cores
#SBATCH -t 120  #Runtime in minutes
#SBATCH -p general  #Partition to submit to
#SBATCH --nodes=1 #Ensure we only run on one node
#SBATCH --mem-per-cpu=1000  #Memory per cpu in MB (see also --mem)

# Don't run this script directly; submit it to the queue with
# sbatch slurm_script.sh
# Check this directory for a file named slurm-jobid.out, which should
# contain stderr and stdout

echo "started" > /n/home04/lgarrison/abacus_tmp/started
~/abacus/Tests/Spiral/SpiralTest.py

