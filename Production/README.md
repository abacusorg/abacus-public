# Production

This directory contains example Abacus simulation configurations and job scripts. To get started with with your own simulation, copy the `Example` simulation, edit the `abacus.par2` configuration file, and follow the example scripts to create derivatives and ICs, run the simulation, and post-process it:

```console
$ cd Production
$ cp -r Example MyCoolSim
```

Edit `MyCoolSim/abacus.par2` with the settings for your simulation. To run the simulation interactively (i.e. without Slurm), invoke the Python driver script:

```console
$ ./run_sim.py MyCoolSim
```

This will create derivatives, ICs, and run the simulation.

To start the simulation from scratch (discarding any state and restarting from the ICs):
```console
$ ./run_sim.py MyCoolSim --clean
```

This will never delete the ICs. Use `--erase-ics` to do that.

If using Slurm, modify the `sim.sbatch` script with the appropriate parameters for your simulation. Then submit it to the queue:

```console
$ sbatch sim.sbatch
```

This is just a thin wrapper to `run_sim.py`, so Abacus will likewise create derivatives and ICs as necessary. However, for big sims, you may want these steps to be done beforehand in a non-GPU allocation. There are helper scripts for these steps too:

```console
$ sbatch derivatives.sbatch
$ sbatch ic.sbatch
```

Abacus should already be compiled, and your environment activated, before submitting any of these scripts.

### Cleaning
Most of the post-processing pipeline is fairly mature, with the exception of the cleaning. We've started to migrate the cleaning to a more modern structure in `Abacus/post`, but it's not ready yet. The cleaning can still be run, but the Slurm scripts depend on the original [abacus_mergertree](https://github.com/sownakbose/abacus_mergertree/) repo. See [Slurm Scripts](#slurm-scripts) below.

## Parameter Files
The `abacus.par2` parameter files are "second level" parameter files that get processed into `abacus.par` files before being passed to Abacus.  The `par2` files support Python math, variable substitutions, and file includes.

The parameters themselves are somewhat irregularly documented. The `include/Parameters.cpp` file contains the parameter registrations and some comments, and `doc/user_guide.pdf` also has fairly good parameter documentation. However, searching the source code for the variable corresponding to the parameter may be necessary in some cases.

## Site Files
"Site files" are parameter files that contain recommended parameter settings for a given computer/cluster, like perlmutter, summit, hal, etc. They are stored in the `Production/site_files` directory. Generally one would put generic parameters here, like performance tunings specific to a node architecture, not anything specific to a certain sim.

## Slurm Scripts
The Slurm scripts in `Example` can be used for end-to-end management of a simulation: making derivatives, making ICs, running the sim, analyzing it, compressing it, and archiving it. The scripts are designed to be short and readable. Some of them, like `sim.sbatch` don't take any arguments and use the current working directory:

```console
$ sbatch sim.sbatch
```

Others, like the post-processing scripts, take a directory as an argument:

```console
$ sbatch 01-post-groups.sbatch $SCRATCH/MyCoolSim/
```

You should open up each script before you run it to examine its usage, change any Slurm parameters, and set any configuration variables.

The post-processing scripts are numbered `01-post-groups.sbatch`, `02-post-associations.sbatch`, etc. This is because some of the steps depend on previous ones. Some steps can be run concurrently, however; one can examine the dependencies listed at the top of the scripts to determine which ones.

The scripts were written on Perlmutter, but are designed to be fairly general and adaptable to other sites. Once processed, the data products are designed to be used with abacusutils.

For processing a simulation locally without job scheduler, one might prefer to run many of these steps manually, possibly bypassing Slurm and/or disBatch. This is fine; the scripts can at least serve as an example.

## Usage of disBatch
Many of the post-processing Slurm scripts use [disBatch](https://github.com/flatironinstitute/disBatch/). disBatch runs a list of command-line tasks, dynamically distributing them to task slots within a Slurm resource allocation. This naturally helps with load-balancing balancing a set of heterogeneous tasks (e.g. cleaning late-time outputs will take much longer than early-time outputs).

The basic pattern of our disBatch Slurm jobs is to run a Python script that writes the list of disBatch tasks, and then run disBatch on that list. The Python scripts themselves live in `Abacus/post`.

disBatch also produces a status file that allows one to track task success or failure. One can re-run failed tasks from this status file with the disBatch flags `-R -r status.txt` (see the disBatch docs). If a Slurm job fails, users should carefully inspect the state of the data on disk and decide whether to resubmit the Slurm job, such that a new task file is generated, or whether to resume from the old status file. While we have tried to design the processing such that deletions of original data are conditional on successful processing, users should be extra cautious if the processing has gotten into an unexpected state, e.g. some tasks have failed.
