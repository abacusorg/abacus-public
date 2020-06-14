# Overview
This directory contains the "assembly line" scripts used to run the AbacusSummit simulations.

The `logs` directory is used to store logs from LSF and SLURM job scripts running on summit and rhea.

The AbacusSummit scripts used to live in `Production/scripts/`, but it seems best to aggregate AbacusSummit-related materials in one directory,
as these are fairly specialized just for AbacusSummit.

## Post-processing

Abacus outputs several data products as it runs.  The primary ones are halo catalogs (with particle subsamples)
and lightcones.  The halo catalogs are written per-slab and need to be concatenated in chunks of \~50 for efficient
data movement.  The light cones are written per-node, but are relatively small to begin with and thus also need to
be concatenated across timestep.  Both halo catalogs and light cones also need to undergo compression.

The data volumes demand that this post-processing be parallelized so we can clear the disk quickly for new sims.  The
following is a high-level summary of the parallelization.

### Lightcones
For the lightcones, we have several thousand independent tasks: (# LC, 1 to 3) x (# file types = 3 [rv,heal,pid]) x (# time step \~ 1000).
For each task we're concatenating anywhere from a few to 60 node files.  Note that not all time steps have all light cones.

Aside from the heterogeneity of the individual tasks, each of the obvious axes over which to parallelize the problem are also
heterogeneous.  We'd like to just launch a giant job array and be done with it, but SLURM on rhea only support job arrays up
to 100 tasks, not 10000.

Flatiron's [disBatch](https://github.com/flatironinstitute/disBatch) package offers a nice solution, in which we specify a list of independent jobs that are dynamically dispatched
to any number of nodes.  It queries SLURM to figure out how many slots we want on each node, then fills slots with jobs from
the list, replacing them with new jobs as they finish.  This immediately resolves the heterogeneity concern, and automatically
scales up and down to any number of nodes.

Concerning the checksums, with the lightcones, each concatenated file gets its own checksum file.  We issue a disBatch BARRIER
to synchronize the nodes after all concatenation tasks, then farm out a final, lightweight job in each directory that removes
the raw data checksum file and concatenates the new checksum files into a single `checksums.crc32`.

### Halo catalogs
With halo catalogs, there's a little more preamble work for each catalog directory.  It makes the chunks a little less than
embarrassingly parallel.  So we have one Python script per catalog, called `convert_raw_groups_to_asdf.py`.  This script
uses Python's multiprocessing to launch several workers

We're also changing the


TODO: time slices

### Compression
Within each compression task, we run two blosc compression threads.  That lets us run ~8 tasks per node, which should saturate
the IO to each node while still using nodes efficiently in the non-blosc parts.
