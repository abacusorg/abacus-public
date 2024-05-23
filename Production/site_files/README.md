# site.def

Each machine we run Abacus on requires some performance
tuning for the specific architecture.  It's annoying to
have to change the main parameter file just to update
these tuning parameters, especially since they can usually
be shared across sims and do not affect the result of the
simulations.  So, we collect the tuning parameters
in this `site_files` directory (one file per machine),
and update the `site.def` symlink to point to the current
machine from the `install.py` script.  Thus, the parameter
files can just import `site.def` without knowing anything
else about the current machine/architecture.

The other place that contains "site" information, of course,
is the environment variables.  But those relate to global
directory structures and not the compute architecture,
which is what the site files focus on.
