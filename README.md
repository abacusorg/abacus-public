# abacus

A high-performance N-body code for cosmological simulations.

See installation instructions in `INSTALL`, and an example of how to run a simulation in `Production/Example/`.

------

## Code directories:

* **singlestep**: the primary abacus timestepping code

* **python**: drivers for the code and related utilies

* **Test**: code for small test cases

* **Convolution**: code to compute the multipoles, Taylors, and convolutions

* **Derivatives**: code to generate derivatives used by the far force

* **ParseHeader**: library to read/write headers on each output file

* **include**: small utility functions (timing, threevectors, files, cosmology)

* **zeldovich**: the Zel'dovich initial condition generator

* **Analysis**: code for generating power spectra, running halo finders, etc

* **clibs**: libraries that codes may want to link against for various analysis tasks

* **Configuration**: files related to the autoconf build system

* **modulefiles**: modules for easy loading/unloading of abacus environmental vars

* **Production**: various configurations used for major production runs

* **util**: utilities to deal with the pack14 output format

------

## Documentation directories:

* **doc**: Code documentation

* **notes**: Musings about the algorithms, computers, etc.

* **papers**: Journal & conference papers (make subdirectories for each one)

------

## Development:

Developers are encouraged to work in a development branch, rather than a forked repository.  This encourages transparency in areas like new code interfaces and helps avoid accidental duplication of work.  Developers are encouraged to merge their work into `master` often, even if it is not entirely stable.  Stable versions will be git `tag`ged with version numbers.

As always, github issues are highly encouraged to keep track of problems and/or feature development in the code.  Pull requests are also an excellent tool to use when asking for feedback on new code.  Specifically, instead of merging directly into master, you may wish to create a pull request of your branch into master, which invites others to review your code.  This doesn't need to be done for every merge (since we want to encourage frequent merging into master) but is a good way to get collaborative feedback.
