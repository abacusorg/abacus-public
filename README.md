# abacus
======

A high-performance N-body code for cosmological simulations.

See installation instructions in INSTALL.

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
