abacus
======

abacus

------

Code directories:

* abacus -- The primary abacus timestepping and convolve code

* derivatives -- Code to generate derivatives

* multipoles -- Code to generate the ASM for the multipoles & Taylors

* parseheader -- For the ParseHeader utility

* setup -- Code to turn external particle sets into legal initial conditions.

* spiral -- Code to generate spiral initial conditions, run the
tests, and analyze outputs.

* test -- Code for other small test cases.

* util -- Small utility functions (timing, FFT, threevectors, files)

* zeldovich -- The Zel'dovich initial condition generator

* bin -- Where executables should go

Do we need a separate include directory?  Perhaps not: we're tending
to build the code in one inline format anyways, so how are includes
different than source files?  Plus we have the util directory for little
snippets.

------

Documentation directories:

* doc -- Code documentation

* notes -- Musings about the algorithms, computers, etc.

* papers -- Journal & conference papers (make subdirectories for each one)
