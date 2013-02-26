abacus
======

abacus

------

Code directories:

* singlestep -- The primary abacus timestepping code

* python -- Drivers for the code

* Test -- Code for small test cases.

* Direct -- Code for the direct forces

* Convolution -- Code to compute the multipoles, Taylors, and convolutions

* Derivatives -- Code to generate derivatives

* Multipoles -- Code to generate the ASM for the multipoles & Taylors

* Parseheader -- For the ParseHeader utility

* include -- Small utility functions (timing, FFT, threevectors, files)

* zeldovich -- The Zel'dovich initial condition generator

------

Documentation directories:

* doc -- Code documentation

* notes -- Musings about the algorithms, computers, etc.

* papers -- Journal & conference papers (make subdirectories for each one)

------

Installation:

# Get the $ABACUS shell variable set
cd [top-level directory]
export ABACUS=`pwd`

# The code is set up for gcc, except for python/clibs, which is icc.
cd Derivatives;   make;   
./CreateDerivatives 15 16 2 8
cd ..
cd Multipoles;    make;   cd ..
cd ParseHeader;   make;   cd ..
cd singlestep;    make;   cd ..
cd python/clibs;  make;   cd ../..

# Now can run the Ewald test
cd Tests/Ewald
./runewaldtest.py





