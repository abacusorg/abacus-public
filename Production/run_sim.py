#!/usr/bin/env python
"""
An example simulation script.

Usage
-----
Make a copy of this directory:
$ cp -r Example/ MySim/

Edit the sim parameters in MySim/abacus.par2.
At least change SimName!

Look at the command line options:
$ ./run_sim.py --help

Launch the sim:
$ ./run_sim.py

"""

from Abacus import abacus
import sys

if __name__ == '__main__':
    parser = abacus.default_parser()
    args = parser.parse_args()
    retcode = abacus.run('abacus.par2', **vars(args))
    
    sys.exit(retcode)
