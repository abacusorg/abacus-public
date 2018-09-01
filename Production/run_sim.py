#!/usr/bin/env python
"""
A top-level script to launch a sim.  Usage is:

 You may wish to give a sim its own `run_sim.py` script and run it that
 way, as is done in the sim in the "Example" directory.  This also gives
you the chance to edit environment variables from the Python script.

For suites of simulations, you probably don't need or want a script
for each individual sim.  In those cases, one can use this generic
script to launch the sim by passing the directory containing the
configuration ("info" directory or "abacus.par2" file).
"""

from Abacus import abacus
import sys

if __name__ == '__main__':
    parser = abacus.default_parser()
    parser.add_argument('config_dir', help="The directory containing the 'info/' directory or the 'abacus.par2' file")
    parser.add_argument('parfn', help="The parameter file, relative to CONFIG_DIR", nargs='?', default='abacus.par2')
    args = parser.parse_args()
    args = vars(args)
    retcode = abacus.run(**args)
    sys.exit(retcode)
