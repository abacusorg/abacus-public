#!/usr/bin/env python3

"""
For a given PPD, calculates all the CPD that have small prime factors
and will stay within a reasonable number of PPC.

Usage: ./choose_cpd.py --help

"""

import sys
import numpy as np
from collections import OrderedDict
import itertools as it
import argparse
from functools import reduce

# Generates variable length combinations, where each factor i can be repeated max_repeats[i] times.
def comb(factors, max_repeats):
    N_per_factor_gen = it.product(*[list(range(0,mr+1)) for mr in max_repeats])
    for N_per_factor in N_per_factor_gen:
        yield sum(([f,]*N for f,N in zip(factors, N_per_factor)), [])

        
# Prints CPD values that satisfy the specified constraints
# If using_fftw, trims combinations that would be slow in FFTW.  If max_prime > 13, this is disabled.
def choose_cpd(ppd, max_prime_factor=13, min_ppc=10, max_ppc=200, show_slow=False):
    primes = np.array([3,5,7,11,13,17,19])
    assert max_prime_factor <= primes.max()
    primes = [p for p in primes if p <= max_prime_factor]
    
    max_cpd = np.ceil(ppd/min_ppc**(1./3))
    
    N_max = OrderedDict()  # The max number of times each factor can be used in max_cpd
    
    for factor in primes:
        N_max[factor] = int(round(np.log(max_cpd)/np.log(factor)))
        
    options = []
    for factors in comb(primes, list(N_max.values())):
        if not show_slow:
            if factors.count(11) + factors.count(13) > 1:
                continue
        cpd = reduce(lambda x,y: x*y, factors,1)
        ppc = (float(ppd)/cpd)**3.
        if min_ppc <= ppc <= max_ppc:
            s = 'CPD: {:4d} = {:20s} [ppc = {:7.2f}]'.format(cpd, '*'.join([str(f) for f in factors]), ppc)
            if factors.count(11) + factors.count(13) > 1:
                s += ' [Slow under FFTW!]'
            options += [(cpd,s)]
    
    print('For max_prime_factor = {:d}, min_ppc = {:d}, max_ppc = {:d}:'.format(max_prime_factor, min_ppc, max_ppc))
    for cpd,s in sorted(options, key=lambda x:x[0]):
        print(s)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('PPD', help='Particle-per-dimension of the simulation', type=int)
    parser.add_argument('--max-prime', help='Maximum prime factor to consider', type=int, default=13)
    parser.add_argument('--show-slow', help='Show all possibilities, even those that would be slow under FFTW and are hidden by default', action='store_true')
    parser.add_argument('--max-ppc', help='Maximum particles per cell', type=int, default=400)
    parser.add_argument('--min-ppc', help='Minimum particles per cell', type=int, default=20)
    
    args = parser.parse_args()
        
    choose_cpd(args.PPD, max_prime_factor=args.max_prime, show_slow=args.show_slow, min_ppc=args.min_ppc, max_ppc=args.max_ppc)
