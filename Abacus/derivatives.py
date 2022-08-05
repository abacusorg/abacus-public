'''
Driver functions to invoke the `CreateDerivatives` executable to precompute
derivatives for a given CPD and Order.  Also handles creating derivatives
suitable for the 2D code via reflection of the original files.

This module can be invoked as an executable via:
$ python -m Abacus.derivatives --help
'''

from pathlib import Path
import shutil
import os
import multiprocessing

import numpy as np
import numba
import tqdm

from .Tools import call_subprocess
from .abacus import abacuspath
from . import convert_derivs_to_float32
abacuspath = Path(abacuspath)

maxthreads = len(os.sched_getaffinity(0))


def make_derivatives(param, search_dirs=True, floatprec=False, twoD=False,
                        nworker=None, nthread=None):
    '''
    Main entry point to create derivatives (1D or 2D).

    Look for the derivatives required for param and make them if they don't exist.
    If derivatives for the 2D code are required, will first create the 1D
    derivatives, and then reflect them into 2D.

    Parameters
    ----------
    param : dict
        The parameters dictionary

    search_dirs : bool, path-like, or iterable of path-like (optional)
        The directories to search. `True` (default) searches canonical dirs
        based on environment variables.

    floatprec : bool (optional)
        Make the derivatives available in float32 precision.
        Default: False

    twoD : bool (optional)
        Make the derivs available in the format for the 2D code.
    
    nworker : int or None (optional)
        The number of processes/threads to use for IO

    nthread : int or None (optional)
        The number of threads to use for CPU work
    '''

    # We will attempt to copy derivatives from the archive dir if they aren't found
    if search_dirs == True:
        search_dirs = []
        for var in ('ABACUS_PERSIST', 'ABACUS_SSD', 'ABACUS_DERIVS'):
            if var in os.environ:
                sdir_name = 'Derivatives/2D' if twoD else 'Derivatives'
                search_dirs += [Path(os.environ[var]) / sdir_name]

    if type(search_dirs) is str:
        search_dirs = [search_dirs]
    search_dirs = [Path(d) for d in search_dirs]

    derivsdir = Path(param['DerivativesDirectory'])
    if twoD:
        derivsdir /= '2D'
    derivsdir.mkdir(parents=True, exist_ok=True)

    search_dirs = list(filter(lambda d: not (d.exists() and derivsdir.samefile(d)), search_dirs))

    CPD = param['CPD']
    order = param['Order']
    NFR = param['NearFieldRadius']
    DER = param['DerivativeExpansionRadius']

    stem32 = '_float32' if floatprec else ''
    stem2D = '_2D' if twoD else ''
    zmax = CPD if twoD else (CPD//2 + 1)
    dnames = [ f"fourierspace{stem32}{stem2D}_{CPD}_{order}_{NFR}_{DER}_{i}"
                for i in range(zmax) ]
    dpaths = [ derivsdir / dn for dn in dnames ]

    nworker = nworker or param.get('PyIOConcurrency',1)
    nthread = nthread or maxthreads

    # First check if the derivatives are in the DerivativesDirectory
    if not all( dpath.is_file() for dpath in dpaths ):
        # If not, maybe they exist elsewhere
        for sdir in search_dirs:
            if all( (sdir / dn).is_file() for dn in dnames ):
                print(f'Found derivatives in "{sdir}". Copying to "{derivsdir}"')
                for dn in dnames:
                    shutil.copy(sdir / dn, derivsdir / dn)
                break
        else:
            print(f'Could not find derivatives in "{derivsdir}" or search dirs "{search_dirs}". Creating them...')

            if twoD:
                # Need to make 2D. First ensure 1D exists.
                make_derivatives(param, twoD=False, floatprec=floatprec)
                make_2d(param, floatprec=floatprec, nworker=nworker, nthread=nthread)

            elif floatprec:
                # first make sure the derivatives exist in double format
                dpaths64 = make_derivatives(param, twoD=False, floatprec=False)
                # now make a 32-bit copy
                for dpath64 in tqdm.tqdm(dpaths64, desc='derivs 64-bit to 32-bit', unit='file'):
                    convert_derivs_to_float32.convert(dpath64)

            else:
                # Call the exe to make 1D, double-precision derivs

                create_derivs_cmd = [str(abacuspath / "Derivatives" / "CreateDerivatives"),
                        str(CPD), str(order), str(NFR), str(DER)]

                (derivsdir / 'farderivatives').unlink(missing_ok=True)
                
                for a in list(range(1,9)) + [16]:
                    ADfn = f'AD32_{a:03d}.dat'
                    source_ADfn = abacuspath / "Derivatives" / ADfn
                    if not (derivsdir / ADfn).is_file():
                        shutil.copy(source_ADfn, derivsdir)
                    
                call_subprocess(create_derivs_cmd, cwd=derivsdir)
                assert all( dpath.is_file() for dpath in dpaths )

    return dpaths


# Helpers to create the 2D derivs

@numba.experimental.jitclass([('order', numba.int64),
                              ('rml', numba.int64),
                              ('_rmap', numba.int64[:]),
                              ])
class rmapper:
    def __init__(self, order):
        self.order = order
        self.rml = (order+1)**2
        self._rmap = np.empty(2*(order+1)**2, dtype=np.int64)

        m = 0
        for c in range(2):
            for a in range(order-c+1):
                for b in range(order-c-a+1):
                    self._rmap[self._lmap(a,b,c)] = m
                    m += 1
        #assert m == self.rml

    def get(self, a, b, c):
        return self._rmap[self._lmap(a,b,c)]

    def _lmap(self, a, b, c):
        return 2*(self.order+1)*a + 2*b + c


@numba.njit(nogil=True)
def reflect_halfquad(d, rmap):
    '''reflect the half-quadrant to the full quadrant'''
    rml = d.shape[0]
    order = int(round(np.sqrt(rml)))-1
    cpdp1half = int(-1 + np.sqrt(1 + 8*d.shape[-1])) // 2
    f = np.empty((rml, cpdp1half, cpdp1half), dtype=d.dtype)

    m = 0
    for c in range(2):
        for a in range(order-c+1):
            for b in range(order-c-a+1):
                i = 0
                for kx in range(cpdp1half):
                    for ky in range(kx+1):
                        #assert rmap(a,b,c) == m
                        f[m,kx,ky] = d[m,i]
                        mba = rmap.get(b,a,c)
                        f[m,ky,kx] = d[mba,i]
                        i += 1
                m += 1
    #assert m == rml
    return f


@numba.njit(nogil=True)
def reflect_z(fq):
    '''make a new array with the Hermitian reflection of z'''
    rml, cpdp1half, _ = fq.shape
    order = int(round(np.sqrt(rml)))-1
    newz = np.empty_like(fq)

    m = 0
    for c in range(2):
        for a in range(order-c+1):
            for b in range(order-c-a+1):
                # The odd parity elements are the imaginary ones
                if c & 1:
                    newz[m] = -fq[m]
                else:
                    newz[m] = fq[m]
                m += 1
    #assert m == rml
    return newz


def make_2d(param, floatprec=True, nworker=1, nthread=1):
    assert nthread >= nworker
    
    #nthread_per_worker = np.diff(nthread*np.arange(nworker+1)//nworker)
    nthread_per_worker = nthread//nworker

    cpd = param['CPD']
    order = param['Order']

    rml = (order+1)**2
    halfquadxy = (cpd+1)*(cpd+3)//8
    cpdp1half = (cpd+1)//2

    rmap = rmapper(order)

    derivsdir = Path(param['DerivativesDirectory'])
    twoDdir = derivsdir / '2D'
    twoDdir.mkdir(exist_ok=True)

    f32 = '_float32' if floatprec else ''
    nfr = param['NearFieldRadius']
    der = param['DerivativeExpansionRadius']

    state = dict(inprefix = derivsdir / f'fourierspace{f32}_{cpd}_{order}_{nfr}_{der}',
                 outprefix = twoDdir / f'fourierspace{f32}_2D_{cpd}_{order}_{nfr}_{der}',
                 rml=rml,
                 halfquadxy=halfquadxy,
                 cpdp1half=cpdp1half,
                 rmap=rmap,
                 cpd=cpd,
                 nthread=nthread_per_worker,
                 )

    print(f'Launching {nworker} workers with {nthread_per_worker} threads each')
    with multiprocessing.pool.ThreadPool(nworker) as pool:
        list(tqdm.tqdm(
                pool.imap(lambda z: _make_2d_worker(state, z), range(cpdp1half)),
                desc='derivs 1D to 2D', unit='file', total=cpdp1half,
            ))

def _make_2d_worker(state, z):
    inprefix = state['inprefix']
    outprefix = state['outprefix']
    rml = state['rml']
    halfquadxy = state['halfquadxy']
    cpd = state['cpd']
    rmap = state['rmap']
    nthread = state['nthread']

    fn = inprefix.parent / (inprefix.name + f'_{z}')

    numba.set_num_threads(nthread)

    d = np.fromfile(fn, dtype=np.float32)
    d = d.reshape(rml,halfquadxy)

    fullquad = reflect_halfquad(d, rmap)  # [m, kx, ky]
    zrefl = reflect_z(fullquad)

    fullquad = np.moveaxis(fullquad, (0,1,2),(1,2,0))  # [ky, m, kx]
    zrefl = np.moveaxis(zrefl, (0,1,2),(1,2,0))

    newz = cpd - z
    
    fullquad.tofile(outprefix.parent / (outprefix.name + f'_{z}'))
    if z > 0:
        zrefl.tofile(outprefix.parent / (outprefix.name + f'_{newz}'))


if __name__ == '__main__':
    import argparse

    from .InputFile import InputFile

    parser = argparse.ArgumentParser(description='Create derivatives')
    parser.add_argument('param', help='abacus.par parameter file')
    parser.add_argument('-cpd', help='override CPD', type=int)
    parser.add_argument('-2d', help='2D derivs', action='store_true')
    parser.add_argument('-out', help='DerivativesDirectory')
    parser.add_argument('-nworker', help='Number of IO workers', type=int)
    parser.add_argument('-nthread', help='Number of CPU threads, to be divided among workers', type=int)
    args = vars(parser.parse_args())

    param = dict(InputFile(args.pop('param')))

    if 'cpd' in args:
        param['CPD'] = args.pop('cpd')
    if '2d' in args:
        twoD = args.pop('2d')
    else:
        twoD = param.get('NumZRanks',1) > 1
    if 'out' in args:
        param['DerivativesDirectory'] = args.pop('out')

    make_derivatives(param, **args, floatprec=True, twoD=twoD)
