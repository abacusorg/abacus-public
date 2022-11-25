'''
Driver functions to invoke the `CreateDerivatives` executable to precompute
derivatives for a given CPD and Order.  Also handles creating derivatives
suitable for the 2D code via reflection of the original files.

This module can be invoked as an executable via:
$ python -m Abacus.derivatives --help
'''

import multiprocessing
import os
import shlex
import shutil
from pathlib import Path

import numba
import numpy as np
import tqdm

from . import convert_derivs_to_float32
from .abacus import abacuspath
from .Tools import call_subprocess

abacuspath = Path(abacuspath)

maxthreads = len(os.sched_getaffinity(0))


def make_derivatives(param, floatprec=False, twoD=False,
                        nworker=None, nthread=None, _stage=True):
    '''
    Main entry point to create derivatives (1D or 2D).

    Look for the derivatives required for param and make them if they don't
    exist. If derivatives for the 2D code are required, will first create the
    1D derivatives, and then reflect them into 2D.

    Derivatives are copied from DerivativesSourceDirectory to
    DerivativesDirectory if necessary.

    DerivativesSourceDirectory must be globally accessible in the parallel
    code.

    Parameters
    ----------
    param : dict
        The parameters dictionary

    floatprec : bool (optional)
        Make the derivatives available in float32 precision.
        Default: False

    twoD : bool (optional)
        Make the derivs available in the format for the 2D code.
        Default: False
    
    nworker : int or None (optional)
        The number of processes/threads to use for IO.
        Default: param['PyIOConcurrency'] or 1

    nthread : int or None (optional)
        The number of threads to use for CPU work.
        Default: based on affinity mask

    _stage : bool (optional)
        Copy the result to DerivativesDirectory, which may be node-local.
        This is used internally.
        Default: True.

    Returns
    -------
    source_paths: list of pathlib.Path
        The paths to the derivative files in the global dir
    '''

    source_dir = Path(param['DerivativesSourceDirectory'])
    derivs_dir = Path(param['DerivativesDirectory'])
    if twoD:
        source_dir /= '2D'
        derivs_dir /= '2D'
    source_dir.mkdir(parents=True, exist_ok=True)
    derivs_dir.mkdir(parents=True, exist_ok=True)

    CPD = param['CPD']
    order = param['Order']
    NFR = param['NearFieldRadius']
    DER = param['DerivativeExpansionRadius']

    stem32 = '_float32' if floatprec else ''
    stem2D = '_2D' if twoD else ''
    zmax = CPD if twoD else (CPD//2 + 1)
    dprefix = f"fourierspace{stem32}{stem2D}_{CPD}_{order}_{NFR}_{DER}"
    dnames = [ f"{dprefix}_{i}" for i in range(zmax) ]
    dpaths = [ derivs_dir / dn for dn in dnames ]
    dpaths_source = [ source_dir / dn for dn in dnames ]

    have_derivs = all( dpath.is_file() for dpath in dpaths )
    have_derivs_in_source = all( p.is_file() for p in dpaths_source )

    nworker = nworker or param.get('PyIOConcurrency',1)
    nthread = nthread or maxthreads

    if not have_derivs and not have_derivs_in_source:
        # Going to have to make derivs
        print(f'Could not find derivatives in DerivativesDirectory="{derivs_dir}"\n'
                f'\tor DerivativesSourceDirectory="{source_dir}". Creating them...')

        if twoD:
            # Need to make 2D. First ensure 1D exists.
            source_paths_1d = make_derivatives(param, twoD=False,
                floatprec=floatprec, _stage=False)
            make_2d(source_paths_1d[0].parent, param, floatprec=floatprec, nworker=nworker, nthread=nthread)
            have_derivs_in_source = True

        elif floatprec:
            # first make sure the derivatives exist in double format
            source_paths_64 = make_derivatives(param, twoD=False, floatprec=False, _stage=False)
            # now make a 32-bit copy
            for p in tqdm.tqdm(source_paths_64, desc='derivs 64-bit to 32-bit', unit='file'):
                convert_derivs_to_float32.convert(p)
            have_derivs_in_source = True

        else:
            # Call the exe to make 1D, double-precision derivs

            create_derivs_cmd = [str(abacuspath / "Derivatives" / "CreateDerivatives"),
                    str(CPD), str(order), str(NFR), str(DER)]

            (derivs_dir / 'farderivatives').unlink(missing_ok=True)
            
            for a in list(range(1,9)) + [16]:
                ADfn = f'AD32_{a:03d}.dat'
                source_ADfn = abacuspath / "Derivatives" / ADfn
                if not (derivs_dir / ADfn).is_file():
                    shutil.copy(source_ADfn, derivs_dir)
                
            call_subprocess(create_derivs_cmd, cwd=derivs_dir)
            assert all( dpath.is_file() for dpath in dpaths )

            # could create in either derivs_dir or source_dir,
            # but maybe derivs_dir is faster because we do many appends
            have_derivs = True

    if derivs_dir.samefile(source_dir):
        have_derivs = True
        have_derivs_in_source = True

    if not have_derivs_in_source:
        print(f'Staging out derivatives from local to global')
        # if we have local, stage out to global
        assert have_derivs
        for dn in dnames:
            shutil.copy(derivs_dir / dn, source_dir / dn)

    if _stage:
        stage_derivs(source_dir, derivs_dir, dnames, param, have_derivs)

    return dpaths_source


def stage_derivs(source_dir, derivs_dir, dnames, params, have_derivs):
    '''Copy the derivatives from global disk to local storage.

    This will just be an ordinary copy in the serial code.

    In the parallel code, this will invoke mpirun to have each node copy to
    itself.  Each node will copy all derivs.  In a larger refactoring, this
    could be improved (TODO).

    Parameters
    ----------
    source_dir : pathlib.Path
        The source directory. In the parallel code, all nodes must see this
        directory.

    derivs_dir : pathlib.Path
        Destination directory. May be node-local.

    dnames : iterable of path-like
        File names to copy

    params: dict
        Simulation parameters

    have_derivs: bool
        Whether derivs_dir has the derivatives
    '''

    #CPD = params['CPD']
    parallel = params.get('Parallel', False)
    #twoD = params.get('NumZRanks',1) > 1

    # Note the serial, have_derivs case is a no-op

    if parallel:
        if source_dir.samefile(derivs_dir):
            # probably pure global
            return

        mpirun_cmd = shlex.split(params['mpirun_cmd'])

        # stage in to local
        mkdir_cmd = mpirun_cmd + ['mkdir', '-p', str(derivs_dir)]
        cp_prefix = mpirun_cmd + ['cp', '-t', str(derivs_dir)]
        spaths = [source_dir / dn for dn in dnames]
        cp_cmd = cp_prefix + spaths
        log_cmd = " ".join(cp_prefix) + f' {spaths[0]} ...'
        print(f'Staging derivatives to nodes with command "{log_cmd}"')
        call_subprocess(mkdir_cmd)
        call_subprocess(cp_cmd)
    elif not have_derivs:
        print(f'Copying derivatives from "{source_dir}" to "{derivs_dir}"')
        for dn in dnames:
            shutil.copy(source_dir / dn, derivs_dir / dn)
        

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


def make_2d(derivs_dir, param, floatprec=True, nworker=1, nthread=1):
    assert nthread >= nworker
    
    #nthread_per_worker = np.diff(nthread*np.arange(nworker+1)//nworker)
    nthread_per_worker = nthread//nworker

    cpd = param['CPD']
    order = param['Order']

    rml = (order+1)**2
    halfquadxy = (cpd+1)*(cpd+3)//8
    cpdp1half = (cpd+1)//2

    rmap = rmapper(order)

    derivs_dir = Path(derivs_dir)
    twoDdir = derivs_dir / '2D'
    twoDdir.mkdir(exist_ok=True)

    f32 = '_float32' if floatprec else ''
    nfr = param['NearFieldRadius']
    der = param['DerivativeExpansionRadius']

    state = dict(inprefix = derivs_dir / f'fourierspace{f32}_{cpd}_{order}_{nfr}_{der}',
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
