#!/usr/bin/env python3
'''
Tests that abacus holds the forces on a homogeneous grid to zero.

Slabs are analyzed one at a time to support large problems.

Usage
-----
$ ./check_forces.py --help
'''

import os
import shutil
import sys
from glob import glob
import argparse
from pathlib import Path

import numpy as np
import matplotlib as mpl
from matplotlib import ticker
#mpl.use('Agg')  # no interactive plotting
import matplotlib.pyplot as plt
import asdf

from Abacus import GenParam
from Abacus import InputFile
from Abacus import Tools
    
DEFAULT_ORDER = 8
DEFAULT_PPD = 64
DEFAULT_CPD = 17

def plt_log_hist(ax, a, label=r'$|a-b|/\sqrt{|a|^2 + |b|^2}$'):
    assert (a >= 0).all()
    a_nonz = a[a > 0]
    bins = np.geomspace(a_nonz.min(),a_nonz.max(),100)
    bins = [0.,] + list(bins)  # make the first bin include 0
    ax.hist(a, bins=bins, log=True)
    ax.set_xscale('symlog', linthresh=a_nonz.min(), linscale=0.01)
    ax.set_xlim(xmin=0.)
    #plt.xscale('log')
    if label:
        ax.set_xlabel(label)

    
def plt_err_scale(disp, diff, logy=True, alpha=.2 ):
    mag = np.sqrt((disp**2).sum(axis=-1))
    if logy:
        plt.loglog(mag, diff, '.', markersize = 1., alpha = alpha)
    else:
        plt.semilogx(mag, diff, '.', markersize = 1., alpha = alpha)
    plt.xlabel('Vector magnitude')
    plt.ylabel('Error magnitude')
    plt.tight_layout()


def run_orders(orders, run_kwargs=None, save=None):
    if run_kwargs is None:
        run_kwargs = {}
    
    all_results = []
    for order in orders:
        res = run(order=order, plot=False, **run_kwargs)
        all_results += [res]
        
    if save:
        print(f'Saving results to {save}')
        af = asdf.AsdfFile(dict(results=all_results))
        af.write_to(save)
        
    plot_orders(all_results)
    
    
def plot_orders(all_results, pltfn='lattice_all_orders.pdf'):
    orders = [r['param']['Order'] for r in all_results]
    
    # Descriptive text in corner
    param = all_results[0]['param']  # use params from one result
    NP = param['NP']
    ppd = int(round(NP**(1/3)))
    ppc = NP/param["CPD"]**3
    text = ['Abacus Lattice Test',
            f'NP ${ppd:d}^3$, CPD {param["CPD"]}',
            'Single precision'
           ]
    text = '\n'.join(text)
        
    fig, ax = plt.subplots(figsize=(3.5,2.6))
    ax.set_xlabel('Multipole order $p$')
    ax.set_ylabel('Force error')
    
    all_orders = [res['param']['Order'] for res in all_results]
    all_maxes = [res['Max'] for res in all_results]
    all_rms = [res['RMS'] for res in all_results]
        
    ax.set_yscale('log')
    ax.plot(all_orders, all_maxes, marker='o', c='k', ls='-', label='Max')
    ax.plot(all_orders, all_rms, marker='o', c='k', ls='--', label='RMS')
    
    from matplotlib.lines import Line2D
    legend_lines = [Line2D([0], [0], color='k', ls='-'),
                    Line2D([0], [0], color='k', ls='--'),]
    ax.legend(legend_lines, ['Max', 'RMS'], loc='lower left')
    
    fig.text(.96, .95, text, transform=ax.transAxes, ha='right', va='top',
                bbox=dict(fc='w', alpha=0.5, ec='k'))
        
    ax.xaxis.set_ticks(all_orders[::2])

    fig.tight_layout()
    fig.savefig(pltfn, bbox_inches='tight')


def run(ppd=None, cpd=None, order=DEFAULT_ORDER, dtype=np.float32, force_output_debug=False, no_analyze=False, plot=True, power=False, **config):
    from Abacus import abacus

    if ppd is not None:
        config['NP'] = int(ppd)**3

    if cpd is not None:
        config['CPD'] = int(cpd)
    
    box = 2000/6912*config['NP']**(1/3)
    
    paramfn = "abacus.par"
    cwd = os.getcwd()
    GenParam.makeInput(paramfn, "checkforces.par2",
        ForceOutputDebug=int(force_output_debug),
        StoreForces=int(not force_output_debug),
        Order=order, BoxSize=box,
        strict=False, **config)
    params = InputFile.InputFile(paramfn)

    print('Running CheckForces with PPD {:d}, CPD {:d}'.format(int(round(params['NP']**(1/3))), params['CPD']))
    
    abacus.run(paramfn, maxsteps=2, erase_ic=True, clean=True)
    
    if not no_analyze:
        if force_output_debug:
            diff, fmag, results = analyze_forceoutputdebug(params, dtype)
            os.chdir(cwd)
            if plot:
                plot_forceoutputdebug(diff, fmag)
        else:
            fmag, results = analyze_storeforces(params, dtype)
            os.chdir(cwd)
            if plot:
                plot_storeforces(fmag)
                if not params['Parallel']:
                    plot_cell(params, dtype)
                if power:
                    plot_force_error_pk(params, dtype)
        return results
    return


# StoreForces does not give the near/far split, but doesn't block on NearForce,
# which may be useful for debugging
# We only return a single acc slab for plotting; all the accelerations might be too big!
def analyze_storeforces(params, dtype, slabfns=None, silent=False, raw=False):
    if slabfns is None:
        slabfns = sorted(Path(params.OutputDirectory).glob("acc_*"))

    n1d = int(round(params.NP**(1./3)))
    
    # Apply the Zel'dovich approximation scaling of forces to displacements, in units of the interparticle spacing
    rescale = 1/(1.5*params.Omega_M)*n1d
    
    running_sum = 0.
    running_sum_of_square = 0.
    running_max = 0.
    running_min = np.inf
    running_min_nz = np.inf

    io_timer = Tools.ContextTimer('Read acc', cumulative=True, output=False)
    io_bytes = 0
    compute_timer = Tools.ContextTimer('Compute stats', cumulative=True, output=False)
    compute_np = 0

    for slabfn in slabfns:
        with io_timer, open(slabfn, 'rb') as fp:
            acc = np.fromfile(fp, dtype = dtype)
            io_bytes += fp.tell()

        compute_timer.Start()

        if acc.size == 0:
            continue  # empty slab, maybe a (very) sparse test?

        # Warning: extremely sketchy! How do we tell if accstruct was FLOAT3 or FLOAT3p1?
        if acc.size%4 == 0:
            acc = acc.reshape(-1,4)
        else:
            acc = acc.reshape(-1,3)
        acc = acc[:,:3]
        compute_np += len(acc)
        
        assert(np.isfinite(acc).all())

        acc *= rescale
    
        acc_mag = np.sqrt((acc**2).sum(axis=-1))

        min_nz = acc_mag[acc_mag != 0].min()
        
        running_max = max(running_max, acc_mag.max())
        running_min = min(running_min, acc_mag.min())
        running_min_nz = min(running_min_nz, min_nz)
        running_sum += acc_mag.sum(dtype=np.float64)
        running_sum_of_square += (acc_mag**2).sum(dtype=np.float64)

        compute_timer.stop(report=False)

    io_timer.report(end=''); print(' ({:.3g} MB/s)'.format(io_bytes/1e6/io_timer.elapsed))
    compute_timer.report(end=''); print(' ({:.3g} Mpart/s)'.format(compute_np/1e6/compute_timer.elapsed))
    
    results = {'Max':running_max,
               #'99%':np.sort(mag(tabs))[int(.99*len(tabs))],
               'Sum':running_sum,
               'Min':running_min,
               'Min (nonzero)':running_min_nz,
               'Sum of squares':running_sum_of_square,
               'param':dict(params),
               }

    if not raw:
        finalize_results(results, params.NP)
    
    if not silent:
        print('Force statistics (equivalent ZA displacement in units of interparticle spacing):')
        for k in results:
            if k == 'param':
                continue
            print('{k}: {v}'.format(k=k, v=results[k]))
    
    return acc_mag, results
    

def finalize_results(results, NP):
    results['Mean'] = results.pop('Sum')/NP
    results['RMS'] = np.sqrt(results.pop('Sum of squares')/NP)

    
def analyze_forceoutputdebug(params, dtype):
    n1d = int(round(params.NP**(1./3)))
    
    # Need to do the Kick rescaling to put these forces in the same units as StoreForces
    rescale = 3.0*params.Omega_M/(8.0*np.pi*params.NP)

    # Additionally apply the Zel'dovich approximation scaling of forces to displacements, in units of the interparticle spacing
    rescale *= 1/(1.5*params.Omega_M)*n1d
    
    running_sum = 0.
    running_sum_of_square = 0.
    running_max = 0.
    running_min = np.inf
    running_min_nz = np.inf
    for nearfn, farfn in zip(sorted(Path(params.OutputDirectory).glob("nearacc_*")),
                             sorted(Path(params.OutputDirectory).glob("faracc_*"))):
        nacc = np.fromfile(nearfn, dtype=dtype)
        facc = np.fromfile(farfn, dtype=dtype)
        nacc *= rescale
        facc *= rescale

        facc = facc.reshape(-1,3)
        if nacc.size == facc.size:
            nacc = nacc.reshape(-1,3)
        else:
            nacc = nacc.reshape(-1,4)
            nacc = nacc[:,:3]
        
        assert(np.isfinite(nacc).all())
        assert(np.isfinite(facc).all())

        acc_mag = np.sqrt(((nacc+facc)**2).sum(axis=-1))

        min_nz = acc_mag[acc_mag != 0].min()

        diff = np.sqrt(((nacc-facc)**2).sum(axis=-1))/acc_mag
        
        running_max = max(running_max, acc_mag.max())
        running_min = min(running_min, acc_mag.min())
        running_min_nz = min(running_min_nz, min_nz)
        running_sum += acc_mag.sum(dtype=np.float64)
        running_sum_of_square += (acc_mag**2).sum(dtype=np.float64)
    
    results = {'Max':running_max,
               #'99%':np.sort(mag(tabs))[int(.99*len(tabs))],
               'Mean':running_sum/params.NP,
               'Min':running_min,
               'Min (nonzero)':running_min_nz,
               'RMS':np.sqrt(running_sum_of_square/params.NP)
               }
    
    print('Force statistics (equivalent ZA displacement in units of interparticle spacing):')
    for k in results:
        print('{k}: {v}'.format(k=k, v=results[k]))
        
    return diff, acc_mag, results
    

def plot_storeforces(fmag, figfn='checkforces_storeforces_absolute.png'):
    fig, ax = plt.subplots()
    print(f'Plotting to {figfn}')
    
    #p.plot(force, diff, '.', markersize=1)
    plt_log_hist(ax, fmag, label='Force magnitude')
    ax.set_title('Absolute residual forces')

    fig.tight_layout()
    fig.savefig(figfn)
    

def _get_all_acc(param, dtype):
    
    accfns = sorted(Path(param['OutputDirectory']).glob('acc_*'))
    auxfns = sorted(Path(param.get('ReadStateDirectory', Path(param['WorkingDirectory']) / 'read')).glob('aux*_*'))
    pid_bitmask=0x7fff7fff7fff
    
    NP = param['NP']
    cpd = param['CPD']
    ppd = int(round(NP**(1/3)))
    acc = np.empty((NP,3), dtype=np.float64)
    pid = np.empty(NP, dtype=np.uint64)
    
    # Apply the Zel'dovich approximation scaling of forces to displacements, in units of the interparticle spacing
    rescale = 1/(1.5*param.Omega_M)*ppd
    
    i = 0
    for fn in accfns:
        thisacc = np.fromfile(fn, dtype=dtype).reshape(-1,4).astype(np.float64)
        acc[i:i+len(thisacc)] = thisacc[:,:3]*rescale
        i += len(thisacc)
    assert i == NP
        
    i = 0
    for fn in auxfns:
        thisaux = np.fromfile(fn, dtype=np.uint64)
        pid[i:i+len(thisaux)] = thisaux & pid_bitmask
        i += len(thisaux)
    assert i == NP
    
    # x is in the low bits, z in the high bits
    iord = pid.argsort()
    acc = acc[iord]
    pid = pid[iord]
    
    return acc, pid
    

def plot_cell(param, dtype):
    '''Plot a vector field of the force residuals within a cell
    '''

    NP = param['NP']
    cpd = param['CPD']
    ppd = int(round(NP**(1/3)))
    
    acc, pid = _get_all_acc(param, dtype)
    
    plane = 0
    # acc now in [y,x] order
    acc = acc[plane*ppd**2:(plane+1)*ppd**2].reshape(ppd,ppd,3)
    
    # Get the cell-wrapped order
    pos = np.linspace(0,1,ppd,endpoint=False)
    pos %= 1/cpd
    iord = np.atleast_2d(pos.argsort())
    
    # apply cell-wrapped order to each axis, careful not to transpose
    acc = acc[iord.T,iord]
    
    fig, ax = plt.subplots(figsize=(3.5*1.5,2.6*1.5))
    ax.set_aspect('equal')
    
    # quiver uses [j,i] indexing
    # but our acc is ordered [y,x], so that's perfect!
    ax.quiver(acc[...,0], acc[...,1])
    ax.set_xlabel('$x$ cell offset')
    ax.set_ylabel('$y$ cell offset')
    
    max_disp = np.sqrt((acc[...,[0,1]]**2).sum(axis=-1)).max()
    rms_disp = np.sqrt((acc[...,[0,1]]**2).mean())
    
    def latex_float(f):
        float_str = "{0:.2g}".format(f)
        if "e" in float_str:
            base, exponent = float_str.split("e")
            return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
        else:
            return float_str
    
    # Descriptive text box
    ppc = NP/param["CPD"]**3
    text = ['Cell-folded Lattice Residuals',
            f'NP ${ppd:d}^3$, {{}} precision'.format('Single' if dtype == 'f4' else 'Double'),
            f'CPD {param["CPD"]}, Order {param["Order"]}',
            rf'Max err: ${latex_float(max_disp)}$',
            #rf'RMS err: ${latex_float(rms_disp)}$',
           ]
    text = '\n'.join(text)
    ax.text(.03, .6, text, transform=ax.transAxes, ha='left', va='top',
             bbox=dict(fc='w', alpha=0.9, ec='k'))
    
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(ppd//8))
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(ppd//8))
    
    ax.set_xticks(list(range(0,ppd,8)) + [ppd-1,])
    ax.set_yticks(list(range(0,ppd,8)) + [ppd-1,])
    
    fig.tight_layout()
    
    fig.savefig(f'cell_residuals_{ppd}.pdf', bbox_inches='tight')
    
    
def plot_force_error_pk(param, dtype):
    '''Measure the Lagrangian power spectrum of the error components
    '''
    
    NP = param['NP']
    cpd = param['CPD']
    ppd = int(round(NP**(1/3)))
    
    # Get the lattice acc
    acc, pid = _get_all_acc(param, dtype)
    acc = acc.reshape(ppd,ppd,ppd,3)
    
    # Now generate the ZA acc
    #from Abacus import zeldovich
    #from Abacus import ReadAbacus
    from Abacus import abacus
    #zeldovich.run_override_dirs('abacus.par', Path(param['InitialConditionsDirectory']).parents[1], ICFormat='RVZel')
    param_kwargs = dict(
        StoreForces=1,
        Order=param['Order'],
        BoxSize=param['BoxSize'],
        ICFormat='RVZel',
        NP=NP,
        CPD=cpd,
    )

    abacus.run('checkforces.par2', param_kwargs=param_kwargs, maxsteps=2, clean=True)
    #za = ReadAbacus.from_dir(param['InitialConditionsDirectory'], format='rvzel',
    #                         return_vel=False)
    
    # The zeldovich code outputs in same units as BoxSize, so divide by the particle spacing
    #za = za['pos'] / (param['BoxSize']/ppd)
    #za = za.reshape(ppd,ppd,ppd,3)
    #print('RMS ZA: {}'.format((za[...,0]**2).mean()**.5))
    
    # Apply the Zel'dovich approximation scaling of forces to displacements, in units of the interparticle spacing
    #rescale = 1/(1.5*param.Omega_M)*ppd
    
    za, pid = _get_all_acc(param, dtype)
    za = za.reshape(ppd,ppd,ppd,3)
    
    from Abacus.Analysis.PowerSpectrum import PowerSpectrum
    k, Px, nmodes = PowerSpectrum.FFTAndBin(acc[...,0], 1, normalize_dens=False, window=None, bins=ppd, multipoles=0)
    k, za_Px, nmodes = PowerSpectrum.FFTAndBin(za[...,0], 1, normalize_dens=False, window=None, bins=ppd, multipoles=0)
    #k, Py, nmodes = PowerSpectrum.FFTAndBin(acc[...,1], 1., normalize_dens=False, bins=ppd)
    #k, Pz, nmodes = PowerSpectrum.FFTAndBin(acc[...,2], 1., normalize_dens=False, bins=ppd)
    
    kny = np.pi*ppd
    k /= kny
    
    fig, ax = plt.subplots(figsize=(3.5,2.6))
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.plot(k,za_Px, c='k', ls='--', label='$z=99$ $\Lambda$CDM')
    ax.plot(k,Px, c='k', ls='-', label='Force error')
    #ax.plot(k,Py, c='k', ls='--', label='$y$')
    #ax.plot(k,Pz, c='k', ls=':', label='$z$')
    
    #ax.axvline(np.pi*cpd/kny, label='Cell', ls=':', c='k')
    #ax.axvline(np.pi*cpd/param['NearFieldRadius']/kny, label='Near Field Radius', ls='--', c='k')
    ax.legend()
    
    ax.set_xlabel('$k/k_\mathrm{Ny}$')
    ax.set_ylabel('$P(k)$')
    
    fig.tight_layout()
    fig.savefig('force_err_pk.pdf', bbox_inches='tight')
    

    
def plot_forceoutputdebug(diff, fmag):
    figfn = 'checkforces_forceoutputdebug_near-far_cancellation.png'
    print('Plotting to {}'.format(figfn))
    
    fig, ax = plt.subplots()
    plt_log_hist(ax, diff)
    ax.set_title('Cancellation of near and far field forces')

    fig.tight_layout()
    fig.savefig(figfn)
    
    figfn = 'checkforces_forceoutputdebug_absolute.png'
    print('Plotting to {}'.format(figfn))
    
    fig, ax = plt.subplots()
    plt_log_hist(ax, fmag, label='Force magnitude')
    ax.set_title('Absolute residual forces')

    fig.tight_layout()
    fig.savefig(figfn)
    
    
def load_results(load, order=None, dtype=None):
    with asdf.open(load, lazy_load=False, copy_arrays=True) as af:
        results = af.tree['results']
        
    if order != None:
        for res in results:
            if res['param']['Order'] == order:
                break
        else:
            raise ValueError(order)
        raise NotImplementedError("cell results aren't stored in asdf, can't load")
        if not res['param']['Parallel']:
            plot_cell(res, dtype)
    else:
        plot_orders(results)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('ppd', help='Particles per dimension', nargs='?', default=DEFAULT_PPD)
    parser.add_argument('cpd', help='Cells per dimension', nargs='?', default=DEFAULT_CPD)
    parser.add_argument('--force-output-debug', help='Use ForceOutputDebug. Slow, but you get the near and far field separately',
                        action='store_true', default=False)
    parser.add_argument('--no-analyze', help='Run singlestep so the forces are stored on disk but do not analyze them',
                        action='store_true', default=False)
    parser.add_argument('--dtype', help='The precision Abacus was compiled with (i.e. the precision of the forces)',
                        default='f4', choices=['f4','f8'])
    parser.add_argument('--order', help='Multipole order', default=DEFAULT_ORDER, type=int)
    parser.add_argument('--sweep', help='Run orders 2 through ORDER', action='store_true')
    parser.add_argument('--save', help='Save results to this filename', type=str)
    parser.add_argument('--load', help="Don't run anything and instead load these results", type=str)

    args = parser.parse_args()
    args = vars(args)

    with Tools.chdir(Path(__file__).parent):
        load = args.pop('load')
        order = args.pop('order')
        dtype = args.get('dtype')
        if load:
            load_results(load, order, dtype=dtype)
        else:
            if args.pop('sweep'):
                save = args.pop('save')
                run_orders(orders=range(2,order+1), run_kwargs=args, save=save)
            else:
                run(**args)
