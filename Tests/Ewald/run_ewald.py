#!/usr/bin/env python3
'''
The "New Ewald" test. Compares Abacus forces to the result computed
in quad-double precision with Ewald summation on 65536 particles.

On ForceOutputDebug vs StoreForces: StoreForces should be considered
the more accurate representation of the forces particles feel, since
the near and far field are combined within Abacus then output.  But
if one is changing a parameter and wants to see which component
changes, then one can use -fod to get the components separately.
One might also want this to see how much error arises from catastrophic
cancellation of large near and far forces.
'''

import argparse
from pathlib import Path
import shutil

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import ticker
import matplotlib.text
import matplotlib.pyplot as plt
import asdf

from Abacus import abacus
from Abacus.Tools import scatter_density, add_tick

DEFAULT_CPD = 11
DEFAULT_ORDER = 8
DEFAULT_DTYPE = 'f4'

POS_FN = Path(__file__).parent / 'reference_ewald_pos.double3'

def run(cpd=DEFAULT_CPD, order=DEFAULT_ORDER, dtype=DEFAULT_DTYPE, fod=False, save=None):
    '''Run the simulation and check the results
    '''
        
    dtype = np.dtype(dtype)
    #if dtype != np.float32:
    #    raise NotImplementedError(f'Ewald only supports np.float32, not {dtype}')
        
    sim_params = dict(CPD=cpd, Order=order, ForceOutputDebug=int(fod), StoreForces=int(not fod))
        
    # TODO: no good way to process params before calling abacus.run()
    par2fn = Path(__file__).parent / 'abacus.par2'
    parfn = Path(__file__).parent / 'abacus.par'
    param = abacus.preprocess_params(parfn, par2fn, param_kwargs=sim_params)
    
    make_ic_files(param, erase_ic=True)
    
    abacus.run(par2fn, param_kwargs=sim_params, clean=True, erase_ic=False, maxsteps=2)
    
    if fod:
        results = check_fod(param, dtype=dtype)
    else:
        results = check_storeforces(param, dtype=dtype)
        
    if save:
        print(f'Saving results to {save}')
        af = asdf.AsdfFile(dict(results=[results]))
        af.write_to(save)
        
    return results

def run_orders(orders=[1,2,3,4,5,6,7,8], run_kwargs=None, save=None):
    if run_kwargs is None:
        run_kwargs = {}
        
    all_results = []
    for order in orders:
        res = run(order=order, **run_kwargs)
        all_results += [res]
        
    if save:
        print(f'Saving results to {save}')
        af = asdf.AsdfFile(dict(results=all_results))
        af.write_to(save)

    plot_orders(all_results)
    
    
def plot_orders(all_results, pltfn='ewald_allorders.pdf', style='line'):
    orders = [r['param']['Order'] for r in all_results]
    
    # Descriptive text in corner
    param = all_results[0]['param']  # use params from one result
    NP = param['NP']
    ppc = NP/param["CPD"]**3
    text = ['Abacus Ewald Test',
            f'NP {NP:d}',
            f'CPD {param["CPD"]}, ppc {ppc:.1f}',
            'Single precision'
           ]
    text = '\n'.join(text)
    
    if style == 'viola':
        fig, ax = plt.subplots(figsize=(0.5*len(orders),4))
        #ax.set_yscale('log')
        ax.set_xlabel('Multipole order $p$')
        ax.set_ylabel('Relative force error')

        violaplot = ax.violinplot  # viola life
        for res in all_results:
            o = res['param']['Order']
            #ax.hist(res['frac_diff'], bins=np.geomspace(res['frac_diff'].min(), res['frac_diff'].max(),200), alpha=0.3, label=o)
            parts = violaplot(np.log10(res['frac_diff']), positions=[o], vert=True, widths=[0.75],
                     showextrema=True, showmeans=False, showmedians=True, bw_method=1e-2, points=100)

            for pc in parts['bodies']:
                pc.set_facecolor('black')
                pc.set_edgecolor('black')
                #pc.set_alpha(1)
            for k in ('cmins','cmaxes','cmedians','cbars'):
                parts[k].set_color('k')

        ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        #ax.yaxis.set_ticks([np.log10(x) for p in range(-6,1) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)

        fig.text(.98, .97, text, transform=ax.transAxes, ha='right', va='top',
                 bbox=dict(fc='w', alpha=0.5, ec='k'))

        fig.tight_layout()
        fig.savefig(pltfn.replace('.pdf','.viola.pdf'), bbox_inches='tight')
        
    elif style == 'line':
        fig, ax = plt.subplots(figsize=(3.5,2.6))
        ax.set_xlabel('Multipole order $p$')
        ax.set_ylabel('Relative force error')
        
        all_orders = [res['param']['Order'] for res in all_results]
        all_medians = [res['median'] for res in all_results]
        all_medians = [res['median'] for res in all_results]
        all_1p, all_99p = [], []
        for res in all_results:
            f = res['frac_diff']
            f.sort()
            all_1p += [f[len(f)//100]]
            all_99p += [f[-len(f)//100-1]]
        
        ax.set_yscale('log')
        ax.plot(all_orders, all_medians, marker='o', c='k')
        ax.fill_between(all_orders, all_1p, all_99p, alpha=0.2, color='k')
        
        fig.text(.95, .95, text, transform=ax.transAxes, ha='right', va='top',
                 bbox=dict(fc='w', alpha=0.5, ec='k'))
        #fig.text(0.05, .05, text, transform=ax.transAxes, ha='left', va='bottom',
        #         bbox=dict(fc='w', alpha=0.5, ec='k'))
        
        ax.xaxis.set_ticks(all_orders[::2])
        yticks = np.logspace(-1,-6,6)
        ax.yaxis.set_ticks(yticks)
        ax.yaxis.set_ticklabels([f'$10^{{{np.log10(y):.0f}}}$' for y in yticks])
        #plt.setp(ax.get_yticklabels()[1::2], visible=True)
        ax.tick_params(axis='y', right=True, direction='in')
        #ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))

        fig.tight_layout()
        fig.savefig(pltfn, bbox_inches='tight')
    plt.close(fig)


def load_results(load):
    with asdf.open(load, lazy_load=False, copy_arrays=True) as af:
        results = af.tree['results']
    
    plot_orders(results)
    
    
def check_storeforces(param, dtype=DEFAULT_DTYPE):
    '''Check the results, using StoreForces
    '''
    
    NP = param['NP']
    rescale = 3.0*(param['Omega_M']-param.get('Omega_Smooth',0.))/(8.0*np.pi)
    
    acc = np.empty((NP,3), dtype=dtype)
    i = 0
    for accfn in sorted(Path(param.get('OutputDirectory')).glob('acc_*')):
        # TODO: better detection of acc3 or acc4
        thisacc = np.fromfile(accfn, dtype=dtype).reshape(-1,4)
        acc[i:i+len(thisacc)] = thisacc[:,:3]
        i += len(thisacc)
    assert i == NP

    read = Path(param.get('ReadStateDirectory', Path(param['WorkingDirectory']) / 'read' ))
    auxfns = sorted(read.glob('aux*_*'))
    if not auxfns:
        # aux might be on another node and thus inaccessible, but maybe we did it all on this node
        auxfns = Path(param.get('LocalWorkingDirectory')).parent.glob(param['SimName'] + '*/read/aux*')
        auxfns = sorted(auxfns, key=lambda p:p.name)
    
    pid = np.empty(NP, dtype=np.int64)
    i = 0
    for auxfn in auxfns:
        thisaux = np.fromfile(auxfn, dtype=np.uint64)
        # we spread out the PID bits into three chunks in the aux field
        # but since we only care about the order, we can ignore that and just sort
        pid[i:i+len(thisaux)] = thisaux & np.uint64(0x7fff7fff7fff)
        i += len(thisaux)
    assert i == NP, (i,NP)
        
    # put the acc back in the original order
    acc = acc[pid.argsort()]
    
    acc = acc.astype(np.float64, copy=False)  # let's do any comparison math in float64
    acc /= rescale
    
    ref_acc = np.fromfile(Path(__file__).parent / 'reference_ewald_acc.double3', dtype=np.float64).reshape(-1,3)
    
    diff = acc - ref_acc
    diff_mag = np.sqrt((diff**2).sum(axis=-1))
    ref_mag = np.sqrt((ref_acc**2).sum(axis=-1))
    frac_diff = diff_mag/ref_mag
    
    frac_diff.sort()
    frac_diff = frac_diff[::-1]
    
    res = {
            'max':frac_diff[0],
            '99%':frac_diff[NP//100],
         'median':frac_diff[NP//2],
            'min':frac_diff[-1],
    'min_nonzero':frac_diff[frac_diff>0][-1],
      'frac_diff':frac_diff,
          'param':param,
        }
    
    restxt = f'''
    (New) Ewald Test
    ================
    
    Fractional error
    ----------------
              Max: {res['max']:.4e}
              99%: {res['99%']:.4e}
           Median: {res['median']:.4e}
              Min: {res['min']:.4e}
    Min (nonzero): {res['min_nonzero']:.4e}
    '''
    
    settings = f'''
    Configuration
    -------------
      CPD: {param['CPD']:d}
    Order: {param['Order']:d}
      ppc: {NP/param['CPD']**3:.1f}
    dtype: {dtype}
    '''
    
    print(restxt)
    print(settings)
    
    plot_storeforces(param, acc, ref_acc)
    
    return res


def plot_storeforces(param, acc, ref_acc, pltfn='ewald_storeforces.png', show_perc=True):
    '''Make some plots, like relative err vs absolute force
    '''
    
    def _plot(figsize, pltfn):
        figsize = np.array(figsize)
        print(f'Plotting to {pltfn} ...', flush=True)  # might take a few seconds to do the KDE

        fig, ax = plt.subplots(figsize=figsize)

        ax.set_xscale('log')
        ax.set_yscale('log')

        # We'll just repeat these calculations in case we want to plot something different
        diff = acc - ref_acc
        diff_mag = np.sqrt((diff**2).sum(axis=-1))
        ref_mag = np.sqrt((ref_acc**2).sum(axis=-1))
        frac_diff = diff_mag/ref_mag

        #ax.scatter(ref_mag, frac_diff, c='k', s=2, alpha=0.5)
        scatter_density(ref_mag, frac_diff, ax, size=2., log=True, bw=0.01)
        ax.set_xlabel('Force magnitude')
        ax.set_ylabel('Relative error')

        NP = param['NP']
        ppc = NP/param["CPD"]**3
        text = ['Abacus Ewald Test',
                f'NP {NP:d}, ppc {ppc:.1f}',
                f'CPD {param["CPD"]}, Order {param["Order"]}'
               ]
        text = '\n'.join(text)
        
        ax.text(.965, .95, text, transform=ax.transAxes, ha='right', va='top',
                 bbox=dict(fc='w', alpha=0.9, ec='k'))
        
        if show_perc:
            rms = (frac_diff**2).mean()**.5
            frac_diff.sort()
            frac_diff = frac_diff[::-1]
            median = frac_diff[len(frac_diff)//2]
            #add_tick(ax, rms, 'RMS', which='y')
            ax.axhline(median, ls=':', c='k')
            ax.annotate('Median', (0.015, median*.8), xycoords=('axes fraction','data'), ha='left', va='top')
            
            p99 = frac_diff[len(frac_diff)//100]
            #add_tick(ax, rms, 'RMS', which='y')
            ax.axhline(p99, ls=':', c='k')
            ax.annotate('99%', (0.015, p99*.8), xycoords=('axes fraction','data'), ha='left', va='top')

        fig.tight_layout()
        fig.savefig(pltfn, bbox_inches='tight')
        
    _plot((3.5,2.6), pltfn)
    _plot((7.0,5.2), pltfn.replace('.png','.large.png'))
    

def make_ic_files(param, erase_ic=True):
    '''Divide the Ewald positions file into CPD slab files
    
    Of course, this is not strictly required, but for
    parallel Ewald this is highly convenient.
    
    Not super optimized, will not scale to huge N
    '''
    
    p = np.fromfile(POS_FN, dtype=np.float64).reshape(-1,3)
    assert np.all(np.abs(p) < 0.5)  # sanity check
    
    iord = p[:,0].argsort()
    p = p[iord]  # sort on x
    
    cpd = param['CPD']
    splits = np.linspace(-.5,.5, cpd, endpoint=False)[1:]
    isplit = np.searchsorted(p[:,0], splits)
    chunks = np.array_split(p, isplit)
    ichunks = np.array_split(iord, isplit)
    assert len(chunks) == cpd
    
    ic_dir = Path(param['InitialConditionsDirectory'])
    
    if erase_ic and ic_dir.exists():
        shutil.rmtree(ic_dir)
    ic_dir.mkdir(parents=True, exist_ok=True)
    
    rvpid_dtype = np.dtype([('pos',p.dtype,3),('vel',p.dtype,3),('pid',np.int64)], align=True)
    
    for i in range(cpd):
        slabfn = ic_dir / f'ic_{i:d}'
        p = chunks[i]
        rvpid = np.empty(len(p), dtype=rvpid_dtype)
        rvpid['pos'][:] = p
        rvpid['vel'][:] = 0.
        rvpid['pid'][:] = ichunks[i]
        
        rvpid.tofile(slabfn)
    
    
class ArgParseFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=ArgParseFormatter)
    parser.add_argument('-cpd', help='Cells-per dimension', type=int, default=DEFAULT_CPD)
    parser.add_argument('-order', help='Multipole order', type=int, default=DEFAULT_ORDER)
    parser.add_argument('-sweep', help='Run orders 2 through ORDER', action='store_true')
    parser.add_argument('-dtype', help='Float type', type=str, default='f4', choices=['f4','f8'])
    parser.add_argument('-fod', help='Use ForceOutputDebug to get separate near and far forces', action='store_true')
    parser.add_argument('-save', help='Save results to this filename', type=str)
    parser.add_argument('-load', help="Don't run anything and instead load these results", type=str)
 
    args = parser.parse_args()
    args = vars(args)
    
    if args['fod']:
        raise NotImplementedError('-fod not yet implemented')
    
    load = args.pop('load')
    if load:
        load_results(load)
    else:
        if args.pop('sweep'):
            order = args.pop('order')
            save = args.pop('save')
            run_orders(orders=range(2,order+1), run_kwargs=args, save=save)
        else:
            run(**args)
