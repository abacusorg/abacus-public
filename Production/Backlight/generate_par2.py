#!/usr/bin/env python

"""
Generate AbacusBacklight parameter files for a cosmology grid.
"""

from pathlib import Path

import click
import numpy as np
from astropy.table import Table

TEMPLATE = """SimName = "{SimName}"
#include "AbacusBacklight_base.par2"

SimComment = "Cosmology grid"

Omega_M = {Omega_M}
Omega_DE = 1 - @Omega_M@
H0 = {H0}
ZD_Pk_filename = "{ZD_Pk_filename}"
ZD_Seed = {ZD_Seed}

CAMB_sigma8 = {sigma8}
CAMB_ns = {ns}
CAMB_Omega_b = {Omega_b}
"""

SEED_PH000 = 0x12B1A761
CAMB_DIR = Path('/mnt/ceph/users/wcoulton/nbodysims/abacusBacklight/camb_files')


def seed(cosm):
    return SEED_PH000 + cosm


def write_one(outdir, params):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    cosm = params['cosm']
    fn = f'AbacusBacklight_base_c{cosm:04d}_ph{cosm:04d}.par2'

    template_params = {
        'SimName': f'AbacusBacklight_base_c{cosm:04d}_ph{cosm:04d}',
        'Omega_M': params['Om'],
        'Omega_b': params['Ob'],
        'sigma8': params['s8'],
        'H0': 100 * params['h'],
        'ZD_Pk_filename': CAMB_DIR / f'Pk_m_z=0.000_simID_{cosm}.txt',
        'ZD_Seed': seed(cosm),
        'ns': params['ns'],
    }

    txt = TEMPLATE.format(**template_params)
    # print(txt)
    with open(outdir / fn, 'w') as fp:
        fp.write(txt)


@click.command()
@click.argument(
    'input_fn', type=click.Path(exists=True), default='cosmo_params_final.txt'
)
@click.argument('outdir', type=click.Path(), default='param')
@click.option('-f', '--first', type=int, default=0)
@click.option('-l', '--last', type=int, default=None)
def write_all(input_fn='cosmo_params_final.txt', outdir='param', first=0, last=None):
    params = Table.read(input_fn, format='ascii', names=['Om', 'Ob', 'h', 'ns', 's8'])
    params['cosm'] = np.arange(len(params))

    if last is None:
        last = len(params)
    params = params[first : last + 1]

    for row in params:
        write_one(outdir, row)


if __name__ == '__main__':
    write_all()
