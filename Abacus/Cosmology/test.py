#!/usr/bin/env python3
# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

'''
Compare our cosmology module with Colossus.

TODO: currently Colossus gives the wrong answer for wCDM
cosmologies.  So these tests will not work until that is fixed!
See: https://bitbucket.org/bdiemer/colossus/issues/8/

Run with:
$ python -m Abacus.Cosmology.test

'''

import pytest
from . import AbacusCosmo
from colossus.cosmology import cosmology as coco
import numpy as np

isclose = lambda *args, **kwargs: np.isclose(*args, rtol=1e-3, **kwargs)

@pytest.fixture(params=[(0.3, 0.7, 70., -1.0, 0. ),
                        (0.3, 0.7, 70., -1.3, 0. ),
                        (0.3, 0.7, 70., -0.7, 0.2),
                        #(0.6, 0.2, 50., -0.7, -0.2),  # evolving DE with curvature doesn't seem to agree
                        (1.0, 0.0, 50., -1., 0. ),
                        (1.5, 0.1, 50., -1., 0. ),
                        (0.5, 0.1, 50., -1., 0. ),
                       ])
def setup_cosmology(request):
    Omega_M, Omega_DE, H0, w0, wa = request.param

    params = {'H0':H0, 'w0':w0, 'wa':wa, 'Om0':Omega_M, 'Ob0':0., 'sigma8':1., 'ns':.96, 'relspecies': False}
    if np.isclose(Omega_M + Omega_DE, 1.):
        params['flat'] = True
    else:
        params['flat'] = False
        params['Ode0'] = Omega_DE

    if wa == 0.:
        if w0 == -1:
            params['de_model'] = 'lambda'
        else:
            params['de_model'] = 'w0'
    else:
        params['de_model'] = 'w0wa'

    #params['de_model'] = 'user'
    #params['wz_function'] = lambda z: -1.3
    #print(params)
    cc_cosm = coco.setCosmology('colossus_cosm', params)
    
    my_cosm = AbacusCosmo.MyCosmology()
    my_cosm.Omega_m = Omega_M
    my_cosm.Omega_DE = Omega_DE
    my_cosm.Omega_K = 1 - Omega_M - Omega_DE
    my_cosm.H0 = H0
    my_cosm.w0 = w0
    my_cosm.wa = wa
                
    return my_cosm, cc_cosm

@pytest.mark.parametrize('z',
                        [0., 1., 9., 99.],
                        )
def test_cosmology(setup_cosmology, z):
    my_cosm, cc_cosm = setup_cosmology
    
    a = 1./(1 + z)
    abacus_cosm = AbacusCosmo.Cosmology(a, my_cosm)
    
    # check linear growth factor D(z)
    assert isclose(cc_cosm.growthFactorUnnormalized(z)/cc_cosm.growthFactorUnnormalized(0.),
                    abacus_cosm.current.growth/abacus_cosm.today.growth)

    # check H(z)
    assert isclose(cc_cosm.Hz(z), abacus_cosm.current.H)

    # check evolution integrand
    total = abacus_cosm.current.OmegaHat_X + abacus_cosm.current.OmegaHat_m + abacus_cosm.current.OmegaHat_K
    assert isclose(cc_cosm.Ez(z)**2, total)

    # check omega_X(z)
    assert isclose(cc_cosm.Ode(z), abacus_cosm.current.OmegaHat_X/total)
    assert isclose(cc_cosm.Om(z), abacus_cosm.current.OmegaHat_m/total)

    # check w(z)
    assert isclose(cc_cosm.wz(z), abacus_cosm.current.w)
    
if __name__ == '__main__':
    pytest.main()
