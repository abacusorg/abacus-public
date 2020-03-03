from . import AbacusCosmo

def from_params(params, z='init'):
    if z == 'init':
        z = params['InitialRedshift']
    my_cosmology = AbacusCosmo.MyCosmology()
    my_cosmology.Omega_m = params['Omega_M']
    my_cosmology.Omega_smooth = params.get('Omega_Smooth',0.)
    my_cosmology.Omega_K = params['Omega_K']
    my_cosmology.Omega_DE = params['Omega_DE']
    my_cosmology.H0 = params['H0']
    my_cosmology.w0 = params['w0']
    my_cosmology.wa = params['wa']

    cosmology = AbacusCosmo.Cosmology(1./(1+z), my_cosmology)
    return cosmology
