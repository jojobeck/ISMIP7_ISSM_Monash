"""Assemble the cold_ensemble.nc and warm_ensemble.nc for the meltMIP local
parameterisation tuning (ISMIP7), from the ISSM-computed, gridded ocean-modelling
melt rates produced by the create_BMB_gD_OceanModelling MATLAB step.

Structure expected by parameter_selection_quadratic_example.ipynb:
  dims  : (model, p2, p1, y, x)
  coords: model = ["mathiot", "naughten_ais_1"], p1 = K values, p2 = [1.0]
  var   : melt_rate  (kg/m2/a)

Mirrors the unit handling of creating_bmbMIP_ensmble.py (pd_ensemble): the gridded
netcdf already holds melt_rate = BMB[m/a] * rho_ice, i.e. kg/m2/a, so it is used
as-is. Only the recommended circum-Antarctic subset is used: Mathiot ("mathiot") and
Naughten FESOM-ACCESS ("naughten_ais_1"). No regional masking is needed for these.
"""
import xarray as xr
import numpy as np

# K parameter values (p1), same range as the present-day ensemble
K = np.arange(0.25e-5, 3.025e-4, 0.25e-5)

ncdir = 'Models/ModelNC/BMBOceanModelling'

# model coordinate name -> MATLAB output tag, per ocean state
cold_tags = {'mathiot': 'Mc', 'naughten_ais_1': 'Nc'}
warm_tags = {'mathiot': 'Mw', 'naughten_ais_1': 'Nw'}

models = ['mathiot', 'naughten_ais_1']


def build_ensemble(tag_map):
    model_ensembles = []
    for model in models:
        tag = tag_map[model]
        runs = []
        for i, k in enumerate(K):
            j = i + 1  # MATLAB index starts at 1
            ncfile = f'{ncdir}/bmelt_OceanModelling_{tag}_K_{j}_gridData_8km.nc'
            ds = xr.open_dataset(ncfile)['melt_rate']  # already kg/m2/a
            ds = ds.expand_dims({'p1': 1}).assign_coords({'p1': [k]})
            runs.append(ds)

        ens = xr.concat(runs, dim='p1')
        ens = ens.expand_dims({'p2': np.ones(1)})
        ens = ens.expand_dims({'model': [model]})
        model_ensembles.append(ens)

    result = xr.concat(model_ensembles, dim='model')
    return result.to_dataset(name='melt_rate')


print('Building cold ensemble')
cold_ensemble = build_ensemble(cold_tags)
print(cold_ensemble)

print('Building warm ensemble')
warm_ensemble = build_ensemble(warm_tags)
print(warm_ensemble)

cold_path = f'{ncdir}/cold_ensemble.nc'
warm_path = f'{ncdir}/warm_ensemble.nc'
cold_ensemble.to_netcdf(cold_path)
warm_ensemble.to_netcdf(warm_path)
print(f'Saved: {cold_path}')
print(f'Saved: {warm_path}')
