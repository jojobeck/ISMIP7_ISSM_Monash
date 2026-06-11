"""Assemble obs_ensemble.nc for the meltMIP local parameterisation tuning
(ISMIP7), from the ISSM-computed, gridded observation-forced melt rates produced
by the create_BMB_gD_ObsData MATLAB step.

Structure expected by parameter_selection_quadratic_example.ipynb
(apply_ocean_observations_melt):
  dims  : (year, p2, p1, y, x)
  coords: year = observation years, p1 = K values, p2 = [1.0]
  var   : melt_rate  (kg/m2/a)

The obs_ensemble is indexed only by year (+ p1, p2). The pig/dotson split is a
separate spatial mask (region_label) applied inside calculate_term4, so it is
NOT a dimension here. ISMIP7 README recommends years 2009 and 2012 for Pine Island.

Mirrors the unit handling of creating_bmbMIP_ensmble.py: the gridded netcdf holds
rho_ice * BMB[m/s]; multiplying by yearsinsec gives kg/m2/a.
"""
import xarray as xr
import numpy as np

yearsinsec = 31556926.080000002

# K parameter values (p1), same range as the present-day ensemble
K = np.arange(0.25e-5, 3.025e-4, 0.25e-5)

# Observation years (must match obs_years in the MATLAB steps)
years = [2009, 2012]

ncdir = 'Models/ModelNC/BMBObsData'

year_ensembles = []
for year in years:
    runs = []
    for i, k in enumerate(K):
        j = i + 1  # MATLAB index starts at 1
        ncfile = f'{ncdir}/bmelt_ObsData_{year}_K_{j}_gridData_8km.nc'
        ds = xr.open_dataset(ncfile)['melt_rate'] * yearsinsec
        ds = ds.expand_dims({'p1': 1}).assign_coords({'p1': [k]})
        runs.append(ds)

    ens = xr.concat(runs, dim='p1')
    ens = ens.expand_dims({'p2': np.ones(1)})
    ens = ens.expand_dims({'year': [year]})
    year_ensembles.append(ens)

print('Combining datasets')
obs_ensemble = xr.concat(year_ensembles, dim='year')
obs_ensemble = obs_ensemble.to_dataset(name='melt_rate')
print(obs_ensemble)

out_path = f'{ncdir}/obs_ensemble.nc'
obs_ensemble.to_netcdf(out_path)
print(f'Saved: {out_path}')
