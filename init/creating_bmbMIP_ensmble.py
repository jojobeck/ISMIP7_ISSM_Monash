import xarray as xr
import numpy as np

K = xr.DataArray(
    data=np.arange(0.25e-5, 3.025e-4, 0.25e-5),
    dims=['p1']
).assign_coords({'p1': np.arange(0.25e-5, 3.025e-4, 0.25e-5)})

model_runs = []

for i, k in enumerate(K):
    j = i + 1  # MATLAB index starts at 1

    ncfile = f'Models/ModelNC/BMBPresentDay/bmelt_test_K_{j}_gridData_8km.nc'
    # gridData wrote melt_rate = BMB[m/a] * rho_ice, i.e. already kg/m2/a
    ds = xr.open_dataset(ncfile)['melt_rate']
    ds = ds.expand_dims({'p1': 1}).assign_coords({'p1': [k.values]})
    model_runs.append(ds)

print('Combining datasets')
pd_ensemble = xr.concat(model_runs, dim='p1')
pd_ensemble = pd_ensemble.expand_dims({"p2": np.ones(1)})
pd_ensemble = pd_ensemble.to_dataset(name='melt_rate')

print(pd_ensemble)

# Save
out_path = 'Models/ModelNC/BMBPresentDay/MeltObs_ensemble.nc'
pd_ensemble.to_netcdf(out_path)
print(f'Saved: {out_path}')
