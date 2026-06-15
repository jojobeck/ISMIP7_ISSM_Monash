"""Plot example basal-melt maps for the meltMIP ensembles against their
observed / target melt rates.

Reads the ensemble netcdfs assembled by the creating_bmbMIP_*.py scripts
(melt_rate already in kg/m2/a):

  present day : ../Models/ModelNC/BMBPresentDay/MeltObs_ensemble.nc   (p2, p1, y, x)
  ocean model : ../Models/ModelNC/BMBOceanModelling/warm_ensemble.nc  (model, p2, p1, y, x)
                ../Models/ModelNC/BMBOceanModelling/cold_ensemble.nc

Two figures are produced:

  1. ../figures/figure_presentday_meltmaps.png
     4 rows x 1 col: modelled melt at Kmin, Kmean, Kmax, and the observed
     present-day target (Paolo/Adusumilli, melt_paolo_err_adusumilli_ismip8km.nc).

  2. ../figures/figure_oceanmodelling_meltmaps.png
     2 rows x 4 cols at Kmean. Top row = modelled melt for Nw/Nc/Mw/Mc
     (Naughten/Mathiot warm/cold); bottom row = the matching ocean-model target
     melt (*_m.nc). Each column (scenario) shares a colour scale so model and
     target are directly comparable.

Everything is plotted in kg/m2/a.

Run from BMB_tuning_python/:  python plot_ensemble_meltmaps.py
"""
import os

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

try:
    import cmocean
    CMAP = cmocean.cm.matter
except ImportError:
    CMAP = plt.get_cmap('magma')
CMAP = CMAP.copy()
CMAP.set_bad('white')

# ------------------------------------------------------------------ constants

ens_pd = '../Models/ModelNC/BMBPresentDay/MeltObs_ensemble.nc'
ens_warm = '../Models/ModelNC/BMBOceanModelling/warm_ensemble.nc'
ens_cold = '../Models/ModelNC/BMBOceanModelling/cold_ensemble.nc'

target_root = './../../raw_data/ISMIP7/AIS/parameterisations/ocean'
paolo_file = os.path.join(target_root, 'meltobs',
                          'melt_paolo_err_adusumilli_ismip8km.nc')
om_target_dir = os.path.join(target_root, 'ocean_modelling_data')
floatmask_file = os.path.join(target_root, 'floatingmasks',
                              'floatingmask_ismip8km.nc')

# tag -> (title, ensemble file, ensemble model coord, ocean-model target file)
OM_SCENARIOS = [
    ('Nw', 'Naughten warm', ens_warm, 'naughten_ais_1',
     'Naughten_FESOM_ACCESS_warm_m.nc'),
    ('Nc', 'Naughten cold', ens_cold, 'naughten_ais_1',
     'Naughten_FESOM_ACCESS_cold_m.nc'),
    ('Mw', 'Mathiot warm',  ens_warm, 'mathiot',
     'Mathiot_NEMO_warm_m.nc'),
    ('Mc', 'Mathiot cold',  ens_cold, 'mathiot',
     'Mathiot_NEMO_cold_m.nc'),
]

figure_dir = '../figures'
os.makedirs(figure_dir, exist_ok=True)


# ------------------------------------------------------------------- helpers
# floating-ice (ice-shelf) mask, applied to the model fields by position
_fmask = xr.open_dataset(floatmask_file)['mask'].values > 0.5


def _finite(field):
    return field.where(np.isfinite(field) & (field != 0))


def robust_vmax(*fields, pct=99.0):
    """Shared upper colour limit from the positive melt of the given fields."""
    vals = np.concatenate([f.values[np.isfinite(f.values)] for f in fields])
    vals = vals[vals > 0]
    return float(np.nanpercentile(vals, pct)) if vals.size else 1.0


def show(ax, field, norm):
    extent = [float(field.x.min()), float(field.x.max()),
              float(field.y.min()), float(field.y.max())]
    im = ax.imshow(field.values, origin='lower', extent=extent,
                   cmap=CMAP, norm=norm, interpolation='nearest')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    return im


# --------------------------------------------------------- figure 1: present day
def plot_presentday():
    ds = xr.open_dataset(ens_pd)['melt_rate'].isel(p2=0)
    p1 = ds.p1.values
    i_min, i_mean, i_max = 0, int(np.argmin(np.abs(p1 - p1.mean()))), p1.size - 1
    print('Present-day: Kmin=%.2e  Kmean=%.2e (idx %d)  Kmax=%.2e'
          % (p1[i_min], p1[i_mean], i_mean, p1[i_max]))

    target = xr.open_dataset(paolo_file)['melt_mean']
    rows = [
        ('Kmin (K=%.2e)' % p1[i_min],   _finite(ds.isel(p1=i_min))),
        ('Kmean (K=%.2e)' % p1[i_mean], _finite(ds.isel(p1=i_mean))),
        ('Kmax (K=%.2e)' % p1[i_max],   _finite(ds.isel(p1=i_max))),
        ('Observed (Paolo/Adusumilli)', _finite(target)),
    ]

    # share the colour scale with the observed target so the K rows can be read
    # directly against the observation (the purpose of the tuning).
    norm = Normalize(vmin=0.0, vmax=robust_vmax(rows[-1][1]))

    fig, axes = plt.subplots(4, 1, figsize=(4.0, 13.5), constrained_layout=True)
    im = None
    for ax, (label, field) in zip(axes, rows):
        im = show(ax, field, norm)
        ax.set_ylabel(label, fontsize=10)
    axes[0].set_title('Present-day basal melt', fontsize=12)
    cbar = fig.colorbar(im, ax=axes, shrink=0.6, location='right', extend='max')
    cbar.set_label('Basal melt (kg m$^{-2}$ yr$^{-1}$)')

    out = os.path.join(figure_dir, 'figure_presentday_meltmaps.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print('Saved:', out)


# ----------------------------------------------------- figure 2: ocean modelling
def plot_oceanmodelling():
    # find the Kmean index from one of the ensembles
    p1 = xr.open_dataset(ens_warm)['p1'].values
    i_mean = int(np.argmin(np.abs(p1 - p1.mean())))
    print('Ocean modelling at Kmean=%.2e (idx %d)' % (p1[i_mean], i_mean))

    fig, axes = plt.subplots(2, 4, figsize=(15.0, 8.0),
                             constrained_layout=True)

    # Each scenario (column) shares one colour scale so the model (top) and its
    # ocean-model target (bottom) are directly comparable.
    for col, (tag, label, ens_file, model_name, tgt_file) in enumerate(
            OM_SCENARIOS):
        model = _finite(
            xr.open_dataset(ens_file)['melt_rate']
              .sel(model=model_name).isel(p2=0, p1=i_mean))
        target = _finite(
            xr.open_dataset(os.path.join(om_target_dir, tgt_file))['melt_rate'])

        norm = Normalize(vmin=0.0, vmax=robust_vmax(model, target))
        im = show(axes[0, col], model, norm)
        show(axes[1, col], target, norm)
        axes[0, col].set_title('%s\n(%s)' % (tag, label), fontsize=11)

        cbar = fig.colorbar(im, ax=axes[:, col], shrink=0.55,
                            location='bottom', extend='max', pad=0.02)
        cbar.set_label('kg m$^{-2}$ yr$^{-1}$', fontsize=9)

    axes[0, 0].set_ylabel('Model (Kmean)', fontsize=11)
    axes[1, 0].set_ylabel('Target (ocean model)', fontsize=11)
    fig.suptitle('Ocean-modelling basal melt: model at Kmean (top) vs '
                 'ocean-model target (bottom)', fontsize=13)

    out = os.path.join(figure_dir, 'figure_oceanmodelling_meltmaps.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print('Saved:', out)


if __name__ == '__main__':
    plot_presentday()
    plot_oceanmodelling()
