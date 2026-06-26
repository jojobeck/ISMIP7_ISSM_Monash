"""
plot_projection_runs.py

Summary plots for the ISMIP6 Appendix-2 NetCDF files written by
proj_run_CESM_WACCM_ssp585_2015_2300.m (step 8, WriteISMIP6_NetCDF) into
<CMIP_Model>/<SCENARIO>/ next to this script. Generic across model runs --
change CMIP_Model and SCENARIO below to point at a different run.

Input NetCDF files: <script_dir>/<CMIP_Model>/<SCENARIO>/
Output figures:     <script_dir>/figures/<CMIP_Model>/<SCENARIO>/

Produces:
  - timeseries_<CMIP_Model>_<SCENARIO>.png : every scalar (time-only) variable,
    one subplot per variable (kg/kg s-1 fields shown in Gt/Gt a-1).
  - maps_<group>_<CMIP_Model>_<SCENARIO>.png : one figure per MAP_GROUPS entry,
    each variable shown as start / end / (start-end) of run, with its own
    independent colorbar/min-max per panel.

Variables are discovered automatically from the NetCDF files in
<CMIP_Model>/<SCENARIO>/ (named <variable>_<IS>_<GROUP>_<MODEL>_<EXP>.nc) and
classified as scalar/2D by their dimensions; any 2D variable not listed in
MAP_GROUPS still gets plotted, in a catch-all "other" figure.
"""
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.colors import LogNorm

CMIP_Model = 'CESM2-WACCM'  # CMIP6 model forcing label; used in filenames and folder structure
SCENARIO   = 'ssp585'        # scenario label; change to 'ssp126' etc. as needed

script_dir = os.path.dirname(os.path.abspath(__file__))
indir  = os.path.join(script_dir, CMIP_Model, SCENARIO)
outdir = os.path.join(script_dir, 'figures', CMIP_Model, SCENARIO)

YTS = 365 * 24 * 3600  # seconds per year (365-day calendar, matches ISSM output)

# ISMIP6 short name -> originating ISSM field name (per
# proj_run_CESM_WACCM_ssp585_2015_2300.m step 8's scalars table / native outputs),
# used only for nicer subplot titles; any variable not listed here just
# falls back to its bare ISMIP6 name.
ISSM_NAME = {
    'lim':             'IceVolumeScaled',
    'limnsw':          'IceVolumeAboveFloatationScaled',
    'iareagr':         'GroundedAreaScaled',
    'iareafl':         'FloatingAreaScaled',
    'tendacabf':       'TotalSmbScaled',
    'tendlibmassbf':   'TotalGroundedBmbScaled+TotalFloatingBmbScaled',
    'tendlibmassbffl': 'TotalFloatingBmbScaled',
    'tendlicalvf':     'IcefrontMassFluxLevelset',
    'tendligroundf':   'GroundinglineMassFlux',
}

# 2D map figure groups: (group_key, title, [variable names in plot order]).
# 'vel' is not a NetCDF variable -- it's derived on the fly as
# sqrt(xvelsurf^2 + yvelsurf^2) and plotted on a log color scale.
# Any 2D variable not listed in any group below still gets plotted,
# in a catch-all "other" figure.
MAP_GROUPS = [
    ('states',         'Fig1. State variables',   ['lithk', 'orog', 'base', 'vel']),
    ('masks',          'Fig2. Masks',              ['sftgif', 'sftgrf', 'sftflf']),
    ('fluxes',         'Fig3. Fluxes',             ['acabf', 'libmassbfgr', 'libmassbffl', 'dlithkdt']),
    ('constant',       'Fig4. Constant fields',    ['topg', 'litemptop', 'hfgeoubed']),
    ('boundaryfluxes', 'Fig5. Boundary fluxes',    ['ligroundf', 'licalvf']),
]

# Derived panel name -> the real NetCDF variables it's computed from, so
# the catch-all "other" figure doesn't redundantly re-plot the raw inputs.
DERIVED_SOURCES = {'vel': ['xvelsurf', 'yvelsurf']}


def to_gt(data, units):
    """Convert kg -> Gt or kg s-1 -> Gt a-1 where applicable; otherwise pass through."""
    if units == 'kg':
        return data / 1e12, 'Gt'
    elif units == 'kg s-1':
        return data * YTS / 1e12, 'Gt a-1'
    else:
        return data, units


def find_variable_files(indir):
    """Map <varname> -> filepath for every <varname>_*.nc file in indir."""
    varfiles = {}
    for f in sorted(glob.glob(os.path.join(indir, '*.nc'))):
        varname = os.path.basename(f).split('_')[0]
        varfiles[varname] = f
    return varfiles


def classify(varfiles):
    """Split varfiles into scalar (time-only) and 2D (x, y, time) groups."""
    scalar_vars, field2d_vars = {}, {}
    for varname, f in varfiles.items():
        with xr.open_dataset(f, decode_times=False) as ds:
            dims = set(ds[varname].dims)
        if dims == {'time'}:
            scalar_vars[varname] = f
        elif {'x', 'y', 'time'}.issubset(dims):
            field2d_vars[varname] = f
        else:
            print(f'Skipping {varname}: unrecognized dims {dims}')
    return scalar_vars, field2d_vars


def time_to_year(ds):
    """Convert the CF 'days since YYYY-MM-DD' time coordinate to decimal years."""
    time = ds['time']
    units = time.attrs.get('units', '')
    m = re.search(r'since\s+(\d{4})-(\d{2})-(\d{2})', units)
    ref_year = int(m.group(1)) if m else 0
    return ref_year + time.values / 365


def plot_timeseries(scalar_vars, outpath):
    nvar = len(scalar_vars)
    if nvar == 0:
        print('No scalar/time-series variables found, skipping timeseries figure.')
        return

    ncols = 3
    nrows = int(np.ceil(nvar / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    for ax, (varname, f) in zip(axes_flat, sorted(scalar_vars.items())):
        with xr.open_dataset(f, decode_times=False) as ds:
            da = ds[varname]
            years = time_to_year(ds)
            values, units = to_gt(da.values, da.attrs.get('units', ''))
            ax.plot(years, values, linewidth=0.8)
            issm_name = ISSM_NAME.get(varname)
            ax.set_title(f'{varname} ({issm_name})' if issm_name else varname)
            ax.set_xlabel('year')
            ax.set_ylabel(units)
            ax.grid(True, alpha=0.3)

    for ax in axes_flat[nvar:]:
        ax.axis('off')

    fig.suptitle(f'{CMIP_Model} {SCENARIO} — time series', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f'Saved: {outpath}')


def load_2d_endpoints(varname, field2d_vars):
    """Load only the first and last time step of a 2D variable.
    Avoids loading all 285 time steps (~450 MB per variable) into memory.
    Returns (first, last, x, y, units) where first/last are (x, y) arrays."""
    with xr.open_dataset(field2d_vars[varname], decode_times=False) as ds:
        da = ds[varname]
        units = da.attrs.get('units', '')
        missing = da.attrs.get('missing_value', None)
        first = da.isel(time=0).values.astype(float)   # (y, x) on disk
        last  = da.isel(time=-1).values.astype(float)
        if missing is not None:
            first = np.where(np.isclose(first, missing), np.nan, first)
            last  = np.where(np.isclose(last,  missing), np.nan, last)
        x = ds['x'].values
        y = ds['y'].values
    # NetCDF stored as (time, y, x); transpose to (x, y) for imshow
    return first.T, last.T, x, y, units


def build_panel(varname, field2d_vars):
    """Return (varname, first, last, x, y, units, log_scale) for one map panel,
    or None if inputs aren't available. 'vel' is derived from xvelsurf/yvelsurf."""
    if varname == 'vel':
        if 'xvelsurf' not in field2d_vars or 'yvelsurf' not in field2d_vars:
            print('Skipping vel: xvelsurf/yvelsurf not both available.')
            return None
        vx0, vx1, x, y, vunits = load_2d_endpoints('xvelsurf', field2d_vars)
        vy0, vy1, _,  _, _     = load_2d_endpoints('yvelsurf', field2d_vars)
        return ('vel', np.sqrt(vx0**2 + vy0**2), np.sqrt(vx1**2 + vy1**2),
                x, y, vunits, True)
    if varname not in field2d_vars:
        print(f'Skipping {varname}: not found in {indir}.')
        return None
    first, last, x, y, units = load_2d_endpoints(varname, field2d_vars)
    return (varname, first, last, x, y, units, False)


def plot_2d_panels(panels, title, outpath):
    """panels: list of (varname, first(x,y), last(x,y), x, y, units, log_scale).
    Only the first and last time steps are held in memory — no full 3D array."""
    nvar = len(panels)
    if nvar == 0:
        print(f'No variables available for "{title}", skipping.')
        return

    fig, axes = plt.subplots(nvar, 3, figsize=(12, 3.5 * nvar), squeeze=False)

    for row, (varname, first, last, x, y, units, log_scale) in enumerate(panels):
        extent = [x.min(), x.max(), y.min(), y.max()]
        diff = first - last

        if log_scale:
            positive = first[np.isfinite(first) & (first > 0)]
            vmin = np.percentile(positive, 2) if positive.size else 1e-3
            vmax = np.percentile(positive, 98) if positive.size else 1.0
            norm = LogNorm(vmin=max(vmin, 1e-6), vmax=max(vmax, vmin * 10))
            imshow_kwargs = dict(norm=norm)
        else:
            all_vals = np.concatenate([first[np.isfinite(first)], last[np.isfinite(last)]])
            vmin = np.percentile(all_vals, 2) if all_vals.size else 0
            vmax = np.percentile(all_vals, 98) if all_vals.size else 1
            imshow_kwargs = dict(vmin=vmin, vmax=vmax)

        for ax, data_slice, label in [
            (axes[row, 0], first, '2015'),
            (axes[row, 1], last,  '2299'),
        ]:
            im = ax.imshow(data_slice.T, origin='lower', extent=extent,
                           cmap='viridis', **imshow_kwargs)
            ax.set_title(f'{varname} ({label})')
            ax.set_aspect('equal')
            fig.colorbar(im, ax=ax, label=units, shrink=0.8)

        dmax = np.nanpercentile(np.abs(diff), 98)
        dmax = dmax if np.isfinite(dmax) and dmax > 0 else 1.0
        ax = axes[row, 2]
        im = ax.imshow(diff.T, origin='lower', extent=extent,
                       cmap='RdBu_r', vmin=-dmax, vmax=dmax)
        ax.set_title(f'{varname} (2015 − 2299)')
        ax.set_aspect('equal')
        fig.colorbar(im, ax=ax, label=units, shrink=0.8)

    fig.suptitle(f'{CMIP_Model} {SCENARIO} — {title}', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f'Saved: {outpath}')


def main():
    os.makedirs(outdir, exist_ok=True)

    varfiles = find_variable_files(indir)
    if not varfiles:
        raise FileNotFoundError(f'No NetCDF files found in {indir}')

    scalar_vars, field2d_vars = classify(varfiles)

    plot_timeseries(scalar_vars, os.path.join(outdir, f'timeseries_{CMIP_Model}_{SCENARIO}.png'))

    grouped = set()
    for key, title, varnames in MAP_GROUPS:
        panels = [p for p in (build_panel(v, field2d_vars) for v in varnames) if p is not None]
        for v in varnames:
            grouped.update(DERIVED_SOURCES.get(v, [v]))
        plot_2d_panels(panels, title, os.path.join(outdir, f'maps_{key}_{CMIP_Model}_{SCENARIO}.png'))

    leftover = sorted(set(field2d_vars) - grouped)
    if leftover:
        panels = [p for p in (build_panel(v, field2d_vars) for v in leftover) if p is not None]
        plot_2d_panels(panels, 'Other 2D variables', os.path.join(outdir, f'maps_other_{CMIP_Model}_{SCENARIO}.png'))


if __name__ == '__main__':
    main()
