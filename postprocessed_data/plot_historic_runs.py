"""
plot_historic_runs.py

Summary plots for the ISMIP6 Appendix-2 NetCDF files written by
hist_run_CESM_WACCM_1995_2014.m (step 2, WriteISMIP6_NetCDF) into
<CMIP_Model>/hist/ next to this script. Generic across model runs --
change CMIP_Model and EXP below to point at a different run.

Input NetCDF files: <script_dir>/<CMIP_Model>/hist/
Output figures:     <script_dir>/figures/<CMIP_Model>/<EXP>/

Produces:
  - timeseries_<CMIP_Model>_<EXP>.png : every scalar (time-only) variable,
    one subplot per variable (kg/kg s-1 fields shown in Gt/Gt a-1).
  - maps_<group>_<CMIP_Model>_<EXP>.png : one figure per MAP_GROUPS entry,
    each variable shown as start / end / (start-end) of run, with its own
    independent colorbar/min-max per panel.

Variables are discovered automatically from the NetCDF files in
<CMIP_Model>/hist/ (named <variable>_<IS>_<GROUP>_<MODEL>_<EXP>.nc) and
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
EXP        = 'historical'   # experiment label; used in filenames and as the output subfolder

script_dir = os.path.dirname(os.path.abspath(__file__))
indir  = os.path.join(script_dir, CMIP_Model, 'hist')
outdir = os.path.join(script_dir, 'figures', CMIP_Model, EXP)

YTS = 365.25 * 24 * 3600  # seconds per year, matches ISSM's md.constants.yts

# ISMIP6 short name -> originating ISSM field name (per
# hist_run_CESM_WACCM_1995_2014.m step 2's scalars table / native outputs),
# used only for nicer subplot titles; any variable not listed here just
# falls back to its bare ISMIP6 name, so this never needs updating to stay
# correct -- only to stay maximally descriptive.
ISSM_NAME = {
    'lim': 'IceVolumeScaled',
    'limnsw': 'IceVolumeAboveFloatationScaled',
    'iareagr': 'GroundedAreaScaled',
    'iareafl': 'FloatingAreaScaled',
    'tendacabf': 'TotalSmbScaled',
    'tendlibmassbf': 'TotalGroundedBmbScaled+TotalFloatingBmbScaled',
    'tendlibmassbffl': 'TotalFloatingBmbScaled',
    'tendlicalvf': 'IcefrontMassFluxLevelset',
    'tendligroundf': 'GroundinglineMassFlux',
}

# 2D map figure groups: (group_key, title, [variable names in plot order]).
# 'vel' is not a NetCDF variable -- it's derived on the fly as
# sqrt(xvelsurf^2 + yvelsurf^2) and plotted on a log color scale (ice
# velocities span orders of magnitude). Any 2D variable not listed in any
# group below still gets plotted, in a catch-all "other" figure.
MAP_GROUPS = [
    ('states', 'Fig1. State variables', ['lithk', 'orog', 'base', 'vel']),
    ('masks', 'Fig2. Masks', ['sftgif', 'sftgrf', 'sftflf']),
    ('fluxes', 'Fig3. Fluxes', ['acabf', 'libmassbfgr', 'libmassbffl', 'dlithkdt']),
    ('constant', 'Fig4. Constant fields', ['topg', 'litemptop', 'hfgeoubed']),
    ('boundaryfluxes', 'Fig5. Boundary fluxes', ['ligroundf', 'licalvf']),
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
    return ref_year + time.values / 365.25


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
            ax.plot(years, values, marker='o', markersize=3)
            issm_name = ISSM_NAME.get(varname)
            ax.set_title(f'{varname} ({issm_name})' if issm_name else varname)
            ax.set_xlabel('year')
            ax.set_ylabel(units)
            ax.grid(True, alpha=0.3)

    for ax in axes_flat[nvar:]:
        ax.axis('off')

    fig.suptitle(f'{CMIP_Model} {EXP} — time series', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f'Saved: {outpath}')


def load_2d(varname, field2d_vars):
    """Load a 2D (x, y, time) variable, with missing_value replaced by NaN.
    Explicitly transposed to (x, y, time) by dimension name -- the MATLAB
    writer actually stores these on disk as (time, y, x), so indexing the
    raw array by assumed axis position (instead of by name) silently picks
    the wrong axis."""
    with xr.open_dataset(field2d_vars[varname], decode_times=False) as ds:
        da = ds[varname].transpose('x', 'y', 'time')
        units = da.attrs.get('units', '')
        missing = da.attrs.get('missing_value', None)
        data = da.values.astype(float)  # guaranteed (x, y, time) by the transpose above
        if missing is not None:
            data = np.where(np.isclose(data, missing), np.nan, data)
        x = ds['x'].values
        y = ds['y'].values
    return data, x, y, units


def build_panel(varname, field2d_vars):
    """Return (varname, data, x, y, units, log_scale) for one map panel, or
    None if its inputs aren't available. 'vel' is derived from
    xvelsurf/yvelsurf rather than read directly."""
    if varname == 'vel':
        if 'xvelsurf' not in field2d_vars or 'yvelsurf' not in field2d_vars:
            print('Skipping vel: xvelsurf/yvelsurf not both available.')
            return None
        vx, x, y, vunits = load_2d('xvelsurf', field2d_vars)
        vy, _, _, _ = load_2d('yvelsurf', field2d_vars)
        vel = np.sqrt(vx**2 + vy**2)
        return ('vel', vel, x, y, vunits, True)
    if varname not in field2d_vars:
        print(f'Skipping {varname}: not found in {indir}.')
        return None
    data, x, y, units = load_2d(varname, field2d_vars)
    return (varname, data, x, y, units, False)


def plot_2d_panels(panels, title, outpath):
    """panels: list of (varname, data(x,y,t), x, y, units, log_scale).
    Each variable gets its own independent colorbar/min-max -- not shared
    across variables or groups."""
    nvar = len(panels)
    if nvar == 0:
        print(f'No variables available for "{title}", skipping.')
        return

    fig, axes = plt.subplots(nvar, 3, figsize=(12, 3.5 * nvar), squeeze=False)

    for row, (varname, data, x, y, units, log_scale) in enumerate(panels):
        extent = [x.min(), x.max(), y.min(), y.max()]
        if log_scale:
            positive = data[np.isfinite(data) & (data > 0)]
            vmin = np.percentile(positive, 2) if positive.size else 1e-3
            vmax = np.percentile(positive, 98) if positive.size else 1.0
            norm = LogNorm(vmin=max(vmin, 1e-6), vmax=max(vmax, vmin * 10))
            imshow_kwargs = dict(norm=norm)
        else:
            vmin = np.nanpercentile(data, 2)
            vmax = np.nanpercentile(data, 98)
            imshow_kwargs = dict(vmin=vmin, vmax=vmax)

        for col, (t_idx, label) in enumerate([(0, 'start'), (-1, 'end')]):
            ax = axes[row, col]
            im = ax.imshow(data[:, :, t_idx].T, origin='lower', extent=extent,
                            cmap='viridis', **imshow_kwargs)
            ax.set_title(f'{varname} ({label})')
            ax.set_aspect('equal')
            fig.colorbar(im, ax=ax, label=units, shrink=0.8)

        # Third panel: start - end, always linear/diverging (even for
        # log-scale variables like vel -- a difference can be negative, so
        # log-scaling it doesn't make sense), centered at zero with its own
        # independent colorbar.
        diff = data[:, :, 0] - data[:, :, -1]
        dmax = np.nanpercentile(np.abs(diff), 98)
        dmax = dmax if np.isfinite(dmax) and dmax > 0 else 1.0
        ax = axes[row, 2]
        im = ax.imshow(diff.T, origin='lower', extent=extent,
                        cmap='RdBu_r', vmin=-dmax, vmax=dmax)
        ax.set_title(f'{varname} (start - end)')
        ax.set_aspect('equal')
        fig.colorbar(im, ax=ax, label=units, shrink=0.8)

    fig.suptitle(f'{CMIP_Model} {EXP} — {title}', fontsize=14)
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

    plot_timeseries(scalar_vars, os.path.join(outdir, f'timeseries_{CMIP_Model}_{EXP}.png'))

    grouped = set()
    for key, title, varnames in MAP_GROUPS:
        panels = [p for p in (build_panel(v, field2d_vars) for v in varnames) if p is not None]
        for v in varnames:
            grouped.update(DERIVED_SOURCES.get(v, [v]))
        plot_2d_panels(panels, title, os.path.join(outdir, f'maps_{key}_{CMIP_Model}_{EXP}.png'))

    leftover = sorted(set(field2d_vars) - grouped)
    if leftover:
        panels = [p for p in (build_panel(v, field2d_vars) for v in leftover) if p is not None]
        plot_2d_panels(panels, 'Other 2D variables', os.path.join(outdir, f'maps_other_{CMIP_Model}_{EXP}.png'))


if __name__ == '__main__':
    main()
