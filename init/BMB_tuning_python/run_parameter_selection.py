"""
ISMIP7 parameter selection for the quadratic local basal-melt parameterisation.

Loads the pre-assembled ensemble NetCDFs (output of creating_bmbMIP_*.py),
calculates the four objective-function terms, samples the weighted objective
function, and saves the K-distribution plot to ../figures/.

Run from BMB_tuning_python/:  python run_parameter_selection.py
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.io import savemat

from parameter_selection_toolbox import (
    calculate_term1,
    calculate_term2,
    calculate_term3,
    calculate_term4,
    calculate_objective_function,
)

# ---------------------------------------------------------------------- paths
data_path = os.path.join('..', '..', 'raw_data', 'ISMIP7', 'AIS')
param_path = os.path.join(data_path, 'parameterisations','ocean')

ens_pd   = os.path.join('..', 'Models', 'ModelNC', 'BMBPresentDay',     'MeltObs_ensemble.nc')
ens_warm = os.path.join('..', 'Models', 'ModelNC', 'BMBOceanModelling', 'warm_ensemble.nc')
ens_cold = os.path.join('..', 'Models', 'ModelNC', 'BMBOceanModelling', 'cold_ensemble.nc')
ens_obs  = os.path.join('..', 'Models', 'ModelNC', 'BMBObsData',        'obs_ensemble.nc')

figure_dir = os.path.join('..', 'figures')
os.makedirs(figure_dir, exist_ok=True)

# ------------------------------------------------------------------ constants
reso       = 8000                  # m
cvt_m      = reso**2 / 1e12       # grid-cell area → Gt/a
model_reso = '8km'
sample_size = 10_000

# --------------------------------------------------------- masks and targets
print('Loading masks and targets...')

basins_m = xr.load_dataset(
    os.path.join(data_path, 'parameterisations', 'ocean', 'imbie2',
                 f'basin_numbers_ismip{model_reso}_v2.nc')
).rename({'basinNumber': 'basins'})

bfrn_m = xr.load_dataset(
    os.path.join(param_path, 'bfrns', f'BFRN_ismip{model_reso}_v2.nc')
)

mask_m = xr.load_dataset(
    os.path.join(param_path, 'floatingmasks', f'floatingmask_ismip{model_reso}.nc')
).mask

nBasins = int(basins_m.basins.max())

MeltDataImbie = pd.read_csv(
    os.path.join(param_path, 'meltobs', 'Melt_Paolo_Davison_Adusumilli_imbie2.csv'),
    index_col=0,
)
buttressing_target = xr.load_dataset(
    os.path.join(param_path, 'meltobs', 'melt_target_term2_v3.nc')
)
cold_target = xr.load_dataset(
    os.path.join(param_path, 'ocean_modelling_data', 'melt_cold_target_term3.nc')
)
warm_target = xr.load_dataset(
    os.path.join(param_path, 'ocean_modelling_data', 'melt_warm_target_term3.nc')
)
t4_obs = xr.load_dataset(
    os.path.join(param_path, 'ocean_observations_data', 'melt_observations_target_term4.nc')
)

# ------------------------------------------------------- load ensembles
print('Loading ensembles...')
pd_ensemble   = xr.open_dataset(ens_pd)
cold_ensemble = xr.open_dataset(ens_cold)
warm_ensemble = xr.open_dataset(ens_warm)
obs_ensemble  = xr.open_dataset(ens_obs)   # years: [2009, 2012]

# ------------------------------------------------- region labels (PIG / Dotson)
shelves = xr.load_dataset(
    os.path.join(param_path, 'shelfmask', 'shelf_mask_ismip8km.nc')
).shelf_mask.isel(time=0)

PIG_ID, DOTSON_ID = 110, 97
pig    = (shelves == PIG_ID) & (shelves.x > -1.625e6)   # restrict to main trunk
dotson = (shelves == DOTSON_ID)
pigdotson_mask = pig.astype(int) + dotson.astype(int) * 2
label_map = {1: 'pig', 2: 'dotson'}
labels = np.vectorize(label_map.get)(pigdotson_mask.values)
region_label_m = xr.DataArray(
    labels, dims=pigdotson_mask.dims, coords=pigdotson_mask.coords
)

# --------------------------------------------------------- calculate terms
print('Calculating terms...')

t1_model, t1_obs_mean, t1_obs_sigma = calculate_term1(
    pd_ensemble, mask_m, basins_m['basins'], nBasins, cvt_m, MeltDataImbie
)
t2_model, t2_obs_mean, t2_obs_sigma = calculate_term2(
    pd_ensemble, mask_m, bfrn_m, cvt_m, buttressing_target
)
t3_model, t3_obs_mean, t3_obs_sigma = calculate_term3(
    cold_ensemble, warm_ensemble, cold_target, warm_target, mask_m, basins_m['basins']
)
t4_model, t4_obs_mean, t4_obs_sigma = calculate_term4(
    obs_ensemble, region_label_m, t4_obs, mask_m, cvt_m
)

# --------------------------------------------------------------- weights
# Term 1: equal weight per basin (all basins)
t1_weights = xr.DataArray(
    np.ones(nBasins + 1),
    dims=['basins'],
    coords={'basins': t1_model.basins.values},
)

# Term 2: weighted by BFRN medians
t2_weights = xr.DataArray(
    (bfrn_m['BFRN_medians'] / bfrn_m['BFRN_median']).values,
    dims=['BFRN_bins'],
    coords={'BFRN_bins': t2_obs_mean.BFRN_bins.values},
)

# Term 3: only Mathiot and Naughten AIS-1, equal weight per basin
# shape must match t3_obs_mean (the target file may have more models than the ensemble)
t3_weights = xr.DataArray(
    np.ones(t3_obs_mean.shape),
    dims=['model', 'basins'],
    coords={'model': t3_obs_mean.model.values, 'basins': t3_obs_mean.basins.values},
)
t3_weights = t3_weights.where(
    (t3_weights.model == 'mathiot') | (t3_weights.model == 'naughten_ais_1'), other=0
)

# Term 4: PIG only (obs_ensemble contains 2009 and 2012)
t4_weights = xr.DataArray(
    np.ones(t4_model.isel(p1=0, p2=0).shape),
    dims=['year', 'region'],
    coords={'year': t4_obs_mean.year.values, 'region': t4_obs_mean.region.values},
)
t4_weights = t4_weights.where(t4_weights.region == 'pig', other=0)

# -------------------------------------------------- sample objective function
print(f'Sampling objective function (n={sample_size:,})...')
min_p1, min_p2 = calculate_objective_function(
    sample_size, reso,
    t1_model, t1_obs_mean, t1_obs_sigma, t1_weights,
    t2_model, t2_obs_mean, t2_obs_sigma, t2_weights,
    t3_model, t3_obs_mean, t3_obs_sigma, t3_weights,
    t4_model, t4_obs_mean, t4_obs_sigma, t4_weights,
)

# ----------------------------------------------------------------- results
p1_vals   = pd_ensemble.p1.values
bin_edges = np.append(p1_vals[0] * 0.5, p1_vals + 1e-7)
counts, _ = np.histogram(min_p1, bins=bin_edges)
positions = np.arange(len(p1_vals))

K_5th  = np.percentile(min_p1, 5)
K_50th = np.median(min_p1)
K_95th = np.percentile(min_p1, 95)
K_mode = p1_vals[np.argmax(counts)]

peaks, _ = find_peaks(counts)
modes = p1_vals[peaks]

print('\nResults:')
print(f'  mode   = {K_mode:.4e}')
print(f'  median = {K_50th:.4e}')
print(f'  5th    = {K_5th:.4e}')
print(f'  95th   = {K_95th:.4e}')
if len(modes) > 1:
    print(f'  additional modes = {modes}')

# --------------------------------------------------------------- plot
fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(positions, counts / sample_size, width=0.8)

highlight_vals = [(K_5th, '5th'), (K_50th, '50th'), (K_95th, '95th')]
for val, lab in highlight_vals:
    idx = int(np.digitize(val, bin_edges)) - 1
    ax.patches[idx].set_facecolor('red')
    ax.annotate(
        lab,
        xy=(idx, counts[idx] / sample_size),
        xytext=(idx, (counts[idx] + sample_size / 100) / sample_size),
        ha='center', va='center',
        arrowprops=dict(arrowstyle='->', color='red'),
    )

idx_mode = int(np.digitize(K_mode, bin_edges)) - 1
ax.patches[idx_mode].set_edgecolor('orange')
ax.annotate(
    'mode',
    xy=(idx_mode, counts[idx_mode] / sample_size),
    xytext=(idx_mode, (counts[idx_mode] + sample_size / 100) / sample_size),
    ha='center', va='center',
    arrowprops=dict(arrowstyle='->', color='orange'),
)

tick_labels = np.round(p1_vals * 1e5, 1)
ax.set_xticks(positions[3::4], tick_labels[3::4], rotation=90)
ax.set_xlabel('K (×10⁻⁵)')
ax.set_ylabel('Relative frequency')
ax.set_ylim(0, ax.get_ylim()[1] * 1.2)
fig.tight_layout()

out = os.path.join(figure_dir, 'parameter_selection_K_distribution.png')
fig.savefig(out, dpi=150)
plt.close(fig)
print(f'Saved: {out}')

# ================================================================ validation
print('\n--- Validation against observation targets ---')

K_compare = {
    '5th':    K_5th,
    'median': K_50th,
    '95th':   K_95th,
    'mode':   K_mode,
}
clrs = {'5th': '#2166AC', 'median': '#4DAC26', '95th': '#D01C8B', 'mode': 'orange'}

# ---- AIS-wide melt summary (Term 1)
obs_total = float(np.nansum(t1_obs_mean))
rows = []
for label, K in K_compare.items():
    mod = t1_model.isel(p2=0).sel(p1=K, method='nearest')
    mod_total = float(mod.sum(skipna=True))
    rows.append({
        'K value':         f'{K:.3e}',
        'Modelled (Gt/a)': round(mod_total),
        'Observed (Gt/a)': round(obs_total),
        'Diff (Gt/a)':     round(mod_total - obs_total),
    })
df_ais = pd.DataFrame(rows, index=K_compare.keys())
print('\nAIS-wide melt (Term 1 — Paolo/Davison/Adusumilli):')
print(df_ais.to_string())

out_ais = os.path.join(figure_dir, 'parameter_selection_AIS_summary.txt')
with open(out_ais, 'w') as f:
    f.write('AIS-wide melt (Term 1 — Paolo/Davison/Adusumilli)\n')
    f.write(df_ais.to_string())
    f.write('\n')
print(f'Saved: {out_ais}')

# ---- basin-by-basin Term 1 table with within-uncertainty check
basins_x = t1_model.basins.values
df_t1 = pd.DataFrame({'obs_mean': t1_obs_mean, 'obs_sigma': t1_obs_sigma},
                     index=basins_x)
for label, K in K_compare.items():
    mod = t1_model.isel(p2=0).sel(p1=K, method='nearest').values
    df_t1[f'mod_{label}']    = mod
    df_t1[f'diff_{label}']   = mod - t1_obs_mean
    df_t1[f'within_{label}'] = (
        (mod >= t1_obs_mean - t1_obs_sigma) &
        (mod <= t1_obs_mean + t1_obs_sigma)
    ).astype(int)

# sum row: count of basins within uncertainty for each K
sum_row = {col: (df_t1[col].sum() if col.startswith('within_') else '')
           for col in df_t1.columns}
sum_row['obs_mean']  = 'SUM within σ'
sum_row['obs_sigma'] = ''
df_t1_out = pd.concat([df_t1, pd.DataFrame(sum_row, index=['TOTAL'])])

print('\nBasin-by-basin melt (Gt/a) + within-obs-uncertainty check (1=yes):')
print(df_t1_out.round(1).to_string())

out_t1 = os.path.join(figure_dir, 'parameter_selection_basin_comparison.txt')
with open(out_t1, 'w') as f:
    f.write('Basin-by-basin melt comparison (Gt/a)\n')
    f.write('within_* = 1 if modelled melt is within obs_mean ± obs_sigma\n\n')
    f.write(df_t1_out.round(1).to_string())
    f.write('\n')
print(f'Saved: {out_t1}')

# ---- validation figure: one panel per term
# band = K_5th..K_95th envelope, solid line = K_50th (median), dashed = K_mode
band_color = '#2166AC'
mode_color = 'orange'


def _band_median_mode(model_da):
    """K_5th/K_50th/K_95th/K_mode slices of a (..., p1, p2) DataArray."""
    med  = model_da.isel(p2=0).sel(p1=K_50th, method='nearest').values
    p05  = model_da.isel(p2=0).sel(p1=K_5th,  method='nearest').values
    p95  = model_da.isel(p2=0).sel(p1=K_95th, method='nearest').values
    mode = model_da.isel(p2=0).sel(p1=K_mode, method='nearest').values
    return med, p05, p95, mode


def _plot_term(ax, x, obs, sig, model_da, xlabel, ylabel, title):
    med, p05, p95, mode = _band_median_mode(model_da)
    ax.fill_between(x, p05, p95, color=band_color, alpha=0.3, label='5th–95th')
    ax.plot(x, med, color=band_color, label='median')
    ax.plot(x, mode, color=mode_color, linestyle='--', label='mode')
    ax.errorbar(x, obs, yerr=sig, fmt='o', color='black', capsize=4,
                label='Target', zorder=5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)


fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)

# [0,0] Term 1: basin-integrated melt
_plot_term(
    axes[0, 0], basins_x, t1_obs_mean, t1_obs_sigma, t1_model,
    'Basin', 'Melt rate (Gt/a)',
    'Term 1: basin-integrated melt vs Paolo/Davison/Adusumilli',
)

# [0,1] Term 2: buttressing-bin melt
bins2 = t2_obs_mean.BFRN_bins.values
_plot_term(
    axes[0, 1], bins2, t2_obs_mean.values, t2_obs_sigma.values, t2_model,
    'BFRN bin', 'Melt rate (Gt/a)',
    'Term 2: buttressing-bin melt vs target',
)

# [1,0] Term 3: warm−cold anomaly, model+basin stacked (only Mathiot and
# Naughten AIS-1 are weighted in the objective function, see t3_weights above)
t3_sel = (
    t3_model.sel(model=['mathiot', 'naughten_ais_1'])
    .stack(modelbasin=('model', 'basins'))
    .dropna(dim='modelbasin')
)
obs3 = t3_obs_mean.sel(model=['mathiot', 'naughten_ais_1']).stack(
    modelbasin=('model', 'basins')).sel(modelbasin=t3_sel.modelbasin)
sig3 = t3_obs_sigma.sel(model=['mathiot', 'naughten_ais_1']).stack(
    modelbasin=('model', 'basins')).sel(modelbasin=t3_sel.modelbasin)
x3 = np.arange(len(t3_sel.modelbasin))
_plot_term(
    axes[1, 0], x3, obs3.values, sig3.values, t3_sel,
    None, 'Warm − cold melt (kg m⁻² yr⁻¹)',
    'Term 3: warm − cold anomaly (Mathiot, Naughten)',
)
axes[1, 0].axhline(0, color='grey', linewidth=0.5, linestyle='--')
axes[1, 0].set_xticks(x3)
axes[1, 0].set_xticklabels(
    [f'{m}-{b}' for m, b in t3_sel.modelbasin.values], rotation=90, fontsize=7)

# [1,1] Term 4: PIG/Dotson observed years, region+year stacked
t4_sel = t4_model.stack(region_year=('region', 'year')).dropna(dim='region_year')
obs4 = t4_obs_mean.stack(region_year=('region', 'year')).sel(
    region_year=t4_sel.region_year)
sig4 = t4_obs_sigma.stack(region_year=('region', 'year')).sel(
    region_year=t4_sel.region_year)
x4 = np.arange(len(t4_sel.region_year))
_plot_term(
    axes[1, 1], x4, obs4.values, sig4.values, t4_sel,
    None, 'Melt rate (Gt/a)',
    'Term 4: PIG/Dotson observed-year melt vs target',
)
axes[1, 1].set_xticks(x4)
axes[1, 1].set_xticklabels(
    [f'{r}-{y}' for r, y in t4_sel.region_year.values], rotation=90, fontsize=8)

out_val = os.path.join(figure_dir, 'parameter_selection_validation.png')
fig.savefig(out_val, dpi=150)
plt.close(fig)
print(f'Saved: {out_val}')

# ============================================= save K values for MATLAB step
K_mat = os.path.join('..', '..', 'preprocessed_data', 'Ocean', 'K_selected.mat')
savemat(K_mat, {
    'K_mode': K_mode,
    'K_5th':  K_5th,
    'K_50th': K_50th,
    'K_95th': K_95th,
})
print(f'Saved: {K_mat}  (run meltMip_ensemble step 16 to write gamma0_local.mat)')
