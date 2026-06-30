# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

ISSM (Ice Sheet System Model) MATLAB workflow for initialising the Antarctic Ice
Sheet around 1995 and tuning basal-melt parameterisations for **ISMIP7 / meltMIP**.
Everything runs on the **NCI Gadi** HPC cluster under PBS, project `au88`, against a
local ISSM build at `$ISSM_DIR=/home/565/jb1863/trunk`. There is no build step for
this repo itself — the "code" is MATLAB driver scripts that call into ISSM and submit
solves to the cluster.

## Running things

All work is launched as PBS jobs from `init/`. Submit with `qsub`:

```bash
cd init
qsub submit_init2.sh      # init / spinup / inversion / relaxation pipeline (tuning_func.m)
qsub submit_meltMIP.sh    # basal-melt parameter ensemble pipeline (meltMip_ensemble.m)
```

Each `.sh` purges modules, loads the ISSM toolchain (openmpi/netcdf/hdf5/petsc/matlab
R2021b + matlab_licence/monash), then invokes MATLAB headless:

```bash
matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
  addpath('$ISSM_DIR/lib'); run('meltMip_ensemble($steps, $j, $loadonly)'), quit" > FuncMelt_$j.log
```

Logs land in `init/FuncInit2.log`, `init/FuncMelt_<j>.log`, and PBS `*.outlog`/`*.errlog`
(all gitignored). To debug a failed run, read the `FuncMelt_<j>.log` — ISSM consistency
errors surface there as `model not consistent: field '...'`.

## The step / organizer pattern (critical to understand)

Both `tuning_func.m` and `meltMip_ensemble.m` are single MATLAB functions that gate each
stage behind `if perform(org,'StageName')`. The `organizer` assigns each stage a **number
by its order of appearance** in the file. You select which stages run by passing
`steps=[...]` in the submit script; outputs are saved/loaded by stage name into `Models/`.

`meltMip_ensemble(steps, j, loadonly)` step map:

| # | stage | role |
|---|---|---|
| 1 | `Obs_clim_TF` | build present-day obs thermal-forcing `.mat` |
| 2 | `OceanModelling_clim_TF` | build Naughten/Mathiot cold+warm TF `.mat` |
| 3 | `ObsData_clim_TF` | build Dutrieux (Amundsen obs) TF `.mat` per year |
| 4 | `melt_run` | present-day melt run over K-ensemble |
| 5 | `melt_run_OceanModelling` | melt runs for Nw/Nc/Mw/Mc forcings |
| 6 | `melt_run_ObsData` | melt runs for obs years (2009, 2012) |
| 7–8 | `create_BMB_gD`, `create_BMB_gD_4km` | grid present-day BMB to 8km/4km NetCDF |
| 9–10 | `create_BMB_gD_OceanModelling`, `create_BMB_gD_ObsData` | grid ocean-model / obs BMB |

`tuning_func(steps, loadonly)` is the upstream init chain (`Param` → `InversionB` →
`InversionFriction1` → `Inversion2B` → `CollapseSSA` → `RunInitSSACollapse` →
`Relax_long` → `Relaxed` → `Assign_Basins` → ...). It produces
`init/Models/AIS_ISMIP7_Relaxed.mat`, which is the `inputmodel_relax` every meltMIP melt
run loads.

### `loadonly` semantics (submit vs. gather)

Each melt/inversion stage calls `md=solve(md,...,'loadonly',loadonly)` with
`md.settings.waitonlock=0` (async). This means runs are **two-phase**:

- `loadonly=0` → builds the model, **submits** the solve to Gadi, returns immediately
  (does not wait for results).
- `loadonly=1` → **gathers** the finished results back from the cluster and saves the
  output `.mat` into `Models/`.

So a full ensemble is: launch with `loadonly=0`, wait for PBS jobs to finish, then re-run
the same `steps` with `loadonly=1` to collect, then run the gridding steps.

### Bash gotcha in the submit scripts

The submit scripts contain several plain assignments like:

```bash
steps=[1,2,3]
loadonly=[0]
steps=[4,5,6]
loadonly=[1]
```

Bash keeps **only the last** `steps=` and `loadonly=` — the earlier lines are dead. The
comment blocks describe a manual, one-stage-at-a-time workflow; you must edit these lines
(or split into separate `matlab` invocations) to actually run a given stage. TF builders
(steps 1–3) take no `j` and only need to run **once**; melt runs (4–6) loop over
`j=1..120` (the K parameter index).

## Data flow

```
tuning_func.m ──> Models/AIS_ISMIP7_Relaxed.mat  (relaxed 1995 state)
       │
meltMip_ensemble.m (loads Relaxed as inputmodel_relax)
   TF build (1-3): raw_data/.../ocean/*.nc ──InterpFromGridToMesh──> preprocessed_data/Ocean/Clim/*.mat
   melt runs (4-6): per K (gamma_0 = K_data(j)/gT_to_K), ISMIP6 local quadratic basal melt
   gridding (7-10): gridData.m ──> Models/ModelNC/BMB*/*_gridData_8km.nc  (holds rho_ice * BMB[m/s])
       │
creating_bmbMIP_*.py  (multiply by yearsinsec ──> kg/m2/a; assemble xarray ensembles)
   creating_bmbMIP_ensmble.py              ──> pd_ensemble  (present-day)
   creating_bmbMIP_OceanModelling_ensemble.py ──> cold_ensemble.nc / warm_ensemble.nc
   creating_bmbMIP_ObsData_ensemble.py     ──> obs_ensemble.nc
       │
BMB_tuning_python/parameter_selection_quadratic_example.ipynb  (ISMIP7 parameter selection)
```

The three `creating_bmbMIP_*.py` scripts must match the MATLAB output filenames, tags
(`Nw/Nc/Mw/Mc`), K range (`np.arange(0.25e-5, 3.025e-4, 0.25e-5)`, 120 values), and obs
years exactly — they index the NetCDFs by the MATLAB `j` (1-based). The xarray structures
they emit (`pd_ensemble` / `cold_ensemble` / `warm_ensemble` / `obs_ensemble`) are the
required inputs to the ISMIP7 tuning toolbox; see `BMB_tuning_python/README.md` for the
exact dims/coords each must have.

## Key conventions and constraints

- **Thermal forcing must be ≥ 0.** ISSM's `basalforcingsismip6` checkconsistency
  (`$ISSM_DIR/src/m/classes/basalforcingsismip6.m`) requires every depth slice
  `basalforcings.tf{i} >= 0`. Cold ocean-model and obs TF fields contain small negative
  (cold-cavity) values, so the TF builders clamp with `max(temp_tfdata, 0)` after
  interpolation. Warm forcings have no negatives. If you regenerate TF `.mat` files you
  must rebuild (rerun the relevant step) before the melt runs can pass consistency.
- The `tf` cell array is `cell(1,1,nDepths)` — **3rd dimension is depth**. Each slice
  is written `tf(:,:,i)={[tf_at_vertices ; t]}` with `t` the (single) time row. Writing
  `tf(:,:)=...` instead silently collapses all depths (historical bug — keep the `i`).
- Gridding targets the fixed ISMIP6 8km grid read from
  `/home/565/jb1863/ismip6_2300/masks/af2_el_ismip6_ant_8km.nc` (see `scripts/gridData.m`).
- Cluster config lives in `scripts/set_cluster.m` (`gadi`: np=48, mem=190, queue=normal,
  project=au88).
- `raw_data/`, `preprocessed_data/`, and `init/Models/` are gitignored (large data /
  HPC outputs). Only the driver scripts (`init/*.m`, `init/*.py`, `init/*.sh`,
  `init/scripts/`, `init/Par/`) and the python tuning toolbox are tracked.

## SMB class unit contracts

`SMBgradients.smbref` and `SMBforcing.mass_balance` expect **different units**:

| Class | Field | Expected units | ISSM internal conversion |
|---|---|---|---|
| `SMBgradients` | `smbref` | mm w.e. yr⁻¹ | applies `/1000 * rho_w/rho_i` internally (from C++ `smb_core`) |
| `SMBforcing` | `mass_balance` | m ice yr⁻¹ | none — used directly |

`smb_forcing` built by `ProjSMB` / `hist_run_tune_CESM_WACCM.m` is in **mm w.e. yr⁻¹**.
When switching from `SMBgradients` to `SMBforcing` (e.g. a noSEF variant), convert first:

```matlab
rho_w = 1000.0;
rho_i = md.materials.rho_ice;
smb_mice = smb_forcing;
smb_mice(1:end-1,:) = smb_forcing(1:end-1,:) / 1000.0 * (rho_w / rho_i);
md.smb.mass_balance = smb_mice;
```

The last row of `smb_forcing` is the time axis (in years) and must **not** be converted.
Passing mm w.e. yr⁻¹ directly to `mass_balance` gives values ~917× too large, causing
immediate solver failure ("Recovery solver failed" on all ranks) within the first few
simulated months.

### `SMBgradients.href` reference surface

`href` is the surface elevation at which `smbref` applies zero lapse-rate correction.
It **must always be the 1995 relaxed surface** (`AIS_ISMIP7_Relaxed_CESM_WACCM.mat`),
because `smbref` is built from the RACMO 1995 climatology.

```matlab
md_relax  = loadmodel([init_dir 'Models/AIS_ISMIP7_Relaxed_CESM_WACCM.mat']);
md.smb.href = [md_relax.geometry.surface ; 1995];
clear md_relax
```

Do **not** use `md.geometry.surface` from the segment start state (2015, 2151, …).
Doing so shifts the lapse-rate neutral point forward in time, suppressing or inflating the
accumulated SEF signal. The time value in the single-snapshot `href` column is not used
by the solver (only matters for multi-snapshot interpolation), but should be set to 1995
for documentation consistency.

## Creating variant scripts

When creating a new script that is a variant of an existing one (e.g. a `_noSEF` version
of a projection script), **always read the corresponding section of the original script
before writing** — never generate field names, output lists, configuration blocks, or
function bodies from memory or general knowledge. Specifically:

- Use the Read tool on the original file first, then copy shared sections verbatim.
- For `ismip6_outputs`, `requested_outputs`, helper functions (`write_proj_netcdf_ismip6`,
  `write_ismip6_2d`, `write_ismip6_scalar`, `flux_along_contour_2d`), and any
  configuration block that must match the original exactly, read the source and copy —
  do not paraphrase or reconstruct.
- After creating the new file, diff the shared sections against the original to verify
  they match: `grep -A30 'ismip6_outputs' original.m` vs the new file.

Failure to do this has previously introduced wrong field names (`TotalCalvingFlux`,
`BasalforcingsMassBalance`) and missing fields (`GroundinglineMassFlux`,
`IcefrontMassFluxLevelset`, `TotalSmbScaled`) that caused silent data errors or
runtime failures.

## Layout

- `init/tuning_func.m` — init/spinup/inversion/relaxation pipeline.
- `init/meltMip_ensemble.m` — meltMIP TF-build / melt-run / gridding pipeline.
- `init/iter_cycle_thickTHW.m` — Thwaites thickness iteration driver.
- `init/scripts/` — helpers: `gridData.m`, `gridData_4km.m`, `set_cluster.m`,
  `Obs_ISMIP7melt.m`, plotting (`plot_*.m`, `melt_colormap_*`).
- `init/Par/Antartica_1995_Thwbedmap.par` — ISSM parameterisation file (mesh → model).
- `init/BMB_tuning_python/` — ISMIP7 parameter-selection toolbox + example notebooks.
- `init/Models/` — ISSM model `.mat` outputs and gridded NetCDF (`Models/ModelNC/BMB*`).
