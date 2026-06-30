# ISMIP7 ISSM Monash — AIS projection workflow

Antarctic Ice Sheet transient simulations for ISMIP7 / meltMIP using
[ISSM](https://issm.jpl.nasa.gov/), run on NCI Gadi (project `au88`).
CMIP6 model: **CESM2-WACCM**.  Scenarios: historical 1995–2014, ssp585 2015–2300
(other SSPs by copying the projection script and changing `SCENARIO`).

ISSM build: `$ISSM_DIR=/home/565/jb1863/trunk`.  All jobs are submitted from the
relevant subdirectory via `qsub`.

---

## Repository layout

```
init/                   initialisation, inversion, relaxation, history matching
  tuning_func.m           mesh → relaxed 1995 state
  meltMip_ensemble.m      basal-melt parameter ensemble for (meltMIP)
  init/BMB_tuning_python/ python script for determining optimal K
  hist_run_tune_CESM_WACCM.m   forcing preprocessing + SMB correction tuning
  scripts/                helpers (gridData, set_cluster, plot_*)
  Par/                    ISSM parameterisation file
  Models/                 ISSM .mat outputs (gitignored)

hist_runs/CESM2-WACCM/
  hist_run_CESM_WACCM_1995_2014.m   historical transient + ISMIP6 NetCDF output

proj_runs/CESM2-WACCM/ssp585/
  proj_run_CESM_WACCM_ssp585_2015_2300.m        main projection (with SEF)
  proj_run_CESM_WACCM_ssp585_2015_2300_noSEF.m  sensitivity: no surface-elevation feedback

preprocessed_data/      interpolated forcing .mat files (gitignored)
raw_data/               ISMIP7 NetCDF forcing files (gitignored)
postprocessed_data/     ISMIP6 output NetCDFs + figures
```

---

## 1. Initialisation (`init/`)

### 1a. Mesh, inversion, relaxation — `tuning_func.m`

Builds the relaxed 1995 AIS state used as the starting point for all runs.
Key stages (run via `qsub submit_init2.sh`):

| Stage | Purpose |
|---|---|
| `Param` | Load BedMap3 topography, build SSA mesh |
| `InversionB` | Invert for basal friction coefficient |
| `InversionFriction1` | Refine friction |
| `Inversion2B` | Second-pass inversion |
| `CollapseSSA` | Collapse to SSA (no higher-order thermal) |
_____________________________ started from here:
| `RunInitSSACollapse` | Short diagnostic run |
| `Relax_long` | Long free-surface relaxationi 50y |
| `Relaxed` | Save relaxed state (20y) → `Models/AIS_ISMIP7_Relaxed.mat` |
| `Assign_Basins` | Attach IMBIE2 basin IDs to the mesh |

Output: `init/Models/AIS_ISMIP7_Relaxed.mat` — the common upstream model.

### 1b. Forcing preprocessing and SMB tuning — `hist_run_tune_CESM_WACCM.m`

Builds all forcing `.mat` files consumed by the historical and projection runs,
and tunes the RACMO SMB climatology against Otosaka et al. (2023) GRACE/GRACE-FO
mass-change observations trends for WAIS / EAIS / Peninsula.

**Ocean TF (thermal forcing)**
Read from `raw_data/.../historical/ocean/tf/v3/` (annual NetCDF, °C).
No SSP126 TF was available; 2014 TF is repeated for 2015–2020 in this tuning run.
Thermal forcing is clamped to ≥ 0 after interpolation (required by
`basalforcingsismip6.checkconsistency`).

**Surface mass balance (SMB)**
Annual fields constructed as:
```
smb_yr = p_vert * smb_racmo + (cesm_yr - cesm_hist_mean)
```
where:
- `smb_racmo` — RACMO2.3p2 1995–2014 climatology (annual mean, mm w.e. yr⁻¹)
- `cesm_yr` — CESM2-WACCM annual SMB (`acabf`, kg m⁻² s⁻¹ × sec\_to\_year)
- `cesm_hist_mean` — CESM 1995–2014 annual mean
- `p_vert` — per-vertex region scaling factor (from step 9 below) 

CESM monthly files are averaged to **annual values** before interpolation to the
mesh (equal weight across 12 months).  The same annual convention is used for the
SMB lapse-rate gradient (`dacabfdz`, mm w.e. yr⁻¹ m⁻¹).

**SMB region correction (step 9 `HistRun_CorrectSMB`)**
After comparing modelled against observed mass change, a per-region scale factor
`p_vert` is applied to the RACMO climatology component:

| Region | Correction |
|---|---|
| WAIS | +5.5 % |
| EAIS | +3.6 % |
| Peninsula | −2.8 % |

Corrections are capped at ±25 %.  `p_vert` and the corrected `smb_forcing` matrix
are saved to `preprocessed_data/Atmosphere/Hist/` and reused by all downstream runs.

**Outputs saved to `preprocessed_data/`:**
- `Atmosphere/Clim/<model>/CESM_WACCM_SMB_clim_1995_2014.mat` (`smb_racmo`, `cesm_mean`)
- `Atmosphere/Hist/CESM_WACCM_SMB_corrected_1995_2020.mat` (`smb_forcing`, `bgrad_forcing`, `p_vert`)
- `Ocean/Hist/CESM_WACCM_TF_1995_2020.mat`
- `Ocean/Hist/Greene_spclevelset_1995_2020.mat`

---

## 2. Historical run (`hist_runs/CESM2-WACCM/`)

**Script:** `hist_run_CESM_WACCM_1995_2014.m`
**Period:** 1995–2014  **Time step:** 1/12 yr (monthly)  **Output frequency:** annual

Starts from `init/Models/AIS_ISMIP7_Relaxed_CESM_WACCM.mat` (one-year CESM-forced
relaxation of `AIS_ISMIP7_Relaxed`).

**Basal melt parameterisation**
Local quadratic ISMIP6 melt (`basalforcingsismip6`, `islocal = 1`).
Basin-specific thermal forcing from CESM2-WACCM historical ocean TF.
No thermal forcing offset correction: `delta_t = zeros(1, num_basins)`.
Melt coefficient `gamma_0` from `preprocessed_data/Ocean/gamma0_local.mat`.

**SMB**
`SMBgradients` with the corrected annual `smb_forcing` (mm w.e. yr⁻¹) plus
lapse-rate gradient `bgrad_forcing` (mm w.e. yr⁻¹ m⁻¹).  Reference surface
`href` set to the 1995 geometry at `start_year = 1995`.

**Calving front**
Annual observed ice-shelf mask from Greene (AntarcticaObsISMIP7-v1.2.nc),
converted to a level-set `spclevelset` time series.

**Step map:**

| Step | Stage | Action |
|---|---|---|
| 1 | `HistRun_1995_2014` | transient solve (submit / gather) |
| 2 | `WriteISMIP6_NetCDF` | grid to ISMIP6 8 km; write NetCDFs to `postprocessed_data/CESM2-WACCM/hist/` |
| 3 | `AIS_state_2015` | save end-state geometry/masks (input for projections) |
| 4–7 | diagnostics | grounding-line flux sanity checks and orientation diagnostics |

**CF time convention (Appendix A2.3.2)**
Base time `1995-01-01`, 365-day calendar (`calendar = "365_day"`).
ST variables (snapshots): `time = 0, 365, 730, ..., 7300` days, no `time_bnds`.
FL variables (annual fluxes): `time = 182.5, 547.5, ...` days (year midpoint),
`time_bnds` = `[0,365], [365,730], ...`.

---

## 3. Projection run (`proj_runs/CESM2-WACCM/ssp585/`)

**Script:** `proj_run_CESM_WACCM_ssp585_2015_2300.m`
**Period:** 2015–2300  **Time step:** 1/12 yr  **Output frequency:** annual

Starts from `hist_runs/CESM2-WACCM/Models/AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat`
(end-state of the historical run, step 3).

Split into two segments (2015–2150 and 2151–2300) to keep `.mat` files manageable;
end-state is saved at 2151 and used as restart.

**Ocean TF**
Annual TF fields from `raw_data/.../ssp585/ocean/tf/v3/` (decade NetCDF chunks).
Last available year (2299) repeated for 2300.  Clamped ≥ 0 after interpolation.

**Basal melt parameterisation**
Local quadratic ISMIP6 melt (`basalforcingsismip6`, `islocal = 1`).
No thermal forcing offset correction: `delta_t = zeros(1, num_basins)`.

**SMB (with surface-elevation feedback, SEF)**
`SMBgradients` using the corrected annual forcing:
```
smb_yr = p_vert * smb_racmo + (cesm_yr_ssp - cesm_hist_mean_1995_2014)
```
All fields in mm w.e. yr⁻¹ as required by `SMBgradients.smbref`.
CESM monthly `acabf` files averaged to **annual values**.
Lapse-rate gradient `bgrad_forcing` (from `dacabfdz`) applied via `b_pos / b_neg`.
Reference surface `href` set to start-of-segment geometry.

**Calving front**
Annual ISMIP7 ice-shelf-collapse mask converted to cumulative `spclevelset`
(once a shelf collapses it remains constrained).

**ISMIP6 NetCDF output (step 8 `WriteISMIP6_NetCDF`)**
Grids combined TransientSolution to ISMIP6 8 km AIS grid and writes one file per
Appendix-2 / Table A1 variable into `postprocessed_data/CESM2-WACCM/proj/` (EXP = `C007`).
Same CF time convention as historical: base time `1995-01-01`, 365-day calendar,
ST/FL distinction, time in days.

**Step map:**

| Step | Stage | Action |
|---|---|---|
| 1 | `ProjTF` | build TF .mat 2015–2300 |
| 2 | `ProjSMB` | build SMB/grad .mat 2015–2300 |
| 3 | `ProjLevelset` | build cumulative collapse spclevelset |
| 4 | `ProjRun_2015_2150` | transient 2015–2150 (submit / gather) |
| 5 | `AIS_state_2151` | save restart state |
| 6 | `ProjRun_2151_2300` | transient 2151–2300 (submit / gather) |
| 7 | `VAFContinuityCheck` | check VAF junction at 2151; save figure |
| 8 | `WriteISMIP6_NetCDF` | write ISMIP6 NetCDFs |

### 3a. noSEF sensitivity — `proj_run_CESM_WACCM_ssp585_2015_2300_noSEF.m`

Identical to the main projection except SMB is applied as a plain time-varying
field with no lapse-rate correction (`SMBforcing.mass_balance`).  The same
corrected `smb_forcing` matrix is used, but converted to m ice yr⁻¹ before
assignment (`/ 1000 * rho_w / rho_i`) because `SMBforcing` does not apply the
internal `/1000 * rho_w/rho_i` conversion that `SMBgradients` does.
Output goes to `postprocessed_data/CESM2-WACCM/proj_noSEF/` (EXP = `C007_noSEF`).

---

## Submitting jobs

```bash
# Historical
cd hist_runs/CESM2-WACCM
qsub submit_hist_1995_2014.sh          # steps=[1], loadonly=0  -> submit
qsub submit_hist_1995_2014.sh          # steps=[1], loadonly=1  -> gather

# Projection (main)
cd proj_runs/CESM2-WACCM/ssp585
qsub submit_proj_ssp585.sh             # edit steps/loadonly in the script

# Projection (noSEF)
qsub submit_proj_ssp585_noSEF.sh
```

Edit `steps=[...]` and `loadonly=[0|1]` directly in each submit script before
queuing.  Forcing-preprocessing steps (TF / SMB / Levelset) are single-node,
fast (~minutes); transient solve steps use 48 cores, 190 GB RAM, ~40 h walltime.

---

## Key references

- ISMIP7 protocol and Appendix A2 output conventions: see `CLAUDE.md` (this repo)
  and the ISMIP7 forcing repository.
- Basal melt parameterisation: ISMIP6 local quadratic scheme
  (Jourdain et al. 2020, *The Cryosphere*).
- SMB: RACMO2.3p2 climatology (van Wessem et al.) + CESM2-WACCM anomaly (CMIP6).
- Observed mass change for SMB tuning: Otosaka et al. 2023, *Earth Syst. Sci. Data*.
- Ice-shelf calving observations: Greene et al., AntarcticaObsISMIP7-v1.2.
