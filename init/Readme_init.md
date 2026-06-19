# ISSM ISMIP7 — init/ pipeline

ISSM ISMIP7 setup, initialized around 1995. This describes the order of
operations carried out in this `init/` folder to go from an existing,
already-initialized old project through to the historical CESM2-WACCM
validation run.

## 1. Start from an old, already-initialized project

Rather than re-running the full inversion/thermal spin-up from scratch, this
project's `tuning_func.m` (`'Param'` step) loads the mesh and inversion
result of an older project as its starting point:

```
inputmodel_path = './../AIS_1850/Models/AIS1850_thTHW_CollapseSSA.mat'
```

i.e. velocity inversion, friction, and the thermal field are already
initialized/converged from that earlier run — they are not redone here.

## 2. Geometry: thicker Thwaites from Bedmap on top of BedMachine

The base geometry comes from BedMachine Antarctica v3
(`BedMachineAntarcticav3_plusSECschroeder1992_plusSmithfloat.nc`), via
`Par/Antartica_1995_Thwbedmap.par`. Over the Thwaites region specifically,
BedMachine's thickness/surface is overridden by a separate Bedmap2 product
(`bedmap2_1km_justTHW_thk_mask_usurf.nc`), which gives a thicker Thwaites
than BedMachine alone. Base/thickness/mask are then re-derived from the
combined surface/thickness/bed fields.Also thickness form smith is added.

## 3. Relaxation: prescribed basal melt + RACMO SMB, ~20 years

Before any melt-rate tuning, the model is relaxed under steady forcing:

- Basal melt: prescribed from observed ISMIP7 melt-per-basin
  (`Obs_ISMIP7melt`), constant in time (`'constant_0.5cmeanBMB_SMB_from_Collapse'`
  step).
- SMB: static RACMO 1995-2014 climatology
  (`racmo_rec_smb_2km_1995_2014_mean.mat`).
- Run for 50 years (`'Relax_long'` step, `Relax50`), but only the first
  ~20 years are used as the relaxed state (`i_t = 20` in the `'Relaxed'`
  step) — chosen because drift in thickness/VAF/grounded area has largely
  plateaued by then.
- The `'Relaxed'` step folds that snapshot's `Base`/`Thickness`/`Surface`/
  `MaskOceanLevelset` back into `md.geometry`/`md.mask` and saves it as
  `./Models/AIS_ISMIP7_Relaxed.mat`.

This `Relaxed` model is the common starting point (`inputmodel_relax`) for
everything downstream: `meltMip_ensemble.m` and
`hist_run_tune_CESM_WACCM.m` both load it directly.

## 4. MeltMIP K tuning

From the relaxed state, `meltMip_ensemble.m` builds the basal-melt ensemble
across the heat-exchange coefficient `K` parameter space (quadratic local
melt parameterisation), gridding ISSM output to NetCDF
(`Models/ModelNC/BMB*/`). The Python toolbox in `BMB_tuning_python/`
(`parameter_selection_quadratic_example.ipynb` /
`parameter_selection_toolbox.py`) then selects the optimal `K` by minimising
the four-term objective function against observational and ocean-model
melt targets (see the top-level `CLAUDE.md` for details on that toolbox).

## 5. Historical matching: CESM2-WACCM

`hist_run_tune_CESM_WACCM.m` takes the relaxed model (now also passed
through a short 1995-1996 `Relaxed_CESM_WACCM` step to absorb the initial
jump from switching on real time-varying forcing) and runs the 1995-2020
historical transient under CESM2-WACCM ocean thermal forcing, SMB, and a
Greene-levelset-based calving front, validating modelled regional
(WAIS/EAIS/Peninsula) mass change against Otosaka et al. See the step map
at the top of `hist_run_tune_CESM_WACCM.m` for the full step-by-step
breakdown of this part of the pipeline.
