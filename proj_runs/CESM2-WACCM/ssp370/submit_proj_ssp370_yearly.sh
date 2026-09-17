#!/bin/bash
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=22:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitProjYearly_ssp370.outlog
#PBS -e SubmitProjYearly_ssp370.errlog
#PBS -l storage=gdata/au88

# Submits proj_run_CESM_WACCM_ssp370_2015_2100_yearly.m -- NOT
# old_scripts/proj_run_CESM_WACCM_ssp370_2015_2100.m /
# old_scripts/submit_proj_ssp370.sh, both left untouched. This ONE PBS job
# runs the ENTIRE requested start_year:end_year loop start to finish, one
# year at a time, inside a single continuous MATLAB session: each year is
# submitted, waited-on (blocking, via the shared ./../../local_matlab/waitonlock.m
# override -- see that file and the main script's header comment for why
# the stock ISSM waitonlock.m fails when run from inside a Gadi compute
# node), and gathered before the next year is built -- no manual
# re-submission or separate gather step needed.
#
# Mirrors submit_proj_ssp585_yearly.sh / submit_proj_ssp126_yearly.sh
# exactly (same PBS resources, same workflow). NOTE: end_year default
# below is 2100, not 2300 -- ssp370 only has real forcing data through
# 2100 (confirmed on disk).

export ISSM_DIR=/home/565/jb1863/trunk
source $ISSM_DIR/etc/environment.sh
module purge
module load openmpi/4.1.3
module load netcdf/4.8.0p
module load hdf5/1.10.7p
module load petsc/3.17.4
module load matlab/R2021b
module load matlab_licence/monash

# ---- Test configuration: edit these three lines between test stages ----
# Stage 1 (1-year test):  start_year=2015  end_year=2015  waitonlock_seconds=3600
start_year=2060
end_year=2100
waitonlock_seconds=3600
step=[4]

# proj_run_CESM_WACCM_ssp370_2015_2100_yearly.m signature: steps is the
# FIRST argument (steps=[4] runs just the annual loop below, assuming
# steps 1-3 -- ProjTF/ProjSMB/ProjLevelsetV21 -- were already run
# separately; steps=[1 2 3] builds all three, steps=[1 2 3 4] does
# everything in one go). See that script's own header for full step map
# and why step 3 uses the v2.1 collapse mask (not v2, matching ssp126).
matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
addpath('$ISSM_DIR/lib'); \
proj_run_CESM_WACCM_ssp370_2015_2100_yearly($step, $start_year, $end_year, $waitonlock_seconds), quit" \
> FuncProjYearly_ssp370.log
