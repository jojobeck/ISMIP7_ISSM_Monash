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
#PBS -o SubmitProjYearly_ssp126.outlog
#PBS -e SubmitProjYearly_ssp126.errlog
#PBS -l storage=gdata/au88

# Submits proj_run_CESM_WACCM_ssp126_2015_2300_yearly.m -- NOT the original
# submit_proj_ssp126.sh / proj_run_CESM_WACCM_ssp126_2015_2300.m, both left
# untouched. This ONE PBS job runs the ENTIRE requested start_year:end_year
# loop start to finish, one year at a time, inside a single continuous
# MATLAB session: each year is submitted, waited-on (blocking, via the
# local ./local_matlab/waitonlock.m override -- see that file and the main
# script's header comment for why the stock ISSM waitonlock.m fails when
# run from inside a Gadi compute node), and gathered before the next year
# is built -- no manual re-submission or separate gather step needed.
#
# Mirrors submit_proj_ssp585_yearly.sh exactly (same PBS resources, same
# workflow). ncpus/mem here match the original submit_proj_ssp126.sh's
# outer-job request.

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
start_year=2015
end_year=2300
waitonlock_seconds=3600
step=[4]

# proj_run_CESM_WACCM_ssp126_2015_2300_yearly.m signature: steps is the
# FIRST argument (steps=[4] runs just the annual loop below, assuming
# steps 1-3 -- ProjTF/ProjSMB/ProjLevelsetV21 -- were already run
# separately; steps=[1 2 3] builds all three, steps=[1 2 3 4] does
# everything in one go). See that script's own header for full step map
# and why step 3 uses the v2.1 collapse mask (not v2, unlike ssp585).
matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
addpath('$ISSM_DIR/lib'); \
proj_run_CESM_WACCM_ssp126_2015_2300_yearly($step, $start_year, $end_year, $waitonlock_seconds), quit" \
> FuncProjYearly_ssp126.log
