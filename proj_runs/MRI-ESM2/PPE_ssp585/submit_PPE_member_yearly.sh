#!/bin/bash
#PBS -P au88
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=22:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -o SubmitPPEyearly.outlog
#PBS -e SubmitPPEyearly.errlog
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -l storage=gdata/au88

# Per-member submission TEMPLATE for the loop (yearly-restart) PPE
# members (P024-P029) -- X is NOT hardcoded here. It must be passed in at
# qsub time via -v, e.g.:
#   qsub -v X=p024 -N PPE_p024 submit_PPE_member_yearly.sh
# One script file now covers every member, submitted once per member as
# its OWN independent PBS job. This replaces the previous pattern (one
# submit_proj_PPE_yearly_ssp585.sh with X hardcoded, manually edited and
# resubmitted once per member) -- see launch_all_PPE_yearly.sh, which
# submits each yearly member via its own explicit, individually
# commentable qsub call to THIS script (not CSV-driven -- ppe_members.csv
# is just the tracking table of what each P-number is). The dispatcher
# itself finishes in seconds (it only calls qsub repeatedly, once per
# uncommented line); it does NOT wait for any submitted job to finish,
# unlike each individual member's own yearly-restart job below, which
# does need to run for its own full duration (it loops over years
# internally, submit-wait-gather, via the shared
# ./../local_matlab/waitonlock.m override).
#
# Walltime kept at 22h, same as the ssp370 version, even though this run
# spans 2015-2300 (285 years) instead of ssp370's 2015-2100 (85 years):
# one completed ssp370 yearly member's own PBS log showed ~2.5h real
# elapsed for the full 85-year loop, so scaling that to 285 years should
# still land well inside 22h -- revisit this if a real ssp585 member
# run ever comes close to it.
#
# start_year/end_year/waitonlock_seconds are optional overrides, also via
# qsub -v (e.g. -v X=p026,start_year=2080 to resume p026 from a crash at
# 2080) -- default to 2015/2300/3600 if not supplied (end_year default
# matches ssp585's own 2015-2300 span, NOT ssp370's 2015-2100), matching
# proj_run_PPE_yearly_ssp585.m's own defaults.

export ISSM_DIR=/home/565/jb1863/trunk
source $ISSM_DIR/etc/environment.sh
module purge
module load openmpi/4.1.3
module load netcdf/4.8.0p
module load hdf5/1.10.7p
module load petsc/3.17.4
module load matlab/R2021b
module load matlab_licence/monash

if [ -z "$X" ]; then
    echo "ERROR: X not set. Submit with: qsub -v X=p024 submit_PPE_member_yearly.sh" >&2
    exit 1
fi

start_year=${start_year:-2015}
end_year=${end_year:-2300}
waitonlock_seconds=${waitonlock_seconds:-3600}

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
addpath('$ISSM_DIR/lib'); \
proj_run_PPE_yearly_ssp585('$X', $start_year, $end_year, $waitonlock_seconds), quit" \
> FuncPPEyearly_$X.log
