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
# members (P004-P009) -- X is NOT hardcoded here. It must be passed in at
# qsub time via -v, e.g.:
#   qsub -v X=p004 -N PPE_p004 submit_PPE_member_yearly.sh
# One script file now covers every member, submitted once per member as
# its OWN independent PBS job. This replaces the previous pattern (one
# submit_proj_PPE_yearly_ssp370.sh with X hardcoded, manually edited and
# resubmitted once per member) -- see launch_all_PPE_yearly.sh, which
# submits each yearly member via its own explicit, individually
# commentable qsub call to THIS script (not CSV-driven -- ppe_members.csv
# is just the tracking table of what each P-number is). The dispatcher
# itself finishes in seconds (it only calls qsub repeatedly, once per
# uncommented line); it does NOT wait for any submitted job to finish,
# unlike each individual member's own ~22-hour yearly-restart job below,
# which does need to run for its own full duration (it loops over years
# internally, submit-wait-gather, via the shared
# ./../local_matlab/waitonlock.m override).
#
# start_year/end_year/waitonlock_seconds are optional overrides, also via
# qsub -v (e.g. -v X=p006,start_year=2080 to resume p006 from a crash at
# 2080) -- default to 2015/2100/3600 if not supplied, matching
# proj_run_PPE_yearly_ssp370.m's own defaults.

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
    echo "ERROR: X not set. Submit with: qsub -v X=p004 submit_PPE_member_yearly.sh" >&2
    exit 1
fi

start_year=${start_year:-2015}
end_year=${end_year:-2100}
waitonlock_seconds=${waitonlock_seconds:-3600}

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
addpath('$ISSM_DIR/lib'); \
proj_run_PPE_yearly_ssp370('$X', $start_year, $end_year, $waitonlock_seconds), quit" \
> FuncPPEyearly_$X.log
