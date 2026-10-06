#!/bin/bash
#PBS -P au88
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=22:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitPPEfrozenfront.outlog
#PBS -e SubmitPPEfrozenfront.errlog
#PBS -l storage=gdata/au88

# Submits proj_run_PPE_frozenfront_ssp585.m -- covers PPE members p011
# (K=mode), p012 (K=95th percentile), p013 (K=5th percentile), each with
# the calving front FROZEN at 2015 (no collapse mask, no floating-gating
# loop). NOT the core ssp585 script or its own submit script -- both
# left untouched.
#
# No separate P-number is reserved for the core/baseline ssp585 run --
# it has no P-number at all. P-numbers continue on from the ssp370 PPE
# family (P001-P009) rather than restarting at P001, so numbering here
# starts at P011.
#
# RESTART-SPLIT at 2150/2151 (NOT a single continuous 2015-2300 transient
# like ../ctrl/submit_proj_ctrl.sh): a single-job attempt at the full
# 285-year span got OOM-killed on this same 48cpu/190GB node (PBS:
# "exceeded memory allocation") -- holding the whole run's
# TransientSolution in RAM at once was too much. Split into two ~140-year
# segments instead, same technique as this project's other 300-year-class
# runs (see ../ssp585/old_scripts/proj_run_CESM_WACCM_ssp585_2015_2300.m).
# See proj_run_PPE_frozenfront_ssp585.m's own module docstring for the
# full 4-step map (ProjRun_2015_2150 / AIS_state_2151 / ProjRun_2151_2300
# / WriteNetCDF).
#
# proj_run_PPE_frozenfront_ssp585.m signature is (X, steps, loadonly) --
# X selects the PPE member (internally resolves K variant + output
# P-number); steps picks which organizer step to run; loadonly is the
# standard two-phase submit(0)/gather(1) flag, used only within the solve
# steps (1 and 3). loadonly=1 ONLY gathers + saves the model -- it does
# NOT write the NetCDF output.
#
# Each PPE member's model is saved to its OWN ./Models_P<NNN>/ folder
# (e.g. ./Models_P011/ for p011) -- distinct per member.
#
# X/steps/loadonly can also be passed in at qsub time via -v, e.g.:
#   qsub -v X=p012,steps='[1]',loadonly=0 -N PPE_frozenfront_p012_submit1 submit_proj_PPE_frozenfront_ssp585.sh
#   qsub -v X=p012,steps='[1]',loadonly=1 -N PPE_frozenfront_p012_gather1 submit_proj_PPE_frozenfront_ssp585.sh
#   qsub -v X=p012,steps='[2]'            -N PPE_frozenfront_p012_restart submit_proj_PPE_frozenfront_ssp585.sh
#   qsub -v X=p012,steps='[3]',loadonly=0 -N PPE_frozenfront_p012_submit2 submit_proj_PPE_frozenfront_ssp585.sh
#   qsub -v X=p012,steps='[3]',loadonly=1 -N PPE_frozenfront_p012_gather2 submit_proj_PPE_frozenfront_ssp585.sh
#   qsub -v X=p012,steps='[4]'            -N PPE_frozenfront_p012_netcdf  submit_proj_PPE_frozenfront_ssp585.sh
# (used by launch_all_PPE_frozenfront.sh's six MODE values to run this
# whole sequence for all three members at once) -- if not supplied, the
# X=/steps=/loadonly= defaults below are used, so editing those lines and
# resubmitting directly still works exactly as before for a single
# one-off run.

export ISSM_DIR=/home/565/jb1863/trunk
source $ISSM_DIR/etc/environment.sh
module purge
module load openmpi/4.1.3
module load netcdf/4.8.0p
module load hdf5/1.10.7p
module load petsc/3.17.4
module load matlab/R2021b
module load matlab_licence/monash

# ---- Edit these lines to pick which PPE member / step / phase to run --
# (ignored if X/steps/loadonly were already set via qsub -v)
# X: 'p011' (K=mode) | 'p012' (K=95th) | 'p013' (K=5th)
# steps: '[1]'=ProjRun_2015_2150 | '[2]'=AIS_state_2151 | '[3]'=ProjRun_2151_2300 | '[4]'=WriteNetCDF
# loadonly: 0=submit, 1=gather (save model only, no NetCDF) -- only used by steps 1/3
X=${X:-p011}
steps=${steps:-'[1]'}
loadonly=${loadonly:-0}

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
addpath('$ISSM_DIR/lib'); \
proj_run_PPE_frozenfront_ssp585('$X', $steps, $loadonly), quit" \
> FuncPPEfrozenfront_$X.log
