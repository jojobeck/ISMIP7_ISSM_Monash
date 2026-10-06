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

# Submits proj_run_PPE_frozenfront_ssp370.m -- covers PPE members p001
# (K=mode), p002 (K=95th percentile), p003 (K=5th percentile), each with
# the calving front FROZEN at 2015 (no collapse mask, no floating-gating
# loop) -- a single continuous 2015-2100 transient, same shape as
# ../ctrl/submit_proj_ctrl.sh. NOT the core ssp370 script or its own
# submit script -- both left untouched.
#
# No separate P-number is reserved for the core/baseline ssp370 run --
# it has no P-number at all -- so PPE numbering starts at P001 here.
#
# proj_run_PPE_frozenfront_ssp370.m signature is (X, steps, loadonly) --
# X selects the PPE member (internally resolves K variant + output
# P-number); steps picks which organizer step to run (1=ProjRun solve,
# 2=WriteNetCDF); loadonly is the standard two-phase submit(0)/gather(1)
# flag, used only within step 1. loadonly=1 ONLY gathers + saves the
# model -- it does NOT write the NetCDF output. Run steps=2 separately
# (after step 1's gather has finished and saved) to write the NetCDF.
#
# Each PPE member's model is saved to its OWN ./Models_P<NNN>/ folder
# (e.g. ./Models_P001/ for p001) -- distinct per member.
#
# X/steps/loadonly can also be passed in at qsub time via -v, e.g.:
#   qsub -v X=p002,steps='[1]',loadonly=0 -N PPE_frozenfront_p002 submit_proj_PPE_frozenfront_ssp370.sh
#   qsub -v X=p002,steps='[2]'            -N PPE_frozenfront_p002_netcdf submit_proj_PPE_frozenfront_ssp370.sh
# (used by launch_all_PPE_frozenfront.sh to submit/gather/write-NetCDF for
# all three members at once) -- if not supplied, the X=/steps=/loadonly=
# defaults below are used, so editing those lines and resubmitting
# directly still works exactly as before for a single one-off run.

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
# X: 'p001' (K=mode) | 'p002' (K=95th) | 'p003' (K=5th)
# steps: '[1]'=ProjRun solve | '[2]'=WriteNetCDF (run after step 1's gather)
# loadonly: 0=submit, 1=gather (save model only, no NetCDF) -- only used by step 1
X=${X:-p001}
steps=${steps:-'[1]'}
loadonly=${loadonly:-0}

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
addpath('$ISSM_DIR/lib'); \
proj_run_PPE_frozenfront_ssp370('$X', $steps, $loadonly), quit" \
> FuncPPEfrozenfront_$X.log
