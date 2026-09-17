#!/bin/bash
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=12:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitAnalyzeYearly_ssp370.outlog
#PBS -e SubmitAnalyzeYearly_ssp370.errlog
#PBS -l storage=gdata/au88

# Submits analyze_proj_ssp370_yearly.m -- pure post-processing (VAF/SLE
# check + ISMIP6 NetCDF write + calving-front plot) over whatever
# Models_yearly/*.mat files proj_run_CESM_WACCM_ssp370_2015_2100_yearly.m
# has produced so far. No PBS solve is submitted from within this MATLAB
# session at all (unlike submit_proj_ssp370_yearly.sh) -- this just
# loads, concatenates, and writes, so waitonlock/the local waitonlock.m
# override are irrelevant here and are not on the path.
#
# Mirrors submit_analyze_proj_ssp585_yearly.sh / submit_analyze_proj_ssp126_yearly.sh
# exactly (same PBS resources, same workflow). This scenario's own run
# only covers 2015-2100 (86 years, not 286), so this should finish
# noticeably faster than ssp585/ssp126's equivalents, especially for
# steps 1 and 3.

export ISSM_DIR=/home/565/jb1863/trunk
source $ISSM_DIR/etc/environment.sh
module purge
module load openmpi/4.1.3
module load netcdf/4.8.0p
module load hdf5/1.10.7p
module load petsc/3.17.4
module load matlab/R2021b
module load matlab_licence/monash

# ---- steps: [1]=VAF check only, [2]=NetCDF write only, [3]=calving
# front plot only, [1 2 3]=all (default) -----
steps=[3]

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
  addpath('$ISSM_DIR/lib'); \
  analyze_proj_ssp370_yearly($steps), quit" \
  > FuncAnalyzeYearly_ssp370.log
