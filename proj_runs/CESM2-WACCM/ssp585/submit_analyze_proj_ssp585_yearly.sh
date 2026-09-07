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
#PBS -o SubmitAnalyzeYearly_ssp585.outlog
#PBS -e SubmitAnalyzeYearly_ssp585.errlog
#PBS -l storage=gdata/au88

# Submits analyze_proj_ssp585_yearly.m -- pure post-processing (VAF/SLE
# check + ISMIP6 NetCDF write) over whatever Models_yearly/*.mat files
# proj_run_CESM_WACCM_ssp585_2015_2300_yearly.m has produced so far. No
# PBS solve is submitted from within this MATLAB session at all (unlike
# submit_proj_ssp585_yearly.sh) -- this just loads, concatenates, and
# writes, so waitonlock/the local waitonlock.m override are irrelevant
# here and are not on the path.
#
# ncpus/mem match submit_proj_ssp585.sh's own NetCDF-writing invocation
# (same write_ismip7_2d_projection / write_ismip7_scalar_projection
# helpers, gridding the same full 2015-2300 span onto the 8km ISMIP6
# grid) -- if anything this script's memory footprint could be a bit
# higher, since it reads up to ~285 separate per-year .mat files
# (repeated mesh/geometry data each) rather than 2 large segment files.
# Not yet profiled/tightened specifically for this script.
#
# walltime=12:00:00 is a conservative starting guess (no PBS-job waiting
# is involved, only file I/O + gridding + NetCDF writes) -- adjust down
# once you've seen how long a real run actually takes.

export ISSM_DIR=/home/565/jb1863/trunk
source $ISSM_DIR/etc/environment.sh
module purge
module load openmpi/4.1.3
module load netcdf/4.8.0p
module load hdf5/1.10.7p
module load petsc/3.17.4
module load matlab/R2021b
module load matlab_licence/monash

# ---- steps: [1]=VAF check only, [2]=NetCDF write only, [1 2]=both -----
steps=[2]

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
  addpath('$ISSM_DIR/lib'); \
  analyze_proj_ssp585_yearly($steps), quit" \
  > FuncAnalyzeYearly_ssp585.log
