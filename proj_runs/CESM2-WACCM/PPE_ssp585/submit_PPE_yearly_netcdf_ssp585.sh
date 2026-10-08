#!/bin/bash
#PBS -P au88
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=12:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitPPEyearlyNetCDF.outlog
#PBS -e SubmitPPEyearlyNetCDF.errlog
#PBS -l storage=gdata/au88

# Per-member submission TEMPLATE for writing the ISMIP7 NetCDF of a loop
# (yearly-restart) PPE member (P014-P019) via
# write_PPE_yearly_netcdf_ssp585.m -- pure post-processing, no PBS solve.
# X is passed in at qsub time, e.g.:
#   qsub -v X=p014 -N PPE_p014_netcdf submit_PPE_yearly_netcdf_ssp585.sh
# (launch_all_PPE_yearly_netcdf.sh does this per uncommented member).
#
# ncpus/mem/walltime match ../ssp585/submit_analyze_proj_ssp585_yearly.sh,
# which does the same load-concatenate-grid of ~285 yearly files.

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
    echo "ERROR: X not set. Submit with: qsub -v X=p014 submit_PPE_yearly_netcdf_ssp585.sh" >&2
    exit 1
fi

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
addpath('$ISSM_DIR/lib'); \
write_PPE_yearly_netcdf_ssp585('$X'), quit" \
> FuncPPEyearlyNetCDF_$X.log
