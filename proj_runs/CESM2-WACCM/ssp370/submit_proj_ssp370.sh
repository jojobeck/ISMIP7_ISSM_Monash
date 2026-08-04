#!/bin/bash
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=24:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitProj_ssp370.outlog
#PBS -e SubmitProj_ssp370.errlog
#PBS -l storage=gdata/au88

export ISSM_DIR=/home/565/jb1863/trunk
source $ISSM_DIR/etc/environment.sh
module purge
module load openmpi/4.1.3
module load netcdf/4.8.0p
module load hdf5/1.10.7p
module load petsc/3.17.4
module load matlab/R2021b
module load matlab_licence/monash

# Workflow (ssp370 only runs 2015-2100 -- single continuous transient, no
# restart split, so this is a shorter step map than ssp126/ssp585):
#   Steps 1-3 (forcing preprocessing): run once, no PBS solve involved.
#     steps=[1] loadonly=[1]    -> build TF forcing
#     steps=[2] loadonly=[1]    -> build SMB forcing
#     steps=[3] loadonly=[1]    -> build collapse levelset
#
#   Step 4 (2015-2100): two-stage PBS solve
#     steps=[4] loadonly=[0]    -> submit PBS job (writes .bin/.queue/.toolkits)
#     steps=[4] loadonly=[1]    -> gather results after PBS job completes
#
#   Step 5: write ISMIP6 NetCDFs test (instant-ish, no PBS solve)
#     steps=[5] loadonly=[1]
#
#   Step 6: write ISMIP6 NetCDFs (slow, no PBS solve; may need longer walltime)
#     steps=[6] loadonly=[1]
#
# -softwareopengl: kept from hist submit for consistency; does not fully
# prevent MATLAB graphics crashes on this node stack, but may reduce frequency.

steps=[7]
loadonly=[1]

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; proj_run_CESM_WACCM_ssp370_2015_2100($steps, $loadonly), quit" >FuncProj_ssp370.log
