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
#PBS -o SubmitProj_ssp585.outlog
#PBS -e SubmitProj_ssp585.errlog
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

# Workflow:
#   Steps 1-3 (forcing preprocessing): run once, no PBS solve involved.
#     steps=[1] loadonly=[1]    -> build TF forcing
#     steps=[2] loadonly=[1]    -> build SMB forcing
#     steps=[3] loadonly=[1]    -> build collapse levelset
#
#   Step 4 phase 1 (2015-2150): two-stage PBS solve
#     steps=[4] loadonly=[0]    -> submit PBS job (writes .bin/.queue/.toolkits)
#     steps=[4] loadonly=[1]    -> gather results after PBS job completes
#
#   Step 5: save AIS_state_2151 (instant, depends on step 4 gathered)
#     steps=[5] loadonly=[1]
#
#   Step 6 phase 2 (2151-2300): two-stage PBS solve
#     steps=[6] loadonly=[0]    -> submit PBS job
#     steps=[6] loadonly=[1]    -> gather results
#
#   Step 7: VAF continuity check (instant, no PBS solve)
#     steps=[7] loadonly=[1]
#
#   Step 8: write ISMIP6 NetCDFs (slow, no PBS solve; may need longer walltime)
#     steps=[8] loadonly=[1]
#
# -softwareopengl: kept from hist submit for consistency; does not fully
# prevent MATLAB graphics crashes on this node stack, but may reduce frequency.

steps=[8]
loadonly=[1]

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; proj_run_CESM_WACCM_ssp585_2015_2300($steps, $loadonly), quit" >FuncProj_ssp585.log
