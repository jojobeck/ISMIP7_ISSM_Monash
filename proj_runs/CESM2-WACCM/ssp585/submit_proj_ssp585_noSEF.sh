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
#PBS -o SubmitProj_ssp585_noSEF.outlog
#PBS -e SubmitProj_ssp585_noSEF.errlog
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

# Prerequisite: main script steps 1-3 must have been run first to generate
# forcing mats in preprocessed_data/Ocean/Proj/ and Atmosphere/Proj/.
#
# Workflow:
#   Step 1 phase 1 (2015-2150): two-stage PBS solve
#     steps=[1] loadonly=[0]    -> submit PBS job
#     steps=[1] loadonly=[1]    -> gather results after PBS job completes
#
#   Step 2: save AIS_state_2151_noSEF (instant, depends on step 1 gathered)
#     steps=[2] loadonly=[1]
#
#   Step 3 phase 2 (2151-2300): two-stage PBS solve
#     steps=[3] loadonly=[0]    -> submit PBS job
#     steps=[3] loadonly=[1]    -> gather results
#
#   Step 4: VAF continuity check + SEF vs noSEF comparison figure (instant)
#     steps=[4] loadonly=[1]
#
#   Step 5: write ISMIP6 NetCDFs (slow, no PBS solve; may need longer walltime)
#     steps=[5] loadonly=[1]

steps=[3,4]
loadonly=[1]

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; proj_run_CESM_WACCM_ssp585_2015_2300_noSEF($steps, $loadonly), quit" >FuncProj_ssp585_noSEF.log
