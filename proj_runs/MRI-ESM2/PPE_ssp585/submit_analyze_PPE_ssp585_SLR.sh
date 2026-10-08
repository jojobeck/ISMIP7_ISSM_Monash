#!/bin/bash
#PBS -P au88
#PBS -q normal
#PBS -l ncpus=4
#PBS -l walltime=01:00:00
#PBS -l mem=16GB
#PBS -l jobfs=20GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitAnalyzePPESLR.outlog
#PBS -e SubmitAnalyzePPESLR.errlog
#PBS -l storage=gdata/au88

# Submits analyze_PPE_ssp585_SLR.m -- pure post-processing (VAF/SLE
# extraction + plot) over the MRI-ESM2-0 historical run, the core ssp585
# run, and every PPE member listed in ppe_members.csv (this folder). No
# PBS solve is submitted from within this MATLAB session -- just loads
# each run's already-finished/in-progress TransientSolution output and
# writes/plots, so this is fast and lightweight (ncpus/mem well below the
# PPE solve jobs themselves).
#
# Gracefully skips any member/run not available yet (or an unfinished
# yearly member is included up to however far it's got) -- see the
# script's own module docstring.

export ISSM_DIR=/home/565/jb1863/trunk
source $ISSM_DIR/etc/environment.sh
module purge
module load openmpi/4.1.3
module load netcdf/4.8.0p
module load hdf5/1.10.7p
module load petsc/3.17.4
module load matlab/R2021b
module load matlab_licence/monash

# ---- steps: [1]=compute+save SLR .mat, [2]=plot only, [1,2]=both -------
steps=[1,2]

matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
  addpath('$ISSM_DIR/lib'); \
  analyze_PPE_ssp585_SLR($steps), quit" \
  > FuncAnalyzePPESLR.log
