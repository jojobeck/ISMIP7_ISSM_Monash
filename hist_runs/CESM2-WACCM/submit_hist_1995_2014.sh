
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=20:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitHist1995_2014.outlog
#PBS -e SubmitHist1995_2014.errlog
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

steps=[1]
loadonly=[1]
# hist_run_CESM_WACCM_1995_2014: submit step 1 (HistRun_1995_2014), loadonly=0
# Once the PBS solve finishes, rerun with steps=[1] loadonly=[1] to gather,
# then steps=[2] loadonly=[1] (WriteISMIP6_NetCDF) to write the
# postprocessed_data/CESM2-WACCM/hist/ NetCDF files.
# steps=[4] (GLFluxSanityCheck) computes the grounding-line sanity-check
# numbers and saves them to hist/tables/ -- no plotting, so it's safe from
# the figure()/saveas() segfault below. Run steps=[5] (GLFluxSanityCheckPlot)
# separately to attempt the plot.
# steps=[6] (GLFluxOrientationDiagnostic) computes the with/without-
# orientation-correction diagnostic and saves to hist/tables/; run
# steps=[7] (GLFluxOrientationDiagnosticPlot) separately to plot it.
# -softwareopengl forces software rendering for figure()/saveas() calls --
# tried as a fix for MATLAB's graphics engine segfaulting on this compute
# node (no display/GPU), but it did NOT resolve the crash (still segfaults
# even with -nodisplay, figure('visible','off'), and this flag); kept in
# case it still helps reduce crash frequency, but isn't a confirmed fix.
matlab -nodisplay -nosplash -softwareopengl -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; hist_run_CESM_WACCM_1995_2014($steps, $loadonly), quit" >FuncHist1995_2014.log
