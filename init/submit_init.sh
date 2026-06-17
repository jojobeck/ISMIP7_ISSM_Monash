
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=20:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitInit2.outlog
#PBS -e SubmitInit2.errlog
#PBS -l storage=gdata/au88

export ISSM_DIR=/home/565/jb1863/trunk
source $ISSM_DIR/etc/environment.sh
module purge
module load openmpi/4.1.3
module load netcdf/4.8.0p
module load hdf5/1.10.7p
module load petsc/3.17.4
module load matlab/R2021b
# module load python3/3.11.0
module load matlab_licence/monash
# source $ISSM_DIR/scripts/startup.sh

steps=[4,5]
loadonly=[1]
# tuning_func (init / inversion / relaxation pipeline):
# matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; tuning_func($steps, $loadonly), quit" >FuncInit2.log
# hist_run_tune_CESM_WACCM (historical validation run):
matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; hist_run_tune_CESM_WACCM($steps, $loadonly), quit" >FuncInit2.log
# runs=[31]

