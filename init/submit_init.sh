
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=20:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitInit.outlog
#PBS -e SubmitInit.errlog
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

steps=[4]
loadonly=[1]
# matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; run('tuning_func($steps, $loadonly)') , quit" >FuncInit.log

# matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; run('meltMip_ensemble($steps,1, $loadonly)') , quit" >FuncInit.log
# for j in $(seq 1 120); do
# for j in $(seq 2 120); do
for j in $(seq 1 1); do
    matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; run('meltMip_ensemble($steps, $j, $loadonly)'), quit" > FuncInit_${j}.log
done
# runs=[31]

