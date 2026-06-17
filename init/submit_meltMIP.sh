
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=20:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o SubmitMeltMIP.outlog
#PBS -e SubmitMeltMIP.errlog
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

# Step map (meltMip_ensemble) -- RENUMBERED: TF builders first, then runs, then gridding.
#   TF builders : 1 Obs_clim_TF   2 OceanModelling_clim_TF   3 ObsData_clim_TF
#   melt runs   : 4 melt_run       5 melt_run_OceanModelling (Nw,Nc,Mw,Mc)   6 melt_run_ObsData
#   gridding    : 7 create_BMB_gD  8 create_BMB_gD_4km
#                 9 create_BMB_gD_OceanModelling   10 create_BMB_gD_ObsData
# NOTE: the standalone melt_run_OCeanModelling_Nw step was folded into step 5.
# TF builders fixed: tf cell now written per depth slice (tf(:,:,i)), not collapsed.
#
# Full (re)build order after the TF depth fix:
#   1) TF mats   : steps=[1] ; steps=[2] ; steps=[3]   (loadonly=[0], run once, j ignored)
#   2) launch    : steps=[4] / steps=[5] / steps=[6]   loadonly=[0]   (overwrites old runs)
#   3) gather    : steps=[4] / steps=[5] / steps=[6]   loadonly=[1]
#   4) grid      : steps=[7] / steps=[8] / steps=[9] / steps=[10]   loadonly=[0]
#   5) ensembles :
#        python3 creating_bmbMIP_ensmble.py                      # pd_ensemble (8km)
#        python3 creating_bmbMIP_OceanModelling_ensemble.py      # cold/warm ensembles
#        python3 creating_bmbMIP_ObsData_ensemble.py             # obs_ensemble
steps=[1,2,3]
loadonly=[0]
steps=[7,9,10]
steps=[11]
# loadonly=[1]
# j = [1];

matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; run('meltMip_ensemble($steps,1, $loadonly)') , quit" >FunMeltMIP.log
# for j in $(seq 1 120); do
# matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; addpath('$ISSM_DIR/lib'); outputDir='$PBS_JOBFS'; numberOfWorkers=$PBS_NCPUS; run('meltMip_ensemble($steps, $j, $loadonly)'), quit" > FuncMelt_${j}.log
# done

