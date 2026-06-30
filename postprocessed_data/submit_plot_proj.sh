#!/bin/bash
#PBS -N plot_proj
#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/au88+scratch/au88
#PBS -P au88
#PBS -q normal
#PBS -j oe
#PBS -o /home/565/jb1863/SAEF/ISMIP7_ISSM_Monash/postprocessed_data/plot_proj.log

module purge
module load python3/3.11.7

cd /home/565/jb1863/SAEF/ISMIP7_ISSM_Monash/postprocessed_data

python3 plot_projection_runs.py
