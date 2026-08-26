#!/bin/bash
#PBS -q normal
#PBS -l ncpus=16
#PBS -l walltime=04:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -m ae
#PBS -l wd
#PBS -o SubmitGreeneRetreat.outlog
#PBS -e SubmitGreeneRetreat.errlog
#PBS -l storage=gdata/au88+gdata/k10

# Runs greene_icemask_retreat_rate.py then plot_greene_icemask_retreat.py.
#
# Memory: the script no longer builds a dense (23, 12161, 12161) stack --
# an earlier version did (~13.6 GB), and nanmedian/nanpercentile on an
# array that size need a full internal working copy to sort, which OOM-
# killed a 64 GB job. It's rebuilt on sparse retreat points now (a few
# thousand to tens of thousands per interval, not 147 million), so actual
# peak usage is more like a few GB. mem=190GB/ncpus=16 below is far more
# than needed at this point -- left as-is rather than silently changed;
# safe to size back down (e.g. mem=16GB, ncpus=4) for faster queueing.
#
# Adjust `conda activate base` below to your actual environment name if
# it isn't `base` -- this just needs scipy/xarray/pandas/numpy/matplotlib
# (cmocean optional, plotting falls back to matplotlib's Reds if absent).

source /g/data/k10/jb1863/mambaforge/etc/profile.d/conda.sh
conda activate base

cd "$PBS_O_WORKDIR" || exit 1

# python greene_icemask_retreat_rate.py
python plot_greene_icemask_retreat.py
