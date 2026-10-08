#!/bin/bash
# Dispatcher for writing the ISMIP7 NetCDF of the loop (yearly-restart)
# PPE members (P014-P019, see ppe_members.csv) -- the yearly counterpart
# of `./launch_all_PPE_frozenfront.sh netcdf`. One explicit qsub line per
# member: comment out a line (prefix with '#') to skip that member.
# Each member is its own independent PBS job (submit_PPE_yearly_netcdf_ssp585.sh),
# writing to postprocessed_data/CESM2-WACCM/ssp585_PPE_P<NNN>/ with the
# P-number as set_counter, e.g.
#   lithk_AIS_Monash_ISSM_m001_CESM2-WACCM_f001_ssp585_P014_2015-2300.nc
#
# Finishes in seconds -- run directly on a login node.
#
# Usage (run from this folder):
#   ./launch_all_PPE_yearly_netcdf.sh

set -e
cd "$(dirname "$0")"

submit_member() {
    local member="$1"
    echo "Submitting PPE_${member}_netcdf"
    qsub -v "X=${member}" \
         -N "PPE_${member}_netcdf" \
         -o "SubmitPPEyearlyNetCDF_${member}.outlog" \
         -e "SubmitPPEyearlyNetCDF_${member}.errlog" \
         submit_PPE_yearly_netcdf_ssp585.sh
}

# ---- Comment out any line below to skip that member on this run -------
submit_member p014   # migration_max = Greene maximum, no eigencalving
# submit_member p015   # migration_max = Greene minimum, no eigencalving
submit_member p016   # migration_max = Greene mean,    + eigencalving
submit_member p017   # migration_max = Greene maximum, + eigencalving
# submit_member p018   # migration_max = Greene minimum, + eigencalving
submit_member p019   # migration_max = Greene mean,    SMB lapserate off
# P020 NetCDF is submitted on its own (once its yearly run is complete), from this folder:
#   qsub -v X=p020 -N PPE_p020_netcdf -o SubmitPPEyearlyNetCDF_p020.outlog -e SubmitPPEyearlyNetCDF_p020.errlog submit_PPE_yearly_netcdf_ssp585.sh

echo ""
echo "Done. Check with: qstat -u \$USER"
