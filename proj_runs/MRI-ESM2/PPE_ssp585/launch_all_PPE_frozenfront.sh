#!/bin/bash
# Dispatcher for the frozen-front PPE members ONLY (P021-P023, see
# ppe_members.csv for what each one is) -- submitted via
# submit_proj_PPE_frozenfront_ssp585.sh. NOT the loop (yearly-restart)
# members (P024-P029) -- those are launched separately via
# launch_all_PPE_yearly.sh.
#
# No separate P-number is reserved for the core/baseline ssp585 run --
# it has no P-number at all. P-numbers continue on from the ssp370 PPE
# families (P001-P009 ssp370, P011-P019 CESM2-WACCM ssp585), so numbering here
# starts at P021.
#
# RESTART-SPLIT at 2150/2151 -- UNLIKE the ssp370 version of this
# dispatcher (which runs a single continuous 2015-2100 transient), this
# one's 2015-2300 span got OOM-killed as a single job, so
# proj_run_PPE_frozenfront_ssp585.m now runs in two ~140-year segments
# (see that script's own module docstring for the full rationale). This
# dispatcher's MODE argument picks which of the SIX phases to submit:
#   submit1  -> steps=[1], loadonly=0  (solve 2015-2150, run on PBS)
#   gather1  -> steps=[1], loadonly=1  (save the solved segment-1 model)
#   restart  -> steps=[2]              (save the 2151 restart state -- needs gather1 done)
#   submit2  -> steps=[3], loadonly=0  (solve 2151-2300, run on PBS -- needs restart done)
#   gather2  -> steps=[3], loadonly=1  (save the solved segment-2 model)
#   netcdf   -> steps=[4]              (write the ISMIP6 NetCDF -- needs gather2 done)
# Each member must go through all six phases in that order before its
# NetCDF is written -- run this dispatcher once per phase, left to right.
#
# One explicit qsub line per member, NOT a CSV-driven loop -- so you can
# comment out individual members (prefix the line with '#') to resubmit
# only some of them, without touching the others. Each line submits ONE
# independent PBS job via qsub -v X=<id>,steps=<...>,loadonly=<0|1> to
# submit_proj_PPE_frozenfront_ssp585.sh, with -o/-e overridden per member
# so concurrent jobs don't clobber each other's PBS log files.
#
# qsub calls are always independent -- there is no blocking between them.
# Each qsub returns its own new job ID immediately and the job queues on
# PBS; submitting PPE_frozenfront_p021 does not wait for, or get blocked
# by, any earlier or later qsub call for a different (or even the same)
# X, unless you explicitly add a `-W depend=...` dependency (not used
# here) -- so e.g. calling this dispatcher with 'restart' before every
# member's own 'submit1'/'gather1' job has actually finished will just
# fail for that member (loadmodel can't find the step-1 output yet), not
# block waiting for it.
#
# This script itself finishes in seconds: it only calls qsub for the
# lines left uncommented and does not wait for any of them to complete,
# so it does NOT need to stay running -- run it directly on a login node
# (no need to qsub this dispatcher itself).
#
# Usage (run from this folder):
#   ./launch_all_PPE_frozenfront.sh submit1    # segment 1: solve every uncommented member
#   ./launch_all_PPE_frozenfront.sh gather1    # segment 1: gather + save (once segment-1 solves finish)
#   ./launch_all_PPE_frozenfront.sh restart    # save each member's own 2151 restart state
#   ./launch_all_PPE_frozenfront.sh submit2    # segment 2: solve (once restart states are saved)
#   ./launch_all_PPE_frozenfront.sh gather2    # segment 2: gather + save (once segment-2 solves finish)
#   ./launch_all_PPE_frozenfront.sh netcdf     # write NetCDF (once both gathers are done)

set -e
cd "$(dirname "$0")"

mode=${1:?"ERROR: mode required -- one of submit1|gather1|restart|submit2|gather2|netcdf"}
case "$mode" in
    submit1) steps='[1]'; loadonly=0 ;;
    gather1) steps='[1]'; loadonly=1 ;;
    restart) steps='[2]'; loadonly=0 ;;   # loadonly unused by step 2, passed for a consistent qsub call
    submit2) steps='[3]'; loadonly=0 ;;
    gather2) steps='[3]'; loadonly=1 ;;
    netcdf)  steps='[4]'; loadonly=0 ;;   # loadonly unused by step 4
    *) echo "ERROR: unknown mode '$mode' -- expected submit1|gather1|restart|submit2|gather2|netcdf" >&2; exit 1 ;;
esac

submit_member() {
    local member="$1"
    echo "Submitting PPE_frozenfront_${member} (mode=${mode}, steps=${steps}, loadonly=${loadonly})"
    qsub -v "X=${member},steps=${steps},loadonly=${loadonly}" \
         -N "PPE_frozenfront_${member}_${mode}" \
         -o "SubmitPPEfrozenfront_${member}_${mode}.outlog" \
         -e "SubmitPPEfrozenfront_${member}_${mode}.errlog" \
         submit_proj_PPE_frozenfront_ssp585.sh
}

# ---- Comment out any line below to skip that member on this run -------
submit_member p021   # K = mode (same K as core run)
submit_member p022   # K = 95th percentile
submit_member p023   # K = 5th percentile

echo ""
echo "Done. Check with: qstat -u \$USER"
