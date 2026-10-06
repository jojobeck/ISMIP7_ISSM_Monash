#!/bin/bash
# Dispatcher for the frozen-front PPE members ONLY (P001-P003, see
# ppe_members.csv for what each one is) -- each a single, continuous
# 2015-2100 transient (no collapse-mask loop), submitted via
# submit_proj_PPE_frozenfront_ssp370.sh. NOT the loop (yearly-restart)
# members (P004-P009) -- those are launched separately via
# launch_all_PPE_yearly.sh.
#
# No separate P-number is reserved for the core/baseline ssp370 run --
# it has no P-number at all -- so PPE numbering starts at P001 here.
#
# proj_run_PPE_frozenfront_ssp370.m now has the NetCDF-writing step split
# OUT of the gather (loadonly=1) step -- gathering only saves the solved
# model, and WriteNetCDF (organizer step 2) is its own separate run, to be
# submitted only after step 1's gather has finished and saved. This
# dispatcher's MODE argument picks which of the three phases to submit:
#   submit  -> steps=[1], loadonly=0  (solve, run on PBS)
#   gather  -> steps=[1], loadonly=1  (save the solved model -- NO NetCDF)
#   netcdf  -> steps=[2]              (write the ISMIP6 NetCDF, needs gather done)
#
# One explicit qsub line per member, NOT a CSV-driven loop -- so you can
# comment out individual members (prefix the line with '#') to resubmit
# only some of them, without touching the others. Each line submits ONE
# independent PBS job via qsub -v X=<id>,steps=<...>,loadonly=<0|1> to
# submit_proj_PPE_frozenfront_ssp370.sh, with -o/-e overridden per member
# so concurrent jobs don't clobber each other's PBS log files.
#
# qsub calls are always independent -- there is no blocking between them.
# Each qsub returns its own new job ID immediately and the job queues on
# PBS; submitting PPE_frozenfront_p001 does not wait for, or get blocked
# by, any earlier or later qsub call for a different (or even the same)
# X, unless you explicitly add a `-W depend=...` dependency (not used
# here).
#
# This script itself finishes in seconds: it only calls qsub for the
# lines left uncommented and does not wait for any of them to complete,
# so it does NOT need to stay running -- run it directly on a login node
# (no need to qsub this dispatcher itself).
#
# Usage (run from this folder):
#   ./launch_all_PPE_frozenfront.sh            # = submit: solve every uncommented member
#   ./launch_all_PPE_frozenfront.sh submit     # same as above, explicit
#   ./launch_all_PPE_frozenfront.sh gather     # gather + save every uncommented member (no NetCDF)
#   ./launch_all_PPE_frozenfront.sh netcdf     # write NetCDF for every uncommented member
#                                               #   (run only once that member's own gather above has finished)

set -e
cd "$(dirname "$0")"

mode=${1:-submit}
case "$mode" in
    submit) steps='[1]'; loadonly=0 ;;
    gather) steps='[1]'; loadonly=1 ;;
    netcdf) steps='[2]'; loadonly=0 ;;   # loadonly unused by step 2, passed for a consistent qsub call
    *) echo "ERROR: unknown mode '$mode' -- expected submit|gather|netcdf" >&2; exit 1 ;;
esac

submit_member() {
    local member="$1"
    echo "Submitting PPE_frozenfront_${member} (mode=${mode}, steps=${steps}, loadonly=${loadonly})"
    qsub -v "X=${member},steps=${steps},loadonly=${loadonly}" \
         -N "PPE_frozenfront_${member}_${mode}" \
         -o "SubmitPPEfrozenfront_${member}_${mode}.outlog" \
         -e "SubmitPPEfrozenfront_${member}_${mode}.errlog" \
         submit_proj_PPE_frozenfront_ssp370.sh
}

# ---- Comment out any line below to skip that member on this run -------
submit_member p001   # K = mode (same K as core run)
submit_member p002   # K = 95th percentile
submit_member p003   # K = 5th percentile

echo ""
echo "Done. Check with: qstat -u \$USER"
