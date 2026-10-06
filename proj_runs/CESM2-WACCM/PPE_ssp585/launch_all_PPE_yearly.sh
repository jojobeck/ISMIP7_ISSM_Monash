#!/bin/bash
# Dispatcher for the loop (yearly-restart) PPE members ONLY (P014-P019,
# see ppe_members.csv for what each one is). Frozen-front members
# (P011-P013) are each a single, much shorter PBS solve and are still
# submitted directly via submit_proj_PPE_frozenfront_ssp585.sh / its own
# launch_all_PPE_frozenfront.sh -- they were not part of this
# restructuring, since the "stays running the whole time" problem this
# script solves is specific to the yearly-restart members' own
# many-submit-wait-gather-cycles-in-one-session design.
#
# No separate P-number is reserved for the core/baseline ssp585 run --
# it has no P-number at all. P-numbers continue on from the ssp370 PPE
# family (P001-P009) rather than restarting at P004, so numbering here
# starts at P014.
#
# One explicit qsub line per member, NOT a CSV-driven loop -- so you can
# comment out individual members (prefix the line with '#') to resubmit
# only some of them, without touching the others. Each line submits ONE
# independent PBS job via qsub -v X=<id> to submit_PPE_member_yearly.sh
# (see that script's own header for why X is passed this way rather than
# hardcoded -- the same template file is reused for every member), with
# -o/-e overridden per member so concurrent jobs don't clobber each
# other's PBS log files (the ssp370 version of this dispatcher does NOT
# do this yet -- its members all share one SubmitPPEyearly.outlog/.errlog
# pair, so running several of its members at once will clobber that
# shared pair; worth the same fix there too if you run several of them
# concurrently).
#
# qsub calls are always independent -- there is no blocking between them.
# Each qsub returns its own new job ID immediately and the job queues on
# PBS; submitting PPE_p014 does not wait for, or get blocked by, any
# earlier or later qsub call for a different (or even the same) X, unless
# you explicitly add a `-W depend=...` dependency (not used here).
#
# This script itself finishes in seconds: it only calls qsub for the
# lines left uncommented and does not wait for any of them to complete,
# so it does NOT need to stay running -- run it directly on a login node
# (no need to qsub this dispatcher itself).
#
# Usage (run from this folder):
#   ./launch_all_PPE_yearly.sh                       # all uncommented members, defaults (2015-2300, 1hr waitonlock)
#   ./launch_all_PPE_yearly.sh 2015 2300 3600         # explicit defaults
#   ./launch_all_PPE_yearly.sh 2080 2300 21600        # e.g. resume every uncommented member from a crash at 2080
#     (resuming assumes every member actually crashed/stopped at the same
#     year -- for resuming just ONE member from its own specific year,
#     submit that one directly instead:
#       qsub -v X=p016,start_year=2080 -N PPE_p016 submit_PPE_member_yearly.sh)

set -e
cd "$(dirname "$0")"

start_year=${1:-2015}
end_year=${2:-2300}
waitonlock_seconds=${3:-3600}

submit_member() {
    local member="$1"
    echo "Submitting PPE_${member} (years ${start_year}-${end_year})"
    qsub -v "X=${member},start_year=${start_year},end_year=${end_year},waitonlock_seconds=${waitonlock_seconds}" \
         -N "PPE_${member}" \
         -o "SubmitPPEyearly_${member}.outlog" \
         -e "SubmitPPEyearly_${member}.errlog" \
         submit_PPE_member_yearly.sh
}

# ---- Comment out any line below to skip that member on this run -------
submit_member p014   # migration_max = Greene maximum (96521.91 m/yr), no eigencalving
submit_member p015   # migration_max = Greene minimum (203.82 m/yr),   no eigencalving
submit_member p016   # migration_max = Greene mean (core value, 2863.78 m/yr), + eigencalving
submit_member p017   # migration_max = Greene maximum,                          + eigencalving
submit_member p018   # migration_max = Greene minimum,                          + eigencalving
submit_member p019   # migration_max = Greene mean (core value), no eigencalving, SMB feedback OFF

echo ""
echo "Done. Check with: qstat -u \$USER"
