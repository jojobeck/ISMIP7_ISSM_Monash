function ispresent = waitonlock(md)
%WAITONLOCK - local override, shadows ISSM trunk's src/m/solve/waitonlock.m
%
%   Logic copied verbatim (function renamed only) from ISSM's own
%   src/m/os/waitonlock_ongadi.m -- an existing, already-written fix for
%   exactly this problem that is never actually called by anything else
%   in the trunk (solve() unconditionally calls the generic waitonlock(),
%   not this one, so the fix sits unused unless something makes MATLAB
%   resolve the name 'waitonlock' to this file instead).
%
%   Why this is needed: the trunk's generic src/m/solve/waitonlock.m polls
%   for job completion over SSH whenever oshostname() ~= cluster.name --
%   which is exactly the case when MATLAB is already running inside a PBS
%   job on a Gadi compute node. Self-SSH from a Gadi compute/interactive
%   node back to itself returns Signal 127, so the wait fails outright.
%   This version polls the lock/outlog files directly on the local
%   filesystem (isfile()) instead, with no SSH involved at all.
%
%   HOW THE SHADOWING WORKS: this file must live in a directory added to
%   the MATLAB path AFTER devpath has already run (devpath adds the whole
%   $ISSM_DIR/src/m tree, including the trunk's own src/m/solve/waitonlock.m).
%   addpath() prepends by default, so calling addpath() on this file's
%   directory from inside proj_run_CESM_WACCM_ssp585_2015_2300_yearly.m
%   (which executes after devpath has already run in the calling shell)
%   places this copy earlier in the search path, so plain calls to
%   waitonlock(md) -- including solve()'s own internal call -- resolve to
%   THIS file instead of the trunk's broken one. Nothing in the shared
%   trunk install is modified.
%
%   UNITS: md.settings.waitonlock is SECONDS here (matching
%   waitonlock_ongadi.m's convention) -- NOT minutes, which is what the
%   trunk's generic waitonlock.m uses for the same field. Do not conflate
%   the two when setting md.settings.waitonlock in calling scripts.
%
%   Usage: ispresent = waitonlock(md)

executionpath = md.cluster.executionpath;
timelimit     = md.settings.waitonlock;   % seconds

% Same path logic as the trunk's waitonlock.m / waitonlock_ongadi.m
if isfield(md, 'private') && isfield(md.private, 'runtimename') && ~isempty(md.private.runtimename)
    dirname = md.private.runtimename;
else
    dirname = md.miscellaneous.name;
end
lockfilename = [executionpath '/' dirname '/' md.miscellaneous.name '.lock'];
logfilename  = [executionpath '/' dirname '/' md.miscellaneous.name '.outlog'];

disp('waitonlock (local Gadi override): waiting for job completion');
disp(['waiting for ' lockfilename ' (Ctrl+C to exit)']);

elapsedtime = 0;
ispresent   = 0;
starttime   = clock;

while ~ispresent && elapsedtime < timelimit
    pause(30);
    ispresent   = isfile(lockfilename) && isfile(logfilename);
    elapsedtime = etime(clock, starttime);
    if ~ispresent
        fprintf('\rchecking for job completion (time: %i min %i sec)', ...
                floor(elapsedtime/60), floor(mod(elapsedtime,60)));
    end
end
fprintf('\n');

if ~ispresent
    disp('Time limit exceeded. Increase md.settings.waitonlock (seconds).');
    disp('Load results manually with md=loadresultsfromcluster(md).');
    error('waitonlock error message: time limit exceeded');
end
