function [time, vaf] = load_yearly_vaf_series(modeldir, SCENARIO)
% LOAD_YEARLY_VAF_SERIES  Scalar-only loader for the VAF/SLE check.
%
%   modeldir  : directory containing AIS_ISMIP7_Proj_<SCENARIO>_yearly_<year>.mat
%   SCENARIO  : e.g. 'ssp585', 'ssp126'.
%
%   Returns time (1xN, ISSM-internal time) and vaf (1xN,
%   IceVolumeAboveFloatationScaled), ONE entry per yearly file (the real
%   annual snapshot, safety step dropped -- see below), concatenated in
%   year order.
%
%   Same file discovery as load_and_concatenate_yearly_models
%   (discover_yearly_files), but extracts ONLY these two scalar fields
%   from each file's TransientSolution and discards everything else
%   (mesh, geometry, every other results field) immediately after each
%   load. A VAF/SLE curve genuinely needs to visit every year (skipping
%   years would leave gaps in the plotted line), but it never needed the
%   full per-year spatial fields (Thickness, Vel, Vx, Vy, Base, Surface,
%   Bed, ...) that load_and_concatenate_yearly_models holds concatenated
%   across every single year for that purpose alone.
%
%   Each yearly loop iteration's solve produces TWO TransientSolution
%   entries, not one: a sub-annual "safety" step a fraction of a timestep
%   after start_time, plus the real annual snapshot at final_time
%   (year+1) -- same reason write_ismip7_2d.m / write_ismip7_scalar.m
%   filter their own (much larger, multi-year) concatenated time array
%   with keep = abs(t_raw-round(t_raw))<0.05 before writing: the safety
%   step's time is not (near-)integer, the real snapshot's is. This
%   function applies the SAME filter per-file, so a run of N years
%   returns N entries here, not 2N.
%
%   NOTE: each yearly file's md is still an ISSM classdef object, and
%   MATLAB's matfile() cannot partially load a classdef object -- so each
%   file is still read fully off disk (this doesn't reduce per-file I/O,
%   which would require changing how the yearly run script itself saves
%   its output). What this avoids is the O(years) growing full-field
%   concatenation this step never used, which was the larger, easily
%   avoidable cost.

    [years, names] = discover_yearly_files(modeldir, SCENARIO);
    fprintf('Found %d yearly model files: %d -> %d\n', numel(years), years(1), years(end));

    time = [];
    vaf  = [];
    for i = 1:numel(names)
        loaded = load([modeldir names{i}], 'md');
        ts    = loaded.md.results.TransientSolution;
        t_raw = [ts.time];
        keep  = abs(t_raw - round(t_raw)) < 0.05;   % drop the sub-annual safety step
        vaf_raw = [ts.IceVolumeAboveFloatationScaled];
        time = [time, t_raw(keep)]; %#ok<AGROW>
        vaf  = [vaf,  vaf_raw(keep)]; %#ok<AGROW>
        clear loaded ts t_raw keep vaf_raw;
    end
end
