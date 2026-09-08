function [md, loaded_years] = load_yearly_at_years(modeldir, SCENARIO, requested_years)
% LOAD_YEARLY_AT_YEARS  Load only SPECIFIC years of a yearly-restart run
%   (not the whole run) into one md, for post-processing that only needs
%   a handful of sparse snapshots -- e.g. a calving-front isoline every 30
%   years out of a run that can span up to 286 -- rather than every year.
%
%   modeldir        : directory containing AIS_ISMIP7_Proj_<SCENARIO>_yearly_<year>.mat
%   SCENARIO        : e.g. 'ssp585', 'ssp126'.
%   requested_years : vector of desired years, e.g. 2015:30:2300.
%
%   requested_years is snapped to the nearest ACTUALLY AVAILABLE year
%   (discovered from disk via discover_yearly_files, not loaded) and
%   de-duplicated, so callers can ask for round marks without needing to
%   know exactly which files exist or whether the run is complete/partial.
%
%   Returns md with TransientSolution containing ONLY the loaded years
%   (one entry per loaded year -- each yearly loop iteration's own file
%   only ever spans one year, so that file's final TransientSolution
%   entry IS simply that solve's own end state; there is no cross-file
%   nearest-time search needed), in year order, and loaded_years (1xN,
%   the actual years loaded -- use these for plot labels, not the
%   original requested_years, since a request can snap to a different
%   available year). md.mask/geometry/etc. reflect whatever the EARLIEST
%   loaded file's own md struct had, matching
%   load_and_concatenate_yearly_models' convention (if the earliest
%   requested year is the run's true first year, this is the true initial
%   state, e.g. AIS_state_2015).
%
%   Loads ONLY the requested files -- for a handful of sparse years, this
%   avoids reading and holding the other 90%+ of years a step like
%   CalvingFrontEvolution_yearly never needed.

    [avail_years, avail_names] = discover_yearly_files(modeldir, SCENARIO);

    target_years = unique(requested_years(:)', 'stable');
    snapped = zeros(size(target_years));
    for i = 1:numel(target_years)
        [~, idx] = min(abs(avail_years - target_years(i)));
        snapped(i) = avail_years(idx);
    end
    loaded_years = unique(sort(snapped));

    fprintf('Loading %d of %d available yearly model files (sparse: %s)\n', ...
            numel(loaded_years), numel(avail_years), mat2str(loaded_years));

    md = [];
    for i = 1:numel(loaded_years)
        idx = find(avail_years == loaded_years(i), 1);
        loaded = load([modeldir avail_names{idx}], 'md');
        this_ts = loaded.md.results.TransientSolution(end);   % this file's own final state
        if isempty(md)
            md = loaded.md;
            md.results.TransientSolution = this_ts;
        else
            md.results.TransientSolution = [md.results.TransientSolution, this_ts];
        end
        clear loaded this_ts;
    end
end
