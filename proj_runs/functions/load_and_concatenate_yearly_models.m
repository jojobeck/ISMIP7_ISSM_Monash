function md = load_and_concatenate_yearly_models(modeldir, SCENARIO)
% LOAD_AND_CONCATENATE_YEARLY_MODELS  Discover and concatenate ALL
%   per-year model .mat files written by a *_yearly.m projection script
%   (e.g. proj_run_CESM_WACCM_ssp585_2015_2300_yearly.m step 4), with the
%   FULL TransientSolution (every spatial field) for every year.
%
%   modeldir  : directory containing AIS_ISMIP7_Proj_<SCENARIO>_yearly_<year>.mat
%               files (e.g. './Models_yearly/').
%   SCENARIO  : scenario string, e.g. 'ssp585', 'ssp126'.
%
%   Returns md with TransientSolution equal to the concatenation of every
%   discovered year's TransientSolution, in year order -- exactly like the
%   original (non-yearly) projection scripts' own concatenation of their
%   fixed 2-segment TransientSolution arrays, just generalised to N
%   arbitrary years.
%
%   This is the EXPENSIVE loader (reads every yearly file in full, holds
%   every spatial field for every year in memory) -- only use it for
%   post-processing that genuinely needs the complete per-year field data,
%   e.g. gridding+writing the ISMIP6 NetCDF output. For anything that only
%   needs a couple of scalars per year (e.g. a VAF/SLE time series), use
%   the much lighter load_yearly_vaf_series instead. For anything that
%   only needs a handful of sparse years (e.g. a calving-front snapshot
%   every 30 years), use load_yearly_at_years instead -- both avoid
%   reading/holding years or fields this function's caller never actually
%   touches.

    [years, names] = discover_yearly_files(modeldir, SCENARIO);
    fprintf('Found %d yearly model files: %d -> %d\n', numel(years), years(1), years(end));

    md = [];
    for i = 1:numel(names)
        loaded = load([modeldir names{i}], 'md');
        if isempty(md)
            md = loaded.md;
        else
            md.results.TransientSolution = [md.results.TransientSolution, ...
                                             loaded.md.results.TransientSolution];
        end
        clear loaded;
    end
end
