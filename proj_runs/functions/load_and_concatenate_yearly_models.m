function md = load_and_concatenate_yearly_models(modeldir, SCENARIO)
% LOAD_AND_CONCATENATE_YEARLY_MODELS  Discover and concatenate all
%   per-year model .mat files written by a *_yearly.m projection script
%   (e.g. proj_run_CESM_WACCM_ssp585_2015_2300_yearly.m step 4).
%
%   modeldir  : directory containing AIS_ISMIP7_Proj_<SCENARIO>_yearly_<year>.mat
%               files (e.g. './Models_yearly/').
%   SCENARIO  : scenario string, e.g. 'ssp585', 'ssp126'.
%
%   Returns md with TransientSolution equal to the concatenation of every
%   discovered year's TransientSolution, in year order -- exactly like the
%   original (non-yearly) projection scripts' own concatenation of their
%   fixed 2-segment TransientSolution arrays, just generalised to N
%   arbitrary years. Works whether the yearly run is complete or only
%   partially run so far (e.g. stopped at some intermediate year) -- years
%   are discovered from disk, not assumed.
%
%   Scans for AIS_ISMIP7_Proj_<SCENARIO>_yearly_*.mat, skipping the small
%   companion appliedcollapse_<year>.mat files saved alongside them, and
%   sorts NUMERICALLY by year (not alphabetically -- not actually
%   ambiguous while all years are 4 digits, but explicit regardless).

    pattern = [modeldir 'AIS_ISMIP7_Proj_' SCENARIO '_yearly_*.mat'];
    files = dir(pattern);
    if isempty(files)
        error('No yearly model files found matching %s', pattern);
    end

    years = [];
    names = {};
    for i = 1:numel(files)
        fname = files(i).name;
        if contains(fname, 'appliedcollapse')
            continue;
        end
        tok = regexp(fname, '_yearly_(\d+)\.mat$', 'tokens', 'once');
        if isempty(tok)
            continue;
        end
        years(end+1) = str2double(tok{1}); %#ok<AGROW>
        names{end+1} = fname; %#ok<AGROW>
    end
    if isempty(years)
        error('No yearly model files (excluding appliedcollapse_*) found matching %s', pattern);
    end

    [years, order] = sort(years);
    names = names(order);
    fprintf('Found %d yearly model files: %d -> %d\n', numel(years), years(1), years(end));

    expected = years(1):years(end);
    missing  = setdiff(expected, years);
    if ~isempty(missing)
        warning('Yearly model files are missing for year(s): %s -- concatenated output will have a gap there.', ...
                mat2str(missing));
    end

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
