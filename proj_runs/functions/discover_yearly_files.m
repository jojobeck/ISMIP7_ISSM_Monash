function [years, names] = discover_yearly_files(modeldir, SCENARIO)
% DISCOVER_YEARLY_FILES  List every yearly model .mat file for a
%   yearly-restart run, sorted numerically by year, WITHOUT loading any
%   of them (fast -- only reads directory entries/filenames).
%
%   modeldir  : directory containing AIS_ISMIP7_Proj_<SCENARIO>_yearly_<year>.mat
%   SCENARIO  : e.g. 'ssp585', 'ssp126'.
%
%   Returns years (1xN, sorted ascending) and names (1xN cell of
%   filenames, same order) -- the small companion
%   appliedcollapse_<year>.mat files are excluded.
%
%   Shared by load_and_concatenate_yearly_models, load_yearly_vaf_series,
%   and load_yearly_at_years so the discovery/parsing/sorting logic lives
%   in exactly one place.

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

    expected = years(1):years(end);
    missing  = setdiff(expected, years);
    if ~isempty(missing)
        warning('Yearly model files are missing for year(s): %s -- some downstream steps may have a gap there.', ...
                mat2str(missing));
    end
end
