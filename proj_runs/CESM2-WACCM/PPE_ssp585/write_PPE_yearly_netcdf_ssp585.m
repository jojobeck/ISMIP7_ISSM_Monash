function md = write_PPE_yearly_netcdf_ssp585(X)
% Writes the ISMIP7 NetCDF output for ONE loop (yearly-restart) PPE member
% (P014-P020, see ppe_members.csv) of proj_run_PPE_yearly_ssp585.m.
%
% Yearly counterpart of proj_run_PPE_frozenfront_ssp585.m's own step 4
% (WriteNetCDF): same output directory pattern, same meta struct, same
% write_ismip7_2d_projection / write_ismip7_scalar_projection helpers --
% the P-number is used as meta.set_counter, e.g.
%   libmassbffl_AIS_Monash_ISSM_m001_CESM2-WACCM_f001_ssp585_P014_2015-2300.nc
% written to postprocessed_data/CESM2-WACCM/ssp585_PPE_P014/.
%
% Loading mirrors ../ssp585/analyze_proj_ssp585_yearly.m step 2
% (WriteISMIP6_NetCDF_yearly): every per-year model file is loaded and
% their TransientSolution concatenated in year order, and time_range is
% computed from the actual data. The shared
% proj_runs/functions/load_and_concatenate_yearly_models.m cannot be reused
% directly because its discovery pattern is hardcoded to the core run's
% AIS_ISMIP7_Proj_<SCENARIO>_yearly_<year>.mat names, whereas this folder
% uses AIS_ISMIP7_PPE_<P_number>_yearly_<year>.mat -- hence the local
% helper below (same pattern as analyze_PPE_ssp585_SLR.m's own
% load_ppe_yearly_vaf_series).
%
% Usage (run from proj_runs/CESM2-WACCM/PPE_ssp585/, ISSM on the path):
%   write_PPE_yearly_netcdf_ssp585('p014')
% or for several members at once via launch_all_PPE_yearly_netcdf.sh.

    if nargin < 1 || isempty(X)
        error('X is required: ''p014''|''p015''|''p016''|''p017''|''p018''|''p019''|''p020''.');
    end

    % ---- PPE member lookup: X -> P_number -----------------------------------
    switch X
        case {'p014', 'p015', 'p016', 'p017', 'p018', 'p019', 'p020'}
            P_number = upper(X);
        otherwise
            error('Unrecognised X ''%s'' -- expected ''p014''..''p020''.', X);
    end
    % -----------------------------------------------------------------------

    CMIP_MODEL = 'CESM2-WACCM';
    SCENARIO   = 'ssp585';
    proj_root  = './../../../';
    modeldir   = ['./Models_' P_number '/'];

    addpath('./../../functions');
    addpath('./../../../init/scripts');

    md = load_and_concatenate_ppe_yearly_models(modeldir, P_number);

    outdir = [proj_root 'postprocessed_data/' CMIP_MODEL '/' SCENARIO '_PPE_' P_number '/'];
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    % time_range from the ACTUAL concatenated data -- same filter/rounding
    % convention as analyze_proj_ssp585_yearly.m step 2.
    t_raw    = [md.results.TransientSolution.time];
    keep     = abs(t_raw - round(t_raw)) < 0.05;
    t_annual = round(t_raw(keep));
    time_yr  = t_annual - 1;

    meta                    = struct();
    meta.experiment_id      = SCENARIO;
    meta.set_counter        = P_number;
    meta.time_range         = sprintf('%d-%d', min(time_yr), max(time_yr));
    meta.ESM_id             = CMIP_MODEL;
    meta.forcing_member_id  = 'f001';
    meta.ISM_member_id      = 'm001';

    [cfflux_tot, glflux_tot] = write_ismip7_2d_projection(md, outdir, meta);
    write_ismip7_scalar_projection(md, outdir, meta, cfflux_tot, glflux_tot);
    fprintf('NetCDF output for %s (%s) written to %s\n', P_number, meta.time_range, outdir);
end

function md = load_and_concatenate_ppe_yearly_models(modeldir, P_number)
% Same as proj_runs/functions/load_and_concatenate_yearly_models.m /
% discover_yearly_files.m, but globbing AIS_ISMIP7_PPE_<P_number>_yearly_*.mat.
% Unlike the VAF plot loader in analyze_PPE_ssp585_SLR.m, this one warns
% about missing years -- a gap would end up in the submitted NetCDF.
    pattern = [modeldir 'AIS_ISMIP7_PPE_' P_number '_yearly_*.mat'];
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

    missing = setdiff(years(1):years(end), years);
    if ~isempty(missing)
        warning('Yearly model files are missing for year(s): %s -- the NetCDF will have a gap there.', ...
                mat2str(missing));
    end
    fprintf('Found %d yearly model files for %s: %d -> %d\n', numel(years), P_number, years(1), years(end));

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
