function md = analyze_proj_ssp585_yearly(steps)
% Post-processing / output-analysis for the yearly-restart ssp585 run
% (proj_run_CESM_WACCM_ssp585_2015_2300_yearly.m).
%
% NEW FILE -- neither the original proj_run_CESM_WACCM_ssp585_2015_2300.m
% nor the yearly run script itself are modified.
%
% Mirrors steps 7 (VAFContinuityCheck) and 9 (WriteISMIP6_NetCDF) of the
% ORIGINAL script exactly -- same VAF/SLE formula, same direct (unfiltered)
% concatenation of TransientSolution arrays, same NetCDF-writing calls --
% just generalised from 2 fixed segments to however many individual
% Models_yearly/*.mat files actually exist. load_and_concatenate_yearly_models
% (proj_runs/functions/) discovers them automatically by scanning that
% directory (not hardcoded to a specific start/end year), so this works
% whether the yearly run is complete (full 2015-2300) or only partially run
% so far (e.g. stopped at 2100) -- other scenarios sharing this same
% analyze-script pattern may only ever run to a shorter final year, and
% nothing here assumes otherwise (see step 2's time_range, computed from
% the actual data rather than hardcoded, for the one place this used to
% matter).
%
% Step map:
%   1  VAFCheck_yearly           load+concatenate all yearly models,
%                                 compute the VAF/SLE curve (identical
%                                 formula to the original's step 7), save
%                                 figure to postprocessed_data/figures/...
%   2  WriteISMIP6_NetCDF_yearly same concatenation, writes ISMIP6 NetCDFs
%                                 via the SAME (untouched)
%                                 proj_runs/functions/write_ismip7_2d_projection.m
%                                 / write_ismip7_scalar_projection.m used
%                                 by the original script's step 9 --
%                                 output goes to the SAME directory the
%                                 original used
%                                 (postprocessed_data/CESM2-WACCM/ssp585/)
%                                 using experiment_id='ssp585' and
%                                 set_counter='C007' -- SAME labels as the
%                                 original (ungated) run, deliberately: the
%                                 yearly floating-gated run is intended to
%                                 SUBSTITUTE that data as the real
%                                 ssp585/C007 submission going forward.
%                                 The original run's own NetCDF output was
%                                 moved to
%                                 postprocessed_data/CESM2-WACCM/ssp585_UnconstrainedCollapse/
%                                 (matching this project's existing
%                                 ssp585_test / ssp585_dT0 sibling-naming
%                                 convention) before this script's step 2
%                                 was first run, so writing into ssp585/
%                                 here does not destroy that prior data --
%                                 it's preserved under the renamed path.
%   3  CalvingFrontEvolution_yearly  mesh-native calving-front isoline
%                                 (MaskIceLevelset==0, via isoline() --
%                                 same technique as
%                                 plot_calvingfront_model_ssp585.m's own
%                                 plot_calvingfront_evolution), every 30
%                                 years across the whole combined run
%                                 (plus the final available year, even if
%                                 it doesn't fall on an exact 30-year mark
%                                 from the first year), colour-coded by
%                                 year over the initial ice-mask
%                                 background. Saves
%                                 CalvingFront_isoline_yearly_....png to
%                                 postprocessed_data/figures/...
%
% Usage (run from proj_runs/CESM2-WACCM/ssp585/, ISSM already on the
% MATLAB path -- pure post-processing, no PBS submission/waitonlock
% involved at all):
%   analyze_proj_ssp585_yearly([1])       % VAF check only
%   analyze_proj_ssp585_yearly([2])       % NetCDF write only
%   analyze_proj_ssp585_yearly([3])       % calving-front evolution plot only
%   analyze_proj_ssp585_yearly([1 2 3])   % all three (default)

    if nargin < 1 || isempty(steps)
        steps = [1 2 3];
    end

    CMIP_MODEL = 'CESM2-WACCM';
    SCENARIO   = 'ssp585';
    proj_root  = './../../../';
    modeldir   = './Models_yearly/';

    addpath('./../../functions');
    addpath('./../../../init/scripts');

    org = organizer('repository', modeldir, ...
                    'prefix', ['AIS_ISMIP7_Proj_' SCENARIO '_yearly_'], ...
                    'steps', steps, 'color', '34;47;2');
    clear steps;

    % load_and_concatenate_yearly_models (proj_runs/functions/) -- shared
    % across every scenario's own analyze_proj_<scenario>_yearly.m, since
    % this discovery/concatenation logic has nothing scenario-specific in
    % it beyond the SCENARIO string passed in. Loaded once here,
    % unconditionally, regardless of which steps were requested below --
    % all three steps need it and re-loading/re-concatenating potentially
    % hundreds of individual yearly .mat files per step would be wasteful.
    md = load_and_concatenate_yearly_models(modeldir, SCENARIO);

    % ================================================================= Step 1
    if perform(org, 'VAFCheck_yearly') % {{{
        % Identical formula to the original script's step 7
        % (VAFContinuityCheck), applied to the single combined (all-years)
        % TransientSolution instead of 2 segments -- no junction-matching
        % plot needed since there's only one continuous array here, not
        % two independently-loaded segments to check against each other.
        time = [md.results.TransientSolution.time];
        vaf  = [md.results.TransientSolution.IceVolumeAboveFloatationScaled];

        % IceVolumeAboveFloatationScaled is in m^3 (ice volume).
        % SLE [m] = -delta_VAF [m^3] * rho_ice / rho_sw / A_ocean
        rho_ice = 917;       % kg/m^3
        rho_sw  = 1028;      % kg/m^3
        A_ocean = 3.625e14;  % m^2 (ISMIP6 standard ocean area)
        vaf_ref = vaf(1);
        sle = -(vaf - vaf_ref) * rho_ice / rho_sw / A_ocean;

        fig = figure('visible', 'off', 'Position', [0 0 900 400]);
        plot(time, sle, 'b-', 'LineWidth', 1.2);
        xlabel('Year');
        ylabel('Sea level contribution (m SLE)');
        title(sprintf('%s %s (yearly-restart run) -- VAF SLE, %d years', ...
                      CMIP_MODEL, SCENARIO, numel(time)), 'Interpreter', 'none');
        grid on;

        figdir = [proj_root 'postprocessed_data/figures/' CMIP_MODEL '/' SCENARIO '/'];
        if ~exist(figdir, 'dir'), mkdir(figdir); end
        figname = fullfile(figdir, sprintf('VAF_continuity_yearly_%s_%s.png', CMIP_MODEL, SCENARIO));
        saveas(fig, figname);
        close(fig);
        fprintf('Saved: %s\n', figname);
    end % }}}

    % ================================================================= Step 2
    if perform(org, 'WriteISMIP6_NetCDF_yearly') % {{{
        % Identical to the original script's step 9 (WriteISMIP6_NetCDF):
        % passes the concatenated TransientSolution straight to the SAME
        % (untouched) write_ismip7_2d_projection / write_ismip7_scalar_projection
        % helpers, with no pre-filtering -- those functions already handle
        % whatever sub-annual "safety step" entries ISSM may add, exactly
        % as they already must for the original's own 2-segment
        % concatenation.
        %
        % Writes into the SAME path/labels as the original script's own
        % output (postprocessed_data/CESM2-WACCM/ssp585/, experiment_id
        % 'ssp585', set_counter 'C007') -- deliberately substituting that
        % data, per instruction. The original run's own output was already
        % moved to ssp585_UnconstrainedCollapse/ (see module docstring)
        % before this was first run, so nothing is destroyed -- but
        % running this again WILL overwrite whatever is currently in
        % ssp585/, including results from any earlier run of this same
        % yearly script.
        outdir = [proj_root 'postprocessed_data/CESM2-WACCM/' SCENARIO '/'];
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        % time_range is computed from the ACTUAL concatenated data (not
        % hardcoded) -- write_ismip7_2d.m / write_ismip7_scalar.m bake
        % meta.time_range into both the output filename and a NetCDF
        % global attribute, so a hardcoded '2015-2300' would be wrong for
        % any run that doesn't cover the full span (e.g. a scenario that
        % only runs to 2100). Same filter/rounding convention those
        % functions already use internally.
        t_raw    = [md.results.TransientSolution.time];
        keep     = abs(t_raw - round(t_raw)) < 0.05;
        t_annual = round(t_raw(keep));
        time_yr  = t_annual - 1;

        meta                    = struct();
        meta.experiment_id      = SCENARIO;   % 'ssp585' -- same as the original run
        meta.set_counter        = 'C007';     % same as the original run
        meta.time_range         = sprintf('%d-%d', min(time_yr), max(time_yr));
        meta.ESM_id             = CMIP_MODEL;
        meta.forcing_member_id  = 'f001';
        meta.ISM_member_id      = 'm001';

        [cfflux_tot, glflux_tot] = write_ismip7_2d_projection(md, outdir, meta);
        write_ismip7_scalar_projection(md, outdir, meta, cfflux_tot, glflux_tot);
        fprintf('NetCDF output written to %s\n', outdir);
    end % }}}

    % ================================================================= Step 3
    if perform(org, 'CalvingFrontEvolution_yearly') % {{{
        % Mesh-native calving-front isoline (MaskIceLevelset==0, the
        % ice/no-ice boundary -- NOT MaskOceanLevelset, the grounding
        % line) every 30 years, colour-coded by year, over a grayscale
        % background of the initial (first-loaded-file) ice mask. Same
        % isoline() technique as plot_calvingfront_model_ssp585.m's own
        % plot_calvingfront_evolution (isoline() gives the exact mesh-edge
        % crossing, not a raster-interpolated contour), just spaced every
        % 30 years instead of every 10, and driven off the concatenated
        % yearly-restart TransientSolution instead of the original's 2
        % fixed segments.
        time      = [md.results.TransientSolution.time];
        nearest_t = @(yr) find(abs(time - yr) == min(abs(time - yr)), 1);

        start_yr = round(min(time));
        end_yr   = round(max(time));
        cf_years = start_yr:30:end_yr;
        if cf_years(end) ~= end_yr
            cf_years(end+1) = end_yr;   % always show the final available year too
        end
        n_cf    = length(cf_years);
        cmap_cf = jet(n_cf);

        hi_mod_all = cell(n_cf, 1);
        for k = 1:n_cf
            ti = nearest_t(cf_years(k));
            hi_mod_all{k} = isoline(md, md.results.TransientSolution(ti).MaskIceLevelset, ...
                                     'value', 0, 'output', 'matrix');
        end

        figure('visible', 'off', 'Position', [0 0 1000 800]);
        plotmodel(md, 'figure', gcf, 'visible', 'off', 'data', md.mask.ice_levelset, ...
                  'colormap', gray, 'caxis', [-1 1], 'title', '', 'colorbar', 0);
        hold on;
        for k = 1:n_cf
            col = cmap_cf(k, :);
            plot(hi_mod_all{k}(:,1), hi_mod_all{k}(:,2), '-', 'Color', col, 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('%d', cf_years(k)));
        end
        axis equal tight off;
        legend('Location', 'eastoutside', 'FontSize', 7, 'NumColumns', 1);
        title(sprintf('%s %s (yearly-restart run) -- calving front (ice-mask isoline), every 30 yr', ...
                      CMIP_MODEL, SCENARIO), 'Interpreter', 'none');

        figdir = [proj_root 'postprocessed_data/figures/' CMIP_MODEL '/' SCENARIO '/'];
        if ~exist(figdir, 'dir'), mkdir(figdir); end
        figname = fullfile(figdir, sprintf('CalvingFront_isoline_yearly_%s_%s.png', CMIP_MODEL, SCENARIO));
        saveas(gcf, figname);
        close(gcf);
        fprintf('Saved: %s\n', figname);
    end % }}}
end
