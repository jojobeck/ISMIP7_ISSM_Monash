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
% Models_yearly/*.mat files actually exist. Each step below loads ONLY
% what it needs (proj_runs/functions/), not a single shared up-front load:
%   - VAFCheck_yearly uses load_yearly_vaf_series -- visits every year
%     (a VAF/SLE curve can't skip years without leaving a gap) but pulls
%     only 2 scalars per year, not the full per-year field set.
%   - WriteISMIP6_NetCDF_yearly uses load_and_concatenate_yearly_models --
%     the one step that genuinely needs every field for every year (the
%     ISMIP6 gridded output requires it), and correspondingly the slow
%     one -- this is inherent, not something this refactor changes.
%   - CalvingFrontEvolution_yearly uses load_yearly_at_years -- loads only
%     the ~10 sparse years (every 30) it actually plots, not the whole run.
% Each loader discovers available years by scanning Models_yearly/ (not
% hardcoded to a specific start/end year), so this works whether the
% yearly run is complete (full 2015-2300) or only partially run so far
% (e.g. stopped at 2100) -- other scenarios sharing this same
% analyze-script pattern may only ever run to a shorter final year, and
% nothing here assumes otherwise (see step 2's time_range, computed from
% the actual data rather than hardcoded, for the one place this used to
% matter). Because each step is now independent, running e.g. just [1] or
% just [3] no longer pays the cost of the full multi-hundred-file
% concatenation that WriteISMIP6_NetCDF_yearly alone actually needs.
%
% Step map:
%   1  VAFCheck_yearly           compute the VAF/SLE curve (identical
%                                 formula to the original's step 7), save
%                                 figure to postprocessed_data/figures/...
%   2  WriteISMIP6_NetCDF_yearly writes ISMIP6 NetCDFs via the SAME
%                                 (untouched)
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
%                                 years across the whole run (plus the
%                                 final available year, even if it doesn't
%                                 fall on an exact 30-year mark from the
%                                 first year), colour-coded by year over
%                                 the initial ice-mask background. Saves
%                                 CalvingFront_isoline_yearly_....png to
%                                 postprocessed_data/figures/...
%
% NOTE: this function declares an md output for consistency with this
% project's other scripts, but only step 2 actually assigns it -- fine
% for the usage pattern below (nargout==0), but a caller that runs only
% step 1 or 3 AND captures the output would hit an "output not assigned"
% error.
%
% Usage (run from proj_runs/CESM2-WACCM/ssp585/, ISSM already on the
% MATLAB path -- pure post-processing, no PBS submission/waitonlock
% involved at all):
%   analyze_proj_ssp585_yearly([1])       % VAF check only (fast)
%   analyze_proj_ssp585_yearly([2])       % NetCDF write only (slow -- needs every field/year)
%   analyze_proj_ssp585_yearly([3])       % calving-front evolution plot only (fast)
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

    % ================================================================= Step 1
    if perform(org, 'VAFCheck_yearly') % {{{
        % load_yearly_vaf_series (proj_runs/functions/) -- scalar-only
        % loader, visits every year but never builds the full concatenated
        % TransientSolution this step doesn't need. Identical VAF/SLE
        % formula to the original script's step 7 (VAFContinuityCheck).
        [time, vaf] = load_yearly_vaf_series(modeldir, SCENARIO);

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
        % load_and_concatenate_yearly_models (proj_runs/functions/) -- the
        % ONE step that genuinely needs the full per-year field set (the
        % ISMIP6 gridded output requires it), so this is the expensive
        % load, done only here.
        %
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
        md = load_and_concatenate_yearly_models(modeldir, SCENARIO);

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
        % load_yearly_at_years (proj_runs/functions/) -- loads only the
        % sparse years this plot actually needs (every 30, plus the final
        % available year), not the whole run. Mesh-native calving-front
        % isoline (MaskIceLevelset==0, the ice/no-ice boundary -- NOT
        % MaskOceanLevelset, the grounding line), colour-coded by year,
        % over a grayscale background of the earliest-loaded-year ice
        % mask. Same isoline() technique as
        % plot_calvingfront_model_ssp585.m's own plot_calvingfront_evolution
        % (isoline() gives the exact mesh-edge crossing, not a
        % raster-interpolated contour).
        [avail_years, ~] = discover_yearly_files(modeldir, SCENARIO);
        start_yr = min(avail_years);
        end_yr   = max(avail_years);
        cf_years = start_yr:30:end_yr;
        if cf_years(end) ~= end_yr
            cf_years(end+1) = end_yr;   % always show the final available year too
        end

        % loaded_years is the ACTUAL year loaded per snapshot (snapped to
        % whatever files exist -- may differ from cf_years if a requested
        % year isn't available) -- used for the legend labels below so
        % they always match what was actually plotted.
        [md_cf, loaded_years] = load_yearly_at_years(modeldir, SCENARIO, cf_years);
        n_cf    = numel(loaded_years);
        cmap_cf = jet(n_cf);

        hi_mod_all = cell(n_cf, 1);
        for k = 1:n_cf
            hi_mod_all{k} = isoline(md_cf, md_cf.results.TransientSolution(k).MaskIceLevelset, ...
                                     'value', 0, 'output', 'matrix');
        end

        figure('visible', 'off', 'Position', [0 0 1000 800]);
        plotmodel(md_cf, 'figure', gcf, 'visible', 'off', 'data', md_cf.mask.ice_levelset, ...
                  'colormap', gray, 'caxis', [-1 1], 'title', '', 'colorbar', 0);
        hold on;
        for k = 1:n_cf
            col = cmap_cf(k, :);
            plot(hi_mod_all{k}(:,1), hi_mod_all{k}(:,2), '-', 'Color', col, 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('%d', loaded_years(k)));
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
