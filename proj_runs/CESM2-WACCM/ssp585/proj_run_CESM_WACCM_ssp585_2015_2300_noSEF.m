function md = proj_run_CESM_WACCM_ssp585_2015_2300_noSEF(steps, loadonly)
% noSEF sensitivity: projection 2015-2300 without surface-elevation feedback.
% Uses SMBforcing() with mass_balance instead of SMBgradients(), same corrected SMB
% forcing as the main SEF run (proj_run_CESM_WACCM_ssp585_2015_2300.m).
%
% Run AFTER the main script steps 1-3 (forcing preprocessing), which writes
% the mat files consumed here.
%
% Uses the same organizer (repository, prefix) as the main script so that
% both the SEF and noSEF model files live in the same Models/ directory and
% can be loaded by name for the VAF comparison in step 4.
%
% Step map:
%   1  ProjRun_2015_2150_noSEF  transient 2015-2150, no SEF
%                               (loadonly=0 submit, =1 gather)
%   2  AIS_state_2151_noSEF     save end-state at 2151 (restart point)
%   3  ProjRun_2151_2300_noSEF  transient 2151-2300, no SEF
%                               (loadonly=0 submit, =1 gather)
%   4  VAFContinuityCheck_noSEF noSEF junction check + comparison vs SEF run;
%                               saves combined figure to
%                               postprocessed_data/figures/CESM2-WACCM/ssp585/
%   5  SMBComparison            total SMB timeline SEF vs noSEF (2015-2299);
%                               2-panel figure: absolute + difference (SEF-noSEF);
%                               saves to postprocessed_data/figures/CESM2-WACCM/ssp585/
%   6  WriteISMIP6_NetCDF_noSEF grids combined noSEF solution and writes
%                               ISMIP6 NetCDFs to
%                               postprocessed_data/CESM2-WACCM/proj_noSEF/
%                               (EXP='C007_noSEF')

    % ---- swappable scenario label ----------------------------------------
    SCENARIO   = 'ssp585';   % ssp126 | ssp370 | ssp534-over | ssp585
    CMIP_MODEL = 'CESM2-WACCM';
    % -----------------------------------------------------------------------

    if ~exist('loadonly','var'), loadonly = 0; end
    addpath('./../../../init/scripts');

    org = organizer('repository', './Models/', ...
                    'prefix', ['AIS_ISMIP7_Proj_' SCENARIO '_'], ...
                    'steps', steps, 'color', '34;47;2');
    clear steps;

    % ------------------------------------------------------------------ paths
    proj_root  = './../../../';
    init_dir   = [proj_root 'init/'];

    inputmodel_2015 = [proj_root 'hist_runs/CESM2-WACCM/Models/' ...
                       'AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat'];

    preproc_ocean      = [proj_root 'preprocessed_data/Ocean/'];
    preproc_proj_ocean = [proj_root 'preprocessed_data/Ocean/Proj/' CMIP_MODEL '/' SCENARIO '/'];
    preproc_proj_atmo  = [proj_root 'preprocessed_data/Atmosphere/Proj/' CMIP_MODEL '/' SCENARIO '/'];

    % ------------------------------------------------------------------ time
    start_year  = 2015;
    mid_year    = 2150;
    end_year    = 2299;

    % ------------------------------------------------------------------ ISMIP6 outputs
    ismip6_outputs = { ...
        'default', ...
        'Thickness', 'Surface', 'Base', 'Bed', ...
        'MaskOceanLevelset', 'MaskIceLevelset', ...
        'SmbMassBalance', ...
        'BasalforcingsFloatingiceMeltingRate', 'BasalforcingsGroundediceMeltingRate', ...
        'Vel', 'Vx', 'Vy', ...
        'GroundinglineMassFlux', 'IcefrontMassFluxLevelset', ...
        'IceVolumeScaled', 'IceVolumeAboveFloatationScaled', 'GroundedAreaScaled', 'FloatingAreaScaled', ...
        'TotalSmbScaled', 'TotalGroundedBmbScaled', 'TotalFloatingBmbScaled'};


    % ================================================================= Step 1
    if perform(org, 'ProjRun_2015_2150_noSEF') % {{{
        % As main script step 4 but without surface-elevation feedback:
        % SMB is assigned as a plain time-varying field (SMBforcing class, no gradient
        % terms). Ocean TF, levelset, and timestepping are identical to step 4.

        md = loadmodel(inputmodel_2015);
        m  = ((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1) ...
              ./ (1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor = (1./m).^2;

        md.inversion.iscontrol       = 0;
        md.transient.isthermal       = 0;
        md.transient.isgroundingline = 1;
        md.transient.ismasstransport = 1;
        md.transient.isstressbalance = 1;
        md.masstransport.spcthickness   = NaN*ones(md.mesh.numberofvertices, 1);
        md.outputdefinition.definitions = {};
        md.timestepping.interp_forcing  = 1;

        md.timestepping.start_time = start_year;
        md.timestepping.final_time = mid_year + 1;
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;

        md.transient.requested_outputs = ismip6_outputs;

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean ---
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_proj_ocean 'CESM_WACCM_TF_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        load([preproc_ocean 'gamma0_local.mat']);

        for di = 1:numel(tf_proj)
            c = tf_proj{di}; t_c = c(end,:);
            tf_proj{di} = c(:, t_c <= mid_year + 1);
        end

        unique_basinid = unique(basinid);
        delta_t        = zeros(1, length(unique_basinid));

        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique_basinid);
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = tf_proj;
        md.basalforcings.islocal    = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_local;

        % --- SMB (no SEF: plain SMB(), corrected forcing, no gradient terms) ---
        % smb_forcing is in mm w.e. yr-1 (same units as SMBgradients.smbref).
        % SMBforcing.mass_balance expects m ice yr-1 -- ISSM does NOT apply the
        % /1000*rho_w/rho_i conversion that SMBgradients does internally, so we
        % must do it here.  Only the data rows are converted; the time row (last
        % row) is left in years.
        load([preproc_proj_atmo 'CESM_WACCM_SMB_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        keep4         = smb_forcing(end,:) <= mid_year + 1;
        smb_forcing   = smb_forcing(:, keep4);
        rho_w = 1000.0;
        rho_i = md.materials.rho_ice;
        smb_mice      = smb_forcing;
        smb_mice(1:end-1,:) = smb_forcing(1:end-1,:) / 1000.0 * (rho_w / rho_i);

        md.smb              = SMBforcing();
        md.smb.mass_balance = smb_mice;

        % --- Calving front ---
        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 1;
        load([preproc_proj_ocean 'CESM_WACCM_levelset_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        keep4            = proj_spclevelset(end,:) <= mid_year + 1;
        proj_spclevelset = proj_spclevelset(:, keep4);
        md.levelset.spclevelset = proj_spclevelset;

        md.miscellaneous.name = ['ProjRun_' CMIP_MODEL '_' SCENARIO '_' ...
                                 num2str(start_year) '_' num2str(mid_year) '_noSEF'];
        clustername = 'gadi';
        md.cluster  = set_cluster(clustername);
        md.settings.waitonlock = 0;
        md.verbose = verbose('solution', true, 'module', true, 'convergence', true);
        md = solve(md, 'tr', 'runtimename', false, 'loadonly', loadonly);
        if loadonly
            savemodel(org, md);
        end
    end % }}}

    % ================================================================= Step 2
    if perform(org, 'AIS_state_2151_noSEF') % {{{
        md = loadmodel(org, 'ProjRun_2015_2150_noSEF');
        md_in = md;
        md.geometry.thickness  = md_in.results.TransientSolution(end).Thickness;
        md.geometry.surface    = md_in.results.TransientSolution(end).Surface;
        md.geometry.base       = md_in.results.TransientSolution(end).Base;
        md.mask.ocean_levelset = md_in.results.TransientSolution(end).MaskOceanLevelset;
        md.mask.ice_levelset   = md_in.results.TransientSolution(end).MaskIceLevelset;
        md.results.TransientSolution = [];
        savemodel(org, md);
    end % }}}

    % ================================================================= Step 3
    if perform(org, 'ProjRun_2151_2300_noSEF') % {{{
        % As main script step 6 but without surface-elevation feedback.

        md = loadmodel(org, 'AIS_state_2151_noSEF');
        m  = ((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1) ...
              ./ (1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor = (1./m).^2;

        md.inversion.iscontrol       = 0;
        md.transient.isthermal       = 0;
        md.transient.isgroundingline = 1;
        md.transient.ismasstransport = 1;
        md.transient.isstressbalance = 1;
        md.masstransport.spcthickness   = NaN*ones(md.mesh.numberofvertices, 1);
        md.outputdefinition.definitions = {};
        md.timestepping.interp_forcing  = 1;

        md.timestepping.start_time = mid_year + 1;
        md.timestepping.final_time = end_year + 1;
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;

        md.transient.requested_outputs = ismip6_outputs;

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean ---
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_proj_ocean 'CESM_WACCM_TF_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        load([preproc_ocean 'gamma0_local.mat']);

        for di = 1:numel(tf_proj)
            c = tf_proj{di}; t_c = c(end,:);
            tf_proj{di} = c(:, t_c >= mid_year + 1);
        end

        unique_basinid = unique(basinid);
        delta_t        = zeros(1, length(unique_basinid));

        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique_basinid);
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = tf_proj;
        md.basalforcings.islocal    = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_local;

        % --- SMB (no SEF) ---
        load([preproc_proj_atmo 'CESM_WACCM_SMB_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        keep6         = smb_forcing(end,:) >= mid_year + 1;
        smb_forcing   = smb_forcing(:, keep6);
        rho_w = 1000.0;
        rho_i = md.materials.rho_ice;
        smb_mice      = smb_forcing;
        smb_mice(1:end-1,:) = smb_forcing(1:end-1,:) / 1000.0 * (rho_w / rho_i);

        md.smb              = SMBforcing();
        md.smb.mass_balance = smb_mice;

        % --- Calving front ---
        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 1;
        load([preproc_proj_ocean 'CESM_WACCM_levelset_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        keep6            = proj_spclevelset(end,:) >= mid_year + 1;
        proj_spclevelset = proj_spclevelset(:, keep6);
        md.levelset.spclevelset = proj_spclevelset;

        md.miscellaneous.name = ['ProjRun_' CMIP_MODEL '_' SCENARIO '_' ...
                                 num2str(mid_year+1) '_' num2str(end_year) '_noSEF'];
        clustername = 'gadi';
        md.cluster  = set_cluster(clustername);
        md.settings.waitonlock = 0;
        md.verbose = verbose('solution', true, 'module', true, 'convergence', true);
        md = solve(md, 'tr', 'runtimename', false, 'loadonly', loadonly);
        if loadonly
            savemodel(org, md);
        end
    end % }}}

    % ================================================================= Step 4
    if perform(org, 'VAFContinuityCheck_noSEF') % {{{
        % noSEF junction continuity check + comparison with SEF run.
        % Loads all four segment models (noSEF and SEF) and plots both full
        % 2015-2300 time series on one axes.

        md1_nsef = loadmodel(org, 'ProjRun_2015_2150_noSEF');
        md2_nsef = loadmodel(org, 'ProjRun_2151_2300_noSEF');
        md1_sef  = loadmodel(org, 'ProjRun_2015_2150');
        md2_sef  = loadmodel(org, 'ProjRun_2151_2300');

        time1_nsef = [md1_nsef.results.TransientSolution.time];
        time2_nsef = [md2_nsef.results.TransientSolution.time];
        vaf1_nsef  = [md1_nsef.results.TransientSolution.IceVolumeAboveFloatationScaled];
        vaf2_nsef  = [md2_nsef.results.TransientSolution.IceVolumeAboveFloatationScaled];

        time1_sef  = [md1_sef.results.TransientSolution.time];
        time2_sef  = [md2_sef.results.TransientSolution.time];
        vaf1_sef   = [md1_sef.results.TransientSolution.IceVolumeAboveFloatationScaled];
        vaf2_sef   = [md2_sef.results.TransientSolution.IceVolumeAboveFloatationScaled];

        rho_ice = 917; rho_sw = 1028; A_ocean = 3.625e14;
        vaf_ref = vaf1_nsef(1);   % 2015 baseline; both runs start from the same state

        sle1_nsef = -(vaf1_nsef - vaf_ref) * rho_ice / rho_sw / A_ocean;
        sle2_nsef = -(vaf2_nsef - vaf_ref) * rho_ice / rho_sw / A_ocean;
        sle1_sef  = -(vaf1_sef  - vaf_ref) * rho_ice / rho_sw / A_ocean;
        sle2_sef  = -(vaf2_sef  - vaf_ref) * rho_ice / rho_sw / A_ocean;

        gap_m = sle2_nsef(1) - sle1_nsef(end);
        fprintf('VAF continuity check (noSEF):\n');
        fprintf('  Segment 1 end   (t=%.1f): %.4f m SLE\n', time1_nsef(end), sle1_nsef(end));
        fprintf('  Segment 2 start (t=%.1f): %.4f m SLE\n', time2_nsef(1),   sle2_nsef(1));
        fprintf('  Junction gap: %.6f m SLE\n', gap_m);
        if abs(gap_m) > 0.001
            warning('VAF junction gap > 1 mm SLE (%.6f m) -- check AIS_state_2151_noSEF restart.', gap_m);
        end

        fig = figure('visible', 'off');
        c_sef  = [0.00 0.45 0.70];   % blue  — SEF
        c_nsef = [0.80 0.40 0.00];   % orange — noSEF
        plot([time1_sef  time2_sef],  [sle1_sef  sle2_sef],  '-',  'Color', c_sef,  'LineWidth', 1.5); hold on;
        plot([time1_nsef time2_nsef], [sle1_nsef sle2_nsef], '--', 'Color', c_nsef, 'LineWidth', 1.5);
        xline(mid_year + 1, 'k:', 'LineWidth', 1, 'Label', '2151 restart', ...
              'LabelVerticalAlignment', 'bottom');
        xlabel('Year');
        ylabel('Sea level contribution (m SLE)');
        title(sprintf('%s %s — SEF vs noSEF', CMIP_MODEL, SCENARIO), ...
              'Interpreter', 'none');
        legend('SEF (with lapse-rate correction)', 'noSEF (plain SMB)', ...
               'Location', 'northwest');
        grid on;

        figdir = [proj_root 'postprocessed_data/figures/' CMIP_MODEL '/' SCENARIO '/'];
        if ~exist(figdir, 'dir'), mkdir(figdir); end
        figname = fullfile(figdir, sprintf('VAF_SEF_vs_noSEF_%s_%s.png', CMIP_MODEL, SCENARIO));
        saveas(fig, figname);
        close(fig);
        fprintf('Saved: %s\n', figname);
    end % }}}

    % ================================================================= Step 5
    if perform(org, 'SMBComparison') % {{{
        % Compare total SMB timelines (SEF vs noSEF) over 2015-2299.
        % TotalSmbScaled from ISSM is the mesh-integral of SMB, in Gt yr-1.
        % Plots: (top) absolute timelines; (bottom) difference SEF minus noSEF.

        md1_nsef = loadmodel(org, 'ProjRun_2015_2150_noSEF');
        md2_nsef = loadmodel(org, 'ProjRun_2151_2300_noSEF');
        md1_sef  = loadmodel(org, 'ProjRun_2015_2150');
        md2_sef  = loadmodel(org, 'ProjRun_2151_2300');

        % Concatenate segments
        ts_nsef = [md1_nsef.results.TransientSolution, md2_nsef.results.TransientSolution];
        ts_sef  = [md1_sef.results.TransientSolution,  md2_sef.results.TransientSolution];
        clear md1_nsef md2_nsef md1_sef md2_sef

        % Keep only annual outputs (same filter as write functions)
        t_nsef_raw = [ts_nsef.time];
        t_sef_raw  = [ts_sef.time];
        keep_nsef  = abs(t_nsef_raw - round(t_nsef_raw)) < 0.05;
        keep_sef   = abs(t_sef_raw  - round(t_sef_raw))  < 0.05;

        t_nsef   = t_nsef_raw(keep_nsef) - 1;   % nominal years (2015…2299)
        t_sef    = t_sef_raw(keep_sef)   - 1;
        smb_nsef = [ts_nsef(keep_nsef).TotalSmbScaled];  % Gt yr-1
        smb_sef  = [ts_sef(keep_sef).TotalSmbScaled];

        % Interpolate noSEF onto SEF time axis for difference
        smb_nsef_i = interp1(t_nsef, smb_nsef, t_sef, 'linear', NaN);
        delta_smb  = smb_sef - smb_nsef_i;   % positive = SEF has more SMB

        c_sef  = [0.00 0.45 0.70];   % blue
        c_nsef = [0.80 0.40 0.00];   % orange

        fig = figure('visible', 'off', 'Position', [100 100 900 600]);

        subplot(2, 1, 1);
        plot(t_sef,  smb_sef,  '-',  'Color', c_sef,  'LineWidth', 1.5); hold on;
        plot(t_nsef, smb_nsef, '--', 'Color', c_nsef, 'LineWidth', 1.5);
        xline(mid_year, 'k:', 'LineWidth', 1);
        ylabel('Total SMB (Gt yr^{-1})');
        title(sprintf('%s %s — Total SMB: SEF vs noSEF', CMIP_MODEL, SCENARIO), ...
              'Interpreter', 'none');
        legend('SEF', 'noSEF', 'Location', 'best');
        grid on; box on;

        subplot(2, 1, 2);
        plot(t_sef, delta_smb, 'k-', 'LineWidth', 1.5);
        yline(0, 'k--', 'LineWidth', 0.8);
        xline(mid_year, 'k:', 'LineWidth', 1);
        xlabel('Year');
        ylabel('\Delta SMB  SEF \minus noSEF (Gt yr^{-1})');
        title('Difference (SEF minus noSEF)');
        grid on; box on;

        figdir = [proj_root 'postprocessed_data/figures/' CMIP_MODEL '/' SCENARIO '/'];
        if ~exist(figdir, 'dir'), mkdir(figdir); end
        figname = fullfile(figdir, sprintf('SMB_SEF_vs_noSEF_%s_%s.png', CMIP_MODEL, SCENARIO));
        saveas(fig, figname);
        close(fig);
        fprintf('Saved: %s\n', figname);
    end % }}}

    % ================================================================= Step 6
    if perform(org, 'WriteISMIP6_NetCDF_noSEF') % {{{
        % Combines the two noSEF TransientSolution arrays (steps 1 + 3) and
        % grids them onto the standard ISMIP7 761×761 8 km AIS grid.
        addpath('./../../functions');

        md1 = loadmodel(org, 'ProjRun_2015_2150_noSEF');
        md2 = loadmodel(org, 'ProjRun_2151_2300_noSEF');

        md = md1;
        md.results.TransientSolution = [md1.results.TransientSolution, ...
                                        md2.results.TransientSolution];
        clear md1 md2;

        outdir = [proj_root 'postprocessed_data/CESM2-WACCM/proj_noSEF/'];
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        meta                    = struct();
        meta.experiment_id      = SCENARIO;        % 'ssp585'
        meta.set_counter        = 'C007';
        meta.time_range         = '2015-2299';
        meta.ESM_id             = CMIP_MODEL;      % 'CESM2-WACCM'
        meta.forcing_member_id  = 'f001';
        meta.ISM_member_id      = 'm001';

        [cfflux_tot, glflux_tot] = write_ismip7_2d_projection(md, outdir, meta);
        write_ismip7_scalar_projection(md, outdir, meta, cfflux_tot, glflux_tot);
    end % }}}

end

% NetCDF writing delegated to proj_runs/functions/write_ismip7_2d_projection.m
% and proj_runs/functions/write_ismip7_scalar_projection.m.

