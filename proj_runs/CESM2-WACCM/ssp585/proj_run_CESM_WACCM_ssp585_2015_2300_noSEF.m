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
%   5  WriteISMIP6_NetCDF_noSEF grids combined noSEF solution and writes
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
    if perform(org, 'WriteISMIP6_NetCDF_noSEF') % {{{
        % Combines the two noSEF TransientSolution arrays (steps 1 + 3) and
        % grids them onto the ISMIP6 8 km AIS grid. EXP tag is 'C007_noSEF'.

        md1 = loadmodel(org, 'ProjRun_2015_2150_noSEF');
        md2 = loadmodel(org, 'ProjRun_2151_2300_noSEF');

        md = md1;
        md.results.TransientSolution = [md1.results.TransientSolution, ...
                                        md2.results.TransientSolution];
        clear md1 md2;

        outdir = [proj_root 'postprocessed_data/CESM2-WACCM/proj_noSEF/'];
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        write_proj_netcdf_ismip6(md, outdir, 'C007_noSEF', start_year);
    end % }}}

end

% ------------------------------------------------------------------ helpers
% (Identical to proj_run_CESM_WACCM_ssp585_2015_2300.m -- keep in sync.)

function write_proj_netcdf_ismip6(md, outdir, EXP, start_year)
% Write all ISMIP6 Appendix-2 NetCDF files for a projection run.
    IS        = 'AIS';
    GROUP     = 'MONASH';
    MODELNAME = 'ISSM';

    rhoi = md.materials.rho_ice;
    yts  = md.constants.yts;

    % Filter to annual outputs; ISSM emits one sub-annual safety step at the
    % start of each segment (e.g. t≈2015.1, t≈2151.1).
    t_raw = [md.results.TransientSolution.time];
    keep  = abs(t_raw - round(t_raw)) < 0.05;
    md.results.TransientSolution = md.results.TransientSolution(keep);
    nT    = sum(keep);

    % ISSM t=Y represents the state after running from start_year to Y, i.e.
    % the end-of-year (Y-1) state.  Subtract 1 so labels read 2015–2299.
    time_yr = round(t_raw(keep)) - 1;   % [2015, 2016, ..., 2299]

    time_units = 'days since 1995-01-01';
    % ST: snapshot at end of year (no time_bnds)
    time_st      = (time_yr - 1994) * 365;       % [7665, 8030, ..., 111325]
    % FL: midpoint of year, with time_bnds
    lb_fl        = (time_yr - 1995) * 365;        % [7300, 7665, ..., 110960]
    ub_fl        = lb_fl + 365;
    time_fl      = lb_fl + 182.5;                 % [7482.5, 7847.5, ...]
    time_bnds_fl = [lb_fl; ub_fl]';

    missing_value = single(1e20);

    [~, xGrid, yGrid] = gridData(md, md.results.TransientSolution(1).Thickness);
    nx = length(xGrid); ny = length(yGrid);

    % ---- 2D snapshot/rate fields ----
    vars2d = { ...
        'lithk',       'Thickness',                            1,        'land_ice_thickness',                           'm',         'ST'; ...
        'orog',        'Surface',                              1,        'surface_altitude',                             'm',         'ST'; ...
        'base',        'Base',                                 1,        'base_altitude',                                'm',         'ST'; ...
        'topg',        'Bed',                                  1,        'bedrock_altitude',                             'm',         'ST'; ...
        'xvelsurf',    'Vx',                                   1/yts,    'land_ice_surface_x_velocity',                  'm s-1',     'ST'; ...
        'yvelsurf',    'Vy',                                   1/yts,    'land_ice_surface_y_velocity',                  'm s-1',     'ST'; ...
        'acabf',       'SmbMassBalance',                       rhoi/yts, 'land_ice_surface_specific_mass_balance_flux',  'kg m-2 s-1','FL'; ...
        'libmassbfgr', 'BasalforcingsGroundediceMeltingRate', -rhoi/yts, 'land_ice_basal_specific_mass_balance_flux',   'kg m-2 s-1','FL'; ...
        'libmassbffl', 'BasalforcingsFloatingiceMeltingRate', -rhoi/yts, 'land_ice_basal_specific_mass_balance_flux',   'kg m-2 s-1','FL'; ...
    };
    for k = 1:size(vars2d, 1)
        varname = vars2d{k,1}; issmfield = vars2d{k,2}; scale = vars2d{k,3};
        long_name = vars2d{k,4}; units = vars2d{k,5}; vtype = vars2d{k,6};
        if strcmp(vtype, 'ST'), t_vec = time_st; t_bnds = [];
        else,                   t_vec = time_fl; t_bnds = time_bnds_fl; end
        data3d = zeros(ny, nx, nT, 'single');
        for t = 1:nT
            data3d(:,:,t) = single(gridData(md, md.results.TransientSolution(t).(issmfield)) * scale);
        end
        write_ismip6_2d(outdir, varname, IS, GROUP, MODELNAME, EXP, data3d, ...
            xGrid, yGrid, t_vec, t_bnds, time_units, long_name, units, missing_value);
    end

    % ---- Static fields ----
    hfgeoubed2d = single(gridData(md, md.basalforcings.geothermalflux));
    write_ismip6_2d(outdir, 'hfgeoubed', IS, GROUP, MODELNAME, EXP, repmat(hfgeoubed2d,1,1,nT), ...
        xGrid, yGrid, time_fl, time_bnds_fl, time_units, 'upward_geothermal_heat_flux_at_ground_level', 'W m-2', missing_value);
    temp2d = single(gridData(md, md.initialization.temperature));
    write_ismip6_2d(outdir, 'litemptop', IS, GROUP, MODELNAME, EXP, repmat(temp2d,1,1,nT), ...
        xGrid, yGrid, time_st, [], time_units, 'temperature_at_top_of_ice_sheet_model', 'K', missing_value);

    % ---- dlithkdt (FL) ----
    H_grid = zeros(ny, nx, nT);
    for t = 1:nT, H_grid(:,:,t) = gridData(md, md.results.TransientSolution(t).Thickness); end
    dlithkdt3d = zeros(ny, nx, nT, 'single');
    for t = 1:nT
        if t < nT
            dlithkdt3d(:,:,t) = single((H_grid(:,:,t+1) - H_grid(:,:,t)) / yts);
        else
            dlithkdt3d(:,:,t) = single((H_grid(:,:,t) - H_grid(:,:,t-1)) / yts);
        end
    end
    write_ismip6_2d(outdir, 'dlithkdt', IS, GROUP, MODELNAME, EXP, dlithkdt3d, ...
        xGrid, yGrid, time_fl, time_bnds_fl, time_units, 'tendency_of_land_ice_thickness', 'm s-1', missing_value);

    % ---- Masks (ST) ----
    sftgif3d = zeros(ny, nx, nT, 'single');
    sftgrf3d = zeros(ny, nx, nT, 'single');
    sftflf3d = zeros(ny, nx, nT, 'single');
    for t = 1:nT
        ice_ls   = md.results.TransientSolution(t).MaskIceLevelset;
        ocean_ls = md.results.TransientSolution(t).MaskOceanLevelset;
        sftgif3d(:,:,t) = single(gridData(md, double(ice_ls < 0)));
        sftgrf3d(:,:,t) = single(gridData(md, double(ice_ls < 0 & ocean_ls >= 0)));
        sftflf3d(:,:,t) = single(gridData(md, double(ice_ls < 0 & ocean_ls < 0)));
    end
    write_ismip6_2d(outdir, 'sftgif', IS, GROUP, MODELNAME, EXP, sftgif3d, ...
        xGrid, yGrid, time_st, [], time_units, 'land_ice_area_fraction', '1', missing_value);
    write_ismip6_2d(outdir, 'sftgrf', IS, GROUP, MODELNAME, EXP, sftgrf3d, ...
        xGrid, yGrid, time_st, [], time_units, 'grounded_ice_sheet_area_fraction', '1', missing_value);
    write_ismip6_2d(outdir, 'sftflf', IS, GROUP, MODELNAME, EXP, sftflf3d, ...
        xGrid, yGrid, time_st, [], time_units, 'floating_ice_shelf_area_fraction', '1', missing_value);

    % ---- licalvf / ligroundf (FL) ----
    cellsize   = xGrid(2) - xGrid(1);
    glflux2d   = zeros(ny, nx, nT, 'single'); glflux_tot = zeros(1, nT);
    cfflux2d   = zeros(ny, nx, nT, 'single'); cfflux_tot = zeros(1, nT);
    for t = 1:nT
        sol = md.results.TransientSolution(t);
        [gl2d, gltot] = flux_along_contour_2d(md, sol, sol.MaskOceanLevelset, true,  xGrid, yGrid, cellsize);
        [cf2d, cftot] = flux_along_contour_2d(md, sol, sol.MaskIceLevelset,   false, xGrid, yGrid, cellsize);
        glflux2d(:,:,t) = gl2d; glflux_tot(t) = gltot;
        cfflux2d(:,:,t) = cf2d; cfflux_tot(t) = cftot;
    end
    write_ismip6_2d(outdir, 'ligroundf', IS, GROUP, MODELNAME, EXP, glflux2d, ...
        xGrid, yGrid, time_fl, time_bnds_fl, time_units, 'land_ice_specific_mass_flux_at_grounding_line', 'kg m-2 s-1', missing_value);
    write_ismip6_2d(outdir, 'licalvf', IS, GROUP, MODELNAME, EXP, cfflux2d, ...
        xGrid, yGrid, time_fl, time_bnds_fl, time_units, 'land_ice_specific_mass_flux_due_to_calving', 'kg m-2 s-1', missing_value);

    % ---- Scalars ----
    scalars = { ...
        'lim',      'IceVolumeScaled',               rhoi,     'land_ice_mass',                                         'kg',    'ST'; ...
        'limnsw',   'IceVolumeAboveFloatationScaled', rhoi,     'land_ice_mass_not_displacing_sea_water',                'kg',    'ST'; ...
        'iareagr',  'GroundedAreaScaled',             1,        'grounded_land_ice_area',                                'm2',    'ST'; ...
        'iareafl',  'FloatingAreaScaled',             1,        'floating_ice_shelf_area',                               'm2',    'ST'; ...
        'tendacabf','TotalSmbScaled',                 1e12/yts, 'tendency_of_land_ice_mass_due_to_surface_mass_balance', 'kg s-1','FL'; ...
    };
    for k = 1:size(scalars, 1)
        varname = scalars{k,1}; issmfield = scalars{k,2}; scale = scalars{k,3};
        long_name = scalars{k,4}; units = scalars{k,5}; vtype = scalars{k,6};
        if strcmp(vtype, 'ST'), t_vec = time_st; t_bnds = [];
        else,                   t_vec = time_fl; t_bnds = time_bnds_fl; end
        data1d = zeros(1, nT);
        for t = 1:nT, data1d(t) = md.results.TransientSolution(t).(issmfield) * scale; end
        write_ismip6_scalar(outdir, varname, IS, GROUP, MODELNAME, EXP, data1d, ...
            t_vec, t_bnds, time_units, long_name, units, missing_value);
    end

    tgbmb = zeros(1, nT); tfbmb = zeros(1, nT);
    for t = 1:nT
        tgbmb(t) = md.results.TransientSolution(t).TotalGroundedBmbScaled * 1e12/yts;
        tfbmb(t) = md.results.TransientSolution(t).TotalFloatingBmbScaled * 1e12/yts;
    end
    write_ismip6_scalar(outdir, 'tendlibmassbf', IS, GROUP, MODELNAME, EXP, -(tgbmb + tfbmb), ...
        time_fl, time_bnds_fl, time_units, 'tendency_of_land_ice_mass_due_to_basal_mass_balance', 'kg s-1', missing_value);
    write_ismip6_scalar(outdir, 'tendlibmassbffl', IS, GROUP, MODELNAME, EXP, -tfbmb, ...
        time_fl, time_bnds_fl, time_units, 'tendency_of_land_ice_mass_due_to_basal_mass_balance_at_ice_shelf_base', 'kg s-1', missing_value);

    cf_native = zeros(1, nT); gl_native = zeros(1, nT);
    for t = 1:nT
        cf_native(t) = md.results.TransientSolution(t).IcefrontMassFluxLevelset * 1e12/yts;
        gl_native(t) = md.results.TransientSolution(t).GroundinglineMassFlux    * 1e12/yts;
    end
    if sign(sum(cf_native)) ~= sign(sum(cfflux_tot)) && any(cfflux_tot)
        warning('IcefrontMassFluxLevelset sign disagrees with cfflux_tot -- flipping to match.');
        cf_native = -cf_native;
    end
    if sign(sum(gl_native)) ~= sign(sum(glflux_tot)) && any(glflux_tot)
        warning('GroundinglineMassFlux sign disagrees with glflux_tot -- flipping to match.');
        gl_native = -gl_native;
    end
    write_ismip6_scalar(outdir, 'tendlicalvf', IS, GROUP, MODELNAME, EXP, cf_native, ...
        time_fl, time_bnds_fl, time_units, 'tendency_of_land_ice_mass_due_to_calving', 'kg s-1', missing_value);
    write_ismip6_scalar(outdir, 'tendligroundf', IS, GROUP, MODELNAME, EXP, gl_native, ...
        time_fl, time_bnds_fl, time_units, 'tendency_of_land_ice_mass_due_to_flux_at_grounding_line', 'kg s-1', missing_value);
end

function fname = ismip6_filename(outdir, varname, IS, GROUP, MODELNAME, EXP)
    fname = fullfile(outdir, sprintf('%s_%s_%s_%s_%s.nc', ...
                     varname, IS, GROUP, MODELNAME, EXP));
end

function write_ismip6_2d(outdir, varname, IS, GROUP, MODELNAME, EXP, data3d, ...
        xGrid, yGrid, time_vec, time_bnds, time_units, long_name, units, missing_value)
    fname = ismip6_filename(outdir, varname, IS, GROUP, MODELNAME, EXP);
    if exist(fname, 'file'), delete(fname); end
    nx = length(xGrid); ny = length(yGrid);

    data3d(isnan(data3d)) = missing_value;

    nccreate(fname, 'x', 'Dimensions', {'x', nx}, 'Datatype', 'single');
    ncwrite(fname, 'x', single(xGrid));
    ncwriteatt(fname, 'x', 'units', 'm');
    ncwriteatt(fname, 'x', 'standard_name', 'projection_x_coordinate');

    nccreate(fname, 'y', 'Dimensions', {'y', ny}, 'Datatype', 'single');
    ncwrite(fname, 'y', single(yGrid));
    ncwriteatt(fname, 'y', 'units', 'm');
    ncwriteatt(fname, 'y', 'standard_name', 'projection_y_coordinate');

    nccreate(fname, 'time', 'Dimensions', {'time', Inf}, 'Datatype', 'double');
    ncwrite(fname, 'time', time_vec);
    ncwriteatt(fname, 'time', 'units', time_units);
    ncwriteatt(fname, 'time', 'long_name', 'time');
    ncwriteatt(fname, 'time', 'standard_name', 'time');
    ncwriteatt(fname, 'time', 'axis', 'T');
    ncwriteatt(fname, 'time', 'calendar', '365_day');
    if ~isempty(time_bnds)
        ncwriteatt(fname, 'time', 'bounds', 'time_bnds');
        nccreate(fname, 'time_bnds', 'Dimensions', {'bnds', 2, 'time', Inf}, 'Datatype', 'double');
        ncwrite(fname, 'time_bnds', time_bnds');
    end

    nccreate(fname, varname, 'Dimensions', {'x', nx, 'y', ny, 'time', Inf}, ...
             'Datatype', 'single', 'FillValue', missing_value);
    ncwrite(fname, varname, permute(data3d, [2 1 3]));
    ncwriteatt(fname, varname, 'long_name', long_name);
    ncwriteatt(fname, varname, 'units', units);
    ncwriteatt(fname, varname, 'missing_value', missing_value);

    fprintf('Saved: %s\n', fname);
end

function write_ismip6_scalar(outdir, varname, IS, GROUP, MODELNAME, EXP, data1d, ...
        time_vec, time_bnds, time_units, long_name, units, missing_value)
    fname = ismip6_filename(outdir, varname, IS, GROUP, MODELNAME, EXP);
    if exist(fname, 'file'), delete(fname); end

    nccreate(fname, 'time', 'Dimensions', {'time', Inf}, 'Datatype', 'double');
    ncwrite(fname, 'time', time_vec);
    ncwriteatt(fname, 'time', 'units', time_units);
    ncwriteatt(fname, 'time', 'long_name', 'time');
    ncwriteatt(fname, 'time', 'standard_name', 'time');
    ncwriteatt(fname, 'time', 'axis', 'T');
    ncwriteatt(fname, 'time', 'calendar', '365_day');
    if ~isempty(time_bnds)
        ncwriteatt(fname, 'time', 'bounds', 'time_bnds');
        nccreate(fname, 'time_bnds', 'Dimensions', {'bnds', 2, 'time', Inf}, 'Datatype', 'double');
        ncwrite(fname, 'time_bnds', time_bnds');
    end

    nccreate(fname, varname, 'Dimensions', {'time', Inf}, 'Datatype', 'single');
    ncwrite(fname, varname, single(data1d));
    ncwriteatt(fname, varname, 'long_name', long_name);
    ncwriteatt(fname, varname, 'units', units);
    ncwriteatt(fname, varname, 'missing_value', missing_value);

    fprintf('Saved: %s\n', fname);
end

function flux_total_kgs = gl_flux_native_mesh(md, sol, apply_correction)
    if nargin < 3, apply_correction = true; end
    flux_total_kgs = 0;
    elems = md.mesh.elements; x = md.mesh.x; y = md.mesh.y;

    gl1 = isoline(md, sol.MaskOceanLevelset, 'value', 0, 'output', 'matrix');
    if isempty(gl1) || size(gl1,1) < 2, return; end

    valid = ~any(isnan(gl1), 2);
    idx = find(valid(1:end-1) & valid(2:end));
    if isempty(idx), return; end

    vx = nan(size(gl1,1),1); vy = nan(size(gl1,1),1); h = nan(size(gl1,1),1);
    vx(valid) = InterpFromMesh2d(elems, x, y, sol.Vx,        gl1(valid,1), gl1(valid,2));
    vy(valid) = InterpFromMesh2d(elems, x, y, sol.Vy,        gl1(valid,1), gl1(valid,2));
    h(valid)  = InterpFromMesh2d(elems, x, y, sol.Thickness, gl1(valid,1), gl1(valid,2));

    x1 = gl1(idx,1);   y1 = gl1(idx,2);
    x2 = gl1(idx+1,1); y2 = gl1(idx+1,2);
    dx = x2 - x1; dy = y2 - y1; L = hypot(dx, dy);
    good = (L > 0) & ~isnan(L);
    if ~any(good), return; end
    idx = idx(good); x1 = x1(good); y1 = y1(good); x2 = x2(good); y2 = y2(good); L = L(good);
    dx = dx(good); dy = dy(good);

    Nx = -dy ./ L; Ny = dx ./ L;
    Vx1 = vx(idx);   Vy1 = vy(idx);   H1 = h(idx);
    Vx2 = vx(idx+1); Vy2 = vy(idx+1); H2 = h(idx+1);
    good2 = ~any(isnan([Vx1 Vy1 H1 Vx2 Vy2 H2]), 2);
    if ~any(good2), return; end
    x1 = x1(good2); y1 = y1(good2); x2 = x2(good2); y2 = y2(good2); L = L(good2);
    Nx = Nx(good2); Ny = Ny(good2);
    Vx1 = Vx1(good2); Vy1 = Vy1(good2); H1 = H1(good2);
    Vx2 = Vx2(good2); Vy2 = Vy2(good2); H2 = H2(good2);

    xm = 0.5*(x1+x2); ym = 0.5*(y1+y2);
    if apply_correction
        eps_n = 1000;
        phi_plus  = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm+eps_n*Nx, ym+eps_n*Ny);
        phi_minus = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm-eps_n*Nx, ym-eps_n*Ny);
        flip = (~isnan(phi_plus) & ~isnan(phi_minus) & (phi_plus > phi_minus));
        Nx(flip) = -Nx(flip); Ny(flip) = -Ny(flip);
    end

    Vxm = InterpFromMesh2d(elems, x, y, sol.Vx,        xm, ym);
    Vym = InterpFromMesh2d(elems, x, y, sol.Vy,        xm, ym);
    Hm  = InterpFromMesh2d(elems, x, y, sol.Thickness, xm, ym);
    good3 = ~any(isnan([Vxm Vym Hm]), 2);

    f1 = H1 .* (Vx1.*Nx + Vy1.*Ny);
    f2 = H2 .* (Vx2.*Nx + Vy2.*Ny);
    fm = Hm .* (Vxm.*Nx + Vym.*Ny);
    fm(~good3) = 0.5*(f1(~good3) + f2(~good3));

    secflux_m3_yr  = L/6 .* (f1 + 4*fm + f2);
    flux_total_kgs = sum(secflux_m3_yr, 'omitnan') * md.materials.rho_ice / md.constants.yts;
end

function [flux2d_kgm2s, flux_total_kgs] = flux_along_contour_2d(md, sol, contour_field, signed_normal, xGrid, yGrid, cellsize)
    nx = length(xGrid); ny = length(yGrid);
    flux2d_kgm2s   = NaN(ny, nx, 'single');
    flux_total_kgs = 0;

    elems = md.mesh.elements; x = md.mesh.x; y = md.mesh.y;

    gl1 = isoline(md, contour_field, 'value', 0, 'output', 'matrix');
    if isempty(gl1) || size(gl1,1) < 2, return; end

    valid = ~any(isnan(gl1), 2);
    idx = find(valid(1:end-1) & valid(2:end));
    if isempty(idx), return; end

    vx = nan(size(gl1,1),1); vy = nan(size(gl1,1),1); h = nan(size(gl1,1),1);
    vx(valid) = InterpFromMesh2d(elems, x, y, sol.Vx,        gl1(valid,1), gl1(valid,2));
    vy(valid) = InterpFromMesh2d(elems, x, y, sol.Vy,        gl1(valid,1), gl1(valid,2));
    h(valid)  = InterpFromMesh2d(elems, x, y, sol.Thickness, gl1(valid,1), gl1(valid,2));

    x1 = gl1(idx,1);   y1 = gl1(idx,2);
    x2 = gl1(idx+1,1); y2 = gl1(idx+1,2);
    dx = x2 - x1; dy = y2 - y1; L = hypot(dx, dy);
    good = (L > 0) & ~isnan(L);
    if ~any(good), return; end
    idx = idx(good); x1 = x1(good); y1 = y1(good); x2 = x2(good); y2 = y2(good); L = L(good);
    dx = dx(good); dy = dy(good);

    Nx = -dy ./ L; Ny = dx ./ L;
    Vx1 = vx(idx);   Vy1 = vy(idx);   H1 = h(idx);
    Vx2 = vx(idx+1); Vy2 = vy(idx+1); H2 = h(idx+1);
    good2 = ~any(isnan([Vx1 Vy1 H1 Vx2 Vy2 H2]), 2);
    if ~any(good2), return; end
    x1 = x1(good2); y1 = y1(good2); x2 = x2(good2); y2 = y2(good2); L = L(good2);
    Nx = Nx(good2); Ny = Ny(good2);
    Vx1 = Vx1(good2); Vy1 = Vy1(good2); H1 = H1(good2);
    Vx2 = Vx2(good2); Vy2 = Vy2(good2); H2 = H2(good2);

    xm = 0.5*(x1+x2); ym = 0.5*(y1+y2);
    if signed_normal
        eps_n = 1000;
        phi_plus  = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm+eps_n*Nx, ym+eps_n*Ny);
        phi_minus = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm-eps_n*Nx, ym-eps_n*Ny);
        flip = (~isnan(phi_plus) & ~isnan(phi_minus) & (phi_plus > phi_minus));
        Nx(flip) = -Nx(flip); Ny(flip) = -Ny(flip);
    end

    Vxm = InterpFromMesh2d(elems, x, y, sol.Vx,        xm, ym);
    Vym = InterpFromMesh2d(elems, x, y, sol.Vy,        xm, ym);
    Hm  = InterpFromMesh2d(elems, x, y, sol.Thickness, xm, ym);
    good3 = ~any(isnan([Vxm Vym Hm]), 2);

    f1 = H1 .* (Vx1.*Nx + Vy1.*Ny);
    f2 = H2 .* (Vx2.*Nx + Vy2.*Ny);
    fm = Hm .* (Vxm.*Nx + Vym.*Ny);
    fm(~good3) = 0.5*(f1(~good3) + f2(~good3));

    secflux_m3_yr = L/6 .* (f1 + 4*fm + f2);
    if ~signed_normal
        secflux_m3_yr = abs(secflux_m3_yr);
    end
    secflux_kgs = secflux_m3_yr * md.materials.rho_ice / md.constants.yts;

    ix = round((xm - xGrid(1)) / cellsize) + 1;
    iy = round((ym - yGrid(1)) / cellsize) + 1;
    inb = ix >= 1 & ix <= nx & iy >= 1 & iy <= ny;

    accum   = zeros(ny, nx, 'single');
    touched = false(ny, nx);
    for s = find(inb)'
        accum(iy(s), ix(s))   = accum(iy(s), ix(s)) + single(secflux_kgs(s));
        touched(iy(s), ix(s)) = true;
    end
    flux_total_kgs = sum(secflux_kgs(inb));
    accum = accum / single(cellsize^2);
    flux2d_kgm2s(touched) = accum(touched);
end
