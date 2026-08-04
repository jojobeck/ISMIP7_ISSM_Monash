function md = proj_run_CESM_WACCM_ctrl_2015_2300(steps, loadonly)
% Control run 2015-2300: CESM2-WACCM ssp585's real 2015 forcing (ocean TF,
% SMB, SMB gradient), held FIXED for the entire run -- no forcing evolves
% after 2015. Used as the ISMIP7 control-run baseline to isolate the ice
% sheet model's own drift from any climate-forcing signal, for comparison
% against the transient ssp585 run (../ssp585/proj_run_CESM_WACCM_ssp585_2015_2300.m).
%
% Calving front: NOT prescribed via a levelset time series here (that
% approach was problematic for this run). Instead md.transient.ismovingfront
% = 0, so the level-set equation is never solved and the ice front stays
% frozen at whatever it is in the AIS_state_2015 starting state for the
% entire run -- no ice-shelf collapse mask is read or used at all.
%
% Only a single (2015) forcing snapshot is built for each of TF/SMB --
% with interp_forcing=0 (step-function forcing, set throughout this
% pipeline), ISSM holds that one value constant for the entire 2015-2300
% run, so no other years need to be read or interpolated. The whole
% 286-year run is a single continuous transient (step 3) -- unlike the
% ssp585 script, there is no need for a mid-run restart split, since the
% forcing never changes and there's nothing for a restart to help with here.
%
% Forcing sources (all under raw_data/ISMIP7/AIS/CESM2-WACCM/ssp585/, since
% there is no separate "CTRL" raw data -- the control run reuses ssp585's
% real 2015 initial fields):
%   Ocean TF   : ocean/tf/v3/tf_AIS_CESM2-WACCM_ssp585_ocean_v3_2015-2024.nc (year 2015 slice)
%   SMB        : SDBN1-2000m/acabf/v2/acabf_AIS_CESM2-WACCM_ssp585_SDBN1-2000m_v2_2015.nc
%   SMB grad   : SDBN1-2000m/dacabfdz/v2/dacabfdz_AIS_CESM2-WACCM_ssp585_SDBN1-2000m_v2_2015.nc
%   SMB anomaly baseline: CESM historical 1995-2014 mean (same as hist_run_tune_CESM_WACCM)
%   RACMO climatology: raw_data/nc_orig/Atmosphere/smb_rec.mean.1995-2014...
%
% Starting state: hist_runs/CESM2-WACCM/Models/AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat
%   (end-state of hist_run_CESM_WACCM_1995_2014.m step 3 -- same starting
%   point as the transient ssp585 run, so any divergence between the two
%   is attributable to the forcing, not the initial state)
%
% Step map:
%   1  ProjTF            build a single 2015 TF snapshot, save to
%                        preprocessed_data/Ocean/Ctrl/CESM2-WACCM/ctrl/
%   2  ProjSMB           build a single 2015 SMB+gradient snapshot, save to
%                        preprocessed_data/Atmosphere/Ctrl/CESM2-WACCM/ctrl/
%   3  ProjRun_2015_2300 single continuous transient 2015-2300 starting from
%                        AIS_state_2015, forcing held at the 2015 snapshot
%                        throughout, calving front frozen (ismovingfront=0)
%                        (loadonly=0 submit, =1 gather)
%   4  WriteISMIP6_NetCDF_test writes 4 variables to ctrl_test/ for quick
%                        compliance checking (sftgif, iareagr, tendlibmassbfgr,
%                        dlithkdt)
%   5  WriteISMIP6_NetCDF grids the TransientSolution onto the ISMIP6 8 km
%                        AIS grid and writes one NetCDF per Appendix-2
%                        variable into postprocessed_data/CESM2-WACCM/ctrl/
%                        (EXP='C009')

    % ---- swappable scenario label ------------------------------------------
    SCENARIO_SRC = 'ssp585';    % real CESM2-WACCM scenario whose 2015 forcing
                                 % is reused (source-data paths only)
    EXP_LABEL    = 'ctrl';      % this control experiment's own label, used
                                 % for all output paths/naming/NetCDF metadata
    CMIP_MODEL   = 'CESM2-WACCM';
    % -------------------------------------------------------------------------

    if ~exist('loadonly','var'), loadonly = 0; end
    addpath('./../../../init/scripts');

    org = organizer('repository', './Models/', ...
                    'prefix', ['AIS_ISMIP7_Proj_' EXP_LABEL '_'], ...
                    'steps', steps, 'color', '34;47;2');
    clear steps;

    % ------------------------------------------------------------------ paths
    proj_root  = './../../../';
    init_dir   = [proj_root 'init/'];

    inputmodel_2015 = [proj_root 'hist_runs/CESM2-WACCM/Models/' ...
                       'AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat'];

    raw_ssp     = [proj_root 'raw_data/ISMIP7/AIS/' CMIP_MODEL '/' SCENARIO_SRC '/'];
    raw_hist    = [proj_root 'raw_data/ISMIP7/AIS/' CMIP_MODEL '/historical/'];

    tf_dir      = [raw_ssp  'ocean/tf/v3/'];
    smb_ssp_dir = [raw_ssp  'SDBN1-2000m/acabf/v2/'];
    grad_ssp_dir= [raw_ssp  'SDBN1-2000m/dacabfdz/v2/'];
    smb_hist_dir= [raw_hist 'SDBN1-2000m/acabf/v2/'];

    preproc_ocean      = [proj_root 'preprocessed_data/Ocean/'];
    preproc_hist_atmo  = [proj_root 'preprocessed_data/Atmosphere/Hist/'];
    preproc_clim       = [proj_root 'preprocessed_data/Atmosphere/Clim/' CMIP_MODEL '/'];
    preproc_proj_ocean = [proj_root 'preprocessed_data/Ocean/Ctrl/' CMIP_MODEL '/' EXP_LABEL '/'];
    preproc_proj_atmo  = [proj_root 'preprocessed_data/Atmosphere/Ctrl/' CMIP_MODEL '/' EXP_LABEL '/'];

    % ------------------------------------------------------------------ time
    fixed_year  = 2015;   % the single forcing year held constant throughout
    start_year  = 2015;
    end_year    = 2300;   % final snapshot lands at t=end_year+1=2301 → nominal
                          % year 2300 (matches the transient ssp585 run's duration)
    sec_to_year = 31556926;   % consistent with hist_run_tune_CESM_WACCM

    % ------------------------------------------------------------------ ISMIP6 outputs (same as hist run)
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
    if perform(org, 'ProjTF') % {{{
        % Build a SINGLE 2015 TF snapshot (not a full time series) -- this is
        % the control run's whole ocean forcing; interp_forcing=0 holds it
        % constant for the entire 2015-2300 run. Only the decade chunk that
        % actually contains 2015 is read (2015-2024), not the full series.
        md = loadmodel(inputmodel_2015);

        if ~exist(preproc_proj_ocean, 'dir'), mkdir(preproc_proj_ocean); end

        tf_files = dir([tf_dir 'tf_AIS_*.nc']);
        tf_files = sort({tf_files.name});

        fpath = '';
        for fi = 1:length(tf_files)
            [~, fname, ~] = fileparts(tf_files{fi});
            parts  = strsplit(fname, '_');
            decade = strsplit(parts{end}, '-');
            yr0 = str2double(decade{1}); yr1 = str2double(decade{2});
            if fixed_year >= yr0 && fixed_year <= yr1
                fpath = [tf_dir tf_files{fi}];
                ti    = fixed_year - yr0 + 1;
                break;
            end
        end
        if isempty(fpath)
            error('No TF file found containing year %d.', fixed_year);
        end

        z_data  = double(ncread(fpath, 'z'));
        nDepths = length(z_data);

        x_n = double(ncread(fpath, 'x'));
        y_n = double(ncread(fpath, 'y'));

        tf_data_all = double(ncread(fpath, 'tf', [1 1 1 ti], [Inf Inf Inf 1]));  % [x,y,depth]

        tf_proj = cell(1, 1, nDepths);
        for i = 1:nDepths
            v = InterpFromGridToMesh(x_n, y_n, tf_data_all(:,:,i)', md.mesh.x, md.mesh.y, 0);
            tf_proj{1,1,i} = [max(v, 0) ; fixed_year];
        end

        save([preproc_proj_ocean 'CESM_WACCM_TF_' EXP_LABEL '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat'], ...
             'tf_proj', 'z_data', '-v7.3');
        fprintf('Saved TF: single %d snapshot, %d depths.\n', fixed_year, nDepths);
    end % }}}

    % ================================================================= Step 2
    if perform(org, 'ProjSMB') % {{{
        % Build a SINGLE 2015 SMB+gradient snapshot -- the control run's
        % whole SMB forcing; interp_forcing=0 holds it constant for the
        % entire 2015-2300 run.
        % Convention (identical to hist_run_tune_CESM_WACCM step 2):
        %   smb_yr = smb_racmo + (cesm_yr_ssp - cesm_hist_mean_1995_2014)
        % Units: mm w.e. yr-1 for smbref; mm w.e. yr-1 m-1 for b_pos/b_neg.
        md = loadmodel(inputmodel_2015);

        if ~exist(preproc_proj_atmo, 'dir'), mkdir(preproc_proj_atmo); end

        % Load smb_racmo (RACMO clim on mesh), cesm_mean (CESM 1995-2014 mean
        % on mesh), and p_vert (region scaling factor) from mats saved by
        % hist_run_tune_CESM_WACCM.m.
        clim_yr0 = 1995; clim_yr1 = 2014;
        hist_smb_yr0 = 1995; hist_smb_yr1 = 2020;
        load([preproc_clim 'CESM_WACCM_SMB_clim_' ...
              num2str(clim_yr0) '_' num2str(clim_yr1) '.mat'], 'smb_racmo', 'cesm_mean');
        load([preproc_hist_atmo 'CESM_WACCM_SMB_corrected_' ...
              num2str(hist_smb_yr0) '_' num2str(hist_smb_yr1) '.mat'], 'p_vert');
        cesm_hist_mean = cesm_mean;

        nc_smb  = [smb_ssp_dir sprintf('acabf_AIS_%s_%s_SDBN1-2000m_v2_%d.nc', CMIP_MODEL, SCENARIO_SRC, fixed_year)];
        nc_grad = [grad_ssp_dir sprintf('dacabfdz_AIS_%s_%s_SDBN1-2000m_v2_%d.nc', CMIP_MODEL, SCENARIO_SRC, fixed_year)];

        % Cheap check: an empty placeholder file has 'time' UNLIMITED with 0
        % records currently (seen for some later years in this series --
        % confirm 2015 itself is real before reading the full spatial field).
        if isempty(ncread(nc_smb, 'time')) || isempty(ncread(nc_grad, 'time'))
            error('acabf/dacabfdz for %d has no data (empty placeholder).', fixed_year);
        end

        x_s = double(ncread(nc_smb, 'x'));
        y_s = double(ncread(nc_smb, 'y'));

        am      = mean(double(ncread(nc_smb, 'acabf')), 3);
        cesm_yr = InterpFromGridToMesh(x_s, y_s, am', md.mesh.x, md.mesh.y, 0) * sec_to_year;
        smb_col = p_vert .* smb_racmo + (cesm_yr - cesm_hist_mean);

        g_raw = squeeze(double(ncread(nc_grad, 'dacabfdz')));
        if ndims(g_raw) == 3
            g_raw = mean(g_raw, 3);
        end
        bgrad_col = InterpFromGridToMesh(x_s, y_s, g_raw', md.mesh.x, md.mesh.y, 0) * sec_to_year;

        smb_forcing   = [smb_col   ; fixed_year];
        bgrad_forcing = [bgrad_col ; fixed_year];

        save([preproc_proj_atmo 'CESM_WACCM_SMB_' EXP_LABEL '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat'], ...
             'smb_forcing', 'bgrad_forcing', '-v7.3');
        fprintf('Saved SMB forcing: single %d snapshot.\n', fixed_year);
    end % }}}

    % ================================================================= Step 3
    if perform(org, 'ProjRun_2015_2300') % {{{
        % Single continuous transient 2015-2300. No restart split is needed
        % (unlike the ssp585 script's ProjRun_2015_2150 / ProjRun_2151_2300
        % pair): forcing is held at the single 2015 snapshot for the whole
        % run, so there's nothing a restart would help with here.
        % Starting from the AIS_state_2015 end-state of the historical run
        % (hist_run_CESM_WACCM_1995_2014 step 3) -- same starting point as
        % the transient ssp585 run. loadonly=0 submits to PBS, =1 gathers.

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
        md.timestepping.interp_forcing  = 0;

        md.timestepping.start_time = start_year;
        md.timestepping.final_time = end_year + 1;   % last snapshot at 2300
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;

        md.transient.requested_outputs = ismip6_outputs;

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean (single fixed-2015 TF, held constant throughout) ---
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_proj_ocean 'CESM_WACCM_TF_' EXP_LABEL '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        load([preproc_ocean 'gamma0_local.mat']);

        unique_basinid = unique(basinid);
        tmp     = load([preproc_ocean 'dT_correction.mat'], 'dT_correction');
        delta_t = tmp.dT_correction;   % 1 x nBasins, from meltMip_ensemble.m step 8 (get_dT_iterate_BMB_j)

        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique_basinid);
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = tf_proj;
        md.basalforcings.islocal    = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_local;

        % --- SMB (single fixed-2015 SMB, held constant throughout) ---
        load([preproc_proj_atmo 'CESM_WACCM_SMB_' EXP_LABEL '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);

        % href must be the 1995 relaxed surface, not the 2015 start state.
        % smbref is calibrated against the 1995 RACMO climatology, so the
        % lapse-rate correction b*(surface-href) must be zero at 1995 conditions.
        md_relax = loadmodel([init_dir 'Models/AIS_ISMIP7_Relaxed_CESM_WACCM.mat']);
        surf_1995 = md_relax.geometry.surface;
        clear md_relax

        md.smb        = SMBgradients();
        md.smb.smbref = smb_forcing;
        md.smb.b_pos  = bgrad_forcing;
        md.smb.b_neg  = bgrad_forcing;
        md.smb.href   = [surf_1995 ; 1995];

        % --- Calving front: frozen at the AIS_state_2015 starting position ---
        % ismovingfront=0 means the level-set equation is never solved, so
        % the ice front/grounding line just evolve passively with the ice
        % dynamics -- no ice-shelf collapse mask or spclevelset time series
        % is used at all.
        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 0;

        md.miscellaneous.name = ['ProjRun_' CMIP_MODEL '_' EXP_LABEL '_' ...
                                 num2str(start_year) '_' num2str(end_year)];
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
    if perform(org, 'VAFContinuityCheck') % {{{
        % Load both projection segments and verify the VAF time series joins
        % continuously at the 2151 restart boundary. Plots sea-level-equivalent
        % contribution relative to 2015, coloured by segment.

        md1 = loadmodel(org, 'ProjRun_2015_2300');

        time1 = [md1.results.TransientSolution.time];
        vaf1  = [md1.results.TransientSolution.IceVolumeAboveFloatationScaled];

        % IceVolumeAboveFloatationScaled is in m³ (ice volume).
        % SLE [m] = -ΔVAF [m³] * rho_ice / rho_sw / A_ocean
        rho_ice = 917;       % kg/m³
        rho_sw  = 1028;      % kg/m³
        A_ocean = 3.625e14;  % m² (ISMIP6 standard ocean area)
        vaf_ref = vaf1(1);
        sle1 = -(vaf1 - vaf_ref) * rho_ice / rho_sw / A_ocean;

        fig = figure('visible', 'off');
        plot(time1, sle1, 'b-', 'LineWidth', 1.5); hold on;
        xlabel('Year');
        ylabel('Sea level contribution (m SLE)');
        title(sprintf('%s %s — VAF SLE', CMIP_MODEL,EXP_LABEL), ...
              'Interpreter', 'none');
        grid on;

        figdir = [proj_root 'postprocessed_data/figures/' CMIP_MODEL '/' EXP_LABEL '/'];
        if ~exist(figdir, 'dir'), mkdir(figdir); end
        figname = fullfile(figdir, sprintf('VAF_continuity_%s_%s.png', CMIP_MODEL,EXP_LABEL));
        saveas(fig, figname);
        close(fig);
        fprintf('Saved: %s\n', figname);
    end % }}}
    % ================================================================= Step 5
    if perform(org, 'WriteISMIP6_NetCDF_test') % {{{
        % Quick 4-variable test: sftgif (ST 2D), iareagr (ST scalar),
        % tendlibmassbfgr (FL scalar), dlithkdt (FL 2D).
        % Output goes to postprocessed_data/CESM2-WACCM/ctrl_test/
        addpath('./../../functions');

        md = loadmodel(org, 'ProjRun_2015_2300');

        outdir_test = [proj_root 'postprocessed_data/CESM2-WACCM/' EXP_LABEL '_test/'];
        if ~exist(outdir_test, 'dir'), mkdir(outdir_test); end

        meta_t                    = struct();
        meta_t.experiment_id      = EXP_LABEL;
        meta_t.set_counter        = 'C009';
        meta_t.time_range         = '2015-2300';
        meta_t.ESM_id             = CMIP_MODEL;
        meta_t.forcing_member_id  = 'f001';
        meta_t.ISM_member_id      = 'm001';

        % ---- shared time setup (mirrors write_ismip7_2d_projection) ----------
        t_raw       = [md.results.TransientSolution.time];
        keep        = abs(t_raw - round(t_raw)) < 0.05;
        ts          = md.results.TransientSolution(keep);
        nT          = sum(keep);
        t_annual     = round(t_raw(keep));   % ISSM time [2016, ..., 2300]
        time_yr      = t_annual - 1;         % calendar year [2015, ..., 2299]
        ref_dn       = datenum(1850, 1, 1);
        time_st      = zeros(1, nT);
        lb_fl        = zeros(1, nT); ub_fl = zeros(1, nT); time_fl = zeros(1, nT);
        for i = 1:nT
            time_st(i) = datenum(t_annual(i), 1, 1) - ref_dn;
            lb_fl(i)   = datenum(time_yr(i),  1, 1) - ref_dn;
            ub_fl(i)   = datenum(time_yr(i)+1,1, 1) - ref_dn;
            time_fl(i) = datenum(time_yr(i),  7, 1) - ref_dn;
        end
        time_bnds_fl = [lb_fl; ub_fl]';

        [~, xGrid, yGrid] = gridData(md, ts(1).Thickness, ...
            'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]);
        nx = length(xGrid); ny = length(yGrid);

        % ---- sftgif (ST, 2D) -------------------------------------------------
        sftgif3d = zeros(ny, nx, nT, 'single');
        for t = 1:nT
            ice = gridData(md, double(ts(t).MaskIceLevelset <= 0), ...
                'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]);
            sftgif3d(:,:,t) = single(max(0, min(1, ice)));
        end
        write_ismip7_2d(outdir_test, 'sftgif', sftgif3d, xGrid, yGrid, ...
            time_st, [], 'land_ice_area_fraction', 'Land ice area fraction', '1', meta_t);

        % ---- iareagr (ST, scalar) --------------------------------------------
        iareagr_v = [ts.GroundedAreaScaled];
        write_ismip7_scalar(outdir_test, 'iareagr', iareagr_v, time_st, [], ...
            'grounded_ice_sheet_area', 'Grounded ice area', 'm^2', meta_t);

        % ---- tendlibmassbfgr (FL, scalar) ------------------------------------
        rhoi = md.materials.rho_ice;
        tendlibmassbfgr_v = -[ts.TotalGroundedBmbScaled] * rhoi / md.constants.yts;
        write_ismip7_scalar(outdir_test, 'tendlibmassbfgr', tendlibmassbfgr_v, ...
            time_fl, time_bnds_fl, ...
            'tendency_of_land_ice_mass_due_to_basal_mass_balance', ...
            'Grounded basal mass balance flux', 'kg s-1', meta_t);

        % ---- dlithkdt (FL, 2D) -----------------------------------------------
        H_grid = zeros(ny, nx, nT);
        for t = 1:nT
            H_grid(:,:,t) = gridData(md, ts(t).Thickness, ...
                'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]);
        end
        dlithkdt3d = zeros(ny, nx, nT, 'single');
        for t = 1:nT
            if t < nT
                dlithkdt3d(:,:,t) = single((H_grid(:,:,t+1) - H_grid(:,:,t)) / md.constants.yts);
            else
                dlithkdt3d(:,:,t) = single((H_grid(:,:,t) - H_grid(:,:,t-1)) / md.constants.yts);
            end
        end
        safe_max = single(1e-4);
        if double(safe_max) > 1e-4, safe_max = safe_max - eps(safe_max); end
        dlithkdt3d = max(single(-1e-4), min(safe_max, dlithkdt3d));
        write_ismip7_2d(outdir_test, 'dlithkdt', dlithkdt3d, xGrid, yGrid, ...
            time_fl, time_bnds_fl, 'tendency_of_land_ice_thickness', ...
            'Ice thickness tendency', 'm s-1', meta_t);

        fprintf('Test output written to %s\n', outdir_test);
    end % }}}

    % ================================================================= Step 6
    if perform(org, 'WriteISMIP6_NetCDF') % {{{
        % Grids the TransientSolution (step 3) onto the standard ISMIP7
        % 761×761 8 km AIS grid. See proj_runs/functions/ for the full
        % implementation.
        addpath('./../../functions');

        md = loadmodel(org, 'ProjRun_2015_2300');

        outdir = [proj_root 'postprocessed_data/CESM2-WACCM/' EXP_LABEL '/'];
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        meta                    = struct();
        meta.experiment_id      = EXP_LABEL;        % 'ctrl'
        meta.set_counter        = 'C009';
        meta.time_range         = '2015-2300';
        meta.ESM_id             = CMIP_MODEL;      % 'CESM2-WACCM'
        meta.forcing_member_id  = 'f001';
        meta.ISM_member_id      = 'm001';

        [cfflux_tot, glflux_tot] = write_ismip7_2d_projection(md, outdir, meta);
        write_ismip7_scalar_projection(md, outdir, meta, cfflux_tot, glflux_tot);
    end % }}}

end

% NetCDF writing delegated to proj_runs/functions/write_ismip7_2d_projection.m
% and proj_runs/functions/write_ismip7_scalar_projection.m.
