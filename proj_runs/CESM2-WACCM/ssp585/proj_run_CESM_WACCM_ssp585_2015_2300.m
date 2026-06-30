function md = proj_run_CESM_WACCM_ssp585_2015_2300(steps, loadonly)
% Projection run 2015-2300 using CESM2-WACCM ssp585 forcing.
%
% To run a different emission scenario, copy this file, rename the function
% and file to match, and change SCENARIO below:
%   ssp126 | ssp370 | ssp534-over | ssp585
%
% Forcing sources (all under raw_data/ISMIP7/AIS/CESM2-WACCM/<SCENARIO>/):
%   Ocean TF   : ocean/tf/v3/  (decade NetCDF chunks 2015-2299; 2299 repeated for 2300)
%   SMB        : SDBN1-2000m/acabf/v2/    (annual NetCDF, kg m-2 s-1, 2015-2300)
%   SMB grad   : SDBN1-2000m/dacabfdz/v2/ (annual NetCDF, kg m-2 s-1 m-1, 2015-2300)
%   Levelset   : fracture/ice_shelf_collapse_mask_<model>_<scenario>_ismip7_8km.nc
%                Annual collapse flag (0=sustainable, 1=collapsed); made cumulative
%                here so once a shelf collapses it stays constrained.
%   SMB anomaly baseline: CESM historical 1995-2014 mean (same as hist_run_tune_CESM_WACCM)
%   RACMO climatology: raw_data/nc_orig/Atmosphere/smb_rec.mean.1995-2014...
%
% Starting state: hist_runs/CESM2-WACCM/Models/AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat
%   (end-state of hist_run_CESM_WACCM_1995_2014.m step 3)
%
% Step map:
%   1  ProjTF            build annual TF mat 2015-2300, save to
%                        preprocessed_data/Ocean/Proj/CESM2-WACCM/ssp585/
%   2  ProjSMB           build annual SMB mat 2015-2300 (CESM anomaly
%                        relative to 1995-2014 historical mean + RACMO clim),
%                        save to preprocessed_data/Atmosphere/Proj/CESM2-WACCM/ssp585/
%   3  ProjLevelset      build cumulative ice-shelf-collapse spclevelset 2015-2300,
%                        save to preprocessed_data/Ocean/Proj/CESM2-WACCM/ssp585/
%   4  ProjRun_2015_2150 transient 2015-2150 starting from AIS_state_2015
%                        (loadonly=0 submit, =1 gather)
%   5  AIS_state_2151    save end-state geometry/masks at 2151 with
%                        TransientSolution cleared (restart point for step 6)
%   6  ProjRun_2151_2300 transient 2151-2300 starting from AIS_state_2151
%                        (loadonly=0 submit, =1 gather)
%   7  VAFContinuityCheck plot VAF (m SLE) from both segments on one axes;
%                        checks junction at 2151 is continuous; saves figure to
%                        postprocessed_data/figures/CESM2-WACCM/ssp585/
%   8  WriteISMIP6_NetCDF grids the combined TransientSolution (steps 4+6)
%                        onto the ISMIP6 8 km AIS grid and writes one NetCDF
%                        per Appendix-2 variable into
%                        postprocessed_data/CESM2-WACCM/proj/  (EXP='C007')
%
%  noSEF sensitivity run lives in a separate script:
%    proj_run_CESM_WACCM_ssp585_2015_2300_noSEF.m  (submit via submit_proj_ssp585_noSEF.sh)

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

    raw_ssp     = [proj_root 'raw_data/ISMIP7/AIS/' CMIP_MODEL '/' SCENARIO '/'];
    raw_hist    = [proj_root 'raw_data/ISMIP7/AIS/' CMIP_MODEL '/historical/'];

    tf_dir      = [raw_ssp  'ocean/tf/v3/'];
    smb_ssp_dir = [raw_ssp  'SDBN1-2000m/acabf/v2/'];
    grad_ssp_dir= [raw_ssp  'SDBN1-2000m/dacabfdz/v2/'];
    smb_hist_dir= [raw_hist 'SDBN1-2000m/acabf/v2/'];
    collapse_nc = [raw_ssp  'fracture/ice_shelf_collapse_mask_' ...
                   lower(strrep(CMIP_MODEL,'-','')) '_' SCENARIO '_ismip7_8km.nc'];

    preproc_ocean      = [proj_root 'preprocessed_data/Ocean/'];
    preproc_hist_atmo  = [proj_root 'preprocessed_data/Atmosphere/Hist/'];
    preproc_clim       = [proj_root 'preprocessed_data/Atmosphere/Clim/' CMIP_MODEL '/'];
    preproc_proj_ocean = [proj_root 'preprocessed_data/Ocean/Proj/' CMIP_MODEL '/' SCENARIO '/'];
    preproc_proj_atmo  = [proj_root 'preprocessed_data/Atmosphere/Proj/' CMIP_MODEL '/' SCENARIO '/'];

    % ------------------------------------------------------------------ time
    start_year  = 2015;
    mid_year    = 2150;   % phase 1 ends here (state saved at 2151)
    end_year    = 2299;   % last year with available forcing; final snapshot
                          % lands at t=end_year+1=2300 (same convention as
                          % hist_run: final_time = end_year+1)
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
        % Build annual TF matrix 2015-2300 from decade NetCDF chunks.
        % TF files cover 2015-2299; year 2300 is filled by repeating 2299.
        % Fields are clamped >= 0 (basalforcingsismip6 consistency check).
        md = loadmodel(inputmodel_2015);

        if ~exist(preproc_proj_ocean, 'dir'), mkdir(preproc_proj_ocean); end

        tf_files = dir([tf_dir 'tf_AIS_*.nc']);
        tf_files = sort({tf_files.name});

        z_data  = double(ncread([tf_dir tf_files{1}], 'z'));
        nDepths = length(z_data);
        nVerts  = md.mesh.numberofvertices;

        % Forcing array spans start_year:end_year+1 so the final_time=end_year+1
        % time step has a valid TF value; end_year+1 is filled by repeating end_year.
        years   = start_year : end_year + 1;
        nyears  = length(years);
        tf_mat  = zeros(nVerts, nyears, nDepths);
        t_vec   = years;
        tf_last = [];

        for fi = 1:length(tf_files)
            fpath = [tf_dir tf_files{fi}];
            [~, fname, ~] = fileparts(fpath);
            parts      = strsplit(fname, '_');
            decade     = strsplit(parts{end}, '-');
            file_years = str2double(decade{1}) : str2double(decade{2});
            x_n = double(ncread(fpath, 'x'));
            y_n = double(ncread(fpath, 'y'));
            for ti = 1:length(file_years)
                yr = file_years(ti);
                if yr < start_year || yr > end_year, continue; end
                ki      = yr - start_year + 1;
                tf_data = double(ncread(fpath, 'tf', [1 1 1 ti], [Inf Inf Inf 1]));
                for i = 1:nDepths
                    v = InterpFromGridToMesh(x_n, y_n, tf_data(:,:,i)', ...
                                            md.mesh.x, md.mesh.y, 0);
                    tf_mat(:, ki, i) = max(v, 0);
                end
                if yr == end_year
                    tf_last = squeeze(tf_mat(:, ki, :));
                end
            end
        end

        % Fill end_year+1 by repeating end_year (no TF data beyond end_year)
        if ~isempty(tf_last)
            tf_mat(:, end, :) = tf_last;
            fprintf('[INFO] TF %d filled by repeating %d.\n', end_year+1, end_year);
        else
            warning('TF for end_year (%d) not found; end_year+1 entry will be zero.', end_year);
        end

        tf_proj = cell(1, 1, nDepths);
        for i = 1:nDepths
            slice = squeeze(tf_mat(:,:,i));
            tf_proj{1,1,i} = [slice ; t_vec];
        end

        save([preproc_proj_ocean 'CESM_WACCM_TF_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat'], ...
             'tf_proj', 'z_data', 't_vec', '-v7.3');
        fprintf('Saved TF: %d years (%d-%d), %d depths.\n', nyears, start_year, end_year+1, nDepths);
    end % }}}

    % ================================================================= Step 2
    if perform(org, 'ProjSMB') % {{{
        % Build annual SMB matrix 2015-2300.
        % Convention (identical to hist_run_tune_CESM_WACCM step 2):
        %   smb_yr = smb_racmo + (cesm_yr_ssp - cesm_hist_mean_1995_2014)
        % Units: mm w.e. yr-1 for smbref; mm w.e. yr-1 m-1 for b_pos/b_neg.
        % CESM acabf / dacabfdz are monthly (12 per year); collapsed to annual
        % mean (equal-weight across 12 months).
        md = loadmodel(inputmodel_2015);

        if ~exist(preproc_proj_atmo, 'dir'), mkdir(preproc_proj_atmo); end

        nVerts = md.mesh.numberofvertices;

        % Load smb_racmo (RACMO clim on mesh), cesm_mean (CESM 1995-2014 mean
        % on mesh), and p_vert (region scaling factor) from mats saved by
        % hist_run_tune_CESM_WACCM.m -- no need to re-read the RACMO or 20
        % years of historical CESM files here.
        clim_yr0 = 1995; clim_yr1 = 2014;
        hist_smb_yr0 = 1995; hist_smb_yr1 = 2020;
        load([preproc_clim 'CESM_WACCM_SMB_clim_' ...
              num2str(clim_yr0) '_' num2str(clim_yr1) '.mat'], 'smb_racmo', 'cesm_mean');
        load([preproc_hist_atmo 'CESM_WACCM_SMB_corrected_' ...
              num2str(hist_smb_yr0) '_' num2str(hist_smb_yr1) '.mat'], 'p_vert');
        cesm_hist_mean = cesm_mean;

        % Read x/y grid coordinates from the first projection year file
        % (same CESM grid as historical; needed for InterpFromGridToMesh in loop)
        nc_first = [smb_ssp_dir sprintf('acabf_AIS_%s_%s_SDBN1-2000m_v2_%d.nc', CMIP_MODEL, SCENARIO, start_year)];
        x_s = double(ncread(nc_first, 'x'));
        y_s = double(ncread(nc_first, 'y'));

        % Build annual projection matrices 2015-2299
        years        = start_year : end_year;
        nyears       = length(years);
        smb_matrix   = zeros(nVerts, nyears);
        bgrad_matrix = zeros(nVerts, nyears);
        t_smb        = years;

        for k = 1:nyears
            yr = years(k);

            nc_smb = [smb_ssp_dir sprintf('acabf_AIS_%s_%s_SDBN1-2000m_v2_%d.nc', CMIP_MODEL, SCENARIO, yr)];
            am     = mean(double(ncread(nc_smb, 'acabf')), 3);
            cesm_yr = InterpFromGridToMesh(x_s, y_s, am', md.mesh.x, md.mesh.y, 0) * sec_to_year;
            smb_matrix(:,k) = p_vert .* smb_racmo + (cesm_yr - cesm_hist_mean);

            nc_grad  = [grad_ssp_dir sprintf('dacabfdz_AIS_%s_%s_SDBN1-2000m_v2_%d.nc', CMIP_MODEL, SCENARIO, yr)];
            g_raw    = squeeze(double(ncread(nc_grad, 'dacabfdz')));
            if ndims(g_raw) == 3
                g_raw = mean(g_raw, 3);
            end
            bgrad_matrix(:,k) = InterpFromGridToMesh(x_s, y_s, g_raw', md.mesh.x, md.mesh.y, 0) * sec_to_year;

            if mod(yr, 10) == 0
                fprintf('  SMB+grad year %d\n', yr);
            end
        end

        smb_forcing   = [smb_matrix  ; t_smb];
        bgrad_forcing = [bgrad_matrix ; t_smb];

        save([preproc_proj_atmo 'CESM_WACCM_SMB_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat'], ...
             'smb_forcing', 'bgrad_forcing', 't_smb', '-v7.3');
        fprintf('Saved SMB forcing.\n');
    end % }}}

    % ================================================================= Step 3
    if perform(org, 'ProjLevelset') % {{{
        % Build a spclevelset time series from the annual ice-shelf-collapse mask.
        % mask=0: ice shelf sustainable (no constraint -> NaN in spclevelset).
        % mask=1: ice shelf collapsed (levelset = +1 = open ocean).
        % Made cumulative: once a grid cell collapses in year Y, it stays
        % constrained as ocean for all Y' > Y (prevents spurious re-advance
        % since the physics alone may not fully suppress it).
        % Advance-suppression: vertices with no ice in the initial state are
        % not forced to NaN -- they start as ocean and the constraint is
        % irrelevant there.

        md = loadmodel(inputmodel_2015);
        nVerts = md.mesh.numberofvertices;

        x_c = double(ncread(collapse_nc, 'x'));
        y_c = double(ncread(collapse_nc, 'y'));
        time_c = double(ncread(collapse_nc, 'time'));  % integer years

        model_years       = start_year : end_year + 1;   % include final_time year
        spclevelset_mat   = NaN(nVerts, length(model_years));
        collapsed_so_far  = false(nVerts, 1);  % cumulative collapse state

        for k = 1:length(model_years)
            yr = model_years(k);

            % Find closest available year in the mask file (covers 1950-2299)
            yr_clamped = min(yr, max(time_c));
            [~, gi] = min(abs(time_c - yr_clamped));

            mask_raw = double(ncread(collapse_nc, 'mask', [1 1 gi], [Inf Inf 1]));
            % mask is stored (y,x) in file; x ascending, y ascending
            % Interpolate to mesh: treat 1 (collapsed) as positive
            v = InterpFromGridToMesh(x_c, y_c, mask_raw', md.mesh.x, md.mesh.y, 0);

            % Update cumulative collapse flag
            collapsed_this_yr = (v >= 0.5);
            collapsed_so_far  = collapsed_so_far | collapsed_this_yr;

            % Constrain: collapsed cells -> ocean (levelset = +1)
            % Uncollapsed cells -> no constraint (NaN)
            lset = NaN(nVerts, 1);
            lset(collapsed_so_far) = 1;
            spclevelset_mat(:, k) = lset;

            if mod(yr, 50) == 0 || yr == start_year
                fprintf('  Levelset year %d: %d cells newly collapsed, %d total constrained\n', ...
                        yr, sum(collapsed_this_yr), sum(collapsed_so_far));
            end
        end

        % Advance-suppression: vertices that already had no ice in the 2015
        % state cannot be further "forced" to no-ice -- those NaN constraints
        % are fine; we suppress the reverse case (trying to force +1 at a
        % vertex that starts with ice, which the collapse mask should always
        % want for retreating fronts, so this mainly prevents data artefacts).
        no_ice0      = md.mask.ice_levelset > 0;
        already_ocean = repmat(no_ice0, 1, length(model_years)) & isnan(spclevelset_mat);
        % no action needed -- NaN at vertices that start as ocean is correct

        proj_spclevelset = [spclevelset_mat ; model_years];
        save([preproc_proj_ocean 'CESM_WACCM_levelset_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat'], ...
             'proj_spclevelset', '-v7.3');
        fprintf('Saved collapse levelset mat.\n');
        clear already_ocean no_ice0;
    end % }}}

    % ================================================================= Step 4
    if perform(org, 'ProjRun_2015_2150') % {{{
        % Transient 2015-2150. Starting from the AIS_state_2015 end-state of
        % the historical run (hist_run_CESM_WACCM_1995_2014 step 3).
        % loadonly=0 submits to PBS; loadonly=1 gathers results.

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
        md.timestepping.final_time = mid_year + 1;   % last snapshot at 2151
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

        % Trim forcing to [start_year, mid_year+1] — no need to carry 2300
        % data into a run that ends at 2151; halves the .bin size.
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

        % --- SMB ---
        load([preproc_proj_atmo 'CESM_WACCM_SMB_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        keep4         = smb_forcing(end,:) <= mid_year + 1;
        smb_forcing   = smb_forcing(:, keep4);
        bgrad_forcing = bgrad_forcing(:, keep4);

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

        % --- Calving front (collapse levelset) ---
        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 1;
        load([preproc_proj_ocean 'CESM_WACCM_levelset_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        keep4            = proj_spclevelset(end,:) <= mid_year + 1;
        proj_spclevelset = proj_spclevelset(:, keep4);
        md.levelset.spclevelset = proj_spclevelset;

        md.miscellaneous.name = ['ProjRun_' CMIP_MODEL '_' SCENARIO '_' ...
                                 num2str(start_year) '_' num2str(mid_year)];
        clustername = 'gadi';
        md.cluster  = set_cluster(clustername);
        md.settings.waitonlock = 0;
        md.verbose = verbose('solution', true, 'module', true, 'convergence', true);
        md = solve(md, 'tr', 'runtimename', false, 'loadonly', loadonly);
        if loadonly
            savemodel(org, md);
        end
    end % }}}

    % ================================================================= Step 5
    if perform(org, 'AIS_state_2151') % {{{
        % Save the end-state at 2151 (final snapshot of ProjRun_2015_2150)
        % as a clean restart for ProjRun_2151_2300, with TransientSolution
        % cleared to keep file size manageable.
        md = loadmodel(org, 'ProjRun_2015_2150');
        md_in = md;
        md.geometry.thickness = md_in.results.TransientSolution(end).Thickness;
        md.geometry.surface   = md_in.results.TransientSolution(end).Surface;
        md.geometry.base      = md_in.results.TransientSolution(end).Base;
        md.mask.ocean_levelset = md_in.results.TransientSolution(end).MaskOceanLevelset;
        md.mask.ice_levelset   = md_in.results.TransientSolution(end).MaskIceLevelset;
        md.results.TransientSolution = [];
        savemodel(org, md);
    end % }}}

    % ================================================================= Step 6
    if perform(org, 'ProjRun_2151_2300') % {{{
        % Transient 2151-2300. Continues from AIS_state_2151.
        % Uses the same full 2015-2300 forcing arrays loaded in step 4 --
        % ISSM truncates to [start_time, final_time] itself, so providing
        % the full series here is intentional and not wasteful of memory
        % (ISSM only evaluates the forcing at the current time step).

        md = loadmodel(org, 'AIS_state_2151');
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

        md.timestepping.start_time = mid_year + 1;    % 2151
        md.timestepping.final_time = end_year + 1;   % last snapshot at 2300
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

        % Trim forcing to [mid_year+1, end_year+1] only.
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

        % --- SMB ---
        load([preproc_proj_atmo 'CESM_WACCM_SMB_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);
        keep6         = smb_forcing(end,:) >= mid_year + 1;
        smb_forcing   = smb_forcing(:, keep6);
        bgrad_forcing = bgrad_forcing(:, keep6);

        % href = 1995 reference surface (same as segment 1 and the hist run).
        md_relax = loadmodel([init_dir 'Models/AIS_ISMIP7_Relaxed_CESM_WACCM.mat']);
        surf_1995 = md_relax.geometry.surface;
        clear md_relax

        md.smb        = SMBgradients();
        md.smb.smbref = smb_forcing;
        md.smb.b_pos  = bgrad_forcing;
        md.smb.b_neg  = bgrad_forcing;
        md.smb.href   = [surf_1995 ; 1995];

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
                                 num2str(mid_year+1) '_' num2str(end_year)];
        clustername = 'gadi';
        md.cluster  = set_cluster(clustername);
        md.settings.waitonlock = 0;
        md.verbose = verbose('solution', true, 'module', true, 'convergence', true);
        md = solve(md, 'tr', 'runtimename', false, 'loadonly', loadonly);
        if loadonly
            savemodel(org, md);
        end
    end % }}}

    % ================================================================= Step 7
    if perform(org, 'VAFContinuityCheck') % {{{
        % Load both projection segments and verify the VAF time series joins
        % continuously at the 2151 restart boundary. Plots sea-level-equivalent
        % contribution relative to 2015, coloured by segment.

        md1 = loadmodel(org, 'ProjRun_2015_2150');
        md2 = loadmodel(org, 'ProjRun_2151_2300');

        time1 = [md1.results.TransientSolution.time];
        time2 = [md2.results.TransientSolution.time];
        vaf1  = [md1.results.TransientSolution.IceVolumeAboveFloatationScaled];
        vaf2  = [md2.results.TransientSolution.IceVolumeAboveFloatationScaled];

        % IceVolumeAboveFloatationScaled is in m³ (ice volume).
        % SLE [m] = -ΔVAF [m³] * rho_ice / rho_sw / A_ocean
        rho_ice = 917;       % kg/m³
        rho_sw  = 1028;      % kg/m³
        A_ocean = 3.625e14;  % m² (ISMIP6 standard ocean area)
        vaf_ref = vaf1(1);
        sle1 = -(vaf1 - vaf_ref) * rho_ice / rho_sw / A_ocean;
        sle2 = -(vaf2 - vaf_ref) * rho_ice / rho_sw / A_ocean;

        gap_m = sle2(1) - sle1(end);
        fprintf('VAF continuity check:\n');
        fprintf('  Segment 1 end   (t=%.1f): %.4f m SLE\n', time1(end), sle1(end));
        fprintf('  Segment 2 start (t=%.1f): %.4f m SLE\n', time2(1),   sle2(1));
        fprintf('  Junction gap: %.6f m SLE\n', gap_m);
        if abs(gap_m) > 0.001
            warning('VAF junction gap > 1 mm SLE (%.6f m) -- check AIS_state_2151 restart.', gap_m);
        end

        fig = figure('visible', 'off');
        plot(time1, sle1, 'b-', 'LineWidth', 1.5); hold on;
        plot(time2, sle2, 'r-', 'LineWidth', 1.5);
        xline(mid_year + 1, 'k--', 'LineWidth', 1, 'Label', '2151 restart', ...
              'LabelVerticalAlignment', 'bottom');
        xlabel('Year');
        ylabel('Sea level contribution (m SLE)');
        title(sprintf('%s %s — VAF continuity check', CMIP_MODEL, SCENARIO), ...
              'Interpreter', 'none');
        legend('2015–2150', '2151–2300', 'Location', 'northwest');
        grid on;

        figdir = [proj_root 'postprocessed_data/figures/' CMIP_MODEL '/' SCENARIO '/'];
        if ~exist(figdir, 'dir'), mkdir(figdir); end
        figname = fullfile(figdir, sprintf('VAF_continuity_%s_%s.png', CMIP_MODEL, SCENARIO));
        saveas(fig, figname);
        close(fig);
        fprintf('Saved: %s\n', figname);
    end % }}}

    % ================================================================= Step 8
    if perform(org, 'WriteISMIP6_NetCDF') % {{{
        % Combines the two TransientSolution arrays (steps 4 + 6), drops the
        % sub-annual safety step from each segment, and grids onto the standard
        % ISMIP7 761×761 8 km AIS grid.  See proj_runs/functions/ for the
        % full implementation.
        addpath('./../../functions');

        md1 = loadmodel(org, 'ProjRun_2015_2150');
        md2 = loadmodel(org, 'ProjRun_2151_2300');
        md = md1;
        md.results.TransientSolution = [md1.results.TransientSolution, ...
                                         md2.results.TransientSolution];
        clear md1 md2;

        outdir = [proj_root 'postprocessed_data/CESM2-WACCM/' SCENARIO '/'];
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

