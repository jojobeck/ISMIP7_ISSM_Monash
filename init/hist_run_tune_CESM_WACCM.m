function md = hist_run_tune_CESM_WACCM(steps, loadonly)
% Historical validation run 1995-2020 using CESM2-WACCM forcing.
% Compares modelled mass loss (WAIS / EAIS / Peninsula) against Otosaka et al.
%
% Forcing strategy:
%   Ocean TF  : CESM2-WACCM historical 1995-2014 (raw_data/.../historical/ocean/tf/v3)
%               No SSP126 ocean TF found -> 2014 TF repeated for 2015-2020 [FLAGGED]
%   SMB       : CESM2-WACCM anomaly + RACMO 1995-2014 climatology
%               historical/SDBN1-2000m/acabf  for 1995-2014
%               ssp126/SDBN1-2000m/acabf       for 2015-2020
%               Same split for dacabfdz gradient
%   Calving   : icemask_greene from AntarcticaObsISMIP7-v1.2.nc (annual, 1997-2021)
%               via md.levelset.spclevelset (level-set calving front)
%
% Step map:
%   1  CESM_TF_hist        build annual TF mat 1995-2020 (2015-2020 = 2014 repeated)
%   2  CESM_SMB_hist       build annual SMB mat (CESM anomaly + RACMO clim, SMBgradients)
%   3  Greene_levelset     build spclevelset time series from icemask_greene
%   4  Relaxed_CESM_WACCM  short 1-year run (1995->1996) under real CESM_WACCM
%                          forcing, used only to relax inputmodel_relax's
%                          geometry/ocean_levelset before HistRun starts
%   5  Assign_Regions                 build & save WAIS/EAIS/Peninsula mask (vertices+elements)
%                                     + modelled mass-change-per-region trend figure
%   6  HistRun             transient 1995-2020 (loadonly=0 submit, =1 gather)
%   7  HistRun_Validation             SMB timeseries + ice mask maps + BMB maps vs obs
%   8  HistRun_Assessment             mass loss by region vs Otosaka
%   9  HistRun_CorrectSMB             rerun transient with region-scaled SMB
%                                     (1+p_calc_corr from step 8, capped +/-25%)
%   10 HistRun_CorrectSMB_Assessment  mass loss by region vs Otosaka, after SMB correction
%
% TODO (step 7): confirm that regions in sectors_8km.nc are 1=WAIS 2=EAIS 3=Peninsula

    if ~exist('loadonly','var'), loadonly = 1; end
    addpath('./scripts');

    org = organizer('repository', './Models/', 'prefix', 'AIS_ISMIP7_', ...
                    'steps', steps, 'color', '34;47;2');
    clear steps;

    % ------------------------------------------------------------------ paths
    inputmodel_relax = './Models/AIS_ISMIP7_Relaxed.mat';

    raw_hist  = './../raw_data/ISMIP7/AIS/CESM2-WACCM/historical';
    raw_ssp   = './../raw_data/ISMIP7/AIS/CESM2-WACCM/ssp126';

    tf_dir        = [raw_hist '/ocean/tf/v3/'];
    smb_hist_dir  = [raw_hist '/SDBN1-2000m/acabf/v2/'];
    grad_hist_dir = [raw_hist '/SDBN1-2000m/dacabfdz/v2/'];
    smb_ssp_dir   = [raw_ssp  '/SDBN1-2000m/acabf/v2/'];
    grad_ssp_dir  = [raw_ssp  '/SDBN1-2000m/dacabfdz/v2/'];

    racmo_clim_nc = './../raw_data/nc_orig/Atmosphere/smb_rec.mean.1995-2014.RACMO2.3p2_ANT27_ERA5-3h.AIS.2km.YY.nc';
    mipkit_nc     = './../raw_data/ISMIP7/AIS/obs/mipkit/AntarcticaObsISMIP7-v1.2.nc';
    sectors_nc    = '/home/565/jb1863/ismip6_2300/masks/sectors_8km.nc';

    preproc_ocean = './../preprocessed_data/Ocean/';
    preproc_atmo  = './../preprocessed_data/Atmosphere/Hist/';
    preproc_clim  = './../preprocessed_data/Atmosphere/Clim/CESM2-WACCM/';
    preproc_front = './../preprocessed_data/Ocean/Hist/';
    for d = {preproc_atmo, preproc_clim, preproc_front}
        if ~exist(d{1}, 'dir'), mkdir(d{1}); end
    end

    % Time range
    start_year  = 1995;
    end_year    = 2020;
    years       = start_year:end_year;
    nyears      = length(years);
    hist_end    = 2014;   % CMIP6 historical boundary
    sec_to_year = 31556926;

    % ================================================================= Step 1
    if perform(org, 'CESM_TF_hist') % {{{
        % Annual TF, 30 depth levels, 8km grid -> interpolate to mesh.
        % 1995-2014: historical files (decade chunks).
        % 2015-2020: no SSP126 TF available -> REPEAT 2014 TF.
        % [FLAG] Replace 2015-2020 with SSP126 ocean TF once available.

        md = loadmodel(inputmodel_relax);
        m  = ((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1) ...
              ./ (1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor = (1./m).^2;

        tf_files = {
            [tf_dir 'tf_AIS_CESM2-WACCM_historical_ocean_v3_1990-1999.nc'], ...
            [tf_dir 'tf_AIS_CESM2-WACCM_historical_ocean_v3_2000-2009.nc'], ...
            [tf_dir 'tf_AIS_CESM2-WACCM_historical_ocean_v3_2010-2014.nc']  ...
        };

        z_data  = double(ncread(tf_files{1}, 'z'));
        nDepths = length(z_data);
        nVerts  = md.mesh.numberofvertices;

        tf_mat = zeros(nVerts, nyears, nDepths);
        t_vec  = years;
        tf2014 = [];

        for f = 1:length(tf_files)
            if ~isfile(tf_files{f}), continue; end
            [~, fname, ~] = fileparts(tf_files{f});
            parts      = strsplit(fname, '_');
            decade     = strsplit(parts{end}, '-');
            file_years = str2double(decade{1}) : str2double(decade{2});
            x_n = double(ncread(tf_files{f}, 'x'));
            y_n = double(ncread(tf_files{f}, 'y'));
            for ti = 1:length(file_years)
                yr = file_years(ti);
                if yr < start_year || yr > hist_end, continue; end
                ki      = yr - start_year + 1;
                tf_data = double(ncread(tf_files{f}, 'tf', [1 1 1 ti], [Inf Inf Inf 1]));
                for i = 1:nDepths
                    v = InterpFromGridToMesh(x_n, y_n, tf_data(:,:,i)', ...
                                            md.mesh.x, md.mesh.y, 0);
                    tf_mat(:, ki, i) = max(v, 0);
                end
                if yr == hist_end
                    tf2014 = squeeze(tf_mat(:, ki, :));
                end
            end
        end

        % Extend 2015-2020 with 2014 TF
        % [FLAG] Replace with SSP126 ocean TF when available
        fprintf('[FLAG] No SSP126 ocean TF: repeating 2014 for %d-%d\n', hist_end+1, end_year);
        for yr = (hist_end+1):end_year
            ki = yr - start_year + 1;
            tf_mat(:, ki, :) = tf2014;
        end

        % Build cell array (3rd dim = depth) with time row appended
        tf_hist = cell(1, 1, nDepths);
        for i = 1:nDepths
            slice = squeeze(tf_mat(:,:,i));
            tf_hist{1,1,i} = [slice ; t_vec];
        end

        save([preproc_front 'CESM_WACCM_TF_' num2str(start_year) '_' num2str(end_year) '.mat'], ...
             'tf_hist', 'z_data', 't_vec', '-v7.3');
        fprintf('Saved TF: %d years, %d depths.\n', nyears, nDepths);
    end % }}}

    % ================================================================= Step 2
    if perform(org, 'CESM_SMB_hist') % {{{
        % SMBgradients units (from ISSM C++ smb_core: smb/1000*rho_w/rho_i):
        %   smbref  must be in mm w.e. yr-1   -> acabf [kg m-2 s-1] * sec_to_year = mm w.e. yr-1
        %   b_pos/b_neg must be mm w.e. yr-1 m-1 -> dacabfdz [kg m-2 s-1 m-1] * sec_to_year
        %   RACMO smb_rec already in mm w.e. yr-1 -> use as-is (no rhoi division)
        %   href is in metres (no scaling)
        %
        % CESM acabf source files are MONTHLY (12 time steps/year, confirmed via
        % ncdump). We collapse each year to a single annual mean (equally-weighted
        % across the 12 months) rather than keeping the real monthly resolution,
        % because the RACMO baseline (smb_racmo) is only a static 1995-2014 ANNUAL
        % climatology -- there is no monthly RACMO climatology to anchor a true
        % seasonal anomaly against. Equal weighting (not day-in-month weighting)
        % matches md.timestepping.time_step = 1/12, which also treats every month
        % as exactly 1/12 year regardless of its real length. Forcing is stored one
        % value per year (t_smb = integer years); with interp_forcing = 1, ISSM
        % linearly interpolates between consecutive annual values for each monthly
        % sub-step, so the seasonal cycle is not resolved -- only the inter-annual
        % trend is.

        md   = loadmodel(inputmodel_relax);
        m    = ((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1) ...
               ./ (1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor = (1./m).^2;
        nVerts = md.mesh.numberofvertices;

        % RACMO climatology [mm w.e. yr-1] — keep as-is
        x_r = double(ncread(racmo_clim_nc, 'x'));
        y_r = double(ncread(racmo_clim_nc, 'y'));
        smb_racmo = InterpFromGridToMesh(x_r, y_r, ...
                        double(ncread(racmo_clim_nc, 'smb_rec'))', ...
                        md.mesh.x, md.mesh.y, 0);   % mm w.e. yr-1

        x_s = []; y_s = [];

        % CESM 1995-2014 mean on mesh [mm w.e. yr-1] for anomaly baseline
        fprintf('Computing CESM SMB mean 1995-%d...\n', hist_end);
        cesm_sum = zeros(nVerts, 1);
        for yr = start_year:hist_end
            nc = smb_nc_path(smb_hist_dir, smb_ssp_dir, hist_end, yr, 'acabf');
            if isempty(x_s)
                x_s = double(ncread(nc, 'x'));
                y_s = double(ncread(nc, 'y'));
            end
            am = mean(double(ncread(nc, 'acabf')), 3);   % kg m-2 s-1
            cesm_sum = cesm_sum + InterpFromGridToMesh(x_s, y_s, am', md.mesh.x, md.mesh.y, 0);
        end
        cesm_mean = cesm_sum / (hist_end - start_year + 1) * sec_to_year;  % mm w.e. yr-1

        % Build annual matrices
        smb_matrix   = zeros(nVerts, nyears);
        bgrad_matrix = zeros(nVerts, nyears);
        t_smb        = years;

        for k = 1:nyears
            yr = years(k);

            nc      = smb_nc_path(smb_hist_dir, smb_ssp_dir, hist_end, yr, 'acabf');
            am      = mean(double(ncread(nc, 'acabf')), 3);   % kg m-2 s-1
            cesm_yr = InterpFromGridToMesh(x_s, y_s, am', md.mesh.x, md.mesh.y, 0) * sec_to_year;  % mm w.e. yr-1
            smb_matrix(:,k) = smb_racmo + (cesm_yr - cesm_mean);   % mm w.e. yr-1

            nc_g  = smb_nc_path(grad_hist_dir, grad_ssp_dir, hist_end, yr, 'dacabfdz');
            g_raw = squeeze(double(ncread(nc_g, 'dacabfdz')));      % kg m-2 s-1 m-1
            bgrad_matrix(:,k) = InterpFromGridToMesh(x_s, y_s, g_raw', md.mesh.x, md.mesh.y, 0) ...
                                 * sec_to_year;                      % mm w.e. yr-1 m-1

            fprintf('  SMB+grad year %d\n', yr);
        end

        smb_forcing   = [smb_matrix  ; t_smb];
        bgrad_forcing = [bgrad_matrix ; t_smb];

        save([preproc_atmo 'CESM_WACCM_SMB_' num2str(start_year) '_' num2str(end_year) '.mat'], ...
             'smb_forcing', 'bgrad_forcing', 't_smb', '-v7.3');
        save([preproc_clim 'CESM_WACCM_SMB_clim_' num2str(start_year) '_' num2str(hist_end) '.mat'], ...
             'smb_racmo', 'cesm_mean', '-v7.3');
        fprintf('Saved SMB forcing and climatology.\n');
    end % }}}

    % ================================================================= Step 3
    if perform(org, 'Greene_levelset') % {{{
        % icemask_greene: int8, 500m grid, (greene_mask_time=24, y=12161, x=12161)
        %   1 = ice -> levelset -1  |  0 = ocean -> levelset +1  |  -128 = fill -> NaN
        % Available from 1997-10; 1997 snapshot used for 1995-1996.

        md     = loadmodel(inputmodel_relax);
        nVerts = md.mesh.numberofvertices;

        x_g = double(ncread(mipkit_nc, 'x'));   % ascending, 500 m
        y_g = double(ncread(mipkit_nc, 'y'));   % descending in file -> flip

        time_raw    = double(ncread(mipkit_nc, 'greene_mask_time'));
        greene_yrs  = 1900 + time_raw / 365.25;

        model_years       = start_year : end_year;
        spclevelset_mat   = NaN(nVerts, length(model_years));

        for k = 1:length(model_years)
            yr = model_years(k);
            [~, gi] = min(abs(greene_yrs - yr));

            raw      = double(ncread(mipkit_nc, 'icemask_greene', [1 1 gi], [Inf Inf 1]));
            raw_flip = flipud(raw');       % (y_asc, x)
            y_asc    = flipud(y_g);

            lset = NaN(size(raw_flip));
            lset(raw_flip == 1) = -1;     % ice
            lset(raw_flip == 0) =  1;     % ocean

            v = InterpFromGridToMesh(x_g, y_asc, lset, md.mesh.x, md.mesh.y, NaN);
            spclevelset_mat(:, k) = v;
            fprintf('  Greene levelset year %d -> snapshot %.1f\n', yr, greene_yrs(gi));
        end

        % Only allow Greene to prescribe ice-front RETREAT relative to the
        % model's initial md.mask.ice_levelset. At vertices that start with
        % no ice (ice_levelset > 0), some Greene years spuriously show ice
        % (advance) -- do not apply those; leave unconstrained (NaN) instead
        % of forcing ice into a region the model never had ice.
        no_ice0      = md.mask.ice_levelset > 0;                  % nVerts x 1
        advance_mask = repmat(no_ice0, 1, length(model_years)) & (spclevelset_mat < 0);
        fprintf('  Suppressing %d/%d Greene advance entries (no ice -> ice) relative to initial mask.\n', ...
                nnz(advance_mask), numel(advance_mask));
        spclevelset_mat(advance_mask) = NaN;

        greene_spclevelset = [spclevelset_mat ; model_years];
        save([preproc_front 'Greene_spclevelset_' num2str(start_year) '_' num2str(end_year) '.mat'], ...
             'greene_spclevelset', '-v7.3');
        fprintf('Saved Greene levelset mat.\n');
    end % }}}

    % ================================================================= Step 4
    if perform(org, 'Relaxed_CESM_WACCM') % {{{
        % Short 1-year run (1995 -> 1996) using the real CESM_WACCM ocean/SMB/
        % calving forcing, used only to relax the static inputmodel_relax
        % geometry into one that is dynamically consistent with this forcing.
        % inputmodel_relax was relaxed under different (steady/inversion)
        % conditions, so switching on the real time-varying forcing at the
        % start of HistRun produces a sharp initial jump; absorbing that
        % adjustment here first, then starting HistRun from the result,
        % removes it from the recorded 1995-2020 series.
        %
        % Mirrors tuning_func.m's 'Relaxed' step: take the final state of a
        % transient and fold it back into geometry/mask as the new initial
        % condition (Base, Thickness, Surface, ocean_levelset, ice_levelset).
        % All of these are folded back together so the new initial state
        % stays internally consistent -- updating ice_levelset on its own,
        % without the matching geometry/ocean_levelset, is what introduces
        % inconsistency, not updating it as part of this joint fold-back.

        md = loadmodel(inputmodel_relax);
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
        md.timestepping.final_time = start_year + 1;   % single year: 1995 -> 1996
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;   % one snapshot, at the end

        md.transient.requested_outputs = { ...
            'default', ...
            'IceVolume', 'IceVolumeAboveFloatation', ...
            'GroundedArea', 'FloatingArea', ...
            'TotalSmb', 'TotalGroundedBmb', 'TotalFloatingBmb', ...
            'IceVolumeAboveFloatationScaled', ...
            'Thickness', ...
            'MaskOceanLevelset', 'MaskIceLevelset', ...
            'SmbMassBalance', ...
            'BasalforcingsFloatingiceMeltingRate'};

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean (same as HistRun) ---
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_front 'CESM_WACCM_TF_' num2str(start_year) '_' num2str(end_year) '.mat']);
        load([preproc_ocean 'gamma0_local.mat']);

        unique_basinid = unique(basinid);
        delta_t        = zeros(1, length(unique_basinid));

        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique_basinid);
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = tf_hist;
        md.basalforcings.islocal    = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_local;

        % --- SMB (same as HistRun) ---
        load([preproc_atmo 'CESM_WACCM_SMB_' num2str(start_year) '_' num2str(end_year) '.mat']);
        md.smb        = SMBgradients();
        md.smb.smbref = smb_forcing;
        md.smb.b_pos  = bgrad_forcing;
        md.smb.b_neg  = bgrad_forcing;
        md.smb.href   = [md.geometry.surface ; start_year];

        % --- Calving front (same as HistRun) ---
        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 1;
        load([preproc_front 'Greene_spclevelset_' num2str(start_year) '_' num2str(end_year) '.mat']);
        md.levelset.spclevelset = greene_spclevelset;

        md.miscellaneous.name = 'Relaxed_CESM_WACCM';
        clustername = 'gadi';
        md.cluster  = set_cluster(clustername);
        md.settings.waitonlock = 0;
        md.verbose = verbose('solution', true, 'module', true, 'convergence', false);
        md = solve(md, 'tr', 'runtimename', false, 'loadonly', loadonly);

        if loadonly
            md_in=md;
            base = md_in.results.TransientSolution(end).Base;
            thickness = md_in.results.TransientSolution(end).Thickness;
            surface = md_in.results.TransientSolution(end).Surface;
            md.geometry.thickness = thickness;
            md.geometry.surface = surface;
            md.geometry.base = base;
            md.mask.ocean_levelset= md_in.results.TransientSolution(end).MaskOceanLevelset; 
            md.mask.ice_levelset= md_in.results.TransientSolution(end).MaskIceLevelset; 
            md.results.TransientSolution = [];
            savemodel(org,md);

        end
    end % }}}

    % ================================================================= Step 5
    if perform(org, 'Assign_Regions') % {{{
        % Build & save the WAIS / EAIS / Peninsula macro-region mask
        % (vertices + elements), analogous to tuning_func.m's 'Assign_Basins'
        % step (which builds the IMBIE2 basin-on-elements/vertices masks
        % already loaded in steps 4/8 from Imbie2_extrap_2km_BasinOnElements.mat).
        % Saved once here and reused by HistRun_Assessment / HistRun_CorrectSMB
        % instead of re-interpolating sectors_8km.nc every time.
        %   regions = 1 -> WAIS  |  2 -> EAIS  |  3 -> Peninsula
        % TODO: confirm these region values with the file on Gadi before trusting output.

        md = loadmodel(org, 'HistRun');

        x_sec    = double(ncread(sectors_nc, 'x'));
        y_sec    = double(ncread(sectors_nc, 'y'));
        reg_grid = double(ncread(sectors_nc, 'regions'));
        reg_vert = round(InterpFromGridToMesh(x_sec, y_sec, reg_grid', md.mesh.x, md.mesh.y, 0));
        reg_elem = mode(reg_vert(md.mesh.elements), 2);

        region_names = {'WAIS', 'EAIS', 'Peninsula'};
        region_ids   = [1, 2, 3];

        region_mask_file = [preproc_ocean 'Basins/HistRun_Regions.mat'];
        if ~exist([preproc_ocean 'Basins'], 'dir'), mkdir([preproc_ocean 'Basins']); end
        save(region_mask_file, 'reg_vert', 'reg_elem', 'region_names', 'region_ids', '-v7.3');
        fprintf('Saved: %s\n', region_mask_file);

        % --- Diagnostic: modelled cumulative mass change per region, with fitted trend ---
        rhoi = md.materials.rho_ice;
        nR   = length(region_ids);
        nT   = length(md.results.TransientSolution);
        time = zeros(1, nT);
        dMass = zeros(nR, nT);

        x1 = md.mesh.x(md.mesh.elements(:,1));  y1 = md.mesh.y(md.mesh.elements(:,1));
        x2 = md.mesh.x(md.mesh.elements(:,2));  y2 = md.mesh.y(md.mesh.elements(:,2));
        x3 = md.mesh.x(md.mesh.elements(:,3));  y3 = md.mesh.y(md.mesh.elements(:,3));
        elem_area = 0.5 * abs((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1));

        H0_vert = md.results.TransientSolution(1).Thickness;
        for t = 1:nT
            time(t) = md.results.TransientSolution(t).time;
            H_vert  = md.results.TransientSolution(t).Thickness;
            dH_vert = H_vert - H0_vert;
            dH_elem = mean(dH_vert(md.mesh.elements), 2);
            for r = 1:nR
                mask = (reg_elem == region_ids(r));
                dMass(r, t) = sum(dH_elem(mask) .* elem_area(mask)) * rhoi / 1e12;  % Gt
            end
        end

        figure('visible','off');
        cols = {'b','r','g'};
        hold on;
        for r = 1:nR
            plot(time, dMass(r,:), cols{r}, 'LineWidth', 1.5, 'DisplayName', region_names{r});
            p = polyfit(time, dMass(r,:), 1);
            plot(time, polyval(p, time), [cols{r} '--'], 'LineWidth', 1, ...
                 'DisplayName', sprintf('%s fit (%.1f Gt/yr)', region_names{r}, p(1)));
        end
        plot(time, sum(dMass,1), 'k-', 'LineWidth', 1.5, 'DisplayName', 'AIS total');
        xlabel('Year'); ylabel('\DeltaMass (Gt)');
        title('Modelled cumulative mass change per region with fitted trend');
        legend('Location','southwest'); grid on;
        if ~exist('./figures','dir'), mkdir('./figures'); end
        saveas(gcf, './figures/Assign_Regions_mass_trend.png');
        fprintf('Saved: figures/Assign_Regions_mass_trend.png\n');
    end % }}}
    % ================================================================= Step 6
    if perform(org, 'HistRun') % {{{

        md = loadmodel(org, 'Relaxed_CESM_WACCM');
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
        md.timestepping.final_time = end_year + 1;
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;   % annual snapshots

        md.transient.requested_outputs = { ...
            'default', ...
            'IceVolume', 'IceVolumeAboveFloatation', ...
            'GroundedArea', 'FloatingArea', ...
            'TotalSmb', 'TotalGroundedBmb', 'TotalFloatingBmb', ...
            'IceVolumeAboveFloatationScaled', ...
            'Thickness', ...
            'MaskOceanLevelset', 'MaskIceLevelset', ...
            'SmbMassBalance', ...
            'BasalforcingsFloatingiceMeltingRate'};

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean ---
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_front 'CESM_WACCM_TF_' num2str(start_year) '_' num2str(end_year) '.mat']);
        load([preproc_ocean 'gamma0_local.mat']);

        unique_basinid = unique(basinid);
        delta_t        = zeros(1, length(unique_basinid));

        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique_basinid);
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = tf_hist;
        md.basalforcings.islocal    = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_local;

        % --- SMB ---
        load([preproc_atmo 'CESM_WACCM_SMB_' num2str(start_year) '_' num2str(end_year) '.mat']);
        md.smb        = SMBgradients();
        md.smb.smbref = smb_forcing;           % RACMO clim + CESM anomaly [m ice yr-1], timeseries
        md.smb.b_pos  = bgrad_forcing;         % dacabfdz > 0 in accumulation zone [m ice yr-1 m-1]
        md.smb.b_neg  = bgrad_forcing;         % dacabfdz > 0 in ablation zone too (SMB always increases with elev.)
        md.smb.href   = [md.geometry.surface ; start_year];  % static reference surface [m], single snapshot

        % --- Calving front ---
        % ismovingfront = 1: activates level-set solver so spclevelset is applied.
        % Without it the level-set equation is never solved and the front is frozen.

        % calving.calvingrate defaults to NaN and is only checked once ismovingfront=1.
        % Front position is fully prescribed by spclevelset, so no extra calving needed.
        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront =1;
        load([preproc_front 'Greene_spclevelset_' num2str(start_year) '_' num2str(end_year) '.mat']);
        md.levelset.spclevelset = greene_spclevelset;

        md.miscellaneous.name = ['HistRun_CESM_WACCM_' num2str(start_year) '_' num2str(end_year)];
        clustername = 'gadi';
        md.cluster  = set_cluster(clustername);
        md.settings.waitonlock = 0;
        md.verbose = verbose('solution', true, 'module', true, 'convergence',true);
        md = solve(md, 'tr', 'runtimename', false, 'loadonly', loadonly);
        if loadonly
            savemodel(org, md);
        end
    end % }}}

    % ================================================================= Step 7
    if perform(org, 'HistRun_Validation') % {{{
        % Three validation figures before Otosaka comparison:
        %   Fig 1: Total AIS SMB over time  (model TotalSmb vs integrated smb_forcing)
        %   Fig 2: Ice mask maps at 1996/2000/2010 vs Greene icemask_greene
        %   Fig 3: BMB melt maps at 1996/2000/2010 vs Paolo/Adusumilli (static obs)

        md = loadmodel(org, 'HistRun');
        rhoi = md.materials.rho_ice;
        if ~exist('./figures', 'dir'), mkdir('./figures'); end

        % --- Common: time vector and element/vertex areas ---
        nT   = length(md.results.TransientSolution);
        time = zeros(1, nT);
        for t = 1:nT, time(t) = md.results.TransientSolution(t).time; end

        x1 = md.mesh.x(md.mesh.elements(:,1)); y1 = md.mesh.y(md.mesh.elements(:,1));
        x2 = md.mesh.x(md.mesh.elements(:,2)); y2 = md.mesh.y(md.mesh.elements(:,2));
        x3 = md.mesh.x(md.mesh.elements(:,3)); y3 = md.mesh.y(md.mesh.elements(:,3));
        elem_area = 0.5 * abs((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1));
        vert_area = accumarray(md.mesh.elements(:), repmat(elem_area/3, 3, 1), ...
                               [md.mesh.numberofvertices, 1]);

        % Helper: find nearest timestep index for a target year
        nearest_t = @(yr) find(abs(time - yr) == min(abs(time - yr)), 1);

        % ---- Figure 1: SMB time series ----
        % TotalSmb is already in Gt/yr
        total_smb_model = zeros(1, nT);
        for t = 1:nT
            total_smb_model(t) = md.results.TransientSolution(t).TotalSmb;
        end
        total_smb_model_Gt = total_smb_model;

        % Forcing: smb_forcing [mm w.e./yr] integrated over mesh -> Gt/yr
        % 1 mm w.e./yr * 1 m2 = 1e-3 m w.e./yr * m2 = 1e-3 m3 w.e./yr -> *1000 kg/m3 -> kg/yr -> /1e12 Gt/yr
        % Simplifies to: sum(smb_mm * area) / 1e12
        load([preproc_atmo 'CESM_WACCM_SMB_' num2str(start_year) '_' num2str(end_year) '.mat']);
        total_smb_forcing_Gt = zeros(1, nyears);
        for k = 1:nyears
            smb_v = smb_forcing(1:end-1, k);   % mm w.e./yr at each vertex
            total_smb_forcing_Gt(k) = sum(smb_v .* vert_area) / 1e12;
        end

        figure('visible','off','Position',[0 0 900 400]);
        plot(time, total_smb_model_Gt, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Model TotalSmb');
        hold on;
        plot(years, total_smb_forcing_Gt, 'r--o', 'LineWidth', 1.5, 'MarkerSize', 4, ...
             'DisplayName', 'Integrated smb\_forcing');
        xlabel('Year'); ylabel('Total SMB (Gt yr^{-1})');
        title('AIS total SMB: model output vs. CESM+RACMO forcing');
        legend('Location','best'); grid on;
        saveas(gcf, './figures/Validation_SMB_timeseries.png');
        close;
        fprintf('Saved: Validation_SMB_timeseries.png\n');

        % ---- Figure 2: Calving front evolution – PIG & Thwaites zoom ----
        % Background via plotmodel (native mesh patch); fronts via isoline (exact
        % mesh-edge contour, no raster interpolation). Greene levelset is the
        % mesh-interpolated field from step 3, so isoline works on it directly too.
        % Union of PIG ('xlim',[-1.7e6,-1.5e6],'ylim',[-3.8e5,-2e5]) and
        % THW ('xlim',[-1.62e6,-1.48e6],'ylim',[-4.9e5,-3.9e5]) boxes.
        xl_cf = [-1.7e6, -1.48e6];
        yl_cf = [-4.9e5, -2e5];

        % plotmodel's default 'axis tight equal' is computed from the FULL mesh
        % extent before our xlim/ylim are applied; with an equal aspect ratio and
        % a figure box that doesn't match the zoom box's shape, MATLAB stretches
        % xlim back out to fill the figure (showing the whole continent). Size the
        % zoom figures to match the zoom box's aspect ratio so no stretch is needed.
        fig_h_cf = 800;
        fig_w_cf = round(fig_h_cf * diff(xl_cf) / diff(yl_cf));

        load([preproc_front 'Greene_spclevelset_' num2str(start_year) '_' num2str(end_year) '.mat']);
        greene_model_yrs = greene_spclevelset(end, :);      % 1995:2020
        greene_lset_mat  = greene_spclevelset(1:end-1, :);  % nVerts x nyears

        cf_years = unique([start_year+1, 2000, 2005, 2010, 2015, end_year]);
        n_cf     = length(cf_years);
        cmap_cf  = jet(n_cf);

        % Compute calving-front isolines once, reuse for both the full-AIS
        % and PIG/Thwaites-zoom figures below.
        hi_mod_all = cell(n_cf, 1);
        hi_gre_all = cell(n_cf, 1);
        for k = 1:n_cf
            yr = cf_years(k);

            ti = nearest_t(yr);
            hi_mod_all{k} = isoline(md, md.results.TransientSolution(ti).MaskIceLevelset, ...
                                     'value', 0, 'output', 'matrix');

            [~, gk] = min(abs(greene_model_yrs - yr));
            hi_gre_all{k} = isoline(md, double(greene_lset_mat(:, gk)), ...
                                     'value', 0, 'output', 'matrix');
        end

        % ---- Figure 2a: full Antarctic Ice Sheet ----
        figure('visible','off','Position',[0 0 900 750]);
        plotmodel(md, 'figure', gcf, 'visible', 'off', 'data', md.mask.ice_levelset, ...
                  'colormap', gray, 'caxis', [-1 1], 'title', '', 'colorbar', 0);
        hold on;
        for k = 1:n_cf
            yr  = cf_years(k);
            col = cmap_cf(k, :);
            plot(hi_mod_all{k}(:,1), hi_mod_all{k}(:,2), '-', 'Color', col, 'LineWidth', 2, ...
                 'DisplayName', sprintf('%.0f model', yr));
            plot(hi_gre_all{k}(:,1), hi_gre_all{k}(:,2), '--', 'Color', col, 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('%.0f Greene', yr));
        end
        axis equal tight off;
        legend('Location', 'southeast', 'FontSize', 8);
        title(sprintf('Calving front %.0f–%.0f  |  Antarctica', ...
                      cf_years(1), cf_years(end)));
        saveas(gcf, './figures/Validation_CalvingFront_AIS.png');
        close;
        fprintf('Saved: Validation_CalvingFront_AIS.png\n');

        % ---- Figure 2b: PIG & Thwaites zoom ----
        figure('visible','off','Position',[0 0 fig_w_cf fig_h_cf]);
        plotmodel(md, 'figure', gcf, 'visible', 'off', 'data', md.mask.ice_levelset, ...
                  'colormap', gray, 'caxis', [-1 1], ...
                  'xlim', xl_cf, 'ylim', yl_cf, 'title', '', 'colorbar', 0);
        hold on;
        for k = 1:n_cf
            yr  = cf_years(k);
            col = cmap_cf(k, :);
            plot(hi_mod_all{k}(:,1), hi_mod_all{k}(:,2), '-', 'Color', col, 'LineWidth', 2, ...
                 'DisplayName', sprintf('%.0f model', yr));
            plot(hi_gre_all{k}(:,1), hi_gre_all{k}(:,2), '--', 'Color', col, 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('%.0f Greene', yr));
        end
        axis equal off;
        xlim(xl_cf); ylim(yl_cf);
        legend('Location', 'southeast', 'FontSize', 8);
        title(sprintf('Calving front %.0f–%.0f  |  PIG & Thwaites', ...
                      cf_years(1), cf_years(end)));
        saveas(gcf, './figures/Validation_CalvingFront_PIGTHW.png');
        close;
        fprintf('Saved: Validation_CalvingFront_PIGTHW.png\n');

        % ---- Figure 2c: initial ice/ocean mask check – PIG & Thwaites zoom ----
        % Diagnostic: overlay the zero contours of the initial md.mask.ice_levelset
        % (ice front) and md.mask.ocean_levelset (grounding line) to check for
        % inconsistencies between the two masks, e.g. near Thwaites.
        hi_ice0   = isoline(md, md.mask.ice_levelset,   'value', 0, 'output', 'matrix');
        hi_ocean0 = isoline(md, md.mask.ocean_levelset, 'value', 0, 'output', 'matrix');

        figure('visible','off','Position',[0 0 fig_w_cf fig_h_cf]);
        plotmodel(md, 'figure', gcf, 'visible', 'off', 'data', md.mask.ice_levelset, ...
                  'colormap', gray, 'caxis', [-1 1], ...
                  'xlim', xl_cf, 'ylim', yl_cf, 'title', '', 'colorbar', 0);
        hold on;
        plot(hi_ice0(:,1), hi_ice0(:,2), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2, ...
             'DisplayName', 'ice\_levelset (ice front)');
        plot(hi_ocean0(:,1), hi_ocean0(:,2), '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, ...
             'DisplayName', 'ocean\_levelset (grounding line)');
        axis equal off;
        xlim(xl_cf); ylim(yl_cf);
        legend('Location', 'southeast', 'FontSize', 8);
        title('Initial mask check  |  PIG & Thwaites');
        saveas(gcf, './figures/Validation_MaskCheck_PIGTHW.png');
        close;
        fprintf('Saved: Validation_MaskCheck_PIGTHW.png\n');

        % ---- Figure 3: BMB melt maps vs Paolo/Adusumilli ----
        plot_years = [start_year+1, 2000, 2010];
        paolo_nc  = ['./../raw_data/ISMIP7/AIS/parameterisations/ocean/meltobs/' ...
                     'melt_paolo_err_adusumilli_ismip8km.nc'];
        x_p   = double(ncread(paolo_nc, 'x'));
        y_p   = double(ncread(paolo_nc, 'y'));
        paolo = double(ncread(paolo_nc, 'melt_mean'));  % kg/m2/yr

        bmb_clim = [0, 50];   % colour scale kg/m2/yr (floating ice only)

        figure('visible','off','Position',[0 0 1200 700]);
        for col = 1:3
            yr = plot_years(col);
            ti = nearest_t(yr);

            % Model BMB [m/yr ice] -> kg/m2/yr: * rhoi
            bmb_model = md.results.TransientSolution(ti).BasalforcingsFloatingiceMeltingRate;
            bmb_model_kgm2 = bmb_model * rhoi;
            % Mask to floating only
            ocean_lset = md.results.TransientSolution(ti).MaskOceanLevelset;
            bmb_model_kgm2(ocean_lset >= 0) = NaN;   % only floating (ocean_lset < 0)

            subplot(2, 3, col);
            patch('Faces', md.mesh.elements, ...
                  'Vertices', [md.mesh.x md.mesh.y], ...
                  'FaceVertexCData', bmb_model_kgm2, ...
                  'FaceColor', 'flat', 'EdgeColor', 'none');
            axis equal off; colormap(gca, parula); caxis(bmb_clim);
            title(sprintf('Model BMB %.0f', time(ti)));
            colorbar;

            % Paolo obs (static, same for all columns)
            subplot(2, 3, col+3);
            imagesc(x_p, y_p, paolo');
            set(gca, 'YDir', 'normal');
            axis equal tight off; colormap(gca, parula); caxis(bmb_clim);
            title('Paolo/Adusumilli obs (mean)');
            colorbar;
        end
        sgtitle('Basal melt (kg m^{-2} yr^{-1}): model (top) vs Paolo/Adusumilli (bottom)');
        saveas(gcf, './figures/Validation_BMB.png');
        close;
        fprintf('Saved: Validation_BMB.png\n');

    end % }}}


    % ================================================================= Step 8
    if perform(org, 'HistRun_Assessment') % {{{
        md = loadmodel(org, 'HistRun');
        region_mask_file = [preproc_ocean 'Basins/HistRun_Regions.mat'];
        assess_vs_otosaka(md, region_mask_file, start_year, ...
            './Data/Tables/HistRun_otosaka_assessment_CESM_WACCM.csv', ...
            './figures/HistRun_regional_mass_change.png', ...
            'Cumulative mass change vs Otosaka et al.  (CESM2-WACCM hist+ssp126)');
    end % }}}

    % ================================================================= Step 9
    if perform(org, 'HistRun_CorrectSMB') % {{{
        % Rerun the full 1995-2020 transient with the RACMO climatology
        % component of the SMB forcing scaled per region by 1+p_calc_corr
        % from step 8's Otosaka assessment. The CESM year-to-year anomaly is
        % left untouched -- only the static climatological level is corrected.
        % Mirrors tuning_func_justine.m's
        % 'submit_from_Collapserelax20y_correctSMB_withC0.5'.
        %
        % p_calc_corr can be poorly constrained (especially for the
        % Peninsula, which is hard to get right), so the scale factor is
        % capped to +/-25% (i.e. clamped to [0.75, 1.25]) for every region.

        md = loadmodel(org, 'Relaxed_CESM_WACCM');
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
        md.timestepping.final_time = end_year + 1;
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;   % annual snapshots

        md.transient.requested_outputs = { ...
            'default', ...
            'IceVolume', 'IceVolumeAboveFloatation', ...
            'GroundedArea', 'FloatingArea', ...
            'TotalSmb', 'TotalGroundedBmb', 'TotalFloatingBmb', ...
            'IceVolumeAboveFloatationScaled', ...
            'Thickness', ...
            'MaskOceanLevelset', 'MaskIceLevelset', ...
            'SmbMassBalance', ...
            'BasalforcingsFloatingiceMeltingRate'};

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean (unchanged from step 5) ---
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_front 'CESM_WACCM_TF_' num2str(start_year) '_' num2str(end_year) '.mat']);
        load([preproc_ocean 'gamma0_local.mat']);

        unique_basinid = unique(basinid);
        delta_t        = zeros(1, length(unique_basinid));

        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique_basinid);
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = tf_hist;
        md.basalforcings.islocal    = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_local;

        % --- SMB: scale RACMO climatology per region, cap correction to +/-25% ---
        load([preproc_atmo 'CESM_WACCM_SMB_' num2str(start_year) '_' num2str(end_year) '.mat']);
        smb_matrix = smb_forcing(1:end-1, :);   % mm w.e. yr-1, [nVerts x nyears]

        load([preproc_clim 'CESM_WACCM_SMB_clim_' num2str(start_year) '_' num2str(hist_end) '.mat'], 'smb_racmo');

        load([preproc_ocean 'Basins/HistRun_Regions.mat'], 'reg_vert', 'region_names', 'region_ids');

        otosaka_assessment = readtable('./Data/Tables/HistRun_otosaka_assessment_CESM_WACCM.csv');
        p_vert = ones(md.mesh.numberofvertices, 1);
        for r = 1:length(region_ids)
            j     = find(strcmpi(otosaka_assessment.model_region, region_names{r}));
            p_raw = 1 + otosaka_assessment.p_calc_corr(j);
            p_cap = min(max(p_raw, 0.75), 1.25);   % cap correction to +/-25%
            fprintf('  %s: p_calc_corr=%.4f -> scale=%.4f (capped %.4f)\n', ...
                    region_names{r}, otosaka_assessment.p_calc_corr(j), p_raw, p_cap);
            p_vert(reg_vert == region_ids(r)) = p_cap;
        end

        smb_corr_matrix  = p_vert .* smb_racmo + (smb_matrix - smb_racmo);
        smb_forcing_corr = [smb_corr_matrix ; smb_forcing(end, :)];

        % Save corrected forcing so hist_run_CESM_WACCM_1995_2014 and the
        % projection scripts can load the same bias-corrected SMB without
        % re-deriving p_vert. smb_forcing is overwritten here with the
        % corrected version so the saved variable name matches the uncorrected
        % mat -- downstream scripts only need a filename change, not a
        % variable rename. p_vert is also saved so the projection script can
        % load it to apply the same correction to the SSP anomaly series.
        smb_forcing = smb_forcing_corr;
        save([preproc_atmo 'CESM_WACCM_SMB_corrected_' num2str(start_year) '_' num2str(end_year) '.mat'], ...
             'smb_forcing', 'bgrad_forcing', 't_smb', 'p_vert', '-v7.3');
        fprintf('Saved corrected SMB: %s\n', [preproc_atmo 'CESM_WACCM_SMB_corrected_' num2str(start_year) '_' num2str(end_year) '.mat']);

        md.smb        = SMBgradients();
        md.smb.smbref = smb_forcing_corr;      % region-corrected RACMO clim + CESM anomaly
        md.smb.b_pos  = bgrad_forcing;
        md.smb.b_neg  = bgrad_forcing;
        md.smb.href   = [md.geometry.surface ; start_year];

        % --- Calving front (unchanged from step 5) ---
        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 1;
        load([preproc_front 'Greene_spclevelset_' num2str(start_year) '_' num2str(end_year) '.mat']);
        md.levelset.spclevelset = greene_spclevelset;

        md.miscellaneous.name = ['HistRun_CorrectSMB_CESM_WACCM_' num2str(start_year) '_' num2str(end_year)];
        clustername = 'gadi';
        md.cluster  = set_cluster(clustername);
        md.settings.waitonlock = 0;
        md.verbose = verbose('solution', true, 'module', true, 'convergence',true);
        md = solve(md, 'tr', 'runtimename', false, 'loadonly', loadonly);
        if loadonly
            savemodel(org, md);
        end
    end % }}}

    % ================================================================= Step 10
    if perform(org, 'HistRun_CorrectSMB_Assessment') % {{{
        md = loadmodel(org, 'HistRun_CorrectSMB');
        region_mask_file = [preproc_ocean 'Basins/HistRun_Regions.mat'];
        assess_vs_otosaka(md, region_mask_file, start_year, ...
            './Data/Tables/HistRun_otosaka_assessment_CESM_WACCM_corrected.csv', ...
            './figures/HistRun_regional_mass_change_corrected.png', ...
            'Cumulative mass change vs Otosaka et al.  (corrected SMB, CESM2-WACCM hist+ssp126)');
    end % }}}

end

% ------------------------------------------------------------------ helpers
function nc = smb_nc_path(hist_dir, ssp_dir, hist_end, yr, var)
    if yr <= hist_end
        nc = [hist_dir sprintf('%s_AIS_CESM2-WACCM_historical_SDBN1-2000m_v2_%d.nc', var, yr)];
    else
        nc = [ssp_dir  sprintf('%s_AIS_CESM2-WACCM_ssp126_SDBN1-2000m_v2_%d.nc',    var, yr)];
    end
end

function assess_vs_otosaka(md, region_mask_file, start_year, csv_out, png_out, plot_title)
    % Regional mass loss using the region mask saved by 'Assign_Regions':
    %   regions = 1 -> WAIS  |  2 -> EAIS  |  3 -> Peninsula
    %
    % Method: dMass_region(t) = rho_ice * sum_elem[ dH_elem(t) * area_elem ]
    %   where sum is over elements whose centroid region == target.
    %
    % Trend vs Otosaka et al. 1995-2020 (Gt/yr) + SMB correction factor:
    % the model trend (slope of cumulative mass change) is compared against
    % Otosaka's min/mean/max trend per region, and p_calc_corr gives the
    % fractional SMB scaling that would close the gap between the model
    % trend and Otosaka's mean trend over the run, given the region's
    % actual cumulative SMB-driven mass.
    rhoi = md.materials.rho_ice;

    load(region_mask_file, 'reg_vert', 'reg_elem', 'region_names', 'region_ids');

    % --- Element areas (triangles) ---
    x1 = md.mesh.x(md.mesh.elements(:,1));  y1 = md.mesh.y(md.mesh.elements(:,1));
    x2 = md.mesh.x(md.mesh.elements(:,2));  y2 = md.mesh.y(md.mesh.elements(:,2));
    x3 = md.mesh.x(md.mesh.elements(:,3));  y3 = md.mesh.y(md.mesh.elements(:,3));
    elem_area = 0.5 * abs((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1));

    nR            = length(region_ids);
    nT            = length(md.results.TransientSolution);
    time          = zeros(1, nT);
    dMass         = zeros(nR, nT);   % Gt

    H0_vert = md.results.TransientSolution(1).Thickness;

    % --- Vertex areas (nodal, 1/3 of each adjacent triangle) for SMB integration ---
    vert_area = accumarray(md.mesh.elements(:), repmat(elem_area/3, 3, 1), ...
                           [md.mesh.numberofvertices, 1]);

    smb_rate = zeros(nR, nT);   % SMB-driven mass rate per region, Gt/yr

    for t = 1:nT
        time(t) = md.results.TransientSolution(t).time;
        H_vert  = md.results.TransientSolution(t).Thickness;
        dH_vert = H_vert - H0_vert;
        % Average dH to elements
        dH_elem = mean(dH_vert(md.mesh.elements), 2);

        smb_vert = md.results.TransientSolution(t).SmbMassBalance;  % m ice eq/yr
        for r = 1:nR
            mask = (reg_elem == region_ids(r));
            dMass(r, t) = sum(dH_elem(mask) .* elem_area(mask)) * rhoi / 1e12;  % Gt

            vmask = (reg_vert == region_ids(r));
            smb_rate(r, t) = sum(smb_vert(vmask) .* vert_area(vmask)) * rhoi / 1e12;  % Gt/yr
        end
    end

    % --- Print table ---
    fprintf('\n=== Cumulative mass change (Gt) relative to %d ===\n', start_year);
    fprintf('%-6s', 'Year');
    for r = 1:nR, fprintf('%12s', region_names{r}); end
    fprintf('%12s\n', 'AIS total');
    for t = 1:12:nT
        fprintf('%-6.0f', time(t));
        for r = 1:nR, fprintf('%12.1f', dMass(r,t)); end
        fprintf('%12.1f\n', sum(dMass(:,t)));
    end

    % --- Trend vs Otosaka et al. 1995-2020 (Gt/yr) + SMB correction factor ---
    otosaka_csv  = './Data/Tables/Otosaka_mass_change_GTy_min_mean_max_1995_2020.csv';
    otosaka_name = {'West', 'East', 'Peninsula'};   % matches region_names order

    opts = detectImportOptions(otosaka_csv);
    opts = setvartype(opts, 'region', 'string');
    n    = readtable(otosaka_csv, opts);
    n.model_region = strings(height(n), 1);
    n.issm         = NaN(height(n), 1);
    n.in_range     = NaN(height(n), 1);
    n.tot_smb      = NaN(height(n), 1);
    n.p_calc_corr  = NaN(height(n), 1);

    nyrs_run = time(end) - time(1) + 1;   % 26 for 1995-2020

    issm_fit = zeros(nR, 2);   % [slope, intercept] per region, for plotting

    for r = 1:nR
        p         = polyfit(time, dMass(r,:), 1);
        rate_issm = p(1);   % model trend, Gt/yr
        issm_fit(r,:) = p;

        j = find(strcmpi(n.region, otosaka_name{r}));

        n.model_region(j) = region_names{r};
        n.issm(j)         = rate_issm;
        n.in_range(j)      = double(rate_issm >= n.min(j) & rate_issm <= n.max(j));

        y0smb           = cumsum(smb_rate(r,:));
        n.tot_smb(j)      = y0smb(end) - y0smb(1);   % Gt, relative to first year
        n.p_calc_corr(j)  = (n.mean(j)*nyrs_run - n.issm(j)*nyrs_run) / n.tot_smb(j);
    end

    fprintf('\n=== Trend vs Otosaka et al. 1995-2020 (Gt/yr) ===\n');
    fprintf('%-10s%10s%10s%10s%10s%8s%12s%12s\n', ...
            'Region', 'Min', 'Mean', 'Max', 'ISSM', 'InRng', 'TotSMB(Gt)', 'p_corr');
    for j = 1:height(n)
        fprintf('%-10s%10.2f%10.2f%10.2f%10.2f%8d%12.1f%12.4f\n', ...
                n.model_region(j), n.min(j), n.mean(j), n.max(j), n.issm(j), ...
                n.in_range(j), n.tot_smb(j), n.p_calc_corr(j));
    end

    [csv_dir, ~, ~] = fileparts(csv_out);
    if ~isempty(csv_dir) && ~exist(csv_dir, 'dir'), mkdir(csv_dir); end
    writetable(n, csv_out);
    fprintf('Saved: %s\n', csv_out);

    % --- Plot: cumulative mass change vs Otosaka mean-trend reference lines ---
    figure('visible','off');
    cols = {'b','r','g'};
    hold on;
    for r = 1:nR
        plot(time, dMass(r,:), cols{r}, 'LineWidth', 1.5, 'DisplayName', region_names{r});
        j = find(strcmpi(n.region, otosaka_name{r}));
        otosaka_lin = n.mean(j) * (time - time(1));
        plot(time, otosaka_lin, [cols{r} '--'], 'LineWidth', 1, ...
             'DisplayName', sprintf('%s Otosaka mean', region_names{r}));
        issm_lin = polyval(issm_fit(r,:), time);   % true regression fit, not anchored to 0
        plot(time, issm_lin, [cols{r} ':'], 'LineWidth', 1.5, ...
             'DisplayName', sprintf('%s ISSM fit (%.1f Gt/yr)', region_names{r}, issm_fit(r,1)));
    end
    plot(time, sum(dMass,1), 'k--', 'LineWidth', 2, 'DisplayName', 'AIS total');
    xlabel('Year'); ylabel('\DeltaMass (Gt)');
    title(plot_title);
    legend('Location','southwest'); grid on;
    [png_dir, ~, ~] = fileparts(png_out);
    if ~isempty(png_dir) && ~exist(png_dir, 'dir'), mkdir(png_dir); end
    saveas(gcf, png_out);
    fprintf('Saved: %s\n', png_out);
end
