function md = hist_run_CESM_WACCM_1995_2014(steps, loadonly)
% Historical validation run 1995-2014 (CMIP6 "historical" period only, no
% SSP126 extension) using CESM2-WACCM forcing, with the full ISMIP6
% Appendix-2 (Table A1) variable set requested so the saved
% md.results.TransientSolution has everything needed to write the
% <variable>_AIS_<GROUP>_<MODEL>_<EXP>.nc submission files directly.
%
% This does NOT rebuild ocean TF / SMB / calving-front forcing or the
% WAIS/EAIS/Peninsula region mask -- those are mesh- and forcing-only
% artefacts already built by ../../init/hist_run_tune_CESM_WACCM.m for the
% 1995-2020 run and are reused as-is:
%   ../../init/Models/AIS_ISMIP7_Relaxed_CESM_WACCM.mat   (relaxed start state)
%   ../../preprocessed_data/Ocean/Hist/CESM_WACCM_TF_1995_2020.mat
%   ../../preprocessed_data/Atmosphere/Hist/CESM_WACCM_SMB_1995_2020.mat
%   ../../preprocessed_data/Ocean/Hist/Greene_spclevelset_1995_2020.mat
%   ../../preprocessed_data/Ocean/Basins/HistRun_Regions.mat
% Each forcing series already spans 1995-2020; ISSM only reads the part of
% the time series up to md.timestepping.final_time, so truncating the run
% to 2014 needs no change to the forcing files themselves.
%
% Step map:
%   1  HistRun_1995_2014       transient 1995-2014 (loadonly=0 submit, =1 gather)
%   2  WriteISMIP6_NetCDF      grids TransientSolution onto the ISMIP6 8 km AIS
%                              grid and writes one NetCDF file per Appendix-2
%                              variable into
%                              ../../postprocessed_data/CESM2-WACCM/hist/
%   3  AIS_state_2015          saves the run's end-state geometry/masks as a
%                              fresh starting point (e.g. for an SSP126
%                              projection), with TransientSolution cleared
%   4  GLFluxSanityCheck       lightweight grounding-line-flux-only check,
%                              split out of step 2 since that step's full
%                              2D gridding is slow; computes ISSM's native
%                              GroundinglineMassFlux against our own
%                              native-mesh total and saves both to
%                              ../../postprocessed_data/CESM2-WACCM/hist/tables/
%                              (no plotting here -- see step 5)
%   5  GLFluxSanityCheckPlot   reads step 4's table and plots it into
%                              ../../postprocessed_data/CESM2-WACCM/hist/figures/
%                              kept separate from step 4 since figure()/
%                              saveas() segfaults in this MATLAB graphics
%                              stack even with -softwareopengl, so a crash
%                              here doesn't lose step 4's numbers
%   6  GLFluxOrientationDiagnostic   diagnostic for the ~6% gap remaining
%                              after step 4's Simpson's-rule fix: computes
%                              our total with vs without the eps_normal
%                              orientation correction (gl_flux_native_mesh's
%                              apply_correction arg), since ISSM's native
%                              GroundinglineMassFlux has no such correction
%                              at all; saves to hist/tables/ (no plotting)
%   7  GLFluxOrientationDiagnosticPlot   reads step 6's table and plots it
%                              into hist/figures/, kept separate for the
%                              same reason step 5 is separate from step 4
%
% Note (step 2): licalvf/ligroundf 2D maps -- ISSM's native
% requested_outputs only exposes these as domain-wide scalar reductions
% (FemModel::IcefrontMassFluxx / GroundinglineMassFluxx in the trunk C++
% source), so step 2 instead rasterises the same segment-wise V.N*H*L
% contour-integral used in
% ../../../ISSM-projects-cluster/AIS_1850_trendtuning/scripts/
% calc_GroundingLineFLux_transient_corr2.m (grounding line, normal-
% orientation corrected) and calc_CalvingFrontFLux_transient.m (ice
% front), onto the ISMIP6 grid cell containing each contour segment, in
% place of summing all segments into one AIS-wide number. lifmassbf/
% tendlifmassbf are NOT written: md.frontalforcings.meltingrate is zero
% everywhere in this run, so they would be an exact duplicate of
% licalvf/tendlicalvf (no separate frontal-melt component to distinguish
% them).

    if ~exist('loadonly','var'), loadonly = 0; end
    addpath('./../../init/scripts');

    org = organizer('repository', './Models/', 'prefix', 'AIS_ISMIP7_Hist1995_2014_', ...
                    'steps', steps, 'color', '34;47;2');
    clear steps;

    % ------------------------------------------------------------------ paths
    init_dir   = './../../init/';
    proj_root  = './../../';

    inputmodel_relax_cesm = [init_dir 'Models/AIS_ISMIP7_Relaxed_CESM_WACCM.mat'];

    preproc_ocean = [proj_root 'preprocessed_data/Ocean/'];
    preproc_atmo  = [proj_root 'preprocessed_data/Atmosphere/Hist/'];
    preproc_front = [proj_root 'preprocessed_data/Ocean/Hist/'];

    % Time range -- historical period only, no SSP126 splice
    start_year = 1995;
    end_year   = 2014;     % CMIP6 historical boundary -- run stops here

    % Forcing .mat files were built for 1995-2020 by init/hist_run_tune_CESM_WACCM.m;
    % they are loaded and used unmodified -- ISSM truncates to final_time itself.
    forcing_start_year = 1995;
    forcing_end_year   = 2020;

    % ------------------------------------------------------------ ISMIP6 outputs
    % Full Appendix-2 / Table A1 variable set this run can produce, mapped to
    % ISSM transient.requested_outputs names (verified against
    % $ISSM_DIR/src/c/shared/Enum/StringToEnumx.cpp on Gadi):
    %
    %   ISMIP6 name   ISSM requested_outputs name(s)
    %   -----------   -------------------------------
    %   lithk         Thickness
    %   orog          Surface
    %   base          Base
    %   topg          Bed
    %   hfgeoubed     (derive: static md.basalforcings.geothermalflux, NOT a
    %     requested_output -- requesting 'BasalforcingsGeothermalflux' as a
    %     TransientSolution output crashes the solve, since that Input is
    %     only marshalled by ThermalAnalysis, which is gated off by
    %     isthermal=0 in this run; the field itself is a real, non-NaN,
    %     per-vertex static map set far upstream in init/Par/
    %     Antartica_1995_Thwbedmap.par and carried through every step since,
    %     so step 2 grids it directly off the model object instead.)
    %   acabf         SmbMassBalance
    %   libmassbfgr   BasalforcingsGroundediceMeltingRate
    %   libmassbffl   BasalforcingsFloatingiceMeltingRate
    %   dlithkdt      (derive: finite-difference consecutive Thickness snapshots)
    %   xvelsurf/yvelsurf   Vx / Vy (this is a pure-SSA run --
    %     StressbalanceAnalysis::InputUpdateFromSolutionSSA assigns the same
    %     solved velocity to Vx/VxSurface/VxBase identically, no vertical
    %     shear profile, so surface/base/depth-mean velocity are all the
    %     same field; only the surface variant is written, xvelbase/
    %     yvelbase/xvelmean/yvelmean are omitted as redundant)
    %   zvelsurf/zvelbase   NOT WRITTEN -- 'Vz' is not a valid requested_output
    %     for this run: StressbalanceAnalysis::InputUpdateFromSolutionSSA
    %     (src/c/analyses/StressbalanceAnalysis.cpp on Gadi) never calls
    %     AddInput(VzEnum,...) for pure SSA; VzEnum is only populated by
    %     Full-Stokes/HOFS/SSAHO analyses, none of which apply here. Vertical
    %     velocity is not a solved quantity for this SSA continental run, so
    %     zvelsurf/zvelbase are omitted (acceptable/standard for 2D/SSA
    %     ISMIP6 submissions).
    %   litemptop   (derive: static md.initialization.temperature, NOT a
    %     requested_output -- this run has isthermal=0 so the thermal field
    %     never evolves; gridded directly off the model object like
    %     hfgeoubed above, rather than requested via TransientSolution.
    %     md.thermal.spctemperature is NOT used here: that is the thermal
    %     solve's Dirichlet boundary constraint array (NaN except where a
    %     constraint is imposed), not a representative temperature field.
    %     litempbotgr/litempbotfl (basal temperature) are NOT written --
    %     this 2D/SSA model has no vertically-resolved temperature profile,
    %     so they would just duplicate litemptop under a misleading label.)
    %   strbasemag    NOT WRITTEN -- 'BasalStress' is dead/unimplemented in
    %     this ISSM version: BasalStressEnum is registered in the Enum
    %     tables (EnumDefinitions.h/StringToEnumx.cpp/EnumToStringx.cpp) but
    %     never referenced anywhere else; Tria::ComputeBasalStress() is just
    %     "_error_(Not Implemented yet)" and is never called by any
    %     dispatcher (contrast with DeviatoricStresseffective, which IS
    %     wired via Element.cpp's case-switch). Requesting it crashes the
    %     solve the same way BasalforcingsGeothermalflux/Vz did. Computing
    %     basal shear stress for the Budd friction law used here
    %     (Sigma_b = coefficient^2 * Neff^(q/p) * |u_b|^(1/p), per
    %     friction.m's disp() docstring) would need to be done separately in
    %     post-processing from md.friction.coefficient/p/q + velocity --
    %     not implemented yet, omitted for now.
    %   licalvf / ligroundf (2D only)   derived in step 2 by
    %     rasterising/summing the V.N*H*L contour integral along the
    %     MaskIceLevelset=0 / MaskOceanLevelset=0 zero-contours (needs Vx, Vy,
    %     Thickness, MaskIceLevelset, MaskOceanLevelset, all requested below)
    %     -- ISSM has no native 2D output for these. lifmassbf NOT written
    %     (see note at top of file -- md.frontalforcings.meltingrate=0
    %     means it would be an exact duplicate of licalvf).
    %   sftgif/sftgrf/sftflf   (derive from MaskIceLevelset/MaskOceanLevelset
    %                           during NetCDF gridding -- no separate ISSM request)
    %   lim/limnsw    (derive: IceVolumeScaled / IceVolumeAboveFloatationScaled
    %     * rho_ice -- the *Scaled variants apply md.mesh.scale_factor, the
    %     polar-stereographic area-distortion correction set in step 1, so
    %     they are the physically correct totals across this continental
    %     mesh; plain IceVolume/IceVolumeAboveFloatation are not)
    %   iareagr/iareafl   GroundedAreaScaled / FloatingAreaScaled (same
    %     area-distortion correction as lim/limnsw above)
    %   tendacabf     TotalSmbScaled
    %   tendlibmassbf (derive: TotalGroundedBmbScaled + TotalFloatingBmbScaled)
    %   tendlibmassbffl   TotalFloatingBmbScaled
    %   tendlicalvf   IcefrontMassFluxLevelset (the more general
    %     zero-crossing test vs IcefrontMassFlux's narrower IsIcefront()
    %     criterion, see Tria.cpp -- sign cross-checked in step 2 against
    %     our own contour-integral cfflux_tot). tendlifmassbf NOT written,
    %     same reasoning as lifmassbf above.
    %   tendligroundf   GroundinglineMassFlux (sign cross-checked against
    %     our own glflux_tot, same reasoning as above)
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
    if perform(org, 'HistRun_1995_2014') % {{{

        md = loadmodel(inputmodel_relax_cesm);
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
        md.timestepping.final_time = end_year + 1;   % stop after the 2014 snapshot
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;   % annual snapshots

        md.transient.requested_outputs = ismip6_outputs;

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean (reuse the existing 1995-2020 series, ISSM truncates to final_time) ---
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_front 'CESM_WACCM_TF_' num2str(forcing_start_year) '_' num2str(forcing_end_year) '.mat']);
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

        % --- SMB (reuse the existing 1995-2020 series) ---
        load([preproc_atmo 'CESM_WACCM_SMB_corrected_' num2str(forcing_start_year) '_' num2str(forcing_end_year) '.mat']);
        md.smb        = SMBgradients();
        md.smb.smbref = smb_forcing;
        md.smb.b_pos  = bgrad_forcing;
        md.smb.b_neg  = bgrad_forcing;
        md.smb.href   = [md.geometry.surface ; start_year];

        % --- Calving front (reuse the existing 1995-2020 series) ---
        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 1;
        load([preproc_front 'Greene_spclevelset_' num2str(forcing_start_year) '_' num2str(forcing_end_year) '.mat']);
        md.levelset.spclevelset = greene_spclevelset;

        md.miscellaneous.name = ['HistRun_CESM_WACCM_' num2str(start_year) '_' num2str(end_year)];
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
    if perform(org, 'WriteISMIP6_NetCDF') % {{{
        % Grids md.results.TransientSolution onto the fixed ISMIP6 8 km AIS
        % grid (init/scripts/gridData.m) and writes one NetCDF file per
        % Appendix-2 variable, named <variable>_<IS>_<GROUP>_<MODEL>_<EXP>.nc.
        % licalvf/ligroundf 2D maps and their AIS-wide scalar
        % totals are derived via flux_along_contour_2d (see note at top of
        % file and the helper function at the bottom of this file).

        md   = loadmodel(org, 'HistRun_1995_2014');
        rhoi = md.materials.rho_ice;
        yts  = md.constants.yts;

        % --- Appendix-2 filename metadata -- TODO confirm official ISMIP7 codes
        IS        = 'AIS';
        GROUP     = 'MONASH';
        MODELNAME = 'ISSM';
        EXP       = 'C001';

        outdir = [proj_root 'postprocessed_data/CESM2-WACCM/hist/'];
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        % Initial state: use the first safety step (t≈1995.08), which has all
        % computed fields (BMB etc.).  Annual outputs: integer t=1996..2015.
        t_raw     = [md.results.TransientSolution.time];
        is_safety = abs(t_raw - round(t_raw)) > 0.01 & t_raw < 1996;
        idx_init  = find(is_safety, 1, 'first');
        if isempty(idx_init), idx_init = 1; end

        keep_annual = abs(t_raw - round(t_raw)) < 0.05 & round(t_raw) >= 1996;

        md.results.TransientSolution = [md.results.TransientSolution(idx_init), ...
                                          md.results.TransientSolution(keep_annual)];
        nT = 1 + sum(keep_annual);

        % Year labels for the 20 annual outputs
        t_annual = round(t_raw(keep_annual));   % [1996, 1997, ..., 2015]
        year_lbl = t_annual - 1;               % [1995, 1996, ..., 2014]

        time_units = 'days since 1995-01-01';
        % 365-day calendar (no leap years): days = (year - 1995) * 365.
        %
        % ST (state variables): snapshot time = end-of-year (Jan 1 of next year).
        %   Includes the initial state at day 0. No time_bnds.
        time_st = [0, (year_lbl - 1994) * 365];   % [0, 365, 730, ..., 7300]  1 x nT
        %
        % FL (flux variables): time = midpoint of year, with time_bnds.
        %   Annual entries only (no initial-state entry).
        nT_fl        = nT - 1;
        lb_fl        = (year_lbl - 1995) * 365;   % [0, 365, ..., 6935]
        ub_fl        = lb_fl + 365;                % [365, 730, ..., 7300]
        time_fl      = lb_fl + 182.5;             % [182.5, 547.5, ..., 7117.5]
        time_bnds_fl = [lb_fl; ub_fl]';           % nT_fl x 2

        missing_value = single(1e20);

        [~, xGrid, yGrid] = gridData(md, md.results.TransientSolution(1).Thickness);
        nx = length(xGrid); ny = length(yGrid);

        % ---- 2D variables: snapshot/rate fields. Units below follow ISSM's
        % MATLAB-side results loader (parseresultsfromdisk.m on Gadi), which
        % rescales select fields by yts on load -- NOT raw internal SI:
        %   Vx/Vy -> m/yr (multiplied by yts), so a 1/yts factor here
        %     converts to the 'm s-1' ISMIP6 wants. Pure SSA: Vx/Vy are the
        %     same field ISSM stores under VxSurface/VxBase/VxAverage, so
        %     only the surface (xvelsurf/yvelsurf) variant is written here;
        %     xvelbase/yvelbase/xvelmean/yvelmean are omitted as redundant.
        %   SmbMassBalance, BasalforcingsGroundediceMeltingRate,
        %     BasalforcingsFloatingiceMeltingRate -> m ice/yr (also
        %     multiplied by yts), so rho_ice/yts converts to kg m-2 s-1.
        %   Thickness/Surface/Base/Bed are NOT in that rescale list and stay
        %     raw SI (m). hfgeoubed/litemptop are handled separately below
        %     (static fields off the model object, not requested_outputs).
        % Sign on libmassbfgr/libmassbffl is flipped: ISSM melting rate is
        % positive=melting=mass loss, ISMIP6 wants positive=mass added.
        % 6th column (ST/FL) retained for documentation but all outputs now use
        % time_fl / time_bnds_fl (uniform annual bounds for every variable).
        vars2d = { ...
            'lithk',       'Thickness',                            1,        'land_ice_thickness',                                'm',         'ST'; ...
            'orog',        'Surface',                               1,        'surface_altitude',                                  'm',         'ST'; ...
            'base',        'Base',                                  1,        'base_altitude',                                     'm',         'ST'; ...
            'topg',        'Bed',                                   1,        'bedrock_altitude',                                  'm',         'ST'; ...
            'xvelsurf',    'Vx',                                    1/yts,    'land_ice_surface_x_velocity',                       'm s-1',     'ST'; ...
            'yvelsurf',    'Vy',                                    1/yts,    'land_ice_surface_y_velocity',                       'm s-1',     'ST'; ...
            'acabf',       'SmbMassBalance',                        rhoi/yts, 'land_ice_surface_specific_mass_balance_flux',       'kg m-2 s-1','FL'; ...
            'libmassbfgr', 'BasalforcingsGroundediceMeltingRate',  -rhoi/yts, 'land_ice_basal_specific_mass_balance_flux',         'kg m-2 s-1','FL'; ...
            'libmassbffl', 'BasalforcingsFloatingiceMeltingRate',  -rhoi/yts, 'land_ice_basal_specific_mass_balance_flux',         'kg m-2 s-1','FL'; ...
        };
        % Note: zvelsurf/zvelbase and strbasemag are omitted -- see the note
        % at the top of this file (Vz/BasalStress are not populated for
        % pure-SSA models).

        for k = 1:size(vars2d, 1)
            varname   = vars2d{k,1};
            issmfield = vars2d{k,2};
            scale     = vars2d{k,3};
            long_name = vars2d{k,4};
            units     = vars2d{k,5};
            vtype     = vars2d{k,6};
            if strcmp(vtype, 'ST')
                t_idx  = 1:nT;   t_vec = time_st;   t_bnds = [];
            else
                t_idx  = 2:nT;   t_vec = time_fl;   t_bnds = time_bnds_fl;
            end
            nTv    = length(t_idx);
            data3d = zeros(ny, nx, nTv, 'single');
            for ti = 1:nTv
                data3d(:,:,ti) = single(gridData(md, md.results.TransientSolution(t_idx(ti)).(issmfield)) * scale);
            end
            write_ismip6_2d(outdir, varname, IS, GROUP, MODELNAME, EXP, data3d, ...
                xGrid, yGrid, t_vec, t_bnds, time_units, long_name, units, missing_value);
        end

        % ---- hfgeoubed: static per-vertex field, NOT a requested_output --
        % see the note at the top of this file (requesting
        % 'BasalforcingsGeothermalflux' as a TransientSolution output
        % crashes the solve under isthermal=0). Gridded once off the model
        % object and replicated across all nT snapshots.
        hfgeoubed2d = single(gridData(md, md.basalforcings.geothermalflux));
        hfgeoubed3d = repmat(hfgeoubed2d, 1, 1, nT_fl);   % FL type
        write_ismip6_2d(outdir, 'hfgeoubed', IS, GROUP, MODELNAME, EXP, hfgeoubed3d, ...
            xGrid, yGrid, time_fl, time_bnds_fl, time_units, 'upward_geothermal_heat_flux_at_ground_level', 'W m-2', missing_value);

        % ---- litemptop: static per-vertex field, NOT a requested_output.
        % md.initialization.temperature ("temperature [K]", per
        % initialization.m's fielddisplay) is the ice temperature carried by
        % the model, inherited unchanged since isthermal=0 means it is never
        % solved transiently. md.thermal.spctemperature is NOT used here --
        % that field is the thermal solve's Dirichlet boundary constraint
        % array ("temperature constraints (NaN means no constraint) [K]",
        % per thermal.m), NaN everywhere a constraint isn't explicitly
        % imposed, not a representative temperature field. Gridded once off
        % the model object and replicated across all nT snapshots, same
        % pattern as hfgeoubed above. Only the surface variant (litemptop)
        % is written -- this 2D/SSA model has no vertically-resolved
        % temperature profile, so litempbotgr/litempbotfl (basal
        % temperature) would just be the same single field mislabelled as a
        % distinct quantity; omitted rather than written misleadingly.
        temp2d = single(gridData(md, md.initialization.temperature));
        temp3d = repmat(temp2d, 1, 1, nT);   % ST type: all nT snapshots
        write_ismip6_2d(outdir, 'litemptop', IS, GROUP, MODELNAME, EXP, temp3d, ...
            xGrid, yGrid, time_st, [], time_units, 'temperature_at_top_of_ice_sheet_model', 'K', missing_value);

        % ---- dlithkdt: finite-difference of Thickness between snapshots ----
        % dlithkdt is FL: use annual snapshots only (indices 2:nT, skip initial state)
        H_grid = zeros(ny, nx, nT_fl);
        for t = 1:nT_fl, H_grid(:,:,t) = gridData(md, md.results.TransientSolution(t+1).Thickness); end
        dlithkdt3d = zeros(ny, nx, nT_fl, 'single');
        for t = 1:nT_fl
            if t < nT_fl
                dlithkdt3d(:,:,t) = single((H_grid(:,:,t+1) - H_grid(:,:,t)) / yts);
            else
                dlithkdt3d(:,:,t) = single((H_grid(:,:,t) - H_grid(:,:,t-1)) / yts);
            end
        end
        write_ismip6_2d(outdir, 'dlithkdt', IS, GROUP, MODELNAME, EXP, dlithkdt3d, ...
            xGrid, yGrid, time_fl, time_bnds_fl, time_units, 'tendency_of_land_ice_thickness', 'm s-1', missing_value);

        % ---- sftgif/sftgrf/sftflf: area fractions derived from the ice/ocean levelsets ----
        % mask.m's fielddisplay documents the sign convention precisely:
        % "presence of ocean if < 0, coastline/grounding line if = 0, no
        % ocean if > 0" -- i.e. ocean_ls < 0 means OCEAN PRESENT (floating),
        % ocean_ls >= 0 means NO OCEAN (grounded). (Confirmed against an
        % earlier sftgrf/sftflf swap caught by visually inspecting
        % plot_historic_runs.py's mask figure -- grounded interior was
        % showing up labelled "floating" and vice versa.)
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

        % ---- licalvf/ligroundf: rasterised contour-integral maps ----
        % Same V.N*H*L segment method as calc_GroundingLineFLux_transient_corr2.m
        % (grounding line, orientation-corrected so grounded->floating is
        % positive, i.e. positive = loss of grounded mass) and
        % calc_CalvingFrontFLux_transient.m (ice front, unsigned magnitude),
        % but binned onto the ISMIP6 grid cell containing each segment
        % instead of summed into one AIS-wide number.
        % FL type: use annual entries only (t+1 skips the initial safety-step entry)
        cellsize = xGrid(2) - xGrid(1);
        glflux2d   = zeros(ny, nx, nT_fl, 'single');
        glflux_tot = zeros(1, nT_fl);
        cfflux2d   = zeros(ny, nx, nT_fl, 'single');
        cfflux_tot = zeros(1, nT_fl);
        for t = 1:nT_fl
            sol = md.results.TransientSolution(t+1);
            [gl2d, gltot] = flux_along_contour_2d(md, sol, sol.MaskOceanLevelset, true,  xGrid, yGrid, cellsize);
            [cf2d, cftot] = flux_along_contour_2d(md, sol, sol.MaskIceLevelset,   false, xGrid, yGrid, cellsize);
            glflux2d(:,:,t) = gl2d; glflux_tot(t) = gltot;
            cfflux2d(:,:,t) = cf2d; cfflux_tot(t) = cftot;
        end
        write_ismip6_2d(outdir, 'ligroundf', IS, GROUP, MODELNAME, EXP, glflux2d, ...
            xGrid, yGrid, time_fl, time_bnds_fl, time_units, 'land_ice_specific_mass_flux_at_grounding_line', 'kg m-2 s-1', missing_value);
        write_ismip6_2d(outdir, 'licalvf', IS, GROUP, MODELNAME, EXP, cfflux2d, ...
            xGrid, yGrid, time_fl, time_bnds_fl, time_units, 'land_ice_specific_mass_flux_due_to_calving', 'kg m-2 s-1', missing_value);
        % lifmassbf NOT written: md.frontalforcings.meltingrate is zero
        % throughout this run, so it would be an exact duplicate of licalvf
        % (no separate frontal-melt component to distinguish them).

        % ---- scalar (t) variables ----
        % *Scaled variants apply md.mesh.scale_factor (the polar-
        % stereographic area-distortion correction set in step 1,
        % FemModel::IceVolumex/GroundedAreax/etc with scaled=true) -- the
        % physically correct totals across this continental mesh.
        % IceVolumeScaled/IceVolumeAboveFloatationScaled/GroundedAreaScaled/
        % FloatingAreaScaled are NOT in parseresultsfromdisk.m's yts-rescale
        % list -- raw SI (m3, m2). TotalSmbScaled IS in that list, rescaled
        % to Gt/yr (/1e12*yts), so 1e12/yts converts back to kg s-1.
        scalars = { ...
            'lim',     'IceVolumeScaled',                rhoi,     'land_ice_mass',                                          'kg',    'ST'; ...
            'limnsw',  'IceVolumeAboveFloatationScaled',  rhoi,     'land_ice_mass_not_displacing_sea_water',                 'kg',    'ST'; ...
            'iareagr', 'GroundedAreaScaled',              1,        'grounded_land_ice_area',                                 'm2',    'ST'; ...
            'iareafl', 'FloatingAreaScaled',              1,        'floating_ice_shelf_area',                                'm2',    'ST'; ...
            'tendacabf','TotalSmbScaled',                 1e12/yts, 'tendency_of_land_ice_mass_due_to_surface_mass_balance',  'kg s-1','FL'; ...
        };
        for k = 1:size(scalars, 1)
            varname   = scalars{k,1};
            issmfield = scalars{k,2};
            scale     = scalars{k,3};
            long_name = scalars{k,4};
            units     = scalars{k,5};
            vtype     = scalars{k,6};
            if strcmp(vtype, 'ST')
                t_idx  = 1:nT;   t_vec = time_st;   t_bnds = [];
            else
                t_idx  = 2:nT;   t_vec = time_fl;   t_bnds = time_bnds_fl;
            end
            data1d = zeros(1, length(t_idx));
            for ti = 1:length(t_idx)
                data1d(ti) = md.results.TransientSolution(t_idx(ti)).(issmfield) * scale;
            end
            write_ismip6_scalar(outdir, varname, IS, GROUP, MODELNAME, EXP, data1d, ...
                t_vec, t_bnds, time_units, long_name, units, missing_value);
        end

        % tendlibmassbf / tendlibmassbffl: TotalGroundedBmbScaled/
        % TotalFloatingBmbScaled (area-distortion corrected, like the
        % scalars above) are also rescaled to Gt/yr by
        % parseresultsfromdisk.m (/1e12*yts), so 1e12/yts converts back to
        % kg s-1; sign flipped to match ISMIP6's positive=mass-added
        % convention, mirroring libmassbfgr/libmassbffl above.
        tgbmb = zeros(1, nT_fl); tfbmb = zeros(1, nT_fl);
        for t = 1:nT_fl
            tgbmb(t) = md.results.TransientSolution(t+1).TotalGroundedBmbScaled * 1e12/yts;
            tfbmb(t) = md.results.TransientSolution(t+1).TotalFloatingBmbScaled * 1e12/yts;
        end
        write_ismip6_scalar(outdir, 'tendlibmassbf', IS, GROUP, MODELNAME, EXP, -(tgbmb + tfbmb), ...
            time_fl, time_bnds_fl, time_units, 'tendency_of_land_ice_mass_due_to_basal_mass_balance', 'kg s-1', missing_value);
        write_ismip6_scalar(outdir, 'tendlibmassbffl', IS, GROUP, MODELNAME, EXP, -tfbmb, ...
            time_fl, time_bnds_fl, time_units, 'tendency_of_land_ice_mass_due_to_basal_mass_balance_at_ice_shelf_base', 'kg s-1', missing_value);

        % tendlicalvf / tendligroundf: use ISSM's native
        % FemModel::IcefrontMassFluxLevelsetx / GroundinglineMassFluxx
        % (requested above) rather than our own contour-integral totals
        % (cfflux_tot/glflux_tot, still used for the 2D maps above, since
        % ISSM has no native 2D output for these). Both are rescaled to
        % Gt/yr by parseresultsfromdisk.m (/1e12*yts, same as
        % TotalSmbScaled etc above), so 1e12/yts converts back to kg s-1.
        %
        % Sign convention: IcefrontMassFluxLevelset/GroundinglineMassFlux's
        % internal normal orientation (Tria.cpp, NormalSection) isn't
        % side-by-side documented against the Table A1 sign convention
        % (positive = mass loss at grounding line / calving), and reverse-
        % engineering it from the geometry code is error-prone -- so instead
        % of asserting a sign, it is checked empirically here against our
        % own already-validated cfflux_tot/glflux_tot (same physical
        % quantity, contour-integral method, known sign per
        % flux_along_contour_2d's docstring) and flipped to match if the
        % two disagree in sign.
        cf_native = zeros(1, nT_fl); gl_native = zeros(1, nT_fl);
        for t = 1:nT_fl
            cf_native(t) = md.results.TransientSolution(t+1).IcefrontMassFluxLevelset * 1e12/yts;
            gl_native(t) = md.results.TransientSolution(t+1).GroundinglineMassFlux   * 1e12/yts;
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
        % tendlifmassbf NOT written -- same reasoning as lifmassbf above.
        write_ismip6_scalar(outdir, 'tendligroundf', IS, GROUP, MODELNAME, EXP, gl_native, ...
            time_fl, time_bnds_fl, time_units, 'tendency_of_land_ice_mass_due_to_flux_at_grounding_line', 'kg s-1', missing_value);
    end % }}}

    % ================================================================= Step 3
    if perform(org, 'AIS_state_2015') % {{{

        md   = loadmodel(org, 'HistRun_1995_2014');
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
    end % }}}

    % ================================================================= Step 4
    if perform(org, 'GLFluxSanityCheck') % {{{
        % Lightweight sanity check, split out of step 2 (WriteISMIP6_NetCDF)
        % since that step's full per-variable 2D gridding is slow and this
        % check only needs the grounding-line contour integral + the native
        % ISSM scalar, not the rest of step 2's output. Computes ISSM's
        % native GroundinglineMassFlux against our own total, BEFORE any
        % sign-correction, so the sign/magnitude agreement that step 2's
        % automatic warning-and-flip relies on can be checked. Saved to a
        % table/CSV here (not plotted) -- figure()/saveas() segfaults on
        % this compute node's MATLAB graphics stack even with
        % -softwareopengl (see SubmitHist1995_2014.errlog), so the
        % numbers are persisted first in case step 5's plot attempt
        % crashes again.
        %
        % Our own total here is computed by gl_flux_native_mesh (below),
        % a direct port of calc_GroundingLineFLux_transient_corr2.m's
        % per-segment V.N*H*L contour integral summed over the ENTIRE
        % native ISSM mesh -- deliberately NOT flux_along_contour_2d (used
        % in step 2 for the 2D maps), which interpolates/bins onto the
        % fixed ISMIP6 grid first and silently drops any segment whose
        % midpoint falls outside that grid's bounding box from its total.
        % That grid-bounds filtering is exactly the kind of discrepancy
        % this sanity check is meant to catch, so it must not be present
        % on either side of the comparison.
        md  = loadmodel(org, 'HistRun_1995_2014');
        yts = md.constants.yts;

        EXP = 'C001';

        nT = length(md.results.TransientSolution);
        time_yr = zeros(1, nT);
        for t = 1:nT, time_yr(t) = md.results.TransientSolution(t).time; end

        glflux_tot = zeros(1, nT);
        gl_native  = zeros(1, nT);
        for t = 1:nT
            sol = md.results.TransientSolution(t);
            glflux_tot(t) = gl_flux_native_mesh(md, sol);
            gl_native(t)  = sol.GroundinglineMassFlux * 1e12/yts;
        end

        tabledir = [proj_root 'postprocessed_data/tables/CESM2-WACCM/historical/'];
        if ~exist(tabledir, 'dir'), mkdir(tabledir); end
        T = table(time_yr', glflux_tot' * yts/1e12, gl_native' * yts/1e12, ...
            'VariableNames', {'year', 'glflux_native_mesh_Gt_per_yr', 'GroundinglineMassFlux_Gt_per_yr'});
        tablepath = [tabledir 'groundinglineflux_sanitycheck_' EXP '.csv'];
        writetable(T, tablepath);
        fprintf('Saved: %s\n', tablepath);
    end % }}}

    % ================================================================= Step 5
    if perform(org, 'GLFluxSanityCheckPlot') % {{{
        % Reads back the table saved by GLFluxSanityCheck and plots it, in
        % its own step so a graphics crash here (see note above) doesn't
        % risk re-running -- or losing -- step 4's already-computed numbers.
        EXP = 'C001';
        tabledir  = [proj_root 'postprocessed_data/tables/CESM2-WACCM/historical/'];
        tablepath = [tabledir 'groundinglineflux_sanitycheck_' EXP '.csv'];
        T = readtable(tablepath);

        figdir = [proj_root 'postprocessed_data/figures/CESM2-WACCM/historical/'];
        if ~exist(figdir, 'dir'), mkdir(figdir); end
        figure('visible', 'off');
        plot(T.year, T.glflux_native_mesh_Gt_per_yr, 'b-o', 'LineWidth', 1.5); hold on;
        plot(T.year, T.GroundinglineMassFlux_Gt_per_yr, 'r-o', 'LineWidth', 1.5);
        xlabel('year'); ylabel('Grounding line flux [Gt/yr]');
        legend('our native-mesh total (gl\_flux\_native\_mesh)', 'ISSM native (GroundinglineMassFlux)', 'Location', 'best');
        title('Sanity check: grounding line flux, native vs contour-integral (pre-sign-correction)');
        grid on;
        saveas(gcf, [figdir 'sanitycheck_groundinglineflux_' EXP '.png']);
        close(gcf);
    end % }}}

    % ================================================================= Step 6
    if perform(org, 'GLFluxOrientationDiagnostic') % {{{
        % Diagnostic for the ~6% gap remaining after the Simpson's-rule fix
        % in gl_flux_native_mesh: computes our total BOTH with (corrected)
        % and without (uncorrected) the eps_normal grounded->floating
        % orientation correction, alongside ISSM's native
        % GroundinglineMassFlux (which has no such correction at all in
        % Tria::GroundinglineMassFlux). If "uncorrected" lands closer to
        % the native value than "corrected" does, that confirms the
        % remaining gap is this orientation-correction asymmetry rather
        % than residual quadrature error. Saved to a table (not plotted)
        % for the same reason as step 4 -- figure()/saveas() segfaults on
        % this compute node even with -softwareopengl.
        md  = loadmodel(org, 'HistRun_1995_2014');
        yts = md.constants.yts;

        EXP = 'C001';

        nT = length(md.results.TransientSolution);
        time_yr = zeros(1, nT);
        for t = 1:nT, time_yr(t) = md.results.TransientSolution(t).time; end

        gl_corrected   = zeros(1, nT);
        gl_uncorrected = zeros(1, nT);
        gl_native      = zeros(1, nT);
        for t = 1:nT
            sol = md.results.TransientSolution(t);
            gl_corrected(t)   = gl_flux_native_mesh(md, sol, true);
            gl_uncorrected(t) = gl_flux_native_mesh(md, sol, false);
            gl_native(t)      = sol.GroundinglineMassFlux * 1e12/yts;
        end

        tabledir = [proj_root 'postprocessed_data/tables/CESM2-WACCM/historical/'];
        if ~exist(tabledir, 'dir'), mkdir(tabledir); end
        T = table(time_yr', gl_corrected' * yts/1e12, gl_uncorrected' * yts/1e12, gl_native' * yts/1e12, ...
            'VariableNames', {'year', 'corrected_Gt_per_yr', 'uncorrected_Gt_per_yr', 'GroundinglineMassFlux_Gt_per_yr'});
        tablepath = [tabledir 'groundinglineflux_orientation_diagnostic_' EXP '.csv'];
        writetable(T, tablepath);
        fprintf('Saved: %s\n', tablepath);
    end % }}}

    % ================================================================= Step 7
    if perform(org, 'GLFluxOrientationDiagnosticPlot') % {{{
        % Reads back the table saved by GLFluxOrientationDiagnostic and
        % plots it, in its own step for the same reason step 5 is separate
        % from step 4 -- a graphics crash here doesn't risk step 6's
        % already-computed numbers.
        EXP = 'C001';
        tabledir  = [proj_root 'postprocessed_data/tables/CESM2-WACCM/historical/'];
        tablepath = [tabledir 'groundinglineflux_orientation_diagnostic_' EXP '.csv'];
        T = readtable(tablepath);

        figdir = [proj_root 'postprocessed_data/figures/CESM2-WACCM/historical/'];
        if ~exist(figdir, 'dir'), mkdir(figdir); end
        figure('visible', 'off');
        plot(T.year, T.corrected_Gt_per_yr, 'b-o', 'LineWidth', 1.5); hold on;
        plot(T.year, T.uncorrected_Gt_per_yr, 'g-s', 'LineWidth', 1.5);
        plot(T.year, T.GroundinglineMassFlux_Gt_per_yr, 'r-o', 'LineWidth', 1.5);
        xlabel('year'); ylabel('Grounding line flux [Gt/yr]');
        legend('ours, orientation-corrected', 'ours, uncorrected', 'ISSM native (GroundinglineMassFlux)', 'Location', 'best');
        title('Diagnostic: does dropping the orientation correction explain the gap vs native?');
        grid on;
        saveas(gcf, [figdir 'orientation_diagnostic_groundinglineflux_' EXP '.png']);
        close(gcf);
    end % }}}
end

% ------------------------------------------------------------------ helpers
function fname = ismip6_filename(outdir, varname, IS, GROUP, MODELNAME, EXP)
    fname = fullfile(outdir, sprintf('%s_%s_%s_%s_%s.nc', varname, IS, GROUP, MODELNAME, EXP));
end

function write_ismip6_2d(outdir, varname, IS, GROUP, MODELNAME, EXP, data3d, ...
        xGrid, yGrid, time_vec, time_bnds, time_units, long_name, units, missing_value)
    % data3d is (ny, nx, nt), as returned by gridData/InterpFromMeshToGrid;
    % permuted to (nx, ny, nt) on write to match the {'x','y','time'} dimension order.
    % time_bnds is [] for ST (state) variables (snapshot time only, no
    % bounds), or an nt-by-2 matrix of [year-start, year-end] days-since-
    % basetime for FL (flux) variables, per Appendix-2 A2.3.2.
    fname = ismip6_filename(outdir, varname, IS, GROUP, MODELNAME, EXP);
    if exist(fname, 'file'), delete(fname); end
    nx = length(xGrid); ny = length(yGrid); nt = length(time_vec);

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
        % MATLAB's nccreate stores Dimensions in reverse order in the file
        % (column-major vs netCDF row-major), so {'bnds',2,'time',Inf} here
        % becomes time_bnds(time,bnds) in the file, matching the Appendix-2
        % A2.3.2 example exactly; the array passed to ncwrite must match
        % that declared order too, hence the transpose.
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
    % time_bnds is [] for ST variables, or an nt-by-2 matrix for FL
    % variables -- see write_ismip6_2d.
    fname = ismip6_filename(outdir, varname, IS, GROUP, MODELNAME, EXP);
    if exist(fname, 'file'), delete(fname); end
    nt = length(time_vec);

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
    % Direct port of calc_GroundingLineFLux_transient_corr2.m's per-segment
    % V.N*H*L contour integral (grounding line, MaskOceanLevelset=0,
    % normal oriented grounded->floating), summed over the ENTIRE native
    % ISSM mesh -- no gridData()/ISMIP6-grid interpolation or grid-bounds
    % filtering involved, unlike flux_along_contour_2d (used for the 2D
    % maps in step 2), which can silently drop segments outside the fixed
    % ISMIP6 grid's bounding box from its scalar total. Used by the
    % GLFluxSanityCheck step, to compare against the native
    % GroundinglineMassFlux output on equal footing.
    %
    % Integration: Simpson's rule (endpoints + midpoint), not the simple
    % two-point average the original port used. Thickness/Vx/Vy are
    % piecewise-linear (P1) FEM fields, so H.(V.N) is an exact quadratic
    % along any straight segment within one element; Simpson's rule
    % integrates quadratics exactly, matching ISSM's own 3-point Gauss
    % quadrature in Tria::GroundinglineMassFlux (the two-point-average
    % version under-integrated this and was the source of the ~10% low
    % bias seen against the native GroundinglineMassFlux comparison).
    %
    % apply_correction (default true): whether to apply the eps_normal
    % grounded->floating orientation correction (see below). Tria::
    % GroundinglineMassFlux has NO equivalent correction at all -- it just
    % uses NormalSection's raw per-segment normal. Pass false here to
    % reproduce that uncorrected behaviour, as a diagnostic for whether
    % this orientation-correction asymmetry (rather than residual
    % quadrature error) explains the remaining gap against the native
    % value after the Simpson's-rule fix above.
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

    xm = 0.5 * (x1 + x2); ym = 0.5 * (y1 + y2);
    if apply_correction
        eps_n = 1000;  % meters, per calc_GroundingLineFLux_transient_corr2.m
        phi_plus  = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm + eps_n*Nx, ym + eps_n*Ny);
        phi_minus = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm - eps_n*Nx, ym - eps_n*Ny);
        flip = (~isnan(phi_plus) & ~isnan(phi_minus) & (phi_plus > phi_minus));
        Nx(flip) = -Nx(flip); Ny(flip) = -Ny(flip);
    end

    % Midpoint sample for Simpson's rule -- interpolated directly AT the
    % segment midpoint, not averaged from the two endpoint values.
    Vxm = InterpFromMesh2d(elems, x, y, sol.Vx,        xm, ym);
    Vym = InterpFromMesh2d(elems, x, y, sol.Vy,        xm, ym);
    Hm  = InterpFromMesh2d(elems, x, y, sol.Thickness, xm, ym);
    good3 = ~any(isnan([Vxm Vym Hm]), 2);

    f1 = H1 .* (Vx1.*Nx + Vy1.*Ny);
    f2 = H2 .* (Vx2.*Nx + Vy2.*Ny);
    fm = Hm .* (Vxm.*Nx + Vym.*Ny);
    % Fall back to the two-point average (trapezoidal) for the rare
    % segment whose midpoint sample itself is NaN (e.g. falls just outside
    % the mesh), so those segments still contribute instead of being
    % silently dropped.
    fm(~good3) = 0.5 * (f1(~good3) + f2(~good3));

    secflux_m3_yr  = L/6 .* (f1 + 4*fm + f2);
    flux_total_kgs = sum(secflux_m3_yr, 'omitnan') * md.materials.rho_ice / md.constants.yts;
end

function [flux2d_kgm2s, flux_total_kgs] = flux_along_contour_2d(md, sol, contour_field, signed_normal, xGrid, yGrid, cellsize)
    % Rasterises a V.N*H*L flux along the zero-contour of contour_field
    % (MaskOceanLevelset for the grounding line, MaskIceLevelset for the ice
    % front) onto the ISMIP6 grid cell containing each contour segment's
    % midpoint, instead of summing all segments into one AIS-wide number.
    % Segment-flux method matches
    % ../../../ISSM-projects-cluster/AIS_1850_trendtuning/scripts/
    % calc_GroundingLineFLux_transient_corr2.m (signed_normal=true: normal
    % oriented grounded->floating, positive = loss of grounded mass) and
    % calc_CalvingFrontFLux_transient.m (signed_normal=false: unsigned
    % magnitude, as that script takes abs(flux)).
    % Cells the contour never touches are NaN (matching the missing_value
    % convention every other 2D field uses via write_ismip6_2d's
    % data3d(isnan(data3d))=missing_value), not 0 -- consistent with Table
    % A1's "only for grid cells in contact with grounding line/ocean".
    %
    % Integration: Simpson's rule (endpoints + midpoint), matching the fix
    % in gl_flux_native_mesh -- Thickness/Vx/Vy are piecewise-linear (P1)
    % FEM fields, so H.(V.N) is an exact quadratic along any straight
    % segment within one element, and Simpson's rule integrates quadratics
    % exactly (the previous two-point-average version under-integrated
    % this, the source of the ~10% low bias found via the
    % GLFluxSanityCheck/GLFluxSanityCheckPlot steps).
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

    Nx = -dy(good) ./ L; Ny = dx(good) ./ L;
    Vx1 = vx(idx);   Vy1 = vy(idx);   H1 = h(idx);
    Vx2 = vx(idx+1); Vy2 = vy(idx+1); H2 = h(idx+1);
    good2 = ~any(isnan([Vx1 Vy1 H1 Vx2 Vy2 H2]), 2);
    if ~any(good2), return; end
    x1 = x1(good2); y1 = y1(good2); x2 = x2(good2); y2 = y2(good2); L = L(good2);
    Nx = Nx(good2); Ny = Ny(good2);
    Vx1 = Vx1(good2); Vy1 = Vy1(good2); H1 = H1(good2);
    Vx2 = Vx2(good2); Vy2 = Vy2(good2); H2 = H2(good2);

    xm = 0.5 * (x1 + x2); ym = 0.5 * (y1 + y2);
    if signed_normal
        eps_n = 1000;  % meters, per calc_GroundingLineFLux_transient_corr2.m
        phi_plus  = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm + eps_n*Nx, ym + eps_n*Ny);
        phi_minus = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm - eps_n*Nx, ym - eps_n*Ny);
        flip = (~isnan(phi_plus) & ~isnan(phi_minus) & (phi_plus > phi_minus));
        Nx(flip) = -Nx(flip); Ny(flip) = -Ny(flip);
    end

    % Midpoint sample for Simpson's rule -- interpolated directly AT the
    % segment midpoint, not averaged from the two endpoint values.
    Vxm = InterpFromMesh2d(elems, x, y, sol.Vx,        xm, ym);
    Vym = InterpFromMesh2d(elems, x, y, sol.Vy,        xm, ym);
    Hm  = InterpFromMesh2d(elems, x, y, sol.Thickness, xm, ym);
    good3 = ~any(isnan([Vxm Vym Hm]), 2);

    f1 = H1 .* (Vx1.*Nx + Vy1.*Ny);
    f2 = H2 .* (Vx2.*Nx + Vy2.*Ny);
    fm = Hm .* (Vxm.*Nx + Vym.*Ny);
    % Fall back to the two-point average (trapezoidal) for the rare
    % segment whose midpoint sample itself is NaN (e.g. falls just outside
    % the mesh), so those segments still contribute instead of being
    % silently dropped.
    fm(~good3) = 0.5 * (f1(~good3) + f2(~good3));

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
