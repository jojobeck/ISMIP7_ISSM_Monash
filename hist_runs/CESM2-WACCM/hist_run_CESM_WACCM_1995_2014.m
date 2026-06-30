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
%                              ../../postprocessed_data/CESM2-WACCM/historical/
%   2b WriteISMIP6_NetCDF_test same as step 2 but truncated to the first 2
%                              annual outputs (1995-1996) for a quick
%                              compliance-checker run; output goes to
%                              ../../postprocessed_data/CESM2-WACCM/historical_test/
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
% front), onto the ISMIP6 grid cell containing each contour segment.
% lifmassbf/tendlifmassbf are written as copies of licalvf/tendlicalvf:
% md.frontalforcings.meltingrate is zero everywhere, so frontal melt = 0
% and the combined flux equals calving alone.

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
        % Grids and writes all ISMIP7 NetCDF files for the historical run using
        % the run-type-specific external functions in hist_runs/functions/.
        % See write_ismip7_2d_historical.m and write_ismip7_scalar_historical.m
        % for the full implementation (time filtering, grid setup, sign fixes).
        addpath('./../functions');

        md = loadmodel(org, 'HistRun_1995_2014');

        meta                    = struct();
        meta.experiment_id      = 'historical';
        meta.set_counter        = 'C001';
        meta.time_range         = '1995-2014';
        meta.ESM_id             = 'CESM2-WACCM';
        meta.forcing_member_id  = 'f001';
        meta.ISM_member_id      = 'm001';

        outdir = [proj_root 'postprocessed_data/CESM2-WACCM/historical/'];
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        [cfflux_tot, glflux_tot] = write_ismip7_2d_historical(md, outdir, meta);
        write_ismip7_scalar_historical(md, outdir, meta, cfflux_tot, glflux_tot);
    end % }}}

    % ================================================================= Step 2b
    if perform(org, 'WriteISMIP6_NetCDF_test') % {{{
        % Quick 2-year compliance test: writes only the first 2 annual outputs
        % (nominal years 1995-1996) to historical_test/ so the checker can be
        % run in seconds without processing the full 20-year solution.
        addpath('./../functions');

        md = loadmodel(org, 'HistRun_1995_2014');

        % Keep only the first 2 annual outputs (after removing the safety step).
        t_raw       = [md.results.TransientSolution.time];
        keep_annual = abs(t_raw - round(t_raw)) < 0.05;
        idx_annual  = find(keep_annual);
        md.results.TransientSolution = md.results.TransientSolution(idx_annual(1:2));

        meta                    = struct();
        meta.experiment_id      = 'historical';
        meta.set_counter        = 'C001';
        meta.time_range         = '1995-1996';
        meta.ESM_id             = 'CESM2-WACCM';
        meta.forcing_member_id  = 'f001';
        meta.ISM_member_id      = 'm001';

        outdir = [proj_root 'postprocessed_data/CESM2-WACCM/historical_test/'];
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        [cfflux_tot, glflux_tot] = write_ismip7_2d_historical(md, outdir, meta);
        write_ismip7_scalar_historical(md, outdir, meta, cfflux_tot, glflux_tot);
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
% (NetCDF writing helpers removed -- replaced by external functions
%  init/scripts/write_ismip7_2d.m and init/scripts/write_ismip7_scalar.m)

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
