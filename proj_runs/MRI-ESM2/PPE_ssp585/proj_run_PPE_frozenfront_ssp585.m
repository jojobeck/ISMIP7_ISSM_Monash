function md = proj_run_PPE_frozenfront_ssp585(X, steps, loadonly)
% PPE (Perturbed-Parameter Ensemble) run: level-set FROZEN at present-day
% (2015) ice extent, to isolate BMB-driven retreat from calving-front
% evolution entirely -- no ice-shelf collapse mask, no spclevelset time
% series, no floating-only gating. MRI-ESM2-0 copy of
% ../../CESM2-WACCM/PPE_ssp585/proj_run_PPE_frozenfront_ssp585.m (P011-P013)
% -- identical mechanism, only the model labels/paths differ. P-numbers
% continue from that family (P001-P009 ssp370, P011-P019 CESM2-WACCM
% ssp585), so every PPE member has its own unique P-number:
%   X='p021'  K = mode              (gamma0_local)      -- SAME K as core ssp585 run
%   X='p022'  K = 95th percentile   (gamma0_local_95th)
%   X='p023'  K = 5th percentile    (gamma0_local_5th)
% (gamma0_local.mat's 4 variants -- mode/5th/50th/95th -- were built by
% init/meltMip_ensemble.m's 'save_gamma0_local' step from
% run_parameter_selection.py's own K selection; nothing here re-derives
% them -- same file for every scenario, not ssp585-specific.) Each X
% internally resolves both the K variant AND this PPE member's output
% label (P021/P022/P023 respectively).
%
% Mechanism: md.transient.ismovingfront = 0, matching
% ../ctrl/proj_run_MRI_ESM2_ctrl_2015_2300.m's own frozen-front
% approach exactly (see that script's header for the full rationale) --
% the level-set equation is never solved, so the ice front/grounding line
% stay exactly where they are in the AIS_state_2015 starting state for
% the whole run. UNLIKE ctrl (which additionally uses a fixed 2000-2029
% CLIMATOLOGICAL forcing, to isolate ice-sheet-model drift from climate
% forcing), this script uses the REAL ssp585 TRANSIENT TF/SMB forcing
% (same as the core ssp585 run) -- only the calving front is frozen, so
% any mass change here is attributable to BMB (via the varying K) and SMB
% alone, with calving-front retreat/advance eliminated as a factor.
%
% RESTART-SPLIT at 2150/2151: a single continuous 2015-2300 solve (what
% ../../CESM2-WACCM/PPE_ssp370/proj_run_PPE_frozenfront_ssp370.m does for its own
% shorter 2015-2100 span) was tried here FIRST and got OOM-killed on a
% 48cpu/190GB 'normal'-queue node (PBS: "Job ... has exceeded memory
% allocation") -- holding the full ~285-year TransientSolution (every
% annual snapshot of ~15-20 full-mesh ISMIP6 output fields) in RAM at
% once is simply too much for a run this long. Split into two ~140-year
% segments instead, exactly the same restart-split technique already
% used for this project's other 300-year-class runs (see
% ../../CESM2-WACCM/ssp585/old_scripts/proj_run_CESM_WACCM_ssp585_2015_2300.m, steps
% 4-6: ProjRun_2015_2150 / AIS_state_2151 / ProjRun_2151_2300) -- each
% segment's own accumulated TransientSolution is roughly half the size,
% comfortably within the 190GB budget (confirmed: the ssp370 version's
% single-segment 85-year run completed fine, and each segment here is
% not much longer than that). The TF/SMB forcing arrays are also TRIMMED
% to each segment's own time window before being assigned to
% md.basalforcings.tf / md.smb.smbref (same technique as that precedent)
% -- this roughly halves the .bin file ISSM writes for the solve itself,
% on top of halving the accumulated output.
%
% NO PREPROCESSING -- loads the ALREADY-BUILT TF/SMB .mat files the core
% ssp585 run produced (confirmed present on disk), same full files for
% BOTH segments (just trimmed differently after loading):
%   preprocessed_data/Ocean/Proj/MRI-ESM2-0/ssp585/MRI_ESM2_TF_ssp585_2015_2300.mat
%   preprocessed_data/Atmosphere/Proj/MRI-ESM2-0/ssp585/MRI_ESM2_SMB_ssp585_2015_2300.mat
%
% Calling convention: (X, steps, loadonly) -- X picks the PPE member,
% steps picks which organizer step(s) of THIS run to execute, loadonly is
% the standard two-phase submit(0)/gather(1) flag used ONLY within the
% two solve steps (1 and 3).
%
% Step map (organizer, per PPE member):
%   1  ProjRun_2015_2150   solve the 2015-2150 transient (loadonly=0
%                           submit, =1 gather + save the model)
%   2  AIS_state_2151      load step 1's gathered output, save its 2151
%                           end-state as a clean restart point with
%                           TransientSolution cleared (no loadonly --
%                           this step is cheap/local, not a PBS solve)
%   3  ProjRun_2151_2300   solve the 2151-2300 transient, continuing from
%                           AIS_state_2151 (loadonly=0 submit, =1 gather
%                           + save the model)
%   4  WriteNetCDF          loads BOTH gathered segments (steps 1 and 3),
%                           concatenates their TransientSolution, and
%                           writes the ISMIP6 NetCDF output -- run this
%                           ONLY after both segments' gathers have
%                           completed and saved successfully
%
% Output: this member's P-number (e.g. 'P021') is used as meta.set_counter,
% matching this group's PPE filename convention directly, e.g.
%   libmassbffl_AIS_Monash_ISSM_m001_MRI-ESM2-0_f001_ssp585_P021_2015-2300.nc
% (ISM_member_id stays fixed at 'm001' -- same physical ice-sheet model
% for every PPE member, an UNRELATED fixed ISMIP7-filename-convention
% field, not this script's own X/P-number -- the P-number is what varies
% per experiment, not the ISM_member_id or the ESM_id).
%
% Each PPE member's solved model is saved to its OWN dedicated
% ./Models_P<NNN>/ repository (e.g. ./Models_P021/ for p021), NOT a
% shared Models/ or Models_yearly/ folder -- keeps every PPE variant's
% output completely separate on disk.
%
% Usage:
%   proj_run_PPE_frozenfront_ssp585('p021', [1], 0)   % submit p021 segment 1 (2015-2150)
%   proj_run_PPE_frozenfront_ssp585('p021', [1], 1)   % gather p021 segment 1
%   proj_run_PPE_frozenfront_ssp585('p021', [2])      % save p021's 2151 restart state (after segment 1 gather)
%   proj_run_PPE_frozenfront_ssp585('p021', [3], 0)   % submit p021 segment 2 (2151-2300)
%   proj_run_PPE_frozenfront_ssp585('p021', [3], 1)   % gather p021 segment 2
%   proj_run_PPE_frozenfront_ssp585('p021', [4])      % write NetCDF for p021 (after both gathers)
%   proj_run_PPE_frozenfront_ssp585('p022', [1], 0)   % submit p022 (K=95th) segment 1
%   proj_run_PPE_frozenfront_ssp585('p023', [1], 0)   % submit p023 (K=5th) segment 1
%
% Run with ISSM already on the MATLAB path (devpath already called), e.g.:
%   matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
%     addpath('$ISSM_DIR/lib'); proj_run_PPE_frozenfront_ssp585('p021',[1],0), quit"

    if nargin < 1 || isempty(X)
        error('X is required: ''p021'' | ''p022'' | ''p023''.');
    end
    if nargin < 2 || isempty(steps)
        steps = [1];
    end
    if nargin < 3 || isempty(loadonly)
        loadonly = 0;
    end

    % ---- PPE member lookup: X -> {K_variant, P_number} ---------------------
    switch X
        case 'p021'
            K_variant = 'mode';
            P_number  = 'P021';
        case 'p022'
            K_variant = '95th';
            P_number  = 'P022';
        case 'p023'
            K_variant = '5th';
            P_number  = 'P023';
        otherwise
            error('Unrecognised X ''%s'' -- expected ''p021''|''p022''|''p023''.', X);
    end
    % -----------------------------------------------------------------------

    % ---- swappable scenario/model label -- change these lines to adapt
    % this script to a different scenario or forcing model (matching this
    % project's existing swappable-label convention). FORCING_PREFIX is
    % the filename prefix that scenario's own ProjTF/ProjSMB steps used
    % when saving (CESM_WACCM_TF_*.mat etc.) AND (confirmed against both
    % the CESM2-WACCM and MRI-ESM2 production hist_run scripts) the
    % suffix on that model's own relaxed-state file
    % (AIS_ISMIP7_Relaxed_<FORCING_PREFIX>.mat) -- NOT a simple transform
    % of CMIP_MODEL. RELAX_MODEL_DIR is where that relaxed-state file
    % actually lives -- also model-specific.
    CMIP_MODEL      = 'MRI-ESM2-0';
    SCENARIO        = 'ssp585';
    FORCING_PREFIX  = 'MRI_ESM2';
    RELAX_MODEL_DIR = 'Models_MRIESM2';   % relative to init/ -- see note above
    HIST_DIR        = 'MRI-ESM2';   % hist_runs/ folder name -- NOT CMIP_MODEL ('MRI-ESM2-0') for this model
    % -----------------------------------------------------------------------

    addpath('./../../../init/scripts');
    addpath('./../../functions');

    modeldir = ['./Models_' P_number '/'];
    if ~exist(modeldir, 'dir'), mkdir(modeldir); end

    org = organizer('repository', modeldir, ...
                    'prefix', ['AIS_ISMIP7_PPE_' P_number '_'], ...
                    'steps', steps, 'color', '34;47;2');
    clear steps;

    % ------------------------------------------------------------------ paths
    proj_root  = './../../../';

    inputmodel_2015 = [proj_root 'hist_runs/' HIST_DIR '/Models/' ...
                       'AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat'];

    preproc_ocean      = [proj_root 'preprocessed_data/Ocean/'];
    preproc_proj_ocean = [proj_root 'preprocessed_data/Ocean/Proj/' CMIP_MODEL '/' SCENARIO '/'];
    preproc_proj_atmo  = [proj_root 'preprocessed_data/Atmosphere/Proj/' CMIP_MODEL '/' SCENARIO '/'];

    % ------------------------------------------------------------------ time
    start_year = 2015;
    mid_year   = 2150;   % phase 1 ends here (state saved at 2151) -- same
                          % split point as ../ssp585/old_scripts/
                          % proj_run_CESM_WACCM_ssp585_2015_2300.m
    end_year   = 2300;

    % ------------------------------------------------------------------ ISMIP6 outputs (same as core run)
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
    if perform(org, 'ProjRun_2015_2150') % {{{
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
        md.timestepping.final_time = mid_year + 1;   % last snapshot at 2151
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;

        md.transient.requested_outputs = ismip6_outputs;

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean: real ssp585 transient TF (NOT a climatology, unlike
        % ctrl), loaded directly from the core run's already-built file --
        % no rebuilding. K (gamma_0) is the ONE thing that varies here.
        % Trimmed to [start_year, mid_year+1] -- halves the .bin size
        % (see module docstring). ---
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_proj_ocean FORCING_PREFIX '_TF_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);   % -> tf_proj, z_data

        for di = 1:numel(tf_proj)
            c = tf_proj{di}; t_c = c(end, :);
            tf_proj{di} = c(:, t_c <= mid_year + 1);
        end

        unique_basinid = unique(basinid);
        tmp     = load([preproc_ocean 'dT_correction.mat'], 'dT_correction');
        delta_t = tmp.dT_correction;   % 1 x nBasins, from meltMip_ensemble.m step 8 (get_dT_iterate_BMB_j)

        % --- K (gamma_0): the PPE-varying parameter for p021/p022/p023.
        % All 4 percentile variants live in the SAME gamma0_local.mat,
        % built once by init/meltMip_ensemble.m's 'save_gamma0_local' step
        % from run_parameter_selection.py's own K selection -- nothing
        % re-derived here, just selecting which field to use. ---
        loaded_gamma0 = load([preproc_ocean 'gamma0_local.mat']);
        switch K_variant
            case 'mode'
                gamma0_this = loaded_gamma0.gamma0_local;        % SAME K as the core ssp585 run
            case '95th'
                gamma0_this = loaded_gamma0.gamma0_local_95th;
            case '5th'
                gamma0_this = loaded_gamma0.gamma0_local_5th;
        end
        clear loaded_gamma0;

        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique_basinid);
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = tf_proj;
        md.basalforcings.islocal    = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_this;

        % --- SMB: real ssp585 transient SMB, loaded directly from the
        % core run's already-built file -- no rebuilding, no PPE variation
        % here. Trimmed to [start_year, mid_year+1], same as TF above. ---
        load([preproc_proj_atmo FORCING_PREFIX '_SMB_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);   % -> smb_forcing, bgrad_forcing
        keep1         = smb_forcing(end, :) <= mid_year + 1;
        smb_forcing   = smb_forcing(:, keep1);
        bgrad_forcing = bgrad_forcing(:, keep1);

        md_relax  = loadmodel([proj_root 'init/' RELAX_MODEL_DIR '/AIS_ISMIP7_Relaxed_' FORCING_PREFIX '.mat']);
        surf_1995 = md_relax.geometry.surface;
        clear md_relax

        md.smb        = SMBgradients();
        md.smb.smbref = smb_forcing;
        md.smb.b_pos  = bgrad_forcing;
        md.smb.b_neg  = bgrad_forcing;
        md.smb.href   = [surf_1995 ; 1995];

        % --- Calving front: FROZEN at the AIS_state_2015 starting position
        % -- ismovingfront=0 means the level-set equation is never solved,
        % so ice front/grounding line evolve passively with ice dynamics
        % only. No ice-shelf collapse mask, no spclevelset at all (matches
        % ../ctrl/'s own frozen-front mechanism exactly). ---
        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 0;

        md.miscellaneous.name = ['ProjRun_PPE_' P_number '_' CMIP_MODEL '_' SCENARIO '_' ...
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

    % ================================================================= Step 2
    if perform(org, 'AIS_state_2151') % {{{
        % Save the end-state at 2151 (final snapshot of ProjRun_2015_2150)
        % as a clean restart for ProjRun_2151_2300, with TransientSolution
        % cleared to keep file size manageable -- same technique as
        % ../../CESM2-WACCM/ssp585/old_scripts/proj_run_CESM_WACCM_ssp585_2015_2300.m's
        % own step 5.
        md = loadmodel(org, 'ProjRun_2015_2150');
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
    if perform(org, 'ProjRun_2151_2300') % {{{
        % Transient 2151-2300, continuing from AIS_state_2151. Reloads the
        % same full TF/SMB files as step 1, but trims to the SECOND half
        % of the span this time.
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
        md.timestepping.interp_forcing  = 0;

        md.timestepping.start_time = mid_year + 1;   % 2151
        md.timestepping.final_time = end_year + 1;   % last snapshot at 2300
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;

        md.transient.requested_outputs = ismip6_outputs;

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_proj_ocean FORCING_PREFIX '_TF_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);   % -> tf_proj, z_data

        for di = 1:numel(tf_proj)
            c = tf_proj{di}; t_c = c(end, :);
            tf_proj{di} = c(:, t_c >= mid_year + 1);
        end

        unique_basinid = unique(basinid);
        tmp     = load([preproc_ocean 'dT_correction.mat'], 'dT_correction');
        delta_t = tmp.dT_correction;

        loaded_gamma0 = load([preproc_ocean 'gamma0_local.mat']);
        switch K_variant
            case 'mode'
                gamma0_this = loaded_gamma0.gamma0_local;
            case '95th'
                gamma0_this = loaded_gamma0.gamma0_local_95th;
            case '5th'
                gamma0_this = loaded_gamma0.gamma0_local_5th;
        end
        clear loaded_gamma0;

        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique_basinid);
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = tf_proj;
        md.basalforcings.islocal    = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_this;

        load([preproc_proj_atmo FORCING_PREFIX '_SMB_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);   % -> smb_forcing, bgrad_forcing
        keep3         = smb_forcing(end, :) >= mid_year + 1;
        smb_forcing   = smb_forcing(:, keep3);
        bgrad_forcing = bgrad_forcing(:, keep3);

        md_relax  = loadmodel([proj_root 'init/' RELAX_MODEL_DIR '/AIS_ISMIP7_Relaxed_' FORCING_PREFIX '.mat']);
        surf_1995 = md_relax.geometry.surface;   % same 1995 reference surface as segment 1
        clear md_relax

        md.smb        = SMBgradients();
        md.smb.smbref = smb_forcing;
        md.smb.b_pos  = bgrad_forcing;
        md.smb.b_neg  = bgrad_forcing;
        md.smb.href   = [surf_1995 ; 1995];

        md.calving.calvingrate         = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.transient.ismovingfront = 0;

        md.miscellaneous.name = ['ProjRun_PPE_' P_number '_' CMIP_MODEL '_' SCENARIO '_' ...
                                 num2str(mid_year + 1) '_' num2str(end_year)];
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
    if perform(org, 'WriteNetCDF') % {{{
        % Loads BOTH gathered segments and concatenates their
        % TransientSolution before gridding -- same technique as
        % ../../CESM2-WACCM/ssp585/old_scripts/proj_run_CESM_WACCM_ssp585_2015_2300.m's
        % own WriteISMIP6_NetCDF step. meta.set_counter = P_number
        % ('P021'/'P022'/'P023'), matching this group's PPE filename
        % convention directly (NOT the 'CXXX' counter used by ../ctrl/'s
        % own WriteISMIP6_NetCDF step -- that's a different experiment
        % family).
        addpath('./../../functions');

        md1 = loadmodel(org, 'ProjRun_2015_2150');
        md2 = loadmodel(org, 'ProjRun_2151_2300');
        md = md1;
        md.results.TransientSolution = [md1.results.TransientSolution, ...
                                         md2.results.TransientSolution];
        clear md1 md2;

        outdir = [proj_root 'postprocessed_data/' CMIP_MODEL '/' SCENARIO '_PPE_' P_number '/'];
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        meta                    = struct();
        meta.experiment_id      = SCENARIO;
        meta.set_counter        = P_number;
        meta.time_range         = [num2str(start_year) '-' num2str(end_year)];
        meta.ESM_id             = CMIP_MODEL;
        meta.forcing_member_id  = 'f001';
        meta.ISM_member_id      = 'm001';

        [cfflux_tot, glflux_tot] = write_ismip7_2d_projection(md, outdir, meta);
        write_ismip7_scalar_projection(md, outdir, meta, cfflux_tot, glflux_tot);
    end % }}}

end
