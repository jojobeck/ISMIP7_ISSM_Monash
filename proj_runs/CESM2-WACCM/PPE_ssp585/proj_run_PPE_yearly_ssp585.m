function proj_run_PPE_yearly_ssp585(X, start_year, end_year, waitonlock_seconds)
% PPE (Perturbed-Parameter Ensemble) run: yearly-restart loop, mirroring
% proj_run_CESM_WACCM_ssp585_2015_2300_yearly.m's own yearly-loop step
% exactly (same floating-only gating mechanism, same 2015-no-advance
% seed, same resume-from-crash design -- see that script's header for the
% full rationale, not repeated here). ssp585 counterpart of
% ../PPE_ssp370/proj_run_PPE_yearly_ssp370.m -- same mechanism/axes,
% different scenario/run span. Covers PPE members that need the per-year
% floating-gated collapse-mask loop (unlike p011/p013's frozen front,
% these all evolve the calving front). P-numbers continue on from the
% ssp370 PPE family (P001-P009) rather than restarting at P004, so every
% PPE member across every scenario folder has its own unique P-number --
% this script's own members are P014-P019:
%
%   X='p014'  migration_max = Greene MAXIMUM (96521.91 m/yr), no eigencalving
%   X='p015'  migration_max = Greene MINIMUM (203.82 m/yr),   no eigencalving
%   X='p016'  migration_max = Greene MEAN (2863.78 m/yr, core value), + eigencalving
%   X='p017'  migration_max = Greene MAXIMUM,                          + eigencalving
%   X='p018'  migration_max = Greene MINIMUM,                          + eigencalving
%   X='p019'  migration_max = Greene MEAN (core value), no eigencalving,
%             SMB surface-elevation feedback SWITCHED OFF (b_pos=b_neg=0)
%
% Each X internally resolves {migration_max, eigencalving on/off,
% smb_feedback on/off} AND this PPE member's own P-number (P014-P019).
%
% EIGENCALVING (p016-p018): md.calving is REPLACED with calvinglevermann
% (coeff below), on top of the SAME floating-gated collapse-mask forcing
% used everywhere else -- these are independent, additive mechanisms in
% ISSM's level-set solver (front velocity w = v - c - m, with c now the
% eigencalving-law-driven rate instead of zero; spclevelset is a separate
% hard override at specific nodes/times, unaffected by the calving law).
% coeff = 6.59e15 m*s, derived from Wilner et al. (2023)'s own
% Levermann-eigencalving calibration across 10 Antarctic ice shelves
% (K in units of 1e8 m*yr each; mean EXCLUDING one outlier shelf at
% K=300e8 -- 9 remaining values 0.5,0.3,1.5,3,4.0,3,2,2,2.5, sum=18.8,
% mean=18.8/9=2.0889e8 m*yr):
%   coeff[m*s] = K_wilner[m*yr] * sec_to_year[s/yr]
%              = 2.0889e8 * 31556926 ~= 6.592e15
% (verified directly against ISSM's own C++ implementation,
% Tria::CalvingRateLevermann in src/c/classes/Elements/Tria.cpp:
% calvingrate = coeff * strainparallel * strainperpendicular, no further
% internal unit scaling -- confirming coeff must be in m*s to convert the
% strain-rate-squared product [s^-2] into a velocity [m/s]). CAVEAT:
% Wilner et al. calibrated at ~1 km resolution at the calving front; this
% project's mesh is Antarctic-wide and coarser, so peak strain rates at
% shear margins/rifts are systematically smoothed/under-estimated here --
% applying this coefficient verbatim likely UNDER-produces the calving
% rate Wilner et al. intended. Treat p017-p019 as "sensitivity to adding
% eigencalving at a literature-representative coefficient," NOT as a
% coefficient calibrated for this mesh.
%
% NO PREPROCESSING -- loads the ALREADY-BUILT TF/SMB/levelset .mat files
% the core ssp585 yearly run produced (confirmed present on disk):
%   preprocessed_data/Ocean/Proj/CESM2-WACCM/ssp585/CESM_WACCM_TF_ssp585_2015_2300.mat
%   preprocessed_data/Atmosphere/Proj/CESM2-WACCM/ssp585/CESM_WACCM_SMB_ssp585_2015_2300.mat
%   preprocessed_data/Ocean/Proj/CESM2-WACCM/ssp585/CESM_WACCM_levelset_v2_ssp585_2015_2300.mat
% NOTE the mask VERSION differs from the ssp370 PPE script: the ssp585
% core run's own collapse mask is 'v2' (variable proj_spclevelset_v2,
% confirmed in proj_run_CESM_WACCM_ssp585_2015_2300_yearly.m), NOT 'v21'
% like ssp370's -- this is a genuine per-scenario difference, not a typo.
%
% Each PPE member's per-year models are saved to its OWN dedicated
% ./Models_P<NNN>/ repository (e.g. ./Models_P014/ for p014), NOT the
% core run's shared Models_yearly/ folder.
%
% start_year need not be 2015 -- if start_year>2015, this RESUMES a prior
% run of this same PPE member from that year (e.g. after a crash at 2080:
% rerun with the SAME X and start_year=2080), exactly like the core
% yearly script's own resume mechanism.
%
% Usage (start_year/end_year/waitonlock_seconds default to 2015/2300/3600
% -- end_year default matches ssp585's own 2015-2300 span, NOT ssp370's
% 2015-2100):
%   proj_run_PPE_yearly_ssp585('p014')                       % full run, defaults
%   proj_run_PPE_yearly_ssp585('p017', 2015, 2015, 3600)      % 1-yr test
%   proj_run_PPE_yearly_ssp585('p018', 2080, 2300, 6*3600)     % resume after a crash at 2080
%
% Run with ISSM already on the MATLAB path (devpath already called), e.g.:
%   matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
%     addpath('$ISSM_DIR/lib'); proj_run_PPE_yearly_ssp585('p014',2015,2015,3600), quit"

    if nargin < 1 || isempty(X)
        error('X is required: ''p014''|''p015''|''p016''|''p017''|''p018''|''p019''.');
    end
    if nargin < 2 || isempty(start_year)
        start_year = 2015;
    end
    if nargin < 3 || isempty(end_year)
        end_year = 2300;
    end
    if nargin < 4 || isempty(waitonlock_seconds)
        waitonlock_seconds = 3600;  % default 1 hour
    end

    % ---- PPE member lookup: X -> {P_number, migration_max, eigencalving, smb_feedback} ---
    MIG_MAX_GREENE_MEAN = 2863.78;    % core value, same as the core ssp585 run
    MIG_MAX_GREENE_MAX  = 96521.91;
    MIG_MAX_GREENE_MIN  = 203.82;
    EIGENCALVING_COEFF  = 6.592e15;   % m*s -- see module docstring for derivation

    switch X
        case 'p014'
            P_number = 'P014'; migration_max = MIG_MAX_GREENE_MAX;  use_eigencalving = false; smb_feedback_on = true;
        case 'p015'
            P_number = 'P015'; migration_max = MIG_MAX_GREENE_MIN;  use_eigencalving = false; smb_feedback_on = true;
        case 'p016'
            P_number = 'P016'; migration_max = MIG_MAX_GREENE_MEAN; use_eigencalving = true;  smb_feedback_on = true;
        case 'p017'
            P_number = 'P017'; migration_max = MIG_MAX_GREENE_MAX;  use_eigencalving = true;  smb_feedback_on = true;
        case 'p018'
            P_number = 'P018'; migration_max = MIG_MAX_GREENE_MIN;  use_eigencalving = true;  smb_feedback_on = true;
        case 'p019'
            P_number = 'P019'; migration_max = MIG_MAX_GREENE_MEAN; use_eigencalving = false; smb_feedback_on = false;
        otherwise
            error('Unrecognised X ''%s'' -- expected ''p014''|''p015''|''p016''|''p017''|''p018''|''p019''.', X);
    end
    % -----------------------------------------------------------------------

    % ---- swappable scenario/model label -- see proj_run_PPE_frozenfront_ssp585.m's
    % own header note for why FORCING_PREFIX/RELAX_MODEL_DIR are NOT simple
    % transforms of CMIP_MODEL and must be set explicitly per model.
    CMIP_MODEL      = 'CESM2-WACCM';
    SCENARIO        = 'ssp585';
    FORCING_PREFIX  = 'CESM_WACCM';
    RELAX_MODEL_DIR = 'Models';   % relative to init/
    % -----------------------------------------------------------------------

    proj_root  = './../../../';
    init_dir   = [proj_root 'init/'];

    addpath('./../../../init/scripts');
    addpath('./../../local_matlab');   % shared waitonlock.m override -- see the
                                        % core yearly script's own header for why

    inputmodel_2015 = [proj_root 'hist_runs/' CMIP_MODEL '/Models/' ...
                       'AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat'];

    preproc_ocean      = [proj_root 'preprocessed_data/Ocean/'];
    preproc_proj_ocean = [proj_root 'preprocessed_data/Ocean/Proj/' CMIP_MODEL '/' SCENARIO '/'];
    preproc_proj_atmo  = [proj_root 'preprocessed_data/Atmosphere/Proj/' CMIP_MODEL '/' SCENARIO '/'];

    modeldir = ['./Models_' P_number '/'];
    if ~exist(modeldir, 'dir'), mkdir(modeldir); end

    % Same requested_outputs list as the core ssp585 yearly script.
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

    % ---- Load the FULL pre-built forcing/basin/mask arrays ONCE (no
    % preprocessing -- these are the core ssp585 run's own already-built
    % files). ----
    load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);  % -> basinid
    load([preproc_ocean 'tf_depths.mat']);                                 % -> tf_depths
    load([preproc_ocean 'gamma0_local.mat']);                              % -> gamma0_local (mode -- PPE K is not varied in this script)
    tmp = load([preproc_ocean 'dT_correction.mat'], 'dT_correction');
    delta_t = tmp.dT_correction;
    unique_basinid = unique(basinid);

    load([preproc_proj_ocean FORCING_PREFIX '_TF_' SCENARIO '_2015_2300.mat']);         % -> tf_proj, z_data
    load([preproc_proj_atmo  FORCING_PREFIX '_SMB_' SCENARIO '_2015_2300.mat']);        % -> smb_forcing, bgrad_forcing

    % NOTE: 'v2', not 'v21' -- see module docstring.
    loaded_lset = load([preproc_proj_ocean FORCING_PREFIX '_levelset_v2_' SCENARIO '_2015_2300.mat']);
    proj_spclevelset = loaded_lset.proj_spclevelset_v2;
    clear loaded_lset;

    md_relax  = loadmodel([init_dir RELAX_MODEL_DIR '/AIS_ISMIP7_Relaxed_' FORCING_PREFIX '.mat']);
    surf_1995 = md_relax.geometry.surface;
    clear md_relax

    % ---- Restart state (identical mechanism to the core yearly script). ----
    if start_year == 2015
        md = loadmodel(inputmodel_2015);
        applied_collapse = md.mask.ice_levelset > 0;
        fprintf('Seeded applied_collapse with %d vertices already ice-free at 2015 (of %d total).\n', ...
                sum(applied_collapse), md.mesh.numberofvertices);
    else
        resume_model_fname = [modeldir 'AIS_ISMIP7_PPE_' P_number '_yearly_' num2str(start_year) '.mat'];
        resume_ac_fname    = [modeldir 'AIS_ISMIP7_PPE_' P_number '_yearly_appliedcollapse_' num2str(start_year) '.mat'];
        if ~isfile(resume_model_fname)
            error('Resume requested from year %d but %s does not exist.', start_year, resume_model_fname);
        end
        if ~isfile(resume_ac_fname)
            error('Resume requested from year %d but %s does not exist.', start_year, resume_ac_fname);
        end

        fprintf('Resuming from year %d: %s\n', start_year, resume_model_fname);
        loaded = load(resume_model_fname, 'md');
        md = loaded.md;
        clear loaded;

        md.geometry.thickness   = md.results.TransientSolution(end).Thickness;
        md.geometry.surface     = md.results.TransientSolution(end).Surface;
        md.geometry.base        = md.results.TransientSolution(end).Base;
        md.mask.ocean_levelset  = md.results.TransientSolution(end).MaskOceanLevelset;
        md.mask.ice_levelset    = md.results.TransientSolution(end).MaskIceLevelset;
        md.results.TransientSolution = [];

        loaded_ac = load(resume_ac_fname, 'applied_collapse');
        applied_collapse = loaded_ac.applied_collapse;
        clear loaded_ac;
        fprintf('Restored applied_collapse: %d vertices already forced to ocean.\n', sum(applied_collapse));
    end

    nVerts = md.mesh.numberofvertices;

    for year = start_year:end_year
        fprintf('\n=== [%s] Year %d -> %d ===\n', X, year, year + 1);

        m = ((1 + sin(71*pi/180)) * ones(nVerts, 1) ...
             ./ (1 + sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor = (1./m).^2;

        md.inversion.iscontrol       = 0;
        md.transient.isthermal       = 0;
        md.transient.isgroundingline = 1;
        md.transient.ismasstransport = 1;
        md.transient.isstressbalance = 1;
        md.masstransport.spcthickness   = NaN*ones(nVerts, 1);
        md.outputdefinition.definitions = {};
        md.timestepping.interp_forcing  = 0;

        md.timestepping.start_time   = year;
        md.timestepping.final_time   = year + 1;
        md.timestepping.time_step    = 1/12.;
        md.settings.output_frequency = 12;

        md.transient.requested_outputs = ismip6_outputs;

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean: slice the full TF series to this 1-year window (same
        % backward-looking "most recent at-or-before year" logic as the
        % core yearly script). K (gamma_0) is NOT varied in this script --
        % that is p011-p013's own axis (proj_run_PPE_frozenfront_ssp585.m). ---
        tf_year = tf_proj;
        for di = 1:numel(tf_year)
            c = tf_year{di};
            t_c = c(end, :);
            idx = find(t_c <= year, 1, 'last');
            tf_year{di} = c(:, idx);
        end

        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique_basinid);
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = tf_year;
        md.basalforcings.islocal    = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_local;

        % --- SMB: same "most recent at-or-before year" logic as TF.
        % smb_feedback_on=false (p019 only): elevation-feedback terms
        % (b_pos/b_neg) zeroed -- smbref alone determines mass balance,
        % with no lapse-rate correction for (surface - href). ---
        idx_smb    = find(smb_forcing(end, :) <= year, 1, 'last');
        smb_year   = smb_forcing(:, idx_smb);
        if smb_feedback_on
            bgrad_year = bgrad_forcing(:, idx_smb);
        else
            bgrad_year = zeros(nVerts, 1);
        end

        md.smb        = SMBgradients();
        md.smb.smbref = smb_year;
        md.smb.b_pos  = bgrad_year;
        md.smb.b_neg  = bgrad_year;
        md.smb.href   = [surf_1995 ; 1995];

        % --- Calving front: collapse-mask forcing GATED to floating ice
        % only (identical mechanism to the core yearly script -- see that
        % script's in-loop comments for the full rationale). ---
        is_floating_now = (md.mask.ice_levelset < 0) & (md.mask.ocean_levelset < 0);

        keep_lset = (proj_spclevelset(end, :) >= year) & (proj_spclevelset(end, :) <= year + 1);
        raw_spclevelset_year = proj_spclevelset(:, keep_lset);
        time_row_year = raw_spclevelset_year(end, :);
        raw_forced_this_year = any(raw_spclevelset_year(1:end-1, :) == 1, 2);

        applied_collapse = applied_collapse | (raw_forced_this_year & is_floating_now);

        gated_col = NaN(nVerts + 1, 1);
        gated_col(applied_collapse) = 1;
        md.levelset.spclevelset = repmat(gated_col, 1, size(raw_spclevelset_year, 2));
        md.levelset.spclevelset(end, :) = time_row_year;

        % --- Calving law: eigencalving (p016-p018) ADDS a real
        % strain-rate-driven calving rate ON TOP of the same collapse-mask
        % forcing above -- these are independent, additive mechanisms (see
        % module docstring). use_eigencalving=false (p014/p015/p019):
        % same zero-calvingrate no-op as the core yearly script -- the
        % collapse mask is the only thing moving the front. ---
        if use_eigencalving
            md.calving       = calvinglevermann();
            md.calving.coeff = EIGENCALVING_COEFF * ones(nVerts, 1);
        else
            md.calving.calvingrate = zeros(nVerts, 1);
        end
        md.frontalforcings.meltingrate = zeros(nVerts, 1);
        md.transient.ismovingfront     = 1;
        md.levelset.migration_max      = migration_max;   % m/yr -- this PPE member's own value

        md.miscellaneous.name = ['ProjRunYearly_PPE_' P_number '_' CMIP_MODEL '_' SCENARIO '_' ...
                                 num2str(year) '_' num2str(year + 1)];
        md.cluster             = set_cluster('gadi');
        md.settings.waitonlock = waitonlock_seconds;   % SECONDS -- local waitonlock.m override
        md.verbose = verbose('solution', true, 'module', true, 'convergence', true);

        try
            md = solve(md, 'tr', 'runtimename', false, 'loadonly', 0);
        catch ME
            fprintf('\n*** [%s] Year %d -> %d FAILED: %s\n', X, year, year + 1, ME.message);
            fprintf('*** Stopping loop. Investigate, then resume from year %d.\n', year);
            rethrow(ME);
        end

        yearly_fname = [modeldir 'AIS_ISMIP7_PPE_' P_number '_yearly_' num2str(year + 1) '.mat'];
        save(yearly_fname, 'md', '-v7.3');
        fprintf('Saved: %s\n', yearly_fname);

        ac_fname = [modeldir 'AIS_ISMIP7_PPE_' P_number '_yearly_appliedcollapse_' num2str(year + 1) '.mat'];
        save(ac_fname, 'applied_collapse', '-v7.3');
        fprintf('Saved: %s\n', ac_fname);

        md_solved = md;
        md.geometry.thickness   = md_solved.results.TransientSolution(end).Thickness;
        md.geometry.surface     = md_solved.results.TransientSolution(end).Surface;
        md.geometry.base        = md_solved.results.TransientSolution(end).Base;
        md.mask.ocean_levelset  = md_solved.results.TransientSolution(end).MaskOceanLevelset;
        md.mask.ice_levelset    = md_solved.results.TransientSolution(end).MaskIceLevelset;
        md.results.TransientSolution = [];
        clear md_solved;
    end

    fprintf('\n=== [%s] complete: %d -> %d ===\n', X, start_year, end_year);
end
