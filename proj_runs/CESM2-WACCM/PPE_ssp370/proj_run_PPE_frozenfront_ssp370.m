function md = proj_run_PPE_frozenfront_ssp370(X, steps, loadonly)
% PPE (Perturbed-Parameter Ensemble) run: level-set FROZEN at present-day
% (2015) ice extent, to isolate BMB-driven retreat from calving-front
% evolution entirely -- no ice-shelf collapse mask, no spclevelset time
% series, no floating-only gating. Covers PPE members (renumbered
% 2026-10-05: there is no separate P000 for the core/baseline run -- it
% has no P-number at all -- so numbering starts at P001, not P002; the
% 'mNNN' input-argument prefix is also retired in favour of 'pNNN',
% matching the yearly script's own convention -- there is no longer any
% 'mNNN'/'MNNN' naming anywhere in this PPE family):
%   X='p001'  K = mode              (gamma0_local)      -- SAME K as core ssp370 run
%   X='p002'  K = 95th percentile   (gamma0_local_95th)
%   X='p003'  K = 5th percentile    (gamma0_local_5th)
% (gamma0_local.mat's 4 variants -- mode/5th/50th/95th -- were built by
% init/meltMip_ensemble.m's 'save_gamma0_local' step from
% run_parameter_selection.py's own K selection; nothing here re-derives
% them.) Each X internally resolves both the K variant AND this PPE
% member's output label (P001/P002/P003 respectively).
%
% Mechanism: md.transient.ismovingfront = 0, matching
% ../ctrl/proj_run_CESM_WACCM_ctrl_2015_2300.m's own frozen-front
% approach exactly (see that script's header for the full rationale) --
% the level-set equation is never solved, so the ice front/grounding line
% stay exactly where they are in the AIS_state_2015 starting state for
% the whole run. UNLIKE ctrl (which additionally uses a fixed 2000-2029
% CLIMATOLOGICAL forcing, to isolate ice-sheet-model drift from climate
% forcing), this script uses the REAL ssp370 TRANSIENT TF/SMB forcing
% (same as the core ssp370 run) -- only the calving front is frozen, so
% any mass change here is attributable to BMB (via the varying K) and SMB
% alone, with calving-front retreat/advance eliminated as a factor.
%
% NO PREPROCESSING -- loads the ALREADY-BUILT TF/SMB .mat files the core
% ssp370 run produced (confirmed present on disk):
%   preprocessed_data/Ocean/Proj/CESM2-WACCM/ssp370/CESM_WACCM_TF_ssp370_2015_2100.mat
%   preprocessed_data/Atmosphere/Proj/CESM2-WACCM/ssp370/CESM_WACCM_SMB_ssp370_2015_2100.mat
% Single continuous 2015-2100 transient (no restart split needed -- like
% ctrl, there's nothing a restart would help with here since the front
% never moves and the loaded forcing already spans the full run).
%
% Calling convention: (X, steps, loadonly) -- X picks the PPE member,
% steps picks which organizer step(s) of THIS run to execute, loadonly is
% the standard two-phase submit(0)/gather(1) flag used ONLY within step 1.
%
% Step map (organizer, per PPE member -- same shape as
% ../ctrl/proj_run_CESM_WACCM_ctrl_2015_2300.m's own ProjRun/WriteISMIP6_NetCDF
% split):
%   1  ProjRun       solve the 2015-2100 transient (loadonly=0 submit,
%                     =1 gather + save the model -- NO NetCDF writing here)
%   2  WriteNetCDF    loads the already-saved model from step 1 and writes
%                     the ISMIP6 NetCDF output -- run this ONLY after
%                     step 1's gather has completed and saved successfully
%
% Output: this member's P-number (e.g. 'P001') is used as meta.set_counter,
% matching this group's PPE filename convention directly, e.g.
%   libmassbffl_AIS_Monash_ISSM_m001_CESM2-WACCM_f001_ssp370_P001_2015-2100.nc
% (ISM_member_id stays fixed at 'm001' -- same physical ice-sheet model
% for every PPE member, an UNRELATED fixed ISMIP7-filename-convention
% field, not this script's own X/P-number -- the P-number is what varies
% per experiment, not the ISM_member_id or the ESM_id).
%
% Each PPE member's solved model is saved to its OWN dedicated
% ./Models_P<NNN>/ repository (e.g. ./Models_P001/ for p001), NOT a
% shared Models/ or Models_yearly/ folder -- keeps every PPE variant's
% output completely separate on disk.
%
% Usage:
%   proj_run_PPE_frozenfront_ssp370('p001', [1], 0)   % submit p001 (K=mode)
%   proj_run_PPE_frozenfront_ssp370('p001', [1], 1)   % gather p001 (save model only, no NetCDF)
%   proj_run_PPE_frozenfront_ssp370('p001', [2])      % write NetCDF for p001 (after gather above)
%   proj_run_PPE_frozenfront_ssp370('p002', [1], 0)   % submit p002 (K=95th)
%   proj_run_PPE_frozenfront_ssp370('p003', [1], 0)   % submit p003 (K=5th)
%
% Run with ISSM already on the MATLAB path (devpath already called), e.g.:
%   matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
%     addpath('$ISSM_DIR/lib'); proj_run_PPE_frozenfront_ssp370('p001',[1],0), quit"

    if nargin < 1 || isempty(X)
        error('X is required: ''p001'' | ''p002'' | ''p003''.');
    end
    if nargin < 2 || isempty(steps)
        steps = [1];
    end
    if nargin < 3 || isempty(loadonly)
        loadonly = 0;
    end

    % ---- PPE member lookup: X -> {K_variant, P_number} ---------------------
    switch X
        case 'p001'
            K_variant = 'mode';
            P_number  = 'P001';
        case 'p002'
            K_variant = '95th';
            P_number  = 'P002';
        case 'p003'
            K_variant = '5th';
            P_number  = 'P003';
        otherwise
            error('Unrecognised X ''%s'' -- expected ''p001''|''p002''|''p003''.', X);
    end
    % -----------------------------------------------------------------------

    % ---- swappable scenario/model label -- change these lines to adapt
    % this script to a different scenario or forcing model (matching this
    % project's existing swappable-label convention). end_year must match
    % that scenario's own real forcing span (2300 for ssp585/ssp126, 2100
    % for ssp370 -- confirmed on disk per each scenario's own proj_run
    % script header). FORCING_PREFIX is the filename prefix that
    % scenario's own ProjTF/ProjSMB steps used when saving
    % (CESM_WACCM_TF_*.mat etc.) AND (confirmed against both the
    % CESM2-WACCM and MRI-ESM2 production hist_run scripts) the suffix on
    % that model's own relaxed-state file (AIS_ISMIP7_Relaxed_<FORCING_PREFIX>.mat)
    % -- NOT a simple transform of CMIP_MODEL (e.g. MRI-ESM2-0 uses
    % 'MRI_ESM2' for both, not 'MRI-ESM2-0' or 'MRI_ESM2_0'). RELAX_MODEL_DIR
    % is where that relaxed-state file actually lives -- also model-specific
    % (CESM2-WACCM's sits in the generic init/Models/; MRI-ESM2's own sits
    % in a per-model init/Models_MRIESM2/ -- confirmed different, not
    % derivable from CMIP_MODEL either).
    CMIP_MODEL      = 'CESM2-WACCM';
    SCENARIO        = 'ssp370';
    FORCING_PREFIX  = 'CESM_WACCM';
    RELAX_MODEL_DIR = 'Models';   % relative to init/ -- see note above
    end_year        = 2100;
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

    inputmodel_2015 = [proj_root 'hist_runs/' CMIP_MODEL '/Models/' ...
                       'AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat'];

    preproc_ocean      = [proj_root 'preprocessed_data/Ocean/'];
    preproc_proj_ocean = [proj_root 'preprocessed_data/Ocean/Proj/' CMIP_MODEL '/' SCENARIO '/'];
    preproc_proj_atmo  = [proj_root 'preprocessed_data/Atmosphere/Proj/' CMIP_MODEL '/' SCENARIO '/'];

    start_year = 2015;

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
    if perform(org, 'ProjRun') % {{{
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
        md.timestepping.final_time = end_year + 1;
        md.timestepping.time_step  = 1/12.;
        md.settings.output_frequency = 12;

        md.transient.requested_outputs = ismip6_outputs;

        md.groundingline.migration              = 'SubelementMigration';
        md.groundingline.friction_interpolation = 'SubelementFriction1';
        md.groundingline.melt_interpolation     = 'SubelementMelt1';

        % --- Ocean: real ssp370 transient TF (NOT a climatology, unlike
        % ctrl), loaded directly from the core run's already-built file --
        % no rebuilding. K (gamma_0) is the ONE thing that varies here. ---
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);
        load([preproc_ocean 'tf_depths.mat']);
        load([preproc_proj_ocean FORCING_PREFIX '_TF_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);   % -> tf_proj, z_data

        unique_basinid = unique(basinid);
        tmp     = load([preproc_ocean 'dT_correction.mat'], 'dT_correction');
        delta_t = tmp.dT_correction;   % 1 x nBasins, from meltMip_ensemble.m step 8 (get_dT_iterate_BMB_j)

        % --- K (gamma_0): the PPE-varying parameter for p001/p002/p003.
        % All 4 percentile variants live in the SAME gamma0_local.mat,
        % built once by init/meltMip_ensemble.m's 'save_gamma0_local' step
        % from run_parameter_selection.py's own K selection -- nothing
        % re-derived here, just selecting which field to use. ---
        loaded_gamma0 = load([preproc_ocean 'gamma0_local.mat']);
        switch K_variant
            case 'mode'
                gamma0_this = loaded_gamma0.gamma0_local;        % SAME K as the core ssp370 run
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

        % --- SMB: real ssp370 transient SMB, loaded directly from the
        % core run's already-built file -- no rebuilding, no PPE variation
        % here (SMB elevation-feedback is a SEPARATE PPE axis, M010 in the
        % yearly-loop script, not this frozen-front script). ---
        load([preproc_proj_atmo FORCING_PREFIX '_SMB_' SCENARIO '_' ...
              num2str(start_year) '_' num2str(end_year) '.mat']);   % -> smb_forcing, bgrad_forcing

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

    % ================================================================= Step 2
    if perform(org, 'WriteNetCDF') % {{{
        % Loads the already-gathered model from step 1 (ProjRun) and writes
        % the ISMIP6 NetCDF output -- kept as its own separate step (NOT run
        % automatically inside step 1's loadonly=1 gather) so the gather can
        % be checked/rerun independently of the NetCDF writing. Grids the
        % TransientSolution onto the standard ISMIP7 761x761 8 km AIS grid;
        % meta.set_counter = P_number ('P001'/'P002'/'P003'), matching this
        % group's PPE filename convention directly (NOT the 'CXXX' counter
        % used by ../ctrl/'s own WriteISMIP6_NetCDF step -- that's a
        % different experiment family).
        addpath('./../../functions');

        md = loadmodel(org, 'ProjRun');

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
