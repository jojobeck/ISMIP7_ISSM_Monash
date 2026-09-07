function md = proj_run_CESM_WACCM_ssp585_2015_2300_yearly(steps, start_year, end_year, waitonlock_seconds)
% Standalone yearly-restart variant of proj_run_CESM_WACCM_ssp585_2015_2300.m.
%
% NEW FILE -- the original script is NOT modified. Unlike the original
% (which the annual loop depends on for its TF/SMB/levelset preprocessing
% .mat files), THIS script can build all of its own inputs, so it does not
% require the original script to have been run first.
%
% Step map (steps=[4] by default -- i.e. assumes 1-3 already run/available):
%   1  ProjTF          -- build annual TF .mat 2015-2300, copied verbatim
%                          from the original script's own ProjTF step (read
%                          fresh from source before writing this).
%   2  ProjSMB         -- build annual SMB+gradient .mat 2015-2300, copied
%                          verbatim from the original script's own ProjSMB
%                          step.
%   3  ProjLevelsetV2  -- build the cumulative collapse-mask spclevelset,
%                          but from the v2 collapse-mask NetCDF (not v1) --
%                          same build_spclevelset() logic already written
%                          and verified in build_and_compare_levelset_v2_ssp585.m,
%                          copied here so this script needs nothing external.
%                          Produces the SAME output file
%                          (CESM_WACCM_levelset_v2_ssp585_2015_2300.mat) that
%                          script's own step 1 produces -- interchangeable,
%                          not a second copy with a different name.
%   4  ProjRunYearly   -- the annual loop -- UNCHANGED logic from the
%                          previous version of this script, except it now
%                          loads the v2 levelset file (step 3's output)
%                          instead of v1's.
%
% Steps 1-3 ALWAYS build the FULL 2015-2300 span regardless of what
% start_year/end_year this invocation's loop (step 4) will actually use --
% they are meant to be run ONCE and reused across many step-4 invocations
% (different test ranges, resumes after a crash, etc.), exactly like the
% original script's own preprocessing steps. Steps 1-2 are NOT changed
% from the original at all (TF/SMB have nothing to do with v1/v2 -- only
% the collapse mask does); only step 3 differs, and only in which raw
% NetCDF it reads.
%
% start_year need not be 2015 for step 4. If start_year>2015, this RESUMES
% a prior run of this same script's step 4 from that year, e.g. after a
% crash at 2100: rerun with start_year=2100. Resuming loads both the
% solved model AND applied_collapse (the persistent floating-only-gating
% record, see below) from what a prior run already saved for that year --
% see the restart-state block for exactly what's required to be on disk.
%
% Calving-front forcing (step 4): the raw v2 collapse mask (proj_spclevelset_v2,
% built by step 3) is GATED every year to floating ice only, using the
% model's own current simulated MaskIceLevelset/MaskOceanLevelset -- a
% vertex is only forced to ocean once it is both proposed by the raw mask
% AND actually floating in the model's own state at that point. See the
% in-loop comments for the full mechanism and why the gating decision must
% persist (applied_collapse) rather than being re-evaluated fresh every year.
%
% applied_collapse is also SEEDED at the very start (start_year==2015)
% with every vertex that already had no ice at 2015 -- since
% md.calving.calvingrate and md.frontalforcings.meltingrate are both zero
% everywhere (calving fully prescribed by the mask), nothing opposes
% forward ice-velocity advection at the front outside the mask's own
% footprint, so without this seed the front can advance indefinitely into
% open ocean over 285 years (observed happening, e.g. in the Ross sector).
% This holds the 2015 extent -- the one point in the run backed by a real
% observation (the historical run's Greene-driven end-state) -- as a
% permanent ceiling the front can never cross. IMPORTANT: this seed only
% takes effect for a run that actually STARTS at 2015 -- if some years
% were already solved by an earlier run of this script without the seed,
% resuming from a later year inherits whatever advance already happened
% in those completed years (baked into their geometry, not just an
% unconstrained flag) -- a full correction requires restarting from 2015.
%
% WHY THIS NEEDS A LOCAL waitonlock.m OVERRIDE:
%   ISSM's generic src/m/solve/waitonlock.m polls for job completion over
%   SSH whenever oshostname() ~= cluster.name -- which is exactly the case
%   when this script is itself already running inside a PBS job on a Gadi
%   compute node. Self-SSH from a Gadi compute node back to itself returns
%   Signal 127, so the wait fails. ISSM already has an unused fix for this
%   (src/m/os/waitonlock_ongadi.m, isfile()-based, no SSH) but nothing
%   calls it -- solve() always calls the generic, broken one by name. This
%   script instead uses a LOCAL copy (proj_runs/local_matlab/waitonlock.m,
%   shared across all models/scenarios, same logic, renamed to shadow the
%   trunk version) added to the path AFTER devpath has already run, so
%   plain calls to waitonlock(md) -- including solve()'s own internal call
%   -- resolve to the working version. Nothing in the shared trunk install
%   is modified. Only relevant to step 4.
%
%   md.settings.waitonlock is in SECONDS with this override (matching
%   waitonlock_ongadi.m's convention) -- NOT minutes, which is what the
%   trunk's generic waitonlock.m uses for the same field. Kept FINITE
%   (never Inf) deliberately: if a submitted job crashes or stalls, the
%   wait throws a real MATLAB error after waitonlock_seconds rather than
%   hanging the whole multi-year loop forever. That error is caught below
%   and the loop stops with the failing year clearly logged, rather than
%   silently continuing or retrying.
%
% MIGRATION_MAX (step 4): set to 2863.78 m/yr -- the front-filtered (see
% greene_icemask_retreat_rate.py's front_mask()) pooled MEAN retreat rate
% from the Greene 1997-2021 observed ice-front record, deliberately not
% the raw/unfiltered pooled max (~96,522 m/yr, dominated by a single
% catastrophic ice-shelf-disintegration event and, more importantly, not
% representative of a cap meant to bind potentially every year given the
% known chronic grounding-line forcing issue this scenario has) nor the
% front-filtered max (~96,522 m/yr also, since that one genuine collapse
% event survived the front filter -- see conversation). ISSM's
% levelset.m marshal step scales this by 1/yts internally, so this is a
% direct m/yr assignment, no manual unit conversion needed.
%
% end_year is INCLUSIVE -- the last calendar year to simulate, matching
% the ORIGINAL script's own end_year=2300 meaning exactly (there,
% final_time=end_year+1 produces the calendar-year-end_year output
% snapshot). A full run therefore still uses end_year=2300, producing 286
% iterations (2015 through 2300 inclusive) -- the same total simulated
% span as the original's two segments combined (final_time=2301).
%
% Real TF/SMB forcing data only exists through 2299 (decade/annual source
% files stop there -- see the original script's ProjTF/ProjSMB comments).
% Every iteration takes the single most recent available forcing entry
% at-or-before that iteration's own start year, so iterations at or
% beyond 2299 correctly keep reusing 2299's real value (interp_forcing=0
% step-function hold) rather than going empty -- see in-loop comments.
%
% Usage (start_year/end_year/waitonlock_seconds default to 2015/2300/3600):
%   proj_run_CESM_WACCM_ssp585_2015_2300_yearly([1])                  % build TF only
%   proj_run_CESM_WACCM_ssp585_2015_2300_yearly([2])                  % build SMB only
%   proj_run_CESM_WACCM_ssp585_2015_2300_yearly([3])                  % build v2 levelset only
%   proj_run_CESM_WACCM_ssp585_2015_2300_yearly([1 2 3])              % build all three
%   proj_run_CESM_WACCM_ssp585_2015_2300_yearly([4], 2015, 2015, 3600)   % 1-yr test, 1 hr timeout
%   proj_run_CESM_WACCM_ssp585_2015_2300_yearly([4], 2015, 2300, 6*3600) % full run, 6 hr timeout/yr
%   proj_run_CESM_WACCM_ssp585_2015_2300_yearly([4], 2100, 2300, 6*3600) % resume after a crash at 2100
%
% Run with ISSM already on the MATLAB path (devpath already called), e.g.:
%   matlab -nodisplay -nosplash -r "addpath('$ISSM_DIR/src/m/dev'); devpath; \
%     addpath('$ISSM_DIR/lib'); proj_run_CESM_WACCM_ssp585_2015_2300_yearly([4],2015,2015,3600), quit"

    if nargin < 1 || isempty(steps)
        steps = [4];
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

    CMIP_MODEL = 'CESM2-WACCM';
    SCENARIO   = 'ssp585';
    proj_root  = './../../../';
    init_dir   = [proj_root 'init/'];

    addpath('./../../../init/scripts');
    addpath('./../../local_matlab');   % shared across all models/scenarios (like
                                        % ./../../functions); must be added AFTER
                                        % devpath (already run by the calling shell
                                        % before this function executes) so this
                                        % shadows the trunk's broken waitonlock.m --
                                        % only used by step 4.

    inputmodel_2015 = [proj_root 'hist_runs/CESM2-WACCM/Models/' ...
                       'AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat'];

    % Raw data paths needed by the build steps (1-3). Same paths as the
    % original script's own header, read fresh from source before writing
    % this.
    raw_ssp     = [proj_root 'raw_data/ISMIP7/AIS/' CMIP_MODEL '/' SCENARIO '/'];
    tf_dir      = [raw_ssp  'ocean/tf/v3/'];
    smb_ssp_dir = [raw_ssp  'SDBN1-2000m/acabf/v2/'];
    grad_ssp_dir= [raw_ssp  'SDBN1-2000m/dacabfdz/v2/'];
    collapse_nc_v2 = [raw_ssp 'fracture/v2/ice_shelf_collapse_mask_' ...
                      lower(strrep(CMIP_MODEL,'-','')) '_' SCENARIO '_ismip7_8km-v2.nc'];
    sec_to_year = 31556926;   % consistent with hist_run_tune_CESM_WACCM

    preproc_ocean      = [proj_root 'preprocessed_data/Ocean/'];
    preproc_hist_atmo  = [proj_root 'preprocessed_data/Atmosphere/Hist/'];
    preproc_clim       = [proj_root 'preprocessed_data/Atmosphere/Clim/' CMIP_MODEL '/'];
    preproc_proj_ocean = [proj_root 'preprocessed_data/Ocean/Proj/' CMIP_MODEL '/' SCENARIO '/'];
    preproc_proj_atmo  = [proj_root 'preprocessed_data/Atmosphere/Proj/' CMIP_MODEL '/' SCENARIO '/'];

    modeldir = './Models_yearly/';
    if ~exist(modeldir, 'dir'), mkdir(modeldir); end

    org = organizer('repository', modeldir, ...
                    'prefix', ['AIS_ISMIP7_Proj_' SCENARIO '_yearly_'], ...
                    'steps', steps, 'color', '34;47;2');
    clear steps;

    % Build steps (1-3) always cover the FULL 2015-2300 span, independent
    % of whichever start_year/end_year this invocation's step 4 loop will
    % actually use.
    build_start_year = 2015;
    build_end_year   = 2300;

    % ================================================================= Step 1
    if perform(org, 'ProjTF') % {{{
        % Copied verbatim from the original script's own ProjTF step (read
        % fresh from source before writing this).
        md = loadmodel(inputmodel_2015);

        if ~exist(preproc_proj_ocean, 'dir'), mkdir(preproc_proj_ocean); end

        tf_files = dir([tf_dir 'tf_AIS_*.nc']);
        tf_files = sort({tf_files.name});

        z_data  = double(ncread([tf_dir tf_files{1}], 'z'));
        nDepths = length(z_data);
        nVerts  = md.mesh.numberofvertices;

        years        = build_start_year : build_end_year;
        tf_mat       = zeros(nVerts, length(years), nDepths);
        last_real_yr = -Inf;

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
                if yr < build_start_year || yr > build_end_year, continue; end
                ki      = yr - build_start_year + 1;
                tf_data = double(ncread(fpath, 'tf', [1 1 1 ti], [Inf Inf Inf 1]));
                for i = 1:nDepths
                    v = InterpFromGridToMesh(x_n, y_n, tf_data(:,:,i)', ...
                                            md.mesh.x, md.mesh.y, 0);
                    tf_mat(:, ki, i) = max(v, 0);
                end
                last_real_yr = max(last_real_yr, yr);
            end
        end

        if isinf(last_real_yr)
            error('No real TF data found in [%d, %d].', build_start_year, build_end_year);
        end

        nkeep  = last_real_yr - build_start_year + 1;
        years  = build_start_year : last_real_yr;
        nyears = length(years);
        tf_mat = tf_mat(:, 1:nkeep, :);
        t_vec  = years;
        if last_real_yr < build_end_year
            fprintf(['[INFO] TF data only found through %d (< end_year=%d); ' ...
                     'forcing series ends there -- relying on ISSM to hold ' ...
                     'the last value for later times (interp_forcing=0).\n'], ...
                    last_real_yr, build_end_year);
        end

        tf_proj = cell(1, 1, nDepths);
        for i = 1:nDepths
            slice = squeeze(tf_mat(:,:,i));
            tf_proj{1,1,i} = [slice ; t_vec];
        end

        save([preproc_proj_ocean 'CESM_WACCM_TF_' SCENARIO '_' ...
              num2str(build_start_year) '_' num2str(build_end_year) '.mat'], ...
             'tf_proj', 'z_data', 't_vec', '-v7.3');
        fprintf('Saved TF: %d years (%d-%d), %d depths.\n', nyears, build_start_year, last_real_yr, nDepths);
    end % }}}

    % ================================================================= Step 2
    if perform(org, 'ProjSMB') % {{{
        % Copied verbatim from the original script's own ProjSMB step (read
        % fresh from source before writing this).
        md = loadmodel(inputmodel_2015);

        if ~exist(preproc_proj_atmo, 'dir'), mkdir(preproc_proj_atmo); end

        nVerts = md.mesh.numberofvertices;

        clim_yr0 = 1995; clim_yr1 = 2014;
        hist_smb_yr0 = 1995; hist_smb_yr1 = 2020;
        load([preproc_clim 'CESM_WACCM_SMB_clim_' ...
              num2str(clim_yr0) '_' num2str(clim_yr1) '.mat'], 'smb_racmo', 'cesm_mean');
        load([preproc_hist_atmo 'CESM_WACCM_SMB_corrected_' ...
              num2str(hist_smb_yr0) '_' num2str(hist_smb_yr1) '.mat'], 'p_vert');
        cesm_hist_mean = cesm_mean;

        nc_first = [smb_ssp_dir sprintf('acabf_AIS_%s_%s_SDBN1-2000m_v2_%d.nc', CMIP_MODEL, SCENARIO, build_start_year)];
        x_s = double(ncread(nc_first, 'x'));
        y_s = double(ncread(nc_first, 'y'));

        years        = build_start_year : build_end_year;
        nyears_max   = length(years);
        smb_matrix   = zeros(nVerts, nyears_max);
        bgrad_matrix = zeros(nVerts, nyears_max);
        last_real_yr = build_start_year - 1;

        for k = 1:nyears_max
            yr = years(k);

            nc_smb  = [smb_ssp_dir sprintf('acabf_AIS_%s_%s_SDBN1-2000m_v2_%d.nc', CMIP_MODEL, SCENARIO, yr)];
            nc_grad = [grad_ssp_dir sprintf('dacabfdz_AIS_%s_%s_SDBN1-2000m_v2_%d.nc', CMIP_MODEL, SCENARIO, yr)];

            if isempty(ncread(nc_smb, 'time')) || isempty(ncread(nc_grad, 'time'))
                fprintf('[INFO] %d: acabf/dacabfdz has no data (empty placeholder) -- stopping SMB series here.\n', yr);
                break;
            end

            am     = mean(double(ncread(nc_smb, 'acabf')), 3);
            cesm_yr = InterpFromGridToMesh(x_s, y_s, am', md.mesh.x, md.mesh.y, 0) * sec_to_year;
            smb_matrix(:,k) = p_vert .* smb_racmo + (cesm_yr - cesm_hist_mean);

            g_raw    = squeeze(double(ncread(nc_grad, 'dacabfdz')));
            if ndims(g_raw) == 3
                g_raw = mean(g_raw, 3);
            end
            bgrad_matrix(:,k) = InterpFromGridToMesh(x_s, y_s, g_raw', md.mesh.x, md.mesh.y, 0) * sec_to_year;

            last_real_yr = yr;

            if mod(yr, 10) == 0
                fprintf('  SMB+grad year %d\n', yr);
            end
        end

        if last_real_yr < build_start_year
            error('No real SMB data found in [%d, %d].', build_start_year, build_end_year);
        end

        nkeep        = last_real_yr - build_start_year + 1;
        years        = build_start_year : last_real_yr;
        nyears       = length(years);
        smb_matrix   = smb_matrix(:, 1:nkeep);
        bgrad_matrix = bgrad_matrix(:, 1:nkeep);
        t_smb        = years;

        smb_forcing   = [smb_matrix  ; t_smb];
        bgrad_forcing = [bgrad_matrix ; t_smb];

        save([preproc_proj_atmo 'CESM_WACCM_SMB_' SCENARIO '_' ...
              num2str(build_start_year) '_' num2str(build_end_year) '.mat'], ...
             'smb_forcing', 'bgrad_forcing', 't_smb', '-v7.3');
        fprintf('Saved SMB forcing: %d years (%d-%d).\n', nyears, build_start_year, last_real_yr);
    end % }}}

    % ================================================================= Step 3
    if perform(org, 'ProjLevelsetV2') % {{{
        % Same build_spclevelset() logic already written and verified in
        % build_and_compare_levelset_v2_ssp585.m, copied here (see that
        % helper function below) so this script needs nothing external.
        % Produces the SAME output file that script's own step 1 produces
        % -- interchangeable, not a second copy.
        proj_spclevelset_v2 = build_spclevelset(inputmodel_2015, collapse_nc_v2, build_start_year, build_end_year);
        v2_fname = [preproc_proj_ocean 'CESM_WACCM_levelset_v2_' SCENARIO '_' ...
                   num2str(build_start_year) '_' num2str(build_end_year) '.mat'];
        save(v2_fname, 'proj_spclevelset_v2', '-v7.3');
        fprintf('Saved: %s\n', v2_fname);
    end % }}}

    % ================================================================= Step 4
    if perform(org, 'ProjRunYearly') % {{{
        % The annual loop -- UNCHANGED from the previous version of this
        % script, except the levelset load now points at v2 (step 3's
        % output) instead of v1.
        fprintf('=== Step 4: annual loop, %d -> %d ===\n', start_year, end_year);

        % Same requested_outputs list as the original ssp585 script.
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

        % ---- Load the FULL pre-built forcing/basin/mask arrays ONCE.
        % These come from steps 1-3 above (run previously, or in this same
        % invocation if requested together) -- nothing here re-derives
        % anything from raw NetCDFs. ----
        load([preproc_ocean 'Basins/Imbie2_extrap_2km_BasinOnElements.mat']);  % -> basinid
        load([preproc_ocean 'tf_depths.mat']);                                 % -> tf_depths
        load([preproc_ocean 'gamma0_local.mat']);                              % -> gamma0_local
        tmp = load([preproc_ocean 'dT_correction.mat'], 'dT_correction');
        delta_t = tmp.dT_correction;
        unique_basinid = unique(basinid);

        load([preproc_proj_ocean 'CESM_WACCM_TF_' SCENARIO '_2015_2300.mat']);         % -> tf_proj, z_data
        load([preproc_proj_atmo  'CESM_WACCM_SMB_' SCENARIO '_2015_2300.mat']);        % -> smb_forcing, bgrad_forcing

        % v2 levelset (step 3's output) -- aliased to proj_spclevelset so
        % the rest of this loop's logic (unchanged from before) doesn't
        % need to reference a different variable name.
        loaded_lset = load([preproc_proj_ocean 'CESM_WACCM_levelset_v2_' SCENARIO '_2015_2300.mat']);   % -> proj_spclevelset_v2
        proj_spclevelset = loaded_lset.proj_spclevelset_v2;
        clear loaded_lset;

        % href must be the 1995 relaxed surface (same for every year --
        % smbref is calibrated against the 1995 RACMO climatology).
        md_relax  = loadmodel([init_dir 'Models/AIS_ISMIP7_Relaxed_CESM_WACCM.mat']);
        surf_1995 = md_relax.geometry.surface;
        clear md_relax

        % ---- Restart state ----
        % start_year==2015: fresh start from the historical run's 2015
        % end-state. applied_collapse is SEEDED (not empty) with every
        % vertex that already had NO ice at 2015 (md.mask.ice_levelset>0)
        % -- this permanently forces the calving front to never ADVANCE
        % beyond the 2015 extent for the rest of the run. Reason:
        % md.calving.calvingrate and md.frontalforcings.meltingrate are
        % both zero everywhere in this run (calving is fully prescribed
        % by the collapse mask, per design), so outside the collapse
        % mask's own footprint there is NOTHING opposing forward
        % ice-velocity advection at the front -- it can only ever push
        % forward, capped in RATE by migration_max but with no DISTANCE
        % limit, letting it advance indefinitely into open ocean over 285
        % years (observed as spurious advance in e.g. the Ross sector).
        % 2015 is also the one point in this whole run backed by a real
        % observation (the historical run's own Greene-driven end-state)
        % -- everything from 2015 onward is collapse-mask/free-physics
        % only, so holding that boundary as a permanent ceiling protects
        % the one piece of ground truth available, exactly like
        % applied_collapse already does for cells the collapse mask
        % itself forces (same mechanism, just seeded from a different
        % source).
        %
        % start_year>2015: RESUME from a previous run of this same
        % script's step 4 from that year (e.g. after a crash at year
        % 2100 -- rerun with start_year=2100). Loads two things saved by
        % a prior run for that year:
        %   1. the solved model AIS_ISMIP7_Proj_ssp585_yearly_<start_year>.mat
        %      (same file every iteration below already saves), restart-
        %      extracted exactly like the in-loop restart-prep block does.
        %   2. applied_collapse itself, from the small companion file
        %      AIS_ISMIP7_Proj_ssp585_yearly_appliedcollapse_<start_year>.mat.
        %      This MUST be restored from disk, not reinitialised to
        %      empty -- it is the permanent record of every vertex
        %      already gated-and-forced to ocean in earlier years. Losing
        %      it on resume would silently revert all of that history to
        %      NaN (unconstrained), risking exactly the
        %      spurious-reglaciation artifact the persistence was built
        %      to prevent (see conversation).
        if start_year == 2015
            md = loadmodel(inputmodel_2015);
            applied_collapse = md.mask.ice_levelset > 0;   % seed: no ice at 2015 -> never ice again
            fprintf('Seeded applied_collapse with %d vertices already ice-free at 2015 (of %d total).\n', ...
                    sum(applied_collapse), md.mesh.numberofvertices);
        else
            resume_model_fname = [modeldir 'AIS_ISMIP7_Proj_' SCENARIO '_yearly_' num2str(start_year) '.mat'];
            resume_ac_fname    = [modeldir 'AIS_ISMIP7_Proj_' SCENARIO '_yearly_appliedcollapse_' num2str(start_year) '.mat'];
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

        % end_year is INCLUSIVE -- the last calendar year you want
        % simulated, matching the ORIGINAL script's own end_year=2300
        % meaning exactly (there, final_time=end_year+1 produces the
        % calendar-year-end_year output snapshot -- see
        % WriteISMIP6_NetCDF_test's time_yr=t_annual-1 convention). The
        % last iteration below is therefore year=end_year itself
        % (start_time=end_year, final_time=end_year+1), NOT end_year-1
        % -- a full run needs end_year=2300 to reproduce the original's
        % true 286 simulated years (2015 through 2301 in ISSM-internal
        % time).
        for year = start_year:end_year
            fprintf('\n=== Year %d -> %d ===\n', year, year + 1);

            m = ((1 + sin(71*pi/180)) * ones(md.mesh.numberofvertices, 1) ...
                 ./ (1 + sin(abs(md.mesh.lat)*pi/180)));
            md.mesh.scale_factor = (1./m).^2;

            md.inversion.iscontrol       = 0;
            md.transient.isthermal       = 0;
            md.transient.isgroundingline = 1;
            md.transient.ismasstransport = 1;
            md.transient.isstressbalance = 1;
            md.masstransport.spcthickness   = NaN*ones(md.mesh.numberofvertices, 1);
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

            % --- Ocean: slice the full TF series to this 1-year window -----
            % ISSM's interp_forcing=0 (step-function) forcing looks BACKWARD
            % at each internal solver sub-step: it uses the value tagged at
            % the largest time <= the CURRENT simulated time, held constant
            % until the next tag is reached. For a solve spanning [year,
            % year+1), the value applicable for that WHOLE interval is
            % therefore the entry tagged at <= year (not year+1 -- a column
            % tagged year+1 would leave every sub-step strictly before
            % year+1 with no backward-applicable tag at all). Using year also
            % correctly reproduces "hold the last real value" once real data
            % has run out: original script's own ProjTF comment notes real TF
            % data only exists through 2299 (decade chunks 2015-2299, no 2300
            % file); for any iteration year >= 2299, find(...,<=year,'last')
            % still resolves to the 2299 entry, exactly matching how the
            % original's own (lower-bound-only trimmed, multi-year) segment-2
            % array naturally holds 2299's value for every later sub-step.
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

            % --- SMB: same "most recent at-or-before year" logic as TF ---
            % (per-year acabf/dacabfdz files can also run out before 2300 --
            % see the original script's ProjSMB comment on empty placeholder
            % files -- a range filter would go empty the same way TF's would).
            idx_smb    = find(smb_forcing(end, :) <= year, 1, 'last');
            smb_year   = smb_forcing(:, idx_smb);
            bgrad_year = bgrad_forcing(:, idx_smb);

            md.smb        = SMBgradients();
            md.smb.smbref = smb_year;
            md.smb.b_pos  = bgrad_year;
            md.smb.b_neg  = bgrad_year;
            md.smb.href   = [surf_1995 ; 1995];

            % --- Calving front: collapse-mask forcing GATED to floating ice
            % only, using the model's OWN current simulated state -- this is
            % the actual fix, not just a static 2015 snapshot.
            %
            % md.mask.ice_levelset / md.mask.ocean_levelset, at this point in
            % the loop, still hold whatever the PREVIOUS iteration's restart
            % block wrote (or the AIS_state_2015 initial condition, for the
            % very first year) -- i.e. the model's true, just-solved state as
            % of the END of the previous year, not yet touched by this year's
            % own solve. That is exactly the "current" floating extent to
            % gate against: ice_levelset<0 (has ice) & ocean_levelset<0
            % (floating, not grounded -- MaskOceanLevelset>=0 means grounded,
            % see write_ismip7_2d_projection.m's sftgrf construction).
            is_floating_now = (md.mask.ice_levelset < 0) & (md.mask.ocean_levelset < 0);

            % raw_spclevelset_year: the SAME columns the original (ungated)
            % logic would have used for this 1-year window, straight from
            % the v2 ProjLevelsetV2 output (step 3) -- 1 means "the raw v2
            % collapse mask says this vertex has collapsed by this point"
            % (already cumulative from that step's own monotonic build),
            % NaN means "no collapse-mask constraint proposed here."
            keep_lset = (proj_spclevelset(end, :) >= year) & (proj_spclevelset(end, :) <= year + 1);
            raw_spclevelset_year = proj_spclevelset(:, keep_lset);
            time_row_year = raw_spclevelset_year(end, :);
            raw_forced_this_year = any(raw_spclevelset_year(1:end-1, :) == 1, 2);

            % Gate: only newly ADD a vertex to the permanently-forced set if
            % the raw mask proposes it AND it is floating right now. A vertex
            % the raw mask proposes while still grounded is simply not forced
            % this year -- it gets re-evaluated next year against the model's
            % then-current state, so it can still be gated in later once the
            % model's own grounding line genuinely retreats past it.
            applied_collapse = applied_collapse | (raw_forced_this_year & is_floating_now);

            % gated_col has nVerts+1 rows (matching proj_spclevelset's own
            % [values; time] shape) so the trailing row is a real time row,
            % not an accidental overwrite of the last vertex's value.
            gated_col = NaN(md.mesh.numberofvertices + 1, 1);
            gated_col(applied_collapse) = 1;   % logical index (nVerts long) only ever touches the first nVerts rows
            md.levelset.spclevelset = repmat(gated_col, 1, size(raw_spclevelset_year, 2));
            md.levelset.spclevelset(end, :) = time_row_year;   % restore the real time tags

            md.calving.calvingrate         = zeros(md.mesh.numberofvertices, 1);
            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);
            md.transient.ismovingfront     = 1;
            md.levelset.migration_max      = 2863.78;   % m/yr

            md.miscellaneous.name = ['ProjRunYearly_' CMIP_MODEL '_' SCENARIO '_' ...
                                     num2str(year) '_' num2str(year + 1)];
            md.cluster             = set_cluster('gadi');
            md.settings.waitonlock = waitonlock_seconds;   % SECONDS -- local waitonlock.m override
            md.verbose = verbose('solution', true, 'module', true, 'convergence', true);

            try
                % loadonly=0 with waitonlock>0: submits, BLOCKS until the
                % local override detects the job's .lock+.outlog files, then
                % loads results -- all in this one call (matches the
                % colleague's runme_transient.m usage pattern).
                md = solve(md, 'tr', 'runtimename', false, 'loadonly', 0);
            catch ME
                fprintf('\n*** Year %d -> %d FAILED: %s\n', year, year + 1, ME.message);
                fprintf('*** Stopping loop. Investigate, then resume from year %d.\n', year);
                rethrow(ME);
            end

            % Save this year's solved model, keeping TransientSolution intact
            % (needed later to concatenate all years for the final ISMIP6
            % NetCDF write, same as the original script's 2-segment
            % concatenation, just longer).
            yearly_fname = [modeldir 'AIS_ISMIP7_Proj_' SCENARIO '_yearly_' num2str(year + 1) '.mat'];
            save(yearly_fname, 'md', '-v7.3');
            fprintf('Saved: %s\n', yearly_fname);

            % Save applied_collapse alongside it -- required to correctly
            % resume a later run with start_year=year+1 (see the restart-state
            % block above for why this must be restored, not reinitialised).
            ac_fname = [modeldir 'AIS_ISMIP7_Proj_' SCENARIO '_yearly_appliedcollapse_' num2str(year + 1) '.mat'];
            save(ac_fname, 'applied_collapse', '-v7.3');
            fprintf('Saved: %s\n', ac_fname);

            % Prepare restart state for the next iteration.
            md_solved = md;
            md.geometry.thickness   = md_solved.results.TransientSolution(end).Thickness;
            md.geometry.surface     = md_solved.results.TransientSolution(end).Surface;
            md.geometry.base        = md_solved.results.TransientSolution(end).Base;
            md.mask.ocean_levelset  = md_solved.results.TransientSolution(end).MaskOceanLevelset;
            md.mask.ice_levelset    = md_solved.results.TransientSolution(end).MaskIceLevelset;
            md.results.TransientSolution = [];
            clear md_solved;
        end

        fprintf('\n=== Step 4 complete: %d -> %d ===\n', start_year, end_year);
    end % }}}
end


function proj_spclevelset = build_spclevelset(inputmodel_2015, collapse_nc, start_year, end_year)
% Copied verbatim from build_and_compare_levelset_v2_ssp585.m's own
% build_spclevelset() (itself copied verbatim from the original script's
% ProjLevelset step) -- only collapse_nc is parameterised, so this same
% logic can be pointed at v1 or v2 data. See that step's own comments for
% the full rationale (cumulative collapse, advance-suppression).

    md = loadmodel(inputmodel_2015);
    nVerts = md.mesh.numberofvertices;

    x_c = double(ncread(collapse_nc, 'x'));
    y_c = double(ncread(collapse_nc, 'y'));
    time_c = double(ncread(collapse_nc, 'time'));  % integer years

    model_years      = start_year : end_year + 1;   % include final_time year
    spclevelset_mat  = NaN(nVerts, length(model_years));
    collapsed_so_far = false(nVerts, 1);  % cumulative collapse state

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
    no_ice0       = md.mask.ice_levelset > 0;
    already_ocean = repmat(no_ice0, 1, length(model_years)) & isnan(spclevelset_mat); %#ok<NASGU>
    % no action needed -- NaN at vertices that start as ocean is correct

    proj_spclevelset = [spclevelset_mat ; model_years];
end
