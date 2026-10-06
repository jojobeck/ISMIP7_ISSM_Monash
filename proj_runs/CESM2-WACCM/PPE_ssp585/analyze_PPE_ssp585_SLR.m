function analyze_PPE_ssp585_SLR(steps)
% ANALYZE_PPE_SSP585_SLR  Sea-level-rise (SLE) comparison across every
%   CESM2-WACCM ssp585 PPE (Perturbed-Parameter-Ensemble) member, plus the
%   core (non-PPE) ssp585 run and the CESM2-WACCM historical run, all on one
%   figure -- scoped to THIS folder only (analysis/analyze_SLR_all_experiments.m
%   stays the cross-model/cross-scenario comparison; this script is the
%   PPE-family-specific one). ssp585 counterpart of
%   ../PPE_ssp370/analyze_PPE_ssp370_SLR.m -- same approach, scenario/span
%   swapped (2015-2300 here, not 2015-2100).
%
% Gracefully includes every PPE yearly member's progress SO FAR -- a member
% whose loop has only reached e.g. 2035 out of its eventual 2300 still gets
% plotted up to 2035, it is NOT skipped just because it is unfinished (see
% load_ppe_yearly_vaf_series below, which simply globs whatever yearly .mat
% files already exist). A member with NO output files yet, or the core run
% / historical run not having been run, is skipped with a printed message
% rather than erroring out (same graceful-partial-rerun convention as
% analysis/analyze_SLR_all_experiments.m).
%
% FROZEN-FRONT MEMBERS (P011-P013) ARE SPLIT ACROSS TWO MODEL FILES, not
% one: proj_run_PPE_frozenfront_ssp585.m runs 2015-2300 as two restart
% segments (2015-2150, then 2151-2300 -- a single continuous 285-year
% solve got OOM-killed, same fix as this project's other 300-year-class
% runs). This script loads BOTH segment files and concatenates their VAF
% series; if only the first segment has been gathered so far, that
% member is plotted up to 2150 rather than being skipped entirely.
%
% SLE [m] = -(VAF(t) - VAF_ref) * rho_ice / rho_sw / A_ocean. VAF_ref is the
% CESM2-WACCM historical run's own FIRST annual-snapshot VAF (nominal 1995)
% -- the core run and every PPE member are all seeded from that same
% historical run's end-state (hist_runs/CESM2-WACCM/Models/
% AIS_ISMIP7_Hist1995_2014_AIS_state_2015.mat), so sharing one reference
% makes the historical curve and every branch join into one continuous
% trajectory from 1995, instead of each branch resetting to zero at 2015.
%
% Each PPE member's legend label is "P<NNN>: <difference_from_main>",
% reading <difference_from_main> straight out of ppe_members.csv (this
% folder) -- so the one-line description of what that member actually
% varies (e.g. "K = 95th percentile; ice front frozen at 2015 extent",
% "+ eigencalving (coeff = 6.592e15 m*s)", "SMB lapserate off") always
% matches the tracking table, with nothing duplicated/hardcoded here that
% could drift out of sync with it.
%
% Line style marks the MECHANISM, since there's only one model/scenario
% here (unlike the cross-experiment script, which uses line style for
% model): solid = frozen-front (ice front held fixed at 2015 extent),
% dashed = every other PPE member (moving calving front, collapse-mask-
% gated yearly loop) -- run_type in ppe_members.csv. Historical + the
% core ssp585 run are both drawn in black solid (continuous trunk); every
% PPE member gets its own colour (Okabe-Ito-extended qualitative palette,
% assigned by that member's row order in ppe_members.csv so it stays
% stable under future renumbering).
%
% Step map:
%   1  ComputeSLR   loads historical + core ssp585 + every ppe_members.csv
%                    row's model output, computes its SLE curve, saves to
%                    postprocessed_data/SLR/SLR_PPE_ssp585_CESM2-WACCM.mat
%   2  PlotSLR       reads that .mat back and plots every curve. Saves
%                    THREE figures to postprocessed_data/figures/SLR/:
%                      PPE_ssp585_CESM2-WACCM.png                 (full span)
%                      PPE_ssp585_CESM2-WACCM_zoom_2015_2040.png  (early zoom)
%                      PPE_ssp585_CESM2-WACCM_zoom_2280_2300.png  (end-of-run zoom)
%
% Usage (run from anywhere; paths are resolved relative to this script's
% own location, not the current directory):
%   analyze_PPE_ssp585_SLR([1])     % compute + save only
%   analyze_PPE_ssp585_SLR([2])     % plot only (needs step 1 already run)
%   analyze_PPE_ssp585_SLR([1 2])   % both (default)

    if nargin < 1 || isempty(steps)
        steps = [1 2];
    end

    here = fileparts(mfilename('fullpath'));   % this folder: proj_runs/CESM2-WACCM/PPE_ssp585/
    root = [here '/../../../'];                % -> repo root

    org = organizer('repository', [root 'postprocessed_data/SLR/'], ...
                    'prefix', 'SLR_PPE_ssp585_', 'steps', steps, 'color', '34;47;2');
    clear steps;

    matfile_out = [root 'postprocessed_data/SLR/SLR_PPE_ssp585_CESM2-WACCM.mat'];
    csv_file    = [here '/ppe_members.csv'];

    % ================================================================= Step 1
    if perform(org, 'ComputeSLR') % {{{
        addpath([root 'proj_runs/functions']);

        rho_ice = 917;       % kg/m^3
        rho_sw  = 1028;      % kg/m^3
        A_ocean = 3.625e14;  % m^2 (ISMIP6 standard ocean area)

        SLR = struct('p_number', {}, 'label', {}, 'run_type', {}, 'time', {}, 'sle', {});

        % ------------------------------------------------------------ Historical
        hist_file = [root 'hist_runs/CESM2-WACCM/Models/AIS_ISMIP7_Hist1995_2014_HistRun_1995_2014.mat'];
        fprintf('Loading CESM2-WACCM historical ...\n');
        vaf_ref = [];
        try
            [time, vaf] = load_vaf_from_matfile(hist_file);
            vaf_ref = vaf(1);
            sle = -(vaf - vaf(1)) * rho_ice / rho_sw / A_ocean;
            SLR(end+1) = struct('p_number', '', 'label', 'CESM2-WACCM historical', ...
                                 'run_type', 'historical', 'time', time, 'sle', sle); %#ok<AGROW>
            fprintf('  %d years (%d-%d), final SLE = %.4f m\n', numel(time), time(1), time(end), sle(end));
        catch ME
            fprintf('  SKIP historical -- not available (%s): %s\n', hist_file, ME.message);
        end

        % ------------------------------------------------------------ Core ssp585
        if ~isempty(vaf_ref)
            fprintf('Loading core ssp585 ...\n');
            core_modeldir = [root 'proj_runs/CESM2-WACCM/ssp585/Models_yearly/'];
            try
                [t_raw, vaf] = load_yearly_vaf_series(core_modeldir, 'ssp585');
                time = round(t_raw) - 1;
                sle  = -(vaf - vaf_ref) * rho_ice / rho_sw / A_ocean;
                SLR(end+1) = struct('p_number', '', 'label', 'Core ssp585 (baseline)', ...
                                     'run_type', 'core', 'time', time, 'sle', sle); %#ok<AGROW>
                fprintf('  %d years (%d-%d), final SLE = %.4f m\n', numel(time), time(1), time(end), sle(end));
            catch ME
                fprintf('  SKIP core ssp585 -- not available yet (%s): %s\n', core_modeldir, ME.message);
            end
        else
            fprintf('SKIP core ssp585 -- no historical VAF_ref available.\n');
        end

        % ------------------------------------------------------------ PPE members
        if ~isempty(vaf_ref)
            ppe_tbl = readtable(csv_file);
            for i = 1:size(ppe_tbl, 1)
                % char(string(...)) tolerates readtable returning either a
                % cell-of-char or a string-array column for text fields
                % (differs by MATLAB version/import-options defaults).
                p_number = char(strtrim(string(ppe_tbl.p_number(i))));
                run_type = char(strtrim(string(ppe_tbl.run_type(i))));
                diff_str = char(strtrim(string(ppe_tbl.difference_from_main(i))));
                label    = sprintf('%s: %s', p_number, diff_str);

                fprintf('Loading %s (%s) ...\n', p_number, run_type);
                try
                    switch run_type
                        case 'frozenfront'
                            % Restart-split at 2150/2151 (see
                            % proj_run_PPE_frozenfront_ssp585.m's own
                            % module docstring for why) -- TWO gathered
                            % model files per member, not one. Segment 1
                            % is required (its own try/catch -- the outer
                            % loop's -- skips the member entirely if it's
                            % missing); segment 2 is appended if present,
                            % so a member whose second segment hasn't
                            % finished/gathered yet still plots its first
                            % segment's progress rather than being skipped.
                            matfile1 = [here '/Models_' p_number '/AIS_ISMIP7_PPE_' p_number '_ProjRun_2015_2150.mat'];
                            matfile2 = [here '/Models_' p_number '/AIS_ISMIP7_PPE_' p_number '_ProjRun_2151_2300.mat'];
                            [time, vaf] = load_vaf_from_matfile(matfile1);
                            try
                                [time2, vaf2] = load_vaf_from_matfile(matfile2);
                                time = [time, time2];
                                vaf  = [vaf,  vaf2];
                            catch
                                % segment 2 not finished/gathered yet -- plot segment 1 alone
                            end
                        case 'yearly'
                            modeldir = [here '/Models_' p_number '/'];
                            [t_raw, vaf] = load_ppe_yearly_vaf_series(modeldir, p_number);
                            time = round(t_raw) - 1;
                        otherwise
                            error('Unrecognised run_type ''%s'' in %s row %d.', run_type, csv_file, i);
                    end
                    sle = -(vaf - vaf_ref) * rho_ice / rho_sw / A_ocean;
                    SLR(end+1) = struct('p_number', p_number, 'label', label, ...
                                         'run_type', run_type, 'time', time, 'sle', sle); %#ok<AGROW>
                    fprintf('  %d years (%d-%d), final SLE = %.4f m\n', numel(time), time(1), time(end), sle(end));
                catch ME
                    fprintf('  SKIP %s -- not available yet: %s\n', p_number, ME.message);
                end
            end
        else
            fprintf('SKIP all PPE members -- no historical VAF_ref available.\n');
        end

        outdir = [root 'postprocessed_data/SLR/'];
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        save(matfile_out, 'SLR', '-v7.3');
        fprintf('\nSaved %d curves to %s\n', numel(SLR), matfile_out);
    end % }}}

    % ================================================================= Step 2
    if perform(org, 'PlotSLR') % {{{
        loaded = load(matfile_out, 'SLR');
        SLR = loaded.SLR;

        % Okabe-Ito-extended qualitative palette, one colour per PPE member,
        % assigned by that member's position among the PPE (non-historical,
        % non-core) entries -- stable under future CSV renumbering.
        palette = { [0.9020 0.6235 0],    [0.3373 0.7059 0.9137], [0 0.6196 0.4510], ...
                    [0.8000 0.7000 0.00], [0 0.4470 0.7410],      [0.8353 0.3686 0], ...
                    [0.8 0.4745 0.6549],  [0.4660 0.6740 0.1880], [0.4940 0.1840 0.5560] };

        figdir = [root 'postprocessed_data/figures/SLR/'];
        if ~exist(figdir, 'dir'), mkdir(figdir); end

        idx_hist = find(strcmp({SLR.run_type}, 'historical'), 1);
        idx_core = find(strcmp({SLR.run_type}, 'core'), 1);
        idx_ppe  = find(strcmp({SLR.run_type}, 'frozenfront') | strcmp({SLR.run_type}, 'yearly'));

        % Build every line's (t,y,style) ONCE -- shared across all three
        % panels below, so the zoomed panels can compute their own tight
        % y-limits from exactly the same data they plot.
        lines = struct('t', {}, 'y', {}, 'color', {}, 'ls', {}, 'lw', {}, 'label', {});
        if ~isempty(idx_hist)
            h = SLR(idx_hist);
            lines(end+1) = struct('t', h.time, 'y', h.sle, ...
                'color', [0 0 0], 'ls', '-', 'lw', 2, 'label', h.label); %#ok<AGROW>
        end
        if ~isempty(idx_core) && ~isempty(idx_hist)
            c = SLR(idx_core);
            lines(end+1) = struct('t', [h.time(end), c.time], 'y', [h.sle(end), c.sle], ...
                'color', [0 0 0], 'ls', '-', 'lw', 2.5, 'label', c.label); %#ok<AGROW>
        end
        for k = 1:numel(idx_ppe)
            m = SLR(idx_ppe(k));
            col = palette{mod(k-1, numel(palette)) + 1};
            % Solid = frozen ice front (ismovingfront=0); dashed = every
            % other (moving-front) PPE member -- see module docstring.
            if strcmp(m.run_type, 'frozenfront'), ls = '-'; else, ls = '--'; end
            if ~isempty(idx_hist)
                t = [h.time(end), m.time]; y = [h.sle(end), m.sle];
            else
                t = m.time; y = m.sle;
            end
            lines(end+1) = struct('t', t, 'y', y, 'color', col, 'ls', ls, 'lw', 1.5, 'label', m.label); %#ok<AGROW>
        end

        % Three views of the same data: full span, early zoom (branch
        % point detail), end-of-run zoom (final-state spread between
        % members, once they've run far enough to show one).
        fig_specs = struct( ...
            'suffix', {'',               'zoom_2015_2040', 'zoom_2280_2300'}, ...
            'title',  {'sea-level contribution', '2015-2040 (zoom)', '2280-2300 (zoom)'}, ...
            'xlim',   {[], [2015 2040], [2280 2300]}, ...
            'vline2015', {false, true, false});

        for which_fig = 1:numel(fig_specs)
            spec = fig_specs(which_fig);
            fig = figure('visible', 'off', 'Position', [0 0 1200 600]);
            hold on;

            for li = 1:numel(lines)
                L = lines(li);
                plot(L.t, L.y, 'Color', L.color, 'LineStyle', L.ls, 'LineWidth', L.lw, ...
                     'DisplayName', L.label);
            end

            xlabel('Year');
            ylabel('Sea level contribution (m SLE)');
            title(['CESM2-WACCM ssp585 PPE -- ' spec.title]);
            if ~isempty(spec.xlim)
                xlim(spec.xlim);

                % Tight y-limits from ONLY the data actually visible in
                % this x-window, not each line's global min/max -- MATLAB's
                % default autoscale otherwise stays pinned to the full
                % 1995/2015-2300 range of every curve (e.g. the historical
                % run's much larger excursion), which is exactly why the
                % non-frozen-front members would otherwise look squashed
                % together near 2300: the y-axis was never actually
                % zoomed, only x.
                y_in_window = [];
                for li = 1:numel(lines)
                    L = lines(li);
                    in_win = L.t >= spec.xlim(1) & L.t <= spec.xlim(2);
                    y_in_window = [y_in_window, L.y(in_win)]; %#ok<AGROW>
                end
                if ~isempty(y_in_window)
                    y_lo = min(y_in_window); y_hi = max(y_in_window);
                    pad  = max(0.05 * (y_hi - y_lo), 1e-6);
                    ylim([y_lo - pad, y_hi + pad]);
                end

                if spec.vline2015
                    yl = ylim;
                    plot([2015 2015], yl, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
                    ylim(yl);
                end
            end
            legend('Location', 'eastoutside', 'FontSize', 8, 'Interpreter', 'none');
            grid on;
            box on;

            if isempty(spec.suffix)
                figname = fullfile(figdir, 'PPE_ssp585_CESM2-WACCM.png');
            else
                figname = fullfile(figdir, ['PPE_ssp585_CESM2-WACCM_' spec.suffix '.png']);
            end
            saveas(fig, figname);
            close(fig);
            fprintf('Saved: %s\n', figname);
        end
    end % }}}

end

% ---------------------------------------------------------------- helpers
function [time, vaf] = load_vaf_from_matfile(matfile)
% Load a single (non-yearly-restart) model .mat file's TransientSolution
% and return its annual-snapshot nominal years + IceVolumeAboveFloatationScaled,
% dropping the sub-annual "safety step" -- same convention as
% analysis/analyze_SLR_all_experiments.m's own local helper of the same name.
    loaded = load(matfile, 'md');
    ts     = loaded.md.results.TransientSolution;
    t_raw  = [ts.time];
    keep   = abs(t_raw - round(t_raw)) < 0.05;
    time   = round(t_raw(keep)) - 1;
    vaf    = [ts.IceVolumeAboveFloatationScaled];
    vaf    = vaf(keep);
end

function [time, vaf] = load_ppe_yearly_vaf_series(modeldir, p_number)
% Scalar-only loader for a PPE yearly-restart member's VAF, mirroring
% proj_runs/functions/load_yearly_vaf_series.m but matching THIS folder's
% own filename convention (AIS_ISMIP7_PPE_<p_number>_yearly_<year>.mat,
% NOT the core run's AIS_ISMIP7_Proj_<SCENARIO>_yearly_<year>.mat -- the
% prefix differs, so that shared helper can't be reused directly here).
%
% Deliberately tolerant of an UNFINISHED run: simply globs and loads
% whatever yearly .mat files already exist, in year order -- a member
% stopped partway through still returns its progress so far, with no gap
% check and no error (unlike discover_yearly_files.m's stricter
% contiguity warning, not needed for this plot-only use).
    pattern = [modeldir 'AIS_ISMIP7_PPE_' p_number '_yearly_*.mat'];
    files = dir(pattern);
    if isempty(files)
        error('No yearly model files found matching %s', pattern);
    end

    years = [];
    names = {};
    for i = 1:numel(files)
        fname = files(i).name;
        if contains(fname, 'appliedcollapse')
            continue;
        end
        tok = regexp(fname, '_yearly_(\d+)\.mat$', 'tokens', 'once');
        if isempty(tok)
            continue;
        end
        years(end+1) = str2double(tok{1}); %#ok<AGROW>
        names{end+1} = fname; %#ok<AGROW>
    end
    if isempty(years)
        error('No yearly model files (excluding appliedcollapse_*) found matching %s', pattern);
    end
    [~, order] = sort(years);
    names = names(order);

    time = [];
    vaf  = [];
    for i = 1:numel(names)
        loaded = load([modeldir names{i}], 'md');
        ts      = loaded.md.results.TransientSolution;
        t_raw   = [ts.time];
        keep    = abs(t_raw - round(t_raw)) < 0.05;
        vaf_raw = [ts.IceVolumeAboveFloatationScaled];
        time = [time, t_raw(keep)]; %#ok<AGROW>
        vaf  = [vaf,  vaf_raw(keep)]; %#ok<AGROW>
        clear loaded ts t_raw keep vaf_raw;
    end
end
