function [cfflux_tot, glflux_tot] = write_ismip7_2d_projection(md, outdir, meta)
% WRITE_ISMIP7_2D_PROJECTION  Grid and write all 2D ISMIP7 NetCDF files for a
%   projection run (e.g. 2015–2300, two concatenated segments).
%
%   md      : ISSM model object.  md.results.TransientSolution must be the
%             concatenation of all run segments (including safety steps).
%   outdir  : output directory (must already exist).
%   meta    : struct with fields experiment_id, set_counter, time_range,
%             ESM_id, forcing_member_id, ISM_member_id (+ optional defaults).
%
%   Returns cfflux_tot and glflux_tot (1×nT kg s-1 AIS-wide totals from
%   flux_along_contour_2d) so write_ismip7_scalar_projection can use them for
%   the IcefrontMassFluxLevelset / GroundinglineMassFlux sign check.
%
%   Requires on MATLAB path (via addpath in calling script):
%     init/scripts/write_ismip7_2d.m
%     init/scripts/flux_along_contour_2d.m

    rhoi     = md.materials.rho_ice;
    yts      = md.constants.yts;
    cellsize = 8000;   % ISMIP7 8 km grid

    % ---- filter TransientSolution -----------------------------------------
    % ISSM emits one sub-annual safety step at the start of each segment
    % (e.g. t≈2015.1, t≈2151.1).  Keep only the annual outputs.
    t_raw = [md.results.TransientSolution.time];
    keep  = abs(t_raw - round(t_raw)) < 0.05;
    md.results.TransientSolution = md.results.TransientSolution(keep);
    nT    = sum(keep);

    % ISSM convention: t=Y means end of year Y-1.
    t_annual = round(t_raw(keep));   % ISSM time [2016, ..., 2300]
    time_yr  = t_annual - 1;         % calendar year [2015, ..., 2299]

    % ---- time vectors (days since 1850-01-01, standard Gregorian calendar) --
    % ST: Jan 1 of following year = date(yr+1,1,1) per ISMIP7 spec.
    % FL: July 1 of labelled year (midpoint), bounds = [Jan 1 yr, Jan 1 yr+1].
    ref_dn  = datenum(1850, 1, 1);
    time_st = zeros(1, nT);
    lb_fl   = zeros(1, nT); ub_fl = zeros(1, nT); time_fl = zeros(1, nT);
    for i = 1:nT
        time_st(i) = datenum(t_annual(i), 1,  1) - ref_dn;
        lb_fl(i)   = datenum(time_yr(i),  1,  1) - ref_dn;
        ub_fl(i)   = datenum(time_yr(i)+1,1,  1) - ref_dn;
        time_fl(i) = datenum(time_yr(i),  7,  1) - ref_dn;
    end
    time_bnds_fl = [lb_fl; ub_fl]';   % nT × 2

    % ---- grid setup -------------------------------------------------------
    [~, xGrid, yGrid] = gridData(md, md.results.TransientSolution(1).Thickness, ...
        'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]);
    nx = length(xGrid); ny = length(yGrid);

    % ---- standard loop variables ------------------------------------------
    % cols: varname | ISSM field | scale | standard_name | long_name | units | ST/FL | clamp_min | clamp_max
    % clamp limits from ISMIP7_variable_request.csv (AIS column).
    % In pure SSA, xvelmean/yvelmean = xvelsurf/yvelsurf (no vertical shear).
    vars2d = { ...
        'lithk',       'Thickness',                           1,        'land_ice_thickness',                          'Ice thickness',                          'm',         'ST',  0,       5000;    ...
        'orog',        'Surface',                             1,        'surface_altitude',                            'Surface elevation',                      'm',         'ST',  0,       4500;    ...
        'base',        'Base',                                1,        '',                                            'Ice base elevation',                     'm',         'ST', -4000,   4000;    ...
        'topg',        'Bed',                                 1,        'bedrock_altitude',                            'Bedrock elevation',                      'm',         'ST', -7000,   4000;    ...
        'xvelsurf',    'Vx',                                  1/yts,    'land_ice_surface_x_velocity',                 'Surface velocity in x',                  'm s-1',     'ST', -4e-4,   4e-4;   ...
        'yvelsurf',    'Vy',                                  1/yts,    'land_ice_surface_y_velocity',                 'Surface velocity in y',                  'm s-1',     'ST', -4e-4,   4e-4;   ...
        'xvelmean',    'Vx',                                  1/yts,    'land_ice_vertical_mean_x_velocity',           'Vertically averaged velocity in x',       'm s-1',     'ST', -4e-4,   4e-4;   ...
        'yvelmean',    'Vy',                                  1/yts,    'land_ice_vertical_mean_y_velocity',           'Vertically averaged velocity in y',       'm s-1',     'ST', -4e-4,   4e-4;   ...
        'acabf',       'SmbMassBalance',                      rhoi/yts, 'land_ice_surface_specific_mass_balance_flux', 'Surface mass balance flux',               'kg m-2 s-1','FL', -6e-4,   1e-3;   ...
        'libmassbfgr', 'BasalforcingsGroundediceMeltingRate', -rhoi/yts,'land_ice_basal_specific_mass_balance_flux',  'Basal mass balance flux beneath grounded ice', 'kg m-2 s-1','FL', -3e-4, 1e-4; ...
        'libmassbffl', 'BasalforcingsFloatingiceMeltingRate', -rhoi/yts,'land_ice_basal_specific_mass_balance_flux',  'Basal mass balance flux beneath floating ice',  'kg m-2 s-1','FL', -8e-3, 1e-3; ...
    };
    for k = 1:size(vars2d, 1)
        varname       = vars2d{k,1};
        issmfield     = vars2d{k,2};
        scale         = vars2d{k,3};
        standard_name = vars2d{k,4};
        long_name     = vars2d{k,5};
        units         = vars2d{k,6};
        vtype         = vars2d{k,7};
        clamp_min     = vars2d{k,8};
        clamp_max     = vars2d{k,9};
        if strcmp(vtype, 'ST'), t_vec = time_st; t_bnds = [];
        else,                   t_vec = time_fl; t_bnds = time_bnds_fl; end
        data3d = zeros(ny, nx, nT, 'single');
        for t = 1:nT
            data3d(:,:,t) = single(gridData(md, ...
                md.results.TransientSolution(t).(issmfield), ...
                'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]) * scale);
        end
        safe_max = single(clamp_max);
        if double(safe_max) > clamp_max, safe_max = safe_max - eps(safe_max); end
        safe_min = single(clamp_min);
        if double(safe_min) < clamp_min, safe_min = safe_min + eps(safe_min); end
        data3d = max(safe_min, min(safe_max, data3d));
        write_ismip7_2d(outdir, varname, data3d, xGrid, yGrid, ...
            t_vec, t_bnds, standard_name, long_name, units, meta);
    end

    % ---- static fields ----------------------------------------------------
    % hfgeoubed: clamp [0, 0.3] W m-2.  single(0.3) rounds to 0.30000001192...
    % in float64, so subtract one ULP to stay strictly ≤ 0.3.
    hfgeoubed2d = single(gridData(md, md.basalforcings.geothermalflux, ...
        'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]));
    hfgeoubed2d = max(single(0), min(single(0.3) - eps(single(0.3)), hfgeoubed2d));
    write_ismip7_2d(outdir, 'hfgeoubed', repmat(hfgeoubed2d, 1, 1, nT), xGrid, yGrid, ...
        time_fl, time_bnds_fl, 'upward_geothermal_heat_flux_in_land_ice', ...
        'Geothermal heat flux', 'W m-2', meta);

    % litemptop: clamp [183, 290] K
    temp2d = single(gridData(md, md.initialization.temperature, ...
        'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]));
    temp2d = max(single(183), min(single(290), temp2d));
    write_ismip7_2d(outdir, 'litemptop', repmat(temp2d, 1, 1, nT), xGrid, yGrid, ...
        time_st, [], 'temperature_at_top_of_ice_sheet_model', ...
        'Surface temperature', 'K', meta);

    % ---- dlithkdt (FL: annual-mean thickness tendency) --------------------
    % clamp [-1e-4, 1e-4] m s-1
    H_grid = zeros(ny, nx, nT);
    for t = 1:nT
        H_grid(:,:,t) = gridData(md, md.results.TransientSolution(t).Thickness, ...
            'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]);
    end
    dlithkdt3d = zeros(ny, nx, nT, 'single');
    for t = 1:nT
        if t < nT
            dlithkdt3d(:,:,t) = single((H_grid(:,:,t+1) - H_grid(:,:,t)) / yts);
        else
            dlithkdt3d(:,:,t) = single((H_grid(:,:,t) - H_grid(:,:,t-1)) / yts);
        end
    end
    dlithkdt3d = max(single(-1e-4), min(single(1e-4), dlithkdt3d));
    write_ismip7_2d(outdir, 'dlithkdt', dlithkdt3d, xGrid, yGrid, ...
        time_fl, time_bnds_fl, 'tendency_of_land_ice_thickness', ...
        'Ice thickness tendency', 'm s-1', meta);

    % ---- ice-area fraction fields (ST) ------------------------------------
    sftgif3d = zeros(ny, nx, nT, 'single');
    sftgrf3d = zeros(ny, nx, nT, 'single');
    sftflf3d = zeros(ny, nx, nT, 'single');
    for t = 1:nT
        sol  = md.results.TransientSolution(t);
        ice  = gridData(md, double(sol.MaskIceLevelset  < 0), ...
            'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]);
        grnd = gridData(md, double(sol.MaskOceanLevelset >= 0), ...
            'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]);
        sftgif3d(:,:,t) = single(ice);
        sftgrf3d(:,:,t) = single(ice .* grnd);
        sftflf3d(:,:,t) = single(ice .* (1 - grnd));
    end
    sftgif3d = max(single(0), min(single(1), sftgif3d));
    sftgrf3d = max(single(0), min(single(1), sftgrf3d));
    sftflf3d = max(single(0), min(single(1), sftflf3d));
    write_ismip7_2d(outdir, 'sftgif', sftgif3d, xGrid, yGrid, ...
        time_st, [], 'land_ice_area_fraction', 'Land ice area fraction', '1', meta);
    write_ismip7_2d(outdir, 'sftgrf', sftgrf3d, xGrid, yGrid, ...
        time_st, [], 'grounded_ice_sheet_area_fraction', 'Grounded ice sheet area fraction', '1', meta);
    write_ismip7_2d(outdir, 'sftflf', sftflf3d, xGrid, yGrid, ...
        time_st, [], 'floating_ice_shelf_area_fraction', 'Floating ice shelf area fraction', '1', meta);

    % ---- strbasemag (ST): Budd p=1, q=1 → r=1, s=1 -----------------------
    % τ_b = C² × Neff × |u_b|  (Pa)
    % Neff = max(ρ_i·g·H + ρ_w·g·min(bed, 0), 0)  [hydrostatic, no subglacial hydrology]
    % Zero under floating ice; uses static md.friction.coefficient.
    rho_w    = md.materials.rho_water;
    C_fric   = md.friction.coefficient;     % nodal, static
    strbasemag3d = zeros(ny, nx, nT, 'single');
    for t = 1:nT
        sol   = md.results.TransientSolution(t);
        H_n   = sol.Thickness;
        bed_n = sol.Bed;
        vx_n  = sol.Vx;
        vy_n  = sol.Vy;
        mask_n = sol.MaskOceanLevelset;
        Neff_n = max(rhoi * md.constants.g * H_n + rho_w * md.constants.g * min(bed_n, 0), 0);
        ub_n   = sqrt(vx_n.^2 + vy_n.^2) / yts;   % m s-1
        str_n  = C_fric.^2 .* Neff_n .* ub_n;
        str_n(mask_n < 0) = 0;                         % zero under floating ice
        str_n  = max(0, min(3e5, str_n));             % clamp [0, 300 kPa]
        strbasemag3d(:,:,t) = single(gridData(md, str_n, ...
            'xRange', [-3040000, 3040000], 'yRange', [-3040000, 3040000]));
    end
    write_ismip7_2d(outdir, 'strbasemag', strbasemag3d, xGrid, yGrid, ...
        time_st, [], 'land_ice_basal_drag', 'Magnitude of basal shear stress', 'Pa', meta);

    % ---- licalvf / ligroundf (FL, contour-integral maps) ------------------
    glflux2d   = zeros(ny, nx, nT, 'single'); glflux_tot = zeros(1, nT);
    cfflux2d   = zeros(ny, nx, nT, 'single'); cfflux_tot = zeros(1, nT);
    for t = 1:nT
        sol = md.results.TransientSolution(t);
        [gl2d, gltot] = flux_along_contour_2d(md, sol, sol.MaskOceanLevelset, true,  xGrid, yGrid, cellsize);
        [cf2d, cftot] = flux_along_contour_2d(md, sol, sol.MaskIceLevelset,   false, xGrid, yGrid, cellsize);
        glflux2d(:,:,t) = gl2d; glflux_tot(t) = gltot;
        cfflux2d(:,:,t) = cf2d; cfflux_tot(t) = cftot;
    end
    write_ismip7_2d(outdir, 'ligroundf', max(single(0), min(single(1e11), glflux2d)), xGrid, yGrid, ...
        time_fl, time_bnds_fl, 'land_ice_specific_grounding_line_flux', ...
        'Grounding line flux', 'kg m-2 s-1', meta);
    write_ismip7_2d(outdir, 'licalvf', max(single(-1e11), min(single(0), -cfflux2d)), xGrid, yGrid, ...
        time_fl, time_bnds_fl, 'land_ice_specific_mass_flux_due_to_calving', ...
        'Calving flux', 'kg m-2 s-1', meta);
    write_ismip7_2d(outdir, 'lifmassbf', max(single(-1e11), min(single(0), -cfflux2d)), xGrid, yGrid, ...
        time_fl, time_bnds_fl, 'land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting', ...
        'Frontal mass flux (calving + melt; no frontal melt in this run)', 'kg m-2 s-1', meta);
end
