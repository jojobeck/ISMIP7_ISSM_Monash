function write_ismip7_scalar_projection(md, outdir, meta, cfflux_tot, glflux_tot)
% WRITE_ISMIP7_SCALAR_PROJECTION  Write all ISMIP7 scalar (t) NetCDF files for
%   a projection run (e.g. 2015–2300, two concatenated segments).
%
%   md          : ISSM model object.  md.results.TransientSolution must be the
%                 concatenation of all run segments (including safety steps).
%   outdir      : output directory (must already exist).
%   meta        : struct with fields experiment_id, set_counter, time_range,
%                 ESM_id, forcing_member_id, ISM_member_id (+ optional defaults).
%   cfflux_tot  : 1×nT calving-front flux totals (kg s-1) from
%                 write_ismip7_2d_projection — used to cross-check ISSM's
%                 IcefrontMassFluxLevelset sign convention.
%   glflux_tot  : 1×nT grounding-line flux totals (kg s-1), same source.
%
%   Requires on MATLAB path (via addpath in calling script):
%     init/scripts/write_ismip7_scalar.m

    yts  = md.constants.yts;
    rhoi = md.materials.rho_ice;

    % ---- filter TransientSolution (same logic as write_ismip7_2d_projection) -
    t_raw = [md.results.TransientSolution.time];
    keep  = abs(t_raw - round(t_raw)) < 0.05;
    md.results.TransientSolution = md.results.TransientSolution(keep);
    nT    = sum(keep);

    t_annual = round(t_raw(keep));   % ISSM time [2016, ..., 2300]
    time_yr  = t_annual - 1;         % calendar year [2015, ..., 2299]

    % ---- time vectors (days since 1850-01-01, standard Gregorian calendar) --
    ref_dn  = datenum(1850, 1, 1);
    time_st = zeros(1, nT);
    lb_fl   = zeros(1, nT); ub_fl = zeros(1, nT); time_fl = zeros(1, nT);
    for i = 1:nT
        time_st(i) = datenum(t_annual(i), 1, 1) - ref_dn;
        lb_fl(i)   = datenum(time_yr(i),  1, 1) - ref_dn;
        ub_fl(i)   = datenum(time_yr(i)+1,1, 1) - ref_dn;
        time_fl(i) = datenum(time_yr(i),  7, 1) - ref_dn;
    end
    time_bnds_fl = [lb_fl; ub_fl]';

    % ---- standard scalar variables ----------------------------------------
    scalars = { ...
        'lim',      'IceVolumeScaled',                rhoi,     'land_ice_mass',                                         'Total ice mass',         'kg',    'ST'; ...
        'limnsw',   'IceVolumeAboveFloatationScaled', rhoi,     'land_ice_mass_not_displacing_sea_water',                'Mass above floatation',  'kg',    'ST'; ...
        'iareagr',  'GroundedAreaScaled',              1,        'grounded_ice_sheet_area',                               'Grounded ice area',      'm^2',   'ST'; ...
        'iareafl',  'FloatingAreaScaled',              1,        'floating_ice_shelf_area',                               'Floating ice area',      'm^2',   'ST'; ...
        'tendacabf','TotalSmbScaled',                  1e12/yts, 'tendency_of_land_ice_mass_due_to_surface_mass_balance', 'Total SMB flux',         'kg s-1','FL'; ...
    };
    for k = 1:size(scalars, 1)
        varname       = scalars{k,1};
        issmfield     = scalars{k,2};
        scale         = scalars{k,3};
        standard_name = scalars{k,4};
        long_name     = scalars{k,5};
        units         = scalars{k,6};
        vtype         = scalars{k,7};
        if strcmp(vtype, 'ST'), t_vec = time_st; t_bnds = [];
        else,                   t_vec = time_fl; t_bnds = time_bnds_fl; end
        data1d = zeros(1, nT);
        for t = 1:nT
            data1d(t) = md.results.TransientSolution(t).(issmfield) * scale;
        end
        write_ismip7_scalar(outdir, varname, data1d, ...
            t_vec, t_bnds, standard_name, long_name, units, meta);
    end

    % ---- tendlibmassbfgr / tendlibmassbffl --------------------------------
    tgbmb = zeros(1, nT); tfbmb = zeros(1, nT);
    for t = 1:nT
        tgbmb(t) = md.results.TransientSolution(t).TotalGroundedBmbScaled * 1e12/yts;
        tfbmb(t) = md.results.TransientSolution(t).TotalFloatingBmbScaled * 1e12/yts;
    end
    write_ismip7_scalar(outdir, 'tendlibmassbfgr', -tgbmb, time_fl, time_bnds_fl, ...
        'tendency_of_land_ice_mass_due_to_basal_mass_balance', ...
        'Total BMB flux beneath grounded ice', 'kg s-1', meta);
    write_ismip7_scalar(outdir, 'tendlibmassbffl', -tfbmb, time_fl, time_bnds_fl, ...
        'tendency_of_land_ice_mass_due_to_basal_mass_balance', ...
        'Total BMB flux beneath floating ice', 'kg s-1', meta);

    % ---- tendlicalvf / tendligroundf --------------------------------------
    cf_native = zeros(1, nT); gl_native = zeros(1, nT);
    for t = 1:nT
        cf_native(t) = md.results.TransientSolution(t).IcefrontMassFluxLevelset * 1e12/yts;
        gl_native(t) = md.results.TransientSolution(t).GroundinglineMassFlux    * 1e12/yts;
    end
    if sign(sum(cf_native)) ~= sign(sum(cfflux_tot)) && any(cfflux_tot)
        warning('IcefrontMassFluxLevelset sign disagrees with cfflux_tot -- flipping to match.');
        cf_native = -cf_native;
    end
    if sign(sum(gl_native)) ~= sign(sum(glflux_tot)) && any(glflux_tot)
        warning('GroundinglineMassFlux sign disagrees with glflux_tot -- flipping to match.');
        gl_native = -gl_native;
    end
    write_ismip7_scalar(outdir, 'tendlicalvf', -cf_native, time_fl, time_bnds_fl, ...
        'tendency_of_land_ice_mass_due_to_calving', ...
        'Total calving flux', 'kg s-1', meta);
    write_ismip7_scalar(outdir, 'tendlifmassbf', -cf_native, time_fl, time_bnds_fl, ...
        'tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting', ...
        'Total frontal mass flux (calving + melt; no frontal melt in this run)', 'kg s-1', meta);
    write_ismip7_scalar(outdir, 'tendligroundf', -gl_native, time_fl, time_bnds_fl, ...
        'tendency_of_land_ice_mass_due_to_flux_at_grounding_line', ...
        'Total grounding line flux', 'kg s-1', meta);
end
