function write_ismip7_2d(outdir, varname, data3d, xGrid, yGrid, ...
        time_vec, time_bnds, standard_name, long_name, units, meta)
% WRITE_ISMIP7_2D  Write a compliant ISMIP7 (x,y,t) NetCDF file.
%
%   write_ismip7_2d(outdir, varname, data3d, xGrid, yGrid, ...
%       time_vec, time_bnds, standard_name, long_name, units, meta)
%
%   data3d        : (ny, nx, nt) float array; NaN → ISMIP7 fill value
%   xGrid         : 1×nx x-coordinates in metres (must be ISMIP7 standard
%                   -3040000:8000:3040000, 761 pts, for compliance)
%   yGrid         : 1×ny y-coordinates in metres (same range)
%   time_vec      : 1×nt days since 1850-01-01 (standard / Gregorian calendar)
%   time_bnds     : nt×2 [lb ub] in the same units for FL variables, or []
%   standard_name : CF standard_name string, or '' if none assigned
%   long_name     : human-readable description
%   units         : CF units string
%   meta          : struct with experiment metadata — see make_ismip7_meta.m
%                   Required fields: experiment_id, set_counter, time_range,
%                     ESM_id, forcing_member_id, ISM_member_id
%                   Optional (defaults applied): source_id, ism_id,
%                     domain_id, group, model, contact_name, contact_email

    FILL_VALUE = single(9.969209968386869e+36);   % netCDF4 float32 default

    meta  = apply_meta_defaults(meta);
    fname = ismip7_filename(outdir, varname, meta);
    if exist(fname, 'file'), delete(fname); end

    nx = length(xGrid);
    ny = length(yGrid);
    nt = length(time_vec);

    data3d          = single(data3d);
    data3d(isnan(data3d)) = FILL_VALUE;

    % x
    nccreate(fname, 'x', 'Dimensions', {'x', nx}, 'Datatype', 'single');
    ncwrite(fname, 'x', single(xGrid(:)'));
    ncwriteatt(fname, 'x', 'units',         'm');
    ncwriteatt(fname, 'x', 'standard_name', 'projection_x_coordinate');
    ncwriteatt(fname, 'x', 'axis',          'X');

    % y
    nccreate(fname, 'y', 'Dimensions', {'y', ny}, 'Datatype', 'single');
    ncwrite(fname, 'y', single(yGrid(:)'));
    ncwriteatt(fname, 'y', 'units',         'm');
    ncwriteatt(fname, 'y', 'standard_name', 'projection_y_coordinate');
    ncwriteatt(fname, 'y', 'axis',          'Y');

    % time  (float32 as required by ISMIP7 checker)
    nccreate(fname, 'time', 'Dimensions', {'time', Inf}, 'Datatype', 'single');
    ncwrite(fname, 'time', single(time_vec));
    ncwriteatt(fname, 'time', 'units',         'days since 1850-01-01');
    ncwriteatt(fname, 'time', 'long_name',     'time');
    ncwriteatt(fname, 'time', 'standard_name', 'time');
    ncwriteatt(fname, 'time', 'axis',          'T');
    ncwriteatt(fname, 'time', 'calendar',      'standard');

    % time_bnds (FL variables only)
    if ~isempty(time_bnds)
        ncwriteatt(fname, 'time', 'bounds', 'time_bnds');
        % MATLAB nccreate reverses dimension order in the file, so
        % {'bnds',2,'time',Inf} → time_bnds(time,bnds) on disk — correct.
        nccreate(fname, 'time_bnds', ...
            'Dimensions', {'bnds', 2, 'time', Inf}, 'Datatype', 'single');
        ncwrite(fname, 'time_bnds', single(time_bnds'));
    end

    % data variable — dimensions (x, y, time) on disk
    nccreate(fname, varname, ...
        'Dimensions', {'x', nx, 'y', ny, 'time', Inf}, ...
        'Datatype',   'single', ...
        'FillValue',  FILL_VALUE);
    ncwrite(fname, varname, permute(data3d, [2 1 3]));   % (ny,nx,nt) → (nx,ny,nt)
    if ~isempty(standard_name)
        ncwriteatt(fname, varname, 'standard_name', standard_name);
    end
    ncwriteatt(fname, varname, 'long_name',     long_name);
    ncwriteatt(fname, varname, 'units',         units);
    ncwriteatt(fname, varname, 'missing_value', FILL_VALUE);
    ncwriteatt(fname, varname, 'coordinates',   'time y x');

    write_global_atts(fname, meta, varname);

    fprintf('Saved: %s\n', fname);
end

% ---------------------------------------------------------------- helpers
function meta = apply_meta_defaults(meta)
    if ~isfield(meta, 'source_id'),     meta.source_id     = 'Monash'; end
    if ~isfield(meta, 'ism_id'),        meta.ism_id        = 'ISSM';   end
    if ~isfield(meta, 'domain_id'),     meta.domain_id     = 'AIS';    end
    if ~isfield(meta, 'group'),         meta.group         = 'Monash'; end
    if ~isfield(meta, 'model'),         meta.model         = 'ISSM';   end
    if ~isfield(meta, 'contact_name'),  meta.contact_name  = 'Johanna Beckmann'; end
    if ~isfield(meta, 'contact_email'), meta.contact_email = 'johanna.beckmann@monash.edu'; end
end

function fname = ismip7_filename(outdir, varname, meta)
    fname = fullfile(outdir, sprintf('%s_%s_%s_%s_%s_%s_%s_%s_%s_%s.nc', ...
        varname, meta.domain_id, meta.source_id, meta.ism_id, ...
        meta.ISM_member_id, meta.ESM_id, meta.forcing_member_id, ...
        meta.experiment_id, meta.set_counter, meta.time_range));
end

function write_global_atts(fname, meta, varname)
    ncwriteatt(fname, '/', 'Conventions',    'CF-1.7 ISMIP7');
    ncwriteatt(fname, '/', 'institution',    meta.source_id);
    ncwriteatt(fname, '/', 'source',         [meta.model ' (' meta.ism_id ')']);
    ncwriteatt(fname, '/', 'group',          meta.group);
    ncwriteatt(fname, '/', 'model',          meta.model);
    ncwriteatt(fname, '/', 'contact_name',   meta.contact_name);
    ncwriteatt(fname, '/', 'contact_email',  meta.contact_email);
    ncwriteatt(fname, '/', 'experiment_id',  meta.experiment_id);
    ncwriteatt(fname, '/', 'ESM_id',         meta.ESM_id);
    ncwriteatt(fname, '/', 'forcing_member_id', meta.forcing_member_id);
    ncwriteatt(fname, '/', 'ISM_member_id',  meta.ISM_member_id);
    ncwriteatt(fname, '/', 'set_counter',    meta.set_counter);
    ncwriteatt(fname, '/', 'time_range',     meta.time_range);
    ncwriteatt(fname, '/', 'variable_id',    varname);
    ncwriteatt(fname, '/', 'crs',            'epsg:3031');
    ncwriteatt(fname, '/', 'history', ...
        ['Created ' datestr(now, 'yyyy-mm-dd') ' by write_ismip7_2d.m']);
end
