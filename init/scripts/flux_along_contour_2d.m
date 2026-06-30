function [flux2d_kgm2s, flux_total_kgs] = flux_along_contour_2d(md, sol, contour_field, signed_normal, xGrid, yGrid, cellsize)
% FLUX_ALONG_CONTOUR_2D  Rasterise a V.N*H*L flux along the zero-contour of
%   contour_field onto the ISMIP6 8 km grid.
%
%   contour_field : MaskOceanLevelset (grounding line) or MaskIceLevelset
%                   (ice front) — zero-contour is the feature of interest.
%   signed_normal : true  → normal oriented grounded→floating (ligroundf);
%                   false → unsigned magnitude (licalvf).
%   xGrid, yGrid  : 1×nx and 1×ny grid coordinate vectors (metres).
%   cellsize      : grid spacing in metres (8000 for ISMIP7 8 km grid).
%
%   Returns flux2d_kgm2s (ny×nx, NaN where contour absent) and
%   flux_total_kgs (scalar AIS-wide integral).
%
%   Integration uses Simpson's rule (endpoints + true midpoint sample).
%   Thickness/Vx/Vy are P1 FEM fields, so H.(V.N) is an exact quadratic
%   along any straight segment; Simpson's rule integrates quadratics exactly.
%   Cells the contour never touches are NaN, consistent with ISMIP7 Table A1.

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
    idx = idx(good); x1 = x1(good); y1 = y1(good); x2 = x2(good); y2 = y2(good);
    dx = dx(good); dy = dy(good); L = L(good);

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
    if signed_normal
        % Orient normal grounded→floating so positive flux = mass crossing GL.
        eps_n = 1000;  % metres, per calc_GroundingLineFLux_transient_corr2.m
        phi_plus  = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm + eps_n*Nx, ym + eps_n*Ny);
        phi_minus = InterpFromMesh2d(elems, x, y, sol.MaskOceanLevelset, xm - eps_n*Nx, ym - eps_n*Ny);
        flip = (~isnan(phi_plus) & ~isnan(phi_minus) & (phi_plus > phi_minus));
        Nx(flip) = -Nx(flip); Ny(flip) = -Ny(flip);
    end

    % Midpoint sample for Simpson's rule (interpolated at true midpoint, not averaged).
    Vxm = InterpFromMesh2d(elems, x, y, sol.Vx,        xm, ym);
    Vym = InterpFromMesh2d(elems, x, y, sol.Vy,        xm, ym);
    Hm  = InterpFromMesh2d(elems, x, y, sol.Thickness, xm, ym);
    good3 = ~any(isnan([Vxm Vym Hm]), 2);

    f1 = H1 .* (Vx1.*Nx + Vy1.*Ny);
    f2 = H2 .* (Vx2.*Nx + Vy2.*Ny);
    fm = Hm .* (Vxm.*Nx + Vym.*Ny);
    % Fall back to trapezoidal for rare midpoint NaN (just outside mesh).
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
