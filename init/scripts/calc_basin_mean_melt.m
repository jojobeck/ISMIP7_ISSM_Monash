function [melt_basin_kgm2a, unique_basins, total_bmb_Gtyr] = calc_basin_mean_melt(melt_v, md, basinid, rho_ice)
% CALC_BASIN_MEAN_MELT  Area-weighted mean basal melt rate [kg m-2 a-1] and
%   total basal mass balance [Gt yr-1] per IMBIE2 basin from per-vertex
%   ISSM melt output.
%
%   Only floating elements (ocean_levelset < 0 for at least 2 of 3 vertices)
%   are included.  Basins with no floating elements return 0.
%
%   Inputs:
%     melt_v  : nVertices×1  BasalforcingsFloatingiceMeltingRate [m_ice/yr]
%     md      : ISSM model   (needs mesh.x/y/elements, mask.ocean_levelset)
%     basinid : nElements×1  IMBIE2 basin ID, MATLAB 1-based (element-level)
%     rho_ice : scalar        ice density [kg m-3]
%
%   Outputs:
%     melt_basin_kgm2a : nBasins×1  area-weighted mean melt [kg m-2 a-1]
%     unique_basins    : nBasins×1  MATLAB basin IDs (1-based)
%     total_bmb_Gtyr   : nBasins×1  total basal mass loss [Gt yr-1]

    unique_basins = unique(basinid);
    nBasins       = length(unique_basins);

    % --- triangle element areas -------------------------------------------
    x = md.mesh.x;  y = md.mesh.y;
    e = md.mesh.elements;   % nElements × 3
    x1=x(e(:,1)); y1=y(e(:,1));
    x2=x(e(:,2)); y2=y(e(:,2));
    x3=x(e(:,3)); y3=y(e(:,3));
    elem_area = 0.5 * abs(x1.*(y2-y3) + x2.*(y3-y1) + x3.*(y1-y2));

    % --- floating mask (element level: 2 of 3 vertices floating) ----------
    float_v   = md.mask.ocean_levelset < 0;
    float_e   = (float_v(e(:,1)) + float_v(e(:,2)) + float_v(e(:,3))) >= 2;

    % --- element-mean melt from vertex values [m_ice/yr] ------------------
    melt_e = (melt_v(e(:,1)) + melt_v(e(:,2)) + melt_v(e(:,3))) / 3;

    % --- per-basin stats --------------------------------------------------
    melt_basin_kgm2a = zeros(nBasins, 1);
    total_bmb_Gtyr   = zeros(nBasins, 1);
    for b = 1:nBasins
        basin_b   = unique_basins(b);
        elem_mask = (basinid == basin_b) & float_e;
        if any(elem_mask)
            A_b = elem_area(elem_mask);
            m_b = melt_e(elem_mask);
            flux = sum(m_b .* A_b) * rho_ice;   % kg/yr
            melt_basin_kgm2a(b) = flux / sum(A_b);   % kg m-2 a-1
            total_bmb_Gtyr(b)   = flux / 1e12;        % Gt yr-1
        end
        % basins with no floating elements stay at 0
    end
end
