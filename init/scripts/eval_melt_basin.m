function melt_kgm2a = eval_melt_basin(TF_base_e, area_e, dT, gamma0, rho_ice)
% EVAL_MELT_BASIN  Analytically compute area-weighted mean basal melt rate
%   [kg m-2 a-1] for a set of floating ice elements, given their effective
%   base thermal forcing (back-calculated from an ISSM dT=+2 run) and a
%   scalar dT correction to test.
%
%   Uses the ISMIP6 local quadratic formula:
%       melt [m_ice/yr] = gamma0 [m/yr/°C²] * max(0, TF_base + dT)^2
%
%   Inputs:
%     TF_base_e : nElem×1  effective base TF per element [°C]
%                 (back-calculated as sqrt(melt_at_dT2/gamma0) - 2)
%     area_e    : nElem×1  triangle element areas [m²]
%     dT        : scalar   temperature correction to evaluate [°C]
%     gamma0    : scalar   heat-exchange coefficient [m/yr/°C²]
%     rho_ice   : scalar   ice density [kg/m³]
%
%   Output:
%     melt_kgm2a : scalar  area-weighted mean melt rate [kg m-2 a-1]
%                  (positive = melting)

    TF_eff     = TF_base_e + dT;                    % [°C], nElem×1
    melt_myr   = gamma0 * max(0, TF_eff).^2;        % [m_ice/yr], nElem×1
    melt_kgm2a = sum(melt_myr .* area_e) / sum(area_e) * rho_ice;
end
