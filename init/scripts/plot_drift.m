function plot_drift(md)
%PLOT_DRIFT Plot relaxation drift diagnostics:
%  1) Thickness change between saved timesteps
%  2) VAF drift
%  3) Grounded area drift

if ~exist('./Figs','dir')
    mkdir('./Figs');
end

sol = md.results.TransientSolution;
nsteps = length(sol);

time = NaN(nsteps,1);
vaf  = NaN(nsteps,1);
ga   = NaN(nsteps,1);

for i = 1:nsteps
    time(i) = sol(i).time;

    if isfield(sol(i),'IceVolumeAboveFloatation')
        vaf(i) = sol(i).IceVolumeAboveFloatation;
    end

    if isfield(sol(i),'GroundedArea')
        ga(i) = sol(i).GroundedArea;
    end
end

% Thickness change between saved outputs
dH = NaN(nsteps-1,1);
time_dH = NaN(nsteps-1,1);

for i = 2:nsteps
    Hnow  = sol(i).Thickness;
    Hprev = sol(i-1).Thickness;

    dH(i-1) = mean(abs(Hnow - Hprev),'omitnan');
    time_dH(i-1) = time(i);
end

% Relative drift from first saved output
vaf_drift = vaf - vaf(1);
ga_drift  = ga  - ga(1);

% Convert units
vaf_drift_km3 = vaf_drift / 1e9;  % m3 to km3
ga_drift_km2  = ga_drift  / 1e6;  % m2 to km2


% ============================================================
% 1) Thickness adjustment
% ============================================================
fig = figure('Visible','off');
plot(time_dH,dH,'LineWidth',2);
xlabel('Time [yr]');
ylabel('Mean DelatH between outputs [m]');
title('Relaxation thickness adjustment');
grid on;

outfile = './Figs/relaxation_dH.jpeg';
print(fig,outfile,'-djpeg','-r300');
fprintf('Saved figure: %s\n', outfile);
close(fig);


% ============================================================
% 2) VAF drift
% ============================================================
fig = figure('Visible','off');
plot(time,vaf_drift_km3,'LineWidth',2);
xlabel('Time [yr]');
ylabel('Delta VAF [km^3]');
title('Volume above floatation drift');
grid on;

outfile = './Figs/relaxation_VAF_drift.jpeg';
print(fig,outfile,'-djpeg','-r300');
fprintf('Saved figure: %s\n', outfile);
close(fig);


% ============================================================
% 3) Grounded area drift
% ============================================================
fig = figure('Visible','off');
plot(time,ga_drift_km2,'LineWidth',2);
xlabel('Time [yr]');
ylabel('Delta grounded area [km^2]');
title('Grounded area drift');
grid on;

outfile = './Figs/relaxation_grounded_area_drift.jpeg';
print(fig,outfile,'-djpeg','-r300');
fprintf('Saved figure: %s\n', outfile);
close(fig);

end
