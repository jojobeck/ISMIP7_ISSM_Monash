function plot_init_BmbSmbC(md, name)
%PLOT_INIT_BMBSMBC Plot and save initial fields as separate figures
%
% Usage:
%   plot_init_BmbSmbC(md)
%   plot_init_BmbSmbC(md, 'my_figure')

if nargin < 2
    name = 'init_field';
end

if ~exist('./Figs','dir')
    mkdir('./Figs');
end

% 1) SMB mass balance
fig = figure('Visible','off');
plotmodel(md, ...
    'figure', fig, ...
    'data', md.smb.mass_balance, ...
    'title', 'SMB mass balance', ...
    'colorbar', 'on', ...
    'colorbartitle', 'm/yr');

outfile = ['./Figs/' name '_SMB_mass_balance.jpeg'];
print(fig, outfile, '-djpeg', '-r300');
fprintf('Saved figure: %s\n', outfile);
close(fig);

% 2) Friction coefficient
fig = figure('Visible','off');
plotmodel(md, ...
    'figure', fig, ...
    'data', md.friction.coefficient, ...
    'title', 'Friction coefficient', ...
    'colorbar', 'on');

outfile = ['./Figs/' name '_friction_coefficient.jpeg'];
print(fig, outfile, '-djpeg', '-r300');
fprintf('Saved figure: %s\n', outfile);
close(fig);

% 3) Floating ice melt rate
data = md.basalforcings.floatingice_melting_rate;

fig = figure('Visible','off');
plotmodel(md, ...
    'figure', fig, ...
    'data', md.basalforcings.floatingice_melting_rate, ...
    'title', 'Floating ice melt rate', ...
    'caxis', [0 3], ...
    'colorbar', 'on', ...
    'colorbartitle', 'm/yr');

outfile = ['./Figs/' name '_floatingice_melting_rate.jpeg'];
print(fig, outfile, '-djpeg', '-r300');
fprintf('Saved figure: %s\n', outfile);
close(fig);

end



function cmap = melt_colormap_red_white_blue()

n = 256;
n2 = n/2;

% Blue -> White
blue = [ ...
    linspace(0,1,n2)', ...
    linspace(0.2,1,n2)', ...
    ones(n2,1)];

% White -> Red
red = [ ...
    ones(n2,1), ...
    linspace(1,0,n2)', ...
    linspace(1,0,n2)'];

cmap = [blue; red];

end
