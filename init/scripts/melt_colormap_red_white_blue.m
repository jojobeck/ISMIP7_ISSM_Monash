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
