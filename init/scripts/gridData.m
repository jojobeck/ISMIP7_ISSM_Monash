function [data, xGrid, yGrid] = gridData(md, field, varargin)

    %Process options
    options = pairoptions(varargin{:});
    
    % Setup default mask (entire domain)
    xRange = getfieldvalue(options, 'xRange', [min(md.mesh.x(:)) max(md.mesh.x(:))]);
    yRange = getfieldvalue(options, 'yRange', [min(md.mesh.y(:)) max(md.mesh.y(:))]);
    
    % Get grid information
    bedmachineFile = '/home/565/jb1863/ismip6_2300/masks/af2_el_ismip6_ant_8km.nc';
    xGrid = double(ncread(bedmachineFile, 'x'));
    yGrid = double(ncread(bedmachineFile, 'y'));
    
    % Get X/Y Indices of model extent
    offset = 1;
    xmin = xRange(1); xmax = xRange(2);
    posx = find(xGrid <= xmax);
    if isempty(posx), posx=numel(xGrid); end
    id1x=max(1,find(xGrid>=xmin,1)-offset);
    id2x=min(numel(xGrid),posx(end)+offset);
    
    ymin = yRange(1); ymax = yRange(2);
    posy = find(yGrid >= ymin);
    if isempty(posy), posy=numel(yGrid); end
    id1y=max(1,find(yGrid<=ymax,1)-offset);
    id2y=min(numel(yGrid),posy(end)+offset);
    
    xGrid=xGrid(id1x:id2x);
    yGrid=yGrid(id1y:id2y);
    
    % Interpolare data to grid
    data = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y, field, xGrid, yGrid, NaN);
