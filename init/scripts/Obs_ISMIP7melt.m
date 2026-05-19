function output = Obs_ISMIP7melt(X,Y,string),
 
% switch (oshostname()),
%   case {'ronne'}
%       rignotmelt='/home/ModelData/Antarctica/RignotMeltingrate/Ant_MeltingRate.nc';
%   case {'totten'}
%       rignotmelt='/totten_1/ModelData/Antarctica/RignotMeltingrate/Ant_MeltingRate.nc';
%   case {'thwaites','murdo','astrid'}
%       rignotmelt=['/home/seroussi/Data/Ant_MeltingRate.nc'];
%   otherwise
%       error('hostname not supported yet');
% end
pth_ismip7='/g/data/au88/ismip6/2300/forcings/ISMIP7/';
melt_data = 'AIS/parameterisations/ocean/meltobs/melt_paolo_err_adusumilli_ismip8km.nc';
melt_nc=fullfile(pth_ismip7, melt_data);
 
if nargin==2,
    string = 'melt_mean';
end
 
disp(['   -- ISMIP7 local melt rate: loading ' string]);
xdata = double(ncread(melt_nc,'x'));
ydata = double(ncread(melt_nc,'y'));
 
disp(['   -- ISMIP7 melt rate: loading' string]);
data  = double(ncread(melt_nc,string));
% Remove crazy fill values / NaNs
data(isnan(data)) = 0;                  
disp(['   -- ISMIP7 lmelt rate interpolating' string]);
output = InterpFromGrid(xdata,ydata,data,X(:),Y(:));
output = reshape(output,size(X,1),size(X,2));
