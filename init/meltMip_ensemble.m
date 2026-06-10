
function md=meltMIp(steps,j,loadonly)
    if ~exist('loadonly','var')
     % loadonly parameter does not exist, so default it to something
      loadonly = 0;
    end


    % loadonly = 1;
    addpath('./scripts');


    sec_to_year = 31556926;
    % org=organizer('repository',['./Models'],'prefix',['ISMIP6Antarctica_'],'steps',steps, 'color', '34;47;2'); 
    org=organizer('repository',['./Models/'],'prefix',['AIS_ISMIP7_'],'steps',steps, 'color', '34;47;2'); 
    clear steps;
    % loadonly=1;
    % InterpFromMeshToMesh2d

    % clim runs 10-12

 %1-2ka
 %%%%%%% tuning for basal melt
  %%%%%%%%%% create forcing file ocean
  inputmodel_path = './../AIS_1850/Models/AIS1850_thTHW_CollapseSSA.mat';
  basin_shelf_mask = './../Data/Ocean/basinid_iceshelf_extrap_davision_interpnearest.mat';
  inputmodel_relax = './Models/AIS_ISMIP7_Relaxed.mat';
  path_dir = '/g/data/au88/jb1863/SAEF/ISMIP7_ISSM_Monash/init/./../scripts '
  directory = 'Data/Tables';
  K_data = 0.25e-5 : 0.25e-5 : 3.0e-4
  sin_alpha=0.0029 %np.arscin(2.9e-3);
  sec_to_year = 31556926;
  %S0 = 34.5
  %gT_to_K = 2*abs(f_coriolis)*rho_sw/(rho_i*g*beta_coeff_lazero*S0*sin_alpha*sec_to_year)
  % calculated in pyhton from multimelt.constants import *
  gT_to_K= 1.2893727023667522e-08;
  if ~exist(directory, 'dir')
      mkdir(directory);
  end
    if perform(org,'melt_run')% {{{
        md=loadmodel(inputmodel_relax);
        m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor=(1./m).^2;
        
        % what is done her with the 1 ?
        md.transient.ismasstransport=1;
        md.transient.isstressbalance=0;
     
        md.inversion.iscontrol=0;
        md.transient.isgroundingline=0;
        % what is done her with the 1 ?                                                                                                                                                                                   
        md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
        md.outputdefinition.definitions={};
        md.timestepping.interp_forcing=0;
        md.transient.isthermal=0;   
        md.transient.issmb =0;
     
        md.timestepping.final_time=1;
        md.timestepping.time_step=1.;
        md.settings.output_frequency=1;
     
        md.transient.requested_outputs={'default','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate'};
     
        load './../preprocessed_data/Ocean/Basins/Imbie2_extrap_2km_BasinOnElements.mat';
        load './../preprocessed_data/Ocean/tf_depths.mat';
        load('./../preprocessed_data/Ocean/Clim/Clim_obs_TF.mat')
         %Set ISMIP6 basal melt rate parameters
        gamma0_median = K_data(j)/gT_to_K;
        unique_basinid = unique(basinid);
        delta_t = 0 * ones(1, length(unique_basinid)); %no correction in dT
     
        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique(basinid));
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = clim_obs_tf;
        md.basalforcings.islocal = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_median;
     
        md.miscellaneous.name=['bmelt_test' '_K_'  num2str(j)];
        clustername = 'gadi';
        cluster = set_cluster(clustername);
        md.cluster=cluster;
        md.settings.waitonlock=0; 
        md.verbose=verbose('solution',true,'module',true,'convergence',false);
        md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
        if loadonly
            savepth = ['Models/MeltMIP' md.miscellaneous.name];
            save(savepth,'md');
        end

        

end %}}}
    if perform(org,'create_BMB_gD'),% {{{
        name=['bmelt_test' '_K_'  num2str(j)];
        save_p = ['Models/MeltMIP' name];

        disp('loading model ...');
        disp(save_p);
        %get the model run
        md=loadmodel(save_p);
        % Create the folder path
        base_folder = 'Models/ModelNC'
        folder_name = sprintf('BMBPresentDay');
        full_folder_path = fullfile(base_folder, folder_name);
        % Check if the folder exists, if not, create it
        if ~exist(full_folder_path, 'dir')
            mkdir(full_folder_path);
        end
        % Extract ice BMB
        BMB= md.results.TransientSolution(end).BasalforcingsFloatingiceMeltingRate;%m/a
        % Load full 761x761 grid extent from bedmachine file
        bedmachineFile = '/home/565/jb1863/ismip6_2300/masks/af2_el_ismip6_ant_8km.nc';
        xFull = double(ncread(bedmachineFile, 'x'));
        yFull = double(ncread(bedmachineFile, 'y'));

        % Interpolate to 8kmx8km grid
        [dsBMB, xGrid, yGrid] = gridData(md, BMB, ...
            'xRange', [min(xFull) max(xFull)], ...
            'yRange', [min(yFull) max(yFull)]);
        dsBMB = dsBMB *md.materials.rho_ice;


        ncfile = fullfile(full_folder_path ,name +"_gridData_8km.nc");
        %check if file exists, delete
        if isfile(ncfile)
            delete(ncfile);
        end

        % Create full variables
        nccreate(ncfile,"x","Dimensions",{"x",length(xGrid)},"FillValue",NaN);
        nccreate(ncfile,"y","Dimensions",{"y",length(yGrid)},"FillValue",NaN);
        nccreate(ncfile, "melt_rate","Dimensions",{"x",length(xGrid),"y",length(yGrid)},"FillValue",NaN)
        % Write grid data to NetCDF
        ncwrite(ncfile,"y",yGrid);
        ncwrite(ncfile,"x",xGrid);

        ncwrite(ncfile,"melt_rate",transpose(dsBMB));
    end% }}}
    if perform(org,'create_BMB_gD_4km'),% {{{
        name=['bmelt_test' '_K_'  num2str(j)];
        save_p = ['Models/MeltMIP' name];

        disp('loading model ...');
        disp(save_p);
        %get the model run
        md=loadmodel(save_p);
        % Create the folder path
        base_folder = 'Models/ModelNC'
        folder_name = sprintf('BMBPresentDay_4km');
        full_folder_path = fullfile(base_folder, folder_name);
        % Check if the folder exists, if not, create it
        if ~exist(full_folder_path, 'dir')
            mkdir(full_folder_path);
        end
        % Extract ice BMB
        BMB= md.results.TransientSolution(end).BasalforcingsFloatingiceMeltingRate;%m/a
        % Interpolate to 4kmx4km grid
        [dsBMB, xGrid, yGrid] = gridData_4km(md, BMB);
        dsBMB = dsBMB *md.materials.rho_ice;


        ncfile = fullfile(full_folder_path ,name +"_gridData_4km.nc");
        %check if file exists, delete
        if isfile(ncfile)
            delete(ncfile);
        end

        % Create full variables
        nccreate(ncfile,"x","Dimensions",{"x",length(xGrid)},"FillValue",NaN);
        nccreate(ncfile,"y","Dimensions",{"y",length(yGrid)},"FillValue",NaN);
        nccreate(ncfile, "melt_rate","Dimensions",{"x",length(xGrid),"y",length(yGrid)},"FillValue",NaN)
        % Write grid data to NetCDF
        ncwrite(ncfile,"y",yGrid);
        ncwrite(ncfile,"x",xGrid);

        ncwrite(ncfile,"melt_rate",transpose(dsBMB));
    end% }}}
    if perform(org,'melt_run_OCeanModelling_Nw')% {{{
        md=loadmodel(inputmodel_relax);
        m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor=(1./m).^2;
        
        % what is done her with the 1 ?
        md.transient.ismasstransport=1;
        md.transient.isstressbalance=0;
     
        md.inversion.iscontrol=0;
        md.transient.isgroundingline=0;
        % what is done her with the 1 ?                                                                                                                                                                                   
        md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
        md.outputdefinition.definitions={};
        md.timestepping.interp_forcing=0;
        md.transient.isthermal=0;   
        md.transient.issmb =0;
     
        md.timestepping.final_time=1;
        md.timestepping.time_step=1.;
        md.settings.output_frequency=1;
     
        md.transient.requested_outputs={'default','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate'};
     
        load './../preprocessed_data/Ocean/Basins/Imbie2_extrap_2km_BasinOnElements.mat';
        load './../preprocessed_data/Ocean/tf_depths.mat';
        load('./../preprocessed_data/Ocean/Clim/Naughten_FESOM_ACCESS_warm_TF.mat')
         %Set ISMIP6 basal melt rate parameters
        gamma0_median = K_data(j)/gT_to_K;
        unique_basinid = unique(basinid);
        delta_t = 0 * ones(1, length(unique_basinid)); %no correction in dT
     
        md.basalforcings            = basalforcingsismip6(md.basalforcings);
        md.basalforcings.basin_id   = basinid;
        md.basalforcings.num_basins = length(unique(basinid));
        md.basalforcings.tf_depths  = tf_depths;
        md.basalforcings.tf         = clim_ocean_modelling_tf;
        md.basalforcings.islocal = 1;
        md.basalforcings.delta_t    = delta_t;
        md.basalforcings.gamma_0    = gamma0_median;
     
        md.miscellaneous.name=['bmelt_OceanModelling_Nw' '_K_'  num2str(j)];
        clustername = 'gadi';
        cluster = set_cluster(clustername);
        md.cluster=cluster;
        md.settings.waitonlock=0; 
        md.verbose=verbose('solution',true,'module',true,'convergence',false);
        md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
        if loadonly
            savepth = ['Models/MeltMIP' md.miscellaneous.name];
            save(savepth,'md');
        end

        

end %}}}

