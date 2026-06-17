
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

  % =====================================================================
  % Step map (run TF builders first, then melt runs, then gridding):
  %   TF builders : 1 Obs_clim_TF  2 OceanModelling_clim_TF  3 ObsData_clim_TF
  %   melt runs   : 4 melt_run  5 melt_run_OceanModelling  6 melt_run_ObsData
  %   gridding    : 7 create_BMB_gD  8 create_BMB_gD_4km
  %                 9 create_BMB_gD_OceanModelling  10 create_BMB_gD_ObsData
  %   gamma_0     : 11 save_gamma0_local  (run AFTER run_parameter_selection.py)
  % The tf cell arrays are size [1,1,numel(tf_depths)] (3rd dim = DEPTH);
  % each tf{i} = [tf_at_vertices ; t] with t the (single) time row.
  % =====================================================================

    %% ---------------------------------------------------------------
    %% TF builders
    %% ---------------------------------------------------------------
    if perform(org,'Obs_clim_TF')% {{{
        md=loadmodel(inputmodel_relax);
        m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor=(1./m).^2;


        ocean_climTF = './../raw_data/ISMIP7/AIS/meltMIP/OI_Climatology_ismip8km_60m_tf_extrap.nc'
        tfnc           = [ocean_climTF];
        x_n            = double(ncread(tfnc,'x'));
        y_n            = double(ncread(tfnc,'y'));
        tf_data        = double(ncread(tfnc,'tf'));
        z_data         = double(ncread(tfnc,'z'));

        %Build tf cell array (3rd dim = depth)
        tf = cell(1,1,size(tf_data,3));
        t=1;

        for i=1:size(tf_data,3)  %Iterate over depths
          temp_matrix=[];
          temp_tfdata=InterpFromGridToMesh(x_n,y_n,tf_data(:,:,i)',md.mesh.x,md.mesh.y,0);
          temp_tfdata = max(temp_tfdata, 0);  % ISSM basalforcingsismip6 requires tf>=0 (clamp cold-cavity TF<0)
          temp_matrix = [temp_matrix temp_tfdata];
          temp_matrix = [temp_matrix ; t];%need time axis , here =1
          tf(:,:,i)={temp_matrix};  % write depth slice i (was tf(:,:) -> overwrote all depths)
        end
        clim_obs_tf = tf;
        save('./../preprocessed_data/Ocean/Clim/Clim_obs_TF.mat','clim_obs_tf','-v7.3');

    end %}}}
    if perform(org, 'OceanModelling_clim_TF')% {{{

        md = loadmodel(inputmodel_relax);
        m = ((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor = (1./m).^2;

        ocean_files = {
            'Naughten_FESOM_ACCESS_cold_TF.nc', ...
            'Naughten_FESOM_ACCESS_warm_TF.nc',  ...
            'Mathiot_NEMO_cold_v2_TF.nc',       ...
            'Mathiot_NEMO_warm_v2_TF.nc'        ...
        };

        base_path  = './../raw_data/ISMIP7/AIS/parameterisations/ocean/ocean_modelling_data/';
        save_path  = './../preprocessed_data/Ocean/Clim/';

        for f = 1:length(ocean_files)

            ocean_climTF = [base_path ocean_files{f}];

            x_n     = double(ncread(ocean_climTF, 'x'));
            y_n     = double(ncread(ocean_climTF, 'y'));
            tf_data = double(ncread(ocean_climTF, 'thermal_forcing'));
            z_data  = double(ncread(ocean_climTF, 'z'));

            % Build tf cell array (3rd dim = depth)
            clim_ocean_modelling_tf = cell(1,1,size(tf_data,3));
            t = 1;

            for i = 1:size(tf_data,3)  % Iterate over depths
                temp_matrix = [];
                temp_tfdata = InterpFromGridToMesh(x_n, y_n, tf_data(:,:,i)', md.mesh.x, md.mesh.y, 0);
                temp_tfdata = max(temp_tfdata, 0);  % ISSM basalforcingsismip6 requires tf>=0 (clamp cold-cavity TF<0)
                temp_matrix = [temp_matrix temp_tfdata];
                temp_matrix = [temp_matrix ; t];  % need time axis, here =1
                clim_ocean_modelling_tf(:,:,i) = {temp_matrix};  % write depth slice i
            end

            % Get filename without extension and save
            [~, fname, ~] = fileparts(ocean_climTF);
            save([save_path fname '.mat'], 'clim_ocean_modelling_tf', '-v7.3');

            fprintf('Saved: %s.mat\n', fname);
        end

    end% }}}
    if perform(org,'ObsData_clim_TF')% {{{
        % Build thermal-forcing cell arrays from the Amundsen Sea ocean
        % OBSERVATIONS (Dutrieux), one mat file per observation year.
        % Mirrors Obs_clim_TF. ISMIP7 README recommends years 2009 and 2012
        % for Pine Island (obs_ensemble is indexed by year; the pig/dotson
        % split is a spatial mask applied later in the tuning notebook).
        obs_years = [2009, 2012];

        base_path = './../raw_data/ISMIP7/AIS/parameterisations/ocean/ocean_observations_data/';
        save_path = './../preprocessed_data/Ocean/Clim/';

        md = loadmodel(inputmodel_relax);
        m = ((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor = (1./m).^2;

        for y = 1:length(obs_years)
            year = obs_years(y);
            ocean_obsTF = [base_path 'Dutrieux_ismip8km_60m_tf_' num2str(year) '.nc'];

            x_n     = double(ncread(ocean_obsTF, 'x'));
            y_n     = double(ncread(ocean_obsTF, 'y'));
            tf_data = double(ncread(ocean_obsTF, 'tf'));

            % Build tf cell array (3rd dim = depth)
            clim_obs_data_tf = cell(1,1,size(tf_data,3));
            t = 1;
            for i = 1:size(tf_data,3)  % Iterate over depths
                temp_matrix = [];
                temp_tfdata = InterpFromGridToMesh(x_n, y_n, tf_data(:,:,i)', md.mesh.x, md.mesh.y, 0);
                temp_tfdata = max(temp_tfdata, 0);  % ISSM basalforcingsismip6 requires tf>=0 (clamp cold-cavity TF<0)
                temp_matrix = [temp_matrix temp_tfdata];
                temp_matrix = [temp_matrix ; t];  % need time axis, here =1
                clim_obs_data_tf(:,:,i) = {temp_matrix};  % write depth slice i
            end

            save([save_path 'Dutrieux_ismip8km_60m_tf_' num2str(year) '.mat'], 'clim_obs_data_tf', '-v7.3');
            fprintf('Saved: Dutrieux_ismip8km_60m_tf_%d.mat\n', year);
        end
    end% }}}

    %% ---------------------------------------------------------------
    %% Melt runs
    %% ---------------------------------------------------------------
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
    if perform(org,'melt_run_OceanModelling')% {{{
        % Generic ocean-modelling melt run over all four forcings.
        % Outputs feed the python step that builds cold/warm_ensemble.nc.
        % Re-runs overwrite existing model outputs (no skip guard).
        forcings = {
            struct('tf','Naughten_FESOM_ACCESS_warm_TF.mat','tag','Nw'), ...
            struct('tf','Naughten_FESOM_ACCESS_cold_TF.mat','tag','Nc'), ...
            struct('tf','Mathiot_NEMO_warm_v2_TF.mat',      'tag','Mw'), ...
            struct('tf','Mathiot_NEMO_cold_v2_TF.mat',      'tag','Mc')  ...
        };

        for f = 1:length(forcings)
            tag     = forcings{f}.tag;
            savepth = ['Models/MeltMIPbmelt_OceanModelling_' tag '_K_' num2str(j)];

            md=loadmodel(inputmodel_relax);
            m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
            md.mesh.scale_factor=(1./m).^2;

            md.transient.ismasstransport=1;
            md.transient.isstressbalance=0;
            md.inversion.iscontrol=0;
            md.transient.isgroundingline=0;
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
            load(['./../preprocessed_data/Ocean/Clim/' forcings{f}.tf]);  % loads clim_ocean_modelling_tf
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

            md.miscellaneous.name=['bmelt_OceanModelling_' tag '_K_'  num2str(j)];
            clustername = 'gadi';
            cluster = set_cluster(clustername);
            md.cluster=cluster;
            md.settings.waitonlock=0;
            md.verbose=verbose('solution',true,'module',true,'convergence',false);
            md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
            if loadonly
                save(savepth,'md');
            end
        end

    end %}}}
    if perform(org,'melt_run_ObsData')% {{{
        % Melt runs forced by Amundsen Sea ocean observations, per year, per K.
        % Re-runs overwrite existing model outputs (no skip guard).
        obs_years = [2009, 2012];

        for y = 1:length(obs_years)
            year    = obs_years(y);
            savepth = ['Models/MeltMIPbmelt_ObsData_' num2str(year) '_K_' num2str(j)];

            md=loadmodel(inputmodel_relax);
            m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
            md.mesh.scale_factor=(1./m).^2;

            md.transient.ismasstransport=1;
            md.transient.isstressbalance=0;
            md.inversion.iscontrol=0;
            md.transient.isgroundingline=0;
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
            load(['./../preprocessed_data/Ocean/Clim/Dutrieux_ismip8km_60m_tf_' num2str(year) '.mat']);  % loads clim_obs_data_tf
            %Set ISMIP6 basal melt rate parameters
            gamma0_median = K_data(j)/gT_to_K;
            unique_basinid = unique(basinid);
            delta_t = 0 * ones(1, length(unique_basinid)); %no correction in dT

            md.basalforcings            = basalforcingsismip6(md.basalforcings);
            md.basalforcings.basin_id   = basinid;
            md.basalforcings.num_basins = length(unique(basinid));
            md.basalforcings.tf_depths  = tf_depths;
            md.basalforcings.tf         = clim_obs_data_tf;
            md.basalforcings.islocal = 1;
            md.basalforcings.delta_t    = delta_t;
            md.basalforcings.gamma_0    = gamma0_median;

            md.miscellaneous.name=['bmelt_ObsData_' num2str(year) '_K_'  num2str(j)];
            clustername = 'gadi';
            cluster = set_cluster(clustername);
            md.cluster=cluster;
            md.settings.waitonlock=0;
            md.verbose=verbose('solution',true,'module',true,'convergence',false);
            md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
            if loadonly
                save(savepth,'md');
            end
        end
    end% }}}

    %% ---------------------------------------------------------------
    %% Gridding
    %% ---------------------------------------------------------------
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
    if perform(org,'create_BMB_gD_OceanModelling'),% {{{
        % Grid the ocean-modelling melt runs (all four forcings) to the 8km
        % ISMIP grid, mirroring create_BMB_gD. Outputs feed the python step
        % that builds cold_ensemble.nc / warm_ensemble.nc.
        tags = {'Nw','Nc','Mw','Mc'};

        base_folder = 'Models/ModelNC';
        folder_name = sprintf('BMBOceanModelling');
        full_folder_path = fullfile(base_folder, folder_name);
        if ~exist(full_folder_path, 'dir')
            mkdir(full_folder_path);
        end

        % Load full 761x761 grid extent from bedmachine file
        bedmachineFile = '/home/565/jb1863/ismip6_2300/masks/af2_el_ismip6_ant_8km.nc';
        xFull = double(ncread(bedmachineFile, 'x'));
        yFull = double(ncread(bedmachineFile, 'y'));

        for t = 1:length(tags)
            tag    = tags{t};
            name   = ['bmelt_OceanModelling_' tag '_K_'  num2str(j)];
            save_p = ['Models/MeltMIP' name];

            if ~isfile([save_p '.mat'])
                fprintf('Missing %s.mat, skipping gridding\n', save_p);
                continue;
            end

            ncfile = fullfile(full_folder_path, name + "_gridData_8km.nc");
            if isfile(ncfile)
                delete(ncfile);
            end

            disp('loading model ...'); disp(save_p);
            md = loadmodel(save_p);

            % Extract ice BMB
            BMB = md.results.TransientSolution(end).BasalforcingsFloatingiceMeltingRate;

            % Interpolate to 8kmx8km grid
            [dsBMB, xGrid, yGrid] = gridData(md, BMB, ...
                'xRange', [min(xFull) max(xFull)], ...
                'yRange', [min(yFull) max(yFull)]);
            dsBMB = dsBMB * md.materials.rho_ice;

            % Create variables
            nccreate(ncfile,"x","Dimensions",{"x",length(xGrid)},"FillValue",NaN);
            nccreate(ncfile,"y","Dimensions",{"y",length(yGrid)},"FillValue",NaN);
            nccreate(ncfile, "melt_rate","Dimensions",{"x",length(xGrid),"y",length(yGrid)},"FillValue",NaN);
            ncwrite(ncfile,"y",yGrid);
            ncwrite(ncfile,"x",xGrid);
            ncwrite(ncfile,"melt_rate",transpose(dsBMB));
        end
    end% }}}
    if perform(org,'create_BMB_gD_ObsData'),% {{{
        % Grid the observation-forced melt runs to the 8km ISMIP grid, per year.
        % Outputs feed creating_bmbMIP_ObsData_ensemble.py -> obs_ensemble.nc.
        obs_years = [2009, 2012];

        base_folder = 'Models/ModelNC';
        folder_name = sprintf('BMBObsData');
        full_folder_path = fullfile(base_folder, folder_name);
        if ~exist(full_folder_path, 'dir')
            mkdir(full_folder_path);
        end

        bedmachineFile = '/home/565/jb1863/ismip6_2300/masks/af2_el_ismip6_ant_8km.nc';
        xFull = double(ncread(bedmachineFile, 'x'));
        yFull = double(ncread(bedmachineFile, 'y'));

        for y = 1:length(obs_years)
            year   = obs_years(y);
            name   = ['bmelt_ObsData_' num2str(year) '_K_'  num2str(j)];
            save_p = ['Models/MeltMIP' name];

            if ~isfile([save_p '.mat'])
                fprintf('Missing %s.mat, skipping gridding\n', save_p);
                continue;
            end

            ncfile = fullfile(full_folder_path, name + "_gridData_8km.nc");
            if isfile(ncfile)
                delete(ncfile);
            end

            disp('loading model ...'); disp(save_p);
            md = loadmodel(save_p);

            BMB = md.results.TransientSolution(end).BasalforcingsFloatingiceMeltingRate;

            [dsBMB, xGrid, yGrid] = gridData(md, BMB, ...
                'xRange', [min(xFull) max(xFull)], ...
                'yRange', [min(yFull) max(yFull)]);
            dsBMB = dsBMB * md.materials.rho_ice;

            nccreate(ncfile,"x","Dimensions",{"x",length(xGrid)},"FillValue",NaN);
            nccreate(ncfile,"y","Dimensions",{"y",length(yGrid)},"FillValue",NaN);
            nccreate(ncfile, "melt_rate","Dimensions",{"x",length(xGrid),"y",length(yGrid)},"FillValue",NaN);
            ncwrite(ncfile,"y",yGrid);
            ncwrite(ncfile,"x",xGrid);
            ncwrite(ncfile,"melt_rate",transpose(dsBMB));
        end
    end% }}}
    if perform(org,'save_gamma0_local'),% {{{
        % Convert K values selected by run_parameter_selection.py to gamma_0
        % and save for use in ISSM projection runs.
        % Requires K_selected.mat written by BMB_tuning_python/run_parameter_selection.py.
        load('./../preprocessed_data/Ocean/K_selected.mat');  % K_mode, K_5th, K_50th, K_95th

        gamma0_local      = K_mode / gT_to_K;
        gamma0_local_5th  = K_5th  / gT_to_K;
        gamma0_local_50th = K_50th / gT_to_K;
        gamma0_local_95th = K_95th / gT_to_K;

        save_path = './../preprocessed_data/Ocean/gamma0_local.mat';
        save(save_path, 'gamma0_local', 'gamma0_local_5th', 'gamma0_local_50th', 'gamma0_local_95th', '-v7.3');
        fprintf('Saved: %s\n', save_path);
        fprintf('  gamma0_local (mode) = %.4e  (K_mode = %.4e)\n', gamma0_local,      K_mode);
        fprintf('  gamma0_local_5th    = %.4e  (K_5th  = %.4e)\n', gamma0_local_5th,  K_5th);
        fprintf('  gamma0_local_50th   = %.4e  (K_50th = %.4e)\n', gamma0_local_50th, K_50th);
        fprintf('  gamma0_local_95th   = %.4e  (K_95th = %.4e)\n', gamma0_local_95th, K_95th);
    end% }}}

