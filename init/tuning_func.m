
function md=run_func(steps,loadonly)
    if ~exist('loadonly','var')
     % loadonly parameter does not exist, so default it to something
      loadonly = 0;
    end


    % loadonly = 1;
    addpath('./scripts');


    sec_to_year = 31556926;
    % org=organizer('repository',['./Models'],'prefix',['ISMIP6Antarctica_'],'steps',steps, 'color', '34;47;2'); 
    org=organizer('repository',['./Models'],'prefix',['AIS_ISMIP7_'],'steps',steps, 'color', '34;47;2'); 
    clear steps;
    % loadonly=1;
    % InterpFromMeshToMesh2d

    % clim runs 10-12

 %1-2ka
 %%%%%%% tuning for basal melt
  %%%%%%%%%% create forcing file ocean
  inputmodel_path = './../AIS_1850/Models/AIS1850_thTHW_CollapseSSA.mat';
  basin_shelf_mask = './../Data/Ocean/basinid_iceshelf_extrap_davision_interpnearest.mat';
  inputmodel_relax = '/Volumes/CrucialX8/computing/data/antarctica/AIS_1850_with_Justine/Models/Testfirst_test_Relax_20y_ocean_local_param.mat';
  path_dir = '/g/data/au88/jb1863/SAEF/ISMIP7_ISSM_Monash/init/./../scripts '
  directory = 'Data/Tables';
  if ~exist(directory, 'dir')
      mkdir(directory);
  end
  %OLD procedure,new new input data
    if perform(org,'Param'),% {{{

        md=loadmodel('./../ISMIP6/Models/ISMIP6Antarctica_Mesh.mat');
        md=setflowequation(md,'SSA','all');
        
        md=parameterize(md,'Par/Antartica_1995_Thwbedmap.par');

        pos=find(md.inversion.vel_obs>5e3);%setting too high velocities to 0, wont be used in inverision then!
        md.inversion.vel_obs(pos)=0;
        md.inversion.vx_obs(pos)=0;
        md.inversion.vy_obs(pos)=0;

        pos=find(md.mask.ocean_levelset<0); % floating ice needs no friction coefficient
        md.friction.coefficient(pos) = 0;
        md.materials.rheology_law = 'Paterson';

        savemodel(org,md);
    end%}}}
    if perform(org,'InversionB')% {{{

        md=loadmodel(org,'Param');

        md=setflowequation(md,'SSA','all');

        md.stressbalance.restol=0.01;
        md.stressbalance.reltol=0.1;
        md.stressbalance.abstol=NaN;

        md.inversion=m1qn3inversion(md.inversion);
        md.inversion.iscontrol=1;
        md.inversion.maxsteps=500;
        md.inversion.maxiter=10*100;
        %Cost functions
        md.inversion.cost_functions=[101 103 502];
        md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
        % md.inversion.cost_functions_coefficients(:,1)=2000;
        % md.inversion.cost_functions_coefficients(:,2)=40;
        % md.inversion.cost_functions_coefficients(:,3)=1e-16;
        % pos=find(md.inversion.vel_obs==0);
        % md.inversion.cost_functions_coefficients(pos,1)=0;
        % md.inversion.cost_functions_coefficients(:,2)=60;
        % pos=find(md.inversion.vel_obs==0);
        % md.inversion.cost_functions_coefficients(pos,2)=0;
        % or:
        md.inversion.cost_functions_coefficients(:,1)=500;
        md.inversion.cost_functions_coefficients(:,2)=60; 
        % md.inversion.cost_functions_coefficients(:,3)=1e-16;%good with 1km
        md.inversion.cost_functions_coefficients(:,3)=1e-18;%good with 1km
        % md.inversion.cost_functions_coefficients(:,1)=1;
        % md.inversion.cost_functions_coefficients(:,2)=1; 
        % md.inversion.cost_functions_coefficients(:,3)=1e-20;
        pos=find(md.inversion.vel_obs==0);
        md.inversion.cost_functions_coefficients(pos,1)=0;
        pos=find(md.inversion.vel_obs==0);
        md.inversion.cost_functions_coefficients(pos,2)=0;
        %make less weight for mask 19992
        mask_1992 =interpMeasureschange_withnans_mask1992(md,md.mesh.x,md.mesh.y,'nearest');
        my_isoline(md,mask_1992,'value',1,'output','./../Exp/vel1992.exp');
        pos1992 = find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/vel1992.exp',2));


%         %% PIG correction
%         posPig=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/PIG.exp',2));
%         md.inversion.cost_functions_coefficients(posPig,1)=15000;%good for 1km
%         md.inversion.cost_functions_coefficients(posPig,2)=1500;
%         %% posPigno1992 = setdiff(posPig,pos1992);
%         %% md.inversion.cost_functions_coefficients(posPigno1992,1)=0;
%         %% md.inversion.cost_functions_coefficients(posPigno1992,2)=0;

%         posPig2=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/THW.exp',2));
%         md.inversion.cost_functions_coefficients(posPig2,1)=15000;
%         md.inversion.cost_functions_coefficients(posPig2,2)=1500;
        % posPig2no1992 = setdiff(posPig2,pos1992);
        % md.inversion.cost_functions_coefficients(posPig2no1992,1)=0;
        % md.inversion.cost_functions_coefficients(posPig2no1992,2)=0;
        md.inversion.control_parameters={'MaterialsRheologyBbar'};
        md.inversion.min_parameters=paterson(273.15-1)*ones(md.mesh.numberofvertices,1);
        md.inversion.max_parameters=paterson(273.15-60)*ones(md.mesh.numberofvertices,1);

        md.groundingline.migration='SubelementMigration'; %VERY important at this resolution
        md.groundingline.friction_interpolation='SubelementFriction1';
        md.groundingline.melt_interpolation='SubelementMelt1';
        md.verbose=verbose('solution',false,'control',true);
        clustername = 'oshostname';
        cluster = set_cluster(clustername);
        md.cluster=cluster;
        mds=extract(md,md.mask.ocean_levelset<0);
        md.stressbalance.restol=0.003;
        mds.cluster = cluster;
        % md.cluster = cluster;
        % mds.settings.waitonlock=0; 
        mds=solve(mds,'sb');

        md.materials.rheology_B(mds.mesh.extractedvertices)=mds.results.StressbalanceSolution.MaterialsRheologyBbar;
        plot_after_tuning(4,mds);
        clear mds;

        savemodel(org,md);
    end% }}}
    if perform(org,'output_for_fitting_L_from_observation'),% {{{

        md =loadmodel('./../AIS_1850/Models/AIS1850_thTHW_InversionB.mat');
        p_vel = '/Users/jbec0008/SAEF/datasets/ICESat1_ICESat2_mass_change_updated_2_2021/';
        vel_set ='Smithgrounded_projbedmachine_mask_change.nc';
        bedm = fullfile(p_vel, vel_set);
        dh= interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'dh','nearest',bedm);

        resultTable = table( md.mask.ocean_levelset, md.inversion.vel_obs,md.geometry.bed, dh,...
        'VariableNames', {'thk_a_fl','vel','bed','dh_smith'});
        filename = fullfile(directory, sprintf('output_for_tunig_K_with_oversational_data.csv'));
     
        % Save the table to a CSV file with the basin name and number
        writetable(resultTable, filename);


    end% }}}
    if perform(org,'InversionFriction1'),% {{{

        md =loadmodel('./../AIS_1850/Models/AIS1850_thTHW_InversionB.mat');
        p_vel = '/Users/jbec0008/SAEF/datasets/ICESat1_ICESat2_mass_change_updated_2_2021/';
        vel_set ='Smithgrounded_projbedmachine_mask_change.nc';
        bedm = fullfile(p_vel, vel_set);
        mask_weight = interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'mask_change','nearest',bedm);

        %Remove some parts of the domain
        pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/IceAroundShelves.exp',2));
        md.mask.ice_levelset(pos)=-1;

        %Control general
        md.inversion=m1qn3inversion(md.inversion);
        md.inversion.iscontrol=1;
        md.verbose=verbose('solution',false,'control',true);

        %Find all elements that include ice and constrain some of them
        % pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
        % flags=zeros(md.mesh.numberofvertices,1); flags(md.mesh.elements(pos_e,:))=1;
        % pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/Constrain.exp',2) & flags);
        % md.stressbalance.spcvx(pos)=md.inversion.vx_obs(pos);
        % md.stressbalance.spcvy(pos)=md.inversion.vy_obs(pos);
        % md.stressbalance.spcvz(pos)=0;

        %Cost functions
        md.inversion.cost_functions=[101 103 501];
        md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
        % md.inversion.cost_functions_coefficients(:,1)=2000;
        % md.inversion.cost_functions_coefficients(:,2)=5;
        % md.inversion.cost_functions_coefficients(:,1)=1;
        % md.inversion.cost_functions_coefficients(:,2)=1;
        md.inversion.cost_functions_coefficients(:,1)=2000;
        md.inversion.cost_functions_coefficients(:,2)=40;
        md.inversion.cost_functions_coefficients(:,3)=2.5*1e-3;
        pos=find(md.mask.ice_levelset>0 | md.inversion.vel_obs==0);% no ice or with velocity 0
        md.inversion.cost_functions_coefficients(pos,1:2)=0;
        % posPig=find(mask_weight ==1);
        % md.inversion.cost_functions_coefficients(posPig,1)=300;
        % md.inversion.cost_functions_coefficients(posPig,2)=1000;

        %Correction for area aorund observed thickness change (smith et al)
        % n = readtable('./Data/Tables/determined_K_for_tunig_weighting_in_inversion.csv');
        n = readtable('./Data/Tables/determined_K_from_smithchanges_for_tuning_weighting_in_inversion.csv');
        md.inversion.cost_functions_coefficients(:,2)=960*n.k+40;

        md.inversion.cost_functions_coefficients(:,1)=2000-n.k*700;

        % md.inversion.cost_functions_coefficients(posPig,1)=0.5;
        % md.inversion.cost_functions_coefficients(posPig,2)=2000;

        %% PIG correction
        % posPig=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/PIG.exp',2));
        % md.inversion.cost_functions_coefficients(posPig,1)=300;
        % md.inversion.cost_functions_coefficients(posPig,2)=1000;
        % % % % %% md.inversion.cost_functions_coefficients(posPig,3)=.2*50^-2;
        % posPig2=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/THW.exp',2));
        % md.inversion.cost_functions_coefficients(posPig2,1)=300;
        % md.inversion.cost_functions_coefficients(posPig2,2)=1000;
        %% md.inversion.cost_functions_coefficients(posPig2,3)=.2*50^-2;

        %Controls
        md.groundingline.migration='SubelementMigration'; %VERY important at this resolution
        md.groundingline.friction_interpolation='SubelementFriction1';
        md.groundingline.melt_interpolation='SubelementMelt1';
        md.inversion.control_parameters={'FrictionCoefficient'};
        md.inversion.maxsteps=300;
        md.inversion.maxiter =300;
        md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
        md.inversion.max_parameters=500*ones(md.mesh.numberofvertices,1);
        % md.inversion.max_parameters=400*ones(md.mesh.numberofvertices,1);
        % md.inversion.max_parameters=1500*ones(md.mesh.numberofvertices,1);
        % md.inversion.max_parameters=500*ones(md.mesh.numberofvertices,1);
        md.inversion.control_scaling_factors=1;

        %Go solve/read_clusterjobfile
        md.miscellaneous.name='SSA_fritction_cmax_velonly';
        clustername = 'gadi';
        cluster = set_cluster(clustername);
        md.cluster=cluster;
        md.settings.waitonlock=0;
        md.cluster.time = 4*60;
        md=solve(md,'sb','runtimename', false,'loadonly',loadonly);
        md.cluster.time = 4*60;
        if loadonly,


            md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
            savemodel(org,md);
            plot_after_tuning(1,md);
        end
    end%}}}
    if perform(org,'Inversion2B'),% {{{

        md=loadmodel(org,'InversionFriction1');
        % use steadty state tmperature
        %Put rheology back
        % md.materials.rheology_law = 'Cuffey';
        % md.materials.rheology_law = 'Paterson';
        % md.materials.rheology_law = 'buddjacka';
        % pos = find(md.mask.ice_levelset<0);%everwhere where ice is
        % md.materials.rheology_B(pos)=cuffey(md.results.ThermalSolution.Temperature(pos));
        % pos = find(md.mask.ice_levelset<0);%everwhere where ice is
        % md.materials.rheology_B(pos)=paterson(md.results.ThermalSolution.Temperature(pos));

        %Control general

        md.stressbalance.restol=0.01;
        md.stressbalance.reltol=0.1;
        md.stressbalance.abstol=NaN;

        md.inversion=m1qn3inversion(md.inversion);
        md.inversion.iscontrol=1;
        md.inversion.maxsteps=500;
        md.inversion.maxiter=10*100;
        %Cost functions

        md.inversion.cost_functions=[101 103 502];
        md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
        % Helene:
        % md.inversion.cost_functions_coefficients(:,1)=2000;
        % md.inversion.cost_functions_coefficients(:,2)=40;
        % md.inversion.cost_functions_coefficients(:,3)=1e-16;
        % pos=find(md.inversion.vel_obs==0);
        % md.inversion.cost_functions_coefficients(pos,1)=0;
        % md.inversion.cost_functions_coefficients(:,2)=60;
        %pos=find(md.inversion.vel_obs==0);
        md.inversion.cost_functions_coefficients(:,1)=500;
        md.inversion.cost_functions_coefficients(:,2)=60; 
        % md.inversion.cost_functions_coefficients(:,3)=1e-16;%good with 1km
        md.inversion.cost_functions_coefficients(:,3)=1e-18;%good with 1km
        % md.inversion.cost_functions_coefficients(:,1)=1;
        % md.inversion.cost_functions_coefficients(:,2)=1; 
        % md.inversion.cost_functions_coefficients(:,3)=1e-20;
        pos=find(md.inversion.vel_obs==0);
        md.inversion.cost_functions_coefficients(pos,1)=0;
        pos=find(md.inversion.vel_obs==0);
        md.inversion.cost_functions_coefficients(pos,2)=0;

        %smith correction
        % n = readtable('./Data/Tables/determined_B_from_smithchanges_for_tuning_weighting_in_inversion.csv');
        % md.inversion.cost_functions_coefficients(:,1)=1000*n.k+500;
        % md.inversion.cost_functions_coefficients(:,2)=440*n.k+60;

        % PIG correction
        % posPig=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/PIG.exp',2));
        % % % md.inversion.cost_functions_coefficients(posPig,1)=1000;
        % % % md.inversion.cost_functions_coefficients(posPig,2)=1000;
        % % % md.inversion.cost_functions_coefficients(posPig,1)=5000;%good for 1km
        % md.inversion.cost_functions_coefficients(posPig,1)=15000;%good for 1km
        % % % md.inversion.cost_functions_coefficients(posPig2,3)=1e-25;%good with 1km
        % md.inversion.cost_functions_coefficients(posPig,2)=500;%good for 1km
        % % md.inversion.cost_functions_coefficients(posPig,2)=1500;
        % posPig2=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/THW.exp',2));
        % md.inversion.cost_functions_coefficients(posPig2,2)=50;
        % md.inversion.cost_functions_coefficients(posPig2,1)=15000;
        %md.inversion.control_parameters={'MaterialsRheologyBbar'};
        %% md.materials.rheology_B = buddjacka(273.15 - 10)*ones(md.mesh.numberofvertices,1);
        %% md.materials.rheology_law = 'BuddJacka';
        %% md.inversion.min_parameters=buddjacka(273.15-5)*ones(md.mesh.numberofvertices,1);
        %% md.inversion.max_parameters=buddjacka(273.15-50)*ones(md.mesh.numberofvertices,1);
        md.inversion.control_parameters={'MaterialsRheologyBbar'};
        md.inversion.min_parameters=paterson(273.15-1)*ones(md.mesh.numberofvertices,1);
        md.inversion.max_parameters=paterson(273.15-60)*ones(md.mesh.numberofvertices,1);

        md.groundingline.migration='SubelementMigration'; %VERY important at this resolution
        md.groundingline.friction_interpolation='SubelementFriction1';
        md.groundingline.melt_interpolation='SubelementMelt1';
        md.verbose=verbose('solution',false,'control',true);

       
        md.stressbalance.restol=0.003;
        % md.cluster = cluster;
        % mds.settings.waitonlock=0; 

        md.miscellaneous.name='SSA_InversionB2';
        %Go solve
        md.inversion.iscontrol=1;
        clustername = 'gadi';
        cluster = set_cluster(clustername);
        md.cluster=cluster;
        md.cluster.time = 4*60;
        md.settings.waitonlock=Inf;
        mds=extract(md,md.mask.ocean_levelset<0);
        % spc to initial verl not oberved
        % md.initialization.vx=md.results.StressbalanceSolution.Vx;
        % md.initialization.vy = md.results.StressbalanceSolution.Vy;
        % md.initialization.vel=md.results.StressbalanceSolution.Vel;
        pos_e = find(min(mds.mask.ice_levelset(mds.mesh.elements),[],2)<0);
        flags=zeros(mds.mesh.numberofvertices,1); flags(mds.mesh.elements(pos_e,:))=1;
        md2=extract(mds,flags); %only floating ice
        pos_icemargin=find(md2.mesh.vertexonboundary);%grounding line and ice front
        pos_icemargin_mds = md2.mesh.extractedvertices(pos_icemargin);%get pos for right mds mesh with ocean
        pos_boundary = find(mds.mesh.vertexonboundary);%grounding line and outer ocean grid cell
        [pos_gl, ia, ib] = intersect(pos_icemargin_mds,pos_boundary);% shared grounding line cells
        mds.stressbalance.spcvx(pos_gl) = mds.results.StressbalanceSolution.Vx(pos_gl);
        mds.stressbalance.spcvy(pos_gl) = mds.results.StressbalanceSolution.Vy(pos_gl);

        
        % mds.cluster = cluster;
        mds=solve(mds,'sb','runtimename',false);
        % mds=solve(mds,'sb','runtimename',false,'loadonly',1);

        md.materials.rheology_B(mds.mesh.extractedvertices)=mds.results.StressbalanceSolution.MaterialsRheologyBbar;
        savemodel(org,md);
        plot_after_tuning(4,mds);
    end% }}}
    if perform(org,'TESTRunInitSSA')% {{{

        % what is done her with the 1 ?
        if 1, md=loadmodel(org,'Inversion2B'); end
        % if 1, md=loadmodel(org,'change_thk_vel'); end
        % if 1, md=loadmodel(org,'InversionB2'); end
        
        % what is done her with the 1 ?
        m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor=(1./m).^2;

        md.inversion.iscontrol=0;
        md.transient.isgroundingline=1;
        % what is done her with the 1 ?
        md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
        md.outputdefinition.definitions={};
        md.timestepping.interp_forcing=0;
        md.transient.isthermal=0;	

        md.timestepping.final_time=1/12;
        md.timestepping.time_step=1/12.;
        md.settings.output_frequency=1;

        md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
            'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled'};

        %NoMeltOnPartiallyFloating results in blown up vertex in Wilkes Basin
        md.groundingline.migration = 'SubelementMigration';
        md.groundingline.friction_interpolation='SubelementFriction1';
        md.groundingline.melt_interpolation='SubelementMelt1';
        
        clustername = 'gadi';
        cluster = set_cluster(clustername);
        md.cluster=cluster;
        md.settings.waitonlock=Inf; 
        md.verbose=verbose('solution',true,'module',true,'convergence',true);
        md=solve(md,'tr');

        savemodel(org,md);
        plot_after_tuning(7,md);
    end %}}}
%Create ocn forcing    
%basins and melt data

if perform(org,'Extrusion'),% {{{

	% md=loadmodel(org,'InversionFriction1');
   md=loadmodel(org, 'flowlineextrap_vertices_processed_interp_field_B2');


	%update initial velocities
	md.initialization.vx=md.results.StressbalanceSolution.Vx;
	md.initialization.vy=md.results.StressbalanceSolution.Vy;

	%Here we want to spc the velocities that are on the inflow boundary because
	% that messes up the thermal model. To use outflow however, we need the segments
	% the easiest way to do that is to extract the model, get the inflow boundary
	% and spc the velocity there.
	pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
	flags=zeros(md.mesh.numberofvertices,1); flags(md.mesh.elements(pos_e,:))=1;
	md2=extract(md,flags);
	pos=find(md2.mesh.vertexonboundary & ~outflow(md2));
	md.stressbalance.spcvx(md2.mesh.extractedvertices(pos)) = 0;
	md.initialization.vx(  md2.mesh.extractedvertices(pos)) = 0;
	md.stressbalance.spcvy(md2.mesh.extractedvertices(pos)) = 0;
	md.initialization.vy(  md2.mesh.extractedvertices(pos)) = 0;
	clear md2;

	%FIXME: don't need if parameterize now
	md.basalforcings.floatingice_melting_rate = md.basalforcings.floatingice_melting_rate(:);

	%Extrude and set flow model as HO now
	md=extrude(md,10,1.1);
	md=setflowequation(md,'HO','all');

	savemodel(org,md);
end% }}}
if perform(org,'HO'),% {{{

	md=loadmodel(org,'Extrusion');

	md.inversion.iscontrol=0;
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md.miscellaneous.name=org.steps(org.currentstep).string;
    md.settings.waitonlock=0;
    md.cluster.time = 4*60;
    md=solve(md,'sb','runtimename',false,'loadonly',loadonly);
    if loadonly,
        md.initialization.vx=md.results.StressbalanceSolution.Vx;
        md.initialization.vy=md.results.StressbalanceSolution.Vy;
        md.initialization.vz=md.results.StressbalanceSolution.Vz;
        md.initialization.vel=md.results.StressbalanceSolution.Vel;
        md.initialization.pressure=md.results.StressbalanceSolution.Pressure;

        savemodel(org,md);

    end

end% }}}
if perform(org,'Thermal'),% {{{

	md=loadmodel(org,'HO');
    pos = md.basalforcings.geothermalflux <0;                                                                                                                                                  
    md.basalforcings.geothermalflux(pos) =0; 

	%Steadystate options
	md.timestepping.time_step = 0;
	md.inversion.iscontrol = 0;
	md.thermal.isenthalpy=1;

	%Spc some temps
	% pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/ThermalSpc.exp',2) & md.mask.ice_levelset>0);
	% md.thermal.spctemperature(pos)=250;


	%NOT NEEDED
	%pos = find(md.mask.groundedice_levelset<0 & md.mesh.vertexonbase);
	%md.thermal.spctemperature(pos)=270;
	%pos = find(md.mesh.vertexonbase);
	%md.thermal.spctemperature(pos)=md.materials.meltingpoint-md.materials.beta*md.materials.rho_ice*md.constants.g*md.geometry.thickness(pos);

    % clustername = 'oshostname';
    % cluster = set_cluster(clustername);
	% md.cluster=cluster;
    % md.settings.waitonlock=Inf; 
	% md=solve(md,'thermal');
    clustername = 'gadi';
    cluster = set_cluster(clustername);
    md.cluster=cluster;
    md.settings.waitonlock=Inf;
    md.cluster.time = 4*60;

    md=solve(md,'thermal','runtimename',false);
	md.initialization.temperature=md.results.ThermalSolution.Temperature;
	%md.basalforcings.melting_rate=md.results.SteadystateSolution.BasalforcingsMeltingRate;

	savemodel(org,md);
end% }}}
    if perform(org,'HO_output_for_fitting_L_from_observation'),% {{{

        md=loadmodel(org,'HO');

        p_vel = '/Users/jbec0008/SAEF/datasets/ICESat1_ICESat2_mass_change_updated_2_2021/';
        vel_set ='Smithgrounded_projbedmachine_mask_change.nc';
        bedm = fullfile(p_vel, vel_set);
        dh= interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'dh','nearest',bedm);

        resultTable = table( md.mask.ocean_levelset, md.geometry.bed, dh,...
        'VariableNames', {'thk_a_fl','bed','dh_smith'});
        filename = fullfile(directory, sprintf('output_for_tunig_K_with_smithchanges_data_HO.csv'));
     

        % Save the table to a CSV file with the basin name and number
        writetable(resultTable, filename);

    end% }}}

if perform(org,'InversionHO_'),% {{{

	md=loadmodel(org,'Thermal');

	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);
    	%Put rheology back
	md.materials.rheology_law = 'Cuffey';
	pos = find(md.mask.ocean_levelset>0);
	md.materials.rheology_B(pos)=cuffey(md.results.ThermalSolution.Temperature(pos));

	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=2000;
	md.inversion.cost_functions_coefficients(:,2)=40;
	pos=find(md.mask.ice_levelset>0 | md.inversion.vel_obs==0);% no ice or with velocity 0
	md.inversion.cost_functions_coefficients(pos,1:2)=0;


	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=200;
	md.inversion.maxiter =200;
	% md.inversion.maxsteps=20;
	% md.inversion.maxiter =20;
	md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=500*ones(md.mesh.numberofvertices,1);
	md.inversion.control_scaling_factors=1;
    n = readtable('./Data/Tables/determined_K_from_smithchanges_for_tuning_weighting_in_inversion_HO.csv');
    % md.inversion.cost_functions_coefficients(:,2)=1995*n.k+5;
    % md.inversion.cost_functions_coefficients(:,1)=1999.5-n.k;
    md.inversion.cost_functions_coefficients(:,2)=960*n.k+40;

    md.inversion.cost_functions_coefficients(:,1)=2000-n.k*700;


	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	% delta=[1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1];
	delta=[ 1e-4];
	%Go solve
	md.inversion.iscontrol=1;
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    if loadonly,
        disp('copy and plot');
        J= zeros(length(delta),length(md.inversion.cost_functions)+1);
        Jo=zeros(length(delta),1);
        Jr=zeros(length(delta),1);
        for num=1:length(delta)
            md.inversion.cost_functions_coefficients(:,3)=2.5*delta(num);
            % md_name=['InversionHO_L-CURVE-bestguess'];
            md_name=['Antarctica_HOdrag-Lcurve-' num2str(delta(num))];
            md.miscellaneous.name=md_name;
            mds=solve(md,'sb','runtimename',false,'loadonly',1);
            pvel = './Models';
            save_p= fullfile(pvel, md.miscellaneous.name);
            save(save_p,'mds'),
            plot_after_tuning(10,md);
        end
    else
        for num=1:length(delta)
            md.inversion.cost_functions_coefficients(:,3)=2.5*delta(num);
            md_name=['Antarctica_HOdrag-Lcurve-' num2str(delta(num))];
            md.miscellaneous.name=md_name;
            md.settings.waitonlock=0;
            mds=solve(md,'sb','runtimename', false);
        end
    end

end% }}}
if perform(org,'HO_deltabest'),% {{{


	md=loadmodel('./Models/Antarctica_HOdrag-Lcurve-0.0001');
    md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
    savemodel(org,md);

end% }}}
%%%% run iteration HO
if perform(org,'CollapseSSA'),% {{{

	% md=loadmodel('./Models/Cycle_thkTHW3_RunInitHO_inv2B.mat');
	md=loadmodel('./Models/Cycle_thkTHWcmax500_3_Inversion2B.mat');
    md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
    mds= md.collapse();
	mds=setflowequation(mds,'SSA','all');
    md1 =loadmodel(org,'Inversion2B.mat');
    mds.mesh =md1.mesh;
    savemodel(org,mds);
    % mds.mesh =md1.mesh;
    savemodel(org,mds);
end% }}}
%new c friciton interplation, and assigning SMB ans BMB IMSIP7 observartion
    if perform(org,'RunInitSSACollapse')% {{


        
        md=loadmodel('./Models/AIS_pd_cmaxCollapseSSA.mat');
        % what is done her with the 1 ?
        m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
        md.mesh.scale_factor=(1./m).^2;

        md.inversion.iscontrol=0;
        md.transient.isgroundingline=1;
        % what is done her with the 1 ?
        md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
        md.outputdefinition.definitions={};
        md.timestepping.interp_forcing=0;
        md.transient.isthermal=0;	

        md.timestepping.final_time=1/12;
        md.timestepping.time_step=1/12.;
        md.settings.output_frequency=1;

        md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
            'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled'};

        %NoMeltOnPartiallyFloating results in blown up vertex in Wilkes Basin
        md.groundingline.migration = 'SubelementMigration';
        md.groundingline.friction_interpolation='SubelementFriction1';
        md.groundingline.melt_interpolation='SubelementMelt1';
        md.miscellaneous.name = 'TESTbeforeISMIP7prep'
        
        clustername = 'gadi';
        cluster = set_cluster(clustername);
        md.cluster=cluster;
        md.settings.waitonlock=0; 
        md.settings.waitonlock=0;
        md.verbose=verbose('solution',true,'module',true,'convergence',true);
        md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
        if loadonly 
            savemodel(org,md);
        end
% plot_after_tuning(7,md);
    end %}}}
    %%%%%%%% ISMIP 7 input
    if perform(org,'constant_0.5cmean_from_Collapse'),% {{{



       md=loadmodel('./Models/AIS_pd_cmaxCollapseSSA.mat');
       p = 0.5;
       c_ground =mean(md.friction.coefficient(md.mask.ocean_levelset>0));
       md.friction.coefficient(md.mask.ocean_levelset<0)=p*c_ground;
       disp(c_ground);
       md.miscellaneous.name = 'CollapseSSACmean0.5'
       savemodel(org,md);

    end% }}}
    %Interpolate Forcing files
    %SMB mean 1995-2015
    %SMB hist 1995-2015
    %BMB obs
    %TF hist 1995 -2015
    % TF mean 1995 -2015
    %relax for 2 years with SMB mean and BMB obs

    %run BMB tuning iteration with tune_BMB_allK.m
    %determine settings K mean ,5 and 95 percentile
    %with settings K mean ,TF and SMB hist run hist
    %Correkt SMB for mass changes after Otosaka

