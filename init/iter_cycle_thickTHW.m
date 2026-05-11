
% loadonly = 0;
addpath('./../scripts');
% i=0;
cycle =3;
steps= [1,2],loadonly=0;
steps =[2],loadonly=1;
steps = [3,4,5];loadonly =0;
steps = [5,6];loadonly = 1;
% steps =[6]
% % steps=[7]
% % % steps =[1];
disp(loadonly);

pref = ['Cycle_thkTHWcmax500_', num2str(cycle),'_'];

org=organizer('repository',['./Models'],'prefix',pref,'steps',steps, 'color', '34;47;2'); 
clear steps;



%Cycle
if perform(org,'Relax')% {{{
    md = loadmodel('./Models/AIS_pd_cmaxHO_deltabest.mat');
    if cycle ==1,
        disp('start iterarion cycle');
	    load './../Data/Atmosphere/MAR_racmo_smb_1995_2014_mean.mat';
        md_B =loadmodel('Models/AIS_pd_cmaxInversion2B.mat') ; 
        md_BMB = loadmodel('Models/bmelt_test_dt_iter_7.mat');% use low melt WAIS
        md_B.smb.mass_balance = MAR_racmo_smb_1995_2014_mean;
        md_B.basalforcings.floatingice_melting_rate =md_BMB.results.TransientSolution(end).BasalforcingsFloatingiceMeltingRate; 
        md_B_extrude =extrude(md_B,10,1.1);
        md.basalforcings.floatingice_melting_rate = md_B_extrude.basalforcings.floatingice_melting_rate;
        md.smb.mass_balance = md_B_extrude.smb.mass_balance;

    else,
        disp('get temperature and friction coefficnet from last iteration');
        num=cycle - 1;
        modelnm=['./Models/Cycle_thkTHWcmax500_', num2str(num),'_Inversion2B.mat'];
        % modelnm=['./Models/Cycle', num2str(num),'_Inversion2HO.mat'];
        md2=loadmodel(modelnm);
        pos =md.mask.ice_levelset <0;
        md.materials.rheology_B(pos)=md2.materials.rheology_B(pos);
        md.initialization.temperature=md2.initialization.temperature;
        md.friction.coefficient=md2.results.StressbalanceSolution.FrictionCoefficient;
    end;

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
	md.outputdefinition.definitions={};
	md.timestepping.interp_forcing=0;

	md.timestepping.final_time=0.2;
	md.timestepping.time_step=0.01;
	md.settings.output_frequency=1;

	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
		'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled'};

	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';


    %Set SMB Forcing Parameters
	md.miscellaneous.name=[pref,'Antarctica_relax'];
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
	md.verbose=verbose('solution',true,'module',true,'convergence',true);
	md.settings.waitonlock=Inf;
	md=solve(md,'tr','runtimename',false);
	% md=solve(md,'tr','runtimename',false,'loadonly',1);
	savemodel(org,md);

end %}}}
if perform(org,'InversionHO'),% {{{

	md=loadmodel(org,'Relax');
% Update geometries from previous relaxation step
    base = md.results.TransientSolution(end).Base;
    thickness = md.results.TransientSolution(end).Thickness;
    surface = md.results.TransientSolution(end).Surface;
    md.mesh.z = base + thickness ./ md.geometry.thickness .* (md.mesh.z - md.geometry.base);
    md.geometry.thickness = thickness;
    md.geometry.surface = surface;
    md.geometry.base = base;
    md.mask.ocean_levelset= md.results.TransientSolution(end).MaskOceanLevelset;

    % clear surface,thickness,base;
	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);

	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=2000;
	md.inversion.cost_functions_coefficients(:,2)=40;
	md.inversion.cost_functions_coefficients(:,3)=2.5*1e-4;% from L curve analysis
	pos=find(md.mask.ice_levelset>0 | md.inversion.vel_obs<0.1);% no ice or with velocity 0
	md.inversion.cost_functions_coefficients(pos,1:2)=0;

	md.miscellaneous.name=[pref,'Antarctica_HOfriciton'];

	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=200;
	md.inversion.maxiter =200;
	% md.inversion.maxsteps=2;
	% md.inversion.maxiter =2;
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

	%Go solve
	md.inversion.iscontrol=1;
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md.settings.waitonlock=0;
	% md.settings.waitonlock=Inf;
    md=solve(md,'sb','runtimename',false,'loadonly',loadonly);
    if loadonly,
        savemodel(org,md);
    end
end% }}}
if perform(org,'Thermal'),% {{{

	md=loadmodel(org,'InversionHO');
    md.initialization.vx=md.results.StressbalanceSolution.Vx;
    md.initialization.vy = md.results.StressbalanceSolution.Vy;
    md.initialization.vz = md.results.StressbalanceSolution.Vz;
    md.initialization.vel=md.results.StressbalanceSolution.Vel;
    % md.initialization.vel=md.results.StressbalanceSolution.Pr;
    pos = md.basalforcings.geothermalflux <0;                                                                                                                                                  
    md.basalforcings.geothermalflux(pos) =0; 

	%Steadystate options
	md.timestepping.time_step = 0;
	md.inversion.iscontrol = 0;
	md.thermal.isenthalpy=1;

    clustername = 'oshostname';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md.settings.waitonlock=Inf; 
	md=solve(md,'thermal');

	md.initialization.temperature=md.results.ThermalSolution.Temperature;
	%md.basalforcings.melting_rate=md.results.SteadystateSolution.BasalforcingsMeltingRate;

	savemodel(org,md);
end% }}}
if perform(org,'InversionB'),% {{{

	md=loadmodel(org,'Thermal');
    % use steadty state tmperature
	%Put rheology back
	md.materials.rheology_law = 'Cuffey';
	pos = find(md.mask.ice_levelset<0);%everwhere where ice is
	md.materials.rheology_B(pos)=cuffey(md.results.ThermalSolution.Temperature(pos));

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
    md.inversion.cost_functions_coefficients(:,1)=2000;
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


	% PIG correction
    posPig=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/PIG.exp',2));
    % md.inversion.cost_functions_coefficients(posPig,1)=1000;
    % md.inversion.cost_functions_coefficients(posPig,2)=1000;
    % md.inversion.cost_functions_coefficients(posPig,1)=5000;%good for 1km
    md.inversion.cost_functions_coefficients(posPig,1)=15000;%good for 1km
	% md.inversion.cost_functions_coefficients(posPig2,3)=1e-25;%good with 1km
    % md.inversion.cost_functions_coefficients(posPig,2)=500;%good for 1km
    md.inversion.cost_functions_coefficients(posPig,2)=1500;
    posPig2=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/THW.exp',2));
    md.inversion.cost_functions_coefficients(posPig2,1)=15000;
    md.inversion.cost_functions_coefficients(posPig2,2)=1500;
	md.inversion.control_parameters={'MaterialsRheologyBbar'};
	md.inversion.min_parameters=paterson(273.15-1)*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=paterson(273.15-50)*ones(md.mesh.numberofvertices,1);

	md.groundingline.migration='SubelementMigration'; %VERY important at this resolution
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
	md.verbose=verbose('solution',false,'control',true);

   
	md.stressbalance.restol=0.003;
    % md.cluster = cluster;
    % mds.settings.waitonlock=0; 

	md.miscellaneous.name=[pref,'Antarctica_HOInversionB'];
	%Go solve
	md.inversion.iscontrol=1;
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md.settings.waitonlock=Inf;
	mds=extract(md,md.mask.ocean_levelset<0);
    % mds.cluster = cluster;
	mds=solve(mds,'sb','runtimename',false);
	% mds=solve(mds,'sb','runtimename',false,'loadonly',1);

	md.materials.rheology_B(mds.mesh.extractedvertices)=mds.results.StressbalanceSolution.MaterialsRheologyBbar;
    savemodel(org,md);
end% }}}
if perform(org,'Inversion2HO'),% {{{
    md=loadmodel(org,'InversionB');
    % use steadty state tmperature

	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);

	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=2000;
	md.inversion.cost_functions_coefficients(:,2)=40;
	md.inversion.cost_functions_coefficients(:,3)=2.5*1e-4;% from L curve analysis
	pos=find(md.mask.ice_levelset>0 | md.inversion.vel_obs<0.1);% no ice or with velocity 0
	md.inversion.cost_functions_coefficients(pos,1:2)=0;

	md.miscellaneous.name=[pref,'Antarctica_HOfriciton2'];

	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=200;
	md.inversion.maxiter =200;
	% md.inversion.maxsteps=2;
	% md.inversion.maxiter =2;
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

	%Go solve
	md.inversion.iscontrol=1;
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md.settings.waitonlock=0;
    % md.settings.waitonlock=Inf;
    md=solve(md,'sb','runtimename',false,'loadonly',loadonly);
    if loadonly,
        savemodel(org,md);
    end
end% }}}
if perform(org,'Inversion2B'),% {{{

	md=loadmodel(org,'Inversion2HO');
    md.initialization.vx=md.results.StressbalanceSolution.Vx;
    md.initialization.vy = md.results.StressbalanceSolution.Vy;
    md.initialization.vz = md.results.StressbalanceSolution.Vz;
    md.initialization.vel=md.results.StressbalanceSolution.Vel;
    % use steadty state tmperature
	%Put rheology back
	md.materials.rheology_law = 'Cuffey';
	% pos = find(md.mask.ice_levelset<0);%everwhere where ice is
	% md.materials.rheology_B(pos)=cuffey(md.results.ThermalSolution.Temperature(pos));

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
    md.inversion.cost_functions_coefficients(:,1)=2000;
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


	% PIG correction
    posPig=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/PIG.exp',2));
    % md.inversion.cost_functions_coefficients(posPig,1)=1000;
    % md.inversion.cost_functions_coefficients(posPig,2)=1000;
    % md.inversion.cost_functions_coefficients(posPig,1)=5000;%good for 1km
    md.inversion.cost_functions_coefficients(posPig,1)=15000;%good for 1km
	% md.inversion.cost_functions_coefficients(posPig2,3)=1e-25;%good with 1km
    % md.inversion.cost_functions_coefficients(posPig,2)=500;%good for 1km
    md.inversion.cost_functions_coefficients(posPig,2)=1500;
    posPig2=find(ContourToNodes(md.mesh.x,md.mesh.y,'./../Exp/THW.exp',2));
    md.inversion.cost_functions_coefficients(posPig2,1)=15000;
    md.inversion.cost_functions_coefficients(posPig2,2)=1500;
	md.inversion.control_parameters={'MaterialsRheologyBbar'};
	md.inversion.min_parameters=paterson(273.15-1)*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=paterson(273.15-50)*ones(md.mesh.numberofvertices,1);

	md.groundingline.migration='SubelementMigration'; %VERY important at this resolution
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
	md.verbose=verbose('solution',false,'control',true);

   
	md.stressbalance.restol=0.003;
    % md.cluster = cluster;
    % mds.settings.waitonlock=0; 

	md.miscellaneous.name=[pref,'Antarctica_HOInversion2B'];
	%Go solve
	md.inversion.iscontrol=1;
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md.settings.waitonlock=Inf;
	mds=extract(md,md.mask.ocean_levelset<0);
    % spc to initial verl not oberved
    pos_e = find(min(mds.mask.ice_levelset(mds.mesh.elements),[],2)<0);
    flags=zeros(mds.mesh.numberofvertices,1); flags(mds.mesh.elements(pos_e,:))=1;
    md2=extract(mds,flags); %only floating ice
    pos_icemargin=find(md2.mesh.vertexonboundary);%grounding line and ice front
    pos_icemargin_mds = md2.mesh.extractedvertices(pos_icemargin);%get pos for right mds mesh with ocean
    pos_boundary = find(mds.mesh.vertexonboundary);%grounding line and outer ocean grid cell
    [pos_gl, ia, ib] = intersect(pos_icemargin_mds,pos_boundary);% shared grounding line cells
    mds.stressbalance.spcvx(pos_gl) = mds.initialization.vx(pos_gl);
    mds.stressbalance.spcvy(pos_gl) = mds.initialization.vy(pos_gl);
    mds.stressbalance.spcvz(pos_gl) = mds.initialization.vz(pos_gl);

    
    % mds.cluster = cluster;
	mds=solve(mds,'sb','runtimename',false);
	% mds=solve(mds,'sb','runtimename',false,'loadonly',1);

	md.materials.rheology_B(mds.mesh.extractedvertices)=mds.results.StressbalanceSolution.MaterialsRheologyBbar;
    savemodel(org,md);
end% }}}
%--------------------------test init:
if perform(org,'RunInitHO_inv2B')% {{{
	if 1, md=loadmodel(org,'Inversion2B'); end
    % mdt = loadmodel('./Models/Cycle1_Relax.mat');
    md.mask.ocean_levelset= md.results.TransientSolution(end).MaskOceanLevelset;
	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;

	md.inversion.iscontrol=0;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
	md.outputdefinition.definitions={};
	md.timestepping.interp_forcing=0;

	md.timestepping.final_time=1/24;
	md.timestepping.time_step=1/24.;
	md.settings.output_frequency=1;

	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
		'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled'};

	%Set friction and melt interpolation methods
	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
	
	
    clustername = 'oshostname';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md.settings.waitonlock=Inf; 
	md.verbose=verbose('solution',true,'module',true,'convergence',true);
	md=solve(md,'tr');

	savemodel(org,md);
end %}}}


