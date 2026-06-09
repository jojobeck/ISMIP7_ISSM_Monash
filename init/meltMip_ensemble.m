
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
    if perform(org,'Assign_Basins')% {{{
    md=loadmodel(inputmodel_relax);
   
    m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
    md.mesh.scale_factor=(1./m).^2;
    
    md.inversion.iscontrol=0;
    md.transient.isthermal=0;
    md.transient.ismasstransport=1;
    md.transient.isstressbalance=1;
    md.transient.isgroundingline=1;
    md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
   
    %load basin
    data_imbie='./../raw_data/ISMIP7/AIS/parameterisations/ocean/imbie2/basin_numbers_ismip2km_v2.nc'
    basin_datanc=[data_imbie];
% basinids                 = double(ncread(basin_datanc,'basinNumber'));%Beware starts with 0!
    %takenfrom interpbedmachine
   
    xdata            = double(ncread(basin_datanc,'x'));
    ydata            = double(ncread(basin_datanc,'y'));
   
    X=md.mesh.x;
    Y=md.mesh.y;
    offset=2;
   
    %%%%%%%%%%%%%%extrapolated ice shelf mask
    % offset=1;
      
    xmin=min(X(:)); xmax=max(X(:));
    posx=find(xdata<=xmax);
    if isempty(posx), posx=numel(xdata); end
    id1x=max(1,find(xdata>=xmin,1)-offset);
    id2x=min(numel(xdata),posx(end)+offset);
      
    ymin=min(Y(:)); ymax=max(Y(:));
    posy=find(ydata>=ymin);
    if isempty(posy), posy=numel(ydata); end
    id1y=max(1,find(ydata<=ymax,1)-offset);
    id2y=min(numel(ydata),posy(end)+offset);
      
    data  = double(ncread(basin_datanc,'basinNumber',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
    xdata=xdata(id1x:id2x);
    ydata=ydata(id1y:id2y);
    basinid_vetices = InterpFromGrid(xdata,ydata,data,double(X),double(Y),'nearest');
    %now get only values of nearest values,no means;
    bvu =unique(basinid_vetices);
    disp(size(bvu));
    % test_why_nogood_shelves    
    nbvu = bvu*0;
    for i= 1:size(bvu)
        bu = bvu(i);
        msk = bu ==basinid_vetices;
        nbvu(i)=sum(msk);
    end
    BasinoOnElements = basinid_vetices(md.mesh.elements);
    Basin_element =BasinoOnElements(1:end,1);
    Basin_element1 =BasinoOnElements(1:end,1);
    Basin_element2 =BasinoOnElements(1:end,2);
    Basin_element3 =BasinoOnElements(1:end,3);
    % disp(num2str( size ( unique(Basin_element1) )));
    % disp(num2str( size ( unique(Basin_element2) )));
    % disp(num2str( size ( unique(Basin_element3) )));
    for elem=1:md.mesh.numberofelements,
        be= BasinoOnElements(elem,1:end); %three corner values
        m = be(1) == be; %mask are all the same 
        m2 = be(2) == be; %mask are all the same 
        if sum(m)==3, %all the same values
            Basin_element(elem)= be(1);
        elseif sum(m)==2, % another similar to be(1)
            Basin_element(elem) = be(1);
        else,%check for 2 and 3
            if sum(m2)==2,% 2 and 3 are similar 
                Basin_element(elem)=be(2);
            else,%m2 =1 and m=1, all three are differnt
                %let us take the less represented ,but not one of the noshelf_davision!
                % bu1 =nbvu( be(1));
                % bu2 = nbvu(be(2));
                % bu3 = nbvu(be(3));
                bu1 = be(1);
                bu2 = be(2);
                bu3 = be(3);
                % Combine the extracted values into an array
                bu_values = [bu1, bu2, bu3];
                %check if there are in no shelfs
                [min_bu, min_idx] = min(bu_values);
      
                 % Get the associated be value
                associated_be = bu_values(min_idx);
      
                Basin_element(elem)=associated_be;%take smallest shelf value
            end
        end
    end
    basinid = Basin_element;
    uniq_bi = unique(Basin_element);
    disp('check');
    disp(size(uniq_bi));
    disp(uniq_bi);
    %original should for from 0 to 15   
    full_set = 1:16;
    
    
    % Display the results
    id_replace = full_set;
    for i =1:size(uniq_bi),
        val = uniq_bi(i);
        m = basinid==val;
        basinid(m)=i;
        id_replace(i)=val;
    end
    
    
    basin_vertices= basinid_vetices;
    
    
    
    save('./../preprocessed_data/Ocean/Basins/Imbie2_extrap_2km_BasinOnElements','basinid');%shifted by 1
    save('./../preprocessed_data/Ocean/Basins/Imbie2_extrap_2km_BasinOnVertices','basin_vertices');%not shiftef
    save('./../preprocessed_data/Ocean/Basins/Imbie2_id_replace','id_replace');%how to shift back
    




end %}}}
    %SMB hist 1995-2015:wq

