function HYLANDS_out = HYLANDS(DEM,T,p,varargin)

% TTLEM_HyLands: Hybrid Landscape evolution model to simulate the impact of
% landslides and landslide-derived sediment
%
% Syntax
% 
%     HYLANDS_out = HYLANDS(DEM,T,p)
%     HYLANDS_out = HYLANDS(DEM,T,p,varargin)
%
% Description
%
%     HYLANDS runs the HyLands landscape evolution model. Parameters that
%     control the LEM as well as output options are supplied to the
%     function in the structure array p, which is generated with
%     HYLANDS_set function (see documentation of HYLANDS_set, or type help
%     HYLANDS_set).
%
% Input arguments: required
%
%     DEM     initial surface (digital elevation model) (GRIDobj)
%     T       structure, containing uplift data:
%             *spatial: GRIDobj of same size as DEM containing spatially
%             variable uplift rates [m/y] (GRIDobj)
%             *other properties can be added in future releases (e.g.
%             temporal uplift variability)
%     p       structure array with parameter definitions (see HYLANDS_set)
% Input arguments: optional
%     iniSed  grid of same size as DEM containing spatially
%             variable sediment thickness [m] (GRIDobj)
%
% Output arguments:
%     HYLANDS_out.I_Domain	mask indicating the extent of the modelled
%                           domain (not including boundary nodes)
%                           HYLANDS_out.bedrock=bedrock;
%     HYLANDS_out.DA        Drainage area
%     HYLANDS_out.sediment	sediment thickness [m]
%     HYLANDS_out.bedrock	bedrock elevation [m]
%     HYLANDS_out.BORDER	flow_BORDER
%     HYLANDS_out.runTime   current time
%     HYLANDS_out.S       final flow object, representing the river network 
%     HYLANDS_out.DSY_out_Riv dissolved and suspended sediment yield per
%     time step
%     HYLANDS_out.DSY_out_Hill dissolved and suspended sediment yield 
%     originating from landsliding per time step
%   
%     if bedrock landslding takes place, the following will be output: 
%         HYLANDS_out.LS_prop_bed_store: structure containing frequency,
%         magnitude and volume of all landslides
%         HYLANDS_out.slidePlane_Bed: GRIDobj indicating landslide areas
%         HYLANDS_out.LS_Total_E: GRIDobj indicating total amount of
%         landslide erosion
%
% Example: 
%     %% Clear environment
%     clearvars
%     clc
%     close all
% 
%     %% Run mode from existing DEM
%     load('DEMs\DEM_Yarlung_Mini.mat','DEM');
%     iniSed=GRIDobj(DEM);
%     T.spatial  = GRIDobj(DEM);
%     p.TimeSpan = 100;
%     p.TimeStep = 10;
%     p.ploteach = 1;
%     p.FlowBC='open';
%     p.BC_Type='Bed_Sed_Open';
%     p.LS_Bed   = true;
%     p.t_LS=5e2;
%     p.Sc_fixed = .8;
%     p   = HYLANDS_set(p);
%     HYLANDS_Out = HYLANDS(DEM,T,p,'iniSed',iniSed);
% 
%     %% Plot Landslide characterstics
%     dx2=DEM.cellsize*DEM.cellsize;
%     bins=logspace(1,5,15);
%     figure; histogram([HYLANDS_Out.LS_prop_bed_store.Size]*dx2,bins);
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
%     xlabel('\bfArea of landslides, m^2')
%     ylabel('\bfNumber of landslides (#)')
%
% See also: HYLANDS_set
%
% =========================================================================
% Papers to cite when using HyLands:
%
% * HyLands: Campforts B., Shobe M.C., et al. : HyLands 1.0: a Hybrid
% Landscape evolution model to simulate the impact of landslides and
% landslide-derived sediment on landscape evolution. Discussion paper in
% Geoscientific Model Development,
% https://geoscientific-model-development.net
%
% * SPACE: Shobe, C. M., Tucker, G. E., & Barnhart, K. R. (2017). The SPACE 1.0
% model: a Landlab component for 2-D calculation of sediment transport,
% bedrock erosion, and landscape evolution. Geoscientific Model
% Development, 10(12), 4577–4604. https://doi.org/10.5194/gmd-10-4577-2017
%
% * TTLEM: Campforts, B., Schwanghart, W., & Govers, G. (2017). 
% Accurate simulation of transient landscape evolution
% by eliminating numerical diffusion: the TTLEM 1.0 model. Earth Surface
% Dynamics, 5(1), 47–66. https://doi.org/10.5194/esurf-5-47-2017
%
% Other relevant references:
%
% * Carretier, S., Martinod, P., Reich, M., & Godderis, Y. (2016).
% Modelling sediment clasts transport during landscape evolution. Earth
% Surface Dynamics, 4(1), 237–251. https://doi.org/10.5194/esurf-4-237-2016
%
% * Densmore, A. L., Ellis, M. A., & Anderson, R. S. (1998). Landsliding
% and the evolution of normal-fault-bounded mountains. Journal of
% Geophysical Research: Solid Earth, 103(B7), 15203–15219.
% https://doi.org/10.1029/98JB00510
%
% =========================================================================
%
% Author:   Benjamin Campforts (benjamin.campforts@gfz-potsdam.de)
%
% Date:     15. March, 2020

cpu_runTime=tic;
if p.verbose
    Progress.iter=0;
    Progress.stage='Initialize';
    Progress.time=toc(cpu_runTime);
    Progress=struct2table(Progress);
    writetable(Progress,[p.resultsdir  'Progress_Table_' p.fileprefix '.csv']);
end


%% Parse inputs
inp = inputParser;
inp.FunctionName = 'HYLANDS';

% required
addRequired(inp,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(inp,'T');
addRequired(inp,'p',@(x) isstruct(x));
% optional
addParameter(inp,'iniSed',[],@(x)isa(x,'GRIDobj'));
parse(inp,DEM,T,p,varargin{:});
clear varargin

%% Fetch data and add boundary nodes
% During model runs, the domain will be surrounded by a boundary nodes wich
% can be open or closed .

%DEM
DEM = inp.Results.DEM;
DEM.Z = double(DEM.Z);
% Check if input DEM is double precision
dx=DEM.cellsize;
dx2=dx*dx;
H1=pad(DEM);
H1.Z=extendBoundary(DEM.Z);
I_Domain=GRIDobj(H1)+1;
I_Domain.Z(getBoundaryNodes(I_Domain.Z))=0;
I_Domain.Z=logical(I_Domain.Z);

% Distance grid
[X,Y] = meshgrid((1:size(H1.Z,2))*dx,(1:size(H1.Z,1))*dx);
X=single(X);
Y=single(Y);

%Uplift
T = inp.Results.T;

%Parameters
p = inp.Results.p;

% Check input Parameters
dt       = p.TimeStep;
tspan    = p.TimeSpan;
rr = rem(tspan,dt);
if rr ~= 0
    tspan = tspan + (dt-rem(tspan,dt));
    warning('TopoToolbox:HYLANDS',...
        ['Total simulation time was increased from %.1f to %.1f \n' ...
        'to be a multiple of the time step'],p.TimeStep,tspan);
    p.ts   = tspan;
end

% number of iterations in the main loop
nriter = round(tspan/dt);

% Check entered Boundary conditions bedrock
if ~(isscalar(p.BC_BedDirVal)||isequal(size(p.BC_BedDirVal),size(getBoundaryNodes(H1.Z))))
    error('Wrong input size of p.BC_upl_Dir')
end

% Check entered Boundary conditions sediment
if ~(isscalar(p.BC_SedDirVal)||isequal(size(p.BC_SedDirVal),size(getBoundaryNodes(H1.Z))))
    error('Wrong input size of p.BC_upl_Dir')
end

% Bedrock
bedrock=H1;
bedrock.Z(bedrock.Z<0)=0;

sediment=GRIDobj(H1);
if ~isempty(inp.Results.iniSed)
    sediment.Z=double(extendBoundary(inp.Results.iniSed.Z));
end

Sc_Map=GRIDobj(H1)+p.Sc_fixed;

%-----------------------------Clear input memory---------------------------
clear inp


%% Additional data structures required during model run

%------------------------Slide plains for landslides-----------------------
%Deepseated
if p.LS_Bed
    slidePlane_Bed=uint16(zeros(H1.size));
    LS_Total_E=uint16(zeros(H1.size));
end
% Flow boundary conditions
%if empty; all boundary nodes are open; otherwise, indicate the boundary
%nodes which have a closed boundary and where no water can leave the area.

BC.BedDirVal=p.BC_BedDirVal;
BC.SedDirVal=p.BC_SedDirVal;
BC.type=p.BC_Type;
[flow_BORDER,open_BC_nodes] = getflowBC(H1,p,'');

BC.boundaryNodesAll=getBoundaryNodes(H1.Z);
BC.open_nodes=open_BC_nodes;
BC.flow_BORDER=flow_BORDER;

%% Calculate parameters being used during model runs
%--------------------------River incision parametres-----------------------
K_sed_L=p.K_sed;
K_bed_L=p.K_bed;
K=max(K_bed_L,K_sed_L);
pureIncision=(p.Ff==1)&&(p.Ff_Hill==1);

%% visualize output
if p.plotOut
    if p.plotSed
        %%
        load color_DEMPoster.mat        
        load color_heat.mat
        H1_Domain=crop(H1,I_Domain);
        Sed_Domain=crop(sediment,I_Domain);
        Bed_Domain=crop(bedrock,I_Domain);
        
        scrsz = get(groot,'ScreenSize');
        sedplot=figure('Position',[5 70 scrsz(3)/1.01 scrsz(4)/3],'color','white');
        sax1=subplot(1,4,1);
        cm1=color_DEMPoster;
        imagesc(H1_Domain)
        colorbar
        colormap(sax1,cm1)
        sax2=subplot(1,4,2);
        cm2=flipud(color_heat);
        imagesc((Sed_Domain))
        colormap(sax2,cm2)
        colorbar
        hold off
        sax2=subplot(1,4,3:4);
        FD=FLOWobj(H1_Domain,'preprocess','carve','mex',true);
        S=STREAMobj(FD,'minarea',100);
        S1=klargestconncomps(S);
        St=trunk(S1);
        hGUI.stopbutton = uicontrol('Style', 'togglebutton', 'String', 'Stop',...
            'units','normalized',...
            'Position', [0.02 .87 .15 .1],...
            'Min',false,'Max',true,'value',false);
        hGUI.timedisp   = uicontrol('Style', 'text', 'String', '0',...
            'units','normalized','backgroundcolor','white',...
            'Position', [0.3 .85 .15 .1]);
        drawnow
        
    else
        H1_Domain=crop(H1,I_Domain);
        hGUI  = preparegui(H1_Domain);
        set(hGUI.ax.time,'xlim',[0 tspan]);
        drawnow
    end
end
hGUIContinue=true;

%% save directory
% does directory exist?
if ~exist(p.resultsdir,'dir')
    mkdir(p.resultsdir)
end
resultsdir = fullfile(p.resultsdir,filesep);

%% Initiate simulation
% keep track of total time
t        = 0;
iter     = 0;

%% Update/Load neighbourhood matrix for DEM
if~exist([resultsdir p.FlowDir '_' p.FlowDirHill '_' ...
        'res_' num2str(round(H1.cellsize)) '_dim' num2str(H1.size(1)) ...
        'X' num2str(H1.size(2)) '_m_neighbours.mat'],'file')
    if p.verbose
        disp('Update neighbourhood matrix for new DEM...')
    end
    
    if strcmp(p.FlowDir,'single')&& strcmp(p.FlowDirHill,'single')
        [nb.M] =neighbourIndices(H1.Z,'M');
    else
        [nb.M,nb.ic,nb.icd] =neighbourIndices(H1.Z,'M&ind');
    end
    save([resultsdir p.FlowDir '_' p.FlowDirHill '_' ...
        'res_' num2str(round(H1.cellsize)) '_dim' num2str(H1.size(1)) ...
        'X' num2str(H1.size(2)) '_m_neighbours.mat'],'nb','-v7.3');
else
    if p.verbose
        disp('Load neighbourhood matrix ...')
    end
    load([resultsdir p.FlowDir '_' p.FlowDirHill '_' ...
        'res_' num2str(round(H1.cellsize)) '_dim' num2str(H1.size(1)) ...
        'X' num2str(H1.size(2)) '_m_neighbours.mat'],'nb');
end

%% Mass balance
% Make row for Mb with ca. the size of the model run
MB_Riv_record=[];
MB_row=zeros(p.TimeSpan/p.TimeStep,1);

%% Soil equation, BOTH PO and alpha in cm!!!!!, output in cm!
mod_Domain_M=logical(I_Domain.Z);
mod_Domain=I_Domain.Z(:);
locationsD=find(mod_Domain_M); %indices of modelled domain

%% ************************* Save *************************
%Save initial bedrock and sediment elevation if storeage is required
if mod(iter,p.saveeach)==0
    save([resultsdir p.fileprefix '_res_' num2str(round(H1.cellsize)) 'm_bedrock_' num2str(round(t)) '.mat'],'bedrock');
    if p.K_sed~=0
        save([resultsdir p.fileprefix '_res_' num2str(round(H1.cellsize)) 'm_sediment_' num2str(round(t)) '.mat'],'sediment');
    end
    
end

%% ***************************** Qs output  *******************************
Qs_out_Riv=nan(nriter,1);
DSY_out_Riv=nan(nriter,1);
Qs_out_Hill=nan(nriter,1);
DSY_out_Hill=nan(nriter,1);
totalRelief=nan(nriter,1);

%% ************************* Model Run *************************
while iter < nriter && hGUIContinue
   
    
    %% *************************** Verbose ***************************
    if p.verbose
        clc
        disp(['HYLANDS: iteration ' num2str(iter) ' of ' num2str(nriter+1)])
        Progress.iter=iter;
        Progress.stage='start';
        Progress.time=toc(cpu_runTime);
        writetable(Progress,[p.resultsdir  'Progress_Table_' p.fileprefix '.csv']);
        
    end
    %% *************************** Update elevation ***************************
    H1=bedrock+sediment;
    
    %% *************************** Track time ***************************
    t = t + dt;
    iter = iter + 1;
    
    
    %% *************************** Mass balance ***************************
    if p.checkMB||  p.checkMB
        bed_Before=bedrock.Z(mod_Domain);
        sed_Before=sediment.Z(mod_Domain);
        if p.checkMB
            checkMB(iter,bed_Before,sed_Before,0,dt,bedrock.Z, sediment.Z,...
                [],[],mod_Domain,p.phi,p.MB_lim,'initial step',p.dispMB);
        end
    end
    
       
    
    %% *************************** Vertical Uplift ***************************
    switch p.U_type
        case {'spatialVar';'uniform'}
            upl=double(T.spatial.Z);
        otherwise
            error('Uplift not correctly inserted!')
    end
    
    %-----------------------uplift in boundary nodes-----------------------
    upl=addNanBorder(upl);
    upl(getBoundaryNodes(upl))=0;
    
    
    % DO Uplift
    bedrock.Z = bedrock.Z +dt.*upl;    
    H1=bedrock+sediment;
    
    %% **************************** Runoff ************************************
    R_myr=p.R_base;
        

    %%  *************************** Update drainage direction ***************************
    
    if iter==1||strcmp(p.DrainDir,'variable')            
            if p.verbose
                disp('Identify river cells...')
                Progress.stage='Identify river cells...';
                Progress.time=toc(cpu_runTime);
                writetable(Progress,[p.resultsdir  'Progress_Table_' p.fileprefix '.csv']);
            end
            
            HT=H1;
            if p.FlowBC
                HT = HT+flow_BORDER;
            end
            
            HT=fillsinks(HT);
            Lakes=HT>H1;
            FD=FLOWobj(HT,'mex',true,'preprocess','none');
            W0=GRIDobj(H1)+dx2;
            W0.Z(BC.boundaryNodesAll)=0;
            W0.Z(BC.open_nodes)=dx2;
            DA  = flowacc(FD,W0);
            dx_ik = double(sqrt((X(FD.ix)-X(FD.ixc)).^2 + (Y(FD.ix)-Y(FD.ixc)).^2));
            
            switch p.FlowDirHill
                case 'single'
                    % OK
                case {'multi','DInf'}
                    FD_single=FD;
                    [FD_Hill,~,DA_Hill,dx_ik_Hill] ...
                        = updateDrainDir_hylands(p.FlowDirHill,HT,X,Y,BC.boundaryNodesAll,nb);
            end            
        clear HT;
    end
        
    %% *************************** Hybrid river Incision ***************************
    if p.verbose
        disp('Hybrid river Incision...')
        Progress.stage='Hybrid river Incision...';
        Progress.time=toc(cpu_runTime);
        writetable(Progress,[p.resultsdir  'Progress_Table_' p.fileprefix '.csv']);
    end
    
    if p.K_bed~=0||p.K_sed~=0
        [bedrock,sediment,lvSedTot,tot_DSY] =hybrid_RivInc(...
            bedrock, sediment,Lakes, DA.Z,R_myr,...
            dx2,K,K_bed_L,K_sed_L,p.m,p.n,...
            FD.ix,FD.ixc,dx_ik,BC.boundaryNodesAll,mod_Domain,...
            p.TimeStep,p.phi,p.Ff,p.V,p.V_Lakes,p.H_star,pureIncision,...
            p.omega_cr,p.omega_cs);
    end
    
    if p.checkMB
        hybridDef=dx2*(sum(bed_Before-bedrock.Z(mod_Domain))+...
            sum(sed_Before-sediment.Z(mod_Domain))*(1-p.phi)+...
            sum(dt.*upl(mod_Domain)))-(tot_DSY+lvSedTot);
        MB_R=(tot_DSY+lvSedTot+hybridDef)/dx2;
        MB_Riv_record.hybridDef(iter)=hybridDef;
        MB_Riv_record.lvSedTot(iter)=lvSedTot;
        MB_Riv_record.tot_DSY(iter)=tot_DSY;
        
        checkMB(iter,bed_Before,sed_Before,upl,dt,bedrock.Z, sediment.Z,...
            MB_R,[],mod_Domain,p.phi,p.MB_lim,'Hybrid river incision',p.dispMB);
        
    end
    Qs_out_Riv(iter)=lvSedTot/p.TimeStep;
    DSY_out_Riv(iter)=tot_DSY/p.TimeStep;   

    
        Qs_hill=zeros(size(bedrock.Z));
        % Calcualte L_hill in case landslides occur bu no linear hillslope
        % erosion
        if p.LS_Bed>0
            switch p.FlowDirHill
                case 'single'
                    L_hill_i=dx_ik./(1-min((max(0,H1.Z(FD.ix)-H1.Z(FD.ixc)).*mod_Domain_M(FD.ix)./dx_ik)./Sc_Map.Z(FD.ix),0.99).^2);
                case {'multi','DInf'}
                    L_hill_i=dx_ik_Hill./(1-min((max(0,H1.Z(FD_Hill.ix)-H1.Z(FD_Hill.ixc)).*mod_Domain_M(FD_Hill.ix)./dx_ik_Hill)./Sc_Map.Z(FD_Hill.ix),0.99).^2);
                    multi_Lakes=Lakes.Z(FD_Hill.ix);
                    L_hill_i(multi_Lakes)=dx_ik_Hill(multi_Lakes);
            
            end
        end
        
    
    %% *************************** Landsliding ****************************
    %----------------------------Bedrock sliding---------------------------
    if p.LS_Bed
        if p.verbose
            disp('Bedrock Landsliding, cullman theory...')
            Progress.stage='Bedrock Landsliding';
            Progress.time=toc(cpu_runTime);
            writetable(Progress,[p.resultsdir  'Progress_Table_' p.fileprefix '.csv']);
        end
        
        [sediment, bedrock,Qs_bed_LS,SSY_LS_Bed,...
            LS_prop_Bed_store(iter),slidePlane_Bed,LSloc{iter},LS_Total_E] =...
            LS_cullman(...
            sediment, bedrock,nb.M,...
            locationsD,X,Y,Sc_Map,mod_Domain_M,slidePlane_Bed,FD,dx_ik,...
            dt,dx2,iter,p.maxLS_Size,p.phi,p.Ff_Hill,p.C_eff,...
            p.grav_earth,p.rho_rock,p.t_LS,...
            p.save_LS_Data,p.resultsdir, p.fileprefix,LS_Total_E,p.verbose_LS,p.saveeach);
    else
        LS_prop_Bed_store(iter)=0;
        Qs_bed_LS=0;
        SSY_LS_Bed=0;
    end
    
    %----------------------------Output landslides-------------------------
    Qs_LS=Qs_bed_LS;
    SSY=SSY_LS_Bed;
    Qs_out=Qs_hill+Qs_LS;
    SSY_m=SSY/dx2;
    
    MB_H_Dom=sum(dt*Qs_out((mod_Domain))/dx2);
    MB_H_OutDom=sum(dt*Qs_out((~mod_Domain))/dx2);
    MB_H=(MB_H_Dom+MB_H_OutDom)+SSY_m;
    
    
    if p.checkMB
        checkMB(iter,bed_Before,sed_Before,upl,dt,bedrock.Z, sediment.Z,...
            MB_R,MB_H,mod_Domain,p.phi,p.MB_lim,'Landsliding',p.dispMB);
    end
   
    %% ********************** Non-Linear deposition ***********************
    if ~pureIncision
        plotInterm=false;
        if plotInterm || p.save_LS_D
            SED1=sediment;
        end
        if p.LS_Bed~=0
            if p.verbose
                disp('Non-linear hillslope deposition...')
                Progress.stage='Non-linear hillslope deposition';
                Progress.time=toc(cpu_runTime);
                writetable(Progress,[p.resultsdir  'Progress_Table_' p.fileprefix '.csv']);
            end
            Qs_leavingDom=sum(Qs_out(mod_Domain==0));
            Qs_out(mod_Domain==0)=0;
            
            switch p.FlowDirHill
                case 'single'
                    [bedrock.Z, sediment.Z,Qs_Hill_OutD] = D_Hill(...
                        bedrock.Z, sediment.Z,L_hill_i,mod_Domain,...
                        FD.ix,FD.ixc,ones(size(FD.ix)),dx_ik,dt,dx,dx2,Qs_out,...
                        p.phi,p.min_HillGrad);
                case {'multi','DInf'}
                    
                     [bedrock.Z, sediment.Z,Qs_Hill_OutD] = D_Hill(...
                        bedrock.Z, sediment.Z,L_hill_i,mod_Domain,...
                        FD_Hill.ix,FD_Hill.ixc,FD_Hill.fraction,dx_ik_Hill,dt,dx,dx2,Qs_out,...
                        p.phi,p.min_HillGrad);

            end
            
            MB_H=(Qs_leavingDom+Qs_Hill_OutD)*dt/dx2+SSY_m;
            if p.checkMB
                checkMB(iter,bed_Before,sed_Before,upl,dt,bedrock.Z, sediment.Z,...
                    MB_R,MB_H,mod_Domain,p.phi,p.MB_lim,'Hillslope deposition',p.dispMB);
            end
            Qs_out_Hill(iter)=Qs_leavingDom/p.TimeStep;
            DSY_out_Hill(iter)=Qs_Hill_OutD/p.TimeStep;
        end
        
        if p.LS_Bed&&p.save_LS_D&& mod(iter,p.saveeach)==0
            LS_D=sediment-SED1;
            save([p.resultsdir p.fileprefix '_res_' num2str(round(bedrock.cellsize)) '_LS_D_' num2str(round(iter*dt)) '.mat'],'LS_D');
        end
        
        if plotInterm
            figure(11);imagesc(log10(LS_D));colorbar;
            [LS_x,LS_y]=ind2coord(H1,LSloc{iter});
            hold on
            sc=scatter(LS_x,LS_y,'filled');
            sc.MarkerFaceColor='r';
        end
     H1=bedrock+sediment;   
        
    end
    %% ********************* General Mass balance *************************
    if p.checkMB
        checkMB(iter,bed_Before,sed_Before,upl,dt,bedrock.Z, sediment.Z,...
            MB_R,MB_H,mod_Domain,p.phi,p.MB_lim,'full iteration',p.dispMB);
    end
    
    %% ************************* Boundary conditions *************************
    if strcmp(BC.type, 'OpenDir')||strcmp(BC.type, 'set_VaropenNodes')||strcmp(BC.type, 'Bed_Sed_Open')
        % Resolve the issue of having open BC conditions
        
        out=streampoi(FD,DA>0,'outlets','ix');
        % Get neighbours of open nodes
        nb_openBC=nb.M(out,:);
        nb_DA=double(nb_openBC);
        nb_DA(nb_DA~=0)=DA.Z(nb_DA(nb_DA~=0));
        nb_el=double(nb_openBC);
        nb_el(nb_el~=0)=double(H1.Z(nb_el(nb_el~=0)));
        [~,mix] = max(nb_DA,[],2);
        mix_C=sub2ind(size(nb_openBC),(1:numel(mix))', mix);
        BC.Giv_To_BC = nb_openBC(mix_C);
        BC.OpenNodes=out;
    end
    [bedrock, sediment]=setBC_HyLands(bedrock, sediment,BC);
    
    %% ************************* Update elevation *************************
    H1=bedrock+sediment;
    
    totalRelief(iter)=max(H1.Z(I_Domain.Z));   
    
    %% ************************* Save *************************
    if mod(iter,p.saveeach)==0
        save([resultsdir p.fileprefix '_res_' num2str(round(H1.cellsize)) 'm_bedrock_' num2str(round(t)) '.mat'],'bedrock');
        if p.K_sed~=0
            save([resultsdir p.fileprefix '_res_' num2str(round(H1.cellsize)) 'm_sediment_' num2str(round(t)) '.mat'],'sediment');
        end
    end
    
    %% ************************* Plot *************************
    if p.plotOut
        hGUIContinue=~get(hGUI.stopbutton,'value');
        if mod(iter,p.ploteach)==0
            H1_Domain=crop(H1,I_Domain);
            if p.plotSed
                Sed_Domain=crop(sediment,I_Domain);
                Bed_Domain=crop(bedrock,I_Domain);
                
                FD_p=FLOWobj(H1_Domain,'preprocess','carve','mex',true);
                S_Plot=STREAMobj(FD_p,'minarea',1e4/(FD.cellsize.^2));
                
                tspanstr = num2str(tspan);
                set(hGUI.timedisp,'String',[num2str(round(t)) ' / ' tspanstr ' years (' num2str(t/tspan*100,'%3.1f') '%)'] );
                
                figure(sedplot)
                sax1=subplot(1,4,1);
                cm1=color_DEMPoster;
                imagesc(H1_Domain)
                colorbar
                colormap(sax1,cm1)
                sax2=subplot(1,4,2);
                cm2=(color_heat);
                imagesc((Sed_Domain))
                colormap(sax2,cm2)
                colorbar
                caxis([0.9*min(Sed_Domain.Z(:)) 1.1*max(Sed_Domain.Z(:))]);
                hold off
                
                St=trunk((S_Plot));
                try
                    hold on; plot(St,'linewidth',2);hold off
                    sax2=subplot(1,4,3:4);
                    plotdz(St,DEM,'color',[0.6 0.6 0.6])
                    hold on
                    plotdz(St,Bed_Domain,'color','black')
                    plotdz(St,Bed_Domain+Sed_Domain,'color','red');
                    legend('Ini Topo','bedrock','Current Topo','location', 'northwest')
                    hold off
                    ylim([min(Bed_Domain.Z(:)) max(Bed_Domain.Z(:)+Sed_Domain.Z(:))]);
                catch
                    disp('Hold')
                end
                drawnow
            else
                tspanstr = num2str(tspan);
                set(hGUI.hIm,'Cdata',H1_Domain.Z);
                set(hGUI.timedisp,'String',[num2str(round(t)) ' / ' tspanstr ' years (' num2str(t/tspan*100,'%3.1f') '%)'] );
                clim  = get(hGUI.ax.main,'clim');
                tt    = get(hGUI.tsmin,'Xdata');
                ydmin = get(hGUI.tsmin,'Ydata');
                ydmax = get(hGUI.tsmax,'Ydata');
                ydmean= get(hGUI.tsmean,'Ydata');
                tt    = [tt(:); t];
                set(hGUI.tsmin,'XData',tt,'YData',[ydmin(:);clim(1)]);
                set(hGUI.tsmax,'XData',tt,'YData',[ydmax(:);clim(2)]);
                set(hGUI.tsmean,'XData',tt,'YData',[ydmean(:);mean(H1_Domain.Z(~isnan(H1_Domain.Z)))]);
                drawnow
            end            
        end
    end
end


%% ************************* Generate output *************************
HYLANDS_out.I_Domain=I_Domain;
HYLANDS_out.bedrock=bedrock;
HYLANDS_out.DA=DA;
HYLANDS_out.sediment=sediment;

HYLANDS_out.runTime=iter*dt;

if p.plotSed
    HYLANDS_out.S=S;
end
HYLANDS_out.DSY_out_Riv=DSY_out_Riv;
HYLANDS_out.DSY_out_Hill=DSY_out_Hill;
if p.LS_Bed~=0
    HYLANDS_out.LS_prop_bed_store=LS_prop_Bed_store;
    HYLANDS_out.slidePlane_Bed=slidePlane_Bed;
    HYLANDS_out.LS_Total_E=LS_Total_E;
end
end


