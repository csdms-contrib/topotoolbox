%% User Guide to TTLEM 3 - Synthetic model run with complex geological configuration
%% Content
% One of the main advantages of LEMs is their potential to explore
% different geo-tectonic configurations and corresponding hypothesis at a
% very low cost. As an example, we simulated landscape evolution over a
% time span of 30 Myr. Successively, three spatial and temporally dependent
% uplift scenarios are imposed to the model. On top of this three
% scenarios, a horizontal shortening field is imposed over the simulated
% domain with high lateral displacements in the left upper corner, linearly
% decreasing in both x and y directions to the Southeastern bottom of the
% simulated domain. 
%
% See also:     ttlemset, ttlem
%
% * TTLEM:
% Campforts B., Schwanghart W., Govers G.: TTLEM 1.0 : a numerical package
%           for accurate simulation of tansient landscape evolution in
%           MATLAB. Discussion paper in GMD.
%
% * TopoToolbox: Schwanghart, W. and Scherler, D.: Short Communication:
%           TopoToolbox 2 – MATLAB-based software for
%           topographic analysis and modeling in Earth surface sciences,
%           Earth Surf. Dyn., 2(1), 1–7,doi:10.5194/esurf-2-1-2014, 2014.
%           <https://www.researchgate.net/publication/259706134_Short_Communication_TopoToolbox_2_-_MATLAB-based_software_for_topographic_analysis_and_modeling_in_Earth_surface_sciences>
% 
% * TVD-FVM: Campforts, B. and Govers, G.:
%           Keeping the edge: A numerical method that avoids knickpoint
%           smearing when solving the stream power law, J. Geophys. Res.
%           Earth Surf., 120(7), 1189–1205, doi:10.1002/2014JF003376, 2015.
%           <https://www.researchgate.net/publication/279957445_Campforts_and_Govers-2015-JGR_ES_Keeping_the_edge_SI>
%
% Authors:   Benjamin Campforts (benjamin.campforts[at]kuleuven.be)
%            Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%
% Date:     8. July, 2016

%% Clear environment
clearvars
clc
close all


%% Initial Surface
% Generate random initial surface of 0m ± 50m

dx=400;%m
Lx=200e3;
Ly=100e3;
x=dx:dx:Lx;
y=dx:dx:Ly;
Z=zeros(numel(y),numel(x));
Z(2:end-1,2:end-1) =rand(numel(y)-2,numel(x)-2)*50;
H1=GRIDobj(x,y,Z);
figure('units','normalized','outerposition',[0.1 0.1 .5 .5],'color','white');
imageschs(H1);


%% Temporal domain
p.TimeSpan=30e6;
p.TimeStep=1e5;

%% Vertical uplift
% Vertical uplift rates are considered here as a sequence of three
% tectonic configurations. 
% Vertical uplift is inserted in TTLEM as an instance of GRIDobj
U = GRIDobj(H1);

%% First uplift phase
% In a first phase, the spatial domain is separated by a tectonic fault. In
% the eastern part, the surface experience a surface uplift of 3 mm/y
% whereas in the western part, uplift equals 1 mm/y. 

fault_x=round(size(U.Z,2)/2);
spat_u_1=ones(size((U.Z)))*3e-3;
spat_u_1(:,fault_x+1:end)=1e-3;
spat_u_1(1,:)=0;spat_u_1(end,:)=0;spat_u_1(:,1)=0;spat_u_1(:,end)=0;

% This first pattern is simulated during 10 Myr. The temporal dimention of
% the inserted uplift is provided in the 3th dimension of the uplift
% matrix. 
nbOfSteps1=10e6/p.TimeStep;
U.Z(2:end-1,2:end-1,1:nbOfSteps1) = repmat(spat_u_1(2:end-1,2:end-1),[1,1,nbOfSteps1]);

%Display first spatial uplift pattern
figure('units','normalized','outerposition',[0.1 0.1 .5 .5],'color','white');
imagesc(spat_u_1)
colorbar
title('First vertical uplift pattern (0 - 10 Myr)')


%% Second uplift phase
% A constant uplift over 10 Myr over the entire spatial domain at 1 mm/y
spat_u_2=ones(size(U.Z,1),size(U.Z,2));
spat_u_2=spat_u_2*1e-3;
spat_u_2(1,:)=0;spat_u_2(end,:)=0;spat_u_2(:,1)=0;spat_u_2(:,end)=0;
nbOfSteps2=10e6/p.TimeStep;
U.Z(2:end-1,2:end-1,nbOfSteps1+1:nbOfSteps1+nbOfSteps2) = ...
    repmat(spat_u_2(2:end-1,2:end-1),[1,1,nbOfSteps2]);

%Display first spatial uplift pattern
figure('units','normalized','outerposition',[0.1 0.1 .5 .5],'color','white');
imagesc(spat_u_2)
colorbar
title('Second vertical uplift pattern (10 - 20 Myr)')


%% Third uplift phase 
% An  uplift gradient of 10 Myr over the entire spatial domain at 1-5 mm/y
spat_u_3=repmat(linspace(1,0.5,size(U.Z,1))',1,size(U.Z,2))*5*1e-3;
spat_u_3(1,:)=0;spat_u_3(end,:)=0;spat_u_3(:,1)=0;spat_u_3(:,end)=0;
nbOfSteps3=10e6/p.TimeStep;
U.Z(2:end-1,2:end-1,nbOfSteps1+nbOfSteps2+1:nbOfSteps1+nbOfSteps2+nbOfSteps3) = ...
    repmat(spat_u_3(2:end-1,2:end-1),[1,1,nbOfSteps3]);

%Display first spatial uplift pattern
figure('units','normalized','outerposition',[0.1 0.1 .5 .5],'color','white');
imagesc(spat_u_3)
colorbar
title('Third vertical uplift pattern (20 - 30 Myr)')

%% Lateral tectonic displacement
% We insert a tectonic shortening  field operational in two directions and
% temporally constant over the entire model simulation
p.shortening=true;
% p.shortening_meth='Upwind_FD';
p.shortening_meth='Upwind_TVD';

x_Speed_c=3e-3;%m/y
y_Speed_c=3e-3;%m/y

%Gradient
grad_x=repmat(linspace(1,0,size(U.Z,2)),size(U.Z,1),1);
grad_y=repmat(linspace(1,0,size(U.Z,1))',1,size(U.Z,2));
p.short_x=grad_x*x_Speed_c;%m/y
p.short_y=grad_y*y_Speed_c;

figure('units','normalized','outerposition',[0.1 0.1 .5 .5],'color','white');
resQv=25;
[xq, yq] =meshgrid(resQv:resQv:size(p.short_x,2),resQv:resQv:size(p.short_x,1));
xq=xq.*dx; yq=yq.*dx;
hq =quiver(xq, yq,p.short_x(resQv:resQv:size(p.short_x,1),resQv:resQv:size(p.short_x,2)),p.short_y(resQv:resQv:size(p.short_y,1),resQv:resQv:size(p.short_y,2)));
hq.Color='k';
set(gca,'ydir','reverse');
xlim([min(xq(:)) max(xq(:))])
ylim([min(yq(:)) max(yq(:))])
title('Lateral tectonic displacement');

%% Diffusion parameters
p.D=0.02;
% p.diffScheme = 'only_sc';
p.diffScheme = 'imp_lin_sc';
% p.diffScheme = 'imp_nonlin_sc';
p.DiffToRiv=false; %If false, river elevations will only be altered by river incision
p.DiffTol=1e-4;

%% Boundary condition
p.BC_nbGhost=2;
p.BC_Type='Dirichlet';
p.BC_dir_value=100;
p.BC_dir_DistSites='tl';
p.BC_dir_Dist_Value=.5;

%% River incision parameters
p.Kw =7.5e-6;
p.m=0.5;
p.n=1;
p.AreaThresh=5e5; % channel contributing area threshold [m^2]
p.NormPrecip=[]; %effective rainfall grid

%% Numerics river incision
% p.riverInc = 'TVD_FVM';
p.riverInc = 'implicit_FDM';
% p.riverInc = 'explicit_FDM';

%% Boundary conditions
% Default

%% Threshold slopes
p.Sc=.7;
p.Sc_unit='tangent';

%% Output
p.ploteach=1;
p.saveeach=1;
p.fileprefix = 'GeologicalSimulation';

%% Initialize parameter structure
% By making p an instance of ttlemset, the user ensures parameter values
% are set in the right way
p   = ttlemset(p);

%% Model run
% TTLEM can be manually interrupted by pushing the 'Stop' Bottom. The
% current model run will quit without losing the information calculated so
% far. 
ttlem_out = ttlem(H1,U,p);

%% Movie and overview figure
% Create simulation movie of the results and generate an overview figure
% illustration of the results

overviewFignbs=floor(linspace(1,nbOfSteps1+nbOfSteps2+nbOfSteps3,9));
scrsz = get(0,'ScreenSize');
movieTopo = figure('OuterPosition',[0.2*scrsz(4) 0.2*scrsz(4) scrsz(4) .5*scrsz(4)],'Name','Topo','color','white');
overv=0;
Z_prev=H1.Z;
ui=0;
OverviewT = figure('units','normalized','outerposition',[0.1 0.1 .9 .9],'color','white');

for t=p.TimeStep:p.TimeStep:p.TimeSpan
    ui=ui+1;
    figure(movieTopo)
    load([p.resultsdir p.fileprefix num2str(t) '.mat'],'H1');
    FD  = FLOWobj(H1,'preprocess','c');
    S = STREAMobj(FD,flowacc(FD)>(p.AreaThresh/FD.cellsize^2));    
    imageschs(H1,[],'ticksToKm',true,'colorBarLabel','\bfm');      
    hold on
    S.x=S.x*1e-3;
    S.y=S.y*1e-3;
    plot(S,'k','linewidth',1.5);
    hold off
    xlabel('\bfX Coordinate, km')
    ylabel('\bfY Coordinate, km')
    title([num2str(round(t*1e-6)) ' Myr']);
    movie1(ui)=getframe(gcf); %#ok<SAGROW>
    if any(ui==overviewFignbs)       
        figure(OverviewT)
        overv=overv+1;
        subplot(3,3,overv)        
        imageschs(H1,[],'ticksToKm',true,'colorBarLabel','\bfm');        
        xlabel('\bfX Coordinate, km')
        ylabel('\bfY Coordinate,k m')
        title([num2str(round(t*1e-6)) ' Myr']);
        drawnow        
    end
end

%% History
%
% This user guide was updated last: June 12, 2016.