%% User Guide to TTLEM 2 - Synthetic model run
%% Content
% In this tutorial, we show how TTLEM can be used to simulate synthetic
% landscape evolution. We illustrate how the user can change between
% different algorithms to simulate hillslope response and how different
% numerical schemes can be set. Finally, we also show how modelled data can
% be processed and converted in a movie of an evolving landscape. Detailed
% information of parameters, their default values and units can be found in
% the help section of ttlemset.
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
close all
clc

%% Initial Surface
% Generate random initial surface of 0m ± 50m
dx=75;%m
Lx=50e3; %m
Ly=50e3;
x=dx:dx:Lx;
y=dx:dx:Ly;
Z=rand(numel(y),numel(x))*50; %Randomized initial condition
DEM=GRIDobj(x,y,Z);

% Display initial DEM
figure
imageschs(DEM)
[~,X_DEM,Y_DEM] = GRIDobj2mat(DEM);
[X_DEM,Y_DEM]   = meshgrid(X_DEM,Y_DEM);

%% Uplift
% Uplift is inserted in TTLEM as an instance of GRIDobj
U    = GRIDobj(DEM);
grad = ones(DEM.size);
upl_rate=1e-3; %Aboslute vertical uplift rate in m/yr
U.Z(2:end-1,2:end-1) = grad(2:end-1,2:end-1)*upl_rate;

%% Temporal domain
p.TimeSpan=20e6;
p.TimeStep=5e4;

%% Hillslope processes, parameters values
%Set diffusivity parameter
p.D=0.03;

%Choose hillslope response scheme
% p.diffScheme = 'imp_lin';
p.diffScheme = 'imp_lin_sc';
% p.diffScheme = 'imp_nonlin_sc';
% p.diffScheme = 'only_sc';


%% Set numerical accuracy
p.DiffTol = 1e-4;

%% River incision parameters
p.Kw = 3e-6;
p.m  = 0.5;
p.n  = 1;
p.AreaThresh = 2e5; % channel contributing area threshold, m²

% (un)comment to get the drainage development type of choice
% Fixed or variable drainage netwrok through time
p.DrainDir='variable';
% p.DrainDir='fixed';

% Insert variability on m value, following Grimaldi et al. (2005)
p.m_var=.3;

% Insert variability on K value. This is not compatible with the nonlinear
% Q-imp diffusion algorithm. Comment the next lines out in the case that
% nonlinear hillslope diffusion is required.
K_weight   = randn(size(DEM.Z));
K_weight   = (K_weight-min(K_weight(:)))./(max(K_weight(:))-min(K_weight(:)));
p.K_weight = GRIDobj(DEM); p.K_weight.Z=K_weight; 
% Correct K so that average K remains equal 
avgK       = mean(K_weight(:)); 
p.Kw       =p.Kw*1/avgK;

%% Numerics river incision
% (un)comment to get the numerical scheme of choice
% p.riverInc = 'TVD_FVM';
p.riverInc = 'implicit_FDM';
% p.riverInc = 'explicit_FDM';
p.cfls=0.95;

%% Threshold slopes
p.Sc=1; %%21°)
p.Sc_unit='tangent';

%% Boundary conditions
% Set to default.

%% Steady state
% Evaluate whether the model evolves towards a dynamic equilibrium where
% uplift balances erosion and absolute elevation differences between model
% time steps are minimal.
p.steadyState  = true;
p.SS_Value     = 1e-6*size(DEM.Z,1)*size(DEM.Z,2);%Max 1mm elevation change per cell allowed in steady state

%% Output
p.ploteach=1;
p.saveeach=1;
% Specify the location where the results can be stored (e.g.
% p.resultsdir='C:\...\'); The default folder is the result file where the
% main model structure is stored.
% p.resultsdir='C:\...';
p.fileprefix='syntheticRun';

%% Initialize parameter structure.
% By making p an instance of ttlemset, the user ensures parameter values
% are set in the right way
p   = ttlemset(p);

%% Model run
% TTLEM can be manually interrupted by pushing the 'Stop' Bottom. The
% current model run will quit without losing the information calculated so
% far. 
ttlem_out = ttlem(DEM,U,p);

%% Show model output
figure('units','normalized','outerposition',[0 0 1 1])
H1=ttlem_out.H1;
imageschs(H1,[],'ticksToKm',true,'colorBarLabel','\bfm');
xlabel('\bfX Coordinate, km');
ylabel('\bfY Coordinate, km');


%% Create simulation movie
scrsz = get(0,'ScreenSize');
hFig  = figure('OuterPosition',[0.1*scrsz(4) 0.1*scrsz(4) .7*scrsz(4) .7*scrsz(4)]);

if strcmp(p.diffScheme,'imp_nonlin_sc')
    parameter=1e-10;
    inc=max(max(U.Z(:)),1e-3);
    dt_max=DEM.cellsize*DEM.cellsize/inc*max(p.Sc(:))*parameter/(p.D*p.Kw);
    dt_max=min(dt_max,5000);
    nbSteps=ceil(p.TimeSpan/dt_max);
    if (p.TimeSpan/nbSteps)<p.TimeStep
        p.TimeStep=p.TimeSpan/nbSteps;
        p.TimeStep       = p.TimeStep;
        disp(['Timestep is decreased to ' num2str(round(p.TimeStep))...
            ' yr in order to keep the non linear solution stable, see Perron 2011, JGR']);
    end
end

figure
files = dir([p.resultsdir p.fileprefix '*.mat']);
[~,ix] = sort([files.datenum]);
files = files(ix);

for r = 1:numel(files);
    t = p.TimeStep*r;
    
    load([p.resultsdir files(r).name],'H1');
    imageschs(H1,[],'tickstokm',true,'colorBarLabel','\bfm');
    hold on
    FD  = FLOWobj(H1,'preprocess','c');
    S   = STREAMobj(FD,flowacc(FD)>(p.AreaThresh/FD.cellsize^2));
    S.x = S.x*1e-3;
    S.y = S.y*1e-3;
    plot(S,'k','linewidth',1)
    hold off
    xlabel('\bfX Coordinate, km')
    ylabel('\bfY Coordinate, km')
    title(['Standard Run ' num2str(round(t*1e-3)) 'k years']);
    drawnow
    
end

%% History
%
% This user guide was updated last: June 12, 2016.