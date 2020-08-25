function p = ttlemset(varargin)
% set options and parameters for TTLEM
%
% Syntax
%
%     p = ttlemset(pn,pv,pn,pv,...)
%     p = ttlemset(p)
%     p = ttlemset(p, pn,pv,pn,pv,...)
%
% Description
%
%     set parameters for ttlem
%
% Input arguments (parameter name-value pairs) {default}
%
% Computational parameters
%
%     'TimeSpan'     simulated years {200000}
%     'TimeStep_outer'     
%                    time step at which flow directions and diffusion are
%                    calculated {2000}.
%     'dt_inner_max'     
%                    maximum timestep used as an inner timestep for river incsion {3000}.
%                    When river incision is calculated with an implicit
%                    scheme, this only applies when the implCFL is set to
%                    true
%
%     'AreaThresh'   area threshold for channel initiation {0}. Set to
%                    values larger than zero (e.g. 1e6) to obtain a
%                    distinct river network.
%     'BC_Type'      General Boundary conditions. Options:
%                           * 'Dirichlet'
%                           * 'Neumann', currently only with 0 gradiient
%                           * 'Periodic'
%                           * 'Dirichlet_Matrix'
%                           * {'Dirichlet_Matrix_Ini'}: use intital
%                           values of border as fixed edges
%    'BC_nbGhost'    Number of Ghost cells for boundary calcuulations:
%                    {1} or 2
%    'BC_dir_value'  Dirichlet fixed numeric value {0}. Can be either a
%                    single value or a (sparse) matrix containing the fixes 
%                    cells at the edges, respecting the number of ghost 
%                    cells. 
%    'BC_dir_DistSites' 
%                    Random disturbance at the edges. Possible combinations:
%                           * {''}
%                           * 't','b','l','r'
%                           * 'tb','tl','bl','br','lr','tr'
%                           * 'tbl','tbr'
%                           * 'tblr'
%     'FlowBC'       flow boundary conditions. Forces specified borders
%                    inward directed flow. The default is '', meaning to
%                    flow boundary conditions. Setting 'b' forces the
%                    bottom boundary to flow inward. 'rl' forces left and
%                    right boundaries to flow inward. FlowBC can take on
%                    any combination of 'b' (bottom), 't' (top), 'r'
%                    (right) and 'l' (left). Combinations should be entered
%                    as:
%                           * {''}
%                           * 't','b','l','r'
%                           * 'tb','tl','bl','br','lr','tr'
%                           * 'tbl','tbr'
%                           * 'tblr'
%     'DrainDir'     {'variable'} or 'fixed' drainage directions. Setting
%                    'fixed' will decrease computation time but will keep
%                    the drainage system at the same locations.
%
% Model parameters
%
%     'D'            diffusivity {1} in [m^2 y^-1]
%     'Kw'           scaling of the incision equation dh/dt = Kw*A^m*S^n
%     'm'            area exponent
%     'n'            slope exponent
%     'm_var'        introduce variability adapted from Grimaldi et al. 2005 
%     'massWasting_river' {false}
%                    rivers cannot exceed a threshold slope if set to true. 
%                    Only needed in case rivers are not fixed. 
%     'NormPrecip'   effective rainfall grid (GRIDobj). Weights flow
%                    accumulation and should have values ranging between 0 
%                    and 1.
%     'K_weight'     effective Kw grid (GRIDobj). Weights Kw and should
%                    have values ranging between 0 and 1.  Eg. lithological
%                    strength
%     'rho_rs'       1.3. ratio between rho_r(rock) and rhos (soil)
%     'Sc'           0.8. Threshold angle (provided in unit Sc_unit
%                    (tangens by default) for hillslope adjustment to
%                    oversteepening
%     'Sc_unit'      unit of 'Sc'. {'tangent'},'degree', 'radian', 
%                    'percent', or 'sine'.
%
% Output related parameters
%
%     'ploteach'     {10}. Plot every nth time step (e.g. at times = 
%                    n*TimeStep). No plotting if n=inf.                   
%     'saveeach'     {inf}. Saves the modelled DEM as GRIDobj each nth time
%                    step.
%     'fileprefix'   {'res_'}. Includes this prefix before each numbered
%                    .mat file when saving the data to the harddrive.
%     'resultsdir'   output directory {['results' filesep]}
%
% Numerics
%
%     'diffScheme'   'imp_lin'. Linear diffusion solved by implicit method
%                    (see Pelletier 2008).
%                    'imp_lin_Sc'. Linear diffusion with critical slope
%                    angle.
%                    'imp_nonlin_sc'
%                    'only_Sc'
%     'diffToRiv'    true or {false}. Determines whether river cells
%                    experience diffusion. Setting to false makes only 
%                    sense if 'AreaThresh' is set to a value >0.
%     'riverInc'     {'implicit_FDM'}, implicit erosion scheme by Braun and
%                    Willett (2013)
%                    'explicit_FDM', explicit finite difference method. 
%                    'TVD_FVM', total variation diminishing finite volume
%                    method by Campforts and Govers (2015)
%     'cfls'         courant-friedrich-lewis value {.7} used for explicit
%                    scheme.
%     'parallel'     parallelization where possible (requires the Parallel
%                    Processing Toolbox, still under construction).
%     'shortening'   {false} tectonic shortening at the surface
%     'short_x'      {nan} horizontal shortening in x direction
%     'short_x'      {nan} horizontal shortening in y direction
%     'steadyState'  logical to check for steady state {false} 
%     'SS_value'     value untill SS is reached in summed meter over 
%                    modelled range {''}
%
% See also: ttlem
%
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
% Authors:   Benjamin Campforts (benjamin.campforts@kuleuven.be)
%           Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%
% Date:     8. July, 2016
%% Parse inputs
p = inputParser;
p.CaseSensitive = true;
p.FunctionName = 'ttlemset';

%% Spatial and temporal domain

% Time span in [years]
addParameter(p,'TimeSpan',2000000,@(x) isscalar(x) && x>0); 
% Time step for outer loop [years]
addParameter(p,'TimeStep',20000,@(x) isscalar(x) && x>0);

%% Boundary conditions
%Main model boundary conditions
addParameter(p,'BC_Type','Dirichlet_Matrix_Ini',@(x) ischar(validatestring(x,{'Dirichlet','Neumann','Periodic','Dirichlet_Matrix','Dirichlet_Matrix_Ini'})));
addParameter(p,'BC_nbGhost',1,@(x) x==1 || x==2);
addParameter(p,'BC_dir_value',0,@(x) ismatrix(x));
addParameter(p,'BC_dir_DistSites','',@(x) ischar(validatestring(x,{'','t','b','l','r','tb','tl','tr','tbl','tbr','tblr','bl','br','lr'})));%Edges where randomisation can occur
%Flow boundary conditions
addParameter(p,'FlowBC','',@(x) ischar(validatestring(x,{'','t','b','l','r','tb','tl','tr','tbl','tbr','tblr','bl','br','lr'})));
addParameter(p,'BC_dir_Dist_Value',1,@(x) isscalar(x) && x>0 && x<10); 


%% River incision
addParameter(p,'Kw',3e-6,@(x) isscalar(x) && x>=0);
addParameter(p,'m',.5,@(x) isscalar(x) && x>0);
addParameter(p,'n',1,@(x) isscalar(x) && x>0);
addParameter(p,'K_weight',[],@(x) isempty(x) || isa(x,'GRIDobj'));
addParameter(p,'m_var',0,@(x) isscalar(x) && x>=0 && x<=1);  
% drainage direction 
addParameter(p,'DrainDir','variable',@(x) ischar(validatestring(x,{'variable','fixed'})));
addParameter(p,'AreaThresh',0,@(x) isscalar(x) && x>=0);
addParameter(p,'NormPrecip',[],@(x) isempty(x) || isa(x,'GRIDobj'));
% Numerics river incision
addParameter(p,'riverInc','implicit_FDM',@(x) ischar(validatestring(x,{'implicit_FDM','explicit_FDM','TVD_FVM'})));
addParameter(p,'implCFL',false,@(x) islogical(x));
addParameter(p,'parallel',false,@(x) isscalar(x));


%% Hillslope processes
addParameter(p,'diffScheme','imp_lin_sc',@(x) ischar(validatestring(x,{'imp_lin','imp_lin_sc','imp_nonlin','imp_nonlin_sc','only_sc'})));
addParameter(p,'DiffToRiv',false,@(x) islogical(x));
addParameter(p,'D',.01,@(x) isscalar(x) && x>=0);

%% Hillslope adjustment to oversteepening
addParameter(p,'Sc',1.2,@(x) isscalar(x) && x>=0.1);
addParameter(p,'Sc_unit','tangent',@(x) ischar(validatestring(x,{'tangent','degree', 'radian', 'percent', 'sine'})));
addParameter(p,'rho_rs',1.3,@(x) isscalar(x) && x>=1);   % ratio between rhor(bedrock)and rhos (soil)

%% Tectonics
%Surface shortening
addParameter(p,'shortening',false,@(x) islogical(x));
addParameter(p,'shortening_meth','Upwind_TVD',@(x) ischar(x));
addParameter(p,'short_x',[],@(x) ismatrix(x)); 
addParameter(p,'short_y',[],@(x) ismatrix(x)); 

%% Numerical properties
addParameter(p,'cfls',.7,@(x) isscalar(x) && x>0 && x<=1);
addParameter(p,'DiffTol',1e-6,@(x) isscalar(x) && x>0); % DiffTol applied in nonlinDiff during pcg 

%% Run to steady state
addParameter(p,'steadyState',false,@(x) islogical(x));
addParameter(p,'SS_Value',0,@(x) isscalar(x) && x>=0);

%% Output
% plot each xth iteration
addParameter(p,'ploteach',1,@(x) isscalar(x) && x>0);
% save DEM at times...
addParameter(p,'saveeach',inf,@(x) isscalar(x) && x>0);
addParameter(p,'fileprefix','res_');
addParameter(p,'resultsdir',['results' filesep]);

%% parse input
parse(p,varargin{:});
% and return
p = p.Results;




