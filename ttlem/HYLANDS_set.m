function p = HYLANDS_set(varargin)
% set options and parameters for HyLands
% Syntax
%
%     p = HYLANDS_set(pn,pv,pn,pv,...)
%     p = HYLANDS_set(p)
% 
% Description
%
%     set parameters for HYLANDS
%
% Input arguments (parameter name-value pairs) {default}
%
% Computational parameters
%
%     'TimeSpan'     simulated years {2e6}
%     'TimeStep'     model timestep, dt {50}.     
%
% Output related parameters
%
%     'plotOut'     produce graphical output {true}
%     'ploteach'     {100}. Plot every nth time step (e.g. at times = 
%                    n*TimeStep). No plotting if n=inf.     
%     'plotSed'     Kind of plot: false for regular TTLEM output, true for
%                    sediment plotting {true}                       
%     'saveeach'     {inf}. Saves the modelled DEM as GRIDobj each nth time
%                    step.
%     'save_LS_D'   Boolean indicating to save or not to save
%                   deposition patterns produced by landslides during
%                   every iteration {false}-
%     'save_LS_Data'Boolean indicating to save or not to save
%                   erosion patterns produced by  landslides during
%                   every iteration {false}
%     'fileprefix'  Includes this prefix before each numbered.mat file when 
%                   saving the data to the harddrive{'HyLands_'}. 
%     'resultsdir'  output directory {['results' filesep]}
%     'verbose'     Print verbose {true}
%     'verbose_LS'  Print verbose landslides {true}
%
% Boundary conditions
%
%     'BC_Type'      General Boundary conditions. Options:
%                           * {'set_VaropenNodes'}: % Bedrock elevation at
%                               open nodes is set to a fixed value,
%                               provided by BC_BedDirVal.
%                               Sediment thickness varies through time as a
%                               function of the SPACE mathematics. Thus, at
%                               open nodes, the sediment thickness is based
%                               on the sediment thickness of the upstream
%                               river cell.
%                           * 'set_openNodes':
%                               fix both the elevation of the bedrok and
%                               the depth of the seiment layer at the
%                               outlet nodes to a predifined value (to the
%                               value of respectively BC_BedDirVal and
%                               BC_SedDirVal). For scenarios where river
%                               incision should be simulated as a hybrid
%                               process, preferable use set_VaropenNodes
%                           * 'Bed_Sed_Open', similar set_VaropenNodes; all
%                               nodes at the boundaries are considered as
%                               open nodes.
%                           * 'none'
%    'BC_nbGhost'    Number of Ghost cells during model run: in HyLands, we
%                       add one row and collumn of grid nodes at the
%                       boundaries which act as open nodes. {1}
%    'BC_BedDirVal'  Dirichlet fixed numeric value {0}. Can be either a
%                    single value or a (sparse) matrix containing the fixes 
%                    cells at the edges, respecting the number of ghost 
%                    cells. 
%    'BC_SedDirVal'  Dirichlet fixed numeric value {0}. Can be either a
%                    single value or a (sparse) matrix containing the fixes 
%                    cells at the edges, respecting the number of ghost 
%                    cells. 
%    'FlowBC'        Flow boundary conditions. At specified borders
%                    flow is forced to direct inward. The default is 'open'
%                    Setting 'b' forces the bottom boundary to flow inward.
%                    'rl' forces left and right boundaries to flow inward.
%                    FlowBC can take on any combination of 'b' (bottom),
%                    't' (top), 'r' (right) and 'l' (left). If a corner is
%                    defined ('ul_cor','ll_cor','ur_cor','lr_cor'), the
%                    corner will be set as the only open node over the
%                    modelled domain.Combinations should be entered as:
%                           * {'open'}
%                           * 't','b','l','r'
%                           * 'tb','tl','bl','br','lr','tr','bl','br','lr'
%                           * 'tbl','tbr'
%                           * 'tblr'
%                           * 'ul_cor','ll_cor','ur_cor','lr_cor'
%     'DrainDir'     {'variable'} 
%                    HyLands is only tested with variable drainage
%                    networks (that is: the drainage network being updated
%                    every iteration)
%     'FlowDir'     in this version of HyLands, only single flow directions
%                   for river incision are supported
%     'FlowDir'     Flow directions used tio redistribute hillslope
%                   sediment {'multi','single','Dinf'}:
%                           * {'multi'}
%                           * 'single'
%                           * 'Dinf'
% Mass balance
%     'checkMB'     Check mass balance during model run {false}
%     'MB_lim'      Absolute number indicating how much the mass balance
%                   might deviate from theoretical value  (units vary){10}
%     'dispMB'      Print mass balance during model run {false}
%
% Model parameters
%
% *General paramters*
%     'U_type'       Uplift pattern:
%                           * {'uniform'}
%                           * 'spatialVar': spatially variable uplift rate
%                           provided as GRIDobj in T 
%     'R_base'       constant runoff rate, in [m y^-1]  
%
% *River incision: SPACE*
%     'K_bed'       Erodibility for bedrock (units vary){5e5}
%     'K_sed'       Erodibility for sediment (units vary) {1e-4}
%     'm'           Drainage area exponent (units vary){.45}
%     'n'           Slope exponent (units vary) {1}
%     'phi'         Sediment porosity [-].{0}
%     'H_star'      Sediment thickness required for full entrainment [L].
%                   {1}
%     'omega_cr'    Critical stream power to erode rock [E/(TL^2)] {0}
%     'omega_cs'    Critical stream power to erode sediment [E/(TL^2)]{0}
%     'Ff'          Fraction of permanently suspendable fines in bedrock 
%                   [-]{0}
%     'V'           Effective settling velocity for chosen grain size 
%                   metric [L/T].{5}
%     'V_Lakes'     Effective settling velocity in flooded cells (lakes) 
%                   for chosen grain size metric [L/T].{5}
% *Landslides
%     'LS_Bed'      boolean indicating presence of bedrock landsliding 
%                   {false}
%     'Sc_fixed'	internal friction angle [%] {1}
%     'rho_rock'    density bedrock, [kg m^-3] {2700}
%     'rho_soil'    density sediment layer, [kg m^-3] {1350}
%     'rho_water'	density water, [kg m^-3] {1000}
%     't_LS'        landslide return time [m²] {500}
%     'maxLS_Size'  maximum landslide size [year] {1e9}
%     'C_eff'       effective cohesian [Pa] {1e4} 
%     'min_HillGrad'minimum sediment spreading slope [%] {0}
%     'Ff_Hill'     Fraction of permanently suspendable fines released 
%                   during landsliding [-]{0}
%
% *Constant paramters
%     's_year'      31557600 s
%     'grav_earth'	9.807 ms^2
%
%
% See also: HYLANDS
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
%
% =========================================================================
%
% Author:   Benjamin Campforts (benjamin.campforts@gfz-potsdam.de)
%
% Date:     15. March, 2020


%% Parse inputs
p = inputParser;
p.CaseSensitive = true;
p.FunctionName = 'HYLANDS_set';

%% Spatial and temporal domain
addParameter(p,'TimeSpan',2e6,@(x) isscalar(x) && x>0); 
addParameter(p,'TimeStep',50,@(x) isnan(x) || (isscalar(x) && x>0));

%% Output
addParameter(p,'plotOut',true,@(x) islogical(x));
addParameter(p,'ploteach',100,@(x) isscalar(x) && x>0);
addParameter(p,'plotSed',true,@(x) islogical(x));
addParameter(p,'saveeach',inf,@(x) isscalar(x) && x>0);
addParameter(p,'save_LS_D',false,@(x) islogical(x));
addParameter(p,'save_LS_Data',false,@(x) islogical(x));
addParameter(p,'fileprefix','HyLands_');
addParameter(p,'resultsdir',['HyLands_results' filesep]);
addParameter(p,'verbose',true,@(x) islogical(x));
addParameter(p,'verbose_LS',true,@(x) islogical(x));

%% Boundary conditions
addParameter(p,'BC_Type','set_VaropenNodes',@(x) ...
    ischar(validatestring(x,{'set_openNodes','set_VaropenNodes','Bed_Sed_Open','none'})));
addParameter(p,'FlowBC','open',@(x) ischar(validatestring(x,{'','t','b','l','r','tb','tl','tr','tbl','tbr','tblr','bl','br','lr',...
    'ul_cor','ll_cor','ur_cor','lr_cor','open'})));
addParameter(p,'BC_nbGhost',1,@(x) x==1);
addParameter(p,'BC_SedDirVal',0,@(x) ismatrix(x));
addParameter(p,'BC_BedDirVal',0,@(x) ismatrix(x));

%% drainage direction 
addParameter(p,'DrainDir','variable',@(x) ischar(validatestring(x,{'variable'})));
addParameter(p,'FlowDir','single',@(x) ischar(validatestring(x,{'single'})));
addParameter(p,'FlowDirHill','multi',@(x) ischar(validatestring(x,{'multi','single','Dinf'})));

%% check Mass balance
addParameter(p,'checkMB',false,@(x) islogical(x));
addParameter(p,'MB_lim',10,@(x) isscalar(x) && x>=0);
addParameter(p,'dispMB',false,@(x) islogical(x));

%% Uplift
addParameter(p,'U_type','uniform',@(x) ischar(validatestring(x,{'spatialVar','uniform'})));

%% Runoff
addParameter(p,'R_base',1,@(x) isscalar(x) && x>=0.01 && x<=10);

%% River incision: SPACE
addParameter(p,'K_bed',5e-5,@(x) ismatrix(x) && all(x(:)>=0));
addParameter(p,'K_sed',1e-4,@(x) ismatrix(x) && all(x(:)>=0));
addParameter(p,'m',.45,@(x) isscalar(x) && x>0);
addParameter(p,'n',1,@(x) isscalar(x) && x>0);
addParameter(p,'phi',0,@(x) isscalar(x) && x>=0 && x<=10);
addParameter(p,'H_star',1,@(x) isscalar(x) && x>=0 && x<=10);
addParameter(p,'omega_cr',0,@(x) isscalar(x) && x>=0 && x<=10);
addParameter(p,'omega_cs',0,@(x) isscalar(x) && x>=0 && x<=10);
addParameter(p,'Ff',0,@(x) isscalar(x) && x>=0 && x<=1);
addParameter(p,'V',5,@(x) isscalar(x) && x>=0 && x<=100);
addParameter(p,'V_Lakes',5,@(x) isscalar(x) && x>=0 && x<=100);

%% Landslides
addParameter(p,'Sc_fixed',1,@(x) (isscalar(x) && x>=0.1) || isempty(x));
addParameter(p,'rho_rock',2700,@(x) isscalar(x) && x>=1000);    
addParameter(p,'rho_soil',1350,@(x) isscalar(x) && x>=1000); 
addParameter(p,'rho_water',1000,@(x) isscalar(x) && x==1000);   
addParameter(p,'LS_Bed',false,@(x) islogical(x));
addParameter(p,'t_LS',500,@(x) isscalar(x) && x>=0);
addParameter(p,'maxLS_Size',1e7,@(x) isscalar(x) && x>=0);
addParameter(p,'C_eff',1e4,@(x) isscalar(x) && x>0); 
addParameter(p,'min_HillGrad',0,@(x) isscalar(x) && x>=0);
addParameter(p,'Ff_Hill',0,@(x) isscalar(x) && x>=0 && x<=1);

%% Constants
addParameter(p,'s_year',31557600,@(x) isscalar(x) && x==31557600);
addParameter(p,'grav_earth',9.807,@(x) isscalar(x) && x==9.807);

%% parse input
parse(p,varargin{:});
% and return
p = p.Results;

