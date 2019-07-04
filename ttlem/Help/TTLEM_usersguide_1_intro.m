%% User Guide to TTLEM - Introduction
%% What is TTLEM?
% TTLEM is a grid based object oriented Landscape Evolution Model based on
% the TopoToolbox library. The aim of TTLEM is to provide an accurate
% numerical tool to simulate landscape evolution in eroding regions. TTLEM
% is multi-process based and hence allows a modular use of the code
% depending on the requirements of the end user.  The LEM has a good
% overall computational performance and obtains a high numerical accuracy
% in solving the partial differential equations involved when simulating
% landscape evolution.
%
%% Who should use TTLEM?
% TTLEM is developed to be easy accessible and adaptable and requires
% knowledge of the Matlab syntax. Users should be able to readily run TTLEM
% with both synthetic land evolution models and natural simulations,
% starting from an existing digital elevation model. With TTLEM, we target
% the same users than those who are using TopoToolpox and want to explore
% landscape functioning beyond present day topographic indices.
%% Installation
% TTLEM comes as a part of the newest release of TopoToolbox and can be
% downloaded from the GitHub repository. Extract the package to a folder on
% your personal computer and add the TopoToolbox folder with its subfolders
% to the search path in the Matlab environment.
%% Where can I find a documentation?
% TTLEM comes with a basic manual with a couple of examples on how to
% initialize model runs, change input source files and adapt model
% parameters. Nearly all functions are documented with a help section where
% the syntax of the function and the background of the used algorithm is
% described. We intentionally kept the user manual short and focus on the
% documentation of the source code in the models functions to allow for
% streamlined updates concerning both the algorithms and the documentation
% that comes with it. Further, the functionality and versatility of the
% model are presented in this discussion paper about TTLEM, published in
% GMD. All scripts used to produce the results and figures shown in this
% paper can be found here.
%
%% Where can I report bugs or suggestions for enhancement?
% In case of bugs in the source code, any suggestions for model
% improvements or elaboration, you can contact us on the further indicated
% mail adresses. 
%% How can I stay informed about updates?
% Small updates will be posted straight on the GitHub repository. Major
% software updates will always be announced at this blog:
% topotoolbox.wordpress.com.
%% Content
% This manual consists out of 3 tutorials. All of them can be executed in
% Matlab. E.g., for the first example just enter the comment:
% 'TTLEM_usersguide_1_intro' in the command window.
%
% * TTLEM_usersguide_1_intro.m Short introduction to the use of TTLEM
%
% * TTLEM_usersguide_2_Synthetic_model_run.m In this tutorial, we show how
% TTLEM can be used to simulate synthetic landscape evolution. We
% illustrate how the user can change between different algorithms to
% simulate hillslope response and how different numerical schemes can be
% set. Finally, we also show how modelled data can be processed and
% converted in a movie of an evolving landscape. Detailed information of
% parameters, their default values and units can be found in the help
% section of ttlemset.
% 
% * TTLEM_usersguide_3_Synthetic_Geological_Configuration One of the main
% advantages of LEMs is their potential to explore different geo-tectonic
% configurations and corresponding hypothesis at a  very low cost. As an
% example, we simulated landscape evolution over a time span of 30 Myr.
% Successively, three spatial and temporally dependent uplift scenarios are
% imposed to the model. On top of this three scenarios, a horizontal
% shortening field is imposed over the simulated domain with high lateral
% displacements in the left upper corner, linearly decreasing in both x and
% y directions to the Southeastern bottom of the simulated domain.
%
%% A note on TopoToolbox
% TopoToolbox provides a set of Matlab functions that support the analysis
% of relief and flow pathways in digital elevation models. The major aim of
% TopoToolbox is to offer helpful analytical GIS utilities in a non-GIS
% environment in order to support the simultaneous application of
% GIS-specific and other quantitative methods. TopoToolbox is written in
% the Matlab language and requires the Image Processing Toolbox for various
% functions.
% We refer to the user guides of TopoToolbox for details on this software
% package. The user guide for TT can be opened by executing this comment: 
% usersguide_1_intro
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


%% Clean workspace 
clc
close all force
clearvars
%% Load a DEM into Matlab
%
% TopoToolbox 2 reads the ESRI ascii grid format and single band geotiffs
% into an instance of GRIDobj. 
% *Note that, throughout the use of TopoToolbox and TTLEM, it is assumed 
% that the DEM has a projected coordinate system (e.g. UTM WGS84) and that
% elevation and horizontal coordinates are in meter units.*
%
DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
imageschs(DEM,[],'ticksToKm',true,'colorBarLabel','\bfm');
xlabel('\bfX Coordinate, km');
ylabel('\bfY Coordinate, km');

%% Set uplift
%
% Uplift is provided to TTLEM as an instance of GRIDobj.
% Here, we illustrate the example of uniform uplift in space and time of 1
% mm per year. 
U = GRIDobj(DEM);
U.Z(2:end-1,2:end-1) = 1e-3;

%% Set model parameter values
% Parameters for the model run are set as an instance of the class
% ttlemset. If no arguments are provided, default parameter values are
% used. 

p = ttlemset %#ok show parameters<NOPTS>

%% Run model 
% TTLEM can be manually interrupted by pushing the 'Stop' Bottom. The
% current model run will quit without losing the information calculated so
% far. 

output = ttlem(DEM,U,p);

%% View results
figure('units','normalized','outerposition',[0 0 1 1])
H1=output.H1;
imageschs(H1,[],'ticksToKm',true,'colorBarLabel','\bfm');
xlabel('\bfX Coordinate, km');
ylabel('\bfY Coordinate, km');

%% History
%
% This user guide was updated last: June 12, 2016.













