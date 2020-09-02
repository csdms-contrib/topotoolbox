# TTLEM and HyLands

<img src= "https://github.com/BCampforts/topotoolbox/blob/HyLands/ttlem/Fig_5%20Syntehtic%20DEM.jpg" align=" center" >

## What is TTLEM
TTLEM is a grid based object oriented Landscape Evolution Model based on the TopoToolbox library. The aim of TTLEM is to provide an accurate numerical tool to simulate landscape evolution in eroding regions. TTLEM is multi-process based and hence allows a modular use of the code depending on the requirements of the end user. The LEM has a good overall computational performance and obtains a high numerical accuracy in solving the partial differential equations involved when simulating landscape evolution.

## Getting started
TTLEM and HyLands are a part of the newest release of TopoToolbox and can be downloaded from the GitHub repository (https://github.com/wschwanghart/topotoolbox). Extract the package to a folder on your personal computer and add the TopoToolbox folder with its subfolders to the search path in the Matlab environment.

Before working with TopoToolbox the functions need to be on the search 
path of Matlab. Enter following code into the command line:
        
        addpath(genpath(['C:\path\to\wherever\you\installed\this\topotoolbox-v2.4-HyLands-v1.0']));

## Who should use TTLEM?
TTLEM is developed to be easy accessible and adaptable and requires knowledge of the Matlab syntax. Users should be able to readily run TTLEM with both synthetic land evolution models and natural simulations, starting from an existing digital elevation model. With TTLEM, we target the same users than those who are using TopoToolpox and want to explore landscape functioning beyond present day topographic indices.

## Where can I report bugs or suggestions for enhancement?
In case of bugs in the source code, any suggestions for model improvements or elaboration, please contact us on one of these mail addresses:

Benjamin Campforts:     benjamin.campforts[at]gfz;potsdam.de

Wolfgang Schwanghart:   w.schwanghart[at]geo.uni-potsdam.de

## How can I stay informed about updates?
Small updates will be posted straight on the GitHub repository. Major software updates will always be announced at this blog: topotoolbox.wordpress.com.

## User manuals TTLEM
TTLEM comes with 3 tutorials. All of them can be executed in Matlab. E.g., for the first example just enter the comment: 'TTLEM_usersguide_1_intro' in the command window.
- TTLEM_usersguide_1_intro.m Short introduction to the use of TTLEM
- TTLEM_usersguide_2_Synthetic_model_run.m In this tutorial, we show how TTLEM can be used to simulate synthetic landscape evolution. We illustrate how the user can change between different algorithms to simulate hillslope response and how different numerical schemes can be set. Finally, we also show how modelled data can be processed and converted in a movie of an evolving landscape. Detailed information of parameters, their default values and units can be found in the help section of ttlemset.
- TTLEM_usersguide_3_Synthetic_Geological_Configuration One of the main advantages of LEMs is their potential to explore different geo-tectonic configurations and corresponding hypothesis at a  very low cost. As an example, we simulated landscape evolution over a time span of 30 Myr. Successively, three spatial and temporally dependent uplift scenarios are imposed to the model. On top of this three scenarios, a horizontal shortening field is imposed over the simulated domain with high lateral displacements in the left upper corner, linearly decreasing in both x and y directions to the Southeastern bottom of the simulated domain.

## User manual HyLands

<img src= "https://github.com/BCampforts/topotoolbox/blob/HyLands/ttlem/HyLands.jpg" align=" center" >

Installation: 
 - Download and extract topotoolbox from: https://github.com/wschwanghart/topotoolbox (the HyLands branch, which is the default branch on the linked repository so nothing has to be adjusted, downloaded fodler will be named: topotoolbox-HyLands by default)
 - Before working with HyLands the directories and functions must be on the search path of Matlab. Enter following code into the command line: addpath(genpath(['C:\path\to\wherever\you\extracted\this\topotoolbox-HyLands']));
 - to verify installation: enter: doc HYLANDS (info on the model) or doc HYLANDS_set (info on the parameter values) in the command window.
 
Documentation: 
 - For documentation, and after adding the topotoolbox folder to the path, enter: doc HYLANDS (info on the model) or doc HYLANDS_set (info on the parameter values) in the command window. 
        

HyLands comes with 1 tutorial (the example provided on the documentation page: doc HYLANDS), and 7 sripts, described in full length in the GMD discussion paper: Campforts B., Shobe M.C., et al. : HyLands 1.0: a Hybrid Landscape evolution model to simulate the impact of landslides and landslide-derived sediment on landscape evolution. Discussion paper in Geoscientific Model Development, https://geoscientific-model-development.net

All the scripts can be downloaded from https://github.com/BCampforts/pub_hylands_campforts_etal_GMD and executed in Matlab. 
- HyLands_NoLS_DL.m;    https://doi.org/10.5446/45969
- HyLands_NoLS_TL.m;    https://doi.org/10.5446/45967
- HyLands_NoLS_Mixed.m; https://doi.org/10.5446/45968
- HyLands_LS_NB.m;      https://doi.org/10.5446/45973
- HyLands_LS_B_LS.m;    https://doi.org/10.5446/45970
- HyLands_LS_LS.m;      https://doi.org/10.5446/45971
- HyLands_LS_A_LS.m;    https://doi.org/10.5446/45972


## Requirements

TTLEM/HyLands is plat-form independent and requires Matlab 2018a or higher 

## References

When you use HyLAnds/TTLEM, please reference to this publication:

 
  HyLands: Campforts B., Shobe M.C., et al. : HyLands 1.0: a Hybrid
  Landscape evolution model to simulate the impact of landslides and
  landslide-derived sediment on landscape evolution. Discussion paper in
  Geoscientific Model Development  
 
  TTLEM: Campforts, B., Schwanghart, W., & Govers, G. (2017). 
  Accurate simulation of transient landscape evolution
  by eliminating numerical diffusion: the TTLEM 1.0 model. Earth Surface
  Dynamics, 5(1), 47–66. https://doi.org/10.5194/esurf-5-47-2017

Other relevant publications: 

Campforts, B. and Govers, G. (2015): Keeping the edge: A numerical method that avoids knickpoint smearing when solving the stream power law, J. Geophys. Res. Earth Surf., 120(7), 1189–1205, [doi:10.1002/2014JF003376] 
https://www.researchgate.net/publication/279957445_Campforts_and_Govers-2015-JGR_ES_Keeping_the_edge_SI

Schwanghart, W. & Scherler, D. (2014): TopoToolbox 2 – MATLAB-based 
software for topographic analysis and modeling in Earth surface sciences. 
Earth Surface Dynamics, 2, 1-7. [DOI: 10.5194/esurf-2-1-2014]
https://www.researchgate.net/publication/259706134_Short_Communication_TopoToolbox_2_-_MATLAB-based_software_for_topographic_analysis_and_modeling_in_Earth_surface_sciences
  
Schwanghart, W., Kuhn, N. J. (2010): TopoToolbox: a set of Matlab 
functions for topographic analysis. Environmental Modelling & Software, 
25, 770-781. [DOI: 10.1016/j.envsoft.2009.12.002]
https://www.researchgate.net/publication/223084856_TopoToolbox_A_set_of_MATLAB_functions_for_topographic_analysis

## Version History

HyLands V1.0 --- 15. March 2020 
- HyLands added to ttlem

V1.0 --- 25. June 2016 
- release
