# TTLEM - TopoToolbox Landscape Evolution Model

<img src= "https://github.com/BCampforts/topotoolbox/blob/master/ttlem/Fig_5%20Syntehtic%20DEM.jpg" align=" center" >

## What is TTLEM
TTLEM is a grid based object oriented Landscape Evolution Model based on the TopoToolbox library. The aim of TTLEM is to provide an accurate numerical tool to simulate landscape evolution in eroding regions. TTLEM is multi-process based and hence allows a modular use of the code depending on the requirements of the end user. The LEM has a good overall computational performance and obtains a high numerical accuracy in solving the partial differential equations involved when simulating landscape evolution.

## Installation
TTLEM comes as a part of the newest release of TopoToolbox and can be downloaded from the GitHub repository. Extract the package to a folder on your personal computer and add the TopoToolbox folder with its subfolders to the search path in the Matlab environment.
For further instructions:
https://github.com/wschwanghart/topotoolbox/blob/master/readme.md

## Who should use TTLEM?
TTLEM is developed to be easy accessible and adaptable and requires knowledge of the Matlab syntax. Users should be able to readily run TTLEM with both synthetic land evolution models and natural simulations, starting from an existing digital elevation model. With TTLEM, we target the same users than those who are using TopoToolpox and want to explore landscape functioning beyond present day topographic indices.

## Where can I report bugs or suggestions for enhancement?
In case of bugs in the source code, any suggestions for model improvements or elaboration, please contact us on one of these mail addresses:

Benjamin Campforts:     benjamin.campforts[at]kuleuven.be

Wolfgang Schwanghart:   w.schwanghart[at]geo.uni-potsdam.de

## How can I stay informed about updates?
Small updates will be posted straight on the GitHub repository. Major software updates will always be announced at this blog: topotoolbox.wordpress.com.

## User manuals
TTLEM comes with 3 tutorials. All of them can be executed in Matlab. E.g., for the first example just enter the comment: 'TTLEM_usersguide_1_intro' in the command window.
- TTLEM_usersguide_1_intro.m Short introduction to the use of TTLEM
- TTLEM_usersguide_2_Synthetic_model_run.m In this tutorial, we show how TTLEM can be used to simulate synthetic landscape evolution. We illustrate how the user can change between different algorithms to simulate hillslope response and how different numerical schemes can be set. Finally, we also show how modelled data can be processed and converted in a movie of an evolving landscape. Detailed information of parameters, their default values and units can be found in the help section of ttlemset.
- TTLEM_usersguide_3_Synthetic_Geological_Configuration One of the main advantages of LEMs is their potential to explore different geo-tectonic configurations and corresponding hypothesis at a  very low cost. As an example, we simulated landscape evolution over a time span of 30 Myr. Successively, three spatial and temporally dependent uplift scenarios are imposed to the model. On top of this three scenarios, a horizontal shortening field is imposed over the simulated domain with high lateral displacements in the left upper corner, linearly decreasing in both x and y directions to the Southeastern bottom of the simulated domain.

## Requirements

SPLM is plat-form independent and requires Matlab 2014a or higher 

## References

When you use SPLM, please reference to this publication:

Campforts B., Schwanghart W., Govers G.: TTLEM 1.0 : a numerical package for accurate simulation of transient landscape evolution in MATLAB. Discussion paper

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

## Links
- 1D Stream Power Law Models: https://github.com/BCampforts/SPLM

- Benjamin: https://www.researchgate.net/profile/Benjamin_Campforts

- Wolfgang: https://www.researchgate.net/profile/Wolfgang_Schwanghart2

- Blog: https://topotoolbox.wordpress.com/

## Version History

V1.0 --- 25. June 2016 
- release
