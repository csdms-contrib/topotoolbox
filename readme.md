# TopoToolbox - a set of Matlab functions for topographic analysis

<img src="https://github.com/wschwanghart/topotoolbox/blob/master/topotoolbox.jpg" align="center" height="100">


TopoToolbox provides a set of Matlab functions that support the analysis
of relief and flow pathways in digital elevation models. The major 
aim of TopoToolbox is to offer helpful analytical GIS utilities in a 
non-GIS environment in order to support the simultaneous application
of GIS-specific and other quantitative methods.

If you have any questions or remarks, please contact the authors:

Wolfgang Schwanghart
w.schwanghart[at]geo.uni-potsdam.de

Dirk Scherler
scherler[at]gfz-potsdam.de

## Requirements

TopoToolbox is plat-form independent and requires
Matlab 2014b or higher and the Image Processing Toolbox.

## References

When you use TopoToolbox in your work, please reference one of these 
publications:

Schwanghart, W. & Scherler, D. (2014): TopoToolbox 2 – MATLAB-based 
software for topographic analysis and modeling in Earth surface sciences. 
Earth Surface Dynamics, 2, 1-7. [DOI: 10.5194/esurf-2-1-2014]

Schwanghart, W., Kuhn, N. J. (2010): TopoToolbox: a set of Matlab 
functions for topographic analysis. Environmental Modelling & Software, 
25, 770-781. [DOI: 10.1016/j.envsoft.2009.12.002]

## Getting started

Before working with TopoToolbox the functions need to be on the search 
path of Matlab. Enter following code into the command line:

        addpath C:\path\to\wherever\you\installed\this\TopoToolbox-2
        addpath C:\path\to\wherever\you\installed\this\TopoToolbox-2\utilities
        addpath C:\path\to\wherever\you\installed\this\TopoToolbox-2\topoapp
        addpath C:\path\to\wherever\you\installed\this\TopoToolbox-2\DEMdata

You may then either look at the user's guide to get an idea of how the 
toolbox works or you may run some of the many examples in the help block
of each function (help function). 


## Version History

2.2 pre
- TTLEM is part of TopoToolbox
- modification: GRIDobj way to store referencing information was changed
- modification: FLOWobj now supports multiple flow directions and Dinf.
- modification: several new options for imageschs
- new function: FLOWobj/dbentropy 
- new function: FLOWobj/updatetoposort
- new function: GRIDobj/line2GRIDobj
- new function: GRIDobj/GRIDobj2pm
- new function: GRIDobj/zscore
- new function: GRIDobj/pad
- new function: STREAMobj/aggregate
- new function: STREAMobj/labelreach
- new function: STREAMobj/distance
- new function: STREAMobj/plotc
- new function: STREAMobj/plotdzshaded
- new function: STREAMobj/meanupstream
- new function: STREAMobj/plot3
- new function: STREAMobj/chitransform
- new function: STREAMobj/cumtrapz
- modification: STREAMobj/modify includes option rmconncomps
- new function: STREAMobj/conncomps
- new function: STREAMobj/isnal
- new function: STREAMobj/info
- new function: STREAMobj/plotstreamorder
- new function: STREAMobj/split
- new function: STREAMobj/networksegment
- new function: STREAMobj/maplateral
- new function: STREAMobj/plotsegmentgeometry
- new function: STREAMobj/randlocs
- modification: STREAMobj/streamorder plotting option removed
- modification: STREAMobj/plotdz includes custom distance option
- modification: STREAMobj/distance includes option to derive distance from different
  STREAMobj
- modification: STREAMobj/STREAMobj2cell
- modification: STREAMobj/STREAMobj2mapstruct
- new function: STREAMobj/transformcoords
- new function: FLOWobj/FLOWobj2cell
- update to several FLOWobj methods to avoid speed loss for MATLAB versions newer
  than R2015b
- removed bug in GRIDobj/curvature

***

2.1
- new function: GRIDobj/excesstopography
- new function: GRIDobj/GRIDobj2polygon
- new function: STREAMobj/getnal
- new function: STREAMobj/sidebranching
- new function: STREAMobj/mincosthydrocon
- new function: STREAMobj/intersectlocs
- new function: STREAMobj/densify
- new function: STREAMobj/plot3d
- new function: STREAMobj/widenstream
- new function: demo_modifystreamnet
- modification of the function slopearea
- better performance of FLOWobj when converting from flow direction matrix 
  by using dmperm to perform topological sort  
- new function: GRIDobj/toposhielding
- new function: demarea (incorporation of Juernjakob Dugge's function: 
  http://www.mathworks.com/matlabcentral/fileexchange/42204-calculate-dem-surface-area )
- new function: GRIDobj/getoutline
- removed a bug when some functions such log, log10 were called with integer 
  GRIDobjs
- additional overloading of built-in functions for GRIDobjs. We added matrix
  arithmetics, which, however, perform element-wise operations (e.g. mtimes can
  be used with GRIDobj now, but performs times)
- the scope of the function GRIDobj/localtopography was enhanced (min, max, 
  percentiles, etc in a disk-shaped neighborhood)
- FLOWobj/streampoi and STREAMobj/streampoi now return 'bconfluences', 
  e.g. stream pixels that in downstream direction are located immediately 
  before confluences.
- new function: STREAMobj/imposemin - limits downstream minima imposition
  to the stream network
- several bug fixes
- demo_modifystreamnet.m
- preprocessapp was removed

***

V2.0.1
- removed bug in GRIDobj
- removed bug with case-sensitivity in some functions
- removed bug with internal drainage option in FLOWobj

***

V2.0
- new functions STREAMobj/intersect, STREAMobj/union
- new interactive tools in STREAMobj/modify
- new interactive tool GRIDobj/measure

***

V2.0beta ---
- V2.0 introduces an object oriented approach towards grid representation, flow direction, stream 
networks and swath objects. The performance of various, inparticular flow related, functions was 
increased. Mex-files have been written to increase the speed of some functions and are delivered as 
64bit Windows binaries. They have been compiled on Windows 7 with an Intel processor, so they should be 
compiled before using them, if your system differs. However, compiling is not mandatory, since m 
versions are available, too, which are a little slower. Please refer to the Contents.m file for a 
complete list of functions.

***

V1.06 --- 11. November 2011
- new function: acv
- new function: cropmat
- new function: dbentropy
- new function: deminpaint
- new function: exaggerate
- new function: label2poly
- new function: routegeodesic (as optimal method for routing through flats and 
  pits). Requires Matlab 2011b and will be made the default routing algorithm
  in future releases.
- new function: upslopesidelength
- new function: upslopestats
- new function: getextent and setextent
- function enhancements: rasterread and rasterwrite. 
- function enhancements: roughness
- new users guide on processing flat areas
- the baranja_hill.mat dataset was added. It was obtained from 
  here: http://geomorphometry.org/content/baranja-hill

***

V1.05 --- 15. September 2010
- some of the functions now employ the function validateattributes to check
  input arguments. Note that this might return a warning on older versions than
  2009a.
- a bug in routeflats was removed
- new function: M2UV
- new function: twi
- new function: aspect
- new function: ismulti
- new function: postprocflats
- new function: demsobel
- new function: hypscurve
- new function: roughness
- removed bug in flowacc. flowacc returned an error when called with three
  input arguments
- gradient8 allows you to output different angular units, major speed 
  increase when called with gradient as sole output argument.
- functional enhancements to flowdistance (see help flowdistance)
- functional enhancements to identifyflats
- rasterread and rasterwrite now feature dialog boxes for reading and
  saving files if no filenames are supplied to the function
- flowdistanceds can now calculate the maximum downward flowpath distance 
  for each cell at one step. 

***

V1.04 --- 5. January 2010 - *first public release*
- a bug in flowacc_lm was removed. When a weight matrix W0
  was supplied as additional input argument, W0 was set to 
  dem.
- new function: imageschs
- minor changes to hillshade were made. The algorithm is now
  based on method proposed by Katzil and Doytsher, 2003.
- flowacc allows for another input argument (runoff ratio). 

***

V1.03 --- 5. November 2009
- sub-basin analysis has been added as new functionality 
  (see sbstruct, sbplot, sbprops)
- new function: adjustgauges

***

V1.02 --- 30. October 2009  
- major speed enhancement for fillsinks with maximum
  fill depth

***

V1.01 --- 16. September 2009  
- hillshade plots the hillshade matrix when no output
  arguments are defined.
- flowacc_lm was optimized, so that large, flat areas
  are handled more memory efficiently.
- a bug in routeflats was removed.
- crossflats was updated to run most efficient on 
  Matlab R2009a.
- drainagebasins functionality was enhanced to allow
  for the delineation of drainage basins of specified 
  order. The former support of multiple flow direction
  has been removed.
- new function: drainagedensity
- new function: shufflelabel
- influencemap got a second output Mstreams

***

V1.0 --- 15. March 2009 
- release
