________________________________________________________________________

TopoToolbox - a set of Matlab functions for topographic analysis
________________________________________________________________________


TopoToolbox provides a set of Matlab functions that support the analysis
of relief and flow pathways in digital elevation models. The major 
aim of TopoToolbox is to offer helpful analytical GIS utilities in a 
non-GIS environment in order to support the simultaneous application
of GIS-specific and other quantitative methods.

If you have any questions or remarks, please contact the authors:

Wolfgang Schwanghart
w.schwanghart[at]geo.uni-potsdam.de

Dirk Scherler
scherler[at]caltech.edu

** Requirements

TopoToolbox is plat-form independent and requires
Matlab 2011b or higher and the Image Processing Toolbox.

** References

When you use TopoToolbox in your work, please refer to this publication:
Schwanghart, W., Kuhn, N. J. (2010): TopoToolbox: a set of Matlab 
functions for topographic analysis. Environmental Modelling & Software, 
25, 770-781. [DOI: 10.1016/j.envsoft.2009.12.002]

** Getting started

Before working with TopoToolbox the functions need to be on the search 
path of Matlab. Enter following code into the command line:

addpath C:\path\to\wherever\you\installed\this\TopoToolbox-2
addpath C:\path\to\wherever\you\installed\this\TopoToolbox-2\tools_and_more
addpath C:\path\to\wherever\you\installed\this\TopoToolbox-2\topoapp

You may than either look at the user's guide to get an idea of how the 
toolbox works or you may run some of the many examples in the help block
of each function (help function). 

** Known issues

- preprocessapp may have some unexpected behavior here and then...


** Version History
V2.0.1
- removed bug in GRIDobj
- removed bug with case-sensitivity in some functions
- removed bug with internal drainage option in FLOWobj

V2.0
- new functions STREAMobj/intersect, STREAMobj/union
- new interactive tools in STREAMobj/modify
- new interactive tool GRIDobj/measure



V2.0beta ---
- V2.0 introduces an object oriented approach towards grid representation, flow direction, stream 
networks and swath objects. The performance of various, inparticular flow related, functions was 
increased. Mex-files have been written to increase the speed of some functions and are delivered as 
64bit Windows binaries. They have been compiled on Windows 7 with an Intel processor, so they should be 
compiled before using them, if your system differs. However, compiling is not mandatory, since m 
versions are available, too, which are a little slower. Please refer to the Contents.m file for a 
complete list of functions.


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

V1.04 --- 5. January 2010 - *first public release*
- a bug in flowacc_lm was removed. When a weight matrix W0
  was supplied as additional input argument, W0 was set to 
  dem.
- new function: imageschs
- minor changes to hillshade were made. The algorithm is now
  based on method proposed by Katzil and Doytsher, 2003.
- flowacc allows for another input argument (runoff ratio). 

V1.03 --- 5. November 2009
- sub-basin analysis has been added as new functionality 
  (see sbstruct, sbplot, sbprops)
- new function: adjustgauges

V1.02 --- 30. October 2009  
- major speed enhancement for fillsinks with maximum
  fill depth

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

V1.0 --- 15. March 2009 
- release

# This is my README
