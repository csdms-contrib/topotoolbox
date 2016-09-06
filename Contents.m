% TopoToolbox
% Version 2.2 (Prerelease) 27-Jul-2016
%
% TopoToolbox provides a set of Matlab functions that support the analysis
% of relief and flow pathways in digital elevation models. The major 
% aim of TopoToolbox is to offer helpful analytical GIS utilities in a 
% non-GIS environment in order to support the simultaneous application
% of GIS-specific and other quantitative methods.
%
% If you have any questions or remarks, please contact the authors:
%
% Wolfgang Schwanghart
% w.schwanghart[at]geo.uni-potsdam.de
%
% Dirk Scherler
% scherler[at]caltech.edu
% 
% When you use TopoToolbox in your work, please reference one of these 
% publications:
% 
% Schwanghart, W. & Scherler, D. (2014): TopoToolbox 2 – MATLAB-based 
% software for topographic analysis and modeling in Earth surface sciences. 
% Earth Surface Dynamics, 2, 1-7. [DOI: 10.5194/esurf-2-1-2014]
% 
% Schwanghart, W., Kuhn, N. J. (2010): TopoToolbox: a set of Matlab 
% functions for topographic analysis. Environmental Modelling & Software, 
% 25, 770-781. [DOI: 10.1016/j.envsoft.2009.12.002]
%
%
% Objects
%
%     GRIDobj         - object for gridded, geospatial data
%     FLOWobj         - object for flow direction
%     STREAMobj       - object for stream (channel) networks
%
% Graphical user interfaces
%
%     flowpathapp     - Map, visualize and export flowpaths that start at manually set channelheads                     
%     preprocessapp   - interactive tool for hydrological conditioning
%     slopeareatool   - Interactively create slope area plots and fit power laws                    
%     topoapp         - Create instance of a topoapp
%
% GRIDobj methods
%
%     GRIDOBJ           : Create instance of a GRIDobj
%     GRIDOBJ2ASCII     : write/export GRIDobj to ESRI ArcGIS ASCII file
%     GRIDOBJ2GEOTIFF   : Exports an instance of GRIDobj to a geotiff file
%     GRIDOBJ2MAT       : convert GRIDobj to matrix and coordinate vectors
%     GRIDOBJ2PM        : combine several GRIDobj into a predictor matrix
%     GRIDOBJ2POLYGON   : Conversion from drainage basin grid to polygon or polyline
%     ACV               : Anisotropic coefficient of variation (ACV) 
%     ASPECT            : aspect (angle of exposition) from a digital elevation model (GRIDobj)
%     CASTSHADOW        : cast shadow
%     CONTOUR           : contour plot of an instance of GRIDobj
%     COORD2IND         : convert x and y coordinates to linear index
%     COORD2SUB         : convert x and y coordinates to subscripts into a GRIDobj
%     CROP              : crop an instance of GRIDobj with axis-aligned minimum bounding box
%     CURVATURE         : 8-connected neighborhood curvature of a digital elevation model 
%     DEMAREA           : Calculate the corrected surface area of a DEM
%     DEMPROFILE        : get profile along path
%     DILATE            : morphological dilation
%     DISTANCE          : distance transform
%     ELEVATEMINIMA     : elevate regional minima in a DEM to their lowest neighbor
%     ERODE             : morphological erosion
%     EXCESSTOPOGRAPHY  : difference between actual elevations and elevations with threshold slope
%     FILLSINKS         : fill/remove pits, sinks or topographic depressions
%     FILTER            : edge detection using sobel filter 
%     GETCOORDINATES    : get coordinate vectors of an instance of GRIDobj
%     GETOUTLINE        : get or plot extent of GRIDobj
%     GRADIENT8         : 8-connected neighborhood gradient and aspect of a digital elevation model
%     GRIDDEDCONTOUR    : plot contours on grid
%     HILLSHADE         : create hillshading from a digital elevation model (GRIDobj)
%     HYPSCURVE         : plot hypsometric curve of a digital elevation model
%     IDENTIFYFLATS     : identify flat terrain in a digital elevation model
%     IMAGESC           : Scale data in GRIDobj and display as image object
%     IMAGESCHS         : plot hillshade image with overlay
%     IND2COORD         : convert linear index to x and y coordinates
%     INFO              : detailed information on GRIDobj instance
%     INPAINTNANS       : interpolate missing values in a grid (GRIDobj)
%     INTERP            : interpolate to query locations
%     INTERP2GRIDOBJ    : Interpolate scattered data to GRIDobj
%     ISNAN             : returns array elements that are NaNs as logical grid
%     LINE2GRIDOBJ      : convert line to a grid
%     LOCALTOPOGRAPHY   : Local topography
%     MEASURE           : take interactive measurements along a polyline
%     MPOWER            : overloaded power for GRIDobj
%     MRDIVIDE          : overloaded right division for GRIDobj
%     MTIMES            : overloaded multiplication for GRIDobj
%     PAD               : add or remove a border of pixels around a GRIDobj
%     POSTPROCFLATS     : postprocess flat terrain for visualization purpose
%     RECLASSIFY        : generate univariate class intervals for an instance of GRIDobj
%     REPROJECT2UTM     : Reproject DEM with WGS84 coordinate system to UTM-WGS84 
%     RESAMPLE          : resample grid to alter spatial resolution
%     ROUGHNESS         : terrain ruggedness, position and roughness indices of DEMs
%     SHUFFLELABEL      : shufflelabel randomly relabels a label matrix
%     SNAP2STREAM       : snap gauges or pour points to stream raster
%     SUB2COORD         : convert subscripts to x and y coordinates
%     SURF              : surface plot for GRIDobj
%     TOPOSHIELDING     : topographic shielding from cosmic rays
%     VALIDATEALIGNMENT : validates whether instances of GRIDobj are spatially aligned
%     ZSCORE            : standardized z-scores for GRIDobj
% 
% FLOWobj methods
% 
%     FLOWOBJ             : Create flow direction object
%     FLOWOBJ2GRIDOBJ     : create ESRI ArcGIS flow direction grid from FLOWobj
%     FLOWOBJ2M           : convert instance of FLOWobj to flow direction matrix 
%     FLOWOBJ2CELL        : Return cell array of FLOWobjs for individual drainage basins
%     FLOWOBJ2GRADIENT    : gradient along flow direction
%     DEPENDENCEMAP       : upslope area for specific locations in a digital elevation model
%     DRAINAGEBASINS      : drainage basin delineation/catchments
%     DRAINAGEBASINSTATS  : Zonal statistics on drainage basins
%     FIND                : find indices and values of edges in the flow direction graph
%     FLOWACC             : flow accumulation (upslope area, contributing area)
%     FLOWCONVERGENCE     : compute flow convergence of a digital elevation model
%     FLOWDISTANCE        : flow distance in upstream and downstream direction
%     FLOWPATHEXTRACT     : extract linear indices of a single flowpath in a DEM
%     IMPOSEMIN           : minima imposition (carving) along drainage network
%     IND2COORD           : convert linear index to x and y coordinates
%     INFLUENCEMAP        : downslope area for specific locations in a digital elevation model
%     ISMULTI             : check if FD is multi or single flow direction
%     MULTI2SINGLE        : converts multiple to single flow direction
%     SAVEOBJ             : Create flow direction object
%     STREAMORDER         : calculate Strahler Stream Order Grid from FLOWobj
%     STREAMPOI           : stream points of interest
%     UPSLOPESTATS        : upslope statistics of a variable based on the flow direction matrix
%     VALIDATEALIGNMENT   : validates whether instances of FLOWobj and GRIDobj are spatially aligned
%     VERTDISTANCE2STREAM : vertical distance to streams
% 
% STREAMobj methods
% 
%     STREAMOBJ           : Create stream object (STREAMobj)
%     STREAMOBJ2GRIDOBJ   : convert an instance of STREAMobj to an instance of GRIDobj
%     STREAMOBJ2SWATHOBJ  : Create swath profile (SWATHobj) from stream network
%     STREAMOBJ2XY        : convert instance of STREAMobj to NaN-separated X and Y coordinates
%     STREAMOBJ2CELL      : convert instance of STREAMobj to cell array of stream objects
%     STREAMOBJ2LATLON    : convert instance of STREAMobj to NaN-separated geographic coordinates
%     STREAMOBJ2MAPSTRUCT : convert instance of STREAMobj to mapstruct
%     CHITRANSFORM        : Coordinate transformation using the integral approach
%     CONNCOMPS           : labels of connected components (individual trees) in a stream network
%     CUMTRAPZ            : Cumulative trapezoidal numerical integration along a stream network
%     CURVATURE           : curvature or 2nd derivative of a STREAMobj
%     DENSIFY             : Increase number of vertices in stream network using splines
%     DISTANCE            : return node attribute list with distances along the stream network
%     EXTRACTCONNCOMPS    : interactive stream network selection
%     GETNAL              : get node attribute list
%     GRADIENT            : stream gradient
%     IMPOSEMIN           : minima imposition (carving) along stream network
%     INFO                : meta information about STREAMobj
%     INPAINTNANS         : inpaint missing values (nans) in a node attribute list
%     INTERSECT           : intersect different instances of STREAMobj 
%     INTERSECTLOCS       : Derive locations where two STREAMobj start to have a common network
%     ISNAL               : test whether a vector is a node attribute list of a STREAMobj
%     KLARGESTCONNCOMPS   : retain k largest connected components in an instance of STREAMobj
%     MINCOSTHYDROCON     : minimum cost hydrological conditioning
%     MODIFY              : modify instance of STREAMobj to meet user-defined criteria
%     NETWORKSEGMENT      : Identify river segments and compute segment geometry
%     PLOT                : plot instance of STREAMobj
%     PLOT3               : 3d-line plot of a STREAMobj
%     PLOT3D              : 3D plot of a stream network
%     PLOTC               : plot a colored stream network
%     PLOTDZ              : plot upstream distance version elevation of a stream network
%     PLOTSEGMENTGEOMETRY : Plot segment geometry obtained using the function networksegment
%     PLOTSTREAMORDER     : calculate Strahler Stream Order from STREAMobj
%     RANDLOCS            : Random locations along the stream network
%     REMOVESHORTSTREAMS  : Remove first order streams with a length less than specified
%     SIDEBRANCHING       : side branching classification according to Tokunaga (1978)
%     SNAP2STREAM         : snap locations to nearest stream location
%     SPLIT               : split drainage network at predefined locations
%     STREAMORDER         : calculate Strahler Stream Order from STREAMobj
%     STREAMPOI           : stream points of interest
%     TRUNK               : extract trunk stream (longest stream) 
%     UNION               : merge different instances of STREAMobj into a new instance
%     VALIDATEALIGNMENT   : is an instance of STREAMobj is spatially aligned with another object of TopoToolbox
%     WIDENSTREAM         : level elevations adjacent to the stream network
% 
% SWATHobj methods
% 
%     SWATHOBJ2GRIDOBJ : SWATHOBJ2GRIDOBJ creates a GRIDobj with Swath-specific information
%     SWATHOBJ2GDS     : SWATHOBJ2GDS creates a geographic data structure from a SWATHobj
%     CONVERT2LATLON   : CONVERT2LATLON converts spatial fields in SWATHobj to lat,lon
%     MAPSWATH         : MAPSWATH obtains Z values along a SWATHobj from an arbitrary GRIDobj
%     PLOT             : plot instance of SWATHobj
%     PLOTDZ           : PLOTDZ creates distance-elevation plot of SWATHobj
%     PLOTDZM          : PLOTDZM creates a color-coded distance-elevation plot from SWATHobj and
%     PROFILES         : PROFILES obtains profiles from a SWATHobj at distinct positions
%     SWATH2LATLON     : convert spatial fields in SWATHobj to geographic coordinates
%     TIDY             : TIDYS removes overlapping points from SWATHobj
%
% TTLEM (TopoToolbox Landscape Evolution Model)
%
%     TTLEM            : TopoToolbox Landscape Evolution Model: TTLEM
%     TTLEMSET         : set options and parameters for TTLEM
%
% Other tools (utitilies)
%
%     compilemexfiles - compile mex-functions that come with TopoToolbox 2
%     coord2ind       - convert xy coordinates to linear index
%     exaggerate      - elevation exaggeration in a 3D surface plot
%     getdistance     - cumulative distance along path defined by vertice coordinates
%     getextent       - get current axis extent
%     interpline      - distribute vertices evenly along line with irregularly spaced vertices
%     ixneighbors     - neighbor indexing in a matrix
%     label2poly      - plot region outlines with polyline
%     landcolor       - colormap for the display of DEMs
%     setextent       - set current axis extent    
%     showmethods     - displays class method names and H1 lines in the command line
%
% Documentation
%
%   usersguide_1_intro      - introductory script to the toolbox
%   usersguide_2_mfd        - how to calculate Multiple Flow Direction in
%                             TopoToolbox 2
%   usersguide_3_ksn        - calculate ksn values with TopoToolbox 2
%
%
%   ---------------------------------------------------------------------
%     TopoToolbox is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

