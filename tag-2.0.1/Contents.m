% TopoToolbox
% Version 2.0 (R2012a) 13-December-2013
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
%     GRIDOBJ           Create instance of a GRIDobj
%     GRIDOBJ2ASCII     write/export GRIDobj to ESRI ArcGIS ASCII file
%     GRIDOBJ2GEOTIFF   Exports an instance of GRIDobj to a geotiff file
%     GRIDOBJ2MAT       convert GRIDobj to matrix and coordinate vectors
%     ACV               Anisotropic coefficient of variation (ACV) 
%     ASPECT            aspect (angle of exposition) from a digital elevation model (GRIDobj)
%     CASTSHADOW        cast shadow
%     CONTOUR           contour plot of an instance of GRIDobj
%     COORD2IND         convert x and y coordinates to linear index
%     COORD2SUB         convert x and y coordinates to subscripts into a GRIDobj
%     CROP              crop an instance of GRIDobj with axis-aligned minimum bounding box
%     CURVATURE         8-connected neighborhood curvature of a digital elevation model 
%     DEMPROFILE        get profile along path
%     DILATE            morphological dilation
%     ELEVATEMINIMA     elevate regional minima in a DEM to their lowest neighbor
%     ERODE             morphological erosion
%     FILLSINKS         fill/remove pits, sinks or topographic depressions
%     FILTER            edge detection using sobel filter 
%     GETCOORDINATES    get coordinate vectors of an instance of GRIDobj
%     GRADIENT8         8-connected neighborhood gradient and aspect of a digital elevation model
%     GRIDDATA          Use different techniques to spatially interpolate values at grid locations
%     GRIDDEDCONTOUR    plot contours on grid
%     HILLSHADE         create hillshading from a digital elevation model (GRIDobj)
%     HYPSCURVE         plot hypsometric curve of a digital elevation model
%     IDENTIFYFLATS     identify flat terrain in a digital elevation model
%     IMAGESC           Scale data in GRIDobj and display as image object
%     IMAGESCHS         plot hillshade image with overlay
%     IND2COORD         convert linear index to x and y coordinates
%     INFO              detailed information on GRIDobj instance
%     INPAINTNANS       interpolate missing values in a grid (GRIDobj)
%     INTERP            interpolate to query locations
%     ISNAN             returns array elements that are NaNs as logical grid
%     LOCALTOPOGRAPHY   Local topography
%     MEASURE           take interactive measurements along a polyline
%     POSTPROCFLATS     postprocess flat terrain for visualization purpose
%     RECLASSIFY        generate univariate class intervals for an instance of GRIDobj
%     RESAMPLE          resample grid to alter spatial resolution
%     ROUGHNESS         terrain ruggedness, position and roughness indices of DEMs
%     SHUFFLELABEL      shufflelabel randomly relabels a label matrix
%     SNAP2STREAM       snap gauges or pour points to stream raster
%     SUB2COORD         convert subscripts to x and y coordinates
%     SURF              surface plot for GRIDobj
%     VALIDATEALIGNMENT validates whether instances of GRIDobj are spatially aligned
% 
% FLOWobj methods
% 
%     FLOWOBJ             Create flow direction object
%     FLOWOBJ2GRIDOBJ     create ESRI ArcGIS flow direction grid from FLOWobj
%     FLOWOBJ2M           convert instance of FLOWobj to flow direction matrix 
%     FLOWOBJ2GRADIENT    gradient along flow direction
%     DEPENDENCEMAP       upslope area for specific locations in a digital elevation model
%     DRAINAGEBASINS      drainage basin delineation/catchments
%     FIND                find indices and values of edges in the flow direction graph
%     FLOWACC             flow accumulation (upslope area, contributing area)
%     FLOWCONVERGENCE     compute flow convergence of a digital elevation model
%     FLOWDISTANCE        flow distance in upstream and downstream direction
%     FLOWPATHEXTRACT     extract linear indices of a single flowpath in a DEM
%     IMPOSEMIN           minima imposition (carving) along drainage network
%     IND2COORD           convert linear index to x and y coordinates
%     INFLUENCEMAP        downslope area for specific locations in a digital elevation model
%     ISMULTI             check if FD is multi or single flow direction
%     SAVEOBJ             Create flow direction object
%     STREAMLINKS         FD is reduced to contain flow directions in channel locations only
%     STREAMORDER         calculate Strahler Stream Order Grid from FLOWobj
%     STREAMPOI           stream points of interest
%     UPSLOPESTATS        upslope statistics of a variable based on the flow direction matrix
%     VALIDATEALIGNMENT   validates whether instances of FLOWobj and GRIDobj are spatially aligned
%     VERTDISTANCE2STREAM vertical distance to streams
% 
% STREAMobj methods
% 
%     STREAMOBJ           Create stream object (STREAMobj)
%     STREAMOBJ2GRIDOBJ   convert an instance of STREAMobj to an instance of GRIDobj
%     STREAMOBJ2XY        convert instance of STREAMobj to NaN-separated X and Y coordinates
%     STREAMOBJ2CELL      convert instance of STREAMobj to cell array of stream objects
%     STREAMOBJ2MAPSTRUCT convert instance of STREAMobj to mapstruct
%     DISTANCE            handle the simple case which is implement as dynamic property of S
%     EXTRACTCONNCOMPS    interactive stream network selection
%     GRADIENT            stream gradient
%     INTERSECT           intersect different instances of STREAMobj 
%     KLARGESTCONNCOMPS   retain k largest connected components in an instance of STREAMobj
%     MODIFY              modify instance of STREAMobj to meet user-defined criteria
%     PLOT                plot instance of STREAMobj
%     PLOTDZ              plot upstream distance version elevation of a stream network
%     REMOVESHORTSTREAMS  Remove first order streams with a length less than specified
%     SNAP2STREAM         snap locations to nearest stream location
%     STREAMMETRICS       metrics
%     STREAMORDER         calculate Strahler Stream Order from STREAMobj
%     STREAMPOI           stream points of interest
%     TRUNK               extract trunk stream (longest stream) 
%     UNION               merge different instances of STREAMobj into a new instance
%     VALIDATEALIGNMENT   is an instance of STREAMobj is spatially aligned with another object of TopoToolbox
% 
% SWATHobj methods
% 
%     SWATHOBJ           Create swath profile object (SWATHobj)
%     SWATHOBJ2GRIDOBJ   convert SWATHobj to GRIDobj with swath-specific information
%     SWATHOBJ2MAPSTRUCT convert instance of SWATHobj to mappstruct
%     BASELEVEL          creates a SWATHobj with the local across-swath minimum
%     CHOPSWATH          chops a SWATHobj into equally long segments
%     MAPSWATH           obtains Z values along a SWATHobj from an arbitrary GRIDobj
%     MODIFY             sets fields in a SWATHobj to NaN according to different rules
%     PLOT               plot instance of SWATHobj in map view
%     PLOTDZ             plot instance of SWATHobj in profile view
%     PLOTDZM            plot instance of SWATHobj in profile view, color-coded based on GRIDobj
%     PLOTZ              PLOTZ creates plot of statistics calculated tranverse to SW(i)ath profile
%     SWATH2LATLON       convert spatial fields in SWATHobj to lat,lon
%     TIDYSWATH          removes overlapping points from SWATHobj
%
% Other tools
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

