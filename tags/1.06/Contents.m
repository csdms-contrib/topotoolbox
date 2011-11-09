% TopoToolbox
% Version 1.06 (R2010a) 05-December-2010
%
% TopoToolbox provides a set of Matlab functions that support the analysis
% of relief and flow pathways in digital elevation models. The major 
% aim of TopoToolbox is to offer helpful analytical GIS utilities in a 
% non-GIS environment in order to support the simultaneous application
% of GIS-specific and other quantitative methods.
%
% If you have any questions or remarks, please contact the author:
%
% Wolfgang Schwanghart
% w.schwanghart[at]unibas.ch
% http://www.physiogeo.unibas.ch/
%
%
% I/O functions
%
%   rasterread      - read/import ESRI ASCII grid
%   rasterwrite     - write/export data to ESRI ASCII file
%
% Digital Elevation Model correction/manipulation/filtering
%
%   acv             - anisotropic coefficient of variation (ACV)
%   crossflats      - cross flats of a digital elevation model
%   deminpaint      - interpolate missing values in a raster
%   demsmooth       - moving average filter for digital elevation models
%   demsobel        - edge detection using sobel filter 
%   fillsinks       - fill/remove pits, sinks or topographic depressions
%   routeflats      - route through flats of a digital elevation model
%   routegeodesic   - use least cost path to route through flats 
%
% Terrain primary and secondary attributes
%
%   aspect          - aspect from a digital elevation model in degrees
%   curvature       - 8-connected neighborhood curvature of a digital 
%                     elevation model 
%   dbentropy       - entropy of drainage basin delineation
%   dependencemap   - drainage area for specific locations in a digital 
%                     elevation model
%   drainagebasins  - segment a digital elevation model in drain basins/
%                     catchments
%   drainagedensity - calculate drainage density of individual drainage
%                     basins
%   entropy         - entropy of a DEM based on the multiple flow direction
%                     matrix
%   ezflowacc       - easy to use flow accumulation algorithm for Digital 
%                     Elevation Models
%   flowacc         - calculate flow accumulation/upslope area from flow 
%                     direction matrix
%   flowacc_lm      - flow accumulation (upslope area) for LARGE Digital 
%                     Elevation Models
%   flowconvergence - compute flow convergence of a digital elevation model
%   flowdir         - multiple and single flow direction algorithm 
%   flowdir_single  - single flow direction algorithm
%   flowdistance    - compute flow distance to cell or to catchment outlet
%   flowdistanceds  - compute maximum downstream flow distance from each 
%                     cell in the DEM 
%   gradient8       - 8-connected neighborhood gradient and aspect of a 
%                     digital elevation model
%   hillshade       - create hillshading from a digital elevation model
%   identifyflats   - identify flat terrain in a digital elevation model
%   influencemap    - downslope area for specific locations in a digital 
%                     elevation model
%   multi2single    - convert multiple to single flow direction matrix
%   multiremfracs   - remove small fractions in a multiple flow direction 
%                     matrix
%   streamorder     - calculate Strahler Stream Order from flow direction 
%                     matrix
%   upslopestats    - upslope statistics of a variable based on the flow 
%                     direction matrix
%
% Terrain indices
%
%   twi             - calculate different topographic wetness indices
%   roughness       - terrain ruggedness, position and roughness indices
%
% Subbasin analysis
%
%   sbplot          - plot subbasins and connectivity between subbasin
%                     gauges
%   sbprops         - measure properties of sub-basins
%   sbstruct        - create structure array for sub-basin analysis
%
% Flow path analysis
%
%   flowpathbuffer  - create a buffer around a stream and extract indices 
%                     of cells
%   flowpathextract - extract linear indices of a single flowpath in a DEM
%
% Other tools
%
%   adjustgauges    - snap gauges or pour points to stream raster
%   coord2ind       - convert xy coordinates to linear index
%   cropmat         - crop a subset of arrays with axis-aligned minimum
%                     bounding box
%   exaggerate      - elevation exaggeration in a 3D surface plot
%   hypscurve       - plot hypsometric curve of a digital elevation model
%   ixneighbors     - neighbor indexing in a matrix
%   ismulti         - determine whether matrix is a multiple flow direction
%                     matrix
%   label2poly      - plot region outlines with polyline
%   M2UV            - calculate horizontal, directional components from
%                     flow direction matrix
%   postprocflats   - postprocess flat terrain for visualization purpose
%   shufflelabel    - shufflelabel randomly relabels a label matrix
%   upslopesidelength - calculate the upslope side length and connectivity
%
% Documentation
%   usersguide_1_intro      - introductory script to the toolbox
%   usersguide_2_flats      - introduction to different flat processing
%                             algorithms
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

