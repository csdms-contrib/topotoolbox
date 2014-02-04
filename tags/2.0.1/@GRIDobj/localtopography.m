function OUT = localtopography(DEM,radius)

% Local topography
%
% Syntax
%
%     H = localtopography(DEM,radius);
%
% Description
%
%     localtopography quantifies local relief, e.g. the elevation range
%     within a specific radius.
%
% Input arguments
%
%     DEM    digital elevation model (GRIDobj)
%     radius radius of the moving window filter in map units. The default
%            value is 5000 (m).
%
% Output arguments
%
%     H      local topography grid (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     H = localtopography(DEM,500);
%     imageschs(DEM,H)
%
% See also: IMDILATE, IMERODE
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 28. January, 2013
narginchk(1,2)
if nargin == 2;
    validateattributes(radius,{'double'},{'scalar'},'localtopography','radius',2);
else
    radius = 5000;
end

dem = DEM.Z;
cs  = DEM.cellsize;

% any nans
INAN = isnan(dem);
flaginan = any(INAN(:));

% structuring element
radiuspx = ceil(radius/cs);
SE = strel('disk',radiuspx,0);

% Maximum filter
if flaginan;
    dem(INAN) = -inf;
end
MA = imdilate(dem,SE);

% Minimum filter
if flaginan;
    dem(INAN) = inf;
end
MI = imerode(dem,SE);

% Difference
H  = MA-MI;
if flaginan;
    H(INAN) = nan;
end

% prepare output
OUT = DEM;
OUT.Z = H;
OUT.name = 'local topography';

% Relation between long-term erosion rate E and local relief H
% E0: erosion rate due to chemical weathering
% K: rate constant
% Hc: limiting relief

% values from Montgomery and Brandon (2003)
% E0 = 0.01 /1000; % m/yr
% K  = 2.5e-4 /1000; % m/yr
% Hc = 1500; % m
% E  = E0 + K*H./(1-(H./Hc).^2);

% Refs: 
% Montgomery, D. & Brandon, M. Topographic controls on erosion rates
% in tectonically active mountain ranges Earth and Planetary Science
% Letters, 2002, 201, 481-489.

% Korup, O.; Clague, J.; Hermanns, R.; Hewitt, K.; Strom, A. & Weidinger,
% J. Giant landslides, topography, and erosion Earth and Planetary Science
% Letters, 2007, 26, 578-589

