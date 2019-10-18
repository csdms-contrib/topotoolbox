function CA = cellarea(DEM,unit)

%CELLAREA calculate cell areas of a GRIDobj in geographic coordinate system
%
% Syntax
%
%     CA = cellarea(DEM)
%     CA = cellarea(DEM,unit)
%
% Description
%
%     cellarea returns the area for each cell in the DEM with a geographic 
%     coordinate system (requires the mapping toolbox).
%
% Input arguments
%
%     DEM     GRIDobj
%     unit    'm' (default) or 'km' (returns either m^2 or km^2)
%
% Output arguments
%
%     CA      cell areas (GRIDobj)
%
% 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 20. May, 2016

if nargin == 1
    unit = 'm';
end

try 
	GST = DEM.georef.SpatialRef.CoordinateSystemType;
catch
	GST = 'geographic';
end

switch GST
    case 'geographic'
    otherwise
        error('DEM has a projected coordinate system')
end


switch lower(unit)
    case 'm'
        scale = 1e6;
    case 'km'
        scale = 1;
    otherwise
        error('unknown unit')
end

total_surface_area = 510072000; %km^2

[~,ca] = areamat(true(DEM.size),DEM.refmat);
ca = total_surface_area * ca * scale;

CA   = DEM;
CA.Z = repmat(ca,1,DEM.size(2));




