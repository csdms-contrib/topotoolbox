function [OUT] = convert2latlon(SW)
%CONVERT2LATLON converts spatial fields in SWATHobj to lat,lon
%
% Requires the Mapping Toolbox
%
% Syntax
%
%     OUT = convert2latlon(SW)
%
% Description
%
%     CONVERT2LATLON uses the mapping toolbox to convert all spatial fields
%     in the SWATHobj SW to geographic coordinates (latitude, longitude).
%     Note that the field 'georef' is not updated to the geographic
%     coordinate system, but simply deleted. This still needs to be done...
%
% Input arguments
%
%     SW     instance of SWATHobj
%
% Output arguments
%
%     OUT    instance of SWATHobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM)
%     SW2 = convert2latlon(SW)
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: May, 2015


OUT = SW;

for i = 1 : length(SW.xy0)
    % convert swath data to lat,lon using projinv
    [lat0,lon0] = minvtran(SW.georef.mstruct,SW.xy0(:,1),SW.xy0(:,2));
    OUT.xy0 = [lon0,lat0];
    
    [lat,lon] = minvtran(SW.georef.mstruct,SW.xy(:,1),SW.xy(:,2));
    OUT.xy = [lon,lat];
    
    ix = find(SW.X);
    [LAT,LON] = minvtran(SW.georef.mstruct,SW.X(ix),SW.Y(ix));
    OUT.X(ix) = LON(ix);
    OUT.Y(ix) = LAT(ix);
end

% still need to change georef!
OUT.georef = [];

