function [OUT] = swath2latlon(SW)
% convert spatial fields in SWATHobj to geographic coordinates
%
% Syntax
%
%     OUT = swath2latlon(SW)
%
% Description
%
%     SWATH2LATLON uses the mapping toolbox to convert all spatial fields
%     in the SWATHobj SW to geographic coordinates (latitude, longitude).
%     Requires Mapping toolbox.
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
%     SW2 = swath2latlon(SW)
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013



if strcmp(SW.georef.ModelType,'ModelTypeGeographic')
    warning('Warning: ModelType of SWATHobj data appears to be in geographic coordinates already. Nothing done.')
    
else
    
    OUT = SW;
    
    % coordinate vectors for SWATHobj
    nrrows = SW.georef.Height;
    nrcols = SW.georef.Width;
    R = SW.georef.RefMatrix;
    x = [ones(nrcols,1) (1:nrcols)' ones(nrcols,1)]*R;
    x = x(:,1)';
    y = [(1:nrrows)' ones(nrrows,2)]*R;
    y = y(:,2);
    
    for i = 1 : length(SW.xy0)
        % convert swath data to lat,lon
        [lat0,lon0] = projinv(SW.georef,SW.xy0{i}(:,1),SW.xy0{i}(:,2));
        OUT.xy0{i} = [lon0,lat0];

        [lat,lon] = projinv(SW.georef,SW.xy{i}(:,1),SW.xy{i}(:,2));
        OUT.xy{i} = [lon,lat];

        ix = find(SW.X{i});
        [LAT,LON] = projinv(SW.georef,SW.X{i}(ix),SW.Y{i}(ix));
        OUT.X{i}(ix) = LON(ix);
        OUT.Y{i}(ix) = LAT(ix);
    end
    
    % still need to change georef!
    OUT.georef = [];
    
end
    