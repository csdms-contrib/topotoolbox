function arclen = sph_distance(lat1,lon1,lat2,lon2,ellipsoid,usemap)

%SPH_DISTANCE Calculate distances on a sphere
%
% Syntax
%
%     d = sph_distance(lat1,lon1,lat2,lon2,ellipsoid)
%     d = sph_distance(...,usemap)
%
% Description
%
%     sph_distance calculates the distance between points on a sphere or
%     ellipsoid.
%
% Input arguments
%
%     lat1, lon1  Coordinate pairs of points
%     lat2, lon2  Coordinate pairs of points
%     ellipsoid   referenceEllipsoid or referenceSphere or two-element
%                 vector with  [SemimajorAxis Eccentricity].
%     usemap      {true} or false. If true, then the Mapping Toolbox will 
%                 used. If false, distance calculations rely on the
%                 equations reported in Florinsky (2017).
%
% Output arguments
% 
%     d           distance in m for each lat1/lon1 - lat2/lon2 pair.
%
% Reference
%
%     Florinsky (2017): Spheroidal equal angular DEMs: The specificity of
%     morphometric treatment. Transactions in GIS, 21, 1115-1129. [DOI:
%     10.1111/tgis.12269].
%
% See also: distance, referenceEllipsoid
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 2. March, 2020


if nargin <= 4
    ellipsoid = [6378137 6.356752314245179e+06];
end

if nargin <= 5
    usemap = true;
end

if usemap
    arclen = distance(lat1,lon1,lat2,lon2,ellipsoid,'degrees');
else
    % ellipsoid is provided as two-element vector 
    % [semimajor_axis exccentricity]
    
    % Major axis
    A = ellipsoid(1);
    % Eccentricity
    eo = ellipsoid(2);
    
    % Minor axis
    B = sqrt((1-eo^2)*A^2);
    
    % Second eccentricity
    es = sqrt(A^2 / B^2 - 1);
    C  = A^2 / B;
    
    % mean latitude
    lat1   = deg2rad(lat1);
    lat2   = deg2rad(lat2);
    lon1   = deg2rad(lon1);
    lon2   = deg2rad(lon2);
    latm = .5 * (lat1 + lat2);
    num2 = (es^2).*cos(latm).^2;
    Nm = C./sqrt(1+num2);
    Mm = Nm./(1+num2);
    
    latdiff = lat1-lat2;
    londiff = lon1-lon2;
    
    Q = latdiff .* Mm .* (1 - (es^2 - 2*num2).* latdiff/8 - ...
                         (1+num2).* (londiff*cos(latm)).^2 / 12 - ...
                         (londiff*sin(latm)).^2 / 8);
    P = londiff .* cos(latm) .* Nm * (1 + (1-9*es^2 + 8*num2) .* latdiff.^2 / 24 - ...
                         (londiff.*sin(latm)).^2 / 24);
                     
    arclen = sqrt(Q.^2 + P.^2);
end
    
    