function ext = getextent(DEM,latlonout)

%GETEXTENT return extent of a GRIDobj
%
% Syntax
%
%     ext = getextent(DEM,latlonout)
%
% Description
%
%     getextent returns the horizontal extent of the GRIDobj DEM. The
%     function will return the maximum extent in WGS84 geographical
%     coordinates if latlonout is true.
%
% Input arguments
%
%     DEM        GRIDobj
%     latlonout  false (default) or true. If true, getextent will return
%                the extent in geographical coordinates. This is, however,
%                only possible if DEM has a known projected coordinate 
%                system (DEM.georef.mstruct must be set), and if the
%                mapping toolbox is available. 
%
% Output arguments
%
%     ext        four element row vector with following format 
%                [west east south north]
%
%
% See also: GRIDobj, readopentopo, GRIDobj/getoutline
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 9. October, 2018


if nargin == 1
    latlonout = false;
end

if ~latlonout
    % This is easy. Just return the extent
    [x,y] = getcoordinates(DEM);
    ext   = [min(x) max(x) min(y) max(y)];
else
    mstructavailable = isfield(DEM.georef,'mstruct');
    if mstructavailable && ~isempty(DEM.georef.mstruct)
        % mstruct is available. Great.
        [x,y] = getcoordinates(DEM);
        % maximum projected extent
        extp  = [min(x) max(x) min(y) max(y)];
        [lat,lon] = minvtran(DEM.georef.mstruct,...
                            [extp(1) extp(1) extp(2) extp(2)]',...
                            [extp(3) extp(4) extp(3) extp(4)]');
        ext   = [min(lon) max(lon) min(lat) max(lat)];
    else
        [x,y] = getcoordinates(DEM);
        ext   = [min(x) max(x) min(y) max(y)]; 
        if any(ext(1:2) > 180) || any(ext(1:2) < -180) || ...
           any(ext(3:4) > 90)  || any(ext(3:4) < -90)
       
            error('TopoToolbox:getextent',...
                ['GRIDobj does not have a map projection structure.\n' ...
                 'In addition, coordinates are out of range to be in \n'...
                 'geographic coordinates.'])
        end
    end
end
    
    

