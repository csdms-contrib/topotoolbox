function varargout = points(P,form)

%POINTS extract a list of points from the point pattern
%
% Syntax
%
%     ix = points(P)
%     [x,y] = points(P)
%     ix = points(P,'IXgrid')
%     xy = points(P,'xy');
%     [x,y] = points(P,'xy');
%     latlon = points(P,'latlon');
%     [lat,lon] = points(P,'latlon');
%     nal = points(P,'nal')
%
% Description
%
%     points extracts a list of points from the point pattern. 
%
% Input arguments
%
%     P            instance of PPS
%     outputtype   'IXgrid' linear index into GRIDobj
%                  'xy' coordinate vectors
%                  'latlon' coordinate vectors (geographic coordinates)
%                  'nal' logical node attribute list
%     
% Output
%
%     depending on outputtype
% 
% 
% See also: PPS, PPS/npoints
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019


if nargin == 1
    if nargout == 2
        form = 'xy';
    else
        form = 'ixgrid';
    end
end

switch lower(form)
    case 'ixgrid'
        varargout{1} = P.S.IXgrid(P.PP);
    case 'xy'
        if nargout == 1
            varargout{1} = P.ppxy;
        elseif nargout > 1
            xy = P.ppxy;
            varargout{1} = xy(:,1);
            varargout{2} = xy(:,2);
        end
    case 'latlon'
        [x,y] = points(P);
        [lat,lon] = minvtran(P.S.georef.mstruct,x,y);
        if nargout == 1
            varargout{1} = [lat lon];
        elseif nargout > 1
            varargout{1} = lat;
            varargout{2} = lon;
        end
    case 'nal'    
        
        nal = false(size(P.S.x));
        nal(P.PP) = true;
        varargout{1} = nal;
        
end