function [dn,z,x,y] = demprofile(DEM,n,x,y)

% get profile along path
%
% Syntax
%
%     [dn,z] = demprofile(DEM,n)
%     [dn,z,x,y] = demprofile(DEM,n,x,y)
%
% Description
%
%     demprofile enables to interactively or programmetically derive
%     elevation profiles from DEMs.
%
% Input arguments
%
%     DEM    digital elevation model (class: GRIDobj)
%     n      number of points along profile
%     x,y    coordinate vectors of profile vertices
%
% Output arguments
%
%     dn     distance along profile
%     z      elevation values interpolated (linear) onto profile
%     x,y    coordinate vectors of profile 
% 
% See also: GRIDobj/interp, GRIDobj/measure
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 1. November, 2013


% interactive
if nargin <= 2;
    h = imagesc(DEM);
    ax = get(h,'Parent');

    [x, y] = getline(ax);
    do = getdistance(x,y);
    
    if nargin == 1;
        n = ceil(do(end)/DEM.cellsize)*2;
    end
    
    dn = linspace(0,do(end),n);
    dn = dn(:);
    xy = interp1(do,[x y],dn,'linear');
    
    x  = xy(:,1);
    y  = xy(:,2);
    
else
    if n ~= numel(x);
        do = getdistance(x,y);
        dn = linspace(0,do(end),n);
        dn = dn(:);
        xy = interp1(do,[x y],dn,'linear');
    
        x  = xy(:,1);
        y  = xy(:,2);
    else
    
        dn = getdistance(x,y);
    end
end

z = interp(DEM,x,y);

end
    


function cumdxy = getdistance(x,y)

% cumulative distance along path defined by vertice coordinates
%
% Syntax
%
%      D = getdistance(x,y)
%
%
% Author: Dirk Scherler (scherler[at]@caltech.edu)
% Date: June 2013
    x = x(:);
    y = y(:);

    dx = diff(x);
    dy = diff(y);
    dxy = hypot(dx, dy); % square root of sum of squares
    cumdxy = [0; cumsum(dxy,1)];
end % 