function C = griddedcontour(DEM,level)

% plot contours on grid
%
% Syntax
%
%     C = griddedcontour(DEM,n)
%     C = griddedcontour(DEM,levels)
%
% Description
%
%     griddedcontour wraps the function contour and approximates the lines
%     on a grid.
%
% Input arguments
%
%     DEM       grid (class: GRIDobj)
%     n         number of contour levels
%     levels    vector with levels (if only one specific level should be 
%               returned, use [level level]).
%  
% Output arguments
%
%     C         instance of GRIDobj
%
% Example
%
%     C = griddedcontour(DEM,[100:100:2000]);
%     imageschs(DEM,dilate(C,ones(3)));
%
% See also: GRIDobj, GRIDobj/imageschs, GRIDobj/contour
%
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. August, 2013

[x,y] = contour(DEM,level);

C = DEM;
C.Z = false(DEM.size);
C.name = 'griddedcontour';

I = ~isnan(x);

C.Z(coord2ind(DEM,x(I),y(I))) = true;