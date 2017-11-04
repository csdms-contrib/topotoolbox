function [D] = dist2line(DEM,x0,y0,alpha)
%DIST2LINE labels pixels in a GRIDobj by their distance to a straight line
%
% Syntax
%     D = dist2line(DEM,x0,y0,alpha)
%
% Description
%
%     dist2line(DEM,x0,y0,alpha) computes the orthogonal distance of each 
%     pixel in a GRIDobj to a line that goes through the point(s) at 
%     (x0,y0) and has an angle alpha from north (or the y-axis). Distances
%     are computed using the function 'distancePointLine' by David Legland,
%     which is part of the geom2d library, available at Matlab Central file
%     exchange. The function is included here.
%
% Input
%
%     DEM       GRIDobj
%     x0        x coordinate of point (scalar)
%     y0        y coordinate of point (scalar)
%     alpha     angle in degrees from north, clockwise (scalar)
%
% Output
%
%     D         GRIDobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [x,y] = getcoordinates(DEM);
%     D = dist2line(DEM,mean(x),mean(y),45);
%     imagesc(D), colorbar
%
% See also: GRIDobj
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 25. June, 2014, updated 30. March, 2017;

[x,y] = getcoordinates(DEM);
[X,Y] = meshgrid(x,y);
x = X(:);
y = Y(:);

D = DEM;
D.Z(:) = 0;

% Create line 
dx = 1;
dy = dx/tand(alpha);
LINE = [x0 y0 dx dy];

D.Z(:) = distancePointLine([x,y],LINE);

end


function dist = distancePointLine(point, line)
%DISTANCEPOINTLINE Minimum distance between a point and a line
%
%   D = distancePointLine(POINT, LINE)
%   Return the euclidean distance between line LINE and point POINT. 
%
%   LINE has the form : [x0 y0 dx dy], and POINT is [x y].
%
%   If LINE is N-by-4 array, result is N-by-1 array computes for each line.
%
%   If POINT is N-by-2, then result is computed for each point.
%
%   If both POINT and LINE are array, result is N-by-1, computed for each
%   corresponding point and line.
%
%
%   See also:
%   lines2d, points2d, distancePoints, distancePointEdge
%
%   
%   ---------
%   author : David Legland 
%   INRA - CEPIA URPOI - MIA MathCell
%   created the 24/06/2005
%

%   HISTORY :


if size(line, 1)==1 && size(point, 1)>1
    line = repmat(line, [size(point, 1) 1]);
end

if size(point, 1)==1 && size(line, 1)>1
    point = repmat(point, [size(line, 1) 1]);
end

dx = line(:, 3);
dy = line(:, 4);

% compute position of points projected on line
tp = ((point(:, 2) - line(:, 2)).*dy + (point(:, 1) - line(:, 1)).*dx) ./ (dx.*dx+dy.*dy);
p0 = line(:, 1:2) + [tp tp].*[dx dy];


% compute distances between points and their projections
dx = point - p0;
dist  = sqrt(sum(dx.*dx, 2));

end


