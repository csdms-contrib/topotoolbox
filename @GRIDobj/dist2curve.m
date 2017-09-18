function [OUT] = dist2curve(DEM,x0,y0,alpha,mindist,maxdist)
%DIST2CURVE labels pixels in a GRIDobj by their directed distance to a curved line
%
% Syntax
%
%     D = dist2curve(DEM,x0,y0,alpha,mindist,maxdist)
%
% Description
%
%     dist2curve(DEM,x0,y0,alpha) computes the distance of pixels in a 
%     GRIDobj to an irregular line defined by the points (x0,y0) along the 
%     direction defined by the angle alpha (in degrees from north, or the 
%     y-axis). mindist and maxdist give the minimum and maximum distance
%     from the curve.
%     Note that in the direction of measuring the distance, the curve
%     should have unique values, otherwise results will be erroneous as the
%     function uses interpolation.
%
% Input
%
%     DEM       GRIDobj
%     x0        x coordinate of points on curve (scalars)
%     y0        y coordinate of points on curve (scalars)
%     alpha     angle in degrees from north, clockwise (scalar)
%     mindist   minimum distance from the curve in map units (scalar)
%     maxdist   maximum distance from the curve in map units (scalar)
%
% Output
%
%     D         GRIDobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [x,y] = getcoordinates(DEM);
%     FD = FLOWobj(DEM,'preprocess','carve');
%     [IXc,d,x0,y0] = flowpathextract(FD,629813);
%     D = dist2curve(DEM,x0,y0,0,-1e4,1e4);
%     imageschs(DEM,D), colorbar, hold on
%     plot(x0,y0,'k-'), hold off
%
% See also: GRIDobj
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 25. June, 2014, updated 30. March, 2017;


f = 10;
cs = DEM.cellsize;
step = f*cs;

[dx,dy] = pol2cart(deg2rad(-alpha+90),step);

d0 = getdistance(x0,y0);
[xt,yt,~] = interpline(x0,y0,d0,step);

disty = (mindist-step)/step:(maxdist+step)/step;

ny = length(disty);
nx = length(xt);

allX = repmat(disty',1,nx) .* repmat(dx.*ones(1,nx),ny,1) + repmat(xt',ny,1);
allY = repmat(disty',1,nx) .* repmat(dy.*ones(1,nx),ny,1) + repmat(yt',ny,1);
allD = repmat(disty',1,nx) .* step;

F = scatteredInterpolant(allX(:),allY(:),allD(:),'linear','none');
[x,y] = getcoordinates(DEM);
[X,Y] = meshgrid(x,y);
qz = F(X(:),Y(:));

OUT = DEM;
OUT.Z(:) = nan;
OUT.Z(:) = qz;
OUT.Z(OUT.Z>maxdist | OUT.Z<mindist) = nan;

