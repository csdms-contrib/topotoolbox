function [L,IX] = reclabel(DEM,dx,dy)

%RECLABEL labels GRIDobj by rectangular fields
%
% Syntax
%     L = reclabel(DEM,dx,dy)
%     [L,IX] = reclabel(DEM,dx,dy)
%
% Description
%
%     reclabel(DEM,dx,dy) creates a new GRIDobj L the same size as DEM,
%     where all cells have been grouped into rectangular fields and 
%     labelled from top left to bottom right. dx and dy are the width and 
%     height of the fields in map units of the DEM. 
%     Optional output IX is a new GRIDobj, in which each field represents 
%     one cell. Only entire fields are considered and dx=dy is required.
%
% Input
%
%     DEM       GRIDobj
%     dx        width of rectangular fields (positive scalar)
%     dy        height of rectangular fields (positive scalar)
%
% Output
%
%     L         labeled GRIDobj
%     IX        new GRIDobj, in which each field is one cell
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     L = reclabel(DEM,5e3,5e3);
%     imageschs(DEM,shufflelabel(L))
%
% See also: GRIDobj, checkerboard
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 25. June, 2014

toggle = 1;
if nargout>1 && dx~=dy
    warning('TopoToolbox2: dx and dy are not equal. Output GRIDobj IX will not be assigned.')
    toggle = 0;
end

dr = round(dy./DEM.cellsize);
dc = round(dx./DEM.cellsize);

LABROWS = zeros(DEM.size(1),1,'uint32');
LABROWS(1:dr:DEM.size(1)) = uint32(1);
LABROWS = cumsum(LABROWS);

LABCOLS = zeros(1,DEM.size(2),'uint32');
LABCOLS(1:dr:DEM.size(2)) = uint32(1);
LABCOLS = (cumsum(LABCOLS)-1)*max(LABROWS);

L = DEM;
L.name = 'labelled grid';

try 
    L.Z = LABROWS + LABCOLS;
catch
    % implicit expansion doesn't work. Use bsxfun
    L.Z = bsxfun(@plus,LABROWS,LABCOLS);
end


if nargout>1 && toggle
    
    nrows = DEM.size(1);
    ncols = DEM.size(2);
    nr = length(1:dr:nrows);
    nc = length(1:dc:ncols);
    IX = nan(nr,nc);
    IX(:) = 1:numel(IX);

    [x,y] = getcoordinates(DEM);
    newx = x(1:dc:end);
    newy = y(1:dr:end);

    IX = GRIDobj(newx,newy,IX);

else
    IX = [];
end









