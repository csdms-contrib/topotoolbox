function [IXgn,d] = snap2stream(I,IXg,plotit)

%SNAP2STREAM snap gauges or pour points to stream raster
%
% Syntax
%
%     [IXn,d] = snap2stream(W,IX)
%     [IXn,d] = snap2stream(W,coords)
%     [IXn,d] = snap2stream(...,plotit)
%
% Description
%
%     Often gauge locations are not precisely associated with 
%     streams derived from a digital elevation model. Calculation of 
%     drainage basins may thus return wrong basin delineations. One method 
%     is to manually adjust their location to a derived stream raster.
%     snap2stream does this automatically by finding the locations'
%     shortest distance to a stream. 
%
%     Note that the function STREAMobj/snap2stream is more elaborate and
%     lets you define various options.
%
% Input
%
%     W         stream raster (GRIDobj with logical stream raster), e.g.
%               obtained from thresholding a flow accumulation grid
%     IX        [m x 1] vector with linear indices of gauge locations
%     coords    [m x 2] matrix with x and y coordinates
%     plotit    false (default) or true. If true, the results of the
%               adjustment are plotted 
% 
% Output
%
%     IXn      adjusted linear index of gauge locations
%     d         euclidean pixel distance between old and adjusted
%               gauges (multiply by cellsize) to obtain true distances
%
%
%     
%
% See also: COORD2IND, BWDIST, STREAMobj/snap2stream
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 7. March, 2018

narginchk(2,3)
if nargin == 2;
    plotit = false;
end

if ~isa(I.Z,'logical');
    error('TopoToolbox:GRIDobj','The instance of GRIDobj must contain a logical stream raster.')
end

if size(IXg,2) == 2
    IXg = coord2ind(I,IXg(:,1),IXg(:,2));
    if isempty(IXg)
        error('TopoToolbox:GRIDobj','All coordinates are outside the grid borders.')
    end
end

% Distance transform
[D,L] = bwdist(I.Z);
D = D.*I.cellsize;

IXg  = IXg(:);
d    = D(IXg);
IXgn = L(IXg);

if plotit
    [io,jo] = ind2sub(I.size,IXg);
    [in,jn] = ind2sub(I.size,IXgn);
    
    imshow(I.Z);
    hold on
    plot(jo,io,'xy',jn,in,'sb',[jo';jn'],[io';in'],'b-');
    hold off
    
    legend('old gauges','adjusted gauges','links','location','southeast');
    
end