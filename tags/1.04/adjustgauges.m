function [IXgn,d] = adjustgauges(IXg,I,plotit)

% snap gauges or pour point to stream raster
%
% Syntax
%
%     [IXgn,d] = adjustgauges(IXg,I)
%     [IXgn,d] = adjustgauges(IXg,I,plotit)
%
% Description
%
%
% Input
%
%     IXg       linear index of gauge locations
%     I         stream raster (logical)
%     plotit    true (default) or false. If true, the results of the
%               adjustment are plotted 
% 
% Output
%
%     IXgn      adjusted linear index of gauge locations
%     d         euclidean pixel distance between old and adjusted
%               gauges (multiply by cellsize) to obtain true distances
%
% Example 
%
%     load exampleDEM;
%     M = flowdir_single(dem);
%     A = flowacc(M,size(dem));
%     W = A>100;
%     IX = [13869 12470 17026 6071];
%     IXn = adjustgauges(IX,W,true);
%
% See also: COORD2IND
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 10. July, 2009

error(nargchk(2, 3, nargin));

if nargin == 2;
    plotit = true;
end

[D,L] = bwdist(I);

IXg  = IXg(:);

d    = D(IXg);
IXgn = L(IXg);

if plotit
    siz     = size(I);
    [io,jo] = ind2sub(siz,IXg);
    [in,jn] = ind2sub(siz,IXgn);
    
    imshow(I);
    hold on
    plot(jo,io,'xy',jn,in,'sb',[jo';jn'],[io';in'],'b-');
    hold off
    
    legend('old gauges','adjusted gauges','links','location','southeast');
    
end








