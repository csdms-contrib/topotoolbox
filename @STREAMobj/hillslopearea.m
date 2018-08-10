function [a,D] = hillslopearea(S,FD)

%HILLSLOPEAREA upslope hillslope area for each stream pixel 
%
% Syntax
%
%     a = hillslopearea(S,FD)
%     [a,D] = hillslopearea(S,FD)
%     
% Description
%
%     hillslopearea returns the upslope hillslope area of each river pixel 
%     in S. Compared to flow accumulation, this function stops accumulation
%     along river pixels so that the accumulated flow calculated for each
%     river pixels includes only the hillslope pixels but not those further
%     upstream along the stream network. 
%
% Input arguments
%
%     S     STREAMobj
%     FD    FLOWobj
%
% Output arguments
%
%     a     node-attribute list with hillslope areas
%     D     GRIDobj with drainage basins for each river pixel
%
%
% See also: FLOWobj/flowacc, FLOWobj/upslopestats
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 10. August, 2018    

A = upslopestats(FD,GRIDobj(FD)+1,'sum',S);
a = getnal(S,A);

if nargout == 2
    D = drainagebasins(FD,S.IXgrid);
end