function dz = diff(S,z)

%DIFF Differences between adjacent pixels in a stream network
%
% Syntax 
% 
%     dz = diff(S,z)
%
% Description
%
%     DIFF calculates the difference of each pixel i and its downstream
%     neighbor j so that dz(i) = dz(i)-dz(j). dz is a node-attribute list.
%     Stream nodes without downstream neighbor (outlets) receive a value of
%     0.
%
% Input arguments
%
%     S     STREAMobj
%     z     node-attribute list
%
% Output arguments
%
%     dz    node-attribute list with differences
%
% Example: Calculate elevation offsets along stream network
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = trunk(klargestconncomps(S));
%     dz = diff(S,DEM);
%     % You can calculate elevation using cumsum
%     zz = cumsum(S,dz,'upstream');
%     
% 
% See also: STREAMobj/gradient, STREAMobj/cumsum
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. November, 2021


% get node attribute list with elevation values
if isa(z,'GRIDobj')
    validatealignment(S,z);
    z = double(getnal(S,z));
elseif isnal(S,z)
    z = double(z);
else
    error('Imcompatible format of second input argument')
end

dz = zeros(size(z));
dz(S.ix) = z(S.ix)-z(S.ixc);