function [z,zbase] = zerobaselevel(S,DEM)

%ZEROBASELEVEL set base level to zero
%
% Syntax
%
%     z0 = zerobaselevel(S,DEM)
%     z0 = zerobaselevel(S,z)
%     [z0,zb] = ...
%
% Description
%
%     ZEROBASELEVEL returns a node-attribute list of elevation values where
%     all values are adjusted so that stream outlet elevations are zero.
%
% Input arguments
%
%     S      STREAMobj
%     DEM    digital elevation model (GRIDobj)
%     z      node-attribute list of elevation values
%
% Output arguments
%
%     z0     node-attribute list of elevation values
%     zb     node-attribute list with base level values. zb is calculated 
%            by zb = z-z0 or zb = getnal(S,DEM)-z0.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S  = STREAMobj(FD,'minarea',1e6,'unit','map');
%     subplot(2,1,1)
%     plotdz(S,DEM);
%     z0 = zerobaselevel(S,DEM);
%     subplot(2,1,2)
%     plotdz(S,z0)
%
% See also: STREAMobj, STREAMobj/plotdz, STREAMobj/getnal
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 12. October, 2017

narginchk(2,2)

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

% set the base level of all streams to zero
zbase = z;
for r = numel(S.ixc):-1:1
    zbase(S.ix(r)) = zbase(S.ixc(r));
end
z = z-zbase;