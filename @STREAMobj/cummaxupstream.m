function z = cummaxupstream(S,DEM)

%CUMMAXUPSTREAM cumulative maximum in upstream direction
%
% Syntax
%
%     zm = cummaxupstream(S,DEM)
%     zm = cummaxupstream(S,z)
%
% Description
%
%     CUMMAXUPSTREAM calculates the cumulative maximum in upstream
%     direction. A value in zm will have the maximum value of its
%     downstream pixels in DEM or z. 
%
% Input arguments
%
%     S     STREAMobj
%     DEM   GRIDobj
%     z     node-attribute list
%
% Output arguments
%
%     zm    node-attribute list
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     zm = cummaxupstream(S,DEM);
%     plotdz(S,DEM)
%     hold on
%     plotdz(S,zm)
%     hold off
%
% See also: STREAMobj/imposemin, cummax
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. May, 2017

narginchk(2,2)

if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

ix = S.ix;
ixc = S.ixc;
for r = numel(ix):-1:1
    z(ix(r)) = max(z(ixc(r)),z(ix(r)));
end

