function [fs,s] = identifyflats(S,DEM)

%IDENTIFYFLATS identify flat sections in a river profile
%
% Syntax
%
%     [flat,sill] =  identifyflats(S,DEM)
%     [flat,sill] =  identifyflats(S,z)
%
% Description
%
%     identifyflats identifies flat sections in a river profile and returns
%     them as a node attribute list. The second output argument provides
%     the sills, i.e. the locations were flat sections spill over.
%
% Input arguments
%
%     S      STREAMobj
%     DEM    digital elevation model (GRIDobj)
%     z      node-attribute list of elevation values
% 
% Output arguments
%
%     flat   logical node-attribute list with flat areas set to true
%     sill   logical node-attribute list with sills set to true
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S  = STREAMobj(FD,A>1000);
%     S  = klargestconncomps(trunk(S));
%     DEM = fillsinks(DEM);
%     flat = identifyflats(S,DEM);
%     plotdz(S,DEM,'color',+flat)
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016


if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM);
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

fs = false(size(S.IXgrid));
s  = zeros(size(S.IXgrid));

for r=1:numel(S.ix);
    if z(S.ixc(r)) == z(S.ix(r));
%         fs(S.ix(r)) = true;
        fs(S.ixc(r)) = true;
    elseif fs(S.ix(r)) && z(S.ixc(r))< z(S.ix(r))
        s(S.ixc(r)) = 1;
    end
end

I = streampoi(S,'channelheads','logical');
fs(I) = false;
I = streampoi(S,'outlets','logical');
fs(I) = false;
s(I) = 0;
    
    
    