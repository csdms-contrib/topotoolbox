function [SW] = STREAMobj2SWATHobj(S,DEM,varargin)

%STREAMOBJ2SWATHOBJ Create swath profile (SWATHobj) from stream network
%
% Syntax
%
%     SW = STREAMobj2SWATHobj(S,DEM,'pn','pv',...)
%
% Description
%
%     STREAMobj2SWATHobj creates a swath profile along individual reaches
%     of a stream network. If the SWATHobj was created from a STREAMobj
%     with multiple channels, the resulting SWATHobj's will be stored in
%     cells.
% 
% Input arguments
%
%     S     STREAMobj
%     DEM   digital elevation model (DEM)
%     pn,pv parameter name value pairs as used in SWATHobj
%
% Output arguments
%
%     SW    SWATHobj or cell array of SWATHobjs
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(trunk(S));
%     SW = STREAMobj2SWATHobj(S,DEM);
%     plotdz(SW);
%
%
% See also: SWATHobj
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 15. May, 2017
            

if ~isa(S,'STREAMobj')
    error('First input needs to be of class STREAMobj')
end

if ~isa(DEM,'GRIDobj') 
    error('Second input needs to be of class GRIDobj')
end


validatealignment(S,DEM);
[x,y,d] = STREAMobj2XY(S,S.distance);
% STREAMobj's are sorted downstream, but we want upstream
x = flipud(x);
y = flipud(y);
d = flipud(d);

ixnan = find(isnan(x));
ix1 = ixnan+1;
ix2 = [ixnan(2:end)-1;length(x)];

SW = cell(length(ix1),1);
for i = 1 : length(ix1)
    x_t = x(ix1(i):ix2(i));
    y_t = y(ix1(i):ix2(i));
    d_t = d(ix1(i):ix2(i));
    SW{i} = SWATHobj(DEM,x_t,y_t,d_t,varargin{:});
end

if length(ix1)==1
    SW = SW{1};
end


