function [SW] = STREAMobj2SWATHobj(S,DEM,varargin)

% Create swath profile (SWATHobj) from stream network
%
% Syntax
%
%     SW = STREAMobj2SWATHobj(S,DEM,'pn','pv',...)
%
% Description
%
%     STREAMobj2SWATHobj creates a swath profile along individual 
%     reaches of a stream network. If the SWATHobj was created from 
%     a STREAMobj with multiple channels, the resulting SWATHobj's 
%     will be stored in cells, which need to be processed individually when
%     using other functions that work on SWATHobj's.
%
%
% Input arguments
%
%     SW      instance of SWATHobj
%
%     DEM     instance of GRIDobj
%
%     Parameter name/value pairs are the same as in the function SWATHobj
%
%
% Output
%
%     SW     swath profile object (SWATHobj)
%
%
% Example 
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,A>100);
%     S = trunk(klargestconncomps(S,1));
%     SW = STREAMobj2SWATHobj(S,DEM,'width',2e3,'smooth',1e3);
%     figure, plot(SW)
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: March, 2016
            

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


