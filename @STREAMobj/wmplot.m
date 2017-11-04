function ht = wmplot(S,varargin)

%WMPLOT plot stream network in the webmap browser
%
% Syntax
%
%     h = wmplot(S)
%
% Description
%
%     wmplot plots the stream network S in MATLAB's webmap browser. This
%     requires the Mapping Toolbox and S must have a valid georeferencing.
%
% Input arguments
%
%     S    STREAMobj
%
% Output arguments
%
%     h    handle to wmline object
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     wmplot(S)
%
% Remark
%
%     wmplot will plot the stream network in a new webmap browser. If you  
%     want to plot the streams into an existent webmap then use following
%     code
%
%     [lat,lon] = STREAMobj2latlon(S);
%     h = wmline(lat,lon,'OverlayName','Stream network');
%
% See also: STREAMobj, STREAMobj/plot
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. June, 2017


% STREAMobj to lat lon
[lat,lon] = STREAMobj2latlon(S);
minlat = min(lat);
maxlat = max(lat);
minlon = min(lon);
maxlon = max(lon);

wm = webmap;
wmlimits(wm,[minlat maxlat],[minlon maxlon]);

h = wmline(lat,lon,'OverlayName','Stream network');

if nargout == 1
    ht = h;
end




