function [hl,hp] = wmplot(P,varargin)

%WMPLOT plot instance of PPS in webmap browswer
%
% Syntax
%
%     wmplot(P)
%     [hl,hp] = wmplot(P)
%
% Description
%
%     wmplot plots the object PPS in MATLAB's webmap browser. This
%     requires the Mapping Toolbox and S must have a valid georeferencing.
%
% Input arguments
%
%     P    instance of PPS
%
% Output arguments
%
%     hl    handle to wmline object
%     hp    handle to wmmarker object
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     P = PPS(S,'rpois',0.001);
%     wmplot(P)
%
% See also: STREAMobj, STREAMobj/wmplot, PPS/points
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. May, 2019

% Check input arguments
p = inputParser;
p.FunctionName = 'PPS/wmplot';
addParameter(p,'linecolor',[]);
addParameter(p,'pointcolor',[]);
addParameter(p,'color','k')
addParameter(p,'linewidth',1);
% Parse
parse(p,varargin{:});

p = p.Results;
if isempty(p.linecolor)
    p.linecolor = p.color;
end
if isempty(p.pointcolor)
    p.pointcolor = p.color;
end



[lat,lon] = STREAMobj2latlon(P.S);
minlat = min(lat);
maxlat = max(lat);
minlon = min(lon);
maxlon = max(lon);

wm = webmap;
wmlimits(wm,[minlat maxlat],[minlon maxlon]);

h = wmline(lat,lon,'OverlayName','Stream network','color',p.linecolor);

[lat,lon] = points(P,'latlon');

% get folder of this function
path = mfilename('fullpath');
[filepath,~,~] = fileparts(path);
switch lower(p.pointcolor) 
    case 'r'
        icon     = 'tt_icon_circle.png';
    otherwise
        icon     = 'tt_icon_circle_4.png';
end

iconpath = [filepath filesep 'private' filesep icon];

h2 = wmmarker(lat,lon,...
            'icon',iconpath,...
            'IconScale',0.5,...
            'Alpha',0.8,...
            'OverlayName','Points on network');

if nargout >= 1
    hl = h;
    hp = h2;
end