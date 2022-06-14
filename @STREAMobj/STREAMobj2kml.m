function kout = STREAMobj2kml(S,varargin)

%STREAMobj2kml Convert STREAMobj to kml (Google Earth)
%
% Syntax
%
%     STREAMobj2kml(S)
%     STREAMobj2kml(S,'pn','pv',...)
%     k = STREAMobj2kml(...)
%
% Description
%
%     STREAMobj2kml converts a stream network stored in a STREAMobj to kml
%     (or kmz) that can be opened in Google Earth and other GIS software.
%     The function requires the kml toolbox which can be found here:
%     https://github.com/rafael-aero/kml-toolbox
%
%     STREAMobj2kml enables plotting of stream networks in Google Earth 
%     that are colored by an attribute (e.g. elevation, gradient, ...).
%
%     By default, and if called without output arguments, STREAMobj2kml
%     saves the stream network to streamnet.kmz in the current folder. With
%     output argument, STREAMobj2kml will not save the kmz file but output
%     a kml object which can be enhanced further using the kml toolbox.
%
% Input arguments
%
%     S     STREAMobj
%
%     Parameter name/value pairs
%
%     'filename'       string/char 
%                      Default is 'streamnet.kmz'
%     'seglength'      Line segment length into which the stream network is
%                      split
%     'attribute'      GRIDobj or node-attribute list that will be used as
%                      attribute for each line segment. Default is
%                      S.distance. May not contain nans or infs.
%     'aggfun'         Anonymous function that will be used to aggregate
%                      node-attribute lists to segment attributes. By
%                      default, this is @mean.
%     'colormap'       string/char, anonymous function or colormap matrix
%                      used for plotting.
%     'nrcolors'       Number of colorvalues. Default is 20.
%     'linewidth'      Linewidth. Default is 2.
%     'name'           String/char. Name of the kmz file as it appears in
%                      Google Earth.
%
% Output arguments
%
%     k        kml object (see kml toolbox)
%
% Example: Plot a chimap in Google Earth
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     c = chitransform(S,A);
%     [~,zb] = zerobaselevel(S,DEM);
%     c = c+zb;    
%     STREAMobj2kml(S,'attribute',c)
%     % Open streamnet.kmz in Google Earth
%
% See also: STREAMobj2mapstruct, STREAMobj2shape, wmplot, kml
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 28. October, 2019

% Check if kml exists
if ~(exist('kml','class')>0)
    error('kml toolbox is not installed or not on the search path.')
end
    
d = S.distance;
defaultseglength = max(max(d)/20,5*S.cellsize);

p = inputParser;
addParameter(p,'filename','streamnet.kmz');
addParameter(p,'seglength',defaultseglength,@(x) isscalar(x) && x>S.cellsize*3);
addParameter(p,'attribute',d)
addParameter(p,'aggfun',@mean)
addParameter(p,'colormap','parula')
addParameter(p,'nrcolors',20);
addParameter(p,'name','Stream network')
addParameter(p,'linewidth',2)
parse(p,varargin{:});


MS = STREAMobj2mapstruct(S,'seglength',p.Results.seglength,...
    'attributes',{'Z' p.Results.attribute p.Results.aggfun});

[lat,lon] = cellfun(@(x,y) minvtran(S.georef.mstruct,x,y),{MS.X},{MS.Y},'UniformOutput',false);
[MS.Lat] = deal(lat{:});
[MS.Lon] = deal(lon{:});
        
MS = rmfield(MS,{'X' 'Y'});

% handle colormap
if ischar(p.Results.colormap)
    cmapfun = str2func(p.Results.colormap);
    map     = cmapfun(p.Results.nrcolors);
else
    cmap = p.Results.colormap;
    map  = interp1(1:size(cmap,1),cmap,linspace(1,size(cmap,1),p.Results.nrcolors));
end

% Gearth wants HEX colors
N = p.Results.nrcolors;        % Number of colors
hexMap = reshape(dec2hex(uint8([255.*ones(N, 1) 255.*flip(map, 2)]).').', 8, []).';

[~,~,bin] = histcounts([MS.Z],p.Results.nrcolors);

if any(bin==0)
    error('Aggregated attributes contain NaNs or Infs')
end

% initiate kml
k = kml(p.Results.name);

% go through line segments in MS
for r = 1:numel(MS)
    k.plot3(MS(r).Lon,MS(r).Lat,zeros(size(MS(r).Lon)),...
        'lineColor',hexMap(bin(r),:),...
        'altitudeMode','relativeToGround',...
        'lineWidth',p.Results.linewidth);
end

% Save to kml
if nargout == 0
    k.save(p.Results.filename);
else
    kout = k;
end
