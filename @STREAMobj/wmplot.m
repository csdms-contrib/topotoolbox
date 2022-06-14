function ht = wmplot(S,varargin)

%WMPLOT plot stream network in the webmap browser
%
% Syntax
%
%     h = wmplot(S)
%     h = wmplot(S,'color',c)
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
%     Parameter name/value pairs
%
%     'color'      three element rgb vector or node attribute list
%     'colormap'   char (default is 'jet')
%     'nrcolors'   number of colors
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     wmplot(S)
%
% Example 2
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     k = ksn(S,DEM,A);
%     wmplot(S,'color',k)
%
% Remark
%
%     Another way to quickly plot a stream object in a network is:
%
%     [lat,lon] = STREAMobj2latlon(S);
%     h = wmline(lat,lon,'OverlayName','Stream network','color','k');
%
% See also: STREAMobj, STREAMobj/plot, STREAMobj/STREAMobj2shape, 
%           STREAMobj/plotc
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. October, 2019


p = inputParser;
addParameter(p,'color','k')
addParameter(p,'colormap','jet')
addParameter(p,'nrcolors',20)
parse(p,varargin{:});


if ~isnal(S,p.Results.color) && ~isa(p.Results.color,'GRIDobj')
    % STREAMobj to lat lon
    switch S.georef.SpatialRef.CoordinateSystemType
        case 'geographic'
            [lon,lat] = STREAMobj2XY(S);
        otherwise
            [lat,lon] = STREAMobj2latlon(S);
    end
    
    minlat = min(lat);
    maxlat = max(lat);
    minlon = min(lon);
    maxlon = max(lon);
    
    % wm = webmap
%     wmlimits(wm,[minlat maxlat],[minlon maxlon]);
    
    h = wmline(lat,lon,'OverlayName','Stream network','Color',p.Results.color);
    
    if nargout == 1
        ht = h;
    end
    
else
    
    [GS,val] = STREAMobj2shape(S,'type','geo',...
        'attributes',{...
        'val' p.Results.color @mean},...
        'summarizeby','val','by',p.Results.nrcolors);
    
    if ischar(p.Results.colormap)
        cmapfun = str2func(p.Results.colormap);
        clr = cmapfun(p.Results.nrcolors);
    else
        cmap = p.Results.colormap;
        clr  = interp1(1:size(cmap,1),cmap,linspace(1,size(cmap,1),p.Results.nrcolors));
    end

    counter = 0;
    for r = 1:numel(val)
        if ~isempty(GS{r})   
            counter = counter + 1;
            h{counter} = wmline(GS{r},'color',clr(r,:));
        end
    end
    if nargout == 1
        ht = h;
    end
end
% wmlimits([minlat maxlat],[minlon maxlon]);


