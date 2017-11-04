function [DEM,MASK] = widenstream(S,DEM,varargin)

%WIDENSTREAM level elevations adjacent to the stream network
%
% Syntax
%
%     DEMw = widenstream(S,DEM,nrpx)
%     DEMw = widenstream(S,DEM,widthnal)
%     DEMw = widenstream(S,DEM,xyw,method)
%     [DEMw,MASK] = ...
%
% Description
%
%     This function sets elevations in a digital elevation model adjacent
%     to the stream network to the same elevations as the nearest stream
%     pixel. The function thus allows to level the valley bottom or stream
%     bed. Latter may be frequently needed when working with
%     high-resolution topographic data (e.g. LiDAR) and numerical,
%     hydrodynamic models (e.g. LISFLOOD FP).
%
%     widenstream(S,DEM,nrpx) levels elevations in DEM within a distance of 
%     nrpx pixels around the stream network S (STREAMobj).
%
%     widenstream(S,DEM,xyw,method) levels elevations in DEM by the width
%     in mapunits measured at the locations xy. The locations and width are
%     must be provided as the nx3 matrix xyw where the first two columns
%     contain x and y coordinates, respectively, and the third column is
%     the width. The function snaps the measured locations to the 
%     stream network S and interpolates between them using
%     one-dimensional interpolation (interp1) with the method (all methods
%     accepted by the function interp1). 
%
%     *Note that widenstream(S,DEM,xyw,method) only works with a STREAMobj
%     of a single river reach, e.g. their must be only one channel head.*
%
% Input arguments
%
%     S        STREAMobj     
%     DEM      Digital Elevation Model (GRIDobj)
%     nrpx     number of pixels (half-widths)
%     widthnal node-attribute list containing the width for each stream
%              network node
%     xyw      n-by-3 matrix with x and y coordinates and river width.
%     method   interpolation method ('linear','spline',...). See interp1 
%              for further options.
%
% Output arguments
%
%     DEMw     Carved DEM.
%     MASK     Stream mask
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     DEM = imposemin(S,DEM);
%     IX = [694708  553319  370738  262825  139925  1801]';
%     w  = [100 150 200 150 50 200]';
%     [x,y] = ind2coord(DEM,IX);
%     xyw = [x y w];
%     DEMw = widenstream(S,DEM,xyw,'linear');
%     imageschs(DEMw)
%
%
% See also: interp1, STREAMobj/imposemin 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 2. June, 2015


narginchk(3,4)

I = false(DEM.size);
I(S.IXgrid) = true;

[D,L] = bwdist(I,'e');

if numel(varargin) == 1 && isscalar(varargin{1});
    I = D<=varargin{1};
elseif numel(varargin) == 1 && isnal(S,varargin{1});
    distnal = varargin{1}/DEM.cellsize;
    HG = zeros(DEM.size);
    HG(S.IXgrid) = distnal;
    I = HG(L) >= D;
else
    if numel(streampoi(S,'Channelheads','ix')) > 1;
        error('TopoToolbox:widenstream',...
            ['widenstream works only with a river reach, i.e., there\n'...
             'must not be more than one channel head.'])
    end
    
    if numel(varargin) == 2;
        method = validatestring(varargin{2},...
            {'linear','nearest','next','previous','spline','pchip','cubic'},...
            'widenstream','method',4);
    else
        method = 'linear';
    end
    
    xyw = varargin{1};
    xyw(:,3) = xyw(:,3)/DEM.cellsize/2;
    
    [~,~,IX] = snap2stream(S,xyw(:,1),xyw(:,2));
    d  = distance(S,'max_from_ch');
    [~,locb] = ismember(IX,S.IXgrid);
    
    nrObsPix = accumarray(locb,1);
    if any(nrObsPix>=2)
        warning('multiple observations per pixel, taking averages');
    end
    w    = accumarray(locb,xyw(:,3),size(S.IXgrid),@mean,nan);
    
    I    = isnan(w);
    w(I) = interp1(d(~I),w(~I),d(I),method,'extrap');
    W    = zeros(DEM.size);
    W(S.IXgrid) = w;
    
    I  = D<=W(L);
end
    
DEM.Z(I) = DEM.Z(L(I));

if nargout == 2;
    MASK = GRIDobj(DEM,'logical');
    MASK.Z = I;
end
