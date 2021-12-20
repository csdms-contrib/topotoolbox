function him = plotdbfringe(FD,varargin)

%PLOTDBFRINGE Plot semitransparent fringe around each drainage basin
%
% Syntax
%
%     h = plotdbfringe(FD)
%     h = plotdbfringe(FD,S,'pn','pv')
%
% Description
%
%     plotdbfringe plots an RGB image with drainage basin outlines which
%     become increasingly transparent towards the center of the drainage
%     basins. 
%
% Input arguments
%
%     FD        FLOWobj
%     S         STREAMobj (optional)
%     
%     Parameter name/value pairs
%
%     'colormap'   colormap for drainage basin outlines {parula}
%     'shuffle'    randomly shuffle drainage basin IDs ({false} or true)
%     'width'      width of the fringe in pixels (default = 10)
%     'maxalpha'   maximum value of alpha (default = 0.9)
%     'type'       'linear' - linear increase of transparency in inward
%                  direction
%                  'uniform' - uniform transparency (defined by maxalpha)
%                  within a fringe defined by 'width'
%                  other options are 'cosine', 'sine', or 'exp'
%     'complementalpha' invert transparency ({false} or true)
%
% Output arguments
%
%     h         handle to image. To extract alpha values, use A = h.
%
% Example: several drainage basins
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = klargestconncomps(S,2);
%     imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false)
%     hold on
%     h = plotdbfringe(FD,S,'colormap',parula,'width',30);
%     hold off
%
% Example: one drainage basin with a specific color
%
%     S   = klargestconncomps(S);
%     imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false)
%     hold on
%     h = plotdbfringe(FD,S,'colormap',[0 .2 .7],'width',30);
%     hold off
%
% See also: FLOWobj/drainagebasins, GRIDobj/imageschs
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 20. December, 2021
    
% parse inputs
p = inputParser;
p.FunctionName = 'plotdbfringe';
addRequired(p,'FD')
addOptional(p,'S',[])
addParameter(p,'colormap',parula)
addParameter(p,'shuffle',false)
addParameter(p,'width',10)
addParameter(p,'maxalpha',0.9,@(x) x> 0 && x <= 1);
addParameter(p,'type','linearinward')
addParameter(p,'complementalpha',false)
parse(p,FD,varargin{:})


FD = p.Results.FD;

% calculate drainage basins 
if isempty(p.Results.S)
    D = drainagebasins(FD);
else
    D = drainagebasins(FD,p.Results.S);
end

% shuffle labels, if required
if p.Results.shuffle
    D = shufflelabel(D);
end

% create colors and RGB image
clr    = p.Results.colormap;
clr    = clr*255;
maxD   = double(max(D));

if maxD > 1 && size(clr,1) > 1
    % interpolate colors if there is more than one color and several
    % drainagebasins
    colors = interp1(linspace(1,maxD,size(clr,1))',clr,1:maxD);
elseif maxD > 1 && size(clr,1) == 1
    % several drainage basins and only one color
    colors = repmat(clr,maxD,1);
else
    % several colors, but only one drainage basin
    colors = clr(1,:);
end

colors = uint8(colors);

iszero  = isnan(D.Z) | D.Z == 0;
RGB     = repmat(uint8(255),prod(FD.size),3);
RGB(~iszero,:)    = colors(D.Z(~iszero(:)),:);
RGB    = reshape(RGB,[FD.size 3]);

% Calculate transparency based on distance from catchment boundaries
BDS    = imerode(D.Z,ones(3)) ~= D.Z | imdilate(D.Z,ones(3)) ~= D.Z; 
DIST   = bwdist(BDS,'quasi-euclidean');

fringetype = validatestring(p.Results.type,{'linearinward','uniform','sine','cosine','exp'});
switch lower(fringetype)
    case 'linearinward'
        ALPHA    = 1 - min(DIST*1/p.Results.width,1);
        ALPHA(iszero) = 0;
    case 'uniform'
        ALPHA = DIST <= p.Results.width;
    case 'sine'
        ALPHA = sin(DIST/p.Results.width * pi);
        ALPHA(DIST > p.Results.width) = 0;      
    case 'cosine'
        ALPHA = (cos(DIST/(p.Results.width) * pi) + 1)*0.5;
        ALPHA(DIST > p.Results.width) = 0;
    case 'exp'
        ALPHA = exp(-DIST/p.Results.width);
end

if p.Results.complementalpha
    ALPHA = 1-ALPHA;
end
ALPHA = ALPHA*p.Results.maxalpha;

% Plot the RGB image
[x,y] = getcoordinates(D);
h = image(x,y,RGB);
h.AlphaData = ALPHA;
axis xy

if nargout == 1
    him = h;
end


