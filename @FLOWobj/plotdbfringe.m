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
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = klargestconncomps(S,2);
%     imageschs(DEM,[],'colormap',[1 1 1])
%     hold on
%     h = plotdbfringe(FD,S,'colormap',parula,'width',30);
%
%
% See also: FLOWobj/drainagebasins, GRIDobj/imageschs
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 16. December, 2021
    

p = inputParser;
addRequired(p,'FD')
addOptional(p,'S',[])
addParameter(p,'colormap',parula)
addParameter(p,'shuffle',false)
addParameter(p,'width',10)
addParameter(p,'maxalpha',0.9,@(x) x> 0 && x <= 1);
parse(p,FD,varargin{:})


FD = p.Results.FD;

if isempty(p.Results.S)
    D = drainagebasins(FD);
else
    D = drainagebasins(FD,p.Results.S);
end

if p.Results.shuffle
    D = shufflelabel(D);
end

clr    = p.Results.colormap;
clr    = clr*255;
maxD   = double(max(D));
colors = interp1(linspace(1,maxD,size(clr,1))',clr,1:maxD);

colors = uint8(colors);

iszero  = isnan(D.Z) | D.Z == 0;
RGB    = repmat(uint8(255),prod(FD.size),3);
RGB(~iszero,:)    = colors(D.Z(~iszero(:)),:);

RGB    = reshape(RGB,[FD.size 3]);

BDS    = imerode(D.Z,ones(3)) ~= D.Z | imdilate(D.Z,ones(3)) ~= D.Z; 
DIST   = bwdist(BDS,'quasi-euclidean');

ALPHA    = 1 - min(DIST*1/p.Results.width,1); 
ALPHA(iszero) = 0;
ALPHA    = ALPHA*p.Results.maxalpha;

[x,y] = getcoordinates(D);
h = image(x,y,RGB);
h.AlphaData = ALPHA;
axis xy

if nargout == 1
    him = h;
end


