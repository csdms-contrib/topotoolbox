function h = convhull(P,varargin)

%CONVHULL Convex hull around points in PPS
%
% Syntax
%
%     psh = convhull(P)
%     psh = convhull(P,'groups',c,'bufferwidth',bw);
%     convhull(P,...)
%
% Description
%
%     convhull returns the convex hull around all or groups of points in P
%     as polyshape object. Without output arguments, convhull plots the
%     convex hull(s).
%
% Input arguments
%
%     P      instance of PPS
%
%     Parameter name/value pairs
%
%     'groups'       vector of group assignments of each point (e.g. as
%                    returned by PPS/cluster.
%     'bufferwidth'  default is the maximum extent of the stream network
%                    divided by 100
%     'useparallel'  true (default) or false. 
%     'text'         false (default) or true (only applicable if plotted)
%     'textcolor'    color of text (if 'text', true)
%     'backgroundcolor background color of text (if 'text', true)
%     'onlynonzero'  encircle only groups with non-zero values 
%
%
% Output arguments
%
%     psh     polyshape object
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',500);
%     S   = klargestconncomps(S,1);
%     P   = PPS(S,'runif',200,'z',DEM);
%     c   = cluster(P,'cutoff',2000);
%     convhull(P,'groups',c,'bufferwidth',200)
%
%
% See also: PPS, PPS/cluster, PPS/aggregate
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2020

bbox = info(P.S,'boundingbox');
maxext  = max(bbox(2)-bbox(1),bbox(4)-bbox(3));

p = inputParser;
p.FunctionName = 'PPS/convhull';
addParameter(p,'groups',ones(npoints(P),1),@(x) (numel(x) == npoints(P)) | isa(x,'GRIDobj'));
addParameter(p,'bufferwidth',maxext/100,@(x) isscalar(x) && x>0);
addParameter(p,'useparallel',true);
addParameter(p,'plotall',true);
addParameter(p,'text',false);
addParameter(p,'textcolor','k');
addParameter(p,'backgroundcolor','w');
addParameter(p,'onlynonzero',false);
parse(p,varargin{:});

groups = p.Results.groups;
if isa(groups,'GRIDobj')
    groups = groups.Z(P.S.IXgrid(P.PP));
end

if p.Results.onlynonzero
    I = groups == 0;
    groups(I) = [];
    xy = P.ppxy(~I,:);
    
else
    xy = P.ppxy;
end

[v,~,ixb] = unique(groups);
n = numel(v);

% assemble points in cell array
ix = accumarray(ixb,(1:numel(ixb))',[n 1],@(x) {x});
xy = cellfun(@(ix) xy(ix,:),ix,'UniformOutput',false);


bw = p.Results.bufferwidth;

warning off
if p.Results.useparallel
    parfor r = 1:numel(xy)
        psh(r) = convhullsub(xy{r},bw);
    end       
else
    psh = cellfun(@(xy) convhullsub(xy,bw),xy,'UniformOutput',true);
end
warning on

if nargout == 0
    if p.Results.plotall
        tf = ishold;
        hold on
        plot(P);
    end
    plot(psh)
    
    if p.Results.text
        [xt,yt] = centroid(psh);
        hold on
        text(xt,yt,num2cell(v),'color',p.Results.textcolor,...
            'EdgeColor','none','Backgroundcolor',p.Results.backgroundcolor)
    end
    
    if ~tf
        hold off
    end
else
    h = psh;
end
end

function h = convhullsub(xy,bw)

% check whether points are unique
xy = unique(xy,'rows');

if size(xy,1) == 1
    h = polybuffer(xy(1,:),'points',bw);
elseif size(xy,1) == 2
    h = polybuffer(xy,'lines',bw);
else
    try
        % try, but may not work if points are collinear
        c = convhull(xy);
        poly = polyshape(xy(c,1),xy(c,2),...
            'SolidBoundaryOrientation','ccw',...
            'Simplify',false,...
            'KeepCollinearPoints',true);
        h = polybuffer(poly,bw);
    catch
        h = polybuffer(xy,'lines',bw);
    end
end
end
    
