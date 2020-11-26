function h = plotpoints(P,varargin)

%PLOTPOINTS plot points of PPS
%
% Syntax
%
%     plotpoints(P)
%     plotpoints(P,'pn',pv,...)
%     h = plotpoints(...)
%
% Description
%
%     plotpoints plots the points of an instance of PPS. The function uses
%     the build-in function scatter. Thus, it accepts all parameter
%     name/value pairs that are accepted by scatter.
%     
% Input arguments
%
%     P    instance of PPS
%
%     Parameter name/value pairs (see scatter)
%
%     'Marker'     Marker type (see scatter)
%     'MarkerEdgeColor'  Marker edge color (see scatter) (default = 'k')
%     'MarkerFaceColor'  Marker face color (see scatter) (default = 'w')
%     'LineWidth'  Width of the marker edge
%     'SizeData'   Default=20. If scalar, then markers will be plotted
%                  using this size (see scatter). If vector, then marker 
%                  size will vary according to the values in this vector. 
%                  Marker size will be automatically scaled to vary within
%                  the bounds set by 'MinSize' and 'MaxSize' if used in
%                  versions of MATLAB older than 2020b.
%     'MinSize'    see 'SizeData' (default = 5)
%     'MaxSize'    see 'SizeData' (default = 75)
%     'MarkerFaceAlpha'  Transparency (default = 0.5). 0 is completely 
%                  transparent and 1 is fully opaque. (see scatter)
%     'ColorData'
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = removeshortstreams(S,100);
%     S = clean(S);
%     P = PPS(S,'rpois',0.001,'z',DEM);
%     plot(S,'k')
%     hold on
%     plotpoints(P,'colordata',DEM,'sizedata',distance(S,'max'))
%
%
% See also: PPS/plot, PPS/plotdz, STREAMobj/plot, STREAMobj/plotc
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 26. November, 2020


p = inputParser;
p.KeepUnmatched = true;
if verLessThan('matlab','9.9')
addParameter(p,'Marker','o');
end
addParameter(p,'MarkerEdgeColor','k');
addParameter(p,'MarkerEdgeAlpha',0.5);
addParameter(p,'MarkerFaceAlpha',0.5);

addParameter(p,'LineWidth',1)
addParameter(p,'SizeData',20);
addParameter(p,'ColorData','w');
addParameter(p,'MinSize',5);
addParameter(p,'MaxSize',75);

% Parse
parse(p,varargin{:});

Results = p.Results;

cols = Results.ColorData;
% Get colordata
if isnal(P.S,cols) || isa(cols,'GRIDobj')
    cols = getmarks(P,cols);
end


if isscalar(Results.SizeData) && ~isa(Results.SizeData,'GRIDobj')
    % all good. All points have the same size
else
    if numel(Results.SizeData) == npoints(P)
    else
        Results.SizeData = getmarks(P,Results.SizeData);
    end
end

if verLessThan('matlab','9.9')
    % scale between 1 and 20
    Results.SizeData = Results.SizeData - min(Results.SizeData);
    Results.SizeData = Results.SizeData./max(Results.SizeData);
    Results.SizeData = Results.SizeData*(Results.MaxSize - Results.MinSize) + Results.MinSize;
else
    sz = Results.SizeData;
    Results = rmfield(Results,'SizeData');
end
        
Results = rmfield(Results,{'MaxSize','MinSize','ColorData'});

xy = P.ppxy;

pnpv = expandstruct(Results);

if verLessThan('matlab','9.9') || isscalar(sz)
    % Use scatter
    ht = scatter(xy(:,1),xy(:,2),[],cols,'filled',pnpv{:});
else
    ht = bubblechart(xy(:,1),xy(:,2),sz,cols,pnpv{:});
    bubblesize([1 10])
end
    

if nargout == 1
    h = ht;
end

end

function pnpv = expandstruct(s)

pn = fieldnames(s);
pv = struct2cell(s);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';
end