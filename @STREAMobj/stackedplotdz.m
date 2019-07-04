function s = stackedplotdz(S,C,varargin)

%STACKEDPLOTDZ plot several stacked variables along the stream networks
%
% Syntax
%
%     s = stackedplotdz(S,C)
%     s = stackedplotdz(S,C,'ylabels',ylabs)
%     s = stackedplotdz(...,pn,pv)
%
% Description
%
%     stackedplotdz takes several GRIDobjs or node-attribute lists (nals) and
%     plots them against the horizontal flow distance of the river network
%     S. C is a cell array that contains the GRIDobjs or nals.
%
%     stackedplotdz relies on the function stackedplot introduced in 2018b.
%     If you are using an older version of MATLAB, stackedplotdz will use
%     a vertical array of subplots.
%
% Input arguments
%
%     S        STREAMobj
%     C        cell array of GRIDobjs or node-attribute lists
%     ylabels  cell array of y-axis labels. The number of elements must
%              be the same as in C
%     pn,pv    additional parameter/value pairs as in the function
%              stackedplot (only for version 2018b)
%
% Output arguments
%
%     s        handle to stackedplot object
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     DEM = imposemin(FD,DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = klargestconncomps(trunk(S));
%     g   = gradient(S,DEM);
%     g   = smooth(S,g,'K',10);
%     stackedplotdz(S,{DEM,g},'ylabels',{'Elevation' 'Gradient (smoothed)'});
%
% See also: STREAMobj/plotdz
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. January, 2019



nals = cell(1,numel(C));
[~,~,d,nals{:}] = STREAMobj2XY(S,S.distance,C{:});

% go through varargin to find 'ylabel'
TF = strcmpi('ylabels',varargin);
ix = find(TF);
if ~isempty(ix)
    ylabels = varargin{ix+1};
    varargin([ix ix+1]) = [];
else
    ylabels = [];
end

% go through varargin to find 'dunit'
TF = strcmpi('dunit',varargin);
ix = find(TF);
if ~isempty(ix)
    dunit = varargin{ix+1};
    varargin([ix ix+1]) = [];
else
    dunit = 'm';
end

switch dunit
    case 'km'
        d = d/1000;
end

if ~verLessThan('matlab','9.5')
h = stackedplot(d,horzcat(nals{:}),varargin{:});
if ~isempty(ylabels)
    h.DisplayLabels = ylabels;
end
else
    n = numel(nals);
    for r = 1:n
        ax(r) = subplot(n,1,r);
        h(r)  = plot(d,nals{r});
    end
    linkaxes(ax,'x');
    if ~isempty(ylabels)
        ylabel(ylabels{r})
    end
        
end

if nargout > 0
    s = h;
end