function varargout = plotdz(P,varargin)

%PLOTDZ plot upstream distance version elevation or covariate of a PPS
%
% Syntax
%
%     plotdz(P)
%
% Description
%
%     plot distance versus elevation of points and stream network in an
%     instance of PPS.
%
% Input arguments
%
%     P       instance of PPS (needs z-property unless parameter z is set
%             to different node-attribute list)
%
%     Parameter name/value pairs
%
%     For the river profile, all parameters that work for STREAMobj/plotdz 
%     can be used except 'color'. Use 'LineColor' instead to define the
%     color of the line. If the line should not be visible, use
%     'LineColor','none'.
%
%     For the points, all parameters that work for scatter or bubblechart 
%     if you have MATLAB 2020b or higher. Note that using bubblechart does
%     not allow you to set the parameter 'marker'.
%
%     For versions older than 2020b, the function uses scatter to plot
%     points. If you supply SizeData, then you can additionally control the
%     range of sizes used for plotting by following arguments.
%
%     'MinSize'    see 'SizeData' (default = 5)
%     'MaxSize'    see 'SizeData' (default = 75)
%
%     Additional arguments are
%
%     'z'         default is the node-attribute of elevation values stored 
%                 in the PPS object. Takes any other node-attribute list.
%     'SizeData'  data (GRIDobj, node-attribute list, or attributes (marks) 
%                 for each point.
%     'ColorData' data (GRIDobj, node-attribute list, or attributes (marks) 
%                 for each point. 
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     P = PPS(S,'runif',30,'z',DEM);
%     plotdz(P,'SizeData',P.S.distance,'ColorData',DEM)
%
% See also: PPS, PPS/npoints, bubblechart, bubblelegend, bubblesize 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 28. September, 2020

% Check version because 2020b onwards will use bubblechart
olderthan2020b = verLessThan('matlab','9.9');
if olderthan2020b
    minsi = 5;
    maxsi = 75;
else
    minsi = 1;
    maxsi = 10;
end

%% Input parameter parsing
p = inputParser;
p.KeepUnmatched = true;

% Custom river elevations and distance
addParameter(p,'z','z');
addParameter(p,'distance',[],@(x) isnal(P.S,x));
addParameter(p,'LineColor',[.7 .7 .7])

% Point properties
if olderthan2020b
addParameter(p,'Marker','o');
end
addParameter(p,'MarkerEdgeColor','k');
addParameter(p,'MarkerFaceAlpha',0.5);
addParameter(p,'MinSize',minsi);
addParameter(p,'MaxSize',maxsi);

addParameter(p,'SizeData',20);    
addParameter(p,'ColorData','w');

% Parse
parse(p,varargin{:});

UnMatched = p.Unmatched;
UnMatched.distance = p.Results.distance;

tf = ishold;

% Get line color
lc = p.Results.LineColor;

% Get elevations
zz  = getcovariate(P,p.Results.z);

if isempty(p.Results.distance)
    d  = P.S.distance;
else
    d  = p.Results.distance;
end

%% Handle additional arguments
Results = p.Results;
sz      = p.Results.SizeData;

% -- SizeData
if isnal(P.S,sz) || isa(sz,'GRIDobj')
    sz = getmarks(P,sz);
elseif numel(sz) == npoints(P)
    % great, go on
end

% If scatter is used, we will adjust the size range of the data between the
% arguments MinSize and MaxSize
if olderthan2020b && ~isscalar(sz)
    % normalize sz between 0 and 1
    sz = sz - min(sz);
    sz = sz./max(sz);
    sz = sz*(Results.MaxSize - Results.MinSize) + Results.MinSize;
else
    maxsz = Results.MaxSize;
    minsz = Results.MinSize;
end

% -- ColorData
cols    = Results.ColorData;
if isnal(P.S,cols) || isa(cols,'GRIDobj')
    cols = getmarks(P,cols);
end

% -- Remove additional arguments from the Results structure
Results = rmfield(Results,{'distance' 'z' 'SizeData' 'ColorData' ...
                           'MaxSize' 'MinSize' 'LineColor'});

% Does the field 'dunit' exist
if isfield(UnMatched,'dunit')
    switch UnMatched.dunit
        case 'km'
            d = d/1000;
    end
end
pnpv    = expandstruct(Results);

%% Plotting
% Plot line
hl = plotdz(P.S,zz,UnMatched,'Color',lc);
hold on

if olderthan2020b || isscalar(sz) 
    hs = scatter(d(P.PP),zz(P.PP),sz,cols,'filled',pnpv{:});
else
    hs = bubblechart(d(P.PP),zz(P.PP),sz,cols,pnpv{:});
    bubblesize([minsz maxsz])
end

if ~tf
    hold off
end

if nargout >= 1
    varargout{1} = hl;
end
if nargout == 2
    varargout{2} = hs;
end
end

function pnpv = expandstruct(s)

pn = fieldnames(s);
pv = struct2cell(s);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';
end
