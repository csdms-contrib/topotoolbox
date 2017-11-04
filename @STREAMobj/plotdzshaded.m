function h = plotdzshaded(S,zz,varargin)

%PLOTDZSHADED plot upstream distance version elevation of a stream network
%
% Syntax
%
%     plotdzshaded(S,[nal1 nal2])
%     plotdzshaded(S,[nal1 nal2],pn,pv,...)
%     h = ...
%
% Description
%
%     Plot stream distance from the outlet versus elevation. PLOTDZSHADED
%     plots a filled area between the higher and lower values for each
%     element contained in two node-attribute lists.
%
% Input arguments
%
%     S      instance of STREAMobj
%     [nal1 nal2]    matrix of two node attribute list (as returned by 
%            various STREAMobj functions). plotdzshaded will plot a filled
%            area between the high and low values in nal1 and nal2,
%            respectively.
%
%     Parameter name/value pairs {default}
%
%     'distance': {S.distance}
%     node attribute list with custom distances (see STREAMobj/distance) or
%     STREAMobj (see function distance(S,S2))
%
%     'dunit': {'m'} 'km'
%     distance unit. plotdz assumes that distance is given in meter. 
%
%     'doffset': {0}
%     add an offset (scalar) to the distance from outlet to shift the 
%     x-axis of the plot.
% 
%     'FaceColor'          see patch properties reference page 
%     'FaceAlpha' {0.3}    see patch properties reference page
%     'EdgeColor' {none}   see patch properties reference page
%     'EdgeAlpha' {1}      see patch properties reference page
%     'LineStyle' {'-'}    see patch properties reference page
%     'LineWidth' {.5}     see patch properties reference page
%
% Output arguments
%
%     h     patch handles.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'mex',true,'preprocess','carve');
%     S  = STREAMobj(FD,flowacc(FD)>1000);
%     S  = klargestconncomps(S);
%     z1  = getnal(S,DEM);
%     z2  = z1+50;
%     plotdzshaded(S,[z1 z2])
%
% See also: STREAMobj, STREAMobj/plotdz
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. March, 2017

nrnodes = numel(S.x);
ax      = gca;
colororderindex = mod(ax.ColorOrderIndex, size(ax.ColorOrder,1));
if colororderindex==0; colororderindex=size(ax.ColorOrder,1); end

% check input
p = inputParser;         
p.FunctionName = 'plotdzshaded';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'zz', @(x) isequal(size(x),[nrnodes 2]));
addParameter(p,'distance',[],@(x) isnal(S,x) || isa(x,'STREAMobj') || ischar(x));
addParameter(p,'dunit','m',@(x) ischar(validatestring(x,{'m' 'km'})));
addParameter(p,'doffset',0,@(x) isscalar(x));

% Faces
addParameter(p,'FaceColor',ax.ColorOrder(colororderindex,:));
addParameter(p,'FaceAlpha',0.3);
% Edges
addParameter(p,'EdgeColor','none');
addParameter(p,'EdgeAlpha',1);
addParameter(p,'LineStyle','-');
addParameter(p,'LineWidth',0.5);


parse(p,S,zz,varargin{:});
S   = p.Results.S;

if isempty(p.Results.distance);
    dist = S.distance;
else
    if isa(p.Results.distance,'STREAMobj');
        dist = distance(S,p.Results.distance);
    elseif ischar(p.Results.distance)
        dist = S.distance;
        dist = dist./max(dist);
    else
        dist = p.Results.distance;
    end
end



switch lower(p.Results.dunit)
    case 'km'
        dist = dist/1000;
end
        

% apply distance offset
dist = dist + p.Results.doffset;

%% get tributaries
[Ctribs,locS] = STREAMobj2cell(S,'tributaries');

%% Plot
keephold = ishold;
if ~keephold
    cla
end

ph = [];
for r = 1:numel(Ctribs)
    ix       = locS{r};
    order    = Ctribs{r}.orderednanlist;
    disttemp = dist(ix);
    zztemp   = zz(ix,:);
    
    I     = ~isnan(order);
    order = order(I);
    
    d     = disttemp(order);
    z1    = zztemp(order,1);
    z2    = zztemp(order,2);

    
    ph = [ph patch([d;d(end:-1:1)],[z1;z2(end:-1:1)],...
                p.Results.FaceColor,...
                'FaceAlpha',p.Results.FaceAlpha,...
                'EdgeColor',p.Results.EdgeColor,...
                'LineStyle',p.Results.LineStyle,...
                'EdgeAlpha',p.Results.EdgeAlpha,...
                'LineWidth',p.Results.LineWidth)];
    hold on
    
end

if ~keephold
    hold off
end

xlabel(['Distance upstream [' lower(p.Results.dunit) ']'])
ylabel('Elevation [m]')


if nargout == 1;
    h = ph;
end

