function h = plotdz(S,DEM,varargin)

% plot upstream distance version elevation of a stream network
%
% Syntax
%
%     plotdz(S,DEM)
%     plotdz(S,DEM,pn,pv,...)
%     plotdz(S,nda,pn,pv,...)
%     h = ...
%
% Description
%
%     Plot stream distance from the outlet versus elevation. 
%
% Input arguments
%
%     S      instance of STREAMobj
%     DEM    digital elevation model (GRIDobj)
%     nda    node attribute vector (as returned by various STREAMobj
%            methods, e.g. STREAMobj/streamorder, STREAMobj/gradient)
%
%     Parameter name/value pairs {default}
%
%     'annotation':  {[]} ix      
%     vector with linear indices of locations into the DEM. The cells
%     referenced by ix must be part of the stream network. Use snap2stream
%     to adjust locations. Annotation is achieved with vertical arrows.
%
%     'annotationtext': cell array of strings
%     if annotated, a cell array of strings can be added to the vertical 
%     arrows
%
%     'color': {'b'}
%     line colors as provided to plot
%
%     'dunit': {'m'} 'km'
%     distance unit. plotdz assumes that distance is given in meter. 
%
%     'doffset': {0}
%     add an offset (scalar) to the distance from outlet to shift the 
%     x-axis of the plot.
%
%     'gradient': {false},true
%     set to true, if you want to plot the gradient instead of elevation.
%
%     'kernelsize'  {11}, odd, scalar integer
%     if option 'smooth' is true, than kernelwidth determines the width of
%     the moving average filter. kernelsize must be odd.
%
%     'smooth'   {false} or true
%     smooth the profiles using a moving average filter (see also the
%     kernelwidth parameter).
%
%
% Output arguments
%
%     h     handle to the line handle
%
% Example
%
% 
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 2. May, 2013

nrnodes = numel(S.x);

% check input
p = inputParser;         
p.FunctionName = 'plotdz';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'DEM', @(x) isa(x,'GRIDobj') || numel(x) == nrnodes);

addParamValue(p,'annotation',[])
addParamValue(p,'color','b');
addParamValue(p,'annotationtext',{});
addParamValue(p,'gradient',false,@(x) isscalar(x));
addParamValue(p,'smooth',false,@(x) isscalar(x));
addParamValue(p,'dunit','m',@(x) ischar(validatestring(x,{'m' 'km'})));
addParamValue(p,'doffset',0,@(x) isscalar(x));
addParamValue(p,'kernelsize',11,@(x) isscalar(x) && mod(x,2)==1 && x>=3);

parse(p,S,DEM,varargin{:});
S   = p.Results.S;
DEM = p.Results.DEM;

if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    
end

% get dynamic properties of S
order    = S.orderednanlist;
distance = S.distance;

switch lower(p.Results.dunit);
    case 'km'
        distance = distance/1000;
end
        

% apply distance offset
distance = distance + p.Results.doffset;

if isa(DEM,'GRIDobj')
    zz    = DEM.Z(S.IXgrid);
else
    zz    = DEM;
    if numel(zz) ~= numel(S.IXgrid);
        error('TopoToolbox:wronginput',...
            ['The number of elements in the node attribute vector \n' ...
             'must equal the number of nodes in the stream network']);
    end
        
end

% smooth
if ~p.Results.smooth;
    I     = ~isnan(order);
    d     = nan(size(order));
    d(I)  = distance(order(I));
    z     = nan(size(order));
    
    if ~p.Results.gradient
        
        z(I)  = zz(order(I));
    else
        
        g     = (zz(S.ix)-zz(S.ixc))./hypot(S.x(S.ix)-S.x(S.ixc),S.y(S.ix)-S.y(S.ixc));
        gg     = zeros(size(S.x));
        gg(S.ix) = g;
        z     = nan(size(order));
        z(I)  = gg(order(I));
    end
        
else
    if ~p.Results.gradient
        [d,z] = smooth(distance,double(zz),order,p.Results.kernelsize);
    else
        g     = (zz(S.ix)-zz(S.ixc))./hypot(S.x(S.ix)-S.x(S.ixc),S.y(S.ix)-S.y(S.ixc));
        gg     = zeros(size(S.x));
        gg(S.ix) = g;
        [d,z] = smooth(distance,double(gg),order,p.Results.kernelsize);
    end
        
end

% plot
ht = plot(d,z,'-','Color',p.Results.color);  
xlabel(['distance upstream [' lower(p.Results.dunit) ']'])
ylabel('elevation [m]')

%% Annotation
if ~isempty(p.Results.annotation);
    ix = p.Results.annotation;
    hold on
    [Lia,Locb] = ismember(ix,S.IXgrid);


    if any(~Lia)
        error('TopoToolbox:STREAMobj',...
            'Some of the annotations are not located on the stream network')
    end
    
    annd = distance(Locb);
    
    if isa(DEM,'GRIDobj')
        annz = DEM.Z(S.IXgrid(Locb));
    else
        annz = zz(Locb);
    end
    
    if ~isempty(p.Results.annotationtext);
        c = p.Results.annotationtext;
        addtext = true;
    else
        addtext = false;
    end
    
    
    for r = 1:numel(ix);
        if addtext
            annotext = [c{r} '\newline\downarrow'];
        else
            annotext = '\downarrow ';
        end
        
        text('Position',[annd(r), annz(r)],...
             'String', annotext,...
             'VerticalAlignment','bottom',...
             'FontWeight','bold');
    end
    hold off
end


if nargout == 1;
    h = ht;
end
end

function [ds,zs] = smooth(d,z,order,ks)
% smooth profiles using an moving average filter
%
% start and endpoints of stream segments are retained or adapted to streams
% that have been filtered before
%

INNAN = ~isnan(order);
nanix1 = 0;
nanix2  = find(~INNAN);
ds     = nan(size(order));
ds(INNAN) = d(order(INNAN));
zs     = nan(size(order));

kernel   = ones(ks,1)/ks;
padwidth = floor(ks/2);

for r = 1:numel(nanix2);
    dd = ds(nanix1+1:nanix2(r)-1);
    zz = z(order(nanix1+1:nanix2(r)-1));
    
    di = linspace(dd(1),dd(end),numel(dd));
    zi = interp1(dd,zz,di,'linear');
    
    zi = [repmat(zi(1),padwidth,1);...
          zi(:);...
          repmat(zi(end),padwidth,1)];
    
    zi = conv(zi,kernel,'valid');
    zh = interp1(di,zi,dd);
    zh(1) = zz(1);
    zh(end) = zz(end);
    
    z(order(nanix1+1:nanix2(r)-1)) = zh;
    zs(nanix1+1:nanix2(r)-1) = zh;
    
    nanix1 = nanix2(r);
    
end
end
