function [st,hs] = plotstreamorder(S,varargin)

%PLOTSTREAMORDER calculate stream order from STREAMobj
%
% Syntax
%
%     [h,hs] = streamorder(S)
%     [h,hs] = streamorder(S,pn,pv,...)
%
% Description
%
%     The Strahler Stream Order is a way to classify rivers based on a
%     hierarchy of tributaries. First-order streams don't have tributaries.
%     When two first-order streams come together, they form a second-order
%     stream. When two second-order streams come together, they form a
%     third-order stream. When two streams of different order confluence,
%     they form a stream of maximum order of both. 
%
%     [h,hs] = plotstreamorder(S)   plots the x- and y-coordinates of the
%     stream network and colors the stream sections according to their
%     stream order and returns a vector of handles to lineseries objects 
%     h and their streamorders s.
%
%     [h,hs] = streamorder(S,'plot',pn,pv,...)   lets you define various
%     parameter name/value pairs (see below).
%     
%
% Input
%   
%     S         STREAMobj 
%
%     parameter name/value pairs (for plotting only) {default}
%     
%     type        {'strahler'} or 'shreve'
%     legend      {true} plots a legend, false doesn't
%     colormap    provide string with name of colormap to be used for 
%                 coloring streams according to streamorder (e.g. {'jet'}, 
%                 'gray', 'hsv', or any other colormap available). If only
%                 one color should be used for all stream orders, provide a
%                 vector with rgb values (e.g. [0 0 0] to plot all lines in
%                 black).
%     parent      handle to parent axis {gca}
%     linewidth   scalar or vector of line widths ({max(so/2,1)} where so
%                 is stream order)
%    
%     
%
% Output
%
%     h         vector of handles to lineseries objects
%     hs        stream order of each line handle
%
% Example:
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     plotstreamorder(S,'colormap',[0 0 0],...
%                       'LineWidth',max([1 2 3 4 5]/2,1));
% 
% See also: STREAMobj, FLOWobj/streamorder, STREAMobj/streamorder
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

p = inputParser;
p.FunctionName = 'STREAMobj/plotstreamorder';
addParamValue(p,'type','strahler',@(x) ischar(validatestring(x,{'strahler','shreve'})));
addParamValue(p,'legend',true,@(x) isscalar(x));
addParamValue(p,'colormap','jet', @(x) ischar(x) || (isnumeric(x) && numel(x)==3));
addParamValue(p,'parent',[],@(x) ishandle(x));
addParamValue(p,'linewidth',[]);

parse(p,varargin{:});

% create mapstruct
MS  = STREAMobj2mapstruct(S,p.Results.type);
un  = unique([MS.streamorder],'sorted');
nrs = max(un);

% get colormap
if ischar(p.Results.colormap)
    cmap = str2func(p.Results.colormap);
    cmap = cmap(nrs);
else
    cmap = repmat(p.Results.colormap(:)',nrs,1);
end

% get linewidth
if isempty(p.Results.linewidth);
    switch p.Results.type
        case 'strahler'
            lw = max((1:nrs)/2,1);
        case 'shreve'
            lw = min(((1:nrs)/nrs)*2,0.5);
    end
else
    lw = p.Results.linewidth;
end

if isempty(p.Results.parent)
    ax = gca;
else
    ax = p.Results.parent;
end
t  = ishold;
fh = zeros(numel(un),1);

s  = [];
hs = 1:nrs;

for r = 1:numel(un);
    so = un(r);
    I = [MS.streamorder]==so;
    h = plot(ax,[MS(I).X],[MS(I).Y],'Color',cmap(so,:),'linewidth',lw(min(r,numel(lw))));
    hold on
    fh(r) = h(1);
    if nargout > 0
        s = [s h];
    end
    
end

if ~t
    hold off
end


if p.Results.legend
    legnames = cellfun(@(x) num2str(x),num2cell(un),'uniformoutput',false);
    legend(fh,legnames);
end


if nargout > 0
    st  = s;
end