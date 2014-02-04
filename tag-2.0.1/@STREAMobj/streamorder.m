function [st,hs] = streamorder(S,pp,varargin)

% calculate Strahler Stream Order from STREAMobj
%
% Syntax
%
%     s = streamorder(S)
%     streamorder(S)
%     [h,hs] = streamorder(S,'plot')
%     [h,hs] = streamorder(S,'plot',pn,pv,...)
%
% Description
%
%     The Strahler Stream Order is a way to classify rivers based on a
%     hierarchy of tributaries. First-order streams don't have tributaries.
%     When two first-order streams come together, they form a second-order
%     stream. When two second-order streams come together, they form a
%     third-order stream. When two streams of different order confluence,
%     they form a stream of maximum order of both. streamorder returns the
%     Strahler Stream Order based on an instance of the STREAMobj (S). 
%
%     s = streamorder(S)     returns a column vector with node attributes 
%     of the directed graph represented by STREAMobj.
%
%     streamorder(S)     plots the x- and y-coordinates of the stream
%     network and colors the stream sections according to their stream
%     order.
%
%     [h,hs] = streamorder(S,'plot')   plots the x- and y-coordinates of the
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
%     s         stream order for each node in STREAMobj
%     h         vector of handles to lineseries objects
%     hs        stream order of each line handle
%
% Example:
%
%
% 
% See also: STREAMobj, FLOWobj/streamorder
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. January, 2013


if nargout == 1 && nargin == 1;
    nrc = numel(S.x);
    s   = zeros(nrc,1);

    for r = 1:numel(S.ix);
        if s(S.ix(r)) == 0;
            s(S.ix(r)) = 1;
        end

        if s(S.ix(r)) < s(S.ixc(r))
        elseif s(S.ix(r)) == s(S.ixc(r))
            s(S.ixc(r)) = s(S.ixc(r))+1;
        else
            s(S.ixc(r)) = s(S.ix(r));
        end
    end
    

else
    
    p = inputParser;
    p.FunctionName = 'STREAMobj/streamorder';
    addParamValue(p,'legend',true,@(x) isscalar(x));
    addParamValue(p,'colormap','jet', @(x) ischar(x) || (isnumeric(x) && numel(x)==3));
    addParamValue(p,'parent',[],@(x) ishandle(x));
    addParamValue(p,'linewidth',[]);
    
    parse(p,varargin{:});
    % streamorder will be plotted
    MS  = STREAMobj2mapstruct(S);
    nrs = max([MS.streamorder]);
    
    % get colormap
    if ischar(p.Results.colormap)
        cmap = str2func(p.Results.colormap);
        cmap = cmap(nrs);
    else
        cmap = repmat(p.Results.colormap(:)',nrs,1);
    end
    
    % get linewidth
    if isempty(p.Results.linewidth);
        lw = max((1:nrs)/2,1);
    else
        lw = p.Results.linewidth;
    end
    
    if isempty(p.Results.parent)
        ax = gca;
    else
        ax = p.Results.parent;
    end
    t  = ishold;
    fh = zeros(nrs,1);
    
    s  = [];
    hs = 1:nrs;
    
    for r = 1:nrs
        I = [MS.streamorder]==r;
        h = plot(ax,[MS(I).X],[MS(I).Y],'Color',cmap(r,:),'linewidth',lw(min(r,numel(lw))));
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
        legnames = cellfun(@(x) num2str(x),num2cell(1:nrs),'uniformoutput',false);
        legend(fh,legnames);
    end
end

if nargout > 0
    st  = s;         
end    
    
    
