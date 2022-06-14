function varargout = plot(P,varargin)

%PLOT plot instance of PPS
%
% Syntax
%
%     plot(P)
%     plot(P,'lineopts',{pn,pv,...},'markeropts',{pn,pv,...})
%     [hl,hm] = plot(...)
%
% Description
%
%     PPS/plot displays the stream network together with the events.
%     
% Input arguments
%
%     P    instance of PPS
%
%     Parameter name/value pairs
%
%     'lineopts'    cell array of parameter name/value pairs accepted by 
%                   plot (e.g. {'Color','k','LineWidth',2})
%     'markeropts'  cell array of parameter name/value pairs accepted for
%                   plotting markers (e.g. {'MarkerSize',20,'Color','b'})
%
% Example
%
%
%
%
% See also: PPS, PPS/plotdz, STREAMobj/plot, STREAMobj/plotc
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019


defaultlineopts = {'LineWidth' 1 ...
    'Color', [.3 .3 .3]};
defaultmarkeropts = {'Marker', '.' ...
    'Color', 'k' ...
    'MarkerFaceColor', 'k',...
    'MarkerSize', 10};

p = inputParser;
p.FunctionName = 'PPS/plot';
addParameter(p,'lineopts',defaultlineopts,@(x) iscell(x));
addParameter(p,'markeropts',defaultmarkeropts,@(x) iscell(x));
parse(p,varargin{:});

pl = inputParser;
pl.KeepUnmatched = true;
addParameter(pl,'MarkerSize', 10);
addParameter(pl,'LineWidth',1);
addParameter(pl,'LineStyle','-');
addParameter(pl,'Marker','none');
addParameter(pl,'Color',[.3 .3 .3]);
parse(pl,p.Results.lineopts{:});

pm = inputParser;
pm.KeepUnmatched = true;
addParameter(pm,'MarkerSize', 10);
addParameter(pm,'LineStyle','none');
addParameter(pm,'Marker','.');
addParameter(pm,'Color','k');
parse(pm,p.Results.markeropts{:});

tf = ishold;
hl = plot(P.S,pl.Results,pl.Unmatched);
hold on
xy = P.ppxy;

if ~isempty(xy)
    hm = plot(xy(:,1),xy(:,2),pm.Results);
else
    hm = [];
end

if ~tf
    hold off
end

if nargout >= 1
    varargout{1} = hl;
end
if nargout == 2
    varargout{2} = hm;
end
end