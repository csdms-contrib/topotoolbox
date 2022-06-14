function niceticks(ax,varargin)

%NICETICKS Makes nice ticks in a 2D plot
%
% Syntax
%
%     niceticks
%     niceticks(ax)
%     niceticks(ax,pn,pv,...)     
%
% Description
%
%     Coordinate values are often quite large numbers which are displayed
%     with an exponent along the graphics axes. This function identifies
%     the best location for axis ticks and labels the first and last one
%     without exponents.
%
% Input arguments
%
%     ax    axis handle (e.g. gca)
%    
%     Parameter name/value pairs
%
%     'degree'      {false} or true. If true, tick labels will be appended
%                   with °E and °N. Note that currently, °W is not
%                   supported.
%     'twoticks'    {true} or false. If true, two ticks per axis.
%     'precision'   Numbers behind the decimal. Default = 0
%     'exponent'    Exponent. Default = 0%
%
% See also: imageschs
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 13. May, 2020

if nargin == 0
    ax = gca;
end

p = inputParser;
p.FunctionName = 'niceticks';
addParameter(p,'degree',false)
addParameter(p,'twoticks',true)
addParameter(p,'precision',0)
addParameter(p,'exponent',0)
addParameter(p,'rotateylabel',true)
parse(p,varargin{:})

prec = ['% .' num2str(p.Results.precision) 'f'];

xticks(ax,'auto')
yticks(ax,'auto')
xticklocs = get(ax,'XTick');
yticklocs = get(ax,'YTick');

nx = numel(xticklocs);
ny = numel(yticklocs);

% xticks
ix = find(round(xticklocs,p.Results.precision) == xticklocs);
if p.Results.twoticks
    if nnz(ix) >= 2    
        xticklocs = xticklocs(ix([1 end]));
    end
else
    if any(ix)
    xticklocs = xticklocs(ix);
    end
end
set(ax,'XTick',xticklocs)

% yticks
ix = find(round(yticklocs,p.Results.precision) == yticklocs);
if p.Results.twoticks
    if nnz(ix) >= 2    
        yticklocs = yticklocs(ix([1 end]));
    end
else
    if any(ix)
    yticklocs = yticklocs(ix);
    end
end
set(ax,'YTick',yticklocs)

if p.Results.degree
    xtickformat([prec '°E'])
    ytickformat([prec '°N'])
else
    xtickformat(prec)
    ytickformat(prec)
end

ax.XAxis.Exponent = p.Results.exponent;
ax.YAxis.Exponent = p.Results.exponent;

% rotate tick labels if matlab 2014b or newer available
if ~verLessThan('matlab','8.4') && p.Results.rotateylabel
    set(ax,'YTickLabelRotation',90);
end

% if p.Results.degree
%     setlabelstodeg(ax)
% end
end

function setlabelstodeg(ax)
xt = xticklabels;
xt = cellfun(@(x) [x '°E'],xt,'UniformOutput',false);
set(ax,'XTickLabel',xt)

yt = yticklabels;
yt = cellfun(@(y) [y '°N'],yt,'UniformOutput',false);
set(ax,'YTickLabel',yt)
end