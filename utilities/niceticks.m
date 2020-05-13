function niceticks(ax)

%NICETICKS Makes nice ticks in a 2D plot
%
% Syntax
%
%     niceticks
%     niceticks(ax)
%
% Description
%
%     Coordinate values are often quite large numbers which are displayed
%     with an exponent along the graphics axes. This function identifies
%     the best location for axis ticks and labels the first and last one
%     without exponents.
%
%
% See also: imageschs
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 13. May, 2020

if nargin == 0
    ax = gca;
end

xticks(ax,'auto')
yticks(ax,'auto')
xticklocs = get(ax,'XTick');
yticklocs = get(ax,'YTick');

set(ax,'XTick',xticklocs([1 end]))
set(ax,'YTick',yticklocs([1 end]))

set(ax,'XTickLabel',num2str(xticklocs([1 end])','%d'));
set(ax,'YTickLabel',num2str(yticklocs([1 end])','%d'));

% rotate tick labels if matlab 2014b or newer available
if ~verLessThan('matlab','8.4')
    set(ax,'YTickLabelRotation',90);
end