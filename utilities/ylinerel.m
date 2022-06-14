function ht = ylinerel(values,rel,varargin)

%YLINEREL Plot vertical lines with constant y values and relative length
%
% Syntax
%
%     ylinerel(y)
%     ylinerel(y, fraction, pn, pv)
%
% Description
%
%     ylinerel plots vertical lines at the values y with length (or
%     positions) relative to the limits of the x-axis. 
%
% Input arguments
%
%     y          values of y
%     fraction   fraction of x-axis limits. If fraction is a scalar, then
%                the function plots the lines starting on the x-axis. If
%                you specify fraction as two-element vector with values
%                between 0 and 1, then the lines will be plotted within
%                these fractions of the x-axis. For example, [0.95 1] will
%                plot lines in the right 5% of the axis. By default,
%                fraction is 0.05.
%     pn,pv      all arguments accepted by the plot function.
%
% Output arguments
%
%     h          handle to line function
%     
% Example
%
%     plot([0 16],[0 5],'k')
%     hold on
%     ylinerel([0:pi/2:5*pi],[.95 1],'r','LineWidth',2);
%
% See also: xline, yline, xlinerel
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 26. August, 2021

if nargin == 1
    rel = 0.05;
end
if isempty(rel)
    rel = 0.05;
end

if isscalar(rel)
    rel = [0; rel];
else
    rel = sort(rel,'ascend');
end

if any(rel > 1 | rel < 0)
    error('The fraction must be within the range of 0 and 1')
end

ax = gca;
xl = xlim(ax);

n  = numel(values);
y  = values(:)';
y  = [y;y];

xlow  = xl(1) + (xl(2)-xl(1))*rel(1);
xhigh = xl(1) + (xl(2)-xl(1))*rel(2);

x  = repmat([xlow; xhigh],1,n);

y  = [y;nan(1,n)];
x  = [x;nan(1,n)];

h  = plot(x(:),y(:),varargin{:});
xlim(ax,xl);

hlist = addlistener(ax,'XLim','PostSet',@(src,evnt)updatedata);
h.DeleteFcn = @(src,evt)upondelete(src,evt,hlist);

if nargout == 1
    ht = h;
end

function updatedata
    xl = xlim(ax);
    xlow  = xl(1) + (xl(2)-xl(1))*rel(1);
    xhigh = xl(1) + (xl(2)-xl(1))*rel(2);
    h.XData(1:3:end) = xlow;
    h.XData(2:3:end) = xhigh;

end
function upondelete(scr,evt,hlist)

        delete(hlist)

    end
end

