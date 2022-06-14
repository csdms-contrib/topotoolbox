function ht = xlinerel(values,rel,varargin)

%XLINEREL Plot vertical lines with constant x values and relative length
%
% Syntax
%
%     xlinerel(x)
%     xlinerel(x, fraction, pn, pv)
%
% Description
%
%     xlinerel plots vertical lines at the values x with length (or
%     positions) relative to the limits of the y-axis. 
%
% Input arguments
%
%     x          values of x
%     fraction   fraction of y-axis limits. If fraction is a scalar, then
%                the function plots the lines starting on the y-axis. If
%                you specify fraction as two-element vector with values
%                between 0 and 1, then the lines will be plotted within
%                these fractions of the y-axis. For example, [0.95 1] will
%                plot lines in the upper 5% of the axis. By default,
%                fraction is 0.05.
%     pn,pv      all arguments accepted by the plot function.
%
% Output arguments
%
%     h          handle to line function
%     
% Example
%
%     plot(linspace(0,5*pi),sin(linspace(0,5*pi)))
%     hold on
%     xlinerel([0:pi/2:5*pi],[.95 1],'r','LineWidth',2);
%
% See also: xline, yline, ylinerel
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
yl = ylim(ax);

n  = numel(values);
x  = values(:)';
x  = [x;x];

ylow  = yl(1) + (yl(2)-yl(1))*rel(1);
yhigh = yl(1) + (yl(2)-yl(1))*rel(2);

y  = repmat([ylow; yhigh],1,n);

x  = [x;nan(1,n)];
y  = [y;nan(1,n)];

h  = plot(x(:),y(:),varargin{:});
ylim(ax,yl);

hlist = addlistener(ax,'YLim','PostSet',@(src,evnt)updatedata);
h.DeleteFcn = @(src,evt)upondelete(src,evt,hlist);

if nargout == 1
    ht = h;
end

function updatedata
    yl = ylim(ax);
    ylow  = yl(1) + (yl(2)-yl(1))*rel(1);
    yhigh = yl(1) + (yl(2)-yl(1))*rel(2);
    h.YData(1:3:end) = ylow;
    h.YData(2:3:end) = yhigh;

end
    function upondelete(scr,evt,hlist)
        delete(hlist)
    end
end

