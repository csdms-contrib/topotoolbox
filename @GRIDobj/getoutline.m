function varargout = getoutline(DEM)

% get or plot extent of GRIDobj
%
% Syntax
%
%     getoutline(DEM)
%     [x,y] = getoutline(DEM)
%
% Description
%
%     getoutline plots the extent of a GRIDobj or returns the coordinate
%     vectors that generate the plot.
%
% Input arguments
%
%     DEM    GRIDobj
%
% Output arguments
%
%     x,y    coordinate vectors that can be used to plot the extent
%            rectangle (plot(x,y))
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. June, 2014


[x,y] = getcoordinates(DEM);
maxx = max(x);
minx = min(x);

maxy = max(y);
miny = min(y);

x = [minx minx maxx maxx minx];
y = [miny maxy maxy miny miny];

if nargout == 0;
    plot(x,y,'LineWidth',3,'Color',[.8 .8 .8]);
    hh = ishold;
    hold on
    plot(x,y,'k--','LineWidth',1);
    if ~hh
        hold off;
    end
else
    varargout{1} = x;
    varargout{2} = y;
end
    
       