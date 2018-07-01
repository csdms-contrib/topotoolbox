function h = plotzobj(app,objdata,dointerp,distadjust)
% PLOTZOBJ draws a distance-elevation plot for simple topoapp objects.
% For STREAMobj, topoapp uses the class-specific function 'plotdz'
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013
x = objdata.x;
y = objdata.y;
d = objdata.distance;
[XX,YY] = getcoordinates(app.DEM);
z = interp2(XX,YY,app.DEM.Z,x,y);

h = plot(d+distadjust,z,'k-');

if dointerp
    set(h,'Linestyle','none','Marker','s','Markeredgecolor','r','Markersize',2)
    d = getdistance(x,y);
    [x,y,d] = interpline(x,y,d,1);
    z = interp2(XX,YY,app.DEM.Z,x,y);
    hold on, plot(d+distadjust,z,'k-'), hold off
end
xlabel('Distance along profile (m)')
ylabel('Elevation (m)')
h = gcf;
end %