function [demn,Xn,Yn] = demresize(X,Y,dem,scale)

cs = (Y(1)-Y(2));

if cs > 0
    Y = flipud(Y);
    dem = flipud(dem);
else
    cs = -cs;
end

csn = cs/scale;

minx = min(X(1,:));
maxx = max(X(1,:));
miny = min(Y(:,1));
maxy = max(Y(:,1));

F = griddedInterpolant(X',Y',dem','cubic');
[Xn,Yn] = ndgrid(minx:csn:maxx,miny:csn:maxy);

demn = F(Xn,Yn);
Xn = Xn';
Yn = Yn';
demn = demn';




